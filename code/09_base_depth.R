# =====================================================================
# 09_base_depth.R -- the depth family.
#
#  MaxDepth-MBD    Lopez-Pintado & Romo (2006, 2009): modified band depth;
#                  classify to the class in which the curve is deepest.
#                  Multivariate extension by Ieva & Paganoni (2013):
#                      MBD(f) = sum_k p_k MBD_k(f_k)
#                  with the p_k left "problem driven".  <-- WE FILL THAT GAP
#                  by plugging in the SFKmL BCSS weights (`sparse = "cv"`).
#                  That is a small but genuine methodological contribution
#                  and it costs nothing: the weights are already computed.
#
#  MaxDepth-MFHD   Claeskens, Hubert, Slaets & Vakili (2014), multivariate
#                  functional halfspace depth (mrfDepth::mfd).  Explicitly
#                  designed to acknowledge amplitude / shape / PHASE
#                  variation, so it is the fair depth competitor in the warp
#                  settings.
#
#  DD-classifier   Li, Cuesta-Albertos & Liu (2012).  Polynomial separator
#                  through the origin of the DD-plot.  Implemented directly
#                  (order chosen by CV; multi-class by majority voting over
#                  pairwise DD-classifiers) with `ddalpha` used when present.
#
#  DistSpace-kNN   Hubert, Rousseeuw & Segaert (2017).  x -> (dist(x,P_1),
#                  ..., dist(x,P_G)) then kNN in that G-dim space.  Their
#                  distance is the bagdistance (halfspace depth based); we
#                  ALSO run the transform with our weighted Frechet distance
#                  to the class prototypes, which is a legitimate instance of
#                  their framework and a strong baseline.
# =====================================================================

# ---- modified band depth, J = 2 (closed form) --------------------------
# In-sample:  MBD(x_i) = (1/T) sum_t [ (r_it-1)(n-r_it) + (n-1) ] / C(n,2)
# Out-of-sample x vs a reference band of n curves:
#             MBD(x)   = (1/T) sum_t  r_t (n - r_t) / C(n,2),
#             r_t = #{ training curves with value <= x(t) }
.mbd_uni <- function(Xref, Xnew = NULL) {
  n <- nrow(Xref); Tn <- ncol(Xref); den <- choose(n, 2)
  if (is.null(Xnew)) {
    R <- apply(Xref, 2, rank, ties.method = "average")
    return(rowMeans(((R - 1) * (n - R) + (n - 1)) / den))
  }
  m <- nrow(Xnew)
  out <- numeric(m)
  for (a in seq_len(m)) {
    r <- vapply(seq_len(Tn), function(t) sum(Xref[, t] <= Xnew[a, t]), 0)
    out[a] <- mean(r * (n - r)) / den
  }
  out
}

# multivariate MBD (Ieva & Paganoni 2013) with weights p_k
.mbd_multi <- function(ref_traj, new_traj = NULL, pk = NULL) {
  p <- dim(ref_traj)[3]
  if (is.null(pk)) pk <- rep(1 / p, p)
  pk <- pk / sum(pk)
  m <- if (is.null(new_traj)) dim(ref_traj)[1] else dim(new_traj)[1]
  out <- numeric(m)
  for (k in seq_len(p)) {
    dk <- .mbd_uni(ref_traj[, , k],
                   if (is.null(new_traj)) NULL else new_traj[, , k])
    out <- out + pk[k] * dk
  }
  out
}

# ---- MaxDepth (MBD) ---------------------------------------------------
fit_maxdepth_mbd <- function(ds, ctx, sparse = c("uniform", "cv", "gap"), ...) {
  sparse <- match.arg(sparse)
  R <- flc_regularize(ds, ctx$grid)
  y <- ds$y; cls <- levels(y)
  pk <- rep(1 / ds$p, ds$p)
  if (sparse != "uniform") {
    Dvar_tr <- ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE]
    W <- flc_learn_weights(Dvar_tr, y,
                           s_select = if (sparse == "cv") "cv" else "gap",
                           seed = ctx$seed)
    pk <- if (sum(W$w) > 0) W$w / sum(W$w) else pk
  }
  structure(list(bands = lapply(cls, function(c) R$traj[y == c, , , drop = FALSE]),
                 levels = cls, pk = pk, grid = R$time, ctx = ctx,
                 method = sprintf("MaxDepth-MBD(%s)", sparse)),
            class = "maxDepthMBD")
}
predict.maxDepthMBD <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  cls <- object$levels
  Dep <- vapply(object$bands, function(b)
    .mbd_multi(b, R$traj, object$pk), numeric(ds_te$n))
  Dep <- matrix(Dep, ds_te$n, length(cls), dimnames = list(NULL, cls))
  P <- pmax(Dep, 0); rs <- rowSums(P); rs[rs <= 0] <- 1; P <- P / rs
  list(prob = P,
       class = factor(cls[max.col(Dep, ties.method = "first")], levels = cls),
       D = -Dep)
}

# ---- MaxDepth (MFHD, Claeskens et al. 2014) ---------------------------
fit_maxdepth_mfhd <- function(ds, ctx, ...) {
  if (!flc_have("mrfDepth")) stop("mrfDepth not installed")
  R <- flc_regularize(ds, ctx$grid)
  y <- ds$y; cls <- levels(y)
  # mrfDepth wants an array [T x n x p]
  toM <- function(A) aperm(A, c(2, 1, 3))
  structure(list(bands = lapply(cls, function(c) toM(R$traj[y == c, , , drop = FALSE])),
                 levels = cls, grid = R$time, ctx = ctx,
                 method = "MaxDepth-MFHD"), class = "maxDepthMFHD")
}
predict.maxDepthMFHD <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  Z <- aperm(R$traj, c(2, 1, 3))
  cls <- object$levels
  Dep <- vapply(object$bands, function(b) {
    r <- try(mrfDepth::mfd(x = b, z = Z, type = "hdepth", diagnostic = FALSE),
             silent = TRUE)
    if (inherits(r, "try-error")) rep(0, dim(Z)[2]) else as.numeric(r$MFDdepthZ)
  }, numeric(ds_te$n))
  Dep <- matrix(Dep, ds_te$n, length(cls), dimnames = list(NULL, cls))
  P <- pmax(Dep, 0); rs <- rowSums(P); rs[rs <= 0] <- 1; P <- P / rs
  list(prob = P,
       class = factor(cls[max.col(Dep, ties.method = "first")], levels = cls),
       D = -Dep)
}

# ---- DD-classifier (Li, Cuesta-Albertos & Liu 2012) -------------------
# Binary core: on the DD-plot (d1, d2) find a polynomial through the origin
#   d2 = a1 d1 + a2 d1^2 + ... + ak d1^k
# minimising misclassification.  We optimise a smooth logistic surrogate from
# many random starts (the paper's own strategy), pick the order by CV, and
# extend to K > 2 by majority voting over all pairs.
.ddplot_depths <- function(bandA, bandB, newtraj, pk) {
  cbind(.mbd_multi(bandA, newtraj, pk), .mbd_multi(bandB, newtraj, pk))
}
.dd_fit_binary <- function(d1, d2, lab, orders = 1:3, ninit = 100L, t_logit = 200) {
  best <- NULL
  for (K in orders) {
    Xp <- outer(d1, seq_len(K), `^`)
    obj <- function(a) {
      z <- t_logit * (Xp %*% a - d2)
      mean(ifelse(lab == 1, log1p(exp(-z)), log1p(exp(z))))
    }
    for (b in seq_len(ninit)) {
      a0 <- stats::rnorm(K, 0, 1)
      r <- try(stats::optim(a0, obj, method = "BFGS",
                            control = list(maxit = 200)), silent = TRUE)
      if (inherits(r, "try-error")) next
      pr <- as.numeric((Xp %*% r$par) > d2)   # above the curve -> class 1
      er <- mean(pr != (lab == 1))
      if (is.null(best) || er < best$err) best <- list(a = r$par, K = K, err = er)
    }
  }
  if (is.null(best)) best <- list(a = 1, K = 1L, err = 0.5)
  best
}
fit_ddclf <- function(ds, ctx, orders = 1:3, ninit = 60L, ...) {
  R <- flc_regularize(ds, ctx$grid)
  y <- ds$y; cls <- levels(y)
  pk <- rep(1 / ds$p, ds$p)
  bands <- lapply(cls, function(c) R$traj[y == c, , , drop = FALSE])
  pairs <- utils::combn(length(cls), 2, simplify = FALSE)
  fits <- lapply(pairs, function(pr) {
    i <- pr[1]; j <- pr[2]
    sel <- which(y %in% cls[c(i, j)])
    dd <- .ddplot_depths(bands[[i]], bands[[j]], R$traj[sel, , , drop = FALSE], pk)
    lab <- ifelse(y[sel] == cls[i], 1L, 0L)
    .dd_fit_binary(dd[, 1], dd[, 2], lab, orders, ninit)
  })
  structure(list(bands = bands, fits = fits, pairs = pairs, levels = cls,
                 pk = pk, grid = R$time, ctx = ctx,
                 method = "DD-classifier (Li et al. 2012)"), class = "ddClf")
}
predict.ddClf <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  cls <- object$levels; n <- ds_te$n
  V <- matrix(0, n, length(cls), dimnames = list(NULL, cls))
  for (q in seq_along(object$pairs)) {
    i <- object$pairs[[q]][1]; j <- object$pairs[[q]][2]
    dd <- .ddplot_depths(object$bands[[i]], object$bands[[j]], R$traj, object$pk)
    a  <- object$fits[[q]]$a; K <- object$fits[[q]]$K
    up <- as.numeric(outer(dd[, 1], seq_len(K), `^`) %*% a) > dd[, 2]
    V[, i] <- V[, i] + as.numeric(up)
    V[, j] <- V[, j] + as.numeric(!up)
  }
  P <- V / rowSums(V)
  list(prob = P,
       class = factor(cls[max.col(V, ties.method = "first")], levels = cls),
       D = -V)
}

# ---- DistSpace (Hubert, Rousseeuw & Segaert 2017) ---------------------
# variant = "bagdist" : their bagdistance (needs mrfDepth), on FPC scores
# variant = "frechet" : distance to each class's Frechet prototype -- an
#                       instance of their framework using our metric
fit_distspace <- function(ds, ctx, variant = c("frechet", "bagdist"),
                          sparse = c("cv", "none"), k = NULL, ...) {
  variant <- match.arg(variant); sparse <- match.arg(sparse)
  y <- ds$y; cls <- levels(y)

  if (variant == "frechet") {
    w <- rep(1, ds$p)
    if (sparse == "cv") {
      W <- flc_learn_weights(ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE], y,
                             s_select = "cv", seed = ctx$seed)
      w <- W$w
    }
    Dw  <- cache_dist_w(list(Dvar = ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE]), w)
    med <- vapply(cls, function(c) {
      idx <- which(y == c); idx[which.min(colSums(Dw[idx, idx, drop = FALSE]))] }, 1L)
    Ztr <- Dw[, med, drop = FALSE]
    obj <- list(kind = "frechet", w = w, med = med, Ztr = Ztr)
  } else {
    if (!flc_have("mrfDepth")) stop("mrfDepth not installed")
    R  <- flc_regularize(ds, ctx$grid)
    Fm <- .fpc_scores_fit(R$traj, ncomp = min(6L, dim(R$traj)[2] - 1L))
    Ztr <- vapply(cls, function(c) {
      Xc <- Fm$scores[y == c, , drop = FALSE]
      r <- try(mrfDepth::bagdistance(x = Xc, z = Fm$scores)$bagdistance,
               silent = TRUE)
      if (inherits(r, "try-error")) rep(NA_real_, nrow(Fm$scores)) else as.numeric(r)
    }, numeric(ds$n))
    Ztr[!is.finite(Ztr)] <- max(Ztr[is.finite(Ztr)], 1)
    obj <- list(kind = "bagdist", fpc = Fm, grid = R$time,
                bands = lapply(cls, function(c) Fm$scores[y == c, , drop = FALSE]),
                Ztr = Ztr)
  }

  if (is.null(k)) {
    kg <- c(1, 3, 5, 7); kg <- kg[kg < min(table(y))]; if (!length(kg)) kg <- 1
    Dl <- as.matrix(stats::dist(obj$Ztr)); diag(Dl) <- Inf
    ord <- t(apply(Dl, 1, order))
    a <- vapply(kg, function(kk) mean(apply(ord[, seq_len(kk), drop = FALSE], 1,
      function(ix) names(which.max(table(factor(y[ix], levels = cls))))) == as.character(y)), 0)
    k <- kg[which.max(a)]
  }
  obj$k <- k; obj$y <- y; obj$levels <- cls; obj$ctx <- ctx
  obj$method <- sprintf("DistSpace-%s-kNN", variant)
  structure(obj, class = "distSpace")
}
predict.distSpace <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  if (object$kind == "frechet") {
    Dc <- ctx$cache$Dvar[ctx$te, ctx$tr, , drop = FALSE]
    Dw <- matrix(0, length(ctx$te), length(ctx$tr))
    for (k in which(object$w > 1e-12)) Dw <- Dw + object$w[k] * Dc[, , k]
    Zte <- Dw[, object$med, drop = FALSE]
  } else {
    R <- flc_regularize(ds_te, object$grid)
    S <- .fpc_scores_apply(object$fpc, R$traj)
    Zte <- vapply(object$bands, function(b) {
      r <- try(mrfDepth::bagdistance(x = b, z = S)$bagdistance, silent = TRUE)
      if (inherits(r, "try-error")) rep(NA_real_, nrow(S)) else as.numeric(r)
    }, numeric(ds_te$n))
    Zte[!is.finite(Zte)] <- max(object$Ztr, na.rm = TRUE)
  }
  D <- sqrt(pmax(outer(rowSums(Zte^2), rowSums(object$Ztr^2), `+`) -
                 2 * Zte %*% t(object$Ztr), 0))
  # BUG FIX: knn_vote() 사용 (01b_helpers.R 주석 참조).
  P <- knn_vote(D, object$y, object$k, cls)
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls), D = D)
}
