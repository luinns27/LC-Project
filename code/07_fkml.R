# =====================================================================
# 07_fkml.R -- the two classifiers this project exists for.
#
# What MFKmL and SFKmL ACTUALLY are (from the fkml source, mfkml.R/sfkml.R,
# not from the abstract):
#
#   MFKmL  = k-MEANS.
#            centre   : weighted FRECHET MEAN (tournament of pairwise means)
#            distance : generalized Frechet on the JOINT curve, i.e. local
#                       cost is the euclidean norm in (lambda*t, gamma_1 x_1,
#                       ..., gamma_p x_p) space
#            weights  : NONE.  All variables enter the joint curve equally
#                       (up to the user-supplied range scaling gamma).
#
#   SFKmL  = k-MEDOIDS.
#            centre   : MEDOID -- an actual observed trajectory
#            distance : SUM OVER VARIABLES of the univariate Frechet
#                       distances, weighted:   d(x,y) = sum_k w_k d_k(x_k,y_k)
#                       *** This is NOT a joint-curve Frechet distance. ***
#            weights  : w from BCSS + optim_w under ||w||_1 <= s
#            missing  : handled natively -- each d_k only uses the times where
#                       variable k was observed
#
# SUPERVISED TRANSCRIPTION.  In both, the k-means / k-medoids ASSIGNMENT step
# is replaced by the class LABEL.  Centre definition, distance definition and
# weight estimation are kept exactly as in the package.  Nothing else changes.
#
#   MFKmL-C : per-class Frechet-mean prototype + joint-curve Frechet distance
#   SFKmL-C : per-class medoid            + weighted sum-over-variables distance
#
# ABLATIONS (to show which design choice does the work):
#   MFKmL-C-medoid   : joint distance, medoid centre       -> mean vs medoid
#   SFKmL-C-mean     : sumvars distance, Frechet-mean centre
#   SFKmL-C-dense    : sumvars distance, medoid, s = sqrt(p)  -> value of sparsity
#   SFKmL-C-kNN/SVM  : same distance, different decision rule
# =====================================================================

# ---------------------------------------------------------------------
# every fit_* here takes (ds_train, cache_train, ...) and returns an object
# with a predict() that takes (ds_test, cross) where `cross` carries the
# test->train distances (also cached).  This keeps ALL distance computation
# in one place -- see 08_engine.R.
# ---------------------------------------------------------------------

# =====================================================================
# MFKmL-C
# =====================================================================
fit_MFKmLC <- function(ds, ctx, ...) {
  y <- ds$y; cls <- levels(y)
  protos <- lapply(cls, function(c)
    flc_frechet_mean(ds, which(y == c), w = rep(1, ds$p),
                     lambda = ctx$lambda, sumOrMax = ctx$sumOrMax))
  structure(list(protos = protos, levels = cls, w = rep(1, ds$p),
                 ctx = ctx, p = ds$p, varNames = ds$varNames,
                 method = "MFKmL-C"),
            class = "MFKmLC")
}

predict.MFKmLC <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  D <- matrix(NA_real_, ds_te$n, length(cls), dimnames = list(NULL, cls))
  for (ci in seq_along(cls)) {
    mu <- flc_curve_as_ds(object$protos[[ci]], object$p, object$varNames)
    D[, ci] <- flc_cross_joint(ds_te, mu, w = object$w,
                               lambda = ctx$lambda, sumOrMax = ctx$sumOrMax)[, 1]
  }
  list(prob = dist_to_prob(D),
       class = factor(cls[max.col(-D, ties.method = "first")], levels = cls),
       D = D)
}

# MFKmL-C with a medoid centre (isolates mean-vs-medoid with NO weighting)
fit_MFKmLC_medoid <- function(ds, ctx, ...) {
  y <- ds$y; cls <- levels(y)
  Dj <- ctx$cache$Djoint[ctx$tr, ctx$tr, drop = FALSE]
  med <- vapply(cls, function(c) {
    idx <- which(y == c)
    if (length(idx) <= 2L) idx[1]
    else idx[which.min(colSums(Dj[idx, idx, drop = FALSE]))]
  }, 1L)
  structure(list(medoid_local = med, levels = cls, w = rep(1, ds$p),
                 ctx = ctx, method = "MFKmL-C-medoid"),
            class = "MFKmLCmedoid")
}

predict.MFKmLCmedoid <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  # cross joint distances (test rows x train cols), straight from cache
  D <- ctx$cache$Djoint[ctx$te, ctx$tr, drop = FALSE][, object$medoid_local, drop = FALSE]
  colnames(D) <- cls
  list(prob = dist_to_prob(D),
       class = factor(cls[max.col(-D, ties.method = "first")], levels = cls),
       D = D)
}

# =====================================================================
# SFKmL-C   (the main contribution)
# =====================================================================
fit_SFKmLC <- function(ds, ctx,
                       s_select = c("cv", "gaplab", "gap", "dense", "fixed"),
                       s = NULL, nperm = 20L, innerK = 5L, rule1se = TRUE,
                       centre = c("medoid", "mean"), ...) {
  s_select <- match.arg(s_select); centre <- match.arg(centre)
  y <- ds$y; cls <- levels(y)
  Dvar_tr <- ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE]

  W <- flc_learn_weights(Dvar_tr, y, s_select = s_select, s = s,
                         nperm = nperm, K = innerK, seed = ctx$seed,
                         rule1se = rule1se)

  obj <- list(w = W$w, w_scaled = W$w_scaled, s = W$s, bcss = W$bcss,
              sel = W$sel, levels = cls, ctx = ctx, centre = centre,
              p = ds$p, varNames = ds$varNames, s_select = s_select,
              method = sprintf("SFKmL-C(%s%s)", centre,
                               if (s_select == "dense") ",dense" else ""))

  if (centre == "medoid") {
    # medoid_k = argmin_{i in class k} sum_{i' in class k} d_w(i, i')
    # 대각(자기 자신, 거리 0)은 모든 후보에 동일하게 0 을 더하므로 argmin 을
    # 바꾸지 않는다.  그래도 정의에 맞게 명시적으로 제외한다.
    Dw <- cache_dist_w(list(Dvar = Dvar_tr), W$w)     # n_tr x n_tr
    obj$medoid_local <- vapply(cls, function(c) {
      idx <- which(y == c)
      if (length(idx) <= 2L) return(idx[1])
      sub <- Dw[idx, idx, drop = FALSE]
      diag(sub) <- 0
      idx[which.min(colSums(sub))]
    }, 1L)

  } else {
    # WELL-DEFINEDNESS FIX (2026-07):
    #   구버전은 flc_frechet_mean(ds, ..., w = W$w_scaled) 를 썼다.  그것은
    #   *joint curve* 위의 가중 Frechet mean 이라, SFKmL 의 거리
    #       d(x,mu) = sum_j w_j FD_j(x_j, mu_j)
    #   에 대한 barycentre 가 아니다.  즉 "다른 거리로 최소화한 점"을 중심으로
    #   쓴 셈이었다.
    #
    #   올바른 barycentre 는 j 에 대해 분리되며 (02_frechet.R 유도 참조),
    #   변수별 1차원 Frechet mean 이고 w 에 의존하지 않는다.
    obj$protos <- lapply(cls, function(c)
      flc_frechet_mean_pervar(ds, which(y == c),
                              lambda = ctx$lambda, sumOrMax = ctx$sumOrMax,
                              grid = ctx$grid))
  }
  structure(obj, class = "SFKmLC")
}

predict.SFKmLC <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  if (object$centre == "medoid") {
    # cached: test x train x p  ->  weighted sum  ->  columns = medoids
    Dcross <- ctx$cache$Dvar[ctx$te, ctx$tr, , drop = FALSE]
    Dw <- matrix(0, length(ctx$te), length(ctx$tr))
    for (k in which(object$w > 1e-12)) Dw <- Dw + object$w[k] * Dcross[, , k]
    D <- Dw[, object$medoid_local, drop = FALSE]
  } else {
    # 중심이 mean 이든 medoid 든 **거리는 항상 SFKmL 의 거리**를 쓴다:
    #     d(x, mu) = sum_j w_j * FD_j(x_j, mu_j)
    # (fit 과 predict 가 같은 거리를 써야 well-defined 하다.)
    D <- matrix(NA_real_, ds_te$n, length(cls))
    for (ci in seq_along(cls)) {
      mu  <- flc_curve_as_ds(object$protos[[ci]], object$p, object$varNames)
      acc <- numeric(ds_te$n)
      for (k in which(object$w > 1e-12)) {
        A <- flc_var_list(ds_te, k)
        B <- flc_var_list(mu, k)
        acc <- acc + object$w[k] *
          frechet_crossmat_cpp(A$X, A$t, B$X, B$t, 1, ctx$lambda,
                               if (ctx$sumOrMax == "sum") 1L else 0L,
                               flc_kernel_threads())[, 1]
      }
      D[, ci] <- acc
    }
  }
  colnames(D) <- cls
  list(prob = dist_to_prob(D),
       class = factor(cls[max.col(-D, ties.method = "first")], levels = cls),
       D = D)
}

# ---- convenience constructors used by the benchmark registry ----------
fit_SFKmLC_cv    <- function(ds, ctx, ...) fit_SFKmLC(ds, ctx, s_select = "cv",     ...)
fit_SFKmLC_gap   <- function(ds, ctx, ...) fit_SFKmLC(ds, ctx, s_select = "gap",    ...)
fit_SFKmLC_gaplab<- function(ds, ctx, ...) fit_SFKmLC(ds, ctx, s_select = "gaplab", ...)
fit_SFKmLC_dense <- function(ds, ctx, ...) fit_SFKmLC(ds, ctx, s_select = "dense",  ...)
fit_SFKmLC_mean  <- function(ds, ctx, ...) fit_SFKmLC(ds, ctx, s_select = "cv",
                                                      centre = "mean", ...)

# =====================================================================
# Same distance, different decision rule.
# The point of these is NOT to beat SFKmL-C but to separate
#   "the Frechet distance is doing the work"  from
#   "the medoid rule is doing the work".
# =====================================================================

# generic distance-based classifier over ANY cached distance
#   metric: "sumvars" (SFKmL) | "joint" (MFKmL) | "dtw" | "euclid"
#   clf   : "knn" | "svm" | "centroid" (centroid only for joint/euclid)
fit_distclf <- function(ds, ctx,
                        metric = c("sumvars", "joint", "dtw", "euclid"),
                        clf    = c("knn", "svm"),
                        sparse = c("none", "cv", "gap"),
                        k = NULL, C = 1, ...) {
  metric <- match.arg(metric); clf <- match.arg(clf); sparse <- match.arg(sparse)
  y <- ds$y; cls <- levels(y)

  w <- rep(1, ds$p)
  if (sparse != "none") {
    Dvar_tr <- ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE]
    W <- flc_learn_weights(Dvar_tr, y,
                           s_select = if (sparse == "cv") "cv" else "gap",
                           seed = ctx$seed)
    w <- W$w
  }

  Dtr <- .train_dist(ctx, metric, w)

  obj <- list(levels = cls, y = y, w = w, metric = metric, clf = clf,
              sparse = sparse, ctx = ctx,
              method = sprintf("%s|%s|%s", metric, sparse, clf))

  if (clf == "knn") {
    # tune k by inner LOOCV on the training distance matrix (cheap: no
    # distance recomputation, just a re-sort)
    if (is.null(k)) {
      kgrid <- c(1, 3, 5, 7, 9, 11)
      kgrid <- kgrid[kgrid < min(table(y))]
      if (!length(kgrid)) kgrid <- 1
      Dl <- Dtr; diag(Dl) <- Inf
      ord <- t(apply(Dl, 1, order))
      accs <- vapply(kgrid, function(kk) {
        nn <- ord[, seq_len(kk), drop = FALSE]
        vote <- apply(nn, 1, function(ix) {
          tb <- table(factor(y[ix], levels = cls)); names(tb)[which.max(tb)] })
        mean(vote == as.character(y))
      }, 0)
      k <- kgrid[which.max(accs)]
    }
    obj$k <- k
  } else if (clf == "svm") {
    if (!flc_have("kernlab")) stop("kernlab not installed")
    med <- stats::median(Dtr[Dtr > 0])
    if (!is.finite(med) || med <= 0) med <- 1
    gamma <- 1 / (2 * med^2)
    K <- exp(-gamma * Dtr^2)
    obj$gamma <- gamma
    obj$svm <- kernlab::ksvm(kernlab::as.kernelMatrix(K), y, type = "C-svc",
                             C = C, prob.model = TRUE)
    obj$sv <- kernlab::SVindex(obj$svm)
  }
  structure(obj, class = "flcDistClf")
}

predict.flcDistClf <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  D <- .cross_dist(ctx, object$metric, object$w)     # n_te x n_tr

  if (object$clf == "knn") {
    # BUG FIX: 예전에는 apply/t()/matrix(byrow=TRUE) 조합으로 확률행렬을
    # 만들었는데, k=1 일 때 행과 열이 뒤섞여 1NN 계열이 chance 수준이 나왔다.
    # knn_vote() 는 명시적으로 행을 채우므로 k 에 상관없이 안전하다.
    P <- knn_vote(D, object$y, object$k, cls)
    pred <- factor(cls[max.col(P, ties.method = "first")], levels = cls)
    return(list(prob = P, class = pred, D = D))
  }

  Kte <- exp(-object$gamma * D^2)[, object$sv, drop = FALSE]
  P <- kernlab::predict(object$svm, kernlab::as.kernelMatrix(Kte), type = "probabilities")
  P <- as.matrix(P)[, cls, drop = FALSE]
  pred <- factor(cls[max.col(P, ties.method = "first")], levels = cls)
  list(prob = P, class = pred, D = D)
}

# ---- cache accessors --------------------------------------------------
.train_dist <- function(ctx, metric, w) {
  switch(metric,
    sumvars = cache_dist_w(list(Dvar = ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE]), w),
    joint   = ctx$cache$Djoint[ctx$tr, ctx$tr, drop = FALSE],
    dtw     = ctx$cache$Ddtw  [ctx$tr, ctx$tr, drop = FALSE],
    euclid  = ctx$cache$Deuc  [ctx$tr, ctx$tr, drop = FALSE])
}
.cross_dist <- function(ctx, metric, w) {
  switch(metric,
    sumvars = {
      Dc <- ctx$cache$Dvar[ctx$te, ctx$tr, , drop = FALSE]
      out <- matrix(0, length(ctx$te), length(ctx$tr))
      for (k in which(w > 1e-12)) out <- out + w[k] * Dc[, , k]
      out
    },
    joint   = ctx$cache$Djoint[ctx$te, ctx$tr, drop = FALSE],
    dtw     = ctx$cache$Ddtw  [ctx$te, ctx$tr, drop = FALSE],
    euclid  = ctx$cache$Deuc  [ctx$te, ctx$tr, drop = FALSE])
}
