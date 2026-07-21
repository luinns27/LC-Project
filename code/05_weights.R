# =====================================================================
# 05_weights.R -- sparse functional-variable weights.
#
# Kang et al. (2023) eq. (5)/(8)/(9), transcribed to the SUPERVISED case:
# the clustering partition G is replaced by the class LABELS.
#
#     a_k  = BCSS of variable k under the class partition
#     w    = argmax_w  w'a   s.t.  ||w||_2 <= 1, ||w||_1 <= s, w >= 0
#          = ST(a_+, Delta) / ||ST(a_+, Delta)||_2                    (eq. 9)
#
# where Delta is DERIVED from the L1 bound s in closed form (Appendix B),
# not a free soft-threshold.  s in [1, sqrt(p)]:
#     s = sqrt(p) -> dense       s -> 1 -> maximally sparse (one variable)
#
# ---------------------------------------------------------------------
# TWO CORRECTIONS / EXTENSIONS OVER A NAIVE PORT
#
# (1) optim_w: the quadratic in Delta has TWO roots.  Only roots that are
#     SELF-CONSISTENT with their truncation level l -- i.e. exactly l entries
#     of a exceed Delta -- are admissible.  Enforcing that (plus Delta >= 0)
#     removes spurious solutions that a "take the minus root, pick the l
#     whose ||w||_1 is closest to s" loop can silently accept.
#
# (2) Choosing s.  The paper uses a GAP STATISTIC over distance-array
#     permutations.  That is the right null for CLUSTERING.  For a
#     CLASSIFIER it is not: (i) the natural null is "this variable does not
#     separate the classes", which is a LABEL permutation, and (ii) we have
#     a supervised loss, so we can just tune s by inner CV on the thing we
#     actually care about.
#
#     Empirically the gap rule is too conservative: in Setting 6, V3 has a
#     real but ~6x smaller BCSS than V1/V2 and gets shrunk to zero (variable
#     TPR 0.67), even though it helps classification.  `s_select = "cv"`
#     fixes that while keeping the noise-variable FPR at 0.
#
#     We implement all three and REPORT ALL THREE, so the comparison is on
#     the table rather than hidden in a default:
#        s_select = "gap"    (paper, distance permutation)
#        s_select = "gaplab" (label permutation -- the supervised null)
#        s_select = "cv"     (inner-CV classification loss)  <- recommended
# =====================================================================

# ---- soft threshold ---------------------------------------------------
soft_threshold <- function(x, c) sign(x) * pmax(abs(x) - c, 0)

# ---- eq. (9): weights under ||w||_2<=1, ||w||_1<=S, w>=0 ---------------
optim_w <- function(a, S, tol = 1e-9) {
  p  <- length(a)
  a  <- pmax(a, 0)                       # a_+ = max(a, 0)
  if (all(a <= tol)) return(rep(1 / sqrt(p), p))

  S  <- max(1, min(S, sqrt(p)))
  l2 <- function(v) sqrt(sum(v^2))

  # Delta = 0 is optimal iff the dense solution already satisfies ||w||_1 <= S
  w0 <- a / l2(a)
  if (sum(w0) <= S + tol) return(w0)

  sa <- sort(a, decreasing = TRUE)
  cand <- list()
  for (l in seq_len(p)) {
    S1 <- sum(sa[seq_len(l)]); S2 <- sum(sa[seq_len(l)]^2)
    den <- l - S^2
    if (abs(den) < 1e-10) next
    disc <- S1^2 / l - (S1^2 - S2 * S^2) / den
    if (!is.finite(disc) || disc < 0) next
    for (sgn in c(-1, 1)) {
      Delta <- S1 / l + sgn * sqrt(disc) / sqrt(l)
      if (!is.finite(Delta) || Delta < -tol) next
      Delta <- max(Delta, 0)
      # SELF-CONSISTENCY: exactly l entries of a must exceed Delta
      nz <- sum(a > Delta + tol)
      if (nz != l) next
      wt <- pmax(a - Delta, 0)
      nn <- l2(wt)
      if (nn <= tol) next
      w <- wt / nn
      cand[[length(cand) + 1L]] <- list(Delta = Delta, w = w, l1 = sum(w))
    }
  }
  if (!length(cand)) {
    # numerically degenerate (ties in a).  Bisect on Delta instead -- ||w||_1
    # is monotone decreasing in Delta, so this always converges.
    f <- function(D) { wt <- pmax(a - D, 0); n <- l2(wt)
                       if (n <= tol) return(1 - S) else return(sum(wt / n) - S) }
    lo <- 0; hi <- max(a)
    for (it in 1:200) { mid <- (lo + hi) / 2
      if (f(mid) > 0) lo <- mid else hi <- mid }
    wt <- pmax(a - (lo + hi) / 2, 0); n <- l2(wt)
    return(if (n <= tol) rep(1 / sqrt(p), p) else wt / n)
  }
  l1 <- vapply(cand, `[[`, 0, "l1")
  cand[[which.min(abs(l1 - S))]]$w
}

# ---- per-variable between-class sum of squares (BCSS) ------------------
# Uses the Huygens identity on a distance array (Witten & Tibshirani 2010,
# eq. 12; the paper's `bcss.feature`).  Works for ANY distance -- which is
# what lets us reuse the cached Frechet array.
#
#   BCSS_k = (1/(2n)) sum_{i,i'} d_k(i,i')  -  sum_g (1/(2 n_g)) sum_{i,i' in g} d_k(i,i')
bcss_feature <- function(Dk, assign) {
  n <- length(assign)
  tot <- sum(Dk) / (2 * n)
  wit <- 0
  for (g in unique(assign)) {
    idx <- which(assign == g)
    wit <- wit + sum(Dk[idx, idx, drop = FALSE]) / (2 * length(idx))
  }
  tot - wit
}

bcss_all <- function(Dvar, assign)
  vapply(seq_len(dim(Dvar)[3]), function(k) bcss_feature(Dvar[, , k], assign), 0)

# ---- permutation nulls ------------------------------------------------
# (paper) permute the off-diagonal of each variable's distance matrix
permute_dist_array <- function(Dvar) {
  Dp <- Dvar
  for (k in seq_len(dim(Dvar)[3])) {
    M <- Dvar[, , k]; lo <- M[lower.tri(M)]
    lo <- sample(lo); M[lower.tri(M)] <- lo
    M <- t(M); M[lower.tri(M)] <- lo
    diag(M) <- 0
    Dp[, , k] <- M
  }
  Dp
}

# ---- gap statistic over the L1 bound s --------------------------------
# perm = "dist"  : paper-faithful (permute distances)
# perm = "label" : supervised null (permute class labels)  -- also gives the
#                  per-variable permutation p-values for free
gap_select_s <- function(Dvar, y, s_grid = NULL, nperm = 20L,
                         perm = c("label", "dist"), seed = 1L) {
  perm <- match.arg(perm)
  set.seed(seed)
  y <- factor(y); p <- dim(Dvar)[3]
  if (is.null(s_grid)) s_grid <- seq(1.02, sqrt(p) * 0.99, length.out = 12)

  a_obs <- bcss_all(Dvar, as.integer(y))
  U <- vapply(s_grid, function(s) sum(optim_w(a_obs, s) * a_obs), 0)

  Ub <- matrix(NA_real_, nperm, length(s_grid))
  Ab <- matrix(NA_real_, nperm, p)                       # for p-values
  for (b in seq_len(nperm)) {
    a_b <- if (perm == "label") bcss_all(Dvar, sample(as.integer(y)))
           else                 bcss_all(permute_dist_array(Dvar), as.integer(y))
    Ab[b, ] <- a_b
    Ub[b, ] <- vapply(s_grid, function(s) sum(optim_w(a_b, s) * a_b), 0)
  }

  gap <- log(pmax(U, 1e-12)) - colMeans(log(pmax(Ub, 1e-12)))
  gap[!is.finite(gap)] <- -Inf
  best <- if (all(!is.finite(gap))) which.max(U) else which.max(gap)

  # per-variable permutation p-value (label null only)
  pval <- if (perm == "label")
    vapply(seq_len(p), function(k) (1 + sum(Ab[, k] >= a_obs[k])) / (nperm + 1), 0)
  else rep(NA_real_, p)

  list(s = s_grid[best], s_grid = s_grid, gap = gap, U = U,
       bcss = a_obs, pval = pval,
       nnz = vapply(s_grid, function(s) sum(optim_w(a_obs, s) > 1e-8), 0L),
       W = t(vapply(s_grid, function(s) optim_w(a_obs, s), numeric(p))))
}

# ---- inner-CV selection of s (RECOMMENDED for a classifier) -----------
# The loss is the actual classification loss of the SFKmL-C rule (medoid +
# weighted sum-over-variables distance), evaluated by inner stratified CV on
# the TRAINING fold only.  Distances come from the cache -> essentially free.
#
# `rule1se = TRUE` prefers the SPARSEST s within one SE of the best -- the
# classifier keeps its accuracy but drops more noise variables.
cv_select_s <- function(Dvar, y, s_grid = NULL, K = 5L, seed = 1L,
                        rule1se = TRUE, centre = c("medoid", "mean_free")) {
  centre <- match.arg(centre)
  set.seed(seed)
  y <- factor(y); n <- length(y); p <- dim(Dvar)[3]
  if (is.null(s_grid)) s_grid <- seq(1.02, sqrt(p) * 0.99, length.out = 12)

  fold <- strat_folds(y, K, seed)
  err  <- matrix(NA_real_, K, length(s_grid))

  for (kk in seq_len(K)) {
    tr <- which(fold != kk); te <- which(fold == kk)
    if (length(unique(y[tr])) < nlevels(y)) next
    a <- bcss_all(Dvar[tr, tr, , drop = FALSE], as.integer(y[tr]))
    for (si in seq_along(s_grid)) {
      w <- optim_w(a, s_grid[si])
      Dtr <- matrix(0, length(tr), length(tr))
      Dte <- matrix(0, length(te), length(tr))
      for (k in which(w > 1e-12)) {
        Dtr <- Dtr + w[k] * Dvar[tr, tr, k]
        Dte <- Dte + w[k] * Dvar[te, tr, k]
      }
      med <- vapply(levels(y), function(c) {
        idx <- which(y[tr] == c)
        if (!length(idx)) return(NA_integer_)
        idx[which.min(colSums(Dtr[idx, idx, drop = FALSE]))]
      }, 1L)
      ok <- !is.na(med)
      Dp   <- Dte[, med[ok], drop = FALSE]
      pred <- factor(levels(y)[which(ok)][max.col(-Dp, ties.method = "first")],
                     levels = levels(y))
      err[kk, si] <- mean(pred != y[te])
    }
  }
  m  <- colMeans(err, na.rm = TRUE)
  se <- apply(err, 2, function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
  se[!is.finite(se)] <- 0
  ibest <- which.min(m)
  if (rule1se) {
    thr <- m[ibest] + se[ibest]
    ok  <- which(m <= thr)
    ibest <- ok[which.min(s_grid[ok])]   # sparsest s within 1 SE
  }
  a_obs <- bcss_all(Dvar, as.integer(y))
  list(s = s_grid[ibest], s_grid = s_grid, cv_err = m, cv_se = se,
       bcss = a_obs,
       nnz = vapply(s_grid, function(s) sum(optim_w(a_obs, s) > 1e-8), 0L),
       W = t(vapply(s_grid, function(s) optim_w(a_obs, s), numeric(p))))
}

# ---- unified front door ----------------------------------------------
flc_learn_weights <- function(Dvar, y,
                              s_select = c("cv", "gaplab", "gap", "fixed", "dense"),
                              s = NULL, nperm = 20L, K = 5L, seed = 1L,
                              rule1se = TRUE) {
  s_select <- match.arg(s_select)
  p <- dim(Dvar)[3]
  a <- bcss_all(Dvar, as.integer(factor(y)))

  sel <- switch(s_select,
    dense  = list(s = sqrt(p), bcss = a),
    fixed  = list(s = s %||% sqrt(p), bcss = a),
    gap    = gap_select_s(Dvar, y, nperm = nperm, perm = "dist",  seed = seed),
    gaplab = gap_select_s(Dvar, y, nperm = nperm, perm = "label", seed = seed),
    cv     = cv_select_s (Dvar, y, K = K, seed = seed, rule1se = rule1se))

  w <- optim_w(a, sel$s)                  # ||w||_2 = 1  (paper convention)
  list(w = w,                             # report THIS
       w_scaled = if (sum(w) > 0) w / sum(w) * p else rep(1, p),  # for distances
       s = sel$s, bcss = a, sel = sel, s_select = s_select)
}
