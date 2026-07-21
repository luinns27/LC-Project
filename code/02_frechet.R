# =====================================================================
# 02_frechet.R -- distances, Frechet mean, and the DISTANCE CACHE.
#
# The single most important performance idea in this project:
#
#   For a given dataset + (lambda, gamma, sumOrMax) we compute the FULL
#   n x n x p variable-wise Frechet array ONCE, plus the n x n joint-curve
#   array, DTW array and Euclidean array.  Every CV fold, every weight
#   vector w, every classifier then just SUBSETS and LINEARLY COMBINES the
#   cache:
#
#         d_w(i,j) = sum_k w_k * Dvar[i,j,k]              (SFKmL distance)
#
#   No distance is ever recomputed for a new w.  This turns the benchmark
#   from O(#methods x #folds) distance computations into O(1).
#
# The ONLY thing that cannot be cached is a weighted JOINT-curve Frechet
# distance (the weights enter inside the sqrt), so `frechet|<sparse>|*`
# combos pay full price -- they are marked `cacheable = FALSE`.
# =====================================================================

# ---- pairwise / cross distance matrices -------------------------------

flc_dist_joint <- function(ds, w = rep(1, ds$p), lambda = 0.1,
                           sumOrMax = "max", idx = seq_len(ds$n),
                           nthreads = flc_kernel_threads()) {
  L <- flc_joint_list(ds, idx)
  frechet_distmat_cpp(L$X, L$t, w, lambda,
                      if (sumOrMax == "sum") 1L else 0L, nthreads)
}

flc_cross_joint <- function(dsA, dsB, w = rep(1, dsA$p), lambda = 0.1,
                            sumOrMax = "max", nthreads = flc_kernel_threads()) {
  A <- flc_joint_list(dsA); B <- flc_joint_list(dsB)
  frechet_crossmat_cpp(A$X, A$t, B$X, B$t, w, lambda,
                       if (sumOrMax == "sum") 1L else 0L, nthreads)
}

# n x n x p : variable-wise Frechet distances (SFKmL's `dist.array`)
flc_dist_var <- function(ds, lambda = 0.1, sumOrMax = "max",
                         idx = seq_len(ds$n), nthreads = flc_kernel_threads()) {
  m <- length(idx); p <- ds$p
  D <- array(0, c(m, m, p), dimnames = list(NULL, NULL, ds$varNames))
  for (k in seq_len(p)) {
    L <- flc_var_list(ds, k, idx)
    D[, , k] <- frechet_distmat_cpp(L$X, L$t, 1, lambda,
                                    if (sumOrMax == "sum") 1L else 0L, nthreads)
  }
  D
}

flc_cross_var <- function(dsA, dsB, lambda = 0.1, sumOrMax = "max",
                          nthreads = flc_kernel_threads()) {
  p <- dsA$p
  D <- array(0, c(dsA$n, dsB$n, p), dimnames = list(NULL, NULL, dsA$varNames))
  for (k in seq_len(p)) {
    A <- flc_var_list(dsA, k); B <- flc_var_list(dsB, k)
    D[, , k] <- frechet_crossmat_cpp(A$X, A$t, B$X, B$t, 1, lambda,
                                     if (sumOrMax == "sum") 1L else 0L, nthreads)
  }
  D
}

flc_dist_dtw <- function(ds, w = rep(1, ds$p), idx = seq_len(ds$n),
                         nthreads = flc_kernel_threads()) {
  L <- flc_joint_list(ds, idx)
  dtw_distmat_cpp(L$X, L$t, w, nthreads)
}
flc_cross_dtw <- function(dsA, dsB, w = rep(1, dsA$p),
                          nthreads = flc_kernel_threads()) {
  A <- flc_joint_list(dsA); B <- flc_joint_list(dsB)
  dtw_crossmat_cpp(A$X, A$t, B$X, B$t, w, nthreads)
}

# Euclidean (requires a common grid -> interpolate; this is the point)
flc_dist_euclid <- function(ds, w = rep(1, ds$p)) {
  R <- flc_regularize(ds)
  Z <- R$traj
  n <- dim(Z)[1]
  F <- t(apply(Z, 1, function(m) as.vector(sweep(matrix(m, dim(Z)[2], dim(Z)[3]),
                                                 2, sqrt(w), `*`))))
  as.matrix(stats::dist(F))
}
flc_cross_euclid <- function(dsA, dsB, w = rep(1, dsA$p), grid = NULL) {
  if (is.null(grid)) grid <- sort(unique(round(unlist(c(dsA$tt, dsB$tt)), 10)))
  if (length(grid) > 60) grid <- seq(min(grid), max(grid), length.out = 40)
  RA <- flc_regularize(dsA, grid); RB <- flc_regularize(dsB, grid)
  fl <- function(Z) t(apply(Z, 1, function(m)
    as.vector(sweep(matrix(m, dim(Z)[2], dim(Z)[3]), 2, sqrt(w), `*`))))
  A <- fl(RA$traj); B <- fl(RB$traj)
  sqrt(pmax(outer(rowSums(A^2), rowSums(B^2), `+`) - 2 * A %*% t(B), 0))
}

# =====================================================================
# Frechet mean of a set of curves (kmlShape-style shape-respecting mean).
#
# Iterative: warp every member onto the current mean via the optimal
# coupling, average the warped values, repeat.  Converges in a few passes.
# Returns list(X = T x p matrix, t = time vector) -- itself an flc curve, so
# it can be fed straight back into the distance kernels.
# =====================================================================
flc_frechet_mean <- function(ds, idx = seq_len(ds$n), w = rep(1, ds$p),
                             lambda = 0.1, sumOrMax = "max", nIter = 3L,
                             grid = NULL) {
  L <- flc_joint_list(ds, idx)
  m <- length(idx)
  if (m == 1L) return(list(X = L$X[[1]], t = L$t[[1]]))

  # initialise on a common grid (pointwise mean after interpolation)
  R  <- flc_regularize(ds[idx], grid)
  mu <- apply(R$traj, c(2, 3), mean)
  tg <- R$time
  Tn <- length(tg); p <- ds$p
  som <- if (sumOrMax == "sum") 1L else 0L

  for (it in seq_len(nIter)) {
    acc <- matrix(0, Tn, p); cnt <- numeric(Tn)
    for (a in seq_len(m)) {
      pth <- frechet_path_cpp(matrix(mu, Tn, p), tg,
                              L$X[[a]], L$t[[a]], w, lambda, som)
      for (r in seq_len(nrow(pth))) {
        ii <- pth[r, 1]; jj <- pth[r, 2]
        acc[ii, ] <- acc[ii, ] + L$X[[a]][jj, ]
        cnt[ii]   <- cnt[ii] + 1
      }
    }
    mu <- acc / pmax(cnt, 1)
  }
  list(X = matrix(mu, Tn, p), t = tg)
}

# DBA -- DTW Barycenter Averaging (Petitjean, Forestier, Webb, Nicholson,
# Chen & Keogh 2014).  Same skeleton, DTW coupling, fixed-length barycentre.
flc_dba_mean <- function(ds, idx = seq_len(ds$n), w = rep(1, ds$p),
                         nIter = 5L, grid = NULL) {
  L <- flc_joint_list(ds, idx)
  m <- length(idx)
  R  <- flc_regularize(ds[idx], grid)
  mu <- apply(R$traj, c(2, 3), mean); tg <- R$time
  Tn <- length(tg); p <- ds$p
  if (m == 1L) return(list(X = L$X[[1]], t = L$t[[1]]))
  for (it in seq_len(nIter)) {
    acc <- matrix(0, Tn, p); cnt <- numeric(Tn)
    for (a in seq_len(m)) {
      pth <- dtw_path_cpp(matrix(mu, Tn, p), tg, L$X[[a]], L$t[[a]], w)
      for (r in seq_len(nrow(pth))) {
        ii <- pth[r, 1]; jj <- pth[r, 2]
        acc[ii, ] <- acc[ii, ] + L$X[[a]][jj, ]; cnt[ii] <- cnt[ii] + 1
      }
    }
    mu <- acc / pmax(cnt, 1)
  }
  list(X = matrix(mu, Tn, p), t = tg)
}

# wrap a single mean curve as a 1-subject flcdata (so cross-distances work)
flc_curve_as_ds <- function(cur, p, varNames) {
  new_flcdata(list(cur$X), list(cur$t), factor("mu"), varNames, "prototype")
}

# =====================================================================
# THE DISTANCE CACHE
# =====================================================================
# Built once per (dataset, lambda, sumOrMax).  Fold-specific views are
# obtained by plain subsetting -- no recomputation.
flc_cache <- function(ds, lambda, sumOrMax = "max",
                      want = c("var", "joint", "dtw", "euclid"),
                      nthreads = flc_kernel_threads(), verbose = FALSE) {
  t0 <- Sys.time()
  ca <- list(lambda = lambda, sumOrMax = sumOrMax, n = ds$n, p = ds$p)
  if ("var"    %in% want) ca$Dvar   <- flc_dist_var(ds, lambda, sumOrMax, nthreads = nthreads)
  if ("joint"  %in% want) ca$Djoint <- flc_dist_joint(ds, rep(1, ds$p), lambda, sumOrMax, nthreads = nthreads)
  if ("dtw"    %in% want) ca$Ddtw   <- flc_dist_dtw(ds, rep(1, ds$p), nthreads = nthreads)
  if ("euclid" %in% want) ca$Deuc   <- flc_dist_euclid(ds, rep(1, ds$p))
  ca$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) message(sprintf("  cache built in %.1fs", ca$secs))
  structure(ca, class = "flccache")
}

# weighted SFKmL distance matrix straight out of the cache -- O(n^2 p), no DP
cache_dist_w <- function(cache, w, rows = NULL, cols = NULL) {
  D <- cache$Dvar
  if (!is.null(rows) || !is.null(cols)) {
    rows <- rows %||% seq_len(dim(D)[1]); cols <- cols %||% seq_len(dim(D)[2])
    D <- D[rows, cols, , drop = FALSE]
  }
  out <- matrix(0, dim(D)[1], dim(D)[2])
  for (k in which(w > 1e-12)) out <- out + w[k] * D[, , k]
  out
}

# =====================================================================
# SFKmL 거리에 대한 올바른 barycentre  (WELL-DEFINEDNESS FIX, 2026-07)
#
# SFKmL 의 거리는
#     d(x, mu) = sum_j w_j * FD_j(x_j, mu_j)          ... (SFKmL, Eq. 6)
# 이다.  즉 변수별 1차원 Frechet 거리의 가중합이지, joint curve 위의
# 가중 Frechet 거리가 *아니다*.
#
# 어떤 거리 d 에 대한 평균(barycentre)은 정의상
#     mu = argmin_mu  sum_i d(x_i, mu)
# 이어야 한다.  위 d 를 대입하면
#     mu = argmin_mu  sum_i sum_j w_j FD_j(x_ij, mu_j)
#        = argmin_mu  sum_j w_j [ sum_i FD_j(x_ij, mu_j) ]
# 이고, 이것은 j 에 대해 **완전히 분리된다**.  따라서
#     mu_j = argmin_{mu_j}  sum_i FD_j(x_ij, mu_j)        (각 j 마다 독립)
# 이며, w_j 는 각 항의 양의 상수배일 뿐이라 argmin 을 바꾸지 못한다.
#   -> w 는 barycentre 계산에 **들어가지 않는다**.
#   -> 변수마다 1차원 Frechet mean 을 따로 구하면 된다.
#
# 구버전은 flc_frechet_mean(ds, idx, w = w_scaled) 를 불렀는데, 그것은
# joint curve 위의 가중 Frechet mean 이라 SFKmL 의 거리에 대한 barycentre
# 가 아니다.  "다른 거리로 최소화한 점"을 중심으로 쓴 셈이라
# 수학적으로 well-defined 가 아니었다.  이 함수가 그것을 바로잡는다.
#
# 반환: flc_curve_as_ds() 가 받는 형식과 동일한 list(t = ..., X = ...)
# =====================================================================
flc_frechet_mean_pervar <- function(ds, idx = seq_len(ds$n),
                                    lambda = 0.1, sumOrMax = "max",
                                    nIter = 3L, grid = NULL) {
  if (is.null(grid)) grid <- flc_common_grid(ds)
  p  <- ds$p
  Tn <- length(grid)
  X  <- matrix(NA_real_, Tn, p)

  for (k in seq_len(p)) {
    # 변수 k 만 담은 1차원 flcdata 를 만들어 그 위에서 Frechet mean 을 구한다.
    sub <- ds
    sub$X  <- lapply(idx, function(i) {
      v  <- ds$X[[i]][, k, drop = FALSE]
      ok <- is.finite(v[, 1])
      if (!any(ok)) { ok[1] <- TRUE; v[1, 1] <- 0 }
      v[ok, , drop = FALSE]
    })
    sub$tt <- lapply(idx, function(i) {
      v  <- ds$X[[i]][, k]
      ok <- is.finite(v)
      if (!any(ok)) ok[1] <- TRUE
      ds$tt[[i]][ok]
    })
    sub$n        <- length(idx)
    sub$p        <- 1L
    sub$varNames <- ds$varNames[k]
    sub$y        <- if (!is.null(ds$y)) ds$y[idx] else NULL

    # w = 1: 위 유도에 따라 가중치는 barycentre 에 영향을 주지 않는다.
    mu <- flc_frechet_mean(sub, seq_len(sub$n), w = 1,
                           lambda = lambda, sumOrMax = sumOrMax,
                           nIter = nIter, grid = grid)
    X[, k] <- stats::approx(mu$t, mu$X[, 1], xout = grid, rule = 2)$y
  }

  list(t = grid, X = X)
}
