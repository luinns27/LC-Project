# =====================================================================
# 08_base_dist.R -- distance-based competitors.
#
#  NC-Euclid      pointwise nearest centroid.  The classification analogue of
#                 KmL3d (Genolini et al. 2013), i.e. the method Kang et al.
#                 compare against.  Needs a common grid -> interpolation.
#  1NN-Euclid     the trivial floor.
#  1NN-DTW        THE benchmark of the TSC literature.  Bagnall et al. (2017)
#                 and Middlehurst, Schaefer & Bagnall (2024) both use it as
#                 the reference every new algorithm must beat.
#  1NN-DDTW       derivative DTW (Keogh & Pazzani) -- shape-sensitive variant.
#  NC-DBA         DTW Barycenter Averaging + nearest centroid.
#                 Petitjean, Forestier, Webb, Nicholson, Chen & Keogh (2014),
#                 "DTW Averaging of Time Series allows Faster and more
#                 Accurate Classification".  This is the closest existing
#                 method to MFKmL-C: a warping-aware AVERAGE used as a class
#                 prototype.  It is the single most important comparator for
#                 the MFKmL-C row.
# =====================================================================

# ---- NC-Euclid (kml3d spirit) -----------------------------------------
fit_nc_euclid <- function(ds, ctx, ...) {
  R <- flc_regularize(ds, ctx$grid)
  y <- ds$y; cls <- levels(y)
  mu <- lapply(cls, function(c) apply(R$traj[y == c, , , drop = FALSE], c(2, 3), mean))
  structure(list(mu = mu, levels = cls, grid = R$time,
                 interpolated = R$interpolated, ctx = ctx,
                 method = "NC-Euclid (kml3d)"), class = "ncEuclid")
}
predict.ncEuclid <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  cls <- object$levels
  D <- vapply(object$mu, function(m)
    apply(R$traj, 1, function(x) sqrt(sum((x - m)^2))), numeric(ds_te$n))
  D <- matrix(D, ds_te$n, length(cls), dimnames = list(NULL, cls))
  list(prob = dist_to_prob(D),
       class = factor(cls[max.col(-D, ties.method = "first")], levels = cls), D = D)
}

# ---- kNN on the Euclidean / DTW caches (uses fit_distclf) --------------
fit_knn_euclid <- function(ds, ctx, ...) fit_distclf(ds, ctx, "euclid", "knn", "none", ...)
fit_knn_dtw    <- function(ds, ctx, ...) fit_distclf(ds, ctx, "dtw",    "knn", "none", ...)
fit_1nn_dtw    <- function(ds, ctx, ...) fit_distclf(ds, ctx, "dtw",    "knn", "none", k = 1, ...)
fit_1nn_euclid <- function(ds, ctx, ...) fit_distclf(ds, ctx, "euclid", "knn", "none", k = 1, ...)
fit_knn_joint  <- function(ds, ctx, ...) fit_distclf(ds, ctx, "joint",  "knn", "none", ...)
fit_knn_sumvars<- function(ds, ctx, ...) fit_distclf(ds, ctx, "sumvars","knn", "cv",   ...)
fit_svm_sumvars<- function(ds, ctx, ...) fit_distclf(ds, ctx, "sumvars","svm", "cv",   ...)
fit_svm_joint  <- function(ds, ctx, ...) fit_distclf(ds, ctx, "joint",  "svm", "none", ...)
fit_svm_dtw    <- function(ds, ctx, ...) fit_distclf(ds, ctx, "dtw",    "svm", "none", ...)
fit_svm_euclid <- function(ds, ctx, ...) fit_distclf(ds, ctx, "euclid", "svm", "none", ...)

# ---- 1NN-DDTW : DTW on the first derivative ---------------------------
.derivative_ds <- function(ds) {
  ds$X <- lapply(seq_len(ds$n), function(i) {
    M <- ds$X[[i]]; t <- ds$tt[[i]]
    D <- M
    for (k in seq_len(ncol(M))) {
      ok <- is.finite(M[, k])
      if (sum(ok) < 3L) { D[, k] <- 0; next }
      xk <- M[ok, k]; tk <- t[ok]; m <- length(xk)
      # Keogh & Pazzani's central estimate
      d <- numeric(m)
      d[2:(m - 1)] <- ((xk[2:(m - 1)] - xk[1:(m - 2)]) +
                       (xk[3:m] - xk[1:(m - 2)]) / 2) / 2
      d[1] <- d[2]; d[m] <- d[m - 1]
      D[, k] <- NA_real_; D[ok, k] <- d
    }
    D
  })
  ds
}
fit_1nn_ddtw <- function(ds, ctx, ...) {
  dsd <- .derivative_ds(ds)
  Dtr <- flc_dist_dtw(dsd, rep(1, ds$p))
  structure(list(ds = dsd, y = ds$y, levels = levels(ds$y), k = 1L,
                 ctx = ctx, method = "1NN-DDTW"), class = "ddtwKNN")
}
predict.ddtwKNN <- function(object, ds_te, ctx = object$ctx, ...) {
  D <- flc_cross_dtw(.derivative_ds(ds_te), object$ds, rep(1, ds_te$p))
  cls <- object$levels
  nn <- apply(D, 1, which.min)
  pred <- factor(as.character(object$y[nn]), levels = cls)
  list(prob = dist_to_prob(vapply(cls, function(c)
         apply(D[, object$y == c, drop = FALSE], 1, min), numeric(nrow(D)))),
       class = pred, D = D)
}

# ---- NC-DBA : DTW Barycenter Averaging nearest centroid ----------------
# Petitjean et al. (2014).  The DTW analogue of MFKmL-C.
fit_nc_dba <- function(ds, ctx, nIter = 5L, ...) {
  y <- ds$y; cls <- levels(y)
  protos <- lapply(cls, function(c)
    flc_dba_mean(ds, which(y == c), w = rep(1, ds$p), nIter = nIter,
                 grid = ctx$grid))
  structure(list(protos = protos, levels = cls, p = ds$p,
                 varNames = ds$varNames, ctx = ctx,
                 method = "NC-DBA (Petitjean 2014)"), class = "ncDBA")
}
predict.ncDBA <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  D <- vapply(object$protos, function(mu) {
    m <- flc_curve_as_ds(mu, object$p, object$varNames)
    flc_cross_dtw(ds_te, m, rep(1, object$p))[, 1]
  }, numeric(ds_te$n))
  D <- matrix(D, ds_te$n, length(cls), dimnames = list(NULL, cls))
  list(prob = dist_to_prob(D),
       class = factor(cls[max.col(-D, ties.method = "first")], levels = cls), D = D)
}
