# =====================================================================
# 01_data.R -- the `flcdata` object
#
# WHY A NEW STRUCTURE.  The whole selling point of the generalized Frechet
# distance in Kang et al. (2023) is that it handles trajectories that are
#   (a) irregularly spaced within a subject,
#   (b) asynchronous across subjects,
#   (c) asynchronous across VARIABLES within a subject,
#   (d) of different lengths,
# without imputing onto a common grid.  An [n x T x p] array cannot express
# any of that.  So `flcdata` is RAGGED by construction:
#
#   X   : list of n matrices, X[[i]] is (T_i x p); NA allowed (variable k not
#         measured at time t_i)
#   tt  : list of n numeric vectors, tt[[i]] is length T_i
#   y   : factor of length n
#
# A regular [n x T x p] array is just a special case (all tt[[i]] equal, no
# NA), and `as_flcdata()` converts either way.  Methods that genuinely need
# a rectangular grid (Euclidean, ROCKET, depth, PCA, RF, ...) call
# `flc_regularize()`, which linearly interpolates onto a common grid and
# RECORDS that it had to.  That is exactly the asymmetry the paper claims
# for SFKmL and that we now measure in a classification setting.
# =====================================================================

new_flcdata <- function(X, tt, y, varNames = NULL, name = "flcdata",
                        informative = NULL, meta = list()) {
  stopifnot(is.list(X), is.list(tt), length(X) == length(tt))
  n <- length(X)
  X <- lapply(X, function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m })
  p <- ncol(X[[1]])
  stopifnot(all(vapply(X, ncol, 1L) == p))
  stopifnot(all(vapply(seq_len(n), function(i) nrow(X[[i]]) == length(tt[[i]]), TRUE)))
  y <- factor(y)
  if (is.null(varNames)) varNames <- paste0("V", seq_len(p))
  if (is.null(informative)) informative <- rep(NA, p)

  structure(list(
    X = X, tt = lapply(tt, as.numeric), y = y,
    n = n, p = p, varNames = varNames, name = name,
    informative = informative,          # logical p-vector: TRUE = signal, FALSE = noise
    meta = meta
  ), class = "flcdata")
}

# ---- convert an [n x T x p] array (the classic layout) -----------------
array_to_flcdata <- function(traj, y, time = seq_len(dim(traj)[2]),
                             varNames = NULL, name = "flcdata",
                             informative = NULL, meta = list()) {
  n <- dim(traj)[1]; Tn <- dim(traj)[2]; p <- dim(traj)[3]
  X  <- lapply(seq_len(n), function(i) matrix(traj[i, , ], Tn, p))
  tt <- rep(list(as.numeric(time)), n)
  new_flcdata(X, tt, y, varNames, name, informative, meta)
}

# ---- back to an array, interpolating if necessary ----------------------
# `grid` NULL -> use the union of observed times (regular case: unchanged).
# Records ds$meta$interpolated = TRUE when interpolation actually happened.
flc_regularize <- function(ds, grid = NULL, rule = 2) {
  same <- length(unique(lapply(ds$tt, function(v) round(v, 10)))) == 1L
  anyNA_ <- any(vapply(ds$X, function(m) anyNA(m), TRUE))
  if (is.null(grid)) {
    grid <- if (same) ds$tt[[1]] else sort(unique(round(unlist(ds$tt), 10)))
    # cap the grid size for the union case
    if (length(grid) > 60) grid <- seq(min(grid), max(grid), length.out = 40)
  }
  Tn <- length(grid); p <- ds$p; n <- ds$n
  traj <- array(NA_real_, c(n, Tn, p))
  for (i in seq_len(n)) {
    ti <- ds$tt[[i]]; Xi <- ds$X[[i]]
    for (k in seq_len(p)) {
      ok <- is.finite(Xi[, k])
      if (sum(ok) == 0L) { traj[i, , k] <- 0; next }
      if (sum(ok) == 1L) { traj[i, , k] <- Xi[ok, k]; next }
      traj[i, , k] <- stats::approx(ti[ok], Xi[ok, k], xout = grid,
                                    rule = rule, ties = mean)$y
    }
  }
  traj[!is.finite(traj)] <- 0
  list(traj = traj, time = grid, y = ds$y,
       interpolated = (!same) || anyNA_)
}

# ---- subsetting -------------------------------------------------------
`[.flcdata` <- function(ds, i) {
  out <- ds
  out$X  <- ds$X[i]
  out$tt <- ds$tt[i]
  out$y  <- droplevels_keep(ds$y[i], levels(ds$y))
  out$n  <- length(out$X)
  out
}
droplevels_keep <- function(y, lv) factor(as.character(y), levels = lv)

print.flcdata <- function(x, ...) {
  Ts <- vapply(x$X, nrow, 1L)
  cat(sprintf("<flcdata> %s\n", x$name))
  cat(sprintf("  n = %d, p = %d, classes = %s\n", x$n, x$p,
              paste(levels(x$y), table(x$y), sep = ":", collapse = " ")))
  cat(sprintf("  time points per subject: %s\n",
              if (length(unique(Ts)) == 1) sprintf("%d (regular)", Ts[1])
              else sprintf("%d-%d (RAGGED)", min(Ts), max(Ts))))
  na <- mean(vapply(x$X, function(m) mean(is.na(m)), 0))
  cat(sprintf("  missing entries: %.1f%%\n", 100 * na))
  if (!all(is.na(x$informative)))
    cat(sprintf("  informative vars: %s | noise vars: %s\n",
                paste(x$varNames[which(x$informative)], collapse = ","),
                paste(x$varNames[which(!x$informative)], collapse = ",")))
  invisible(x)
}

# =====================================================================
# scale parameters (paper Sect. 2.2.1 & 4)
#   lambda (time-scale)  : A : (t, x) -> (lambda t, x)
#   gamma  (range-scale) : R : (t, x) -> (t, gamma_k x_k),  gamma_k = c/range_k
# BOTH must be estimated on the TRAINING fold only -- otherwise the test
# curves leak their range into the metric.
# =====================================================================

flc_gamma <- function(ds, c = 100) {
  rng <- vapply(seq_len(ds$p), function(k) {
    v <- unlist(lapply(ds$X, function(m) m[, k]))
    r <- diff(range(v, na.rm = TRUE))
    if (!is.finite(r) || r <= 0) 1 else r
  }, 0)
  c / rng
}

# paper's rule of thumb: lambda = range(variable) / range(time) * 0.1
# applied AFTER range-rescaling, all variables have range c, so this is
# just c / range(time) * 0.1.
flc_lambda_rule <- function(ds, gamma = NULL, c = 100) {
  rt <- diff(range(unlist(ds$tt)))
  if (!is.finite(rt) || rt <= 0) rt <- 1
  if (is.null(gamma)) c / rt * 0.1 else c / rt * 0.1
}

# apply the range scaling to the data itself (equivalent to carrying gamma
# in the metric, but keeps every downstream method on the same footing)
flc_rescale <- function(ds, gamma) {
  ds$X <- lapply(ds$X, function(m) sweep(m, 2, gamma, `*`))
  ds$meta$gamma <- gamma
  ds
}

# =====================================================================
# helpers that feed the C++ kernels
# =====================================================================

# joint-curve view: rows with ANY NA are dropped (a joint curve is undefined
# where a coordinate is missing).  This is precisely why MFKmL needs
# imputation and SFKmL does not.
flc_joint_list <- function(ds, idx = seq_len(ds$n)) {
  Xs <- vector("list", length(idx)); ts <- vector("list", length(idx))
  for (a in seq_along(idx)) {
    i <- idx[a]; m <- ds$X[[i]]; ok <- stats::complete.cases(m)
    if (!any(ok)) { ok[1] <- TRUE; m[1, ][is.na(m[1, ])] <- 0 }
    Xs[[a]] <- matrix(m[ok, , drop = FALSE], ncol = ds$p)
    ts[[a]] <- ds$tt[[i]][ok]
  }
  list(X = Xs, t = ts)
}

# variable-k view: only the times where variable k was actually observed.
# Different variables may have completely different time supports.
flc_var_list <- function(ds, k, idx = seq_len(ds$n)) {
  Xs <- vector("list", length(idx)); ts <- vector("list", length(idx))
  for (a in seq_along(idx)) {
    i <- idx[a]; v <- ds$X[[i]][, k]; ok <- is.finite(v)
    if (!any(ok)) { ok[1] <- TRUE; v[1] <- 0 }
    Xs[[a]] <- matrix(v[ok], ncol = 1L)
    ts[[a]] <- ds$tt[[i]][ok]
  }
  list(X = Xs, t = ts)
}

# how many curves have at least one incomplete row?  (reporting)
flc_missing_report <- function(ds) {
  inc <- vapply(ds$X, function(m) any(!stats::complete.cases(m)), TRUE)
  list(frac_curves_incomplete = mean(inc),
       frac_cells_missing = mean(vapply(ds$X, function(m) mean(is.na(m)), 0)))
}
