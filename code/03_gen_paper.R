# =====================================================================
# 03_gen_paper.R
# Kang, Choi, Yoon, Park, Kwon & Park (2023), Stat. Comput. 33:75,
# Appendix C -- Settings 1, 2, 2a, 2b, 3, 4, 5, 6, reproduced exactly.
#
#   n = 90 (3 groups x 30), K = 3.
#   phi_X(t) = N(X, 1) density,  Phi_X(t) = N(X, 1) cdf.
#   Sample trajectories randomise the peak location X ~ N(U(.), 1) and the
#   amplitude Z ~ N(.,.).
#
# NOTE for the paper we are writing.  In ALL of these settings the classes
# are separated by AMPLITUDE (the Z's differ across groups) while the shape
# family is shared.  That is a regime where a plain pointwise / summary-
# statistic method is already close to optimal -- which is why `tsfeat-RF`
# hits 1.000 here.  These settings are therefore the FIDELITY baseline
# (does our classifier reproduce the clustering paper's conditions?), NOT
# the evidence for Frechet.  That evidence lives in 04_gen_new.R.
# =====================================================================

.phi <- function(t, X) stats::dnorm(t, mean = X, sd = 1)
.Phi <- function(t, X) stats::pnorm(t, mean = X, sd = 1)

# ---------------- Setting 1 -------------------------------------------
gen_setting1 <- function(nEach = 30, seed = 1) {
  set.seed(seed); t <- 1:20; Tn <- length(t); n <- nEach * 3
  traj <- array(NA_real_, c(n, Tn, 2)); y <- factor(rep(1:3, each = nEach))
  r <- 0L
  for (g in 1:3) for (s in seq_len(nEach)) {
    r <- r + 1L
    Z1 <- stats::rnorm(Tn, stats::runif(1, -1, 1), 2)
    traj[r, , 1] <- switch(g,
      (1 / 19) * (16 * t + 22) + Z1,
      rep(10, Tn) + Z1,
      -(1 / 19) * (10 * t - 295) + Z1)
    if (g == 1L) {
      X  <- stats::rnorm(1, stats::runif(1, 2, 8), 1)
      Z1a <- stats::rnorm(1, 5, 3); Z2 <- stats::rnorm(1, 10, 2)
      traj[r, , 2] <- Z1a * .phi(t, X) + Z2
    } else if (g == 2L) {
      X1 <- stats::rnorm(1, stats::runif(1, 2, 8), 1)
      X2 <- stats::rnorm(1, stats::runif(1, 7, 13), 1)
      X3 <- stats::rnorm(1, stats::runif(1, 12, 18), 1)
      Z  <- stats::rnorm(1, 30, 3)
      traj[r, , 2] <- Z * (.phi(t, X1) + .phi(t, X2) + .phi(t, X3))
    } else {
      X1 <- stats::rnorm(1, stats::runif(1, 2, 8), 1)
      X2 <- stats::rnorm(1, stats::runif(1, 12, 18), 1)
      traj[r, , 2] <- stats::rnorm(1, 45, 3) * .phi(t, X1) +
                      stats::rnorm(1, 35, 3) * .phi(t, X2)
    }
  }
  array_to_flcdata(traj, y, t, c("V1", "V2"), "Setting1",
                   informative = c(TRUE, TRUE))
}

# ---------------- Setting 2 / 2a / 2b ---------------------------------
gen_setting2 <- function(nEach = 30, seed = 1, variant = c("2", "2a", "2b")) {
  variant <- match.arg(variant)
  set.seed(seed); t <- 1:10; Tn <- length(t); n <- nEach * 3
  traj <- array(NA_real_, c(n, Tn, 2)); y <- factor(rep(1:3, each = nEach))
  v1mean <- c(5, 10, 15);  v1sd <- 2
  v2mean <- c(15, 30, 45); v2sd <- 2
  if (variant == "2a") { v1mean <- c(10, 20, 30);  v1sd <- 6  }
  if (variant == "2b") { v2mean <- c(40, 80, 120); v2sd <- 10 }
  r <- 0L
  for (g in 1:3) for (s in seq_len(nEach)) {
    r <- r + 1L
    Xv1 <- stats::rnorm(1, stats::runif(1, 3, 7), 1)
    traj[r, , 1] <- (1 - .Phi(t, Xv1)) * stats::rnorm(1, v1mean[g], v1sd)
    Xv2 <- stats::rnorm(1, stats::runif(1, 2, 8), 1)
    traj[r, , 2] <- .phi(t, Xv2) * stats::rnorm(1, v2mean[g], v2sd)
  }
  array_to_flcdata(traj, y, t, c("V1", "V2"), paste0("Setting", variant),
                   informative = c(TRUE, TRUE))
}

# ---------------- Setting 3 : Setting 2 + Var3 ------------------------
gen_setting3 <- function(nEach = 30, seed = 1) {
  base <- gen_setting2(nEach, seed, "2")
  R <- flc_regularize(base); t <- R$time; Tn <- length(t); n <- base$n
  set.seed(seed + 1000L)
  v3 <- matrix(NA_real_, n, Tn); r <- 0L
  for (g in 1:3) for (s in seq_len(nEach)) {
    r <- r + 1L
    if (g == 1L) {
      X <- stats::rnorm(1, stats::runif(1, 3, 5), 1); Z <- stats::rnorm(1, 35, 2)
    } else if (g == 2L) {
      X <- stats::rnorm(1, stats::runif(1, 6, 8), 1); Z <- stats::rnorm(1, 45, 2)
    } else {
      U <- stats::runif(1, -1, 1)
      X <- if (stats::runif(1) < 0.4) stats::rnorm(1, 3 + U, 1) else stats::rnorm(1, 7 + U, 1)
      Z <- stats::rnorm(1, 40, 2)
    }
    v3[r, ] <- Z * .phi(t, X)
  }
  traj <- array(NA_real_, c(n, Tn, 3))
  traj[, , 1:2] <- R$traj; traj[, , 3] <- v3
  array_to_flcdata(traj, base$y, t, c("V1", "V2", "V3"), "Setting3",
                   informative = c(TRUE, TRUE, TRUE))
}

# ---------------- Setting 4 : two groups share a model (hardest) ------
gen_setting4 <- function(nEach = 30, seed = 1) {
  set.seed(seed); t <- 1:10; Tn <- length(t); n <- nEach * 3
  traj <- array(NA_real_, c(n, Tn, 2)); y <- factor(rep(1:3, each = nEach))
  v1mean <- c(10, 15, 15); v2mean <- c(45, 45, 35); r <- 0L
  for (g in 1:3) for (s in seq_len(nEach)) {
    r <- r + 1L
    Xv1 <- stats::rnorm(1, stats::runif(1, 3, 7), 1)
    traj[r, , 1] <- (1 - .Phi(t, Xv1)) * stats::rnorm(1, v1mean[g], 1)
    Xv2 <- stats::rnorm(1, stats::runif(1, 2, 8), 1)
    traj[r, , 2] <- .phi(t, Xv2) * stats::rnorm(1, v2mean[g], 1)
  }
  array_to_flcdata(traj, y, t, c("V1", "V2"), "Setting4",
                   informative = c(TRUE, TRUE))
}

# ---------------- Setting 5 / 6 : + pure-noise variables ---------------
# noise ~ N(7.5, 2.8^2), so its range is ~15, matching the signal variables
.add_noise_vars <- function(ds, nNoise, seed, name) {
  R <- flc_regularize(ds); t <- R$time; Tn <- length(t); n <- ds$n
  set.seed(seed)
  d <- ds$p + nNoise
  traj <- array(NA_real_, c(n, Tn, d))
  traj[, , seq_len(ds$p)] <- R$traj
  for (k in seq_len(nNoise))
    traj[, , ds$p + k] <- matrix(stats::rnorm(n * Tn, 7.5, 2.8), n, Tn)
  array_to_flcdata(traj, ds$y, t,
                   c(ds$varNames, paste0("N", seq_len(nNoise))), name,
                   informative = c(rep(TRUE, ds$p), rep(FALSE, nNoise)))
}

gen_setting5 <- function(nEach = 30, seed = 1, nNoise = 3)
  .add_noise_vars(gen_setting2(nEach, seed, "2"), nNoise, seed + 2000L, "Setting5")

gen_setting6 <- function(nEach = 30, seed = 1, nNoise = 7)
  .add_noise_vars(gen_setting3(nEach, seed), nNoise, seed + 3000L, "Setting6")

# ---------------- registry --------------------------------------------
FLC_PAPER_SETTINGS <- list(
  Setting1  = function(seed, ...) gen_setting1(seed = seed),
  Setting2  = function(seed, ...) gen_setting2(seed = seed, variant = "2"),
  Setting2a = function(seed, ...) gen_setting2(seed = seed, variant = "2a"),
  Setting2b = function(seed, ...) gen_setting2(seed = seed, variant = "2b"),
  Setting3  = function(seed, ...) gen_setting3(seed = seed),
  Setting4  = function(seed, ...) gen_setting4(seed = seed),
  Setting5  = function(seed, ...) gen_setting5(seed = seed),
  Setting6  = function(seed, ...) gen_setting6(seed = seed)
)
