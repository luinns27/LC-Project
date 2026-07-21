# =====================================================================
# 00_setup.R -- packages, C++ compilation, parallel backend
# =====================================================================

FLC_PKGS_REQUIRED <- c("Rcpp", "foreach", "doParallel", "iterators")

# optional; each unlocks a family of baselines.  Missing ones are skipped
# gracefully and recorded as NA in the results, never as a crash.
FLC_PKGS_OPTIONAL <- c(
  kernlab      = "distance-substitution kernel SVM",
  randomForest = "flat-RF / catch22-RF / TSF",
  glmnet       = "ROCKET ridge classifier, functional GLM",
  class        = "flat kNN",
  e1071        = "flat SVM",
  MASS         = "LDA on FPC scores",
  fda.usc      = "classif.knn / kernel / glm / depth / DD (functional)",
  fda          = "B-spline bases for MFPCA / MFPLS",
  mrfDepth     = "MFHD, bagdistance -> DistSpace (Hubert et al. 2017)",
  ddalpha      = "DD-classifier (Li, Cuesta-Albertos & Liu 2012)",
  Rcatch22     = "catch22 features (bake-off-redux 'feature' category)",
  pls          = "functional PLS-DA (Preda, Saporta & Levedere 2007)"
)

flc_have <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))

flc_setup <- function(src_dir = file.path(getOption("FLC_ROOT", "."), "src"),
                      nthreads = max(1L, parallel::detectCores() - 1L),
                      verbose  = TRUE) {

  miss <- FLC_PKGS_REQUIRED[!vapply(FLC_PKGS_REQUIRED, flc_have, logical(1))]
  if (length(miss))
    stop("Missing required packages: ", paste(miss, collapse = ", "),
         "\n  install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")

  suppressPackageStartupMessages({
    library(Rcpp); library(foreach); library(doParallel)
  })

  cpp <- file.path(src_dir, "flc_core.cpp")
  if (!file.exists(cpp)) stop("cannot find ", cpp,
                              "\n  set options(FLC_ROOT = '/path/to/flc') first.")
  Rcpp::sourceCpp(cpp, rebuild = FALSE)

  omp <- flc_omp_threads()
  .flc$omp      <- omp
  .flc$nthreads <- if (omp > 0) min(nthreads, omp) else 1L

  if (verbose) {
    message("flc: C++ core compiled.  OpenMP threads available: ",
            if (omp > 0) omp else "0 (OpenMP OFF -- see note below)")
    if (omp == 0)
      message("     -> on macOS install libomp; on Windows use Rtools. ",
              "Everything still works, just single-threaded.")
    opt <- names(FLC_PKGS_OPTIONAL)
    ok  <- vapply(opt, flc_have, logical(1))
    if (any(!ok))
      message("flc: optional packages missing (those baselines will be skipped): ",
              paste(opt[!ok], collapse = ", "))
  }
  invisible(TRUE)
}

# small mutable environment for global state
.flc <- new.env(parent = emptyenv())
.flc$omp      <- 0L
.flc$nthreads <- 1L

flc_nthreads <- function() as.integer(.flc$nthreads %||% 1L)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- parallel backend for the OUTER loop (replications x settings) -----
# We deliberately do NOT nest foreach inside OpenMP.  Pick one:
#   level = "outer": foreach over jobs, OpenMP forced to 1 inside workers
#   level = "inner": sequential jobs, OpenMP inside the distance kernels
# "outer" is faster when there are many jobs (the usual case).
flc_parallel_start <- function(ncores = max(1L, parallel::detectCores() - 1L),
                               level  = c("outer", "inner")) {
  level <- match.arg(level)
  if (level == "inner" || ncores <= 1L) {
    foreach::registerDoSEQ()
    .flc$par_level <- "inner"
    return(invisible(NULL))
  }
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  .flc$cl        <- cl
  .flc$par_level <- "outer"
  invisible(cl)
}

flc_parallel_stop <- function() {
  if (!is.null(.flc$cl)) { try(parallel::stopCluster(.flc$cl), silent = TRUE); .flc$cl <- NULL }
  foreach::registerDoSEQ()
  invisible(NULL)
}

# threads to hand to the C++ kernels *right now*
flc_kernel_threads <- function() {
  if (identical(.flc$par_level, "outer")) 1L else flc_nthreads()
}
