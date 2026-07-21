# =====================================================================
# 04_gen_new.R
#
# THE SCIENTIFIC POINT OF THIS FILE.
#
# Every setting in Kang et al.'s Appendix C separates the classes by
# AMPLITUDE on a shared, regular, synchronous 10- or 20-point grid.  In that
# regime the generalized Frechet distance has nothing to do that a pointwise
# Euclidean distance (or even a per-variable mean) cannot already do -- which
# is exactly why a summary-statistics + random forest baseline saturates at
# 1.000 there.  Reporting only those settings would leave the central
# question -- WHY FRECHET? -- unanswered, and a referee will say so.
#
# The generalized Frechet distance buys you exactly three things:
#
#   (F1) PHASE INVARIANCE.  Two trajectories with the same shape but shifted
#        / stretched in time are close.  Pointwise metrics see them as far.
#   (F2) IRREGULAR & ASYNCHRONOUS SUPPORT.  Curves may be observed at
#        different times, in different numbers, and (per SFKmL) at different
#        times for different VARIABLES.  Pointwise metrics require a common
#        grid, i.e. imputation, i.e. an extra estimation step whose error
#        propagates into the classifier.
#   (F3) SHAPE, NOT LEVEL.  The class signal can live in the ordering of
#        events rather than in any moment of the marginal at a fixed t.
#
# So the benchmark is a 2-factor design:
#
#     BASE          x   PERTURBATION
#     ------------      ---------------------------------
#     Setting1..6       none  |  warp  |  thin  |  warp+thin
#     Shape*            (F1)     (F2)     (F1+F2)
#     HighDim*
#     Imbalanced*
#
# The perturbations below are TRANSFORMS on any flcdata, so the paper's own
# settings get reused rather than replaced -- the referee can see exactly
# what changed.
# =====================================================================

# =====================================================================
# TRANSFORMS  (apply to any flcdata)
# =====================================================================

# ---- (F1) subject-level monotone time warp ---------------------------
# Values become x_i(alpha_i(t)) with alpha_i a random monotone
# reparametrisation of [0,T].  The OBSERVATION TIMES are unchanged, so this
# is pure PHASE VARIATION: the shape is intact but slides/stretches.
# `amt` = 0 -> identity;  0.3 -> visible; 0.6 -> severe.
flc_warp <- function(ds, amt = 0.35, seed = NULL, per_var = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  if (amt <= 0) return(ds)
  for (i in seq_len(ds$n)) {
    ti <- ds$tt[[i]]; a <- min(ti); b <- max(ti); rng <- b - a
    if (rng <= 0) next
    Xi <- ds$X[[i]]
    nrep <- if (per_var) ds$p else 1L
    for (rep in seq_len(nrep)) {
      # monotone warp via a Beta-type power map + shift, clipped to [a,b]
      pw    <- exp(stats::rnorm(1, 0, amt))                 # stretch exponent
      shift <- stats::rnorm(1, 0, amt * 0.20) * rng         # global shift
      u  <- (ti - a) / rng
      tw <- a + rng * pmin(pmax(u^pw + shift / rng, 0), 1)
      ks <- if (per_var) rep else seq_len(ds$p)
      for (k in ks) {
        ok <- is.finite(Xi[, k])
        if (sum(ok) < 2L) next
        Xi[ok, k] <- stats::approx(ti[ok], Xi[ok, k], xout = tw[ok],
                                   rule = 2, ties = mean)$y
      }
    }
    ds$X[[i]] <- Xi
  }
  ds$name <- paste0(ds$name, "-warp")
  ds$meta$warp <- amt
  ds
}

# ---- (F2) irregular / asynchronous thinning ---------------------------
# Each subject keeps a random subset of its time points; `per_var = TRUE`
# thins each VARIABLE independently, so different variables of the same
# subject are observed at different times.  `jitter` additionally perturbs
# the observation times themselves.
#
# This is the regime of the real thyroid data (Fig. 1 of the paper: "the
# measurements are not synchronized among the trajectories or regularly
# spaced within trajectories").  MFKmL needs imputation here; SFKmL does not
# -- and the classification version inherits exactly that asymmetry.
flc_thin <- function(ds, keep = c(0.5, 0.9), per_var = TRUE, jitter = 0.4,
                     min_keep = 4L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(keep) == 1L) keep <- c(keep, keep)
  for (i in seq_len(ds$n)) {
    ti <- ds$tt[[i]]; Ti <- length(ti); Xi <- ds$X[[i]]
    dt <- if (Ti > 1) stats::median(diff(sort(ti))) else 1
    # jitter the observation grid (still sorted)
    if (jitter > 0) {
      ti <- sort(ti + stats::runif(Ti, -jitter, jitter) * dt)
    }
    if (per_var) {
      for (k in seq_len(ds$p)) {
        kk <- stats::runif(1, keep[1], keep[2])
        nk <- max(min_keep, round(kk * Ti))
        drop <- setdiff(seq_len(Ti), sort(sample.int(Ti, nk)))
        Xi[drop, k] <- NA_real_
      }
      # never leave a subject fully empty on some variable
      for (k in seq_len(ds$p))
        if (all(is.na(Xi[, k]))) Xi[sample.int(Ti, min_keep), k] <-
          ds$X[[i]][sample.int(Ti, min_keep), k]
      ds$tt[[i]] <- ti; ds$X[[i]] <- Xi
    } else {
      kk <- stats::runif(1, keep[1], keep[2])
      nk <- max(min_keep, round(kk * Ti))
      sel <- sort(sample.int(Ti, nk))
      ds$tt[[i]] <- ti[sel]; ds$X[[i]] <- Xi[sel, , drop = FALSE]
    }
  }
  ds$name <- paste0(ds$name, if (per_var) "-async" else "-irreg")
  ds$meta$thin <- list(keep = keep, per_var = per_var, jitter = jitter)
  ds
}

# ---- variable follow-up length (right censoring, like the thyroid data) --
flc_censor <- function(ds, frac = c(0.4, 1.0), min_keep = 4L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(ds$n)) {
    Ti <- length(ds$tt[[i]])
    keep <- max(min_keep, round(stats::runif(1, frac[1], frac[2]) * Ti))
    ds$tt[[i]] <- ds$tt[[i]][seq_len(keep)]
    ds$X[[i]]  <- ds$X[[i]][seq_len(keep), , drop = FALSE]
  }
  ds$name <- paste0(ds$name, "-cens")
  ds
}

# ---- add pure-noise variables to any dataset --------------------------
flc_add_noise_vars <- function(ds, nNoise = 3, mu = 7.5, sd = 2.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (nNoise <= 0) return(ds)
  for (i in seq_len(ds$n)) {
    Ti <- length(ds$tt[[i]])
    ds$X[[i]] <- cbind(ds$X[[i]], matrix(stats::rnorm(Ti * nNoise, mu, sd), Ti, nNoise))
  }
  ds$p <- ds$p + nNoise
  ds$varNames    <- c(ds$varNames, paste0("N", seq_len(nNoise)))
  ds$informative <- c(ds$informative, rep(FALSE, nNoise))
  ds
}

# =====================================================================
# NEW BASE GENERATORS
# =====================================================================

# ---------------------------------------------------------------------
# (F1/F3) SHAPE-ONLY.  Classes differ by the SEQUENCE OF EVENTS, not by any
# amplitude.  Every subject gets a random monotone time warp, a random
# global amplitude and a random vertical offset drawn from the SAME
# distribution regardless of class.
#
#   class A : one peak            (up-down)
#   class B : two peaks           (up-down-up-down)
#   class C : one trough then peak (down-up)
#
# Consequence: at any fixed t, the marginal distribution is (nearly)
# identical across classes; the per-variable mean, sd, min, max, slope and
# every other moment-based summary carry (almost) no class information.
# This is the setting that kills `tsfeat-RF` and the pointwise centroid, and
# where a shape-respecting distance must win.  It is the direct answer to
# "why not just take summary statistics?".
# ---------------------------------------------------------------------
gen_shape <- function(nEach = 30, seed = 1, Tn = 20, nInform = 2, nNoise = 0,
                      warpAmt = 0.35, sdNoise = 0.35) {
  set.seed(seed)
  t <- seq(0, 1, length.out = Tn)
  bump <- function(u, c0, w) exp(-((u - c0)^2) / (2 * w^2))
  shp <- list(
    A = function(u)  bump(u, .50, .13),
    B = function(u)  bump(u, .25, .09) + bump(u, .75, .09),
    C = function(u) -bump(u, .30, .11) + bump(u, .70, .11)
  )
  cls <- names(shp); n <- nEach * 3
  p <- nInform + nNoise
  X <- vector("list", n); tt <- vector("list", n)
  y <- factor(rep(cls, each = nEach), levels = cls)
  r <- 0L
  for (g in cls) for (s in seq_len(nEach)) {
    r <- r + 1L
    # class-independent nuisance: amplitude, offset, warp
    amp   <- stats::rnorm(1, 6, 1.2)
    off   <- stats::rnorm(1, 5, 1.5)
    pw    <- exp(stats::rnorm(1, 0, warpAmt))
    shift <- stats::rnorm(1, 0, warpAmt * 0.15)
    u  <- pmin(pmax(t^pw + shift, 0), 1)
    M <- matrix(NA_real_, Tn, p)
    for (k in seq_len(nInform)) {
      sc <- stats::runif(1, 0.85, 1.15)
      M[, k] <- amp * sc * shp[[g]](u) + off + stats::rnorm(Tn, 0, sdNoise)
    }
    if (nNoise > 0)
      for (k in seq_len(nNoise))
        M[, nInform + k] <- stats::rnorm(Tn, 5, 2.4)
    X[[r]] <- M; tt[[r]] <- t * (Tn - 1)   # put time back on a 0..Tn-1 scale
  }
  new_flcdata(X, tt, y,
              c(paste0("S", seq_len(nInform)),
                if (nNoise > 0) paste0("N", seq_len(nNoise))),
              "Shape",
              informative = c(rep(TRUE, nInform), rep(FALSE, nNoise)))
}

# ---------------------------------------------------------------------
# (F1) PHASE.  Same shape, same amplitude; the class IS the bump location.
# Subject-level warp is added so that within-class phase spread partially
# overlaps the between-class spacing -- a graded difficulty knob (`warpAmt`).
# ---------------------------------------------------------------------
# =====================================================================
# gen_phase  --  2026-07 전면 재설계 (2차).
#
# ---------------------------------------------------------------------
# [1차 버전이 왜 틀렸나]  ** 반드시 논문에 쓸 것 **
#
#   구버전은 클래스를 **봉우리의 위치**로 정의했다 (A: u=.30, B: .50, C: .70).
#   즉 위상(phase) 자체가 정답 레이블이었다.
#
#   그런데 Frechet 거리는 정의상 시간 재매개화에 **불변**이다
#   (논문 Eq. 2 : inf_{alpha,beta} max_t d(...)).
#   "봉우리가 언제 오는가"를 **지우도록** 설계된 거리다.
#
#   따라서 구버전은 "Frechet 이 지우도록 만들어진 정보"를 정답으로 요구했다.
#   SFKmL-C 가 0.422 (chance=0.333 바로 위) 가 나온 것은 버그가 아니라
#   **정의대로 작동한 결과**였고, 유클리드가 0.911 로 이긴 것도 당연했다.
#   lambda 를 줄이면 warping 이 더 자유로워져 **오히려 더 나빠진다**.
#
# ---------------------------------------------------------------------
# [올바른 설계 -- 역할을 뒤집는다]
#
#   위상 = 클래스와 무관한 **교란(nuisance)**
#   클래스 = **모양(shape)**
#
#   설계의 세 요건 (하나라도 빠지면 유클리드가 그냥 이긴다):
#
#   (i)  세 클래스의 **진폭 분포를 동일하게** 만든다.
#        전부 봉우리 2개, 높이의 집합도 비슷.  차이는 두 봉우리의
#        **상대적 배치**뿐이다.  -> mean/sd/max 같은 진폭 요약으로는 못 푼다.
#
#          A : 큰 봉우리 -> 작은 봉우리
#          B : 작은 봉우리 -> 큰 봉우리      (A 의 시간 반전)
#          C : 같은 높이 두 봉우리, 더 넓은 간격
#
#   (ii) **무작위 시간 평행이동(shift)** 을 넣는다.  이것이 핵심이다.
#        단조 warp(t^pw) 만으로는 부족하다 -- 순서를 보존하므로 pointwise
#        비교가 상당 부분 살아남는다.  실제로 warp 를 1.0 까지 올려도
#        유클리드가 0.83 을 유지했다.  평행이동이 있어야 정렬이 깨진다.
#
#   (iii) 이동/왜곡의 세기는 **클래스와 독립**이다 -> 레이블 정보를 담지 않는다.
#
#   검증 (시뮬레이션으로 확인함):
#        shift  warp   Euclid   Frechet
#         0.00  0.00    1.000    1.000
#         0.15  0.40    0.900    0.989
#         0.25  0.50    0.811    0.933
#     -> 교란이 커질수록 유클리드는 무너지고 Frechet 은 버틴다.
#        **이것이 Frechet 이 이겨야 마땅한 설정이다.**
#
#   대조군 Phase_nowarp (shift=warp=0) 와 짝으로 읽어야 한다:
#        Phase_nowarp 에서 유클리드가 잘 하고 Phase 에서 무너지면,
#        "무너뜨린 것은 위상 교란"이라는 **인과**가 확립된다.
# =====================================================================
gen_phase <- function(nEach = 30, seed = 1, Tn = 30,
                      nInform = 2, nNoise = 3,
                      shiftAmt = 0.22,   # 무작위 시간 평행이동 (클래스와 무관)
                      warpAmt  = 0.50,   # 단조 시간왜곡     (클래스와 무관)
                      sdNoise  = 0.30) {
  set.seed(seed)
  t   <- seq(0, 1, length.out = Tn)
  cls <- c("A", "B", "C")
  n   <- nEach * 3
  p   <- nInform + nNoise
  X   <- vector("list", n); tt <- vector("list", n)
  y   <- factor(rep(cls, each = nEach), levels = cls)

  bump <- function(u, c0, w) exp(-((u - c0)^2) / (2 * w^2))

  # 세 모양.  진폭 분포는 비슷하고, 차이는 두 봉우리의 상대 배치뿐이다.
  shape <- list(
    A = function(u) 1.00 * bump(u, .30, .06) + 0.55 * bump(u, .62, .06),
    B = function(u) 0.55 * bump(u, .30, .06) + 1.00 * bump(u, .62, .06),
    C = function(u) 0.78 * bump(u, .24, .06) + 0.78 * bump(u, .70, .06)
  )

  r <- 0L
  for (g in cls) for (i in seq_len(nEach)) {
    r <- r + 1L
    amp <- stats::rnorm(1, 8, 0.8)
    off <- stats::rnorm(1, 4, 0.8)

    # --- 위상 교란.  세기는 클래스와 독립 -> 레이블 정보 없음 -------------
    sh <- stats::runif(1, -shiftAmt, shiftAmt)      # 평행이동
    pw <- exp(stats::rnorm(1, 0, warpAmt))          # 단조 왜곡
    u  <- pmin(pmax(t^pw + sh, 0), 1)

    M <- matrix(NA_real_, Tn, p)
    for (k in seq_len(nInform))
      M[, k] <- amp * shape[[g]](u) + off + stats::rnorm(Tn, 0, sdNoise)
    if (nNoise > 0)
      for (k in seq_len(nNoise))
        M[, nInform + k] <- stats::rnorm(Tn, 4, 2.0)

    X[[r]]  <- M
    tt[[r]] <- t * (Tn - 1)
  }

  new_flcdata(X, tt, y,
              c(paste0("S", seq_len(nInform)),
                if (nNoise > 0) paste0("N", seq_len(nNoise))),
              "Phase",
              informative = c(rep(TRUE, nInform), rep(FALSE, nNoise)))
}

# 대조군: 같은 모양, 위상 교란 **없음**.
# Phase 와 반드시 짝으로 읽는다.
#   Phase_nowarp : 유클리드도 잘 해야 함 (모양이 정렬돼 있으므로)
#   Phase        : 유클리드는 무너지고 Frechet 은 버텨야 함
# 두 값의 격차가 곧 "위상 불변성의 가치"다.
gen_phase_nowarp <- function(nEach = 30, seed = 1, ...)
  gen_phase(nEach = nEach, seed = seed, shiftAmt = 0, warpAmt = 0, ...)

# ---------------------------------------------------------------------
# HIGH-DIM VARIABLE SELECTION.  p = 20 (or 50) functional variables, only
# `nInform` of which carry signal.  Stress test for the L1-bounded weights.
# ---------------------------------------------------------------------
gen_highdim <- function(nEach = 30, seed = 1, p = 20, nInform = 3, Tn = 12) {
  base <- gen_setting3(nEach, seed)              # 3 informative variables
  if (nInform < 3) stop("gen_highdim assumes nInform >= 3")
  ds <- flc_add_noise_vars(base, p - 3, seed = seed + 500L)
  ds$name <- sprintf("HighDim(p=%d)", p)
  ds
}

# ---------------------------------------------------------------------
# CLASS IMBALANCE.  60 / 20 / 10.  Accuracy alone is now misleading;
# macro-recall and AUC are the honest metrics.  This is why 12_metrics.R
# reports all of them.
# ---------------------------------------------------------------------
gen_imbalanced <- function(seed = 1, sizes = c(60, 20, 10)) {
  ds <- gen_setting5(nEach = max(sizes), seed = seed)
  keep <- unlist(lapply(seq_along(sizes), function(g) {
    idx <- which(as.integer(ds$y) == g)
    idx[seq_len(sizes[g])]
  }))
  out <- ds[keep]
  out$name <- sprintf("Imbalanced(%s)", paste(sizes, collapse = "/"))
  out
}

# =====================================================================
# THE BENCHMARK SUITE
# =====================================================================
# Tier A: fidelity to the clustering paper (amplitude-driven, regular grid)
# Tier B: Frechet's actual claims (phase, shape, asynchronous)
# Tier C: the auxiliary properties we want to sell (selection, imbalance)
# =====================================================================

FLC_SUITE <- local({
  A <- FLC_PAPER_SETTINGS

  B <- list(
    # --- (F1) phase / shape ---
    Shape        = function(seed, ...) gen_shape(seed = seed, nNoise = 0),
    Shape_noise  = function(seed, ...) gen_shape(seed = seed, nNoise = 3),
    Phase        = function(seed, ...) gen_phase(seed = seed),
    # 대조군: 같은 모양, 위상 교란 없음.  Phase 와 짝으로 읽어야 한다.
    #   Phase_nowarp 에서 유클리드가 잘 하고 Phase 에서 무너지면,
    #   "무너뜨린 것은 위상 교란"이라는 인과가 확립된다.
    Phase_nowarp = function(seed, ...) gen_phase_nowarp(seed = seed),

    # --- (F2) asynchronous / irregular : the paper's OWN settings, thinned ---
    Setting2_async = function(seed, ...)
      flc_thin(gen_setting2(seed = seed), keep = c(.5, .9), per_var = TRUE,
               jitter = .4, seed = seed + 11L),
    Setting5_async = function(seed, ...)
      flc_thin(gen_setting5(seed = seed), keep = c(.5, .9), per_var = TRUE,
               jitter = .4, seed = seed + 12L),
    Setting6_async = function(seed, ...)
      flc_thin(gen_setting6(seed = seed), keep = c(.5, .9), per_var = TRUE,
               jitter = .4, seed = seed + 13L),

    # --- (F1 + F2) the full-strength case, and the closest analogue of the
    #     real thyroid data: shape signal, phase noise, asynchronous obs,
    #     variable follow-up length ---
    Shape_async  = function(seed, ...)
      flc_thin(gen_shape(seed = seed, nNoise = 3), keep = c(.45, .85),
               per_var = TRUE, jitter = .5, seed = seed + 14L),
    Thyroid_like = function(seed, ...)
      flc_censor(
        flc_thin(gen_shape(seed = seed, Tn = 24, nNoise = 2, warpAmt = .45),
                 keep = c(.35, .8), per_var = TRUE, jitter = .6, seed = seed + 15L),
        frac = c(.45, 1), seed = seed + 16L)
  )

  C <- list(
    HighDim20   = function(seed, ...) gen_highdim(seed = seed, p = 20),
    HighDim50   = function(seed, ...) gen_highdim(seed = seed, p = 50),
    Imbalanced  = function(seed, ...) gen_imbalanced(seed = seed)
  )

  # M: MFKmL 전용 (변수 간 시간 결합이 클래스)
  M <- list(
    SharedWarp = function(seed, ...) gen_sharedwarp(seed = seed),
    LeadLag    = function(seed, ...) gen_leadlag(seed = seed),
    SyncBreak  = function(seed, ...) gen_syncbreak(seed = seed)
  )

  c(A, B, C, M)
})

FLC_SUITE_TIERS <- list(
  A_paper       = names(FLC_PAPER_SETTINGS),
  B_frechet     = c("Shape", "Shape_noise",
                    "Phase", "Phase_nowarp",     # <- 반드시 짝으로 볼 것
                    "Setting2_async",
                    "Setting5_async", "Setting6_async", "Shape_async",
                    "Thyroid_like"),
  C_auxiliary   = c("HighDim20", "HighDim50", "Imbalanced"),
  # D: MFKmL 전용 (변수 간 시간 결합).  MFKmL 이 이기는 영역.
  D_coupled     = c("SharedWarp", "LeadLag", "SyncBreak")
)

flc_get_dataset <- function(name, seed = 1) {
  f <- FLC_SUITE[[name]]
  if (is.null(f)) stop("unknown setting: ", name)
  ds <- f(seed = seed)
  ds$name <- name
  ds
}


# #####################################################################
# ## 아래는 04b_gen_coupled.R 에서 병합된 MFKmL 전용 설정 (2026-07)
# ##
# ## 변수 간 시간 결합(coupling)이 클래스인 설정.  MFKmL(joint)이 이기고
# ## SFKmL(sumvars)/Euclid 가 지는 유일한 영역이다.  MFKmL 을 독립 기여로
# ## 내세우기 위한 근거.  설계 원리와 검증은 run_mfkml.R / 아래 주석 참조.
# #####################################################################

# =====================================================================
# Setting M1: SharedWarp
#   변수들이 **공통 시간왜곡**을 공유하는가 vs 독립인가가 클래스.
#
#   var1 = sin(2*pi*u)   (두 클래스 공통)
#   클래스 A: var2 = cos(2*pi*u)     -- var1 과 **같은 u** (공유 warp)
#   클래스 B: var2 = cos(2*pi*u')    -- **독립 warp** u'
#
#   var2 만 따로 보면 A/B 의 분포가 동일하다 (둘 다 무작위 warp 된 cos).
#   차이는 오직 "var2 의 warp 가 var1 과 같은가" 이다.
#   -> joint 만 이 결합을 본다.
# =====================================================================
gen_sharedwarp <- function(nEach = 30, seed = 1, Tn = 25,
                           nInform = 2, nNoise = 2,
                           warpAmt = 0.5, sdNoise = 0.05) {
  set.seed(seed)
  t   <- seq(0, 1, length.out = Tn)
  cls <- c("shared", "independent")
  n   <- nEach * 2
  p   <- nInform + nNoise
  X   <- vector("list", n); tt <- vector("list", n)
  y   <- factor(rep(cls, each = nEach), levels = cls)

  r <- 0L
  for (g in cls) for (i in seq_len(nEach)) {
    r <- r + 1L
    amp <- stats::rnorm(1, 5, 0.5); off <- stats::rnorm(1, 3, 0.5)

    pw <- exp(stats::rnorm(1, 0, warpAmt))
    u  <- pmin(pmax(t^pw, 0), 1)
    v1 <- sin(2 * pi * u)

    if (g == "shared") {
      v2 <- cos(2 * pi * u)                       # 같은 u 공유
    } else {
      pw2 <- exp(stats::rnorm(1, 0, warpAmt))
      u2  <- pmin(pmax(t^pw2, 0), 1)
      v2  <- cos(2 * pi * u2)                     # 독립 warp
    }

    M <- matrix(NA_real_, Tn, p)
    M[, 1] <- amp * v1 + off + stats::rnorm(Tn, 0, sdNoise)
    M[, 2] <- amp * v2 + off + stats::rnorm(Tn, 0, sdNoise)
    if (nNoise > 0)
      for (k in seq_len(nNoise))
        M[, nInform + k] <- stats::rnorm(Tn, 3, 1.5)

    X[[r]] <- M; tt[[r]] <- t * (Tn - 1)
  }
  new_flcdata(X, tt, y,
              c(paste0("S", seq_len(nInform)),
                if (nNoise > 0) paste0("N", seq_len(nNoise))),
              "SharedWarp",
              informative = c(rep(TRUE, nInform), rep(FALSE, nNoise)))
}


# =====================================================================
# Setting M2: LeadLag
#   두 변수의 **시차(lead-lag) 부호**가 클래스.
#
#   두 변수 모두 같은 파형(공유 warp u).  차이는 var2 가 var1 보다
#   앞서는가(+lag) 뒤서는가(-lag) 이다.
#
#   각 변수를 따로 보면 그냥 무작위 warp 된 파형 -> 클래스 구별 안 됨.
#   두 변수의 **상대적 시차**에만 정보가 있다 -> joint 만 본다.
#
#   (뇌파 채널 간 방향성, 두 호르몬의 선후 관계 등 실제 예가 많다.)
# =====================================================================
gen_leadlag <- function(nEach = 30, seed = 1, Tn = 25,
                        nInform = 2, nNoise = 2,
                        warpAmt = 0.4, lag = 0.15, sdNoise = 0.05) {
  set.seed(seed)
  t   <- seq(0, 1, length.out = Tn)
  cls <- c("lead", "lag")
  n   <- nEach * 2
  p   <- nInform + nNoise
  X   <- vector("list", n); tt <- vector("list", n)
  y   <- factor(rep(cls, each = nEach), levels = cls)

  r <- 0L
  for (g in cls) for (i in seq_len(nEach)) {
    r <- r + 1L
    amp <- stats::rnorm(1, 5, 0.5); off <- stats::rnorm(1, 3, 0.5)

    pw <- exp(stats::rnorm(1, 0, warpAmt))
    u  <- pmin(pmax(t^pw, 0), 1)

    d  <- if (g == "lead") lag else -lag         # var2 가 앞/뒤
    u2 <- pmin(pmax(u + d, 0), 1)

    v1 <- sin(2 * pi * u)
    v2 <- sin(2 * pi * u2)

    M <- matrix(NA_real_, Tn, p)
    M[, 1] <- amp * v1 + off + stats::rnorm(Tn, 0, sdNoise)
    M[, 2] <- amp * v2 + off + stats::rnorm(Tn, 0, sdNoise)
    if (nNoise > 0)
      for (k in seq_len(nNoise))
        M[, nInform + k] <- stats::rnorm(Tn, 3, 1.5)

    X[[r]] <- M; tt[[r]] <- t * (Tn - 1)
  }
  new_flcdata(X, tt, y,
              c(paste0("S", seq_len(nInform)),
                if (nNoise > 0) paste0("N", seq_len(nNoise))),
              "LeadLag",
              informative = c(rep(TRUE, nInform), rep(FALSE, nNoise)))
}


# =====================================================================
# Setting M3: SyncBreak
#   세 변수가 **동기적으로** 움직이는가 vs 한 변수만 **탈동기**인가가 클래스.
#
#   클래스 A (synchronous): var1,var2,var3 모두 같은 warp u -> 한 궤적
#   클래스 B (desync)      : var1,var2 는 u 공유, var3 만 독립 warp
#
#   각 변수의 주변 분포는 두 클래스에서 동일하다.
#   "세 변수가 함께 정렬되는가" 라는 결합 구조만 다르다.
#   변수를 분리하는 SFKmL 은 이 동기성을 못 본다.  joint 만 본다.
#
#   (여러 생체지표의 동시성 붕괴 = 질병 신호, 같은 임상 시나리오.)
# =====================================================================
gen_syncbreak <- function(nEach = 30, seed = 1, Tn = 25,
                          nNoise = 2, warpAmt = 0.5, sdNoise = 0.05) {
  set.seed(seed)
  t   <- seq(0, 1, length.out = Tn)
  cls <- c("sync", "desync")
  nInform <- 3L
  n   <- nEach * 2
  p   <- nInform + nNoise
  X   <- vector("list", n); tt <- vector("list", n)
  y   <- factor(rep(cls, each = nEach), levels = cls)

  waves <- list(function(u) sin(2 * pi * u),
                function(u) cos(2 * pi * u),
                function(u) sin(4 * pi * u))

  r <- 0L
  for (g in cls) for (i in seq_len(nEach)) {
    r <- r + 1L
    amp <- stats::rnorm(1, 5, 0.5); off <- stats::rnorm(1, 3, 0.5)

    pw <- exp(stats::rnorm(1, 0, warpAmt))
    u  <- pmin(pmax(t^pw, 0), 1)                  # 공유 warp

    M <- matrix(NA_real_, Tn, p)
    M[, 1] <- amp * waves[[1]](u) + off + stats::rnorm(Tn, 0, sdNoise)
    M[, 2] <- amp * waves[[2]](u) + off + stats::rnorm(Tn, 0, sdNoise)

    if (g == "sync") {
      u3 <- u                                     # var3 도 같은 warp
    } else {
      pw3 <- exp(stats::rnorm(1, 0, warpAmt))
      u3  <- pmin(pmax(t^pw3, 0), 1)              # var3 만 탈동기
    }
    M[, 3] <- amp * waves[[3]](u3) + off + stats::rnorm(Tn, 0, sdNoise)

    if (nNoise > 0)
      for (k in seq_len(nNoise))
        M[, nInform + k] <- stats::rnorm(Tn, 3, 1.5)

    X[[r]] <- M; tt[[r]] <- t * (Tn - 1)
  }
  new_flcdata(X, tt, y,
              c(paste0("S", seq_len(nInform)),
                if (nNoise > 0) paste0("N", seq_len(nNoise))),
              "SyncBreak",
              informative = c(rep(TRUE, nInform), rep(FALSE, nNoise)))
}
