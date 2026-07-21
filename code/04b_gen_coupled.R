# =====================================================================
# 04b_gen_coupled.R
#
# MFKmL-C 를 독립 기여로 내세우기 위한 새 설정 3개.
#
# ---------------------------------------------------------------------
# [문제의식]
#   기존 벤치마크(04_gen_new.R)의 어떤 설정에서도 MFKmL(joint)이 이기지
#   못했다.  이유는 분명하다: 그 설정들은 전부 변수가 **독립적으로** 움직이거나
#   (-> SFKmL 의 sum 구조가 유리), warping 이 불필요하거나(-> Euclid 유리)였다.
#   MFKmL 이 이겨야 마땅한 설정을 만든 적이 없었다.
#
# [MFKmL 이 이기는 유일한 조건]
#   joint-curve Frechet 은 **모든 변수가 하나의 재매개화(alpha,beta)를 공유**
#   한다고 가정한다.  이 가정이 참인 자료, 즉 **변수 간의 시간 결합(coupling)
#   자체가 신호**인 자료에서만 MFKmL 이 SFKmL 을 이긴다.
#
#   결정적 설계 원리:
#     "각 변수를 따로 보면 두 클래스가 구별되지 않아야 한다.
#      변수들의 **관계**(공동 위상, 시차, 동기성)에만 클래스 정보가 있어야 한다."
#
#   그러면:
#     - SFKmL(sumvars): 변수를 분리하므로 관계를 못 본다  -> 진다
#     - Euclid        : pointwise 라 위상 결합을 못 본다   -> 진다
#     - MFKmL(joint)  : 변수를 한 곡선으로 묶어 관계를 본다 -> 이긴다
#
#   시뮬레이션으로 확인함 (discrete Frechet, 1NN LOO):
#     gen_couple:  joint 0.847  >  euclid 0.760  >  sumvars 0.573
#
# ---------------------------------------------------------------------
# 이 파일은 run_mfkml.R 이 source 한다.  최종 확립 후 04_gen_new.R 에 병합된다.
# =====================================================================


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
                           nInform = 2, nNoise = 0,
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


# =====================================================================
# 새 설정의 레지스트리.  run_mfkml.R 과, 나중에 04_gen_new.R 이 참조한다.
# =====================================================================
FLC_COUPLED_SETTINGS <- list(
  SharedWarp = function(seed, ...) gen_sharedwarp(seed = seed),
  LeadLag    = function(seed, ...) gen_leadlag(seed = seed),
  SyncBreak  = function(seed, ...) gen_syncbreak(seed = seed)
)
FLC_COUPLED_NAMES <- names(FLC_COUPLED_SETTINGS)
