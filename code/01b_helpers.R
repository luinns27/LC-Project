# =====================================================================
# 01b_helpers.R -- 다른 파일들이 참조하는 보조 함수 모음.
# 파일명이 01b 인 이유: run_master.R 이 R/ 아래를 알파벳순으로 source 하므로
# 01_data.R 직후, 02_frechet.R 이전에 로드되어야 한다.
# =====================================================================

# ---------------------------------------------------------------------
# 공통 시간격자.  보간이 필요한 방법들(euclid, depth, ROCKET, FPCA 등)이
# 전부 이걸 쓴다.  ragged 데이터의 관측 시각 전체를 훑어 [min, max]를 잡고
# 등간격 len개로 자른다.
#
# 중요: 이 격자는 "보간이 필요한 방법"에만 쓰인다.
#       Frechet/DTW 계열은 원본 ragged 구조를 그대로 쓰므로 이걸 통과하지
#       않는다.  이 비대칭이 비동기 setting에서의 성능차의 원인이며,
#       측정하려는 대상 그 자체다.
# ---------------------------------------------------------------------
flc_common_grid <- function(ds, len = 100L) {
  tr <- range(unlist(ds$tt), na.rm = TRUE)
  if (!all(is.finite(tr)) || tr[1] == tr[2]) tr <- c(0, 1)
  seq(tr[1], tr[2], length.out = len)
}

# NULL 안전 기본값 (00_setup.R 에도 있지만 방어적으로 재정의)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------
# OpenMP 스레드 수.  C++ 쪽 flc_omp_threads()가 컴파일되지 않았을 때를 대비.
# ---------------------------------------------------------------------
flc_omp_available <- function() {
  exists("flc_omp_threads", mode = "function") &&
    isTRUE(try(flc_omp_threads() > 1L, silent = TRUE))
}

# ---------------------------------------------------------------------
# 목록 -> 안전한 rbind (열 이름이 달라도 죽지 않음)
# ---------------------------------------------------------------------
rbind_fill <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(NULL)
  nm <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    miss <- setdiff(nm, names(d))
    for (m in miss) d[[m]] <- NA
    d[, nm, drop = FALSE]
  }))
}

# ---------------------------------------------------------------------
# 진행 상황 출력
# ---------------------------------------------------------------------
flc_msg <- function(fmt, ...) cat(sprintf(fmt, ...))


# =====================================================================
# kNN 다수결 확률 행렬.  (BUG FIX 2026-07)
#
# 원래 코드는 이랬다:
#     P <- t(apply(ord[, 1:kk, drop=FALSE], 1, function(ix) table(...)/kk))
#     if (kk == 1L) P <- matrix(P, ncol = length(cls), byrow = TRUE)   # <-- 버그
#
# apply(X, 1, f) 에서 f 가 길이 K 벡터를 돌려주면 결과는 K x n_te 행렬이다
# (결과가 열에 쌓인다).  거기에 t() 를 하면 n_te x K 로 이미 올바르다.
# 그런데 두 번째 줄이 그 올바른 행렬을 다시 column-major 로 읽어
# byrow=TRUE 로 재조립한다 -> 행과 열이 뒤섞인다.
#
# 그 결과 k=1 인 모든 분류기(1NN-DTW, 1NN-Euclid)가 chance 수준(0.33)이
# 나왔다.  1NN-DTW 는 Bagnall et al.(2017) 의 표준 기준선인데 chance 가
# 나올 리 없다 -- 그게 버그의 신호였다.
#
# 아래 함수는 apply/t() 를 아예 쓰지 않고 명시적으로 행을 채운다.
# k=1 이든 k>1 이든 같은 경로를 타므로 특수분기가 없다.
#
#   D    : n_te x n_tr 거리행렬
#   ytr  : 길이 n_tr 인 훈련 레이블 (factor)
#   k    : 이웃 수
#   cls  : 클래스 레벨 (levels(ytr))
# 반환 : n_te x K 확률행렬 (열이름 = cls)
# =====================================================================
knn_vote <- function(D, ytr, k, cls = levels(ytr)) {
  D   <- as.matrix(D)
  n   <- nrow(D)
  K   <- length(cls)
  ytr <- factor(ytr, levels = cls)
  k   <- max(1L, min(as.integer(k), ncol(D)))

  P <- matrix(0, n, K, dimnames = list(NULL, cls))
  for (i in seq_len(n)) {
    di <- D[i, ]
    di[!is.finite(di)] <- Inf
    nn <- order(di)[seq_len(k)]           # 가장 가까운 k개의 train 인덱스
    tb <- tabulate(as.integer(ytr[nn]), nbins = K)
    P[i, ] <- tb / k
  }
  P
}
