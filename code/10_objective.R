# =====================================================================
# 10_objective.R
#
# SFKmL의 목적함수를 classification 문제로 옮기는 부분.
# 이 파일이 논문의 핵심 기여를 코드로 옮긴 지점이다.
#
# ---------------------------------------------------------------------
# [1] 원 논문(Kang et al. 2023, Eq. 5)의 SFKmL 목적함수
# ---------------------------------------------------------------------
#
#   maximize_{G, w}   sum_{i,i'} FD*(x_i, x_i')
#                   - sum_{k=1}^{K} (1/n_k) sum_{i,i' in G_k} FD*(x_i, x_i')
#
#   subject to   ||w||_2 <= 1,  ||w||_1 <= s,  w_j >= 0
#   where        FD*(x_i, x_i') = sum_{j=1}^{p} w_j * FD_{lambda,gamma_j}(x_ij, x_i'j)
#
#   - 첫 항  : 전체 pairwise distance 총합 (TOTAL)
#   - 둘째 항: cluster 내부 pairwise distance 총합 (WITHIN)
#   - 따라서 목적함수 = BCSS(w, G)  -- between-cluster 성분을 최대화
#   - 미지수가 (G, w) 두 개이므로 논문은 alternating optimization을 쓴다:
#         (a) w 고정 -> G에 대해 최소화 (7)   [K-medoids 배정]
#         (b) G 고정 -> w에 대해 최대화 (8)   [linear programming, Eq. 9]
#
# ---------------------------------------------------------------------
# [2] classification으로의 전환 (본 연구)
# ---------------------------------------------------------------------
#
# 핵심 관찰: 목적함수에서 G가 등장하는 곳은 "둘째 항"뿐이다.
# 지도학습에서는 G가 미지수가 아니라 **관측된 클래스 레이블 y**다.
# 즉 G := {G_1, ..., G_K} = {i : y_i = c_k}  로 고정된다.
#
# 결과적으로:
#   (a) 단계 [G에 대한 최적화]는 **소멸한다**.   <- alternating이 불필요
#   (b) 단계 [w에 대한 최적화]만 남고, 이는 CLOSED FORM (Eq. 9)이다.
#
#   따라서 supervised SFKmL의 목적함수는
#
#       maximize_w   a(y)' w
#       s.t.         ||w||_2 <= 1, ||w||_1 <= s, w_j >= 0
#
#       a_j(y) = sum_{i,i'} FD_{lambda,gamma_j}(x_ij, x_i'j)
#                - sum_{k} (1/n_k) sum_{i,i' in class k} FD_{lambda,gamma_j}(x_ij, x_i'j)
#
#   = "변수 j 하나만 썼을 때의 between-class 분리력".
#   해는 w = ST(a_+, Delta) / ||ST(a_+, Delta)||_2   (Eq. 9, optim_w).
#
#   이것은 클러스터링의 iterative 문제가 supervised에서는
#   **단일 볼록 문제(one-shot convex problem)** 로 붕괴함을 의미한다.
#   -> 논문에서 반드시 언급할 만한 구조적 사실이다.
#   -> 계산량도 O(iterations) 만큼 줄어든다.
#
# ---------------------------------------------------------------------
# [3] s의 선택: gap statistic -> classification loss
# ---------------------------------------------------------------------
#
# 클러스터링에는 정답이 없어서 논문은 permutation 기반 gap statistic으로
# s를 골랐다 (Eq. 10). 분류에는 정답(y)이 있다. 따라서 s는
#   (i)  gap      : 논문 방식 (거리행렬 permutation) -- 재현용
#   (ii) gaplab   : label permutation (supervised null) -- 변수별 p-value 부산물
#   (iii) cv      : inner-CV classification loss로 직접 선택  <- 권장
# 세 가지를 모두 구현하고 모두 보고한다 (05_weights.R).
#
# 원 논문 Setting 6에서 정보변수 V3가 탈락하던 문제(TPR 0.67)는
# gap이 보수적인 s*를 고르기 때문이며, cv 기준이 이를 완화한다.
# 이 대비 자체가 논문의 결과 하나가 된다.
#
# ---------------------------------------------------------------------
# [4] lambda 와 range-scale gamma  (논문 Sec. 2.2, Eq. 1, Eq. 11)
# ---------------------------------------------------------------------
#
# lambda (time-scale): affine transform A:(t,x) -> (lambda t, x).
#   Frechet distance는 시간축과 값축의 단위 비율에 민감하다.
#   논문 rule of thumb (Genolini et al. 2016, Sec. 3):
#       lambda = range(variable) / (range(time) * 0.1)
#   실데이터에서는 gap statistic을 lambda까지 확장해서 고름 (Eq. 11):
#       lambda_hat = argmax_l Gap( s_hat(lambda_l) )
#   본 코드에서는 lambda를 다음 세 가지로 지원한다:
#       "rule"  : 위 rule of thumb           (flc_lambda_rule)
#       "gap"   : 논문 Eq. 11 확장           (flc_tune_lambda, method="gap")
#       "cv"    : inner-CV 분류손실 최소화   (flc_tune_lambda, method="cv")  <- 권장
#   *** lambda는 반드시 TRAIN fold에서만 추정한다 (leakage 방지). ***
#
# gamma (range-scale): R:(t,x_1..x_p) -> (t, gamma_1 x_1, ..., gamma_p x_p),
#   gamma_j = c / range(variable j),  c = 100 (논문 Sec. 4).
#   진폭이 큰 변수가 Frechet distance를 지배하는 것을 막는다.
#   *** gamma도 TRAIN fold의 range로만 계산한다. ***
#   *** 중요: gamma는 w와 다른 역할이다.
#       gamma = 단위 통일 (전처리, 항상 적용)
#       w     = 변수 중요도 (모형 파라미터, 학습됨)
#       둘을 혼동하면 sparsity 해석이 무너진다. ***
# =====================================================================


# ---------------------------------------------------------------------
# supervised BCSS 벡터 a(y).  거리 캐시를 그대로 재활용한다.
# Dvar : n x n x p  (변수별 pairwise Frechet distance)
# 반환  : 길이 p 벡터. a_j > 0 이면 변수 j가 클래스를 분리하는 방향.
# ---------------------------------------------------------------------
sfkml_objective_a <- function(Dvar, y) {
  y <- factor(y); n <- dim(Dvar)[1]; p <- dim(Dvar)[3]
  a <- numeric(p)
  for (j in seq_len(p)) {
    Dj    <- Dvar[, , j]
    total <- sum(Dj)                              # sum_{i,i'}
    within <- 0
    for (c in levels(y)) {
      idx <- which(y == c); nk <- length(idx)
      if (nk == 0L) next
      within <- within + sum(Dj[idx, idx, drop = FALSE]) / nk   # (1/n_k) sum_{i,i' in G_k}
    }
    a[j] <- total - within                        # == BCSS_j
  }
  a
}

# 목적함수 값 자체 (보고용).  U(w) = a' w
sfkml_objective_value <- function(Dvar, y, w) sum(sfkml_objective_a(Dvar, y) * w)

# ---------------------------------------------------------------------
# supervised SFKmL 전체 해 (one-shot, no alternating)
#   1. a <- sfkml_objective_a(Dvar, y)
#   2. w <- optim_w(a, s)                 [Eq. 9]
#   3. medoid_k <- argmin_{i in class k} sum_{i' in class k} sum_j w_j D_j(i,i')
# ---------------------------------------------------------------------
sfkml_solve_supervised <- function(Dvar, y, s) {
  a <- sfkml_objective_a(Dvar, y)
  w <- optim_w(a, s)                     # 05_weights.R, 논문 Eq.9의 closed form
  y <- factor(y)
  Dw <- matrix(0, dim(Dvar)[1], dim(Dvar)[2])
  for (j in which(w > 1e-12)) Dw <- Dw + w[j] * Dvar[, , j]
  medoid <- vapply(levels(y), function(c) {
    idx <- which(y == c)
    if (length(idx) <= 2L) return(idx[1])
    sub <- Dw[idx, idx, drop = FALSE]
    diag(sub) <- 0
    idx[which.min(colSums(sub))]
  }, 1L)
  list(a = a, w = w, s = s, medoid = medoid, objective = sum(a * w),
       Dw = Dw)
}


# =====================================================================
# lambda 튜닝.  전부 TRAIN fold 안에서만 호출된다.
# =====================================================================

# 후보 lambda 격자.  rule of thumb 값을 중심으로 log-scale로 퍼뜨린다.
# lambda 후보 격자.
#
# lambda 는 시간축과 값축의 "환율"이다 (논문 Eq. 1, A:(t,x)->(lambda t, x)).
#   lambda 大 -> 시간축 이동이 비싸다 -> warping 을 안 한다 -> Euclid 에 가까워짐
#   lambda 小 -> 시간축 이동이 싸다   -> 자유롭게 warping -> 위상 불변에 가까워짐
#
# 첫 실행에서 Phase 설정(클래스=bump 위치)의 SFKmL-C 가 0.422 로 chance 수준이
# 나왔다.  rule of thumb 의 lambda 가 위상 문제에는 너무 커서, Frechet 이
# 위상 차이를 "거리"로 벌점화해 버린 것으로 의심된다.  그래서 격자를 아래쪽으로
# 넓게 (0.05배까지) 잡는다.  이것이 논문 Eq. 11 이 lambda 를 튜닝하라고 한 이유다.
# 격자를 위/아래로 모두 넓게 잡는다.
#
#   아래쪽 (0.05배) : Phase 처럼 위상 교란이 있는 자료.  lambda 가 작아야
#                     시간축 이동이 싸져서 warping 이 자유로워진다.
#                     실제로 Phase 에서 0.017~0.17 이 선택되었다.
#
#   위쪽 (16배)     : joint 거리(MFKmL) 용.  joint 의 local cost 는
#                       sqrt( (lambda*dt)^2 + sum_k gamma_k^2 (x_k-y_k)^2 )
#                     이라 값 쪽 항이 **p 개 누적**된다.  시간축이 상대적으로
#                     묻히므로 최적 lambda 가 sumvars 보다 커야 한다.
#                     구 격자(최대 8배)가 그 지점에 못 닿았을 수 있다.
flc_lambda_grid <- function(ds,
                            mult = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)) {
  base <- flc_lambda_rule(ds)
  sort(unique(base * mult))
}

# --- (a) 논문 Eq. 11 방식: lambda를 gap statistic으로 -------------------
flc_tune_lambda_gap <- function(ds, lambda_grid = NULL, sumOrMax = "max",
                                nperm = 20L, seed = 1L, ncores = 1L) {
  if (is.null(lambda_grid)) lambda_grid <- flc_lambda_grid(ds)
  res <- lapply(lambda_grid, function(lam) {
    Dvar <- flc_dist_var(ds, lambda = lam, sumOrMax = sumOrMax)
    g <- gap_select_s(Dvar, ds$y, perm = "dist", nperm = nperm, seed = seed)
    data.frame(lambda = lam, s = g$s, gap = g$gap_max)
  })
  res <- do.call(rbind, res)
  list(lambda = res$lambda[which.max(res$gap)],
       s = res$s[which.max(res$gap)], table = res)
}

# --- (b) 권장: inner-CV 분류손실로 lambda 선택 -------------------------
# 분류 문제이므로 gap 대신 직접 목적(정확도)을 최적화하는 것이 정직하다.
flc_tune_lambda_cv <- function(ds, lambda_grid = NULL, sumOrMax = "max",
                               K = 5L, seed = 1L, s_select = "cv") {
  if (is.null(lambda_grid)) lambda_grid <- flc_lambda_grid(ds)
  y <- ds$y; fold <- strat_folds(y, K, seed)
  res <- lapply(lambda_grid, function(lam) {
    Dvar <- flc_dist_var(ds, lambda = lam, sumOrMax = sumOrMax)
    acc <- numeric(K)
    for (f in seq_len(K)) {
      tr <- which(fold != f); te <- which(fold == f)
      sol <- NULL
      Wf <- flc_learn_weights(Dvar[tr, tr, , drop = FALSE], y[tr],
                              s_select = s_select, seed = seed)
      Dc <- Dvar[te, tr, , drop = FALSE]
      Dw <- matrix(0, length(te), length(tr))
      for (j in which(Wf$w > 1e-12)) Dw <- Dw + Wf$w[j] * Dc[, , j]
      Dtr <- matrix(0, length(tr), length(tr))
      for (j in which(Wf$w > 1e-12))
        Dtr <- Dtr + Wf$w[j] * Dvar[tr, tr, j]
      med <- vapply(levels(y), function(c) {
        idx <- which(y[tr] == c)
        if (length(idx) <= 2L) idx[1]
        else idx[which.min(colSums(Dtr[idx, idx, drop = FALSE]))] }, 1L)
      pr <- levels(y)[max.col(-Dw[, med, drop = FALSE], ties.method = "first")]
      acc[f] <- mean(pr == as.character(y[te]))
    }
    data.frame(lambda = lam, acc = mean(acc), se = stats::sd(acc) / sqrt(K))
  })
  res <- do.call(rbind, res)
  list(lambda = res$lambda[which.max(res$acc)], table = res)
}

flc_tune_lambda <- function(ds, method = c("rule", "cv", "gap"), ...) {
  method <- match.arg(method)
  switch(method,
    rule = list(lambda = flc_lambda_rule(ds), table = NULL),
    cv   = flc_tune_lambda_cv(ds, ...),
    gap  = flc_tune_lambda_gap(ds, ...))
}


# =====================================================================
# ctx : 하나의 (dataset, fold, lambda, gamma, cache) 실행 문맥.
# 모든 fit_*/predict_* 가 이것 하나만 받는다.  누수 방지의 단일 지점.
# =====================================================================
flc_make_ctx <- function(ds, tr, te, cache, lambda, gamma, sumOrMax = "max",
                         grid = NULL, seed = 1L) {
  if (is.null(grid)) grid <- flc_common_grid(ds)
  list(tr = tr, te = te, cache = cache, lambda = lambda, gamma = gamma,
       sumOrMax = sumOrMax, grid = grid, seed = seed, p = ds$p)
}
