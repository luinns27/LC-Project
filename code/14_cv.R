# =====================================================================
# 14_cv.R -- 교차검증 하네스.
#
# 설계 원칙 (전부 논문에서 방어 가능해야 함):
#
#  [1] LEAKAGE 방지.  fold마다 TRAIN에서만 추정하는 것:
#        - gamma (range-scale) : TRAIN 변수 range
#        - lambda (time-scale) : TRAIN 위에서 rule / cv / gap
#        - w (sparse weights), s
#        - k (kNN), C/gamma (SVM), q (subspace), 표준화 상수
#      TEST는 오직 "이미 고정된 파라미터로 변환 + 예측"만 한다.
#
#  [2] 거리 캐시.  n x n x p 변수별 Frechet 배열 + joint/dtw/euclid를
#      데이터셋 전체에 대해 **한 번** 계산하고, fold는 인덱싱만 한다.
#      정당성: 거리는 (x_i, x_j)의 함수일 뿐 y를 쓰지 않는다.
#      즉 unsupervised transform이므로 fold 간 정보 누출이 없다.
#      (표준화/가중치처럼 y나 fold-특이 통계를 쓰는 것은 캐시하지 않는다.)
#
#      gamma가 fold마다 달라지면 거리도 달라지는데, 현재는 전체 데이터의
#      range로 gamma를 한 번 고정한다("global" 모드).  range는 y를 쓰지
#      않는 통계이므로 실무상 허용되고 선행 문헌 대부분이 이렇게 한다.
#      엄격한 fold별 재계산이 필요하면 flc_prepare()를 fold 안에서 부르면
#      되지만 비용이 크다.  최종 표는 두 방식을 모두 돌려 값이 같은지
#      보고하는 것이 안전하다 (부록 민감도 분석).
#
#  [3] 모든 방법은 실패해도 죽지 않는다.  try() -> error 문자열 + NA metric.
#
#  [4] fit 객체의 필드는 신뢰하지 않는다.  방법마다 $s, $w 가 전혀 다른
#      의미로 존재할 수 있다 (예: tabular 방법의 내부 표준편차 벡터).
#      따라서 스칼라/길이 p 인지 **반드시 검사한 뒤에만** 결과표에 넣는다.
#      (이 검사가 없으면 data.frame 대입에서 죽는다.)
# =====================================================================


# ---------------------------------------------------------------------
# 데이터셋 하나에 대한 캐시 준비
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# lambda 를 어떻게 고르는가 -- 누수 문제 (2026-07 수정)
#
# 문제: 거리 캐시는 lambda 에 의존한다.  그런데 lambda 를 CV 로 고르려면
#       레이블 y 가 필요하다.  fold 를 나누기 *전에* 전체 데이터로 lambda 를
#       고르면 test fold 의 레이블이 lambda 선택에 스며든다 = 누수.
#
# 해결: lambda 후보마다 캐시를 **미리 전부** 계산해 둔다 (거리는 y 를 쓰지
#       않으므로 이 단계는 누수가 아니다).  그런 다음 fold 안에서, TRAIN
#       레이블만으로 최적 lambda 를 골라 해당 캐시를 꺼내 쓴다.
#
#       lambda_method = "rule" : 후보 1개 (Genolini rule of thumb).  누수 없음.
#                       "cv"   : 후보 여러 개.  fold 안에서 TRAIN 으로 선택.
#                       "gap"  : 후보 여러 개.  gap statistic (레이블 미사용,
#                                논문 Eq. 11) 이므로 전역 선택도 무방.
#
# 비용: 후보 L 개면 캐시 계산이 L 배.  그래서 기본 격자를 6개로 제한했다.
# ---------------------------------------------------------------------
flc_prepare <- function(ds, lambda_method = c("rule", "cv", "gap", "align"),
                        sumOrMax = "max", grid_len = 100L,
                        lambda_grid = NULL, verbose = TRUE) {
  lambda_method <- match.arg(lambda_method)
  t0 <- Sys.time()

  gamma <- flc_gamma(ds)                 # c / range_j   (c = 100)
  dsr   <- flc_rescale(ds, gamma)        # gamma 적용된 데이터

  # 공통격자.  보간이 필요한 방법(euclid/depth/FPCA 등)만 이걸 쓴다.
  # Frechet/DTW 계열은 원본 ragged 구조를 그대로 쓰므로 통과하지 않는다.
  grid <- flc_common_grid(dsr, len = grid_len)

  # ---- lambda 후보 ----------------------------------------------------
  if (identical(lambda_method, "rule")) {
    lam_grid <- flc_lambda_rule(dsr)                 # 단일 값
  } else {
    lam_grid <- if (is.null(lambda_grid)) flc_lambda_grid(dsr) else lambda_grid
  }

  if (verbose) cat(sprintf("  [prepare] gamma=%s  lambda=%s  ... ",
                           paste(sprintf("%.3g", gamma), collapse = ","),
                           if (length(lam_grid) == 1L) sprintf("%.4g", lam_grid)
                           else sprintf("%d cands", length(lam_grid))))

  # ---- 후보마다 캐시를 미리 계산 (y 를 쓰지 않으므로 누수 아님) --------
  caches <- lapply(lam_grid, function(lam)
    flc_cache(dsr, lambda = lam, sumOrMax = sumOrMax))
  names(caches) <- sprintf("%.6g", lam_grid)

  # gap 은 레이블을 쓰지 않으므로 (논문 Eq. 11) 전역 선택이 정당하다.
  # cv 는 레이블을 쓰므로 fold 안에서 골라야 한다 -> 여기서는 고르지 않고,
  #    일단 rule 값을 기본으로 두고 flc_make_ctx() 가 fold 안에서 재선택한다.
  sel <- 1L
  gap_tab <- NULL
  if (identical(lambda_method, "gap") && length(lam_grid) > 1L) {
    g <- vapply(seq_along(lam_grid), function(l) {
      gs <- try(gap_select_s(caches[[l]]$Dvar, dsr$y, perm = "dist",
                             nperm = 20L, seed = 1L), silent = TRUE)
      if (inherits(gs, "try-error")) -Inf else gs$gap_max
    }, 0)
    gap_tab <- data.frame(lambda = lam_grid, gap = g)
    sel <- which.max(g)
  }

  if (verbose) cat(sprintf("done (%.1fs)\n",
                           as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  list(ds = dsr, gamma = gamma, grid = grid, sumOrMax = sumOrMax,
       lambda_method = lambda_method,
       lambda_grid = lam_grid,
       caches = caches,
       lambda = lam_grid[sel],            # 기본값 (rule/gap 에서 확정)
       cache  = caches[[sel]],            # 기본 캐시 (하위호환)
       gap_table = gap_tab,
       prep_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}


# ---------------------------------------------------------------------
# fold 안에서 TRAIN 레이블만으로 lambda 를 고른다 (lambda_method="cv" 일 때).
#
# 선택 기준: inner-CV 분류손실.  분류 문제이므로 gap 대신 실제 목적(정확도)을
# 직접 최적화하는 것이 정직하다.  분류기는 SFKmL-C 의 규칙(medoid) 을 쓴다.
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# fold 안에서 TRAIN 레이블만으로 lambda 를 고른다.
#
# *** 왜 거리마다 따로 골라야 하는가 (2026-07) ***
#
# 구버전은 fold 당 lambda 를 **하나만** 골라 모든 방법이 공유했고, 그 기준은
# SFKmL 의 규칙(sumvars + medoid)이었다.  이것은 MFKmL 에 불공정하다.
#
#   sumvars : lambda 가 각 1차원 곡선의 시간축에만 작용한다.
#             local cost = sqrt( (lambda*dt)^2 + (x_j - y_j)^2 )   [변수마다 따로]
#
#   joint   : lambda 가 (p+1)차원 공간의 **한 축**이고, 나머지 p 개 축은 gamma 로
#             스케일된 변수들이다.
#             local cost = sqrt( (lambda*dt)^2 + sum_k gamma_k^2 (x_k - y_k)^2 )
#             -> 값 쪽 항이 p 개 누적되므로 시간축과의 균형점이 **다르다**.
#                p 가 커질수록 최적 lambda 도 커져야 한다.
#
# 즉 두 거리의 최적 lambda 가 같을 이유가 없다.  같은 값을 강제하면
# "MFKmL 이 진다"는 결론이 구현 아티팩트일 수 있다.  리뷰어가 반드시 묻는다.
#
# 그래서 criterion 별로 따로 고른다:
#   "sumvars" : SFKmL 규칙(가중합 거리 + medoid)으로 평가
#   "joint"   : MFKmL 규칙(joint 거리 + medoid)으로 평가
# 두 lambda 를 모두 저장하고, 방법마다 자기 거리에 맞는 것을 쓴다.
# ---------------------------------------------------------------------
.pick_lambda_train <- function(prep, tr, seed = 1L, innerK = 3L,
                               criterion = c("sumvars", "joint"),
                               method = c("cv", "align")) {
  criterion <- match.arg(criterion)
  method    <- match.arg(method)
  L <- length(prep$lambda_grid)
  if (L == 1L) return(1L)

  y    <- prep$ds$y[tr]
  cls  <- levels(y)

  # ===================================================================
  # 층위 B: supervised 거리 파라미터 선택 (거리 공식은 불변)
  #
  # method = "align":  거리행렬 target alignment 를 최대화하는 lambda.
  #
  #     A(l) = <-D_l, T>_F / ( ||D_l||_F * ||T||_F )
  #     T_ij = +1 (같은 클래스), -1 (다른 클래스)
  #
  #   거리행렬 D_l 이 레이블 구조 T 와 정렬될수록 A 가 크다: 같은 클래스
  #   쌍의 거리가 작고 다른 클래스 쌍의 거리가 크면 <-D,T> 가 커진다.
  #   Fréchet 공식(원 논문 정의)은 전혀 바뀌지 않는다.  단지 원 논문이
  #   gap statistic(unsupervised)으로 고르던 lambda 를, 분류에서는 레이블을
  #   써서 supervised 로 고를 뿐이다.  metric 성질은 모든 lambda>0 에서
  #   원 논문 Theorem 1 로 이미 보장되므로 재증명이 불필요하다.
  #
  #   누출 없음: tr(=TRAIN) 안의 거리행렬만 쓴다.  align 은 held-out 이
  #   아니라 TRAIN 전체의 거리-레이블 정렬을 보므로 inner-fold 도 불필요.
  # ===================================================================
  if (method == "align") {
    ytr <- as.integer(y)
    n   <- length(ytr)
    T   <- matrix(-1, n, n)
    for (ci in seq_along(cls)) { ix <- which(ytr == ci); T[ix, ix] <- 1 }
    Tn  <- sqrt(sum(T * T))

    ali <- vapply(seq_len(L), function(l) {
      ca <- prep$caches[[l]]
      if (criterion == "sumvars") {
        # w 를 TRAIN 안에서 학습한 뒤 가중합 거리행렬을 만든다.
        Dvar <- ca$Dvar[tr, tr, , drop = FALSE]
        W <- try(flc_learn_weights(Dvar, y, s_select = "cv", seed = seed),
                 silent = TRUE)
        if (inherits(W, "try-error")) return(-Inf)
        D <- matrix(0, n, n)
        for (k in which(W$w > 1e-12)) D <- D + W$w[k] * Dvar[, , k]
      } else {
        D <- ca$Djoint[tr, tr, drop = FALSE]
      }
      Dn <- sqrt(sum(D * D))
      if (!is.finite(Dn) || Dn <= 0) return(-Inf)
      sum(-D * T) / (Dn * Tn)       # <-D, T>_F / (||D|| ||T||)
    }, 0)

    if (all(!is.finite(ali))) return(1L)
    return(which.max(ali))
  }

  # ===================================================================
  # method = "cv" (기존): inner-CV 분류손실로 lambda 선택.
  # ===================================================================
  fold <- strat_folds(y, min(innerK, min(table(y))), seed)

  acc <- vapply(seq_len(L), function(l) {
    ca <- prep$caches[[l]]

    a <- vapply(sort(unique(fold)), function(f) {
      i1 <- which(fold != f); i2 <- which(fold == f)

      if (criterion == "sumvars") {
        # SFKmL: 변수별 거리의 가중합.  w 도 TRAIN 안에서 학습.
        Dvar <- ca$Dvar[tr, tr, , drop = FALSE]
        W <- try(flc_learn_weights(Dvar[i1, i1, , drop = FALSE], y[i1],
                                   s_select = "cv", seed = seed), silent = TRUE)
        if (inherits(W, "try-error")) return(NA_real_)
        Dtr <- matrix(0, length(i1), length(i1))
        Dte <- matrix(0, length(i2), length(i1))
        for (k in which(W$w > 1e-12)) {
          Dtr <- Dtr + W$w[k] * Dvar[i1, i1, k]
          Dte <- Dte + W$w[k] * Dvar[i2, i1, k]
        }
      } else {
        # MFKmL: joint curve 거리.  가중치 없음 (구조상 희소화 불가).
        Dj  <- ca$Djoint[tr, tr, drop = FALSE]
        Dtr <- Dj[i1, i1, drop = FALSE]
        Dte <- Dj[i2, i1, drop = FALSE]
      }

      med <- vapply(cls, function(c) {
        ii <- which(y[i1] == c)
        if (length(ii) <= 2L) return(ii[1])
        sub <- Dtr[ii, ii, drop = FALSE]; diag(sub) <- 0
        ii[which.min(colSums(sub))]
      }, 1L)
      pr <- cls[max.col(-Dte[, med, drop = FALSE], ties.method = "first")]
      mean(pr == as.character(y[i2]))
    }, 0)

    mean(a, na.rm = TRUE)
  }, 0)

  if (all(!is.finite(acc))) return(1L)
  which.max(acc)
}


# 방법이 어떤 거리를 쓰는가 -> 어떤 lambda 를 줘야 하는가.
# joint 계열(MFKmL)만 "joint" lambda 를 쓰고, 나머지는 "sumvars" lambda 를 쓴다.
# (DTW/Euclid 는 lambda 를 아예 안 쓰므로 무엇을 주든 무관하다.)
.lambda_criterion_for <- function(m) {
  if (grepl("^MFKmL", m$name)) return("joint")
  "sumvars"
}

# ---------------------------------------------------------------------
# fit 객체에서 안전하게 스칼라 / 길이-p 벡터를 꺼낸다.  (복원 2026-07)
#
# 왜 필요한가: 방법마다 $s, $w 라는 이름의 필드가 완전히 다른 것을 담는다.
#   - SFKmL-C   : $s = L1 bound(스칼라),  $w = 변수 가중치 (길이 p)
#   - tabular   : 내부 표준편차 벡터 등 (길이 34)
#   - kNN/SVM   : 아예 없음
# 검사 없이 data.frame 열에 대입하면 "replacement has N rows, data has 1" 로
# 죽는다.  아래 두 함수가 그것을 막는다.
#
# (이 정의 블록이 파일 재작성 과정에서 유실되어, 호출만 남고 정의가 사라진
#  적이 있었다.  ".safe_wvec 를 찾을 수 없습니다" 의 원인.  다시 넣는다.)
# ---------------------------------------------------------------------
.safe_scalar <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (!is.numeric(x) || length(x) != 1L) return(NA_real_)
  if (!is.finite(x)) return(NA_real_)
  as.numeric(x)
}

.safe_wvec <- function(x, p) {
  if (is.null(x) || !is.numeric(x) || length(x) != p) return(rep(NA_real_, p))
  as.numeric(x)
}


# ---------------------------------------------------------------------
# 한 (방법 x fold) 실행
# ---------------------------------------------------------------------
flc_run_one <- function(prep, m, tr, te, seed = 1L, args = list(),
                        lam_idx = NULL) {
  ds <- prep$ds

  # lambda 선택.  "cv" 이면 TRAIN 레이블만으로 fold 안에서 고른다 (누수 방지).
  # lam_idx 를 인자로 받으면 재계산하지 않는다 (fold 안의 모든 방법이 같은
  # lambda 를 공유해야 공정하고, 매 방법마다 다시 고르면 낭비다).
  if (is.null(lam_idx)) {
    lam_idx <- if (identical(prep$lambda_method, "cv"))
                 .pick_lambda_train(prep, tr, seed,
                                    criterion = .lambda_criterion_for(m)) else 1L
  }
  lam_idx <- max(1L, min(lam_idx, length(prep$lambda_grid)))

  ctx <- flc_make_ctx(ds, tr, te,
                      prep$caches[[lam_idx]],
                      prep$lambda_grid[lam_idx],
                      prep$gamma, prep$sumOrMax, prep$grid, seed)
  t0 <- Sys.time()

  miss <- m$needs[!vapply(m$needs, flc_have, TRUE)]
  if (length(miss))
    return(list(ok = FALSE,
                error = paste("missing:", paste(miss, collapse = ",")),
                seconds = 0, w = rep(NA_real_, ds$p), s = NA_real_,
                varsel = NULL))

  r <- try({
    fit <- do.call(m$fit, c(list(ds[tr], ctx), args))
    pr  <- stats::predict(fit, ds[te], ctx)
    list(fit = fit, pr = pr)
  }, silent = TRUE)

  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(r, "try-error"))
    return(list(ok = FALSE, error = trimws(as.character(r)), seconds = sec,
                w = rep(NA_real_, ds$p), s = NA_real_, varsel = NULL))

  met <- flc_metrics(ds$y[te], r$pr$class, r$pr$prob, levels(ds$y))

  # ---- 변수선택 부산물 (반드시 방어적으로 꺼낼 것) --------------------
  w <- .safe_wvec(r$fit$w, ds$p)
  s <- .safe_scalar(r$fit$s)

  vs <- NULL
  if (!is.null(ds$informative) && !all(is.na(w)))
    vs <- varsel_scores(w, ds$informative)

  list(ok = TRUE, metrics = met, w = w, s = s, varsel = vs,
       seconds = sec, error = "")
}


# 변수선택 정확도 (TPR / FPR / exact recovery)
varsel_scores <- function(w, informative, tol = 1e-8) {
  if (length(w) != length(informative)) return(NULL)
  sel <- w > tol
  sel[is.na(sel)] <- FALSE
  tp <- sum(sel & informative);  fn <- sum(!sel & informative)
  fp <- sum(sel & !informative); tn <- sum(!sel & !informative)
  list(TPR   = if (tp + fn) tp / (tp + fn) else NA_real_,
       FPR   = if (fp + tn) fp / (fp + tn) else NA_real_,
       exact = all(sel == informative),
       n_sel = sum(sel))
}


# ---------------------------------------------------------------------
# 전역 진행률 카운터.
#
# flc_benchmark() 가 총 fit 수(datasets x seeds x methods x folds)를 미리
# 계산해 세팅하고, flc_cv() 안의 runner 가 fit 하나 끝날 때마다 1씩 올린다.
# 데이터셋마다 n, p 가 달라 속도가 다르므로 ETA 는 초반에 부정확하다.
# (HighDim50 같은 뒤쪽 설정이 훨씬 느리므로 실제로는 ETA 보다 오래 걸린다.)
# ---------------------------------------------------------------------
.flc_prog <- new.env(parent = emptyenv())
.flc_prog$done  <- 0L
.flc_prog$total <- 0L
.flc_prog$t0    <- Sys.time()

flc_progress_init <- function(total) {
  .flc_prog$done  <- 0L
  .flc_prog$total <- as.integer(total)
  .flc_prog$t0    <- Sys.time()
  invisible(NULL)
}

# "  123/5605 ( 2.2%) eta 41m" 형태의 문자열을 돌려준다.
flc_progress_tick <- function() {
  .flc_prog$done <- .flc_prog$done + 1L
  d <- .flc_prog$done
  n <- .flc_prog$total
  if (n <= 0L) return(sprintf("%d", d))

  el  <- as.numeric(difftime(Sys.time(), .flc_prog$t0, units = "secs"))
  eta <- if (d > 0L) el / d * (n - d) else NA_real_

  etastr <- if (!is.finite(eta)) {
    "--"
  } else if (eta < 90) {
    sprintf("%.0fs", eta)
  } else if (eta < 5400) {
    sprintf("%.0fm", eta / 60)
  } else {
    sprintf("%.1fh", eta / 3600)
  }

  w <- nchar(as.character(n))
  sprintf("%*d/%d (%4.1f%%) eta %-5s", w, d, n, 100 * d / n, etastr)
}


# ---------------------------------------------------------------------
# 결과 한 줄의 뼈대.  모든 열을 미리 올바른 타입으로 만들어 둔다.
# (나중에 $<- 로 채울 때 타입 충돌이 나지 않도록.)
# ---------------------------------------------------------------------
.empty_row <- function(dataset, m, fold, seed) {
  data.frame(
    dataset = as.character(dataset),
    method  = as.character(m$name),
    family  = as.character(m$family),
    source  = as.character(m$source),
    needs_grid = as.logical(m$grid),
    fold = as.integer(fold),
    seed = as.integer(seed),
    acc = NA_real_, bal_acc = NA_real_, macro_prec = NA_real_,
    macro_rec = NA_real_, macro_f1 = NA_real_, auc = NA_real_,
    kappa = NA_real_, logloss = NA_real_, brier = NA_real_,
    s = NA_real_, n_sel = NA_integer_, lambda = NA_real_,
    TPR = NA_real_, FPR = NA_real_, exact = NA,
    seconds = NA_real_, ok = FALSE, error = "",
    stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------
# 데이터셋 하나 x 방법 집합 전체에 대한 K-fold CV
# ---------------------------------------------------------------------
flc_cv <- function(ds, methods = flc_method_set("core"), K = 5L, seed = 1L,
                   lambda_method = "rule", sumOrMax = "max",
                   dataset = "unnamed", parallel = FALSE, verbose = TRUE) {

  if (verbose) cat(sprintf("\n== %s  (n=%d, p=%d, K=%d classes) ==\n",
                           dataset, ds$n, ds$p, nlevels(ds$y)))

  prep <- flc_prepare(ds, lambda_method, sumOrMax, verbose = verbose)
  fold <- strat_folds(ds$y, K, seed)

  jobs <- expand.grid(mi = seq_along(methods), f = seq_len(K))

  # fold 마다 lambda 를 고른다 (TRAIN 레이블만 사용).
  #
  # *** 거리마다 따로 고른다 ***
  # joint(MFKmL) 와 sumvars(SFKmL) 는 최적 lambda 가 다를 이유가 충분하다.
  # (.pick_lambda_train 주석 참조)  같은 값을 강제하면 MFKmL 이 불공정하게
  # 불리해지고, "MFKmL 이 진다"는 결론이 구현 아티팩트가 된다.
  #
  # lambda_method:
  #   "cv"    inner-CV 분류손실로 선택 (기존)
  #   "align" 거리행렬 target alignment 로 선택 (층위 B, supervised 거리파라미터)
  #   그 외   fold 내 선택 안 함 (rule/gap 은 전역)
  .lam_method <- if (prep$lambda_method %in% c("cv", "align"))
                   prep$lambda_method else NA_character_
  lam_sum <- vapply(seq_len(K), function(f) {
    if (is.na(.lam_method)) return(1L)
    .pick_lambda_train(prep, which(fold != f), seed = seed * 100L + f,
                       criterion = "sumvars", method = .lam_method)
  }, 1L)

  lam_joint <- vapply(seq_len(K), function(f) {
    if (is.na(.lam_method)) return(1L)
    .pick_lambda_train(prep, which(fold != f), seed = seed * 100L + f,
                       criterion = "joint", method = .lam_method)
  }, 1L)

  if (verbose && is.na(.lam_method))
    invisible(NULL)
  else if (verbose)
    cat(sprintf("  [lambda:%s] sumvars: %s   joint: %s\n", .lam_method,
                paste(sprintf("%.3g", prep$lambda_grid[lam_sum]), collapse = ","),
                paste(sprintf("%.3g", prep$lambda_grid[lam_joint]), collapse = ",")))

  lam_for <- function(m, f) {
    if (identical(.lambda_criterion_for(m), "joint")) lam_joint[f] else lam_sum[f]
  }

  runner <- function(j) {
    m  <- methods[[jobs$mi[j]]]
    f  <- jobs$f[j]
    tr <- which(fold != f)
    te <- which(fold == f)

    li <- lam_for(m, f)
    r  <- flc_run_one(prep, m, tr, te, seed = seed * 100L + f, args = m$args,
                      lam_idx = li)

    row <- .empty_row(dataset, m, f, seed)
    row$lambda  <- .safe_scalar(prep$lambda_grid[li])
    row$seconds <- .safe_scalar(r$seconds)
    row$ok      <- isTRUE(r$ok)
    row$error   <- substr(as.character(r$error %||% ""), 1, 300)

    if (isTRUE(r$ok)) {
      mr <- flc_metrics_row(r$metrics)
      for (nm in names(mr)) row[[nm]] <- .safe_scalar(mr[[nm]])
      row$s     <- .safe_scalar(r$s)
      row$n_sel <- as.integer(sum(r$w > 1e-8, na.rm = TRUE))
      if (!is.null(r$varsel)) {
        row$TPR   <- .safe_scalar(r$varsel$TPR)
        row$FPR   <- .safe_scalar(r$varsel$FPR)
        row$exact <- isTRUE(r$varsel$exact)
      }
    }

    list(row = row, w = r$w, method = m$name, fold = f)
  }

  res <- if (parallel && flc_have("foreach") && !is.null(.flc$cl)) {
    foreach::foreach(j = seq_len(nrow(jobs)),
                     .errorhandling = "pass") %dopar% runner(j)
  } else {
    lapply(seq_len(nrow(jobs)), function(j) {
      r <- runner(j)
      if (verbose) {
        cat(sprintf("  %s  %-26s f%d  acc=%s bal=%s %5.1fs%s\n",
          flc_progress_tick(),
          substr(r$row$method, 1, 26), r$row$fold,
          ifelse(is.na(r$row$acc),     "  -  ", sprintf("%.3f", r$row$acc)),
          ifelse(is.na(r$row$bal_acc), "  -  ", sprintf("%.3f", r$row$bal_acc)),
          r$row$seconds %||% 0,
          ifelse(r$row$ok, "", paste0(" [", substr(r$row$error, 1, 45), "]"))))
      }
      r })
  }

  # foreach(.errorhandling="pass") 는 실패 시 condition 객체를 돌려준다.
  # 그런 항목은 버린다 (이미 runner 안에서 try()로 잡히므로 드물다).
  res <- Filter(function(x) is.list(x) && !is.null(x$row), res)

  rows <- do.call(rbind, lapply(res, `[[`, "row"))

  Wmat <- do.call(rbind, lapply(res, function(r) {
    data.frame(
      dataset  = dataset,
      method   = r$method,
      fold     = r$fold,
      variable = ds$varNames,
      w        = .safe_wvec(r$w, ds$p),
      informative = if (!is.null(ds$informative)) ds$informative
                    else rep(NA, ds$p),
      stringsAsFactors = FALSE)
  }))

  list(results = rows, weights = Wmat, prep = prep, fold = fold)
}


# ---------------------------------------------------------------------
# 여러 시드 x 여러 데이터셋
#
# 매 (데이터셋, 시드) 마다 CSV 를 덮어쓰며 중간 저장한다.
# 오래 걸리는 실행이 중간에 죽어도 그때까지의 결과는 남는다.
# ---------------------------------------------------------------------
flc_benchmark <- function(datasets, methods = flc_method_set("core"),
                          seeds = 1:5, K = 5L, lambda_method = "rule",
                          sumOrMax = "max", outdir = "results",
                          verbose = TRUE) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  R <- list(); W <- list()

  # 총 fit 수를 미리 알려줘야 진행률 표시가 "123/5605 (2.2%)" 로 나온다.
  flc_progress_init(length(datasets) * length(seeds) * length(methods) * K)

  for (dn in datasets) {
    for (sd in seeds) {

      ds <- try(flc_get_dataset(dn, seed = sd), silent = TRUE)
      if (inherits(ds, "try-error")) {
        warning("dataset 생성 실패: ", dn, " (seed ", sd, ")")
        next
      }

      cv <- try(flc_cv(ds, methods, K = K, seed = sd,
                       lambda_method = lambda_method, sumOrMax = sumOrMax,
                       dataset = dn, verbose = verbose), silent = TRUE)
      if (inherits(cv, "try-error")) {
        warning("CV 실패: ", dn, " (seed ", sd, ")\n", as.character(cv))
        next
      }

      R[[length(R) + 1L]] <- cv$results
      W[[length(W) + 1L]] <- cbind(cv$weights, seed = sd)

      utils::write.csv(do.call(rbind, R),
                       file.path(outdir, "results.csv"), row.names = FALSE)
      utils::write.csv(do.call(rbind, W),
                       file.path(outdir, "weights.csv"), row.names = FALSE)
    }
  }

  if (!length(R)) stop("성공한 실행이 하나도 없습니다.")

  list(results = do.call(rbind, R), weights = do.call(rbind, W))
}


# ---------------------------------------------------------------------
# 요약 표
# ---------------------------------------------------------------------
flc_summarise <- function(res, by = c("dataset", "method", "family")) {
  ok <- res[res$ok, , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  cols <- c("acc", "bal_acc", "macro_f1", "auc", "kappa", "n_sel", "seconds")
  agg <- stats::aggregate(ok[, cols],
                          by = as.list(ok[, by, drop = FALSE]),
                          FUN = function(v) mean(v, na.rm = TRUE))
  sdv <- stats::aggregate(ok[, cols],
                          by = as.list(ok[, by, drop = FALSE]),
                          FUN = function(v) stats::sd(v, na.rm = TRUE))
  names(sdv)[-seq_along(by)] <- paste0(cols, "_sd")
  out <- cbind(agg, sdv[, -seq_along(by), drop = FALSE])
  out[order(out$dataset, -out$bal_acc), ]
}


# 방법별 순위 (bake-off 스타일: 데이터셋 내 평균순위)
# Middlehurst et al. (2024) 가 쓰는 방식.  값 자체보다 순위가 안정적이다.
flc_ranks <- function(res, metric = "bal_acc") {
  ok <- res[res$ok, , drop = FALSE]
  if (!nrow(ok)) return(NULL)
  a <- stats::aggregate(ok[[metric]],
                        by = list(dataset = ok$dataset, method = ok$method),
                        FUN = mean, na.rm = TRUE)
  names(a)[3] <- "v"
  a <- a[is.finite(a$v), , drop = FALSE]
  a$rank <- stats::ave(-a$v, a$dataset,
                       FUN = function(z) rank(z, ties.method = "average"))
  r <- stats::aggregate(rank ~ method, a, mean)
  r[order(r$rank), ]
}
