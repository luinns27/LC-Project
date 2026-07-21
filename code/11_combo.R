# =====================================================================
# 11_combo.R
#
# 요청사항: "MFKmL / SFKmL의 구조를 파악하고, classification으로 전환하기 위해
#            kNN, SVM 등 대중적인 분류기 4~5가지를 모두 조합"
#
# 조합의 축(axis)을 명확히 분리한다.  이 3축 factorial이 논문의 실험 설계다.
#
#   [축 1] DISTANCE (거리)
#     joint    : MFKmL의 거리.  joint curve 위의 generalized Frechet.
#                local cost = sqrt( (lambda*dt)^2 + sum_k gamma_k^2 (x_k - y_k)^2 )
#     sumvars  : SFKmL의 거리.  d(x,y) = sum_j w_j * FD_j(x_j, y_j)
#                (변수별 1차원 Frechet의 가중합. joint Frechet이 아님!)
#     dtw      : multivariate dependent DTW  (경쟁 warping 거리)
#     euclid   : 공통격자 보간 후 pointwise L2  (KmL3d 계열의 거리)
#
#   [축 2] SPARSITY (가중치)
#     dense    : w_j = 1 (or 1/sqrt(p)).  변수 선택 없음.
#     cv       : s를 inner-CV 분류손실로 선택   <- 권장
#     gap      : s를 논문 gap statistic으로 선택 (재현용)
#     gaplab   : s를 label-permutation gap으로 선택 (변수별 p-value 부산물)
#
#   [축 3] CLASSIFIER (결정규칙)  -- 5종
#     centroid : Frechet-mean prototype 최근접   (MFKmL의 k-means 정신)
#     medoid   : 클래스 medoid 최근접            (SFKmL의 k-medoids 정신)
#     knn      : k-nearest neighbour (k는 inner-LOOCV로 선택)
#     svm      : kernel SVM.  K = exp(-gamma D^2)   (distance-substitution kernel)
#     kfda     : kernel Fisher discriminant analysis (동일 커널, 다른 규칙)
#     [+ rf   : proximity/거리기반 랜덤포레스트는 12_base_ml.R 에 별도]
#
# 왜 medoid와 centroid를 둘 다 두는가:
#   원 패키지에서 MFKmL은 Frechet MEAN, SFKmL은 MEDOID를 중심으로 쓴다.
#   이 차이가 성능차를 만드는지, 아니면 거리(joint vs sumvars)가 만드는지
#   분리하려면 두 축을 교차시켜야만 한다.  그래서 격자로 만든다.
#
# 주의: joint 거리 + w 가중은 원리적으로 캐시가 불가능하다
#       (w가 sqrt 안으로 들어가므로 선형결합으로 복원 불가).
#       따라서 sparsity는 sumvars 계열에서만 의미가 있고,
#       joint 계열은 dense로만 돌린다.  이 비대칭 자체가 논문의 논거이다:
#       "SFKmL의 sum-over-variables 구조가 sparsity를 가능케 한 이유".
# =====================================================================


# ---------------------------------------------------------------------
# 통합 분류기.  축 3개를 인자로 받는다.
# ---------------------------------------------------------------------
fit_combo <- function(ds, ctx,
                      distance = c("sumvars", "joint", "dtw", "euclid"),
                      sparsity = c("dense", "cv", "gap", "gaplab"),
                      clf      = c("medoid", "centroid", "knn", "svm", "kfda", "gp"),
                      k = NULL, C = 1, gamma_mult = 1, reg = 1e-3,
                      knn_weighted = FALSE, knn_nested = FALSE, ...) {
  distance <- match.arg(distance)
  sparsity <- match.arg(sparsity)
  clf      <- match.arg(clf)

  y <- ds$y; cls <- levels(y); p <- ds$p

  # ---- 축 2 : 가중치 w 학습 -----------------------------------------
  if (distance == "sumvars" && sparsity != "dense") {
    W <- flc_learn_weights(ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE], y,
                           s_select = sparsity, seed = ctx$seed)
    w <- W$w; s_used <- W$s; sel <- W$sel; bcss <- W$bcss
  } else {
    if (distance != "sumvars" && sparsity != "dense")
      warning("sparsity is only defined for distance='sumvars'; forcing dense")
    W <- NULL
    w <- rep(1, p); s_used <- sqrt(p); sel <- rep(TRUE, p)
    bcss <- if (distance == "sumvars")
      sfkml_objective_a(ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE], y)
      else rep(NA_real_, p)
  }

  # ---- 축 1 : train 거리행렬 (캐시에서 조립) --------------------------
  Dtr <- .train_dist(ctx, distance, w)

  obj <- list(distance = distance, sparsity = sparsity, clf = clf,
              w = w, s = s_used, sel = sel, bcss = bcss, W = W,
              y = y, levels = cls, ctx = ctx, p = p,
              varNames = ds$varNames,
              gamma_mult = gamma_mult, C = C, knn_weighted = knn_weighted,
              method = sprintf("%s|%s|%s", distance, sparsity, clf))

  # ---- 축 3 : 결정규칙 -------------------------------------------------
  if (clf == "medoid") {
    obj$medoid <- vapply(cls, function(c) {
      idx <- which(y == c)
      if (length(idx) <= 2L) return(idx[1])
      sub <- Dtr[idx, idx, drop = FALSE]
      diag(sub) <- 0                       # 자기 자신 제외 (정의상)
      idx[which.min(colSums(sub))]
    }, 1L)

  } else if (clf == "centroid") {
    # 거리에 맞는 barycentre를 쓴다.
    obj$protos <- switch(distance,
      joint   = lapply(cls, function(c)
                  flc_frechet_mean(ds, which(y == c), w = rep(1, p),
                                   lambda = ctx$lambda, sumOrMax = ctx$sumOrMax)),
      # WELL-DEFINEDNESS FIX: sumvars 거리의 barycentre 는 변수별 1차원
      # Frechet mean 이며 w 에 의존하지 않는다 (02_frechet.R 유도 참조).
      sumvars = lapply(cls, function(c)
                  flc_frechet_mean_pervar(ds, which(y == c),
                                          lambda = ctx$lambda,
                                          sumOrMax = ctx$sumOrMax,
                                          grid = ctx$grid)),
      dtw     = lapply(cls, function(c)
                  flc_dba_mean(ds, which(y == c), w = rep(1, p),
                               nIter = 5L, grid = ctx$grid)),
      euclid  = { R <- flc_regularize(ds, ctx$grid)
                  obj$grid <- R$time
                  lapply(cls, function(c)
                    apply(R$traj[y == c, , , drop = FALSE], c(2, 3), mean)) })

  } else if (clf == "knn") {
    # k 선택.  기본은 원래 방식(train resubstitution LOO) -> 기존 결과 재현.
    # knn_nested=TRUE 이면 nested inner-CV (덜 낙관적).  옵션으로만 켠다.
    if (is.null(k)) {
      kg <- c(1, 3, 5, 7, 9, 11)
      kg <- kg[kg < min(table(y))]; if (!length(kg)) kg <- 1
      if (isTRUE(knn_nested)) {
        k <- .knn_select_k_nested(Dtr, y, kg, cls,
                                  weighted = knn_weighted, seed = ctx$seed)
      } else {
        # 원래 방식: train 자기예측 LOO (diag=Inf).
        Dl  <- Dtr; diag(Dl) <- Inf
        ord <- t(apply(Dl, 1, order))
        vote <- if (isTRUE(knn_weighted)) {
          function(kk) vapply(seq_len(nrow(Dl)), function(i) {
            nn <- ord[i, seq_len(kk)]; wt <- 1 / (Dl[i, nn] + 1e-8)
            names(which.max(tapply(wt, factor(y[nn], levels = cls), sum)))
          }, "")
        } else {
          function(kk) apply(ord[, seq_len(kk), drop = FALSE], 1, function(ix)
            names(which.max(table(factor(y[ix], levels = cls)))))
        }
        a <- vapply(kg, function(kk) mean(vote(kk) == as.character(y)), 0)
        k <- kg[which.max(a)]
      }
    }
    obj$k <- k
    obj$knn_weighted <- knn_weighted

  } else if (clf == "svm") {
    if (!flc_have("kernlab")) stop("kernlab required for clf='svm'")
    # distance-substitution kernel: K(x,y) = exp(-g * d(x,y)^2)
    # d가 metric이면 (논문 Appendix A: generalized Frechet는 metric임을 증명)
    # 이 커널은 실용상 잘 동작한다.  PD가 보장되지 않을 수 있으므로
    # 대각 jitter를 추가한다.
    med <- stats::median(Dtr[Dtr > 0]); if (!is.finite(med) || med <= 0) med <- 1
    g <- gamma_mult / (2 * med^2)
    K <- exp(-g * Dtr^2)
    K <- K + diag(1e-8, nrow(K))
    obj$gamma_k <- g
    obj$svm <- kernlab::ksvm(kernlab::as.kernelMatrix(K), y, type = "C-svc",
                             C = C, prob.model = TRUE)
    obj$sv <- kernlab::SVindex(obj$svm)

  } else if (clf == "kfda") {
    # kernel Fisher discriminant.  동일한 distance-substitution kernel 위에서
    # SVM과 다른 결정규칙을 준다 -> "거리가 일하는가, 규칙이 일하는가" 분리용.
    med <- stats::median(Dtr[Dtr > 0]); if (!is.finite(med) || med <= 0) med <- 1
    g <- gamma_mult / (2 * med^2)
    K <- exp(-g * Dtr^2)
    obj$gamma_k <- g
    obj$kfda <- .kfda_fit(K, y, reg = reg)

  } else if (clf == "gp") {
    # ---- (2) 커널 로지스틱 회귀 (GP 분류기의 라플라스 근사와 동형) --------
    # distance-substitution kernel 위에서 kernel (multinomial) logistic
    # regression 을 IRLS 로 적합.  SVM/kFDA 와 같은 커널을 쓰되, 확률적
    # 결정경계를 준다.  GP 분류기(라플라스 근사)는 릿지 벌점을 가진 커널
    # 로지스틱 회귀와 사실상 동일한 예측을 준다 -- 그 형태로 구현한다.
    med <- stats::median(Dtr[Dtr > 0]); if (!is.finite(med) || med <= 0) med <- 1
    g <- gamma_mult / (2 * med^2)
    K <- exp(-g * Dtr^2)
    obj$gamma_k <- g
    obj$gp <- .klr_fit(K, y, lambda = reg)
  }

  structure(obj, class = "flcCombo")
}


predict.flcCombo <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels; clf <- object$clf

  if (clf == "medoid") {
    D <- .cross_dist(ctx, object$distance, object$w)[, object$medoid, drop = FALSE]
    colnames(D) <- cls
    return(list(prob = dist_to_prob(D),
                class = factor(cls[max.col(-D, ties.method = "first")], levels = cls),
                D = D))
  }

  if (clf == "centroid") {
    D <- switch(object$distance,
      joint = vapply(object$protos, function(mu)
                flc_cross_joint(ds_te,
                  flc_curve_as_ds(mu, object$p, object$varNames),
                  w = rep(1, object$p), lambda = ctx$lambda,
                  sumOrMax = ctx$sumOrMax)[, 1], numeric(ds_te$n)),
      sumvars = vapply(object$protos, function(mu) {
                m <- flc_curve_as_ds(mu, object$p, object$varNames)
                acc <- numeric(ds_te$n)
                for (j in which(object$w > 1e-12)) {
                  A <- flc_var_list(ds_te, j); B <- flc_var_list(m, j)
                  acc <- acc + object$w[j] *
                    frechet_crossmat_cpp(A$X, A$t, B$X, B$t, 1, ctx$lambda,
                      if (ctx$sumOrMax == "sum") 1L else 0L,
                      flc_kernel_threads())[, 1]
                }
                acc }, numeric(ds_te$n)),
      dtw = vapply(object$protos, function(mu)
              flc_cross_dtw(ds_te,
                flc_curve_as_ds(mu, object$p, object$varNames),
                rep(1, object$p))[, 1], numeric(ds_te$n)),
      euclid = { R <- flc_regularize(ds_te, object$grid)
                 vapply(object$protos, function(m)
                   apply(R$traj, 1, function(x) sqrt(sum((x - m)^2))),
                   numeric(ds_te$n)) })
    D <- matrix(D, ds_te$n, length(cls), dimnames = list(NULL, cls))
    return(list(prob = dist_to_prob(D),
                class = factor(cls[max.col(-D, ties.method = "first")], levels = cls),
                D = D))
  }

  D <- .cross_dist(ctx, object$distance, object$w)     # n_te x n_tr

  if (clf == "knn") {
    # (3) 거리 가중 투표 옵션.  knn_weighted=TRUE 면 이웃을 1/(d+eps) 로 가중.
    if (isTRUE(object$knn_weighted)) {
      P <- .knn_vote_weighted(D, object$y, object$k, cls)
    } else {
      # BUG FIX: knn_vote() 사용 (01b_helpers.R 주석 참조).
      P <- knn_vote(D, object$y, object$k, cls)
    }
    return(list(prob = P,
                class = factor(cls[max.col(P, ties.method = "first")], levels = cls),
                D = D))
  }

  if (clf == "svm") {
    Kte <- exp(-object$gamma_k * D^2)[, object$sv, drop = FALSE]
    P <- as.matrix(kernlab::predict(object$svm, kernlab::as.kernelMatrix(Kte),
                                    type = "probabilities"))[, cls, drop = FALSE]
    return(list(prob = P,
                class = factor(cls[max.col(P, ties.method = "first")], levels = cls),
                D = D))
  }

  if (clf == "kfda") {
    Kte <- exp(-object$gamma_k * D^2)
    return(.kfda_predict(object$kfda, Kte, cls, D))
  }

  if (clf == "gp") {
    Kte <- exp(-object$gamma_k * D^2)      # n_te x n_tr
    P <- .klr_predict(object$gp, Kte, cls)
    return(list(prob = P,
                class = factor(cls[max.col(P, ties.method = "first")], levels = cls),
                D = D))
  }
}


# =====================================================================
# (4) k 선택: nested inner-CV
#   TRAIN 거리행렬만으로, inner K-fold held-out 정확도를 최대화하는 k.
#   Dtr = ctx$tr 부분행렬 (n_tr x n_tr).  test 안 들어옴 -> 누출 없음.
# =====================================================================
.knn_select_k_nested <- function(Dtr, y, kg, cls, weighted = FALSE,
                                 innerK = 5L, seed = 1L) {
  y  <- factor(y, levels = cls); n <- nrow(Dtr)
  nf <- min(innerK, min(table(y)))
  if (nf < 2L || length(kg) == 1L) return(kg[1])
  fold <- strat_folds(y, nf, seed)

  voter <- if (weighted) .knn_vote_weighted else knn_vote
  acc <- vapply(kg, function(kk) {
    a <- vapply(sort(unique(fold)), function(f) {
      i1 <- which(fold != f); i2 <- which(fold == f)
      D12 <- Dtr[i2, i1, drop = FALSE]         # val x train
      P <- voter(D12, y[i1], kk, cls)
      pr <- cls[max.col(P, ties.method = "first")]
      mean(pr == as.character(y[i2]))
    }, 0)
    mean(a, na.rm = TRUE)
  }, 0)
  kg[which.max(acc)]
}


# =====================================================================
# (3) 거리 가중 kNN 투표
#   균등 투표 대신 이웃을 1/(d+eps) 로 가중.  가까운 이웃이 더 큰 표를 준다.
#   D: n_te x n_tr 거리행렬,  반환: n_te x K 확률행렬.
# =====================================================================
.knn_vote_weighted <- function(D, ytr, k, cls = levels(ytr)) {
  D   <- as.matrix(D); n <- nrow(D); K <- length(cls)
  ytr <- factor(ytr, levels = cls)
  k   <- max(1L, min(as.integer(k), ncol(D)))
  eps <- 1e-8
  P <- matrix(0, n, K, dimnames = list(NULL, cls))
  for (i in seq_len(n)) {
    di <- D[i, ]; di[!is.finite(di)] <- Inf
    nn <- order(di)[seq_len(k)]
    wt <- 1 / (di[nn] + eps)
    for (a in seq_along(nn))
      P[i, as.integer(ytr[nn[a]])] <- P[i, as.integer(ytr[nn[a]])] + wt[a]
    s <- sum(P[i, ]); if (s > 0) P[i, ] <- P[i, ] / s else P[i, ] <- 1 / K
  }
  P
}


# =====================================================================
# (2) kernel (multinomial) logistic regression  -- GP 분류기의 실용 형태
#
#   모형:  f_c = K a_c,   p_c = softmax_c(f),   벌점 lambda * a_c' K a_c
#   IRLS/뉴턴 스텝으로 계수 a (n x K) 를 구한다.  라플라스 근사 GP 분류기와
#   같은 예측족.  이진/다중 모두 처리.
#
#   K   : n x n 훈련 커널 (PSD 가정; 대각 jitter 로 안정화)
#   반환: list(a = n x K 계수, cls, K 참조 불필요 -- 예측은 Kte 로)
# =====================================================================
.klr_fit <- function(K, y, lambda = 1e-2, iters = 25L) {
  y  <- factor(y); cls <- levels(y); n <- nrow(K); Kc <- length(cls)
  Y  <- model.matrix(~ y - 1); colnames(Y) <- cls    # n x K one-hot
  Kj <- K + diag(1e-6, n)

  A <- matrix(0, n, Kc)                               # 계수
  for (it in seq_len(iters)) {
    F  <- Kj %*% A                                    # n x K logits
    F  <- F - apply(F, 1, max)
    Pm <- exp(F); Pm <- Pm / rowSums(Pm)
    # 뉴턴 스텝 (클래스별 근사: 대각 가중 W_c = p_c(1-p_c))
    Anew <- A
    for (c in seq_len(Kc)) {
      w  <- pmax(Pm[, c] * (1 - Pm[, c]), 1e-4)
      z  <- F[, c] + (Y[, c] - Pm[, c]) / w           # 작업 반응
      # (diag(w) K + lambda I) a_c = diag(w) z  형태의 릿지 해
      Kw <- Kj * w                                    # 행 스케일
      sol <- try(solve(Kw + lambda * diag(n), w * z), silent = TRUE)
      if (inherits(sol, "try-error"))
        sol <- solve(Kw + (lambda + 1e-2) * diag(n), w * z)
      Anew[, c] <- sol
    }
    if (max(abs(Anew - A)) < 1e-6) { A <- Anew; break }
    A <- Anew
  }
  list(A = A, cls = cls)
}

.klr_predict <- function(fit, Kte, cls) {
  F <- Kte %*% fit$A                                  # n_te x K logits
  F <- F - apply(F, 1, max)
  P <- exp(F); P <- P / rowSums(P)
  colnames(P) <- fit$cls
  P[, cls, drop = FALSE]
}


# =====================================================================
# kernel Fisher discriminant analysis (multi-class, one-vs-rest projection
# + LDA in the projected space).  Mika et al. (1999); Baudat & Anouar (2000).
# =====================================================================
.kfda_fit <- function(K, y, reg = 1e-3) {
  y <- factor(y); n <- nrow(K); cls <- levels(y)
  # centre in feature space
  one <- matrix(1 / n, n, n)
  Kc <- K - one %*% K - K %*% one + one %*% K %*% one

  M_star <- colMeans(Kc)
  M <- vapply(cls, function(c) colMeans(Kc[y == c, , drop = FALSE]), numeric(n))
  B <- matrix(0, n, n)                      # between
  W <- matrix(0, n, n)                      # within
  for (ci in seq_along(cls)) {
    idx <- which(y == cls[ci]); nk <- length(idx)
    d <- M[, ci] - M_star
    B <- B + nk * tcrossprod(d)
    Kk <- Kc[idx, , drop = FALSE]
    Hk <- Kk - matrix(M[, ci], nk, n, byrow = TRUE)
    W <- W + crossprod(Hk)
  }
  W <- W + diag(reg * mean(diag(W)) + 1e-8, n)
  ev <- try(eigen(solve(W, B), symmetric = FALSE), silent = TRUE)
  if (inherits(ev, "try-error")) stop("kfda: eigen failed")
  nd <- min(length(cls) - 1L, n - 1L)
  A <- Re(ev$vectors[, seq_len(nd), drop = FALSE])
  Ztr <- Kc %*% A
  # cen 은 반드시 (K x nd) 행렬이어야 한다.  BUG FIX (2026-07):
  #   nd = 1 (=2클래스) 일 때 vapply(..., numeric(1)) 은 길이-K 벡터를 주고
  #   t() 를 하면 (1 x K) 가 되어, 뒤에서 cen[ci, ] 가 범위를 벗어난다.
  #   행렬을 명시적으로 K x nd 로 만들어 이 문제를 없앤다.
  cen <- matrix(0, length(cls), nd, dimnames = list(cls, NULL))
  for (ci in seq_along(cls))
    cen[ci, ] <- colMeans(Ztr[y == cls[ci], , drop = FALSE])
  list(A = A, cen = cen, Kmean_col = colMeans(K), Kmean_all = mean(K),
       levels = cls, n = n, nd = nd)
}
.kfda_predict <- function(fit, Kte, cls, D) {
  m <- nrow(Kte)
  Kc <- Kte - matrix(fit$Kmean_col, m, fit$n, byrow = TRUE) -
        matrix(rowMeans(Kte), m, fit$n) + fit$Kmean_all
  Z <- Kc %*% fit$A
  Dz <- vapply(seq_along(cls), function(ci)
    sqrt(rowSums((Z - matrix(fit$cen[ci, ], m, fit$nd, byrow = TRUE))^2)),
    numeric(m))
  Dz <- matrix(Dz, m, length(cls), dimnames = list(NULL, cls))
  list(prob = dist_to_prob(Dz),
       class = factor(cls[max.col(-Dz, ties.method = "first")], levels = cls),
       D = D)
}


# =====================================================================
# 격자 정의.  run_master.R 이 이 표를 그대로 순회한다.
# =====================================================================
flc_combo_grid <- function(full = TRUE) {
  g <- expand.grid(
    distance = c("sumvars", "joint", "dtw", "euclid"),
    sparsity = c("dense", "cv", "gap"),
    clf      = c("medoid", "centroid", "knn", "svm", "kfda"),  # gp 제외 (성능 낮음)
    stringsAsFactors = FALSE)
  # sparsity 는 sumvars 에서만 정의된다 (joint 는 w 가 sqrt 안으로 들어가
  # 선형계획이 되지 않으므로 원리적으로 희소화가 불가능하다).
  g <- g[!(g$distance != "sumvars" & g$sparsity != "dense"), ]
  if (!full) g <- g[g$sparsity %in% c("dense", "cv"), ]

  g$name <- flc_method_label(g$distance, g$sparsity, g$clf)
  rownames(g) <- NULL
  g
}


# =====================================================================
# .tune_kernel_params  --  svm/kfda 의 (C, gamma) 를 fold 안에서 선택.
#
# 누출 방지: ctx$tr (TRAIN) 만 쓴다.  TRAIN 을 다시 inner K-fold 로 쪼개
# (C, gamma_mult) 격자에서 bal_acc 를 최대화하는 조합을 고른다.
#
#   C          : SVM 여유(soft-margin).  {0.25, 1, 4, 16}
#   gamma_mult : 커널 폭 배율.  median heuristic(=1) 주변 {0.25,0.5,1,2,4}
#                실제 gamma = gamma_mult / (2 * median(D)^2)
#
# 거리행렬은 ctx 캐시에서 이미 만들어져 있으므로, 튜닝 비용은 커널 적합뿐이다.
# =====================================================================
.tune_kernel_params <- function(ds, ctx, distance, sparsity, clf,
                                seed = 1L, innerK = 3L) {
  y   <- ds$y
  cls <- levels(y)
  tr  <- ctx$tr
  ytr <- y[tr]

  Cg <- c(0.25, 1, 4, 16)
  Gg <- c(0.25, 0.5, 1, 2, 4)
  if (clf == "kfda") Cg <- 1        # kfda 는 C 무관 -> gamma 만 튜닝

  nf   <- min(innerK, min(table(ytr)))
  if (nf < 2L) return(list(C = 1, gamma_mult = 1))   # 너무 작으면 기본값
  fold <- strat_folds(ytr, nf, seed)

  # TRAIN 안에서 거리행렬을 만들려면 combo 를 직접 적합하는 게 가장 간단.
  # inner fold 마다 fit_combo(tune=FALSE) 를 (C,gamma) 조합으로 돌린다.
  score <- function(Cv, Gv) {
    accs <- vapply(sort(unique(fold)), function(f) {
      i1 <- tr[fold != f]; i2 <- tr[fold == f]
      ctx_in <- ctx; ctx_in$tr <- i1; ctx_in$te <- i2
      fit <- try(fit_combo(ds, ctx_in, distance = distance, sparsity = sparsity,
                           clf = clf, C = Cv, gamma_mult = Gv, tune = FALSE),
                 silent = TRUE)
      if (inherits(fit, "try-error")) return(NA_real_)
      pr <- try(stats::predict(fit, ds[i2], ctx_in), silent = TRUE)
      if (inherits(pr, "try-error")) return(NA_real_)
      mean(pr$class == y[i2])
    }, 0)
    mean(accs, na.rm = TRUE)
  }

  grid <- expand.grid(C = Cg, G = Gg)
  sc <- vapply(seq_len(nrow(grid)),
               function(i) score(grid$C[i], grid$G[i]), 0)
  if (all(!is.finite(sc))) return(list(C = 1, gamma_mult = 1))
  b <- which.max(sc)
  list(C = grid$C[b], gamma_mult = grid$G[b])
}



#
# [왜 바꿨나]
#   구표기는 "sumvars|cv|medoid" 같은 파이프 형식이었다.  문제가 둘 있었다.
#
#   (1) sumvars 가 SFKmL 의 거리라는 것이 이름에서 안 보인다.
#       그래서 "SFKmL-C 가 sumvars|cv|knn 보다 낮다" 처럼, **같은 계열인데
#       다른 방법인 것처럼** 읽히는 곡해가 생겼다.
#       -> 사실은 둘 다 SFKmL-C 이고, 결정규칙만 다르다.
#
#   (2) dense 가 무슨 축인지 안 보인다.  dense 는 "cross-validation 의 반대"가
#       아니라 "sparse 의 반대"다.  즉 희소성 축의 값이다:
#           dense  : 희소성 없음.  w_j = 1 전부 (s = sqrt(p))
#           cv     : s 를 inner-CV 분류손실로 선택   (제안)
#           gap    : s 를 논문의 gap statistic 으로 선택
#           gaplab : s 를 label-permutation gap 으로 선택
#
# [새 표기]
#   SFKmL-C + kNN           SFKmL 거리 + 희소가중(cv) + kNN 규칙
#   SFKmL-C + SVM
#   SFKmL-C + medoid        <- 원 논문(K-medoids)에 가장 충실한 버전
#   SFKmL-C(dense) + kNN    <- 희소성 없는 대조군
#   SFKmL-C(gap) + medoid   <- s 를 논문 방식으로 고른 버전
#   MFKmL-C + kNN           MFKmL 거리(joint curve) + kNN 규칙
#   MFKmL-C + mean          <- 원 논문(K-means)에 가장 충실한 버전
#   DTW + kNN,  Euclid + kNN, ...
#
# [중요한 개념 정리 -- 논문에 반드시 쓸 것]
#   SFKmL-C 는 "하나의 방법"이 아니라 **계열**이다.
#   거리(변수별 Frechet 의 가중합)와 희소 가중치를 지도학습으로 옮긴 것이
#   기여이고, **결정규칙은 그 위에서 자유롭게 고를 수 있다.**
#   medoid 를 쓴 것은 원 논문이 K-medoids 였기 때문일 뿐, 필연이 아니다.
#   실제로 kNN/SVM 이 medoid 보다 낫다 -- 이는 Kang et al. 이 이미 지적한
#   medoid 의 약점(짧은 궤적이 중심으로 뽑히는 문제)이, 분류에서는 보고용
#   그림이 아니라 **결정규칙 자체**로 작동하기 때문이다.
# =====================================================================
flc_method_label <- function(distance, sparsity, clf) {
  base <- c(sumvars = "SFKmL-C", joint = "MFKmL-C",
            dtw = "DTW", euclid = "Euclid")[distance]

  # 희소성: cv 가 기본(제안)이므로 이름에 안 붙인다.  나머지는 괄호로 표시.
  tag <- ifelse(sparsity == "cv", "",
         ifelse(sparsity == "dense", "(dense)", paste0("(", sparsity, ")")))
  # joint/dtw/euclid 는 희소성이 정의되지 않으므로 dense 태그도 생략.
  tag[distance != "sumvars"] <- ""

  rule <- c(medoid = "medoid", centroid = "centroid",
            knn = "kNN", svm = "SVM", kfda = "kFDA", gp = "GP")[clf]
  # 각 거리의 "자연스러운" 중심 규칙에는 원 논문 용어를 쓴다.
  rule[clf == "centroid" & distance == "joint"]   <- "mean"      # MFKmL = K-means
  rule[clf == "centroid" & distance == "sumvars"] <- "mean"
  rule[clf == "centroid" & distance == "dtw"]     <- "DBA"       # Petitjean
  rule[clf == "centroid" & distance == "euclid"]  <- "centroid"  # KmL3d

  paste0(base, tag, " + ", rule)
}
