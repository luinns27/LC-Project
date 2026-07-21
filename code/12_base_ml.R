# =====================================================================
# 12_base_ml.R -- 나머지 벤치마크. 제공된 txt 논문에서 R로 구현/호출 가능한
#                 방법을 전부 끌어온다.
#
# 출처별 정리
# -------------------------------------------------------------------
#  Bagnall et al. (2017) "The Great Time Series Classification Bake Off"
#  Middlehurst, Schaefer & Bagnall (2024) "Bake off redux" (DMKD)
#     -> 이 두 논문이 지정한 "반드시 이겨야 하는" 표준 baseline:
#        1NN-DTW (08_base_dist.R), ROCKET/MultiROCKET, catch22+RF,
#        그리고 flat feature + RF/SVM.
#     -> HIVE-COTE v2 / Hydra 는 R 구현이 없다.  대신 ROCKET(=Hydra의 조상,
#        bake-off redux에서 accuracy 상위권)과 catch22-RF로 대표한다.
#        md 파일에 이 대체 사유를 명시한다.
#
#  Preda, Saporta & Leveder (2007) "PLS classification of functional data"
#     -> PLS-DA on functional data.  다변량 확장 = MFPLS
#        (Statistics and Computing 2024, 34:5).
#
#  Fukuda et al. (2023) "Multivariate functional subspace classification"
#     -> 클래스별 subspace + 정준각(canonical angle) 기반 분류.
#
#  fda.usc (Febrero-Bande & Oviedo de la Fuente)
#     -> classif.knn / classif.kernel / classif.glm / classif.depth
#
# 공통 처리: 이 방법들은 전부 "공통 격자 + 결측 없음"을 요구한다.
#            따라서 flc_regularize()를 통과해야 하고, 그 사실 자체를
#            (interpolated = TRUE) 기록해서 결과표에 남긴다.
#            비동기/불규칙 setting에서 이들이 손해를 보는 것은
#            버그가 아니라 **측정하려는 대상**이다.
# =====================================================================

# ---------------------------------------------------------------------
# 공통: [n x T x p] -> [n x (T*p)] flat matrix
# ---------------------------------------------------------------------
.flatten <- function(traj) matrix(traj, nrow = dim(traj)[1])

# ---------------------------------------------------------------------
# summary-statistic features.  선행연구에서 이 baseline이 논문의 setting을
# 1.000으로 풀어버렸다 -> "왜 Frechet인가"가 성립하지 않는다는 증거.
# 새 setting(Shape/Phase)에서 이것이 무너지는지 보는 것이 실험의 핵심.
# ---------------------------------------------------------------------
.tsfeat <- function(traj) {
  n <- dim(traj)[1]; Tn <- dim(traj)[2]; p <- dim(traj)[3]
  out <- lapply(seq_len(p), function(k) {
    X <- traj[, , k, drop = FALSE][, , 1]
    dX <- t(apply(X, 1, diff))
    if (Tn == 2L) dX <- matrix(dX, n, 1)
    cbind(
      mean  = rowMeans(X),
      sd    = apply(X, 1, stats::sd),
      min   = apply(X, 1, min),
      max   = apply(X, 1, max),
      med   = apply(X, 1, stats::median),
      iqr   = apply(X, 1, stats::IQR),
      q10   = apply(X, 1, stats::quantile, .10),
      q90   = apply(X, 1, stats::quantile, .90),
      slope = apply(X, 1, function(v) stats::coef(stats::lsfit(seq_len(Tn), v))[2]),
      auc   = rowSums(X),
      argmax= apply(X, 1, which.max) / Tn,
      argmin= apply(X, 1, which.min) / Tn,
      dmean = rowMeans(dX),
      dsd   = apply(dX, 1, stats::sd),
      dabs  = rowMeans(abs(dX)),
      ncross= apply(X, 1, function(v) sum(diff(sign(v - mean(v))) != 0)),
      acf1  = apply(X, 1, function(v) {
                if (stats::sd(v) < 1e-12) 0 else stats::acf(v, 1, plot = FALSE)$acf[2] })
    )
  })
  M <- do.call(cbind, out)
  colnames(M) <- paste0("v", rep(seq_len(p), each = ncol(out[[1]])), "_",
                        rep(colnames(out[[1]]), p))
  M
}

# generic "regularize -> features -> tabular learner" wrapper
.fit_tabular <- function(ds, ctx, featfun, learner, ...) {
  R <- flc_regularize(ds, ctx$grid)
  F <- featfun(R$traj)
  ctr <- colMeans(F); sdv <- apply(F, 2, stats::sd); sdv[sdv < 1e-12] <- 1
  Fs <- scale(F, ctr, sdv)
  Fs[!is.finite(Fs)] <- 0
  mdl <- learner(Fs, ds$y, ...)
  structure(list(mdl = mdl, featfun = featfun, ctr = ctr, sdv = sdv,
                 grid = R$time, levels = levels(ds$y), ctx = ctx),
            class = "flcTabular")
}
predict.flcTabular <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  F <- object$featfun(R$traj)
  Fs <- scale(F, object$ctr, object$sdv); Fs[!is.finite(Fs)] <- 0
  P <- object$mdl$predict(Fs)
  cls <- object$levels
  P <- as.matrix(P)[, cls, drop = FALSE]
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls),
       D = -P)
}

# ---- tabular learners -------------------------------------------------
.learn_rf <- function(X, y, ntree = 500L, ...) {
  if (!flc_have("randomForest")) stop("randomForest not installed")
  m <- randomForest::randomForest(X, y, ntree = ntree, ...)
  list(predict = function(Z) stats::predict(m, Z, type = "prob"),
       imp = randomForest::importance(m))
}
.learn_svm <- function(X, y, ...) {
  if (!flc_have("e1071")) stop("e1071 not installed")
  m <- e1071::svm(X, y, kernel = "radial", probability = TRUE, ...)
  list(predict = function(Z) {
    pr <- stats::predict(m, Z, probability = TRUE)
    attr(pr, "probabilities") })
}
.learn_knn <- function(X, y, k = 5L, ...) {
  list(predict = function(Z) {
    D <- as.matrix(stats::dist(rbind(Z, X)))[seq_len(nrow(Z)),
                                             nrow(Z) + seq_len(nrow(X)), drop = FALSE]
    cls <- levels(y)
    ord <- t(apply(D, 1, order))
    P <- t(apply(ord[, seq_len(k), drop = FALSE], 1, function(ix)
      as.numeric(table(factor(y[ix], levels = cls))) / k))
    colnames(P) <- cls; P })
}
.learn_glmnet <- function(X, y, alpha = 0, ...) {
  if (!flc_have("glmnet")) stop("glmnet not installed")
  cv <- glmnet::cv.glmnet(X, y, family = "multinomial", alpha = alpha, nfolds = 5)
  list(predict = function(Z) {
    P <- stats::predict(cv, Z, s = "lambda.min", type = "response")
    matrix(P, nrow(Z), dim(P)[2], dimnames = list(NULL, dimnames(P)[[2]])) },
    cv = cv)
}
.learn_lda <- function(X, y, ...) {
  if (!flc_have("MASS")) stop("MASS not installed")
  keep <- apply(X, 2, stats::sd) > 1e-8
  m <- MASS::lda(X[, keep, drop = FALSE], y)
  list(predict = function(Z) stats::predict(m, Z[, keep, drop = FALSE])$posterior)
}

# =====================================================================
# !! 아래 방법들은 2026-07 에 벤치마크에서 제외되었다 (코드는 남겨둠).
#
# 제외 사유: "궤적 구조를 버리는" 방법이라 본 연구가 보고자 하는 대상이
#            아니다.  곡선을 스칼라 요약이나 순서 없는 벡터로 붕괴시킨다.
#
#   tsfeat-RF / tsfeat-SVM : 17개 요약통계(mean, sd, argmax, ...)로 붕괴.
#       *** 특히 argmax/argmin 이 문제였다.  gen_phase() 는 클래스를 bump
#       위치로 정의하므로 argmax 가 사실상 레이블을 직접 인코딩한다.
#       "요약통계가 강하다"가 아니라 "특성 하나가 레이블을 누설한다"에
#       가깝다.  이 방법이 19개 데이터셋 전부에서 1위였던 이유. ***
#   catch22-RF             : 22개 요약통계.
#   flat-RF/SVM/ridge      : 값을 이어붙인 벡터 (모형이 순서를 활용 안 함).
#   ROCKET-ridge           : 곡선을 버리진 않으나 함께 제외.
#                            rocket_transform_cpp 버그는 미수정 상태.
#
# 되살리려면 13_registry.R 의 FLC_METHODS_EXCLUDED 를 FLC_METHODS 에
# 합치면 된다.  함수 정의는 아래에 그대로 살아 있다.
# =====================================================================

# ---- registered baselines ---------------------------------------------
fit_tsfeat_rf   <- function(ds, ctx, ...) .fit_tabular(ds, ctx, .tsfeat, .learn_rf)
fit_tsfeat_svm  <- function(ds, ctx, ...) .fit_tabular(ds, ctx, .tsfeat, .learn_svm)
fit_flat_rf     <- function(ds, ctx, ...) .fit_tabular(ds, ctx, .flatten, .learn_rf)
fit_flat_svm    <- function(ds, ctx, ...) .fit_tabular(ds, ctx, .flatten, .learn_svm)
fit_flat_ridge  <- function(ds, ctx, ...) .fit_tabular(ds, ctx, .flatten, .learn_glmnet)

# catch22 (Lubba et al.; bake-off redux의 표준 feature set)
fit_catch22_rf <- function(ds, ctx, ...) {
  if (!flc_have("Rcatch22")) stop("Rcatch22 not installed")
  ff <- function(traj) {
    n <- dim(traj)[1]; p <- dim(traj)[3]
    M <- lapply(seq_len(p), function(k)
      t(vapply(seq_len(n), function(i)
        Rcatch22::catch22_all(traj[i, , k])$values, numeric(22))))
    M <- do.call(cbind, M)
    M[!is.finite(M)] <- 0
    M
  }
  .fit_tabular(ds, ctx, ff, .learn_rf)
}

# ---------------------------------------------------------------------
# ROCKET (Dempster, Petitjean & Webb 2020) + ridge.
# bake-off redux에서 정확도/시간 trade-off 최상위.  C++ 커널(rocket_transform_cpp)
# 이 convolution을 처리한다.  MultiROCKET 근사: ppv + max + mean-of-positives.
# ---------------------------------------------------------------------
.rocket_kernels <- function(p, nkernel = 2000L, Tn = 100L, seed = 1L) {
  set.seed(seed)
  lens <- sample(c(7L, 9L, 11L), nkernel, TRUE)
  lapply(seq_len(nkernel), function(i) {
    L <- lens[i]
    W <- matrix(stats::rnorm(L * p), p, L)
    W <- W - rowMeans(W)
    dil <- 2^stats::runif(1, 0, log2(max((Tn - 1) / (L - 1), 1)))
    list(W = W, len = L, bias = stats::runif(1, -1, 1),
         dilation = max(1L, as.integer(dil)),
         padding = if (stats::runif(1) < .5) 0L else as.integer(((L - 1) * max(1L, as.integer(dil))) / 2))
  })
}
fit_rocket <- function(ds, ctx, nkernel = 2000L, ...) {
  R <- flc_regularize(ds, ctx$grid)
  ker <- .rocket_kernels(ds$p, nkernel, length(R$time), seed = ctx$seed)
  Xtr <- rocket_transform_cpp(R$traj, ker, flc_kernel_threads())
  Xtr[!is.finite(Xtr)] <- 0
  ctr <- colMeans(Xtr); sdv <- apply(Xtr, 2, stats::sd); sdv[sdv < 1e-12] <- 1
  Z <- scale(Xtr, ctr, sdv)
  mdl <- .learn_glmnet(Z, ds$y, alpha = 0)      # ridge, per the ROCKET paper
  structure(list(mdl = mdl, ker = ker, ctr = ctr, sdv = sdv, grid = R$time,
                 levels = levels(ds$y), ctx = ctx, method = "ROCKET-ridge"),
            class = "flcRocket")
}
predict.flcRocket <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  X <- rocket_transform_cpp(R$traj, object$ker, flc_kernel_threads())
  X[!is.finite(X)] <- 0
  Z <- scale(X, object$ctr, object$sdv); Z[!is.finite(Z)] <- 0
  P <- object$mdl$predict(Z)
  cls <- object$levels
  P <- as.matrix(P)[, cls, drop = FALSE]
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls), D = -P)
}

# ---------------------------------------------------------------------
# MFPCA (multivariate functional PCA) + LDA / SVM
# ---------------------------------------------------------------------
.fpc_scores_fit <- function(traj, ncomp = 6L) {
  X <- .flatten(traj)
  ctr <- colMeans(X)
  Xc <- sweep(X, 2, ctr)
  sv <- svd(Xc, nu = 0, nv = min(ncomp, ncol(Xc), nrow(Xc) - 1L))
  V <- sv$v
  list(ctr = ctr, V = V, scores = Xc %*% V,
       varexp = sv$d[seq_len(ncol(V))]^2 / sum(sv$d^2))
}
.fpc_scores_apply <- function(fit, traj) sweep(.flatten(traj), 2, fit$ctr) %*% fit$V

fit_mfpca <- function(ds, ctx, ncomp = 6L, learner = c("lda", "svm", "knn"), ...) {
  learner <- match.arg(learner)
  R <- flc_regularize(ds, ctx$grid)
  Fm <- .fpc_scores_fit(R$traj, ncomp)
  L <- switch(learner, lda = .learn_lda, svm = .learn_svm, knn = .learn_knn)
  mdl <- L(Fm$scores, ds$y)
  structure(list(fpc = Fm, mdl = mdl, grid = R$time, levels = levels(ds$y),
                 ctx = ctx, method = paste0("MFPCA-", toupper(learner))),
            class = "flcFPCA")
}
predict.flcFPCA <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  S <- .fpc_scores_apply(object$fpc, R$traj)
  P <- as.matrix(object$mdl$predict(S))
  cls <- object$levels
  P <- P[, cls, drop = FALSE]
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls), D = -P)
}
fit_mfpca_lda <- function(ds, ctx, ...) fit_mfpca(ds, ctx, learner = "lda", ...)
fit_mfpca_svm <- function(ds, ctx, ...) fit_mfpca(ds, ctx, learner = "svm", ...)

# ---------------------------------------------------------------------
# MFPLS-DA.  Preda, Saporta & Leveder (2007) / MFPLS (2024).
# 다중클래스는 클래스 지시행렬 Y (n x K)에 대한 PLS2로 처리하고,
# 잠재점수 위에서 LDA를 적합한다 (PLS-DA의 표준 형태).
# ---------------------------------------------------------------------
fit_mfpls <- function(ds, ctx, ncomp = 5L, ...) {
  if (!flc_have("pls")) stop("pls not installed")
  R <- flc_regularize(ds, ctx$grid)
  X <- .flatten(R$traj)
  y <- ds$y; cls <- levels(y)
  Y <- stats::model.matrix(~ y - 1)
  nc <- min(ncomp, ncol(X) - 1L, nrow(X) - 2L)
  df <- data.frame(row.names = seq_len(nrow(X)))
  df$Y <- Y; df$X <- X
  m <- pls::plsr(Y ~ X, ncomp = nc, data = df, method = "kernelpls")
  S <- m$scores[, seq_len(nc), drop = FALSE]
  lda <- .learn_lda(as.matrix(S), y)
  structure(list(m = m, lda = lda, nc = nc, grid = R$time, levels = cls,
                 ctx = ctx, method = "MFPLS-DA"), class = "flcPLS")
}
predict.flcPLS <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  X <- .flatten(R$traj)
  S <- stats::predict(object$m, newdata = data.frame(X = I(X)),
                      type = "scores")[, seq_len(object$nc), drop = FALSE]
  P <- as.matrix(object$lda$predict(as.matrix(S)))
  cls <- object$levels
  P <- P[, cls, drop = FALSE]
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls), D = -P)
}

# ---------------------------------------------------------------------
# Functional subspace classification.  Fukuda et al. (2023).
# 클래스 c 의 곡선들이 span하는 부분공간 U_c (SVD의 상위 q개 좌특이벡터)를
# 잡고, 새 곡선 x 를 각 U_c 에 정사영해서 그 "정준각/유사도"
#     sim_c(x) = || U_c^T x ||^2 / || x ||^2
# 가 가장 큰 클래스로 배정한다.  q는 inner-CV로 고른다.
# ---------------------------------------------------------------------
fit_subspace <- function(ds, ctx, q = NULL, qgrid = c(1, 2, 3, 5, 8), ...) {
  R <- flc_regularize(ds, ctx$grid)
  X <- .flatten(R$traj)
  ctr <- colMeans(X)
  Xc <- sweep(X, 2, ctr)
  nrm <- sqrt(rowSums(Xc^2)); nrm[nrm < 1e-12] <- 1
  Xn <- Xc / nrm
  y <- ds$y; cls <- levels(y)

  mkU <- function(qq) lapply(cls, function(c) {
    A <- Xn[y == c, , drop = FALSE]
    sv <- svd(t(A), nu = min(qq, nrow(A), ncol(A)), nv = 0)
    sv$u
  })
  simf <- function(U, Z) vapply(U, function(u)
    colSums((t(Z %*% u))^2), numeric(nrow(Z)))

  if (is.null(q)) {
    qgrid <- qgrid[qgrid < min(table(y))]
    if (!length(qgrid)) qgrid <- 1
    fold <- strat_folds(y, min(5L, min(table(y))), ctx$seed)
    accs <- vapply(qgrid, function(qq) {
      a <- vapply(sort(unique(fold)), function(f) {
        tr <- fold != f; te <- !tr
        Ut <- lapply(cls, function(c) {
          A <- Xn[tr & y == c, , drop = FALSE]
          if (nrow(A) < 1) return(matrix(0, ncol(Xn), 1))
          svd(t(A), nu = min(qq, nrow(A)), nv = 0)$u })
        S <- simf(Ut, Xn[te, , drop = FALSE])
        mean(cls[max.col(S, ties.method = "first")] == as.character(y[te]))
      }, 0)
      mean(a) }, 0)
    q <- qgrid[which.max(accs)]
  }
  structure(list(U = mkU(q), q = q, ctr = ctr, grid = R$time, levels = cls,
                 ctx = ctx, method = "Subspace (Fukuda 2023)"),
            class = "flcSubspace")
}
predict.flcSubspace <- function(object, ds_te, ctx = object$ctx, ...) {
  R <- flc_regularize(ds_te, object$grid)
  X <- sweep(.flatten(R$traj), 2, object$ctr)
  nrm <- sqrt(rowSums(X^2)); nrm[nrm < 1e-12] <- 1
  Xn <- X / nrm
  cls <- object$levels
  S <- vapply(object$U, function(u) rowSums((Xn %*% u)^2), numeric(nrow(Xn)))
  S <- matrix(S, nrow(Xn), length(cls), dimnames = list(NULL, cls))
  P <- pmax(S, 0); rs <- rowSums(P); rs[rs <= 0] <- 1
  list(prob = P / rs,
       class = factor(cls[max.col(S, ties.method = "first")], levels = cls), D = -S)
}

# ---------------------------------------------------------------------
# fda.usc wrappers (Febrero-Bande & Oviedo de la Fuente).
# 이 패키지는 1차원 functional data만 받으므로, 변수별로 적합한 뒤
# 클래스 사후확률을 곱(로그합)한다 -- naive multivariate 확장.
# ---------------------------------------------------------------------
.fdausc <- function(ds, ctx, method) {
  if (!flc_have("fda.usc")) stop("fda.usc not installed")
  R <- flc_regularize(ds, ctx$grid)
  y <- ds$y
  fits <- lapply(seq_len(ds$p), function(k) {
    fd <- fda.usc::fdata(R$traj[, , k], argvals = R$time)
    switch(method,
      knn    = fda.usc::classif.knn(y, fd, knn = c(1, 3, 5, 7, 9)),
      kernel = fda.usc::classif.kernel(y, fd),
      glm    = fda.usc::classif.glm(y ~ fd, data = list("df" = data.frame(y = y), "fd" = fd)),
      depth  = fda.usc::classif.depth(y, fd, depth = "FM"))
  })
  structure(list(fits = fits, grid = R$time, levels = levels(y), p = ds$p,
                 ctx = ctx, method = paste0("fda.usc-", method)),
            class = "flcFdausc")
}
predict.flcFdausc <- function(object, ds_te, ctx = object$ctx, ...) {
  # BUG FIX (2026-07):
  #   구버전은 예측이 실패하면 next 로 조용히 건너뛰었다.  fda.usc 의
  #   predict 는 type="probs" 를 안 받는 경우가 있어서 *모든* 변수가
  #   실패했고, 그 결과 LP 가 전부 0 -> 균등확률 -> 항상 클래스 1 예측
  #   -> bal_acc 가 정확히 1/K (=0.333) 로 분산 0.  버그의 명백한 신호.
  #
  #   이제 (1) 여러 반환 형식을 모두 시도하고, (2) 하나도 성공하지
  #   못하면 오류를 던져서 결과표에 NA 로 남게 한다.  조용히 틀린 답을
  #   내는 것보다 실패를 드러내는 편이 낫다.
  R   <- flc_regularize(ds_te, object$grid)
  cls <- object$levels
  K   <- length(cls)
  LP  <- matrix(0, ds_te$n, K, dimnames = list(NULL, cls))
  nok <- 0L

  for (k in seq_len(object$p)) {
    fd <- fda.usc::fdata(R$traj[, , k], argvals = R$time)
    P  <- .fdausc_prob(object$fits[[k]], fd, cls)
    if (is.null(P)) next
    LP  <- LP + log(pmax(P, 1e-9))
    nok <- nok + 1L
  }

  if (nok == 0L)
    stop("fda.usc: 모든 변수에서 예측 실패 (확률 추출 불가)")

  P <- exp(LP - apply(LP, 1, max))
  P <- P / rowSums(P)
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls),
       D = -P)
}

# fda.usc 의 predict 는 객체 종류마다 반환 형식이 다르다.
# 확률행렬을 뽑을 수 있는 모든 경로를 시도하고, 안 되면 NULL 을 돌려준다.
# (하드 레이블만 나오면 one-hot 으로 변환한다 -- AUC 는 나빠지지만
#  accuracy/F1 은 정상적으로 계산된다.)
.fdausc_prob <- function(fit, fd, cls) {
  K <- length(cls)

  # (1) type="probs"
  pr <- try(stats::predict(fit, fd, type = "probs"), silent = TRUE)
  if (!inherits(pr, "try-error")) {
    P <- NULL
    if (is.list(pr) && !is.null(pr$prob.group))      P <- as.matrix(pr$prob.group)
    else if (is.matrix(pr))                           P <- pr
    else if (is.data.frame(pr))                       P <- as.matrix(pr)
    if (!is.null(P) && ncol(P) == K) {
      colnames(P) <- cls
      return(P)
    }
  }

  # (2) 기본 predict -> 하드 레이블 -> one-hot
  pr <- try(stats::predict(fit, fd), silent = TRUE)
  if (!inherits(pr, "try-error")) {
    lab <- if (is.list(pr) && !is.null(pr$group.pred)) pr$group.pred else pr
    lab <- try(factor(as.character(lab), levels = cls), silent = TRUE)
    if (!inherits(lab, "try-error") && length(lab) > 0L && !all(is.na(lab))) {
      P <- matrix(1e-6, length(lab), K, dimnames = list(NULL, cls))
      ok <- !is.na(lab)
      P[cbind(which(ok), as.integer(lab[ok]))] <- 1
      return(P / rowSums(P))
    }
  }

  NULL
}
fit_fdausc_knn    <- function(ds, ctx, ...) .fdausc(ds, ctx, "knn")
fit_fdausc_kernel <- function(ds, ctx, ...) .fdausc(ds, ctx, "kernel")
fit_fdausc_glm    <- function(ds, ctx, ...) .fdausc(ds, ctx, "glm")
fit_fdausc_depth  <- function(ds, ctx, ...) .fdausc(ds, ctx, "depth")

# ---------------------------------------------------------------------
# Proximity Forest (Lucas et al. 2019) -- 간이 구현.
# 각 노드에서 무작위 (거리, exemplar) 를 뽑아 가장 가까운 exemplar로 분기.
# 우리 캐시(joint/dtw/euclid/sumvars)를 거리 후보 풀로 그대로 쓴다.
# 이렇게 하면 Frechet 거리가 forest 안에서 얼마나 선택되는지도 셀 수 있어
# "어느 거리가 유용한가"에 대한 부가 증거가 된다.
# ---------------------------------------------------------------------
fit_proxforest <- function(ds, ctx, ntree = 100L, ncand = 5L,
                           metrics = c("joint", "dtw", "euclid", "sumvars"), ...) {
  y <- ds$y; cls <- levels(y)
  w <- rep(1, ds$p)
  if ("sumvars" %in% metrics) {
    W <- try(flc_learn_weights(ctx$cache$Dvar[ctx$tr, ctx$tr, , drop = FALSE], y,
                               s_select = "cv", seed = ctx$seed), silent = TRUE)
    if (!inherits(W, "try-error")) w <- W$w
  }
  DL <- lapply(metrics, function(m) .train_dist(ctx, m, w))
  names(DL) <- metrics

  grow <- function(idx, depth) {
    yy <- y[idx]
    if (length(unique(yy)) == 1L || length(idx) < 3L || depth > 15L)
      return(list(leaf = TRUE,
                  prob = as.numeric(table(factor(yy, levels = cls))) / length(idx)))
    best <- NULL
    for (b in seq_len(ncand)) {
      m <- sample(metrics, 1)
      ex <- vapply(cls, function(c) {
        cand <- idx[yy == c]
        if (!length(cand)) NA_integer_ else sample(cand, 1) }, 1L)
      ok <- !is.na(ex)
      if (sum(ok) < 2L) next
      br <- max.col(-DL[[m]][idx, ex[ok], drop = FALSE], ties.method = "first")
      # gini gain
      g <- function(v) { t <- table(v); 1 - sum((t / sum(t))^2) }
      gain <- g(yy) - sum(vapply(unique(br), function(bb)
        mean(br == bb) * g(yy[br == bb]), 0))
      if (is.null(best) || gain > best$gain)
        best <- list(gain = gain, metric = m, ex = ex[ok], br = br,
                     exlab = cls[ok])
    }
    if (is.null(best) || best$gain <= 1e-9 || length(unique(best$br)) < 2L)
      return(list(leaf = TRUE,
                  prob = as.numeric(table(factor(yy, levels = cls))) / length(idx)))
    kids <- lapply(sort(unique(best$br)), function(bb)
      grow(idx[best$br == bb], depth + 1L))
    list(leaf = FALSE, metric = best$metric, ex = best$ex,
         branches = sort(unique(best$br)), kids = kids)
  }
  set.seed(ctx$seed)
  trees <- lapply(seq_len(ntree), function(b) {
    bs <- sample(seq_along(y), length(y), TRUE)
    grow(bs, 1L)
  })
  # 어떤 거리가 얼마나 선택됐는지 센다 (부가 결과)
  cnt <- new.env(); for (m in metrics) assign(m, 0L, cnt)
  walk <- function(nd) if (!nd$leaf) {
    assign(nd$metric, get(nd$metric, cnt) + 1L, cnt)
    lapply(nd$kids, walk) }
  invisible(lapply(trees, walk))
  usage <- vapply(metrics, function(m) get(m, cnt), 0L)

  structure(list(trees = trees, metrics = metrics, w = w, levels = cls,
                 usage = usage, ctx = ctx,
                 method = "ProximityForest (Lucas 2019)"), class = "flcPF")
}
predict.flcPF <- function(object, ds_te, ctx = object$ctx, ...) {
  cls <- object$levels
  DC <- lapply(object$metrics, function(m) .cross_dist(ctx, m, object$w))
  names(DC) <- object$metrics
  n <- length(ctx$te)
  P <- matrix(0, n, length(cls), dimnames = list(NULL, cls))
  for (tr in object$trees) {
    for (i in seq_len(n)) {
      nd <- tr
      while (!nd$leaf) {
        d <- DC[[nd$metric]][i, nd$ex]
        b <- which.min(d)
        j <- match(nd$branches[b], nd$branches)
        if (is.na(j) || j > length(nd$kids)) { nd <- list(leaf = TRUE,
          prob = rep(1 / length(cls), length(cls))); break }
        nd <- nd$kids[[j]]
      }
      P[i, ] <- P[i, ] + nd$prob
    }
  }
  P <- P / length(object$trees)
  list(prob = P,
       class = factor(cls[max.col(P, ties.method = "first")], levels = cls), D = -P)
}
