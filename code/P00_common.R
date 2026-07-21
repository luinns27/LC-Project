# =====================================================================
#  P00_common.R  --  모든 P0x 실행 스크립트가 맨 앞에서 source 한다.
#  ---------------------------------------------------------------
#  여기 있는 것:
#    - 환경/경로, flc 로드
#    - FLC_MAX_GRID = 250 패치      (01_data.R / 02_frechet.R 의 하드코딩 대체)
#    - 42개 방법 목록               (fda.usc 2종 제외)
#    - varsel_scores 안전 패치      (informative = NA 일 때 죽지 않게)
#    - 공통 지표/헬퍼
#
#  ! 원본 numbered 스크립트(00~15)는 건드리지 않는다.  전부 전역 재할당.
#  ! 콘솔 통째 복붙 금지.  RStudio Source 로 실행.
# =====================================================================

# ---------------------------------------------------------------------
# 0. 환경
# ---------------------------------------------------------------------
if (!exists("ENV")) ENV <- "main"

FLC_ROOTS <- list(
  main = "C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/LC",
  sub  = ""
)
DATA_ROOT <- "C:/Users/bassc/luinn27/DATA/long"

ROOT <- FLC_ROOTS[[ENV]]
if (is.null(ROOT) || !nzchar(ROOT) || !dir.exists(ROOT)) {
  stop("경로 확인: ", ROOT)
}
setwd(ROOT); options(FLC_ROOT = ROOT)

#  결과 저장 위치.  코드(ROOT)와 분리한다.
#  없으면 만든다.  RESULT_ROOT 를 미리 정의해 두면 그 값을 쓴다.
if (!exists("RESULT_ROOT")) {
  RESULT_ROOT <- file.path(dirname(ROOT), "Longitudinal Classification")
}
dir.create(RESULT_ROOT, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(RESULT_ROOT)) stop("결과 경로를 만들 수 없다: ", RESULT_ROOT)
options(FLC_RESULT_ROOT = RESULT_ROOT)

RDIR <- ROOT
if (dir.exists(file.path(ROOT, "R"))) RDIR <- file.path(ROOT, "R")
rfiles <- list.files(RDIR, pattern = "^[0-9]{2}b?_.*\\.R$", full.names = TRUE)
rfiles <- rfiles[order(basename(rfiles), method = "radix")]
for (f in rfiles) source(f, encoding = "UTF-8")

SRC_DIR <- ROOT
if (dir.exists(file.path(ROOT, "src"))) SRC_DIR <- file.path(ROOT, "src")
flc_setup(src_dir = SRC_DIR)

stopifnot(exists("flc_get_dataset"), exists("flc_prepare"), exists("fit_combo"),
          exists("flc_metrics"), exists("strat_folds"), exists("new_flcdata"),
          exists("flc_learn_weights"), exists("cache_dist_w"), exists("knn_vote"),
          exists("FLC_SUITE"), exists("FLC_SUITE_TIERS"))


# ---------------------------------------------------------------------
# 1. 공통 실행 파라미터
# ---------------------------------------------------------------------
if (!exists("SEEDS"))   SEEDS   <- 1:20
if (!exists("K_FOLDS")) K_FOLDS <- 5L

#  [중요] 격자 상한.  원본은 union grid > 60 이면 40점으로 압축했다.
#  비동기 설정에서 union grid 는 수천 점이므로 40점 압축은 grid=TRUE
#  비교군을 부당하게 불리하게 만든다 (= 제안법에 유리한 편향).
FLC_MAX_GRID <- 250L

#  kNN 의 k 후보 (홀수)
K_GRID <- c(1, 3, 5, 7, 9, 11)

#  kFDA 의 gamma 후보.
#  gamma 자체는 스케일 의존이므로 median heuristic 을 기준점으로 두고
#  로그 등간격 배율을 쓴다:  gamma = gamma_mult / (2 * median(D)^2)
#  gamma_mult = 1 이 표준 median heuristic.
GAMMA_GRID <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)

#  희소성 s 스윕 격자 크기
S_GRID_LEN <- 8L


# ---------------------------------------------------------------------
# 2. [패치] flc_regularize / flc_cross_euclid 의 격자 상한
# ---------------------------------------------------------------------
if (!exists(".flc_regularize_orig")) {
  .flc_regularize_orig <- flc_regularize
}
flc_regularize <- function(ds, grid = NULL, rule = 2) {
  same   <- length(unique(lapply(ds$tt, function(v) round(v, 10)))) == 1L
  anyNA_ <- any(vapply(ds$X, function(m) anyNA(m), TRUE))
  shrunk <- FALSE
  if (is.null(grid)) {
    grid <- if (same) {
      ds$tt[[1]]
    } else {
      sort(unique(round(unlist(ds$tt), 10)))
    }
    if (length(grid) > FLC_MAX_GRID) {
      grid <- seq(min(grid), max(grid), length.out = FLC_MAX_GRID)
      shrunk <- TRUE
    }
  }
  Tn <- length(grid); p <- ds$p; n <- ds$n
  traj <- array(NA_real_, c(n, Tn, p))
  for (i in seq_len(n)) {
    ti <- ds$tt[[i]]; Xi <- ds$X[[i]]
    for (k in seq_len(p)) {
      ok <- is.finite(Xi[, k])
      if (sum(ok) == 0L) {
        traj[i, , k] <- 0
      } else if (sum(ok) == 1L) {
        traj[i, , k] <- Xi[ok, k]
      } else {
        traj[i, , k] <- stats::approx(ti[ok], Xi[ok, k], xout = grid,
                                      rule = rule, ties = mean)$y
      }
    }
  }
  traj[!is.finite(traj)] <- 0
  list(traj = traj, time = grid, y = ds$y,
       interpolated = (!same) || anyNA_ || shrunk,
       grid_len = Tn, shrunk = shrunk)
}
assign("flc_regularize", flc_regularize, envir = globalenv())

if (exists("flc_cross_euclid") && !exists(".flc_cross_euclid_orig")) {
  .flc_cross_euclid_orig <- flc_cross_euclid
}
flc_cross_euclid <- function(dsA, dsB, w = rep(1, dsA$p), grid = NULL) {
  if (is.null(grid)) {
    grid <- sort(unique(round(unlist(c(dsA$tt, dsB$tt)), 10)))
    if (length(grid) > FLC_MAX_GRID) {
      grid <- seq(min(grid), max(grid), length.out = FLC_MAX_GRID)
    }
  }
  RA <- flc_regularize(dsA, grid); RB <- flc_regularize(dsB, grid)
  fl <- function(Z) t(apply(Z, 1, function(m)
    as.vector(sweep(matrix(m, dim(Z)[2], dim(Z)[3]), 2, sqrt(w), `*`))))
  A <- fl(RA$traj); B <- fl(RB$traj)
  sqrt(pmax(outer(rowSums(A^2), rowSums(B^2), `+`) - 2 * A %*% t(B), 0))
}
assign("flc_cross_euclid", flc_cross_euclid, envir = globalenv())


# ---------------------------------------------------------------------
# 3. [패치] varsel_scores  --  informative 가 NA 면 조용히 NULL
#     원본은 if (tp + fn) 에서 NA 를 만나 죽는데, 이 호출이 try() 밖이라
#     runner 전체가 무너진다.
# ---------------------------------------------------------------------
if (exists("varsel_scores")) {
  if (!exists(".varsel_scores_orig")) {
    .varsel_scores_orig <- varsel_scores
  }
  varsel_scores <- function(w, informative, tol = 1e-8) {
    if (is.null(informative)) return(NULL)
    inf <- suppressWarnings(as.logical(informative))
    if (length(w) != length(inf) || all(is.na(inf)) || anyNA(inf)) return(NULL)
    .varsel_scores_orig(w, inf, tol = tol)
  }
  assign("varsel_scores", varsel_scores, envir = globalenv())
}


# ---------------------------------------------------------------------
# 4. 42개 방법 목록
#     현재 44개에서 fda.usc-knn / fda.usc-kernel 을 제외한다.
#     제외 사유: classif.knn 은 Euclid + kNN 과 역할이 중복되고,
#                classif.kernel 은 fold 당 70초로 비용 대비 정보가 없다.
# ---------------------------------------------------------------------
DROP_METHODS <- c("fda.usc-knn", "fda.usc-kernel")

flc_methods_42 <- function() {
  M <- flc_method_set("all")
  M[!(names(M) %in% DROP_METHODS)]
}

#  보간이 필요한 방법인지 (분석 C: grid 그룹 대조용)
methods_grid_flag <- function(M = flc_methods_42()) {
  data.frame(method = names(M),
             family = vapply(M, function(m) m$family, ""),
             needs_grid = vapply(M, function(m) isTRUE(m$grid), TRUE),
             stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------
# 5. 설정 목록 (23종) + tier
# ---------------------------------------------------------------------
SETTINGS_ALL <- unique(unlist(FLC_SUITE_TIERS))
SETTINGS_ALL <- intersect(SETTINGS_ALL, names(FLC_SUITE))

tier_of <- function(nm) {
  for (tn in names(FLC_SUITE_TIERS)) {
    if (nm %in% FLC_SUITE_TIERS[[tn]]) return(sub("^[A-Z]_", "", tn))
  }
  NA_character_
}
SETTING_TIER <- data.frame(
  dataset = SETTINGS_ALL,
  tier    = vapply(SETTINGS_ALL, tier_of, ""),
  stringsAsFactors = FALSE)

#  변수선택 분석 대상: informative 가 정의된 설정만
SETTINGS_VARSEL <- c("Setting5", "Setting6", "Setting5_async", "Setting6_async",
                     "Shape_noise", "Shape_async", "Phase", "Phase_nowarp",
                     "HighDim20", "HighDim50", "Imbalanced", "Thyroid_like")
SETTINGS_VARSEL <- intersect(SETTINGS_VARSEL, SETTINGS_ALL)

#  시각화에서 뺄 고차원 설정 (분석에는 포함, 그림만 별도)
SETTINGS_HIGHDIM <- c("HighDim20", "HighDim50")

#  짝 대조 (분석 B)
PAIRED_CONTRASTS <- list(
  c("Phase",  "Phase_nowarp"),
  c("Shape",  "Shape_async"),
  c("Setting2", "Setting2_async"),
  c("Setting5", "Setting5_async"),
  c("Setting6", "Setting6_async")
)


# ---------------------------------------------------------------------
# 6. 공통 헬퍼
# ---------------------------------------------------------------------

#  안전한 논리 변환 (informative 가 NA 를 포함할 수 있다)
as_logical_safe <- function(v) {
  z <- suppressWarnings(as.logical(v))
  z[is.na(z)] <- FALSE
  z
}

#  prepare 캐시.  같은 (dataset, seed, lambda_method) 는 한 번만 계산한다.
.PREP_CACHE <- new.env(parent = emptyenv())

prep_for <- function(ds_name, seed = 1L, lambda_method = "rule",
                     sumOrMax = "max", ds = NULL, cache = TRUE) {
  ky <- paste(ds_name, seed, lambda_method, sumOrMax, sep = "|")
  if (cache && !is.null(.PREP_CACHE[[ky]])) return(.PREP_CACHE[[ky]])
  if (is.null(ds)) ds <- flc_get_dataset(ds_name, seed = seed)
  pr <- flc_prepare(ds, lambda_method = lambda_method,
                    sumOrMax = sumOrMax, verbose = FALSE)
  if (cache) .PREP_CACHE[[ky]] <- pr
  pr
}

prep_cache_clear <- function() {
  rm(list = ls(.PREP_CACHE), envir = .PREP_CACHE)
  invisible(gc(verbose = FALSE))
}

#  전체 지표를 한 행으로.  다중클래스면 AUC 를 NA 로 둔다.
metrics_row <- function(truth, pred, prob = NULL, lv = levels(factor(truth))) {
  m <- flc_metrics(truth, pred, prob = prob, lv = lv)
  r <- data.frame(
    acc        = m$acc,
    bal_acc    = m$bal_acc,
    macro_prec = m$macro_prec,
    macro_rec  = m$macro_rec,
    macro_f1   = m$macro_f1,
    micro_f1   = m$acc,
    kappa      = m$kappa,
    auc        = m$auc,
    logloss    = m$logloss,
    brier      = m$brier,
    stringsAsFactors = FALSE)
  #  다중클래스(K>2) 에서는 AUC 를 보고하지 않는다.
  if (length(lv) > 2L) r$auc <- NA_real_
  r
}

#  거리행렬 + 레이블로 kNN 예측 (k 를 주면 고정, 안 주면 train-LOO 로 선택)
knn_predict <- function(Dfull, y, tr, te, k = NULL, kg = K_GRID) {
  cls <- levels(y)
  kk  <- kg[kg < length(tr)]
  if (!length(kk)) kk <- 1L
  if (is.null(k)) {
    Dl <- Dfull[tr, tr, drop = FALSE]; diag(Dl) <- Inf
    ord <- t(apply(Dl, 1, order))
    a <- vapply(kk, function(k2)
      mean(apply(ord[, seq_len(k2), drop = FALSE], 1, function(ix)
        names(which.max(table(factor(y[tr][ix], levels = cls))))) ==
          as.character(y[tr])), 0)
    k <- kk[which.max(a)]
  }
  P  <- knn_vote(Dfull[te, tr, drop = FALSE], y[tr], k, cls)
  pr <- factor(cls[max.col(P, ties.method = "first")], levels = cls)
  list(pred = pr, prob = P, k = k)
}

#  kFDA 예측.  gamma_mult 를 주면 고정, 안 주면 train inner-CV 로 선택.
kfda_predict <- function(Dfull, y, tr, te, gamma_mult = NULL,
                         gg = GAMMA_GRID, seed = 1L, innerK = 3L,
                         reg = 1e-3) {
  cls <- levels(y)
  Dtr <- Dfull[tr, tr, drop = FALSE]
  med <- stats::median(Dtr[Dtr > 0])
  if (!is.finite(med) || med <= 0) med <- 1
  
  .fit_pred <- function(gm, i1, i2) {
    D1 <- Dfull[i1, i1, drop = FALSE]
    m1 <- stats::median(D1[D1 > 0]); if (!is.finite(m1) || m1 <= 0) m1 <- 1
    g  <- gm / (2 * m1^2)
    K1 <- exp(-g * D1^2)
    ft <- try(.kfda_fit(K1, y[i1], reg = reg), silent = TRUE)
    if (inherits(ft, "try-error")) return(NULL)
    K2 <- exp(-g * Dfull[i2, i1, drop = FALSE]^2)
    pr <- try(.kfda_predict(ft, K2, cls, Dfull[i2, i1, drop = FALSE]),
              silent = TRUE)
    if (inherits(pr, "try-error")) return(NULL)
    pr
  }
  
  if (is.null(gamma_mult)) {
    ytr <- y[tr]
    nf  <- min(innerK, min(table(ytr)))
    if (nf < 2L) {
      gamma_mult <- 1
    } else {
      fold <- strat_folds(ytr, nf, seed)
      sc <- vapply(gg, function(gm) {
        a <- vapply(sort(unique(fold)), function(f) {
          i1 <- tr[fold != f]; i2 <- tr[fold == f]
          pr <- .fit_pred(gm, i1, i2)
          if (is.null(pr)) return(NA_real_)
          flc_metrics(y[i2], pr$class)$bal_acc
        }, 0)
        mean(a, na.rm = TRUE)
      }, 0)
      gamma_mult <- if (all(!is.finite(sc))) 1 else gg[which.max(sc)]
    }
  }
  
  pr <- .fit_pred(gamma_mult, tr, te)
  if (is.null(pr)) {
    return(list(pred = factor(rep(cls[1], length(te)), levels = cls),
                prob = NULL, gamma_mult = gamma_mult, ok = FALSE))
  }
  list(pred = pr$class, prob = pr$prob, gamma_mult = gamma_mult, ok = TRUE)
}

#  쌍체 비교 요약 (같은 seed/fold 끼리)
paired_summary <- function(df, a, b, key = c("dataset", "seed", "fold"),
                           value = "bal_acc", method_col = "method") {
  A <- df[df[[method_col]] == a, ]
  B <- df[df[[method_col]] == b, ]
  if (!nrow(A) || !nrow(B)) return(NULL)
  ka <- do.call(paste, c(A[key], sep = "|"))
  kb <- do.call(paste, c(B[key], sep = "|"))
  common <- intersect(ka, kb)
  if (!length(common)) return(NULL)
  va <- A[[value]][match(common, ka)]
  vb <- B[[value]][match(common, kb)]
  d  <- va - vb
  d  <- d[is.finite(d)]
  if (length(d) < 2L) return(NULL)
  s <- stats::sd(d)
  tt <- if (s > 0) mean(d) / (s / sqrt(length(d))) else NA_real_
  wt <- try(stats::wilcox.test(va, vb, paired = TRUE, exact = FALSE),
            silent = TRUE)
  data.frame(method_a = a, method_b = b, n_pairs = length(d),
             mean_a = mean(va, na.rm = TRUE), mean_b = mean(vb, na.rm = TRUE),
             diff = mean(d), sd = s, t = tt,
             wins = sum(d > 0), ties = sum(d == 0),
             p_wilcox = if (inherits(wt, "try-error")) NA_real_ else wt$p.value,
             stringsAsFactors = FALSE)
}

#  Nemenyi critical difference
nemenyi_cd <- function(k, N, alpha = 0.05) {
  q <- c(1.960, 2.344, 2.569, 2.728, 2.850, 2.949, 3.031, 3.102, 3.164,
         3.219, 3.268, 3.313, 3.354, 3.391, 3.426, 3.458, 3.489, 3.517,
         3.544, 3.569)
  qa <- if (k <= 21) q[min(max(k - 1L, 1L), 20L)] else 3.6
  qa * sqrt(k * (k + 1) / (6 * N))
}

#  설정 내 평균순위 (동률 평균)
mean_ranks <- function(df, value = "bal_acc",
                       by = "dataset", method_col = "method") {
  ag <- stats::aggregate(df[[value]],
                         by = list(ds = df[[by]], m = df[[method_col]]),
                         FUN = function(z) mean(z, na.rm = TRUE))
  names(ag)[3] <- "v"
  ag$rank <- stats::ave(-ag$v, ag$ds,
                        FUN = function(z) rank(z, ties.method = "average"))
  R <- stats::aggregate(rank ~ m, ag, mean)
  names(R) <- c("method", "mean_rank")
  R[order(R$mean_rank), ]
}

#  색 (Okabe-Ito)
.PAL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
          "#56B4E9", "#F0E442", "#000000")
.C_SF <- "#0072B2"; .C_MF <- "#D55E00"; .C_NOISE <- "grey55"
.fade <- function(col, a) {
  r <- grDevices::col2rgb(col)
  grDevices::rgb(r[1], r[2], r[3], alpha = a * 255, maxColorValue = 255)
}

#  결과 폴더 만들기
mk_out <- function(nm) {
  d <- file.path(RESULT_ROOT, nm)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(d, "figs"), showWarnings = FALSE, recursive = TRUE)
  d
}

#  CSV 쓰기 (NULL 안전)
wr <- function(x, path) {
  if (is.null(x) || !NROW(x)) return(invisible(FALSE))
  utils::write.csv(x, path, row.names = FALSE)
  invisible(TRUE)
}

#  CSV 이어쓰기 (RESUME 용).  첫 호출이면 헤더를 쓴다.
wr_append <- function(x, path) {
  if (is.null(x) || !NROW(x)) return(invisible(FALSE))
  first <- !file.exists(path)
  utils::write.table(x, path, sep = ",", row.names = FALSE,
                     col.names = first, append = !first, qmethod = "double")
  invisible(TRUE)
}


# ---------------------------------------------------------------------
# 7. 벤치 실행 엔진
#     원본 flc_benchmark() 는 내부에서 flc_prepare 를 다시 부르므로
#     (dataset x seed) 마다 준비 비용이 중복된다.  여기서는 prep 을
#     밖에서 받아 재사용한다.
#
#     반환: list(results = <fold 단위 지표>, weights = <fold 단위 w>)
# ---------------------------------------------------------------------
flc_benchmark_one <- function(ds, ds_name, seed, prep, fold, METHODS,
                              verbose = TRUE) {
  y <- ds$y; lv <- levels(y); p <- ds$p
  rows <- list(); wrows <- list()
  folds <- sort(unique(fold))
  
  for (f in folds) {
    tr <- which(fold != f); te <- which(fold == f)
    
    #  fold 안에서 lambda 를 한 번만 고르고 모든 방법이 공유한다.
    #  (방법마다 다시 고르면 불공정하고 비싸다 -- 14_cv.R 주석과 동일 원칙)
    lam_idx <- if (identical(prep$lambda_method, "cv") &&
                   length(prep$lambda_grid) > 1L) {
      li <- try(.pick_lambda_train(prep, tr, seed), silent = TRUE)
      if (inherits(li, "try-error")) 1L else li
    } else {
      1L
    }
    lam_idx <- max(1L, min(lam_idx, length(prep$lambda_grid)))
    
    ctx <- flc_make_ctx(prep$ds, tr, te,
                        prep$caches[[lam_idx]],
                        prep$lambda_grid[lam_idx],
                        prep$gamma, prep$sumOrMax, prep$grid, seed)
    
    for (mi in seq_along(METHODS)) {
      M <- METHODS[[mi]]; mname <- names(METHODS)[mi]
      
      miss <- M$needs[!vapply(M$needs, flc_have, TRUE)]
      if (length(miss)) {
        rows[[length(rows) + 1L]] <- .bench_row(
          ds_name, mname, M, seed, f, lv, ok = FALSE,
          err = paste("missing:", paste(miss, collapse = ",")),
          secs = 0, lambda = ctx$lambda)
        next
      }
      
      t0 <- Sys.time()
      r <- try({
        fit <- do.call(M$fit, c(list(prep$ds[tr], ctx), M$args))
        pr  <- stats::predict(fit, prep$ds[te], ctx)
        list(fit = fit, pr = pr)
      }, silent = TRUE)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      
      if (inherits(r, "try-error")) {
        rows[[length(rows) + 1L]] <- .bench_row(
          ds_name, mname, M, seed, f, lv, ok = FALSE,
          err = substr(trimws(as.character(r)), 1, 200),
          secs = secs, lambda = ctx$lambda)
        next
      }
      
      mrow <- metrics_row(y[te], r$pr$class, prob = r$pr$prob, lv = lv)
      rows[[length(rows) + 1L]] <- .bench_row(
        ds_name, mname, M, seed, f, lv, ok = TRUE, err = "",
        secs = secs, lambda = ctx$lambda, mrow = mrow, fit = r$fit)
      
      w <- r$fit$w
      if (!is.null(w) && length(w) == p) {
        wrows[[length(wrows) + 1L]] <- data.frame(
          dataset = ds_name, method = mname, seed = seed, fold = f,
          variable = ds$varNames, w = as.numeric(w),
          informative = if (!is.null(ds$informative)) ds$informative else NA,
          stringsAsFactors = FALSE)
      }
    }
    if (verbose) cat(".")
  }
  
  list(results = do.call(rbind, rows),
       weights = if (length(wrows)) do.call(rbind, wrows) else NULL)
}

#  결과 한 행 조립 (성공/실패 공통)
.bench_row <- function(ds_name, mname, M, seed, f, lv, ok, err, secs,
                       lambda = NA_real_, mrow = NULL, fit = NULL) {
  if (is.null(mrow)) {
    mrow <- data.frame(acc = NA_real_, bal_acc = NA_real_,
                       macro_prec = NA_real_, macro_rec = NA_real_,
                       macro_f1 = NA_real_, micro_f1 = NA_real_,
                       kappa = NA_real_, auc = NA_real_,
                       logloss = NA_real_, brier = NA_real_)
  }
  gv <- function(nm) {
    if (is.null(fit) || is.null(fit[[nm]])) return(NA_real_)
    v <- fit[[nm]]
    if (length(v) != 1L) return(NA_real_)
    as.numeric(v)
  }
  # 선택된 변수 수.
  #  주의: flc_learn_weights 의 `sel` 은 논리벡터가 아니라 s 선택 진단
  #  리스트(list(s=, bcss=))다.  변수 개수는 w 에서 직접 센다.
  n_sel <- NA_real_
  if (!is.null(fit)) {
    if (!is.null(fit$w) && is.numeric(fit$w)) {
      n_sel <- sum(fit$w > 1e-6)
    } else if (!is.null(fit$sel) && is.logical(fit$sel)) {
      n_sel <- sum(fit$sel)
    }
  }
  cbind(
    data.frame(dataset = ds_name, method = mname, family = M$family,
               needs_grid = isTRUE(M$grid), seed = seed, fold = f,
               stringsAsFactors = FALSE),
    mrow,
    data.frame(
      k_sel     = gv("k"),
      gamma_sel = gv("gamma_mult"),
      s_sel     = gv("s"),
      n_sel     = n_sel,
      lambda    = lambda,
      seconds   = secs, ok = ok, error = err,
      stringsAsFactors = FALSE))
}

cat("P00_common.R 로드 완료\n")
cat(sprintf("  설정 %d종, 방법 %d개, SEEDS=%s, FLC_MAX_GRID=%d\n",
            length(SETTINGS_ALL), length(flc_methods_42()),
            paste(range(SEEDS), collapse = ":"), FLC_MAX_GRID))
cat("  코드:", ROOT, "\n")
cat("  결과:", RESULT_ROOT, "\n")