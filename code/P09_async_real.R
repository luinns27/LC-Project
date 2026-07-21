# =====================================================================
#  P09_async_real.R
#  ---------------------------------------------------------------
#  P12 / P19  --  진정한 비동기(asynchronous) 실데이터 분석.
#
#  이 자료가 기존 3종(Epilepsy/JV/arabic)과 결정적으로 다른 점:
#    기존 : 프레임 간격 균등, 길이만 다름 (speed-varying)  -> DTW 의 홈그라운드
#    여기 : 개체마다, 그리고 같은 개체의 변수마다 관측 시각이 다름
#           -> Frechet 의 lambda*|t-t'| 항이 기여할 수 있는 유일한 구조
#
#  따라서 방법 선정 원칙이 다르다:
#    (a) 공통 격자를 요구하는 방법(grid=TRUE)은 원리적으로 부적합.
#        없는 격자를 발명해야 하므로 "방법"이 아니라 "보간 어댑터"를 재게 된다.
#        -> 기본 제외.  단 COMPARE_GRID_METHODS=TRUE 로 켜면 부록용으로 포함.
#    (b) Frechet/DTW 평균(mean/DBA/centroid)도 내부적으로 flc_regularize 를
#        호출하므로 제외.
#    (c) 남는 것: medoid / kNN / SVM / kFDA  x  {SFCL, MFCL, DTW}
#                 + 1NN-DTW, 1NN-DDTW, DistSpace-Frechet
#
#  n 이 크므로(P12 ~4000, P19 ~40000) 층화 부분표본을 20회 반복한다.
#  각 반복이 독립된 자료 실현이므로, seed 는 부분표본 추출 + fold 분할을
#  동시에 지배한다.
#
#  산출:  R09_async_<name>/
#    results.csv               원자료 (subsample x fold x method)
#    weights.csv
#    T01_main_<metric>.csv     방법별 전 지표
#    T02_rank.csv
#    T03_paired.csv            핵심 쌍체 대비
#    T04_timing.csv
#    T_ablation_k.csv / T_ablation_gamma.csv
#    A1_weight_profile.csv     변수 가중치 (희소성이 작동하는가)
#    A2_addcurve.csv
#    A3_sparsity.csv
#    A8_async_diag.csv         비동기 정도 진단
#    methods_excluded.csv
# =====================================================================

ENV <- "main"
source("P00_common.R", encoding = "UTF-8")

# ---------------------------------------------------------------------
# 0. 설정
# ---------------------------------------------------------------------
DATA <- "P12"              # "P12" | "P19" | "both"

N_SUBSAMPLE   <- 20L       # 랜덤 부분표본 반복 수
MAX_PER_CLASS <- 250L      # 클래스당 최대 개체 수 -> P12 기준 n=500
K_FOLDS       <- 5L
LAMBDA_METHOD <- "cv"

#  격자를 요구하는 방법을 포함할 것인가.
#  FALSE(기본): 비동기 자료에 원리적으로 부적합하므로 제외.
#  TRUE       : 부록/반론 대응용.  보간 후 강제로 돌린다.
COMPARE_GRID_METHODS <- FALSE

#  P19 는 클래스 불균형이 극심하다(양성 ~4%).  양성을 전부 쓰고
#  음성을 같은 수로 뽑는다.
BALANCE_MODE <- "min_class"   # "min_class" | "cap"

TARGETS <- if (DATA == "both") c("P12", "P19") else DATA


# ---------------------------------------------------------------------
# 1. 로딩
#    long CSV -> flcdata.  *** 보간하지 않는다 ***
#    개체 i 의 시점은 그 개체에서 관측된 시각의 합집합이고,
#    그 시점에 측정되지 않은 변수는 NA 로 남긴다.
#    flc_frechet 계열은 변수별로 is.finite 를 보므로 비동기가 보존된다.
# ---------------------------------------------------------------------
load_async <- function(name, root = DATA_ROOT) {
  f <- file.path(root, paste0(name, "_processed"),
                 paste0(name, "_long.csv"))
  if (!file.exists(f)) {
    stop("CSV 없음 (preprocess_physionet.R 먼저): ", f)
  }
  cat("    읽는 중 ...\n")
  L <- utils::read.csv(f, stringsAsFactors = FALSE)
  cat(sprintf("    %d 행\n", nrow(L)))
  L$id <- as.character(L$id)
  L
}

#  층화 부분표본 -> flcdata
subsample_flcdata <- function(L, name, seed,
                              max_per_class = MAX_PER_CLASS,
                              mode = BALANCE_MODE) {
  key <- L[!duplicated(L$id), c("id", "label")]
  set.seed(seed)

  by_cls <- split(key$id, key$label)
  if (mode == "min_class") {
    m <- min(vapply(by_cls, length, 0L), max_per_class)
  } else {
    m <- max_per_class
  }
  ids <- unlist(lapply(by_cls, function(v) {
    if (length(v) > m) sample(v, m) else v
  }), use.names = FALSE)
  ids <- sample(ids)                       # 순서 섞기

  S <- L[L$id %in% ids, ]
  vnames <- sort(unique(S$variable))
  key2 <- S[!duplicated(S$id), c("id", "label")]
  key2 <- key2[match(ids, key2$id), ]
  stopifnot(identical(as.character(key2$id), as.character(ids)))

  #  개체별 (T_i x p) 행렬 조립.  NA 를 채우지 않는 것이 핵심.
  Ssp <- split(S[, c("time", "variable", "value")],
               factor(S$id, levels = ids))
  X <- vector("list", length(ids)); tt <- vector("list", length(ids))
  for (k in seq_along(ids)) {
    di <- Ssp[[k]]
    ts <- sort(unique(di$time))
    M  <- matrix(NA_real_, length(ts), length(vnames),
                 dimnames = list(NULL, vnames))
    M[cbind(match(di$time, ts), match(di$variable, vnames))] <- di$value
    X[[k]]  <- M
    tt[[k]] <- as.numeric(ts)
  }

  ds <- new_flcdata(X = X, tt = tt, y = factor(key2$label),
                    varNames = vnames, name = name, informative = NULL,
                    meta = list(id = ids, seed = seed))
  ds$informative <- NULL
  ds
}

#  비동기 정도 진단
async_diag <- function(ds) {
  Ts <- vapply(ds$X, nrow, 0L)
  obs <- vapply(ds$X, function(m) sum(is.finite(m)), 0L)
  fill <- obs / (Ts * ds$p)
  #  변수 쌍이 같은 시점에 동시 관측되는 비율
  co <- vapply(ds$X, function(m) {
    f <- is.finite(m)
    if (ncol(f) < 2L) return(NA_real_)
    mean(rowSums(f) >= 2L)
  }, 0)
  data.frame(
    n = ds$n, p = ds$p, K = nlevels(ds$y),
    T_min = min(Ts), T_med = stats::median(Ts), T_max = max(Ts),
    fill_rate = mean(fill),
    cooccur_rate = mean(co, na.rm = TRUE),
    stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------
# 2. 방법 선정
# ---------------------------------------------------------------------
methods_async <- function(ALLM, include_grid = COMPARE_GRID_METHODS) {
  ng <- vapply(ALLM, function(m) isTRUE(m$grid), TRUE)
  im <- grepl("\\+ *mean$|\\+ *centroid$|\\+ *DBA$", names(ALLM))

  reason <- rep(NA_character_, length(ALLM))
  reason[ng] <- "requires a common time grid; not defined under asynchronous sampling"
  reason[im] <- "class barycentre calls flc_regularize (interpolation)"

  drop <- (if (include_grid) im else (ng | im))
  list(use = ALLM[!drop], drop = ALLM[drop], reason = reason[drop])
}


# ---------------------------------------------------------------------
# 3. 자료 하나 처리
# ---------------------------------------------------------------------
run_async_dataset <- function(dn) {
  cat("\n", strrep("=", 64), "\n", sep = "")
  cat("=== ", dn, "  (asynchronous real data) ===\n", sep = "")
  OUT <- mk_out(paste0("R09_async_", dn))
  t_all <- Sys.time()

  L <- load_async(dn)
  key <- L[!duplicated(L$id), c("id", "label")]
  cat(sprintf("  전체 %d명\n", nrow(key)))
  print(table(key$label))

  ALLM <- flc_methods_42()
  MS <- methods_async(ALLM)
  METHODS <- MS$use
  cat(sprintf("\n  방법 %d개 사용, %d개 제외 (COMPARE_GRID_METHODS=%s)\n",
              length(METHODS), length(MS$drop), COMPARE_GRID_METHODS))
  cat("  사용:", paste(names(METHODS), collapse = ", "), "\n")
  if (length(MS$drop)) {
    wr(data.frame(dataset = dn, method = names(MS$drop),
                  family = vapply(MS$drop, function(m) m$family, ""),
                  reason = MS$reason, stringsAsFactors = FALSE),
       file.path(OUT, "methods_excluded.csv"))
  }

  f_res <- file.path(OUT, "results.csv"); if (file.exists(f_res)) file.remove(f_res)
  f_wt  <- file.path(OUT, "weights.csv"); if (file.exists(f_wt)) file.remove(f_wt)
  f_dg  <- file.path(OUT, "A8_async_diag.csv")

  diag_rows <- list(); ab_k <- list(); ab_g <- list()
  a2_rows <- list(); a3_rows <- list()

  for (sd in seq_len(N_SUBSAMPLE)) {
    cat(sprintf("\n  [부분표본 %d/%d] ", sd, N_SUBSAMPLE))
    ds <- subsample_flcdata(L, dn, seed = sd)
    dg <- async_diag(ds); dg$seed <- sd
    diag_rows[[sd]] <- dg
    cat(sprintf("n=%d p=%d  fill=%.1f%% cooccur=%.1f%%\n",
                ds$n, ds$p, 100 * dg$fill_rate, 100 * dg$cooccur_rate))

    prep <- flc_prepare(ds, lambda_method = LAMBDA_METHOD, verbose = FALSE)
    fold <- strat_folds(ds$y, K_FOLDS, sd)

    #  --- 메인 벤치 ---
    B <- flc_benchmark_one(ds, dn, sd, prep, fold, METHODS)
    wr_append(B$results, f_res)
    wr_append(B$weights, f_wt)

    #  --- 부수 분석 (같은 캐시 재사용) ---
    ca <- prep$caches[[1]]
    sub <- .async_side(ds, ca, fold, sd, dn)
    ab_k[[sd]]   <- sub$k
    ab_g[[sd]]   <- sub$g
    a2_rows[[sd]] <- sub$a2
    a3_rows[[sd]] <- sub$a3

    rm(prep, ca, ds); invisible(gc(verbose = FALSE))
    cat("    done\n")
  }

  wr(do.call(rbind, diag_rows), f_dg)
  wr(do.call(rbind, ab_k), file.path(OUT, "T_ablation_k.csv"))
  wr(do.call(rbind, ab_g), file.path(OUT, "T_ablation_gamma.csv"))
  wr(do.call(rbind, a2_rows), file.path(OUT, "A2_addcurve.csv"))
  wr(do.call(rbind, a3_rows), file.path(OUT, "A3_sparsity.csv"))

  # ---- 집계 -------------------------------------------------------
  res <- utils::read.csv(f_res, stringsAsFactors = FALSE)
  res <- res[res$ok %in% c(TRUE, "TRUE"), ]
  cat(sprintf("\n  성공 %d 행\n", nrow(res)))

  METRICS <- c("acc", "bal_acc", "macro_prec", "macro_rec", "macro_f1",
               "kappa", "auc", "logloss", "brier")
  for (mt in METRICS) {
    if (!(mt %in% names(res)) || all(!is.finite(res[[mt]]))) next
    tb <- stats::aggregate(res[[mt]], by = list(method = res$method),
                           FUN = function(z) c(m = mean(z, na.rm = TRUE),
                                               s = stats::sd(z, na.rm = TRUE)))
    T <- data.frame(method = tb$method, mean = tb$x[, "m"],
                    sd = tb$x[, "s"], stringsAsFactors = FALSE)
    wr(T[order(-T$mean), ], file.path(OUT, sprintf("T01_main_%s.csv", mt)))
  }

  wr(mean_ranks(res, "bal_acc", by = "seed"), file.path(OUT, "T02_rank.csv"))

  #  핵심 쌍체: 비동기에서 Frechet 이 DTW 를 이기는가
  PAIRS <- list(
    c("SFKmL-C + kNN",       "DTW + kNN"),
    c("SFKmL-C + kNN",       "1NN-DTW"),
    c("SFKmL-C + kFDA",      "DTW + kFDA"),
    c("SFKmL-C + medoid",    "DTW + medoid"),
    c("SFKmL-C + kNN",       "MFKmL-C + kNN"),
    c("SFKmL-C + kNN",       "SFKmL-C(dense) + kNN"),
    c("SFKmL-C + kFDA",      "SFKmL-C(dense) + kFDA"),
    c("SFKmL-C + kNN",       "DistSpace-Frechet"))
  rows <- lapply(PAIRS, function(pp) paired_summary(res, pp[1], pp[2]))
  P <- do.call(rbind, rows)
  wr(P, file.path(OUT, "T03_paired.csv"))
  if (!is.null(P)) {
    cat("\n  === 핵심 쌍체 대비 ===\n")
    print(P[, c("method_a", "method_b", "diff", "t", "wins", "n_pairs")],
          row.names = FALSE, digits = 3)
  }

  if ("seconds" %in% names(res)) {
    tb <- stats::aggregate(seconds ~ method, res,
                           function(z) c(m = mean(z, na.rm = TRUE),
                                         s = stats::median(z, na.rm = TRUE)))
    wr(data.frame(method = tb$method, sec_mean = tb$seconds[, "m"],
                  sec_median = tb$seconds[, "s"], stringsAsFactors = FALSE),
       file.path(OUT, "T04_timing.csv"))
  }

  #  변수 가중치 프로파일: 희소성이 실제로 작동하는가
  if (file.exists(f_wt)) {
    W <- utils::read.csv(f_wt, stringsAsFactors = FALSE)
    W <- W[W$method == "SFKmL-C + kNN", ]
    if (nrow(W)) {
      tb <- stats::aggregate(w ~ variable, W,
                             function(z) c(m = mean(z), s = stats::sd(z),
                                           z0 = mean(z <= 1e-6)))
      A1 <- data.frame(variable = tb$variable, mean_w = tb$w[, "m"],
                       sd_w = tb$w[, "s"], zero_rate = tb$w[, "z0"],
                       stringsAsFactors = FALSE)
      A1 <- A1[order(-A1$mean_w), ]
      wr(A1, file.path(OUT, "A1_weight_profile.csv"))
      cat("\n  === 변수 가중치 (SFCL + kNN) ===\n")
      print(A1, row.names = FALSE, digits = 3)
      cat(sprintf("  -> %d/%d 변수가 절반 이상의 fold 에서 0\n",
                  sum(A1$zero_rate > 0.5), nrow(A1)))
    }
  }

  D <- do.call(rbind, diag_rows)
  cat(sprintf("\n  === 비동기 진단 (평균) ===\n"))
  cat(sprintf("    격자 채움율   %.1f%%   (100%% = 완전 동기)\n",
              100 * mean(D$fill_rate)))
  cat(sprintf("    변수 동시관측 %.1f%%   (낮을수록 변수간 비동기)\n",
              100 * mean(D$cooccur_rate)))

  cat(sprintf("\n  === %s 완료 (%.1f h) -> %s\n", dn,
              as.numeric(difftime(Sys.time(), t_all, units = "hours")), OUT))
  invisible(TRUE)
}


# ---------------------------------------------------------------------
# 4. 부수 분석 (k/gamma ablation, A2, A3)
# ---------------------------------------------------------------------
.async_side <- function(ds, ca, fold, sd, dn) {
  y <- ds$y; p <- ds$p
  krow <- list(); grow <- list(); a2 <- list(); a3 <- list()
  sg <- exp(seq(log(1.02), log(sqrt(p)), length.out = S_GRID_LEN))

  for (f in sort(unique(fold))) {
    tr <- which(fold != f); te <- which(fold == f)
    W  <- flc_learn_weights(ca$Dvar[tr, tr, , drop = FALSE], y[tr],
                            s_select = "cv", seed = sd)
    Dsf <- cache_dist_w(ca, W$w)
    Dmf <- ca$Djoint
    ordv <- order(W$w, decreasing = TRUE)

    #  k / gamma ablation
    for (nm in c("SFCL", "MFCL")) {
      D <- if (nm == "SFCL") Dsf else Dmf
      for (k in K_GRID) {
        if (k >= length(tr)) next
        kp <- knn_predict(D, y, tr, te, k = k)
        m  <- metrics_row(y[te], kp$pred, kp$prob, levels(y))
        krow[[length(krow) + 1L]] <- cbind(
          data.frame(dataset = dn, seed = sd, fold = f, distance = nm,
                     k = k, stringsAsFactors = FALSE), m)
      }
      for (gm in GAMMA_GRID) {
        kf <- kfda_predict(D, y, tr, te, gamma_mult = gm, seed = sd)
        if (!isTRUE(kf$ok)) next
        m <- metrics_row(y[te], kf$pred, kf$prob, levels(y))
        grow[[length(grow) + 1L]] <- cbind(
          data.frame(dataset = dn, seed = sd, fold = f, distance = nm,
                     gamma_mult = gm, stringsAsFactors = FALSE), m)
      }
    }

    #  A2: 중요도 순 변수 추가
    for (j in seq_len(p)) {
      wj <- numeric(p); wj[ordv[seq_len(j)]] <- 1 / j
      kp <- knn_predict(cache_dist_w(ca, wj), y, tr, te)
      m  <- metrics_row(y[te], kp$pred, kp$prob, levels(y))
      a2[[length(a2) + 1L]] <- cbind(
        data.frame(dataset = dn, seed = sd, fold = f, n_var = j,
                   k_sel = kp$k, stringsAsFactors = FALSE), m)
    }

    #  A3: 희소성 스윕
    for (s in sg) {
      Ws <- flc_learn_weights(ca$Dvar[tr, tr, , drop = FALSE], y[tr],
                              s_select = "fixed", s = s, seed = sd)
      D <- cache_dist_w(ca, Ws$w)
      kp <- knn_predict(D, y, tr, te)
      kf <- kfda_predict(D, y, tr, te, seed = sd)
      mk <- metrics_row(y[te], kp$pred, kp$prob, levels(y))
      mf <- metrics_row(y[te], kf$pred, kf$prob, levels(y))
      a3[[length(a3) + 1L]] <- data.frame(
        dataset = dn, seed = sd, fold = f, s = s,
        n_sel = sum(Ws$w > 1e-6), p = p,
        k_sel = kp$k, gamma_sel = kf$gamma_mult,
        bal_acc_knn = mk$bal_acc, acc_knn = mk$acc, auc_knn = mk$auc,
        bal_acc_kfda = mf$bal_acc, acc_kfda = mf$acc, auc_kfda = mf$auc,
        stringsAsFactors = FALSE)
    }
  }
  list(k = do.call(rbind, krow), g = do.call(rbind, grow),
       a2 = do.call(rbind, a2), a3 = do.call(rbind, a3))
}


# ---------------------------------------------------------------------
# 5. 실행
# ---------------------------------------------------------------------
cat("\n=== P09: 비동기 실데이터 ===\n")
cat("  대상:", paste(TARGETS, collapse = ", "), "\n")
cat(sprintf("  부분표본 %d회 x 클래스당 최대 %d명 x %d-fold\n",
            N_SUBSAMPLE, MAX_PER_CLASS, K_FOLDS))

for (dn in TARGETS) {
  r <- try(run_async_dataset(dn), silent = FALSE)
  if (inherits(r, "try-error")) {
    cat("  (!) 실패:", dn, "\n")
  }
  prep_cache_clear()
}

cat("\n=== P09 완료 ===\n")
