# =====================================================================
#  P01_bench_synthetic.R
#  ---------------------------------------------------------------
#  [1] 42개 방법 x 23개 합성 설정 x SEEDS x 5-fold 전체 벤치마크.
#      전 지표 기록 (acc, bal_acc, macro_prec/rec/f1, kappa, logloss,
#      brier, auc(2클래스만), seconds).
#
#  산출:
#    results.csv              원자료 (fold 단위)
#    weights.csv              변수 가중치 (fold 단위)
#    T01_main_<metric>.csv    방법 x 설정 평균표 (지표별)
#    T02_rank.csv             평균순위 + Nemenyi CD
#    T03_tier.csv             tier 별 요약
#    T04_timing.csv           [분석 A] 계산시간
#    T05_paired_contrast.csv  [분석 B] 짝 대조 (Phase/Shape/async)
#    T06_grid_group.csv       [분석 C] 보간 필요 vs 불필요 그룹 대조
#    T07_tierD_full.csv       [분석 D] 결합설정 전체 순위
#    T08_euclid_vs_frechet.csv 유클리드 vs 프레셰 쌍체 (요청 2)
#
#  ! 가장 오래 걸린다.  RESUME=TRUE 로 이어받기 가능.
# =====================================================================

ENV <- "main"
source("C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification/P00_common.R", encoding = "UTF-8")

SEEDS   <- 1:20
K_FOLDS <- 5L
LAMBDA_METHOD <- "cv"
RESUME  <- TRUE

OUT <- mk_out("R01_bench_synth")
METHODS <- flc_methods_42()

cat("\n=== P01: 합성 벤치마크 ===\n")
cat(sprintf("  설정 %d x 방법 %d x seed %d x fold %d = %d 실행\n",
            length(SETTINGS_ALL), length(METHODS), length(SEEDS), K_FOLDS,
            length(SETTINGS_ALL) * length(METHODS) * length(SEEDS) * K_FOLDS))
cat("  출력:", OUT, "\n\n")


# ---------------------------------------------------------------------
# 1. 실행
# ---------------------------------------------------------------------
f_res <- file.path(OUT, "results.csv")
f_wt  <- file.path(OUT, "weights.csv")

done_key <- character(0)
if (RESUME && file.exists(f_res)) {
  prev <- utils::read.csv(f_res, stringsAsFactors = FALSE)
  done_key <- unique(paste(prev$dataset, prev$seed, sep = "|"))
  cat(sprintf("  RESUME: %d (dataset x seed) 조합 완료됨\n", length(done_key)))
}

t_all <- Sys.time()
for (dn in SETTINGS_ALL) {
  for (sd in SEEDS) {
    ky <- paste(dn, sd, sep = "|")
    if (ky %in% done_key) next

    cat(sprintf("[%s | seed %d] ", dn, sd))
    t0 <- Sys.time()

    ds <- flc_get_dataset(dn, seed = sd)
    prep <- flc_prepare(ds, lambda_method = LAMBDA_METHOD, verbose = FALSE)
    fold <- strat_folds(ds$y, K_FOLDS, sd)

    B <- flc_benchmark_one(ds, dn, sd, prep, fold, METHODS)

    wr_append(B$results, f_res)
    wr_append(B$weights, f_wt)

    rm(prep); invisible(gc(verbose = FALSE))
    cat(sprintf("done (%.1f min)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
cat(sprintf("\n벤치 완료 (%.1f h)\n",
            as.numeric(difftime(Sys.time(), t_all, units = "hours"))))


# ---------------------------------------------------------------------
# 2. 집계
# ---------------------------------------------------------------------
res <- utils::read.csv(f_res, stringsAsFactors = FALSE)
res <- res[res$ok %in% c(TRUE, "TRUE"), ]
res <- merge(res, SETTING_TIER, by = "dataset", all.x = TRUE)


res$needs_grid <- res$needs_grid %in% c(TRUE, "TRUE")

METRICS <- c("acc", "bal_acc", "macro_prec", "macro_rec", "macro_f1",
             "kappa", "logloss", "brier", "auc")

cat("\n[집계] 지표별 방법 x 설정 표 ...\n")
for (mt in METRICS) {
  if (all(!is.finite(res[[mt]]))) next
  tb <- stats::aggregate(res[[mt]],
                         by = list(method = res$method, dataset = res$dataset),
                         FUN = function(z) mean(z, na.rm = TRUE))
  names(tb)[3] <- mt
  wide <- stats::reshape(tb, idvar = "method", timevar = "dataset",
                         direction = "wide")
  names(wide) <- sub(paste0("^", mt, "\\."), "", names(wide))
  ov <- stats::aggregate(res[[mt]], by = list(method = res$method),
                         FUN = function(z) mean(z, na.rm = TRUE))
  names(ov)[2] <- "OVERALL"
  wide <- merge(wide, ov, by = "method")
  wide <- wide[order(-wide$OVERALL), ]
  wr(wide, file.path(OUT, sprintf("T01_main_%s.csv", mt)))
}

# ---- T02: 평균순위 + Nemenyi -----------------------------------------
R <- mean_ranks(res, "bal_acc")
k <- length(unique(res$method)); N <- length(unique(res$dataset))
CD <- nemenyi_cd(k, N)
R$cd <- CD
R$indistinguishable_from_top <- R$mean_rank < (min(R$mean_rank) + CD)
wr(R, file.path(OUT, "T02_rank.csv"))
cat(sprintf("  Nemenyi CD (k=%d, N=%d) = %.2f -> 상위와 무차별한 방법 %d개\n",
            k, N, CD, sum(R$indistinguishable_from_top)))

# ---- T03: tier 요약 --------------------------------------------------
tb <- stats::aggregate(bal_acc ~ method + tier, res,
                       function(z) mean(z, na.rm = TRUE))
T3 <- stats::reshape(tb, idvar = "method", timevar = "tier", direction = "wide")
names(T3) <- sub("^bal_acc\\.", "", names(T3))
ov <- stats::aggregate(bal_acc ~ method, res, function(z) mean(z, na.rm = TRUE))
names(ov)[2] <- "OVERALL"
T3 <- merge(T3, ov, by = "method")
T3 <- T3[order(-T3$OVERALL), ]
wr(T3, file.path(OUT, "T03_tier.csv"))

# ---- T04: [분석 A] 계산시간 -------------------------------------------
if ("seconds" %in% names(res)) {
  T4 <- stats::aggregate(seconds ~ method, res,
                         function(z) c(mean = mean(z, na.rm = TRUE),
                                       median = stats::median(z, na.rm = TRUE),
                                       max = max(z, na.rm = TRUE)))
  T4 <- data.frame(method = T4$method,
                   sec_mean = T4$seconds[, "mean"],
                   sec_median = T4$seconds[, "median"],
                   sec_max = T4$seconds[, "max"],
                   stringsAsFactors = FALSE)
  ba <- stats::aggregate(bal_acc ~ method, res, function(z) mean(z, na.rm = TRUE))
  T4 <- merge(T4, ba, by = "method")
  T4$sec_per_point_balacc <- T4$sec_mean / pmax(T4$bal_acc, 1e-8)
  T4 <- T4[order(-T4$sec_mean), ]
  wr(T4, file.path(OUT, "T04_timing.csv"))
}

# ---- T05: [분석 B] 짝 대조 -------------------------------------------
cat("[집계] 짝 대조 ...\n")
rows <- list()
for (pc in PAIRED_CONTRASTS) {
  a <- pc[1]; b <- pc[2]
  if (!all(c(a, b) %in% res$dataset)) next
  for (m in sort(unique(res$method))) {
    A <- res[res$method == m & res$dataset == a, ]
    B <- res[res$method == m & res$dataset == b, ]
    if (!nrow(A) || !nrow(B)) next
    ka <- paste(A$seed, A$fold, sep = "|")
    kb <- paste(B$seed, B$fold, sep = "|")
    cm <- intersect(ka, kb)
    if (!length(cm)) next
    va <- A$bal_acc[match(cm, ka)]; vb <- B$bal_acc[match(cm, kb)]
    d <- vb - va
    d <- d[is.finite(d)]
    if (length(d) < 2L) next
    s <- stats::sd(d)
    rows[[length(rows) + 1L]] <- data.frame(
      base = a, perturbed = b, method = m, n_pairs = length(d),
      mean_base = mean(va, na.rm = TRUE), mean_pert = mean(vb, na.rm = TRUE),
      loss = -mean(d),
      t = if (s > 0) mean(d) / (s / sqrt(length(d))) else NA_real_,
      stringsAsFactors = FALSE)
  }
}
T5 <- do.call(rbind, rows)
if (!is.null(T5)) {
  T5 <- T5[order(T5$base, T5$loss), ]
  wr(T5, file.path(OUT, "T05_paired_contrast.csv"))
}

# ---- T06: [분석 C] 보간 그룹 대조 ------------------------------------
cat("[집계] 보간 필요/불필요 그룹 ...\n")
res$grid_group <- ifelse(res$needs_grid, "needs_grid", "no_grid")
T6 <- stats::aggregate(bal_acc ~ grid_group + dataset + tier, res,
                       function(z) mean(z, na.rm = TRUE))
T6 <- stats::reshape(T6, idvar = c("dataset", "tier"),
                     timevar = "grid_group", direction = "wide")
names(T6) <- sub("^bal_acc\\.", "", names(T6))
if (all(c("needs_grid", "no_grid") %in% names(T6))) {
  T6$gap <- T6$no_grid - T6$needs_grid
  T6 <- T6[order(-T6$gap), ]
}
wr(T6, file.path(OUT, "T06_grid_group.csv"))

# ---- T07: [분석 D] 결합설정 전체 순위 --------------------------------
cat("[집계] Tier D 전체 순위 ...\n")
dD <- res[res$tier == "coupled", ]
if (nrow(dD)) {
  T7 <- stats::aggregate(bal_acc ~ method + dataset, dD,
                         function(z) mean(z, na.rm = TRUE))
  T7 <- stats::reshape(T7, idvar = "method", timevar = "dataset",
                       direction = "wide")
  names(T7) <- sub("^bal_acc\\.", "", names(T7))
  T7$MEAN <- rowMeans(T7[, -1, drop = FALSE], na.rm = TRUE)
  T7 <- T7[order(-T7$MEAN), ]
  wr(T7, file.path(OUT, "T07_tierD_full.csv"))
  cat("  Tier D 상위 8:\n")
  print(utils::head(T7, 8), row.names = FALSE, digits = 3)
}

# ---- T08: 유클리드 vs 프레셰 (요청 2) --------------------------------
cat("[집계] Euclid vs Frechet ...\n")
PAIRS <- list(
  c("SFKmL-C + kNN",         "Euclid + kNN"),
  c("SFKmL-C(dense) + kNN",  "Euclid + kNN"),
  c("MFKmL-C + kNN",         "Euclid + kNN"),
  c("SFKmL-C + kFDA",        "Euclid + kFDA"),
  c("SFKmL-C + medoid",      "Euclid + medoid"),
  c("SFKmL-C + SVM",         "Euclid + SVM"),
  c("SFKmL-C + kNN",         "DTW + kNN"),
  c("SFKmL-C + kNN",         "1NN-DTW"),
  c("DTW + kNN",             "Euclid + kNN")
)
rows <- list()
for (pp in PAIRS) {
  #  전체
  r <- paired_summary(res, pp[1], pp[2])
  if (!is.null(r)) { r$scope <- "ALL"; rows[[length(rows) + 1L]] <- r }
  #  tier 별
  for (tt in unique(res$tier)) {
    sub <- res[res$tier == tt, ]
    r <- paired_summary(sub, pp[1], pp[2])
    if (!is.null(r)) { r$scope <- tt; rows[[length(rows) + 1L]] <- r }
  }
}
T8 <- do.call(rbind, rows)
wr(T8, file.path(OUT, "T08_euclid_vs_frechet.csv"))

cat("\n=== P01 완료 ===\n")
cat("->", OUT, "\n")
