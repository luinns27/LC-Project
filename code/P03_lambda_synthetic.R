# =====================================================================
#  P03_lambda_synthetic.R
#  ---------------------------------------------------------------
#  [8] lambda 선택 방식 비교 (요청: 기존 그대로 유지).
#
#      rule  : 원 논문 rule-of-thumb (단일 lambda, unsupervised)
#      cv    : inner-CV 분류손실 (레이블 사용)
#      align : 거리행렬 target alignment (레이블 사용)
#
#  Frechet 공식은 건드리지 않는다.  모든 lambda>0 에서 metric 이므로
#  (Kang et al. 2023, Thm.1) 어느 선택이든 거리의 성질은 보존된다.
#
#  대상: lambda 가 성능을 좌우하는 위상/비동기 설정.
#  산출: results_lambda.csv (원자료), T_lambda.csv (요약)
# =====================================================================

ENV <- "main"
source("C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification/P00_common.R", encoding = "UTF-8")

SEEDS   <- 1:20
K_FOLDS <- 5L
OUT <- mk_out("R03_lambda")

LAMBDA_METHODS <- c("rule", "cv", "align")
DATASETS_LAM <- c("Phase", "Phase_nowarp", "Setting2_async",
                  "Setting5_async", "Setting6_async", "Shape_async")
DATASETS_LAM <- intersect(DATASETS_LAM, SETTINGS_ALL)

#  비교할 방법: 제안 방법의 대표 판본만 (전 방법을 돌릴 필요 없다)
METHOD_SUBSET <- c("SFKmL-C + kNN", "SFKmL-C + kFDA", "SFKmL-C + medoid",
                   "MFKmL-C + kNN", "MFKmL-C + kFDA", "MFKmL-C + medoid")

cat("\n=== P03: lambda 선택 비교 ===\n")
cat("  설정", length(DATASETS_LAM), "x 방법", length(METHOD_SUBSET),
    "x lambda", length(LAMBDA_METHODS), "\n")
cat("  출력:", OUT, "\n\n")

ALLM <- flc_methods_42()
METHODS <- ALLM[intersect(METHOD_SUBSET, names(ALLM))]

f_res <- file.path(OUT, "results_lambda.csv")
if (file.exists(f_res)) file.remove(f_res)

t0 <- Sys.time()
for (lm in LAMBDA_METHODS) {
  for (dn in DATASETS_LAM) {
    for (sd in SEEDS) {
      cat(sprintf("[%s | %s | seed %d] ", lm, dn, sd))
      ds   <- flc_get_dataset(dn, seed = sd)
      prep <- flc_prepare(ds, lambda_method = lm, verbose = FALSE)
      fold <- strat_folds(ds$y, K_FOLDS, sd)
      B <- flc_benchmark_one(ds, dn, sd, prep, fold, METHODS, verbose = FALSE)
      if (!is.null(B$results)) {
        B$results$lambda_method <- lm
        wr_append(B$results, f_res)
      }
      rm(prep); invisible(gc(verbose = FALSE))
      cat("done\n")
    }
  }
}

# ---- 요약 -----------------------------------------------------------
res <- utils::read.csv(f_res, stringsAsFactors = FALSE)
res <- res[res$ok %in% c(TRUE, "TRUE"), ]

tb <- stats::aggregate(bal_acc ~ dataset + method + lambda_method, res,
                       function(z) mean(z, na.rm = TRUE))
T <- stats::reshape(tb, idvar = c("dataset", "method"),
                    timevar = "lambda_method", direction = "wide")
names(T) <- sub("^bal_acc\\.", "", names(T))
wr(T, file.path(OUT, "T_lambda.csv"))

#  쌍체: cv - rule,  align - rule
rows <- list()
for (dn in unique(res$dataset)) {
  for (m in unique(res$method)) {
    sub <- res[res$dataset == dn & res$method == m, ]
    if (!nrow(sub)) next
    for (lm in c("cv", "align")) {
      A <- sub[sub$lambda_method == lm, ]
      B <- sub[sub$lambda_method == "rule", ]
      if (!nrow(A) || !nrow(B)) next
      ka <- paste(A$seed, A$fold, sep = "|")
      kb <- paste(B$seed, B$fold, sep = "|")
      cm <- intersect(ka, kb); if (!length(cm)) next
      d <- A$bal_acc[match(cm, ka)] - B$bal_acc[match(cm, kb)]
      d <- d[is.finite(d)]; if (length(d) < 2L) next
      s <- stats::sd(d)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = dn, method = m, contrast = paste(lm, "- rule"),
        diff = mean(d), n = length(d),
        t = if (s > 0) mean(d) / (s / sqrt(length(d))) else NA_real_,
        wins = sum(d > 0), stringsAsFactors = FALSE)
    }
  }
}
wr(do.call(rbind, rows), file.path(OUT, "T_lambda_paired.csv"))

cat(sprintf("\n=== P03 완료 (%.1f h) ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
print(T, row.names = FALSE, digits = 3)
cat("->", OUT, "\n")
