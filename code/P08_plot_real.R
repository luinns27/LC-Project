# =====================================================================
#  P08_plot_real.R   --  P04 (실데이터) 결과 작도
#  ---------------------------------------------------------------
#  입력:  R04_real_<name>/*.csv        (여러 폴더를 자동 탐색)
#  출력:  각 폴더의 figs/*.pdf  +  합본 R04_real_ALL/figs/*.pdf
#
#  자료마다:
#    K01_box            방법별 bal_acc (우연수준 표시)
#    K02_metric_heat    방법 x 지표
#    K03_ablation_k     kNN k-ablation
#    K04_ablation_gamma kFDA gamma-ablation
#    K05_paired         주요 쌍체 대비
#    K06_addcurve       변수 추가 곡선
#    K07_sparsity       희소성 트레이드오프
#    K08_weights        변수 가중치 프로파일
#    K09_confusion      혼동행렬 (SFCL/MFCL)
#
#  합본:
#    K10_all_datasets   자료 x 상위방법 히트맵
# =====================================================================

ROOT <- "C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification"
RESULT_ROOT <- ROOT
setwd(ROOT)
source("P0X_plot_common.R", encoding = "UTF-8")

DIRS <- list.dirs(RESULT_ROOT, recursive = FALSE)
DIRS <- DIRS[grepl("R04_real_", basename(DIRS))]
DIRS <- DIRS[basename(DIRS) != "R04_real_ALL"]

if (!length(DIRS)) stop("R04_real_* 폴더가 없다.  P04 를 먼저 실행할 것.")

cat("\n=== P08: 실데이터 작도 ===\n")
cat("  자료", length(DIRS), "개:",
    paste(sub("^R04_real_", "", basename(DIRS)), collapse = ", "), "\n")

TOP_N <- 18L
ALL <- list()


# =====================================================================
# 자료 하나
# =====================================================================
for (DIR in DIRS) {
  dn <- sub("^R04_real_", "", basename(DIR))
  cat("\n--", dn, "--\n")

  res <- rd(DIR, "results.csv")
  if (is.null(res)) next
  res <- res[res$ok %in% c(TRUE, "TRUE"), ]
  if (!nrow(res)) next

  K <- length(unique(c(res$dataset)))          # placeholder
  chance <- NA_real_
  pr <- rd(DIR, "A6_predictions.csv")
  if (!is.null(pr)) chance <- 1 / length(unique(pr$true))

  # ---- K01: 방법별 분포 ---------------------------------------------
  pdf_open(DIR, "K01_box.pdf", w = 8, h = 6.2)
  plot_method_box(res, "bal_acc", top_n = TOP_N,
                  title = sprintf("%s -- balanced accuracy", dn),
                  sub = sprintf("%d seeds x %d folds",
                                length(unique(res$seed)),
                                length(unique(res$fold))),
                  chance = chance)
  pdf_close()

  # ---- K02: 방법 x 지표 ---------------------------------------------
  METS <- c("bal_acc", "acc", "macro_f1", "macro_prec", "kappa")
  METS <- METS[METS %in% names(res)]
  ag <- stats::aggregate(res[METS], by = list(method = res$method),
                         FUN = function(z) mean(z, na.rm = TRUE))
  ag <- ag[order(-ag$bal_acc), ]
  if (nrow(ag) > TOP_N) ag <- ag[seq_len(TOP_N), ]
  M <- as.matrix(ag[, METS, drop = FALSE])
  rownames(M) <- short_method(ag$method)
  pdf_open(DIR, "K02_metric_heat.pdf", w = 7.4, h = 6.2)
  plot_heat(M, title = sprintf("%s -- method x metric", dn),
            sub = "AUC omitted for multi-class problems")
  pdf_close()

  ALL[[dn]] <- ag

  # ---- K03 / K04: k / gamma ablation ---------------------------------
  ak <- rd(DIR, "T_ablation_k.csv")
  if (!is.null(ak)) {
    a <- stats::aggregate(bal_acc ~ k + distance, ak,
                          function(z) mean(z, na.rm = TRUE))
    ks <- sort(unique(a$k))
    Mk <- matrix(unlist(lapply(c("SFCL", "MFCL"), function(nm) {
      s <- a[a$distance == nm, ]; s$bal_acc[match(ks, s$k)]
    })), ncol = 2, dimnames = list(NULL, c("SFCL", "MFCL")))
    pdf_open(DIR, "K03_ablation_k.pdf", w = 6.6, h = 5.0)
    plot_lines(ks, Mk, cols = c(C_SF, C_MF),
               xlab = "k (number of neighbours)", ylab = "bal_acc",
               title = sprintf("%s -- kNN k-ablation", dn),
               sub = "main table reports the value chosen inside each fold",
               xaxis_at = ks)
    pdf_close()
  }
  ag2 <- rd(DIR, "T_ablation_gamma.csv")
  if (!is.null(ag2)) {
    a <- stats::aggregate(bal_acc ~ gamma_mult + distance, ag2,
                          function(z) mean(z, na.rm = TRUE))
    gs <- sort(unique(a$gamma_mult))
    Mg <- matrix(unlist(lapply(c("SFCL", "MFCL"), function(nm) {
      s <- a[a$distance == nm, ]; s$bal_acc[match(gs, s$gamma_mult)]
    })), ncol = 2, dimnames = list(NULL, c("SFCL", "MFCL")))
    pdf_open(DIR, "K04_ablation_gamma.pdf", w = 6.6, h = 5.0)
    plot_lines(gs, Mg, cols = c(C_SF, C_MF), logx = TRUE,
               xlab = expression(paste("kernel width multiplier  ",
                                       gamma[mult], "   (log scale)")),
               ylab = "bal_acc",
               title = sprintf("%s -- kFDA gamma-ablation", dn),
               sub = "gamma_mult = 1 is the median heuristic",
               xaxis_at = gs, xaxis_lab = formatC(gs, format = "g"))
    graphics::abline(v = 1, lty = 3, col = "grey55")
    pdf_close()
  }

  # ---- K05: 쌍체 대비 -------------------------------------------------
  P5 <- rd(DIR, "T05_paired.csv")
  if (!is.null(P5)) {
    P5$lab <- sprintf("%s\n vs %s", short_method(P5$method_a),
                      short_method(P5$method_b))
    pdf_open(DIR, "K05_paired.pdf", w = 7.2, h = 5.2)
    plot_paired_dots(P5, "lab",
                     title = sprintf("%s -- paired contrasts", dn),
                     sub = "positive = first method better; bar = 95% CI",
                     xlab = "difference in bal_acc")
    pdf_close()
  }

  # ---- K06: 변수 추가 곡선 -------------------------------------------
  A2 <- rd(DIR, "A2_addcurve.csv")
  if (!is.null(A2)) {
    a <- stats::aggregate(bal_acc ~ n_var, A2,
                          function(z) mean(z, na.rm = TRUE))
    a <- a[order(a$n_var), ]
    sdv <- stats::aggregate(bal_acc ~ n_var, A2,
                            function(z) stats::sd(z, na.rm = TRUE))
    sdv <- sdv[order(sdv$n_var), ]
    pdf_open(DIR, "K06_addcurve.pdf", w = 6.8, h = 5.0)
    op <- graphics::par(mar = c(4.2, 4.2, 3.8, 1.0), mgp = c(2.5, 0.7, 0))
    yl <- range(c(a$bal_acc - sdv$bal_acc, a$bal_acc + sdv$bal_acc),
                na.rm = TRUE)
    graphics::plot(NA, xlim = range(a$n_var), ylim = yl,
                   xlab = "number of variables included (importance order)",
                   ylab = "bal_acc", main = "")
    graphics::polygon(c(a$n_var, rev(a$n_var)),
                      c(a$bal_acc - sdv$bal_acc, rev(a$bal_acc + sdv$bal_acc)),
                      col = fade(C_SF, 0.13), border = NA)
    graphics::lines(a$n_var, a$bal_acc, type = "b", col = C_SF, pch = 19,
                    lwd = 2.1, cex = 0.7)
    b <- which.max(a$bal_acc)
    graphics::points(a$n_var[b], a$bal_acc[b], pch = 1, cex = 1.8,
                     col = C_MF, lwd = 1.8)
    graphics::title(main = sprintf("%s -- variable-addition curve", dn),
                    line = 2.4, cex.main = 1.05)
    graphics::mtext(sprintf("peak at %d variables; band = 1 sd across folds",
                            a$n_var[b]), side = 3, line = 0.9, cex = 0.75,
                    col = "grey30")
    graphics::par(op); pdf_close()
  }

  # ---- K07: 희소성 ---------------------------------------------------
  A3 <- rd(DIR, "A3_sparsity.csv")
  if (!is.null(A3)) {
    a <- stats::aggregate(cbind(n_sel, bal_acc_knn, bal_acc_kfda) ~ s, A3,
                          function(z) mean(z, na.rm = TRUE))
    a <- a[order(-a$n_sel), ]
    pdf_open(DIR, "K07_sparsity.pdf", w = 6.8, h = 5.0)
    M <- cbind(kNN = a$bal_acc_knn, kFDA = a$bal_acc_kfda)
    plot_lines(a$n_sel, M, cols = c(C_SF, C_DTW),
               xlab = "number of selected variables",
               ylab = "bal_acc",
               title = sprintf("%s -- sparsity trade-off", dn),
               sub = "no sparsity level beats the dense baseline here")
    pdf_close()
  }

  # ---- K08: 가중치 프로파일 -------------------------------------------
  W <- rd(DIR, "A1_weight_profile.csv")
  if (!is.null(W)) {
    W <- W[order(-W$mean_w), ]
    pdf_open(DIR, "K08_weights.pdf", w = 7.2, h = 5.0)
    op <- graphics::par(mar = c(6.6, 4.4, 3.8, 1.0), mgp = c(2.6, 0.7, 0))
    bp <- graphics::barplot(W$mean_w, names.arg = W$variable, las = 2,
                            cex.names = 0.66, col = fade(C_SF, 0.7),
                            border = "grey40", ylab = "mean learned weight  w",
                            main = "")
    if ("sd_w" %in% names(W)) {
      graphics::arrows(bp, W$mean_w - W$sd_w, bp, W$mean_w + W$sd_w,
                       angle = 90, code = 3, length = 0.02,
                       col = "grey45", lwd = 0.9)
    }
    graphics::title(main = sprintf("%s -- variable weights (SFCL + kNN)", dn),
                    line = 2.4, cex.main = 1.05)
    nz <- sum(W$zero_rate > 0.5)
    graphics::mtext(sprintf("%d of %d variables set to zero in >50%% of folds",
                            nz, nrow(W)),
                    side = 3, line = 0.9, cex = 0.75, col = "grey30")
    graphics::par(op); pdf_close()
  }

  # ---- K09: 혼동행렬 ---------------------------------------------------
  if (!is.null(pr)) {
    mts <- sort(unique(pr$method))
    pdf_open(DIR, "K09_confusion.pdf", w = 5.4, h = 5.0)
    for (mt in mts) {
      d <- pr[pr$method == mt, ]
      lv <- sort(unique(c(d$true, d$pred)))
      cm <- table(factor(d$true, levels = lv), factor(d$pred, levels = lv))
      cmn <- sweep(cm, 1, pmax(rowSums(cm), 1), "/")
      M <- as.matrix(cmn)
      dimnames(M) <- list(lv, lv)
      plot_heat(M, title = sprintf("%s -- %s", dn, mt),
                sub = "row-normalised confusion matrix (recall per class)",
                xlab = "predicted", ylab = "true", digits = 2)
    }
    pdf_close()
  }

  cat("   그림 저장:", file.path(DIR, "figs"), "\n")
}


# =====================================================================
# 합본: 자료 x 방법
# =====================================================================
if (length(ALL) > 1L) {
  OUTA <- file.path(RESULT_ROOT, "R04_real_ALL")
  dir.create(file.path(OUTA, "figs"), showWarnings = FALSE, recursive = TRUE)

  #  모든 자료에 공통으로 있는 방법만
  ms <- Reduce(intersect, lapply(ALL, function(d) d$method))
  if (length(ms) > 2L) {
    M <- do.call(cbind, lapply(ALL, function(d)
      d$bal_acc[match(ms, d$method)]))
    colnames(M) <- names(ALL); rownames(M) <- short_method(ms)
    M <- M[order(-rowMeans(M, na.rm = TRUE)), , drop = FALSE]
    if (nrow(M) > TOP_N) M <- M[seq_len(TOP_N), , drop = FALSE]

    grDevices::pdf(file.path(OUTA, "figs", "K10_all_datasets.pdf"),
                   width = 3.0 + 1.5 * ncol(M), height = 6.4)
    plot_heat(M, title = "Real data: methods common to all datasets",
              sub = "bal_acc; rows sorted by mean across datasets")
    grDevices::dev.off()

    utils::write.csv(data.frame(method = rownames(M), M,
                                check.names = FALSE),
                     file.path(OUTA, "T_all_datasets.csv"), row.names = FALSE)
    cat("\n합본 저장:", file.path(OUTA, "figs"), "\n")
  }
}

cat("\n=== P08 완료 ===\n")
