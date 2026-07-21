# =====================================================================
#  P06_plot_analysis.R   --  P02 (합성 심층분석) 결과 작도
#  ---------------------------------------------------------------
#  입력:  R02_analysis_synth/*.csv
#  출력:  R02_analysis_synth/figs/*.pdf
#
#  G01_A4_k              kNN 의 k-ablation
#  G02_A4_gamma          kFDA 의 gamma-ablation          [요청 3]
#  G03_A5_noise          노이즈 변수 robustness          [요청 4]
#  G04_A2_addcurve       변수 추가 곡선 (HighDim 제외)   [요청 5]
#  G04b_A2_highdim       HighDim20/50 전용               [요청 5]
#  G05_A3_sparsity       희소성 트레이드오프
#  G05b_A3_selected_k    s 별 선택된 k / gamma           [요청 7]
#  G06_A1_importance     변수 중요도 boxplot             [요청 6]
#  G07_A1_summary        연속 지표 요약 (TPR/FPR 없음)   [요청 6]
#  G08_A7_2x2            거리 x 희소성
#
#  실행: RStudio Source
# =====================================================================

ROOT <- "C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification"
RESULT_ROOT <- ROOT
setwd(ROOT)
source("P0X_plot_common.R", encoding = "UTF-8")
IN <- file.path(RESULT_ROOT, "R02_analysis_synth")
HIGHDIM <- c("HighDim20", "HighDim50")   # 요청 5: 시각화에서 분리

cat("\n=== P06: 심층분석 작도 ===\n입력:", IN, "\n")


# ---------------------------------------------------------------------
# G01 / G02.  A4: k-ablation (kNN) + gamma-ablation (kFDA)
#   저장 형식이 long 이다: rule / param_name / param / distance
# ---------------------------------------------------------------------
a4 <- rd(IN, "A4_kablation.csv")
if (!is.null(a4)) {

  # ---- G01: kNN, k ----------------------------------------------------
  d <- a4[a4$rule == "kNN", ]
  if (nrow(d)) {
    ag <- stats::aggregate(bal_acc ~ param + distance, d,
                           function(z) mean(z, na.rm = TRUE))
    ks <- sort(unique(ag$param))
    M <- matrix(unlist(lapply(c("SFCL", "MFCL"), function(nm) {
      s <- ag[ag$distance == nm, ]
      s$bal_acc[match(ks, s$param)]
    })), ncol = 2, dimnames = list(NULL, c("SFCL (sumvars)", "MFCL (joint)")))
    pdf_open(IN, "G01_A4_k.pdf", w = 6.6, h = 5.2)
    plot_lines(ks, M, cols = c(C_SF, C_MF),
               xlab = "k (number of neighbours)", ylab = "bal_acc",
               title = "kNN: sensitivity to k",
               sub = "flat curve = the distance geometry is sound, not just the rule",
               xaxis_at = ks)
    pdf_close()
    cat("  G01 완료\n")
  }

  # ---- G02: kFDA, gamma ----------------------------------------------
  #   gamma 는 로그 스케일이 자연스럽다.
  #   gamma = gamma_mult / (2 * median(D)^2),  gamma_mult=1 이 median heuristic.
  d <- a4[a4$rule == "kFDA", ]
  if (nrow(d)) {
    ag <- stats::aggregate(bal_acc ~ param + distance, d,
                           function(z) mean(z, na.rm = TRUE))
    gs <- sort(unique(ag$param))
    M <- matrix(unlist(lapply(c("SFCL", "MFCL"), function(nm) {
      s <- ag[ag$distance == nm, ]
      s$bal_acc[match(gs, s$param)]
    })), ncol = 2, dimnames = list(NULL, c("SFCL (sumvars)", "MFCL (joint)")))
    pdf_open(IN, "G02_A4_gamma.pdf", w = 6.6, h = 5.2)
    plot_lines(gs, M, cols = c(C_SF, C_MF), logx = TRUE,
               xlab = expression(paste("kernel width multiplier  ",
                                       gamma[mult], "   (log scale)")),
               ylab = "bal_acc",
               title = "kFDA: sensitivity to kernel width",
               sub = paste("gamma = gamma_mult / (2 x median(D)^2);",
                           "gamma_mult = 1 is the median heuristic"),
               xaxis_at = gs, xaxis_lab = formatC(gs, format = "g"))
    graphics::abline(v = 1, lty = 3, col = "grey55")
    pdf_close()
    cat("  G02 완료\n")
  }
}


# ---------------------------------------------------------------------
# G03.  A5: 노이즈 변수 추가에 대한 robustness  [요청 4]
# ---------------------------------------------------------------------
a5 <- rd(IN, "A5_noise.csv")
if (!is.null(a5)) {
  bases <- sort(unique(a5$base))
  pdf_open(IN, "G03_A5_noise.pdf", w = 3.9 * length(bases) + 0.6, h = 5.0)
  op <- graphics::par(mfrow = c(1, length(bases)),
                      mar = c(4.2, 4.0, 3.6, 0.8), oma = c(0, 0, 1.6, 0),
                      mgp = c(2.4, 0.7, 0))
  for (bs in bases) {
    d <- a5[a5$base == bs, ]
    ag <- stats::aggregate(bal_acc ~ n_noise + distance, d,
                           function(z) mean(z, na.rm = TRUE))
    ns <- sort(unique(ag$n_noise))
    graphics::plot(NA, xlim = range(ns), ylim = c(0.3, 1.02),
                   xlab = "# noise variables added", ylab = "bal_acc",
                   main = bs, cex.main = 1.0)
    for (j in seq_along(c("SFCL", "MFCL"))) {
      nm <- c("SFCL", "MFCL")[j]
      s <- ag[ag$distance == nm, ]
      graphics::lines(ns, s$bal_acc[match(ns, s$n_noise)], type = "b",
                      col = c(C_SF, C_MF)[j], pch = c(19, 17)[j], lwd = 2.2)
    }
    if (bs == bases[1]) {
      graphics::legend("bottomleft", c("SFCL (sparse)", "MFCL (dense)"),
                       col = c(C_SF, C_MF), pch = c(19, 17), lwd = 2.2,
                       bty = "n", cex = 0.72)
    }
  }
  graphics::mtext("Robustness to added noise variables", outer = TRUE,
                  line = 0.1, cex = 1.05, font = 2)
  graphics::par(op); pdf_close()
  cat("  G03 완료\n")
}


# ---------------------------------------------------------------------
# G04 / G04b.  A2: 변수 추가 곡선  [요청 5]
#   HighDim20/50 은 변수가 많아 곡선이 뭉개지므로 분리한다.
# ---------------------------------------------------------------------
a2 <- rd(IN, "A2_addcurve.csv")
if (!is.null(a2)) {

  draw_add <- function(d, title, sub) {
    dss <- sort(unique(d$dataset))
    cols <- pal_n(length(dss))
    ag <- stats::aggregate(bal_acc ~ n_var + dataset, d,
                           function(z) mean(z, na.rm = TRUE))
    xr <- range(ag$n_var); yr <- range(ag$bal_acc, na.rm = TRUE)
    op <- graphics::par(mar = c(4.2, 4.2, 3.8, 1.0), mgp = c(2.5, 0.7, 0))
    graphics::plot(NA, xlim = xr, ylim = c(max(0, yr[1] - 0.03), min(1, yr[2] + 0.03)),
                   xlab = "number of variables included (importance order)",
                   ylab = "bal_acc", main = "")
    for (j in seq_along(dss)) {
      s <- ag[ag$dataset == dss[j], ]
      s <- s[order(s$n_var), ]
      graphics::lines(s$n_var, s$bal_acc, type = "b", col = cols[j],
                      pch = 19, cex = 0.6, lwd = 1.9)
      #  각 곡선의 최고점 표시
      b <- which.max(s$bal_acc)
      graphics::points(s$n_var[b], s$bal_acc[b], pch = 1, cex = 1.5,
                       col = cols[j], lwd = 1.6)
    }
    graphics::title(main = title, line = 2.4, cex.main = 1.05)
    graphics::mtext(sub, side = 3, line = 0.9, cex = 0.75, col = "grey30")
    graphics::legend("bottomleft", dss, col = cols, lwd = 1.9, pch = 19,
                     bty = "n", cex = 0.68, ncol = if (length(dss) > 6) 2 else 1)
    graphics::par(op)
  }

  d1 <- a2[!(a2$dataset %in% HIGHDIM), ]
  if (nrow(d1)) {
    pdf_open(IN, "G04_A2_addcurve.pdf", w = 7.0, h = 5.4)
    draw_add(d1, "Variable-addition curve",
             "open circle = peak; accuracy declines once noise variables enter")
    pdf_close()
    cat("  G04 완료\n")
  }
  d2 <- a2[a2$dataset %in% HIGHDIM, ]
  if (nrow(d2)) {
    pdf_open(IN, "G04b_A2_highdim.pdf", w = 7.0, h = 5.4)
    draw_add(d2, "Variable-addition curve: high-dimensional settings",
             "p = 20 and 50; plotted separately because the x-range differs")
    pdf_close()
    cat("  G04b 완료\n")
  }
}


# ---------------------------------------------------------------------
# G05 / G05b.  A3: 희소성  [요청 7]
#   TPR/FPR 은 쓰지 않는다.  선택 변수 수 + 선택된 k / gamma 를 보인다.
# ---------------------------------------------------------------------
a3 <- rd(IN, "A3_sparsity.csv")
if (!is.null(a3)) {
  dss <- sort(unique(a3$dataset))

  # ---- G05: (선택 변수 수, 정확도) 궤적 -------------------------------
  pdf_open(IN, "G05_A3_sparsity.pdf", w = 7.2, h = 5.4)
  ag <- stats::aggregate(cbind(n_sel, bal_acc_knn, bal_acc_kfda) ~ s + dataset,
                         a3, function(z) mean(z, na.rm = TRUE))
  cols <- pal_n(length(dss))
  op <- graphics::par(mar = c(4.2, 4.2, 3.8, 1.0), mgp = c(2.5, 0.7, 0))
  graphics::plot(NA, xlim = rev(range(ag$n_sel)),
                 ylim = range(ag$bal_acc_knn, na.rm = TRUE),
                 xlab = "number of selected variables  (sparser ->)",
                 ylab = "bal_acc (kNN)", main = "")
  for (j in seq_along(dss)) {
    s <- ag[ag$dataset == dss[j], ]; s <- s[order(-s$n_sel), ]
    graphics::lines(s$n_sel, s$bal_acc_knn, type = "b", col = cols[j],
                    pch = 19, cex = 0.6, lwd = 1.9)
  }
  graphics::title(main = "Sparsity trade-off", line = 2.4, cex.main = 1.05)
  graphics::mtext("x-axis reversed: moving right means keeping fewer variables",
                  side = 3, line = 0.9, cex = 0.75, col = "grey30")
  graphics::legend("bottomleft", dss, col = cols, lwd = 1.9, pch = 19,
                   bty = "n", cex = 0.68, ncol = if (length(dss) > 6) 2 else 1)
  graphics::par(op); pdf_close()
  cat("  G05 완료\n")

  # ---- G05b: s 별로 고른 k / gamma ------------------------------------
  if (all(c("k_sel", "gamma_sel") %in% names(a3))) {
    pdf_open(IN, "G05b_A3_selected_hyper.pdf", w = 8.4, h = 4.6)
    op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.6, 0.8),
                        oma = c(0, 0, 1.4, 0), mgp = c(2.4, 0.7, 0))
    ak <- stats::aggregate(k_sel ~ n_sel, a3, function(z) mean(z, na.rm = TRUE))
    ak <- ak[order(ak$n_sel), ]
    graphics::plot(ak$n_sel, ak$k_sel, type = "b", pch = 19, col = C_SF,
                   lwd = 2.0, xlab = "number of selected variables",
                   ylab = "selected k (kNN)", main = "kNN")
    ag2 <- stats::aggregate(gamma_sel ~ n_sel, a3,
                            function(z) mean(z, na.rm = TRUE))
    ag2 <- ag2[order(ag2$n_sel), ]
    graphics::plot(ag2$n_sel, ag2$gamma_sel, type = "b", pch = 17, col = C_MF,
                   lwd = 2.0, log = "y",
                   xlab = "number of selected variables",
                   ylab = expression(paste("selected  ", gamma[mult])),
                   main = "kFDA")
    graphics::abline(h = 1, lty = 3, col = "grey55")
    graphics::mtext("Hyperparameters chosen as sparsity changes",
                    outer = TRUE, line = 0.0, cex = 1.02, font = 2)
    graphics::par(op); pdf_close()
    cat("  G05b 완료\n")
  }
}


# ---------------------------------------------------------------------
# G06.  A1: 변수 중요도 boxplot  [요청 6]
#   informative vs noise 를 색으로 구분.  TPR/FPR 은 쓰지 않고
#   부제에는 연속 지표(AUC-of-w, gap)만 적는다.
# ---------------------------------------------------------------------
a1 <- rd(IN, "A1_importance.csv")
s1 <- rd(IN, "A1_summary.csv")
if (!is.null(a1)) {
  a1$informative <- as.logical(a1$informative)
  a1$informative[is.na(a1$informative)] <- FALSE
  dss <- unique(a1$dataset)
  MAXV <- 24L
  pdf_open(IN, "G06_A1_importance.pdf", w = 8.4, h = 5.0)
  for (dn in dss) {
    d <- a1[a1$dataset == dn, ]
    mw <- stats::aggregate(w ~ variable + informative, d,
                           function(z) mean(z, na.rm = TRUE))
    mw <- mw[order(-mw$w), ]
    shown <- mw$variable
    trimmed <- FALSE
    if (nrow(mw) > MAXV) {
      ki <- mw$variable[mw$informative]
      kn <- utils::head(mw$variable[!mw$informative], MAXV - length(ki))
      shown <- mw$variable[mw$variable %in% c(ki, kn)]
      trimmed <- TRUE
    }
    d2 <- d[d$variable %in% shown, ]
    d2$variable <- factor(d2$variable, levels = shown)
    isinfo <- vapply(levels(d2$variable), function(v)
      any(d2$informative[d2$variable == v]), TRUE)
    cols <- ifelse(isinfo, fade(C_SF, 0.6), fade(C_NOI, 0.65))

    op <- graphics::par(mar = c(6.4, 4.4, 3.8, 1.0), mgp = c(2.7, 0.7, 0))
    graphics::boxplot(w ~ variable, data = d2, col = cols, border = "grey35",
                      las = 2, xlab = "", ylab = "learned weight  w",
                      outpch = 1, outcex = 0.4, whisklty = 1,
                      cex.axis = 0.68, main = "")
    graphics::abline(h = 0, col = "grey80", lty = 3)
    graphics::title(main = sprintf("%s -- variable importance", dn),
                    line = 2.4, cex.main = 1.05)
    sb <- ""
    if (!is.null(s1) && dn %in% s1$dataset) {
      r <- s1[s1$dataset == dn, ][1, ]
      sb <- sprintf("AUC(w) = %.3f    gap = %.3f    info share of total w = %.2f",
                    r$auc_w, r$gap, r$w_ratio)
    }
    if (trimmed) {
      sb <- paste0(sb, sprintf("    [showing %d of %d variables]",
                               length(shown), nrow(mw)))
    }
    graphics::mtext(sb, side = 3, line = 0.9, cex = 0.72, col = "grey30")
    graphics::legend("topright", c("informative", "noise"),
                     fill = c(fade(C_SF, 0.6), fade(C_NOI, 0.65)),
                     border = "grey35", bty = "n", cex = 0.72)
    graphics::par(op)
  }
  pdf_close()
  cat("  G06 완료 (", length(dss), "페이지)\n")
}


# ---------------------------------------------------------------------
# G07.  A1 요약: 연속 지표만  [요청 6]
# ---------------------------------------------------------------------
if (!is.null(s1)) {
  s1 <- s1[order(-s1$auc_w), ]
  pdf_open(IN, "G07_A1_summary.pdf", w = 7.4, h = 5.4)
  op <- graphics::par(mar = c(4.4, mar_left_for(s1$dataset, 0.75), 3.8, 1.2),
                      mgp = c(2.5, 0.7, 0))
  M <- rbind(s1$auc_w, s1$cor_spearman, s1$w_ratio)
  graphics::barplot(M, beside = TRUE, horiz = TRUE, names.arg = s1$dataset,
                    las = 1, cex.names = 0.72, xlim = c(0, 1.05),
                    col = c(fade(C_SF, 0.75), fade(C_DTW, 0.7),
                            fade(C_MF, 0.6)),
                    border = "grey40", xlab = "value", main = "")
  graphics::abline(v = 1, lty = 3, col = "grey65")
  graphics::title(main = "Recovery of ground-truth importance", line = 2.4,
                  cex.main = 1.05)
  graphics::mtext(
    "threshold-free measures only; AUC(w) = 1 means every informative variable outranks every noise variable",
    side = 3, line = 0.9, cex = 0.68, col = "grey30")
  graphics::legend("bottomright",
                   c("AUC of w", "Spearman cor", "info share of total w"),
                   fill = c(fade(C_SF, 0.75), fade(C_DTW, 0.7), fade(C_MF, 0.6)),
                   border = "grey40", bty = "n", cex = 0.72)
  graphics::par(op); pdf_close()
  cat("  G07 완료\n")
}


# ---------------------------------------------------------------------
# G08.  A7: 거리 x 희소성 2x2
# ---------------------------------------------------------------------
a7 <- rd(IN, "A7_ablation2x2.csv")
if (!is.null(a7)) {
  CN <- c("euclid_dense", "euclid_sparse", "frechet_dense", "frechet_sparse")
  CN <- CN[CN %in% names(a7)]
  dss <- sort(unique(a7$dataset))
  ag <- do.call(rbind, lapply(dss, function(dn) {
    d <- a7[a7$dataset == dn, ]
    vapply(CN, function(c) mean(d[[c]], na.rm = TRUE), 0)
  }))
  rownames(ag) <- dss
  cols <- c(fade(C_EU, 0.55), fade(C_EU, 0.9),
            fade(C_SF, 0.55), fade(C_SF, 0.95))[seq_along(CN)]

  pdf_open(IN, "G08_A7_2x2.pdf", w = 8.0, h = 5.4)
  op <- graphics::par(mar = c(6.8, 4.4, 3.8, 1.0), mgp = c(2.6, 0.7, 0))
  bp <- graphics::barplot(t(ag), beside = TRUE, col = cols, border = "grey40",
                          ylim = c(0, 1.06), las = 2, cex.names = 0.72,
                          ylab = "bal_acc", main = "")
  graphics::abline(h = 0.5, lty = 3, col = "grey70")
  graphics::title(main = "Distance x sparsity (2 x 2)", line = 2.4,
                  cex.main = 1.05)
  graphics::mtext("Frechet + sparse is best in almost every setting",
                  side = 3, line = 0.9, cex = 0.75, col = "grey30")
  graphics::legend("bottomright", gsub("_", " ", CN), fill = cols,
                   border = "grey40", bty = "n", cex = 0.7, ncol = 2)
  graphics::par(op); pdf_close()

  # ---- 쌍체 대비 요약 (콘솔 + 그림) -----------------------------------
  cmp <- list(c("frechet_sparse", "euclid_sparse"),
              c("frechet_dense",  "euclid_dense"),
              c("frechet_sparse", "frechet_dense"),
              c("euclid_sparse",  "euclid_dense"))
  rows <- lapply(cmp, function(cc) {
    if (!all(cc %in% names(a7))) return(NULL)
    d <- a7[[cc[1]]] - a7[[cc[2]]]
    d <- d[is.finite(d)]
    s <- stats::sd(d)
    data.frame(contrast = paste(cc, collapse = " - "), diff = mean(d),
               t = if (s > 0) mean(d) / (s / sqrt(length(d))) else NA_real_,
               n_pairs = length(d), wins = sum(d > 0),
               stringsAsFactors = FALSE)
  })
  R7 <- do.call(rbind, rows)
  if (!is.null(R7)) {
    utils::write.csv(R7, file.path(IN, "A7_contrasts.csv"), row.names = FALSE)
    pdf_open(IN, "G08b_A7_contrasts.pdf", w = 6.8, h = 4.2)
    plot_paired_dots(R7, "contrast", n_col = "n_pairs",
                     title = "A7 paired contrasts",
                     sub = "positive favours the first term; bar = 95% CI")
    pdf_close()
    print(R7, row.names = FALSE, digits = 3)
  }
  cat("  G08 완료\n")
}

cat("\n=== P06 완료 ===\n->", file.path(IN, "figs"), "\n")
