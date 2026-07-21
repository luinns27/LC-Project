# =====================================================================
#  P05_plot_bench.R   --  P01 (합성 벤치마크) 결과 작도
#  ---------------------------------------------------------------
#  입력:  R01_bench_synth/*.csv
#  출력:  R01_bench_synth/figs/*.pdf
#
#  F01_overall_box       방법별 bal_acc 분포 (상위 N)
#  F02_metric_heat       방법 x 지표 히트맵
#  F03_tier_heat         방법 x tier 히트맵  (대칭적 실패)
#  F04_rank_cd           평균순위 + Nemenyi CD
#  F05_euclid_frechet    유클리드 vs 프레셰 쌍체
#  F06_paired_contrast   짝 대조 (thinning 손실)  [분석 B]
#  F07_grid_group        보간 필요/불필요 그룹    [분석 C]
#  F08_tierD             결합설정 전체 순위       [분석 D]
#  F09_timing            정확도 vs 계산시간       [분석 A]
#
#  실행: RStudio Source (CSV 만 읽으므로 수 초)
# =====================================================================

ROOT <- "C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification"
RESULT_ROOT <- ROOT
setwd(ROOT)
source("P0X_plot_common.R", encoding = "UTF-8")

IN  <- file.path(RESULT_ROOT, "R01_bench_synth")
TOP_N <- 20L

cat("\n=== P05: 벤치마크 작도 ===\n입력:", IN, "\n")

res <- rd(IN, "results.csv")
if (!is.null(res)) res <- res[res$ok %in% c(TRUE, "TRUE"), ]

# ---------------------------------------------------------------------
# F01. 방법별 bal_acc 분포
# ---------------------------------------------------------------------
if (!is.null(res)) {
  nds <- length(unique(res$dataset))
  pdf_open(IN, "F01_overall_box.pdf", w = 8, h = 6.4)
  plot_method_box(res, "bal_acc", top_n = TOP_N,
                  title = sprintf("Balanced accuracy across %d settings", nds),
                  sub = sprintf("%d seeds x %d folds; box = all runs",
                                length(unique(res$seed)),
                                length(unique(res$fold))))
  pdf_close()
  cat("  F01 완료\n")
}

# ---------------------------------------------------------------------
# F02. 방법 x 지표 히트맵
#   지표마다 스케일이 다르므로 열별로 0-1 정규화한 색을 쓰되
#   칸에는 원값을 적는다.
# ---------------------------------------------------------------------
if (!is.null(res)) {
  METS <- c("bal_acc", "acc", "macro_f1", "macro_prec", "kappa")
  METS <- METS[METS %in% names(res)]
  ag <- stats::aggregate(res[METS], by = list(method = res$method),
                         FUN = function(z) mean(z, na.rm = TRUE))
  ag <- ag[order(-ag$bal_acc), ]
  if (nrow(ag) > TOP_N) ag <- ag[seq_len(TOP_N), ]
  M <- as.matrix(ag[, METS, drop = FALSE])
  rownames(M) <- short_method(ag$method)
  #  색은 열별 정규화, 표시는 원값
  Mn <- apply(M, 2, function(z) {
    r <- range(z, na.rm = TRUE)
    if (diff(r) == 0) return(rep(0.5, length(z)))
    (z - r[1]) / diff(r)
  })
  pdf_open(IN, "F02_metric_heat.pdf", w = 7.6, h = 6.4)
  op <- graphics::par(mar = c(4.6, mar_left_for(rownames(M), 0.78), 3.6, 1.4))
  pal <- grDevices::colorRampPalette(
    c("#F7F7F7", "#FDDBC7", "#92C5DE", "#1B7837"))(64)
  nr <- nrow(M); nc <- ncol(M)
  graphics::image(seq_len(nc), seq_len(nr), t(Mn[nr:1, , drop = FALSE]),
                  col = pal, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(1, seq_len(nc), colnames(M), tick = FALSE, las = 2,
                 cex.axis = 0.78)
  graphics::axis(2, seq_len(nr), rev(rownames(M)), tick = FALSE, las = 1,
                 cex.axis = 0.75)
  graphics::box(col = "grey70")
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      graphics::text(j, i, formatC(M[nr - i + 1, j], format = "f", digits = 3),
                     cex = 0.6,
                     col = if (Mn[nr - i + 1, j] > 0.55) "white" else "grey15")
    }
  }
  graphics::title(main = "Method x metric", line = 2.0, cex.main = 1.05)
  graphics::mtext("colour scaled within each column; cells show raw values",
                  side = 3, line = 0.5, cex = 0.75, col = "grey30")
  graphics::par(op); pdf_close()
  cat("  F02 완료\n")
}

# ---------------------------------------------------------------------
# F03. 방법 x tier 히트맵  (SFCL/MFCL 의 대칭적 실패)
# ---------------------------------------------------------------------
T3 <- rd(IN, "T03_tier.csv")
if (!is.null(T3)) {
  tc <- setdiff(names(T3), c("method", "OVERALL"))
  ord <- c("paper", "frechet", "auxiliary", "coupled")
  tc <- c(intersect(ord, tc), setdiff(tc, ord))
  T3 <- T3[order(-T3$OVERALL), ]
  if (nrow(T3) > TOP_N) T3 <- T3[seq_len(TOP_N), ]
  M <- as.matrix(T3[, c(tc, "OVERALL"), drop = FALSE])
  rownames(M) <- short_method(T3$method)
  pdf_open(IN, "F03_tier_heat.pdf", w = 7.2, h = 6.4)
  plot_heat(M, title = "Method x setting tier (bal_acc)",
            sub = "coupled = cross-variable temporal coupling (MFCL's regime)")
  pdf_close()
  cat("  F03 완료\n")
}

# ---------------------------------------------------------------------
# F04. 평균순위 + Nemenyi CD
# ---------------------------------------------------------------------
R <- rd(IN, "T02_rank.csv")
if (!is.null(R)) {
  R <- R[order(R$mean_rank), ]
  CD <- if ("cd" %in% names(R)) R$cd[1] else NA_real_
  sh <- if (nrow(R) > TOP_N) R[seq_len(TOP_N), ] else R
  lab <- short_method(sh$method)
  pdf_open(IN, "F04_rank_cd.pdf", w = 7.6, h = 6.0)
  op <- graphics::par(mar = c(4.4, mar_left_for(lab, 0.72), 3.8, 1.2),
                      mgp = c(2.5, 0.7, 0))
  xr <- range(sh$mean_rank) + c(-0.6, 0.6)
  graphics::plot(NA, xlim = xr, ylim = c(0.5, nrow(sh) + 0.5), yaxt = "n",
                 xlab = "mean rank (lower is better)", ylab = "", main = "")
  graphics::axis(2, seq_len(nrow(sh)), rev(lab), las = 1, tick = FALSE,
                 cex.axis = 0.72)
  yy <- rev(seq_len(nrow(sh)))
  isp <- grepl("SFCL|MFCL", lab)
  graphics::points(sh$mean_rank, yy, pch = 19, cex = 1.05,
                   col = ifelse(isp, C_SF, "grey45"))
  if (is.finite(CD)) {
    b <- min(sh$mean_rank)
    graphics::abline(v = b + CD, lty = 2, col = C_MF)
    graphics::mtext(sprintf("CD = %.2f", CD), side = 3, at = b + CD,
                    line = 0.2, cex = 0.7, col = C_MF)
    graphics::rect(b, 0.4, b + CD, nrow(sh) + 0.6,
                   col = fade(C_MF, 0.06), border = NA)
  }
  graphics::title(main = "Mean rank with Nemenyi critical difference",
                  line = 2.4, cex.main = 1.05)
  graphics::mtext(
    "methods left of the dashed line are not distinguishable from the best",
    side = 3, line = 0.9, cex = 0.75, col = "grey30")
  graphics::par(op); pdf_close()
  cat("  F04 완료\n")
}

# ---------------------------------------------------------------------
# F05. 유클리드 vs 프레셰 (요청 2)
# ---------------------------------------------------------------------
T8 <- rd(IN, "T08_euclid_vs_frechet.csv")
if (!is.null(T8)) {
  a <- T8[T8$scope == "ALL", ]
  if (nrow(a)) {
    a$lab <- sprintf("%s\n vs %s", short_method(a$method_a),
                     short_method(a$method_b))
    pdf_open(IN, "F05_euclid_frechet.pdf", w = 7.4, h = 5.6)
    plot_paired_dots(a, "lab", title = "Distance contrasts (paired, all settings)",
                     sub = "positive = first method better; bar = 95% CI",
                     xlab = "difference in bal_acc")
    pdf_close()
  }
  #  tier 별 패널
  sc <- setdiff(unique(T8$scope), "ALL")
  if (length(sc)) {
    pdf_open(IN, "F05b_euclid_frechet_by_tier.pdf",
             w = 4.2 * length(sc), h = 5.4)
    op <- graphics::par(mfrow = c(1, length(sc)))
    for (s in sc) {
      b <- T8[T8$scope == s, ]
      if (!nrow(b)) next
      b$lab <- sprintf("%s vs %s", short_method(b$method_a),
                       short_method(b$method_b))
      plot_paired_dots(b, "lab", title = s, xlab = "diff in bal_acc")
    }
    graphics::par(op); pdf_close()
  }
  cat("  F05 완료\n")
}

# ---------------------------------------------------------------------
# F06. 짝 대조 [분석 B]
#   base -> perturbed 로 갈 때 각 방법이 얼마나 잃는가.
# ---------------------------------------------------------------------
T5 <- rd(IN, "T05_paired_contrast.csv")
if (!is.null(T5)) {
  prs <- unique(paste(T5$base, T5$perturbed, sep = " -> "))
  pdf_open(IN, "F06_paired_contrast.pdf", w = 8, h = 5.6)
  for (pr in prs) {
    pp <- strsplit(pr, " -> ", fixed = TRUE)[[1]]
    d <- T5[T5$base == pp[1] & T5$perturbed == pp[2], ]
    if (!nrow(d)) next
    d <- d[order(-d$loss), ]
    if (nrow(d) > TOP_N) {
      d <- rbind(utils::head(d, TOP_N %/% 2), utils::tail(d, TOP_N %/% 2))
    }
    lab <- short_method(d$method)
    op <- graphics::par(mar = c(4.4, mar_left_for(lab, 0.72), 3.8, 1.2),
                        mgp = c(2.5, 0.7, 0))
    isp <- grepl("SFCL|MFCL", lab)
    graphics::barplot(rev(d$loss), horiz = TRUE, names.arg = rev(lab),
                      las = 1, cex.names = 0.7,
                      col = ifelse(rev(isp), fade(C_SF, 0.7),
                                   fade("grey60", 0.6)),
                      border = "grey40",
                      xlab = "loss in bal_acc (base - perturbed)", main = "")
    graphics::abline(v = 0, col = "grey50")
    graphics::title(main = sprintf("%s  ->  %s", pp[1], pp[2]),
                    line = 2.4, cex.main = 1.05)
    graphics::mtext("larger bar = method breaks down under the perturbation",
                    side = 3, line = 0.9, cex = 0.75, col = "grey30")
    graphics::par(op)
  }
  pdf_close()
  cat("  F06 완료 (", length(prs), "페이지)\n")
}

# ---------------------------------------------------------------------
# F07. 보간 필요/불필요 그룹 [분석 C]
# ---------------------------------------------------------------------
T6 <- rd(IN, "T06_grid_group.csv")
if (!is.null(T6) && all(c("needs_grid", "no_grid") %in% names(T6))) {
  T6 <- T6[order(T6$gap), ]
  pdf_open(IN, "F07_grid_group.pdf", w = 7.6, h = 5.6)
  op <- graphics::par(mar = c(4.4, mar_left_for(T6$dataset, 0.75), 3.8, 1.2),
                      mgp = c(2.5, 0.7, 0))
  M <- rbind(T6$needs_grid, T6$no_grid)
  bp <- graphics::barplot(M, beside = TRUE, horiz = TRUE,
                          names.arg = T6$dataset, las = 1, cex.names = 0.72,
                          col = c(fade(C_EU, 0.75), fade(C_SF, 0.72)),
                          border = "grey40", xlim = c(0, 1),
                          xlab = "mean bal_acc", main = "")
  graphics::title(main = "Interpolation-based vs distance-native methods",
                  line = 2.4, cex.main = 1.05)
  graphics::mtext("grouped by the `needs_grid` flag; mean over all methods in the group",
                  side = 3, line = 0.9, cex = 0.75, col = "grey30")
  graphics::legend("bottomright",
                   c("needs common grid (interpolated)", "no grid needed"),
                   fill = c(fade(C_EU, 0.75), fade(C_SF, 0.72)),
                   border = "grey40", bty = "n", cex = 0.75)
  graphics::par(op); pdf_close()
  cat("  F07 완료\n")
}

# ---------------------------------------------------------------------
# F08. Tier D 전체 순위 [분석 D]
# ---------------------------------------------------------------------
T7 <- rd(IN, "T07_tierD_full.csv")
if (!is.null(T7)) {
  T7 <- T7[order(-T7$MEAN), ]
  sh <- if (nrow(T7) > TOP_N) T7[seq_len(TOP_N), ] else T7
  dsc <- setdiff(names(sh), c("method", "MEAN"))
  M <- as.matrix(sh[, c(dsc, "MEAN"), drop = FALSE])
  rownames(M) <- short_method(sh$method)
  pdf_open(IN, "F08_tierD.pdf", w = 6.8, h = 6.2)
  plot_heat(M, title = "Coupled settings: full ranking",
            sub = "every method shown, not a hand-picked subset")
  pdf_close()
  cat("  F08 완료\n")
}

# ---------------------------------------------------------------------
# F09. 정확도 vs 계산시간 [분석 A]
# ---------------------------------------------------------------------
T4 <- rd(IN, "T04_timing.csv")
if (!is.null(T4) && all(c("sec_mean", "bal_acc") %in% names(T4))) {
  T4 <- T4[is.finite(T4$sec_mean) & is.finite(T4$bal_acc), ]
  T4$sec_plot <- pmax(T4$sec_mean, 1e-3)
  pdf_open(IN, "F09_timing.pdf", w = 7.4, h = 5.6)
  op <- graphics::par(mar = c(4.4, 4.4, 3.8, 1.2), mgp = c(2.6, 0.7, 0))
  isp <- grepl("SFKmL|MFKmL", T4$method)
  graphics::plot(T4$sec_plot, T4$bal_acc, log = "x", pch = 19, cex = 1.0,
                 col = ifelse(isp, C_SF, "grey55"),
                 xlab = "seconds per fit (log scale)", ylab = "mean bal_acc",
                 main = "")
  #  라벨은 상위/느린 것만 (겹침 방지)
  showi <- which(T4$bal_acc > stats::quantile(T4$bal_acc, 0.75) |
                 T4$sec_plot > stats::quantile(T4$sec_plot, 0.85))
  graphics::text(T4$sec_plot[showi], T4$bal_acc[showi],
                 short_method(T4$method[showi]), pos = 4, cex = 0.58,
                 col = "grey25", offset = 0.3)
  graphics::title(main = "Accuracy vs computational cost", line = 2.4,
                  cex.main = 1.05)
  graphics::mtext("upper-left is better: high accuracy at low cost",
                  side = 3, line = 0.9, cex = 0.75, col = "grey30")
  graphics::legend("bottomright", c("proposed", "comparators"),
                   pch = 19, col = c(C_SF, "grey55"), bty = "n", cex = 0.75)
  graphics::par(op); pdf_close()
  cat("  F09 완료\n")
}

cat("\n=== P05 완료 ===\n->", file.path(IN, "figs"), "\n")
