# =====================================================================
#  P07_plot_lambda.R   --  P03 (lambda 선택) 결과 작도
#  ---------------------------------------------------------------
#  입력:  R03_lambda/*.csv
#  출력:  R03_lambda/figs/*.pdf
#
#  H01_lambda_heat     설정 x 방법 x lambda 규칙 히트맵
#  H02_lambda_paired   rule 대비 cv / align 의 차이
#  H03_lambda_box      lambda 규칙별 분포
# =====================================================================

ROOT <- "C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification"
RESULT_ROOT <- ROOT
setwd(ROOT)
source("P0X_plot_common.R", encoding = "UTF-8")

IN <- file.path(RESULT_ROOT, "R03_lambda")
cat("\n=== P07: lambda 작도 ===\n입력:", IN, "\n")

LM_COL <- c(rule = "grey55", cv = C_SF, align = C_MF)


# ---------------------------------------------------------------------
# H01. 설정 x 방법 히트맵 (규칙마다 한 페이지)
# ---------------------------------------------------------------------
T <- rd(IN, "T_lambda.csv")
if (!is.null(T)) {
  lms <- intersect(c("rule", "cv", "align"), names(T))
  if (length(lms)) {
    pdf_open(IN, "H01_lambda_heat.pdf", w = 7.2, h = 5.2)
    for (lm in lms) {
      w <- stats::reshape(T[, c("dataset", "method", lm)],
                          idvar = "method", timevar = "dataset",
                          direction = "wide")
      M <- as.matrix(w[, -1, drop = FALSE])
      colnames(M) <- sub("^.*\\.", "", colnames(M))
      rownames(M) <- short_method(w$method)
      plot_heat(M, title = sprintf("lambda selection: %s", lm),
                sub = "bal_acc; same colour scale within each page")
    }
    pdf_close()
    cat("  H01 완료 (", length(lms), "페이지)\n")
  }
}


# ---------------------------------------------------------------------
# H02. rule 대비 차이
# ---------------------------------------------------------------------
P <- rd(IN, "T_lambda_paired.csv")
if (!is.null(P)) {
  P$lab <- sprintf("%s | %s", P$dataset, short_method(P$method))
  ctr <- sort(unique(P$contrast))
  pdf_open(IN, "H02_lambda_paired.pdf", w = 7.4, h = 6.2)
  for (cc in ctr) {
    d <- P[P$contrast == cc, ]
    if (!nrow(d)) next
    plot_paired_dots(d, "lab", diff_col = "diff", t_col = "t", n_col = "n",
                     title = sprintf("lambda: %s", cc),
                     sub = "positive means the supervised rule helps",
                     xlab = "difference in bal_acc")
  }
  pdf_close()
  cat("  H02 완료 (", length(ctr), "페이지)\n")
}


# ---------------------------------------------------------------------
# H03. 규칙별 분포
# ---------------------------------------------------------------------
res <- rd(IN, "results_lambda.csv")
if (!is.null(res)) {
  res <- res[res$ok %in% c(TRUE, "TRUE"), ]
  res$lambda_method <- factor(res$lambda_method,
                              levels = intersect(c("rule", "cv", "align"),
                                                 unique(res$lambda_method)))
  dss <- sort(unique(res$dataset))
  pdf_open(IN, "H03_lambda_box.pdf",
           w = min(12, 2.5 * length(dss) + 1.2), h = 5.0)
  op <- graphics::par(mar = c(6.2, 4.2, 3.8, 1.0), mgp = c(2.6, 0.7, 0))
  graphics::boxplot(bal_acc ~ lambda_method + dataset, data = res,
                    col = LM_COL[levels(res$lambda_method)],
                    border = "grey35", las = 2, cex.axis = 0.62,
                    xlab = "", ylab = "bal_acc", outpch = 1, outcex = 0.4,
                    whisklty = 1, main = "")
  graphics::title(main = "lambda selection rules", line = 2.4, cex.main = 1.05)
  graphics::mtext("differences are small and no rule wins consistently",
                  side = 3, line = 0.9, cex = 0.75, col = "grey30")
  graphics::legend("bottomright", levels(res$lambda_method),
                   fill = LM_COL[levels(res$lambda_method)],
                   border = "grey35", bty = "n", cex = 0.72)
  graphics::par(op); pdf_close()
  cat("  H03 완료\n")
}

cat("\n=== P07 완료 ===\n->", file.path(IN, "figs"), "\n")
