# =====================================================================
#  P0X_plot_common.R  --  모든 작도 스크립트가 source 한다.
#  ---------------------------------------------------------------
#  요청 7 반영: 축 라벨 / 제목 / 범례가 잘리지 않도록
#    - 여백(mar)을 라벨 길이에서 계산한다
#    - 긴 방법명은 줄여 쓰되 원문을 범례에 남긴다
#    - 로그축이 필요한 곳(gamma)은 로그 눈금을 명시한다
#    - 제목은 항상 한 줄, 부제는 mtext 로 분리
#
#  ! 이 파일만 따로 실행하지 않는다.
# =====================================================================

if (!exists("ROOT")) {
  ROOT <- getOption("FLC_ROOT", getwd())
}

# ---- 색 (Okabe-Ito, colorblind-safe) --------------------------------
PAL   <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
           "#56B4E9", "#F0E442", "#000000")
C_SF  <- "#0072B2"     # SFCL / Frechet / sparse
C_MF  <- "#D55E00"     # MFCL / joint
C_EU  <- "#7A7A7A"     # Euclid
C_DTW <- "#009E73"     # DTW
C_NOI <- "grey55"      # noise

fade <- function(col, a) {
  r <- grDevices::col2rgb(col)
  grDevices::rgb(r[1], r[2], r[3], alpha = a * 255, maxColorValue = 255)
}

pal_n <- function(n) {
  if (n <= length(PAL)) return(PAL[seq_len(n)])
  grDevices::hcl.colors(n, "Dark3")
}

# ---- CSV 읽기 (없으면 NULL, 경고만) ---------------------------------
rd <- function(dir, f) {
  p <- file.path(dir, f)
  if (!file.exists(p)) {
    message("  (없음) ", f); return(NULL)
  }
  x <- utils::read.csv(p, stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  x
}

# ---- PDF 열기/닫기 --------------------------------------------------
pdf_open <- function(dir, name, w = 7, h = 5) {
  d <- file.path(dir, "figs")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(file.path(d, name), width = w, height = h)
  invisible(file.path(d, name))
}
pdf_close <- function() grDevices::dev.off()

# ---- 방법명 축약 (축 라벨이 잘리지 않게) -----------------------------
#   원문은 유지하되 표시용으로만 줄인다.
short_method <- function(x) {
  x <- gsub("SFKmL-C", "SFCL", x, fixed = TRUE)
  x <- gsub("MFKmL-C", "MFCL", x, fixed = TRUE)
  x <- gsub("KmL3d-C \\(NC-Euclid\\)", "KmL3d-C", x)
  x <- gsub("DistSpace-", "DS-", x, fixed = TRUE)
  x <- gsub("MaxDepth-", "MD-", x, fixed = TRUE)
  x <- gsub("\\(dense\\)", "(d)", x)
  x <- gsub("\\(gaplab\\)", "(gl)", x)
  x <- gsub("DD-classifier", "DD-clf", x, fixed = TRUE)
  x <- gsub(" \\+ ", "+", x)
  x
}

#  왼쪽 여백을 라벨 길이에서 계산 (문자 폭 기준)
mar_left_for <- function(labels, cex = 0.72, base = 1.2, maxw = 16) {
  w <- max(nchar(labels), 1)
  min(maxw, base + w * cex * 0.52)
}

# ---- 공통 작도: 수평 boxplot (방법 비교) -----------------------------
#   metric 을 방법별로 그린다.  제안 방법은 파랑, 나머지는 회색.
plot_method_box <- function(df, value = "bal_acc", top_n = NULL,
                            title = "", sub = "", chance = NA,
                            highlight = "SFCL|MFCL") {
  if (is.null(df) || !nrow(df)) return(invisible())
  df <- df[is.finite(df[[value]]), ]
  ag <- stats::aggregate(df[[value]], by = list(method = df$method),
                         FUN = function(z) mean(z, na.rm = TRUE))
  names(ag)[2] <- "m"
  ag <- ag[order(-ag$m), ]
  if (!is.null(top_n) && nrow(ag) > top_n) ag <- ag[seq_len(top_n), ]
  d <- df[df$method %in% ag$method, ]
  lab <- short_method(ag$method)
  d$lab <- factor(short_method(d$method), levels = rev(lab))
  isp <- grepl(highlight, levels(d$lab))
  bcol <- ifelse(isp, fade(C_SF, 0.62), fade("grey60", 0.55))

  op <- graphics::par(mar = c(4.2, mar_left_for(lab), 3.4, 1.0),
                      mgp = c(2.5, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)
  graphics::boxplot(d[[value]] ~ d$lab, horizontal = TRUE, las = 1,
                    col = bcol, border = "grey35", cex.axis = 0.72,
                    xlab = value, ylab = "", outpch = 1, outcex = 0.45,
                    whisklty = 1, main = "")
  graphics::title(main = title, line = 2.0, cex.main = 1.05)
  if (nzchar(sub)) {
    graphics::mtext(sub, side = 3, line = 0.5, cex = 0.78, col = "grey30")
  }
  if (is.finite(chance)) {
    graphics::abline(v = chance, lty = 3, col = "grey65")
    graphics::mtext("chance", side = 1, at = chance, line = -1.1,
                    cex = 0.6, col = "grey50")
  }
  graphics::legend("bottomright",
                   c("proposed (SFCL/MFCL)", "comparators"),
                   fill = c(fade(C_SF, 0.62), fade("grey60", 0.55)),
                   border = "grey35", bty = "n", cex = 0.7)
  invisible()
}

# ---- 공통 작도: 선 그래프 (x 별 여러 계열) ---------------------------
plot_lines <- function(x, ymat, cols = NULL, pchs = NULL,
                       xlab = "", ylab = "", title = "", sub = "",
                       legend_pos = "bottomleft", ylim = NULL,
                       xaxis_at = NULL, xaxis_lab = NULL, logx = FALSE) {
  k <- ncol(ymat)
  if (is.null(cols)) cols <- pal_n(k)
  if (is.null(pchs)) pchs <- rep(c(19, 17, 15, 18), length.out = k)
  if (is.null(ylim)) {
    rg <- range(ymat, na.rm = TRUE)
    pad <- diff(rg) * 0.08; if (!is.finite(pad) || pad == 0) pad <- 0.02
    ylim <- c(rg[1] - pad, rg[2] + pad)
  }
  op <- graphics::par(mar = c(4.2, 4.2, 3.4, 1.0), mgp = c(2.5, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(NA, xlim = range(x), ylim = ylim, xlab = xlab, ylab = ylab,
                 main = "", xaxt = if (is.null(xaxis_at)) "s" else "n",
                 log = if (logx) "x" else "")
  if (!is.null(xaxis_at)) {
    graphics::axis(1, at = xaxis_at,
                   labels = if (is.null(xaxis_lab)) xaxis_at else xaxis_lab)
  }
  for (j in seq_len(k)) {
    graphics::lines(x, ymat[, j], type = "b", col = cols[j], pch = pchs[j],
                    lwd = 2.1, cex = 0.8)
  }
  graphics::title(main = title, line = 2.0, cex.main = 1.05)
  if (nzchar(sub)) {
    graphics::mtext(sub, side = 3, line = 0.5, cex = 0.78, col = "grey30")
  }
  if (!is.null(colnames(ymat))) {
    graphics::legend(legend_pos, colnames(ymat), col = cols, pch = pchs,
                     lwd = 2.1, bty = "n", cex = 0.75)
  }
  invisible()
}

# ---- 공통 작도: 히트맵 ------------------------------------------------
plot_heat <- function(M, title = "", sub = "", xlab = "", ylab = "",
                      digits = 3, cex_txt = 0.62) {
  if (is.null(M) || !length(M)) return(invisible())
  nr <- nrow(M); nc <- ncol(M)
  op <- graphics::par(mar = c(4.6, mar_left_for(rownames(M), 0.78),
                              3.6, 1.4), mgp = c(2.6, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)
  pal <- grDevices::colorRampPalette(
    c("#F7F7F7", "#FDDBC7", "#92C5DE", "#1B7837"))(64)
  rg <- range(M, na.rm = TRUE)
  graphics::image(seq_len(nc), seq_len(nr), t(M[nr:1, , drop = FALSE]),
                  col = pal, zlim = rg, axes = FALSE, xlab = xlab, ylab = ylab)
  graphics::axis(1, seq_len(nc), colnames(M), tick = FALSE, las = 2,
                 cex.axis = 0.72)
  graphics::axis(2, seq_len(nr), rev(rownames(M)), tick = FALSE, las = 1,
                 cex.axis = 0.75)
  graphics::box(col = "grey70")
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      v <- M[nr - i + 1, j]
      if (!is.finite(v)) next
      graphics::text(j, i, formatC(v, format = "f", digits = digits),
                     cex = cex_txt,
                     col = if (v > mean(rg, na.rm = TRUE)) "white" else "grey15")
    }
  }
  graphics::title(main = title, line = 2.0, cex.main = 1.05)
  if (nzchar(sub)) {
    graphics::mtext(sub, side = 3, line = 0.5, cex = 0.78, col = "grey30")
  }
  invisible()
}

# ---- 공통 작도: 쌍체 차이 (dot + CI) ---------------------------------
plot_paired_dots <- function(df, label_col, diff_col = "diff",
                             t_col = "t", n_col = "n_pairs",
                             title = "", sub = "", xlab = "difference") {
  if (is.null(df) || !nrow(df)) return(invisible())
  df <- df[is.finite(df[[diff_col]]), ]
  df <- df[order(df[[diff_col]]), ]
  lab <- df[[label_col]]
  se <- rep(NA_real_, nrow(df))
  ok <- is.finite(df[[t_col]]) & df[[t_col]] != 0
  se[ok] <- abs(df[[diff_col]][ok] / df[[t_col]][ok])
  lo <- df[[diff_col]] - 1.96 * se
  hi <- df[[diff_col]] + 1.96 * se

  op <- graphics::par(mar = c(4.2, mar_left_for(lab, 0.72), 3.4, 1.0),
                      mgp = c(2.5, 0.7, 0))
  on.exit(graphics::par(op), add = TRUE)
  xr <- range(c(lo, hi, 0), na.rm = TRUE)
  graphics::plot(NA, xlim = xr, ylim = c(0.5, nrow(df) + 0.5),
                 yaxt = "n", xlab = xlab, ylab = "", main = "")
  graphics::abline(v = 0, lty = 2, col = "grey60")
  graphics::axis(2, seq_len(nrow(df)), lab, las = 1, tick = FALSE,
                 cex.axis = 0.72)
  cols <- ifelse(df[[diff_col]] > 0, C_SF, C_MF)
  graphics::segments(lo, seq_len(nrow(df)), hi, seq_len(nrow(df)),
                     col = fade("grey40", 0.8), lwd = 1.4)
  graphics::points(df[[diff_col]], seq_len(nrow(df)), pch = 19,
                   col = cols, cex = 1.0)
  graphics::title(main = title, line = 2.0, cex.main = 1.05)
  if (nzchar(sub)) {
    graphics::mtext(sub, side = 3, line = 0.5, cex = 0.78, col = "grey30")
  }
  invisible()
}

cat("plot_common 로드 완료\n")
