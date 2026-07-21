# =====================================================================
#  P02_analysis_synthetic.R
#  ---------------------------------------------------------------
#  합성자료 심층분석.  대상은 제안 방법(SFCL/MFCL)의 kNN / kFDA 판본.
#
#  [3] A4  k-ablation (kNN) + gamma-ablation (kFDA)
#  [4] A5  노이즈 변수 추가에 대한 robustness
#  [5] A2  변수 추가 곡선 (HighDim20/50 은 분석 포함, 시각화만 분리)
#  [6] A1  ground-truth variable importance  -- TPR/FPR 제외, 연속 지표만
#  [7] A3  variable-selection: 선택 변수 수 + 선택된 k/gamma 함께 기록
#          (TPR/FPR 제외)
#      A7  거리(Euclid/Frechet) x 희소성(dense/sparse) 2x2
#
#  산출:
#    A1_importance.csv      변수별 w (fold 단위)
#    A1_summary.csv         설정별 연속 지표 (cor / AUC-of-w / L2)
#    A2_addcurve.csv        변수 추가 곡선
#    A3_sparsity.csv        s 스윕 + n_sel + k/gamma
#    A4_kablation.csv       kNN 의 k, kFDA 의 gamma
#    A5_noise.csv           노이즈 변수 추가
#    A7_ablation2x2.csv     2x2
# =====================================================================

ENV <- "main"
source("C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification/P00_common.R", encoding = "UTF-8")

SEEDS   <- 1:20
K_FOLDS <- 5L
OUT <- mk_out("R02_analysis_synth")

cat("\n=== P02: 합성 심층분석 ===\n")
cat("  변수선택 대상 설정:", length(SETTINGS_VARSEL), "\n")
cat("  출력:", OUT, "\n\n")


# ---------------------------------------------------------------------
# 공통: 한 (설정, seed) 의 캐시 + fold
# ---------------------------------------------------------------------
.setup <- function(dn, sd) {
  ds   <- flc_get_dataset(dn, seed = sd)
  prep <- prep_for(dn, sd, "rule", ds = ds, cache = FALSE)
  list(ds = ds, ca = prep$caches[[1]], y = ds$y,
       fold = strat_folds(ds$y, K_FOLDS, sd))
}

.learn_w <- function(ca, y, tr, s_select = "cv", s = NULL, seed = 1L) {
  flc_learn_weights(ca$Dvar[tr, tr, , drop = FALSE], y[tr],
                    s_select = s_select, s = s, seed = seed)
}


# =====================================================================
# A1 : 변수 중요도  (연속 지표만.  TPR/FPR 제외)
#
#   정답은 informative 이진이지만, 학습된 w 는 연속이다.
#   임계값에 의존하는 TPR/FPR 대신 임계값과 무관한 지표를 쓴다:
#     cor_pearson / cor_spearman : w 와 정답 지시자의 상관
#     auc_w   : w 를 점수로, informative 를 정답으로 본 Mann-Whitney AUC
#               (1 이면 모든 정보변수의 w 가 모든 노이즈변수보다 크다)
#     gap     : 정보변수 평균 w  -  노이즈변수 평균 w
#     w_ratio : 정보변수 w 합 / 전체 w 합
# =====================================================================
.auc_w <- function(w, info) {
  x <- w[info]; y <- w[!info]
  if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) /
    (length(x) * length(y))
}

run_A1 <- function(datasets = SETTINGS_VARSEL) {
  rows <- list()
  for (dn in datasets) {
    for (sd in SEEDS) {
      S <- .setup(dn, sd)
      info <- as_logical_safe(S$ds$informative)
      if (!any(info) || all(info)) next
      for (f in sort(unique(S$fold))) {
        tr <- which(S$fold != f)
        W  <- .learn_w(S$ca, S$y, tr, "cv", seed = sd)
        ord <- order(W$w, decreasing = TRUE)
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = dn, seed = sd, fold = f,
          variable = S$ds$varNames, w = W$w,
          rank = match(seq_along(W$w), ord),
          informative = info,
          stringsAsFactors = FALSE)
      }
    }
    cat("  A1", dn, "\n")
  }
  do.call(rbind, rows)
}

summarize_A1 <- function(a1) {
  if (is.null(a1) || !nrow(a1)) return(NULL)
  key <- paste(a1$dataset, a1$seed, a1$fold, sep = "|")
  sp  <- split(a1, key)
  rows <- lapply(sp, function(d) {
    info <- as_logical_safe(d$informative)
    if (!any(info) || all(info)) return(NULL)
    ind <- as.numeric(info)
    cp <- suppressWarnings(stats::cor(d$w, ind))
    cs <- suppressWarnings(stats::cor(d$w, ind, method = "spearman"))
    data.frame(dataset = d$dataset[1], seed = d$seed[1], fold = d$fold[1],
               p = nrow(d), p_info = sum(info),
               cor_pearson = cp, cor_spearman = cs,
               auc_w = .auc_w(d$w, info),
               gap = mean(d$w[info]) - mean(d$w[!info]),
               w_ratio = sum(d$w[info]) / max(sum(d$w), 1e-12),
               stringsAsFactors = FALSE)
  })
  raw <- do.call(rbind, rows)
  if (is.null(raw)) return(NULL)
  ag <- stats::aggregate(
    cbind(cor_pearson, cor_spearman, auc_w, gap, w_ratio) ~ dataset + p + p_info,
    raw, function(z) mean(z, na.rm = TRUE))
  ag[order(-ag$auc_w), ]
}


# =====================================================================
# A2 : 변수 추가 곡선
#   중요도 순으로 변수를 하나씩 넣으며 held-out 성능.
#   HighDim20/50 도 분석에는 포함한다 (시각화에서만 분리).
# =====================================================================
run_A2 <- function(datasets = SETTINGS_VARSEL) {
  rows <- list()
  for (dn in datasets) {
    for (sd in SEEDS) {
      S <- .setup(dn, sd); p <- S$ds$p
      for (f in sort(unique(S$fold))) {
        tr <- which(S$fold != f); te <- which(S$fold == f)
        W   <- .learn_w(S$ca, S$y, tr, "cv", seed = sd)
        ord <- order(W$w, decreasing = TRUE)
        info <- as_logical_safe(S$ds$informative)
        for (j in seq_len(p)) {
          wj <- numeric(p); wj[ord[seq_len(j)]] <- 1 / j
          D  <- cache_dist_w(S$ca, wj)
          kp <- knn_predict(D, S$y, tr, te)
          m  <- metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))
          rows[[length(rows) + 1L]] <- cbind(
            data.frame(dataset = dn, seed = sd, fold = f, n_var = j,
                       n_info_included = sum(info[ord[seq_len(j)]]),
                       k_sel = kp$k, stringsAsFactors = FALSE), m)
        }
      }
    }
    cat("  A2", dn, "\n")
  }
  do.call(rbind, rows)
}


# =====================================================================
# A3 : 희소성 s 스윕
#   요청 7: TPR/FPR 대신 n_sel + 선택된 k / gamma 를 함께 기록.
# =====================================================================
run_A3 <- function(datasets = SETTINGS_VARSEL) {
  rows <- list()
  for (dn in datasets) {
    for (sd in SEEDS) {
      S <- .setup(dn, sd); p <- S$ds$p
      sg <- exp(seq(log(1.02), log(sqrt(p)), length.out = S_GRID_LEN))
      info <- as_logical_safe(S$ds$informative)
      for (f in sort(unique(S$fold))) {
        tr <- which(S$fold != f); te <- which(S$fold == f)
        for (s in sg) {
          W <- .learn_w(S$ca, S$y, tr, "fixed", s = s, seed = sd)
          D <- cache_dist_w(S$ca, W$w)
          sel <- W$w > 1e-6
          kp <- knn_predict(D, S$y, tr, te)
          kf <- kfda_predict(D, S$y, tr, te, seed = sd)
          mk <- metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))
          mf <- metrics_row(S$y[te], kf$pred, kf$prob, levels(S$y))
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = dn, seed = sd, fold = f, s = s,
            n_sel = sum(sel), p = p,
            n_info_selected = sum(sel & info),
            info_share = if (sum(sel)) sum(sel & info) / sum(sel) else NA_real_,
            k_sel = kp$k, gamma_sel = kf$gamma_mult,
            bal_acc_knn = mk$bal_acc, acc_knn = mk$acc,
            macro_f1_knn = mk$macro_f1,
            bal_acc_kfda = mf$bal_acc, acc_kfda = mf$acc,
            macro_f1_kfda = mf$macro_f1,
            stringsAsFactors = FALSE)
        }
      }
    }
    cat("  A3", dn, "\n")
  }
  do.call(rbind, rows)
}


# =====================================================================
# A4 : k-ablation (kNN) + gamma-ablation (kFDA)
#
#   요청 3.  kFDA 에는 k 가 없다 (판별차원은 K-1 고정).  대응되는
#   하이퍼파라미터는 커널 폭 gamma 이므로, median heuristic 을 기준점으로
#   로그 등간격 배율 gamma_mult 를 훑는다:
#       gamma = gamma_mult / (2 * median(D)^2)
#   gamma_mult = 1 이 표준 median heuristic 이다.
# =====================================================================
run_A4 <- function(datasets = SETTINGS_ALL) {
  rows <- list()
  for (dn in datasets) {
    for (sd in SEEDS) {
      S <- .setup(dn, sd)
      for (f in sort(unique(S$fold))) {
        tr <- which(S$fold != f); te <- which(S$fold == f)
        W   <- .learn_w(S$ca, S$y, tr, "cv", seed = sd)
        Dsf <- cache_dist_w(S$ca, W$w)          # SFCL (sparse)
        Dmf <- S$ca$Djoint                       # MFCL (joint)

        #  --- kNN: k 훑기 ---
        for (k in K_GRID) {
          if (k >= length(tr)) next
          for (nm in c("SFCL", "MFCL")) {
            D <- if (nm == "SFCL") Dsf else Dmf
            kp <- knn_predict(D, S$y, tr, te, k = k)
            m  <- metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))
            rows[[length(rows) + 1L]] <- cbind(
              data.frame(dataset = dn, seed = sd, fold = f,
                         distance = nm, rule = "kNN",
                         param_name = "k", param = k,
                         stringsAsFactors = FALSE), m)
          }
        }

        #  --- kFDA: gamma_mult 훑기 ---
        for (gm in GAMMA_GRID) {
          for (nm in c("SFCL", "MFCL")) {
            D <- if (nm == "SFCL") Dsf else Dmf
            kf <- kfda_predict(D, S$y, tr, te, gamma_mult = gm, seed = sd)
            if (!isTRUE(kf$ok)) next
            m <- metrics_row(S$y[te], kf$pred, kf$prob, levels(S$y))
            rows[[length(rows) + 1L]] <- cbind(
              data.frame(dataset = dn, seed = sd, fold = f,
                         distance = nm, rule = "kFDA",
                         param_name = "gamma_mult", param = gm,
                         stringsAsFactors = FALSE), m)
          }
        }
      }
    }
    cat("  A4", dn, "\n")
  }
  do.call(rbind, rows)
}


# =====================================================================
# A5 : 노이즈 변수 추가에 대한 robustness
# =====================================================================
run_A5 <- function(bases = c("Setting1", "Setting3", "Shape"),
                   n_noise = c(0, 2, 4, 8, 16)) {
  rows <- list()
  for (bs in bases) {
    for (sd in SEEDS) {
      ds0 <- flc_get_dataset(bs, seed = sd)
      for (nn in n_noise) {
        ds <- if (nn > 0) flc_add_noise_vars(ds0, nn, seed = sd + 100L) else ds0
        prep <- flc_prepare(ds, lambda_method = "rule", verbose = FALSE)
        ca <- prep$caches[[1]]; y <- ds$y
        fold <- strat_folds(y, K_FOLDS, sd)
        for (f in sort(unique(fold))) {
          tr <- which(fold != f); te <- which(fold == f)
          W   <- .learn_w(ca, y, tr, "cv", seed = sd)
          Dsf <- cache_dist_w(ca, W$w)
          Dmf <- ca$Djoint
          for (nm in c("SFCL", "MFCL")) {
            D <- if (nm == "SFCL") Dsf else Dmf
            kp <- knn_predict(D, y, tr, te)
            m  <- metrics_row(y[te], kp$pred, kp$prob, levels(y))
            rows[[length(rows) + 1L]] <- cbind(
              data.frame(base = bs, seed = sd, fold = f, n_noise = nn,
                         distance = nm, p = ds$p,
                         n_sel = sum(W$w > 1e-6), k_sel = kp$k,
                         stringsAsFactors = FALSE), m)
          }
        }
        rm(prep); invisible(gc(verbose = FALSE))
      }
    }
    cat("  A5", bs, "\n")
  }
  do.call(rbind, rows)
}


# =====================================================================
# A7 : (Euclid vs Frechet) x (dense vs sparse) 2x2
#   비동기 설정은 Euclid 가 정의되지 않으므로 제외한다.
# =====================================================================
run_A7 <- function(datasets = SETTINGS_VARSEL) {
  rows <- list()
  for (dn in datasets) {
    ds1 <- flc_get_dataset(dn, seed = SEEDS[1])
    Ts  <- vapply(ds1$X, nrow, 0L)
    ragged <- length(unique(Ts)) > 1L ||
      any(vapply(ds1$X, function(m) anyNA(m), TRUE))
    if (ragged) {
      cat("  A7 skip", dn, "(비동기/ragged)\n"); next
    }
    for (sd in SEEDS) {
      S <- .setup(dn, sd); p <- S$ds$p
      R <- flc_regularize(S$ds)
      Xf <- t(apply(R$traj, 1, as.vector))
      for (f in sort(unique(S$fold))) {
        tr <- which(S$fold != f); te <- which(S$fold == f)
        W  <- .learn_w(S$ca, S$y, tr, "cv", seed = sd)
        wd <- rep(1, p)

        Dfs <- cache_dist_w(S$ca, W$w)
        Dfd <- cache_dist_w(S$ca, wd)
        Des <- .euclid_dist_w(R$traj, W$w)
        Ded <- .euclid_dist_w(R$traj, wd)

        gv <- function(D) {
          kp <- knn_predict(D, S$y, tr, te)
          metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))$bal_acc
        }
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = dn, seed = sd, fold = f,
          frechet_sparse = gv(Dfs), frechet_dense = gv(Dfd),
          euclid_sparse  = gv(Des), euclid_dense  = gv(Ded),
          n_sel = sum(W$w > 1e-6), p = p,
          stringsAsFactors = FALSE)
      }
    }
    cat("  A7", dn, "\n")
  }
  do.call(rbind, rows)
}

#  가중 유클리드 거리 ([n x T x p] -> n x n)
.euclid_dist_w <- function(traj, w) {
  n <- dim(traj)[1]; Tn <- dim(traj)[2]; p <- dim(traj)[3]
  A <- matrix(0, n, Tn * p)
  for (k in seq_len(p)) {
    A[, ((k - 1) * Tn + 1):(k * Tn)] <- traj[, , k] * sqrt(w[k])
  }
  sqrt(pmax(outer(rowSums(A^2), rowSums(A^2), `+`) - 2 * tcrossprod(A), 0))
}


# =====================================================================
# 실행
# =====================================================================
t0 <- Sys.time()

cat("\n[1/6] A1 변수 중요도 ...\n")
a1 <- run_A1()
wr(a1, file.path(OUT, "A1_importance.csv"))
s1 <- summarize_A1(a1)
wr(s1, file.path(OUT, "A1_summary.csv"))
if (!is.null(s1)) {
  cat("\n  A1 요약 (연속 지표):\n")
  print(s1, row.names = FALSE, digits = 3)
}

cat("\n[2/6] A2 변수 추가 곡선 ...\n")
a2 <- run_A2()
wr(a2, file.path(OUT, "A2_addcurve.csv"))

cat("\n[3/6] A3 희소성 스윕 ...\n")
a3 <- run_A3()
wr(a3, file.path(OUT, "A3_sparsity.csv"))

cat("\n[4/6] A4 k / gamma ablation ...\n")
a4 <- run_A4()
wr(a4, file.path(OUT, "A4_kablation.csv"))

cat("\n[5/6] A5 노이즈 robustness ...\n")
a5 <- run_A5()
wr(a5, file.path(OUT, "A5_noise.csv"))

cat("\n[6/6] A7 2x2 ablation ...\n")
a7 <- run_A7()
wr(a7, file.path(OUT, "A7_ablation2x2.csv"))

if (!is.null(a7) && nrow(a7)) {
  cat("\n  A7 쌍체 대비:\n")
  cn <- c("frechet_sparse", "frechet_dense", "euclid_sparse", "euclid_dense")
  cmp <- list(c("frechet_sparse", "euclid_sparse"),
              c("frechet_dense",  "euclid_dense"),
              c("frechet_sparse", "frechet_dense"),
              c("euclid_sparse",  "euclid_dense"))
  for (cc in cmp) {
    d <- a7[[cc[1]]] - a7[[cc[2]]]
    d <- d[is.finite(d)]
    s <- stats::sd(d)
    cat(sprintf("    %-32s %+.4f (t=%+.2f, win %d/%d)\n",
                paste(cc, collapse = " - "), mean(d),
                if (s > 0) mean(d) / (s / sqrt(length(d))) else NA_real_,
                sum(d > 0), length(d)))
  }
}

cat(sprintf("\n=== P02 완료 (%.1f h) ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
cat("->", OUT, "\n")
