# =====================================================================
#  P04_real.R
#  ---------------------------------------------------------------
#  실데이터 분석.  DATA 변수 하나로 대상을 고른다.
#
#    기존 3종 : Epilepsy, japanese_vowels, spoken_arabic_digit
#    신규 3종 : Heartbeat, FaceDetection_reduced, PEMS-SF
#    DATA <- "all3"  (기존)  |  "new3" (신규)  |  "all" (전부)
#                    |  개별 이름
#
#  합성과 같은 분석을 그대로 수행하되, 실데이터 고유 사항을 반영한다:
#    - informative 정답이 없으므로 A1 은 "w 프로파일 + 안정성"만 본다
#    - 격자 불규칙 자료는 grid=TRUE 방법을 제외한다 (methods_excluded.csv)
#    - 전 지표 기록.  K>2 면 AUC 는 NA.
#    - kNN 은 k, kFDA 는 gamma_mult 를 ablation 하고,
#      메인표에는 fold 안에서 고른 최적값만 올린다 (누출 없음).
#      전체 ablation 은 T_ablation_k.csv / T_ablation_gamma.csv 로 별도.
#
#  산출:
#    results.csv, weights.csv
#    T01_main_<metric>.csv     방법별 전 지표
#    T02_rank.csv
#    T03_timing.csv            [분석 A]
#    T04_grid_group.csv        [분석 C]
#    T05_paired.csv            주요 쌍체 대비
#    T_ablation_k.csv          kNN k-ablation
#    T_ablation_gamma.csv      kFDA gamma-ablation
#    A2_addcurve.csv, A3_sparsity.csv, A6_predictions.csv
#    methods_excluded.csv
# =====================================================================

ENV <- "main"
source("C:/Users/bassc/luinn27/연구/[OnGoing] Longitudinal Classification/Longitudinal Classification/P00_common.R", encoding = "UTF-8")

# ---------------------------------------------------------------------
# 0. 대상
# ---------------------------------------------------------------------
DATA <- "all3"          # "all3" | "new3" | "all" | 개별 이름

SEEDS   <- 1:20
K_FOLDS <- 5L
LAMBDA_METHOD <- "cv"
TIME_MODE <- "raw"       # [0,1] 정규화는 길이 정보를 파괴한다
SUBSAMPLE_SEED <- 2026L

#  자료별 클래스당 상한.  NULL 이면 전체.
#  Frechet 은 O(n^2 T^2 p), Dvar 은 [n x n x p] double 이므로
#  n 을 키우면 메모리가 제곱으로 는다.
MAX_PER_CLASS <- list(
  Epilepsy               = NULL,
  japanese_vowels        = NULL,
  spoken_arabic_digit    = 200L,   # 10 클래스 -> n=2000
  Heartbeat              = NULL,
  FaceDetection_reduced  = NULL,   # 전처리에서 이미 클래스당 250
  `PEMS-SF`              = NULL
)

GROUPS <- list(
  all3 = c("Epilepsy", "japanese_vowels", "spoken_arabic_digit"),
  new3 = c("Heartbeat", "FaceDetection_reduced", "PEMS-SF")
)
GROUPS$all <- c(GROUPS$all3, GROUPS$new3)

TARGETS <- if (DATA %in% names(GROUPS)) GROUPS[[DATA]] else DATA


# ---------------------------------------------------------------------
# 1. 로딩
#   [B1] 반드시 new_flcdata() 를 거친다 (class="flcdata" -> [.flcdata 동작)
#   [B3] key 를 ids 순서에 match() 로 재정렬 (레이블-궤적 대응 보존)
#   [B4] informative=NULL 을 rep(NA,p) 로 바꾸므로 다시 NULL 로
#   [속도] 개체별 스캔 대신 split() 한 번
# ---------------------------------------------------------------------
load_real <- function(name, root = DATA_ROOT,
                      max_per_class = NULL, seed = SUBSAMPLE_SEED) {
  f <- file.path(root, paste0(name, "_processed"), paste0(name, "_long.csv"))
  if (!file.exists(f)) stop("CSV 없음 (preprocess 먼저): ", f)
  cat("    읽는 중 ...\n")
  L <- utils::read.csv(f, stringsAsFactors = FALSE)
  cat(sprintf("    %d 행\n", nrow(L)))

  vnames <- unique(L$variable)
  ids    <- unique(L$id)
  key <- L[!duplicated(L$id), c("id", "split", "label")]
  key <- key[match(ids, key$id), ]

  if (!is.null(max_per_class)) {
    set.seed(seed)
    keep <- unlist(lapply(split(key$id, key$label), function(v) {
      if (length(v) > max_per_class) sample(v, max_per_class) else v
    }))
    ids <- ids[ids %in% keep]
    key <- key[match(ids, key$id), ]
    L   <- L[L$id %in% ids, ]
    cat(sprintf("    부분표본 n=%d\n", length(ids)))
  }
  stopifnot(identical(as.character(key$id), as.character(ids)))

  cat("    행렬 조립 ...\n")
  Lsp <- split(L[, c("time", "variable", "value")], factor(L$id, levels = ids))
  rm(L); invisible(gc(verbose = FALSE))

  X <- vector("list", length(ids)); tt <- vector("list", length(ids))
  for (k in seq_along(ids)) {
    di <- Lsp[[k]]
    ts <- sort(unique(di$time))
    M  <- matrix(NA_real_, length(ts), length(vnames),
                 dimnames = list(NULL, vnames))
    M[cbind(match(di$time, ts), match(di$variable, vnames))] <- di$value
    if (anyNA(M)) {
      for (j in seq_len(ncol(M))) {
        v <- M[, j]
        if (all(is.na(v))) {
          M[, j] <- 0
        } else if (anyNA(v)) {
          M[, j] <- stats::approx(seq_along(v), v, seq_along(v), rule = 2)$y
        }
      }
    }
    X[[k]]  <- M
    tt[[k]] <- as.numeric(ts)
    if (k %% 500L == 0L) cat(sprintf("      %d / %d\n", k, length(ids)))
  }
  rm(Lsp); invisible(gc(verbose = FALSE))

  ds <- new_flcdata(X = X, tt = tt, y = factor(key$label),
                    varNames = vnames, name = name, informative = NULL,
                    meta = list(split = key$split, id = ids))
  ds$informative <- NULL
  ds
}

is_regular <- function(ds, tol = 1e-8) {
  Ts <- vapply(ds$X, nrow, 1L)
  if (length(unique(Ts)) != 1L) return(FALSE)
  t1 <- ds$tt[[1]]
  all(vapply(ds$tt, function(v)
    length(v) == length(t1) && max(abs(v - t1)) < tol, TRUE))
}

methods_for <- function(ds, ALLM) {
  if (is_regular(ds)) {
    return(list(use = ALLM, drop = list(), regular = TRUE, reason = character(0)))
  }
  ng <- vapply(ALLM, function(m) isTRUE(m$grid), TRUE)
  im <- grepl("\\+ *mean$|\\+ *centroid$|\\+ *DBA$", names(ALLM))
  keep <- !ng & !im
  list(use = ALLM[keep], drop = ALLM[!keep], regular = FALSE,
       reason = ifelse(ng[!keep], "needs common grid",
                       "Frechet/DTW mean calls interpolation"))
}


# ---------------------------------------------------------------------
# 2. 자료 하나 처리
# ---------------------------------------------------------------------
run_dataset <- function(dn) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("=== ", dn, " ===\n", sep = "")
  OUT <- mk_out(paste0("R04_real_", gsub("[^A-Za-z0-9]", "_", dn)))
  t_all <- Sys.time()

  mpc <- MAX_PER_CLASS[[dn]]
  ds  <- load_real(dn, max_per_class = mpc)
  stopifnot(inherits(ds, "flcdata"))

  Ts  <- vapply(ds$X, nrow, 0L)
  reg <- is_regular(ds)
  chance <- 1 / nlevels(ds$y)
  tstr <- if (length(unique(Ts)) == 1L) {
    as.character(Ts[1])
  } else {
    sprintf("%d-%d", min(Ts), max(Ts))
  }
  cat(sprintf("  n=%d p=%d K=%d T=%s %s  chance=%.3f\n",
              ds$n, ds$p, nlevels(ds$y), tstr,
              if (reg) "(regular)" else "(irregular)", chance))
  cat(sprintf("  Dvar 예상 %.2f GB / 캐시\n", ds$n^2 * ds$p * 8 / 1e9))

  ALLM <- flc_methods_42()
  MS <- methods_for(ds, ALLM)
  METHODS <- MS$use
  cat(sprintf("  방법 %d개 사용, %d개 제외\n",
              length(METHODS), length(MS$drop)))
  if (length(MS$drop)) {
    wr(data.frame(dataset = dn, method = names(MS$drop),
                  family = vapply(MS$drop, function(m) m$family, ""),
                  reason = MS$reason, stringsAsFactors = FALSE),
       file.path(OUT, "methods_excluded.csv"))
  }

  # ---- 벤치 ---------------------------------------------------------
  cat("\n  [1/6] 벤치 ...\n")
  f_res <- file.path(OUT, "results.csv"); if (file.exists(f_res)) file.remove(f_res)
  f_wt  <- file.path(OUT, "weights.csv"); if (file.exists(f_wt)) file.remove(f_wt)

  for (sd in SEEDS) {
    cat(sprintf("    seed %d ", sd))
    prep <- flc_prepare(ds, lambda_method = LAMBDA_METHOD, verbose = FALSE)
    fold <- strat_folds(ds$y, K_FOLDS, sd)
    B <- flc_benchmark_one(ds, dn, sd, prep, fold, METHODS)
    wr_append(B$results, f_res); wr_append(B$weights, f_wt)
    rm(prep); invisible(gc(verbose = FALSE))
    cat(" ok\n")
  }

  res <- utils::read.csv(f_res, stringsAsFactors = FALSE)
  res <- res[res$ok %in% c(TRUE, "TRUE"), ]
  cat(sprintf("    성공 %d 행\n", nrow(res)))

  METRICS <- c("acc", "bal_acc", "macro_prec", "macro_rec", "macro_f1",
               "kappa", "logloss", "brier", "auc")
  for (mt in METRICS) {
    if (all(!is.finite(res[[mt]]))) next
    tb <- stats::aggregate(res[[mt]], by = list(method = res$method),
                           FUN = function(z) c(m = mean(z, na.rm = TRUE),
                                               s = stats::sd(z, na.rm = TRUE)))
    T <- data.frame(method = tb$method, mean = tb$x[, "m"], sd = tb$x[, "s"],
                    stringsAsFactors = FALSE)
    T <- T[order(-T$mean), ]
    wr(T, file.path(OUT, sprintf("T01_main_%s.csv", mt)))
  }

  R <- mean_ranks(res, "bal_acc", by = "seed")
  wr(R, file.path(OUT, "T02_rank.csv"))

  if ("seconds" %in% names(res)) {
    tb <- stats::aggregate(seconds ~ method, res,
                           function(z) c(m = mean(z, na.rm = TRUE),
                                         s = stats::median(z, na.rm = TRUE)))
    T3 <- data.frame(method = tb$method, sec_mean = tb$seconds[, "m"],
                     sec_median = tb$seconds[, "s"], stringsAsFactors = FALSE)
    T3 <- T3[order(-T3$sec_mean), ]
    wr(T3, file.path(OUT, "T03_timing.csv"))
  }

  if ("needs_grid" %in% names(res)) {
    res$grid_group <- ifelse(res$needs_grid %in% c(TRUE, "TRUE"),
                             "needs_grid", "no_grid")
    T4 <- stats::aggregate(bal_acc ~ grid_group, res,
                           function(z) c(m = mean(z, na.rm = TRUE),
                                         n = length(z)))
    T4 <- data.frame(grid_group = T4$grid_group,
                     bal_acc = T4$bal_acc[, "m"], n = T4$bal_acc[, "n"],
                     stringsAsFactors = FALSE)
    wr(T4, file.path(OUT, "T04_grid_group.csv"))
  }

  PAIRS <- list(
    c("SFKmL-C(dense) + kNN", "SFKmL-C + kNN"),
    c("SFKmL-C + kNN",        "MFKmL-C + kNN"),
    c("SFKmL-C + kNN",        "Euclid + kNN"),
    c("SFKmL-C + kNN",        "1NN-DTW"),
    c("SFKmL-C + kNN",        "DTW + kNN"),
    c("SFKmL-C + kFDA",       "DTW + kFDA"))
  rows <- lapply(PAIRS, function(pp) paired_summary(res, pp[1], pp[2]))
  wr(do.call(rbind, rows), file.path(OUT, "T05_paired.csv"))

  # ---- k / gamma ablation ------------------------------------------
  cat("  [2/6] k / gamma ablation ...\n")
  ab <- run_ablation_real(ds, dn)
  wr(ab$k,  file.path(OUT, "T_ablation_k.csv"))
  wr(ab$g,  file.path(OUT, "T_ablation_gamma.csv"))
  wr(ab$best, file.path(OUT, "T_ablation_best.csv"))

  # ---- A2 / A3 -------------------------------------------------------
  cat("  [3/6] A2 변수 추가 곡선 ...\n")
  wr(run_A2_real(ds, dn), file.path(OUT, "A2_addcurve.csv"))

  cat("  [4/6] A3 희소성 ...\n")
  wr(run_A3_real(ds, dn), file.path(OUT, "A3_sparsity.csv"))

  # ---- A1: w 프로파일 -----------------------------------------------
  cat("  [5/6] 변수 가중치 요약 ...\n")
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
      print(A1, row.names = FALSE, digits = 3)
    }
  }

  # ---- A6: 예측 저장 -------------------------------------------------
  cat("  [6/6] 예측 저장 ...\n")
  wr(run_A6_real(ds, dn), file.path(OUT, "A6_predictions.csv"))

  cat(sprintf("\n  === %s 완료 (%.1f h) -> %s\n", dn,
              as.numeric(difftime(Sys.time(), t_all, units = "hours")), OUT))
  invisible(TRUE)
}


# ---------------------------------------------------------------------
# 3. 실데이터용 하위 분석
# ---------------------------------------------------------------------
.setup_real <- function(ds, sd) {
  prep <- flc_prepare(ds, lambda_method = "rule", verbose = FALSE)
  list(ca = prep$caches[[1]], y = ds$y,
       fold = strat_folds(ds$y, K_FOLDS, sd), prep = prep)
}

run_ablation_real <- function(ds, dn) {
  krows <- list(); grows <- list()
  for (sd in SEEDS) {
    S <- .setup_real(ds, sd)
    for (f in sort(unique(S$fold))) {
      tr <- which(S$fold != f); te <- which(S$fold == f)
      W  <- flc_learn_weights(S$ca$Dvar[tr, tr, , drop = FALSE], S$y[tr],
                              s_select = "cv", seed = sd)
      Dsf <- cache_dist_w(S$ca, W$w)
      Dmf <- S$ca$Djoint
      for (nm in c("SFCL", "MFCL")) {
        D <- if (nm == "SFCL") Dsf else Dmf
        for (k in K_GRID) {
          if (k >= length(tr)) next
          kp <- knn_predict(D, S$y, tr, te, k = k)
          m  <- metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))
          krows[[length(krows) + 1L]] <- cbind(
            data.frame(dataset = dn, seed = sd, fold = f, distance = nm,
                       k = k, stringsAsFactors = FALSE), m)
        }
        for (gm in GAMMA_GRID) {
          kf <- kfda_predict(D, S$y, tr, te, gamma_mult = gm, seed = sd)
          if (!isTRUE(kf$ok)) next
          m <- metrics_row(S$y[te], kf$pred, kf$prob, levels(S$y))
          grows[[length(grows) + 1L]] <- cbind(
            data.frame(dataset = dn, seed = sd, fold = f, distance = nm,
                       gamma_mult = gm, stringsAsFactors = FALSE), m)
        }
      }
    }
    rm(S); invisible(gc(verbose = FALSE))
  }
  kt <- do.call(rbind, krows); gt <- do.call(rbind, grows)

  #  메인표용 최적값 (ablation 평균 기준.  누출을 피하려면 fold 내부 선택을
  #  써야 하므로, 메인 results.csv 의 k_sel/gamma_sel 이 정식 결과이고
  #  여기 best 는 "사후 최적이 어디였나"를 보는 진단용이다.)
  best <- NULL
  if (!is.null(kt)) {
    a <- stats::aggregate(bal_acc ~ distance + k, kt,
                          function(z) mean(z, na.rm = TRUE))
    b <- do.call(rbind, lapply(split(a, a$distance), function(d)
      d[which.max(d$bal_acc), ]))
    b$rule <- "kNN"; names(b)[names(b) == "k"] <- "param"
    best <- rbind(best, b[, c("distance", "rule", "param", "bal_acc")])
  }
  if (!is.null(gt)) {
    a <- stats::aggregate(bal_acc ~ distance + gamma_mult, gt,
                          function(z) mean(z, na.rm = TRUE))
    b <- do.call(rbind, lapply(split(a, a$distance), function(d)
      d[which.max(d$bal_acc), ]))
    b$rule <- "kFDA"; names(b)[names(b) == "gamma_mult"] <- "param"
    best <- rbind(best, b[, c("distance", "rule", "param", "bal_acc")])
  }
  list(k = kt, g = gt, best = best)
}

run_A2_real <- function(ds, dn) {
  rows <- list()
  for (sd in SEEDS) {
    S <- .setup_real(ds, sd); p <- ds$p
    for (f in sort(unique(S$fold))) {
      tr <- which(S$fold != f); te <- which(S$fold == f)
      W   <- flc_learn_weights(S$ca$Dvar[tr, tr, , drop = FALSE], S$y[tr],
                               s_select = "cv", seed = sd)
      ord <- order(W$w, decreasing = TRUE)
      for (j in seq_len(p)) {
        wj <- numeric(p); wj[ord[seq_len(j)]] <- 1 / j
        D  <- cache_dist_w(S$ca, wj)
        kp <- knn_predict(D, S$y, tr, te)
        m  <- metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))
        rows[[length(rows) + 1L]] <- cbind(
          data.frame(dataset = dn, seed = sd, fold = f, n_var = j,
                     k_sel = kp$k, stringsAsFactors = FALSE), m)
      }
    }
    rm(S); invisible(gc(verbose = FALSE))
  }
  do.call(rbind, rows)
}

run_A3_real <- function(ds, dn) {
  rows <- list(); p <- ds$p
  sg <- exp(seq(log(1.02), log(sqrt(p)), length.out = S_GRID_LEN))
  for (sd in SEEDS) {
    S <- .setup_real(ds, sd)
    for (f in sort(unique(S$fold))) {
      tr <- which(S$fold != f); te <- which(S$fold == f)
      for (s in sg) {
        W <- flc_learn_weights(S$ca$Dvar[tr, tr, , drop = FALSE], S$y[tr],
                               s_select = "fixed", s = s, seed = sd)
        D <- cache_dist_w(S$ca, W$w)
        kp <- knn_predict(D, S$y, tr, te)
        kf <- kfda_predict(D, S$y, tr, te, seed = sd)
        mk <- metrics_row(S$y[te], kp$pred, kp$prob, levels(S$y))
        mf <- metrics_row(S$y[te], kf$pred, kf$prob, levels(S$y))
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = dn, seed = sd, fold = f, s = s,
          n_sel = sum(W$w > 1e-6), p = p,
          k_sel = kp$k, gamma_sel = kf$gamma_mult,
          bal_acc_knn = mk$bal_acc, acc_knn = mk$acc, macro_f1_knn = mk$macro_f1,
          bal_acc_kfda = mf$bal_acc, acc_kfda = mf$acc,
          macro_f1_kfda = mf$macro_f1,
          stringsAsFactors = FALSE)
      }
    }
    rm(S); invisible(gc(verbose = FALSE))
  }
  do.call(rbind, rows)
}

run_A6_real <- function(ds, dn) {
  rows <- list(); cls <- levels(ds$y)
  sd <- SEEDS[1]
  S <- .setup_real(ds, sd)
  for (f in sort(unique(S$fold))) {
    tr <- which(S$fold != f); te <- which(S$fold == f)
    W  <- flc_learn_weights(S$ca$Dvar[tr, tr, , drop = FALSE], S$y[tr],
                            s_select = "cv", seed = sd)
    for (nm in c("SFCL", "MFCL")) {
      D <- if (nm == "SFCL") cache_dist_w(S$ca, W$w) else S$ca$Djoint
      kp <- knn_predict(D, S$y, tr, te)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = dn, method = nm, seed = sd, fold = f, sample = te,
        true = as.character(S$y[te]), pred = as.character(kp$pred),
        correct = as.character(S$y[te]) == as.character(kp$pred),
        k = kp$k, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


# ---------------------------------------------------------------------
# 4. 실행
# ---------------------------------------------------------------------
cat("\n=== P04: 실데이터 ===\n")
cat("  대상:", paste(TARGETS, collapse = ", "), "\n")

for (dn in TARGETS) {
  r <- try(run_dataset(dn), silent = FALSE)
  if (inherits(r, "try-error")) {
    cat("  (!) 실패:", dn, "-", conditionMessage(attr(r, "condition")), "\n")
  }
  prep_cache_clear()
}

cat("\n=== P04 완료 ===\n")
