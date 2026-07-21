# =====================================================================
# 06_metrics.R
#
# The previous round reported accuracy and macro-F1 only.  With the
# imbalanced setting in play (60/20/10), accuracy is actively misleading --
# predicting the majority class everywhere already gives 0.67.  So every
# classifier now returns CLASS PROBABILITIES and we score:
#
#   acc          overall accuracy
#   bal_acc      balanced accuracy == macro-recall  (the honest headline)
#   macro_prec / macro_rec / macro_f1
#   micro_f1     (== acc for single-label multiclass; kept for completeness)
#   auc          multiclass AUC, Hand & Till (2001): average pairwise AUC.
#                Reduces to the usual AUC when K = 2.
#   kappa        Cohen's kappa
#   logloss      multiclass cross-entropy (calibration)
#   brier        multiclass Brier score
#
# Methods that are not naturally probabilistic (medoid / centroid rules)
# expose a softmax over negative distances.  That is a monotone map, so it
# does not change acc/F1, and it gives a usable ranking for AUC -- which is
# all AUC needs.  Log-loss/Brier for those are reported but should be read
# as "ranking quality", not calibrated probability.
# =====================================================================

flc_confusion <- function(truth, pred, lv = levels(truth)) {
  table(factor(truth, levels = lv), factor(pred, levels = lv),
        dnn = c("truth", "pred"))
}

# ---- binary AUC (Mann-Whitney U), ties handled by mid-ranks -----------
.auc_binary <- function(score, pos) {
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- Hand & Till (2001) multiclass AUC --------------------------------
# M = (1/(K(K-1))) * sum_{i != j} A(i|j),  where A(i|j) uses P(class i) to
# separate class i from class j.  Symmetrised by averaging A(i|j), A(j|i).
auc_hand_till <- function(truth, prob, lv = levels(truth)) {
  K <- length(lv)
  if (K < 2) return(NA_real_)
  if (is.null(colnames(prob))) colnames(prob) <- lv
  prob <- prob[, lv, drop = FALSE]
  if (K == 2L) {
    return(.auc_binary(prob[, 2], truth == lv[2]))
  }
  acc <- c()
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    sel <- truth %in% c(lv[i], lv[j])
    if (!any(sel)) next
    aij <- .auc_binary(prob[sel, i], truth[sel] == lv[i])   # i vs j via P(i)
    aji <- .auc_binary(prob[sel, j], truth[sel] == lv[j])   # j vs i via P(j)
    acc <- c(acc, mean(c(aij, aji), na.rm = TRUE))
  }
  mean(acc, na.rm = TRUE)
}

flc_metrics <- function(truth, pred, prob = NULL, lv = levels(factor(truth))) {
  truth <- factor(truth, levels = lv)
  pred  <- factor(pred,  levels = lv)
  n <- length(truth); K <- length(lv)

  cm <- flc_confusion(truth, pred, lv)
  acc <- sum(diag(cm)) / n

  prec <- rec <- f1 <- numeric(K)
  for (k in seq_len(K)) {
    tp <- cm[k, k]; fp <- sum(cm[, k]) - tp; fn <- sum(cm[k, ]) - tp
    prec[k] <- if (tp + fp == 0) 0 else tp / (tp + fp)
    rec[k]  <- if (tp + fn == 0) 0 else tp / (tp + fn)
    f1[k]   <- if (prec[k] + rec[k] == 0) 0 else 2 * prec[k] * rec[k] / (prec[k] + rec[k])
  }

  # Cohen's kappa
  pe <- sum(rowSums(cm) * colSums(cm)) / n^2
  kappa <- if (pe == 1) NA_real_ else (acc - pe) / (1 - pe)

  out <- list(
    acc        = acc,
    bal_acc    = mean(rec),          # == macro-recall
    macro_prec = mean(prec),
    macro_rec  = mean(rec),
    macro_f1   = mean(f1),
    micro_f1   = acc,
    kappa      = kappa,
    auc        = NA_real_,
    logloss    = NA_real_,
    brier      = NA_real_,
    per_class_recall = setNames(rec, lv),
    confusion  = cm
  )

  if (!is.null(prob)) {
    prob <- as.matrix(prob)
    if (is.null(colnames(prob)) || !all(lv %in% colnames(prob)))
      colnames(prob) <- lv
    prob <- prob[, lv, drop = FALSE]
    prob[!is.finite(prob)] <- 1 / K
    prob <- pmax(prob, 0)
    rs <- rowSums(prob); rs[rs <= 0] <- 1
    prob <- prob / rs

    out$auc <- auc_hand_till(truth, prob, lv)

    idx <- cbind(seq_len(n), as.integer(truth))
    out$logloss <- -mean(log(pmax(prob[idx], 1e-15)))

    Y <- matrix(0, n, K); Y[idx] <- 1
    out$brier <- mean(rowSums((prob - Y)^2))
  }
  out
}

# flatten to a one-row data.frame (for rbind into the results CSV)
flc_metrics_row <- function(m) {
  data.frame(
    acc = m$acc, bal_acc = m$bal_acc, macro_prec = m$macro_prec,
    macro_rec = m$macro_rec, macro_f1 = m$macro_f1,
    auc = m$auc, kappa = m$kappa, logloss = m$logloss, brier = m$brier,
    stringsAsFactors = FALSE)
}

# softmax over negative distances -- the probability surrogate for
# prototype/medoid rules.  `tau` rescales by the median distance so the
# temperature is scale free.
dist_to_prob <- function(D, tau = 1) {
  s <- stats::median(D[is.finite(D) & D > 0])
  if (!is.finite(s) || s <= 0) s <- 1
  Z <- -D / (tau * s)
  Z <- Z - apply(Z, 1, max)
  E <- exp(Z)
  E / rowSums(E)
}

# ---- stratified folds -------------------------------------------------
strat_folds <- function(y, K = 5L, seed = 1L) {
  set.seed(seed)
  y <- factor(y); fold <- integer(length(y))
  for (c in levels(y)) {
    idx <- which(y == c)
    idx <- sample(idx)
    fold[idx] <- rep(seq_len(K), length.out = length(idx))
  }
  fold
}
