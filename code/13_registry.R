# =====================================================================
# 13_registry.R -- 벤치마크에 들어가는 모든 방법의 단일 목록.
#
# 각 항목:
#   name     결과표에 찍히는 이름
#   family   그룹 (표를 묶는 단위)
#   fit      fit function (ds, ctx, ...) -> object with predict()
#   needs    필요한 optional 패키지 (없으면 NA로 건너뛰고 절대 죽지 않음)
#   source   출처 논문
#   grid     TRUE면 공통격자 보간이 필요 (=> 비동기 setting에서 불리)
#
# ---------------------------------------------------------------------
# *** 개념 정리 -- 논문에 반드시 쓸 것 ***
#
# SFKmL-C / MFKmL-C 는 "하나의 방법"이 아니라 **계열(family)** 이다.
#
#   기여  = (1) 거리 + (2) 희소 가중치를 지도학습으로 옮긴 것
#   자유  = **결정규칙은 그 위에서 자유롭게 고를 수 있다**
#
# 따라서 이들은 flc_combo_methods() 의 격자 안에 전부 들어 있다:
#
#   SFKmL-C + medoid        <- 원 논문(K-medoids)에 가장 충실한 버전
#   SFKmL-C + kNN / SVM / kFDA / mean
#   SFKmL-C(dense) + *      <- 희소성 없는 대조군
#   MFKmL-C + mean          <- 원 논문(K-means)에 가장 충실한 버전
#   MFKmL-C + kNN / SVM / kFDA / medoid
#
# medoid 를 쓴 것은 원 논문이 K-medoids 였기 때문일 뿐 필연이 아니다.
# 실제로 kNN/SVM 이 medoid 보다 낫다 -- 이는 Kang et al. (2023, Sec. 3) 이
# 이미 지적한 medoid 의 약점(짧은 궤적이 중심으로 뽑히는 문제)이, 분류에서는
# 보고용 그림이 아니라 **결정규칙 자체**로 작동하기 때문이다.
#
# => "SFKmL-C 가 SFKmL-C + kNN 보다 낮다" 같은 말은 성립하지 않는다.
#    둘 다 SFKmL-C 이고, 결정규칙만 다르다.
# ---------------------------------------------------------------------
#
# "grid" 열을 명시하는 이유:
#   보간이 필요한 방법은 비동기/불규칙 데이터에서 없는 정보를 만들어내야 한다.
#   이것이 Frechet 계열의 존재 이유이므로, 결과표에서 이 열로 묶어
#   "보간 필요 vs 불필요"를 대비시킨다.
# =====================================================================

.m <- function(name, family, fit, needs = character(0), source = "",
               grid = FALSE, args = list())
  list(name = name, family = family, fit = fit, needs = needs,
       source = source, grid = grid, args = args)

FLC_METHODS <- list(

  # ================= 제안 계열 중 격자로 표현 안 되는 것 ==============
  # (s 를 고르는 *다른 기준*.  격자는 cv 와 dense 만 쓴다.)
  .m("SFKmL-C(gap) + medoid", "proposed", fit_SFKmLC_gap, character(0),
     "s 를 논문의 gap statistic 으로 선택 (Kang et al. 2023, Eq. 10)", FALSE),

  .m("SFKmL-C(gaplab) + medoid", "proposed", fit_SFKmLC_gaplab, character(0),
     "s 를 label-permutation gap 으로 선택 (변수별 p-value 부산물)", FALSE),

  # ================= 거리 기반 경쟁 방법 ==============================
  .m("KmL3d-C (NC-Euclid)", "distance", fit_nc_euclid, character(0),
     "Genolini et al. (2013) KmL3d -> supervised. 논문 자신의 비교대상.", TRUE),

  .m("1NN-Euclid", "distance", fit_1nn_euclid, character(0),
     "trivial floor", TRUE),

  .m("1NN-DTW", "distance", fit_1nn_dtw, character(0),
     "Bagnall et al. (2017): TSC 의 표준 기준선", FALSE),

  .m("1NN-DDTW", "distance", fit_1nn_ddtw, character(0),
     "Keogh & Pazzani derivative DTW", FALSE),

  .m("NC-DBA", "distance", fit_nc_dba, character(0),
     "Petitjean et al. (2014) DTW Barycenter Averaging. MFKmL-C 의 DTW 대응물.",
     FALSE),

  # ================= depth 기반 =======================================
  .m("MaxDepth-MBD", "depth", function(ds, ctx, ...)
       fit_maxdepth_mbd(ds, ctx, sparse = "uniform", ...), character(0),
     "Lopez-Pintado & Romo (2009) MBD; Ieva & Paganoni (2013) 다변량 확장", TRUE),

  .m("MaxDepth-MBD (w)", "depth", function(ds, ctx, ...)
       fit_maxdepth_mbd(ds, ctx, sparse = "cv", ...), character(0),
     "Ieva & Paganoni 가 p_k 를 'problem driven' 으로 남겨둔 자리에 SFKmL 의 BCSS 가중치를 꽂은 것.",
     TRUE),

  .m("MaxDepth-MFHD", "depth", fit_maxdepth_mfhd, "mrfDepth",
     "Claeskens, Hubert, Slaets & Vakili (2014) 다변량 함수형 halfspace depth",
     TRUE),

  .m("DD-classifier", "depth", fit_ddclf, character(0),
     "Li, Cuesta-Albertos & Liu (2012) DD-plot 다항 분리자", TRUE),

  .m("DistSpace-Frechet", "depth", function(ds, ctx, ...)
       fit_distspace(ds, ctx, variant = "frechet", sparse = "cv", ...),
     character(0),
     "Hubert, Rousseeuw & Segaert (2017) 프레임워크 + 우리의 가중 Frechet 거리",
     FALSE),

  .m("DistSpace-bagdist", "depth", function(ds, ctx, ...)
       fit_distspace(ds, ctx, variant = "bagdist", ...), "mrfDepth",
     "Hubert, Rousseeuw & Segaert (2017) bagdistance", TRUE),

  # ================= functional 통계 계열 ==============================
  .m("MFPCA-LDA", "functional", fit_mfpca_lda, "MASS",
     "다변량 FPCA 점수 + LDA", TRUE),

  .m("MFPCA-SVM", "functional", fit_mfpca_svm, "e1071",
     "다변량 FPCA 점수 + RBF SVM", TRUE),

  .m("MFPLS-DA", "functional", fit_mfpls, "pls",
     "Preda, Saporta & Leveder (2007); MFPLS (Stat.Comput. 2024, 34:5)", TRUE),

  .m("Subspace", "functional", fit_subspace, character(0),
     "Fukuda et al. (2023) 다변량 함수형 부분공간 분류", TRUE),

  .m("fda.usc-knn", "functional", fit_fdausc_knn, "fda.usc",
     "Febrero-Bande & Oviedo de la Fuente, classif.knn", TRUE),

  .m("fda.usc-kernel", "functional", fit_fdausc_kernel, "fda.usc",
     "classif.kernel (Nadaraya-Watson)", TRUE)

  # fda.usc-glm / fda.usc-depth 는 제외.
  #   classif.glm  : ldata 객체 규약이 맞지 않아 predict 가 실패한다.
  #   classif.depth: predict 메서드가 없다 (fit/predict 분리 구조와 충돌).
  #   두 방법 모두 이 논문의 논지에 필수적이지 않고, functional GLM 은
  #   MFPLS-DA 와, depth 는 MBD/MFHD/DD/DistSpace 와 역할이 겹친다.
)
names(FLC_METHODS) <- vapply(FLC_METHODS, `[[`, "", "name")


# =====================================================================
# 제외된 방법들 (2026-07).  코드는 살아 있고, 목록만 분리해 두었다.
#
# [제외 사유 1] "궤적 구조를 버리는" 방법.  본 연구가 보고자 하는 대상이 아니다.
#
#   첫 전체 실행에서 tsfeat-RF 가 19개 데이터셋 **전부** 1위였다
#   (평균순위 2.26; SFKmL-C 14.2).  Tier B 에서도 졌다:
#       Phase       : tsfeat-RF 0.856  vs  SFKmL-C 0.422
#       Shape_async : tsfeat-RF 0.967  vs  SFKmL-C 0.811
#
#   원인: .tsfeat() 의 특성 중 argmax / argmin 이 있다.  구버전 gen_phase() 는
#   클래스를 bump 위치로 정의했으므로 argmax 가 레이블을 **직접 인코딩**했다.
#   즉 "요약통계가 강하다"가 아니라 "특성 하나가 레이블을 누설한다"였다.
#   (gen_phase 는 이후 전면 재설계되었다 -- 04_gen_new.R 주석 참조.)
#
# [제외 사유 2] 앙상블.
#   ProximityForest (Lucas et al. 2019) 는 아예 **삭제**했다.  단일 분류기와
#   공정 비교가 되지 않으며, 앙상블까지 포함하면 비교 대상이 한도 끝도 없어진다.
#
# 되살리는 법:  METHODS <- flc_method_set("all_plus_excluded")
# 부록이나 리뷰어 반론 대응에 필요할 수 있으므로 삭제하지 말 것.
# =====================================================================
FLC_METHODS_EXCLUDED <- list(

  .m("ROCKET-ridge", "excluded", fit_rocket, "glmnet",
     "Dempster, Petitjean & Webb (2020). 제외. rocket_transform_cpp 버그 미수정.",
     TRUE),

  .m("catch22-RF", "excluded", fit_catch22_rf, c("Rcatch22", "randomForest"),
     "Lubba et al. (2019). 제외: 22개 요약통계로 곡선을 붕괴시킴.", TRUE),

  .m("tsfeat-RF", "excluded", fit_tsfeat_rf, "randomForest",
     "17개 요약통계 + RF. 제외: 곡선을 버림. argmax 가 phase 레이블을 누설.",
     TRUE),

  .m("tsfeat-SVM", "excluded", fit_tsfeat_svm, "e1071",
     "17개 요약통계 + RBF SVM. 제외: 위와 동일.", TRUE),

  .m("flat-RF", "excluded", fit_flat_rf, "randomForest",
     "값을 이어붙인 벡터 + RF. 제외: 모형이 시간 순서를 활용하지 않음.", TRUE),

  .m("flat-SVM", "excluded", fit_flat_svm, "e1071",
     "값을 이어붙인 벡터 + RBF SVM. 제외: 위와 동일.", TRUE),

  .m("flat-ridge", "excluded", fit_flat_ridge, "glmnet",
     "값을 이어붙인 벡터 + multinomial ridge. 제외: 위와 동일.", TRUE)
)
names(FLC_METHODS_EXCLUDED) <-
  vapply(FLC_METHODS_EXCLUDED, `[[`, "", "name")


# ---------------------------------------------------------------------
# 거리 x 희소성 x 분류기 격자를 레지스트리 항목으로 확장.
# SFKmL-C / MFKmL-C 계열 전체가 여기서 나온다.
# ---------------------------------------------------------------------
flc_combo_methods <- function(full = TRUE) {
  g <- flc_combo_grid(full)
  out <- lapply(seq_len(nrow(g)), function(i) {
    d <- g$distance[i]; sp <- g$sparsity[i]; cl <- g$clf[i]
    .m(name = g$name[i],
       family = if (d %in% c("sumvars", "joint")) "proposed" else "distance",
       fit = local({
         dd <- d; ss <- sp; cc <- cl
         function(ds, ctx, ...) fit_combo(ds, ctx, distance = dd,
                                          sparsity = ss, clf = cc, ...)
       }),
       needs = if (cl == "svm") "kernlab" else character(0),
       source = sprintf("%s distance + %s rule", d, cl),
       grid = (d == "euclid"))
  })
  names(out) <- g$name
  out
}


# ---------------------------------------------------------------------
# 실제로 돌릴 방법 목록을 조립한다.
#   which = "all" | "core" | "proposed" | "combo"
#           | "excluded" | "all_plus_excluded"
#           | 방법 이름 벡터
# ---------------------------------------------------------------------
flc_method_set <- function(which = c("all", "core", "proposed", "combo",
                                     "excluded", "all_plus_excluded"),
                           combo_full = FALSE) {
  .known <- c("all", "core", "proposed", "combo", "excluded", "all_plus_excluded")
  if (length(which) > 1L || !which[1] %in% .known) {
    pool <- c(FLC_METHODS, flc_combo_methods(TRUE), FLC_METHODS_EXCLUDED)
    return(pool[which])
  }
  which <- match.arg(which)
  cm <- flc_combo_methods(full = combo_full)

  switch(which,
    all               = c(FLC_METHODS, cm),
    excluded          = FLC_METHODS_EXCLUDED,
    all_plus_excluded = c(FLC_METHODS, cm, FLC_METHODS_EXCLUDED),
    combo             = cm,

    proposed = c(
      FLC_METHODS[vapply(FLC_METHODS, `[[`, "", "family") == "proposed"],
      cm[grepl("^(SFKmL|MFKmL)", names(cm))]),

    core = c(
      FLC_METHODS[c("SFKmL-C(gap) + medoid",
                    "KmL3d-C (NC-Euclid)", "1NN-DTW", "1NN-DDTW", "NC-DBA",
                    "MaxDepth-MBD", "DD-classifier", "DistSpace-Frechet",
                    "MFPCA-LDA", "MFPLS-DA", "Subspace")],
      cm[intersect(names(cm),
                   c("SFKmL-C + medoid", "SFKmL-C + kNN", "SFKmL-C + SVM",
                     "SFKmL-C(dense) + medoid",
                     "MFKmL-C + mean", "MFKmL-C + kNN",
                     "DTW + kNN", "Euclid + kNN"))])
  )
}


# 사용 가능 여부 점검 (실행 전에 무엇이 빠지는지 미리 보고)
flc_check_methods <- function(methods = flc_method_set("all")) {
  do.call(rbind, lapply(methods, function(m) {
    miss <- m$needs[!vapply(m$needs, flc_have, TRUE)]
    data.frame(name = m$name, family = m$family,
               needs = paste(m$needs, collapse = ","),
               missing = paste(miss, collapse = ","),
               available = length(miss) == 0L,
               grid = m$grid, source = m$source,
               stringsAsFactors = FALSE)
  }))
}
