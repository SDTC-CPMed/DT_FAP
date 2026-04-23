#============================================================
# PROJECT: UK Biobank Proteomics Risk Prediction Analysis Cox model (Machine Learning)
# DESCRIPTION:
# This script performs protein-based risk prediction for Colorectal Cancer using the UK Biobank proteomics data. 
# Similar analysis was performed adding clinical baseline information (age, sex, BMI) and a colorectal cancer-specific polygenic risk score (PRS_CRC) to evaluate the added value of proteomics beyond traditional risk factors.
# Performances were evaluated in the testing data at year 2, 4, 6, 8, 10, 12, 14 
# Cases that diagnosed at year 0, 2 and 4 were removed for reverse-causation sensitivity analyses
# Calibration / Brier / sensitivity were measured at year 2, 10 and 14 for early, late, and very late diagnosis, respectively, with fixed specificity at 0.80.
# Key steps include: data loading, case and control sub-setting, training and testing splitting, KNN imputation, Cox model (alpha=0) modeling in training, 
# KM curve, C-index, time-dependent ROC curve generation in test, all performed 1000 times for robust estimation.
#============================================================

# -----------------------------
# 0) Packages
# -----------------------------
packages <- c(
  "dplyr", "stringr", "tibble", "survival", "Matrix",
  "ggplot2", "data.table", "glmnet", "purrr", "broom",
  "tidyr", "foreach", "doRNG", "doParallel"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

if (!requireNamespace("survivalROC", quietly = TRUE)) install.packages("survivalROC")
library(survivalROC)

if (!requireNamespace("impute", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("impute")
}
library(impute)

if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")
library(survminer)

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

rm(list = setdiff(ls(), "UKBB_Proteomics"))

# -----------------------------
# 1) User parameters ####
# -----------------------------
dis2 <- "Colorectal_cancer"

set.seed(123)
n_bootstrap <- 1000
alpha <- 0
ncores <- 40

eval_years <- c(2, 4, 6, 8, 10, 12, 14) #2, 4, 6, 8, 10, 12, 14, 16
exclude_first_n_years_list <- c(0,2,4)

# fixed horizons for calibration / Brier / sensitivity
fixed_eval_horizons <- c(2,10,14)
fixed_specificity <- 0.80

control_groups <- c("HC")

known_marker <- c("CEACAM5")
signature <- unique(c(
  "LAMB1","LAMC2","TGFBI","SEMA4D","AGRN","SEMA5A","SEMA3F","INHBA",
  "LAMA5","COL7A1","SOX9","CCND1","TGFBI","MET","ID1","PTMA","CLU",
  "TGIF1","EDN1","GDF15"
))

# -----------------------------
# 2) Directories
# -----------------------------
input.dir  <- ".../UKBB/input"
input.dir2 <- ".../UKBB/output"
out.base   <- ".../UKBB/output"
out.dir    <- file.path(out.base, "UKBB_incidentCRC_survival_final")


dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3) Small helpers
# -----------------------------
save_csv_rds <- function(obj, file_stub) {
  saveRDS(obj, paste0(file_stub, ".rds"))
  write.csv(obj, paste0(file_stub, ".csv"), row.names = FALSE)
}

make_analysis_dirs <- function(base_dir) {
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_dir, "descriptive"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_dir, "km"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_dir, "predictions"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_dir, "fixed_metrics"), showWarnings = FALSE, recursive = TRUE)
}

make_analysis_dirs(out.dir)

get_model_colors <- function() {
  c(
    "FDS" = "#F4A261",
    "CEA" = "#5DA5DA",
    "Random" = "grey60",
    "age_sex_bmi_PRS" =    "#ADD8E6", 
    "FDS_age_sex_bmi_PRS" = "#C45A3A" 
  )
}

# -----------------------------
# 4) Load data
# -----------------------------
study_end_date <- as.Date("2022-12-30")
# UKBB_Proteomics = read.table(file = paste(input.dir, '/Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt', sep = ''), sep='\t', header = TRUE, fill = TRUE, row.names = 'PID')

All_clinic <- readRDS(paste0(input.dir2, "/ukb_clinic.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
  dplyr::mutate(across(
    c(date_of_attending_assessment_centre_f53_0_0,
      date_of_death_f40000_0_0,
      date_of_death_f40000_1_0,
      date_lost_to_followup_f191_0_0),
    as.Date
  )) %>%
  dplyr::mutate(
    sex_f31_0_0 = as.factor(sex_f31_0_0),
    ethnic_background_f21000_0_0 = as.factor(ethnic_background_f21000_0_0),
    OSstatus.raw = ifelse(is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0), 0, 1),
    date_of_last_record = dplyr::case_when(
      OSstatus.raw == 1 ~ pmax(date_of_death_f40000_0_0, date_of_death_f40000_1_0, na.rm = TRUE),
      is.na(date_lost_to_followup_f191_0_0) ~ study_end_date,
      TRUE ~ pmin(date_lost_to_followup_f191_0_0, study_end_date)
    ),
    followup_days = as.numeric(date_of_last_record - date_of_attending_assessment_centre_f53_0_0),
    followup_years = followup_days / 365
  ) %>%
  dplyr::filter(
    !is.na(ethnic_background_f21000_0_0),
    !is.na(age_at_recruitment_f21022_0_0),
    !is.na(uk_biobank_assessment_centre_f54_0_0),
    !is.na(sex_f31_0_0)
  ) %>%
  dplyr::mutate(eid = as.character(eid))

dis2_clinic <- readRDS(paste0(input.dir2, "/ukb_", dis2, "_subset_clinic.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
  dplyr::mutate(eid = as.character(eid))

All_cancer_clinic_eid <- readRDS(paste0(input.dir2, "/ukb_AllCancers_subset_clinic.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
  dplyr::pull(eid)

dis2_inci_clinic <- dis2_clinic %>%
  dplyr::filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0)

dis2_inci_eid <- dis2_inci_clinic %>% dplyr::pull(eid)
dis2_eid <- dis2_clinic %>% dplyr::pull(eid)

HC_clinic <- readRDS(paste0(input.dir2, "/ukb_healthy_controls.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics))
HC_clinic <- All_clinic %>% dplyr::filter(eid %in% as.character(HC_clinic$eid))
HC_eid <- HC_clinic %>% dplyr::pull(eid)

PRS_data_path <- '/data/sharedData/Dina_genetics/PRS_directory/Martin_pan_cancer_PRS/Martin_Cancer_PRSs/Matrix_scaled_CombinedPRS_EUROS_Martin_pancamcer_PRSs.csv'
PRS_data <- data.table::fread(PRS_data_path, data.table = FALSE) %>%
  dplyr::mutate(eid = as.character(eid)) %>%
  dplyr::select(eid, CRC_new_SCORE) %>%
  dplyr::rename(PRS_CRC = CRC_new_SCORE)

UKBB_Proteomics <- UKBB_Proteomics[, colSums(is.na(UKBB_Proteomics)) <= 0.25 * nrow(UKBB_Proteomics), drop = FALSE]
UKBB_Proteomics$eid <- rownames(UKBB_Proteomics)

exist_signature_proteins <- intersect(signature, colnames(UKBB_Proteomics))
all_other_proteins <- setdiff(colnames(UKBB_Proteomics), c(signature, "eid"))

# -----------------------------
# 5) Save descriptive outputs
# -----------------------------
save_descriptive_outputs <- function() {
  
  followup_df <- All_clinic %>%
    dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
    dplyr::select(eid, followup_years, sex_f31_0_0, age_at_recruitment_f21022_0_0)
  
  save_csv_rds(followup_df, file.path(out.dir, "descriptive", "followup_distribution_data"))
  
  followup_summary <- followup_df %>%
    dplyr::summarise(
      N = dplyr::n(),
      median_followup = median(followup_years, na.rm = TRUE),
      IQR_followup_L = quantile(followup_years, 0.25, na.rm = TRUE),
      IQR_followup_U = quantile(followup_years, 0.75, na.rm = TRUE),
      min_followup = min(followup_years, na.rm = TRUE),
      max_followup = max(followup_years, na.rm = TRUE)
    )
  
  write.csv(followup_summary, file.path(out.dir, "descriptive", "followup_summary.csv"), row.names = FALSE)
  
  p_followup <- ggplot(followup_df, aes(x = followup_years)) +
    geom_histogram(bins = 40, fill = "#5DA5DA", color = "white") +
    theme_classic() +
    labs(
      title = "Distribution of follow-up time",
      x = "Follow-up time (years)",
      y = "Number of participants"
    )
  
  ggsave(file.path(out.dir, "descriptive", "followup_distribution_histogram.pdf"),
         p_followup, width = 5.5, height = 4.0, dpi = 300)
  ggsave(file.path(out.dir, "descriptive", "followup_distribution_histogram.png"),
         p_followup, width = 5.5, height = 4.0, dpi = 300)
  
  delta_df <- dis2_inci_clinic %>%
    dplyr::mutate(
      delta_days = as.numeric(diagnosis_time - date_of_attending_assessment_centre_f53_0_0),
      delta_years = delta_days / 365
    ) %>%
    dplyr::filter(is.finite(delta_years), delta_years > 0) %>%
    dplyr::select(eid, diagnosis_time, date_of_attending_assessment_centre_f53_0_0, delta_days, delta_years)
  
  save_csv_rds(delta_df, file.path(out.dir, "descriptive", "incident_crc_delta_time_data"))
  
  delta_summary <- delta_df %>%
    dplyr::summarise(
      N = dplyr::n(),
      median_delta_years = median(delta_years, na.rm = TRUE),
      IQR_delta_L = quantile(delta_years, 0.25, na.rm = TRUE),
      IQR_delta_U = quantile(delta_years, 0.75, na.rm = TRUE),
      min_delta_years = min(delta_years, na.rm = TRUE),
      max_delta_years = max(delta_years, na.rm = TRUE),
      n_le_2y = sum(delta_years <= 2, na.rm = TRUE),
      n_le_4y = sum(delta_years <= 4, na.rm = TRUE),
      n_gt_10y = sum(delta_years > 10, na.rm = TRUE)
    )
  
  write.csv(delta_summary, file.path(out.dir, "descriptive", "incident_crc_delta_time_summary.csv"), row.names = FALSE)
  
  p_delta <- ggplot(delta_df, aes(x = delta_years)) +
    geom_histogram(bins = 35, fill = "#F4A261", color = "white") +
    theme_classic() +
    labs(
      title = "Time from blood sampling to CRC diagnosis",
      x = "Years from sampling to diagnosis",
      y = "Number of incident CRC cases"
    )
  
  ggsave(file.path(out.dir, "descriptive", "incident_crc_delta_time_histogram.pdf"),
         p_delta, width = 5.5, height = 4.0, dpi = 300)
  ggsave(file.path(out.dir, "descriptive", "incident_crc_delta_time_histogram.png"),
         p_delta, width = 5.5, height = 4.0, dpi = 300)
}

save_descriptive_outputs()

# -----------------------------
# 6) Model helpers####
# -----------------------------
get_model_spec <- function(model_name) {
  if (model_name == "FDS") {
    list(
      protein_set = intersect(signature, colnames(UKBB_Proteomics)),
      clinical_vars = character(0)
    )
  } else if (model_name == "CEA") {
    list(
      protein_set = intersect("CEACAM5", colnames(UKBB_Proteomics)),
      clinical_vars = character(0)
    )
  } else if (model_name == "Random") {
    length_random <- length(intersect(signature, colnames(UKBB_Proteomics)))
    protein_set_i <- if (length(all_other_proteins) < length_random) {
      all_other_proteins
    } else {
      sample(all_other_proteins, length_random)
    }
    list(
      protein_set = protein_set_i,
      clinical_vars = character(0)
    )
  } else if (model_name == "age_sex_bmi_PRS") {
    list(
      protein_set = character(0),
      clinical_vars = c("age", "sex", "bmi", "PRS_CRC")
    )
  } else if (model_name == "FDS_age_sex_bmi_PRS") {
    list(
      protein_set = intersect(signature, colnames(UKBB_Proteomics)),
      clinical_vars = c("age", "sex", "bmi", "PRS_CRC")
    )
  } else {
    stop(paste("Unknown model_name:", model_name))
  }
}

build_survival_dataset <- function(current_control = "HC",
                                   exclude_first_n_years = 0) {
  
  case_surv <- All_clinic %>%
    dplyr::filter(eid %in% dis2_inci_eid) %>%
    dplyr::left_join(UKBB_Proteomics, by = "eid") %>%
    dplyr::left_join(PRS_data, by = "eid") %>%
    dplyr::left_join(dis2_inci_clinic %>% dplyr::select(eid, diagnosis_time), by = "eid") %>%
    dplyr::mutate(
      delta_days = as.numeric(diagnosis_time - date_of_attending_assessment_centre_f53_0_0),
      observed_time_yrs = delta_days / 365,
      event = 1
    ) %>%
    dplyr::filter(observed_time_yrs > exclude_first_n_years) %>%
    dplyr::mutate(
      time_yrs = observed_time_yrs - exclude_first_n_years,
      binary_group = 1,
      age = as.numeric(age_at_recruitment_f21022_0_0),
      sex = ifelse(sex_f31_0_0 == "Female", 0, ifelse(sex_f31_0_0 == "Male", 1, NA)),
      bmi = as.numeric(body_mass_index_bmi_f21001_0_0)
    ) %>%
    dplyr::select(
      eid, binary_group, age, sex, bmi, PRS_CRC, event, time_yrs,
      dplyr::all_of(colnames(UKBB_Proteomics))
    )
  
  control_ids <- if (current_control == "HC") HC_eid else stop("Only HC implemented.")
  
  control_surv <- All_clinic %>%
    dplyr::filter(eid %in% control_ids) %>%
    dplyr::left_join(UKBB_Proteomics, by = "eid") %>%
    dplyr::left_join(PRS_data, by = "eid") %>%
    dplyr::mutate(
      observed_time_yrs = followup_years,
      event = 0
    ) %>%
    dplyr::filter(observed_time_yrs > exclude_first_n_years) %>%
    dplyr::mutate(
      time_yrs = observed_time_yrs - exclude_first_n_years,
      binary_group = 0,
      age = as.numeric(age_at_recruitment_f21022_0_0),
      sex = ifelse(sex_f31_0_0 == "Female", 0, ifelse(sex_f31_0_0 == "Male", 1, NA)),
      bmi = as.numeric(body_mass_index_bmi_f21001_0_0)
    ) %>%
    dplyr::select(
      eid, binary_group, age, sex, bmi, PRS_CRC, event, time_yrs,
      dplyr::all_of(colnames(UKBB_Proteomics))
    )
  
  dplyr::bind_rows(case_surv, control_surv) %>%
    dplyr::mutate(binary_group = factor(binary_group))
}

make_stratified_split <- function(full_data_all, train_frac = 0.7) {
  event_idx <- which(full_data_all$event == 1)
  nonevent_idx <- which(full_data_all$event == 0)
  
  train_event_idx <- sample(event_idx, round(train_frac * length(event_idx)))
  train_nonevent_idx <- sample(nonevent_idx, round(train_frac * length(nonevent_idx)))
  
  train_idx <- c(train_event_idx, train_nonevent_idx)
  test_idx <- setdiff(seq_len(nrow(full_data_all)), train_idx)
  
  list(
    train_data = full_data_all[train_idx, , drop = FALSE],
    test_data  = full_data_all[test_idx, , drop = FALSE]
  )
}

impute_train_test <- function(train_x_df, test_x_df, protein_set_i, clinical_vars) {
  
  if (length(protein_set_i) > 0) {
    protein_matrix_train <- as.matrix(train_x_df[, protein_set_i, drop = FALSE])
    protein_matrix_train <- apply(protein_matrix_train, 2, as.numeric)
    
    if (is.null(dim(protein_matrix_train))) {
      protein_matrix_train <- matrix(protein_matrix_train, ncol = 1)
      colnames(protein_matrix_train) <- protein_set_i
    }
    
    imputed_train <- tryCatch(
      impute::impute.knn(protein_matrix_train)$data,
      error = function(e) protein_matrix_train
    )
    
    if (is.null(dim(imputed_train))) {
      imputed_train <- matrix(imputed_train, ncol = 1)
      colnames(imputed_train) <- protein_set_i
    }
    
    train_x_df[, protein_set_i] <- imputed_train
    
    for (col in protein_set_i) {
      train_mean <- mean(train_x_df[[col]], na.rm = TRUE)
      if (is.na(train_mean)) train_mean <- 0
      train_x_df[[col]][is.na(train_x_df[[col]])] <- train_mean
      test_x_df[[col]][is.na(test_x_df[[col]])] <- train_mean
    }
  }
  
  if (length(clinical_vars) > 0) {
    for (col in clinical_vars) {
      train_mean <- mean(train_x_df[[col]], na.rm = TRUE)
      if (is.na(train_mean)) train_mean <- 0
      train_x_df[[col]][is.na(train_x_df[[col]])] <- train_mean
      test_x_df[[col]][is.na(test_x_df[[col]])] <- train_mean
    }
  }
  
  list(train_x_df = train_x_df, test_x_df = test_x_df)
}

fit_model_and_predict <- function(x_train, x_test, train_data, test_data, alpha = 0) {
  y_train <- survival::Surv(train_data$time_yrs, train_data$event)
  
  if (ncol(x_train) == 1) {
    train_tmp <- train_data
    train_tmp$marker <- as.numeric(x_train[, 1])
    
    fit <- tryCatch(
      survival::coxph(survival::Surv(time_yrs, event) ~ marker, data = train_tmp),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    
    train_lp <- tryCatch(as.numeric(stats::predict(fit, type = "lp")),
                         error = function(e) rep(NA_real_, nrow(train_data)))
    
    test_tmp <- test_data
    test_tmp$marker <- as.numeric(x_test[, 1])
    test_lp <- tryCatch(as.numeric(stats::predict(fit, newdata = test_tmp, type = "lp")),
                        error = function(e) rep(NA_real_, nrow(test_data)))
    
  } else {
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x = x_train,
        y = y_train,
        family = "cox",
        alpha = alpha,
        standardize = TRUE,
        nfolds = 10
      ),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) return(NULL)
    
    train_lp <- tryCatch(
      as.numeric(predict(cv_fit, newx = x_train, s = "lambda.min", type = "link")),
      error = function(e) rep(NA_real_, nrow(train_data))
    )
    test_lp <- tryCatch(
      as.numeric(predict(cv_fit, newx = x_test, s = "lambda.min", type = "link")),
      error = function(e) rep(NA_real_, nrow(test_data))
    )
  }
  
  if (all(is.na(test_lp)) || all(is.na(train_lp)) || is.na(sd(test_lp, na.rm = TRUE)) || sd(test_lp, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  list(train_lp = train_lp, test_lp = test_lp)
}

calc_cindex <- function(test_data, lp) {
  tryCatch({
    cidx_df <- test_data %>%
      dplyr::mutate(lp = lp) %>%
      dplyr::select(time_yrs, event, lp) %>%
      dplyr::filter(!is.na(time_yrs), !is.na(event), !is.na(lp), is.finite(lp))
    
    if (nrow(cidx_df) < 2 || length(unique(cidx_df$event)) < 2 || sd(cidx_df$lp, na.rm = TRUE) == 0) {
      NA_real_
    } else {
      survival::concordance(
        survival::Surv(time_yrs, event) ~ lp,
        data = cidx_df,
        reverse = TRUE
      )$concordance
    }
  }, error = function(e) NA_real_)
}

# -----------------------------
# 7) Fixed-horizon helpers
# -----------------------------
get_H0_at_t <- function(basehaz_df, t0) {
  if (is.null(basehaz_df) || nrow(basehaz_df) == 0) return(NA_real_)
  if (t0 <= min(basehaz_df$time, na.rm = TRUE)) return(0)
  idx <- which(basehaz_df$time <= t0)
  if (length(idx) == 0) return(0)
  max(basehaz_df$hazard[idx], na.rm = TRUE)
}

get_evaluable_horizon_df <- function(df, risk, t0) {
  df %>%
    dplyr::mutate(pred_risk = risk) %>%
    dplyr::filter(
      !is.na(time_yrs), !is.na(event), !is.na(pred_risk),
      (event == 1 & time_yrs <= t0) | (time_yrs >= t0)
    ) %>%
    dplyr::mutate(event_by_t0 = ifelse(event == 1 & time_yrs <= t0, 1, 0))
}

fit_censor_km <- function(time, event) {
  censor_status <- 1 - event
  tryCatch(
    survival::survfit(survival::Surv(time, censor_status) ~ 1),
    error = function(e) NULL
  )
}

get_Ghat <- function(km_fit, t) {
  if (is.null(km_fit)) return(NA_real_)
  s <- tryCatch(summary(km_fit, times = t, extend = TRUE), error = function(e) NULL)
  if (is.null(s) || length(s$surv) == 0) return(NA_real_)
  g <- as.numeric(s$surv[1])
  if (!is.finite(g)) return(NA_real_)
  g
}

calc_ipcw_brier_score <- function(time, event, pred_risk, t0, eps = 1e-6) {
  keep <- is.finite(time) & is.finite(event) & is.finite(pred_risk)
  time <- time[keep]
  event <- event[keep]
  pred_risk <- pred_risk[keep]
  
  if (length(time) == 0) return(NA_real_)
  
  km_cens <- fit_censor_km(time = time, event = event)
  if (is.null(km_cens)) return(NA_real_)
  
  G_t0 <- get_Ghat(km_cens, t0)
  if (is.na(G_t0) || G_t0 <= eps) return(NA_real_)
  
  event_before_t0 <- (time <= t0 & event == 1)
  G_Tminus <- rep(NA_real_, length(time))
  
  if (any(event_before_t0)) {
    t_left <- pmax(time[event_before_t0] - eps, 0)
    G_Tminus[event_before_t0] <- vapply(t_left, function(tt) get_Ghat(km_cens, tt), numeric(1))
  }
  
  survive_past_t0 <- (time > t0)
  w1 <- rep(0, length(time))
  w2 <- rep(0, length(time))
  
  w1[event_before_t0] <- 1 / pmax(G_Tminus[event_before_t0], eps)
  w2[survive_past_t0] <- 1 / pmax(G_t0, eps)
  
  brier_i <- w1 * (1 - pred_risk)^2 + w2 * (0 - pred_risk)^2
  mean(brier_i, na.rm = TRUE)
}

calc_null_and_scaled_brier <- function(train_time, train_event, test_time, test_event, pred_risk, t0, eps = 1e-6) {
  keep_test <- is.finite(test_time) & is.finite(test_event) & is.finite(pred_risk)
  test_time <- test_time[keep_test]
  test_event <- test_event[keep_test]
  pred_risk <- pred_risk[keep_test]
  
  keep_train <- is.finite(train_time) & is.finite(train_event)
  train_time <- train_time[keep_train]
  train_event <- train_event[keep_train]
  
  if (length(test_time) == 0 || length(train_time) == 0) {
    return(list(brier_model = NA_real_, brier_null = NA_real_, scaled_brier = NA_real_, null_risk = NA_real_))
  }
  
  km_train <- tryCatch(survival::survfit(survival::Surv(train_time, train_event) ~ 1), error = function(e) NULL)
  if (is.null(km_train)) {
    return(list(brier_model = NA_real_, brier_null = NA_real_, scaled_brier = NA_real_, null_risk = NA_real_))
  }
  
  s_train <- tryCatch(summary(km_train, times = t0, extend = TRUE), error = function(e) NULL)
  if (is.null(s_train) || length(s_train$surv) == 0) {
    return(list(brier_model = NA_real_, brier_null = NA_real_, scaled_brier = NA_real_, null_risk = NA_real_))
  }
  
  null_risk <- 1 - as.numeric(s_train$surv[1])
  
  brier_model <- calc_ipcw_brier_score(test_time, test_event, pred_risk, t0, eps)
  brier_null  <- calc_ipcw_brier_score(test_time, test_event, rep(null_risk, length(test_time)), t0, eps)
  
  scaled_brier <- ifelse(
    is.na(brier_model) || is.na(brier_null) || brier_null <= eps,
    NA_real_,
    1 - (brier_model / brier_null)
  )
  
  list(
    brier_model = brier_model,
    brier_null = brier_null,
    scaled_brier = scaled_brier,
    null_risk = null_risk
  )
}

calc_sensitivity_at_fixed_specificity <- function(event_by_t0, pred_risk, fixed_specificity = 0.90) {
  if (length(unique(event_by_t0)) < 2) {
    return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
  }
  
  roc_obj <- tryCatch(
    pROC::roc(response = event_by_t0, predictor = pred_risk, quiet = TRUE, direction = "<"),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) {
    return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
  }
  
  coords_df <- tryCatch(
    pROC::coords(
      roc_obj,
      x = roc_obj$thresholds,
      input = "threshold",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(coords_df) || nrow(coords_df) == 0) {
    return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
  }
  
  coords_df <- as.data.frame(coords_df)
  coords_df$dist <- abs(coords_df$specificity - fixed_specificity)
  best <- coords_df[which.min(coords_df$dist), , drop = FALSE]
  
  list(
    threshold = best$threshold[1],
    sensitivity = best$sensitivity[1],
    specificity = best$specificity[1]
  )
}

make_calibration_df <- function(eval_df, n_bins = 10) {
  if (nrow(eval_df) < n_bins) return(NULL)
  
  eval_df %>%
    dplyr::mutate(risk_bin = dplyr::ntile(pred_risk, n_bins)) %>%
    dplyr::group_by(risk_bin) %>%
    dplyr::summarise(
      mean_pred = mean(pred_risk, na.rm = TRUE),
      obs_rate = mean(event_by_t0, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

# -----------------------------
# 8) Plotting helpers
# -----------------------------
plot_calibration_curve <- function(cal_df, model_name, current_control, horizon, out.dir) {
  if (is.null(cal_df) || nrow(cal_df) == 0) return(NULL)
  
  p <- ggplot(cal_df, aes(x = mean_pred, y = obs_rate)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    coord_equal(
      xlim = c(0, max(0.15, max(cal_df$mean_pred, na.rm = TRUE))),
      ylim = c(0, max(0.15, max(cal_df$obs_rate, na.rm = TRUE)))
    ) +
    labs(
      title = paste0("Calibration at ", horizon, "-year risk"),
      subtitle = paste0(model_name, " | ", current_control, " | Bootstrap 1"),
      x = "Mean predicted risk",
      y = "Observed event rate"
    ) +
    theme_classic()
  
  ggsave(file.path(out.dir, paste0("Calibration_", model_name, "_", current_control, "_", horizon, "y_bootstrap1.pdf")),
         p, width = 4.5, height = 4.2, dpi = 300)
  ggsave(file.path(out.dir, paste0("Calibration_", model_name, "_", current_control, "_", horizon, "y_bootstrap1.png")),
         p, width = 4.5, height = 4.2, dpi = 300)
  
  invisible(p)
}

plot_combined_km_first_run <- function(models = c("FDS", "CEA", "Random"),
                                       current_control = "HC",
                                       out.dir) {
  km_list <- lapply(models, function(m) {
    f <- file.path(out.dir, "km", paste0("KMdata_", m, "_", current_control, "_bootstrap1.rds"))
    if (!file.exists(f)) return(NULL)
    x <- readRDS(f)
    x$model <- m
    x
  })
  km_list <- km_list[!vapply(km_list, is.null, logical(1))]
  if (length(km_list) == 0) return(NULL)
  
  surv_df <- dplyr::bind_rows(lapply(km_list, function(test_data) {
    fit <- survival::survfit(survival::Surv(time_yrs, event) ~ risk_group, data = test_data)
    s <- survminer::surv_summary(fit, data = test_data)
    s$model <- unique(test_data$model)
    s$event_prob <- 1 - s$surv
    s
  }))
  
  surv_df$model <- factor(surv_df$model, levels = models)
  
  p <- ggplot(surv_df, aes(x = time, y = event_prob, color = strata)) +
    geom_step(linewidth = 1) +
    geom_point(data = surv_df %>% dplyr::filter(n.censor > 0), size = 1.1) +
    facet_wrap(~ model, nrow = 1) +
    coord_cartesian(ylim = c(0, 0.30)) +
    scale_color_manual(
      values = c("risk_group=Low risk" = "#1F77B4", "risk_group=High risk" = "#a12728")
    ) +
    labs(
      title = paste0("Held-out test-set KM curves (Bootstrap 1, ", current_control, ")"),
      subtitle = "Risk groups defined by the median training-set linear predictor; points indicate censoring",
      x = "Years since baseline",
      y = "Probability of CRC",
      color = "Predicted risk"
    ) +
    theme_classic()
  
  ggsave(file.path(out.dir, paste0("KM_combined_FDS_CEA_Random_", current_control, "_bootstrap1_withCensor.pdf")),
         p, width = 9.5, height = 3.7, dpi = 300)
  ggsave(file.path(out.dir, paste0("KM_combined_FDS_CEA_Random_", current_control, "_bootstrap1_withCensor.png")),
         p, width = 9.5, height = 3.7, dpi = 300)
  
  invisible(p)
}

# -----------------------------
# 9) Main analysis function
# -----------------------------
analyze_survival_models_incident <- function(
    model_name,
    controls = c("HC"),
    eval_years = c(2,4,6,8,10,12,14),
    n_bootstrap = 100,
    alpha = 0,
    ncores = 20,
    exclude_first_n_years = 0,
    fixed_eval_horizons = c(5, 14),
    fixed_specificity = 0.90,
    out.dir = NULL) {
  
  purrr::map_dfr(controls, function(current_control) {
    
    message(paste("Working: Model:", model_name, "| Control:", current_control,
                  "| Exclude first:", exclude_first_n_years, "years"))
    
    full_data_all <- build_survival_dataset(
      current_control = current_control,
      exclude_first_n_years = exclude_first_n_years
    )
    
    if (nrow(full_data_all) < 50 || sum(full_data_all$event == 1, na.rm = TRUE) < 20) return(NULL)
    
    doParallel::registerDoParallel(cores = ncores)
    
    foreach::foreach(
      i = 1:n_bootstrap, .combine = rbind,
      .packages = c("dplyr","glmnet","survival","impute","survivalROC","tibble","pROC","tidyr","purrr","survminer")
    ) %dopar% {
      
      set.seed(1000 + i)
      
      spec <- get_model_spec(model_name)
      protein_set_i <- spec$protein_set
      clinical_vars <- spec$clinical_vars
      selected_vars <- c(protein_set_i, clinical_vars)
      
      if (length(selected_vars) == 0) return(NULL)
      
      split_obj <- make_stratified_split(full_data_all, train_frac = 0.7)
      train_data <- split_obj$train_data
      test_data  <- split_obj$test_data
      
      if (nrow(train_data) < 20 || nrow(test_data) < 10) return(NULL)
      if (sum(train_data$event == 1, na.rm = TRUE) < 10 || sum(test_data$event == 1, na.rm = TRUE) < 5) return(NULL)
      
      train_x_df <- train_data[, selected_vars, drop = FALSE]
      test_x_df  <- test_data[, selected_vars, drop = FALSE]
      
      imp_obj <- impute_train_test(train_x_df, test_x_df, protein_set_i, clinical_vars)
      train_x_df <- imp_obj$train_x_df
      test_x_df  <- imp_obj$test_x_df
      
      x_train <- as.matrix(train_x_df)
      x_test  <- as.matrix(test_x_df)
      
      if (is.null(dim(x_train))) {
        x_train <- matrix(x_train, ncol = 1)
        colnames(x_train) <- selected_vars
      }
      if (is.null(dim(x_test))) {
        x_test <- matrix(x_test, ncol = 1)
        colnames(x_test) <- selected_vars
      }
      
      mode(x_train) <- "numeric"
      mode(x_test)  <- "numeric"
      
      fit_obj <- fit_model_and_predict(x_train, x_test, train_data, test_data, alpha = alpha)
      if (is.null(fit_obj)) return(NULL)
      
      train_data$lp <- fit_obj$train_lp
      test_data$lp  <- fit_obj$test_lp
      
      # save KM data for bootstrap 1
      if (!is.null(out.dir) && i == 1) {
        km_cutoff <- median(fit_obj$train_lp, na.rm = TRUE)
        km_df <- test_data %>%
          dplyr::mutate(
            risk_group = ifelse(lp >= km_cutoff, "High risk", "Low risk"),
            model = model_name,
            Bootstrap = i,
            ExcludeFirstYears = exclude_first_n_years,
            Control_Group = current_control
          ) %>%
          dplyr::select(
            eid, time_yrs, event, lp, risk_group, model, Bootstrap,
            ExcludeFirstYears, Control_Group
          )
        
        saveRDS(km_df,
                file.path(out.dir, "km", paste0("KMdata_", model_name, "_", current_control, "_bootstrap1.rds")))
      }
      
      c_index <- calc_cindex(test_data, fit_obj$test_lp)
      
      # tdAUC rows
      valid_years <- eval_years[eval_years > exclude_first_n_years]
      auc_rows <- purrr::map_dfr(valid_years, function(tt) {
        horizon_adj <- tt - exclude_first_n_years
        
        roc_obj <- tryCatch(
          survivalROC::survivalROC(
            Stime = test_data$time_yrs,
            status = test_data$event,
            marker = fit_obj$test_lp,
            predict.time = horizon_adj,
            method = "KM"
          ),
          error = function(e) NULL
        )
        
        tibble::tibble(
          Group = "dis2_inci_Cox",
          # N_case = sum(full_data_all$event == 1, na.rm = TRUE),
          N_used_Protein = length(protein_set_i),
          Selected_Proteins = ifelse(length(protein_set_i) == 0, NA, paste(protein_set_i, collapse = ";")),
          Selected_Clinical = ifelse(length(clinical_vars) == 0, NA, paste(clinical_vars, collapse = ";")),
          C_index = c_index,
          tdAUC = if (is.null(roc_obj)) NA_real_ else roc_obj$AUC,
          Model = model_name,
          Bootstrap = i,
          Observation_Years = tt,
          Control_Group = current_control,
          ExcludeFirstYears = exclude_first_n_years
        )
      })
      
      # fixed-horizon outputs
      cal_fit <- tryCatch(
        survival::coxph(survival::Surv(time_yrs, event) ~ lp, data = train_data),
        error = function(e) NULL
      )
      
      if (!is.null(cal_fit)) {
        bh_df <- tryCatch(survival::basehaz(cal_fit, centered = FALSE), error = function(e) NULL)
        
        if (!is.null(bh_df) && nrow(bh_df) > 0) {
          fixed_metric_rows <- list()
          
          for (t0 in fixed_eval_horizons) {
            if (t0 <= exclude_first_n_years) next
            t0_adj <- t0 - exclude_first_n_years
            
            H0_t0 <- get_H0_at_t(bh_df, t0_adj)
            if (is.na(H0_t0)) next
            
            pred_risk_t0 <- 1 - exp(-H0_t0 * exp(fit_obj$test_lp))
            
            eval_df_t0 <- get_evaluable_horizon_df(test_data, pred_risk_t0, t0_adj) %>%
              dplyr::mutate(
                Model = model_name,
                Bootstrap = i,
                HorizonYears = t0,
                HorizonAdjusted = t0_adj,
                ExcludeFirstYears = exclude_first_n_years,
                Control_Group = current_control
              )
            
            brier_obj <- calc_null_and_scaled_brier(
              train_time = train_data$time_yrs,
              train_event = train_data$event,
              test_time = test_data$time_yrs,
              test_event = test_data$event,
              pred_risk = pred_risk_t0,
              t0 = t0_adj
            )
            
            sens_info <- if (nrow(eval_df_t0) > 10 && length(unique(eval_df_t0$event_by_t0)) > 1) {
              calc_sensitivity_at_fixed_specificity(
                event_by_t0 = eval_df_t0$event_by_t0,
                pred_risk = eval_df_t0$pred_risk,
                fixed_specificity = fixed_specificity
              )
            } else {
              list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_)
            }
            
            fixed_metric_rows[[length(fixed_metric_rows) + 1]] <- tibble::tibble(
              Model = model_name,
              Bootstrap = i,
              Control_Group = current_control,
              ExcludeFirstYears = exclude_first_n_years,
              HorizonYears = t0,
              HorizonAdjusted = t0_adj,
              Brier_IPCW = brier_obj$brier_model,
              Brier_Null = brier_obj$brier_null,
              Scaled_Brier = brier_obj$scaled_brier,
              Null_Risk = brier_obj$null_risk,
              Sensitivity_at_fixed_spec = sens_info$sensitivity,
              Achieved_Specificity = sens_info$specificity,
              Threshold_at_fixed_spec = sens_info$threshold,
              Spec_target = fixed_specificity,
              N_eval = nrow(eval_df_t0)
            )
            
            if (!is.null(out.dir) && nrow(eval_df_t0) > 0) {
              saveRDS(
                eval_df_t0 %>%
                  dplyr::select(
                    eid, time_yrs, event, event_by_t0, pred_risk,
                    Model, Bootstrap, HorizonYears, HorizonAdjusted,
                    ExcludeFirstYears, Control_Group
                  ),
                file.path(out.dir, "predictions",
                          paste0("PredRisk_", model_name, "_", current_control,
                                 "_bt", i, "_", t0, "y_excl", exclude_first_n_years, ".rds"))
              )
            }
            
            if (!is.null(out.dir) &&
                model_name == "FDS_age_sex_bmi_PRS" &&
                i == 1 &&
                nrow(eval_df_t0) >= 20 &&
                length(unique(eval_df_t0$event_by_t0)) > 1) {
              
              cal_df <- make_calibration_df(eval_df_t0, n_bins = 10)
              if (!is.null(cal_df)) {
                saveRDS(
                  cal_df,
                  file.path(out.dir, paste0("CalibrationData_", model_name, "_", current_control, "_", t0, "y_bootstrap1.rds"))
                )
              }
            }
          }
          
          if (!is.null(out.dir) && length(fixed_metric_rows) > 0) {
            saveRDS(
              dplyr::bind_rows(fixed_metric_rows),
              file.path(out.dir, "fixed_metrics",
                        paste0("FixedMetrics_", model_name, "_", current_control,
                               "_bt", i, "_excl", exclude_first_n_years, ".rds"))
            )
          }
        }
      }
      
      auc_rows
    }
  })
}

# -----------------------------
# 10) Run analyses ####
# -----------------------------
runs_list <- list(
  list(model_name = "FDS", controls = c("HC"), alpha = alpha),
  list(model_name = "CEA", controls = c("HC"), alpha = alpha),
  list(model_name = "Random", controls = c("HC"), alpha = alpha),
  list(model_name = "age_sex_bmi_PRS", controls = c("HC"), alpha = alpha),
  list(model_name = "FDS_age_sex_bmi_PRS", controls = c("HC"), alpha = alpha)
)

start_time <- Sys.time()
all_intermediate_files <- c()

for (exclude_first_n_years in exclude_first_n_years_list) {
  
  out.dir.sens <- file.path(out.dir, paste0("excludeFirst_", exclude_first_n_years, "y"))
  make_analysis_dirs(out.dir.sens)
  
  intermediate_files <- c()
  
  for (i in seq_along(runs_list)) {
    run <- runs_list[[i]]

    res <- analyze_survival_models_incident(
      model_name = run$model_name,
      controls = run$controls,
      eval_years = eval_years,
      n_bootstrap = n_bootstrap,
      alpha = run$alpha,
      ncores = ncores,
      exclude_first_n_years = exclude_first_n_years,
      fixed_eval_horizons = fixed_eval_horizons,
      fixed_specificity = fixed_specificity,
      out.dir = out.dir.sens
    )

    rds_file <- file.path(
      out.dir.sens,
      paste0(run$model_name, "_survival_excl", exclude_first_n_years, "y_runlist", i, "_bt", n_bootstrap, ".rds")
    )
    saveRDS(res, rds_file)
    intermediate_files <- c(intermediate_files, rds_file)
  }
  
  if (exclude_first_n_years == 0) {
    plot_combined_km_first_run(
      models = c("FDS", "CEA", "Random"),
      current_control = "HC",
      out.dir = out.dir.sens
    )
  }

  sens_results <- dplyr::bind_rows(lapply(intermediate_files, readRDS))
  save_csv_rds(
    sens_results,
    file.path(out.dir.sens, paste0("Survival_results_final_excl", exclude_first_n_years, "y_bt", n_bootstrap))
  )
  
  all_intermediate_files <- c(
    all_intermediate_files,
    paste0(file.path(out.dir.sens, paste0("Survival_results_final_excl", exclude_first_n_years, "y_bt", n_bootstrap)), ".rds")
  )
}

end_time <- Sys.time()
run_time <- end_time - start_time
print(paste0("Total survival bootstrap run time: ", run_time))

# -----------------------------
# 11) Load all saved outputs ####
# -----------------------------
all_results <- dplyr::bind_rows(lapply(all_intermediate_files, readRDS))
save_csv_rds(all_results, file.path(out.dir, paste0("Survival_results_allSensitivity_bt", n_bootstrap)))

all_fixed_metric_files <- unlist(lapply(
  file.path(out.dir, paste0("excludeFirst_", exclude_first_n_years_list, "y")),
  function(x) {
    p <- file.path(x, "fixed_metrics")
    if (!dir.exists(p)) return(character(0))
    list.files(p, pattern = "\\.rds$", full.names = TRUE)
  }
))

all_pred_files <- unlist(lapply(
  file.path(out.dir, paste0("excludeFirst_", exclude_first_n_years_list, "y")),
  function(x) {
    p <- file.path(x, "predictions")
    if (!dir.exists(p)) return(character(0))
    list.files(p, pattern = "\\.rds$", full.names = TRUE)
  }
))

all_fixed_metrics <- if (length(all_fixed_metric_files) > 0) {
  dplyr::bind_rows(lapply(all_fixed_metric_files, readRDS))
} else {
  tibble::tibble()
}

all_predictions_fixedH <- if (length(all_pred_files) > 0) {
  dplyr::bind_rows(lapply(all_pred_files, readRDS))
} else {
  tibble::tibble()
}

save_csv_rds(all_fixed_metrics, file.path(out.dir, "all_fixed_horizon_metrics"))
save_csv_rds(all_predictions_fixedH, file.path(out.dir, "all_fixed_horizon_predictions"))

# -----------------------------
# 12) Summary tables
# -----------------------------
summary_df <- all_results %>%
  dplyr::filter(!is.na(Observation_Years)) %>%
  dplyr::mutate(contrast = paste0(Model, " - ", Group, " vs ", Control_Group)) %>%
  dplyr::group_by(Group, Model, Observation_Years, Control_Group, contrast, ExcludeFirstYears) %>%
  dplyr::summarise(
    # N_case = dplyr::first(N_case),
    tdAUC_median = median(tdAUC, na.rm = TRUE),
    tdAUC_lower  = quantile(tdAUC, 0.025, na.rm = TRUE),
    tdAUC_upper  = quantile(tdAUC, 0.975, na.rm = TRUE),
    C_index_median = median(C_index, na.rm = TRUE),
    C_index_lower  = quantile(C_index, 0.025, na.rm = TRUE),
    C_index_upper  = quantile(C_index, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

save_csv_rds(summary_df, file.path(out.dir, paste0("survival_summary_allSensitivity_bt", n_bootstrap)))
write.csv2(summary_df, file.path(out.dir, paste0("survival_summary_allSensitivity_bt", n_bootstrap, ".csv")), row.names = FALSE)
# read back in with read.csv2 to check formatting:
# summary_df <- read.csv2(file.path(out.dir, paste0("survival_summary_allSensitivity_bt", n_bootstrap, ".csv")))

summary_auc <- all_results %>%
  dplyr::filter(!is.na(Observation_Years)) %>%
  dplyr::group_by(Model, Observation_Years, Control_Group, ExcludeFirstYears) %>%
  dplyr::summarise(
    tdAUC_median = median(tdAUC, na.rm = TRUE),
    tdAUC_lower  = quantile(tdAUC, 0.025, na.rm = TRUE),
    tdAUC_upper  = quantile(tdAUC, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

summary_cindex <- all_results %>%
  dplyr::filter(!is.na(C_index), ExcludeFirstYears == 0) %>%
  dplyr::group_by(Model, Control_Group, Bootstrap) %>%
  dplyr::summarise(C_index = dplyr::first(C_index), .groups = "drop")
summary_cindex$Model <- factor(summary_cindex$Model, levels = c("FDS", "CEA", "Random","age_sex_bmi_PRS","FDS_age_sex_bmi_PRS"))

# -----------------------------
# 13) Main plots ####
# -----------------------------
plot_tdauc_models <- function(summary_auc, models, title, outfile_stub) {
  p <- summary_auc %>%
    dplyr::filter(ExcludeFirstYears == 0, Observation_Years <=14,Model %in% models) %>%
    ggplot(aes(x = Observation_Years, y = tdAUC_median, color = Model, group = Model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = tdAUC_lower, ymax = tdAUC_upper), width = 0.15) +
    facet_wrap(~ Control_Group) +
    # geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = get_model_colors()[names(get_model_colors()) %in% models]) +
     scale_x_continuous(breaks = seq(0, 14, by = 2)) +
    theme_classic() +
    labs(
      title = title,
      subtitle = "Held-out test sets; incident CRC analysis",
      x = "Prediction horizon (years)",
      y = "Time-dependent AUC"
    )
  
  ggsave(file.path(out.dir, paste0(outfile_stub, ".pdf")), p, width = 6.5, height = 3.5, dpi = 300)
  ggsave(file.path(out.dir, paste0(outfile_stub, ".png")), p, width = 6.5, height = 3.5, dpi = 300)
}

plot_cindex_models <- function(summary_cindex, models, title, outfile_stub) {
  plot_df <- summary_cindex %>%
    dplyr::filter(Model %in% models) %>%
    dplyr::group_by(Model, Control_Group) %>%
    dplyr::summarise(
      C_index_median = median(C_index, na.rm = TRUE),
      C_index_lower  = quantile(C_index, 0.025, na.rm = TRUE),
      C_index_upper  = quantile(C_index, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(plot_df, aes(x = Model, y = C_index_median, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = C_index_lower, ymax = C_index_upper), width = 0.2) +
    facet_wrap(~ Control_Group) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = get_model_colors()[names(get_model_colors()) %in% models]) +
    theme_classic() +
    labs(title = title, x = "Model", y = "C-index")
  
  ggsave(file.path(out.dir, paste0(outfile_stub, ".pdf")), p, width = 3, height = 3.4, dpi = 300)
  ggsave(file.path(out.dir, paste0(outfile_stub, ".png")), p, width = 3, height = 3.4, dpi = 300)
}

plot_tdauc_models(
  summary_auc,
  models = c("FDS", "CEA", "Random"),
  title = "Time-dependent AUC: protein-only models",
  outfile_stub = "tdAUC_FDS_vs_CEA_vs_Random"
)

plot_tdauc_models(
  summary_auc,
  models = c("age_sex_bmi_PRS", "FDS_age_sex_bmi_PRS"),
  title = "Time-dependent AUC: incremental value beyond clinical + PRS",
  outfile_stub = "tdAUC_ClinicalPRS_vs_FDSplusClinicalPRS"
)

plot_cindex_models(
  summary_cindex,
  models = c("FDS", "CEA", "Random"),
  title = "C-index: protein-only models",
  outfile_stub = "Cindex_FDS_vs_CEA_vs_Random"
)

plot_cindex_models(
  summary_cindex,
  models = c("age_sex_bmi_PRS", "FDS_age_sex_bmi_PRS"),
  title = "C-index: incremental value beyond clinical + PRS",
  outfile_stub = "Cindex_ClinicalPRS_vs_FDSplusClinicalPRS"
)

# -----------------------------
# 14) Delta AUC / Delta C-index ####
# -----------------------------
main_results <- all_results %>% dplyr::filter(ExcludeFirstYears == 0)

delta_auc_df <- main_results %>%
  dplyr::filter(
    Model %in% c("age_sex_bmi_PRS", "FDS_age_sex_bmi_PRS"),
    !is.na(tdAUC),
    Observation_Years <= 14,
    !is.na(Observation_Years)
  ) %>%
  dplyr::select(Control_Group, Bootstrap, Observation_Years, Model, tdAUC) %>%
  tidyr::pivot_wider(names_from = Model, values_from = tdAUC) %>%
  dplyr::mutate(delta_AUC = FDS_age_sex_bmi_PRS - age_sex_bmi_PRS)

save_csv_rds(delta_auc_df, file.path(out.dir, "deltaAUC_raw_bootstrap_pairs"))

delta_auc_summary <- delta_auc_df %>%
  dplyr::group_by(Control_Group, Observation_Years) %>%
  dplyr::summarise(
    delta_AUC_median = median(delta_AUC, na.rm = TRUE),
    delta_AUC_lower  = quantile(delta_AUC, 0.025, na.rm = TRUE),
    delta_AUC_upper  = quantile(delta_AUC, 0.975, na.rm = TRUE),
    n_boot = sum(!is.na(delta_AUC)),
    .groups = "drop"
  )

write.csv(delta_auc_summary,
          file.path(out.dir, paste0("deltaAUC_FDSplusClinical_vs_Clinical_bt", n_bootstrap, ".csv")),
          row.names = FALSE)

delta_auc_plot <- ggplot(delta_auc_summary, aes(x = Observation_Years, y = delta_AUC_median, group = 1)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = delta_AUC_lower, ymax = delta_AUC_upper), width = 0.15, color = "black") +
   scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  facet_wrap(~ Control_Group) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = expression(Delta*"AUC of FDS beyond clinical + PRS model"),
    subtitle = "(FDS + age + sex + BMI + PRS) minus (age + sex + BMI + PRS)",
    x = "Prediction horizon (years)",
    y = expression(Delta*"AUC")
  ) +
  theme_classic()

ggsave(file.path(out.dir, paste0("deltaAUC_FDSplusClinical_vs_Clinical_bt", n_bootstrap, ".pdf")),
       delta_auc_plot, width = 6.2, height = 3.2, dpi = 300)

delta_cindex_df <- main_results %>%
  dplyr::filter(Model %in% c("age_sex_bmi_PRS", "FDS_age_sex_bmi_PRS")) %>%
  dplyr::group_by(Control_Group, Bootstrap, Model) %>%
  dplyr::summarise(C_index = dplyr::first(C_index), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Model, values_from = C_index) %>%
  dplyr::mutate(delta_Cindex = FDS_age_sex_bmi_PRS - age_sex_bmi_PRS)

save_csv_rds(delta_cindex_df, file.path(out.dir, "deltaCindex_raw_bootstrap_pairs"))

delta_cindex_summary <- delta_cindex_df %>%
  dplyr::group_by(Control_Group) %>%
  dplyr::summarise(
    delta_Cindex_median = median(delta_Cindex, na.rm = TRUE),
    delta_Cindex_lower  = quantile(delta_Cindex, 0.025, na.rm = TRUE),
    delta_Cindex_upper  = quantile(delta_Cindex, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(delta_cindex_summary,
          file.path(out.dir, paste0("deltaCindex_FDSplusClinical_vs_Clinical_bt", n_bootstrap, ".csv")),
          row.names = FALSE)

# -----------------------------
# 15) Fixed-horizon summaries ####
# -----------------------------
fixedH_summary <- all_fixed_metrics %>%
  dplyr::filter(ExcludeFirstYears == 0 ) %>% #& (Model %in% c('FDS','CEA','Random'))
  dplyr::group_by(HorizonYears, Control_Group,Model) %>%
  dplyr::summarise(
    Brier_median = median(Brier_IPCW, na.rm = TRUE),
    Brier_lower = quantile(Brier_IPCW, 0.025, na.rm = TRUE),
    Brier_upper = quantile(Brier_IPCW, 0.975, na.rm = TRUE),
    ScaledBrier_median = median(Scaled_Brier, na.rm = TRUE),
    ScaledBrier_lower = quantile(Scaled_Brier, 0.025, na.rm = TRUE),
    ScaledBrier_upper = quantile(Scaled_Brier, 0.975, na.rm = TRUE),
    Sens_median = median(Sensitivity_at_fixed_spec, na.rm = TRUE),
    Sens_lower = quantile(Sensitivity_at_fixed_spec, 0.025, na.rm = TRUE),
    Sens_upper = quantile(Sensitivity_at_fixed_spec, 0.975, na.rm = TRUE),
    Spec_target = dplyr::first(Spec_target),
    .groups = "drop"
  )

save_csv_rds(fixedH_summary, file.path(out.dir, "FixedH_Brier_ScaledBrier_Sensitivity_FDSplusClinical"))
write.csv2( fixedH_summary, file.path(out.dir, "FixedH_Brier_ScaledBrier_Sensitivity_FDSplusClinical.csv"), row.names = FALSE)
fixedH_summary <- read.csv2(file.path(out.dir, "FixedH_Brier_ScaledBrier_Sensitivity_FDSplusClinical.csv"))

# keep model order consistent
fixedH_summary$Model <- factor(fixedH_summary$Model, levels = c("FDS", "CEA", "Random"))

#  Scaled Brier plot
p_brier = ggplot(fixedH_summary %>% dplyr::filter(Model %in% c("FDS")),
       aes(x = HorizonYears, 
           y = ScaledBrier_median, 
           color = Control_Group, 
           group = Control_Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ScaledBrier_lower, 
                    ymax = ScaledBrier_upper),
                width = 0.2, alpha = 0.6) +
  labs(
    title = "Scaled Brier Score across Prediction Horizons",
    x = "Horizon (Years)",
    y = "Scaled Brier Score",
    color = "Control Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )
ggsave(plot = p_brier, file.path(out.dir, "ScaledBrier_FDSplusClinical_fixedHorizons.pdf"), width = 6.2, height = 3.5, dpi = 300)
ggsave(plot = p_brier, file.path(out.dir, "ScaledBrier_FDSplusClinical_fixedHorizons.png"), width = 6.2, height = 3.5, dpi = 300)


# Sensitivity at fixed specificity plot
spec_lab <- unique(na.omit(fixedH_summary$Spec_target))
spec_lab <- if (length(spec_lab) == 1) paste0(round(spec_lab * 100), "%") else "fixed specificity"

p_sens_at_spec <- ggplot(
  fixedH_summary %>% dplyr::filter(Model %in% c("FDS","CEA","Random")),
  aes(
    x = HorizonYears,
    y = Sens_median,
    color = Model,
    group = Model
  )
) +
  geom_line(linewidth = 1) +
   scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = Sens_lower, ymax = Sens_upper),
    width = 0.2,
    alpha = 0.6
  ) +
  scale_color_manual(values = get_model_colors()[c("FDS", "CEA", "Random")]) +
  facet_wrap(~ Control_Group) +
  labs(
    title = paste0("Sensitivity at ", spec_lab, " Specificity across Prediction Horizons"),
    x = "Horizon (Years)",
    y = paste0("Sensitivity at ", spec_lab, " specificity"),
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(out.dir, "SensitivityAtFixedSpec_byModel_fixedHorizons.pdf"),
  plot = p_sens_at_spec, width = 6.2, height = 3.5, dpi = 300
)
ggsave(
  file.path(out.dir, "SensitivityAtFixedSpec_byModel_fixedHorizons.png"),
  plot = p_sens_at_spec, width = 6.2, height = 3.5, dpi = 300
)

# Calibration plots
for (hh in fixed_eval_horizons) {
  cal_file <- file.path(
    out.dir, "excludeFirst_0y",
    paste0("CalibrationData_FDS_age_sex_bmi_PRS_HC_", hh, "y_bootstrap1.rds")
  )
  if (file.exists(cal_file)) {
    cal_df <- readRDS(cal_file)
    plot_calibration_curve(
      cal_df = cal_df,
      model_name = "FDS_age_sex_bmi_PRS",
      current_control = "HC",
      horizon = hh,
      out.dir = out.dir
    )
  }
}

# -----------------------------
# 16) Reverse-causation sensitivity plot ####
# -----------------------------
auc_sens_plot <- summary_auc %>%
  dplyr::filter(Model == "FDS",
                summary_auc$Observation_Years <= 14) %>%
  ggplot(aes(x = Observation_Years, y = tdAUC_median,
             color = factor(ExcludeFirstYears),
             group = factor(ExcludeFirstYears))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = tdAUC_lower, ymax = tdAUC_upper), width = 0.15) +
   scale_x_continuous(breaks = seq(0, 14, by = 2)) +
   scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
  facet_wrap(~ Control_Group) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("0" = "#F4A261", "2" = "#5DA5DA", "4" = "#C44A3A")) +
  theme_classic() +
  labs(
    title = "Reverse-causation sensitivity analysis",
    subtitle = "FDS after excluding early incident cases",
    x = "Prediction horizon (years)",
    y = "Time-dependent AUC",
    color = "Excluded first years"
  )

ggsave(file.path(out.dir, "ReverseCausation_tdAUC_FDS.pdf"),
       auc_sens_plot, width = 6.6, height = 3.6, dpi = 300)
ggsave(file.path(out.dir, "ReverseCausation_tdAUC_FDS.png"),
       auc_sens_plot, width = 6.6, height = 3.6, dpi = 300)

print(paste0("DONE. Output folder: ", out.dir))
print(paste0("Total survival bootstrap run time: ", run_time))