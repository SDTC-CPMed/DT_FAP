# ==============================================================================
# PROJECT: UK Biobank Proteomics Risk Prediction Analysis (Machine Learning)
# DESCRIPTION:
# This script performs protein-based risk prediction for specified diseases
# (e.g., Colorectal Cancer) using the UK Biobank proteomics data.
# Key steps include: data loading, Healthy Control (HC) subsetting, matching
# cases to HC by age/sex, KNN imputation, Ridge Regression (alpha=0) modeling
# within a bootstrap framework, and final results summary and visualization in 
# the test data.

# ==============================================================================
# --- 0. INITIAL SETUP AND PACKAGE LOADING --- #
rm(list = ls())
packages <- c("dplyr", "stringr", "MatchIt", "tibble", "survival", "Matrix",
              "ggplot2", "patchwork", "reshape2", "readxl", "data.table",
              "glmnet", "reticulate", "VIM", "purrr", "broom", "impute",
              "pROC", "ROCR", "tidyverse", "foreach", "doRNG", "doParallel")
lapply(packages, library, character.only = TRUE)

# Remove all variables except UKBB_Proteomics (if it was passed from another session)

# --- 1. CONFIGURATION AND DIRECTORY SETUP --- #
input.dir <- 'UKBB/input'
input.dir2 <- input.dir 
out.dir <- 'output'
dir.create(out.dir, showWarnings = FALSE)

# Set analysis parameters
dis1 <- 'Colon_premalignant' 
dis2 <- 'Colorectal_cancer' 
cases <- c('dis2_prev', 'dis2_inci') # Disease comparison groups
control <- "HC"
goi <- 'signature' # Gene of interest set to use in the analysis
n_bootstrap <- 1000
set.seed(123)

# --- 2. DATA LOADING AND PREPARATION (UKBB) --- #
if (!exists("UKBB_Proteomics")) {
  UKBB_Proteomics <- read.table(
    file = file.path(input.dir, 'Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt'),
    sep = '\t', header = TRUE, fill = TRUE, row.names = 'PID'
  )
}

# Load and process core clinical data
All_clinic <- readRDS(file.path(input.dir2, 'ukb_clinic.rda')) %>% 
  filter(eid %in% rownames(UKBB_Proteomics)) %>% 
  # Convert date fields to Date format and factors
  mutate(across(c(date_of_attending_assessment_centre_f53_0_0, date_of_death_f40000_0_0, date_of_death_f40000_1_0, date_lost_to_followup_f191_0_0), as.Date)) %>%
  mutate(sex_f31_0_0 = as.factor(sex_f31_0_0),
         ethnic_background_f21000_0_0 = as.factor(ethnic_background_f21000_0_0),
         delta_diag_enroll = NA) %>%
  # Calculate OS status and time
  mutate(OSstatus.raw = ifelse(is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0), 0, 1),
         OStime = ifelse(OSstatus.raw == 1, 
                         (date_of_death_f40000_0_0 - date_of_attending_assessment_centre_f53_0_0)/365,
                         ifelse(is.na(date_lost_to_followup_f191_0_0),
                                (as.Date("2022-12-30") - date_of_attending_assessment_centre_f53_0_0)/365,
                                (pmin(date_lost_to_followup_f191_0_0, as.Date("2022-12-30")) - date_of_attending_assessment_centre_f53_0_0)/365))) %>% 
  # Filter out individuals with missing key demographics
  filter(!is.na(ethnic_background_f21000_0_0),
         !is.na(age_at_recruitment_f21022_0_0),
         !is.na(uk_biobank_assessment_centre_f54_0_0),
         !is.na(sex_f31_0_0)) %>% 
  mutate(eid = as.character(eid))

# Load disease-specific clinical data and filter by proteomics participants
load_disease_data <- function(dis_name) {
  readRDS(file.path(input.dir2, paste0('ukb_', dis_name, '_subset_clinic.rda'))) %>% 
    filter(eid %in% rownames(UKBB_Proteomics)) %>% 
    mutate(eid = as.character(eid))
}

dis1_clinic <- load_disease_data(dis1)
dis2_clinic <- load_disease_data(dis2)
All_cancer_clinic_eid <- readRDS(file.path(input.dir2,'ukb_AllCancers_subset_clinic.rda')) %>% filter(eid %in% rownames(UKBB_Proteomics)) %>% pull(eid)

# Function to extract EIDs for prevalent and incident cases
extract_case_eids <- function(clinic_df, prevalent = TRUE) {
  if (prevalent) {
    clinic_df %>% filter(diagnosis_time <= date_of_attending_assessment_centre_f53_0_0) %>% pull(eid)
  } else {
    clinic_df %>% filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0) %>% pull(eid)
  }
}

# Extract specific EID lists
dis1_eid <- dis1_clinic %>% dplyr::select(eid, diagnosis_time) %>% pull(eid)
dis1_prev_eid <- extract_case_eids(dis1_clinic, prevalent = TRUE)
dis1_inci_eid <- extract_case_eids(dis1_clinic, prevalent = FALSE)
dis1_prev_clinic <- dis1_clinic %>% filter(diagnosis_time <= date_of_attending_assessment_centre_f53_0_0)
dis1_inci_clinic <- dis1_clinic %>% filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0)

dis2_eid <- dis2_clinic %>% dplyr::select(eid, diagnosis_time) %>% pull(eid)
dis2_prev_eid <- extract_case_eids(dis2_clinic, prevalent = TRUE)
dis2_inci_eid <- extract_case_eids(dis2_clinic, prevalent = FALSE)
dis2_prev_clinic <- dis2_clinic %>% filter(diagnosis_time <= date_of_attending_assessment_centre_f53_0_0)
dis2_inci_clinic <- dis2_clinic %>% filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0)

# Load and prepare Healthy Controls (HC)
HC_clinic <- readRDS(file.path(input.dir2,'ukb_healthy_controls.rda')) %>% filter(eid %in% rownames(UKBB_Proteomics)) 
HC_clinic <- All_clinic %>% filter(eid %in% as.character(HC_clinic$eid)) 
HC_eid <- HC_clinic %>% filter(eid %in% rownames(UKBB_Proteomics)) %>% pull(eid)

# --- Subgroup definition for comparing Dis1 vs Dis2 cases (re-filter for pre-malignant) ---
# Goal: Define Dis1 cases that eventually progressed to Dis2 (CRC), and Dis1 cases that did not (Dis1 no Dis2).
All_cancer_clinic_2 <- All_clinic %>% 
  left_join(dis1_clinic %>% select(eid, diagnosis_time.dis1 = diagnosis_time), by = 'eid') %>% 
  left_join(dis2_clinic %>% select(eid, diagnosis_time.dis2 = diagnosis_time), by = 'eid') %>% 
  mutate(across(starts_with("diagnosis_time"), as.Date)) %>%
  # Check if Dis1 diagnosis existed and was before or concurrent with Dis2 diagnosis (or Dis2 is NA)
  mutate(dis1_before_dis2 = ifelse((!is.na(diagnosis_time.dis1) & (is.na(diagnosis_time.dis2) | diagnosis_time.dis1 <= diagnosis_time.dis2)), 1, 0)) %>% 
  mutate(dis1 = ifelse(eid %in% dis1_eid, 1, 0)) %>% 
  mutate(CRC = ifelse(eid %in% dis2_eid, 1, 0))

# Dis1 cases that progressed to Dis2 (CRC)
Dis1andDis2_clinic <- All_cancer_clinic_2 %>% 
  filter(dis1_before_dis2 == 1, CRC == 1) %>% 
  rename(diagnosis_time = diagnosis_time.dis2) # Use CRC diagnosis time for follow-up
Dis1andDis2_eid <- Dis1andDis2_clinic %>% pull(eid)

# Dis1 cases that did NOT progress to Dis2 (No CRC)
Dis1noDis2_clinic <- All_cancer_clinic_2 %>% 
  filter(dis1_before_dis2 == 1, CRC == 0) %>%
  # Need to ensure no CRC diagnosis date exists for this group.
  rename(diagnosis_time = diagnosis_time.dis1) # Use Dis1 diagnosis time for clarity
Dis1noDis2_eid <- Dis1noDis2_clinic %>% pull(eid)


# --- 3. PROTEIN SET DEFINITION --- #
known_marker <- c('CEACAM5') # CEA
signature <- c("LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F", "INHBA", 
               "LAMA5", "COL7A1", "SOX9", "CCND1", "MET", "ID1", "PTMA", "CLU", "TGIF1", 
               "EDN1", "GDF15") 

# Filter proteins for missing values (less than 25% missingness)
UKBB_Proteomics <- UKBB_Proteomics[, colSums(is.na(UKBB_Proteomics)) <= 0.25 * nrow(UKBB_Proteomics)]
exist_signature_proteins <- intersect(get(goi), colnames(UKBB_Proteomics))
all_other_proteins <- setdiff(colnames(UKBB_Proteomics), signature)
UKBB_Proteomics$eid <- rownames(UKBB_Proteomics)

# --- 4. MACHINE LEARNING ANALYSIS FUNCTION (RIDGE REGRESSION) --- #
# This function performs sample matching, bootstrapping, KNN imputation, 
# Ridge regression, and AUC calculation.
analyze_proteins_multiyear_multicontrol <- function(
    protein_set, set_name,
    cases = c('dis2_inci'),
    controls = c('HC'), 
    n_obs_yrs_vec = c(2,4, 6, 8, 10,12,14),
    ratio = 5,
    n_bootstrap = 100,
    alpha = 0, # Alpha=0 forces Ridge Regression
    ncores = 20) {
  
  # --- Define proteins to test ---
  if (grepl("signature", set_name, ignore.case = TRUE) || set_name == 'CEA') {
    all_protein_cols <- intersect(get(set_name), colnames(UKBB_Proteomics))
  } else if (set_name == 'random') {
    all_protein_cols <- all_other_proteins
  } else {
    all_protein_cols <- intersect(protein_set, colnames(UKBB_Proteomics))
  }
  
  map_dfr(cases, function(case) {
    map_dfr(controls, function(current_control) { 
      map_dfr(n_obs_yrs_vec, function(n_obs_yrs) {
        message(paste("Working: Setname:", set_name,
                      "| Case:", case,
                      "| Control:", current_control,
                      "| Years:", n_obs_yrs))
        
        # --- 4.1. Prepare Case Data (Filtering by Observation Window) ---
        
        # Determine the source data for the case group EIDs and clinic
        case_eid_var <- switch(case,
                               'dis1_prev' = "dis1_prev_eid",
                               'dis1_inci' = "dis1_inci_eid",
                               'dis2_prev' = "dis2_prev_eid",
                               'dis2_inci' = "dis2_inci_eid",
                               'Dis1andDis2' = "Dis1andDis2_eid",
                               stop(paste("Unknown case group:", case)))
        
        case_clinic_var <- switch(case,
                                  'dis1_prev' = "dis1_prev_clinic",
                                  'dis1_inci' = "dis1_inci_clinic",
                                  'dis2_prev' = "dis2_prev_clinic",
                                  'dis2_inci' = "dis2_inci_clinic",
                                  'Dis1andDis2' = "Dis1andDis2_clinic",
                                  stop(paste("Unknown case clinic:", case)))
        
        case_data <- All_clinic %>%
          filter(eid %in% get(case_eid_var)) %>%
          left_join(UKBB_Proteomics, by = "eid") %>%
          left_join(get(case_clinic_var) %>% select(eid, diagnosis_time), by = "eid") %>%
          # Filter cases where diagnosis is within the observation window
          filter(abs(diagnosis_time - date_of_attending_assessment_centre_f53_0_0) <= n_obs_yrs*365) %>%
          mutate(binary_group = 1, # Case = 1
                 age = as.numeric(age_at_recruitment_f21022_0_0),
                 sex = ifelse(sex_f31_0_0 == "Female", 0, ifelse(sex_f31_0_0 == "Male", 1, NA))) %>%
          select(eid, binary_group, age, sex, all_of(colnames(UKBB_Proteomics)))
        
        if(nrow(case_data) < 20) return(tibble(
          Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
          N_case = nrow(case_data), N_used_Protein = NA, Selected_Proteins = NA,
          AUC = NA, Protein_Set = set_name, Bootstrap = NA,
          Observation_Years = n_obs_yrs, Control_Group = current_control
        ))
        
        # --- 4.2. Prepare Control Data (Filtering by Follow-up Time) ---
        
        # Note: 'Dis1noDis2' is a specific comparison group, not a standard control.
        if (current_control == 'HC') {
          control_data_clinic <- HC_clinic
        } else if (current_control == 'Dis1noDis2') {
          control_data_clinic <- Dis1noDis2_clinic
        } else {
          # Should not happen based on configuration, but ensures robustness
          stop("Invalid control group specified.")
        }
        
        control_data <- All_clinic %>%
          filter(eid %in% control_data_clinic$eid) %>%
          # Filter out controls that died or were censored before the observation window ends
          filter(OSstatus.raw != 1 | 
                   (OStime * 365) > n_obs_yrs*365) %>%
          left_join(UKBB_Proteomics, by = "eid") %>%
          mutate(binary_group = 0, # Control = 0
                 age = as.numeric(age_at_recruitment_f21022_0_0),
                 sex = ifelse(sex_f31_0_0 == "Female", 0, ifelse(sex_f31_0_0 == "Male", 1, NA))) %>%
          select(eid, binary_group, age, sex, all_of(colnames(UKBB_Proteomics)))
        
        # --- 4.3. Age and Sex Matching (1:Ratio with caliper) ---
        
        # Match cases to controls by age (caliper=2 years) and exactly by sex
        match_obj <- matchit(binary_group ~ age,
                             data = bind_rows(case_data, control_data),
                             method = "nearest", ratio = ratio, replace = TRUE,
                             exact = ~ sex, caliper = c(age = 2))
        
        full_data <- match.data(match_obj) %>%
          select(eid, binary_group, age, sex, all_of(colnames(UKBB_Proteomics))) %>%
          mutate(binary_group = factor(binary_group))
        
        # --- 4.4. Bootstrap Modeling (Parallel Processing) ---
        
        registerDoParallel(cores = ncores)
        
        results_boot <- foreach(i = 1:n_bootstrap, .combine = rbind,
                                .packages = c('dplyr','glmnet','broom','pROC','ROCR','impute','MatchIt','tibble')) %dopar% {
                                  
                                  set.seed(1000 + i) # Bootstrap-specific seed for reproducibility
                                  
                                  # --- Random protein selection per bootstrap (for 'random' set only) ---
                                  if (set_name == 'random') {
                                    length_signature <- length(exist_signature_proteins)
                                    # Sample a random set of proteins of the same size as the main signature
                                    protein_set_i <- sample(all_protein_cols, length_signature)
                                  } else {
                                    # Use the defined protein set for signature/CEA runs
                                    protein_set_i <- intersect(all_protein_cols, colnames(UKBB_Proteomics))
                                  }
                                  
                                  if(length(protein_set_i) == 0) {
                                    # Handle case where no proteins are available for the set
                                    return(tibble(
                                      Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                                      N_case = nrow(case_data), N_used_Protein = 0, Selected_Proteins = NA,
                                      AUC = NA, Protein_Set = set_name, Bootstrap = i,
                                      Observation_Years = n_obs_yrs, Control_Group = current_control
                                    ))
                                  }
                                  
                                  # --- Train/Test Split (70/30) ---
                                  set.seed(2000 + i)
                                  case_idx <- which(full_data$binary_group == 1)
                                  control_idx <- which(full_data$binary_group == 0)
                                  
                                  # Sample 70% of cases and controls for the training set
                                  train_idx <- c(sample(case_idx, round(0.7*length(case_idx))),
                                                 sample(control_idx, round(0.7*length(control_idx))))
                                  test_idx <- setdiff(seq_len(nrow(full_data)), train_idx)
                                  
                                  train_data <- full_data[train_idx, ]
                                  test_data <- full_data[test_idx, ]
                                  
                                  proteins_to_impute <- setdiff(colnames(UKBB_Proteomics), "eid")
                                  proteins_to_keep <- intersect(proteins_to_impute, protein_set_i)
                                  
                                  # --- Imputation using Training Data (impute.knn) ---
                                  protein_matrix_train <- as.matrix(train_data[, proteins_to_impute])
                                  protein_matrix_train <- apply(protein_matrix_train, 2, as.numeric)
                                  # Impute missing values in the training set using k-NN
                                  imputed_train <- impute.knn(protein_matrix_train)$data
                                  train_data[, proteins_to_impute] <- imputed_train
                                  train_data <- train_data[, c('binary_group', proteins_to_keep)] 
                                  
                                  # --- Imputation/Handling of Test Data (Mean Substitution based on Train) ---
                                  for(col in proteins_to_keep) {
                                    test_mean <- mean(train_data[[col]], na.rm = TRUE)
                                    # Substitute NA in test data with the mean from the training data
                                    test_data[is.na(test_data[[col]]), col] <- test_mean
                                  }
                                  test_data <- test_data[, c('binary_group', proteins_to_keep)] 
                                  
                                  # --- Modeling (Ridge Regression via glmnet) ---
                                  if(set_name %in% c("signature", "random")) {
                                    x_train <- as.matrix(train_data[, proteins_to_keep])
                                    y_train <- as.numeric(train_data$binary_group) - 1 # Needs to be 0 or 1
                                    
                                    # Cross-validated glmnet (Ridge)
                                    cv_fit <- cv.glmnet(x_train, y_train, family = "binomial",
                                                        alpha = alpha, type.measure = "auc", standardize = TRUE)
                                    
                                    # Extract coefficients at lambda.min
                                    coefs <- coef(cv_fit, s = "lambda.min")
                                    sel_weights <- as.numeric(coefs[proteins_to_keep, ])
                                    
                                    # Calculate protein risk score on the TEST set
                                    test_data$protein_score <- as.numeric(as.matrix(test_data[, proteins_to_keep, drop=FALSE]) %*% sel_weights)
                                  } else { 
                                    # For single-protein markers (like CEA), score is just the protein level itself
                                    test_data <- test_data %>% mutate(protein_score = get(proteins_to_keep))
                                  }
                                  
                                  # --- Final Logistic Regression and Evaluation on TEST set ---
                                  final_model <- glm(binary_group ~ protein_score, data = test_data, family = binomial)
                                  tidy_model <- broom::tidy(final_model)
                                  score_row <- filter(tidy_model, term == "protein_score")
                                  
                                  # Calculate AUC
                                  pred <- ROCR::prediction(predict(final_model, type="response"), test_data$binary_group)
                                  auc_val <- tryCatch(as.numeric(performance(pred, "auc")@y.values), error = function(e) NA_real_)
                                  
                                  # Compile bootstrap result
                                  tibble(
                                    Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                                    N_case = nrow(case_data),
                                    N_used_Protein = length(proteins_to_keep),
                                    Selected_Proteins = paste(proteins_to_keep, collapse = ";"),
                                    AUC = auc_val,
                                    Protein_Set = set_name,
                                    Bootstrap = i,
                                    Observation_Years = n_obs_yrs,
                                    Control_Group = current_control
                                  )
                                }
        # Stop parallel backend for cleanup
        stopImplicitCluster() 
        return(results_boot)
      })
    })
  })
}

# --- 5. EXECUTE ANALYSES AND MERGE RESULTS --- #
n_bootstrap <- 1000
observation_years <- c(2, 4, 6, 8, 10, 12, 14) 
alpha <- 0
set.seed(123)

runs_list <- list(
  # CRC Incident Cases vs Healthy Controls (HC)
  list(protein_set = signature, set_name = "signature", cases = c('dis2_inci'), controls = c("HC"), alpha = alpha),
  list(protein_set = known_marker, set_name = "CEA", cases = c('dis2_inci'), controls = c("HC"), alpha = alpha),
  list(protein_set = "random", set_name = "random", cases = c('dis2_inci'), controls = c("HC"), alpha = alpha),
  
  # CRC Prevalent Cases vs Healthy Controls (HC)
  list(protein_set = signature, set_name = "signature", cases = c('dis2_prev'), controls = c("HC"), alpha = alpha),
  list(protein_set = known_marker, set_name = "CEA", cases = c('dis2_prev'), controls = c("HC"), alpha = alpha),
  list(protein_set = "random", set_name = "random", cases = c('dis2_prev'), controls = c("HC"), alpha = alpha),
  
  # Dis1 progressed to Dis2 (CRC) vs Dis1 non-progressed (Case vs Case)
  # NOTE: 'Dis1noDis2' is a specific non-progressing case cohort, not a general control.
  list(protein_set = signature, set_name = "signature", cases = 'Dis1andDis2', controls = "Dis1noDis2", alpha = alpha),
  list(protein_set = known_marker, set_name = "CEA", cases = 'Dis1andDis2', controls = "Dis1noDis2", alpha = alpha),
  list(protein_set = "random", set_name = "random", cases = 'Dis1andDis2', controls = "Dis1noDis2", alpha = alpha)
)

# Execute the runs
intermediate_files <- c()
start_time <- Sys.time()

for(i in seq_along(runs_list)){
  run <- runs_list[[i]]
  
  # Execute the analysis function
  res <- analyze_proteins_multiyear_multicontrol(
    protein_set = run$protein_set,
    set_name = run$set_name,
    cases = run$cases,
    controls = run$controls,
    n_obs_yrs_vec = observation_years,
    n_bootstrap = n_bootstrap
  )
  
  # Save intermediate RDS file for safety
  rds_file <- file.path(out.dir, paste0(run$set_name, "_runlist", i, "_bt", n_bootstrap, ".rds"))
  saveRDS(res, rds_file)
  intermediate_files <- c(intermediate_files, rds_file)
}
end_time <- Sys.time()
run_time <- end_time - start_time
print(paste0('Total bootstrap run time: ', round(run_time, 2), " ", units(run_time)))

# Merge and save final results
all_results <- bind_rows(lapply(intermediate_files, readRDS))
saveRDS(all_results, file.path(out.dir, paste0("ML_results_final_bt", n_bootstrap, ".rds")))
write.csv(all_results, file.path(out.dir, paste0("ML_results_final_bt", n_bootstrap, ".csv")), row.names = FALSE)


# --- 6. SUMMARIZE RESULTS (Bootstrap Aggregation) --- #
all_results <- readRDS(file.path(out.dir, paste0("ML_results_final_bt", n_bootstrap, ".rds")))
summary_df <- all_results %>% 
  mutate(contrast = paste0(Protein_Set,' - ',Group, " vs ", Control_Group)) %>%
  group_by(Group, Protein_Set, Observation_Years, Control_Group, contrast) %>% 
  dplyr::summarise(
    N_case = unique(N_case, na.rm = TRUE),
    # Calculate median and 95% confidence intervals (2.5% and 97.5% quantiles)
    AUC_median = median(AUC, na.rm = TRUE),
    AUC_lower = quantile(AUC, 0.025, na.rm = TRUE),
    AUC_upper = quantile(AUC, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  # Factorize Control_Group, removing 'non_cancer' and 'OC' levels
  mutate(Control_Group = factor(Control_Group, 
                                levels = c('HC', 'Dis1noDis2'), # Only keeping HC and the subgroup control
                                labels = c('Healthy Controls', 'Dis1 No Progression'))) 

# Save summary results
write.csv(summary_df, file.path(out.dir, paste0("glm_result_summary_bt",n_bootstrap,"_HC_only.csv")), row.names = FALSE)
saveRDS(summary_df, file.path(out.dir, paste0("glm_result_summary_bt",n_bootstrap,"_HC_only.rds")))


# --- 7. VISUALIZATION (AUC Plots) --- # 
auc_plot <- ggplot(summary_df[summary_df$Control_Group == "Healthy Controls", ], 
                   aes(x = Observation_Years, y = AUC_median, 
                       color = Protein_Set, group = Protein_Set)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper), width = 0.2) +
  facet_grid(. ~ Group, scales = "fixed") +
  labs(title = "Model AUC vs Healthy Controls",
       subtitle = "Median AUC (95% CI) by Observation Period",
       x = "Observation Years Before/At Diagnosis", 
       y = "AUC",
       color = "Protein Set") +
  scale_color_manual(values = c("signature" = "#FC7F0E", "random" = '#7F7F7F',"CEA" ="#2CA02C","signature_EMBO2015" = "pink", "signature_JTM2024" = "lightblue")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = 'grey') + # Reference line at AUC=0.5
  scale_y_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, by = 0.1)) +
  scale_x_continuous(limits = c(1.9, 14.1), breaks = observation_years) +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out.dir, paste0("AUC_plot_bt",n_bootstrap,"_HC_only.pdf")), auc_plot,
       width = 7, height = 3.5, dpi = 300)

print("Visualization complete. Check the output directory for PDF plots and CSV/RDS summaries.")