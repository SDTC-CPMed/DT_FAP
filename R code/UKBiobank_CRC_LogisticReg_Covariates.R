# ==============================================================================
# PROJECT: UK Biobank Proteomics Risk Prediction Analysis (Sum Score Approach)
# DESCRIPTION:
# This script performs protein-based risk prediction using the mean protein level
# (Sum Score) derived from a defined Gene of Interest (GOI) list.
# The analysis involves bootstrapped control resampling,
# logistic regression adjusted for age, sex, BMI, PRS, smoking, alcohol,
# and diabetes status, and calculation of OR.
#
# METHOD: logistic regression.
# ==============================================================================

rm(list = ls())

# --- 0. INITIAL SETUP AND PACKAGE LOADING --- #
packages <- c(
  "dplyr",  "ggplot2", "patchwork",  "stringr",  "data.table",  "glmnet",
  "broom",  "MatchIt", "purrr",  "pROC",  "foreach",  "doParallel", "tidyr")

invisible(lapply(packages, library, character.only = TRUE))

# --- 1. CONFIGURATION AND DIRECTORY SETUP --- #
input.dir <- "~/Documents/fap_ukbb/inputs"
input.dir2 <- "~/Documents/fap_ukbb/subsets"
out.dir <- "output"
dir.create(out.dir, showWarnings = FALSE)

dis1 <- "Colon_premalignant"
dis2 <- "Colorectal"
cases <- "dis2_inci"
control <- "HC"
goi <- "FAP_list"
n_bootstrap <- 1000
observation_years <- c(2, 4, 6, 8, 10, 12, 14)

set.seed(123)
registerDoParallel(cores = 10)

# --- 2. DATA LOADING AND CLINICAL SUBSETTING --- #
UKBB_Proteomics <- read.table(
  file = file.path(input.dir, "Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt"),
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  row.names = "PID",
  check.names = FALSE
)

# Load and process core clinical data
All_clinic <- readRDS(file.path(input.dir, "ukb_clinic.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        date_of_attending_assessment_centre_f53_0_0,
        date_of_death_f40000_0_0,
        date_of_death_f40000_1_0,
        date_lost_to_followup_f191_0_0
      ),
      as.Date
    )
  ) %>%
  dplyr::mutate(
    sex_f31_0_0 = as.factor(sex_f31_0_0),
    ethnic_background_f21000_0_0 = as.factor(ethnic_background_f21000_0_0),
    
    smoking_status_f20116_0_0 = stats::relevel(
      factor(dplyr::na_if(trimws(as.character(smoking_status_f20116_0_0)), "Prefer not to answer")),
      ref = "Never"
    ),
    
    alcohol_intake_frequency_f1558_0_0 = stats::relevel(
      factor(dplyr::na_if(trimws(as.character(alcohol_intake_frequency_f1558_0_0)), "Prefer not to answer")),
      ref = "Never"
    ),
    
    diabetes_diagnosed_by_doctor_f2443_0_0 = stats::relevel(
      factor(
        dplyr::na_if(
          dplyr::na_if(trimws(as.character(diabetes_diagnosed_by_doctor_f2443_0_0)), "Prefer not to answer"),
          "Do not know"
        )
      ),
      ref = "No"
    ),
    
    delta_diag_enroll = NA_real_,
    OSstatus.raw = ifelse(
      is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0),
      0,
      1
    ),
    OStime = ifelse(
      OSstatus.raw == 1,
      (date_of_death_f40000_0_0 - date_of_attending_assessment_centre_f53_0_0) / 365,
      ifelse(
        is.na(date_lost_to_followup_f191_0_0),
        (as.Date("2022-12-30") - date_of_attending_assessment_centre_f53_0_0) / 365,
        (pmin(date_lost_to_followup_f191_0_0, as.Date("2022-12-30")) -
           date_of_attending_assessment_centre_f53_0_0) / 365
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(ethnic_background_f21000_0_0),
    !is.na(age_at_recruitment_f21022_0_0),
    !is.na(uk_biobank_assessment_centre_f54_0_0),
    !is.na(sex_f31_0_0)
  ) %>%
  dplyr::mutate(eid = as.character(eid))

# --- 3A. MERGE PRS --- #
prs_data <- read.csv(
  file.path(input.dir, "Matrix_scaled_CombinedPRS_EUROS_Martin_pancamcer_PRSs.csv"),
  sep = ";"
)

prs_data <- prs_data %>%
  dplyr::rename(
    eid = eid,
    PRS = CRC_new_SCORE
  ) %>%
  dplyr::mutate(
    eid = as.character(eid),
    PRS = as.numeric(PRS)
  ) %>%
  dplyr::select(eid, PRS)

All_clinic <- All_clinic %>%
  dplyr::left_join(prs_data, by = "eid")

# Add simplified covariate aliases used later
All_clinic <- All_clinic %>%
  dplyr::mutate(
    age = age_at_recruitment_f21022_0_0,
    sex = sex_f31_0_0
  )

# --- 3B. LOAD DISEASE-SPECIFIC CLINICAL DATA --- #
load_disease_data <- function(dis_name) {
  readRDS(file.path(input.dir2, paste0("ukb_", dis_name, "_subset_clinic.rda"))) %>%
    dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
    dplyr::mutate(
      eid = as.character(eid),
      diagnosis_time = as.Date(diagnosis_time)
    )
}

print(input.dir2)

dis1_clinic <- load_disease_data(dis1)
dis2_clinic <- load_disease_data(dis2)

All_cancer_clinic_eid <- readRDS(file.path(input.dir2, "ukb_AllCancers_subset_clinic.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
  dplyr::pull(eid)

dis1_eid <- dis1_clinic %>% dplyr::select(eid, diagnosis_time) %>% dplyr::pull(eid)
dis2_eid <- dis2_clinic %>% dplyr::select(eid, diagnosis_time) %>% dplyr::pull(eid)

# Function to extract EIDs and clinic data for prevalent and incident cases
extract_cases <- function(clinic_df, prevalent = TRUE) {
  tmp <- clinic_df %>%
    dplyr::select(-dplyr::any_of("date_of_attending_assessment_centre_f53_0_0")) %>%
    dplyr::left_join(
      All_clinic %>% dplyr::select(eid, date_of_attending_assessment_centre_f53_0_0),
      by = "eid"
    )
  
  if (prevalent) {
    tmp %>% dplyr::filter(diagnosis_time <= date_of_attending_assessment_centre_f53_0_0)
  } else {
    tmp %>% dplyr::filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0)
  }
}

dis1_prev_clinic <- extract_cases(dis1_clinic, prevalent = TRUE)
dis1_inci_clinic <- extract_cases(dis1_clinic, prevalent = FALSE)
dis2_prev_clinic <- extract_cases(dis2_clinic, prevalent = TRUE)
dis2_inci_clinic <- extract_cases(dis2_clinic, prevalent = FALSE)

dis1_prev_eid <- dis1_prev_clinic %>% dplyr::pull(eid)
dis1_inci_eid <- dis1_inci_clinic %>% dplyr::pull(eid)
dis2_prev_eid <- dis2_prev_clinic %>% dplyr::pull(eid)
dis2_inci_eid <- dis2_inci_clinic %>% dplyr::pull(eid)

# --- 3C. LOAD AND PREPARE HEALTHY CONTROLS --- #
HC_clinic <- readRDS(file.path(input.dir2, "ukb_healthy_controls.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics))

HC_clinic <- All_clinic %>%
  dplyr::filter(eid %in% as.character(HC_clinic$eid))

HC_eid <- HC_clinic %>% dplyr::pull(eid)

# --- 3D. SUBGROUP DEFINITION FOR PRE-MALIGNANT PROGRESSION --- #
All_cancer_clinic_2 <- All_clinic %>%
  dplyr::left_join(
    dis1_clinic %>% dplyr::select(eid, diagnosis_time.dis1 = diagnosis_time),
    by = "eid"
  ) %>%
  dplyr::left_join(
    dis2_clinic %>% dplyr::select(eid, diagnosis_time.dis2 = diagnosis_time),
    by = "eid"
  ) %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("diagnosis_time"), as.Date))

table(is.na(All_cancer_clinic_2$diagnosis_time.dis2))

Dis1andDis2_clinic <- All_cancer_clinic_2 %>%
  dplyr::mutate(CRC = ifelse(eid %in% dis2_eid, 1, 0)) %>%
  dplyr::filter(
    !is.na(diagnosis_time.dis1),
    is.na(diagnosis_time.dis2) | diagnosis_time.dis1 <= diagnosis_time.dis2,
    CRC == 1
  ) %>%
  dplyr::rename(diagnosis_time = diagnosis_time.dis2)

Dis1andDis2_eid <- Dis1andDis2_clinic %>% dplyr::pull(eid)

Dis1noDis2_clinic <- All_cancer_clinic_2 %>%
  dplyr::mutate(CRC = ifelse(eid %in% dis2_eid, 1, 0)) %>%
  dplyr::filter(!is.na(diagnosis_time.dis1), CRC == 0) %>%
  dplyr::rename(diagnosis_time = diagnosis_time.dis1)

Dis1noDis2_eid <- Dis1noDis2_clinic %>% dplyr::pull(eid)

# --- 4. PROTEIN SET DEFINITION --- #
known_marker <- c("CEACAM5")

FAP_list <- c(
  "LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F", "INHBA",
  "LAMA5", "COL7A1", "SOX9", "CCND1", "MET", "ID1", "PTMA", "CLU", "TGIF1",
  "EDN1", "GDF15"
)

UKBB_Proteomics <- UKBB_Proteomics[, colSums(is.na(UKBB_Proteomics)) <= 0.25 * nrow(UKBB_Proteomics)]
fap_proteins <- intersect(FAP_list, colnames(UKBB_Proteomics))
all_other_proteins <- setdiff(colnames(UKBB_Proteomics), FAP_list)
UKBB_Proteomics$eid <- rownames(UKBB_Proteomics) %>% as.character()

# --- 5. BASELINE CHARACTERISTICS --- #
case_for_table <- "dis2_inci"

# --- 6. SUM SCORE ANALYSIS FUNCTION --- #
analyze_proteins_multiyear_multicontrol <- function(
    protein_set_name,
    set_name,
    cases = "dis2_inci",
    controls = c("HC"),
    ratio = 5,
    n_obs_yrs_vec = c(2, 4, 6, 8, 10, 12, 14),
    n_bootstrap = 1000) {
  
  purrr::map_dfr(controls, function(current_control) {
    purrr::map_dfr(n_obs_yrs_vec, function(n_obs_yrs) {
      
      foreach(
        i = 1:n_bootstrap,
        .combine = rbind,
        .packages = c("dplyr", "purrr", "tibble", "stringr", "broom", "tidyr")
      ) %dopar% {
        
        set.seed(1234 + i)
        
        if (set_name == "FAP_list") {
          proteins_used <- fap_proteins
        } else if (set_name == "random") {
          proteins_used <- sample(all_other_proteins, length(fap_proteins))
        } else {
          proteins_used <- intersect(protein_set_name, colnames(UKBB_Proteomics))
        }
        
        if (length(proteins_used) == 0) {
          return(NULL)
        }
        
        protein_data <- UKBB_Proteomics[, proteins_used, drop = FALSE]
        protein_data$protein_mean <- rowMeans(protein_data, na.rm = TRUE)
        protein_data$eid <- rownames(protein_data) %>% as.character()
        
        purrr::map_dfr(cases, function(case) {
          
          case_clinic_var <- paste0(case, "_clinic")
          case_eid_var <- paste0(case, "_eid")
          
          case_data <- All_clinic %>%
            dplyr::filter(eid %in% get(case_eid_var)) %>%
            dplyr::left_join(
              protein_data %>% dplyr::select(eid, protein_mean),
              by = "eid"
            ) %>%
            dplyr::left_join(
              get(case_clinic_var) %>% dplyr::select(eid, diagnosis_time),
              by = "eid"
            )
          
          if (grepl("_inci$", case)) {
            case_data <- case_data %>%
              dplyr::filter(
                diagnosis_time > date_of_attending_assessment_centre_f53_0_0,
                diagnosis_time <= date_of_attending_assessment_centre_f53_0_0 + n_obs_yrs * 365
              )
          } else if (grepl("_prev$", case)) {
            case_data <- case_data %>%
              dplyr::filter(
                diagnosis_time <= date_of_attending_assessment_centre_f53_0_0,
                diagnosis_time >= date_of_attending_assessment_centre_f53_0_0 - n_obs_yrs * 365
              )
          } else {
            case_data <- case_data %>%
              dplyr::filter(
                abs(diagnosis_time - date_of_attending_assessment_centre_f53_0_0) <= n_obs_yrs * 365
              )
          }
          
          case_data <- case_data %>%
            dplyr::mutate(binary_group = 1)
          
          if (current_control == "HC") {
            control_eid_var <- "HC_eid"
          } else if (current_control == "Dis1noDis2") {
            control_eid_var <- "Dis1noDis2_eid"
          } else {
            stop("Invalid control group.")
          }
          
          control_data <- All_clinic %>%
            dplyr::filter(eid %in% get(control_eid_var)) %>%
            dplyr::left_join(
              protein_data %>% dplyr::select(eid, protein_mean),
              by = "eid"
            ) %>%
            dplyr::mutate(binary_group = 0)
          
          if (nrow(case_data) < 1 || nrow(control_data) < 1) {
            return(
              tibble::tibble(
                Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                OR = NA_real_,
                P = NA_real_,
                Protein_Set = set_name,
                Bootstrap = i,
                Observation_Years = n_obs_yrs,
                Control_Group = current_control,
                matched_hash = NA_character_,
                n_total = NA_real_,
                n_cases = nrow(case_data),
                n_controls = nrow(control_data),
                num_keep = NA_character_,
                fac_keep = NA_character_,
                all_terms = NA_character_,
                formula_used = NA_character_,
                n_model = NA_real_
              )
            )
          }
          
          sampled_controls <- control_data %>%
            dplyr::sample_n(size = ratio * nrow(case_data), replace = TRUE)
          
          sampled_data <- dplyr::bind_rows(case_data, sampled_controls) %>%
            dplyr::mutate(binary_group = factor(binary_group))
          
          matched_ids <- sort(as.character(sampled_data$eid))
          matched_hash <- paste(matched_ids, collapse = "|")
          n_total <- nrow(sampled_data)
          n_cases <- sum(sampled_data$binary_group == "1", na.rm = TRUE)
          n_controls <- sum(sampled_data$binary_group == "0", na.rm = TRUE)
          
          if (dplyr::n_distinct(sampled_data$binary_group) < 2) {
            return(
              tibble::tibble(
                Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                OR = NA_real_,
                P = NA_real_,
                Protein_Set = set_name,
                Bootstrap = i,
                Observation_Years = n_obs_yrs,
                Control_Group = current_control,
                matched_hash = matched_hash,
                n_total = n_total,
                n_cases = n_cases,
                n_controls = n_controls,
                num_keep = NA_character_,
                fac_keep = NA_character_,
                all_terms = NA_character_,
                formula_used = NA_character_,
                n_model = NA_real_
              )
            )
          }
          
          covars_num <- c("age", "body_mass_index_bmi_f21001_0_0", "PRS")
          covars_fac <- c(
            "sex",
            "smoking_status_f20116_0_0",
            "alcohol_intake_frequency_f1558_0_0",
            "diabetes_diagnosed_by_doctor_f2443_0_0"
          )
          
          covars_num_present <- covars_num[covars_num %in% colnames(sampled_data)]
          covars_fac_present <- covars_fac[covars_fac %in% colnames(sampled_data)]
          
          sampled_data <- sampled_data %>%
            dplyr::mutate(
              dplyr::across(all_of(covars_fac_present), as.factor)
            )
          
          if ("smoking_status_f20116_0_0" %in% covars_fac_present) {
            smoking_levels <- levels(sampled_data$smoking_status_f20116_0_0)
            if ("Never" %in% smoking_levels) {
              sampled_data$smoking_status_f20116_0_0 <- stats::relevel(
                sampled_data$smoking_status_f20116_0_0,
                ref = "Never"
              )
            }
          }
          
          if ("alcohol_intake_frequency_f1558_0_0" %in% covars_fac_present) {
            alcohol_levels <- levels(sampled_data$alcohol_intake_frequency_f1558_0_0)
            if ("Never" %in% alcohol_levels) {
              sampled_data$alcohol_intake_frequency_f1558_0_0 <- stats::relevel(
                sampled_data$alcohol_intake_frequency_f1558_0_0,
                ref = "Never"
              )
            }
          }
          
          if ("diabetes_diagnosed_by_doctor_f2443_0_0" %in% covars_fac_present) {
            diabetes_levels <- levels(sampled_data$diabetes_diagnosed_by_doctor_f2443_0_0)
            if ("No" %in% diabetes_levels) {
              sampled_data$diabetes_diagnosed_by_doctor_f2443_0_0 <- stats::relevel(
                sampled_data$diabetes_diagnosed_by_doctor_f2443_0_0,
                ref = "No"
              )
            }
          }
          
          model_data <- sampled_data %>%
            dplyr::select(
              binary_group,
              protein_mean,
              dplyr::all_of(covars_num_present),
              dplyr::all_of(covars_fac_present)
            ) %>%
            tidyr::drop_na()
          
          if (nrow(model_data) < 2 || dplyr::n_distinct(model_data$binary_group) < 2) {
            return(
              tibble::tibble(
                Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                OR = NA_real_,
                P = NA_real_,
                Protein_Set = set_name,
                Bootstrap = i,
                Observation_Years = n_obs_yrs,
                Control_Group = current_control,
                matched_hash = matched_hash,
                n_total = n_total,
                n_cases = n_cases,
                n_controls = n_controls,
                num_keep = NA_character_,
                fac_keep = NA_character_,
                all_terms = NA_character_,
                formula_used = NA_character_,
                n_model = nrow(model_data)
              )
            )
          }
          
          model_data <- model_data %>%
            dplyr::mutate(
              dplyr::across(dplyr::all_of(covars_fac_present), droplevels)
            )
          
          fac_keep <- covars_fac_present[
            sapply(model_data[covars_fac_present], function(x) nlevels(x) >= 2)
          ]
          
          num_keep <- covars_num_present[
            sapply(model_data[covars_num_present], function(x) dplyr::n_distinct(x) >= 2)
          ]
          
          formula_terms <- c("protein_mean", num_keep, fac_keep)
          formula_str <- paste("binary_group ~", paste(formula_terms, collapse = " + "))
          
          model <- tryCatch(
            glm(
              as.formula(formula_str),
              data = model_data,
              family = binomial
            ),
            error = function(e) NULL
          )
          
          if (is.null(model)) {
            return(
              tibble::tibble(
                Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                OR = NA_real_,
                P = NA_real_,
                Protein_Set = set_name,
                Bootstrap = i,
                Observation_Years = n_obs_yrs,
                Control_Group = current_control,
                matched_hash = matched_hash,
                n_total = n_total,
                n_cases = n_cases,
                n_controls = n_controls,
                num_keep = paste(num_keep, collapse = ";"),
                fac_keep = paste(fac_keep, collapse = ";"),
                all_terms = paste(formula_terms, collapse = ";"),
                formula_used = formula_str,
                n_model = nrow(model_data)
              )
            )
          }
          
          tidy_model <- broom::tidy(model)
          protein_row <- dplyr::filter(tidy_model, term == "protein_mean")
          
          tibble::tibble(
            Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
            OR = ifelse(nrow(protein_row) == 1, exp(protein_row$estimate), NA_real_),
            P = ifelse(nrow(protein_row) == 1, protein_row$p.value, NA_real_),
            Protein_Set = set_name,
            Bootstrap = i,
            Observation_Years = n_obs_yrs,
            Control_Group = current_control,
            matched_hash = matched_hash,
            n_total = n_total,
            n_cases = n_cases,
            n_controls = n_controls,
            num_keep = paste(num_keep, collapse = ";"),
            fac_keep = paste(fac_keep, collapse = ";"),
            all_terms = paste(formula_terms, collapse = ";"),
            formula_used = formula_str,
            n_model = nrow(model_data)
          )
        })
      }
    })
  })
}

# --- 7. EXECUTE ANALYSES AND MERGE RESULTS --- #
runs_list <- list(
  list(p_set = FAP_list, set_name = "FAP_list", controls = c("HC"), cases = "dis2_inci"),
  list(p_set = "random", set_name = "random", controls = c("HC"), cases = "dis2_inci"),
  list(p_set = known_marker, set_name = "CEACAM5", controls = c("HC"), cases = "dis2_inci"),
  list(p_set = FAP_list, set_name = "FAP_list", controls = c("Dis1noDis2"), cases = c("Dis1andDis2")),
  list(p_set = "random", set_name = "random", controls = c("Dis1noDis2"), cases = c("Dis1andDis2")),
  list(p_set = known_marker, set_name = "CEACAM5", controls = c("Dis1noDis2"), cases = c("Dis1andDis2"))
)

start_time <- Sys.time()
results_list <- list()

for (r in runs_list) {
  res <- analyze_proteins_multiyear_multicontrol(
    protein_set_name = r$p_set,
    set_name = r$set_name,
    controls = r$controls,
    cases = r$cases,
    n_obs_yrs_vec = observation_years,
    n_bootstrap = n_bootstrap
  )
  results_list[[length(results_list) + 1]] <- res
}

end_time <- Sys.time()
print(paste("Analysis completed in", round(end_time - start_time, 1), units(end_time - start_time)))

results <- dplyr::bind_rows(results_list)

# --- 8. SUMMARIZE RESULTS (Bootstrap Aggregation) --- #
summary_df <- results %>%
  dplyr::mutate(contrast = paste0(Protein_Set, " - ", Group, " vs ", Control_Group)) %>%
  dplyr::group_by(Group, Protein_Set, Observation_Years, Control_Group, contrast) %>%
  dplyr::summarise(
    n_non_na = sum(!is.na(OR)),
    OR_mean = ifelse(n_non_na > 0, mean(OR, na.rm = TRUE), NA_real_),
    OR_lower = ifelse(n_non_na > 0, quantile(OR, 0.025, na.rm = TRUE), NA_real_),
    OR_upper = ifelse(n_non_na > 0, quantile(OR, 0.975, na.rm = TRUE), NA_real_),
    P_median = ifelse(n_non_na > 0, median(P, na.rm = TRUE), NA_real_),
    n_cases_mean = mean(n_cases, na.rm = TRUE),
    n_controls_mean = mean(n_controls, na.rm = TRUE),
    n_model_mean = mean(n_model, na.rm = TRUE),
    n_boot_valid = sum(!is.na(OR)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Control_Group = factor(
      Control_Group,
      levels = c("HC", "Dis1noDis2"),
      labels = c("Healthy Controls", "Dis1 No Progression")
    )
  )

# --- 9. VISUALIZATION (OR PLOTS) --- #
plot_data_incident <- summary_df %>%
  dplyr::filter(Group == "dis2_inci", Control_Group != "Dis1 No Progression")

or_plot_incident <- ggplot(
  plot_data_incident,
  aes(
    x = Observation_Years,
    y = log2(OR_mean),
    color = Protein_Set,
    group = Protein_Set
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey60") +
  geom_errorbar(
    aes(ymin = log2(OR_lower), ymax = log2(OR_upper)),
    width = 0.18,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  labs(
    x = "Years before diagnosis",
    y = expression(log[2] ~ "odds ratio"),
    color = NULL
  ) +
  scale_color_manual(
    name = "Model:",
    values = c(
      "FAP_list" = "#FC7F0E",
      "random" = "#7F7F7F",
      "CEACAM5" = "#5DADE2"
    ),
    labels = c(
      "FAP_list" = "FdS",
      "random" = "Random set",
      "CEACAM5" = "CEA"
    )
  ) +
  scale_x_continuous(
    limits = c(1.9, 14.1),
    breaks = observation_years
  ) +
  scale_y_continuous(
    limits = c(-2.5, 6),
    breaks = c(-2, 0, 2, 4, 6),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.title = element_text(face = "plain"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black")
  )

or_plot_incident

