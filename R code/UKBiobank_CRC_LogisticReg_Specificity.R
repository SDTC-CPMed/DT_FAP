# ==============================================================================
# Project -> UK Biobank Proteomics Specificity to CRC Analysis
#
# Description:
# This script evaluates proteomic specificity across three comparisons:
#   1) CRC vs Healthy controls
#   2) Other cancers vs Healthy controls
#   3) CRC vs Other cancers
# 
# ==============================================================================

rm(list = ls())

# --- 0. PACKAGE LOADING --- #
packages <- c(
  "dplyr",
  "ggplot2",
  "purrr",
  "broom",
  "foreach",
  "doParallel",
  "tibble"
)

invisible(lapply(packages, library, character.only = TRUE))

# --- 1. CONFIGURATION AND DIRECTORY SETUP --- #
# Set input/output directories and analysis parameters here.
input_data_dir <- "path/to/input_data"
subset_data_dir <- "path/to/subset_data"
output_dir <- "output_specificity"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

target_disease <- "Colorectal"
n_bootstrap <- 1000
observation_years <- c(2, 4, 6, 8, 10, 12, 14)
random_seed <- 123
n_cores <- 10

set.seed(random_seed)
registerDoParallel(cores = n_cores)

# --- 2. DATA LOADING --- #

UKBB_Proteomics <- read.table(
  file = file.path(input_data_dir, "Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt"),
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  row.names = "PID",
  check.names = FALSE)

All_clinic <- readRDS(file.path(input_data_dir, "ukb_clinic.rda")) %>%
  filter(eid %in% rownames(UKBB_Proteomics)) %>%
  mutate(across(
    c(date_of_attending_assessment_centre_f53_0_0,
      date_of_death_f40000_0_0,
      date_of_death_f40000_1_0,
      date_lost_to_followup_f191_0_0),
    as.Date)) %>%
  mutate(
    sex_f31_0_0 = as.factor(sex_f31_0_0),
    ethnic_background_f21000_0_0 = as.factor(ethnic_background_f21000_0_0),
    delta_diag_enroll = NA,
    OSstatus.raw = ifelse(is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0),0, 1),
    OStime = ifelse(OSstatus.raw == 1,
                    (date_of_death_f40000_0_0 - date_of_attending_assessment_centre_f53_0_0) / 365,
                    ifelse(
                      is.na(date_lost_to_followup_f191_0_0),
                      (as.Date("2022-12-30") - date_of_attending_assessment_centre_f53_0_0) / 365,
                      (pmin(date_lost_to_followup_f191_0_0, as.Date("2022-12-30")) -
                         date_of_attending_assessment_centre_f53_0_0) / 365))) %>%
  filter(
    !is.na(ethnic_background_f21000_0_0),
    !is.na(age_at_recruitment_f21022_0_0),
    !is.na(uk_biobank_assessment_centre_f54_0_0),
    !is.na(sex_f31_0_0)) %>%
  mutate(
    eid = as.character(eid),
    age = age_at_recruitment_f21022_0_0,
    sex = sex_f31_0_0)

# Function to load disease-specific clinical data
load_disease_data <- function(dis_name) {
  readRDS(file.path(subset_data_dir, paste0("ukb_", dis_name, "_subset_clinic.rda"))) %>%
    filter(eid %in% rownames(UKBB_Proteomics)) %>%
    mutate(
      eid = as.character(eid),
      diagnosis_time = as.Date(diagnosis_time))
}

disC_clinic <- load_disease_data(target_disease)

AllCancer_clinic <- readRDS(file.path(subset_data_dir, "ukb_AllCancers_subset_clinic.rda")) %>%
  filter(eid %in% rownames(UKBB_Proteomics)) %>%
  mutate(
    eid = as.character(eid),
    diagnosis_time = as.Date(diagnosis_time))

# Keep earliest cancer diagnosis per individual
AllCancer_clinic <- AllCancer_clinic %>%
  group_by(eid) %>%
  slice_min(order_by = diagnosis_time, n = 1, with_ties = FALSE) %>%
  ungroup()

disC_clinic <-  disC_clinic %>%
  group_by(eid) %>%
  slice_min(order_by = diagnosis_time, n = 1, with_ties = FALSE) %>%
  ungroup()

disC_eid <-  disC_clinic %>% pull(eid)

othercancer_clinic <- AllCancer_clinic %>%
  filter(!(eid %in%  disC_eid))

# Function to extract EIDs and clinic data for prevalent and incident cases
extract_cases <- function(clinic_df, prevalent = TRUE) {
  tmp <- clinic_df %>%
    dplyr::select(-any_of("date_of_attending_assessment_centre_f53_0_0")) %>%
    left_join(
      All_clinic %>% dplyr::select(eid, date_of_attending_assessment_centre_f53_0_0),
      by = "eid")
  
  if (prevalent) {
    tmp %>% filter(diagnosis_time <= date_of_attending_assessment_centre_f53_0_0)
  } else {
    tmp %>% filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0)
  }}

disC_prev_clinic <- extract_cases( disC_clinic, prevalent = TRUE)
disC_inci_clinic <- extract_cases( disC_clinic, prevalent = FALSE)

othercancer_prev_clinic <- extract_cases(othercancer_clinic, prevalent = TRUE)
othercancer_inci_clinic <- extract_cases(othercancer_clinic, prevalent = FALSE)

disC_prev_eid <-  disC_prev_clinic %>% pull(eid)
disC_inci_eid <-  disC_inci_clinic %>% pull(eid)

othercancer_prev_eid <- othercancer_prev_clinic %>% pull(eid)
othercancer_inci_eid <- othercancer_inci_clinic %>% pull(eid)

# Load and prepare Healthy Controls (HC) EIDs
HC_clinic_raw <- readRDS(file.path(subset_data_dir, "ukb_healthy_controls.rda")) %>%
  filter(eid %in% rownames(UKBB_Proteomics))

HC_clinic <- All_clinic %>%
  filter(eid %in% as.character(HC_clinic_raw$eid))

HC_eid <- HC_clinic %>% pull(eid)

# --- 3. PROTEIN SETS --- #

known_marker <- c("CEACAM5")

FAP_list <- c(
  "LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F", "INHBA",
  "LAMA5", "COL7A1", "SOX9", "CCND1", "MET", "ID1", "PTMA", "CLU", "TGIF1",
  "EDN1", "GDF15")

UKBB_Proteomics <- UKBB_Proteomics[,colSums(is.na(UKBB_Proteomics)) <= 0.25 * nrow(UKBB_Proteomics)]
fap_proteins <- intersect(FAP_list, colnames(UKBB_Proteomics))
all_other_proteins <- setdiff(colnames(UKBB_Proteomics), FAP_list)
UKBB_Proteomics$eid <- rownames(UKBB_Proteomics) %>% as.character()

# --- 4. CLINICAL OBJECTS --- #

clinic_list <- list(
  disC_prev =  disC_prev_clinic,
  disC_inci =  disC_inci_clinic,
  othercancer_prev = othercancer_prev_clinic,
  othercancer_inci = othercancer_inci_clinic,
  HC = HC_clinic)

eid_list <- list(
  disC_prev =  disC_prev_eid,
  disC_inci =  disC_inci_eid,
  othercancer_prev = othercancer_prev_eid,
  othercancer_inci = othercancer_inci_eid,
  HC = HC_eid)

# --- 5. ANALYSIS FUNCTION --- #
# This function calculates the mean protein score, performs control resampling, and runs simple logistic regression.
analyze_specificity_or  <- function(
    protein_set_name,
    set_name,
    case_control_pairs,
    clinic_list,
    eid_list,
    ratio = 1,
    n_obs_yrs_vec = c(2, 4, 6, 8, 10, 12, 14),
    n_bootstrap = 1000) {
  
  map_dfr(n_obs_yrs_vec, function(n_obs_yrs) {
    
    foreach(
      i = 1:n_bootstrap,
      .combine = rbind,
      .packages = c("dplyr", "purrr", "tibble", "broom", "tidyr")) %dopar% {
        
        set.seed(1234 + i) # Set seed for random processes within bootstrap
        
        # Protein mean calculation
        if (set_name == "FAP_list") {
          proteins_used <- fap_proteins
        } else if (set_name == "random") {
          # Randomly sample proteins of the same size as FAP_list for null comparison
          proteins_used <- sample(all_other_proteins, length(fap_proteins))
        } else {
          proteins_used <- intersect(protein_set_name, colnames(UKBB_Proteomics))
        }
        
        if (length(proteins_used) == 0) {
          return(
            map_dfr(seq_len(nrow(case_control_pairs)), function(k) {
              tibble(
                Group = case_control_pairs$case[k],
                OR = NA_real_,
                P = NA_real_,
                Protein_Set = set_name,
                Bootstrap = i,
                Observation_Years = n_obs_yrs,
                Control_Group = case_control_pairs$control[k],
                Comparison = case_control_pairs$comparison[k])}))
        }
        
        # Calculate row mean score across the selected protein set (Sum Score)
        
        protein_data <- UKBB_Proteomics[, proteins_used, drop = FALSE]
        protein_data$protein_mean <- rowMeans(protein_data, na.rm = TRUE)
        protein_data$eid <- rownames(protein_data) %>% as.character()
        
        
        # Process Case and Control Groups
        map_dfr(seq_len(nrow(case_control_pairs)), function(k) {
          
          case_name <- case_control_pairs$case[k]
          control_name <- case_control_pairs$control[k]
          comparison_name <- case_control_pairs$comparison[k]
          
          case_clinic_df <- clinic_list[[case_name]]
          case_eids <- eid_list[[case_name]]
          
          case_data <- All_clinic %>%
            dplyr::filter(eid %in% case_eids) %>%
            dplyr::left_join(protein_data %>% dplyr::select(eid, protein_mean), by = "eid") %>%
            dplyr::left_join(case_clinic_df %>% dplyr::select(eid, diagnosis_time), by = "eid") %>%
            dplyr::filter(abs(diagnosis_time - date_of_attending_assessment_centre_f53_0_0) <= n_obs_yrs * 365) %>%
            dplyr::mutate(binary_group = 1)
          
          if (nrow(case_data) < 1) {
            return(
              tibble(
                Group = case_name,
                OR = NA_real_,
                P = NA_real_,
                Protein_Set = set_name,
                Bootstrap = i,
                Observation_Years = n_obs_yrs,
                Control_Group = control_name,
                Comparison = comparison_name))
          }
          
          if (control_name == "HC") {
            
            control_data <- All_clinic %>%
              dplyr::filter(eid %in% HC_eid) %>%
              dplyr::left_join(protein_data %>% dplyr::select(eid, protein_mean), by = "eid") %>%
              dplyr::mutate(binary_group = 0)
            
            sampled_data <- bind_rows(
              case_data,
              control_data %>% sample_n(nrow(case_data), replace = TRUE)) %>%
              dplyr::mutate(binary_group = factor(binary_group))
            
          } else {
            
            control_clinic_df <- clinic_list[[control_name]]
            control_eids <- eid_list[[control_name]]
            
            control_data <- All_clinic %>%
              dplyr::filter(eid %in% control_eids) %>%
              dplyr::left_join(protein_data %>% dplyr::select(eid, protein_mean), by = "eid") %>%
              dplyr::left_join(control_clinic_df %>% dplyr::select(eid, diagnosis_time), by = "eid") %>%
              dplyr::filter(abs(diagnosis_time - date_of_attending_assessment_centre_f53_0_0) <= n_obs_yrs * 365) %>%
              dplyr::mutate(binary_group = 0)
            
            if (nrow(control_data) < 1) {
              return(
                tibble(
                  Group = case_name,
                  OR = NA_real_,
                  P = NA_real_,
                  Protein_Set = set_name,
                  Bootstrap = i,
                  Observation_Years = n_obs_yrs,
                  Control_Group = control_name,
                  Comparison = comparison_name))
            }
            
            sampled_data <- bind_rows(
              case_data,
              control_data %>% sample_n(nrow(case_data), replace = TRUE)) %>%
              dplyr::mutate(binary_group = factor(binary_group))
          }
          
          # Fit logistic regression model (adjusted for Age, Sex, BMI, and smoking status)
          if (
            n_distinct(sampled_data$binary_group) == 2 &&
            sum(!is.na(sampled_data$protein_mean)) > 0) {
            model <- glm(
              binary_group ~ protein_mean + sex + age + body_mass_index_bmi_f21001_0_0 +
                smoking_status_f20116_0_0,
              data = sampled_data,
              family = binomial,
              na.action = na.exclude)
            
            tidy_model <- broom::tidy(model)
            protein_row <- dplyr::filter(tidy_model, term == "protein_mean")
            
            # Compile bootstrap result
            tibble(
              Group = case_name,
              OR = ifelse(nrow(protein_row) == 1, exp(protein_row$estimate), NA_real_),
              P = ifelse(nrow(protein_row) == 1, protein_row$p.value, NA_real_),
              Protein_Set = set_name,
              Bootstrap = i,
              Observation_Years = n_obs_yrs,
              Control_Group = control_name,
              Comparison = comparison_name)
          } else {
            tibble(
              Group = case_name,
              OR = NA_real_,
              P = NA_real_,
              Protein_Set = set_name,
              Bootstrap = i,
              Observation_Years = n_obs_yrs,
              Control_Group = control_name,
              Comparison = comparison_name)
          }
        })
      }
  })
}

# --- 6. COMPARISON DEFINITIONS --- #

case_control_pairs <- tibble::tribble(
  ~case,              ~control,             ~comparison,
  "disC_inci",        "HC",                 "CRC_vs_HC",
  "othercancer_inci", "HC",                 "OtherCancers_vs_HC",
  "disC_inci",        "othercancer_inci",   "CRC_vs_OtherCancers")

# --- 7. RUN ANALYSES --- #

runs_list <- list(
  list(p_set = FAP_list,     set_name = "FAP_list"),
  list(p_set = "random",     set_name = "random"),
  list(p_set = known_marker, set_name = "CEACAM5"))

start_time <- Sys.time()
results_list <- list()

for (r in runs_list) {
  res <- analyze_specificity_or (
    protein_set_name = r$p_set,
    set_name = r$set_name,
    case_control_pairs = case_control_pairs,
    clinic_list = clinic_list,
    eid_list = eid_list,
    n_obs_yrs_vec = observation_years,
    n_bootstrap = n_bootstrap)
  results_list[[length(results_list) + 1]] <- res
}

end_time <- Sys.time()
print(paste("Analysis completed in",round(end_time - start_time, 1), units(end_time - start_time)))

results <- bind_rows(results_list)


# --- 8. SUMMARIZE RESULTS --- #

summary_df <- results %>%
  mutate(
    contrast = paste0(Protein_Set, " - ", Group, " vs ", Control_Group),
    Timing = case_when(
      Group %in% c(" disC_prev", "othercancer_prev") ~ "Prevalent",
      Group %in% c(" disC_inci", "othercancer_inci") ~ "Incident",
      TRUE ~ "Other"),
    Case_Group_Label = case_when(
      Group %in% c(" disC_prev", " disC_inci") ~ "CRC",
      Group %in% c("othercancer_prev", "othercancer_inci") ~ "Other cancers",
      TRUE ~ Group),
    Comparison_Label = case_when(
      Comparison == "CRC_vs_HC" ~ "CRC vs HC",
      Comparison == "OtherCancers_vs_HC" ~ "Other cancers vs HC",
      Comparison == "CRC_vs_OtherCancers" ~ "CRC vs Other cancers",
      TRUE ~ Comparison)) %>%
  group_by(
    Comparison, Comparison_Label, Group, Case_Group_Label,
    Protein_Set, Observation_Years, Control_Group, Timing, contrast) %>%
  summarise(
    OR_mean = mean(OR, na.rm = TRUE),
    OR_lower = quantile(OR, 0.025, na.rm = TRUE),
    OR_upper = quantile(OR, 0.975, na.rm = TRUE),
    P_median = median(P, na.rm = TRUE),
    n_boot_valid = sum(!is.na(OR)),
    .groups = "drop") %>%
  mutate(
    log2_OR_mean  = log2(OR_mean),
    log2_OR_lower = log2(OR_lower),
    log2_OR_upper = log2(OR_upper))

mean_or_range_df <- summary_df %>%
  group_by(Comparison, Comparison_Label, Timing, Protein_Set) %>%
  summarise(
    min_OR_mean = min(OR_mean, na.rm = TRUE),
    max_OR_mean = max(OR_mean, na.rm = TRUE),
    min_log2_OR_mean = min(log2_OR_mean, na.rm = TRUE),
    max_log2_OR_mean = max(log2_OR_mean, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(Comparison_Label, Timing, Protein_Set)

print(mean_or_range_df)
print(summary_df)


# --- 9. VISUALIZATION --- #

plot_colors <- c(
  "FAP_list" = "#FC7F0E",
  "random"   = "#7F7F7F",
  "CEACAM5"  = "#5DADE2")

summary_all_panels <- summary_df %>%
  filter(
    Comparison %in% c("CRC_vs_HC", "OtherCancers_vs_HC", "CRC_vs_OtherCancers")) %>%
  mutate(
    Comparison_Label = recode(
      Comparison,
      "CRC_vs_HC" = "CRC vs HC",
      "OtherCancers_vs_HC" = "Other cancers vs HC",
      "CRC_vs_OtherCancers" = "CRC vs Other cancers"),
    Comparison_Label = factor(
      Comparison_Label,
      levels = c("CRC vs HC", "Other cancers vs HC", "CRC vs Other cancers")))

or_plot_3panel <- ggplot(
  summary_all_panels,
  aes(
    x = Observation_Years,
    y = log2_OR_mean,
    color = Protein_Set,
    group = Protein_Set)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log2_OR_lower, ymax = log2_OR_upper), width = 0.2) +
  facet_grid(~Comparison_Label, scales = "fixed") +
  labs(
    x = "Observation years before diagnosis",
    y = expression(log[2] ~ "Odds Ratio"),
    color = "Protein Set") +
  scale_color_manual(
    name = "Model:",
    values = plot_colors,
    labels = c(
      "FAP_list" = "FdS",
      "random"   = "Random set",
      "CEACAM5"  = "CEA")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_continuous(limits = c(1.9, 14.1), breaks = observation_years) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"))

print(or_plot_3panel)


