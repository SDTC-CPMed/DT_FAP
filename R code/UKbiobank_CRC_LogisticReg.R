# ==============================================================================
# PROJECT: UK Biobank Proteomics Risk Prediction Analysis (Sum Score Approach)
# DESCRIPTION:
# This script performs protein-based risk prediction using the mean protein level
# (Sum Score) derived from a defined Gene of Interest (GOI) list.
# The analysis involves parallelized bootstrapping,
# logistic regression adjusted for age and sex, and calculation of OR
#
# METHOD: logistic regression.
# ==============================================================================

# --- 0. INITIAL SETUP AND PACKAGE LOADING --- #
# Define and load necessary packages
packages <- c("STRINGdb", "org.Hs.eg.db", "igraph", "dplyr", "survival", 
              "Matrix", "ggplot2", "patchwork", "RColorBrewer", "reshape2", 
              "ComplexHeatmap", "circlize", "viridis", "readxl", "stringr", 
              "data.table", "glmnet", "reticulate", "broom", "MatchIt", "VIM", 
              "purrr", "pROC", "tidyverse", "foreach", "doParallel", "writexl", 
              "tableone") # Added writexl and tableone for Baseline section

# Load all packages 
lapply(packages, library, character.only = TRUE)

# --- 1. CONFIGURATION AND DIRECTORY SETUP --- #
input.dir <- '/UKBB/input' # Assumed input path
out.dir <- 'output'
dir.create(out.dir, showWarnings = FALSE)

# Set analysis parameters
dis1 <- 'Colon_premalignant' 
dis2 <- 'Colorectal_cancer' 
cases <- c('dis2_prev', 'dis2_inci')
control <- 'HC'  
goi <- 'FAP_list' # Gene of interest set name
n_bootstrap <- 1000
set.seed(123)
# Set up parallel processing (adjust cores as needed)
registerDoParallel(cores = 10) 


# --- 2. DATA LOADING AND CLINICAL SUBSETTING --- #
UKBB_Proteomics = read.table(file = file.path(input.dir, 'Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt'), sep='\t', header = TRUE, fill = TRUE, row.names = 'PID')

# Load and process core clinical data
All_clinic <- readRDS(file.path(input.dir, 'ukb_clinic.rda')) %>% 
  filter(eid %in% rownames(UKBB_Proteomics)) %>% 
  mutate(across(c(date_of_attending_assessment_centre_f53_0_0, date_of_death_f40000_0_0, date_of_death_f40000_1_0, date_lost_to_followup_f191_0_0), as.Date)) %>%
  mutate(sex_f31_0_0 = as.factor(sex_f31_0_0),
         ethnic_background_f21000_0_0 = as.factor(ethnic_background_f21000_0_0),
         delta_diag_enroll = NA,
         OSstatus.raw = ifelse(is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0), 0, 1),
         OStime = ifelse(OSstatus.raw == 1, 
                         (date_of_death_f40000_0_0 - date_of_attending_assessment_centre_f53_0_0)/365,
                         ifelse(is.na(date_lost_to_followup_f191_0_0),
                                (as.Date("2022-12-30") - date_of_attending_assessment_centre_f53_0_0)/365,
                                (pmin(date_lost_to_followup_f191_0_0, as.Date("2022-12-30")) - date_of_attending_assessment_centre_f53_0_0)/365))) %>% 
  filter(!is.na(ethnic_background_f21000_0_0), !is.na(age_at_recruitment_f21022_0_0),
         !is.na(uk_biobank_assessment_centre_f54_0_0), !is.na(sex_f31_0_0)) %>% 
  mutate(eid = as.character(eid))

# Function to load disease-specific clinical data
load_disease_data <- function(dis_name) {
  readRDS(file.path(input.dir, paste0('ukb_', dis_name, '_subset_clinic.rda'))) %>% 
    filter(eid %in% rownames(UKBB_Proteomics)) %>% 
    mutate(eid = as.character(eid))
}

dis1_clinic <- load_disease_data(dis1)
dis2_clinic <- load_disease_data(dis2)
All_cancer_clinic_eid <- readRDS(file.path(input.dir,'ukb_AllCancers_subset_clinic.rda')) %>% filter(eid %in% rownames(UKBB_Proteomics)) %>% pull(eid)
dis1_eid <- dis1_clinic %>% dplyr::select(eid, diagnosis_time) %>% pull(eid)
dis2_eid <- dis2_clinic %>% dplyr::select(eid, diagnosis_time) %>% pull(eid)

# Function to extract EIDs and clinic data for prevalent and incident cases
extract_cases <- function(clinic_df, prevalent = TRUE) {
  if (prevalent) {
    clinic_df %>% filter(diagnosis_time <= date_of_attending_assessment_centre_f53_0_0)
  } else {
    clinic_df %>% filter(diagnosis_time > date_of_attending_assessment_centre_f53_0_0)
  }
}

dis1_prev_clinic <- extract_cases(dis1_clinic, prevalent = TRUE)
dis1_inci_clinic <- extract_cases(dis1_clinic, prevalent = FALSE)
dis2_prev_clinic <- extract_cases(dis2_clinic, prevalent = TRUE)
dis2_inci_clinic <- extract_cases(dis2_clinic, prevalent = FALSE)

dis1_prev_eid <- dis1_prev_clinic %>% pull(eid)
dis1_inci_eid <- dis1_inci_clinic %>% pull(eid)
dis2_prev_eid <- dis2_prev_clinic %>% pull(eid)
dis2_inci_eid <- dis2_inci_clinic %>% pull(eid)

# Load and prepare Healthy Controls (HC) EIDs
HC_clinic <- readRDS(file.path(input.dir,'ukb_healthy_controls.rda')) %>% filter(eid %in% rownames(UKBB_Proteomics))
HC_clinic <- All_clinic %>% filter(eid %in% as.character(HC_clinic$eid))
HC_eid <- HC_clinic %>% filter(eid %in% rownames(UKBB_Proteomics)) %>% pull(eid)

# --- Subgroup definition for pre-malignant progression (Case-vs-Case) ---

# Merge Dis1 and Dis2 diagnosis times for progression analysis
All_cancer_clinic_2 <- All_clinic %>% 
  left_join(dis1_clinic %>% select(eid, diagnosis_time.dis1 = diagnosis_time), by = 'eid') %>% 
  left_join(dis2_clinic %>% select(eid, diagnosis_time.dis2 = diagnosis_time), by = 'eid') %>% 
  mutate(across(starts_with("diagnosis_time"), as.Date))

# Define Dis1 progressing to Dis2 (CRC)
Dis1andDis2_clinic <- All_cancer_clinic_2 %>% 
  mutate(CRC = ifelse(eid %in% dis2_eid, 1, 0)) %>%
  filter(!is.na(diagnosis_time.dis1), 
         is.na(diagnosis_time.dis2) | diagnosis_time.dis1 <= diagnosis_time.dis2, # Dis1 before or concurrent with Dis2
         CRC == 1) %>% 
  rename(diagnosis_time = diagnosis_time.dis2) # Use CRC diagnosis time for follow-up
Dis1andDis2_eid <- Dis1andDis2_clinic %>% pull(eid)

# Define Dis1 NOT progressing to Dis2 (No CRC)
Dis1noDis2_clinic <- All_cancer_clinic_2 %>% 
  mutate(CRC = ifelse(eid %in% dis2_eid, 1, 0)) %>%
  filter(!is.na(diagnosis_time.dis1), 
         CRC == 0) %>% # No CRC diagnosis
  rename(diagnosis_time = diagnosis_time.dis1)
Dis1noDis2_eid <- Dis1noDis2_clinic %>% pull(eid)


# --- 3. PROTEIN SET DEFINITION --- #
known_marker <- c('CEACAM5') 
FAP_list <- c("LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F", "INHBA", 
              "LAMA5", "COL7A1", "SOX9", "CCND1", "MET", "ID1", "PTMA", "CLU", "TGIF1", 
              "EDN1", "GDF15") 

# Prepare protein sets
UKBB_Proteomics <- UKBB_Proteomics[, colSums(is.na(UKBB_Proteomics)) <= 0.25 * nrow(UKBB_Proteomics)]
fap_proteins <- intersect(FAP_list, colnames(UKBB_Proteomics))
all_other_proteins <- setdiff(colnames(UKBB_Proteomics), FAP_list)
UKBB_Proteomics$eid <- rownames(UKBB_Proteomics) %>% as.character()

# --- 4. BASELINE CHARACTERISTICS  --- #
library(tableone)
library(writexl)

case_for_table <- 'dis2' # Example case group for table generation
if (case_for_table %in% c('HC')){
  data <- All_clinic %>% filter(eid %in% get(paste0(case_for_table, '_eid')))
  # Define clinical variables...
} else if (case_for_table %in% c('dis1','dis2')){
  # Prepare data for prev vs inci table...
}
write_xlsx(table_df, paste0(out.dir, "/Baseline_characteristics_",case_for_table,".xlsx"))


# --- 5. SUM SCORE ANALYSIS FUNCTION  --- #

# This function calculates the mean protein score, performs control resampling 
# (instead of MatchIt), and runs simple logistic regression (not penalized).
analyze_proteins_multiyear_multicontrol <- function(
    protein_set_name, set_name,
    cases = c('dis2_inci', 'dis2_prev'),
    controls = c('HC'),
    ratio = 5,
    n_obs_yrs_vec = c(2, 4, 6, 8, 10, 12, 14),
    n_bootstrap = 1000) {
  
  # Run analysis for each control group
  map_dfr(controls, function(current_control) {
    # Run analysis for each observation period
    map_dfr(n_obs_yrs_vec, function(n_obs_yrs) {
      
      # Parallel processing for bootstrapping
      foreach(i = 1:n_bootstrap, .combine = rbind,
              .packages = c('dplyr','purrr','tibble','pROC','stringr')) %dopar% {
                
                set.seed(1234 + i) # Set seed for random processes within bootstrap
                
                # --- 5.1. Protein Mean Calculation ---
                if (set_name == "FAP_list") {
                  proteins_used <- fap_proteins
                } else if (set_name == "random") {
                  # Randomly sample proteins of the same size as FAP_list for null comparison
                  proteins_used <- sample(all_other_proteins, length(fap_proteins))
                } else {
                  proteins_used <- intersect(protein_set_name, colnames(UKBB_Proteomics))
                }
                
                # Calculate row mean score across the selected protein set (Sum Score)
                protein_data <- UKBB_Proteomics[, proteins_used, drop = FALSE]
                protein_data$protein_mean <- rowMeans(protein_data, na.rm = TRUE)
                protein_data$eid <- rownames(protein_data) %>% as.character()
                
                # --- 5.2. Process Case and Control Groups ---
                map_dfr(cases, function(case) {
                  
                  # Retrieve clinic data object name for the current case
                  case_clinic_var <- paste0(case, '_clinic')
                  case_eid_var <- paste0(case, '_eid')
                  
                  case_data <- All_clinic %>% 
                    filter(eid %in% get(case_eid_var)) %>% 
                    left_join(protein_data %>% select(eid, protein_mean), by = "eid") %>% 
                    left_join(get(case_clinic_var) %>% select(eid, diagnosis_time), by = 'eid') %>% 
                    # Filter cases: diagnosis must be within the observation window
                    filter(abs(diagnosis_time - date_of_attending_assessment_centre_f53_0_0) <= n_obs_yrs*365) %>% 
                    mutate(binary_group = 1)
                  
                  # Retrieve clinic data object name for the current control
                  if (current_control == 'HC') {
                    control_eid_var <- 'HC_eid'
                  } else if (current_control == 'Dis1noDis2') {
                    control_eid_var <- 'Dis1noDis2_eid'
                  } else {
                    stop("Invalid control group.") # Should be caught by the run list configuration
                  }
                  
                  control_data <- All_clinic %>% 
                    filter(eid %in% get(control_eid_var)) %>% 
                    # Note: Survival filter is removed here compared to the Ridge analysis, 
                    # but is handled in the overall HC definition in Section 2.
                    left_join(protein_data %>% select(eid, protein_mean), by = "eid") %>% 
                    mutate(binary_group = 0) 
                  
                  # --- 5.3. Resampling and Modeling ---
                  if(nrow(case_data) < 1) return(NULL)
                  match_obj <- matchit(binary_group ~ age,
                                       data = bind_rows(case_data, control_data),
                                       method = "nearest", ratio = ratio, replace = TRUE,
                                       exact = ~ sex, caliper = c(age = 2))
                  
                  sampled_data <- match.data(match_obj) %>%
                    select(eid, binary_group, age, sex, all_of(colnames(UKBB_Proteomics))) %>%
                    mutate(binary_group = factor(binary_group))
                  
                  
                  # Fit logistic regression model (adjusted for Age and Sex)
                  if (n_distinct(sampled_data$binary_group) == 2 && sum(!is.na(sampled_data$protein_mean)) > 0) {
                    model <- glm(binary_group ~ protein_mean + sex_f31_0_0 + age_at_recruitment_f21022_0_0,
                                 data = sampled_data, family = binomial, na.action = na.exclude)
                    
                    tidy_model <- broom::tidy(model)
                    protein_row <- filter(tidy_model, term == "protein_mean")
                    
                    
                    # Compile bootstrap result
                    tibble(
                      Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                      OR = ifelse(nrow(protein_row) == 1, exp(protein_row$estimate), NA),
                      P = ifelse(nrow(protein_row) == 1, protein_row$p.value, NA),
                      Protein_Set = set_name,
                      Bootstrap = i,
                      Observation_Years = n_obs_yrs,
                      Control_Group = current_control
                    )
                  } else {
                    # Return NA if modeling fails due to insufficient data
                    tibble(
                      Group = paste0(strsplit(case, "_")[[1]][1], "_", strsplit(case, "_")[[1]][2]),
                      OR = NA, 
                      Protein_Set = set_name,
                      Bootstrap = i,
                      Observation_Years = n_obs_yrs,
                      Control_Group = current_control
                    )
                  }
                })
              }
    })
  })
}


# --- 6. EXECUTE ANALYSES AND MERGE RESULTS --- #

# Define run parameters
observation_years <- c(2, 4, 6, 8, 10, 12, 14) 

# Run list focuses on HC for primary comparisons, and Dis1noDis2 for subgroup comparison
runs_list <- list(
  # Primary comparisons vs Healthy Controls (HC)
  list(p_set = FAP_list, set_name = "FAP_list", controls = c('HC'), cases = c('dis2_inci', 'dis2_prev')),
  list(p_set = "random", set_name = "random", controls = c('HC'), cases = c('dis2_inci', 'dis2_prev')),
  list(p_set = known_marker, set_name = 'CEACAM5', controls = c('HC'), cases = c('dis2_inci', 'dis2_prev')),
  
  # Subgroup progression analysis (Case vs Case)
  list(p_set = FAP_list, set_name = "FAP_list", controls = c('Dis1noDis2'), cases = c('Dis1andDis2')),
  list(p_set = "random", set_name = "random", controls = c('Dis1noDis2'), cases = c('Dis1andDis2')),
  list(p_set = known_marker, set_name = 'CEACAM5', controls = c('Dis1noDis2'), cases = c('Dis1andDis2'))
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

# Merge and save full results
results <- bind_rows(results_list)
saveRDS(results, file.path(out.dir, paste0("full_results_bt",n_bootstrap,"_sumscore.rds")))
write_csv(results, file.path(out.dir, paste0("full_results_bt",n_bootstrap,"_sumscore.csv")))


# --- 7. SUMMARIZE RESULTS (Bootstrap Aggregation) --- #

# Summarize bootstrap results to get mean/median and 95% CIs
summary_df <- results %>% 
  mutate(contrast = paste0(Protein_Set,' - ',Group, " vs ", Control_Group)) %>%
  group_by(Group, Protein_Set, Observation_Years, Control_Group, contrast) %>% 
  dplyr::summarise(
    OR_mean = mean(OR, na.rm = TRUE),
    OR_lower = quantile(OR, 0.025, na.rm = TRUE),
    OR_upper = quantile(OR, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  # Factorize Control_Group, removing other controls
  mutate(Control_Group = factor(Control_Group, 
                                levels = c('HC', 'Dis1noDis2'),
                                labels = c('Healthy Controls', 'Dis1 No Progression')))

# Save summary
write_csv(summary_df, file.path(out.dir, paste0("summary_bt",n_bootstrap,"_sumscore.csv")))
saveRDS(summary_df, file.path(out.dir, paste0("summary_bt",n_bootstrap,"_sumscore.rds")))


# --- 8. CASE COUNT CALCULATION --- #
# Calculate number of cases included in each time interval for documentation

calculate_case_counts <- function(clinic_data, max_years = 14) {
  time_intervals <- c(2, 4, 6, 8, 10, 12, 14) 
  clinic_data %>%
    mutate(delta_diag_enroll = as.numeric(diagnosis_time - date_of_attending_assessment_centre_f53_0_0) / 365.25) %>%
    map_dfr(time_intervals, function(yrs) {
      filter(., delta_diag_enroll <= yrs) %>%
        count() %>%
        mutate(years = yrs, .before = 1)
    })
}

case_counts_inci <- calculate_case_counts(dis2_inci_clinic)
# Save case counts
write_csv(case_counts_inci, file.path(out.dir, "case_counts_dis2_incident.csv"))

# --- 9. VISUALIZATION (OR Plots) --- #

# Prepare data for plotting
plot_data <- summary_df %>% 
  filter(Control_Group != "Dis1 No Progression") # Focus on HC comparison for primary plots

# Create faceted plots for Odds Ratio (OR)
or_plot <- ggplot(plot_data,
                  aes(x = Observation_Years, y = log2(OR_mean), 
                      color = Protein_Set, group = Protein_Set)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log2(OR_lower), ymax = log2(OR_upper)), width = 0.2) +
  facet_grid(. ~ Group, scales = "fixed") + 
  labs(title = "Protein Sum Score Odds Ratios (vs HC)",
       subtitle = paste("Mean Log2(OR) (95% CI) for", dis2),
       x = "Observation Years Before/At Diagnosis", 
       y = expression(log[2]~"Odds Ratio"),
       color = "Protein Set") +
  scale_color_manual(values = c("FAP_list" = "#FC7F0E", "random" = '#7F7F7F',"CEACAM5" ="#2CA02C")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') + # Reference line at log2(OR)=0 (OR=1)
  scale_y_continuous(limits = c(-2.5, 6), breaks = c(-2, 0, 2, 4, 6)) +
  scale_x_continuous(limits = c(1.9, 14.1), breaks = observation_years) +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out.dir, paste0("OR_plot_bt",n_bootstrap,"_sumscore.pdf")), or_plot, 
       width = 7, height = 3.5, dpi = 300)