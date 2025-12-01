# ==============================================================================
# Analysis of GDF15 IHC (AOD) across FAP/CRC Tissue Types
# DESCRIPTION:
# This script loads data from an Excel file containing sample information for
# Familial Adenomatous Polyposis (FAP) and Colorectal Cancer (CRC) patients.
# It performs patient filtering, data preparation, linear mixed-effects modeling
# (LMM) to account for repeated measures, post-hoc analysis, and generates
# a boxplot visualization of AOD values by tissue type.
# ==============================================================================
# Load libraries####
library(Matrix)      
library(dplyr)        
library(readxl)      
library(ggplot2)     
library(tidyr)       
library(broom)        
library(lme4)         
library(lmerTest)    
library(emmeans)      

# Define output directory
output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 1. DATA LOADING AND PREPROCESSING ####
FAP <- read_excel("data/FAP and CRC.xlsx", sheet = "FAP")
CRC <- read_excel("data/FAP and CRC.xlsx", sheet = "CRC")

# Define Disease type for current analysis (Change as needed: 'FAP' or 'CRC')
Disease <- 'FAP' 
message(paste("Running analysis for:", Disease))

# Select the appropriate dataset based on the 'Disease' variable
if (Disease == 'FAP') {
  data <- FAP
} else if (Disease == 'CRC') {
  data <- CRC
} else {
  stop("Invalid disease type. Choose either 'FAP' or 'CRC'.")
}

# Data Filtering: Ensure each patient has samples from both primary groups####
print(table(data$Patient_ID, data$Tissue_type))
# Filter 1: Require at least one "Normal" or "Inflammation" sample (Control/Baseline)
data <- data %>%
  group_by(Patient_ID) %>%
  filter(any(Tissue_type %in% c("Normal", "Inflammation"))) %>%
  ungroup()

# Filter 2: Require at least one "Neoplasia" or "Adenocarcinoma" sample (Disease/Progression)
data <- data %>%
  group_by(Patient_ID) %>%
  filter(any(Tissue_type %in% c("Low-grade neoplasia", "High-grade neoplasia", "Adenocarcinoma"))) %>%
  ungroup()

# Create derived sample columns
# 'sample': Extracts the second-to-last character from the 'Label' column
data$sample <- substr(data$Label, nchar(data$Label) - 1, nchar(data$Label) - 1)
# 'pt_sample': Unique identifier for Patient-Sample combination
data$pt_sample <- paste(data$Patient_ID, data$sample, sep = "_")
print(table(data$Patient_ID, data$Tissue_type))

# Convert Tissue_type to a factor with specific, desired level order and labels
data$Tissue_type <- factor(
  data$Tissue_type,
  levels = c("Normal", "Inflammation", "Low-grade neoplasia", "High-grade neoplasia", "Adenocarcinoma"),
  labels = c("Normal", "Inflammation", "LGD", "HGD", "Adenocarcinoma")
)

# 2. STATISTICAL ANALYSIS - LINEAR MIXED-EFFECTS MODEL (LMM) ####
# Use LMM to analyze the effect of Tissue_type on AOD, while adjusting for the non-independence of samples taken from the same patient (repeated measures).
# (1 | Patient_ID) includes a random intercept for each patient.
model <- lmer(AOD ~ Tissue_type + (1 | Patient_ID), data = data)
print("Summary of Linear Mixed-Effects Model:")
summary(model)

# Extract fixed effects coefficients and p-values
fixed_effects <- as.data.frame(coef(summary(model)))
colnames(fixed_effects) <- c("Estimate", "Std_Error", "DF", "t_value", "p_value")

# Save the fixed effects results (compares each tissue type against the reference: 'Normal')
write.csv(
  fixed_effects, 
  file = file.path(output_dir, paste0(Disease, "_lmer_fixed_effects_results.csv")), 
  row.names = TRUE # Keep row names which are the model terms
)


# 3. POST-HOC ANALYSIS (PAIRWISE COMPARISONS) ####
# Post-hoc tests compare the estimated marginal means (EMMs) of each Tissue_type group.
# This explicitly tests differences between groups after accounting for the random effects.

# Calculate Estimated Marginal Means (EMMs)
emm <- emmeans(model, specs = ~ Tissue_type)

# Perform Pairwise Comparisons, referencing 'Normal'
posthoc <- pairs(
  emm,  
  adjust = "fdr",      # Adjust p-values using the False Discovery Rate (FDR/Benjamini-Hochberg) method
  ref = "Normal"       # Set 'Normal' tissue as the comparison reference group
)
print("Summary of Post-hoc (Normal vs Other) Comparisons:")
print(summary(posthoc, infer = TRUE))

# Save the post-hoc comparison results
write.csv(
  summary(posthoc, infer = TRUE),
  file = file.path(output_dir, paste0(Disease, "_posthoc_results.csv")),
  row.names = FALSE
)


# 4. DATA VISUALIZATION - BOX PLOT   #####
data$Tissue_type_rev <- factor(data$Tissue_type, 
                               levels = rev(levels(data$Tissue_type)))

pdf(file = file.path(output_dir, paste0(Disease, "_AOD_by_Tissue_Type.pdf")), 
    width = 5, height = 10)
boxplot(
  AOD ~ Tissue_type_rev,
  data = data,
  horizontal = TRUE,      # Orient the plot horizontally
  las = 1,                # Rotate y-axis labels to be readable
  main = paste("AOD Signal by Tissue Type for", Disease),
  xlab = "AOD Value",
  ylab = "Tissue Type (Progression)",
  col = "lightblue",      # Box fill color
  border = "darkblue",
  outline = FALSE         # Hide standard outliers to use custom stripchart points
)
stripchart(
  AOD ~ Tissue_type_rev,
  data = data,
  method = "jitter",
  vertical = FALSE,       # Must be FALSE for horizontal boxplot
  pch = 16,               # Solid circle points
  cex = 0.8,              # Point size
  col = rgb(0, 0, 0, 0.5), # Semi-transparent black color
  add = TRUE              # Add to the existing plot
)

# Close the PDF device to save the file
dev.off()

