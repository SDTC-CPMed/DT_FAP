rm(list = ls())

# =========================
# Packages
# =========================
packages <- c(
  "Seurat", "tidyverse", "tibble", "dplyr", "stringr", "ggplot2",
  "patchwork", "lme4", "lmerTest", "emmeans", "broom.mixed", "forcats"
)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}
install_if_missing(packages)

library(Seurat)
library(tidyverse)
library(tibble)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(forcats)

# =========================
# Input / output
# =========================
inputdir <- "hest_data/"
outdir   <- "output_FdS_mixedmodel/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

samples_of_interest <- c(
  "MISC35","MISC36","MISC43","MISC48","MISC51",
  "MISC62","MISC70","MISC71","MISC72"
)

goi <- c(
  "LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F",
  "INHBA", "LAMA5", "COL7A1", "SOX9", "CCND1", "MET", "ID1",
  "PTMA", "CLU", "TGIF1", "EDN1", "GDF15"
)

# use the stages actually included in current analysis
label_levels <- c("Normal", "LGD", "CRC")

# Name for the signature feature
signature_name <- "FDS_score"

# =========================
# Load processed Seurat list
# =========================
load("output/seurat_list_processed_Epi_subset.Rdata")  # contains seurat_list_processed

# =========================
# 0) Build combined spot-level dataframe####
# =========================
spot_data_list <- list()

for (sample in names(seurat_list_processed)) {
  if (!sample %in% samples_of_interest) next
  
  se <- seurat_list_processed[[sample]]
  
  se$Label <- factor(as.character(se$Label), levels = label_levels)
  
  # keep only genes present
  genes_present <- intersect(goi, rownames(se))
  if (length(genes_present) < 2) {
    message("Skipping sample ", sample, ": too few genes present.")
    next
  }
  
  expr_df <- FetchData(se, vars = genes_present) %>%
    rownames_to_column("Spot")
  
  meta_df <- se@meta.data %>%
    rownames_to_column("Spot") %>%
    mutate(Sample = sample) %>%
    dplyr::select(Spot, Sample, Label, nCount_Spatial, nFeature_Spatial, everything())
  
  tmp <- meta_df %>%
    left_join(expr_df, by = "Spot") %>%
    filter(!is.na(Label))
  
  spot_data_list[[sample]] <- tmp
}

spot_df <- bind_rows(spot_data_list)

# make sure all requested genes exist as columns; if missing, add NA columns
missing_gene_cols <- setdiff(goi, colnames(spot_df))
if (length(missing_gene_cols) > 0) {
  for (g in missing_gene_cols) {
    spot_df[[g]] <- NA_real_
  }
}

# Keep only spots with valid labels
spot_df <- spot_df %>%
  filter(Label %in% label_levels) %>%
  filter(Sample %in% samples_of_interest) %>%
  mutate(
    Label = factor(Label, levels = label_levels),
    log_nCount = log10(nCount_Spatial + 1)
  )

write.csv(
  spot_df[, c("Spot", "Sample", "Label", "nCount_Spatial", "nFeature_Spatial")],
  file.path(outdir, "spot_metadata_used_for_model.csv"),
  row.names = FALSE
)

# =========================
# Add signature score (PC1)
# PCA across spots using GOI expression
# compute the signature score using the first principal component value derived from Principal Component Analysis of the expression levels of the input genes
# =========================
genes_for_signature <- goi[colSums(!is.na(spot_df[, goi, drop = FALSE])) > 0]

if (length(genes_for_signature) < 2) {
  stop("Fewer than 2 genes available to compute PCA signature score.")
}

# impute missing gene expression with gene-wise median for PCA only
expr_for_pca <- spot_df[, genes_for_signature, drop = FALSE]

for (g in genes_for_signature) {
  med_g <- median(expr_for_pca[[g]], na.rm = TRUE)
  if (is.na(med_g)) med_g <- 0
  expr_for_pca[[g]][is.na(expr_for_pca[[g]])] <- med_g
}

# PCA expects rows = spots, cols = genes
pca_res <- prcomp(expr_for_pca, center = TRUE, scale. = TRUE)

spot_df[[signature_name]] <- pca_res$x[, 1]

# orient PC1 so that CRC has higher score than Normal
mean_crc <- mean(spot_df[[signature_name]][spot_df$Label == "CRC"], na.rm = TRUE)
mean_normal <- mean(spot_df[[signature_name]][spot_df$Label == "Normal"], na.rm = TRUE)

if (!is.na(mean_crc) && !is.na(mean_normal) && mean_crc < mean_normal) {
  spot_df[[signature_name]] <- -spot_df[[signature_name]]
}

write.csv(
  data.frame(
    Gene = genes_for_signature,
    PC1_loading = pca_res$rotation[genes_for_signature, 1]
  ),
  file.path(outdir, "FDS_signature_PC1_loadings.csv"),
  row.names = FALSE
)

# =========================
# 1) Summarize stage composition per sample####
# =========================
stage_count_table <- spot_df %>%
  dplyr::count(Sample, Label, name = "n_spots") %>%
  tidyr::pivot_wider(names_from = Label, values_from = n_spots, values_fill = 0)

stage_presence_table <- spot_df %>%
  dplyr::distinct(Sample, Label) %>%
  dplyr::count(Sample, name = "n_stages_present") %>%
  dplyr::arrange(desc(n_stages_present), Sample)

stage_summary <- stage_count_table %>%
  dplyr::left_join(stage_presence_table, by = "Sample") %>%
  dplyr::arrange(desc(n_stages_present), Sample)

write.csv(
  stage_summary,
  file.path(outdir, "sample_stage_composition.csv"),
  row.names = FALSE
)

# Samples with >=2 stages
eligible_samples <- stage_presence_table %>%
  filter(n_stages_present >= 2) %>%
  pull(Sample)

spot_df_model <- spot_df %>%
  filter(Sample %in% eligible_samples)

cat("Number of samples with >=2 stages:", length(unique(spot_df_model$Sample)), "\n")

# =========================
# Helper: sample-stratified permutation test
# =========================
perm_test_within_sample <- function(df_gene, n_perm = 1000, seed = 123) {
  set.seed(seed)
  
  fit_full <- lm(expr ~ Label + log_nCount + Sample, data = df_gene)
  fit_red  <- lm(expr ~ log_nCount + Sample, data = df_gene)
  obs_tab  <- anova(fit_red, fit_full)
  obs_stat <- obs_tab$F[2]
  
  perm_stats <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    df_perm <- df_gene %>%
      group_by(Sample) %>%
      mutate(Label_perm = sample(Label, size = n(), replace = FALSE)) %>%
      ungroup()
    
    fit_full_perm <- lm(expr ~ Label_perm + log_nCount + Sample, data = df_perm)
    fit_red_perm  <- lm(expr ~ log_nCount + Sample, data = df_perm)
    perm_tab      <- anova(fit_red_perm, fit_full_perm)
    perm_stats[b] <- perm_tab$F[2]
  }
  
  emp_p <- (sum(perm_stats >= obs_stat, na.rm = TRUE) + 1) / (n_perm + 1)
  
  tibble(
    perm_F_obs = obs_stat,
    perm_p_empirical = emp_p
  )
}

# =========================
# 2) Mixed model signature####
# =========================
features_to_model <- c( signature_name)

overall_results <- list()
contrast_results <- list()

for (feature in features_to_model) {
  message("Processing feature: ", feature)
  
  if (!feature %in% colnames(spot_df_model)) {
    message("Skipping missing feature: ", feature)
    next
  }
  
  df_gene <- spot_df_model %>%
    dplyr::select(Spot, Sample, Label, log_nCount, expr = all_of(feature)) %>%
    dplyr::filter(!is.na(expr), !is.na(Label), !is.na(log_nCount))
  
  if (nrow(df_gene) < 20 || length(unique(df_gene$Label)) < 2 || sd(df_gene$expr) == 0) {
    message("Skipping low-information feature: ", feature)
    next
  }
  
  fit_full <- lmer(expr ~ Label + log_nCount + (1 | Sample), data = df_gene, REML = FALSE)
  fit_red  <- lmer(expr ~ log_nCount + (1 | Sample), data = df_gene, REML = FALSE)
  
  lr <- anova(fit_red, fit_full)
  
  overall_res <- tibble(
    Gene = feature,
    n_spots = nrow(df_gene),
    n_samples = n_distinct(df_gene$Sample),
    n_stages = n_distinct(df_gene$Label),
    LRT = lr$Chisq[2],
    df = lr$`Chi Df`[2],
    P_value = lr$`Pr(>Chisq)`[2]
  )
  
  emm <- emmeans(fit_full, ~ Label)
  contr <- contrast(emm, method = "trt.vs.ctrl", ref = "Normal") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(Gene = feature)
  
  perm_res <- perm_test_within_sample(
    df_gene,
    n_perm = 1000,
    seed = 100 + which(features_to_model == feature)
  )
  
  overall_results[[feature]] <- overall_res %>% bind_cols(perm_res)
  contrast_results[[feature]] <- contr
}

overall_df <- bind_rows(overall_results) %>%
  mutate(
    P_adj_BH = p.adjust(P_value, method = "BH"),
    perm_p_adj_BH = p.adjust(perm_p_empirical, method = "BH")
  ) %>%
  arrange(P_adj_BH, perm_p_adj_BH)

contrast_df <- bind_rows(contrast_results) %>%
  dplyr::rename(
    Contrast   = contrast,
    EffectSize = estimate,
    CI_lower   = asymp.LCL,
    CI_upper   = asymp.UCL,
    P_value    = p.value
  ) %>%
  dplyr::group_by(Contrast) %>%
  dplyr::mutate(P_adj_BH = p.adjust(P_value, method = "BH")) %>%
  dplyr::ungroup()

write.csv(
  overall_df,
  file.path(outdir, "mixed_model_overall_stage_effects.csv"),
  row.names = FALSE
)

write.csv(
  contrast_df,
  file.path(outdir, "mixed_model_stage_contrasts_vs_normal.csv"),
  row.names = FALSE
)

# separate signature outputs
write.csv(
  overall_df %>% filter(Gene == signature_name),
  file.path(outdir, "signature_score_overall_stage_effects.csv"),
  row.names = FALSE
)

write.csv(
  contrast_df %>% filter(Gene == signature_name),
  file.path(outdir, "signature_score_stage_contrasts_vs_normal.csv"),
  row.names = FALSE
)

# =========================
# 3) Sample-level descriptive summaries####
# now includes signature score
# =========================
plot_df <- spot_df_model %>%
  dplyr::select(Sample, Label, all_of(features_to_model)) %>%
  pivot_longer(cols = all_of(features_to_model), names_to = "Gene", values_to = "expr") %>%
  group_by(Sample, Label, Gene) %>%
  summarise(
    median_expr = median(expr, na.rm = TRUE),
    mean_expr = mean(expr, na.rm = TRUE),
    n_spots = n(),
    .groups = "drop"
  )

write.csv(
  plot_df,
  file.path(outdir, "sample_stage_feature_summary.csv"),
  row.names = FALSE
)

 

# =========================
# 4) Signature summary plot ####
# =========================
signature_plot_df <- plot_df %>%
  filter(Gene == signature_name)

p_signature <- ggplot(signature_plot_df, aes(x = Label, y = median_expr)) +
  geom_violin(aes(fill = Label), color = "grey30", alpha = 0.7, trim = FALSE) +
  scale_fill_manual(values = c("#5268B1", "orange", "darkred")) +
  # geom_line(aes(group = Sample), color = "black", alpha = 0.25) +
  geom_point(color = "black", size = 1.8, alpha = 0.7) +
  coord_cartesian(ylim = c(-4.0, 5.0)) +
  theme_bw(base_size = 12) +
  labs(
    title = "FDS signature score across stages",
    x = "Histopathological stage",
    y = "FdS signature score",
    fill = "Stage"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

pdf(file.path(outdir, "signature_score_per_sample_stage_plot.pdf"), width = 4, height = 3)
print(p_signature)
dev.off()

# =========================
# 5) Overall summary table for manuscript use ####
# includes signature score
# =========================
manuscript_summary <- overall_df %>%
  dplyr::left_join(
    contrast_df %>%
      filter(Contrast == "CRC - Normal") %>%
      dplyr::select(
        Gene,
        EffectSize_crc_vs_normal = EffectSize,
        CI_lower_crc_vs_normal = CI_lower,
        CI_upper_crc_vs_normal = CI_upper,
        P_adj_crc_vs_normal = P_adj_BH
      ),
    by = "Gene"
  ) %>%
  dplyr::left_join(
    contrast_df %>%
      filter(Contrast == "LGD - Normal") %>%
      dplyr::select(
        Gene,
        EffectSize_lgd_vs_normal = EffectSize,
        CI_lower_lgd_vs_normal = CI_lower,
        CI_upper_lgd_vs_normal = CI_upper,
        P_adj_lgd_vs_normal = P_adj_BH
      ),
    by = "Gene"
  ) %>%
  dplyr::arrange(P_adj_BH, perm_p_adj_BH)

write.csv(
  manuscript_summary,
  file.path(outdir, "manuscript_summary_mixed_model_with_signature.csv"),
  row.names = FALSE
)

message("All done.")

