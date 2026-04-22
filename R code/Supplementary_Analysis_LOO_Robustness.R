############################################################
## Supplementary Analyses - LOO Ranking
############################################################

############################################################
# 1. LIBRARIES
############################################################

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(Seurat)
library(tidyverse)
library(zellkonverter)
library(readr)
library(irr)  

############################################################
# 2. UTILITY FUNCTIONS
############################################################

safe_dir <- function(path){
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

save_object <- function(obj, path){
  saveRDS(obj, path)
}

############################################################
# 3. GENE RANKING FUNCTION
############################################################

gene_ranking <- function(sender_receiver_de,
                         prioritization_tables,
                         lr_target_df,
                         cell_abundance_info_tbl,
                         num_cells){
  
  rank_lfc_ligand <- sender_receiver_de %>%
    select(contrast, sender, ligand, lfc_ligand) %>%
    distinct() %>%
    arrange(desc(abs(lfc_ligand)))
  
  rank_lfc_avg <- sender_receiver_de %>%
    select(contrast, sender, ligand, receptor, ligand_receptor_lfc_avg) %>%
    distinct() %>%
    arrange(desc(abs(ligand_receptor_lfc_avg)))
  
  rank_prioritization_score <- prioritization_tables$group_prioritization_table_source %>%
    select(contrast, sender, ligand, prioritization_score) %>%
    distinct() %>%
    arrange(desc(abs(prioritization_score)))
  
  top_target_genes <- lr_target_df %>%
    group_by(group, sender, ligand, receiver) %>%
    summarise(total_vals = n(), .groups = "drop")
  
  top_ligand_target <- lr_target_df %>%
    group_by(sender, ligand, target) %>%
    summarise(total = n(), .groups = "drop") %>%
    arrange(desc(total))
  
  weighted_df <- lr_target_df %>%
    group_by(group, sender, ligand, receiver) %>%
    summarise(total = n(), .groups = "drop") %>%
    left_join(cell_abundance_info_tbl,
              by = c("receiver", "group")) %>%
    mutate(cumulative_effect = total * (Freq / num_cells)) %>%
    arrange(desc(cumulative_effect))
  
  total_score_per_ligand <- weighted_df %>%
    group_by(sender.x, ligand) %>%
    summarise(total_cumulative_per_ligand = sum(cumulative_effect, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(ligand) %>%
    summarise(total_score = sum(total_cumulative_per_ligand, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(total_score))
  
  list(
    rank_lfc_ligand          = rank_lfc_ligand,
    rank_lfc_avg             = rank_lfc_avg,
    rank_prioritization_score = rank_prioritization_score,
    top_target_genes         = top_target_genes,
    top_ligand_target        = top_ligand_target,
    weighted_df              = weighted_df,
    total_score_per_ligand   = total_score_per_ligand
  )
}

############################################################
# 4. LIGAND RANKING ROBUSTNESS FUNCTION
############################################################

ligand_ranking_robustness <- function(ligand_ranking_df, dir.list, top_n = 40, max_k = 40){
  
  # build ranking table from specified columns
  ligand_ranking_df <- ligand_ranking_df %>%
    select(all_of(dir.list)) %>%
    head(top_n)
  
  # convert to long format
  df_long <- ligand_ranking_df %>%
    mutate(rank = row_number()) %>%
    pivot_longer(-rank, names_to = "run", values_to = "ligand")
  
  # ligand x run rank matrix
  rank_table <- df_long %>%
    pivot_wider(names_from = run, values_from = rank)
  
  rank_matrix <- as.data.frame(rank_table)
  rownames(rank_matrix) <- rank_matrix$ligand
  rank_matrix$ligand <- NULL
  
  # Spearman correlation
  cor_matrix    <- cor(rank_matrix, method = "spearman", use = "pairwise.complete.obs")
  mean_spearman <- mean(cor_matrix[upper.tri(cor_matrix)])
  
  # rank variability
  rank_sd <- apply(rank_matrix, 1, sd, na.rm = TRUE)
  
  # Kendall W
  kendall_res <- kendall(rank_matrix)
  
  # ---- Top-k Jaccard ----
  topk_jaccard <- function(df, k){
    runs <- colnames(df)
    n    <- length(runs)
    jac  <- matrix(NA, n, n)
    for(i in 1:n){
      for(j in 1:n){
        A <- df[[runs[i]]][1:k]
        B <- df[[runs[j]]][1:k]
        jac[i, j] <- length(intersect(A, B)) / length(union(A, B))
      }
    }
    mean(jac[upper.tri(jac)])
  }
  
  k_vals        <- 1:max_k
  jaccard_curve <- sapply(k_vals, function(k) topk_jaccard(ligand_ranking_df, k))
  
  list(
    ligand_ranking_df = ligand_ranking_df,
    rank_matrix       = rank_matrix,
    spearman_matrix   = cor_matrix,
    mean_spearman     = mean_spearman,
    rank_sd           = rank_sd,
    kendall           = kendall_res,
    k_vals            = k_vals,
    jaccard_curve     = jaccard_curve
  )
}

############################################################
# 5. MAIN LIGAND PIPELINE (per sample)
############################################################

run_multinichenet_pipeline <- function(disease,
                                       input_dir,
                                       output_dir){
  
  cat("Running ligand analysis for:", disease, "\n")
  safe_dir(output_dir)
  
  # ---- Load data ----
  DE_output_contrastlist_2 <- read_csv(paste0(input_dir, "DE_output_contrastlist_2.csv"),
                                       show_col_types = FALSE)
  
  # ---- Extract objects ----
  prioritization_tables   <- readRDS(paste0(input_dir, "prioritization_tables_1.rds"))
  cell_abundance_info_tbl <- readRDS(paste0(input_dir, "cell_abundance_info_tbl_1.rds"))
  
  lr_target_df <- readRDS(paste0(input_dir, "lr_target_df_1.rds")) %>%
    left_join(
      prioritization_tables$group_prioritization_tbl %>%
        select(id, prioritization_score),
      by = "id"
    )
  
  # ---- Ligand target barplot ----
  ligand_counts <- lr_target_df %>%
    group_by(ligand) %>%
    summarise(n = n(), .groups = "drop")
  
  p1 <- ggplot(ligand_counts %>% head(40),
               aes(x = reorder(ligand, -n), y = n)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(paste0(output_dir, "ligand_target_counts.png"), p1, width = 12, height = 6)
  
  # ---- Cell abundance processing ----
  cell_abundance_info_tbl <- cell_abundance_info_tbl %>%
    filter(group == "FAP_Polyp") %>%
    mutate(sender = receiver,
           prop   = Freq / sum(Freq))
  
  num_cells <- sum(cell_abundance_info_tbl$Freq)
  
  # ---- Prioritization score ----
  prioritization_score_ligand <- lr_target_df %>%
    group_by(group, ligand, receiver) %>%
    summarise(
      total                = n(),
      prioritization_score = mean(prioritization_score, na.rm = TRUE),
      .groups              = "drop"
    ) %>%
    left_join(cell_abundance_info_tbl, by = c("receiver", "group")) %>%
    mutate(cumulative_effect = total * prop) %>%
    arrange(desc(cumulative_effect)) %>%
    group_by(ligand, group) %>%
    mutate(
      total_cumulative_per_ligand = sum(cumulative_effect, na.rm = TRUE),
      mean_prior_score            = mean(prioritization_score),
      max_prior_score             = max(prioritization_score),
      sum_prior_score             = sum(prioritization_score)
    ) %>%
    distinct()
  
  save_object(prioritization_score_ligand,
              paste0(output_dir, "prioritization_score_ligand_", disease, ".rds"))
  
  # ---- Prioritization plot ----
  p2 <- ggplot(prioritization_score_ligand %>% head(40),
               aes(x = total_cumulative_per_ligand, y = max_prior_score)) +
    geom_point(size = 2, color = "#0073C2") +
    geom_text(aes(label = ligand), size = 3) +
    theme_classic()
  
  ggsave(paste0(output_dir, "prioritization_plot_", disease, ".png"), p2,
         width = 12, height = 8)
  
  # ---- Differential expression ----
  DE_FAP_Polyp_vs_Healthy     <- DE_output_contrastlist_2 %>%
    filter(contrast == "FAP_Polyp-Healthy_Normal")
  
  DE_FAP_Unaffected_vs_Healthy <- DE_output_contrastlist_2 %>%
    filter(contrast == "FAP_Unaffected-Healthy_Normal")
  
  DE_wide <- DE_FAP_Unaffected_vs_Healthy %>%
    select(gene, cluster_id, logFC, p_adj) %>%
    rename(logFC_FAP_Unaffected = logFC,
           p_adj_FAP_Unaffected  = p_adj) %>%
    inner_join(
      DE_FAP_Polyp_vs_Healthy %>%
        select(gene, cluster_id, logFC, p_adj) %>%
        rename(logFC_FAP_Polyp = logFC,
               p_adj_FAP_Polyp  = p_adj),
      by = c("gene", "cluster_id")
    ) %>%
    mutate(DEG_FLAG = (logFC_FAP_Unaffected * logFC_FAP_Polyp) > 0) %>%
    filter(DEG_FLAG)
  
  save_object(DE_wide, paste0(output_dir, "DE_wide_", disease, ".rds"))
  
  # ---- Target comparison ----
  celltype_de  <- readRDS(paste0(input_dir, "celltype_de_2.rds"))
  lr_target_df_2 <- readRDS(paste0(input_dir, "lr_target_df_2.rds"))
  
  gene_list <- celltype_de %>%
    filter(contrast == "FAP_Unaffected-Healthy_Normal") %>%
    select(gene) %>%
    distinct()
  
  lrtdf <- lr_target_df_2 %>% filter(group == "FAP_Unaffected")
  
  ligand_target_counts_FAP_Unaffected <- gene_list %>%
    rowwise() %>%
    mutate(n_targets = sum(lrtdf$ligand == gene)) %>%
    ungroup()
  
  ligand_target_counts_FAP_Polyp <- lr_target_df_2 %>%
    filter(group == "FAP_Polyp") %>%
    group_by(ligand) %>%
    summarise(num_targets = n(), .groups = "drop")
  
  ligand_comparison <- ligand_target_counts_FAP_Polyp %>%
    left_join(ligand_target_counts_FAP_Unaffected,
              by = c("ligand" = "gene")) %>%
    mutate(inc_target_vs_prev_stage = num_targets >= n_targets)
  
  filtered_ligand_comparison <- ligand_comparison %>%
    filter(inc_target_vs_prev_stage)
  
  save_object(filtered_ligand_comparison,
              paste0(output_dir, "filtered_ligands_", disease, ".rds"))
  
  # ---- Gene ranking ----
  filtered_df <- readRDS(paste0(input_dir, "lr_target_df_1.rds")) %>%
    filter(ligand %in% filtered_ligand_comparison$ligand)
  
  sender_receiver_de_1 <- readRDS(paste0(input_dir, "sender_receiver_de_1.rds"))
  
  ranking_outputs <- gene_ranking(
    sender_receiver_de   = sender_receiver_de_1,
    prioritization_tables = prioritization_tables,
    lr_target_df         = filtered_df,
    cell_abundance_info_tbl = cell_abundance_info_tbl,
    num_cells            = num_cells
  )
  
  for(name in names(ranking_outputs)){
    save_object(ranking_outputs[[name]],
                paste0(output_dir, name, "_", disease, ".rds"))
  }
  
  cat("Finished:", disease, "\n")
  
  return(list(ligand_ranking = ranking_outputs$total_score_per_ligand))
}

############################################################
# 6. TARGET RANKING FUNCTION
############################################################

target_ranking <- function(input_dir, ligand_list){
  
  lr_target_df_1         <- read_csv(paste0(input_dir, "lr_target_df_1.csv"),
                                     show_col_types = FALSE)
  DE_output_contrastlist_1 <- read_csv(paste0(input_dir, "DE_output_contrastlist_1.csv"),
                                       show_col_types = FALSE)
  
  df <- lr_target_df_1 %>%
    left_join(DE_output_contrastlist_1,
              by = c("target" = "gene", "receiver" = "cluster_id"))
  
  ligand_counts <- df %>%
    filter(ligand %in% ligand_list) %>%
    group_by(target) %>%
    summarise(num_ligands = n_distinct(ligand), .groups = "drop")
  
  df_wide <- df %>%
    select(target, ligand, logFC) %>%
    distinct() %>%
    pivot_wider(
      names_from  = ligand,
      values_from = logFC,
      values_fn   = sum
    ) %>%
    select(target, all_of(ligand_list))
  
  score_df <- ligand_counts %>%
    left_join(df_wide, by = "target") %>%
    mutate(TE = rowSums(across(all_of(ligand_list)), na.rm = TRUE))
  
  ranked_target_score_df <- score_df %>%
    arrange(desc(abs(TE))) %>%
    select(target, num_ligands, TE)
  
  save_dir  <- paste0(input_dir, "Top_", length(ligand_list), "/")
  safe_dir(save_dir)
  save_object(ranked_target_score_df, paste0(save_dir, "ranked_target_score_df.rds"))
}

############################################################
# 7. RUN — LIGAND PIPELINE (LOO loop)
############################################################

output_dir <- "~/DT/FAP v2 Resubmission/"
dir.list   <- c("A001", "A002", "A008", "A010", "A014", "A015", "A018", "A022", "F")

# dir.list <- c("A001","A002")

result_df <- list()

for(subfold in dir.list){
  
  output_dir_loop        <- paste0(output_dir, subfold, "/")
  result_df[[subfold]]   <- run_multinichenet_pipeline(
    disease    = "FAP",
    input_dir  = output_dir_loop,
    output_dir = output_dir_loop
  )
}

saveRDS(result_df, paste0(output_dir, "LOO_ligand_ranking_list.rds"))

ligand_lists <- lapply(result_df, \(x) x$ligand_ranking$ligand)

max_len      <- max(lengths(ligand_lists))
ligand_lists <- lapply(ligand_lists, function(x){ length(x) <- max_len; x })

ligand_ranking_df <- bind_cols(ligand_lists)

saveRDS(ligand_ranking_df, paste0(output_dir, "LOO_ligand_ranking_df.rds"))

# ---- Ligand robustness analysis ----

ligand_ranking_robustness_results <- ligand_ranking_robustness(
  ligand_ranking_df,
  dir.list = dir.list,
  top_n    = 40,
  max_k    = 40
)

saveRDS(ligand_ranking_robustness_results,
        paste0(output_dir, "ligand_ranking_robustness_results.rds"))


############################################################
# 8. RUN — TARGET PIPELINE (LOO loop)
############################################################

for(i in dir.list){
  
  cat("Processing target ranking for:", i, "\n")
  
  master_ligand_list <- readRDS(paste0(output_dir, i, "/total_score_per_ligand_FAP.rds"))
  
  for(n in seq(5, 40, 5)){
    ligand_list <- master_ligand_list %>%
      arrange(desc(total_score)) %>%
      slice_head(n = n) %>%
      pull(ligand)
    
    cat("  Top", n, "ligands — selecting targets\n")
    target_ranking(paste0(output_dir, i, "/"), ligand_list)
  }
}

# ---- Assemble target ranking dataframe (top-10 ligands) ----

target_list <- list()

for(i in dir.list){
  ranked_target_score_df <- readRDS(paste0(output_dir, i, "/Top_10/ranked_target_score_df.rds"))
  target_list[[i]]       <- ranked_target_score_df %>%
    arrange(desc(abs(TE))) %>%
    slice_head(n = 10)
}

saveRDS(target_list, paste0(output_dir, "LOO_target_ranking_list.rds"))

target_lists <- lapply(target_list, \(x) x$target)

max_len      <- max(lengths(target_lists))
target_lists <- lapply(target_lists, function(x){ length(x) <- max_len; x })

target_ranking_df <- bind_cols(target_lists)

saveRDS(target_ranking_df, paste0(output_dir, "LOO_target_ranking_df.rds"))

# ---- Target robustness analysis ----

target_ranking_robustness_results <- ligand_ranking_robustness(
  target_ranking_df,
  dir.list = dir.list,
  top_n    = 10,
  max_k    = 10
)

saveRDS(target_ranking_robustness_results,
        paste0(output_dir, "LOO_target_ranking_robustness_results.rds"))

############################################################
# 9. COMBINED LIGAND + TARGET RANKING & VISUALISATION
############################################################

# NOTE: output_dir and dir.list are inherited from sections 7/8.

# ---- Load and trim LOO ranking outputs ----
ligand_ranking_df <- readRDS(paste0(output_dir, "LOO_ligand_ranking_df.rds")) %>%
  select(-A022) #A022 is HC

target_ranking_df <- readRDS(paste0(output_dir, "LOO_target_ranking_df.rds")) %>%
  select(-A022) #A022 is HC

# ---- Build combined ranking dataframe ----
# Top 10 ligands stacked above top 10 targets, then deduplicated per column
combined_ranking_df <- rbind(ligand_ranking_df[1:10, ], target_ranking_df)

combined_ranking_df <- as.data.frame(
  lapply(combined_ranking_df, function(col){
    unique_vals          <- unique(col)
    length(unique_vals)  <- max(lengths(lapply(combined_ranking_df, unique)))
    unique_vals
  })
)

# ---- Append FdS reference column ----
full_Fds <- c(
  "LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F",
  "INHBA", "LAMA5", "COL7A1", "SOX9", "CCND1", "TGFBI", "MET",
  "ID1", "PTMA", "CLU", "TGIF1", "EDN1", "GDF15"
) %>% unique()  # 19 unique genes after dedup of TGFBI

combined_ranking_df <- cbind(combined_ranking_df, "FdS" = full_Fds)

# ---- Build gene × run rank matrix ----
df_long <- combined_ranking_df %>%
  mutate(rank = row_number()) %>%
  pivot_longer(-rank, names_to = "run", values_to = "gene")

rank_table <- df_long %>%
  pivot_wider(names_from = run, values_from = rank)

rank_matrix <- as.data.frame(rank_table)
rownames(rank_matrix) <- rank_matrix$gene
rank_matrix$gene      <- NULL

# ---- Spearman correlation heatmap ----
cor_matrix <- cor(rank_matrix, method = "spearman", use = "pairwise.complete.obs")

pheatmap(
  cor_matrix,
  main = "Spearman correlation between ranking runs"
)

# ---- Top-k Jaccard (combined) ----
# NOTE: na.omit applied per run so missing genes don't inflate union
topk_jaccard_combined <- function(df, k){
  runs <- colnames(df)
  n    <- length(runs)
  jac  <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      A <- na.omit(df[[runs[i]]][1:k])
      B <- na.omit(df[[runs[j]]][1:k])
      jac[i, j] <- length(intersect(A, B)) / length(union(A, B))
    }
  }
  mean(jac[upper.tri(jac)])
}

k_vals        <- 1:nrow(combined_ranking_df)
jaccard_curve <- sapply(k_vals, function(k) topk_jaccard_combined(combined_ranking_df, k))

plot(
  k_vals, jaccard_curve,
  type = "l", lwd = 3,
  xlab = "Top-k genes",
  ylab = "Mean Jaccard similarity",
  main = "Top-k ranking robustness (combined)"
)

# ---- Top-k overlap matrix ----
topk_overlap_matrix <- function(df, k){
  runs <- colnames(df)
  n    <- length(runs)
  mat  <- matrix(NA, n, n, dimnames = list(runs, runs))
  for(i in 1:n){
    for(j in 1:n){
      A <- df[[runs[i]]][1:k]
      B <- df[[runs[j]]][1:k]
      mat[i, j] <- length(intersect(A, B)) / k
    }
  }
  mat
}

overlap_mat <- topk_overlap_matrix(combined_ranking_df, k = 19)

pheatmap(
  overlap_mat,
  main = "Top-19 gene overlap across runs"
)

# ---- Kendall W (raw and NA-filled) ----
kendall(rank_matrix)

rank_matrix_filled                          <- rank_matrix
rank_matrix_filled[is.na(rank_matrix_filled)] <- max(rank_matrix, na.rm = TRUE) + 1

kendall(rank_matrix_filled)

# ---- Rank stability heatmap ----
pheatmap(
  rank_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main         = "Gene Rank Stability Across LOO Runs"
)

# ---- Rank trajectory plot (missing as X) ----
rank_plot_df <- rank_matrix %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "rank") %>%
  mutate(
    rank_plot = ifelse(is.na(rank), max(rank, na.rm = TRUE) + 1, rank),
    missing   = is.na(rank)
  )

# ORPHANED — initial trajectory plot without missing-value handling; superseded below
# ggplot(rank_plot_df, aes(x = sample, y = rank, group = gene, color = gene)) +
#   geom_line(alpha = 0.6) +
#   geom_point(size = 3) +
#   scale_y_reverse() + theme_minimal() + ...

ggplot(rank_plot_df, aes(sample, rank_plot, group = gene, color = gene)) +
  geom_line(alpha = 0.6) +
  geom_point(data = subset(rank_plot_df, !missing), size = 3) +
  geom_point(
    data   = subset(rank_plot_df, missing),
    shape  = 4,     # X symbol for genes absent in that LOO run
    size   = 3,
    stroke = 1.5
  ) +
  scale_y_reverse() +
  theme_minimal() +
  labs(
    title = "FdS Gene Rank Variation across Leave One Out Runs",
    x     = "Left Out Sample",
    y     = "Rank"
  )
