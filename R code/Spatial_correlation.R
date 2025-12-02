# ==============================================================================
# PROJECT: Spatial Transcriptomics (ST) Analysis - Tissue Heterogeneity & Gene Correlation
# AUTHOR: [Yelin Zhao]
# DESCRIPTION:
# This script processes and analyzes annotated Seurat objects from Spatial
# Transcriptomics experiments (Colon/H&E data). It performs Quality Control (QC),
# dimensionality reduction (PCA/UMAP), and calculates the Spearman correlation
# between selected Genes of Interest (GOI) and tissue pathology labels (Label)
# across multiple patient samples.
# ==============================================================================

#  0. INITIAL SETUP #####
rm(list = ls())
library(Seurat)
library(tidyverse)   # Includes ggplot2, dplyr, tibble, tidyr, stringr
library(tibble)      # Data frame utilities
library(dplyr)       # Data manipulation (part of tidyverse)
library(stringr)     # String manipulation (part of tidyverse)
library(patchwork)   # For combining ggplot plots (e.g., plot1 | plot2)
library(ggplot2)     # Data visualization (part of tidyverse)
library(scater)      # Single-cell analysis tools (often used for QC)
library(reticulate)  # Interface to Python
library(cowplot)     # Additional plot utilities (Used for combining plots)

inputdir <- paste0('data/ST') 
outdir <- paste0('output/')
if (dir.exists(outdir) == FALSE) {
  dir.create(outdir)
}

# Define sample lists and genes of interest (GOI)
samples_of_interest <- c("MISC35","MISC36","MISC43","MISC48","MISC51","MISC62","MISC70","MISC71","MISC72")
disease_status <- c("Cancer","Cancer","Cancer","Cancer","Cancer","Cancer","Cancer","Cancer","Cancer","Cancer")
goi <- c("LAMB1", "LAMC2", "TGFBI", "SEMA4D", "AGRN", "SEMA5A", "SEMA3F", "INHBA", 
         "LAMA5", "COL7A1", "SOX9", "CCND1", "TGFBI", "MET", "ID1", "PTMA", 
         "CLU", "TGIF1", "EDN1", "GDF15")


# 1. QC AND INITIAL SPATIAL PLOTTING (PRE-PROCESSING) ####
seurat_list <- list()
for (i in 1:length(samples_of_interest)) {
  sample <- samples_of_interest[i]
  print(paste("Processing initial QC for sample:", sample))
  load(paste0(inputdir, "/", sample, "_annotated21.Rdata"))
  
  # Check metadata and tissue labels
  print(head(se@meta.data))
  print(table(se$Label, useNA = 'always'))
  
  # QC: Filter out spots with very low read counts (potential empty spots/background)
  se <- subset(se, subset = nCount_Spatial > 100)
  seurat_list[[sample]] <- se
  
  # --- Generate Basic QC and Spatial Plots ---
  # Plot 1: Violin plot of UMI counts per spot
  plot1 <- VlnPlot(se, features = "nCount_Spatial", pt.size = 0.1) + 
    NoLegend() + 
    labs(title = paste(sample, "nCount_Spatial"))
  plot2 <- SpatialFeaturePlot(se, features = "nFeature_Spatial", pt.size.factor = 1.6, crop = TRUE) + 
    theme(legend.position = "right") +
    labs(title = paste(sample, "nFeature_Spatial"))
  
  pdf(paste0(outdir, "/HE_Labels_", sample, ".pdf"), height = 10, width = 12)
  plot3 <- SpatialDimPlot(
    se,
    group.by = 'Label',
    pt.size.factor = 1.8,
    crop = TRUE,
    alpha = 0, # Makes the spots invisible, just showing boundaries/color legend
    cols = c("CRC" = "darkred", "Normal" = "#5268B1", "LGD" = "orange", "HGD" = "pink")
  ) + NoLegend()
  plot4 <- SpatialDimPlot(
    se,
    group.by = 'Label',
    pt.size.factor = 3,
    crop = TRUE,
    alpha = 1,
    cols = c("CRC" = "darkred", "Normal" = "#5268B1", "LGD" = "orange", "HGD" = "pink")
  ) + theme(legend.position = "right")
  print(wrap_plots(plot3, plot4))
  dev.off()
  
  print(paste0(sample, " Initial QC Plots Done!"))
}


#  2. PREPROCESSING AND DIMENSIONALITY REDUCTION ####
seurat_list_processed <- list()
for (i in 1:length(samples_of_interest)) { 
  sample <- samples_of_interest[i]
  print(paste("Processing dimensionality reduction for sample:", sample))
  
  # Load the object again
  load(paste0(inputdir, "/", sample, "_annotated21.Rdata"))
  se <- subset(se, subset = nCount_Spatial > 100)  
  
  # Filter to keep only spots with defined pathology labels
  se <- subset(se, subset = Label == "CRC" | Label == "Normal" | Label == "LGD" | Label == "HGD")
  
  # Standard Seurat Preprocessing Workflow
  se <- NormalizeData(se)                
  se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)  
  se <- ScaleData(se)                    
  se <- RunPCA(se, npcs = 30)            
  se <- FindNeighbors(se, dims = 1:30)   
  se <- FindClusters(se, resolution = 0.5)  
  se <- RunUMAP(se, dims = 1:30)         
  
  # --- Store/Plot Results ---
  p_1 <- DimPlot(se, group.by = "seurat_clusters", reduction = "umap") + 
    ggtitle(paste(sample, "Clusters"))
  p_2 <- DimPlot(
    se, 
    group.by = "Label", 
    reduction = "umap", 
    cols = c("CRC" = "darkred", "Normal" = "#5268B1", "LGD" = "orange", "HGD" = "pink")
  ) + ggtitle(paste(sample, "Labels"))
  
  # Save combined UMAP plots
  pdf(paste0(outdir, "/umap_label_", sample, ".pdf"), height = 5, width = 10)
  print(wrap_plots(p_1 + p_2))
  dev.off()
  
  # Update metadata for sample identification
  se$BC <- rownames(se@meta.data) # Barcode as a column
  se$orig.ident <- sample         # Sample name
  Idents(se) <- se$BC             # Set Ident to Barcode (often useful for merged object)
  
  # Store the processed object
  seurat_list_processed[[sample]] <- se
  
  print(paste0(sample, " Processing Done!"))
}

# Save the final lists of Seurat objects
save(seurat_list, file = file.path(outdir, "seurat_list.Rdata"))
save(seurat_list_processed, file = file.path(outdir, "seurat_list_processed_Epi_subset.Rdata"))


# 3. GENE EXPRESSION CORRELATION ANALYSIS  #####
# Correlate gene expression of GOI with the progression of tissue labels
# The label order is critical for correlation: Normal < LGD < CRC
label_levels <- c("Normal", "LGD", "CRC") 
cor_results <- list()

for (sample in names(seurat_list_processed)) { 
  se <- seurat_list_processed[[sample]]
  print(paste("Calculating correlation for:", sample))
  
  se <- subset(se, subset = Label %in% label_levels)
  se$Label <- factor(se$Label, levels = label_levels)
  label_numeric <- as.numeric(se$Label) 
  expr_data <- FetchData(se, vars = goi)
  
  res <- map_dfr(goi, function(gene) {
    # Check if gene exists in the dataset (FetchData might return all NA if not found)
    if (!gene %in% colnames(expr_data) || all(is.na(expr_data[[gene]]))) {
      return(tibble(Sample = sample, Gene = gene, Correlation = NA_real_, P_value = NA_real_))
    }
    
    expr <- expr_data[[gene]]
    # Spearman correlation test between gene expression and numeric tissue progression
    cor_test <- cor.test(expr, label_numeric, method = "spearman", exact = FALSE) 
    
    tibble(
      Sample = sample,
      Gene = gene,
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value
    )
  }) %>% unique()
  
  cor_results[[sample]] <- res
}

cor_df <- bind_rows(cor_results)

# 4. RESULTS ADJUSTMENT AND SUMMARIZATION ####
# Adjust p-values for multiple testing across all samples for each gene
cor_df <- cor_df %>%
  group_by(Gene) %>%
  mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(Significant = P_adj < 0.05)

# 4.1 Determine correlation direction for each test
cor_df <- cor_df %>%
  mutate(Direction = ifelse(Correlation > 0, "Positive", "Negative"))

# 4.2 Count how many samples support the majority direction per gene
direction_summary <- cor_df %>%
  # Filter out NA correlation results (e.g., if gene was not expressed)
  filter(!is.na(Correlation)) %>% 
  group_by(Gene, Direction) %>%
  summarise(Sample_Count = n(), .groups = "drop") %>%
  group_by(Gene) %>%
  # Keep only the direction with the highest count (the 'majority' direction)
  slice_max(Sample_Count, n = 1, with_ties = FALSE) %>% 
  rename(Majority_Count = Sample_Count) %>%
  select(Gene, Majority_Count) # Only keep Gene and the Count

# 4.3 Join the summary back to the main data frame and reorder Gene factor
cor_df <- cor_df %>%
  left_join(direction_summary, by = "Gene") %>%
  # Reorder genes based on the number of samples supporting the majority direction
  mutate(Gene = fct_reorder(Gene, Majority_Count, .desc = TRUE))

# 5. VISUALIZATION AND FINAL OUTPUT ####
# Define the final order of genes for plotting  
gene_order <- cor_df %>%
  mutate(Direction = ifelse(Correlation > 0, "Positive", "Negative")) %>%
  group_by(Gene, Direction) %>%
  summarise(Sample_Count = n(), .groups = "drop") %>% 
  filter(Sample_Count > 5) %>% 
  arrange(Direction, Sample_Count) %>% 
  distinct() %>% 
  pull(Gene) 

cor_df$Gene <- factor(cor_df$Gene, levels = gene_order)

# Create Dot Plot visualization
dot_plot <- cor_df %>%
  filter(!is.na(Gene)) %>% 
  ggplot(aes(x = Sample, y = Gene, color = Correlation, shape = Significant)) +
  geom_point(size = 5) +
  
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0, 
    limits = c(-0.8, 0.8) # Set fixed limits for consistent color key
  ) +  scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1)) +
  theme_minimal() + 
  labs(
    title = "Spearman Correlation of Gene Expression with Cancer Progression (per Sample)",
    x = "Patient Sample ID",
    y = "Gene of Interest",
    color = "Spearman R", 
    shape = "BH-FDR < 0.05"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate sample names for readability
  )

pdf(file.path(outdir, "correlation_dotplot.pdf"), height = 5, width = 7)
print(dot_plot)
dev.off()

