# UKBB risk prediction####
library(STRINGdb)
library(org.Hs.eg.db)
library(igraph)
library(dplyr)
library(survival)
library(Matrix)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(readxl)
library(stringr)
library(data.table)
library(glmnet)
library(reticulate)
library(broom)
library(MatchIt)
library(VIM)
library(purrr)
library(broom)
library(pROC)
library(tidyverse)
library(DMwR2)
library(foreach)
library(doParallel)

#################
input.dir = 'data'
out.dir = 'output'
dir.create(out.dir, showWarnings = FALSE)
################
FAP_list = c("LAMB1", "LAMC2", "TGFBI", "SEMA4D" ,"AGRN" ,"SEMA5A" ,"SEMA3F" ,"INHBA", "LAMA5" ,"COL7A1", "SOX9"  ,"CCND1" ,"TGFBI" ,"MET"  , "ID1" , "PTMA",  "CLU" ,  "TGIF1", "EDN1" , "GDF15" ) # combined list of top10ligand_to10target from FAP PvsN

################
#load data####
load(file.path("data/Esplin2024NatCan_Proteo_DEG_with_symbols.RData"))

Proteo_DEG_faplist = Proteo_DEG %>% 
  filter(Symbol %in% FAP_list) %>% select(datatype,Contr, Symbol, Log2FC, `-Log10FDR`, up_down)

# plot the distribution of up/down regulated proteins across contrasts ####
# change log2FC to FC
Proteo_DEG_faplist$FC = 2^Proteo_DEG_faplist$Log2FC

ggplot(Proteo_DEG_faplist, 
       aes(x = Log2FC, y = reorder(Symbol, Log2FC), fill = Contr)) +
  geom_bar(stat = "identity") +
  #geom_text(aes(label = format(`-Log10FDR`, digits = 2)), vjust = -1, size = 3) +
  facet_wrap(~ Contr) +
  scale_fill_manual(values = c("M-B" = "#1f77b4", "M-D" = "#ff7f0e")) +
  labs(title = "Top DEGs by Log2FC", 
       x = "Log2 Fold Change", 
       y = "Gene Symbol") +
  theme_minimal()
# save as pdf
ggsave(filename = file.path(out.dir, "Proteomics_DEG_FAPlist.pdf"), 
       width = 8, height = 6)

