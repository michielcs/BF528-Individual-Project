library(Seurat)
library(dplyr)
install.packages('patchwork')
library(patchwork)
library(reticulate)
Install.packages('umap')
library(umap)


val <- load("/projectnb/bf528/users/group6/project_4/output/processed_panc_with_clusters.rda")
#cells <- readRDS("/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")

pbmc.markers <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#filter by log fold change
log_filter = subset(pbmc.markers, pbmc.markers$avg_logFC>=0.58)
write.csv(log_filter,"log_filter.csv")
#filter by p value
p_filter = subset(pbmc.markers, pbmc.markers$p_val_adj<0.05)
write.csv(p_filter,"p_filter.csv")