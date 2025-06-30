# Title: Blastoid Blastocyst Analysis —— Refine pipeline & Delete E4.5 outliers
# Author: Gaozy
# Time: 2025-04-04

# comput172-r_env-R
# 0. Basic settings ----
setwd("/home/qcao02/gaozy/blastocyst/")

.libPaths(c("/usr/lib64/R/library",
            "/usr/share/R/library",
            "/home/lushi02/.local/share/r-miniconda/envs/r-reticulate",
            "/home/lushi02/anaconda3/lib/R/library"))

# 1. library ----
library(Matrix, lib.loc = "/usr/lib64/R/library") 
library(sva, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4") 
library(Seurat, lib.loc = "/usr/lib64/R/library") 
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggcorrplot, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4")
library(Hmisc)
library(dplyr, lib.loc = "/usr/lib64/R/library")
library(reshape2)
library(ggrepel)
library(biomaRt)
library(patchwork)  # For combining plots

if (dir.exists("250404") == F) dir.create('250404')

# 2. load data ----
corrected_seurat_combat_raw <- readRDS("250326/combat.rds")

## 1) delete outliers ----
pca_coords <- Embeddings(corrected_seurat_combat_raw, "pca")
cells_to_keep <- rownames(pca_coords)[!(pca_coords[, "PC_1"] > 20 & pca_coords[, "PC_2"] < -110)]
corrected_seurat_combat <- subset(corrected_seurat_combat_raw, cells = cells_to_keep)

Idents(corrected_seurat_combat) <- "cell_info"
corrected_seurat_combat <- FindVariableFeatures(corrected_seurat_combat, selection.method = "vst", nfeatures = 2000)
corrected_seurat_combat <- ScaleData(corrected_seurat_combat, features = rownames(corrected_seurat_combat)) %>%
  RunPCA(npcs=50, features = rownames(corrected_seurat_combat))

## 2) choose proper PCs ----
corrected_seurat_combat <- JackStraw(corrected_seurat_combat, num.replicate = 50, dims = 50)
corrected_seurat_combat <- ScoreJackStraw(corrected_seurat_combat, dims = 1:50) 

p1 <- JackStrawPlot(corrected_seurat_combat, dims = 1:50)
p2 <- ElbowPlot(corrected_seurat_combat, ndims = 50)

pdf(paste0("250404/corrected_pca_elbow_deleteOuteliers.pdf"), width = 12, height = 8, onefile = F)
p1 / p2
dev.off()

corrected_seurat_combat <- ScaleData(corrected_seurat_combat, features = rownames(corrected_seurat_combat)) %>%
  RunPCA(npcs=30, features = rownames(corrected_seurat_combat))

corrected_seurat_combat$cell_info <- factor(corrected_seurat_combat$cell_info, levels = unique(corrected_seurat_combat$cell_info))
## Get PCA coordinates and metadata
pca_var <- corrected_seurat_combat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( corrected_seurat_combat, reduction = "pca", dims = c(1, 2), group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 8) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )+ guides(color = guide_legend(ncol = 1))

pdf("250404/corrected_seurat_combat_PCA.pdf", width=10, height=9, onefile = F)
print(p)
dev.off()

## 3) UMAP & tSNE ----
corrected_seurat_combat <- RunTSNE(corrected_seurat_combat, dims = 1:30)
corrected_seurat_combat <- RunUMAP(corrected_seurat_combat, dims = 1:30)

p1 <- DimPlot(corrected_seurat_combat, reduction = "umap", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 3)+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(corrected_seurat_combat, reduction = "tsne", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 3)+ guides(color = guide_legend(ncol = 1))

pdf("250404/corrected_seurat_combat_umap_tsne.pdf", width=25, height=12, onefile = F)
print(p1 | p2 )
dev.off()

saveRDS(corrected_seurat_combat, "250404/corrected_seurat_combat.rds")

dataset <- corrected_seurat_combat
output_folder <- "250404"
# 3. Cell info & HVG ----
ident_by <- "cell_info"
use_all_genes <- FALSE

## 1) Pseudobulk data ----
Idents(dataset) <- ident_by

on.exit({
  while(dev.cur() > 1) dev.off()
})

file_names_prefix <- paste0(output_folder, "/", ident_by, "-", ifelse(use_all_genes, "All_genes", "Hvg"), "-")

# Pseudobulk
pb_seurat <- Seurat:::PseudobulkExpression(dataset, return.seurat = T, method = "aggregate", # "average"
                                           normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat[[ident_by]] <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- ident_by

# Variable features
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)

## 2) PCA with selected PCs ----
features_to_use <- if (use_all_genes) rownames(pb_seurat) else VariableFeatures(pb_seurat)
pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 14, features = features_to_use)

pb_seurat <- JackStraw(pb_seurat, num.replicate = 10, dims = 14)
pb_seurat <- ScoreJackStraw(pb_seurat, dims = 1:14) 

p1 <- JackStrawPlot(pb_seurat, dims = 1:14)
p2 <- ElbowPlot(pb_seurat, ndims = 14)

pdf(paste0(file_names_prefix, "pca_elbow.pdf"), width = 12, height = 5, onefile = F)
p1 | p2
dev.off()

pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 6, features = features_to_use)

## Get PCA coordinates and metadata
pca_var <- pb_seurat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( pb_seurat, reduction = "pca", dims = c(1, 2), group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=8, height=6, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=25, height=23, onefile = F)
}
print(p)
dev.off()

## 3) UMAP & tSNE ----
pb_seurat <- RunTSNE(pb_seurat, dims = 1:6, perplexity = 3)
pb_seurat <- RunUMAP(pb_seurat, dims = 1:6, n.neighbors = 2)

p1 <- DimPlot(pb_seurat, reduction = "umap", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(pb_seurat, reduction = "tsne", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=20, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=33, height=15, onefile = F)
}
print(p1 | p2 )
dev.off()

## 4) Clustering ----
pb_seurat <- FindNeighbors(pb_seurat, dims = 1:6)
pb_seurat <- FindClusters(pb_seurat, resolution = 0.5)

## use scale.data
data_for_clustering <- t(GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data"))

## calculate distance
dist_matrix <- dist(data_for_clustering, method = "euclidean")

## hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

## output
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=6, height=6, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=12, height=7, onefile = F)
}
print(plot(hclust_result, main = "Transcriptome-based clustering", xlab = "", sub = "", cex = 0.9, hang = -1))
dev.off()

## 5) Pearson correlation ----
counts_mat <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")

## Extract relevant columns for blastoid and E3.5-6.5 embryo
blastoid_columns <- grep("blastoid", colnames(counts_mat), value = T)  # Select columns containing "blastoid"
embryo_columns <- grep("E[3-6]\\.5|In house-blastocyst", colnames(counts_mat), value = T)  # Select columns for E3.5-6.5 embryos

## Subset the data to focus on these columns (blastoid samples as x-axis, E3.5-6.5 embryos as y-axis)
data_subset <- counts_mat[, c(blastoid_columns, embryo_columns)]

## Perform Pearson correlation on the subsetted data
correlation_matrix <- cor(data_subset[, embryo_columns],
                          data_subset[, blastoid_columns],
                          method = "pearson")

## corr dotplot
corr_long <- melt(correlation_matrix)
colnames(corr_long) <- c("Sample1", "Sample2", "Correlation")

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=6, height=5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=11, height=6, onefile = F)
}
print(  ggplot(corr_long, aes(x = Sample2, y = Sample1, size = (Correlation), color = Correlation)) +
          geom_point() +
          scale_size_continuous(range = c(2, 10)) +  # 调整点的大小范围
          scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)) +  # 设置颜色渐变
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),) +  # 旋转 x 轴标签
          labs(title = "Correlation Dotplot", x = "", y = "", color = "Correlation") +
          guides(size = "none") )
dev.off()

## 6) post-implantation stage-specific genes ----
## Find markers for post-implantation stages
markers_post_implantation <- FindMarkers(
  pb_seurat, 
  slot = "data",
  min.cells.group = 2,
  ident.1 = c( "E5.5", "E6.5"),
  ident.2 = c("E3.5", "In house-blastocyst", "public Blastocyst-SSII"), 
  only.pos = T,
  logfc.threshold = 0.25)

markers_post_implantation <- subset(markers_post_implantation, #p_val < 0.05 & 
                                    avg_log2FC > 0.5)

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_data <- exp_data[intersect(rownames(exp_data), rownames(markers_post_implantation)), ]
colnames(exp_data) <- colnames(pb_seurat)

meta <- pb_seurat@meta.data
exp_scale = exp_data

## 7) hclust by post-implantation stage-specific genes ----
### a. subset genes high expressed in E4.5-E6.5 while low expressed in blastocyst ----
post_implantation_genes <- rownames(exp_scale)

high_exp_post_implantation <- apply(exp_scale[intersect(rownames(exp_data), rownames(markers_post_implantation)),], 1, function(x) {
  # Mean expression in E5.5, E6.5 should be greater than in blastocyst stages
  high_exp_condition <- mean(x[colnames(exp_scale) %in% c("E5.5", "E6.5")]) > 
    mean(x[grep("E3.5|blastocyst|Blastocyst",colnames(exp_scale),value = T)])
  
  return(high_exp_condition )
})

# Filter marker genes with higher expression in post-implantation stages and lower in blastocyst
filtered_markers <- post_implantation_genes[high_exp_post_implantation]

### b. plot heatmap ---- 
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=8, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=10, height=15, onefile = F)
}
print(
  Heatmap(as.matrix(exp_scale[filtered_markers,]),
          column_split = meta[ident_by],
          show_row_names = F,
          show_column_names = T,  # Show column names
          show_column_dend = T,  # Show column dendrogram
          column_title = "Post-implantation stage-specific genes expression",
          use_raster = F,
          col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
          heatmap_legend_param = list(title='')))
dev.off()

## 8) Similarity ----
### a. train data (blastocyst) ----
train <- subset(dataset, sample_info %in% as.character(
  grep("E3.5|E4.5|E5.5|E6.5|blastocyst",dataset$sample_info, value = T)))

RhpcBLASctl::blas_set_num_threads(10)
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data (all) ----
test <- dataset

### c. similarity ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

Test_similarity =
  glm.predict(
    train.data = as.matrix(GetAssayData(train, assay = "RNA", slot = "data")),
    train.group = as.character(train$cell_info),
    # genes.used = rownames(train),
    genes.used = VariableFeatures(train),
    downsample = F,
    test.data = as.matrix(GetAssayData(test, assay = "RNA", slot = "data")),
    test.group = as.character(test$cell_info),
    alpha = 0.99,
    nfolds = 10,
    seed = 123)

heatmap_mat = Test_similarity$heatmap@matrix[,c(
  "E3.5", "E4.5", "E5.5", "E6.5",
  "In house-blastocyst", "public Blastocyst-SSII", 
  "EPSC-public blastoid", "EPSC-public blastoid-ssII", 
  "2MYCP-Blastoid", "2MYCP-Blastoid-bullk", "EPT-blastoid", "ET-blastoid",
  "TBLC-blastoid", "TPSC-blastoid", "In house-blastoid"
)]


if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=7, height=3.5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=10, height=3, onefile = F)
}
print(  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                           color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                           cluster_rows = F))
dev.off()

## 9) Marker genes expr heatmap ----
### a. load marker list ----
marker <- openxlsx::read.xlsx('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/20240221/data/SSII-seq sample 信息及各个lineage marker.xlsx',
                              sheet = 2, colNames = F)
colnames(marker) <- c('type', 'gene')
rownames(marker) <- marker$gene

### b. Marker Expression Heatmap - Scale Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=8, height=8, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=10, height=8, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

### c. Marker Expression Heatmap - Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=8, height=9, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=10, height=8, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

## 10) save data ----
saveRDS(pb_seurat, paste0(file_names_prefix, "data.rds"))
write.csv(markers_post_implantation, paste0(file_names_prefix, "DEG.csv"), quote = F, row.names = T)

# 4. Cell info & All genes ----
ident_by <- "cell_info"
use_all_genes <- TRUE

## 1) Pseudobulk data ----
Idents(dataset) <- ident_by

on.exit({
  while(dev.cur() > 1) dev.off()  # 确保所有图形设备关闭
})

file_names_prefix <- paste0(output_folder, "/", ident_by, "-", ifelse(use_all_genes, "All_genes", "Hvg"), "-")

# Pseudobulk
pb_seurat <- Seurat:::PseudobulkExpression(dataset, return.seurat = T, method = "aggregate", # "average"
                                           normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat[[ident_by]] <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- ident_by

# Variable features
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)

## 2) PCA with selected PCs ----
features_to_use <- if (use_all_genes) rownames(pb_seurat) else VariableFeatures(pb_seurat)
pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 14, features = features_to_use)

pb_seurat <- JackStraw(pb_seurat, num.replicate = 10, dims = 14)
pb_seurat <- ScoreJackStraw(pb_seurat, dims = 1:14) 

p1 <- JackStrawPlot(pb_seurat, dims = 1:14)
p2 <- ElbowPlot(pb_seurat, ndims = 14)

pdf(paste0(file_names_prefix, "pca_elbow.pdf"), width = 12, height = 5, onefile = F)
p1 | p2
dev.off()

pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 11, features = features_to_use)

## Get PCA coordinates and metadata
pca_var <- pb_seurat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( pb_seurat, reduction = "pca", dims = c(1, 2), group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 8) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=9, height=7, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=25, height=23, onefile = F)
}
print(p)
dev.off()

## 3) UMAP & tSNE ----
pb_seurat <- RunTSNE(pb_seurat, dims = 1:11, perplexity = 3)
pb_seurat <- RunUMAP(pb_seurat, dims = 1:11, n.neighbors = 2)

p1 <- DimPlot(pb_seurat, reduction = "umap", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(pb_seurat, reduction = "tsne", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=20, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=33, height=15, onefile = F)
}
print(p1 | p2 )
dev.off()

## 4) Clustering ----
pb_seurat <- FindNeighbors(pb_seurat, dims = 1:11)
pb_seurat <- FindClusters(pb_seurat, resolution = 0.5)

## use scale.data
data_for_clustering <- t(GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data"))

## calculate distance
dist_matrix <- dist(data_for_clustering, method = "euclidean")

## hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

## output
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=6, height=6, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=12, height=7, onefile = F)
}
print(plot(hclust_result, main = "Transcriptome-based clustering", xlab = "", sub = "", cex = 0.9, hang = -1))
dev.off()

## 5) Pearson correlation ----
counts_mat <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")

## Extract relevant columns for blastoid and E3.5-6.5 embryo
blastoid_columns <- grep("blastoid", colnames(counts_mat), value = T)  # Select columns containing "blastoid"
embryo_columns <- grep("E[3-6]\\.5|In house-blastocyst", colnames(counts_mat), value = T)  # Select columns for E3.5-6.5 embryos

## Subset the data to focus on these columns (blastoid samples as x-axis, E3.5-6.5 embryos as y-axis)
data_subset <- counts_mat[, c(blastoid_columns, embryo_columns)]

## Perform Pearson correlation on the subsetted data
correlation_matrix <- cor(data_subset[, embryo_columns],
                          data_subset[, blastoid_columns],
                          method = "pearson")

## corr dotplot
corr_long <- melt(correlation_matrix)
colnames(corr_long) <- c("Sample1", "Sample2", "Correlation")

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=6, height=5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=11, height=6, onefile = F)
}
print(  ggplot(corr_long, aes(x = Sample2, y = Sample1, size = (Correlation), color = Correlation)) +
          geom_point() +
          scale_size_continuous(range = c(2, 10)) +  # 调整点的大小范围
          scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)) +  # 设置颜色渐变
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),) +  # 旋转 x 轴标签
          labs(title = "Correlation Dotplot", x = "", y = "", color = "Correlation") +
          guides(size = "none") )
dev.off()

## 6) post-implantation stage-specific genes ----
## Find markers for post-implantation stages
markers_post_implantation <- FindMarkers(
  pb_seurat, 
  slot = "data",
  min.cells.group = 2,
  ident.1 = c( "E5.5", "E6.5"),
  ident.2 = c("E3.5", "In house-blastocyst", "public Blastocyst-SSII"), 
  only.pos = T,
  logfc.threshold = 0.25)

markers_post_implantation <- subset(markers_post_implantation, #p_val < 0.05 & 
                                    avg_log2FC > 0.5)

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_data <- exp_data[intersect(rownames(exp_data), rownames(markers_post_implantation)), ]
colnames(exp_data) <- colnames(pb_seurat)

meta <- pb_seurat@meta.data
exp_scale = exp_data

## 7) hclust by post-implantation stage-specific genes ----
### a. subset genes high expressed in E4.5-E6.5 while low expressed in blastocyst ----
post_implantation_genes <- rownames(exp_scale)

high_exp_post_implantation <- apply(exp_scale[intersect(rownames(exp_data), rownames(markers_post_implantation)),], 1, function(x) {
  # Mean expression in E5.5, E6.5 should be greater than in blastocyst stages
  high_exp_condition <- mean(x[colnames(exp_scale) %in% c("E5.5", "E6.5")]) > 
    mean(x[grep("E3.5|blastocyst|Blastocyst",colnames(exp_scale),value = T)])
  
  return(high_exp_condition )
})

# Filter marker genes with higher expression in post-implantation stages and lower in blastocyst
filtered_markers <- post_implantation_genes[high_exp_post_implantation]

### b. plot heatmap ---- 
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=8, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=10, height=15, onefile = F)
}
print(
  Heatmap(as.matrix(exp_scale[filtered_markers,]),
          column_split = meta[ident_by],
          show_row_names = F,
          show_column_names = T,  # Show column names
          show_column_dend = T,  # Show column dendrogram
          column_title = "Post-implantation stage-specific genes expression",
          use_raster = F,
          col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
          heatmap_legend_param = list(title='')))
dev.off()

## 8) Similarity ----
### a. train data (blastocyst) ----
train <- subset(dataset, sample_info %in% as.character(
  grep("E3.5|E4.5|E5.5|E6.5|blastocyst",dataset$sample_info, value = T)))

RhpcBLASctl::blas_set_num_threads(10)
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data (all) ----
test <- dataset

### c. similarity ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

Test_similarity =
  glm.predict(
    train.data = as.matrix(GetAssayData(train, assay = "RNA", slot = "data")),
    train.group = as.character(train$cell_info),
    # genes.used = rownames(train),
    genes.used = VariableFeatures(train),
    downsample = F,
    test.data = as.matrix(GetAssayData(test, assay = "RNA", slot = "data")),
    test.group = as.character(test$cell_info),
    alpha = 0.99,
    nfolds = 10,
    seed = 123)

heatmap_mat = Test_similarity$heatmap@matrix[,c(
  "E3.5", "E4.5", "E5.5", "E6.5",
  "In house-blastocyst", "public Blastocyst-SSII",
  "EPSC-public blastoid", "EPSC-public blastoid-ssII", 
  "2MYCP-Blastoid", "2MYCP-Blastoid-bullk", "EPT-blastoid", "ET-blastoid",
  "TBLC-blastoid", "TPSC-blastoid", "In house-blastoid"
)]


if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=7, height=3.5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=10, height=3, onefile = F)
}
print(  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                           color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                           cluster_rows = F))
dev.off()

## 9) Marker genes expr heatmap ----
### a. load marker list ----
marker <- openxlsx::read.xlsx('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/20240221/data/SSII-seq sample 信息及各个lineage marker.xlsx',
                              sheet = 2, colNames = F)
colnames(marker) <- c('type', 'gene')
rownames(marker) <- marker$gene

### b. Marker Expression Heatmap - Scale Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=8, height=8, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=10, height=8, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

### c. Marker Expression Heatmap - Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=8, height=9, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=10, height=8, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

## 10) save data ----
saveRDS(pb_seurat, paste0(file_names_prefix, "data.rds"))
write.csv(markers_post_implantation, paste0(file_names_prefix, "DEG.csv"), quote = F, row.names = T)

# 5. Sample info & HVG ----
ident_by <- "sample_info"
use_all_genes <- FALSE

## 1) Pseudobulk data ----
Idents(dataset) <- ident_by

on.exit({
  while(dev.cur() > 1) dev.off()  # 确保所有图形设备关闭
})

file_names_prefix <- paste0(output_folder, "/", ident_by, "-", ifelse(use_all_genes, "All_genes", "Hvg"), "-")

# Pseudobulk
pb_seurat <- Seurat:::PseudobulkExpression(dataset, return.seurat = T, method = "aggregate", # "average"
                                           normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat[[ident_by]] <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- ident_by

# Variable features
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)

## 2) PCA with selected PCs ----
features_to_use <- if (use_all_genes) rownames(pb_seurat) else VariableFeatures(pb_seurat)
pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 39, features = features_to_use)

pb_seurat <- JackStraw(pb_seurat, num.replicate = 10, dims = 39)
pb_seurat <- ScoreJackStraw(pb_seurat, dims = 1:39) 

p1 <- JackStrawPlot(pb_seurat, dims = 1:39)
p2 <- ElbowPlot(pb_seurat, ndims = 39)

pdf(paste0(file_names_prefix, "pca_elbow.pdf"), width = 12, height = 5, onefile = F)
p1 | p2
dev.off()

pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 26, features = features_to_use)

## Get PCA coordinates and metadata
pca_var <- pb_seurat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( pb_seurat, reduction = "pca", dims = c(1, 2), group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 8) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=9, height=7, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=30, height=28, onefile = F)
}
print(p)
dev.off()

## 3) UMAP & tSNE ----
pb_seurat <- RunTSNE(pb_seurat, dims = 1:26, perplexity = 5)
pb_seurat <- RunUMAP(pb_seurat, dims = 1:26, n.neighbors = 10)

p1 <- DimPlot(pb_seurat, reduction = "umap", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(pb_seurat, reduction = "tsne", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=20, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=33, height=15, onefile = F)
}
print(p1 | p2 )
dev.off()

## 4) Clustering ----
pb_seurat <- FindNeighbors(pb_seurat, dims = 1:26)
pb_seurat <- FindClusters(pb_seurat, resolution = 0.5)

## use scale.data
data_for_clustering <- t(GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data"))

## calculate distance
dist_matrix <- dist(data_for_clustering, method = "euclidean")

## hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

## output
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=6, height=6, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=12, height=7, onefile = F)
}
print(plot(hclust_result, main = "Transcriptome-based clustering", xlab = "", sub = "", cex = 0.9, hang = -1))
dev.off()

## 5) Pearson correlation ----
counts_mat <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")

## Extract relevant columns for blastoid and E3.5-6.5 embryo
blastoid_columns <- grep("blastoid", colnames(counts_mat), value = T)  # Select columns containing "blastoid"
embryo_columns <- grep("E[3-6]\\.5|In house-blastocyst", colnames(counts_mat), value = T)  # Select columns for E3.5-6.5 embryos

## Subset the data to focus on these columns (blastoid samples as x-axis, E3.5-6.5 embryos as y-axis)
data_subset <- counts_mat[, c(blastoid_columns, embryo_columns)]

## Perform Pearson correlation on the subsetted data
library(tibble)
correlation_matrix <- cor(data_subset[, embryo_columns],
                          data_subset[, blastoid_columns],
                          method = "pearson") %>% as.data.frame %>% rownames_to_column("sample")

in_house_mean <- correlation_matrix %>%
  filter(grepl("In house-blastocyst_", sample)) %>%
  dplyr::select(-sample) %>%
  summarise(across(everything(), mean)) %>%
  mutate(sample = "In house-blastocyst")

# 过滤掉原始的 blastocyst 行并添加均值行
new_cor_df <- correlation_matrix %>%
  filter(!grepl("In house-blastocyst_", sample)) %>%
  bind_rows(in_house_mean) %>%
  column_to_rownames("sample") %>% as.matrix

## corr dotplot
corr_long <- melt(new_cor_df)
colnames(corr_long) <- c("Sample1", "Sample2", "Correlation")

corr_long$Sample1 <- factor(corr_long$Sample1, levels = c("In house-blastocyst","E6.5", "E5.5", "E4.5", "E3.5"))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=6, height=5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=11, height=4, onefile = F)
}
print(  ggplot(corr_long, aes(x = Sample2, y = Sample1, size = (Correlation), color = Correlation)) +
          geom_point() +
          scale_size_continuous(range = c(2, 10)) +  # 调整点的大小范围
          scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)) +  # 设置颜色渐变
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),) +  # 旋转 x 轴标签
          labs(title = "Correlation Dotplot", x = "", y = "", color = "Correlation") +
          guides(size = "none") )
dev.off()

## 6) post-implantation stage-specific genes ----
## Find markers for post-implantation stages
markers_post_implantation <- FindMarkers(
  pb_seurat, 
  slot = "data",
  min.cells.group = 2,
  ident.1 = grep("E5.5|E6.5", colnames(pb_seurat), value = T),
  ident.2 = grep("E3.5|In house-blastocyst|public Blastocyst-SSII", colnames(pb_seurat), value = T), 
  only.pos = T,
  logfc.threshold = 0.25)

markers_post_implantation <- subset(markers_post_implantation, #p_val < 0.05 & 
                                    avg_log2FC > 0.5)

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_data <- exp_data[intersect(rownames(exp_data), rownames(markers_post_implantation)), ]
colnames(exp_data) <- colnames(pb_seurat)

meta <- pb_seurat@meta.data
exp_scale = exp_data

## 7) hclust by post-implantation stage-specific genes ----
### a. subset genes high expressed in E4.5-E6.5 while low expressed in blastocyst ----
post_implantation_genes <- rownames(exp_scale)

high_exp_post_implantation <- apply(exp_scale[intersect(rownames(exp_data), rownames(markers_post_implantation)),], 1, function(x) {
  # Mean expression in E5.5, E6.5 should be greater than in blastocyst stages
  high_exp_condition <- mean(x[colnames(exp_scale) %in% c("E5.5", "E6.5")]) > 
    mean(x[grep("E3.5|blastocyst|Blastocyst",colnames(exp_scale),value = T)])
  
  return(high_exp_condition )
})

# Filter marker genes with higher expression in post-implantation stages and lower in blastocyst
filtered_markers <- post_implantation_genes[high_exp_post_implantation]

### b. plot heatmap ---- 
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=8, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=10, height=15, onefile = F)
}
print(
  Heatmap(as.matrix(exp_scale[filtered_markers,]),
          column_split = meta[ident_by],
          show_row_names = F,
          show_column_names = T,  # Show column names
          show_column_dend = T,  # Show column dendrogram
          column_title = "Post-implantation stage-specific genes expression",
          use_raster = F,
          col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
          heatmap_legend_param = list(title='')))
dev.off()

## 8) Similarity ----
### a. train data (blastocyst) ----
train <- subset(dataset, sample_info %in% as.character(
  grep("E3.5|E4.5|E5.5|E6.5|blastocyst",dataset$sample_info, value = T)))

RhpcBLASctl::blas_set_num_threads(10)
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data (all) ----
test <- dataset

### c. similarity ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

if (ident_by == "sample_info") {
  train$simi_info <- train$cell_info
  test$simi_info <- gsub("*-(\\d+)$", "", test$sample_info)
  
  Test_similarity =
    glm.predict(
      train.data = as.matrix(GetAssayData(train, assay = "RNA", slot = "data")),
      train.group = as.character(train$simi_info),
      # genes.used = rownames(train),
      genes.used = VariableFeatures(train),
      downsample = T,
      test.data = as.matrix(GetAssayData(test, assay = "RNA", slot = "data")),
      test.group = as.character(test$simi_info),
      alpha = 0.99,
      nfolds = 10,
      seed = 123)
  
  heatmap_mat = Test_similarity$heatmap@matrix[, c(
    "E3.5", "E4.5", "E5.5", "E6.5",
    
    "In house-blastocyst_1", "In house-blastocyst_2", "In house-blastocyst_3", 
    "In house-blastocyst_4", "In house-blastocyst_5", "In house-blastocyst_6", 
    
    "public Blastocyst-SSII_1", "public Blastocyst-SSII_2", "public Blastocyst-SSII_3", 
    
    "EPSC-public blastoid", "EPSC-public blastoid-ssII_1", "EPSC-public blastoid-ssII_2", 
    "EPSC-public blastoid-ssII_3", "EPSC-public blastoid-ssII_4", "EPSC-public blastoid-ssII_5", 
    "EPSC-public blastoid-ssII_6", "EPSC-public blastoid-ssII_7", "EPSC-public blastoid-ssII_8", 
    "2MYCP-Blastoid", "2MYCP-Blastoid-bullk_1", #"2MYCP-Blastoid-bullk_2", 
    "EPT-blastoid_1", "EPT-blastoid_2", "EPT-blastoid_3", 
    "EPT-blastoid_4", "EPT-blastoid_5", "ET-blastoid", 
    "TBLC-blastoid_1", "TBLC-blastoid_2", "TPSC-blastoid",
    
    "In house-blastoid_1", "In house-blastoid_2", "In house-blastoid_3", 
    "In house-blastoid_4", "In house-blastoid_5", "In house-blastoid_6"
  )]
}

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=7, height=3.5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=10, height=3, onefile = F)
}
print(  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                           color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                           cluster_rows = F))
dev.off()

## 9) Marker genes expr heatmap ----
### a. load marker list ----
marker <- openxlsx::read.xlsx('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/20240221/data/SSII-seq sample 信息及各个lineage marker.xlsx',
                              sheet = 2, colNames = F)
colnames(marker) <- c('type', 'gene')
rownames(marker) <- marker$gene

### b. Marker Expression Heatmap - Scale Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=8, height=8, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=10, height=8, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

### c. Marker Expression Heatmap - Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=8, height=9, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=10, height=10, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

## 10) save data ----
saveRDS(pb_seurat, paste0(file_names_prefix, "data.rds"))
write.csv(markers_post_implantation, paste0(file_names_prefix, "DEG.csv"), quote = F, row.names = T)

# 6. Sample info & All genes ----
ident_by <- "sample_info"
use_all_genes <- TRUE

## 1) Pseudobulk data ----
Idents(dataset) <- ident_by

on.exit({
  while(dev.cur() > 1) dev.off()  # 确保所有图形设备关闭
})

file_names_prefix <- paste0(output_folder, "/", ident_by, "-", ifelse(use_all_genes, "All_genes", "Hvg"), "-")

# Pseudobulk
pb_seurat <- Seurat:::PseudobulkExpression(dataset, return.seurat = T, method = "aggregate", # "average"
                                           normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat[[ident_by]] <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- ident_by

# Variable features
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)

## 2) PCA with selected PCs ----
features_to_use <- if (use_all_genes) rownames(pb_seurat) else VariableFeatures(pb_seurat)
pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 39, features = features_to_use)

pb_seurat <- JackStraw(pb_seurat, num.replicate = 10, dims = 39)
pb_seurat <- ScoreJackStraw(pb_seurat, dims = 1:39) 

p1 <- JackStrawPlot(pb_seurat, dims = 1:39)
p2 <- ElbowPlot(pb_seurat, ndims = 39)

pdf(paste0(file_names_prefix, "pca_elbow.pdf"), width = 12, height = 5, onefile = F)
p1 | p2
dev.off()

pb_seurat <- ScaleData(pb_seurat, features = features_to_use) %>%
  RunPCA(npcs = 34, features = features_to_use)

## Get PCA coordinates and metadata
pca_var <- pb_seurat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( pb_seurat, reduction = "pca", dims = c(1, 2), group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 8) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=9, height=7, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "PCA.pdf"), width=30, height=28, onefile = F)
}
print(p)
dev.off()

## 3) UMAP & tSNE ----
pb_seurat <- RunTSNE(pb_seurat, dims = 1:34, perplexity = 5)
pb_seurat <- RunUMAP(pb_seurat, dims = 1:34, n.neighbors = 10)

p1 <- DimPlot(pb_seurat, reduction = "umap", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))
p2 <- DimPlot(pb_seurat, reduction = "tsne", group.by = ident_by, label = TRUE, repel = TRUE, pt.size = 5)+ guides(color = guide_legend(ncol = 1))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=20, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "umap_tsne.pdf"), width=33, height=15, onefile = F)
}
print(p1 | p2 )
dev.off()

## 4) Clustering ----
pb_seurat <- FindNeighbors(pb_seurat, dims = 1:34)
pb_seurat <- FindClusters(pb_seurat, resolution = 0.5)

## use scale.data
data_for_clustering <- t(GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data"))

## calculate distance
dist_matrix <- dist(data_for_clustering, method = "euclidean")

## hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")

## output
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=6, height=6, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "hierachical_cluster.pdf"), width=12, height=7, onefile = F)
}
print(plot(hclust_result, main = "Transcriptome-based clustering", xlab = "", sub = "", cex = 0.9, hang = -1))
dev.off()

## 5) Pearson correlation ----
counts_mat <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")

## Extract relevant columns for blastoid and E3.5-6.5 embryo
blastoid_columns <- grep("blastoid", colnames(counts_mat), value = T)  # Select columns containing "blastoid"
embryo_columns <- grep("E[3-6]\\.5|In house-blastocyst", colnames(counts_mat), value = T)  # Select columns for E3.5-6.5 embryos

## Subset the data to focus on these columns (blastoid samples as x-axis, E3.5-6.5 embryos as y-axis)
data_subset <- counts_mat[, c(blastoid_columns, embryo_columns)]

## Perform Pearson correlation on the subsetted data
library(tibble)
correlation_matrix <- cor(data_subset[, embryo_columns],
                          data_subset[, blastoid_columns],
                          method = "pearson") %>% as.data.frame %>% rownames_to_column("sample")

in_house_mean <- correlation_matrix %>%
  filter(grepl("In house-blastocyst_", sample)) %>%
  dplyr::select(-sample) %>%
  summarise(across(everything(), mean)) %>%
  mutate(sample = "In house-blastocyst")

# 过滤掉原始的 blastocyst 行并添加均值行
new_cor_df <- correlation_matrix %>%
  filter(!grepl("In house-blastocyst_", sample)) %>%
  bind_rows(in_house_mean) %>%
  column_to_rownames("sample") %>% as.matrix

## corr dotplot
corr_long <- melt(new_cor_df)
colnames(corr_long) <- c("Sample1", "Sample2", "Correlation")

corr_long$Sample1 <- factor(corr_long$Sample1, levels = c("In house-blastocyst","E6.5", "E5.5", "E4.5", "E3.5"))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=6, height=5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "Corr.pdf"), width=11, height=4, onefile = F)
}
print(  ggplot(corr_long, aes(x = Sample2, y = Sample1, size = (Correlation), color = Correlation)) +
          geom_point() +
          scale_size_continuous(range = c(2, 10)) +  # 调整点的大小范围
          scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)) +  # 设置颜色渐变
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),) +  # 旋转 x 轴标签
          labs(title = "Correlation Dotplot", x = "", y = "", color = "Correlation") +
          guides(size = "none") )
dev.off()

## 6) post-implantation stage-specific genes ----
## Find markers for post-implantation stages
markers_post_implantation <- FindMarkers(
  pb_seurat, 
  slot = "data",
  min.cells.group = 2,
  ident.1 = grep("E5.5|E6.5", colnames(pb_seurat), value = T),
  ident.2 = grep("E3.5|In house-blastocyst|public Blastocyst-SSII", colnames(pb_seurat), value = T), 
  only.pos = T,
  logfc.threshold = 0.25)

markers_post_implantation <- subset(markers_post_implantation, #p_val < 0.05 & 
                                    avg_log2FC > 0.5)

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_data <- exp_data[intersect(rownames(exp_data), rownames(markers_post_implantation)), ]
colnames(exp_data) <- colnames(pb_seurat)

meta <- pb_seurat@meta.data
exp_scale = exp_data

## 7) hclust by post-implantation stage-specific genes ----
### a. subset genes high expressed in E4.5-E6.5 while low expressed in blastocyst ----
post_implantation_genes <- rownames(exp_scale)

high_exp_post_implantation <- apply(exp_scale[intersect(rownames(exp_data), rownames(markers_post_implantation)),], 1, function(x) {
  # Mean expression in E5.5, E6.5 should be greater than in blastocyst stages
  high_exp_condition <- mean(x[colnames(exp_scale) %in% c("E5.5", "E6.5")]) > 
    mean(x[grep("E3.5|blastocyst|Blastocyst",colnames(exp_scale),value = T)])
  
  return(high_exp_condition )
})

# Filter marker genes with higher expression in post-implantation stages and lower in blastocyst
filtered_markers <- post_implantation_genes[high_exp_post_implantation]

### b. plot heatmap ---- 
if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=8, height=10, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "filtered_post-implantation stage-specific genes_heatmap.pdf"), width=10, height=15, onefile = F)
}
print(
  Heatmap(as.matrix(exp_scale[filtered_markers,]),
          column_split = meta[ident_by],
          show_row_names = F,
          show_column_names = T,  # Show column names
          show_column_dend = T,  # Show column dendrogram
          column_title = "Post-implantation stage-specific genes expression",
          use_raster = F,
          col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
          heatmap_legend_param = list(title='')))
dev.off()

## 8) Similarity ----
### a. train data (blastocyst) ----
train <- subset(dataset, sample_info %in% as.character(
  grep("E3.5|E4.5|E5.5|E6.5|blastocyst",dataset$sample_info, value = T)))

RhpcBLASctl::blas_set_num_threads(10)
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data (all) ----
test <- dataset

### c. similarity ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

if (ident_by == "sample_info") {
  train$simi_info <- train$cell_info
  test$simi_info <- gsub("*-(\\d+)$", "", test$sample_info)
  
  Test_similarity =
    glm.predict(
      train.data = as.matrix(GetAssayData(train, assay = "RNA", slot = "data")),
      train.group = as.character(train$simi_info),
      # genes.used = rownames(train),
      genes.used = VariableFeatures(train),
      downsample = T,
      test.data = as.matrix(GetAssayData(test, assay = "RNA", slot = "data")),
      test.group = as.character(test$simi_info),
      alpha = 0.99,
      nfolds = 10,
      seed = 123)
  
  heatmap_mat = Test_similarity$heatmap@matrix[, c(
    "E3.5", "E4.5", "E5.5", "E6.5",
    
    "In house-blastocyst_1", "In house-blastocyst_2", "In house-blastocyst_3", 
    "In house-blastocyst_4", "In house-blastocyst_5", "In house-blastocyst_6", 
    
    "public Blastocyst-SSII_1", "public Blastocyst-SSII_2", "public Blastocyst-SSII_3", 
    
    "EPSC-public blastoid", "EPSC-public blastoid-ssII_1", "EPSC-public blastoid-ssII_2", 
    "EPSC-public blastoid-ssII_3", "EPSC-public blastoid-ssII_4", "EPSC-public blastoid-ssII_5", 
    "EPSC-public blastoid-ssII_6", "EPSC-public blastoid-ssII_7", "EPSC-public blastoid-ssII_8", 
    "2MYCP-Blastoid", "2MYCP-Blastoid-bullk_1", #"2MYCP-Blastoid-bullk_2", 
    "EPT-blastoid_1", "EPT-blastoid_2", "EPT-blastoid_3", 
    "EPT-blastoid_4", "EPT-blastoid_5", "ET-blastoid", 
    "TBLC-blastoid_1", "TBLC-blastoid_2", "TPSC-blastoid",
    
    "In house-blastoid_1", "In house-blastoid_2", "In house-blastoid_3", 
    "In house-blastoid_4", "In house-blastoid_5", "In house-blastoid_6"
  )]
}

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=7, height=3.5, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "similarity.pdf"), width=10, height=3, onefile = F)
}
print(  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                           color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                           cluster_rows = F))
dev.off()

## 9) Marker genes expr heatmap ----
### a. load marker list ----
marker <- openxlsx::read.xlsx('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/20240221/data/SSII-seq sample 信息及各个lineage marker.xlsx',
                              sheet = 2, colNames = F)
colnames(marker) <- c('type', 'gene')
rownames(marker) <- marker$gene

### b. Marker Expression Heatmap - Scale Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=8, height=8, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_ScaleData.pdf"), width=10, height=8, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

### c. Marker Expression Heatmap - Data ----
exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "data")
exp_scale <- exp_data[intersect(rownames(exp_data), marker$gene), ]

meta <- pb_seurat@meta.data
#exp_scale <- t(scale(t(exp_data), center = T, scale = T))

if (ident_by == "cell_info") {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=8, height=9, onefile = F)
} else {
  pdf(paste0(file_names_prefix, "marker_expr_heatmap_Data.pdf"), width=10, height=10, onefile = F)
}
if (nrow(as.matrix(exp_scale) %>% t) > 1){
  print(Heatmap(as.matrix(exp_scale),
                row_split = marker[rownames(exp_scale),]$type,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
} else {
  exp_scale <- as.matrix(exp_scale) %>% t
  rownames(exp_scale) <- intersect(rownames(exp_data), marker$gene)
  
  print(Heatmap(exp_scale,
                column_split = meta[ident_by],
                show_column_names = TRUE,  # Show column names
                show_column_dend = TRUE,  # Show column dendrogram
                column_title = "Markers Expression Heatmap",
                col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
                heatmap_legend_param = list(title='')))
}
dev.off()

## 10) save data ----
saveRDS(pb_seurat, paste0(file_names_prefix, "data.rds"))
write.csv(markers_post_implantation, paste0(file_names_prefix, "DEG.csv"), quote = F, row.names = T)
