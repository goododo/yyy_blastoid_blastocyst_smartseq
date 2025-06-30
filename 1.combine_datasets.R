# Title: Blastoid Blastocyst Analysis —— Combine datasets & Batch correction
# Author: Gaozy
# Time: 2025-03-28

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

# 2. load data ----
if (dir.exists("250326") == F) dir.create('250326')

## 1) subset NH0221 & NH0416 ----
## NH0221
nh0221 <- readRDS('rawdata/NH0221.rds')
nh0221$cell_type <- capitalize(nh0221$cell_type)

nh0221 <- subset(nh0221, subset = nFeature_RNA >= 300)

## NH0416
nh0416 <- readRDS('rawdata/NH0416.rds')
nh0416$cell_type <- paste0("0416_",1:10)

nh0416 <- subset(nh0416, subset = nFeature_RNA >= 300)

## contain samples close to "blastocyst"
nh0221_blastoid <- subset(nh0221, cell_type %in% "16#pick 5#-10")
nh0221_blastocyst <- subset(nh0221, cell_type %in% "Blastocyst")
nh0416 <- subset(nh0416, cell_type %in% c("0416_2", "0416_4", "0416_5", "0416_6", "0416_10"))

## add cell info
nh0221_blastoid$cell_info <- "In house-blastoid"
nh0221_blastocyst$cell_info <- "In house-blastocyst"
nh0416$cell_info <- "In house-blastoid"

## add cell info + sample
nh0221_blastoid$sample_info <- paste0("In house-blastoid","_",1)
nh0221_blastocyst$sample_info <- paste0("In house-blastocyst","_",1:6)
nh0416$sample_info <- paste0("In house-blastoid","_",2:6)

## 2) GSM4026211 ----
data <- Read10X(data.dir = "rawdata/GSE135701/GSM4026211/")
#assay_data <- CreateAssayObject(counts = data, min.cells = 3, min.features = 300) # for Seurat V5
GSM4026211 <- CreateSeuratObject(counts = data, project = "GSM4026211", assay = "RNA",
                                 min.cells = 0, min.features = 300)
GSM4026211$cell_info <- "EPSC-public blastoid"
GSM4026211$sample_info <- "EPSC-public blastoid"
GSM4026211[["RNA"]]

## 3) GSE197779 ----
data <- read.delim("rawdata/GSE197779_SC-TBLCs-Blastoids.matrix.txt",
                   header = TRUE,
                   check.names = F,
                   row.names = 1)
#assay_data <- CreateAssayObject(counts = Matrix(as.matrix(data), sparse = TRUE), min.cells = 3, min.features = 300)
GSE197779 <-CreateSeuratObject(counts = data, assay = "RNA",
                               project = "GSE197779", min.cells = 0, min.features = 300)
GSE197779$cell_info <- "TBLC-blastoid"
GSE197779$sample_info <- ifelse (sub("^[A-Za-z]+-(\\d+)$", "\\1", colnames(GSE197779)) == 0,
                                 paste0("TBLC-blastoid","_1"), paste0("TBLC-blastoid","_2"))
GSE197779[["RNA"]]

## 4) GSM6070537-TPS-blastoid ----
data <- Read10X(data.dir = "rawdata/GSM6070537-TPS-blastoid/")
#assay_data <- CreateAssayObject(counts = data, min.cells = 3, min.features = 300)
GSM6070537 <- CreateSeuratObject(counts = data, assay = "RNA",
                                 project = "GSM6070537", min.cells = 0, min.features = 300)
GSM6070537$cell_info <- "TPSC-blastoid"
GSM6070537$sample_info <- "TPSC-blastoid"
GSM6070537[["RNA"]]

## 5) GSM3940220 ----
data_list <- lapply(list.files('rawdata/GSM3940220/'), function(x){
  each_dir <- paste0('rawdata/GSM3940220/', x)
  each_data <- Read10X(data.dir = each_dir)
  #each_assay <- CreateAssayObject(counts = each_data, min.cells = 3, min.features = 300)
  each_obj <- CreateSeuratObject(counts = each_data, assay = "RNA",
                                 project = "GSM3940220", min.cells = 0, min.features = 300)
})
GSM3940220 <- merge(data_list[[1]],
                    c(data_list[[2]], data_list[[3]], data_list[[4]], data_list[[5]]))
GSM3940220$cell_info <- "EPT-blastoid"
GSM3940220$sample_info <- paste0("EPT-blastoid_",sub("^[A-Za-z]+_(\\d+)$", "\\1", colnames(GSM3940220)))
GSM3940220[["RNA"]]

## 6) GSM2916885 ----
data <- read.delim("rawdata/GSM2916885_NR-BLAS1_H77FJBGX3_S3.coutt.csv",
                   header = TRUE,
                   check.names = F,
                   row.names = 1)

data$Gene <- gsub("__chr[0-9]+", "", rownames(data))
data <- data %>% group_by(Gene) %>%
  dplyr::summarise(across(everything(), ~mean(. , na.rm = TRUE))) %>% 
  tibble::column_to_rownames("Gene")

#assay_data <- CreateAssayObject(counts = Matrix(as.matrix(data), sparse = TRUE), min.cells = 3, min.features = 300)
GSM2916885 <-CreateSeuratObject(counts = data, assay = "RNA",
                                project = "GSM2916885", min.cells = 0, min.features = 300)
GSM2916885$cell_info <- "ET-blastoid"
GSM2916885$sample_info <- "ET-blastoid"
GSM2916885[["RNA"]]

## 7) GSM7798459 ----
data <- read.delim("rawdata/GSM7798459/GSM7798459_mTBLC-2MYCP-blastoid_count.txt",
                   header = TRUE,
                   check.names = F)
#assay_data <- CreateAssayObject(counts = Matrix(as.matrix(data), sparse = TRUE), min.cells = 3, min.features = 300)
GSM7798459 <- CreateSeuratObject(counts = data, assay = "RNA",
                                 project = "GSM7798459", min.cells = 0, min.features = 300)
GSM7798459$cell_info <- "2MYCP-Blastoid"
GSM7798459$sample_info <- "2MYCP-Blastoid"
GSM7798459[["RNA"]]

## 8) GSE243925 ----
GSE243925_expr <- read.delim("rawdata/GSE243921_Blastoid.csv",
                             header = TRUE,
                             check.names = F,
                             sep = ',') %>%
  group_by(Gene) %>%
  dplyr::summarise(across(everything(), ~mean(. , na.rm = TRUE))) %>% 
  tibble::column_to_rownames("Gene")
#assay_data <- CreateAssayObject(counts = Matrix(as.matrix(GSE243925_expr), sparse = TRUE), min.cells = 1, min.features = 300)
GSE243925 <- CreateSeuratObject(counts = GSE243925_expr, assay = "RNA",
                                project = "GSE243925", min.cells = 0, min.features = 300)
GSE243925$cell_info <- "2MYCP-Blastoid-bullk"
GSE243925$sample_info <- paste0("2MYCP-Blastoid-bullk_",sub(".*-(\\d+)$", "\\1", colnames(GSE243925)))
GSE243925[["RNA"]]

## 9) GSE135289 ----
GSE135289 <- readRDS("rawdata/GSE135289.rds")
GSE135289$cell_info <- "EPSC-public blastoid-ssII"
GSE135289$sample_info <- paste0("EPSC-public blastoid-ssII_", 1:8)
GSE135289[["RNA"]]

GSE135289 <- subset(GSE135289, subset = nFeature_RNA >= 300)

## 10) GSE87504 ----
GSE87504 <- readRDS("rawdata/GSE87504.rds")
GSE87504$cell_info <- "public Blastocyst-SSII"
GSE87504$sample_info <- paste0("public Blastocyst-SSII_", rep(1:4,each=2))

GSE87504 <- subset(GSE87504, subset = nFeature_RNA >= 300)
GSE87504[["RNA"]]

## 11) stem cells ----
stemcell <- readRDS('/home/fengyan02/celloracle_data/scRNAData/TrajectoryScRNACelloracle/stem_cells_log1pData.rds')
stemcell <- subset(stemcell, cell_type2 %in% c("E3.5", "E4.5" ,"E5.5" ,"E6.5"))
stemcell$cell_info <- stemcell$cell_type2
stemcell$sample_info <- stemcell$cell_type2

stemcell <- subset(stemcell, subset = nFeature_RNA >= 300)

## 12) marker genes ----
marker <- openxlsx::read.xlsx('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/20240221/data/SSII-seq sample 信息及各个lineage marker.xlsx',
                              sheet = 2, colNames = F)
colnames(marker) <- c('type', 'gene')
rownames(marker) <- marker$gene

# 4. merge data ----
merge_data <- merge(nh0221_blastoid,
                    c(nh0221_blastocyst, nh0416, GSE135289, GSE197779, GSE243925)) %>% 
  merge(.,c(GSE87504,GSM2916885, GSM3940220, GSM4026211)) %>% 
  merge(.,c(GSM6070537, GSM7798459, stemcell))

saveRDS(merge_data, "250326/raw_merge_data.rds")

# 5. batch correction ----
method <- "ComBat_seq"
### a. batch correction ----
library(edgeR, lib.loc = "/home/qcao02/R/x86_64-conda-linux-gnu-library/4.4")
table(merge_data$orig.ident, merge_data$cell_info) 

batch_list <- SplitObject(merge_data, split.by = "orig.ident")
#### Normalize for each dataset
for(i in 1:length(batch_list)){
  batch_list[[i]] <- NormalizeData(batch_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

combined_obj <- merge(batch_list[[1]], y=batch_list[2:length(batch_list)], add.cell.ids = names(batch_list), project = "Combined")
expr_mat <- as.matrix(combined_obj@assays$RNA@data)
batch_info <- combined_obj@meta.data$orig.ident
expr_corrected <- ComBat_seq(counts = expm1(expr_mat), batch = batch_info, group = NULL)
expr_corrected_log <- log1p(expr_corrected)

# re-create obj
corrected_seurat_combat <- CreateSeuratObject(counts = expr_corrected_log, meta.data = combined_obj@meta.data)

### b. PCA ----
Idents(corrected_seurat_combat) <- "cell_info"
corrected_seurat_combat <- FindVariableFeatures(corrected_seurat_combat, selection.method = "vst", nfeatures = 2000)
corrected_seurat_combat <- ScaleData(corrected_seurat_combat, features = rownames(corrected_seurat_combat)) %>%
  RunPCA(npcs=50, features = rownames(corrected_seurat_combat))

### choose proper PCs
corrected_seurat_combat <- JackStraw(corrected_seurat_combat, num.replicate = 50, dims = 50)
corrected_seurat_combat <- ScoreJackStraw(corrected_seurat_combat, dims = 1:50) 

p1 <- JackStrawPlot(corrected_seurat_combat, dims = 1:50)
p2 <- ElbowPlot(corrected_seurat_combat, ndims = 50)

pdf(paste0("250326/corrected_pca_elbow_", method, ".pdf"), width = 12, height = 8, onefile = F)
p1 / p2
dev.off()

corrected_seurat_combat <- ScaleData(corrected_seurat_combat, features = rownames(corrected_seurat_combat)) %>%
  RunPCA(npcs=30, features = rownames(corrected_seurat_combat))

corrected_seurat_combat$cell_info <- factor(corrected_seurat_combat$cell_info, levels = unique(corrected_seurat_combat$cell_info))

## Get PCA coordinates and metadata
pca_var <- corrected_seurat_combat[["pca"]]@stdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)

p <- DimPlot( corrected_seurat_combat, reduction = "pca", dims = c(1, 2), group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 1) +
  labs( x = paste0("PC1 (", pca_var_perc[1], "%)"),
        y = paste0("PC2 (", pca_var_perc[2], "%)") )

pdf(paste0("250326/corrected_pca_", method, ".pdf"), width=8, height=6, onefile = F)
p
dev.off()

### c. UMAP & tSNE ----
corrected_seurat_combat <- FindNeighbors(corrected_seurat_combat, dims = 1:30)
corrected_seurat_combat <- FindClusters(corrected_seurat_combat, resolution = 0.5)

corrected_seurat_combat <- RunTSNE(corrected_seurat_combat, dims = 1:30 )
corrected_seurat_combat <- RunUMAP(corrected_seurat_combat, dims = 1:30 )

# UMAP plot
p1 <- DimPlot(corrected_seurat_combat, reduction = "umap", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("UMAP") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# tSNE plot
p2 <- DimPlot(corrected_seurat_combat, reduction = "tsne", group.by = "cell_info", label = TRUE, repel = TRUE, pt.size = 1)  +
  ggtitle("tSNE") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Combine plots
final_plot <- (p1 | p2) + plot_layout(guides = "collect")

pdf(paste0("250326/corrected_umap_tsne_", method, ".pdf"), width=25, height=12, onefile = F)
final_plot
dev.off()

### d. marker genes heatmap ----
Idents(corrected_seurat_combat) <- "cell_info"
options(future.globals.maxSize= 4000*1024^2)

pb_seurat <- Seurat:::PseudobulkExpression(corrected_seurat_combat, return.seurat = T)
pb_seurat@meta.data$cell_info <- rownames(pb_seurat@meta.data)

RhpcBLASctl::blas_set_num_threads(8)
Idents(pb_seurat) <- "cell_info"
pb_seurat <- NormalizeData(pb_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pb_seurat <- FindVariableFeatures(pb_seurat, selection.method = "vst", nfeatures = 2000)
pb_seurat <- ScaleData(pb_seurat, features = VariableFeatures(pb_seurat)) %>%
  RunPCA(npcs=10, features = VariableFeatures(pb_seurat))

pb_seurat$cell_info <- factor(pb_seurat$cell_info, levels = unique(pb_seurat$cell_info))

exp_data <- GetAssayData(pb_seurat, assay = "RNA", slot = "scale.data")
exp_data <- exp_data[intersect(rownames(exp_data), marker$gene), ]
colnames(exp_data) <- pb_seurat$cell_info

meta <- pb_seurat@meta.data
exp_scale <- exp_data

pdf(paste0("250326/corrected_marker_heatmap_", method, ".pdf"), width = 8, height = 8, onefile = F)
Heatmap(as.matrix(exp_scale),
        row_split = marker[rownames(exp_scale),]$type,
        column_split = meta$cell_info,
        show_column_names = TRUE,  # Show column names
        show_column_dend = TRUE,  # Show column dendrogram
        column_title = "Markers Expression Heatmap",
        col = colorRampPalette(rev(brewer.pal(9,"Spectral")))(100),
        heatmap_legend_param = list(title=''))
dev.off()

### e. save data ----
saveRDS(corrected_seurat_combat, "250326/combat.rds")
