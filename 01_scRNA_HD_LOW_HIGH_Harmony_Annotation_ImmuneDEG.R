############################################################
## scRNA-seq analysis: HD / LOW / HIGH
## Datasets: GSE192483 and GSE135893
############################################################

## =========================
## 0. Load packages
## =========================
library(Seurat)
library(harmony)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

setwd("D:/atraining/scGSE192483/GSE192483_RAW")


## =========================
## 1. Read PTB lesion samples
## =========================
HIGH1 <- Read10X("./HIGH1/")
LOW1  <- Read10X("./LOW1/")
HIGH2 <- Read10X("./HIGH2/")
LOW2  <- Read10X("./LOW2/")
HIGH3 <- Read10X("./HIGH3/")
LOW3  <- Read10X("./LOW3/")
HIGH4 <- Read10X("./HIGH4/")
LOW4  <- Read10X("./LOW4/")
HIGH5 <- Read10X("./HIGH5/")
LOW5  <- Read10X("./LOW5/")

HIGH1 <- CreateSeuratObject(HIGH1, project = "HIGH1")
LOW1  <- CreateSeuratObject(LOW1,  project = "LOW1")
HIGH2 <- CreateSeuratObject(HIGH2, project = "HIGH2")
LOW2  <- CreateSeuratObject(LOW2,  project = "LOW2")
HIGH3 <- CreateSeuratObject(HIGH3, project = "HIGH3")
LOW3  <- CreateSeuratObject(LOW3,  project = "LOW3")
HIGH4 <- CreateSeuratObject(HIGH4, project = "HIGH4")
LOW4  <- CreateSeuratObject(LOW4,  project = "LOW4")
HIGH5 <- CreateSeuratObject(HIGH5, project = "HIGH5")
LOW5  <- CreateSeuratObject(LOW5,  project = "LOW5")


## =========================
## 2. Read normal lung control cells from GSE135893
## =========================
setwd("D:/atraining/scGSE192483")

CT1 <- Read10X(data.dir = "./GSE135893/", gene.column = 1)
CT1 <- CreateSeuratObject(CT1, project = "CT1")

normal_ids <- c("HD65", "HD66", "HD67")
CT1_normal <- subset(CT1, subset = orig.ident %in% normal_ids)


## =========================
## 3. Merge datasets and define groups
## =========================
scRNA <- merge(
  CT1_normal,
  y = c(
    HIGH1, HIGH2, HIGH3, HIGH4, HIGH5,
    LOW1, LOW2, LOW3, LOW4, LOW5
  )
)

orig <- scRNA$orig.ident

scRNA$tissue_type <- ifelse(
  orig %in% c("HD65", "HD66", "HD67"), "HD",
  ifelse(grepl("^HIGH", orig), "HIGH",
         ifelse(grepl("^LOW", orig), "LOW", "Other"))
)

table(scRNA$tissue_type)


## =========================
## 4. Quality control
## =========================
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

VlnPlot(
  scRNA,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

scRNA <- subset(
  scRNA,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 10
)

VlnPlot(
  scRNA,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)


## =========================
## 5. Normalization, feature selection, PCA
## =========================
scRNA <- NormalizeData(
  scRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

scRNA <- FindVariableFeatures(
  scRNA,
  selection.method = "vst",
  nfeatures = 2000
)

scale.genes <- VariableFeatures(scRNA)

scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(
  scRNA,
  features = VariableFeatures(scRNA)
)

DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")


## =========================
## 6. Harmony batch correction
## =========================
setwd("D:/atraining/scGSE192483/GSE192483_RAW")

scRNA$orig.ident <- dplyr::recode(
  scRNA$orig.ident,
  "HD65" = "HD1",
  "HD66" = "HD2",
  "HD67" = "HD3"
)

saveRDS(scRNA, file = "scRNA_merge_preharmony.RDS")

scRNA <- readRDS("scRNA_merge_preharmony.RDS")

scRNA_harmony <- RunHarmony(
  scRNA,
  group.by.vars = "orig.ident"
)

DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_harmony, reduction = "harmony")


## =========================
## 7. Clustering and dimensionality reduction
## =========================
scRNA <- FindNeighbors(
  scRNA_harmony,
  dims = 1:15,
  reduction = "harmony"
)

scRNA <- FindClusters(scRNA)

scRNA <- RunUMAP(
  scRNA,
  dims = 1:15,
  reduction = "harmony"
)

scRNA <- RunTSNE(
  scRNA,
  dims = 1:15,
  reduction = "harmony",
  perplexity = 30
)

DimPlot(scRNA, reduction = "umap", group.by = "tissue_type", label = FALSE)
DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", label = FALSE)

saveRDS(scRNA, file = "scRNA_HD_HIGH_LOW.RDS")


## =========================
## 8. Final cell annotation
## =========================
scRNA <- readRDS("scRNA_HD_HIGH_LOW.RDS")

Idents(scRNA) <- "seurat_clusters"

## Exclude cluster 28 according to final analysis
scRNA <- subset(scRNA, idents = setdiff(levels(scRNA), "28"))

scRNA <- RenameIdents(
  scRNA,
  "0"  = "T cell",
  "1"  = "NK cell",
  "2"  = "T cell",
  "3"  = "T cell",
  "4"  = "Macrophage",
  "5"  = "Macrophage",
  "6"  = "Monocyte",
  "7"  = "Monocyte",
  "8"  = "B cell",
  "9"  = "Monocyte",
  "10" = "Epithelial",
  "11" = "Endothelial",
  "12" = "Epithelial",
  "13" = "Macrophage",
  "14" = "T cell",
  "15" = "Macrophage",
  "16" = "NK cell",
  "17" = "T cell",
  "18" = "Smooth muscle",
  "19" = "Epithelial",
  "20" = "T cell",
  "21" = "DC",
  "22" = "NK cell",
  "23" = "T cell",
  "24" = "NK cell",
  "25" = "Monocyte",
  "26" = "Macrophage",
  "27" = "Macrophage"
)

scRNA$celltype <- Idents(scRNA)

table(scRNA$celltype)
table(scRNA$celltype, scRNA$tissue_type)

saveRDS(scRNA, file = "scRNA_HD_HIGH_LOW260220.RDS")

############################################################
## Differential expression analysis in lymphoid and myeloid cells
############################################################

library(Seurat)
library(dplyr)
library(tibble)
library(purrr)

setwd("D:/atraining/scGSE192483/GSE192483_RAW")

scRNA <- readRDS("scRNA_HD_HIGH_LOW260220.RDS")

## =========================
## 1. Define immune compartments
## =========================
Idents(scRNA) <- "seurat_clusters"

lymphoid_clusters <- c(
  0, 2, 3, 14, 17, 20, 23,
  1, 16, 22, 24,
  8
)

myeloid_clusters <- c(
  4, 5, 13, 15, 26, 27,
  6, 7, 9, 25,
  21
)

scRNA$immune_major <- "Other"
scRNA$immune_major[scRNA$seurat_clusters %in% lymphoid_clusters] <- "Lymphoid"
scRNA$immune_major[scRNA$seurat_clusters %in% myeloid_clusters]  <- "Myeloid"

table(scRNA$immune_major, scRNA$tissue_type)


## =========================
## 2. Subset lymphoid and myeloid cells
## =========================
sc_myeloid  <- subset(scRNA, subset = immune_major == "Myeloid")
sc_lymphoid <- subset(scRNA, subset = immune_major == "Lymphoid")

Idents(sc_myeloid)  <- "tissue_type"
Idents(sc_lymphoid) <- "tissue_type"

sc_myeloid <- subset(
  sc_myeloid,
  idents = c("HD", "LOW", "HIGH")
)

sc_lymphoid <- subset(
  sc_lymphoid,
  idents = c("HD", "LOW", "HIGH")
)


## =========================
## 3. Define pairwise comparisons
## =========================
comparisons <- list(
  HD_vs_LOW   = c("HD", "LOW"),
  HD_vs_HIGH  = c("HD", "HIGH"),
  LOW_vs_HIGH = c("LOW", "HIGH")
)


## =========================
## 4. Function for DEG analysis
## =========================
run_deg_all <- function(obj, g1, g2, label, cellclass) {
  
  deg <- FindMarkers(
    obj,
    ident.1 = g1,
    ident.2 = g2,
    min.pct = 0.1,
    logfc.threshold = 0.1,
    test.use = "wilcox"
  )
  
  deg %>%
    rownames_to_column("Gene") %>%
    mutate(
      CellClass  = cellclass,
      Comparison = label
    )
}


## =========================
## 5. Run DEG analysis
## =========================
deg_myeloid_all <- imap_dfr(
  comparisons,
  ~ run_deg_all(
    obj = sc_myeloid,
    g1 = .x[1],
    g2 = .x[2],
    label = .y,
    cellclass = "Myeloid"
  )
)

deg_lymphoid_all <- imap_dfr(
  comparisons,
  ~ run_deg_all(
    obj = sc_lymphoid,
    g1 = .x[1],
    g2 = .x[2],
    label = .y,
    cellclass = "Lymphoid"
  )
)

deg_all <- bind_rows(deg_myeloid_all, deg_lymphoid_all)


## =========================
## 6. Define significant DEGs
## =========================
deg_sig_all <- deg_all %>%
  filter(p_val_adj <= 0.05)

sig_genes_all <- unique(deg_sig_all$Gene)

cat("Number of significant myeloid DEGs:\n")
print(length(unique(deg_myeloid_all$Gene[deg_myeloid_all$p_val_adj <= 0.05])))

cat("Number of significant lymphoid DEGs:\n")
print(length(unique(deg_lymphoid_all$Gene[deg_lymphoid_all$p_val_adj <= 0.05])))

cat("Number of significant DEGs across immune compartments:\n")
print(length(sig_genes_all))


## =========================
## 7. Save results
## =========================
write.csv(
  deg_myeloid_all,
  "deg_myeloid_all_final.csv",
  row.names = FALSE
)

write.csv(
  deg_lymphoid_all,
  "deg_lymphoid_all_final.csv",
  row.names = FALSE
)

write.csv(
  deg_all,
  "deg_all_final.csv",
  row.names = FALSE
)

write.csv(
  deg_sig_all,
  "deg_sig_all_final.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = sig_genes_all),
  "sig_genes_all_final.csv",
  row.names = FALSE
)