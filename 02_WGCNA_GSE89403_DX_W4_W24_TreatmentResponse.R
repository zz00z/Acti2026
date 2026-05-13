############################################################
## 02_WGCNA_GSE89403_DX_W4_W24_TreatmentResponse.R
##
## Weighted gene co-expression network analysis of the
## GSE89403 blood transcriptomic cohort to identify treatment-
## associated gene modules between diagnosis and follow-up
## time points, including DX vs week 4 and DX vs week 24.
##
## Main steps:
## 1. Read and preprocess expression matrices
## 2. Select the top 50% most variable genes by MAD
## 3. Construct WGCNA networks using unsigned Pearson correlation
## 4. Identify module–trait correlations with treatment time points
## 5. Select top five treatment-associated modules
## 6. Visualize top paired differential genes within selected modules
## 7. Intersect DX–W4 and DX–W24 module genes
############################################################


## =========================================================
## 0. Load packages
## =========================================================

library(WGCNA)
library(dplyr)
library(readr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()


## =========================================================
## 1. Set working directory and global parameters
## =========================================================

setwd("D:/atraining/blood/GSE89403/WGCNA")

network_type <- "unsigned"
cor_type <- "pearson"

corFnc <- ifelse(cor_type == "pearson", cor, bicor)
maxPOutliers <- ifelse(cor_type == "pearson", 1, 0.05)
robustY <- ifelse(cor_type == "pearson", TRUE, FALSE)

powers <- c(1:10, seq(from = 12, to = 20, by = 2))


## =========================================================
## 2. Define reusable WGCNA function
## =========================================================

run_wgcna_treatment <- function(
    expr_file,
    followup_suffix,
    followup_label,
    out_dir,
    network_type = "unsigned",
    cor_type = "pearson",
    powers = c(1:10, seq(from = 12, to = 20, by = 2)),
    mad_quantile = 0.5,
    min_module_size = 30,
    merge_cut_height = 0.25,
    top_module_n = 5,
    filter_module_p = TRUE
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("============================================================")
  message("Running WGCNA for: ", followup_label)
  message("Expression file: ", expr_file)
  message("Output directory: ", out_dir)
  message("============================================================")
  
  ## ---------------------------------------------------------
  ## 2.1 Read expression matrix
  ## ---------------------------------------------------------
  
  expr_raw <- read.csv(
    expr_file,
    sep = ",",
    row.names = 1,
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
  
  expr_raw <- expr_raw[!is.na(rownames(expr_raw)) & rownames(expr_raw) != "", ]
  
  ## Collapse duplicated genes by mean expression
  if (any(duplicated(rownames(expr_raw)))) {
    expr_raw <- aggregate(
      expr_raw,
      by = list(rownames(expr_raw)),
      FUN = mean,
      na.rm = TRUE
    )
    rownames(expr_raw) <- expr_raw$Group.1
    expr_raw$Group.1 <- NULL
  }
  
  message("Raw expression matrix dimension: ")
  print(dim(expr_raw))
  
  
  ## ---------------------------------------------------------
  ## 2.2 Select top 50% most variable genes by MAD
  ## ---------------------------------------------------------
  
  gene_mad <- apply(expr_raw, 1, mad)
  mad_cutoff <- quantile(gene_mad, probs = mad_quantile, na.rm = TRUE)
  
  expr_var <- expr_raw[gene_mad > mad_cutoff, ]
  
  message("MAD cutoff: ", mad_cutoff)
  message("Genes retained after MAD filtering: ", nrow(expr_var))
  message("Genes removed after MAD filtering: ", sum(gene_mad <= mad_cutoff))
  
  
  ## ---------------------------------------------------------
  ## 2.3 Prepare WGCNA input: samples as rows, genes as columns
  ## ---------------------------------------------------------
  
  dataExpr <- as.data.frame(t(expr_var))
  
  gsg <- goodSamplesGenes(dataExpr, verbose = 3)
  
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      message("Removing genes: ",
              paste(names(dataExpr)[!gsg$goodGenes], collapse = ", "))
    }
    if (sum(!gsg$goodSamples) > 0) {
      message("Removing samples: ",
              paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", "))
    }
    
    dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  nGenes <- ncol(dataExpr)
  nSamples <- nrow(dataExpr)
  
  message("Final WGCNA input dimension: ")
  print(dim(dataExpr))
  message("Number of genes: ", nGenes)
  message("Number of samples: ", nSamples)
  
  
  ## ---------------------------------------------------------
  ## 2.4 Sample clustering
  ## ---------------------------------------------------------
  
  pdf(file.path(out_dir, paste0("Sample_clustering_", followup_label, ".pdf")),
      width = 10, height = 6)
  sampleTree <- hclust(dist(dataExpr), method = "average")
  plot(
    sampleTree,
    main = paste0("Sample clustering to detect outliers: ", followup_label),
    sub = "",
    xlab = ""
  )
  dev.off()
  
  
  ## ---------------------------------------------------------
  ## 2.5 Select soft-thresholding power
  ## ---------------------------------------------------------
  
  sft <- pickSoftThreshold(
    dataExpr,
    powerVector = powers,
    networkType = network_type,
    verbose = 5
  )
  
  pdf(file.path(out_dir, paste0("Soft_threshold_", followup_label, ".pdf")),
      width = 10, height = 5)
  par(mfrow = c(1, 2))
  
  cex1 <- 0.9
  
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    type = "n",
    main = "Scale independence"
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  abline(h = 0.8, col = "red")
  
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = "Mean connectivity"
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  dev.off()
  
  power <- sft$powerEstimate
  
  if (is.na(power)) {
    power <- ifelse(
      nSamples < 20,
      ifelse(network_type == "unsigned", 9, 18),
      ifelse(
        nSamples < 30,
        ifelse(network_type == "unsigned", 8, 16),
        ifelse(
          nSamples < 40,
          ifelse(network_type == "unsigned", 7, 14),
          ifelse(network_type == "unsigned", 6, 12)
        )
      )
    )
  }
  
  message("Final soft-thresholding power: ", power)
  
  
  ## ---------------------------------------------------------
  ## 2.6 Construct WGCNA network
  ## ---------------------------------------------------------
  
  net <- blockwiseModules(
    dataExpr,
    power = power,
    maxBlockSize = nGenes,
    TOMType = network_type,
    minModuleSize = min_module_size,
    reassignThreshold = 0,
    mergeCutHeight = merge_cut_height,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    corType = cor_type,
    maxPOutliers = maxPOutliers,
    loadTOMs = TRUE,
    saveTOMFileBase = file.path(out_dir, paste0(basename(expr_file), ".tom")),
    verbose = 3
  )
  
  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  
  module_table <- table(moduleColors)
  write.csv(
    as.data.frame(module_table),
    file.path(out_dir, paste0("Module_size_", followup_label, ".csv")),
    row.names = FALSE
  )
  
  pdf(file.path(out_dir, paste0("Gene_dendrogram_modules_", followup_label, ".pdf")),
      width = 12, height = 6)
  plotDendroAndColors(
    net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  
  ## ---------------------------------------------------------
  ## 2.7 Module eigengene network
  ## ---------------------------------------------------------
  
  MEs <- net$MEs
  
  MEs_col <- MEs
  colnames(MEs_col) <- paste0(
    "ME",
    labels2colors(as.numeric(str_replace_all(colnames(MEs), "ME", "")))
  )
  MEs_col <- orderMEs(MEs_col)
  
  pdf(file.path(out_dir, paste0("Eigengene_network_", followup_label, ".pdf")),
      width = 10, height = 8)
  plotEigengeneNetworks(
    MEs_col,
    "Eigengene adjacency heatmap",
    marDendro = c(3, 3, 2, 4),
    marHeatmap = c(3, 4, 2, 2),
    plotDendrograms = TRUE,
    xLabelsAngle = 90
  )
  dev.off()
  
  
  ## ---------------------------------------------------------
  ## 2.8 Build trait matrix
  ## ---------------------------------------------------------
  
  sample_names <- rownames(dataExpr)
  
  timepoints <- sapply(sample_names, function(sample) {
    if (grepl("_DX$", sample)) {
      return("DX")
    } else if (grepl(paste0("_", followup_suffix, "$"), sample)) {
      return(followup_label)
    } else {
      return("unknown")
    }
  })
  
  traitData <- data.frame(
    row.names = sample_names,
    DX = as.numeric(timepoints == "DX"),
    Followup = as.numeric(timepoints == followup_label)
  )
  
  colnames(traitData) <- c("DX", followup_label)
  traitData <- traitData[match(sample_names, rownames(traitData)), ]
  
  write.csv(
    traitData,
    file.path(out_dir, paste0("Trait_matrix_", followup_label, ".csv"))
  )
  
  
  ## ---------------------------------------------------------
  ## 2.9 Module–trait correlation
  ## ---------------------------------------------------------
  
  if (cor_type == "pearson") {
    modTraitCor <- cor(MEs_col, traitData, use = "p")
    modTraitP <- corPvalueStudent(modTraitCor, nSamples)
  } else {
    modTraitCorP <- bicorAndPvalue(MEs_col, traitData, robustY = robustY)
    modTraitCor <- modTraitCorP$bicor
    modTraitP <- modTraitCorP$p
  }
  
  write.csv(
    modTraitCor,
    file.path(out_dir, paste0("Module_trait_correlation_", followup_label, ".csv"))
  )
  
  write.csv(
    modTraitP,
    file.path(out_dir, paste0("Module_trait_pvalue_", followup_label, ".csv"))
  )
  
  textMatrix <- paste(
    signif(modTraitCor, 2),
    "\n(",
    signif(modTraitP, 1),
    ")",
    sep = ""
  )
  dim(textMatrix) <- dim(modTraitCor)
  
  pdf(file.path(out_dir, paste0("Module_trait_heatmap_", followup_label, ".pdf")),
      width = 8, height = 10)
  par(mar = c(6, 8, 4, 2))
  labeledHeatmap(
    Matrix = modTraitCor,
    xLabels = colnames(traitData),
    yLabels = colnames(MEs_col),
    cex.lab = 0.8,
    ySymbols = colnames(MEs_col),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = paste0("Module-trait relationships: ", followup_label)
  )
  dev.off()
  
  
  ## ---------------------------------------------------------
  ## 2.10 Save gene-module information
  ## ---------------------------------------------------------
  
  geneInfo <- data.frame(
    gene = colnames(dataExpr),
    moduleColor = moduleColors,
    moduleLabel = net$colors,
    stringsAsFactors = FALSE
  )
  
  geneInfo$gene <- gsub('"', "", geneInfo$gene)
  
  write.csv(
    geneInfo,
    file.path(out_dir, paste0("All_gene_module_info_", followup_label, ".csv")),
    row.names = FALSE
  )
  
  
  ## ---------------------------------------------------------
  ## 2.11 Select top treatment-associated modules
  ## ---------------------------------------------------------
  
  module_max_cor <- apply(abs(modTraitCor), 1, max)
  module_min_p <- apply(modTraitP, 1, min)
  
  module_rank_df <- data.frame(
    Module = rownames(modTraitCor),
    MaxAbsCor = module_max_cor,
    MinP = module_min_p,
    stringsAsFactors = FALSE
  )
  
  if (filter_module_p) {
    module_rank_df <- module_rank_df %>%
      filter(MinP < 0.05)
  }
  
  module_rank_df <- module_rank_df %>%
    arrange(desc(MaxAbsCor))
  
  top5_modules <- module_rank_df$Module[1:min(top_module_n, nrow(module_rank_df))]
  
  write.csv(
    module_rank_df,
    file.path(out_dir, paste0("Module_rank_by_trait_correlation_", followup_label, ".csv")),
    row.names = FALSE
  )
  
  message("Selected top modules:")
  print(top5_modules)
  
  
  ## ---------------------------------------------------------
  ## 2.12 Export genes from top modules
  ## ---------------------------------------------------------
  
  top5_genes <- data.frame()
  
  for (module in top5_modules) {
    
    module_color <- gsub("ME", "", module)
    module_genes <- geneInfo[geneInfo$moduleColor == module_color, ]
    
    module_genes$moduleName <- module
    module_genes$corWithDX <- modTraitCor[module, "DX"]
    module_genes$pValueDX <- modTraitP[module, "DX"]
    module_genes$corWithFollowup <- modTraitCor[module, followup_label]
    module_genes$pValueFollowup <- modTraitP[module, followup_label]
    module_genes$FollowupLabel <- followup_label
    
    top5_genes <- rbind(top5_genes, module_genes)
    
    write.table(
      module_genes$gene,
      file = file.path(out_dir, paste0("module_", module, "_genes_", followup_label, ".txt")),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
  
  write.csv(
    top5_genes,
    file.path(out_dir, paste0("Top5_modules_combined_genes_", followup_label, ".csv")),
    row.names = FALSE
  )
  
  write.table(
    top5_genes,
    file.path(out_dir, paste0("Top5_modules_combined_genes_", followup_label, ".txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  
  ## ---------------------------------------------------------
  ## 2.13 Paired differential analysis for heatmap genes
  ## ---------------------------------------------------------
  
  dx_cols <- grep("_DX$", colnames(expr_raw), value = TRUE)
  fu_cols <- grep(paste0("_", followup_suffix, "$"), colnames(expr_raw), value = TRUE)
  
  dx_ids <- sub("_DX$", "", dx_cols)
  fu_ids <- sub(paste0("_", followup_suffix, "$"), "", fu_cols)
  
  common_ids <- intersect(dx_ids, fu_ids)
  
  dx_cols <- paste0(common_ids, "_DX")
  fu_cols <- paste0(common_ids, "_", followup_suffix)
  
  expr_dx <- expr_raw[, dx_cols, drop = FALSE]
  expr_fu <- expr_raw[, fu_cols, drop = FALSE]
  
  message("Number of paired samples for ", followup_label, ": ", length(common_ids))
  
  deg <- data.frame(
    gene = rownames(expr_raw),
    pvalue = NA_real_,
    logFC = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(expr_raw))) {
    
    dx <- as.numeric(expr_dx[i, ])
    fu <- as.numeric(expr_fu[i, ])
    
    keep <- complete.cases(dx, fu)
    dx2 <- dx[keep]
    fu2 <- fu[keep]
    
    if (length(dx2) >= 2) {
      deg$pvalue[i] <- tryCatch(
        t.test(dx2, fu2, paired = TRUE)$p.value,
        error = function(e) NA_real_
      )
      
      deg$logFC[i] <- mean(fu2, na.rm = TRUE) - mean(dx2, na.rm = TRUE)
    }
  }
  
  deg$padj <- p.adjust(deg$pvalue, method = "BH")
  
  write.csv(
    deg,
    file.path(out_dir, paste0("Paired_DEG_DX_vs_", followup_label, ".csv")),
    row.names = FALSE
  )
  
  
  ## ---------------------------------------------------------
  ## 2.14 Select top 5 differential genes from each top module
  ## ---------------------------------------------------------
  
  deg2 <- merge(deg, top5_genes, by = "gene")
  deg2 <- deg2 %>% filter(!is.na(pvalue))
  
  top_heatmap_genes <- deg2 %>%
    group_by(moduleColor) %>%
    arrange(pvalue, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    arrange(moduleColor, pvalue)
  
  write.csv(
    top_heatmap_genes,
    file.path(out_dir, paste0("Top25_genes_by_module_for_heatmap_", followup_label, ".csv")),
    row.names = FALSE
  )
  
  
  ## ---------------------------------------------------------
  ## 2.15 Draw heatmap
  ## ---------------------------------------------------------
  
  gene_order <- top_heatmap_genes$gene
  
  heat_expr <- expr_raw[gene_order, c(dx_cols, fu_cols), drop = FALSE]
  rownames(heat_expr) <- gene_order
  
  heat_expr_scaled <- t(scale(t(as.matrix(heat_expr))))
  
  keep_row <- apply(heat_expr_scaled, 1, function(x) all(!is.na(x)))
  heat_expr_scaled <- heat_expr_scaled[keep_row, , drop = FALSE]
  
  top_heatmap_genes2 <- top_heatmap_genes[
    match(rownames(heat_expr_scaled), top_heatmap_genes$gene),
    ,
    drop = FALSE
  ]
  
  annotation_col <- data.frame(
    Group = c(rep("DX", length(dx_cols)), rep(followup_label, length(fu_cols))),
    row.names = c(dx_cols, fu_cols),
    stringsAsFactors = FALSE
  )
  
  annotation_row <- data.frame(
    Module = top_heatmap_genes2$moduleColor,
    row.names = top_heatmap_genes2$gene,
    stringsAsFactors = FALSE
  )
  
  module_levels <- unique(annotation_row$Module)
  
  is_valid_color <- function(x) {
    sapply(x, function(cc) {
      tryCatch({
        col2rgb(cc)
        TRUE
      }, error = function(e) FALSE)
    })
  }
  
  valid_flag <- is_valid_color(module_levels)
  
  module_color_map <- setNames(
    ifelse(
      valid_flag,
      module_levels,
      brewer.pal(max(3, min(8, length(module_levels))), "Set2")[seq_along(module_levels)]
    ),
    module_levels
  )
  
  ann_colors <- list(
    Group = setNames(c("#E64B35", "#4DBBD5"), c("DX", followup_label)),
    Module = module_color_map
  )
  
  module_table <- table(annotation_row$Module)
  gaps_row <- cumsum(as.numeric(module_table))
  gaps_row <- gaps_row[-length(gaps_row)]
  
  pdf(
    file.path(out_dir, paste0("WGCNA_top25_genes_heatmap_", followup_label, ".pdf")),
    width = 10,
    height = 8
  )
  pheatmap(
    heat_expr_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    gaps_row = gaps_row,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 9,
    fontsize_col = 10,
    border_color = NA,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
    main = paste0("Top differential genes from WGCNA modules: DX vs ", followup_label)
  )
  dev.off()
  
  
  ## ---------------------------------------------------------
  ## 2.16 Return results
  ## ---------------------------------------------------------
  
  result <- list(
    dataExpr = dataExpr,
    power = power,
    net = net,
    moduleColors = moduleColors,
    MEs_col = MEs_col,
    traitData = traitData,
    modTraitCor = modTraitCor,
    modTraitP = modTraitP,
    geneInfo = geneInfo,
    moduleRank = module_rank_df,
    topModules = top5_modules,
    topModuleGenes = top5_genes,
    pairedDEG = deg,
    heatmapGenes = top_heatmap_genes
  )
  
  saveRDS(
    result,
    file.path(out_dir, paste0("WGCNA_result_", followup_label, ".rds"))
  )
  
  return(result)
}


## =========================================================
## 3. Run WGCNA for DX vs week 4
## =========================================================

result_w4 <- run_wgcna_treatment(
  expr_file = "expr_filtered_4W_TGRV_notimetonegativity04.csv",
  followup_suffix = "week_4",
  followup_label = "week_4",
  out_dir = "WGCNA_DX_vs_week4",
  network_type = network_type,
  cor_type = cor_type,
  powers = powers,
  mad_quantile = 0.5,
  min_module_size = 30,
  merge_cut_height = 0.25,
  top_module_n = 5,
  filter_module_p = TRUE
)


## =========================================================
## 4. Run WGCNA for DX vs week 24
## =========================================================
## Replace the file name below with your actual DX–week 24
## expression matrix file.

result_w24 <- run_wgcna_treatment(
  expr_file = "expr_filtered_24W_TGRV_notimetonegativity04.csv",
  followup_suffix = "week_24",
  followup_label = "week_24",
  out_dir = "WGCNA_DX_vs_week24",
  network_type = network_type,
  cor_type = cor_type,
  powers = powers,
  mad_quantile = 0.5,
  min_module_size = 30,
  merge_cut_height = 0.25,
  top_module_n = 5,
  filter_module_p = TRUE
)


## =========================================================
## 5. Intersect genes from DX–week 4 and DX–week 24 top modules
## =========================================================

genes_w4 <- unique(result_w4$topModuleGenes$gene)
genes_w24 <- unique(result_w24$topModuleGenes$gene)

genes_overlap <- intersect(genes_w4, genes_w24)

message("Number of genes in DX–week 4 top modules: ", length(genes_w4))
message("Number of genes in DX–week 24 top modules: ", length(genes_w24))
message("Number of overlapping treatment-associated genes: ", length(genes_overlap))

write.csv(
  data.frame(Gene = genes_w4),
  "WGCNA_DX_week4_top_module_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = genes_w24),
  "WGCNA_DX_week24_top_module_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = genes_overlap),
  "WGCNA_DX_week4_week24_overlap_genes.csv",
  row.names = FALSE
)


## =========================================================
## 6. Draw Venn diagram for WGCNA gene overlap
## =========================================================

venn_df <- data.frame(
  Gene = unique(c(genes_w4, genes_w24)),
  DX_week4 = unique(c(genes_w4, genes_w24)) %in% genes_w4,
  DX_week24 = unique(c(genes_w4, genes_w24)) %in% genes_w24
)

write.csv(
  venn_df,
  "WGCNA_DX_week4_week24_venn_input.csv",
  row.names = FALSE
)

## Optional Venn plot using ggplot2 only
venn_counts <- data.frame(
  Category = c("DX_week4", "DX_week24", "Overlap"),
  Count = c(length(genes_w4), length(genes_w24), length(genes_overlap))
)

write.csv(
  venn_counts,
  "WGCNA_DX_week4_week24_venn_counts.csv",
  row.names = FALSE
)

message("WGCNA analysis completed.")