# Single-cell RNA-seq analysis and plots.
# Inputs: 10x matrices, HPCA reference and reviewed annotation map; paths in acti_config().
# Run: Rscript 1_ACTI_single_cell_analysis.R --stage=qc
# source() loads functions without starting an analysis.
# R 4.1.3 / Seurat 4.x; review the cluster annotation before analysis.

acti_config <- function(project_dir = ".", data_dir = file.path(project_dir, "data"), output_dir = file.path(project_dir,
  "results")) {
  list(project_dir = project_dir, data_dir = data_dir, output_dir = output_dir, case_dir = file.path(data_dir,
    "GSE192483_RAW"), control_dir = file.path(data_dir, "GSE135893"), hpca_rds = file.path(data_dir, "reference",
    "HumanPrimaryCellAtlasData.rds"), annotation_map = file.path(data_dir, "annotation_map_FINAL.csv"), annotated_rds = "",
    use_manuscript_annotation = FALSE, stpmg_csv = file.path(data_dir, "bulk", "STPMGs_1174_shared_genes.csv"),
    figure2_mode = "ALL")
}

acti_input_object <- function(cfg) {
  if (nzchar(cfg$annotated_rds)) {
    if (!file.exists(cfg$annotated_rds))
      stop("Annotated checkpoint not found: ", cfg$annotated_rds)
    return(cfg$annotated_rds)
  }
  candidates <- file.path(cfg$output_dir, c("03_lineage_pseudobulk/RDS/scRNA_major_annotation_DE_audited.rds",
    "02_annotation/RDS/scRNA_major_classes_annotated.rds"))
  found <- candidates[file.exists(candidates)]
  if (!length(found))
    stop("Run the annotation stage or supply --annotated-rds=... first.")
  found[1L]
}

acti_candidate_filter <- function(x) {
  required <- c("PValue", "logFC", "same_direction_pairs", "n_complete_pairs_checked")
  if (!all(required %in% names(x)))
    stop("Candidate results lack required fields.")
  x$PValue < 0.05 & abs(x$logFC) >= 0.2 & x$n_complete_pairs_checked >= 5L &
    x$same_direction_pairs/x$n_complete_pairs_checked >= 0.8
}

acti_package_version <- function(package) {
  directory <- find.package(package, quiet = TRUE)
  if (!length(directory))
    return(NA_character_)
  tryCatch(unname(read.dcf(file.path(directory, "DESCRIPTION"), fields = "Version")[1L, 1L]), error = function(e) NA_character_)
}

acti_stage_packages <- list(qc = c("Seurat", "Matrix", "SingleCellExperiment", "scDblFinder", "BiocParallel", "harmony",
  "dplyr", "tidyr", "tibble", "ggplot2"), annotation_evidence = c("Seurat", "SingleR", "SingleCellExperiment",
  "SummarizedExperiment", "MAST", "Matrix", "edgeR", "dplyr", "tidyr", "tibble", "ggplot2"), annotation_audit = c("Seurat",
  "dplyr", "tidyr", "ggplot2"), annotation = c("Seurat", "dplyr", "ggplot2"), lineage = c("Seurat", "Matrix",
  "edgeR", "dplyr", "tidyr", "tibble", "ggplot2", "statmod"), celltype = c("Seurat", "Matrix", "edgeR", "statmod"),
  harmony = c("Seurat", "Matrix", "RANN", "cluster", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork"), figure2 = c("Seurat",
    "ggplot2", "dplyr", "tidyr", "scales", "RColorBrewer", "scattermore", "ggrepel", "edgeR", "svglite", "ragg"),
  dotplot = c("Seurat", "Matrix", "dplyr", "tidyr", "ggplot2", "scales", "svglite"), paired_umap = c("Seurat",
    "ggplot2", "svglite", "MASS", "sp"), merged_umap = c("ggplot2", "svglite"), overlap = c("ggplot2", "svglite",
    "ragg"), validate = character())

acti_check <- function(cfg = acti_config(), stages = names(acti_stage_packages)) {
  if (any(!stages %in% names(acti_stage_packages)))
    stop("Unknown stage in preflight.")
  packages <- unique(unlist(acti_stage_packages[stages], use.names = FALSE))
  versions <- vapply(packages, acti_package_version, character(1))
  report <- data.frame(package = packages, version = unname(versions), stringsAsFactors = FALSE)
  report$installed <- !is.na(report$version)
  cat(R.version.string, "\nTarget: R 4.1.3 / Seurat 4.4.0\n")
  print(report, row.names = FALSE)
  inputs <- c(case_matrices = cfg$case_dir, control_matrix = cfg$control_dir, HPCA = cfg$hpca_rds, reviewed_annotation_map = cfg$annotation_map,
    optional_STPMGs = cfg$stpmg_csv)
  print(data.frame(input = names(inputs), path = unname(inputs), exists = file.exists(inputs)), row.names = FALSE)
  invisible(report)
}

acti_require_environment <- function(stage) {
  packages <- acti_stage_packages[[stage]]
  missing <- packages[!vapply(packages, function(p) length(find.package(p, quiet = TRUE)) > 0L, logical(1))]
  if (length(missing))
    stop("Missing packages for ", stage, ": ", paste(missing, collapse = ", "))
  if ("Seurat" %in% acti_stage_packages[[stage]]) {
    if (numeric_version(acti_package_version("Seurat")) >= "5.0.0")
      stop("This manuscript workflow requires Seurat 4.x; use the recorded environment.")
    if (getRversion() != "4.1.3")
      warning("The recorded analysis used R 4.1.3; current R is ", getRversion(), ". Results may differ.")
  }
}

acti_manuscript_annotation <- function() {
  data.frame(cluster = as.character(0:25), final_celltype = c("T cell", "Macrophage", "NK cell", "Macrophage",
    "T cell", "Monocyte", "Epithelial", "Endothelial", "B cell", "Epithelial", "Macrophage", "Epithelial",
    "T cell", "Smooth muscle", "Macrophage", "DC", "Macrophage", "Endothelial", "Other", "Mast cell", "Epithelial",
    "T cell", "Macrophage", "Epithelial", "Macrophage", "B cell"), lineage = c("LYMPHOID", "MYELOID", "LYMPHOID",
    "MYELOID", "LYMPHOID", "MYELOID", "OTHER", "OTHER", "LYMPHOID", "OTHER", "MYELOID", "OTHER", "LYMPHOID",
    "OTHER", "MYELOID", "MYELOID", "MYELOID", "OTHER", "OTHER", "MYELOID", "OTHER", "LYMPHOID", "MYELOID",
    "OTHER", "MYELOID", "LYMPHOID"), use_for_DE = !(0:25 %in% c(6, 7, 9, 11, 13, 14, 16, 17, 18, 20, 21, 23)),
    review_status = c("final", "final", "final", "final", "final", "corrected", "final", "final", "final",
      "final", "final", "corrected", "final", "corrected", "provisional_excluded", "final", "provisional_excluded",
      "final", "mixed_excluded", "corrected", "corrected", "provisional_excluded", "final", "final", "final",
      "final"), stringsAsFactors = FALSE)
}

acti_stage_qc <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(SingleCellExperiment)
    library(scDblFinder)
    library(BiocParallel)
    library(harmony)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
  })
  set.seed(20260825)
  options(stringsAsFactors = FALSE)
  project_dir <- cfg$output_dir
  case_dir <- cfg$case_dir
  control_dir <- cfg$control_dir
  outdir <- file.path(project_dir, "01_preannotation")
  rds_dir <- file.path(outdir, "RDS")
  table_dir <- file.path(outdir, "tables")
  plot_dir <- file.path(outdir, "plots")
  dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(case_dir))
    stop("Case directory not found: ", case_dir)
  if (!dir.exists(control_dir))
    stop("Control directory not found: ", control_dir)
  if (getRversion() < "4.1.0" || getRversion() >= "4.2.0") {
    warning("This script targets R 4.1.x; current R is ", getRversion())
  }
  if (packageVersion("Seurat") >= "5.0.0") {
    warning("Seurat 5 detected. This script was written for Seurat 4.x.")
  }
  qc_min_features <- 200L
  qc_max_percent_mt <- 25
  qc_max_counts <- 50000L
  n_hvg <- 2000L
  n_pcs <- 30L
  harmony_dims <- 1:30
  cluster_resolution <- 0.5
  neighbor_k <- 20L
  umap_neighbors <- 30L
  umap_min_dist <- 0.3
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_start.txt"))
  case_info <- tibble(sample_id = c("HIGH1", "LOW1", "HIGH2", "LOW2", "HIGH3", "LOW3", "HIGH4", "LOW4", "HIGH5",
    "LOW5"), patient_id = c("SP019", "SP019", "SP021", "SP021", "SP023", "SP023", "SP024", "SP024", "SP025",
    "SP025"), tissue_type = c("HIGH", "LOW", "HIGH", "LOW", "HIGH", "LOW", "HIGH", "LOW", "HIGH", "LOW"), disease_status = "case",
    dataset = "GSE192483", platform = "10x Genomics", sequencing_batch = "GSE192483")
  control_map <- tribble(~library_id, ~sample_id, ~patient_id, "F00409", "VUHD69", "VUHD69", "F01157", "THD0001",
    "THD0001", "F01174", "THD0002", "THD0002", "F01394", "VUHD71", "VUHD71", "HD68", "VUHD68", "VUHD68", "HD70",
    "VUHD70", "VUHD70")
  stopifnot(!anyDuplicated(case_info$sample_id))
  stopifnot(!anyDuplicated(control_map$sample_id))
  stopifnot(!anyDuplicated(control_map$library_id))
  stopifnot(all(case_info$tissue_type %in% c("HIGH", "LOW")))
  pair_check <- xtabs(~patient_id + tissue_type, data = case_info)
  if (!all(rowSums(pair_check == 1L) == 2L)) {
    stop("Each case patient must have exactly one HIGH and one LOW sample.")
  }
  write.csv(case_info, file.path(table_dir, "case_sample_metadata.csv"), row.names = FALSE)
  write.csv(control_map, file.path(table_dir, "control_sample_metadata.csv"), row.names = FALSE)
  read_10x_gene_expression <- function(data_dir, gene.column = 2L) {
    x <- Read10X(data.dir = data_dir, gene.column = gene.column)
    if (is.list(x)) {
      if ("Gene Expression" %in% names(x)) {
        x <- x[["Gene Expression"]]
      }
      else {
        warning("Multiple 10x assays detected; using the first assay: ", names(x)[1])
        x <- x[[1]]
      }
    }
    if (!inherits(x, "Matrix"))
      x <- as(x, "dgCMatrix")
    x
  }
  assign_control_library <- function(barcodes, library_ids) {
    assigned <- rep(NA_character_, length(barcodes))
    ambiguous <- rep(FALSE, length(barcodes))
    for (id in library_ids) {
      pattern <- paste0("(^|[^[:alnum:]])", id, "([^[:alnum:]]|$)")
      hit <- grepl(pattern, barcodes)
      ambiguous[hit & !is.na(assigned)] <- TRUE
      assigned[hit & is.na(assigned)] <- id
    }
    if (any(ambiguous)) {
      stop("Some control barcodes matched more than one library ID. Examples: ", paste(head(barcodes[ambiguous],
        10L), collapse = ", "))
    }
    assigned
  }
  case_list_raw <- list()
  for (i in seq_len(nrow(case_info))) {
    sid <- case_info$sample_id[i]
    data_path <- file.path(case_dir, sid)
    if (!dir.exists(data_path))
      stop("Case sample directory not found: ", data_path)
    message("Reading case sample: ", sid)
    counts <- read_10x_gene_expression(data_path, gene.column = 2L)
    obj <- CreateSeuratObject(counts = counts, project = sid, min.cells = 3, min.features = 0)
    obj <- RenameCells(obj, add.cell.id = sid)
    obj$sample_id <- sid
    obj$patient_id <- case_info$patient_id[i]
    obj$tissue_type <- case_info$tissue_type[i]
    obj$disease_status <- case_info$disease_status[i]
    obj$dataset <- case_info$dataset[i]
    obj$platform <- case_info$platform[i]
    obj$sequencing_batch <- case_info$sequencing_batch[i]
    case_list_raw[[sid]] <- obj
  }
  control_counts <- read_10x_gene_expression(control_dir, gene.column = 1L)
  control_barcodes <- colnames(control_counts)
  control_library <- assign_control_library(control_barcodes, control_map$library_id)
  barcode_audit <- data.frame(barcode = control_barcodes, assigned_library_id = control_library, stringsAsFactors = FALSE)
  write.csv(head(barcode_audit, 500L), file.path(table_dir, "control_barcode_assignment_preview.csv"), row.names = FALSE)
  assignment_summary <- as.data.frame(table(factor(control_library, levels = control_map$library_id), useNA = "always"),
    stringsAsFactors = FALSE)
  colnames(assignment_summary) <- c("library_id", "raw_barcodes")
  write.csv(assignment_summary, file.path(table_dir, "control_barcode_assignment_summary.csv"), row.names = FALSE)
  unassigned_n <- sum(is.na(control_library))
  if (unassigned_n > 0L) {
    write.csv(data.frame(barcode = head(control_barcodes[is.na(control_library)], 500L)), file.path(table_dir,
      "control_unassigned_barcode_preview.csv"), row.names = FALSE)
    message(unassigned_n, " control-matrix barcodes do not belong to the six selected libraries ", "and will be excluded. A preview was saved for audit.")
  }
  missing_libraries <- setdiff(control_map$library_id, unique(control_library))
  if (length(missing_libraries) > 0L) {
    stop("No barcodes were found for control libraries: ", paste(missing_libraries, collapse = ", "))
  }
  control_objects_raw <- list()
  for (i in seq_len(nrow(control_map))) {
    library_id <- control_map$library_id[i]
    sample_id <- control_map$sample_id[i]
    keep <- which(control_library == library_id)
    counts <- control_counts[, keep, drop = FALSE]
    obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 0)
    obj <- RenameCells(obj, add.cell.id = sample_id)
    obj$sample_id <- sample_id
    obj$patient_id <- control_map$patient_id[i]
    obj$tissue_type <- "CONTROL"
    obj$disease_status <- "control"
    obj$dataset <- "GSE135893"
    obj$platform <- "10x Genomics"
    obj$sequencing_batch <- "GSE135893"
    obj$source_library_id <- library_id
    control_objects_raw[[sample_id]] <- obj
  }
  rm(control_counts)
  gc()
  all_raw <- c(case_list_raw, control_objects_raw)
  stopifnot(length(all_raw) == 16L)
  stopifnot(all(vapply(all_raw, ncol, integer(1)) > 0L))
  process_one_sample <- function(obj, qc_min_features = 200L, qc_mt = 25, qc_max_counts = 50000L, seed = 20260825) {
    sid <- unique(obj$sample_id)
    if (length(sid) != 1L) {
      stop("Each object must contain exactly one sample_id.")
    }
    mt_features <- grep("^MT-", rownames(obj), value = TRUE)
    if (length(mt_features) == 0L) {
      stop("No ^MT- genes were detected in ", sid, ". Confirm that gene symbols, rather than Ensembl IDs, are row names.")
    }
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    md_before <- obj[[]]
    before_n <- nrow(md_before)
    keep_primary <- md_before$nFeature_RNA >= qc_min_features & md_before$percent.mt < qc_mt & md_before$nCount_RNA <
      qc_max_counts
    keep_reference_500_7500 <- md_before$nFeature_RNA >= 500 & md_before$nFeature_RNA <= 7500 & md_before$percent.mt <
      25 & md_before$nCount_RNA < 50000
    obj <- subset(obj, cells = rownames(md_before)[keep_primary])
    after_qc_n <- ncol(obj)
    if (after_qc_n < 50L) {
      warning("Sample ", sid, " has fewer than 50 cells after QC.")
    }
    if (after_qc_n < 10L) {
      stop("Too few cells for reliable doublet detection in sample ", sid)
    }
    set.seed(seed)
    sce <- as.SingleCellExperiment(obj, assay = "RNA")
    sce$sample_id <- sid
    sce <- scDblFinder(sce, samples = "sample_id", BPPARAM = SerialParam())
    obj$doublet_score <- colData(sce)$scDblFinder.score
    obj$doublet_class <- colData(sce)$scDblFinder.class
    doublet_n <- sum(obj$doublet_class == "doublet", na.rm = TRUE)
    singlet_n <- sum(obj$doublet_class == "singlet", na.rm = TRUE)
    doublet_rate <- doublet_n/after_qc_n
    summary_row <- data.frame(sample_id = sid, patient_id = unique(obj$patient_id), tissue_type = unique(obj$tissue_type),
      disease_status = unique(obj$disease_status), dataset = unique(obj$dataset), cells_before_qc = before_n,
      cells_after_primary_qc = after_qc_n, cells_passing_reference_500_7500 = sum(keep_reference_500_7500),
      predicted_doublets = doublet_n, predicted_doublet_percent = 100 * doublet_rate, retained_singlets = singlet_n,
      median_nFeature_before = median(md_before$nFeature_RNA), median_nCount_before = median(md_before$nCount_RNA),
      median_percent_mt_before = median(md_before$percent.mt), median_nFeature_after_qc = median(obj$nFeature_RNA),
      median_nCount_after_qc = median(obj$nCount_RNA), median_percent_mt_after_qc = median(obj$percent.mt),
      stringsAsFactors = FALSE)
    obj@misc$processing_info <- as.list(summary_row[1, ])
    obj@misc$processing_info$qc_min_features <- qc_min_features
    obj@misc$processing_info$qc_mt <- qc_mt
    obj@misc$processing_info$qc_max_counts <- qc_max_counts
    obj@misc$processing_info$scDblFinder_version <- as.character(packageVersion("scDblFinder"))
    obj <- subset(obj, subset = doublet_class == "singlet")
    list(object = obj, summary = summary_row)
  }
  processed_results <- lapply(all_raw, process_one_sample, qc_min_features = qc_min_features, qc_mt = qc_max_percent_mt,
    qc_max_counts = qc_max_counts)
  qc_summary <- bind_rows(lapply(processed_results, `[[`, "summary"))
  processed_list <- lapply(processed_results, `[[`, "object")
  rm(processed_results)
  gc()
  write.csv(qc_summary, file.path(table_dir, "QC_scDblFinder_summary_by_sample.csv"), row.names = FALSE)
  qc_total <- qc_summary %>% summarise(n_samples = n(), cells_before_qc = sum(cells_before_qc), cells_after_primary_qc = sum(cells_after_primary_qc),
    predicted_doublets = sum(predicted_doublets), predicted_doublet_percent = 100 * predicted_doublets/cells_after_primary_qc,
    retained_singlets = sum(retained_singlets))
  write.csv(qc_total, file.path(table_dir, "QC_scDblFinder_summary_total.csv"), row.names = FALSE)
  for (sid in names(processed_list)) {
    saveRDS(processed_list[[sid]], file.path(rds_dir, paste0(sid, "_QC_singlet.rds")))
  }
  qc_long <- qc_summary %>% select(sample_id, tissue_type, cells_before_qc, cells_after_primary_qc, retained_singlets) %>%
    pivot_longer(cols = c(cells_before_qc, cells_after_primary_qc, retained_singlets), names_to = "processing_stage",
      values_to = "n_cells")
  p_qc_cells <- ggplot(qc_long, aes(x = sample_id, y = n_cells, fill = processing_stage)) + geom_col(position = "dodge") +
    facet_grid(. ~ tissue_type, scales = "free_x", space = "free_x") + theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45,
    hjust = 1)) + labs(x = NULL, y = "Number of cells", fill = "Stage")
  ggsave(file.path(plot_dir, "QC_cell_retention_by_sample.pdf"), p_qc_cells, width = 14, height = 6)
  p_doublet <- ggplot(qc_summary, aes(x = sample_id, y = predicted_doublet_percent, fill = tissue_type)) + geom_col() +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL, y = "Predicted doublets (%)",
    fill = "Group")
  ggsave(file.path(plot_dir, "scDblFinder_rate_by_sample.pdf"), p_doublet, width = 12, height = 6)
  scRNA <- merge(x = processed_list[[1]], y = processed_list[-1], add.cell.ids = NULL, merge.data = FALSE, project = "GSE192483_GSE135893")
  DefaultAssay(scRNA) <- "RNA"
  scRNA$tissue_type <- factor(as.character(scRNA$tissue_type), levels = c("CONTROL", "LOW", "HIGH"))
  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
  scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA), verbose = FALSE)
  scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA), npcs = n_pcs, seed.use = 20260825, verbose = FALSE)
  saveRDS(scRNA, file.path(rds_dir, "merged_singlets_before_Harmony.rds"))
  p_elbow <- ElbowPlot(scRNA, ndims = n_pcs)
  ggsave(file.path(plot_dir, "PCA_ElbowPlot.pdf"), p_elbow, width = 7, height = 5)
  scRNA <- RunHarmony(object = scRNA, group.by.vars = "sample_id", reduction = "pca", dims.use = 1:n_pcs, assay.use = "RNA",
    plot_convergence = TRUE, verbose = TRUE)
  scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = harmony_dims, k.param = neighbor_k, verbose = FALSE)
  scRNA <- FindClusters(scRNA, resolution = cluster_resolution, algorithm = 1, random.seed = 20260825, verbose = FALSE)
  scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = harmony_dims, n.neighbors = umap_neighbors, min.dist = umap_min_dist,
    seed.use = 20260825, verbose = FALSE)
  cluster_sample_table <- as.data.frame(table(scRNA$seurat_clusters, scRNA$sample_id), stringsAsFactors = FALSE)
  colnames(cluster_sample_table) <- c("cluster", "sample_id", "cell_count")
  write.csv(cluster_sample_table, file.path(table_dir, "cluster_cell_counts_by_sample.csv"), row.names = FALSE)
  plot_umap <- function(group_by, title, filename, label = FALSE, width = 10) {
    p <- DimPlot(scRNA, reduction = "umap", group.by = group_by, label = label, repel = label, raster = FALSE,
      pt.size = 0.15) + ggtitle(title)
    ggsave(file.path(plot_dir, filename), p, width = width, height = 7)
  }
  plot_umap("seurat_clusters", "Harmony-integrated representation by cluster", "UMAP_by_cluster.pdf", label = TRUE,
    width = 10)
  plot_umap("sample_id", "Harmony-integrated representation by sample", "UMAP_by_sample.pdf", width = 13)
  plot_umap("tissue_type", "Harmony-integrated representation by tissue group", "UMAP_by_tissue_group.pdf", width = 9)
  plot_umap("dataset", "Harmony-integrated representation by dataset", "UMAP_by_dataset.pdf", width = 9)
  Idents(scRNA) <- "seurat_clusters"
  scRNA_marker <- subset(scRNA, downsample = 1000)
  cluster_markers <- FindAllMarkers(scRNA_marker, assay = "RNA", only.pos = TRUE, test.use = "wilcox", min.pct = 0.25,
    logfc.threshold = 0.25, return.thresh = 0.05, verbose = FALSE)
  cluster_markers <- cluster_markers %>% mutate(cluster = as.character(cluster)) %>% arrange(cluster, p_val_adj,
    desc(avg_log2FC))
  cluster_markers_ranked <- cluster_markers %>% group_by(cluster) %>% mutate(marker_rank = row_number()) %>%
    ungroup()
  write.csv(cluster_markers_ranked, file.path(table_dir, "cluster_markers_for_annotation.csv"), row.names = FALSE)
  write.csv(filter(cluster_markers_ranked, marker_rank <= 20L), file.path(table_dir, "cluster_top20_markers_for_annotation.csv"),
    row.names = FALSE)
  canonical_markers <- c("PTPRC", "CD3D", "CD3E", "TRAC", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD79A", "CD37",
    "MZB1", "JCHAIN", "CD14", "FCGR3A", "LYZ", "LST1", "FCER1A", "CST3", "MARCO", "FABP4", "PPARG", "EPCAM",
    "KRT8", "KRT18", "SCGB1A1", "SFTPC", "AGER", "PECAM1", "VWF", "COL1A1", "COL1A2", "DCN", "ACTA2", "RGS5",
    "TPSAB1")
  canonical_markers_present <- intersect(canonical_markers, rownames(scRNA))
  canonical_markers_missing <- setdiff(canonical_markers, rownames(scRNA))
  write.csv(data.frame(gene = canonical_markers_missing), file.path(table_dir, "canonical_markers_missing_from_matrix.csv"),
    row.names = FALSE)
  if (length(canonical_markers_present) > 0L) {
    p_canonical <- DotPlot(scRNA, features = canonical_markers_present, group.by = "seurat_clusters") + RotatedAxis()
    write.csv(p_canonical$data, file.path(table_dir, "canonical_marker_DotPlot_values.csv"), row.names = FALSE)
    ggsave(file.path(plot_dir, "canonical_marker_DotPlot.pdf"), p_canonical, width = 16, height = 8)
  }
  saveRDS(scRNA, file.path(rds_dir, "scRNA_Harmony_clustered_PRE_ANNOTATION.rds"))
  write.csv(scRNA[[]], file.path(table_dir, "cell_metadata_PRE_ANNOTATION.csv"), row.names = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_final.txt"))
  message("Pre-annotation workflow completed.")
  message("Output directory: ", outdir)
  message("Next step: review cluster markers and canonical-marker DotPlot, then assign final annotations.")
  invisible(NULL)
}

acti_stage_annotation_evidence <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleR)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(MAST)
    library(Matrix)
    library(edgeR)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
  })
  set.seed(20260825)
  options(stringsAsFactors = FALSE)
  project_dir <- cfg$output_dir
  preannotation_outdir <- file.path(project_dir, "01_preannotation")
  input_rds <- file.path(preannotation_outdir, "RDS", "scRNA_Harmony_clustered_PRE_ANNOTATION.rds")
  outdir <- file.path(project_dir, "02_annotation")
  rds_dir <- file.path(outdir, "RDS")
  table_dir <- file.path(outdir, "tables")
  plot_dir <- file.path(outdir, "plots")
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(input_rds)) {
    stop("Pre-annotation object not found: ", input_rds)
  }
  ref_path <- cfg$hpca_rds
  if (!file.exists(ref_path))
    stop("HPCA reference not found: ", ref_path)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_start.txt"))
  scRNA <- readRDS(input_rds)
  DefaultAssay(scRNA) <- "RNA"
  required_metadata <- c("sample_id", "patient_id", "tissue_type", "disease_status", "dataset", "seurat_clusters")
  missing_metadata <- setdiff(required_metadata, colnames(scRNA[[]]))
  if (length(missing_metadata) > 0L) {
    stop("Missing required metadata: ", paste(missing_metadata, collapse = ", "))
  }
  if (!"harmony" %in% Reductions(scRNA))
    stop("Harmony reduction is missing")
  if (!"umap" %in% Reductions(scRNA))
    stop("UMAP reduction is missing")
  scRNA$tissue_type <- factor(as.character(scRNA$tissue_type), levels = c("CONTROL", "LOW", "HIGH"))
  sample_design <- unique(scRNA[[]][, c("sample_id", "patient_id", "tissue_type", "dataset")])
  sample_design <- sample_design[order(sample_design$tissue_type, sample_design$patient_id), ]
  write.csv(sample_design, file.path(table_dir, "sample_design_used_for_DE.csv"), row.names = FALSE)
  case_design <- sample_design[sample_design$tissue_type %in% c("LOW", "HIGH"), ]
  pair_check <- xtabs(~patient_id + tissue_type, data = case_design)
  if (!all(rowSums(pair_check[, c("LOW", "HIGH"), drop = FALSE] == 1L) == 2L)) {
    stop("HIGH/LOW pairing is incomplete or duplicated in sample metadata")
  }
  refdata <- readRDS(ref_path)
  if (!"label.main" %in% colnames(colData(refdata))) {
    stop("HPCA reference does not contain colData(ref)$label.main")
  }
  if (!"logcounts" %in% assayNames(refdata)) {
    stop("HPCA reference does not contain a logcounts assay")
  }
  rna_data <- GetAssayData(scRNA, assay = "RNA", slot = "data")
  if (nrow(rna_data) == 0L || ncol(rna_data) == 0L) {
    scRNA <- NormalizeData(scRNA, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000,
      verbose = FALSE)
    rna_data <- GetAssayData(scRNA, assay = "RNA", slot = "data")
  }
  cellpred <- SingleR(test = rna_data, ref = refdata, labels = refdata$label.main, clusters = as.character(scRNA$seurat_clusters),
    assay.type.ref = "logcounts", de.method = "wilcox", prune = TRUE)
  cellpred_df <- as.data.frame(cellpred)
  singleR_annotation <- data.frame(cluster = rownames(cellpred_df), SingleR_label = as.character(cellpred_df$labels),
    SingleR_pruned_label = as.character(cellpred_df$pruned.labels), stringsAsFactors = FALSE)
  singleR_annotation$SingleR_pruned_label[is.na(singleR_annotation$SingleR_pruned_label)] <- singleR_annotation$SingleR_label[is.na(singleR_annotation$SingleR_pruned_label)]
  if ("delta.next" %in% colnames(cellpred_df)) {
    singleR_annotation$delta_next <- cellpred_df$delta.next
  }
  else {
    singleR_annotation$delta_next <- NA_real_
  }
  write.csv(singleR_annotation, file.path(table_dir, "SingleR_cluster_annotation.csv"), row.names = FALSE)
  auto_label_map <- setNames(singleR_annotation$SingleR_pruned_label, singleR_annotation$cluster)
  scRNA$celltype_SingleR <- unname(auto_label_map[as.character(scRNA$seurat_clusters)])
  p_singleR <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_SingleR", label = TRUE, repel = TRUE,
    pt.size = 0.15, raster = FALSE) + ggtitle("Provisional cluster-level SingleR annotation")
  ggsave(file.path(plot_dir, "UMAP_provisional_SingleR_annotation.pdf"), p_singleR, width = 13, height = 8)
  Idents(scRNA) <- "seurat_clusters"
  scRNA_marker <- subset(scRNA, downsample = 1000)
  markers_MAST <- FindAllMarkers(scRNA_marker, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", return.thresh = 0.05, verbose = FALSE)
  markers_MAST <- markers_MAST %>% mutate(cluster = as.character(cluster)) %>% group_by(cluster) %>% arrange(desc(avg_log2FC),
    desc(pct.1), p_val_adj, .by_group = TRUE) %>% mutate(marker_rank = row_number()) %>% ungroup()
  write.csv(markers_MAST, file.path(table_dir, "All_cluster_markers_MAST.csv"), row.names = FALSE)
  write.csv(filter(markers_MAST, marker_rank <= 20L), file.path(table_dir, "All_cluster_top20_markers_MAST.csv"),
    row.names = FALSE)
  write.csv(filter(markers_MAST, marker_rank <= 50L), file.path(table_dir, "All_cluster_top50_markers_MAST.csv"),
    row.names = FALSE)
  canonical_markers <- c("PTPRC", "CD3D", "CD3E", "TRAC", "IL7R", "LTB", "CD8A", "CD8B", "NKG7", "GNLY", "KLRD1",
    "MS4A1", "CD79A", "CD37", "CD74", "MZB1", "JCHAIN", "SDC1", "LYZ", "LST1", "TYROBP", "FCER1G", "CTSS",
    "CD14", "FCGR3A", "FCER1A", "CLEC10A", "CST3", "MARCO", "FABP4", "PPARG", "APOC1", "SPP1", "FCN1", "S100A8",
    "S100A9", "EPCAM", "KRT8", "KRT18", "SCGB1A1", "SFTPC", "AGER", "PECAM1", "VWF", "COL1A1", "COL1A2", "DCN",
    "ACTA2", "RGS5")
  marker_panel <- intersect(canonical_markers, rownames(scRNA))
  if (length(marker_panel) > 0L) {
    p_marker <- DotPlot(scRNA, features = marker_panel, group.by = "seurat_clusters") + RotatedAxis()
    write.csv(p_marker$data, file.path(table_dir, "canonical_marker_DotPlot_values.csv"), row.names = FALSE)
    ggsave(file.path(plot_dir, "canonical_marker_DotPlot.pdf"), p_marker, width = 18, height = 9)
  }
  annotation_template <- singleR_annotation
  annotation_template$final_celltype <- ""
  annotation_template$lineage <- ""
  annotation_template$use_for_DE <- FALSE
  annotation_template$review_status <- "REVIEW_REQUIRED"
  write.csv(annotation_template, file.path(outdir, "annotation_map_TO_REVIEW.csv"), row.names = FALSE)
  saveRDS(scRNA, file.path(rds_dir, "scRNA_provisional_SingleR_annotation.rds"))
  message("Annotation evidence saved. Supply a reviewed map to the annotation stage.")
  invisible(NULL)
}

acti_stage_annotation_audit <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })
  project_dir <- cfg$output_dir
  input_candidates <- c(file.path(project_dir, "02_annotation", "RDS", "scRNA_provisional_SingleR_annotation.rds"),
    file.path(project_dir, "01_preannotation", "RDS", "scRNA_Harmony_clustered_PRE_ANNOTATION.rds"))
  input_rds <- input_candidates[file.exists(input_candidates)][1]
  if (is.na(input_rds))
    stop("No provisional/pre-annotation RDS was found")
  outdir <- file.path(project_dir, "02_annotation")
  table_dir <- file.path(outdir, "tables")
  plot_dir <- file.path(outdir, "plots")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  scRNA <- readRDS(input_rds)
  DefaultAssay(scRNA) <- "RNA"
  Idents(scRNA) <- "seurat_clusters"
  target_clusters <- c("8", "10", "14", "16", "18", "21")
  target_clusters <- intersect(target_clusters, unique(as.character(scRNA$seurat_clusters)))
  marker_panel <- c("PTPRC", "LYZ", "LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "C1QA", "C1QB", "C1QC", "MARCO",
    "FABP4", "CHI3L1", "CHIT1", "CD3D", "CD3E", "TRAC", "CD2", "IL7R", "CD8A", "NKG7", "GNLY", "MS4A1", "CD79A",
    "CD79B", "CD19", "JCHAIN", "MZB1", "SDC1", "EPCAM", "KRT8", "PECAM1", "VWF", "COL1A1", "DCN", "MKI67",
    "TOP2A", "STMN1")
  marker_panel <- intersect(marker_panel, rownames(scRNA))
  target_obj <- subset(scRNA, idents = target_clusters)
  p <- DotPlot(target_obj, features = marker_panel, group.by = "seurat_clusters") + RotatedAxis()
  write.csv(p$data, file.path(table_dir, "uncertain_clusters_targeted_DotPlot_values.csv"), row.names = FALSE)
  ggsave(file.path(plot_dir, "uncertain_clusters_targeted_DotPlot.pdf"), p, width = 16, height = 7)
  required_genes <- c("LST1", "TYROBP", "CD3D", "TRAC", "MS4A1", "CD79A", "EPCAM", "KRT8", "PECAM1", "VWF", "MKI67",
    "TOP2A")
  present <- intersect(required_genes, rownames(scRNA))
  expr <- t(as.matrix(GetAssayData(target_obj, assay = "RNA", slot = "data")[present, , drop = FALSE]))
  expr <- as.data.frame(expr)
  expr$cluster <- as.character(target_obj$seurat_clusters)
  is_on <- function(df, genes) {
    genes <- intersect(genes, colnames(df))
    if (length(genes) == 0L)
      return(rep(FALSE, nrow(df)))
    rowSums(df[, genes, drop = FALSE] > 0) == length(genes)
  }
  expr$myeloid_core <- is_on(expr, c("LST1", "TYROBP"))
  expr$t_core <- is_on(expr, c("CD3D", "TRAC"))
  expr$b_core <- is_on(expr, c("MS4A1", "CD79A"))
  expr$epithelial_core <- is_on(expr, c("EPCAM", "KRT8"))
  expr$endothelial_core <- is_on(expr, c("PECAM1", "VWF"))
  expr$cycling_core <- is_on(expr, c("MKI67", "TOP2A"))
  coexpression_summary <- expr %>% group_by(cluster) %>% summarise(n_cells = n(), pct_myeloid_core = 100 * mean(myeloid_core),
    pct_T_core = 100 * mean(t_core), pct_B_core = 100 * mean(b_core), pct_epithelial_core = 100 * mean(epithelial_core),
    pct_endothelial_core = 100 * mean(endothelial_core), pct_cycling_core = 100 * mean(cycling_core), pct_myeloid_and_T = 100 *
      mean(myeloid_core & t_core), pct_myeloid_and_B = 100 * mean(myeloid_core & b_core), .groups = "drop")
  write.csv(coexpression_summary, file.path(table_dir, "uncertain_clusters_cell_level_coexpression.csv"), row.names = FALSE)
  sample_cluster <- scRNA[[]] %>% filter(as.character(seurat_clusters) %in% target_clusters) %>% count(seurat_clusters,
    sample_id, tissue_type, name = "cell_count") %>% group_by(seurat_clusters) %>% mutate(cluster_fraction = cell_count/sum(cell_count)) %>%
    ungroup()
  write.csv(sample_cluster, file.path(table_dir, "uncertain_clusters_sample_composition.csv"), row.names = FALSE)
  qc_by_cluster <- scRNA[[]] %>% filter(as.character(seurat_clusters) %in% target_clusters) %>% group_by(seurat_clusters) %>%
    summarise(n_cells = n(), median_nFeature_RNA = median(nFeature_RNA), median_nCount_RNA = median(nCount_RNA),
      median_percent_mt = median(percent.mt), q90_percent_mt = quantile(percent.mt, 0.9),
      .groups = "drop")
  write.csv(qc_by_cluster, file.path(table_dir, "uncertain_clusters_QC_summary.csv"), row.names = FALSE)
  message("Targeted annotation checks completed: ", outdir)
  invisible(NULL)
}

acti_stage_annotation <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
  })
  project_dir <- cfg$output_dir
  annotation_outdir <- file.path(project_dir, "02_annotation")
  rds_dir <- file.path(annotation_outdir, "RDS")
  table_dir <- file.path(annotation_outdir, "tables")
  plot_dir <- file.path(annotation_outdir, "plots")
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  input_candidates <- c(file.path(rds_dir, "scRNA_provisional_SingleR_annotation.rds"), file.path(project_dir,
    "01_preannotation", "RDS", "scRNA_Harmony_clustered_PRE_ANNOTATION.rds"))
  input_rds <- input_candidates[file.exists(input_candidates)][1]
  if (is.na(input_rds))
    stop("No provisional/pre-annotation RDS was found")
  scRNA <- readRDS(input_rds)
  if (file.exists(cfg$annotation_map)) {
    annotation_map <- read.csv(cfg$annotation_map, check.names = FALSE)
  }
  else if (isTRUE(cfg$use_manuscript_annotation)) {
    annotation_map <- acti_manuscript_annotation()
    warning("Using the historical manuscript cluster map. Confirm the cluster identities against annotation evidence.")
  }
  else {
    write.csv(acti_manuscript_annotation(), file.path(annotation_outdir, "manuscript_annotation_map_REFERENCE_ONLY.csv"),
      row.names = FALSE)
    stop("Manual review checkpoint: supply --annotation-map=... or explicitly use --use-manuscript-annotation for the original manuscript clustering.")
  }
  required_map <- c("cluster", "final_celltype", "lineage", "use_for_DE", "review_status")
  if (!all(required_map %in% names(annotation_map)))
    stop("Annotation map lacks: ", paste(setdiff(required_map, names(annotation_map)), collapse = ", "))
  if (anyDuplicated(annotation_map$cluster) || anyNA(annotation_map[, required_map]))
    stop("Annotation map has duplicate clusters or missing values.")
  if (any(!nzchar(trimws(annotation_map$final_celltype))) || any(annotation_map$review_status == "REVIEW_REQUIRED"))
    stop("Complete manual annotation review first.")
  de_use <- as.logical(annotation_map$use_for_DE)
  if (anyNA(de_use))
    stop("use_for_DE must contain TRUE or FALSE.")
  if (any(de_use & !annotation_map$review_status %in% c("final", "corrected")))
    stop("All cells entering DE require final/corrected review status.")
  annotation_map$cluster <- as.character(annotation_map$cluster)
  annotation_map$use_for_DE <- as.logical(annotation_map$use_for_DE)
  object_clusters <- unique(as.character(scRNA$seurat_clusters))
  if (!setequal(object_clusters, annotation_map$cluster)) {
    stop("Clusters in the annotation map do not match the Seurat object")
  }
  idx <- match(as.character(scRNA$seurat_clusters), annotation_map$cluster)
  scRNA$celltype_major <- annotation_map$final_celltype[idx]
  scRNA$lineage_major <- annotation_map$lineage[idx]
  scRNA$use_for_DE <- annotation_map$use_for_DE[idx]
  scRNA$annotation_review_status <- annotation_map$review_status[idx]
  celltype_levels <- c("T cell", "NK cell", "Macrophage", "Monocyte", "B cell", "DC", "Mast cell", "Epithelial",
    "Endothelial", "Smooth muscle", "Other")
  if (any(!as.character(scRNA$celltype_major) %in% celltype_levels))
    stop("Unsupported major cell type in map.")
  scRNA$celltype_major <- factor(scRNA$celltype_major, levels = celltype_levels)
  scRNA$lineage_major <- factor(scRNA$lineage_major, levels = c("LYMPHOID", "MYELOID", "OTHER"))
  celltype_colours <- c(`T cell` = "#A6CEE3", `NK cell` = "#1F78B4", Macrophage = "#E88C8D", Monocyte = "#33A02C",
    `B cell` = "#A6D17A", DC = "#E31A1C", `Mast cell` = "#6A3D9A", Epithelial = "#FDBF6F", Endothelial = "#FF7F00",
    `Smooth muscle` = "#CAB2D6", Other = "#BDBDBD")
  lineage_colours <- c(LYMPHOID = "#438CC7", MYELOID = "#E83E3E", OTHER = "#C9C9C9")
  p_celltype <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_major", split.by = "tissue_type", cols = celltype_colours,
    raster = FALSE, pt.size = 0.15, ncol = 3) + ggtitle("Broad cell classes")
  ggsave(file.path(plot_dir, "UMAP_major_celltypes_by_group.pdf"), p_celltype, width = 17, height = 7)
  p_lineage <- DimPlot(scRNA, reduction = "umap", group.by = "lineage_major", split.by = "tissue_type", cols = lineage_colours,
    raster = FALSE, pt.size = 0.15, ncol = 3) + ggtitle("Broad immune lineages")
  ggsave(file.path(plot_dir, "UMAP_lymphoid_myeloid_other_by_group.pdf"), p_lineage, width = 17, height = 7)
  annotation_counts <- scRNA[[]] %>% count(sample_id, patient_id, tissue_type, celltype_major, lineage_major,
    use_for_DE, name = "cell_count")
  write.csv(annotation_counts, file.path(table_dir, "major_annotation_cell_counts_by_sample.csv"), row.names = FALSE)
  write.csv(annotation_map, file.path(annotation_outdir, "annotation_map_FINAL.csv"), row.names = FALSE)
  saveRDS(scRNA, file.path(rds_dir, "scRNA_major_classes_annotated.rds"))
  message("Broad annotation completed: ", annotation_outdir)
  invisible(NULL)
}

acti_stage_lineage <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(edgeR)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
  })
  set.seed(20260826)
  options(stringsAsFactors = FALSE)
  project_dir <- cfg$output_dir
  annotation_outdir <- file.path(project_dir, "02_annotation")
  input_rds <- acti_input_object(cfg)
  needs_annotation_fallback <- FALSE
  outdir <- file.path(project_dir, "03_lineage_pseudobulk")
  rds_dir <- file.path(outdir, "RDS")
  table_dir <- file.path(outdir, "tables")
  plot_dir <- file.path(outdir, "plots")
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  min_cells_per_pseudobulk <- 20L
  min_complete_pairs <- 3L
  min_samples_context_group <- 3L
  candidate_fdr <- 0.1
  candidate_abs_log2fc <- 0.25
  conventional_fdr <- 0.05
  conventional_abs_log2fc <- 0.5
  run_contextual_control_comparisons <- TRUE
  run_leave_one_control_out <- TRUE
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_start.txt"))
  scRNA <- readRDS(input_rds)
  DefaultAssay(scRNA) <- "RNA"
  md <- scRNA[[]]
  if (!all(c("celltype_major", "lineage_major", "use_for_DE") %in% colnames(md))) {
    stop("Input object lacks reviewed major annotations. Run the annotation stage first.")
  }
  md <- scRNA[[]]
  if (!"celltype_major" %in% colnames(md) && "final_celltype" %in% colnames(md)) {
    scRNA$celltype_major <- md$final_celltype
  }
  if (!"lineage_major" %in% colnames(md) && "lineage" %in% colnames(md)) {
    scRNA$lineage_major <- md$lineage
  }
  if (!"annotation_review_status" %in% colnames(md)) {
    if ("review_status" %in% colnames(md)) {
      scRNA$annotation_review_status <- md$review_status
    }
    else {
      scRNA$annotation_review_status <- NA_character_
    }
  }
  md <- scRNA[[]]
  required_metadata <- c("sample_id", "patient_id", "tissue_type", "dataset", "seurat_clusters", "celltype_major",
    "lineage_major", "use_for_DE")
  missing_metadata <- setdiff(required_metadata, colnames(md))
  if (length(missing_metadata) > 0L) {
    stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
  }
  md$sample_id <- as.character(md$sample_id)
  md$patient_id <- as.character(md$patient_id)
  md$tissue_type <- toupper(as.character(md$tissue_type))
  md$dataset <- as.character(md$dataset)
  md$seurat_clusters <- as.character(md$seurat_clusters)
  md$celltype_major <- as.character(md$celltype_major)
  md$lineage_major <- toupper(as.character(md$lineage_major))
  md$use_for_DE <- as.logical(md$use_for_DE)
  md$annotation_review_status <- as.character(md$annotation_review_status)
  if (!all(md$tissue_type %in% c("CONTROL", "LOW", "HIGH"))) {
    stop("tissue_type must contain only CONTROL, LOW, and HIGH")
  }
  if (anyNA(md$use_for_DE))
    stop("use_for_DE contains NA")
  if (anyNA(md$celltype_major) || any(trimws(md$celltype_major) == "")) {
    stop("Some cells lack a final major-cell-class annotation")
  }
  if (anyNA(md$lineage_major) || any(trimws(md$lineage_major) == "")) {
    stop("Some cells lack a final lineage annotation")
  }
  scRNA$sample_id <- md$sample_id
  scRNA$patient_id <- md$patient_id
  scRNA$tissue_type <- factor(md$tissue_type, levels = c("CONTROL", "LOW", "HIGH"))
  scRNA$dataset <- md$dataset
  scRNA$celltype_major <- md$celltype_major
  scRNA$lineage_major <- md$lineage_major
  scRNA$use_for_DE <- md$use_for_DE
  scRNA$annotation_review_status <- md$annotation_review_status
  md <- scRNA[[]]
  cluster_annotation_audit <- md %>% transmute(cluster = as.character(seurat_clusters), celltype_major = as.character(celltype_major),
    lineage_major = as.character(lineage_major), use_for_DE = as.logical(use_for_DE), review_status = as.character(annotation_review_status)) %>%
    distinct() %>% arrange(as.numeric(cluster), celltype_major)
  cluster_uniqueness <- cluster_annotation_audit %>% count(cluster, name = "n_distinct_annotation_rows")
  if (any(cluster_uniqueness$n_distinct_annotation_rows != 1L)) {
    bad <- cluster_uniqueness$cluster[cluster_uniqueness$n_distinct_annotation_rows != 1L]
    stop("A cluster maps to more than one annotation: ", paste(bad, collapse = ", "))
  }
  valid_de_lineages <- c("MYELOID", "LYMPHOID")
  bad_de_lineage <- cluster_annotation_audit$use_for_DE & !cluster_annotation_audit$lineage_major %in% valid_de_lineages
  if (any(bad_de_lineage)) {
    stop("Some use_for_DE=TRUE clusters are not MYELOID/LYMPHOID: ", paste(cluster_annotation_audit$cluster[bad_de_lineage],
      collapse = ", "))
  }
  accepted_review_status <- c("final", "corrected")
  bad_de_review <- cluster_annotation_audit$use_for_DE & (is.na(cluster_annotation_audit$review_status) | !cluster_annotation_audit$review_status %in%
    accepted_review_status)
  if (any(bad_de_review)) {
    stop("Some use_for_DE=TRUE clusters do not have final/corrected review status: ", paste(cluster_annotation_audit$cluster[bad_de_review],
      collapse = ", "))
  }
  write.csv(cluster_annotation_audit, file.path(table_dir, "annotation_map_audit.csv"), row.names = FALSE)
  cluster_cell_count <- md %>% count(seurat_clusters, celltype_major, lineage_major, use_for_DE, annotation_review_status,
    name = "n_cells") %>% arrange(as.numeric(as.character(seurat_clusters)))
  write.csv(cluster_cell_count, file.path(table_dir, "annotation_cluster_cell_counts.csv"), row.names = FALSE)
  excluded_annotation_count <- cluster_cell_count %>% filter(!use_for_DE) %>% arrange(desc(n_cells))
  write.csv(excluded_annotation_count, file.path(table_dir, "annotation_clusters_excluded_from_DE.csv"), row.names = FALSE)
  sample_design <- md %>% select(sample_id, patient_id, tissue_type, dataset) %>% distinct()
  if (anyDuplicated(sample_design$sample_id)) {
    stop("A sample_id maps to more than one patient/group/dataset")
  }
  all_lineage_counts <- md %>% filter(lineage_major %in% valid_de_lineages) %>% count(sample_id, lineage_major,
    name = "all_annotated_cells")
  accepted_lineage_counts <- md %>% filter(use_for_DE, lineage_major %in% valid_de_lineages) %>% count(sample_id,
    lineage_major, name = "cells_used_for_DE")
  sample_lineage_audit <- tidyr::crossing(sample_id = sample_design$sample_id, lineage_major = valid_de_lineages) %>%
    left_join(sample_design, by = "sample_id") %>% left_join(all_lineage_counts, by = c("sample_id", "lineage_major")) %>%
    left_join(accepted_lineage_counts, by = c("sample_id", "lineage_major")) %>% mutate(all_annotated_cells = replace_na(all_annotated_cells,
    0L), cells_used_for_DE = replace_na(cells_used_for_DE, 0L), excluded_cells = all_annotated_cells - cells_used_for_DE,
    eligible_pseudobulk = cells_used_for_DE >= min_cells_per_pseudobulk, exclusion_reason = ifelse(eligible_pseudobulk,
      "included", paste0("fewer_than_", min_cells_per_pseudobulk, "_accepted_cells"))) %>% arrange(lineage_major,
    tissue_type, patient_id)
  write.csv(sample_lineage_audit, file.path(table_dir, "sample_lineage_cell_count_audit.csv"), row.names = FALSE)
  pair_audit <- sample_lineage_audit %>% filter(tissue_type %in% c("LOW", "HIGH")) %>% select(lineage_major,
    patient_id, tissue_type, sample_id, cells_used_for_DE, eligible_pseudobulk) %>% pivot_wider(names_from = tissue_type,
    values_from = c(sample_id, cells_used_for_DE, eligible_pseudobulk), names_sep = "_") %>% mutate(complete_pair = !is.na(sample_id_LOW) &
    !is.na(sample_id_HIGH), eligible_pair = complete_pair & eligible_pseudobulk_LOW & eligible_pseudobulk_HIGH,
    pair_status = case_when(!complete_pair ~ "missing_LOW_or_HIGH_sample", !eligible_pseudobulk_LOW ~ "LOW_below_cell_threshold",
      !eligible_pseudobulk_HIGH ~ "HIGH_below_cell_threshold", TRUE ~ "included")) %>% arrange(lineage_major,
    patient_id)
  write.csv(pair_audit, file.path(table_dir, "paired_HIGH_LOW_eligibility_audit.csv"), row.names = FALSE)
  lineage_pair_summary <- pair_audit %>% group_by(lineage_major) %>% summarise(total_matched_patients = sum(complete_pair),
    eligible_pairs = sum(eligible_pair), included_patient_ids = paste(patient_id[eligible_pair], collapse = ";"),
    excluded_patient_ids = paste(patient_id[!eligible_pair], collapse = ";"), .groups = "drop")
  write.csv(lineage_pair_summary, file.path(table_dir, "paired_HIGH_LOW_pair_summary.csv"), row.names = FALSE)
  p_cell_count <- ggplot(sample_lineage_audit, aes(x = sample_id, y = cells_used_for_DE, fill = tissue_type)) +
    geom_col() + geom_hline(yintercept = min_cells_per_pseudobulk, linetype = 2, colour = "black") + facet_wrap(~lineage_major,
    ncol = 1, scales = "free_y") + theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45,
    hjust = 1)) + labs(x = NULL, y = "Accepted cells per pseudobulk library", fill = "Group")
  ggsave(file.path(plot_dir, "sample_lineage_cells_used_for_DE.pdf"), p_cell_count, width = 13, height = 8)
  aggregate_sample_lineage <- function(object) {
    meta <- object[[]]
    keep <- meta$use_for_DE & as.character(meta$lineage_major) %in% valid_de_lineages
    meta <- meta[keep, , drop = FALSE]
    if (nrow(meta) == 0L)
      stop("No accepted MYELOID/LYMPHOID cells")
    meta$lineage_major <- as.character(meta$lineage_major)
    meta$pb_id <- paste(meta$sample_id, meta$lineage_major, sep = "___")
    pb_factor <- factor(meta$pb_id, levels = unique(meta$pb_id))
    aggregation_matrix <- Matrix::sparse.model.matrix(~0 + pb_factor)
    colnames(aggregation_matrix) <- levels(pb_factor)
    raw_counts <- GetAssayData(object, assay = "RNA", slot = "counts")
    if (!all(rownames(meta) %in% colnames(raw_counts))) {
      stop("Cell names in metadata and RNA counts do not match")
    }
    raw_counts <- raw_counts[, rownames(meta), drop = FALSE]
    pb_counts <- raw_counts %*% aggregation_matrix
    pb_meta <- meta %>% mutate(pb_id = as.character(pb_id)) %>% select(pb_id, sample_id, patient_id, tissue_type,
      dataset, lineage_major) %>% distinct()
    if (anyDuplicated(pb_meta$pb_id)) {
      stop("A pseudobulk ID maps to duplicated metadata")
    }
    pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$pb_id), , drop = FALSE]
    pb_meta$n_cells <- as.integer(table(pb_factor)[pb_meta$pb_id])
    rownames(pb_meta) <- pb_meta$pb_id
    if (!identical(colnames(pb_counts), rownames(pb_meta))) {
      stop("Pseudobulk count columns and metadata rows are not aligned")
    }
    list(counts = pb_counts, meta = pb_meta)
  }
  pb <- aggregate_sample_lineage(scRNA)
  write.csv(pb$meta, file.path(table_dir, "pseudobulk_library_metadata_all.csv"), row.names = FALSE)
  saveRDS(pb, file.path(rds_dir, "pseudobulk_raw_counts_and_metadata.rds"))
  add_threshold_columns <- function(tab) {
    tab %>% mutate(candidate = FDR < candidate_fdr & abs(logFC) >= candidate_abs_log2fc, conventional_DEG = FDR <
      conventional_fdr & abs(logFC) >= conventional_abs_log2fc)
  }
  fit_paired_lineage <- function(pb_counts, pb_meta, lineage_name) {
    valid_patients <- pair_audit %>% filter(lineage_major == lineage_name, eligible_pair) %>% pull(patient_id)
    if (length(valid_patients) < min_complete_pairs) {
      warning(lineage_name, ": only ", length(valid_patients), " eligible pairs; paired DE was skipped")
      return(NULL)
    }
    keep <- pb_meta$lineage_major == lineage_name & pb_meta$patient_id %in% valid_patients & pb_meta$tissue_type %in%
      c("LOW", "HIGH") & pb_meta$n_cells >= min_cells_per_pseudobulk
    counts <- as.matrix(pb_counts[, keep, drop = FALSE])
    meta <- pb_meta[keep, , drop = FALSE]
    meta$patient_id <- factor(meta$patient_id, levels = valid_patients)
    meta$tissue_type <- factor(meta$tissue_type, levels = c("LOW", "HIGH"))
    meta <- meta[order(meta$patient_id, meta$tissue_type), , drop = FALSE]
    counts <- counts[, rownames(meta), drop = FALSE]
    check_table <- table(meta$patient_id, meta$tissue_type)
    if (any(check_table != 1L)) {
      stop(lineage_name, ": each patient must have exactly one LOW and one HIGH")
    }
    design <- model.matrix(~patient_id + tissue_type, data = meta)
    if (qr(design)$rank < ncol(design)) {
      stop(lineage_name, ": paired design matrix is not full rank")
    }
    y <- DGEList(counts = counts)
    keep_gene <- filterByExpr(y, design = design)
    if (sum(keep_gene) < 100L) {
      warning(lineage_name, ": fewer than 100 genes passed filterByExpr")
    }
    if (!any(keep_gene))
      return(NULL)
    y <- y[keep_gene, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = "tissue_typeHIGH")
    tab <- topTags(qlf, n = Inf, sort.by = "P")$table %>% rownames_to_column("gene") %>% mutate(lineage = lineage_name,
      contrast = "HIGH_vs_LOW_paired", comparison_role = "PRIMARY_PAIRED", logFC_definition = "HIGH_minus_LOW",
      n_pairs = length(valid_patients), patient_ids = paste(valid_patients, collapse = ";")) %>% add_threshold_columns()
    saveRDS(list(y = y, design = design, fit = fit, qlf = qlf, metadata = meta), file.path(rds_dir, paste0("edgeR_paired_",
      lineage_name, ".rds")))
    logcpm <- cpm(y, log = TRUE, prior.count = 2)
    write.csv(logcpm, file.path(table_dir, paste0("paired_", lineage_name, "_TMM_logCPM.csv")))
    mds <- plotMDS(y, plot = FALSE)
    mds_df <- data.frame(sample_id = meta$sample_id, patient_id = as.character(meta$patient_id), tissue_type = as.character(meta$tissue_type),
      Dim1 = mds$x, Dim2 = mds$y, stringsAsFactors = FALSE)
    p_mds <- ggplot(mds_df, aes(Dim1, Dim2, colour = tissue_type)) + geom_point(size = 3) + geom_text(aes(label = patient_id),
      vjust = -0.7, size = 3) + theme_bw(base_size = 11) + labs(title = paste("Paired pseudobulk MDS:",
      lineage_name), colour = "Group")
    ggsave(file.path(plot_dir, paste0("paired_", lineage_name, "_MDS.pdf")), p_mds, width = 7, height = 6)
    tab
  }
  fit_contextual_lineage <- function(pb_counts, pb_meta, lineage_name, case_group, drop_control = NULL) {
    keep <- pb_meta$lineage_major == lineage_name & pb_meta$tissue_type %in% c("CONTROL", case_group) & pb_meta$n_cells >=
      min_cells_per_pseudobulk
    if (!is.null(drop_control)) {
      keep <- keep & pb_meta$sample_id != drop_control
    }
    counts <- as.matrix(pb_counts[, keep, drop = FALSE])
    meta <- pb_meta[keep, , drop = FALSE]
    group_n <- table(meta$tissue_type)
    required_groups <- c("CONTROL", case_group)
    if (!all(required_groups %in% names(group_n)) || any(group_n[required_groups] < min_samples_context_group)) {
      return(NULL)
    }
    meta$tissue_type <- factor(meta$tissue_type, levels = c("CONTROL", case_group))
    design <- model.matrix(~tissue_type, data = meta)
    y <- DGEList(counts = counts)
    keep_gene <- filterByExpr(y, design = design)
    if (!any(keep_gene))
      return(NULL)
    y <- y[keep_gene, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    coef_name <- paste0("tissue_type", case_group)
    qlf <- glmQLFTest(fit, coef = coef_name)
    tab <- topTags(qlf, n = Inf, sort.by = "P")$table %>% rownames_to_column("gene") %>% mutate(lineage = lineage_name,
      contrast = paste0(case_group, "_vs_CONTROL"), comparison_role = "SECONDARY_CONTEXTUAL_CROSS_DATASET",
      logFC_definition = paste0(case_group, "_minus_CONTROL"), n_case = unname(group_n[case_group]), n_control = unname(group_n["CONTROL"]),
      dropped_control = ifelse(is.null(drop_control), NA_character_, drop_control)) %>% add_threshold_columns()
    tab
  }
  paired_results <- lapply(valid_de_lineages, function(x) fit_paired_lineage(pb$counts, pb$meta, x))
  names(paired_results) <- valid_de_lineages
  paired_results <- paired_results[!vapply(paired_results, is.null, logical(1))]
  de_paired <- bind_rows(paired_results)
  if (nrow(de_paired) == 0L) {
    stop("Neither lineage had enough complete HIGH/LOW pairs. See ", file.path(table_dir, "paired_HIGH_LOW_eligibility_audit.csv"))
  }
  de_paired <- de_paired %>% mutate(direction = case_when(logFC > 0 ~ "HIGH_UP", logFC < 0 ~ "LOW_UP", TRUE ~
    "NO_CHANGE"))
  write.csv(de_paired, file.path(table_dir, "paired_HIGH_vs_LOW_all_tested_genes.csv"), row.names = FALSE)
  write.csv(filter(de_paired, candidate), file.path(table_dir, "paired_HIGH_vs_LOW_candidates_FDR0.10_logFC0.25.csv"),
    row.names = FALSE)
  write.csv(filter(de_paired, conventional_DEG), file.path(table_dir, "paired_HIGH_vs_LOW_conventional_FDR0.05_logFC0.50.csv"),
    row.names = FALSE)
  paired_summary <- de_paired %>% group_by(lineage, contrast, n_pairs, patient_ids) %>% summarise(tested_genes = n(),
    candidate_genes = sum(candidate), candidate_HIGH_up = sum(candidate & logFC > 0), candidate_LOW_up = sum(candidate &
      logFC < 0), conventional_DEGs = sum(conventional_DEG), conventional_HIGH_up = sum(conventional_DEG &
      logFC > 0), conventional_LOW_up = sum(conventional_DEG & logFC < 0), .groups = "drop")
  write.csv(paired_summary, file.path(table_dir, "paired_HIGH_vs_LOW_DE_summary.csv"), row.names = FALSE)
  for (lin in unique(de_paired$lineage)) {
    z <- de_paired %>% filter(lineage == lin) %>% mutate(plot_class = case_when(conventional_DEG & logFC >
      0 ~ "HIGH up (conventional)", conventional_DEG & logFC < 0 ~ "LOW up (conventional)", candidate & logFC >
      0 ~ "HIGH up (candidate)", candidate & logFC < 0 ~ "LOW up (candidate)", TRUE ~ "Not selected"), minus_log10_FDR = -log10(pmax(FDR,
      .Machine$double.xmin)))
    labels <- z %>% filter(candidate) %>% arrange(FDR, desc(abs(logFC))) %>% slice_head(n = 12L)
    p_volcano <- ggplot(z, aes(logFC, minus_log10_FDR, colour = plot_class)) + geom_point(size = 1.1,
      alpha = 0.65) + geom_vline(xintercept = c(-candidate_abs_log2fc, candidate_abs_log2fc),
      linetype = 2) + geom_hline(yintercept = -log10(candidate_fdr), linetype = 2) + geom_text(data = labels,
      aes(label = gene), size = 2.7, check_overlap = TRUE, vjust = -0.5, show.legend = FALSE) +
      scale_colour_manual(values = c(`HIGH up (conventional)` = "#B2182B", `LOW up (conventional)` = "#2166AC",
        `HIGH up (candidate)` = "#EF8A62", `LOW up (candidate)` = "#67A9CF", `Not selected` = "grey75")) +
      theme_bw(base_size = 11) + labs(title = paste("Paired HIGH vs LOW:", lin), x = "log2 fold change (HIGH - LOW)",
      y = expression(-log[10](FDR)), colour = NULL)
    ggsave(file.path(plot_dir, paste0("paired_HIGH_vs_LOW_volcano_", lin, ".pdf")), p_volcano, width = 8, height = 6.5)
  }
  de_contextual <- tibble(lineage = character(), contrast = character(), gene = character(), logFC = numeric(),
    FDR = numeric(), candidate = logical(), conventional_DEG = logical(), comparison_role = character(), n_case = integer(),
    n_control = integer(), direction = character())
  if (run_contextual_control_comparisons) {
    contextual_results <- list()
    for (lin in valid_de_lineages) {
      for (case_group in c("HIGH", "LOW")) {
        key <- paste(lin, case_group, sep = "___")
        contextual_results[[key]] <- fit_contextual_lineage(pb$counts, pb$meta, lin, case_group)
      }
    }
    contextual_results <- contextual_results[!vapply(contextual_results, is.null, logical(1))]
    contextual_bound <- bind_rows(contextual_results)
    if (nrow(contextual_bound) > 0L) {
      de_contextual <- contextual_bound %>% mutate(direction = case_when(contrast == "HIGH_vs_CONTROL" &
        logFC > 0 ~ "HIGH_UP", contrast == "HIGH_vs_CONTROL" & logFC < 0 ~ "CONTROL_UP", contrast == "LOW_vs_CONTROL" &
        logFC > 0 ~ "LOW_UP", contrast == "LOW_vs_CONTROL" & logFC < 0 ~ "CONTROL_UP", TRUE ~ "NO_CHANGE"))
      write.csv(de_contextual, file.path(table_dir, "contextual_HIGH_LOW_vs_CONTROL_all_tested_genes.csv"),
        row.names = FALSE)
      write.csv(filter(de_contextual, candidate), file.path(table_dir, "contextual_vs_CONTROL_candidates_FDR0.10_logFC0.25.csv"),
        row.names = FALSE)
      write.csv(filter(de_contextual, conventional_DEG), file.path(table_dir, "contextual_vs_CONTROL_conventional_FDR0.05_logFC0.50.csv"),
        row.names = FALSE)
      contextual_summary <- de_contextual %>% group_by(lineage, contrast, comparison_role, n_case, n_control) %>%
        summarise(tested_genes = n(), candidate_genes = sum(candidate), conventional_DEGs = sum(conventional_DEG),
         .groups = "drop")
      write.csv(contextual_summary, file.path(table_dir, "contextual_vs_CONTROL_DE_summary.csv"), row.names = FALSE)
    }
    else {
      warning("No contextual comparison met the minimum sample requirements")
    }
  }
  loo_summary <- tibble()
  if (run_contextual_control_comparisons && run_leave_one_control_out && nrow(de_contextual) > 0L) {
    target <- de_contextual %>% filter(candidate) %>% select(lineage, contrast, gene, full_logFC = logFC, full_FDR = FDR)
    loo_rows <- list()
    if (nrow(target) > 0L) {
      for (lin in valid_de_lineages) {
        control_ids <- pb$meta %>% filter(lineage_major == lin, tissue_type == "CONTROL", n_cells >= min_cells_per_pseudobulk) %>%
         pull(sample_id) %>% unique()
        for (case_group in c("HIGH", "LOW")) {
         contrast_name <- paste0(case_group, "_vs_CONTROL")
         target_now <- target %>% filter(lineage == lin, contrast == contrast_name)
         if (nrow(target_now) == 0L)
          next
         for (drop_id in control_ids) {
          z <- fit_contextual_lineage(pb$counts, pb$meta, lin, case_group, drop_control = drop_id)
          if (is.null(z))
           next
          z <- z %>% semi_join(target_now, by = c("lineage", "contrast", "gene"))
          if (nrow(z) > 0L) {
           loo_rows[[paste(lin, case_group, drop_id, sep = "___")]] <- z
          }
         }
        }
      }
    }
    loo_all <- bind_rows(loo_rows)
    if (nrow(loo_all) > 0L) {
      loo_summary <- loo_all %>% inner_join(target, by = c("lineage", "contrast", "gene")) %>% group_by(lineage,
        contrast, gene, full_logFC, full_FDR) %>% summarise(n_leave_one_out_runs = n(), direction_consistency = mean(sign(logFC) ==
        sign(full_logFC)), candidate_retention_fraction = mean(candidate), min_leave_one_out_logFC = min(logFC),
        max_leave_one_out_logFC = max(logFC), max_leave_one_out_FDR = max(FDR), .groups = "drop")
    }
    write.csv(loo_summary, file.path(table_dir, "contextual_leave_one_control_out_sensitivity.csv"), row.names = FALSE)
  }
  all_selected <- bind_rows(de_paired %>% select(lineage, contrast, gene, logFC, FDR, candidate), de_contextual %>%
    select(lineage, contrast, gene, logFC, FDR, candidate)) %>% filter(candidate) %>% mutate(direction_sign = ifelse(logFC >
    0, "UP", "DOWN"))
  cross_comparison_pattern <- all_selected %>% mutate(pattern_item = paste(contrast, direction_sign, sep = ":")) %>%
    group_by(lineage, gene) %>% summarise(n_candidate_comparisons = n_distinct(contrast), comparison_direction_pattern = paste(sort(unique(pattern_item)),
    collapse = ";"), .groups = "drop")
  write.csv(cross_comparison_pattern, file.path(table_dir, "candidate_cross_comparison_patterns.csv"), row.names = FALSE)
  saveRDS(scRNA, file.path(rds_dir, "scRNA_major_annotation_DE_audited.rds"))
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_final.txt"))
  message("Major-lineage annotation audit and pseudobulk DE completed.")
  message("Primary result: paired HIGH versus LOW within patient.")
  message("Positive paired logFC means higher expression in HIGH.")
  message("Control comparisons are contextual because dataset and group are confounded.")
  message("Output directory: ", outdir)
  invisible(NULL)
}

acti_stage_celltype <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(edgeR)
  })
  project_dir <- cfg$output_dir
  input_rds <- acti_input_object(cfg)
  outdir <- file.path(project_dir, "04_celltype_pseudobulk")
  table_dir <- file.path(outdir, "tables")
  rds_dir <- file.path(outdir, "RDS")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  min_cells_per_library <- 20L
  required_n_pairs <- 5L
  nominal_p_cutoff <- 0.05
  abs_log2fc_cutoff <- 0.2
  min_direction_fraction <- 0.8
  if (!file.exists(input_rds))
    stop("Input RDS not found: ", input_rds)
  message("Reading: ", input_rds)
  scRNA <- readRDS(input_rds)
  DefaultAssay(scRNA) <- "RNA"
  md <- scRNA[[]]
  required_columns <- c("sample_id", "patient_id", "tissue_type", "celltype_major", "lineage_major", "use_for_DE")
  missing_columns <- setdiff(required_columns, colnames(md))
  if (length(missing_columns) > 0L) {
    stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
  }
  md$sample_id <- as.character(md$sample_id)
  md$patient_id <- as.character(md$patient_id)
  md$tissue_type <- toupper(as.character(md$tissue_type))
  md$celltype_major <- as.character(md$celltype_major)
  md$lineage_major <- toupper(as.character(md$lineage_major))
  md$use_for_DE <- as.logical(md$use_for_DE)
  if (!all(colnames(scRNA) == rownames(md))) {
    stop("Seurat cell order does not match metadata row order")
  }
  if (!all(c("HIGH", "LOW") %in% unique(md$tissue_type))) {
    stop("Both HIGH and LOW labels are required in tissue_type")
  }
  immune_celltypes <- c("T cell", "NK cell", "Macrophage", "Monocyte", "B cell", "DC", "Mast cell")
  accepted <- !is.na(md$use_for_DE) & md$use_for_DE & md$tissue_type %in% c("HIGH", "LOW") & !is.na(md$celltype_major) &
    nzchar(md$celltype_major) & md$celltype_major %in% immune_celltypes
  audit_cells <- as.data.frame(xtabs(~celltype_major + tissue_type, data = md[accepted, , drop = FALSE]), stringsAsFactors = FALSE)
  colnames(audit_cells)[3] <- "n_cells"
  write.csv(audit_cells, file.path(table_dir, "01_cell_counts_by_celltype_and_group.csv"), row.names = FALSE)
  celltypes <- immune_celltypes[immune_celltypes %in% unique(md$celltype_major[accepted])]
  message("Cell types approved for DE: ", paste(celltypes, collapse = ", "))
  counts <- GetAssayData(scRNA, assay = "RNA", slot = "counts")
  counts <- counts[, accepted, drop = FALSE]
  meta_use <- md[accepted, , drop = FALSE]
  meta_use$pb_id <- paste(meta_use$sample_id, meta_use$celltype_major, sep = "___")
  pb_factor <- factor(meta_use$pb_id, levels = unique(meta_use$pb_id))
  design_sum <- Matrix::sparse.model.matrix(~0 + pb_factor)
  colnames(design_sum) <- levels(pb_factor)
  pb_counts <- counts %*% design_sum
  pb_counts <- as(pb_counts, "dgCMatrix")
  pb_meta <- unique(meta_use[, c("pb_id", "sample_id", "patient_id", "tissue_type", "celltype_major", "lineage_major")])
  pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$pb_id), , drop = FALSE]
  pb_meta$n_cells <- as.integer(table(pb_factor)[pb_meta$pb_id])
  rownames(pb_meta) <- pb_meta$pb_id
  if (anyNA(pb_meta$sample_id) || !identical(colnames(pb_counts), rownames(pb_meta))) {
    stop("Pseudobulk count columns do not match pseudobulk metadata")
  }
  write.csv(pb_meta, file.path(table_dir, "02_pseudobulk_library_metadata_all_celltypes.csv"), row.names = FALSE)
  saveRDS(list(counts = pb_counts, metadata = pb_meta), file.path(rds_dir, "celltype_pseudobulk_raw_counts_and_metadata.rds"))
  safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  analyse_one_celltype <- function(celltype_name) {
    m0 <- pb_meta[pb_meta$celltype_major == celltype_name & pb_meta$tissue_type %in% c("LOW", "HIGH"), , drop = FALSE]
    patient_group <- aggregate(n_cells ~ patient_id + tissue_type, data = m0, FUN = sum)
    patient_group$eligible <- patient_group$n_cells >= min_cells_per_library
    low_ok <- patient_group$patient_id[patient_group$tissue_type == "LOW" & patient_group$eligible]
    high_ok <- patient_group$patient_id[patient_group$tissue_type == "HIGH" & patient_group$eligible]
    paired_patients <- sort(intersect(low_ok, high_ok))
    eligibility_out <- merge(expand.grid(patient_id = sort(unique(m0$patient_id)), tissue_type = c("LOW", "HIGH"),
      stringsAsFactors = FALSE), patient_group, by = c("patient_id", "tissue_type"), all.x = TRUE)
    eligibility_out$n_cells[is.na(eligibility_out$n_cells)] <- 0L
    eligibility_out$eligible[is.na(eligibility_out$eligible)] <- FALSE
    eligibility_out$celltype_major <- celltype_name
    eligibility_out$paired_patient <- eligibility_out$patient_id %in% paired_patients
    status <- data.frame(celltype_major = celltype_name, n_paired_patients = length(paired_patients), paired_patient_ids = paste(paired_patients,
      collapse = ";"), min_cells_per_library = min_cells_per_library, required_n_pairs = required_n_pairs,
      analysis_status = if (length(paired_patients) >= required_n_pairs) {
        "analysed"
      }
      else {
        "not_analysed_insufficient_complete_pairs"
      }, tested_genes = NA_integer_, widest_candidates = NA_integer_, widest_HIGH_up = NA_integer_, widest_LOW_up = NA_integer_,
      conventional_DEGs = NA_integer_, conventional_HIGH_up = NA_integer_, conventional_LOW_up = NA_integer_,
      stringsAsFactors = FALSE)
    if (length(paired_patients) < required_n_pairs) {
      return(list(status = status, eligibility = eligibility_out))
    }
    keep_meta <- m0$patient_id %in% paired_patients & m0$n_cells >= min_cells_per_library
    m <- m0[keep_meta, , drop = FALSE]
    y_counts <- pb_counts[, rownames(m), drop = FALSE]
    pair_table <- table(m$patient_id, m$tissue_type)
    if (!all(pair_table[, c("LOW", "HIGH"), drop = FALSE] == 1L)) {
      stop("Unexpected duplicate/missing library for cell type: ", celltype_name)
    }
    m$patient_id <- factor(m$patient_id)
    m$tissue_type <- factor(m$tissue_type, levels = c("LOW", "HIGH"))
    design <- model.matrix(~patient_id + tissue_type, data = m)
    y <- DGEList(counts = y_counts)
    keep_gene <- filterByExpr(y, design = design)
    y <- y[keep_gene, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    coef_name <- grep("^tissue_typeHIGH$", colnames(design), value = TRUE)
    if (length(coef_name) != 1L) {
      stop("Could not identify HIGH-versus-LOW coefficient for ", celltype_name)
    }
    qlf <- glmQLFTest(fit, coef = match(coef_name, colnames(design)))
    tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
    tab$gene <- rownames(tab)
    rownames(tab) <- NULL
    lcpm <- cpm(y, log = TRUE, prior.count = 2)
    patient_diffs <- sapply(levels(m$patient_id), function(pid) {
      hi <- which(m$patient_id == pid & m$tissue_type == "HIGH")
      lo <- which(m$patient_id == pid & m$tissue_type == "LOW")
      lcpm[, hi] - lcpm[, lo]
    })
    if (is.null(dim(patient_diffs))) {
      patient_diffs <- matrix(patient_diffs, ncol = 1L)
    }
    rownames(patient_diffs) <- rownames(lcpm)
    colnames(patient_diffs) <- levels(m$patient_id)
    idx <- match(tab$gene, rownames(patient_diffs))
    d <- patient_diffs[idx, , drop = FALSE]
    expected_sign <- sign(tab$logFC)
    tab$same_direction_pairs <- rowSums(sign(d) == expected_sign, na.rm = TRUE)
    tab$n_complete_pairs_checked <- ncol(d)
    tab$direction_consistency <- tab$same_direction_pairs/tab$n_complete_pairs_checked
    tab$pairwise_median_logCPM_difference <- apply(d, 1, median, na.rm = TRUE)
    tab$pairwise_min_logCPM_difference <- apply(d, 1, min, na.rm = TRUE)
    tab$pairwise_max_logCPM_difference <- apply(d, 1, max, na.rm = TRUE)
    tab$celltype_major <- celltype_name
    tab$contrast <- "HIGH_vs_LOW_paired"
    tab$logFC_definition <- "HIGH_minus_LOW"
    tab$n_pairs <- length(paired_patients)
    tab$patient_ids <- paste(paired_patients, collapse = ";")
    tab$direction <- ifelse(tab$logFC > 0, "HIGH_UP", "LOW_UP")
    tab$widest_candidate <- tab$PValue < nominal_p_cutoff & abs(tab$logFC) >= abs_log2fc_cutoff & tab$direction_consistency >=
      min_direction_fraction
    tab$conventional_DEG <- tab$FDR < 0.05 & abs(tab$logFC) >= 0.5
    status$tested_genes <- nrow(tab)
    status$widest_candidates <- sum(tab$widest_candidate)
    status$widest_HIGH_up <- sum(tab$widest_candidate & tab$logFC > 0)
    status$widest_LOW_up <- sum(tab$widest_candidate & tab$logFC < 0)
    status$conventional_DEGs <- sum(tab$conventional_DEG)
    status$conventional_HIGH_up <- sum(tab$conventional_DEG & tab$logFC > 0)
    status$conventional_LOW_up <- sum(tab$conventional_DEG & tab$logFC < 0)
    saveRDS(list(celltype = celltype_name, metadata = m, design = design, dge = y, fit = fit, test = qlf, patient_logCPM_differences = patient_diffs,
      results = tab), file.path(rds_dir, paste0("edgeR_paired_", safe_name(celltype_name), ".rds")))
    list(status = status, eligibility = eligibility_out, results = tab)
  }
  results_list <- lapply(celltypes, analyse_one_celltype)
  names(results_list) <- celltypes
  status_all <- do.call(rbind, lapply(results_list, `[[`, "status"))
  eligibility_all <- do.call(rbind, lapply(results_list, `[[`, "eligibility"))
  analysed <- vapply(results_list, function(z) !is.null(z$results), logical(1))
  if (any(analysed)) {
    all_results <- do.call(rbind, lapply(results_list[analysed], `[[`, "results"))
    rownames(all_results) <- NULL
    widest_results <- all_results[all_results$widest_candidate, , drop = FALSE]
    conventional_results <- all_results[all_results$conventional_DEG, , drop = FALSE]
  }
  else {
    all_results <- data.frame()
    widest_results <- data.frame()
    conventional_results <- data.frame()
  }
  write.csv(eligibility_all, file.path(table_dir, "03_celltype_pair_eligibility_audit.csv"), row.names = FALSE)
  write.csv(status_all, file.path(table_dir, "04_celltype_DE_summary.csv"), row.names = FALSE)
  write.csv(all_results, file.path(table_dir, "05_celltype_HIGH_vs_LOW_all_tested_genes.csv"), row.names = FALSE)
  write.csv(widest_results, file.path(table_dir, "06_celltype_WIDEST_candidates.csv"), row.names = FALSE)
  write.csv(conventional_results, file.path(table_dir, "07_celltype_conventional_FDR0.05_logFC0.50.csv"), row.names = FALSE)
  threshold_record <- data.frame(parameter = c("primary_comparison", "experimental_unit", "minimum_cells_per_sample_celltype_library",
    "required_complete_patient_pairs", "nominal_PValue_strictly_less_than", "absolute_log2FC_greater_than_or_equal_to",
    "minimum_pair_direction_consistency", "candidate_interpretation"), value = c("HIGH_vs_LOW_paired", "patient",
    min_cells_per_library, required_n_pairs, nominal_p_cutoff, abs_log2fc_cutoff, min_direction_fraction, "exploratory nominally associated candidates; not FDR-significant DEGs"),
    stringsAsFactors = FALSE)
  write.csv(threshold_record, file.path(table_dir, "00_analysis_threshold_record.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
  message("Completed. Output folder: ", outdir)
  message("First inspect: tables/04_celltype_DE_summary.csv")
  write.csv(data.frame(gene = sort(unique(widest_results$gene))), file.path(table_dir, "08_celltype_candidate_gene_union.csv"),
    row.names = FALSE)
  invisible(NULL)
}

acti_stage_harmony <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(RANN)
    library(cluster)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
  })
  set.seed(20260826)
  options(stringsAsFactors = FALSE)
  project_dir <- cfg$output_dir
  input_rds <- acti_input_object(cfg)
  outdir <- file.path(project_dir, "05_harmony_audit")
  table_dir <- file.path(outdir, "tables")
  plot_dir <- file.path(outdir, "plots")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  n_dims <- 30L
  knn_k <- 30L
  max_cells_per_sample_celltype <- 300L
  min_cells_per_sample_celltype <- 20L
  min_samples_per_stratum <- 2L
  max_cells_per_celltype_for_bio_silhouette <- 200L
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_start.txt"))
  scRNA <- readRDS(input_rds)
  md <- scRNA[[]]
  if (!"celltype_major" %in% colnames(md) && "final_celltype" %in% colnames(md)) {
    scRNA$celltype_major <- md$final_celltype
  }
  if (!"lineage_major" %in% colnames(md) && "lineage" %in% colnames(md)) {
    scRNA$lineage_major <- md$lineage
  }
  if (!"annotation_review_status" %in% colnames(md) && "review_status" %in% colnames(md)) {
    scRNA$annotation_review_status <- md$review_status
  }
  md <- scRNA[[]]
  required_metadata <- c("sample_id", "tissue_type", "celltype_major", "seurat_clusters")
  missing_metadata <- setdiff(required_metadata, colnames(md))
  if (length(missing_metadata) > 0L) {
    stop("Missing metadata: ", paste(missing_metadata, collapse = ", "))
  }
  if (!all(c("pca", "harmony") %in% Reductions(scRNA))) {
    stop("Both pca and harmony reductions must be present in the Seurat object")
  }
  pca_all <- Embeddings(scRNA, "pca")
  harmony_all <- Embeddings(scRNA, "harmony")
  if (ncol(pca_all) < n_dims || ncol(harmony_all) < n_dims) {
    stop("PCA and Harmony must each contain at least ", n_dims, " dimensions")
  }
  if (!identical(rownames(pca_all), rownames(harmony_all))) {
    stop("PCA and Harmony cell order differs")
  }
  md$sample_id <- as.character(md$sample_id)
  md$tissue_type <- toupper(as.character(md$tissue_type))
  md$celltype_major <- as.character(md$celltype_major)
  md$seurat_clusters <- as.character(md$seurat_clusters)
  keep_annotation <- !is.na(md$celltype_major) & trimws(md$celltype_major) != "" & md$celltype_major != "Other"
  if ("annotation_review_status" %in% colnames(md)) {
    review_status <- as.character(md$annotation_review_status)
    if (any(!is.na(review_status))) {
      keep_annotation <- keep_annotation & review_status %in% c("final", "corrected")
    }
  }
  md_keep <- md[keep_annotation, , drop = FALSE]
  if (nrow(md_keep) == 0L)
    stop("No high-confidence annotated cells remain")
  balanced_sample <- function(cell_names, group_labels, maximum_n, minimum_group_n = 1L, seed = 20260826L) {
    stopifnot(length(cell_names) == length(group_labels))
    split_cells <- split(cell_names, group_labels)
    split_cells <- split_cells[lengths(split_cells) >= minimum_group_n]
    set.seed(seed)
    selected <- unlist(lapply(split_cells, function(x) {
      if (length(x) <= maximum_n)
        x
      else sample(x, maximum_n)
    }), use.names = FALSE)
    unique(selected)
  }
  sample_celltype_group <- paste(md_keep$sample_id, md_keep$celltype_major, sep = "___")
  selected_cells <- balanced_sample(rownames(md_keep), sample_celltype_group, maximum_n = max_cells_per_sample_celltype,
    minimum_group_n = min_cells_per_sample_celltype)
  md_eval <- md_keep[selected_cells, , drop = FALSE]
  pca_eval <- pca_all[selected_cells, seq_len(n_dims), drop = FALSE]
  harmony_eval <- harmony_all[selected_cells, seq_len(n_dims), drop = FALSE]
  if (!identical(rownames(md_eval), rownames(pca_eval)) || !identical(rownames(md_eval), rownames(harmony_eval))) {
    stop("Selected metadata and embedding rows are not aligned")
  }
  selection_audit <- md_eval %>% count(sample_id, tissue_type, celltype_major, name = "selected_cells") %>% arrange(celltype_major,
    tissue_type, sample_id)
  write.csv(selection_audit, file.path(table_dir, "balanced_cell_selection_audit.csv"), row.names = FALSE)
  local_sample_metrics <- function(embedding, sample_label, k = 30L) {
    n <- nrow(embedding)
    sample_label <- as.character(sample_label)
    n_samples <- length(unique(sample_label))
    if (n_samples < 2L || n < 3L)
      return(NULL)
    k_use <- min(k, n - 1L)
    nn <- RANN::nn2(embedding, k = k_use + 1L)$nn.idx
    inverse_simpson <- numeric(n)
    other_sample_fraction <- numeric(n)
    for (i in seq_len(n)) {
      neighbours <- nn[i, ]
      neighbours <- neighbours[neighbours != i]
      neighbours <- head(neighbours, k_use)
      neighbour_samples <- sample_label[neighbours]
      proportions <- table(neighbour_samples)/length(neighbour_samples)
      inverse_simpson[i] <- 1/sum(proportions^2)
      other_sample_fraction[i] <- mean(neighbour_samples != sample_label[i])
    }
    normalized_inverse_simpson <- (inverse_simpson - 1)/(n_samples - 1)
    expected_other_sample_fraction <- 1 - sum((table(sample_label)/length(sample_label))^2)
    normalized_other_sample_mixing <- other_sample_fraction/expected_other_sample_fraction
    data.frame(normalized_sample_inverse_simpson = normalized_inverse_simpson, other_sample_neighbour_fraction = other_sample_fraction,
      normalized_other_sample_mixing = normalized_other_sample_mixing, stringsAsFactors = FALSE)
  }
  sample_structure_metrics <- function(embedding, sample_label) {
    sample_label <- factor(sample_label)
    if (nlevels(sample_label) < 2L || nrow(embedding) <= nlevels(sample_label)) {
      return(c(sample_ASW = NA_real_, sample_PCR_R2 = NA_real_))
    }
    sample_asw <- mean(cluster::silhouette(as.integer(sample_label), dist(embedding))[, "sil_width"])
    design <- model.matrix(~0 + sample_label)
    fit <- lm.fit(design, embedding)
    sse <- colSums(fit$residuals^2)
    centred <- sweep(embedding, 2L, colMeans(embedding), "-")
    sst <- colSums(centred^2)
    valid <- is.finite(sst) & sst > 0
    r2 <- 1 - sse[valid]/sst[valid]
    weighted_r2 <- sum(r2 * sst[valid])/sum(sst[valid])
    c(sample_ASW = sample_asw, sample_PCR_R2 = weighted_r2)
  }
  stratum_id <- paste(md_eval$celltype_major, md_eval$tissue_type, sep = "___")
  stratum_cells <- split(rownames(md_eval), stratum_id)
  local_cell_results <- list()
  stratum_results <- list()
  for (st in names(stratum_cells)) {
    cells <- stratum_cells[[st]]
    meta_st <- md_eval[cells, , drop = FALSE]
    n_samples_st <- length(unique(meta_st$sample_id))
    if (n_samples_st < min_samples_per_stratum)
      next
    for (method_name in c("Before Harmony (PCA)", "After Harmony")) {
      embedding_st <- if (method_name == "Before Harmony (PCA)") {
        pca_eval[cells, , drop = FALSE]
      }
      else {
        harmony_eval[cells, , drop = FALSE]
      }
      mixing <- local_sample_metrics(embedding_st, meta_st$sample_id, k = knn_k)
      if (is.null(mixing))
        next
      mixing$cell <- cells
      mixing$sample_id <- meta_st[cells, "sample_id"]
      mixing$tissue_type <- meta_st[cells, "tissue_type"]
      mixing$celltype_major <- meta_st[cells, "celltype_major"]
      mixing$stratum <- st
      mixing$method <- method_name
      local_cell_results[[paste(st, method_name, sep = "___")]] <- mixing
      structure <- sample_structure_metrics(embedding_st, meta_st$sample_id)
      stratum_results[[paste(st, method_name, sep = "___")]] <- data.frame(stratum = st, celltype_major = unique(meta_st$celltype_major),
        tissue_type = unique(meta_st$tissue_type), method = method_name, n_cells = nrow(meta_st), n_samples = n_samples_st,
        median_normalized_sample_inverse_simpson = median(mixing$normalized_sample_inverse_simpson, na.rm = TRUE),
        median_other_sample_neighbour_fraction = median(mixing$other_sample_neighbour_fraction, na.rm = TRUE),
        median_normalized_other_sample_mixing = median(mixing$normalized_other_sample_mixing, na.rm = TRUE),
        sample_ASW = unname(structure["sample_ASW"]), sample_PCR_R2 = unname(structure["sample_PCR_R2"]),
        stringsAsFactors = FALSE)
    }
  }
  local_cell_metrics <- bind_rows(local_cell_results)
  metrics_by_stratum <- bind_rows(stratum_results)
  if (nrow(metrics_by_stratum) == 0L) {
    stop("No celltype x tissue stratum contained at least two eligible samples")
  }
  write.csv(local_cell_metrics, file.path(table_dir, "local_sample_mixing_per_cell.csv"), row.names = FALSE)
  write.csv(metrics_by_stratum, file.path(table_dir, "batch_metrics_by_celltype_tissue_stratum.csv"), row.names = FALSE)
  biological_structure_metrics <- function(embedding, celltype_label, k = 30L, silhouette_max_per_type = 200L) {
    celltype_label <- as.character(celltype_label)
    n <- nrow(embedding)
    k_use <- min(k, n - 1L)
    nn <- RANN::nn2(embedding, k = k_use + 1L)$nn.idx
    celltype_inverse_simpson <- numeric(n)
    same_celltype_neighbour_fraction <- numeric(n)
    for (i in seq_len(n)) {
      neighbours <- nn[i, ]
      neighbours <- neighbours[neighbours != i]
      neighbours <- head(neighbours, k_use)
      labels_now <- celltype_label[neighbours]
      proportions <- table(labels_now)/length(labels_now)
      celltype_inverse_simpson[i] <- 1/sum(proportions^2)
      same_celltype_neighbour_fraction[i] <- mean(labels_now == celltype_label[i])
    }
    silhouette_cells <- balanced_sample(rownames(embedding), celltype_label, maximum_n = silhouette_max_per_type,
      minimum_group_n = 2L, seed = 20260827L)
    embedding_sil <- embedding[silhouette_cells, , drop = FALSE]
    labels_sil <- factor(celltype_label[match(silhouette_cells, rownames(embedding))])
    celltype_asw <- mean(cluster::silhouette(as.integer(labels_sil), dist(embedding_sil))[, "sil_width"])
    list(per_cell = data.frame(cell = rownames(embedding), celltype_major = celltype_label, celltype_inverse_simpson = celltype_inverse_simpson,
      same_celltype_neighbour_fraction = same_celltype_neighbour_fraction, stringsAsFactors = FALSE), summary = c(median_celltype_inverse_simpson = median(celltype_inverse_simpson),
      median_same_celltype_neighbour_fraction = median(same_celltype_neighbour_fraction), celltype_ASW = celltype_asw))
  }
  bio_before <- biological_structure_metrics(pca_eval, md_eval$celltype_major, k = knn_k, silhouette_max_per_type = max_cells_per_celltype_for_bio_silhouette)
  bio_after <- biological_structure_metrics(harmony_eval, md_eval$celltype_major, k = knn_k, silhouette_max_per_type = max_cells_per_celltype_for_bio_silhouette)
  bio_before$per_cell$method <- "Before Harmony (PCA)"
  bio_after$per_cell$method <- "After Harmony"
  biological_per_cell <- bind_rows(bio_before$per_cell, bio_after$per_cell)
  biological_summary <- bind_rows(as.data.frame(as.list(bio_before$summary)) %>% mutate(method = "Before Harmony (PCA)"),
    as.data.frame(as.list(bio_after$summary)) %>% mutate(method = "After Harmony")) %>% select(method, everything())
  write.csv(biological_per_cell, file.path(table_dir, "biological_structure_per_cell.csv"), row.names = FALSE)
  write.csv(biological_summary, file.path(table_dir, "biological_structure_summary.csv"), row.names = FALSE)
  metric_definitions <- tibble(metric = c("median_normalized_sample_inverse_simpson", "median_other_sample_neighbour_fraction",
    "median_normalized_other_sample_mixing", "sample_ASW", "sample_PCR_R2"), favorable_change = c("increase",
    "increase", "increase", "decrease", "decrease"), interpretation = c("Local diversity of sample identities within the same cell type and tissue group",
    "Fraction of neighbours originating from another sample", "Observed other-sample mixing divided by random-mixing expectation",
    "Separation of cells by sample identity; lower is less sample-associated", "Embedding variance explained by sample identity; lower is less sample-associated"))
  metrics_long <- metrics_by_stratum %>% pivot_longer(cols = all_of(metric_definitions$metric), names_to = "metric",
    values_to = "value")
  metric_paired_wide <- metrics_long %>% select(stratum, celltype_major, tissue_type, metric, method, value) %>%
    pivot_wider(names_from = method, values_from = value)
  names(metric_paired_wide)[names(metric_paired_wide) == "Before Harmony (PCA)"] <- "before"
  names(metric_paired_wide)[names(metric_paired_wide) == "After Harmony"] <- "after"
  metric_summary <- metric_paired_wide %>% group_by(metric) %>% summarise(n_strata = sum(is.finite(before) &
    is.finite(after)), before_median = median(before, na.rm = TRUE), after_median = median(after, na.rm = TRUE),
    median_paired_change = median(after - before, na.rm = TRUE), paired_wilcoxon_PValue = ifelse(n_strata >=
      3L, suppressWarnings(wilcox.test(after, before, paired = TRUE, exact = FALSE)$p.value), NA_real_),
    .groups = "drop") %>% left_join(metric_definitions, by = "metric")
  write.csv(metric_paired_wide, file.path(table_dir, "batch_metrics_paired_before_after.csv"), row.names = FALSE)
  write.csv(metric_summary, file.path(table_dir, "Harmony_quantitative_metric_summary.csv"), row.names = FALSE)
  method_levels <- c("Before Harmony (PCA)", "After Harmony")
  metrics_long$method <- factor(metrics_long$method, levels = method_levels)
  palette_method <- c(`Before Harmony (PCA)` = "#9E9E9E", `After Harmony` = "#2C7FB8")
  theme_batch_eval <- function(base_size = 7) {
    theme_classic(base_size = base_size, base_family = "sans") + theme(axis.line = element_line(size = 0.35,
      colour = "black"), axis.ticks = element_line(size = 0.35, colour = "black"), axis.text = element_text(colour = "black"),
      strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
      legend.position = "none", panel.grid = element_blank())
  }
  paired_metric_plot <- function(metric_name, title_text, y_text) {
    z <- metrics_long %>% filter(metric == metric_name)
    ggplot(z, aes(method, value, group = stratum)) + geom_line(colour = "#BDBDBD", size = 0.35,
      alpha = 0.75) + geom_point(aes(colour = method), size = 1.5, alpha = 0.9) + scale_colour_manual(values = palette_method) +
      labs(x = NULL, y = y_text, title = title_text) + theme_batch_eval() + theme(axis.text.x = element_text(angle = 18,
      hjust = 1))
  }
  p_a <- paired_metric_plot("median_normalized_sample_inverse_simpson", "Local sample diversity", "Normalized inverse Simpson")
  p_b <- paired_metric_plot("sample_ASW", "Sample-associated separation", "Sample silhouette width")
  p_c <- paired_metric_plot("sample_PCR_R2", "Variance explained by sample", "Weighted sample R2")
  bio_plot_data <- biological_summary %>% pivot_longer(cols = c(median_celltype_inverse_simpson, median_same_celltype_neighbour_fraction,
    celltype_ASW), names_to = "metric", values_to = "value") %>% mutate(metric = recode(metric, median_celltype_inverse_simpson = "Cell-type inverse Simpson (lower better)",
    median_same_celltype_neighbour_fraction = "Same-cell-type neighbours (higher better)", celltype_ASW = "Cell-type silhouette (higher better)"),
    method = factor(method, levels = method_levels))
  p_d <- ggplot(bio_plot_data, aes(method, value, group = metric)) + geom_line(colour = "#BDBDBD", size = 0.4) +
    geom_point(aes(colour = method), size = 1.8) + facet_wrap(~metric, scales = "free_y", nrow = 1) + scale_colour_manual(values = palette_method) +
    labs(x = NULL, y = "Metric value", title = "Broad cell-class preservation") + theme_batch_eval() + theme(axis.text.x = element_text(angle = 18,
    hjust = 1))
  figure_quantitative <- (p_a | p_b | p_c)/p_d + plot_layout(heights = c(1, 0.9)) + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  figure_base <- file.path(plot_dir, "Harmony_quantitative_batch_evaluation")
  grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = 183/25.4, height = 130/25.4,
    family = "sans")
  print(figure_quantitative)
  dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(figure_base, ".svg"), width = 183/25.4, height = 130/25.4)
    print(figure_quantitative)
    dev.off()
  }
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(paste0(figure_base, ".tiff"), width = 183/25.4, height = 130/25.4,
      units = "in", res = 600, compression = "lzw")
    print(figure_quantitative)
    dev.off()
  }
  methods_record <- c("Harmony batch-correction evaluation", paste0("Input object: ", basename(input_rds)), paste0("Dimensions compared: PCA 1-",
    n_dims, " versus Harmony 1-", n_dims), paste0("Neighbourhood size: k=", knn_k), paste0("Maximum cells per sample x cell type: ",
    max_cells_per_sample_celltype), paste0("Minimum cells per sample x cell type: ", min_cells_per_sample_celltype),
    "Sample-mixing metrics were calculated within broad cell type x tissue group strata.", "HIGH, LOW, and CONTROL cells were not expected or forced to mix across biological groups.",
    "The same selected cells were used before and after Harmony.", "Statistical comparison used paired celltype x tissue strata, not cells as independent replicates.",
    "Lower sample ASW and sample PCR R2 indicate reduced sample-associated structure.", "Higher normalized sample inverse-Simpson diversity indicates improved local sample mixing.",
    "Cell-type ASW, same-cell-type neighbour fraction, and cell-type inverse Simpson assess biological structure preservation.",
    "These metrics support but do not prove complete removal of all technical effects.")
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_final.txt"))
  message("Harmony quantitative evaluation completed.")
  message("Primary table: ", file.path(table_dir, "Harmony_quantitative_metric_summary.csv"))
  message("Figure: ", paste0(figure_base, ".pdf"))
  message("Interpret mixing and biological-conservation metrics together.")
  invisible(NULL)
}

acti_stage_figure2 <- function(cfg) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(RColorBrewer)
    library(scattermore)
    library(ggrepel)
    library(edgeR)
    library(svglite)
    library(ragg)
  })
  set.seed(20261008)
  script_args <- commandArgs(trailingOnly = TRUE)
  figure2_run_mode <- toupper(cfg$figure2_mode)
  only_bc_update <- "--only-bc" %in% script_args || figure2_run_mode == "BC"
  only_abcd_update <- "--only-abcd" %in% script_args || figure2_run_mode == "ABCD"
  only_e_update <- "--only-e" %in% script_args || figure2_run_mode == "E"
  only_f_update <- any(c("--only-f", "--only-h") %in% script_args) || figure2_run_mode %in% c("F", "H")
  panel_d_width_mm <- 220
  panel_d_height_mm <- 90
  panel_e_width_mm <- 90
  panel_e_height_mm <- 90
  revision_dir <- cfg$output_dir
  input_rds <- acti_input_object(cfg)
  celltype_deg_file <- file.path(revision_dir, "04_celltype_pseudobulk", "tables", "06_celltype_WIDEST_candidates.csv")
  lineage_deg_file <- file.path(revision_dir, "03_lineage_pseudobulk", "tables", "paired_HIGH_vs_LOW_all_tested_genes.csv")
  celltype_edger_rds_dir <- file.path(revision_dir, "04_celltype_pseudobulk", "RDS")
  outdir <- file.path(revision_dir, "06_figure2")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "source_data"), showWarnings = FALSE)
  stopifnot(file.exists(input_rds), file.exists(celltype_deg_file))
  group_levels <- c("HD", "High", "Low")
  group_cols <- c(HD = "#4C78A8", High = "#D55E00", Low = "#009E73")
  celltype_levels <- c("T cell", "NK cell", "B cell", "Monocyte", "Macrophage", "DC", "Mast cell", "Epithelial",
    "Endothelial", "Smooth muscle", "Other")
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", , drop = FALSE]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals), SIMPLIFY = FALSE,
    USE.NAMES = FALSE))
  celltype_cols <- setNames(col_vector[seq_along(celltype_levels)], celltype_levels)
  lineage_cols <- c(Other = "#BDBDBD", Lymphoid = "#377EB8", Myeloid = "#E41A1C")
  marker_list <- list(`T cell` = c("CD3D", "TRAC"), `NK cell` = c("NKG7", "GNLY"), `B cell` = c("MS4A1", "CD79A"),
    Monocyte = c("S100A8", "FCN1"), Macrophage = c("C1QB", "MARCO"), DC = c("IDO1", "DNASE1L3"), `Mast cell` = c("TPSAB1",
      "KIT"), Epithelial = c("EPCAM", "KRT19"), Endothelial = c("PECAM1", "CLDN5"), `Smooth muscle` = c("ACTA2",
      "TAGLN"))
  base_theme <- theme_classic(base_size = 9, base_family = "Arial") + theme(plot.title = element_text(size = 10,
    face = "bold", hjust = 0.5), axis.title = element_text(size = 8.5), axis.text = element_text(size = 7.5,
    colour = "black"), legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7),
    strip.text = element_text(size = 8.5, face = "bold"), strip.background = element_rect(fill = "white", colour = "black",
      linewidth = 0.35), plot.margin = margin(4, 4, 4, 4))
  place_export_file <- function(temp_file, requested_file) {
    final_file <- requested_file
    if (file.exists(final_file) && !file.remove(final_file)) {
      suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
      ext <- tools::file_ext(final_file)
      stem <- tools::file_path_sans_ext(basename(final_file))
      final_file <- file.path(dirname(final_file), paste0(stem, "_new_", suffix, ".", ext))
      warning("The existing output appears to be open or locked. Wrote a versioned ", "file instead: ", final_file)
    }
    if (!file.rename(temp_file, final_file)) {
      copied <- file.copy(temp_file, final_file, overwrite = TRUE)
      if (!copied)
        stop("Could not place exported file: ", final_file)
      unlink(temp_file)
    }
    normalizePath(final_file, winslash = "/", mustWork = TRUE)
  }
  close_new_graphics_devices <- function(start_device) {
    while (grDevices::dev.cur() != start_device && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  }
  remove_plot_titles <- function(plot) {
    plot + labs(title = NULL, subtitle = NULL, caption = NULL) + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank())
  }
  clean_previous_panel_outputs <- function(mode = "ALL") invisible(NULL)
  export_panel <- function(plot, stem, width_mm, height_mm) {
    pdf_file <- file.path(outdir, paste0(stem, ".pdf"))
    svg_file <- file.path(outdir, paste0(stem, ".svg"))
    tif_file <- file.path(outdir, paste0(stem, ".tiff"))
    png_file <- file.path(outdir, paste0(stem, ".png"))
    width_in <- width_mm/25.4
    height_in <- height_mm/25.4
    pdf_temp <- tempfile(pattern = paste0(stem, "_"), tmpdir = outdir, fileext = ".pdf")
    start_device <- grDevices::dev.cur()
    pdf_method <- "cairo_pdf"
    cairo_message <- NULL
    cairo_ok <- tryCatch({
      grDevices::cairo_pdf(pdf_temp, width = width_in, height = height_in, family = "Arial", bg = "white")
      print(plot)
      grDevices::dev.off()
      TRUE
    }, error = function(e) {
      cairo_message <<- conditionMessage(e)
      close_new_graphics_devices(start_device)
      FALSE
    })
    if (!cairo_ok) {
      pdf_method <- "base_pdf_fallback"
      warning("cairo_pdf failed; using base PDF fallback. Reason: ", cairo_message)
      fallback_plot <- plot + theme(text = element_text(family = "sans"))
      grDevices::pdf(pdf_temp, width = width_in, height = height_in, family = "sans", useDingbats = FALSE,
        bg = "white")
      print(fallback_plot)
      grDevices::dev.off()
    }
    pdf_file <- place_export_file(pdf_temp, pdf_file)
    svg_temp <- tempfile(pattern = paste0(stem, "_"), tmpdir = outdir, fileext = ".svg")
    ggsave(svg_temp, plot = plot, width = width_mm, height = height_mm, units = "mm", device = svglite::svglite,
      limitsize = FALSE)
    svg_file <- place_export_file(svg_temp, svg_file)
    tif_temp <- tempfile(pattern = paste0(stem, "_"), tmpdir = outdir, fileext = ".tiff")
    ragg::agg_tiff(tif_temp, width = width_mm, height = height_mm, units = "mm", res = 600, compression = "lzw",
      background = "white")
    print(plot)
    grDevices::dev.off()
    tif_file <- place_export_file(tif_temp, tif_file)
    png_temp <- tempfile(pattern = paste0(stem, "_"), tmpdir = outdir, fileext = ".png")
    ragg::agg_png(png_temp, width = width_mm, height = height_mm, units = "mm", res = 300, background = "white")
    print(plot)
    grDevices::dev.off()
    png_file <- place_export_file(png_temp, png_file)
    data.frame(panel = stem, width_mm = width_mm, height_mm = height_mm, pdf_method = pdf_method, pdf = pdf_file,
      svg = svg_file, tiff = tif_file, png = png_file, stringsAsFactors = FALSE)
  }
  umap_base <- function(data, colour_var, colours, title, show_legend = TRUE, facet_var = NULL, point_size = 0.45) {
    p <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, colour = .data[[colour_var]])) + scattermore::geom_scattermore(pointsize = point_size,
      pixels = c(1400, 1100)) + scale_colour_manual(values = colours, drop = FALSE) + coord_equal() + labs(title = title,
      x = "UMAP 1", y = "UMAP 2", colour = NULL) + base_theme + theme(axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35), legend.key.height = unit(3.6,
        "mm"))
    if (!is.null(facet_var)) {
      p <- p + facet_wrap(stats::as.formula(paste("~", facet_var)), nrow = 1)
    }
    if (!show_legend)
      p <- p + theme(legend.position = "none")
    p
  }
  dot_summary <- function(object, features, group_vector, group_levels) {
    features <- features[features %in% rownames(object)]
    if (!length(features))
      stop("None of the requested dot-plot genes were found.")
    mat <- Seurat::GetAssayData(object, assay = "RNA", slot = "data")[features, , drop = FALSE]
    group_vector <- factor(group_vector, levels = group_levels)
    ans <- lapply(group_levels, function(g) {
      idx <- which(group_vector == g)
      if (!length(idx))
        return(NULL)
      submat <- mat[, idx, drop = FALSE]
      data.frame(gene = features, group = g, avg_expression = as.numeric(Matrix::rowMeans(submat)), percent_expressed = as.numeric(Matrix::rowMeans(submat >
        0)) * 100, n_cells = length(idx), stringsAsFactors = FALSE)
    })
    ans <- bind_rows(ans)
    ans <- ans %>% group_by(gene) %>% mutate(scaled_average = if (sd(avg_expression) > 0) {
      as.numeric(scale(avg_expression))
    }
    else 0, scaled_average = pmax(-2, pmin(2, scaled_average))) %>% ungroup()
    ans
  }
  message("Loading fixed audited Seurat object: ", input_rds)
  obj <- readRDS(input_rds)
  DefaultAssay(obj) <- "RNA"
  required_meta <- c("sample_id", "patient_id", "tissue_type", "celltype_major", "lineage_major", "dataset",
    "use_for_DE")
  missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
  if (length(missing_meta))
    stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
  if (!"umap" %in% names(obj@reductions))
    stop("The fixed object has no UMAP reduction.")
  meta <- obj@meta.data
  meta$cell_barcode <- rownames(meta)
  meta$figure_group <- dplyr::recode(meta$tissue_type, CONTROL = "HD", HIGH = "High", LOW = "Low")
  meta$figure_group <- factor(meta$figure_group, levels = group_levels)
  meta$celltype_major <- factor(meta$celltype_major, levels = celltype_levels)
  meta$lineage_major <- factor(meta$lineage_major, levels = c("OTHER", "LYMPHOID", "MYELOID"))
  sample_map <- meta %>% distinct(sample_id, figure_group) %>% arrange(figure_group, sample_id) %>% group_by(figure_group) %>%
    mutate(display_sample = paste0(as.character(figure_group), row_number())) %>% ungroup()
  sample_levels <- sample_map$display_sample
  meta <- meta %>% left_join(sample_map, by = c("sample_id", "figure_group"))
  lesion_patient_ids <- meta %>% filter(figure_group %in% c("High", "Low")) %>% distinct(patient_id) %>% arrange(patient_id) %>%
    pull(patient_id)
  hd_patient_ids <- meta %>% filter(figure_group == "HD") %>% distinct(patient_id) %>% arrange(patient_id) %>%
    pull(patient_id)
  donor_key <- bind_rows(data.frame(patient_id = lesion_patient_ids, display_donor = paste0("Pt", seq_along(lesion_patient_ids)),
    donor_class = "Paired lesion patient", stringsAsFactors = FALSE), data.frame(patient_id = hd_patient_ids,
    display_donor = paste0("HD", seq_along(hd_patient_ids)), donor_class = "Independent HD donor", stringsAsFactors = FALSE))
  meta <- meta %>% left_join(donor_key, by = "patient_id")
  umap <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  umap$cell_barcode <- rownames(umap)
  plot_df <- left_join(umap, meta, by = "cell_barcode")
  if (anyNA(plot_df$figure_group) || anyNA(plot_df$celltype_major)) {
    stop("Group or major cell-type annotation contains missing values.")
  }
  plot_df <- plot_df[sample(seq_len(nrow(plot_df))), , drop = FALSE]
  write.csv(sample_map, file.path(outdir, "source_data", "sample_display_mapping.csv"), row.names = FALSE)
  write.csv(donor_key, file.path(outdir, "source_data", "donor_display_mapping.csv"), row.names = FALSE)
  write.csv(meta %>% count(figure_group, sample_id, display_sample, celltype_major, lineage_major, name = "n_cells"),
    file.path(outdir, "source_data", "cell_counts_by_sample_annotation.csv"), row.names = FALSE)
  pA <- umap_base(plot_df, "figure_group", group_cols, "Sample origin", show_legend = TRUE, point_size = 0.68) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
  donor_levels <- donor_key$display_donor
  donor_cols <- setNames(col_vector[seq_along(donor_levels)], donor_levels)
  plot_df$display_donor <- factor(plot_df$display_donor, levels = donor_levels)
  pB <- umap_base(plot_df, "display_donor", donor_cols, "Donor", show_legend = TRUE, point_size = 1.2) + guides(colour = guide_legend(ncol = 2,
    override.aes = list(size = 2, alpha = 1))) + theme(legend.key.height = unit(2.6, "mm"),
    legend.text = element_text(size = 6.5))
  pC <- umap_base(plot_df, "celltype_major", celltype_cols, "Cell type", show_legend = TRUE, point_size = 1.2) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.2, alpha = 1))) + theme(legend.key.height = unit(3,
    "mm"), legend.text = element_text(size = 6.6))
  dotplot_celltype_levels <- setdiff(celltype_levels, "Other")
  markers <- unique(unlist(marker_list, use.names = FALSE))
  missing_markers <- setdiff(markers, rownames(obj))
  if (length(missing_markers))
    warning("Markers absent from object: ", paste(missing_markers, collapse = ", "))
  markers <- setdiff(markers, missing_markers)
  dotD <- dot_summary(obj, markers, meta$celltype_major, dotplot_celltype_levels)
  marker_order <- unlist(marker_list, use.names = FALSE)
  marker_order <- marker_order[marker_order %in% markers]
  marker_key <- bind_rows(lapply(names(marker_list), function(ct) {
    data.frame(gene = marker_list[[ct]], marker_celltype = ct, stringsAsFactors = FALSE)
  })) %>% filter(gene %in% markers)
  dotD <- dotD %>% left_join(marker_key, by = "gene")
  dotD$gene <- factor(dotD$gene, levels = marker_order)
  dotD$group <- factor(dotD$group, levels = rev(dotplot_celltype_levels))
  dotD$marker_celltype <- factor(dotD$marker_celltype, levels = dotplot_celltype_levels)
  strip_labels <- c(`T cell` = "T\ncell", `NK cell` = "NK\ncell", `B cell` = "B\ncell", Monocyte = "Monocyte",
    Macrophage = "Macrophage", DC = "DC", `Mast cell` = "Mast\ncell", Epithelial = "Epithelial", Endothelial = "Endothelial",
    `Smooth muscle` = "Smooth\nmuscle")
  pD <- ggplot(dotD, aes(x = gene, y = group)) + geom_point(aes(size = percent_expressed, colour = scaled_average)) +
    facet_grid(. ~ marker_celltype, scales = "free_x", space = "free_x", labeller = labeller(marker_celltype = strip_labels)) +
    scale_size_continuous(range = c(0.4, 5.8), limits = c(0, 100), breaks = c(25,
      50, 75, 100), name = "Percent expressed") + scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7",
    high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Scaled average\nexpression") +
    labs(x = NULL, y = NULL) + base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
    size = 8, face = "italic"), axis.text.y = element_text(size = 7.8), axis.line = element_blank(),
    axis.ticks = element_blank(), panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.3),
    panel.grid.major.x = element_blank(), strip.background = element_rect(fill = "white", colour = "#4D4D4D",
      linewidth = 0.45), strip.text.x = element_text(size = 6.8, face = "bold",
      lineheight = 0.95, margin = margin(4, 2, 4, 2)), panel.spacing.x = unit(1.5, "mm"),
    legend.position = "right", legend.key.height = unit(3.5, "mm")) + guides(size = guide_legend(order = 1),
    colour = guide_colourbar(order = 2, barheight = unit(22, "mm"), barwidth = unit(3.2, "mm")))
  composition <- meta %>% count(figure_group, display_sample, celltype_major, name = "n_cells") %>% group_by(display_sample) %>%
    mutate(proportion = n_cells/sum(n_cells)) %>% ungroup()
  composition$figure_group <- factor(composition$figure_group, levels = group_levels)
  composition$display_sample <- factor(composition$display_sample, levels = rev(sample_levels))
  composition$celltype_major <- factor(composition$celltype_major, levels = rev(celltype_levels))
  pE <- ggplot(composition, aes(y = display_sample, x = proportion, fill = celltype_major)) + geom_col(width = 0.82,
    colour = NA) + facet_grid(figure_group ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_fill_manual(values = celltype_cols,
    breaks = celltype_levels, name = "Cell type", drop = FALSE) + scale_x_continuous(labels = scales::percent_format(accuracy = 1),
    breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.01))) + labs(x = "Cell proportion", y = NULL) +
    base_theme + theme(axis.text.x = element_text(size = 6.8), axis.text.y = element_text(size = 6.5),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), strip.placement = "outside", strip.background = element_rect(fill = "white",
      colour = "#4D4D4D", linewidth = 0.4), strip.text.y.left = element_text(angle = 0, face = "bold",
      size = 7.5, margin = margin(2, 3, 2, 3)), panel.spacing.y = unit(1.2, "mm"), legend.position = "bottom",
    legend.title = element_text(size = 7.2, face = "bold"), legend.text = element_text(size = 6.2),
    legend.key.size = unit(2.7, "mm"), plot.margin = margin(4, 4, 4, 4)) + guides(fill = guide_legend(reverse = FALSE,
    ncol = 4, byrow = TRUE, title.position = "top", title.hjust = 0.5, keywidth = unit(3, "mm"), keyheight = unit(2.7,
      "mm")))
  celltype_deg <- read.csv(celltype_deg_file, check.names = FALSE)
  lineage_deg <- if (file.exists(lineage_deg_file))
    read.csv(lineage_deg_file, check.names = FALSE)
  else data.frame(gene = character())
  immune_types <- c("T cell", "NK cell", "B cell", "Monocyte", "Macrophage", "DC")
  celltype_candidate_set <- celltype_deg %>% filter(celltype_major %in% immune_types) %>% distinct(celltype_major,
    gene, .keep_all = TRUE)
  if (any(celltype_candidate_set$PValue >= 0.05, na.rm = TRUE) || any(abs(celltype_candidate_set$logFC) <
    0.2, na.rm = TRUE) || any(celltype_candidate_set$same_direction_pairs < 4L, na.rm = TRUE) ||
    any(celltype_candidate_set$n_complete_pairs_checked != 5L, na.rm = TRUE)) {
    stop("The cell-type candidate table does not uniformly satisfy the stated ", "nominal P, |log2FC| and >=4/5 paired-direction criteria.")
  }
  celltype_candidate_counts <- celltype_candidate_set %>% mutate(celltype_major = factor(celltype_major, levels = immune_types)) %>%
    group_by(celltype_major) %>% summarise(n_candidate_genes = n_distinct(gene), n_high_up = n_distinct(gene[direction ==
    "HIGH_UP"]), n_low_up = n_distinct(gene[direction == "LOW_UP"]), .groups = "drop") %>% arrange(celltype_major)
  total_gene_celltype_associations <- nrow(celltype_candidate_set)
  total_unique_candidate_genes <- n_distinct(celltype_candidate_set$gene)
  celltype_count_labels <- setNames(paste0(as.character(celltype_candidate_counts$celltype_major), "\nCandidates: n=",
    scales::comma(celltype_candidate_counts$n_candidate_genes)), as.character(celltype_candidate_counts$celltype_major))
  representative_exclusion_pattern <- "^(IG[HKL]|TR[ABDG]|MT-|MT[0-9]|RP[SL]|AC[0-9]|AL[0-9]|LINC)"
  eligible_g <- celltype_deg %>% filter(celltype_major %in% immune_types, gene %in% rownames(obj), !grepl(representative_exclusion_pattern,
    gene, ignore.case = TRUE)) %>% mutate(celltype_major = factor(celltype_major, levels = immune_types)) %>%
    arrange(celltype_major, PValue, desc(abs(logFC)))
  top_genes_by_subtype <- eligible_g %>% group_by(celltype_major) %>% slice_head(n = 2) %>% ungroup() %>% arrange(celltype_major)
  if (nrow(top_genes_by_subtype) != 2L * length(immune_types)) {
    stop("Could not select exactly two representative genes for every immune subtype.")
  }
  panel_g_genes <- unique(top_genes_by_subtype$gene)
  if (length(panel_g_genes) != 2L * length(immune_types)) {
    stop("The top-two lists contain duplicated genes across immune subtypes; review is required.")
  }
  panel_g_groups <- c("High", "Low")
  meta_g <- meta %>% filter(celltype_major %in% immune_types, figure_group %in% panel_g_groups)
  meta_g$celltype_group <- paste(meta_g$celltype_major, meta_g$figure_group, sep = " | ")
  celltype_group_levels <- unlist(lapply(immune_types, function(ct) {
    paste(ct, panel_g_groups, sep = " | ")
  }), use.names = FALSE)
  obj_g <- subset(obj, cells = meta_g$cell_barcode)
  meta_g <- meta_g[match(colnames(obj_g), meta_g$cell_barcode), , drop = FALSE]
  stopifnot(identical(meta_g$cell_barcode, colnames(obj_g)))
  dotG <- dot_summary(obj_g, panel_g_genes, meta_g$celltype_group, celltype_group_levels)
  dotG$gene <- factor(dotG$gene, levels = panel_g_genes)
  dotG$group <- factor(dotG$group, levels = rev(celltype_group_levels))
  pG <- ggplot(dotG, aes(x = gene, y = group)) + geom_point(aes(size = percent_expressed, colour = scaled_average)) +
    scale_size_continuous(range = c(0.4, 5.8), limits = c(0, 100), breaks = c(25,
      50, 75, 100), name = "Percent expressed") + scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7",
    high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Scaled average\nexpression") +
    labs(x = NULL, y = NULL) + base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
    size = 8, face = "italic"), axis.text.y = element_text(size = 7.8), axis.line = element_blank(),
    axis.ticks = element_blank(), panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.3),
    panel.grid.major.x = element_blank(), legend.position = "right", legend.key.height = unit(3.5, "mm"), plot.margin = margin(7,
      8, 7, 7)) + guides(size = guide_legend(order = 1), colour = guide_colourbar(order = 2, barheight = unit(22,
    "mm"), barwidth = unit(3.2, "mm")))
  paired_patient_key <- donor_key %>% filter(donor_class == "Paired lesion patient") %>% select(patient_id, display_donor)
  paired_heatmap_list <- lapply(immune_types, function(ct) {
    selected_ct <- top_genes_by_subtype %>% filter(as.character(celltype_major) == ct)
    edge_file <- file.path(celltype_edger_rds_dir, paste0("edgeR_paired_", gsub(" ", "_", ct), ".rds"))
    if (!file.exists(edge_file))
      stop("Missing paired edgeR object: ", edge_file)
    edge_obj <- readRDS(edge_file)
    pair_diff <- edge_obj$patient_logCPM_differences
    if (is.null(pair_diff)) {
      stop("patient_logCPM_differences is absent from: ", edge_file)
    }
    missing_genes <- setdiff(selected_ct$gene, rownames(pair_diff))
    if (length(missing_genes)) {
      stop("Selected genes missing from paired differences for ", ct, ": ", paste(missing_genes, collapse = ", "))
    }
    pair_diff <- as.data.frame(pair_diff[selected_ct$gene, , drop = FALSE], check.names = FALSE)
    pair_diff$gene <- rownames(pair_diff)
    pair_long <- pair_diff %>% pivot_longer(cols = -gene, names_to = "patient_id", values_to = "high_minus_low_logCPM") %>%
      left_join(selected_ct %>% transmute(gene, celltype_major = as.character(celltype_major), model_logFC = logFC,
        nominal_P = PValue, same_direction_pairs, n_complete_pairs_checked, direction), by = "gene") %>%
      left_join(paired_patient_key, by = "patient_id")
    if (anyNA(pair_long$display_donor)) {
      stop("A paired patient in ", ct, " could not be mapped to Pt1-Pt5.")
    }
    direction_check <- pair_long %>% group_by(gene, direction, same_direction_pairs) %>% summarise(recalculated_same_direction_pairs = ifelse(first(direction) ==
      "HIGH_UP", sum(high_minus_low_logCPM > 0, na.rm = TRUE), sum(high_minus_low_logCPM < 0, na.rm = TRUE)),
      .groups = "drop")
    if (any(direction_check$recalculated_same_direction_pairs != direction_check$same_direction_pairs)) {
      stop("Direction-consistency audit failed for immune subtype: ", ct)
    }
    pair_long
  })
  paired_heatmap <- bind_rows(paired_heatmap_list) %>% mutate(display_donor = factor(display_donor, levels = paired_patient_key$display_donor),
    celltype_major = factor(celltype_major, levels = immune_types), direction_symbol = ifelse(direction ==
      "HIGH_UP", "\u2191", "\u2193"), consistency_label = paste0(direction_symbol, " ", same_direction_pairs, "/",
      n_complete_pairs_checked), gene_key = paste(celltype_major, gene, sep = "___"))
  heatmap_gene_order <- top_genes_by_subtype %>% transmute(gene_key = paste(as.character(celltype_major), gene,
    sep = "___")) %>% pull(gene_key)
  paired_heatmap$gene_key <- factor(paired_heatmap$gene_key, levels = rev(heatmap_gene_order))
  paired_heatmap_summary <- paired_heatmap %>% distinct(celltype_major, gene, gene_key, model_logFC, nominal_P,
    same_direction_pairs, n_complete_pairs_checked, direction, direction_symbol, consistency_label) %>% arrange(celltype_major,
    nominal_P)
  paired_heatmap$patient_x <- as.numeric(paired_heatmap$display_donor)
  paired_heatmap$cell_text_colour <- ifelse(abs(paired_heatmap$high_minus_low_logCPM) >= 3, "white", "black")
  consistency_annotation <- paired_heatmap_summary %>% mutate(annotation_x = 6.15)
  heatmap_limit <- max(abs(paired_heatmap$high_minus_low_logCPM), na.rm = TRUE)
  heatmap_limit <- max(1, ceiling(heatmap_limit))
  pF <- ggplot(paired_heatmap, aes(x = patient_x, y = gene_key, fill = high_minus_low_logCPM)) + geom_tile(colour = "white",
    linewidth = 0.55, width = 0.96, height = 0.92) + geom_text(aes(label = sprintf("%.2f",
    high_minus_low_logCPM), colour = cell_text_colour), size = 2.45, family = "Arial", show.legend = FALSE) +
    geom_text(data = consistency_annotation, aes(x = annotation_x, y = gene_key, label = consistency_label),
      inherit.aes = FALSE, size = 3, family = "Arial", fontface = "bold") + scale_colour_identity() + scale_fill_gradient2(low = "#2166AC",
    mid = "white", high = "#B2182B", midpoint = 0, limits = c(-heatmap_limit, heatmap_limit), oob = scales::squish,
    name = "High - Low\nlogCPM") + scale_x_continuous(breaks = c(seq_len(nrow(paired_patient_key)), 6.15),
    labels = c(paired_patient_key$display_donor, "Consistency"), limits = c(0.5, 6.75), expand = c(0, 0), position = "top") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x), drop = TRUE) + facet_grid(celltype_major ~
    ., scales = "free_y", space = "free_y", switch = "y", labeller = labeller(celltype_major = as_labeller(celltype_count_labels))) +
    labs(title = "Patient-paired High-Low expression differences", subtitle = paste0(scales::comma(total_unique_candidate_genes),
      " unique candidates; P < 0.05, |log2FC| >= 0.20, ", "concordant in >=4/5 pairs"), x = NULL, y = NULL) +
    base_theme + theme(plot.title = element_text(size = 10.5, face = "bold", hjust = 0), plot.subtitle = element_text(size = 7.6,
    hjust = 0, colour = "black", margin = margin(b = 5)), axis.text.x = element_text(size = 8.2,
    face = "bold", angle = 0, hjust = 0.5), axis.text.y = element_text(size = 8.2, face = "italic",
    colour = "black"), axis.line = element_blank(), axis.ticks = element_blank(), strip.placement = "outside",
    strip.background = element_rect(fill = "#F2F2F2", colour = "#707070", linewidth = 0.4),
    strip.text.y.left = element_text(size = 7.8, face = "bold", angle = 0, margin = margin(3,
      4, 3, 4)), panel.spacing.y = unit(1.1, "mm"), panel.grid = element_blank(), legend.position = "right",
    legend.title = element_text(size = 7.5), legend.text = element_text(size = 7), plot.margin = margin(8,
      10, 7, 8))
  pA <- remove_plot_titles(pA)
  pB <- remove_plot_titles(pB)
  pC <- remove_plot_titles(pC)
  pD <- remove_plot_titles(pD)
  pE <- remove_plot_titles(pE)
  pF <- remove_plot_titles(pF)
  pG <- remove_plot_titles(pG)
  write.csv(plot_df %>% select(cell_barcode, UMAP_1, UMAP_2, figure_group, sample_id, display_sample, patient_id,
    display_donor, celltype_major, lineage_major), file.path(outdir, "source_data", "panels_A_B_C_UMAP_coordinates.csv"),
    row.names = FALSE)
  write.csv(dotD, file.path(outdir, "source_data", "panel_D_marker_dotplot.csv"), row.names = FALSE)
  write.csv(composition, file.path(outdir, "source_data", "panel_E_celltype_proportions.csv"), row.names = FALSE)
  write.csv(dotG, file.path(outdir, "source_data", "panel_G_expression_dotplot.csv"), row.names = FALSE)
  write.csv(top_genes_by_subtype, file.path(outdir, "source_data", "panel_G_top2_genes_by_immune_subtype.csv"),
    row.names = FALSE)
  write.csv(eligible_g, file.path(outdir, "source_data", "panel_G_eligible_candidates_ranked.csv"), row.names = FALSE)
  write.csv(lineage_deg[lineage_deg$gene %in% panel_g_genes, , drop = FALSE], file.path(outdir, "source_data",
    "panel_G_major_lineage_DE_context.csv"), row.names = FALSE)
  write.csv(paired_heatmap %>% transmute(celltype_major, gene, patient_id, displayed_patient = as.character(display_donor),
    high_minus_low_logCPM, model_logFC, nominal_P, same_direction_pairs, n_complete_pairs_checked, direction,
    consistency_label), file.path(outdir, "source_data", "panel_F_patient_paired_logCPM_differences.csv"),
    row.names = FALSE)
  write.csv(paired_heatmap_summary %>% select(-gene_key), file.path(outdir, "source_data", "panel_F_direction_consistency_summary.csv"),
    row.names = FALSE)
  write.csv(bind_rows(celltype_candidate_counts %>% transmute(scope = as.character(celltype_major), n_unique_genes = n_candidate_genes,
    n_gene_celltype_associations = n_candidate_genes, n_high_up, n_low_up), data.frame(scope = "All six immune cell subtypes (union)",
    n_unique_genes = total_unique_candidate_genes, n_gene_celltype_associations = total_gene_celltype_associations,
    n_high_up = NA_integer_, n_low_up = NA_integer_, stringsAsFactors = FALSE)), file.path(outdir, "source_data",
    "panel_F_candidate_gene_counts.csv"), row.names = FALSE)
  if (only_e_update) {
    clean_previous_panel_outputs("E")
    manifest <- export_panel(pE, "E_celltype_proportions", panel_e_width_mm, panel_e_height_mm)
    write.csv(manifest, file.path(outdir, "Figure2_E_update_manifest.csv"), row.names = FALSE)
  }
  else if (only_f_update) {
    clean_previous_panel_outputs("F")
    manifest <- export_panel(pF, "F_patient_paired_HIGH_minus_LOW_difference_heatmap", 170, 118)
    write.csv(manifest, file.path(outdir, "Figure2_F_update_manifest.csv"), row.names = FALSE)
  }
  else if (only_abcd_update) {
    clean_previous_panel_outputs("ABCD")
    manifest <- bind_rows(export_panel(pA, "A_UMAP_sample_origin", 66, 62), export_panel(pB, "B_UMAP_donor",
      82, 62), export_panel(pC, "C_UMAP_celltype_legend", 86, 62), export_panel(pD, "D_canonical_marker_dotplot_annotated",
      panel_d_width_mm, panel_d_height_mm))
    write.csv(manifest, file.path(outdir, "Figure2_ABCD_update_manifest.csv"), row.names = FALSE)
  }
  else if (only_bc_update) {
    clean_previous_panel_outputs("BC")
    manifest <- bind_rows(export_panel(pB, "B_UMAP_donor_paired_pt1p2", 82, 62), export_panel(pC, "C_UMAP_celltype_brewer_pt1p2",
      76, 62))
    write.csv(manifest, file.path(outdir, "Figure2_BC_update_manifest.csv"), row.names = FALSE)
  }
  else {
    clean_previous_panel_outputs("ALL")
    manifest <- bind_rows(export_panel(pA, "A_UMAP_sample_origin", 66, 62), export_panel(pB, "B_UMAP_donor",
      82, 62), export_panel(pC, "C_UMAP_celltype_legend", 86, 62), export_panel(pD, "D_canonical_marker_dotplot_annotated",
      panel_d_width_mm, panel_d_height_mm), export_panel(pE, "E_celltype_proportions", panel_e_width_mm,
      panel_e_height_mm), export_panel(pF, "F_patient_paired_HIGH_minus_LOW_difference_heatmap", 170, 118),
      export_panel(pG, "G_top2_genes_per_immune_subtype_dotplot", 183, 102))
    write.csv(manifest, file.path(outdir, "Figure2_panel_manifest.csv"), row.names = FALSE)
  }
  audit <- data.frame(item = c("Seurat RDS", "Cells", "Genes", "Groups", "Samples", "Cell types", "Lineages",
    "Formal DE contrast", "Paired patients", "Panel G rule", "Panel G genes", "Panel F evidence", "HD interpretation"),
    value = c(input_rds, ncol(obj), nrow(obj), paste(levels(meta$figure_group), collapse = "; "), length(unique(meta$sample_id)),
      paste(celltype_levels, collapse = "; "), "LYMPHOID; MYELOID; OTHER", "High minus Low, paired pseudobulk",
      paste(sort(unique(top_genes_by_subtype$patient_ids)), collapse = "; "), paste0("Within each immune subtype: widest candidates only; exclude ",
        representative_exclusion_pattern, "; rank by nominal P ascending then |log2FC| descending; retain top 2"),
      paste(panel_g_genes, collapse = "; "), paste0("Exact patient-level High-minus-Low TMM-logCPM differences from the ",
        "saved cell-type edgeR objects; direction consistency independently ", "recalculated and checked against the candidate table"),
      "Descriptive reference only; not used in the formal High-vs-Low DE test"), stringsAsFactors = FALSE)
  write.csv(audit, file.path(outdir, "Figure2_analysis_audit.csv"), row.names = FALSE)
  legend_text <- c("Figure 2 | PET-CT metabolic heterogeneity is associated with immune transcriptional reprogramming in pulmonary tuberculosis lesions.",
    "(A) UMAP of the Harmony-integrated single-cell atlas coloured by sample origin (HD, High and Low).", "(B) The same UMAP coloured by donor. Paired High and Low lesions from the same patient share a colour, whereas each HD donor has a distinct colour; donor display labels are mapped in the accompanying source-data file.",
    "(C) UMAP coloured by the audited major cell-type annotation; cell types are identified in the colour legend.",
    "(D) Dot plot of canonical markers supporting major immune and structural cell identities: CD3D and TRAC for T cells; NKG7 and GNLY for NK cells; MS4A1 and CD79A for B cells; S100A8 and FCN1 for monocytes; C1QB and MARCO for macrophages; IDO1 and DNASE1L3 for dendritic cells; TPSAB1 and KIT for mast cells; EPCAM and KRT19 for epithelial cells; PECAM1 and CLDN5 for endothelial cells; and ACTA2 and TAGLN for smooth-muscle cells. Dot size indicates the percentage of expressing cells and colour indicates scaled average expression.",
    "(E) Horizontal 100% stacked bars showing cell-type proportions for each sample, grouped by HD, High and Low.",
    paste0("(F) Heatmap of patient-level paired expression differences for the representative genes subsequently shown in G. Each cell is the TMM-normalized logCPM difference calculated as High minus Low for one patient; red and blue therefore indicate higher expression in High and Low lesions, respectively. Fractions on the right report the number of five patient pairs sharing the majority direction. Facet labels report the complete number of candidate genes in each immune cell subtype under the prespecified nominal P < 0.05, |log2FC| >= 0.20 and >=4/5 direction-concordance criteria. Across the six immune subtypes, there were ",
      scales::comma(total_gene_celltype_associations), " gene-cell-type associations corresponding to ",
      scales::comma(total_unique_candidate_genes), " unique candidate genes. Arrows indicate the majority direction and do not imply that each individual pair was independently statistically significant."),
    "(G) Expression of the same representative paired High-Low candidates shown in F across the High and Low cell-type groups. Within each subtype, eligible genes were ranked by nominal P value, with absolute log2 fold change used to break ties. Immunoglobulin/TCR, mitochondrial, ribosomal and unannotated non-coding features were excluded from representative-gene selection. Dot size indicates the percentage of expressing cells and colour indicates scaled average expression.",
    "Formal differential-expression testing used the five paired High and Low lesions only. HD cells are displayed as a descriptive reference and were not included in the paired High-Low test.")
  writeLines(legend_text, file.path(outdir, "Figure2_legend_current_fixed.txt"))
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  message("Completed. Output directory: ", outdir)
  invisible(NULL)
}

acti_stage_dotplot <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })
  b <- cfg$output_dir
  out <- file.path(b, "07_ACTI_expression")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  genes <- c("DUSP3", "FCGR1A", "NPC2", "TAP1", "ANKRD22", "GBP5", "VAMP5")
  types <- c("T cell", "NK cell", "B cell", "Monocyte", "Macrophage", "DC")
  obj <- readRDS(acti_input_object(cfg))
  m <- obj@meta.data
  stopifnot(identical(rownames(m), colnames(obj)), all(genes %in% rownames(obj)))
  writeLines(names(m), file.path(out, "metadata_columns.txt"))
  m$Group <- recode(toupper(as.character(m$tissue_type)), CONTROL = "HD", HIGH = "HIGH", LOW = "LOW")
  keep <- m$celltype_major %in% types & m$Group %in% c("HD", "LOW", "HIGH") & !is.na(m$use_for_DE) & as.logical(m$use_for_DE)
  m <- m[keep, , drop = FALSE]
  x <- GetAssayData(obj, assay = "RNA", slot = "data")[genes, rownames(m), drop = FALSE]
  stopifnot(ncol(x) == nrow(m), all(is.finite(x@x)))
  message("Donors by group:")
  print(m %>% distinct(patient_id, Group) %>% count(Group))
  id <- interaction(m$celltype_major, m$patient_id, m$Group, drop = TRUE)
  summ <- lapply(split(seq_len(nrow(m)), id), function(idx) {
    v <- as.matrix(x[, idx, drop = FALSE])
    md <- m[idx[1], ]
    data.frame(Celltype = as.character(md$celltype_major), Donor = as.character(md$patient_id), Group = md$Group,
      Gene = genes, n_cells = length(idx), avg_normalized_expression = rowMeans(expm1(v)), percent_expressed = rowMeans(v >
        0) * 100)
  }) %>% bind_rows()
  audit <- summ %>% distinct(Celltype, Donor, Group, n_cells) %>% mutate(enough_cells = n_cells >= 20)
  pairs <- audit %>% filter(Group %in% c("HIGH", "LOW"), enough_cells) %>% group_by(Celltype, Donor) %>% summarise(n_conditions = n_distinct(Group),
    .groups = "drop") %>% filter(n_conditions == 2)
  audit <- audit %>% left_join(pairs %>% mutate(complete_pair = TRUE) %>% select(-n_conditions), by = c("Celltype",
    "Donor")) %>% mutate(complete_pair = replace_na(complete_pair, FALSE), Included = enough_cells & (Group ==
    "HD" | complete_pair))
  summ <- summ %>% left_join(audit %>% select(Celltype, Donor, Group, Included), by = c("Celltype", "Donor",
    "Group"))
  use <- summ %>% filter(Included)
  stopifnot(all(types %in% unique(use$Celltype)))
  dot <- use %>% group_by(Celltype, Group, Gene) %>% summarise(avg_normalized_expression = mean(avg_normalized_expression),
    percent_expressed = mean(percent_expressed), n_donors = n_distinct(Donor), .groups = "drop") %>% mutate(log1p_average = log1p(avg_normalized_expression)) %>%
    group_by(Gene) %>% mutate(scaled_average = if (sd(log1p_average) > 0)
    as.numeric(scale(log1p_average))
  else 0) %>% ungroup()
  stopifnot(nrow(dot) == length(types) * 3 * length(genes))
  dot$Gene <- factor(dot$Gene, levels = genes)
  group_display <- c(HD = "HD", LOW = "Low", HIGH = "High")
  dot$row <- paste(dot$Celltype, unname(group_display[dot$Group]), sep = " | ")
  dot$Compartment <- ifelse(dot$Celltype %in% types[1:3], "Lymphoid", "Myeloid")
  make_plot <- function(dd, ct) {
    rows <- unlist(lapply(ct, function(z) paste(z, c("HD", "Low", "High"), sep = " | ")))
    dd$row <- factor(dd$row, levels = rev(rows))
    ggplot(dd, aes(Gene, row)) + geom_point(aes(size = percent_expressed, colour = scaled_average)) + scale_size_area(max_size = 6,
      limits = c(0, 100), breaks = c(25, 50, 75, 100), name = "Percent expressed") + scale_colour_gradient2(low = "#2166AC",
      mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Scaled average\nexpression") +
      theme_bw(base_size = 12) + labs(x = NULL, y = NULL, title = NULL) + theme(panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.35),
      panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1.2), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic", size = 12, colour = "black"),
      axis.text.y = element_text(size = 11, colour = "black"), legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10))
  }
  plots <- list(All_immune = list(data = dot, ct = types, height = 8), Lymphoid = list(data = filter(dot, Compartment ==
    "Lymphoid"), ct = types[1:3], height = 4.8), Myeloid = list(data = filter(dot, Compartment ==
    "Myeloid"), ct = types[4:6], height = 4.8))
  for (nm in names(plots)) {
    a <- plots[[nm]]
    p <- make_plot(a$data, a$ct)
    prefix <- file.path(out, paste0("ACTI7_", nm, "_HD_Low_High_dotplot_labelcase"))
    ggsave(paste0(prefix, ".pdf"), p, width = 8, height = a$height, device = cairo_pdf)
    ggsave(paste0(prefix, ".png"), p, width = 8, height = a$height, dpi = 600, bg = "white")
    ggsave(paste0(prefix, ".svg"), p, width = 8, height = a$height, device = svglite::svglite)
  }
  delta <- use %>% filter(Group %in% c("HIGH", "LOW")) %>% mutate(log1p_expression = log1p(avg_normalized_expression)) %>%
    select(Celltype, Donor, Group, Gene, log1p_expression) %>% pivot_wider(names_from = Group, values_from = log1p_expression) %>%
    mutate(HIGH_minus_LOW = HIGH - LOW)
  write.csv(dot, file.path(out, "Dotplot_source_data.csv"), row.names = FALSE)
  write.csv(summ, file.path(out, "Donor_level_expression_and_detection.csv"), row.names = FALSE)
  write.csv(audit, file.path(out, "Donor_pair_and_cell_count_audit.csv"), row.names = FALSE)
  write.csv(delta, file.path(out, "Patient_paired_expression_differences.csv"), row.names = FALSE)
  print(audit %>% filter(Included) %>% count(Celltype, Group, name = "n_donors"))
  cat("OUTPUT:", out, "\n")
  invisible(NULL)
}

acti_stage_paired_umap <- function(cfg) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(svglite)
  })
  b <- cfg$output_dir
  out <- file.path(b, "08_paired_expression")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(acti_input_object(cfg))
  md <- obj@meta.data
  u <- Embeddings(obj, reduction = "umap")
  stopifnot(identical(rownames(md), rownames(u)))
  genes <- c("NPC2", "VAMP5")
  rna <- GetAssayData(obj, assay = "RNA", slot = "data")[genes, , drop = FALSE]
  keep <- toupper(as.character(md$tissue_type)) %in% c("HIGH", "LOW")
  d <- data.frame(Cell = rownames(md)[keep], patient_id = as.character(md$patient_id[keep]), sample_id = as.character(md$sample_id[keep]),
    Group = toupper(as.character(md$tissue_type[keep])), Celltype = as.character(md$celltype_major[keep]),
    UMAP_1 = u[keep, 1], UMAP_2 = u[keep, 2], NPC2 = as.numeric(rna["NPC2", keep]), VAMP5 = as.numeric(rna["VAMP5",
      keep]))
  pt <- sort(unique(d$patient_id))
  stopifnot(length(pt) == 5)
  stopifnot(all(table(d$patient_id, d$Group) > 0))
  d$Patient <- factor(paste0("Pt", match(d$patient_id, pt)), levels = paste0("Pt", 1:5))
  d$Group <- factor(d$Group, levels = c("HIGH", "LOW"))
  limits <- c(0, max(as.matrix(d[, genes]), na.rm = TRUE))
  xr <- range(u[, 1])
  yr <- range(u[, 2])
  write.csv(data.frame(Patient = paste0("Pt", 1:5), patient_id = pt), file.path(out, "Patient_display_mapping.csv"),
    row.names = FALSE)
  write.csv(as.data.frame(table(d$Patient, d$Group)), file.path(out, "Cell_counts_by_patient_group.csv"), row.names = FALSE)
  write.csv(d, file.path(out, "NPC2_VAMP5_fixed_UMAP_source_data.csv"), row.names = FALSE)
  de <- read.csv(file.path(b, "04_celltype_pseudobulk/tables/06_celltype_WIDEST_candidates.csv"), check.names = FALSE)
  de <- de[de$gene %in% genes, , drop = FALSE]
  stopifnot(nrow(de) == 3)
  write.csv(de, file.path(out, "Outlined_paired_candidate_celltypes.csv"), row.names = FALSE)
  ref_keep <- keep & !is.na(md$use_for_DE) & as.logical(md$use_for_DE)
  paths <- list()
  contour_audit <- list()
  min_component_cell_fraction <- 0.05
  for (ct in unique(de$celltype_major)) {
    idx <- which(ref_keep & as.character(md$celltype_major) == ct)
    z <- MASS::kde2d(u[idx, 1], u[idx, 2], n = 160, lims = c(xr, yr))
    zs <- sort(as.vector(z$z), decreasing = TRUE)
    threshold <- zs[which(cumsum(zs)/sum(zs) >= 0.9)[1]]
    cs <- contourLines(z$x, z$y, z$z, levels = threshold)
    component_n <- vapply(cs, function(cc) sum(sp::point.in.polygon(u[idx, 1], u[idx, 2], cc$x, cc$y) > 0),
      integer(1))
    component_fraction <- component_n/length(idx)
    retain <- component_fraction >= min_component_cell_fraction
    contour_audit[[ct]] <- data.frame(Celltype = ct, Component = seq_along(cs), Total_reference_cells = length(idx),
      Cells_inside_component = component_n, Cell_fraction = component_fraction, Minimum_cell_fraction = min_component_cell_fraction,
      Retained = retain)
    stopifnot(any(retain))
    paths[[ct]] <- do.call(rbind, lapply(which(retain), function(i) data.frame(UMAP_1 = cs[[i]]$x, UMAP_2 = cs[[i]]$y,
      Celltype = ct, Contour = paste(ct, i, sep = "_"))))
  }
  path_data <- do.call(rbind, paths)
  write.csv(do.call(rbind, contour_audit), file.path(out, "Contour_component_filter_audit.csv"), row.names = FALSE)
  write.csv(path_data, file.path(out, "Reference_contour_coordinates.csv"), row.names = FALSE)
  line_cols <- c(Macrophage = "#164B63", `T cell` = "#6B3F80")
  line_types <- c(Macrophage = "solid", `T cell` = "dashed")
  for (g in genes) {
    dd <- d[order(d[[g]]), ]
    dd$Expression <- dd[[g]]
    p <- ggplot(dd, aes(UMAP_1, UMAP_2, colour = Expression)) + geom_point(size = 0.25, alpha = 1, stroke = 0) +
      facet_grid(Patient ~ Group, switch = "y") + scale_colour_gradientn(colours = c("#E8E8E8", "#F8D9D5",
      "#E78B82", "#B2182B"), limits = limits, name = "Normalized\nexpression") + coord_fixed(xlim = xr, ylim = yr) +
      labs(x = NULL, y = NULL, title = NULL) + theme_bw(base_size = 9) + theme(panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey30", fill = NA, size = 1), axis.text = element_blank(), axis.ticks = element_blank(),
      strip.background = element_blank(), strip.text.x = element_text(size = 11, face = "bold"), strip.text.y.left = element_text(size = 10,
        face = "bold", angle = 0), strip.placement = "outside", panel.spacing = grid::unit(2, "mm"), legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8), plot.margin = margin(3,
        3, 3, 3, "mm")) + guides(colour = guide_colourbar(direction = "horizontal", barwidth = grid::unit(35,
      "mm"), barheight = grid::unit(3, "mm"), title.position = "left"))
    chosen <- unique(de$celltype_major[de$gene == g])
    for (ct in chosen) {
      cc <- path_data[path_data$Celltype == ct, ]
      p <- p + geom_path(data = cc, aes(UMAP_1, UMAP_2, group = Contour, linetype = Celltype), inherit.aes = FALSE,
        colour = line_cols[ct], size = 0.45, show.legend = TRUE)
    }
    p <- p + scale_linetype_manual(values = line_types, breaks = chosen, name = "Outlined cell types") + guides(linetype = guide_legend(nrow = 1,
      override.aes = list(colour = unname(line_cols[chosen]), size = 0.6)), colour = guide_colourbar(direction = "horizontal",
      barwidth = grid::unit(30, "mm"), barheight = grid::unit(3, "mm"), title.position = "left")) + theme(legend.box = "vertical",
      legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), legend.spacing.y = grid::unit(0,
        "mm"))
    prefix <- file.path(out, paste0(g, "_paired_HIGH_LOW_UMAP_vertical_outlined"))
    ggsave(paste0(prefix, ".pdf"), p, width = 110, height = 235, units = "mm", device = cairo_pdf)
    ggsave(paste0(prefix, ".svg"), p, width = 110, height = 235, units = "mm", device = svglite)
    ggsave(paste0(prefix, ".png"), p, width = 110, height = 235, units = "mm", dpi = 600, bg = "white")
  }
  cat("OUTPUT:", out, "\n")
  print(table(d$Patient, d$Group))
  print(limits)
  invisible(NULL)
}

acti_stage_merged_umap <- function(cfg) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(svglite)
  })
  out <- file.path(cfg$output_dir, "08_paired_expression")
  d <- read.csv(file.path(out, "NPC2_VAMP5_fixed_UMAP_source_data.csv"), check.names = FALSE)
  paths <- read.csv(file.path(out, "Reference_contour_coordinates.csv"), check.names = FALSE)
  de <- read.csv(file.path(out, "Outlined_paired_candidate_celltypes.csv"), check.names = FALSE)
  genes <- c("NPC2", "VAMP5")
  ld <- do.call(rbind, lapply(genes, function(g) {
    z <- d
    z$Gene <- g
    z$Expression <- z[[g]]
    z
  }))
  ld$Gene <- factor(ld$Gene, levels = genes)
  ld$Group <- factor(ld$Group, levels = c("HIGH", "LOW"))
  ld$Patient <- factor(ld$Patient, levels = paste0("Pt", 1:5))
  ld <- ld[order(ld$Expression), ]
  cs <- do.call(rbind, lapply(genes, function(g) {
    z <- paths[paths$Celltype %in% de$celltype_major[de$gene == g], ]
    z$Gene <- factor(g, levels = genes)
    z
  }))
  cols <- c(Macrophage = "#164B63", `T cell` = "#6B3F80")
  types <- c(Macrophage = "solid", `T cell` = "dashed")
  limits <- c(0, max(ld$Expression))
  p <- ggplot(ld, aes(UMAP_1, UMAP_2, colour = Expression)) + geom_point(size = 0.3, alpha = 1,
    stroke = 0) + facet_grid(Patient ~ Gene + Group, switch = "y", labeller = labeller(Group = c(HIGH = "High",
    LOW = "Low"))) + scale_colour_gradientn(colours = c("#E8E8E8", "#F8D9D5", "#E78B82", "#B2182B"), limits = limits,
    name = "Normalized expression") + coord_fixed(ratio = 1/1.5, xlim = range(d$UMAP_1), ylim = range(d$UMAP_2)) +
    labs(title = NULL, x = NULL, y = NULL) + theme_bw(base_size = 11) + theme(panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey30", fill = NA, size = 1), axis.text = element_blank(), axis.ticks = element_blank(),
    strip.background = element_blank(), strip.text.x = element_text(size = 13, face = "bold"), strip.text.y.left = element_text(size = 12,
      face = "bold", angle = 0), strip.placement = "outside", panel.spacing.x = grid::unit(2, "mm"), panel.spacing.y = grid::unit(2,
      "mm"), legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 10,
      face = "bold"), legend.text = element_text(size = 10), plot.margin = margin(3, 3, 3, 3, "mm"))
  for (ct in names(cols)) {
    cc <- cs[cs$Celltype == ct, ]
    p <- p + geom_path(data = cc, aes(UMAP_1, UMAP_2, group = Contour, linetype = Celltype), inherit.aes = FALSE,
      colour = cols[ct], size = 0.45)
  }
  p <- p + scale_linetype_manual(values = types, breaks = names(cols), name = "Outlined cell types") + guides(linetype = guide_legend(order = 1,
    nrow = 1, override.aes = list(colour = unname(cols), size = 0.6)), colour = guide_colourbar(order = 2,
    direction = "horizontal", barwidth = grid::unit(35, "mm"), barheight = grid::unit(3, "mm"), title.position = "top"))
  prefix <- file.path(out, "NPC2_VAMP5_paired_High_Low_UMAP_merged_5x4_wide1p5_labelcase")
  ggsave(paste0(prefix, ".pdf"), p, width = 360, height = 260, units = "mm", device = cairo_pdf)
  ggsave(paste0(prefix, ".svg"), p, width = 360, height = 260, units = "mm", device = svglite)
  ggsave(paste0(prefix, ".png"), p, width = 360, height = 260, units = "mm", dpi = 600, bg = "white")
  write.csv(ld, file.path(out, "Merged_UMAP_source_data.csv"), row.names = FALSE)
  cat(prefix, "\n")
  invisible(NULL)
}

acti_stage_overlap <- function(cfg) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(svglite)
    library(ragg)
  })
  sc_file <- file.path(cfg$output_dir, "04_celltype_pseudobulk", "tables", "06_celltype_WIDEST_candidates.csv")
  stpmg_file <- cfg$stpmg_csv
  outdir <- file.path(cfg$output_dir, "09_integrated_overlap")
  source_dir <- file.path(outdir, "source_data")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(sc_file), file.exists(stpmg_file))
  sc_raw <- read.csv(sc_file, check.names = FALSE)
  stpmg_raw <- read.csv(stpmg_file, check.names = FALSE)
  stopifnot("gene" %in% names(sc_raw), "gene" %in% names(stpmg_raw))
  clean_gene <- function(x) {
    x <- toupper(trimws(as.character(x)))
    sort(unique(x[!is.na(x) & nzchar(x)]))
  }
  sc_genes <- clean_gene(sc_raw$gene)
  stpmg_genes <- clean_gene(stpmg_raw$gene)
  shared_genes <- intersect(sc_genes, stpmg_genes)
  sc_only <- setdiff(sc_genes, stpmg_genes)
  stpmg_only <- setdiff(stpmg_genes, sc_genes)
  stopifnot(length(sc_genes) == 4155L)
  stopifnot(length(stpmg_genes) == 1174L)
  stopifnot(length(shared_genes) == 338L)
  all_genes <- sort(union(sc_genes, stpmg_genes))
  membership <- data.frame(gene = all_genes, in_scRNA_candidate_set = all_genes %in% sc_genes, in_STPMGs = all_genes %in%
    stpmg_genes, category = ifelse(all_genes %in% shared_genes, "Integrated candidate genes (ICGs)", ifelse(all_genes %in%
    sc_only, "scRNA-seq candidates only", "STPMGs only")), stringsAsFactors = FALSE)
  summary_tab <- data.frame(set = c("scRNA-seq candidates (total)", "STPMGs (total)", "Integrated candidate genes (intersection)",
    "scRNA-seq candidates only", "STPMGs only"), n_genes = c(length(sc_genes), length(stpmg_genes), length(shared_genes),
    length(sc_only), length(stpmg_only)), stringsAsFactors = FALSE)
  write.csv(data.frame(gene = shared_genes, set_name = "ICGs", full_name = "Integrated candidate genes", stringsAsFactors = FALSE),
    file.path(source_dir, "A_integrated_candidate_genes_338.csv"), row.names = FALSE)
  write.csv(membership, file.path(source_dir, "A_gene_set_overlap_membership.csv"), row.names = FALSE)
  write.csv(summary_tab, file.path(source_dir, "A_gene_set_overlap_summary.csv"), row.names = FALSE)
  sc_shared_evidence <- sc_raw[toupper(trimws(sc_raw$gene)) %in% shared_genes, , drop = FALSE]
  write.csv(sc_shared_evidence, file.path(source_dir, "A_integrated_candidate_genes_celltype_evidence.csv"),
    row.names = FALSE)
  circle_xy <- function(cx, cy, radius, group, n = 720L) {
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(x = cx + radius * cos(theta), y = cy + radius * sin(theta), group = group)
  }
  circles <- rbind(circle_xy(-0.78, 0, 1.22, "scRNA-seq candidates"), circle_xy(0.78,
    0, 1.22, "STPMGs"))
  fill_cols <- c(`scRNA-seq candidates` = "#8FB8DE", STPMGs = "#E8A09A")
  p <- ggplot(circles, aes(x = x, y = y, group = group, fill = group, colour = group)) + geom_polygon(alpha = 0.58,
    linewidth = 1.25) + scale_fill_manual(values = fill_cols) + scale_colour_manual(values = c(`scRNA-seq candidates` = "#3F78A8",
    STPMGs = "#B75550")) + annotate("text", x = -1.18, y = 0.42, label = "scRNA-seq\ncandidates",
    family = "Arial", fontface = "bold", size = 5.2, lineheight = 0.92) + annotate("text",
    x = -1.18, y = -0.28, label = "Total n = 4,155\ngenes", family = "Arial",
    fontface = "bold", size = 4.1, lineheight = 0.95) + annotate("text", x = 1.18,
    y = 0.42, label = "STPMGs", family = "Arial", fontface = "bold", size = 5.2) +
    annotate("text", x = 1.18, y = -0.28, label = "Total n = 1,174\ngenes", family = "Arial",
      fontface = "bold", size = 4.1, lineheight = 0.95) + annotate("text",
    x = 0, y = 0.27, label = "ICGs", family = "Arial", fontface = "bold", size = 5.2) +
    annotate("text", x = 0, y = -0.28, label = "n = 338\ngenes", family = "Arial", fontface = "bold",
      size = 4.3, lineheight = 0.95) + coord_fixed(xlim = c(-2.18,
    2.18), ylim = c(-1.35, 1.35), clip = "off") + theme_void(base_family = "Arial") +
    theme(legend.position = "none", plot.margin = margin(5, 5, 5, 5, unit = "mm"))
  stem <- file.path(outdir, "A_scRNA_candidates_STPMGs_venn")
  svglite::svglite(paste0(stem, ".svg"), width = 5.4, height = 3.45, bg = "white")
  print(p)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 5.4, height = 3.45, bg = "white",
    family = "Arial")
  print(p)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = 137, height = 88, units = "mm", res = 600, background = "white",
    compression = "lzw")
  print(p)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = 137, height = 88, units = "mm", res = 300, background = "white")
  print(p)
  grDevices::dev.off()
  capture.output(sessionInfo(), file = file.path(outdir, "A_sessionInfo.txt"))
  message("Completed Figure 4A: ", paste0(stem, ".pdf"))
  message("Counts: scRNA-seq=", length(sc_genes), "; STPMGs=", length(stpmg_genes), "; overlap=", length(shared_genes))
  invisible(NULL)
}

acti_stage_validate <- function(cfg) {
  f <- file.path(cfg$output_dir, "04_celltype_pseudobulk", "tables", "06_celltype_WIDEST_candidates.csv")
  if (!file.exists(f))
    stop("Run celltype first; candidate table not found: ", f)
  x <- read.csv(f, check.names = FALSE)
  stopifnot(nrow(x) > 0L, all(acti_candidate_filter(x)), all(x$n_complete_pairs_checked == 5L))
  expected <- c(`T cell` = 1049L, `NK cell` = 1317L, `B cell` = 385L, Monocyte = 713L, Macrophage = 1805L, DC = 372L)
  actual <- table(x$celltype_major)
  counts <- as.integer(actual[names(expected)])
  report <- data.frame(celltype = names(expected), observed = counts, manuscript = unname(expected), matches = counts ==
    unname(expected))
  union_n <- length(unique(x$gene))
  print(report, row.names = FALSE)
  cat("Gene-cell-type associations:", nrow(x), "; unique genes:", union_n, "\n")
  write.csv(report, file.path(cfg$output_dir, "manuscript_count_validation.csv"), row.names = FALSE)
  if (anyNA(report$matches) || !all(report$matches) || nrow(x) != 5641L || union_n != 4155L)
    stop("Results differ from the manuscript reference counts; inspect data, environment and annotation. No counts were forced.")
  qf <- file.path(cfg$output_dir, "01_preannotation", "tables", "QC_scDblFinder_summary_total.csv")
  if (file.exists(qf)) {
    q <- read.csv(qf)
    stopifnot(q$n_samples == 16L, q$cells_before_qc == 109472L, q$cells_after_primary_qc == 104322L, q$predicted_doublets ==
      10174L, q$retained_singlets == 94148L)
  }
  invisible(report)
}

run_acti_scrna <- function(stages, cfg = acti_config()) {
  if (missing(stages) || !length(stages))
    stop("Choose an explicit analysis stage.")
  if (identical(stages, "all")) {
    stages <- setdiff(names(acti_stage_packages), "overlap")
    message("all stops at manual annotation if no reviewed map/explicit manuscript-map flag is supplied.")
  }
  if (any(!stages %in% names(acti_stage_packages)))
    stop("Unknown stages: ", paste(setdiff(stages, names(acti_stage_packages)), collapse = ", "))
  normalize <- function(x) tolower(normalizePath(x, winslash = "/", mustWork = FALSE))
  data_root <- sub("/+$", "", normalize(cfg$data_dir))
  out_root <- sub("/+$", "", normalize(cfg$output_dir))
  if (identical(data_root, out_root) || startsWith(paste0(out_root, "/"), paste0(data_root, "/")))
    stop("output_dir must be outside data_dir.")
  for (stage in stages) {
    acti_require_environment(stage)
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
    message("Running stage: ", stage)
    f <- get(paste0("acti_stage_", stage), mode = "function")
    f(cfg)
    gc(verbose = FALSE)
  }
  invisible(cfg)
}

acti_help <- function() cat("Rscript 1_ACTI_single_cell_analysis.R --stage=NAME\nStages: qc, annotation_evidence, annotation_audit, annotation, lineage, celltype, harmony, figure2, dotplot, paired_umap, merged_umap, overlap, validate\n--check checks dependencies; edit input/output paths in acti_config().\nAnnotation requires a reviewed map or --use-manuscript-annotation.\n")

acti_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- c("--help", "--check", "--use-manuscript-annotation")
  valued <- c("project-dir", "data-dir", "output-dir", "stage", "annotation-map", "annotated-rds", "hpca-rds",
    "case-dir", "control-dir", "stpmg-csv", "figure2-mode")
  for (arg in args) {
    if (arg %in% flags)
      next
    if (!grepl("^--[^=]+=.+$", arg) || !sub("^--([^=]+)=.*$", "\\1", arg) %in% valued)
      stop("Unknown/empty argument: ", arg, ". Use --help.")
  }
  value <- function(key, default = "") {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (length(hit) > 1L)
      stop("Duplicate argument: --", key)
    if (length(hit))
      substring(hit, nchar(key) + 4L)
    else default
  }
  if (!length(args) || "--help" %in% args) {
    acti_help()
    return(invisible(NULL))
  }
  root <- value("project-dir", ".")
  cfg <- acti_config(root, value("data-dir", file.path(root, "data")), value("output-dir", file.path(root, "results")))
  for (key in c("annotation-map", "annotated-rds", "hpca-rds", "case-dir", "control-dir", "stpmg-csv", "figure2-mode")) {
    v <- value(key)
    if (nzchar(v))
      cfg[[gsub("-", "_", key)]] <- v
  }
  cfg$figure2_mode <- toupper(cfg$figure2_mode)
  if (!cfg$figure2_mode %in% c("ALL", "ABCD", "BC", "E", "F"))
    stop("Unsupported Figure 2 mode.")
  cfg$use_manuscript_annotation <- "--use-manuscript-annotation" %in% args
  stage_arg <- value("stage")
  stages <- if (nzchar(stage_arg))
    strsplit(stage_arg, ",", fixed = TRUE)[[1L]]
  else character()
  if ("--check" %in% args) {
    acti_check(cfg, if (length(stages) && !identical(stages, "all"))
      stages
    else names(acti_stage_packages))
    return(invisible(cfg))
  }
  if (!length(stages)) {
    acti_help()
    return(invisible(cfg))
  }
  run_acti_scrna(stages, cfg)
}

if (sys.nframe() == 0L) acti_main()

