# Paired WGCNA and shared module genes.
# Inputs: GSE89403 expression and phenotype tables; paths in acti_wgcna_config().
# Run: Rscript 2_WGCNA_analysis.R --stage week4
# source() loads functions without starting an analysis.

acti_wgcna_config <- function(data_dir = file.path("data", "GSE89403"), output_dir = file.path("results", "WGCNA"),
  expr_file = NULL, pheno_file = NULL, week4_dir = NULL, week24_dir = NULL, overwrite = FALSE) {
  list(expr_file = if (is.null(expr_file)) file.path(data_dir, "GSE89403_log2ExpGeneNames.csv") else expr_file,
    pheno_file = if (is.null(pheno_file)) file.path(data_dir, "GSE89403_matrix_pheno.csv") else pheno_file,
    output_dir = output_dir, week4_dir = if (is.null(week4_dir)) file.path(output_dir, "01_DX_week4") else week4_dir,
    week24_dir = if (is.null(week24_dir)) file.path(output_dir, "02_DX_week24") else week24_dir, overwrite = isTRUE(overwrite))
}

acti_wgcna_check <- function(cfg = acti_wgcna_config()) {
  packages <- c("WGCNA", "dplyr", "pheatmap", "RColorBrewer", "ggplot2")
  available <- vapply(packages, function(p) {
    nzchar(system.file(package = p))
  }, logical(1))
  inputs <- c(expression = cfg$expr_file, phenotype = cfg$pheno_file)
  print(data.frame(package = packages, available = available), row.names = FALSE)
  print(data.frame(input = names(inputs), path = unname(inputs), exists = file.exists(inputs)), row.names = FALSE)
  invisible(all(available) && all(file.exists(inputs)))
}

acti_wgcna_week4 <- function(cfg = acti_wgcna_config()) {
  suppressPackageStartupMessages({
    library(WGCNA)
  })
  options(stringsAsFactors = FALSE)
  allowWGCNAThreads()
  set.seed(1234)
  project_dir <- cfg$output_dir
  expr_file <- cfg$expr_file
  pheno_file <- cfg$pheno_file
  outdir <- cfg$week4_dir
  overwrite_existing_output <- cfg$overwrite
  if (dir.exists(outdir) && length(list.files(outdir, all.files = FALSE)) > 0L && !overwrite_existing_output) {
    stop("Output directory is not empty: ", outdir, ". Use a new directory or explicitly set overwrite_existing_output <- TRUE.")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  network_type <- "signed"
  tom_type <- "signed"
  cor_type <- "bicor"
  max_p_outliers <- 0.05
  pearson_fallback <- "individual"
  powers <- c(1:10, seq(12, 20, 2))
  min_module_size <- 30
  merge_cut_height <- 0.25
  mad_keep_fraction <- 0.5
  scale_free_R2_cutoff <- 0.85
  expected_n_pairs <- 67L
  expected_time1_only <- c("S134", "S144", "S85", "S86")
  expected_time2_only <- c("S130", "S4", "S5")
  figure_cex <- 1.15
  figure_cex_main <- 1.45
  figure_cex_lab <- 1.3
  figure_cex_axis <- 1.15
  figure_cex_sub <- 1.1
  sample_label_cex <- 0.72
  annotation_text_cex <- 0.9
  module_label_cex <- 1.15
  remove_outlier_pairs <- character(0)
  run_analysis <- function() {
    log_file <- file.path(outdir, "analysis_log.txt")
    sink(log_file, split = TRUE)
    on.exit({
      sink()
    }, add = TRUE)
    cat("Analysis started:", format(Sys.time()), "\n")
    cat("Expression file:", expr_file, "\n")
    cat("Phenotype file:", pheno_file, "\n")
    cat("Output directory:", outdir, "\n\n")
    cat("Correlation:", cor_type, "\n")
    cat("maxPOutliers:", max_p_outliers, "\n")
    cat("Pearson fallback:", pearson_fallback, "\n\n")
    if (!file.exists(expr_file))
      stop("Expression file not found: ", expr_file)
    if (!file.exists(pheno_file))
      stop("Phenotype file not found: ", pheno_file)
    expr_raw <- read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(expr_raw) < 3L) {
      stop("Expression file must contain Ensembl ID, symbol, and sample columns.")
    }
    ensembl_id <- trimws(as.character(expr_raw[[1]]))
    gene_symbol <- trimws(as.character(expr_raw[[2]]))
    has_symbol <- !is.na(gene_symbol) & nzchar(gene_symbol) & gene_symbol != "NA"
    gene_id <- ifelse(has_symbol, gene_symbol, ensembl_id)
    if (any(is.na(gene_id) | !nzchar(gene_id))) {
      stop("Some expression rows lack both gene symbol and Ensembl ID.")
    }
    gene_identifier_map <- data.frame(ensembl_id = ensembl_id, original_symbol = gene_symbol, analysis_gene_id = gene_id,
      identifier_source = ifelse(has_symbol, "symbol", "Ensembl_ID"), stringsAsFactors = FALSE)
    write.csv(gene_identifier_map, file.path(outdir, "gene_identifier_mapping.csv"), row.names = FALSE, quote = FALSE,
      na = "")
    expr <- as.matrix(expr_raw[, -(1:2), drop = FALSE])
    storage.mode(expr) <- "double"
    rownames(expr) <- gene_id
    rm(expr_raw)
    if (any(!is.finite(expr))) {
      stop("Expression matrix contains NA, NaN, or infinite values.")
    }
    cat("Article-supplied log2 expression range:", paste(range(expr), collapse = " to "), "\n")
    cat("Negative log2 values retained without offset; WGCNA correlation does not require nonnegative input.\n")
    valid_gene <- !is.na(rownames(expr)) & nzchar(trimws(rownames(expr)))
    expr <- expr[valid_gene, , drop = FALSE]
    date_like_gene_ids <- grep("^([0-9]{1,2}-[A-Za-z]{3}|[A-Za-z]{3}-[0-9]{1,2})$", rownames(expr), value = TRUE)
    if (length(date_like_gene_ids) > 0L) {
      warning(length(date_like_gene_ids), " gene identifiers look like spreadsheet-converted dates; see ",
        file.path(outdir, "date_like_gene_identifiers_to_check.txt"))
      writeLines(unique(date_like_gene_ids), file.path(outdir, "date_like_gene_identifiers_to_check.txt"))
    }
    if (anyDuplicated(rownames(expr))) {
      duplicate_counts <- table(rownames(expr))
      expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
      expr <- expr/as.numeric(duplicate_counts[rownames(expr)])
    }
    cat("Input matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
    pheno <- read.csv(pheno_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_pheno <- c("subject", "treatmentresult", "timetonegativity")
    if (!all(required_pheno %in% colnames(pheno))) {
      stop("Phenotype file lacks columns: ", paste(setdiff(required_pheno, colnames(pheno)), collapse = ", "))
    }
    valid_subject <- !is.na(pheno$subject) & nzchar(trimws(as.character(pheno$subject)))
    pheno <- pheno[valid_subject, required_pheno, drop = FALSE]
    pheno$subject <- as.character(pheno$subject)
    pheno$treatmentresult <- trimws(pheno$treatmentresult)
    clinical_rows <- unique(pheno)
    clinical_split <- split(clinical_rows, clinical_rows$subject, drop = TRUE)
    inconsistent <- names(clinical_split)[vapply(clinical_split, function(z) {
      x <- unique(z$treatmentresult[!is.na(z$treatmentresult) & nzchar(z$treatmentresult)])
      length(x) > 1L
    }, logical(1))]
    if (length(inconsistent) > 0L) {
      stop("Inconsistent treatment result within subjects: ", paste(inconsistent, collapse = ", "))
    }
    subject_clinical <- do.call(rbind, lapply(clinical_split, function(z) {
      treatment_result <- unique(z$treatmentresult[!is.na(z$treatmentresult) & nzchar(z$treatmentresult)])
      data.frame(subject = paste0("S", z$subject[1]), treatmentresult = if (length(treatment_result) == 0L) {
        NA_character_
      }
      else {
        treatment_result[1]
      }, timetonegativity = {
        x <- unique(z$timetonegativity[!is.na(z$timetonegativity)])
        if (length(x) == 0L)
         NA_character_
        else x[1]
      }, stringsAsFactors = FALSE)
    }))
    rownames(subject_clinical) <- NULL
    definite_cure_subjects <- subject_clinical$subject[!is.na(subject_clinical$treatmentresult) & subject_clinical$treatmentresult ==
      "Definite Cure"]
    all_sample_names <- colnames(expr)
    is_dx_or_week4 <- grepl("_(DX|week_4)$", all_sample_names)
    sample_subject <- sub("_(DX|week_4)$", "", all_sample_names)
    keep_sample <- is_dx_or_week4 & sample_subject %in% definite_cure_subjects
    expr <- expr[, keep_sample, drop = FALSE]
    cat("Definite Cure subjects in phenotype table:", length(definite_cure_subjects), "\n")
    cat("Definite Cure DX/week4 samples before pairing:", ncol(expr), "\n")
    sample_name <- colnames(expr)
    recognized <- grepl("_(DX|week_4)$", sample_name)
    if (!all(recognized)) {
      stop("Unrecognized sample-name suffixes: ", paste(sample_name[!recognized], collapse = ", "))
    }
    sample_info <- data.frame(sample = sample_name, subject = sub("_(DX|week_4)$", "", sample_name), timepoint = ifelse(grepl("_DX$",
      sample_name), "DX", "week4"), stringsAsFactors = FALSE)
    pair_counts <- table(sample_info$subject, sample_info$timepoint)
    if (!all(c("DX", "week4") %in% colnames(pair_counts))) {
      stop("Both DX and week4 labels were not found.")
    }
    complete_subjects <- rownames(pair_counts)[pair_counts[, "DX"] == 1L & pair_counts[, "week4"] == 1L]
    dx_only <- rownames(pair_counts)[pair_counts[, "DX"] == 1L & pair_counts[, "week4"] == 0L]
    week4_only <- rownames(pair_counts)[pair_counts[, "DX"] == 0L & pair_counts[, "week4"] == 1L]
    if (length(complete_subjects) != expected_n_pairs) {
      stop("Expected ", expected_n_pairs, " complete pairs, but found ", length(complete_subjects), ".")
    }
    if (!setequal(dx_only, expected_time1_only)) {
      stop("DX-only subjects differ from expectation. Detected: ", paste(sort(dx_only), collapse = ", "))
    }
    if (!setequal(week4_only, expected_time2_only)) {
      stop("week4-only subjects differ from expectation. Detected: ", paste(sort(week4_only), collapse = ", "))
    }
    if (length(remove_outlier_pairs) > 0L) {
      invalid_remove <- setdiff(remove_outlier_pairs, complete_subjects)
      if (length(invalid_remove) > 0L) {
        stop("Unknown outlier-pair IDs: ", paste(invalid_remove, collapse = ", "))
      }
      complete_subjects <- setdiff(complete_subjects, remove_outlier_pairs)
    }
    complete_subjects <- sort(complete_subjects)
    paired_samples <- as.vector(rbind(paste0(complete_subjects, "_DX"), paste0(complete_subjects, "_week_4")))
    stopifnot(all(paired_samples %in% colnames(expr)))
    sample_info_paired <- data.frame(sample = paired_samples, subject = rep(complete_subjects, each = 2L),
      timepoint = rep(c("DX", "week4"), times = length(complete_subjects)), time_numeric = rep(c(0, 1), times = length(complete_subjects)),
      stringsAsFactors = FALSE)
    clinical_index <- match(sample_info_paired$subject, subject_clinical$subject)
    sample_info_paired$treatmentresult <- subject_clinical$treatmentresult[clinical_index]
    sample_info_paired$timetonegativity <- subject_clinical$timetonegativity[clinical_index]
    stopifnot(all(sample_info_paired$treatmentresult == "Definite Cure"))
    rownames(sample_info_paired) <- sample_info_paired$sample
    expr <- expr[, paired_samples, drop = FALSE]
    cat("Complete pairs retained:", length(complete_subjects), "\n")
    cat("Samples retained:", ncol(expr), "\n")
    cat("Treatment result:", paste(names(table(sample_info_paired$treatmentresult)), as.integer(table(sample_info_paired$treatmentresult)),
      collapse = "; "), "\n")
    cat("Time to negativity among paired subjects:\n")
    print(table(sample_info_paired$timetonegativity[!duplicated(sample_info_paired$subject)], useNA = "ifany"))
    cat("DX-only excluded:", paste(sort(dx_only), collapse = ", "), "\n")
    cat("week4-only excluded:", paste(sort(week4_only), collapse = ", "), "\n\n")
    write.csv(sample_info_paired, file.path(outdir, "paired_sample_manifest.csv"), row.names = FALSE, quote = FALSE)
    write.csv(data.frame(subject = c(dx_only, week4_only), available_timepoint = c(rep("DX", length(dx_only)),
      rep("week4", length(week4_only)))), file.path(outdir, "excluded_unpaired_samples.csv"), row.names = FALSE,
      quote = FALSE)
    gene_mad <- apply(expr, 1, mad, na.rm = TRUE)
    mad_cutoff <- unname(quantile(gene_mad[is.finite(gene_mad)], probs = 1 - mad_keep_fraction, na.rm = TRUE))
    keep_gene <- is.finite(gene_mad) & gene_mad > mad_cutoff
    expr_filtered <- expr[keep_gene, , drop = FALSE]
    cat("MAD cutoff:", mad_cutoff, "\n")
    cat("Genes retained after top-50% MAD filter:", nrow(expr_filtered), "\n")
    datExpr <- as.data.frame(t(expr_filtered))
    gsg <- goodSamplesGenes(datExpr, verbose = 3)
    if (!gsg$allOK) {
      removed_genes <- colnames(datExpr)[!gsg$goodGenes]
      removed_samples <- rownames(datExpr)[!gsg$goodSamples]
      if (length(removed_samples) > 0L) {
        stop("goodSamplesGenes would remove individual samples and break pairing: ", paste(removed_samples,
         collapse = ", "), ". Investigate and remove complete subject pairs explicitly.")
      }
      datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
      writeLines(removed_genes, file.path(outdir, "genes_removed_by_WGCNA_QC.txt"))
    }
    stopifnot(identical(rownames(datExpr), paired_samples))
    nSamples <- nrow(datExpr)
    nGenes <- ncol(datExpr)
    cat("Final WGCNA matrix:", nSamples, "samples x", nGenes, "genes\n\n")
    sample_cor <- bicor(t(as.matrix(datExpr)), use = "pairwise.complete.obs", maxPOutliers = max_p_outliers,
      pearsonFallback = pearson_fallback)
    diag(sample_cor) <- 0
    sample_connectivity <- rowSums(abs(sample_cor), na.rm = TRUE)
    sample_connectivity_Z <- as.numeric(scale(sample_connectivity))
    sample_connectivity_qc <- data.frame(sample = names(sample_connectivity), subject = sub("_(DX|week_4)$",
      "", names(sample_connectivity)), timepoint = ifelse(grepl("_DX$", names(sample_connectivity)), "DX",
      "week4"), connectivity = sample_connectivity, connectivity_Z = sample_connectivity_Z, flag_Z_below_minus2_5 = sample_connectivity_Z <
      -2.5, stringsAsFactors = FALSE)
    write.csv(sample_connectivity_qc, file.path(outdir, "sample_connectivity_QC.csv"), row.names = FALSE, quote = FALSE)
    flagged_samples <- sample_connectivity_qc$sample[sample_connectivity_qc$flag_Z_below_minus2_5]
    cat("Sample-connectivity flags (Z < -2.5; no automatic removal):", if (length(flagged_samples))
      paste(flagged_samples, collapse = ", ")
    else "none", "\n\n")
    pdf(file.path(outdir, "01_sample_clustering.pdf"), width = 16, height = 9)
    par(cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis,
      cex.sub = figure_cex_sub, mar = c(6, 5, 5, 2))
    sampleTree <- hclust(dist(datExpr), method = "average")
    plot(sampleTree, main = paste0("Sample clustering: ", length(complete_subjects), " complete DX-week4 pairs"),
      sub = "No automatic outlier removal", xlab = "", cex = sample_label_cex)
    dev.off()
    sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = scale_free_R2_cutoff, networkType = network_type,
      corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = max_p_outliers, pearsonFallback = pearson_fallback),
      verbose = 5)
    power <- sft$powerEstimate
    if (is.na(power)) {
      power <- if (network_type == "unsigned")
        6
      else 12
    }
    cat("Automatic soft-threshold estimate:", sft$powerEstimate, "\n")
    cat("Final soft-threshold power:", power, "\n\n")
    write.csv(sft$fitIndices, file.path(outdir, "soft_threshold_fit_indices.csv"), row.names = FALSE)
    pdf(file.path(outdir, "02_soft_threshold_selection.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2), cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis,
      mar = c(5.5, 5.5, 4.5, 2))
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft-threshold power",
      ylab = "Signed scale-free topology fit (R^2)", type = "n", main = "Scale independence")
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = sft$fitIndices[, 1],
      col = "red", cex = annotation_text_cex)
    abline(h = scale_free_R2_cutoff, col = "red", lty = 2)
    abline(v = power, col = "blue", lty = 2)
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft-threshold power", ylab = "Mean connectivity",
      type = "n", main = "Mean connectivity")
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = sft$fitIndices[, 1], col = "red", cex = annotation_text_cex)
    abline(v = power, col = "blue", lty = 2)
    dev.off()
    net <- blockwiseModules(datExpr, power = power, maxBlockSize = nGenes, TOMType = tom_type, networkType = network_type,
      minModuleSize = min_module_size, reassignThreshold = 0, mergeCutHeight = merge_cut_height, numericLabels = TRUE,
      pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = file.path(outdir, "paired_DX_week4_bicor_signed_merge025_TOM"),
      corType = cor_type, maxPOutliers = max_p_outliers, pearsonFallback = pearson_fallback, useCorOptionsThroughout = TRUE,
      verbose = 3)
    moduleLabels <- net$colors
    moduleColors <- labels2colors(moduleLabels)
    MEs <- orderMEs(moduleEigengenes(datExpr, colors = moduleColors)$eigengenes)
    module_size <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
    colnames(module_size) <- c("moduleColor", "nGenes")
    write.csv(module_size, file.path(outdir, "module_sizes.csv"), row.names = FALSE)
    pdf(file.path(outdir, "03_gene_dendrogram_and_modules.pdf"), width = 16, height = 9)
    par(cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis)
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE,
      cex.colorLabels = module_label_cex, cex.dendroLabels = figure_cex_axis, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05)
    dev.off()
    pdf(file.path(outdir, "04_eigengene_network.pdf"), width = 14, height = 12)
    par(cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis)
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marDendro = c(3, 3, 2, 4), marHeatmap = c(3,
      4, 2, 2), plotDendrograms = TRUE, xLabelsAngle = 90, cex.adjacency = annotation_text_cex, cex.preservation = annotation_text_cex)
    dev.off()
    MEs <- MEs[paired_samples, , drop = FALSE]
    dx_rows <- paste0(complete_subjects, "_DX")
    w4_rows <- paste0(complete_subjects, "_week_4")
    ME_delta <- as.matrix(MEs[w4_rows, , drop = FALSE]) - as.matrix(MEs[dx_rows, , drop = FALSE])
    rownames(ME_delta) <- complete_subjects
    paired_results <- do.call(rbind, lapply(seq_len(ncol(ME_delta)), function(j) {
      delta <- ME_delta[, j]
      delta_ok <- delta[is.finite(delta)]
      if (length(delta_ok) >= 3L && sd(delta_ok) > 0) {
        tt <- t.test(delta_ok, mu = 0)
        t_stat <- unname(tt$statistic)
        df_value <- unname(tt$parameter)
        p_value <- tt$p.value
        ci_low <- unname(tt$conf.int[1])
        ci_high <- unname(tt$conf.int[2])
      }
      else {
        t_stat <- df_value <- p_value <- ci_low <- ci_high <- NA_real_
      }
      data.frame(module = colnames(ME_delta)[j], moduleColor = sub("^ME", "", colnames(ME_delta)[j]), nPairs = sum(is.finite(delta)),
        mean_DX = mean(MEs[dx_rows, j], na.rm = TRUE), mean_week4 = mean(MEs[w4_rows, j], na.rm = TRUE),
        mean_week4_minus_DX = mean(delta, na.rm = TRUE), sd_difference = sd(delta, na.rm = TRUE), cohen_dz = mean(delta,
         na.rm = TRUE)/sd(delta, na.rm = TRUE), t = t_stat, df = df_value, pValue = p_value, ci95_low = ci_low,
        ci95_high = ci_high, stringsAsFactors = FALSE)
    }))
    paired_results$FDR <- p.adjust(paired_results$pValue, method = "BH")
    wilcox_p <- apply(ME_delta, 2, function(delta) {
      delta <- delta[is.finite(delta)]
      if (length(delta) < 3L || all(delta == 0))
        return(NA_real_)
      suppressWarnings(wilcox.test(delta, mu = 0, exact = FALSE)$p.value)
    })
    paired_results$wilcox_pValue <- wilcox_p[paired_results$module]
    paired_results$wilcox_FDR <- p.adjust(paired_results$wilcox_pValue, method = "BH")
    paired_results$direction <- ifelse(paired_results$mean_week4_minus_DX > 0, "higher_at_week4", "lower_at_week4")
    paired_results <- paired_results[order(paired_results$FDR, paired_results$pValue), ]
    write.csv(paired_results, file.path(outdir, "paired_module_time_association_all.csv"), row.names = FALSE)
    write.csv(subset(paired_results, FDR < 0.05), file.path(outdir, "paired_module_time_association_FDR05.csv"),
      row.names = FALSE)
    write.csv(data.frame(subject = rownames(ME_delta), ME_delta, check.names = FALSE), file.path(outdir, "module_eigengene_week4_minus_DX_by_subject.csv"),
      row.names = FALSE)
    non_grey_results <- paired_results[paired_results$moduleColor != "grey", ]
    top5_modules <- head(non_grey_results$module, 5L)
    write.csv(head(non_grey_results, 5L), file.path(outdir, "top5_modules_by_paired_test.csv"), row.names = FALSE)
    heat_order <- paired_results$module
    heat_effect <- matrix(paired_results$cohen_dz, ncol = 1, dimnames = list(heat_order, "week4 - DX"))
    heat_text <- matrix(paste0("dz=", sprintf("%.2f", paired_results$cohen_dz), "\nP=", format.pval(paired_results$pValue,
      digits = 2), "\nFDR=", format.pval(paired_results$FDR, digits = 2)), ncol = 1, dimnames = dimnames(heat_effect))
    heat_limit <- max(1, max(abs(heat_effect), na.rm = TRUE))
    pdf(file.path(outdir, "05_paired_module_time_association_heatmap.pdf"), width = 9, height = max(8, 0.7 *
      nrow(heat_effect) + 2.5))
    par(mar = c(8, 11, 5, 4), cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis)
    labeledHeatmap(Matrix = heat_effect, xLabels = colnames(heat_effect), yLabels = rownames(heat_effect),
      ySymbols = rownames(heat_effect), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = heat_text,
      setStdMargins = FALSE, cex.text = annotation_text_cex, cex.lab = module_label_cex, cex.legendLabel = figure_cex_axis,
      zlim = c(-heat_limit, heat_limit), main = paste0("Paired module-time association (", length(complete_subjects),
        " pairs)"))
    dev.off()
    kME <- signedKME(datExpr, MEs, corFnc = "bicor", corOptions = paste0("use = 'p', maxPOutliers = ", max_p_outliers,
      ", pearsonFallback = '", pearson_fallback, "'"))
    own_kME <- vapply(seq_len(nGenes), function(i) {
      kme_name <- paste0("kME", moduleColors[i])
      if (kme_name %in% colnames(kME))
        kME[i, kme_name]
      else NA_real_
    }, numeric(1))
    geneInfo <- data.frame(gene = colnames(datExpr), moduleColor = moduleColors, moduleLabel = moduleLabels,
      kME_own_module = own_kME, stringsAsFactors = FALSE)
    geneInfo <- merge(geneInfo, paired_results[, c("moduleColor", "mean_week4_minus_DX", "cohen_dz", "pValue",
      "FDR", "wilcox_pValue", "wilcox_FDR", "direction")], by = "moduleColor", all.x = TRUE, sort = FALSE)
    geneInfo <- geneInfo[order(geneInfo$moduleLabel, -abs(geneInfo$kME_own_module)), ]
    write.csv(geneInfo, file.path(outdir, "all_gene_module_information.csv"), row.names = FALSE)
    top5_gene_info <- geneInfo[paste0("ME", geneInfo$moduleColor) %in% top5_modules, ]
    write.csv(top5_gene_info, file.path(outdir, "top5_paired_modules_combined_genes.csv"), row.names = FALSE)
    for (module in top5_modules) {
      color <- sub("^ME", "", module)
      genes <- geneInfo$gene[geneInfo$moduleColor == color]
      writeLines(genes, file.path(outdir, paste0("module_", module, "_genes.txt")))
    }
    target_genes <- c("FCGR1A", "FCGR1B", "ANKRD22", "DUSP3", "GBP5", "VAMP5", "PSTPIP2", "GBP1", "NPC2", "TIFA",
      "GBP2", "LHFPL2", "GYG1", "FLVCR2", "SECTM1", "POMP", "FBXO6", "GADD45B", "RHBDF2", "CASP4")
    target_gene_results <- geneInfo[geneInfo$gene %in% target_genes, c("gene", "moduleColor", "moduleLabel",
      "kME_own_module", "mean_week4_minus_DX", "cohen_dz", "pValue", "FDR", "wilcox_pValue", "wilcox_FDR",
      "direction")]
    target_gene_results$in_top5_paired_modules <- paste0("ME", target_gene_results$moduleColor) %in% top5_modules
    write.csv(target_gene_results, file.path(outdir, "target_gene_module_lookup.csv"), row.names = FALSE)
    saveRDS(list(parameters = list(network_type = network_type, cor_type = cor_type, max_p_outliers = max_p_outliers,
      pearson_fallback = pearson_fallback, power = power, min_module_size = min_module_size, merge_cut_height = merge_cut_height,
      mad_keep_fraction = mad_keep_fraction, scale_free_R2_cutoff = scale_free_R2_cutoff), sample_info = sample_info_paired,
      datExpr = datExpr, net = net, MEs = MEs, ME_delta = ME_delta, paired_results = paired_results, geneInfo = geneInfo),
      file.path(outdir, "paired_WGCNA_complete_results.rds"))
    writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
    cat("\nAnalysis completed:", format(Sys.time()), "\n")
    cat("Number of subjects:", length(complete_subjects), "\n")
    cat("Number of samples:", nSamples, "\n")
    cat("Number of genes:", nGenes, "\n")
    cat("Number of non-grey modules:", sum(unique(moduleColors) != "grey"), "\n")
    cat("FDR < 0.05 paired time-associated modules:", sum(paired_results$FDR < 0.05), "\n")
    cat("Top five exploratory modules:", paste(top5_modules, collapse = ", "), "\n")
    cat("Results directory:", outdir, "\n")
  }
  run_analysis()
}

acti_wgcna_week24 <- function(cfg = acti_wgcna_config()) {
  suppressPackageStartupMessages({
    library(WGCNA)
  })
  options(stringsAsFactors = FALSE)
  allowWGCNAThreads()
  set.seed(1234)
  project_dir <- cfg$output_dir
  expr_file <- cfg$expr_file
  pheno_file <- cfg$pheno_file
  outdir <- cfg$week24_dir
  overwrite_existing_output <- cfg$overwrite
  if (dir.exists(outdir) && length(list.files(outdir, all.files = FALSE)) > 0L && !overwrite_existing_output) {
    stop("Output directory is not empty: ", outdir, ". Use a new directory or explicitly set overwrite_existing_output <- TRUE.")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  network_type <- "signed"
  tom_type <- "signed"
  cor_type <- "bicor"
  max_p_outliers <- 0.05
  pearson_fallback <- "individual"
  powers <- c(1:10, seq(12, 20, 2))
  min_module_size <- 30
  merge_cut_height <- 0.25
  mad_keep_fraction <- 0.5
  scale_free_R2_cutoff <- 0.85
  expected_n_pairs <- 70L
  expected_time1_only <- c("S152")
  expected_time2_only <- c("S130", "S4", "S5")
  figure_cex <- 1.15
  figure_cex_main <- 1.45
  figure_cex_lab <- 1.3
  figure_cex_axis <- 1.15
  figure_cex_sub <- 1.1
  sample_label_cex <- 0.72
  annotation_text_cex <- 0.9
  module_label_cex <- 1.15
  remove_outlier_pairs <- character(0)
  run_analysis <- function() {
    log_file <- file.path(outdir, "analysis_log.txt")
    sink(log_file, split = TRUE)
    on.exit({
      sink()
    }, add = TRUE)
    cat("Analysis started:", format(Sys.time()), "\n")
    cat("Expression file:", expr_file, "\n")
    cat("Phenotype file:", pheno_file, "\n")
    cat("Output directory:", outdir, "\n\n")
    cat("Correlation:", cor_type, "\n")
    cat("maxPOutliers:", max_p_outliers, "\n")
    cat("Pearson fallback:", pearson_fallback, "\n\n")
    if (!file.exists(expr_file))
      stop("Expression file not found: ", expr_file)
    if (!file.exists(pheno_file))
      stop("Phenotype file not found: ", pheno_file)
    expr_raw <- read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(expr_raw) < 3L) {
      stop("Expression file must contain Ensembl ID, symbol, and sample columns.")
    }
    ensembl_id <- trimws(as.character(expr_raw[[1]]))
    gene_symbol <- trimws(as.character(expr_raw[[2]]))
    has_symbol <- !is.na(gene_symbol) & nzchar(gene_symbol) & gene_symbol != "NA"
    gene_id <- ifelse(has_symbol, gene_symbol, ensembl_id)
    if (any(is.na(gene_id) | !nzchar(gene_id))) {
      stop("Some expression rows lack both gene symbol and Ensembl ID.")
    }
    gene_identifier_map <- data.frame(ensembl_id = ensembl_id, original_symbol = gene_symbol, analysis_gene_id = gene_id,
      identifier_source = ifelse(has_symbol, "symbol", "Ensembl_ID"), stringsAsFactors = FALSE)
    write.csv(gene_identifier_map, file.path(outdir, "gene_identifier_mapping.csv"), row.names = FALSE, quote = FALSE,
      na = "")
    expr <- as.matrix(expr_raw[, -(1:2), drop = FALSE])
    storage.mode(expr) <- "double"
    rownames(expr) <- gene_id
    rm(expr_raw)
    if (any(!is.finite(expr))) {
      stop("Expression matrix contains NA, NaN, or infinite values.")
    }
    cat("Article-supplied log2 expression range:", paste(range(expr), collapse = " to "), "\n")
    cat("Negative log2 values retained without offset; WGCNA correlation does not require nonnegative input.\n")
    valid_gene <- !is.na(rownames(expr)) & nzchar(trimws(rownames(expr)))
    expr <- expr[valid_gene, , drop = FALSE]
    date_like_gene_ids <- grep("^([0-9]{1,2}-[A-Za-z]{3}|[A-Za-z]{3}-[0-9]{1,2})$", rownames(expr), value = TRUE)
    if (length(date_like_gene_ids) > 0L) {
      warning(length(date_like_gene_ids), " gene identifiers look like spreadsheet-converted dates; see ",
        file.path(outdir, "date_like_gene_identifiers_to_check.txt"))
      writeLines(unique(date_like_gene_ids), file.path(outdir, "date_like_gene_identifiers_to_check.txt"))
    }
    if (anyDuplicated(rownames(expr))) {
      duplicate_counts <- table(rownames(expr))
      expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
      expr <- expr/as.numeric(duplicate_counts[rownames(expr)])
    }
    cat("Input matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
    pheno <- read.csv(pheno_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_pheno <- c("subject", "treatmentresult", "timetonegativity")
    if (!all(required_pheno %in% colnames(pheno))) {
      stop("Phenotype file lacks columns: ", paste(setdiff(required_pheno, colnames(pheno)), collapse = ", "))
    }
    valid_subject <- !is.na(pheno$subject) & nzchar(trimws(as.character(pheno$subject)))
    pheno <- pheno[valid_subject, required_pheno, drop = FALSE]
    pheno$subject <- as.character(pheno$subject)
    pheno$treatmentresult <- trimws(pheno$treatmentresult)
    clinical_rows <- unique(pheno)
    clinical_split <- split(clinical_rows, clinical_rows$subject, drop = TRUE)
    inconsistent <- names(clinical_split)[vapply(clinical_split, function(z) {
      x <- unique(z$treatmentresult[!is.na(z$treatmentresult) & nzchar(z$treatmentresult)])
      length(x) > 1L
    }, logical(1))]
    if (length(inconsistent) > 0L) {
      stop("Inconsistent treatment result within subjects: ", paste(inconsistent, collapse = ", "))
    }
    subject_clinical <- do.call(rbind, lapply(clinical_split, function(z) {
      treatment_result <- unique(z$treatmentresult[!is.na(z$treatmentresult) & nzchar(z$treatmentresult)])
      data.frame(subject = paste0("S", z$subject[1]), treatmentresult = if (length(treatment_result) == 0L) {
        NA_character_
      }
      else {
        treatment_result[1]
      }, timetonegativity = {
        x <- unique(z$timetonegativity[!is.na(z$timetonegativity)])
        if (length(x) == 0L)
         NA_character_
        else x[1]
      }, stringsAsFactors = FALSE)
    }))
    rownames(subject_clinical) <- NULL
    definite_cure_subjects <- subject_clinical$subject[!is.na(subject_clinical$treatmentresult) & subject_clinical$treatmentresult ==
      "Definite Cure"]
    all_sample_names <- colnames(expr)
    is_dx_or_week24 <- grepl("_(DX|week_24)$", all_sample_names)
    sample_subject <- sub("_(DX|week_24)$", "", all_sample_names)
    keep_sample <- is_dx_or_week24 & sample_subject %in% definite_cure_subjects
    expr <- expr[, keep_sample, drop = FALSE]
    cat("Definite Cure subjects in phenotype table:", length(definite_cure_subjects), "\n")
    cat("Definite Cure DX/week24 samples before pairing:", ncol(expr), "\n")
    sample_name <- colnames(expr)
    recognized <- grepl("_(DX|week_24)$", sample_name)
    if (!all(recognized)) {
      stop("Unrecognized sample-name suffixes: ", paste(sample_name[!recognized], collapse = ", "))
    }
    sample_info <- data.frame(sample = sample_name, subject = sub("_(DX|week_24)$", "", sample_name), timepoint = ifelse(grepl("_DX$",
      sample_name), "DX", "week24"), stringsAsFactors = FALSE)
    pair_counts <- table(sample_info$subject, sample_info$timepoint)
    if (!all(c("DX", "week24") %in% colnames(pair_counts))) {
      stop("Both DX and week24 labels were not found.")
    }
    complete_subjects <- rownames(pair_counts)[pair_counts[, "DX"] == 1L & pair_counts[, "week24"] == 1L]
    dx_only <- rownames(pair_counts)[pair_counts[, "DX"] == 1L & pair_counts[, "week24"] == 0L]
    week24_only <- rownames(pair_counts)[pair_counts[, "DX"] == 0L & pair_counts[, "week24"] == 1L]
    if (length(complete_subjects) != expected_n_pairs) {
      stop("Expected ", expected_n_pairs, " complete pairs, but found ", length(complete_subjects), ".")
    }
    if (!setequal(dx_only, expected_time1_only)) {
      stop("DX-only subjects differ from expectation. Detected: ", paste(sort(dx_only), collapse = ", "))
    }
    if (!setequal(week24_only, expected_time2_only)) {
      stop("week24-only subjects differ from expectation. Detected: ", paste(sort(week24_only), collapse = ", "))
    }
    if (length(remove_outlier_pairs) > 0L) {
      invalid_remove <- setdiff(remove_outlier_pairs, complete_subjects)
      if (length(invalid_remove) > 0L) {
        stop("Unknown outlier-pair IDs: ", paste(invalid_remove, collapse = ", "))
      }
      complete_subjects <- setdiff(complete_subjects, remove_outlier_pairs)
    }
    complete_subjects <- sort(complete_subjects)
    paired_samples <- as.vector(rbind(paste0(complete_subjects, "_DX"), paste0(complete_subjects, "_week_24")))
    stopifnot(all(paired_samples %in% colnames(expr)))
    sample_info_paired <- data.frame(sample = paired_samples, subject = rep(complete_subjects, each = 2L),
      timepoint = rep(c("DX", "week24"), times = length(complete_subjects)), time_numeric = rep(c(0, 1),
        times = length(complete_subjects)), stringsAsFactors = FALSE)
    clinical_index <- match(sample_info_paired$subject, subject_clinical$subject)
    sample_info_paired$treatmentresult <- subject_clinical$treatmentresult[clinical_index]
    sample_info_paired$timetonegativity <- subject_clinical$timetonegativity[clinical_index]
    stopifnot(all(sample_info_paired$treatmentresult == "Definite Cure"))
    rownames(sample_info_paired) <- sample_info_paired$sample
    expr <- expr[, paired_samples, drop = FALSE]
    cat("Complete pairs retained:", length(complete_subjects), "\n")
    cat("Samples retained:", ncol(expr), "\n")
    cat("Treatment result:", paste(names(table(sample_info_paired$treatmentresult)), as.integer(table(sample_info_paired$treatmentresult)),
      collapse = "; "), "\n")
    cat("Time to negativity among paired subjects:\n")
    print(table(sample_info_paired$timetonegativity[!duplicated(sample_info_paired$subject)], useNA = "ifany"))
    cat("DX-only excluded:", paste(sort(dx_only), collapse = ", "), "\n")
    cat("week24-only excluded:", paste(sort(week24_only), collapse = ", "), "\n\n")
    write.csv(sample_info_paired, file.path(outdir, "paired_sample_manifest.csv"), row.names = FALSE, quote = FALSE)
    write.csv(data.frame(subject = c(dx_only, week24_only), available_timepoint = c(rep("DX", length(dx_only)),
      rep("week24", length(week24_only)))), file.path(outdir, "excluded_unpaired_samples.csv"), row.names = FALSE,
      quote = FALSE)
    gene_mad <- apply(expr, 1, mad, na.rm = TRUE)
    mad_cutoff <- unname(quantile(gene_mad[is.finite(gene_mad)], probs = 1 - mad_keep_fraction, na.rm = TRUE))
    keep_gene <- is.finite(gene_mad) & gene_mad > mad_cutoff
    expr_filtered <- expr[keep_gene, , drop = FALSE]
    cat("MAD cutoff:", mad_cutoff, "\n")
    cat("Genes retained after top-50% MAD filter:", nrow(expr_filtered), "\n")
    datExpr <- as.data.frame(t(expr_filtered))
    gsg <- goodSamplesGenes(datExpr, verbose = 3)
    if (!gsg$allOK) {
      removed_genes <- colnames(datExpr)[!gsg$goodGenes]
      removed_samples <- rownames(datExpr)[!gsg$goodSamples]
      if (length(removed_samples) > 0L) {
        stop("goodSamplesGenes would remove individual samples and break pairing: ", paste(removed_samples,
         collapse = ", "), ". Investigate and remove complete subject pairs explicitly.")
      }
      datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
      writeLines(removed_genes, file.path(outdir, "genes_removed_by_WGCNA_QC.txt"))
    }
    stopifnot(identical(rownames(datExpr), paired_samples))
    nSamples <- nrow(datExpr)
    nGenes <- ncol(datExpr)
    cat("Final WGCNA matrix:", nSamples, "samples x", nGenes, "genes\n\n")
    sample_cor <- bicor(t(as.matrix(datExpr)), use = "pairwise.complete.obs", maxPOutliers = max_p_outliers,
      pearsonFallback = pearson_fallback)
    diag(sample_cor) <- 0
    sample_connectivity <- rowSums(abs(sample_cor), na.rm = TRUE)
    sample_connectivity_Z <- as.numeric(scale(sample_connectivity))
    sample_connectivity_qc <- data.frame(sample = names(sample_connectivity), subject = sub("_(DX|week_24)$",
      "", names(sample_connectivity)), timepoint = ifelse(grepl("_DX$", names(sample_connectivity)), "DX",
      "week24"), connectivity = sample_connectivity, connectivity_Z = sample_connectivity_Z, flag_Z_below_minus2_5 = sample_connectivity_Z <
      -2.5, stringsAsFactors = FALSE)
    write.csv(sample_connectivity_qc, file.path(outdir, "sample_connectivity_QC.csv"), row.names = FALSE, quote = FALSE)
    flagged_samples <- sample_connectivity_qc$sample[sample_connectivity_qc$flag_Z_below_minus2_5]
    cat("Sample-connectivity flags (Z < -2.5; no automatic removal):", if (length(flagged_samples))
      paste(flagged_samples, collapse = ", ")
    else "none", "\n\n")
    pdf(file.path(outdir, "01_sample_clustering.pdf"), width = 16, height = 9)
    par(cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis,
      cex.sub = figure_cex_sub, mar = c(6, 5, 5, 2))
    sampleTree <- hclust(dist(datExpr), method = "average")
    plot(sampleTree, main = paste0("Sample clustering: ", length(complete_subjects), " complete DX-week24 pairs"),
      sub = "No automatic outlier removal", xlab = "", cex = sample_label_cex)
    dev.off()
    sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = scale_free_R2_cutoff, networkType = network_type,
      corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = max_p_outliers, pearsonFallback = pearson_fallback),
      verbose = 5)
    power <- sft$powerEstimate
    if (is.na(power)) {
      power <- if (network_type == "unsigned")
        6
      else 12
    }
    cat("Automatic soft-threshold estimate:", sft$powerEstimate, "\n")
    cat("Final soft-threshold power:", power, "\n\n")
    write.csv(sft$fitIndices, file.path(outdir, "soft_threshold_fit_indices.csv"), row.names = FALSE)
    pdf(file.path(outdir, "02_soft_threshold_selection.pdf"), width = 14, height = 7)
    par(mfrow = c(1, 2), cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis,
      mar = c(5.5, 5.5, 4.5, 2))
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft-threshold power",
      ylab = "Signed scale-free topology fit (R^2)", type = "n", main = "Scale independence")
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = sft$fitIndices[, 1],
      col = "red", cex = annotation_text_cex)
    abline(h = scale_free_R2_cutoff, col = "red", lty = 2)
    abline(v = power, col = "blue", lty = 2)
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft-threshold power", ylab = "Mean connectivity",
      type = "n", main = "Mean connectivity")
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = sft$fitIndices[, 1], col = "red", cex = annotation_text_cex)
    abline(v = power, col = "blue", lty = 2)
    dev.off()
    net <- blockwiseModules(datExpr, power = power, maxBlockSize = nGenes, TOMType = tom_type, networkType = network_type,
      minModuleSize = min_module_size, reassignThreshold = 0, mergeCutHeight = merge_cut_height, numericLabels = TRUE,
      pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = file.path(outdir, "paired_DX_week24_bicor_signed_merge025_TOM"),
      corType = cor_type, maxPOutliers = max_p_outliers, pearsonFallback = pearson_fallback, useCorOptionsThroughout = TRUE,
      verbose = 3)
    moduleLabels <- net$colors
    moduleColors <- labels2colors(moduleLabels)
    MEs <- orderMEs(moduleEigengenes(datExpr, colors = moduleColors)$eigengenes)
    module_size <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
    colnames(module_size) <- c("moduleColor", "nGenes")
    write.csv(module_size, file.path(outdir, "module_sizes.csv"), row.names = FALSE)
    pdf(file.path(outdir, "03_gene_dendrogram_and_modules.pdf"), width = 16, height = 9)
    par(cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis)
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE,
      cex.colorLabels = module_label_cex, cex.dendroLabels = figure_cex_axis, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05)
    dev.off()
    pdf(file.path(outdir, "04_eigengene_network.pdf"), width = 14, height = 12)
    par(cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis)
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marDendro = c(3, 3, 2, 4), marHeatmap = c(3,
      4, 2, 2), plotDendrograms = TRUE, xLabelsAngle = 90, cex.adjacency = annotation_text_cex, cex.preservation = annotation_text_cex)
    dev.off()
    MEs <- MEs[paired_samples, , drop = FALSE]
    dx_rows <- paste0(complete_subjects, "_DX")
    w24_rows <- paste0(complete_subjects, "_week_24")
    ME_delta <- as.matrix(MEs[w24_rows, , drop = FALSE]) - as.matrix(MEs[dx_rows, , drop = FALSE])
    rownames(ME_delta) <- complete_subjects
    paired_results <- do.call(rbind, lapply(seq_len(ncol(ME_delta)), function(j) {
      delta <- ME_delta[, j]
      delta_ok <- delta[is.finite(delta)]
      if (length(delta_ok) >= 3L && sd(delta_ok) > 0) {
        tt <- t.test(delta_ok, mu = 0)
        t_stat <- unname(tt$statistic)
        df_value <- unname(tt$parameter)
        p_value <- tt$p.value
        ci_low <- unname(tt$conf.int[1])
        ci_high <- unname(tt$conf.int[2])
      }
      else {
        t_stat <- df_value <- p_value <- ci_low <- ci_high <- NA_real_
      }
      data.frame(module = colnames(ME_delta)[j], moduleColor = sub("^ME", "", colnames(ME_delta)[j]), nPairs = sum(is.finite(delta)),
        mean_DX = mean(MEs[dx_rows, j], na.rm = TRUE), mean_week24 = mean(MEs[w24_rows, j], na.rm = TRUE),
        mean_week24_minus_DX = mean(delta, na.rm = TRUE), sd_difference = sd(delta, na.rm = TRUE), cohen_dz = mean(delta,
         na.rm = TRUE)/sd(delta, na.rm = TRUE), t = t_stat, df = df_value, pValue = p_value, ci95_low = ci_low,
        ci95_high = ci_high, stringsAsFactors = FALSE)
    }))
    paired_results$FDR <- p.adjust(paired_results$pValue, method = "BH")
    wilcox_p <- apply(ME_delta, 2, function(delta) {
      delta <- delta[is.finite(delta)]
      if (length(delta) < 3L || all(delta == 0))
        return(NA_real_)
      suppressWarnings(wilcox.test(delta, mu = 0, exact = FALSE)$p.value)
    })
    paired_results$wilcox_pValue <- wilcox_p[paired_results$module]
    paired_results$wilcox_FDR <- p.adjust(paired_results$wilcox_pValue, method = "BH")
    paired_results$direction <- ifelse(paired_results$mean_week24_minus_DX > 0, "higher_at_week24", "lower_at_week24")
    paired_results <- paired_results[order(paired_results$FDR, paired_results$pValue), ]
    write.csv(paired_results, file.path(outdir, "paired_module_time_association_all.csv"), row.names = FALSE)
    write.csv(subset(paired_results, FDR < 0.05), file.path(outdir, "paired_module_time_association_FDR05.csv"),
      row.names = FALSE)
    write.csv(data.frame(subject = rownames(ME_delta), ME_delta, check.names = FALSE), file.path(outdir, "module_eigengene_week24_minus_DX_by_subject.csv"),
      row.names = FALSE)
    non_grey_results <- paired_results[paired_results$moduleColor != "grey", ]
    top5_modules <- head(non_grey_results$module, 5L)
    write.csv(head(non_grey_results, 5L), file.path(outdir, "top5_modules_by_paired_test.csv"), row.names = FALSE)
    heat_order <- paired_results$module
    heat_effect <- matrix(paired_results$cohen_dz, ncol = 1, dimnames = list(heat_order, "week24 - DX"))
    heat_text <- matrix(paste0("dz=", sprintf("%.2f", paired_results$cohen_dz), "\nP=", format.pval(paired_results$pValue,
      digits = 2), "\nFDR=", format.pval(paired_results$FDR, digits = 2)), ncol = 1, dimnames = dimnames(heat_effect))
    heat_limit <- max(1, max(abs(heat_effect), na.rm = TRUE))
    pdf(file.path(outdir, "05_paired_module_time_association_heatmap.pdf"), width = 9, height = max(8, 0.7 *
      nrow(heat_effect) + 2.5))
    par(mar = c(8, 11, 5, 4), cex = figure_cex, cex.main = figure_cex_main, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis)
    labeledHeatmap(Matrix = heat_effect, xLabels = colnames(heat_effect), yLabels = rownames(heat_effect),
      ySymbols = rownames(heat_effect), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = heat_text,
      setStdMargins = FALSE, cex.text = annotation_text_cex, cex.lab = module_label_cex, cex.legendLabel = figure_cex_axis,
      zlim = c(-heat_limit, heat_limit), main = paste0("Paired module-time association (", length(complete_subjects),
        " pairs)"))
    dev.off()
    kME <- signedKME(datExpr, MEs, corFnc = "bicor", corOptions = paste0("use = 'p', maxPOutliers = ", max_p_outliers,
      ", pearsonFallback = '", pearson_fallback, "'"))
    own_kME <- vapply(seq_len(nGenes), function(i) {
      kme_name <- paste0("kME", moduleColors[i])
      if (kme_name %in% colnames(kME))
        kME[i, kme_name]
      else NA_real_
    }, numeric(1))
    geneInfo <- data.frame(gene = colnames(datExpr), moduleColor = moduleColors, moduleLabel = moduleLabels,
      kME_own_module = own_kME, stringsAsFactors = FALSE)
    geneInfo <- merge(geneInfo, paired_results[, c("moduleColor", "mean_week24_minus_DX", "cohen_dz", "pValue",
      "FDR", "wilcox_pValue", "wilcox_FDR", "direction")], by = "moduleColor", all.x = TRUE, sort = FALSE)
    geneInfo <- geneInfo[order(geneInfo$moduleLabel, -abs(geneInfo$kME_own_module)), ]
    write.csv(geneInfo, file.path(outdir, "all_gene_module_information.csv"), row.names = FALSE)
    top5_gene_info <- geneInfo[paste0("ME", geneInfo$moduleColor) %in% top5_modules, ]
    write.csv(top5_gene_info, file.path(outdir, "top5_paired_modules_combined_genes.csv"), row.names = FALSE)
    for (module in top5_modules) {
      color <- sub("^ME", "", module)
      genes <- geneInfo$gene[geneInfo$moduleColor == color]
      writeLines(genes, file.path(outdir, paste0("module_", module, "_genes.txt")))
    }
    target_genes <- c("FCGR1A", "FCGR1B", "ANKRD22", "DUSP3", "GBP5", "VAMP5", "PSTPIP2", "GBP1", "NPC2", "TIFA",
      "GBP2", "LHFPL2", "GYG1", "FLVCR2", "SECTM1", "POMP", "FBXO6", "GADD45B", "RHBDF2", "CASP4")
    target_gene_results <- geneInfo[geneInfo$gene %in% target_genes, c("gene", "moduleColor", "moduleLabel",
      "kME_own_module", "mean_week24_minus_DX", "cohen_dz", "pValue", "FDR", "wilcox_pValue", "wilcox_FDR",
      "direction")]
    target_gene_results$in_top5_paired_modules <- paste0("ME", target_gene_results$moduleColor) %in% top5_modules
    write.csv(target_gene_results, file.path(outdir, "target_gene_module_lookup.csv"), row.names = FALSE)
    saveRDS(list(parameters = list(network_type = network_type, cor_type = cor_type, max_p_outliers = max_p_outliers,
      pearson_fallback = pearson_fallback, power = power, min_module_size = min_module_size, merge_cut_height = merge_cut_height,
      mad_keep_fraction = mad_keep_fraction, scale_free_R2_cutoff = scale_free_R2_cutoff), sample_info = sample_info_paired,
      datExpr = datExpr, net = net, MEs = MEs, ME_delta = ME_delta, paired_results = paired_results, geneInfo = geneInfo),
      file.path(outdir, "paired_WGCNA_complete_results.rds"))
    writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
    cat("\nAnalysis completed:", format(Sys.time()), "\n")
    cat("Number of subjects:", length(complete_subjects), "\n")
    cat("Number of samples:", nSamples, "\n")
    cat("Number of genes:", nGenes, "\n")
    cat("Number of non-grey modules:", sum(unique(moduleColors) != "grey"), "\n")
    cat("FDR < 0.05 paired time-associated modules:", sum(paired_results$FDR < 0.05), "\n")
    cat("Top five exploratory modules:", paste(top5_modules, collapse = ", "), "\n")
    cat("Results directory:", outdir, "\n")
  }
  run_analysis()
}

acti_wgcna_figures <- function(cfg = acti_wgcna_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages(library(WGCNA))
  allowWGCNAThreads()
  project_dir <- cfg$output_dir
  input_dirs <- c(DX_week4 = cfg$week4_dir, DX_week24 = cfg$week24_dir)
  outdir <- file.path(project_dir, "03_figure3")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  figure_cex <- 1.44
  figure_cex_lab <- 1.58
  figure_cex_axis <- 1.42
  annotation_text_cex <- 1.16
  module_label_cex <- 1.42
  figure_line_width <- 2.6
  reference_line_width <- 2.4
  save_base_panel <- function(draw_fun, stem, width, height, dpi = 600) {
    grDevices::cairo_pdf(file.path(outdir, paste0(stem, ".pdf")), width = width, height = height, family = "sans")
    draw_fun()
    grDevices::dev.off()
    grDevices::svg(file.path(outdir, paste0(stem, ".svg")), width = width, height = height, family = "sans")
    draw_fun()
    grDevices::dev.off()
    grDevices::tiff(file.path(outdir, paste0(stem, ".tiff")), width = width, height = height, units = "in",
      res = dpi, compression = "lzw", family = "sans")
    draw_fun()
    grDevices::dev.off()
    grDevices::png(file.path(outdir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
    draw_fun()
    grDevices::dev.off()
  }
  draw_soft_threshold <- function(x, input_dir) {
    sft <- read.csv(file.path(input_dir, "soft_threshold_fit_indices.csv"), check.names = FALSE)
    power <- x$parameters$power
    r2_cutoff <- x$parameters$scale_free_R2_cutoff
    par(mfrow = c(1, 2), family = "sans", cex = figure_cex, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis,
      font = 2, font.lab = 2, font.axis = 2, lwd = figure_line_width, lend = "round", mar = c(5.5, 5.8,
        1.2, 1.5))
    plot(sft[, 1], -sign(sft[, 3]) * sft[, 2], xlab = "Soft-threshold power", ylab = expression("Signed scale-free topology fit (" *
      R^2 * ")"), type = "n", main = "")
    text(sft[, 1], -sign(sft[, 3]) * sft[, 2], labels = sft[, 1], col = "red", cex = annotation_text_cex, font = 2)
    box(lwd = figure_line_width)
    abline(h = r2_cutoff, col = "red", lty = 2, lwd = reference_line_width)
    abline(v = power, col = "blue", lty = 2, lwd = reference_line_width)
    plot(sft[, 1], sft[, 5], xlab = "Soft-threshold power", ylab = "Mean connectivity", type = "n", main = "")
    text(sft[, 1], sft[, 5], labels = sft[, 1], col = "red", cex = annotation_text_cex, font = 2)
    box(lwd = figure_line_width)
    abline(v = power, col = "blue", lty = 2, lwd = reference_line_width)
  }
  draw_gene_dendrogram <- function(x) {
    module_colors <- labels2colors(x$net$colors)
    par(family = "sans", cex = figure_cex, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis, font = 2,
      font.lab = 2, font.axis = 2, lwd = figure_line_width, lend = "round", mar = c(4, 5.6,
        1, 1.5))
    plotDendroAndColors(x$net$dendrograms[[1]], module_colors[x$net$blockGenes[[1]]], "Module colors", main = "",
      dendroLabels = FALSE, cex.colorLabels = module_label_cex, cex.dendroLabels = figure_cex_axis, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05, lwd = figure_line_width)
  }
  draw_eigengene_network <- function(x) {
    par(family = "sans", cex = figure_cex, cex.lab = figure_cex_lab, cex.axis = figure_cex_axis, font = 2,
      font.lab = 2, font.axis = 2, lwd = figure_line_width, lend = "round")
    plotEigengeneNetworks(x$MEs, "", marDendro = c(3.8, 4, 0.7, 4.5), marHeatmap = c(4.2,
      5, 0.7, 2.5), plotDendrograms = TRUE, xLabelsAngle = 90, cex.adjacency = annotation_text_cex,
      cex.preservation = annotation_text_cex, lwd = figure_line_width)
  }
  draw_module_association <- function(x, x_label) {
    paired <- x$paired_results
    heat_effect <- matrix(paired$cohen_dz, ncol = 1, dimnames = list(paired$module, x_label))
    heat_text <- matrix(paste0("dz=", sprintf("%.2f", paired$cohen_dz), "\nP=", format.pval(paired$pValue,
      digits = 2), "\nFDR=", format.pval(paired$FDR, digits = 2)), ncol = 1, dimnames = dimnames(heat_effect))
    heat_limit <- max(1, max(abs(heat_effect), na.rm = TRUE))
    par(family = "sans", mar = c(6.5, 11.5, 1, 4.5), cex = 1.3, cex.lab = 1.42, cex.axis = 1.26,
      font = 1, font.lab = 1, font.axis = 1, lwd = 1)
    labeledHeatmap(Matrix = heat_effect, xLabels = colnames(heat_effect), yLabels = rownames(heat_effect),
      ySymbols = rownames(heat_effect), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = heat_text,
      setStdMargins = FALSE, xLabelsAngle = 0, xLabelsAdj = 0.5, cex.text = 1.02, cex.lab = 1.28, font.lab.x = 1,
      font.lab.y = 1, cex.legendLabel = 1.26, verticalSeparator.lwd = 1, horizontalSeparator.lwd = 1, zlim = c(-heat_limit,
        heat_limit), main = "")
  }
  audit_rows <- list()
  for (comparison in names(input_dirs)) {
    input_dir <- input_dirs[[comparison]]
    rds_file <- file.path(input_dir, "paired_WGCNA_complete_results.rds")
    stopifnot(file.exists(rds_file))
    x <- readRDS(rds_file)
    stopifnot(identical(x$parameters$network_type, "signed"), identical(x$parameters$cor_type, "bicor"), isTRUE(all.equal(x$parameters$merge_cut_height,
      0.25)))
    late_label <- if (comparison == "DX_week4")
      "week4 - DX"
    else "week24 - DX"
    prefix <- if (comparison == "DX_week4")
      "DX_week4"
    else "DX_week24"
    save_base_panel(function() draw_soft_threshold(x, input_dir), paste0(prefix, "_02_soft_threshold_selection"),
      14, 7)
    save_base_panel(function() draw_gene_dendrogram(x), paste0(prefix, "_03_gene_dendrogram_and_modules"),
      16, 9)
    save_base_panel(function() draw_eigengene_network(x), paste0(prefix, "_04_eigengene_network"), 14, 12)
    association_height <- max(8, 0.7 * nrow(x$paired_results) + 2.5)
    save_base_panel(function() draw_module_association(x, late_label), paste0(prefix, "_05_paired_module_time_association_heatmap"),
      9, association_height)
    audit_rows[[comparison]] <- data.frame(comparison = gsub("_", "-", comparison), source_rds = normalizePath(rds_file,
      winslash = "/"), network_type = x$parameters$network_type, correlation = x$parameters$cor_type, merge_cut_height = x$parameters$merge_cut_height,
      soft_power = x$parameters$power, samples = nrow(x$datExpr), paired_patients = length(unique(x$sample_info$subject)),
      genes = ncol(x$datExpr), modules_including_grey = nrow(x$paired_results), stringsAsFactors = FALSE)
  }
  audit <- do.call(rbind, audit_rows)
  write.csv(audit, file.path(outdir, "Figure3_redraw_source_audit.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  message("Completed original-style Figure 3 redraw: ", outdir)
}

acti_wgcna_heatmaps <- function(cfg = acti_wgcna_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(dplyr)
    library(pheatmap)
    library(RColorBrewer)
    library(grid)
  })
  project_dir <- cfg$output_dir
  figure3_dir <- file.path(project_dir, "03_figure3")
  outdir <- file.path(figure3_dir, "top3_modules")
  source_dir <- file.path(outdir, "source_data")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  analyses <- list(DX_week4 = list(input_dir = cfg$week4_dir, late_timepoint = "week4", late_display = "week 4",
    group_name = "week_4", stem = "G_DX_week4_top3_modules_top5_genes_heatmap"), DX_week24 = list(input_dir = cfg$week24_dir,
    late_timepoint = "week24", late_display = "week 24", group_name = "week_24", stem = "H_DX_week24_top3_modules_top5_genes_heatmap"))
  is_valid_colour <- function(x) {
    vapply(x, function(cc) {
      tryCatch({
        grDevices::col2rgb(cc)
        TRUE
      }, error = function(e) FALSE)
    }, logical(1))
  }
  paired_gene_test <- function(expr_gene_by_sample, sample_info, late_timepoint) {
    subjects <- intersect(sample_info$subject[sample_info$timepoint == "DX"], sample_info$subject[sample_info$timepoint ==
      late_timepoint])
    subjects <- sort(unique(subjects))
    dx_samples <- vapply(subjects, function(id) {
      sample_info$sample[sample_info$subject == id & sample_info$timepoint == "DX"][1]
    }, character(1))
    late_samples <- vapply(subjects, function(id) {
      sample_info$sample[sample_info$subject == id & sample_info$timepoint == late_timepoint][1]
    }, character(1))
    stopifnot(all(dx_samples %in% colnames(expr_gene_by_sample)), all(late_samples %in% colnames(expr_gene_by_sample)),
      length(dx_samples) == length(late_samples))
    dx <- expr_gene_by_sample[, dx_samples, drop = FALSE]
    late <- expr_gene_by_sample[, late_samples, drop = FALSE]
    delta <- late - dx
    n <- rowSums(is.finite(delta))
    delta_mean <- rowMeans(delta, na.rm = TRUE)
    delta_sd <- apply(delta, 1, sd, na.rm = TRUE)
    t_stat <- delta_mean/(delta_sd/sqrt(n))
    pvalue <- 2 * stats::pt(-abs(t_stat), df = pmax(n - 1, 1))
    pvalue[n < 2 | !is.finite(t_stat)] <- NA_real_
    data.frame(gene = rownames(expr_gene_by_sample), n_pairs = n, logFC = delta_mean, t = t_stat, pvalue = pvalue,
      padj = p.adjust(pvalue, method = "BH"), stringsAsFactors = FALSE)
  }
  draw_saved_heatmap <- function(ph, file, device_fun, width, height, ...) {
    device_fun(file, width = width, height = height, ...)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    grDevices::dev.off()
  }
  make_heatmap <- function(cfg, comparison_name) {
    rds_file <- file.path(cfg$input_dir, "paired_WGCNA_complete_results.rds")
    top_module_file <- file.path(cfg$input_dir, "top5_modules_by_paired_test.csv")
    stopifnot(file.exists(rds_file), file.exists(top_module_file))
    x <- readRDS(rds_file)
    stopifnot(identical(x$parameters$network_type, "signed"), identical(x$parameters$cor_type, "bicor"), isTRUE(all.equal(x$parameters$merge_cut_height,
      0.25)))
    top_table <- read.csv(top_module_file, check.names = FALSE)
    top3_modules <- head(as.character(top_table$module), 3L)
    top3_colours <- sub("^ME", "", top3_modules)
    expr <- t(as.matrix(x$datExpr))
    gene_stats <- paired_gene_test(expr, x$sample_info, cfg$late_timepoint)
    module_map <- x$geneInfo[, c("gene", "moduleColor"), drop = FALSE]
    module_map <- module_map[!duplicated(module_map$gene), , drop = FALSE]
    ranked <- merge(gene_stats, module_map, by = "gene") %>% filter(moduleColor %in% top3_colours, !is.na(pvalue)) %>%
      mutate(module_rank = match(moduleColor, top3_colours)) %>% group_by(moduleColor) %>% arrange(pvalue,
      desc(abs(logFC)), gene, .by_group = TRUE) %>% slice_head(n = 5L) %>% ungroup() %>% arrange(module_rank,
      pvalue, desc(abs(logFC)), gene)
    if (nrow(ranked) != 15L) {
      stop(comparison_name, ": expected 15 genes but selected ", nrow(ranked), ".")
    }
    subjects <- sort(unique(x$sample_info$subject))
    dx_samples <- vapply(subjects, function(id) {
      x$sample_info$sample[x$sample_info$subject == id & x$sample_info$timepoint == "DX"][1]
    }, character(1))
    late_samples <- vapply(subjects, function(id) {
      x$sample_info$sample[x$sample_info$subject == id & x$sample_info$timepoint == cfg$late_timepoint][1]
    }, character(1))
    column_order <- c(dx_samples, late_samples)
    heat_expr <- expr[ranked$gene, column_order, drop = FALSE]
    heat_scaled <- t(scale(t(heat_expr)))
    keep <- apply(heat_scaled, 1, function(z) all(is.finite(z)))
    heat_scaled <- heat_scaled[keep, , drop = FALSE]
    ranked <- ranked[match(rownames(heat_scaled), ranked$gene), , drop = FALSE]
    if (nrow(heat_scaled) != 15L)
      stop("A selected gene failed row scaling.")
    annotation_col <- data.frame(Group = c(rep("DX", length(dx_samples)), rep(cfg$group_name, length(late_samples))),
      row.names = column_order, stringsAsFactors = FALSE)
    annotation_row <- data.frame(Module = factor(ranked$moduleColor, levels = top3_colours), row.names = ranked$gene,
      stringsAsFactors = FALSE)
    valid <- is_valid_colour(top3_colours)
    fallback <- RColorBrewer::brewer.pal(3, "Set2")
    module_colours <- ifelse(valid, top3_colours, fallback[seq_along(top3_colours)])
    names(module_colours) <- top3_colours
    group_colours <- c(DX = "#E64B35", week_4 = "#4DBBD5", week_24 = "#4DBBD5")
    ann_colours <- list(Group = group_colours[c("DX", cfg$group_name)], Module = module_colours)
    module_sizes <- table(ranked$moduleColor)
    module_sizes <- as.integer(module_sizes[top3_colours])
    gaps_row <- cumsum(module_sizes)
    gaps_row <- gaps_row[-length(gaps_row)]
    ph <- pheatmap(heat_scaled, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation_col,
      annotation_row = annotation_row, annotation_colors = ann_colours, gaps_row = gaps_row, show_colnames = FALSE,
      show_rownames = TRUE, fontsize = 14, fontsize_row = 14, fontsize_col = 13, fontface = "bold", border_color = NA,
      color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100), breaks = seq(-3, 3, length.out = 101),
      silent = TRUE)
    module_bar_columns <- unique(ph$gtable$layout$l[ph$gtable$layout$name == "row_annotation"])
    for (j in module_bar_columns) {
      ph$gtable$widths[j] <- ph$gtable$widths[j] + grid::unit(5, "mm")
    }
    ph$gtable <- gtable::gtable_add_cols(ph$gtable, grid::unit(24, "mm"), pos = -1)
    matrix_columns <- unique(ph$gtable$layout$l[ph$gtable$layout$name == "matrix"])
    for (j in matrix_columns) {
      ph$gtable$widths[j] <- ph$gtable$widths[j] - grid::unit(29, "mm")
    }
    width <- 13.5
    height <- 6.4
    draw_saved_heatmap(ph, file.path(outdir, paste0(cfg$stem, ".pdf")), grDevices::cairo_pdf, width, height,
      family = "sans")
    draw_saved_heatmap(ph, file.path(outdir, paste0(cfg$stem, ".svg")), grDevices::svg, width, height, family = "sans")
    draw_saved_heatmap(ph, file.path(outdir, paste0(cfg$stem, ".tiff")), grDevices::tiff, width, height, units = "in",
      res = 600, compression = "lzw", family = "sans")
    draw_saved_heatmap(ph, file.path(outdir, paste0(cfg$stem, ".png")), grDevices::png, width, height, units = "in",
      res = 300)
    write.csv(ranked, file.path(source_dir, paste0(cfg$stem, "_selected_genes.csv")), row.names = FALSE)
    write.csv(data.frame(gene = rownames(heat_scaled), heat_scaled, check.names = FALSE), file.path(source_dir,
      paste0(cfg$stem, "_row_zscores.csv")), row.names = FALSE)
    write.csv(data.frame(sample = rownames(annotation_col), annotation_col, check.names = FALSE), file.path(source_dir,
      paste0(cfg$stem, "_sample_annotation.csv")), row.names = FALSE)
    data.frame(comparison = comparison_name, source_rds = normalizePath(rds_file, winslash = "/"), paired_patients = length(subjects),
      top3_modules = paste(top3_modules, collapse = "; "), selected_genes = nrow(ranked), selection_rule = "Five smallest paired t-test P values per top-three module",
      stringsAsFactors = FALSE)
  }
  audit <- do.call(rbind, lapply(names(analyses), function(nm) {
    make_heatmap(analyses[[nm]], nm)
  }))
  write.csv(audit, file.path(outdir, "top3_module_heatmap_audit.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  message("Completed top-three module heatmaps: ", outdir)
}

acti_wgcna_stpmgs <- function(cfg = acti_wgcna_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(ggplot2)
  })
  project_dir <- cfg$output_dir
  week4_dir <- cfg$week4_dir
  week24_dir <- cfg$week24_dir
  outdir <- file.path(project_dir, "03_figure3", "top3_modules")
  source_dir <- file.path(outdir, "source_data")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  read_top3_set <- function(input_dir) {
    rds_file <- file.path(input_dir, "paired_WGCNA_complete_results.rds")
    top_file <- file.path(input_dir, "top5_modules_by_paired_test.csv")
    stopifnot(file.exists(rds_file), file.exists(top_file))
    x <- readRDS(rds_file)
    top3_modules <- head(as.character(read.csv(top_file)$module), 3L)
    top3_colours <- sub("^ME", "", top3_modules)
    module_genes <- lapply(top3_colours, function(colour) {
      sort(unique(as.character(x$geneInfo$gene[x$geneInfo$moduleColor == colour])))
    })
    names(module_genes) <- top3_modules
    list(source_rds = normalizePath(rds_file, winslash = "/"), modules = top3_modules, module_genes = module_genes,
      genes = sort(unique(unlist(module_genes, use.names = FALSE))), parameters = x$parameters)
  }
  w4 <- read_top3_set(week4_dir)
  w24 <- read_top3_set(week24_dir)
  stopifnot(identical(w4$parameters$network_type, "signed"), identical(w24$parameters$network_type, "signed"),
    isTRUE(all.equal(w4$parameters$merge_cut_height, 0.25)), isTRUE(all.equal(w24$parameters$merge_cut_height,
      0.25)))
  shared <- sort(intersect(w4$genes, w24$genes))
  w4_only <- sort(setdiff(w4$genes, w24$genes))
  w24_only <- sort(setdiff(w24$genes, w4$genes))
  all_genes <- sort(union(w4$genes, w24$genes))
  stopifnot(length(w4$genes) == 1315L, length(w24$genes) == 1738L, length(shared) == 1174L)
  circle_polygon <- function(x0, y0, radius, group) {
    angle <- seq(0, 2 * pi, length.out = 721L)
    data.frame(x = x0 + radius * cos(angle), y = y0 + radius * sin(angle), group = group)
  }
  circle_data <- rbind(circle_polygon(-0.72, 0, 1.18, "DX-week4"), circle_polygon(0.72,
    0, 1.18, "DX-week24"))
  p_venn <- ggplot(circle_data, aes(x = x, y = y, group = group, fill = group, colour = group)) + geom_polygon(alpha = 0.42,
    linewidth = 1.05) + scale_fill_manual(values = c(`DX-week4` = "#4C78A8", `DX-week24` = "#E45756")) + scale_colour_manual(values = c(`DX-week4` = "#2D5F8B",
    `DX-week24` = "#B33B3B")) + annotate("text", x = 0, y = 0.26, label = "STPMGs", size = 9.6,
    fontface = "bold", family = "sans") + annotate("text", x = 0, y = -0.22, label = paste0("n = ", format(length(shared),
    big.mark = ","), "\ngenes"), size = 6.1, lineheight = 0.95, fontface = "bold",
    family = "sans") + annotate("text", x = -1.22, y = 0.28, label = "DX-week 4", size = 10,
    fontface = "bold", family = "sans") + annotate("text", x = -1.22, y = -0.22, label = paste0("Total n = ",
    format(length(w4$genes), big.mark = ","), "\ngenes"), size = 5.7, lineheight = 0.95,
    fontface = "bold", family = "sans") + annotate("text", x = 1.22, y = 0.28, label = "DX-week 24",
    size = 10, fontface = "bold", family = "sans") + annotate("text", x = 1.22, y = -0.22, label = paste0("Total n = ",
    format(length(w24$genes), big.mark = ","), "\ngenes"), size = 5.7, lineheight = 0.95,
    fontface = "bold", family = "sans") + coord_equal(xlim = c(-2.18, 2.18), ylim = c(-1.3,
    1.3), clip = "off") + theme_void(base_family = "sans") + theme(legend.position = "none", plot.margin = margin(8,
    12, 8, 12, unit = "pt"))
  stem <- file.path(outdir, "I_top3_module_gene_overlap_STPMGs_venn")
  width <- 9
  height <- 6
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "sans")
  print(p_venn)
  grDevices::dev.off()
  grDevices::svg(paste0(stem, ".svg"), width = width, height = height, family = "sans")
  print(p_venn)
  grDevices::dev.off()
  grDevices::tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = 600, compression = "lzw",
    family = "sans")
  print(p_venn)
  grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300, bg = "white")
  print(p_venn)
  grDevices::dev.off()
  write.csv(data.frame(gene = shared, set_name = "STPMGs", full_name = "Shared treatment-progression module genes"),
    file.path(source_dir, "STPMGs_1174_shared_genes.csv"), row.names = FALSE)
  write.csv(data.frame(gene = w4$genes), file.path(source_dir, "DX_week4_top3_module_genes_1315.csv"), row.names = FALSE)
  write.csv(data.frame(gene = w24$genes), file.path(source_dir, "DX_week24_top3_module_genes_1738.csv"), row.names = FALSE)
  write.csv(data.frame(gene = all_genes, in_DX_week4_top3 = all_genes %in% w4$genes, in_DX_week24_top3 = all_genes %in%
    w24$genes, set = ifelse(all_genes %in% shared, "STPMGs", ifelse(all_genes %in% w4_only, "DX-week4 only",
    "DX-week24 only"))), file.path(source_dir, "top3_module_gene_overlap_membership.csv"), row.names = FALSE)
  module_audit <- rbind(data.frame(comparison = "DX-week4", module = names(w4$module_genes), n_genes = lengths(w4$module_genes),
    source_rds = w4$source_rds, stringsAsFactors = FALSE), data.frame(comparison = "DX-week24", module = names(w24$module_genes),
    n_genes = lengths(w24$module_genes), source_rds = w24$source_rds, stringsAsFactors = FALSE))
  write.csv(module_audit, file.path(source_dir, "top3_module_gene_counts_audit.csv"), row.names = FALSE)
  summary_table <- data.frame(set = c("DX-week4 top3", "STPMGs overlap", "DX-week24 top3", "DX-week4 only", "DX-week24 only",
    "Union"), n_genes = c(length(w4$genes), length(shared), length(w24$genes), length(w4_only), length(w24_only),
    length(all_genes)), stringsAsFactors = FALSE)
  write.csv(summary_table, file.path(source_dir, "top3_module_gene_overlap_summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "I_STPMGs_sessionInfo.txt"))
  message("Completed STPMGs Venn diagram: ", paste0(stem, ".pdf"))
}

acti_wgcna_validate <- function(cfg = acti_wgcna_config()) {
  expected <- list(week4 = list(n = 67L, power = 6, genes = 1315L), week24 = list(n = 70L, power = 7, genes = 1738L))
  sets <- list()
  audit <- list()
  for (stage in names(expected)) {
    input_dir <- cfg[[paste0(stage, "_dir")]]
    x <- readRDS(file.path(input_dir, "paired_WGCNA_complete_results.rds"))
    e <- expected[[stage]]
    p <- x$parameters
    stopifnot(p$network_type == "signed", p$cor_type == "bicor", p$pearson_fallback == "individual", p$power ==
      e$power, p$merge_cut_height == 0.25, p$min_module_size == 30, p$max_p_outliers == 0.05,
      p$mad_keep_fraction == 0.5, p$scale_free_R2_cutoff == 0.85, nrow(x$datExpr) == 2L *
        e$n, ncol(x$datExpr) == 8069L, nrow(x$ME_delta) == e$n)
    dx <- x$sample_info[x$sample_info$timepoint == "DX", , drop = FALSE]
    late <- x$sample_info[x$sample_info$timepoint == stage, , drop = FALSE]
    late <- late[match(dx$subject, late$subject), , drop = FALSE]
    stopifnot(nrow(dx) == e$n, !anyNA(late$subject), !anyDuplicated(dx$subject))
    sample_field <- intersect(c("sample", "sample_id", "sample_name"), names(dx))
    stopifnot(length(sample_field) >= 1L)
    sample_field <- sample_field[[1L]]
    delta <- as.matrix(x$MEs[as.character(late[[sample_field]]), , drop = FALSE]) - as.matrix(x$MEs[as.character(dx[[sample_field]]),
      , drop = FALSE])
    rownames(delta) <- as.character(dx$subject)
    saved_delta <- as.matrix(x$ME_delta[match(rownames(delta), rownames(x$ME_delta)), colnames(delta), drop = FALSE])
    stopifnot(isTRUE(all.equal(unname(delta), unname(saved_delta), tolerance = 1e-10)))
    pv <- apply(delta, 2L, function(d) t.test(d)$p.value)
    dz <- colMeans(delta)/apply(delta, 2L, sd)
    pr <- x$paired_results[match(names(pv), x$paired_results$module), , drop = FALSE]
    stopifnot(isTRUE(all.equal(unname(pv), pr$pValue, tolerance = 1e-10)), isTRUE(all.equal(unname(dz), pr$cohen_dz,
      tolerance = 1e-10)), isTRUE(all.equal(unname(p.adjust(pv, "BH")), pr$FDR, tolerance = 1e-10)))
    top_file <- read.csv(file.path(input_dir, "top5_modules_by_paired_test.csv"))
    eligible <- x$paired_results[x$paired_results$module != "MEgrey", , drop = FALSE]
    selected <- head(as.character(top_file$module), 3L)
    stopifnot(setequal(selected, head(eligible$module[order(eligible$pValue)], 3L)))
    sets[[stage]] <- sort(unique(as.character(x$geneInfo$gene[x$geneInfo$moduleColor %in% sub("^ME", "", selected)])))
    stopifnot(length(sets[[stage]]) == e$genes)
    audit[[stage]] <- data.frame(comparison = stage, pairs = e$n, genes = ncol(x$datExpr), power = p$power,
      selected_modules = paste(selected, collapse = ";"), selected_module_genes = length(sets[[stage]]))
  }
  shared <- sort(intersect(sets$week4, sets$week24))
  stopifnot(length(shared) == 1174L)
  result <- do.call(rbind, audit)
  print(result, row.names = FALSE)
  message("Validated saved statistics and top-three module membership: 1174 STPMGs.")
  invisible(list(audit = result, STPMGs = shared, gene_sets = sets))
}

acti_wgcna_run <- function(stage, cfg = acti_wgcna_config()) {
  stages <- c("week4", "week24", "validate", "figures", "heatmaps", "stpmgs")
  if (length(stage) != 1L || !(stage %in% c("all", stages))) {
    stop("Choose stage: all, ", paste(stages, collapse = ", "))
  }
  selected <- if (stage == "all")
    stages
  else stage
  for (s in selected) {
    message("Running ACTI WGCNA stage: ", s)
    get(paste0("acti_wgcna_", s), mode = "function")(cfg)
  }
  invisible(cfg)
}

acti_wgcna_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- "Rscript 2_WGCNA_analysis.R --stage NAME\nStages: all, week4, week24, validate, figures, heatmaps, stpmgs\nUse --check to check inputs; see the configuration function for paths.\n"
  if (!length(args) || "--help" %in% args) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  values <- list()
  stage <- NULL
  check <- FALSE
  i <- 1L
  value_flags <- c("--data-dir", "--output-dir", "--expr-file", "--pheno-file", "--week4-dir", "--week24-dir",
    "--stage")
  while (i <= length(args)) {
    flag <- args[[i]]
    if (flag == "--overwrite") {
      values$overwrite <- TRUE
    }
    else if (flag == "--check") {
      check <- TRUE
    }
    else if (flag %in% value_flags) {
      if (i == length(args) || startsWith(args[[i + 1L]], "--"))
        stop("Missing value: ", flag)
      i <- i + 1L
      if (flag == "--stage")
        stage <- args[[i]]
      else {
        key <- gsub("-", "_", substring(flag, 3L), fixed = TRUE)
        values[[key]] <- args[[i]]
      }
    }
    else stop("Unknown argument: ", flag)
    i <- i + 1L
  }
  cfg <- do.call(acti_wgcna_config, values)
  if (check) {
    ok <- acti_wgcna_check(cfg)
    if (!ok)
      stop("Preflight incomplete: supply the missing inputs/packages.")
    return(invisible(cfg))
  }
  if (is.null(stage))
    stop("Supply --stage or --check.\n", help)
  acti_wgcna_run(stage, cfg)
}

if (sys.nframe() == 0L) acti_wgcna_main()

