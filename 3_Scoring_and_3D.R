# Candidate-gene scoring and 3D plots.
# Inputs: Single-cell/WGCNA tables, GSE89403 expression and imaging data; acti_fig4_config().
# Run: Rscript 3_scoring_and_3D.R --stage score
# source() loads functions without starting an analysis.

acti_fig4_config <- function(data_dir = file.path("data", "GSE89403"), output_dir = file.path("results", "Figure4"),
  sc_results_dir = file.path("results", "04_celltype_pseudobulk"), wgcna_dir = file.path("results", "WGCNA"),
  week4_dir = NULL, week24_dir = NULL, score_dir = NULL, expression_file = NULL, phenotype_file = NULL, pet_file = NULL,
  suv_file = NULL, sc_all_tested_file = NULL, sc_candidates_file = NULL, stpmg_file = NULL, overwrite = FALSE) {
  default <- function(value, fallback) if (is.null(value))
    fallback
  else value
  list(expression_file = default(expression_file, file.path(data_dir, "GSE89403_log2ExpGeneNames.csv")), phenotype_file = default(phenotype_file,
    file.path(data_dir, "GSE89403_matrix_pheno.csv")), pet_file = default(pet_file, file.path(data_dir, "pet_data.csv")),
    suv_file = default(suv_file, file.path(data_dir, "SUVmax_data.csv")), sc_all_tested_file = default(sc_all_tested_file,
      file.path(sc_results_dir, "tables", "05_celltype_HIGH_vs_LOW_all_tested_genes.csv")), sc_candidates_file = default(sc_candidates_file,
      file.path(sc_results_dir, "tables", "06_celltype_WIDEST_candidates.csv")), week4_dir = default(week4_dir,
      file.path(wgcna_dir, "01_DX_week4")), week24_dir = default(week24_dir, file.path(wgcna_dir, "02_DX_week24")),
    stpmg_file = default(stpmg_file, file.path(wgcna_dir, "03_figure3", "top3_modules", "source_data", "STPMGs_1174_shared_genes.csv")),
    score_dir = default(score_dir, file.path(output_dir, "01_scores")), figure_dir = file.path(output_dir,
      "02_figures"), overwrite = isTRUE(overwrite))
}

acti_fig4_stages <- function(stage) {
  if (length(stage) != 1L || !(stage %in% c("all", "score", "validate", "plot3d", "heatmap", "overlap"))) {
    stop("Choose stage: all, score, validate, plot3d, heatmap, overlap")
  }
  if (stage == "all")
    c("score", "validate", "plot3d", "heatmap")
  else stage
}

acti_fig4_check <- function(cfg = acti_fig4_config(), stage = "all") {
  stages <- acti_fig4_stages(stage)
  packages <- list(score = "data.table", validate = character(), plot3d = c("dplyr", "plotly", "viridisLite",
    "htmlwidgets", "scales"), heatmap = c("ggplot2", "dplyr", "tidyr", "scales", "svglite", "ragg"), overlap = c("ggplot2",
    "svglite", "ragg"))
  pkgs <- unique(unlist(packages[stages], use.names = FALSE))
  available <- vapply(pkgs, function(p) nzchar(system.file(package = p)), logical(1))
  if (length(pkgs))
    print(data.frame(package = pkgs, available = available), row.names = FALSE)
  files <- character()
  if ("score" %in% stages)
    files <- c(files, cfg$expression_file, cfg$phenotype_file, cfg$pet_file, cfg$suv_file, cfg$sc_all_tested_file,
      file.path(cfg$week4_dir, c("top5_modules_by_paired_test.csv", "all_gene_module_information.csv")),
      file.path(cfg$week24_dir, c("top5_modules_by_paired_test.csv", "all_gene_module_information.csv")))
  if (!("score" %in% stages) && any(c("validate", "plot3d", "heatmap") %in% stages)) {
    files <- c(files, file.path(cfg$score_dir, c("04_final_signed_Top3_scRNA_candidates_ranked.csv", "05_Top20_direction_unrestricted_equal_weight.csv")))
  }
  if ("overlap" %in% stages)
    files <- c(files, cfg$sc_candidates_file, cfg$stpmg_file)
  print(data.frame(path = files, exists = file.exists(files)), row.names = FALSE)
  invisible(all(available) && all(file.exists(files)))
}

acti_fig4_score <- function(cfg = acti_fig4_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages(library(data.table))
  base_wgcna <- NULL
  dx_w4_modules_file <- file.path(cfg$week4_dir, "top5_modules_by_paired_test.csv")
  dx_w4_genes_file <- file.path(cfg$week4_dir, "all_gene_module_information.csv")
  dx_w24_modules_file <- file.path(cfg$week24_dir, "top5_modules_by_paired_test.csv")
  dx_w24_genes_file <- file.path(cfg$week24_dir, "all_gene_module_information.csv")
  expression_file <- cfg$expression_file
  phenotype_file <- cfg$phenotype_file
  pet_file <- cfg$pet_file
  suv_file <- cfg$suv_file
  sc_all_tested_file <- cfg$sc_all_tested_file
  output_dir <- cfg$score_dir
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = FALSE)) > 0L && !cfg$overwrite) {
    stop("Scoring output directory is not empty: ", output_dir, ". Choose a new directory or explicitly use --overwrite.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  required_files <- c(dx_w4_modules_file, dx_w4_genes_file, dx_w24_modules_file, dx_w24_genes_file, expression_file,
    phenotype_file, pet_file, suv_file, sc_all_tested_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
  }
  EXPECTED <- list(TGRV = c(n_observations = 214L, n_patients = 74L), SUVmax = c(n_observations = 140L, n_patients = 70L),
    Time = c(n_observations = 214L, n_patients = 74L))
  normalize_sample_code <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("day_?7$", "day_7", x, ignore.case = TRUE)
    x <- gsub("week_?4$", "week_4", x, ignore.case = TRUE)
    x <- gsub("week_?24$", "week_24", x, ignore.case = TRUE)
    x
  }
  assert_unique <- function(x, label) {
    duplicated_values <- unique(x[duplicated(x)])
    if (length(duplicated_values) > 0L) {
      stop(label, " contains duplicated identifiers: ", paste(head(duplicated_values, 20L), collapse = ", "))
    }
  }
  assert_no_metadata_conflict <- function(pheno, columns) {
    for (column in columns) {
      conflict <- pheno[, {
        values <- unique(trimws(as.character(get(column))))
        values <- values[!is.na(values) & nzchar(values) & values != "NA"]
        .(n_unique = length(values))
      }, by = sample_code][n_unique > 1L]
      if (nrow(conflict) > 0L) {
        stop("Conflicting phenotype values for column '", column, "' in sample_code(s): ", paste(head(conflict$sample_code,
         20L), collapse = ", "))
      }
    }
  }
  rm_rank_correlation <- function(x, y, patient) {
    ok <- is.finite(x) & is.finite(y) & !is.na(patient)
    x <- as.numeric(x[ok])
    y <- as.numeric(y[ok])
    patient <- as.character(patient[ok])
    patient_n <- table(patient)
    keep <- patient_n[patient] >= 2L
    x <- x[keep]
    y <- y[keep]
    patient <- patient[keep]
    n_obs <- length(x)
    n_patients <- length(unique(patient))
    df <- n_obs - n_patients - 1L
    if (n_obs < 4L || n_patients < 2L || df < 1L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
      return(c(rho = NA_real_, p = NA_real_, n_observations = n_obs, n_patients = n_patients))
    }
    x_rank <- rank(x, ties.method = "average")
    y_rank <- rank(y, ties.method = "average")
    x_within <- x_rank - ave(x_rank, patient, FUN = mean)
    y_within <- y_rank - ave(y_rank, patient, FUN = mean)
    denominator <- sqrt(sum(x_within^2) * sum(y_within^2))
    if (!is.finite(denominator) || denominator <= 0) {
      return(c(rho = NA_real_, p = NA_real_, n_observations = n_obs, n_patients = n_patients))
    }
    rho <- sum(x_within * y_within)/denominator
    rho <- max(min(rho, 1), -1)
    t_stat <- rho * sqrt(df/max(1 - rho^2, .Machine$double.eps))
    p_value <- 2 * pt(-abs(t_stat), df = df)
    c(rho = rho, p = p_value, n_observations = n_obs, n_patients = n_patients)
  }
  calculate_endpoint <- function(expression_matrix, endpoint_data, value_column) {
    sample_ids <- endpoint_data$sample_code
    if (!all(sample_ids %in% colnames(expression_matrix))) {
      stop("Some endpoint samples are absent from the expression matrix")
    }
    endpoint_data <- endpoint_data[match(sample_ids, endpoint_data$sample_code), , drop = FALSE]
    stopifnot(identical(sample_ids, endpoint_data$sample_code))
    result_matrix <- t(vapply(rownames(expression_matrix), function(gene) {
      rm_rank_correlation(x = as.numeric(expression_matrix[gene, sample_ids]), y = endpoint_data[[value_column]],
        patient = endpoint_data$subject)
    }, FUN.VALUE = c(rho = 0, p = 0, n_observations = 0, n_patients = 0)))
    data.table(Gene = rownames(expression_matrix), rho = result_matrix[, "rho"], p = result_matrix[, "p"],
      n_observations = as.integer(result_matrix[, "n_observations"]), n_patients = as.integer(result_matrix[,
        "n_patients"]))
  }
  add_endpoint_suffix <- function(result, endpoint_name) {
    old_names <- setdiff(names(result), "Gene")
    setnames(result, old_names, paste0(old_names, "_", endpoint_name))
    result
  }
  percentile_evidence <- function(x) {
    output <- rep(NA_real_, length(x))
    ok <- is.finite(x)
    n_ok <- sum(ok)
    if (n_ok == 1L) {
      output[ok] <- 1
    }
    else if (n_ok > 1L) {
      output[ok] <- (rank(x[ok], ties.method = "average") - 1)/(n_ok - 1)
    }
    output
  }
  check_expected_count <- function(data, endpoint_name) {
    actual <- c(n_observations = nrow(data), n_patients = uniqueN(data$subject))
    expected <- EXPECTED[[endpoint_name]]
    if (!identical(as.integer(actual), as.integer(expected))) {
      stop(endpoint_name, " sample audit failed. Expected ", expected["n_observations"], " observations / ",
        expected["n_patients"], " patients, but found ", actual["n_observations"], " observations / ",
        actual["n_patients"], " patients")
    }
  }
  dx_w4_modules <- fread(dx_w4_modules_file, check.names = FALSE)
  dx_w24_modules <- fread(dx_w24_modules_file, check.names = FALSE)
  dx_w4_all <- fread(dx_w4_genes_file, check.names = FALSE)
  dx_w24_all <- fread(dx_w24_genes_file, check.names = FALSE)
  for (x_name in c("dx_w4_modules", "dx_w24_modules")) {
    x <- get(x_name)
    if (!all(c("module", "moduleColor") %in% names(x))) {
      stop(x_name, " lacks module/moduleColor columns")
    }
    if (nrow(x) < 3L)
      stop(x_name, " contains fewer than three modules")
  }
  for (x_name in c("dx_w4_all", "dx_w24_all")) {
    x <- get(x_name)
    if (!all(c("gene", "moduleColor") %in% names(x))) {
      stop(x_name, " lacks gene/moduleColor columns")
    }
    x[, `:=`(gene, trimws(as.character(gene)))]
    x <- x[!is.na(gene) & nzchar(gene)]
    assign(x_name, x)
  }
  dx_w4_top3_colors <- dx_w4_modules$moduleColor[1:3]
  dx_w24_top3_colors <- dx_w24_modules$moduleColor[1:3]
  dx_w4 <- dx_w4_all[moduleColor %chin% dx_w4_top3_colors]
  dx_w24 <- dx_w24_all[moduleColor %chin% dx_w24_top3_colors]
  wgcna_intersection_genes <- sort(intersect(dx_w4$gene, dx_w24$gene))
  if (length(wgcna_intersection_genes) == 0L) {
    stop("Signed-WGCNA Top3 intersection is empty")
  }
  sc_for_candidates <- fread(sc_all_tested_file, check.names = FALSE)
  required_sc_candidate_columns <- c("gene", "celltype_major", "logFC", "PValue", "direction_consistency", "direction")
  if (!all(required_sc_candidate_columns %in% names(sc_for_candidates))) {
    stop("Single-cell all-tested table lacks required columns")
  }
  sc_selected <- sc_for_candidates[PValue < 0.05 & abs(logFC) >= 0.2 & direction_consistency >=
    0.8]
  sc_summary <- sc_selected[, .(sc_n_qualifying_celltypes = uniqueN(celltype_major), sc_qualifying_celltypes = paste(sort(unique(celltype_major)),
    collapse = ";"), sc_directions = paste(sort(unique(direction)), collapse = ";"), sc_best_nominal_P = min(PValue,
    na.rm = TRUE), sc_max_abs_log2FC = max(abs(logFC), na.rm = TRUE), sc_max_direction_consistency = max(direction_consistency,
    na.rm = TRUE)), by = gene]
  intersection_genes <- sort(intersect(wgcna_intersection_genes, sc_summary$gene))
  if (length(intersection_genes) == 0L) {
    stop("WGCNA Top3 and scRNA candidate intersection is empty")
  }
  dx_w4_info <- copy(dx_w4)
  setnames(dx_w4_info, setdiff(names(dx_w4_info), "gene"), paste0("DX_week4_", setdiff(names(dx_w4_info), "gene")))
  dx_w24_info <- copy(dx_w24)
  setnames(dx_w24_info, setdiff(names(dx_w24_info), "gene"), paste0("DX_week24_", setdiff(names(dx_w24_info),
    "gene")))
  intersection_detail <- merge(dx_w4_info[gene %chin% intersection_genes], dx_w24_info[gene %chin% intersection_genes],
    by = "gene", all = FALSE, sort = TRUE)
  stopifnot(nrow(intersection_detail) == length(intersection_genes))
  intersection_detail <- merge(intersection_detail, sc_summary[gene %chin% intersection_genes], by = "gene",
    all.x = TRUE, sort = TRUE)
  fwrite(intersection_detail, file.path(output_dir, "01_signed_Top3_scRNA_intersection_gene_details.csv"))
  writeLines(intersection_genes, file.path(output_dir, "02_signed_Top3_scRNA_intersection_gene_symbols.txt"))
  expression_table <- fread(expression_file, check.names = FALSE, showProgress = FALSE)
  setnames(expression_table, 1:2, c("ensembl_id", "symbol"))
  expression_table[, `:=`(ensembl_id, trimws(as.character(ensembl_id)))]
  expression_table[, `:=`(symbol, trimws(as.character(symbol)))]
  expression_table[, `:=`(analysis_gene_id, fifelse(!is.na(symbol) & nzchar(symbol) & symbol != "NA", symbol,
    ensembl_id))]
  expression_table <- expression_table[!is.na(analysis_gene_id) & nzchar(analysis_gene_id) & analysis_gene_id %chin%
    intersection_genes]
  sample_columns <- setdiff(names(expression_table), c("ensembl_id", "symbol", "analysis_gene_id"))
  expression_numeric <- as.matrix(expression_table[, ..sample_columns])
  storage.mode(expression_numeric) <- "double"
  expression_sum <- rowsum(expression_numeric, group = expression_table$analysis_gene_id, reorder = FALSE)
  symbol_count <- table(expression_table$analysis_gene_id)
  expression_matrix <- expression_sum/as.numeric(symbol_count[rownames(expression_sum)])
  missing_expression_genes <- setdiff(intersection_genes, rownames(expression_matrix))
  if (length(missing_expression_genes) > 0L) {
    stop("WGCNA-intersection genes absent from expression matrix: ", paste(missing_expression_genes, collapse = ", "))
  }
  expression_matrix <- expression_matrix[intersection_genes, , drop = FALSE]
  assert_unique(colnames(expression_matrix), "Expression sample columns")
  phenotype_raw <- fread(phenotype_file, check.names = FALSE)
  required_pheno_columns <- c("sample_code", "subject", "disease.state", "treatmentresult", "time", "timetonegativity")
  missing_pheno <- setdiff(required_pheno_columns, names(phenotype_raw))
  if (length(missing_pheno) > 0L) {
    stop("Phenotype columns missing: ", paste(missing_pheno, collapse = ", "))
  }
  phenotype_raw[, `:=`(sample_code, normalize_sample_code(sample_code))]
  assert_no_metadata_conflict(phenotype_raw, c("subject", "disease.state", "treatmentresult", "time", "timetonegativity"))
  phenotype <- unique(phenotype_raw[, ..required_pheno_columns], by = "sample_code")
  assert_unique(phenotype$sample_code, "De-duplicated phenotype sample_code")
  phenotype[, `:=`(subject, trimws(as.character(subject)))]
  phenotype[, `:=`(time, trimws(as.character(time)))]
  phenotype[, `:=`(treatmentresult, trimws(as.character(treatmentresult)))]
  phenotype[, `:=`(disease.state, trimws(as.character(disease.state)))]
  phenotype_analysis <- phenotype[toupper(treatmentresult) == "DEFINITE CURE" & toupper(disease.state) == "TB SUBJECTS" &
    time %chin% c("DX", "week_4", "week_24") & sample_code %chin% colnames(expression_matrix)]
  pet <- fread(pet_file, check.names = FALSE)
  suv <- fread(suv_file, check.names = FALSE)
  pet[, `:=`(sample_code, normalize_sample_code(sample_code))]
  suv[, `:=`(sample_code, normalize_sample_code(sample_code))]
  assert_unique(pet$sample_code, "PET sample_code")
  assert_unique(suv$sample_code, "SUVmax sample_code")
  analysis_metadata <- phenotype_analysis[, .(sample_code, subject, disease.state, treatmentresult, time, timetonegativity)]
  pet_analysis <- merge(analysis_metadata, pet[, .(sample_code, TGRV_log10)], by = "sample_code", all = FALSE,
    sort = FALSE)
  pet_analysis <- pet_analysis[is.finite(TGRV_log10)]
  pet_analysis <- pet_analysis[, if (.N >= 2L)
    .SD, by = subject]
  setorder(pet_analysis, subject, time)
  suv_analysis <- merge(analysis_metadata, suv[, .(sample_code, SUVmax)], by = "sample_code", all = FALSE, sort = FALSE)
  suv_analysis <- suv_analysis[is.finite(SUVmax)]
  suv_analysis <- suv_analysis[, if (.N >= 2L)
    .SD, by = subject]
  setorder(suv_analysis, subject, time)
  time_analysis <- copy(analysis_metadata)
  time_analysis[, `:=`(time_week, fcase(time == "DX", 0, time == "week_4", 4, time == "week_24", 24, default = NA_real_))]
  time_analysis <- time_analysis[is.finite(time_week)]
  time_analysis <- time_analysis[, if (.N >= 2L)
    .SD, by = subject]
  setorder(time_analysis, subject, time_week)
  assert_unique(pet_analysis$sample_code, "Eligible PET sample_code")
  assert_unique(suv_analysis$sample_code, "Eligible SUVmax sample_code")
  assert_unique(time_analysis$sample_code, "Eligible time sample_code")
  check_expected_count(pet_analysis, "TGRV")
  check_expected_count(suv_analysis, "SUVmax")
  check_expected_count(time_analysis, "Time")
  sample_audit <- data.table(Endpoint = c("TGRV", "SUVmax", "Treatment_stage"), Observations = c(nrow(pet_analysis),
    nrow(suv_analysis), nrow(time_analysis)), Patients = c(uniqueN(pet_analysis$subject), uniqueN(suv_analysis$subject),
    uniqueN(time_analysis$subject)))
  fwrite(sample_audit, file.path(output_dir, "03_sample_audit_summary.csv"))
  message("Calculating TGRV correlations for ", nrow(expression_matrix), " genes")
  result_tgrv <- add_endpoint_suffix(calculate_endpoint(expression_matrix, pet_analysis, "TGRV_log10"), "TGRV")
  message("Calculating SUVmax correlations")
  result_suv <- add_endpoint_suffix(calculate_endpoint(expression_matrix, suv_analysis, "SUVmax"), "SUVmax")
  message("Calculating treatment-stage correlations")
  result_time <- add_endpoint_suffix(calculate_endpoint(expression_matrix, time_analysis, "time_week"), "Time")
  result <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE, sort = FALSE), list(result_tgrv, result_suv,
    result_time))
  result[, `:=`(candidate_order, match(Gene, intersection_genes))]
  setorder(result, candidate_order)
  for (endpoint in c("TGRV", "SUVmax", "Time")) {
    p_column <- paste0("p_", endpoint)
    result[, `:=`((paste0("FDR_", endpoint)), p.adjust(get(p_column), method = "BH"))]
    result[, `:=`((paste0("evidence_", endpoint)), -log10(pmax(get(p_column), .Machine$double.xmin)))]
    result[, `:=`((paste0("percentile_", endpoint)), percentile_evidence(get(paste0("evidence_", endpoint))))]
  }
  result[, `:=`(score_raw_logp_sum, evidence_TGRV + evidence_SUVmax + evidence_Time)]
  result[, `:=`(score_equal_weight, rowMeans(.SD)), .SDcols = c("percentile_TGRV", "percentile_SUVmax", "percentile_Time")]
  result[, `:=`(Rank_final_candidates, as.integer(rank(-score_equal_weight, ties.method = "min", na.last = "keep")))]
  setnames(intersection_detail, "gene", "Gene")
  result <- merge(result, intersection_detail, by = "Gene", all.x = TRUE, sort = FALSE)
  result[, `:=`(candidate_order, match(Gene, intersection_genes))]
  setorder(result, Rank_final_candidates, Gene)
  fwrite(result, file.path(output_dir, "04_final_signed_Top3_scRNA_candidates_ranked.csv"))
  fwrite(result[order(Rank_final_candidates, Gene)][1:min(.N, 20L)], file.path(output_dir, "05_Top20_direction_unrestricted_equal_weight.csv"))
  overlap_020 <- copy(result)
  setorder(overlap_020, Rank_final_candidates, Gene)
  fwrite(overlap_020, file.path(output_dir, "06_primary_final_candidates_log2FC0.20.csv"))
  sc_050_genes <- unique(sc_for_candidates[PValue < 0.05 & abs(logFC) >= 0.5 & direction_consistency >=
    0.8, gene])
  overlap_050 <- result[Gene %chin% sc_050_genes]
  overlap_050[, `:=`(Rank_within_log2FC0.50_subset, as.integer(rank(-score_equal_weight, ties.method = "min",
    na.last = "keep")))]
  setorder(overlap_050, Rank_final_candidates, Gene)
  fwrite(overlap_050, file.path(output_dir, "06_sensitivity_candidates_log2FC0.50.csv"))
  panel20 <- c("FCGR1A", "FCGR1B", "ANKRD22", "DUSP3", "GBP5", "VAMP5", "PSTPIP2", "GBP1", "NPC2", "TIFA", "GBP2",
    "LHFPL2", "GYG1", "FLVCR2", "SECTM1", "POMP", "FBXO6", "GADD45B", "RHBDF2", "CASP4")
  panel61 <- c("TUBGCP6", "MYADM", "NR4A1", "KLF2", "AHNAK", "TNFAIP3", "FOS", "NFKBIA", "BTG2", "ACTB", "LMNA",
    "TRMT2A", "EIF5", "SDR39U1", "GBP1", "HSP90AA1", "GLUD1", "SOCS1", "MCL1", "TAGLN2", "GADD45B", "CD83",
    "LHFPL2", "ANKRD22", "RHBDF2", "FCGR1A", "HLA-F", "FCGR1B", "GBP5", "CRIP1", "FBXO6", "GZMB", "FGFBP2",
    "CASP4", "PSMB9", "NPC2", "DUSP3", "MMP9", "SLC18A2", "CALM1", "SECTM1", "FLVCR2", "VAMP5", "PSTPIP2",
    "GBP4", "GYG1", "IRF1", "TIFA", "ETV7", "STAT1", "POMP", "C1QC", "S100A10", "PFN1", "SERPING1", "IL2RG",
    "TAP1", "APOL6", "GBP2", "APOL1", "GSTP1")
  audit_panel <- function(panel, panel_name) {
    panel_table <- data.table(Gene = panel, panel_input_order = seq_along(panel))
    panel_table[, `:=`(in_DX_week4_signed_Top3, Gene %chin% dx_w4$gene)]
    panel_table[, `:=`(in_DX_week24_signed_Top3, Gene %chin% dx_w24$gene)]
    panel_table[, `:=`(in_both_signed_Top3, Gene %chin% wgcna_intersection_genes)]
    panel_table[, `:=`(in_single_cell_4155, Gene %chin% sc_summary$gene)]
    audit <- merge(panel_table, overlap_020, by = "Gene", all.x = TRUE, sort = FALSE)
    audit[, `:=`(in_final_candidate_set, !is.na(Rank_final_candidates))]
    setorder(audit, panel_input_order)
    fwrite(audit, file.path(output_dir, paste0("08_", panel_name, "_rank_in_final_set.csv")))
    audit
  }
  panel20_audit <- audit_panel(panel20, "panel20")
  panel61_audit <- audit_panel(panel61, "panel61")
  summary_table <- data.table(item = c("DX_week4_signed_top3_unique_genes", "DX_week24_signed_top3_unique_genes",
    "signed_top3_WGCNA_intersection_unique_genes", "single_cell_log2FC0.20_unique_candidates", "final_WGCNA_scRNA_intersection_candidates",
    "final_candidates_meeting_log2FC0.50_sensitivity", "panel20_genes_in_final_candidate_set", "panel61_genes_in_final_candidate_set"),
    value = c(uniqueN(dx_w4$gene), uniqueN(dx_w24$gene), length(wgcna_intersection_genes), uniqueN(sc_summary$gene),
      uniqueN(overlap_020$Gene), uniqueN(overlap_050$Gene), sum(panel20_audit$in_final_candidate_set), sum(panel61_audit$in_final_candidate_set)))
  fwrite(summary_table, file.path(output_dir, "00_analysis_summary.csv"))
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  message("Completed: ", output_dir)
  print(summary_table)
}

acti_fig4_plot3d <- function(cfg = acti_fig4_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(dplyr)
    library(plotly)
    library(viridisLite)
    library(htmlwidgets)
    library(scales)
  })
  analysis_dir <- cfg$score_dir
  all_file <- file.path(analysis_dir, "04_final_signed_Top3_scRNA_candidates_ranked.csv")
  top20_file <- file.path(analysis_dir, "05_Top20_direction_unrestricted_equal_weight.csv")
  outdir <- cfg$figure_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "source_data"), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(all_file), file.exists(top20_file))
  all338 <- read.csv(all_file, check.names = FALSE)
  top20 <- read.csv(top20_file, check.names = FALSE)
  required_cols <- c("Gene", "rho_TGRV", "p_TGRV", "rho_SUVmax", "p_SUVmax", "rho_Time", "p_Time", "score_equal_weight",
    "Rank_final_candidates")
  if (!all(required_cols %in% names(all338))) {
    stop("The 338-gene table lacks required columns: ", paste(setdiff(required_cols, names(all338)), collapse = ", "))
  }
  if (!all(required_cols %in% names(top20))) {
    stop("The Top20 table lacks required columns: ", paste(setdiff(required_cols, names(top20)), collapse = ", "))
  }
  all338 <- all338 %>% filter(is.finite(rho_TGRV), is.finite(rho_SUVmax), is.finite(rho_Time), is.finite(score_equal_weight)) %>%
    arrange(Rank_final_candidates)
  top20 <- top20 %>% filter(is.finite(rho_TGRV), is.finite(rho_SUVmax), is.finite(rho_Time), is.finite(score_equal_weight)) %>%
    arrange(Rank_final_candidates)
  stopifnot(nrow(all338) == 338L)
  stopifnot(nrow(top20) == 20L)
  stopifnot(setequal(top20$Gene, all338$Gene[all338$Rank_final_candidates <= 20]))
  all_score_limits <- range(all338$score_equal_weight, na.rm = TRUE)
  score_to_size <- function(x) {
    score_span <- diff(all_score_limits)
    if (!is.finite(score_span) || score_span <= 0)
      return(rep(7, length(x)))
    score_01 <- pmin(1, pmax(0, (x - all_score_limits[1])/score_span))
    3.5 + score_01 * (10.5 - 3.5)
  }
  all338$display_size <- score_to_size(all338$score_equal_weight)
  top20$display_size <- score_to_size(top20$score_equal_weight)
  make_hover <- function(dat) {
    paste0("<b>", dat$Gene, "</b>", "<br>Final rank: ", dat$Rank_final_candidates, "<br>Equal-weight score: ",
      sprintf("%.3f", dat$score_equal_weight), "<br>TGRV: rho = ", sprintf("%.3f", dat$rho_TGRV), ", P = ",
      format(dat$p_TGRV, digits = 3, scientific = TRUE), "<br>SUVmax: rho = ", sprintf("%.3f", dat$rho_SUVmax),
      ", P = ", format(dat$p_SUVmax, digits = 3, scientific = TRUE), "<br>Treatment time: rho = ", sprintf("%.3f",
        dat$rho_Time), ", P = ", format(dat$p_Time, digits = 3, scientific = TRUE))
  }
  all338$hover_text <- make_hover(all338)
  top20$hover_text <- make_hover(top20)
  camera_eye <- list(x = 1.55, y = 1.55, z = 1.1)
  camera_up <- list(x = 0, y = 0, z = 1)
  camera_center <- list(x = 0, y = 0, z = 0)
  axis_title_size <- 18
  axis_tick_size <- 13
  top20_gene_label_size <- 13
  top20_gene_label_family <- "Arial Black"
  leader_line_width <- 2.2
  leader_line_alpha <- 0.35
  score_palette <- viridisLite::viridis(256, option = "D", direction = 1)
  plotly_score_scale <- cbind(seq(0, 1, length.out = length(score_palette)), score_palette)
  scene_style <- list(camera = list(eye = camera_eye, up = camera_up, center = camera_center), aspectmode = "cube",
    xaxis = list(title = list(text = "Correlation with TGRV", font = list(size = axis_title_size)), tickfont = list(size = axis_tick_size),
      showline = TRUE, linewidth = 3, linecolor = "black", gridcolor = "#D9D9D9", zerolinecolor = "#8C8C8C",
      backgroundcolor = "white"), yaxis = list(title = list(text = "Correlation with SUVmax", font = list(size = axis_title_size)),
      tickfont = list(size = axis_tick_size), showline = TRUE, linewidth = 3, linecolor = "black", gridcolor = "#D9D9D9",
      zerolinecolor = "#8C8C8C", backgroundcolor = "white"), zaxis = list(title = list(text = "Correlation with treatment time",
      font = list(size = axis_title_size)), tickfont = list(size = axis_tick_size), showline = TRUE, linewidth = 3,
      linecolor = "black", gridcolor = "#D9D9D9", zerolinecolor = "#8C8C8C", backgroundcolor = "white"))
  common_layout <- list(scene = scene_style, paper_bgcolor = "white", plot_bgcolor = "white", font = list(family = "Arial",
    size = 14, color = "black"), margin = list(l = 20, r = 20, b = 20, t = 10), showlegend = FALSE)
  fig_338 <- plot_ly(data = all338, x = ~rho_TGRV, y = ~rho_SUVmax, z = ~rho_Time, type = "scatter3d", mode = "markers",
    text = ~hover_text, hoverinfo = "text", marker = list(color = all338$score_equal_weight, colorscale = plotly_score_scale,
      cmin = all_score_limits[1], cmax = all_score_limits[2], showscale = TRUE, size = all338$display_size,
      opacity = 0.78, line = list(color = "rgba(30,30,30,0.35)", width = 0.7),
      colorbar = list(title = list(text = "Association evidence<br>score", font = list(size = 15)), tickfont = list(size = 12),
        tickvals = c(0, 0.25, 0.5, 0.75, 1), tickformat = ".2f", thickness = 16, len = 0.66,
        x = 1.03, xanchor = "left"))) %>% layout(scene = common_layout$scene, paper_bgcolor = common_layout$paper_bgcolor,
    plot_bgcolor = common_layout$plot_bgcolor, font = common_layout$font, margin = common_layout$margin, showlegend = FALSE) %>%
    config(displaylogo = FALSE, toImageButtonOptions = list(format = "svg", filename = "Figure4_338_ICGs_3D",
      width = 1800, height = 1500, scale = 1))
  score_span <- diff(all_score_limits)
  if (!is.finite(score_span) || score_span <= 0)
    score_span <- 1
  top20$score_01 <- pmin(1, pmax(0, (top20$score_equal_weight - all_score_limits[1])/score_span))
  top20$palette_index <- 1L + round(top20$score_01 * (length(score_palette) - 1L))
  top20$point_colour <- score_palette[top20$palette_index]
  cross3 <- function(a, b) {
    c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
  }
  unit3 <- function(v) v/sqrt(sum(v^2))
  xyz_scaled <- cbind(scales::rescale(top20$rho_TGRV, to = c(-1, 1)), scales::rescale(top20$rho_SUVmax, to = c(-1,
    1)), scales::rescale(top20$rho_Time, to = c(-1, 1)))
  eye_vec <- c(camera_eye$x, camera_eye$y, camera_eye$z)
  up_vec <- c(camera_up$x, camera_up$y, camera_up$z)
  view_forward <- unit3(-eye_vec)
  view_right <- unit3(cross3(view_forward, up_vec))
  view_up <- unit3(cross3(view_right, view_forward))
  projected_x <- drop(xyz_scaled %*% view_right)
  projected_y <- drop(xyz_scaled %*% view_up)
  natural_angles <- atan2(projected_y - median(projected_y), projected_x - median(projected_x))
  natural_angles <- (natural_angles + 2 * pi)%%(2 * pi)
  n_labels <- nrow(top20)
  angle_order <- order(natural_angles)
  angle_grid <- seq(0, 2 * pi, length.out = n_labels + 1L)[seq_len(n_labels)]
  phase_shift <- atan2(mean(sin(natural_angles[angle_order] - angle_grid)), mean(cos(natural_angles[angle_order] -
    angle_grid)))
  label_angles <- numeric(n_labels)
  label_angles[angle_order] <- (angle_grid + phase_shift)%%(2 * pi)
  label_radius_px <- rep(c(62, 76, 90, 68, 82), length.out = n_labels)
  top20$label_ax <- round(label_radius_px * cos(label_angles))
  top20$label_ay <- -round(label_radius_px * sin(label_angles))
  label_annotations <- lapply(seq_len(nrow(top20)), function(i) {
    list(x = top20$rho_TGRV[i], y = top20$rho_SUVmax[i], z = top20$rho_Time[i], text = paste0("<b>", top20$Gene[i],
      "</b>"), showarrow = TRUE, ax = top20$label_ax[i], ay = top20$label_ay[i], xanchor = "center", yanchor = "middle",
      arrowhead = 0, arrowsize = 1, arrowwidth = leader_line_width, arrowcolor = scales::alpha(top20$point_colour[i],
        leader_line_alpha), font = list(family = top20_gene_label_family, size = top20_gene_label_size,
        color = "black"), bgcolor = "rgba(255,255,255,0.55)", borderpad = 1, opacity = 1)
  })
  top20_scene <- scene_style
  top20_scene$annotations <- label_annotations
  fig_top20 <- plot_ly() %>% add_trace(data = top20, x = ~rho_TGRV, y = ~rho_SUVmax, z = ~rho_Time, type = "scatter3d",
    mode = "markers", text = ~hover_text, hoverinfo = "text", marker = list(color = top20$score_equal_weight,
      colorscale = plotly_score_scale, cmin = all_score_limits[1], cmax = all_score_limits[2], showscale = TRUE,
      size = top20$display_size, opacity = 0.92, line = list(color = "black", width = 1.2),
      colorbar = list(title = list(text = "Association evidence<br>score", font = list(size = 15)), tickfont = list(size = 12),
        tickvals = c(0, 0.25, 0.5, 0.75, 1), tickformat = ".2f", thickness = 16, len = 0.66,
        x = 1.03, xanchor = "left")), showlegend = FALSE, name = "") %>% layout(scene = top20_scene, paper_bgcolor = common_layout$paper_bgcolor,
    plot_bgcolor = common_layout$plot_bgcolor, font = common_layout$font, margin = list(l = 30, r = 30, b = 30,
      t = 25), showlegend = FALSE) %>% config(displaylogo = FALSE, toImageButtonOptions = list(format = "svg",
    filename = "Figure4_Top20_ICGs_3D_labeled", width = 2000, height = 1700, scale = 1))
  fig_338
  fig_top20
  htmlwidgets::saveWidget(fig_338, file.path(outdir, "B_338_ICGs_3D_interactive.html"), selfcontained = FALSE)
  htmlwidgets::saveWidget(fig_top20, file.path(outdir, "F_Top20_ICGs_3D_labeled_interactive.html"), selfcontained = FALSE)
  write.csv(all338, file.path(outdir, "source_data", "B_338_ICGs_3D_plot_data.csv"), row.names = FALSE)
  write.csv(top20, file.path(outdir, "source_data", "F_Top20_ICGs_3D_plot_data.csv"), row.names = FALSE)
  message("Created interactive objects: fig_338 and fig_top20")
  message("Adjust camera_eye near the top of the script to change the viewing angle.")
  invisible(list(fig_338 = fig_338, fig_top20 = fig_top20))
}

acti_fig4_heatmap <- function(cfg = acti_fig4_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(svglite)
    library(ragg)
  })
  input_file <- file.path(cfg$score_dir, "05_Top20_direction_unrestricted_equal_weight.csv")
  outdir <- cfg$figure_dir
  source_dir <- file.path(outdir, "source_data")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(input_file))
  top20 <- read.csv(input_file, check.names = FALSE)
  required_columns <- c("Gene", "Rank_final_candidates", "percentile_TGRV", "percentile_SUVmax", "percentile_Time",
    "score_equal_weight")
  if (!all(required_columns %in% names(top20))) {
    stop("Missing required columns: ", paste(setdiff(required_columns, names(top20)), collapse = ", "))
  }
  top20 <- top20 %>% arrange(Rank_final_candidates)
  stopifnot(nrow(top20) == 20L)
  stopifnot(identical(top20$Rank_final_candidates, seq_len(20L)))
  recalculated_score <- rowMeans(top20[, c("percentile_TGRV", "percentile_SUVmax", "percentile_Time")], na.rm = FALSE)
  if (!isTRUE(all.equal(recalculated_score, top20$score_equal_weight, tolerance = 1e-12))) {
    stop("The stored composite score is not the equal-weight mean of the three endpoint scores.")
  }
  top20$gene_label <- top20$Gene
  heatmap_data <- top20 %>% select(Gene, Rank_final_candidates, gene_label, percentile_TGRV, percentile_SUVmax,
    percentile_Time, score_equal_weight) %>% pivot_longer(cols = c(percentile_TGRV, percentile_SUVmax, percentile_Time,
    score_equal_weight), names_to = "score_component", values_to = "score") %>% mutate(component_label = recode(score_component,
    percentile_TGRV = "TGRV", percentile_SUVmax = "SUVmax", percentile_Time = "Treatment\ntime", score_equal_weight = "Composite"),
    x_position = recode(score_component, percentile_TGRV = 1, percentile_SUVmax = 2.3, percentile_Time = 3.6,
      score_equal_weight = 5.15), gene_label = factor(gene_label, levels = rev(top20$gene_label)),
    score_text = sprintf("%.2f", score), text_colour = ifelse(score >= 0.955, "white", "black"),
    is_composite = score_component == "score_equal_weight")
  component_breaks <- c(1, 2.3, 3.6, 5.15)
  component_labels <- c("TGRV", "SUVmax", "Treatment\ntime", "Composite")
  p <- ggplot(heatmap_data, aes(x = x_position, y = gene_label, fill = score)) + geom_tile(data = subset(heatmap_data,
    !is_composite), width = 1.15, height = 0.9, colour = "white", linewidth = 0.55) +
    geom_tile(data = subset(heatmap_data, is_composite), width = 1.15, height = 0.9,
      colour = "black", linewidth = 0.72) + geom_text(aes(label = score_text, colour = text_colour),
    family = "Arial", fontface = "bold", size = 3.15, show.legend = FALSE) + scale_colour_identity() +
    scale_fill_gradientn(colours = c("#FBE1DE", "#F4B9B1", "#E67F73", "#CA493D", "#861C1C"), values = scales::rescale(c(0.9,
      0.925, 0.95, 0.975, 1), from = c(0.9, 1)),
      limits = c(0.9, 1), breaks = seq(0.9, 1, by = 0.02), labels = number_format(accuracy = 0.01),
      oob = squish, name = "Association\nevidence score") + scale_x_continuous(breaks = component_breaks,
    labels = component_labels, limits = c(0.35, 5.84), expand = c(0, 0), position = "top") +
    scale_y_discrete(drop = FALSE) + geom_vline(xintercept = 4.37, linewidth = 0.8,
    colour = "#4D4D4D") + coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + theme_classic(base_size = 10,
    base_family = "Arial") + theme(plot.title = element_blank(), plot.subtitle = element_blank(), axis.line = element_blank(),
    axis.ticks = element_blank(), axis.text.x = element_text(size = 10.2, face = "bold", colour = "black",
      margin = margin(b = 5)), axis.text.y = element_text(size = 9.2, face = "bold", colour = "black",
      margin = margin(r = 5)), panel.grid = element_blank(), legend.position = "right", legend.title = element_text(size = 9.2,
      face = "bold", lineheight = 0.95), legend.text = element_text(size = 8.4),
    legend.key.height = grid::unit(18, "mm"), legend.key.width = grid::unit(4.2, "mm"), plot.margin = margin(7,
      8, 6, 6, unit = "mm")) + guides(fill = guide_colourbar(title.position = "top", title.hjust = 0, barheight = grid::unit(46,
    "mm"), barwidth = grid::unit(4.2, "mm"), ticks.colour = "black", frame.colour = "black"))
  stem <- file.path(outdir, "F_Top20_equal_weight_score_heatmap")
  width_mm <- 145
  height_mm <- 132
  svglite::svglite(paste0(stem, ".svg"), width = width_mm/25.4, height = height_mm/25.4,
    bg = "white")
  print(p)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm/25.4, height = height_mm/25.4,
    bg = "white", family = "Arial")
  print(p)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm", res = 600, background = "white",
    compression = "lzw")
  print(p)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm", res = 300, background = "white")
  print(p)
  grDevices::dev.off()
  write.csv(heatmap_data %>% select(Gene, Rank_final_candidates, score_component, component_label, score), file.path(source_dir,
    "F_Top20_equal_weight_score_heatmap_data.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(outdir, "F_sessionInfo.txt"))
  message("Completed Figure 4F: ", paste0(stem, ".pdf"))
}

acti_fig4_overlap <- function(cfg = acti_fig4_config()) {
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(svglite)
    library(ragg)
  })
  sc_file <- cfg$sc_candidates_file
  stpmg_file <- cfg$stpmg_file
  outdir <- cfg$figure_dir
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
}

acti_fig4_validate <- function(cfg = acti_fig4_config()) {
  all <- read.csv(file.path(cfg$score_dir, "04_final_signed_Top3_scRNA_candidates_ranked.csv"), check.names = FALSE)
  top <- read.csv(file.path(cfg$score_dir, "05_Top20_direction_unrestricted_equal_weight.csv"), check.names = FALSE)
  stopifnot(nrow(all) == 338L, nrow(top) == 20L, !anyDuplicated(all$Gene))
  near <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-12))
  components <- matrix(NA_real_, nrow(all), 3L, dimnames = list(NULL, c("TGRV", "SUVmax", "Time")))
  for (endpoint in c("TGRV", "SUVmax", "Time")) {
    p <- all[[paste0("p_", endpoint)]]
    stopifnot(all(is.finite(p) & p >= 0 & p <= 1))
    evidence <- -log10(pmax(p, .Machine$double.xmin))
    percentile <- (rank(evidence, ties.method = "average") - 1)/(length(evidence) - 1)
    components[, endpoint] <- percentile
    stopifnot(near(evidence, all[[paste0("evidence_", endpoint)]]), near(percentile, all[[paste0("percentile_",
      endpoint)]]), near(p.adjust(p, "BH"), all[[paste0("FDR_", endpoint)]]))
    obs <- if (endpoint == "SUVmax")
      140L
    else 214L
    patients <- if (endpoint == "SUVmax")
      70L
    else 74L
    stopifnot(all(all[[paste0("n_observations_", endpoint)]] == obs), all(all[[paste0("n_patients_", endpoint)]] ==
      patients))
  }
  score <- rowMeans(components)
  ranks <- as.integer(rank(-score, ties.method = "min", na.last = "keep"))
  expected_top <- all[order(ranks, all$Gene), ][seq_len(20L), ]
  stopifnot(near(score, all$score_equal_weight), identical(ranks, all$Rank_final_candidates), identical(as.character(top$Gene),
    as.character(expected_top$Gene)), near(top$score_equal_weight, expected_top$score_equal_weight))
  print(top[, c("Gene", "Rank_final_candidates", "score_equal_weight")], row.names = FALSE)
  message("Validated: 338 candidates; three endpoint percentiles; equal-weight scores; Top20 ordering.")
  invisible(list(all = all, top20 = top))
}

acti_fig4_run <- function(stage, cfg = acti_fig4_config()) {
  result <- NULL
  for (s in acti_fig4_stages(stage)) {
    message("Running Figure 4 stage: ", s)
    result <- get(paste0("acti_fig4_", s), mode = "function")(cfg)
  }
  invisible(result)
}

acti_fig4_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- "Rscript 3_scoring_and_3D.R --stage NAME\nStages: all, score, validate, plot3d, heatmap, overlap\nUse --check to check inputs; see the configuration function for paths.\n"
  if (!length(args) || "--help" %in% args) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  values <- list()
  stage <- "all"
  check <- FALSE
  i <- 1L
  flags <- c("--data-dir", "--output-dir", "--sc-results-dir", "--wgcna-dir", "--week4-dir", "--week24-dir",
    "--score-dir", "--expression-file", "--phenotype-file", "--pet-file", "--suv-file", "--sc-all-tested-file",
    "--sc-candidates-file", "--stpmg-file", "--stage")
  while (i <= length(args)) {
    flag <- args[[i]]
    if (flag == "--check")
      check <- TRUE
    else if (flag == "--overwrite") {
      values$overwrite <- TRUE
    }
    else if (flag %in% flags) {
      if (i == length(args) || startsWith(args[[i + 1L]], "--"))
        stop("Missing value: ", flag)
      i <- i + 1L
      if (flag == "--stage")
        stage <- args[[i]]
      else {
        values[[gsub("-", "_", substring(flag, 3L), fixed = TRUE)]] <- args[[i]]
      }
    }
    else stop("Unknown argument: ", flag)
    i <- i + 1L
  }
  cfg <- do.call(acti_fig4_config, values)
  if (check) {
    if (!acti_fig4_check(cfg, stage))
      stop("Preflight incomplete: supply missing inputs/packages.")
    return(invisible(cfg))
  }
  acti_fig4_run(stage, cfg)
}

if (sys.nframe() == 0L) acti_fig4_main()

