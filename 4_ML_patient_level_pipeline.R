# Patient-level data splitting and machine learning.
# Inputs: Three expression parts, phenotype table and candidate genes; acti_ml_config().
# Run: Rscript 4_ML_patient_level_pipeline.R --stage prepare --split-seed 1234
# source() loads functions without starting an analysis.
# Type 0: DX/day7/week4; Type 1: week24. Split by patient.
# --split-seed reproduces a split; the model seed remains 1234.

acti_ml_config <- function(output_dir = file.path("results", "ML"), source_part_dir = file.path("data", "ML", "source_parts"),
  pheno_file = file.path("data", "GSE89403", "GSE89403_matrix_pheno.csv"), gene_file = file.path("data", "ML",
    "intersect.txt"), part_dir = NULL, input_dir = NULL, model_dir = NULL, split_seed = NULL, overwrite = FALSE) {
  default <- function(value, fallback) if (is.null(value))
    fallback
  else value
  list(source_part_dir = source_part_dir, pheno_file = pheno_file, gene_file = gene_file, part_dir = default(part_dir,
    file.path(output_dir, "01_parts")), input_dir = default(input_dir, file.path(output_dir, "02_input")),
    model_dir = default(model_dir, file.path(output_dir, "03_models")), split_seed = split_seed, overwrite = isTRUE(overwrite))
}

acti_ml_packages <- function() c("randomForestSRC", "glmnet", "plsRglm", "gbm", "mboost", "e1071", "MASS", "xgboost",
  "pROC")

acti_ml_load_packages <- function() {
  packages <- acti_ml_packages()
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!available))
    stop("Missing R packages: ", paste(packages[!available], collapse = ", "))
  suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))
}

acti_ml_stages <- function(stage) {
  if (length(stage) != 1L || !(stage %in% c("all", "prepare", "train", "validate"))) {
    stop("Choose stage: all, prepare, train, validate")
  }
  if (stage == "all")
    c("prepare", "validate", "train")
  else stage
}

acti_ml_check <- function(cfg = acti_ml_config(), stage = "all") {
  stages <- acti_ml_stages(stage)
  files <- if ("prepare" %in% stages)
    c(file.path(cfg$source_part_dir, paste0("GSE89403part", 1:3, ".normalize.txt")), cfg$pheno_file, cfg$gene_file)
  else file.path(cfg$input_dir, c("data.train.txt", "data.test.txt"))
  print(data.frame(path = files, exists = file.exists(files)), row.names = FALSE)
  ok <- all(file.exists(files))
  if ("train" %in% stages) {
    packages <- acti_ml_packages()
    found <- vapply(packages, function(p) nzchar(system.file(package = p)), logical(1))
    print(data.frame(package = packages, available = found), row.names = FALSE)
    ok <- ok && all(found)
  }
  invisible(ok)
}

acti_ml_prepare <- function(cfg = acti_ml_config()) {
  absolute <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)
  outputs <- c(cfg$part_dir, cfg$input_dir)
  for (path in outputs) {
    if (tolower(absolute(path)) == tolower(absolute(cfg$source_part_dir))) {
      stop("Preparation output must differ from immutable source_part_dir.")
    }
    if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) && !cfg$overwrite) {
      stop("Preparation output directory is not empty: ", path, ". Use a new directory or explicitly supply --overwrite.")
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  options(stringsAsFactors = FALSE)
  part_dir <- cfg$part_dir
  ml_out_dir <- cfg$input_dir
  source_part_dir <- cfg$source_part_dir
  pheno_file <- cfg$pheno_file
  gene_file <- cfg$gene_file
  part_files <- file.path(part_dir, c("GSE89403part1.normalize.txt", "GSE89403part2.normalize.txt", "GSE89403part3.normalize.txt"))
  source_part_files <- file.path(source_part_dir, basename(part_files))
  merge_file <- file.path(part_dir, "merge.normalize.txt")
  train_file <- file.path(ml_out_dir, "data.train.txt")
  test_file <- file.path(ml_out_dir, "data.test.txt")
  manifest_file <- file.path(ml_out_dir, "patient_level_train_test_split_manifest.csv")
  summary_file <- file.path(ml_out_dir, "patient_level_train_test_split_summary.txt")
  excluded_subject_file <- file.path(ml_out_dir, "excluded_subjects_without_week24_and_earlier_samples.csv")
  train_fraction <- 0.7
  seed_override <- if (is.null(cfg$split_seed)) {
    trimws(Sys.getenv("ML_SPLIT_SEED", unset = ""))
  }
  else trimws(as.character(cfg$split_seed))
  if (nzchar(seed_override)) {
    run_seed <- suppressWarnings(as.integer(seed_override))
    if (is.na(run_seed) || run_seed < 1L) {
      stop("\u73af\u5883\u53d8\u91cfML_SPLIT_SEED\u5fc5\u987b\u662f1\u81f32147483647\u4e4b\u95f4\u7684\u6574\u6570")
    }
  }
  else {
    run_seed <- as.integer((as.numeric(Sys.time()) * 1000000 + Sys.getpid())%%(.Machine$integer.max - 1)) +
      1L
  }
  part_seed <- run_seed
  split_seed <- if (run_seed == .Machine$integer.max)
    1L
  else run_seed + 1L
  cat("\u672c\u6b21\u968f\u673aRun seed\uff1a", run_seed, "\n", sep = "")
  cat("\u5982\u9700\u590d\u73b0\uff1aSys.setenv(ML_SPLIT_SEED = \"", run_seed, "\")\n", sep = "")
  if (!dir.exists(ml_out_dir)) {
    dir.create(ml_out_dir, recursive = TRUE)
  }
  required_files <- c(source_part_files, pheno_file, gene_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("\u4ee5\u4e0b\u8f93\u5165\u6587\u4ef6\u4e0d\u5b58\u5728\uff1a\n", paste(missing_files, collapse = "\n"))
  }
  pheno <- read.csv(pheno_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  required_pheno_columns <- c("subject", "treatmentresult", "timetonegativity")
  if (!all(required_pheno_columns %in% colnames(pheno))) {
    stop("pheno\u7f3a\u5c11\u4ee5\u4e0b\u5b57\u6bb5\uff1a", paste(setdiff(required_pheno_columns, colnames(pheno)), collapse = ", "))
  }
  valid_subject <- !is.na(pheno$subject) & grepl("^[0-9]+$", as.character(pheno$subject))
  pheno <- pheno[valid_subject, , drop = FALSE]
  pheno$subject_id <- paste0("S", pheno$subject)
  pheno$neg_week <- suppressWarnings(as.numeric(sub("^Week0*", "", pheno$timetonegativity)))
  subject_pheno <- unique(pheno[, c("subject_id", "treatmentresult", "timetonegativity", "neg_week")])
  duplicated_subjects <- subject_pheno$subject_id[duplicated(subject_pheno$subject_id)]
  if (length(duplicated_subjects) > 0) {
    stop("\u4ee5\u4e0b\u60a3\u8005\u5728pheno\u4e2d\u5b58\u5728\u4e0d\u4e00\u81f4\u7684\u7ed3\u5c40\u6216\u8f6c\u9634\u65f6\u95f4\uff1a", paste(unique(duplicated_subjects),
      collapse = ", "))
  }
  rownames(subject_pheno) <- subject_pheno$subject_id
  read_part <- function(path) {
    x <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
    colnames(x)[1] <- "id"
    if (anyDuplicated(x$id)) {
      stop("\u6587\u4ef6\u4e2d\u5b58\u5728\u91cd\u590d\u57fa\u56e0\u540d\uff1a", basename(path))
    }
    x
  }
  parts <- lapply(source_part_files, read_part)
  reference_genes <- parts[[1]]$id
  for (i in seq_along(parts)) {
    if (!identical(parts[[i]]$id, reference_genes)) {
      stop("\u4e09\u4e2apart\u6587\u4ef6\u7684\u57fa\u56e0\u884c\u6216\u57fa\u56e0\u987a\u5e8f\u4e0d\u4e00\u81f4\uff1a", basename(source_part_files[i]))
    }
  }
  for (i in seq_along(parts)) {
    colnames(parts[[i]])[-1] <- sub("^GSE89403part[123]_", "", colnames(parts[[i]])[-1])
  }
  raw_sample_names <- unlist(lapply(parts, function(x) colnames(x)[-1]), use.names = FALSE)
  if (anyDuplicated(raw_sample_names)) {
    stop("\u6c47\u603b\u4e09\u4e2apart\u540e\u53d1\u73b0\u91cd\u590d\u6837\u672c\u540d\uff1a", paste(unique(raw_sample_names[duplicated(raw_sample_names)]),
      collapse = ", "))
  }
  expr_pool <- parts[[1]]
  for (i in 2:length(parts)) {
    expr_pool <- cbind(expr_pool, parts[[i]][, -1, drop = FALSE])
  }
  canonical_sample_order <- order(raw_sample_names)
  raw_sample_names <- raw_sample_names[canonical_sample_order]
  expr_pool <- cbind(expr_pool[, 1, drop = FALSE], expr_pool[, canonical_sample_order + 1L, drop = FALSE])
  source_sample_pattern <- paste0("^GSMS([0-9]+)", "(DX|day7|week4|week24)_(Control|Treat)$")
  valid_raw_name <- grepl(source_sample_pattern, raw_sample_names)
  if (!all(valid_raw_name)) {
    stop("\u4ee5\u4e0b\u539f\u59cb\u6837\u672c\u540d\u65e0\u6cd5\u89e3\u6790\uff1a\n", paste(raw_sample_names[!valid_raw_name], collapse = "\n"))
  }
  raw_subject_id <- sub(source_sample_pattern, "S\\1", raw_sample_names)
  raw_timepoint <- sub(source_sample_pattern, "\\2", raw_sample_names)
  raw_label_in_source <- sub(source_sample_pattern, "\\3", raw_sample_names)
  all_part_subjects <- sort(unique(raw_subject_id))
  missing_part_pheno <- setdiff(all_part_subjects, rownames(subject_pheno))
  if (length(missing_part_pheno) > 0) {
    stop("\u4ee5\u4e0b\u60a3\u8005\u5728pheno\u4e2d\u4e0d\u5b58\u5728\uff1a", paste(missing_part_pheno, collapse = ", "))
  }
  raw_treatment_result <- subject_pheno[raw_subject_id, "treatmentresult"]
  if (any(is.na(raw_treatment_result)) || any(raw_treatment_result != "Definite Cure")) {
    stop("\u539f\u59cbpart\u5305\u542b\u975eDefinite Cure\u60a3\u8005\u6216\u7ed3\u5c40\u7f3a\u5931\uff1a", paste(unique(raw_subject_id[is.na(raw_treatment_result) |
      raw_treatment_result != "Definite Cure"]), collapse = ", "))
  }
  has_week24 <- tapply(raw_timepoint == "week24", raw_subject_id, any)
  has_earlier <- tapply(raw_timepoint %in% c("DX", "day7", "week4"), raw_subject_id, any)
  eligible_subjects <- intersect(names(has_earlier)[has_earlier], names(has_week24)[has_week24])
  excluded_subjects <- setdiff(all_part_subjects, eligible_subjects)
  excluded_table <- data.frame(subject_id = excluded_subjects, reason = ifelse(!has_earlier[excluded_subjects] &
    !has_week24[excluded_subjects], "No earlier sample and no week24 sample", ifelse(!has_earlier[excluded_subjects],
    "No DX/day7/week4 sample", "No week24 sample")), stringsAsFactors = FALSE)
  write.csv(excluded_table, file = excluded_subject_file, row.names = FALSE, quote = TRUE)
  keep_sample <- raw_subject_id %in% eligible_subjects & raw_timepoint %in% c("DX", "day7", "week4", "week24")
  raw_sample_names <- raw_sample_names[keep_sample]
  raw_subject_id <- raw_subject_id[keep_sample]
  raw_timepoint <- raw_timepoint[keep_sample]
  raw_label_in_source <- raw_label_in_source[keep_sample]
  expr_pool <- cbind(expr_pool[, 1, drop = FALSE], expr_pool[, which(keep_sample) + 1L, drop = FALSE])
  all_part_subjects <- sort(unique(raw_subject_id))
  sample_key <- paste(raw_subject_id, raw_timepoint, sep = "__")
  if (anyDuplicated(sample_key)) {
    stop("\u540c\u4e00\u60a3\u8005\u540c\u4e00\u65f6\u95f4\u70b9\u5b58\u5728\u91cd\u590d\u8868\u8fbe\u6837\u672c")
  }
  check_has_week24 <- tapply(raw_timepoint == "week24", raw_subject_id, any)
  check_has_earlier <- tapply(raw_timepoint %in% c("DX", "day7", "week4"), raw_subject_id, any)
  if (any(!check_has_week24) || any(!check_has_earlier)) {
    stop("\u5b8c\u6574\u6027\u68c0\u67e5\u5931\u8d25\uff1a\u6bcf\u540d\u60a3\u8005\u5fc5\u987b\u6709week24\u53ca\u81f3\u5c11\u4e00\u4e2a\u66f4\u65e9\u65f6\u95f4\u70b9")
  }
  raw_corrected_label <- ifelse(raw_timepoint == "week24", "Endpoint", "NonEndpoint")
  raw_corrected_sample_names <- paste0("GSMS", sub("^S", "", raw_subject_id), raw_timepoint, "_", raw_corrected_label)
  if (anyDuplicated(raw_corrected_sample_names)) {
    stop("\u6309pheno\u4fee\u6b63\u6807\u7b7e\u540e\u4ea7\u751f\u91cd\u590d\u6837\u672c\u540d")
  }
  source_label_changes <- sum(raw_label_in_source != raw_corrected_label)
  raw_sample_names <- raw_corrected_sample_names
  colnames(expr_pool)[-1] <- raw_sample_names
  analysis_sample_pattern <- paste0("^GSMS([0-9]+)", "(DX|day7|week4|week24)_(NonEndpoint|Endpoint)$")
  cat("\u4ece\u5b8c\u6574\u5907\u4efd\u8bfb\u53d6\uff1a", length(raw_sample_names), "\u4e2a\u6837\u672c\uff1b\u7eb3\u5165", length(all_part_subjects),
    "\u540d\u5177\u6709week24\u53ca\u81f3\u5c11\u4e00\u4e2a\u66f4\u65e9\u65f6\u95f4\u70b9\u7684Definite Cure\u60a3\u8005\uff1b\u6392\u9664", length(excluded_subjects),
    "\u540d\u4e0d\u5b8c\u6574\u60a3\u8005\uff1b", "\u5ffd\u7565\u5e76\u66ff\u6362", source_label_changes, "\u4e2a\u65e7Control/Treat\u6807\u7b7e\n", sep = "")
  set.seed(part_seed)
  subject_part_number <- setNames(integer(length(all_part_subjects)), all_part_subjects)
  part_strata <- subject_pheno[all_part_subjects, "timetonegativity"]
  for (stratum in sort(unique(part_strata))) {
    stratum_subjects <- all_part_subjects[part_strata == stratum]
    stratum_subjects <- sample(stratum_subjects, length(stratum_subjects))
    random_start <- sample(0:2, size = 1)
    assigned_part <- ((seq_along(stratum_subjects) - 1L + random_start)%%3L) + 1L
    subject_part_number[stratum_subjects] <- assigned_part
  }
  sample_part_number <- unname(subject_part_number[raw_subject_id])
  parts <- lapply(seq_along(part_files), function(i) {
    sample_index <- which(sample_part_number == i)
    sample_index <- sample_index[sample.int(length(sample_index))]
    cbind(data.frame(id = expr_pool$id, check.names = FALSE), expr_pool[, sample_index + 1L, drop = FALSE])
  })
  redistributed_names <- unlist(lapply(parts, function(x) colnames(x)[-1]), use.names = FALSE)
  if (!setequal(raw_sample_names, redistributed_names) || length(raw_sample_names) != length(redistributed_names)) {
    stop("\u91cd\u65b0\u5206\u914dpart\u540e\u6837\u672c\u96c6\u5408\u53d1\u751f\u6539\u53d8")
  }
  subject_part_check <- unique(data.frame(subject_id = sub(analysis_sample_pattern, "S\\1", redistributed_names),
    part = rep(paste0("part", seq_along(parts)), lengths(lapply(parts, function(x) colnames(x)[-1]))), stringsAsFactors = FALSE))
  if (anyDuplicated(subject_part_check$subject_id)) {
    stop("\u91cd\u65b0\u5206\u914d\u540e\u53d1\u73b0\u540c\u4e00\u60a3\u8005\u51fa\u73b0\u5728\u591a\u4e2apart\u4e2d")
  }
  for (i in seq_along(parts)) {
    write.table(parts[[i]], file = part_files[i], sep = "\t", quote = FALSE, row.names = FALSE)
  }
  for (i in seq_along(parts)) {
    colnames(parts[[i]])[-1] <- paste0("GSE89403part", i, "_", colnames(parts[[i]])[-1])
  }
  expr_all <- parts[[1]]
  for (i in 2:length(parts)) {
    expr_all <- cbind(expr_all, parts[[i]][, -1, drop = FALSE])
  }
  colnames(expr_all)[1] <- "geneNames"
  if (anyDuplicated(colnames(expr_all)[-1])) {
    stop("\u5408\u5e76\u540e\u51fa\u73b0\u91cd\u590d\u6837\u672c\u540d")
  }
  write.table(expr_all, file = merge_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\u5408\u5e76\u8868\u8fbe\u77e9\u9635\uff1a", nrow(expr_all), "\u4e2a\u57fa\u56e0 \u00d7", ncol(expr_all) - 1, "\u4e2a\u6837\u672c\n")
  sample_names <- colnames(expr_all)[-1]
  sample_pattern <- paste0("^GSE89403part([123])_GSMS([0-9]+)", "(DX|day7|week4|week24)_(NonEndpoint|Endpoint)$")
  valid_sample_name <- grepl(sample_pattern, sample_names)
  if (!all(valid_sample_name)) {
    stop("\u4ee5\u4e0b\u6837\u672c\u540d\u65e0\u6cd5\u89e3\u6790\uff1a\n", paste(sample_names[!valid_sample_name], collapse = "\n"))
  }
  part_id <- sub(sample_pattern, "part\\1", sample_names)
  subject_id <- sub(sample_pattern, "S\\2", sample_names)
  timepoint <- sub(sample_pattern, "\\3", sample_names)
  label_in_column <- sub(sample_pattern, "\\4", sample_names)
  time_to_week <- c(DX = 0, day7 = 1, week4 = 4, week24 = 24)
  sample_week <- unname(time_to_week[timepoint])
  sample_meta <- data.frame(sample_id = sample_names, part = part_id, subject_id = subject_id, timepoint = timepoint,
    sample_week = sample_week, label_in_column = label_in_column, stringsAsFactors = FALSE)
  sample_meta$treatmentresult <- subject_pheno[sample_meta$subject_id, "treatmentresult"]
  sample_meta$timetonegativity <- subject_pheno[sample_meta$subject_id, "timetonegativity"]
  sample_meta$neg_week <- subject_pheno[sample_meta$subject_id, "neg_week"]
  if (any(is.na(sample_meta$neg_week))) {
    stop("\u4ee5\u4e0b\u6837\u672c\u7f3a\u5c11\u53ef\u7528\u7684timetonegativity\uff1a\n", paste(sample_meta$sample_id[is.na(sample_meta$neg_week)],
      collapse = "\n"))
  }
  if (any(sample_meta$treatmentresult != "Definite Cure")) {
    stop("\u5f53\u524d\u811a\u672c\u9884\u671f\u4ec5\u4f7f\u7528Definite Cure\u60a3\u8005\uff0c\u4f46\u53d1\u73b0\u5176\u4ed6\u7ed3\u5c40\uff1a", paste(unique(sample_meta$treatmentresult[sample_meta$treatmentresult !=
      "Definite Cure"]), collapse = ", "))
  }
  sample_meta$group <- ifelse(sample_meta$timepoint == "week24", "Endpoint", "NonEndpoint")
  label_mismatch <- sample_meta$label_in_column != sample_meta$group
  if (any(label_mismatch)) {
    stop("part\u6587\u4ef6\u5217\u540d\u4ecd\u6709", sum(label_mismatch), "\u4e2a\u6807\u7b7e\u4e0epheno\u89c4\u5219\u4e0d\u4e00\u81f4\u3002\u8bf7\u5148\u4fee\u6b63\u4ee5\u4e0b\u6837\u672c\uff1a\n",
      paste(sample_meta$sample_id[label_mismatch], collapse = "\n"))
  }
  candidate_genes <- unique(read.table(gene_file, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)[,
    1])
  available_genes <- intersect(candidate_genes, expr_all$geneNames)
  missing_candidate_genes <- setdiff(candidate_genes, expr_all$geneNames)
  if (length(available_genes) == 0) {
    stop("intersect.txt\u4e2d\u7684\u5019\u9009\u57fa\u56e0\u5747\u672a\u5728\u8868\u8fbe\u77e9\u9635\u4e2d\u627e\u5230")
  }
  if (length(missing_candidate_genes) > 0) {
    warning("\u4ee5\u4e0b\u5019\u9009\u57fa\u56e0\u672a\u5728\u8868\u8fbe\u77e9\u9635\u4e2d\u627e\u5230\uff1a", paste(missing_candidate_genes, collapse = ", "))
  }
  gene_index <- match(available_genes, expr_all$geneNames)
  model_expr <- t(as.matrix(expr_all[gene_index, -1, drop = FALSE]))
  storage.mode(model_expr) <- "numeric"
  colnames(model_expr) <- available_genes
  rownames(model_expr) <- sample_names
  if (!identical(rownames(model_expr), sample_meta$sample_id)) {
    stop("\u8868\u8fbe\u77e9\u9635\u6837\u672c\u987a\u5e8f\u4e0esample_meta\u4e0d\u4e00\u81f4")
  }
  subject_split_table <- unique(sample_meta[, c("subject_id", "timetonegativity", "neg_week")])
  if (anyDuplicated(subject_split_table$subject_id)) {
    stop("\u60a3\u8005\u5bf9\u5e94\u4e86\u591a\u4e2a\u4e0d\u540c\u7684timetonegativity")
  }
  set.seed(split_seed)
  train_subjects <- unlist(lapply(split(subject_split_table$subject_id, subject_split_table$timetonegativity),
    function(ids) {
      ids <- sort(unique(ids))
      if (length(ids) == 1) {
        return(ids)
      }
      n_train <- round(length(ids) * train_fraction)
      n_train <- max(1, min(length(ids) - 1, n_train))
      sample(ids, size = n_train, replace = FALSE)
    }), use.names = FALSE)
  all_subjects <- sort(unique(subject_split_table$subject_id))
  train_subjects <- sort(unique(train_subjects))
  test_subjects <- setdiff(all_subjects, train_subjects)
  sample_meta$dataset <- ifelse(sample_meta$subject_id %in% train_subjects, "Training", "Internal_test")
  train_index <- sample_meta$dataset == "Training"
  test_index <- sample_meta$dataset == "Internal_test"
  subject_overlap <- intersect(unique(sample_meta$subject_id[train_index]), unique(sample_meta$subject_id[test_index]))
  if (length(subject_overlap) > 0) {
    stop("\u53d1\u73b0\u60a3\u8005\u7ea7\u6570\u636e\u6cc4\u6f0f\uff1a", paste(subject_overlap, collapse = ", "))
  }
  sample_overlap <- intersect(sample_meta$sample_id[train_index], sample_meta$sample_id[test_index])
  if (length(sample_overlap) > 0) {
    stop("\u53d1\u73b0\u6837\u672c\u7ea7\u6570\u636e\u6cc4\u6f0f\uff1a", paste(sample_overlap, collapse = ", "))
  }
  if (length(unique(sample_meta$group[train_index])) < 2) {
    stop("\u8bad\u7ec3\u96c6\u6ca1\u6709\u540c\u65f6\u5305\u542bNonEndpoint\u548cEndpoint")
  }
  if (length(unique(sample_meta$group[test_index])) < 2) {
    stop("\u5185\u90e8\u6d4b\u8bd5\u96c6\u6ca1\u6709\u540c\u65f6\u5305\u542bNonEndpoint\u548cEndpoint")
  }
  train_expr <- model_expr[train_index, , drop = FALSE]
  test_expr <- model_expr[test_index, , drop = FALSE]
  train_type <- ifelse(sample_meta$group[train_index] == "Endpoint", 1, 0)
  test_type <- ifelse(sample_meta$group[test_index] == "Endpoint", 1, 0)
  train_output <- data.frame(id = paste0("Train.", rownames(train_expr)), train_expr, Type = train_type, check.names = FALSE)
  test_output <- data.frame(id = rownames(test_expr), test_expr, Type = test_type, check.names = FALSE)
  write.table(train_output, file = train_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(test_output, file = test_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  sample_meta$Type <- ifelse(sample_meta$group == "Endpoint", 1, 0)
  sample_meta$label_definition <- ifelse(sample_meta$Type == 1, "Week24 end-of-treatment blood draw in Definite Cure patients",
    "DX/day7/week4 pre-/mid-treatment blood draw in Definite Cure patients")
  sample_meta$run_seed <- run_seed
  sample_meta$part_seed <- part_seed
  sample_meta$split_seed <- split_seed
  sample_meta$train_fraction <- train_fraction
  write.csv(sample_meta, file = manifest_file, quote = TRUE, row.names = FALSE)
  train_check <- read.table(train_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  test_check <- read.table(test_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  train_original_ids <- sub("^Train\\.", "", train_check$id)
  if (length(intersect(train_original_ids, test_check$id)) > 0) {
    stop("\u5199\u51fa\u6587\u4ef6\u540e\u68c0\u67e5\u5230\u8bad\u7ec3\u96c6\u4e0e\u6d4b\u8bd5\u96c6\u5b58\u5728\u91cd\u590d\u6837\u672c")
  }
  train_check_subject <- sub(sample_pattern, "S\\2", train_original_ids)
  test_check_subject <- sub(sample_pattern, "S\\2", test_check$id)
  if (length(intersect(train_check_subject, test_check_subject)) > 0) {
    stop("\u5199\u51fa\u6587\u4ef6\u540e\u68c0\u67e5\u5230\u540c\u4e00\u60a3\u8005\u540c\u65f6\u8fdb\u5165\u8bad\u7ec3\u96c6\u548c\u6d4b\u8bd5\u96c6")
  }
  summary_lines <- c(paste0("Run seed: ", run_seed), paste0("Part-allocation seed: ", part_seed), paste0("Train/test split seed: ",
    split_seed), paste0("Reproduce with: Sys.setenv(ML_SPLIT_SEED = \"", run_seed, "\")"), paste0("Training fraction: ",
    train_fraction), "Label 0: NonEndpoint = DX/day7/week4 pre-/mid-treatment blood draws in Definite Cure patients",
    "Label 1: Endpoint = week24/M6 end-of-treatment blood draw in Definite Cure patients", "timetonegativity is NOT used to define Type; it is used only for stratified patient allocation",
    paste0("Excluded subjects without week24 and at least one earlier sample: ", length(excluded_subjects)),
    paste0("Excluded subject IDs: ", ifelse(length(excluded_subjects), paste(excluded_subjects, collapse = ", "),
      "None")), paste0("Total subjects: ", length(all_subjects)), paste0("Training subjects: ", length(train_subjects)),
    paste0("Internal-test subjects: ", length(test_subjects)), paste0("Total samples: ", nrow(sample_meta)),
    paste0("DX samples: ", sum(sample_meta$timepoint == "DX")), paste0("Day7 samples: ", sum(sample_meta$timepoint ==
      "day7")), paste0("Week4 samples: ", sum(sample_meta$timepoint == "week4")), paste0("Week24 samples: ",
      sum(sample_meta$timepoint == "week24")), paste0("Training samples: ", nrow(train_output)), paste0("Internal-test samples: ",
      nrow(test_output)), paste0("Candidate genes supplied: ", length(candidate_genes)), paste0("Candidate genes retained: ",
      length(available_genes)), paste0("Patient overlap: ", length(subject_overlap)), paste0("Sample overlap: ",
      length(sample_overlap)), "Caution: This split checks patient/sample overlap only. If candidate genes were selected using all GSE89403 subjects, internal-test AUC is not fully independent.",
    "", "Randomized part subject counts:", capture.output(print(table(unique(sample_meta[, c("subject_id",
      "part")])$part))), "", "Randomized part sample counts:", capture.output(print(table(sample_meta$part))),
    "", "Randomized part subjects by timetonegativity:", capture.output(print(table(unique(sample_meta[, c("subject_id",
      "part", "timetonegativity")])$part, unique(sample_meta[, c("subject_id", "part", "timetonegativity")])$timetonegativity))),
    "", "Training label counts:", capture.output(print(table(sample_meta$group[train_index]))), "", "Internal-test label counts:",
    capture.output(print(table(sample_meta$group[test_index]))), "", "Training timepoint by label:", capture.output(print(table(sample_meta$timepoint[train_index],
      sample_meta$group[train_index]))), "", "Internal-test timepoint by label:", capture.output(print(table(sample_meta$timepoint[test_index],
      sample_meta$group[test_index]))), "", "Subjects by timetonegativity and dataset:", capture.output(print(table(subject_split_table$timetonegativity,
      ifelse(subject_split_table$subject_id %in% train_subjects, "Training", "Internal_test")))))
  writeLines(summary_lines, summary_file)
  cat(paste(summary_lines, collapse = "\n"), "\n")
  cat("\n\u5df2\u751f\u6210\uff1a\n")
  cat(train_file, "\n")
  cat(test_file, "\n")
  cat(manifest_file, "\n")
  cat(summary_file, "\n")
  cat(excluded_subject_file, "\n")
  cat("\n\u6ce8\u610f\uff1adata.test.txt\u53ea\u80fd\u7528\u4e8e\u6700\u7ec8\u5185\u90e8\u9a8c\u8bc1\uff0c", "\u4e0d\u5f97\u7528\u4e8e\u4ece\u5019\u9009\u7b97\u6cd5\u4e2d\u9009\u62e9\u6700\u4f73\u6a21\u578b\u3002\n",
    sep = "")
  invisible(list(input_dir = ml_out_dir, manifest = sample_meta, seed = run_seed, train_file = train_file, test_file = test_file))
}

acti_ml_models <- function(cfg = acti_ml_config(), functions_only = FALSE) {
  if (!functions_only)
    acti_ml_load_packages()
  RunML <- function(method, Train_set, Train_label, mode = "Model", classVar) {
    method = gsub(" ", "", method)
    method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
    method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)
    method_param = switch(EXPR = method_name, Enet = list(alpha = as.numeric(gsub("alpha=", "", method_param))),
      Stepglm = list(direction = method_param), NULL)
    message("Run ", method_name, " algorithm for ", mode, "; ", method_param, ";", " using ", ncol(Train_set),
      " Variables")
    args = list(Train_set = Train_set, Train_label = Train_label, mode = mode, classVar = classVar)
    args = c(args, method_param)
    obj <- do.call(what = paste0("Run", method_name), args = args)
    if (mode == "Variable") {
      message(length(obj), " Variables retained;\n")
    }
    else {
      message("\n")
    }
    return(obj)
  }
  RunLasso <- function(Train_set, Train_label, mode, classVar) {
    RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
  }
  RunRidge <- function(Train_set, Train_label, mode, classVar) {
    RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
  }
  RunStepglm <- function(Train_set, Train_label, mode, classVar, direction) {
    fit <- step(glm(formula = Train_label[[classVar]] ~ ., family = "binomial", data = as.data.frame(Train_set)),
      direction = direction, trace = 0)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model")
      return(fit)
    if (mode == "Variable")
      return(ExtractVar(fit))
  }
  RunRF <- function(Train_set, Train_label, mode, classVar) {
    rf_nodesize = 5
    Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
    fit <- rfsrc(formula = formula(paste0(classVar, "~.")), data = cbind(Train_set, Train_label[classVar]),
      ntree = 1000, nodesize = rf_nodesize, importance = T, proximity = T, forest = T)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model")
      return(fit)
    if (mode == "Variable")
      return(ExtractVar(fit))
  }
  RunNaiveBayes <- function(Train_set, Train_label, mode, classVar) {
    data <- cbind(Train_set, Train_label[classVar])
    data[[classVar]] <- as.factor(data[[classVar]])
    fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))), data = data)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model")
      return(fit)
    if (mode == "Variable")
      return(ExtractVar(fit))
  }
  quiet <- function(..., messages = FALSE, cat = FALSE) {
    if (!cat) {
      sink(tempfile())
      on.exit(sink())
    }
    out <- if (messages)
      eval(...)
    else suppressMessages(eval(...))
    out
  }
  ExtractVar <- function(fit) {
    Feature <- quiet(switch(EXPR = class(fit)[1], lognet = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
      glm = names(coef(fit)), svm.formula = fit$subFeature, train = fit$coefnames, glmboost = names(coef(fit)[abs(coef(fit)) >
        0]), plsRglmmodel = rownames(fit$Coeffs)[fit$Coeffs != 0], rfsrc = var.select(fit, verbose = F)$topvars,
      gbm = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf > 0], xgb.Booster = fit$subFeature,
      naiveBayes = fit$subFeature))
    Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
    return(Feature)
  }
  input_dir <- cfg$input_dir
  output_root <- cfg$model_dir
  ml_seed <- 1234L
  outer_k <- 5L
  inner_k <- 5L
  min_retained_features <- 1L
  gbm_max_trees <- 10000L
  gbm_tree_grid <- seq.int(100L, gbm_max_trees, 100L)
  xgb_max_rounds <- 10L
  fit_scaler <- function(x) {
    s <- apply(x, 2, sd)
    keep <- names(s)[is.finite(s) & s > 1e-8]
    if (!length(keep))
      stop("\u8bad\u7ec3\u6298\u65e0\u975e\u96f6\u65b9\u5dee\u57fa\u56e0")
    list(features = keep, center = colMeans(x[, keep, drop = FALSE]), scale = s[keep])
  }
  apply_scaler <- function(x, s) {
    if (!all(s$features %in% colnames(x)))
      stop("\u9884\u6d4b\u6570\u636e\u7f3a\u5c11\u57fa\u56e0")
    sweep(sweep(x[, s$features, drop = FALSE], 2, s$center, "-"), 2, s$scale, "/")
  }
  make_group_folds <- function(ids, y, k = 5L, seed = ml_seed) {
    if (length(ids) != length(y) || anyNA(ids) || any(!nzchar(ids)))
      stop("\u60a3\u8005ID\u65e0\u6548")
    ss <- sort(unique(ids))
    strata <- vapply(ss, function(s) paste(sort(unique(y[ids == s])), collapse = ""), character(1))
    if (length(ss) < 4L)
      stop("\u60a3\u8005\u6570\u8fc7\u5c11")
    for (kk in seq.int(min(k, length(ss)), 2L)) for (attempt in seq_len(100L)) {
      set.seed(seed + kk * 100L + attempt)
      assignment <- setNames(integer(length(ss)), ss)
      counts <- integer(kk)
      for (stratum in unique(strata)) {
        zz <- ss[strata == stratum]
        zz <- zz[sample.int(length(zz))]
        for (s in zz) {
         candidates <- which(counts == min(counts))
         f <- candidates[sample.int(length(candidates), 1L)]
         assignment[s] <- f
         counts[f] <- counts[f] + 1L
        }
      }
      fold <- unname(assignment[ids])
      ok <- all(vapply(seq_len(kk), function(f) length(unique(y[fold == f])) == 2L && length(unique(y[fold !=
        f])) == 2L, logical(1)))
      if (ok)
        return(fold)
    }
    stop("\u65e0\u6cd5\u6784\u5efa\u4e24\u7c7b\u9f50\u5168\u7684\u60a3\u8005\u6298\uff1b\u4e0d\u5f97\u56de\u9000\u6837\u672c\u7ea7CV")
  }
  fixed_auc <- function(y, p) {
    if (length(p) != length(y) || any(!is.finite(p)) || length(unique(y)) != 2L)
      return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)))
  }
  patient_balanced_auc <- function(y, p, subject) {
    if (length(y) != length(p) || length(y) != length(subject) || any(!is.finite(p)) || anyNA(subject))
      return(NA_real_)
    d <- data.frame(subject = as.character(subject), Type = as.integer(y), probability = as.numeric(p), stringsAsFactors = FALSE)
    a <- aggregate(probability ~ subject + Type, data = d, FUN = mean)
    if (length(unique(a$Type)) != 2L)
      return(NA_real_)
    fixed_auc(a$Type, a$probability)
  }
  logloss <- function(y, p) {
    p <- pmin(pmax(as.numeric(p), 1e-8), 1 - 1e-8)
    if (length(p) != length(y) || any(!is.finite(p)))
      stop("CV\u6982\u7387\u65e0\u6548")
    mean(-y * log(p) - (1 - y) * log(1 - p))
  }
  current_subjects <- NULL
  inner_foldid <- function(x, y) make_group_folds(unname(current_subjects[rownames(x)]), y, inner_k)
  OriginalExtractVar <- ExtractVar
  ExtractVar <- function(fit) {
    if (!is.null(fit$selectedFeature))
      return(fit$selectedFeature)
    f <- switch(class(fit)[1], lda = fit$subFeature, glm = names(coef(fit))[is.finite(coef(fit))], gbm = {
      z <- summary(fit, plotit = FALSE, n.trees = fit$best_trees)
      as.character(z$var[z$rel.inf > 0])
    }, OriginalExtractVar(fit))
    intersect(unique(setdiff(f, c("(Intercept)", "Intercept"))), fit$subFeature)
  }
  CalPredictScore <- function(fit, new_data, type = "response") {
    x <- new_data[, fit$subFeature, drop = FALSE]
    p <- switch(class(fit)[1], lognet = predict(fit, newx = as.matrix(x), type = "response"), glm = predict(fit,
      newdata = as.data.frame(x), type = "response"), svm.formula = {
      pp <- attr(predict(fit, as.data.frame(x), probability = TRUE), "probabilities")
      if (is.null(pp) || !"1" %in% colnames(pp)) stop("SVM\u672a\u8fd4\u56de\u7c7b\u522b1\u6982\u7387")
      pp[, "1"]
    }, lda = predict(fit, as.data.frame(x))$posterior[, "1"], glmboost = predict(fit, newdata = as.data.frame(x),
      type = "response"), plsRglmmodel = predict(fit, newdata = as.data.frame(x), type = "response"), rfsrc = predict(fit,
      as.data.frame(x))$predicted[, "1"], gbm = predict(fit, newdata = as.data.frame(x), n.trees = fit$best_trees,
      type = "response"), xgb.Booster = predict(fit, as.matrix(x)), naiveBayes = predict(fit, newdata = as.data.frame(x),
      type = "raw")[, "1"], stop("\u672a\u652f\u6301\u6a21\u578b\u7c7b\u522b\uff1a", class(fit)[1]))
    p <- as.numeric(p)
    if (length(p) != nrow(x) || any(!is.finite(p)) || any(p < -1e-8 | p > 1 + 1e-8))
      stop("\u6982\u7387\u65e0\u6548")
    setNames(pmin(pmax(p, 0), 1), rownames(x))
  }
  PredictClass <- function(fit, new_data) ifelse(CalPredictScore(fit, new_data) >= 0.5, "1", "0")
  RunEnet <- function(Train_set, Train_label, mode, classVar, alpha) {
    y <- Train_label[[classVar]]
    fold <- inner_foldid(Train_set, y)
    cvfit <- glmnet::cv.glmnet(Train_set, y, family = "binomial", alpha = alpha, foldid = fold, type.measure = "deviance")
    fit <- glmnet::glmnet(Train_set, y, family = "binomial", alpha = alpha, lambda = cvfit$lambda.min)
    fit$tuning <- list(lambda = cvfit$lambda.min, folds = max(fold), alpha = alpha)
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  RunSVM <- function(Train_set, Train_label, mode, classVar) {
    d <- as.data.frame(Train_set)
    d[[classVar]] <- factor(Train_label[[classVar]], levels = c(0, 1))
    fit <- e1071::svm(reformulate(colnames(Train_set), classVar), data = d, probability = TRUE, scale = FALSE)
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  RunLDA <- function(Train_set, Train_label, mode, classVar) {
    d <- as.data.frame(Train_set)
    d[[classVar]] <- factor(Train_label[[classVar]], levels = c(0, 1))
    fit <- MASS::lda(reformulate(colnames(Train_set), classVar), data = d)
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  RunglmBoost <- function(Train_set, Train_label, mode, classVar) {
    d <- as.data.frame(Train_set)
    d[[classVar]] <- factor(Train_label[[classVar]], levels = c(0, 1))
    form <- reformulate(colnames(Train_set), classVar)
    fit <- mboost::glmboost(form, data = d, family = mboost::Binomial())
    fold <- inner_foldid(Train_set, Train_label[[classVar]])
    weights <- model.weights(fit)
    folds <- vapply(seq_len(max(fold)), function(f) weights * (fold != f), numeric(nrow(d)))
    cvfit <- mboost::cvrisk(fit, folds = folds, papply = lapply)
    best <- mboost::mstop(cvfit)
    fit <- mboost::glmboost(form, data = d, family = mboost::Binomial(), control = mboost::boost_control(mstop = best))
    fit$tuning <- list(mstop = best, folds = max(fold))
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  RunplsRglm <- function(Train_set, Train_label, mode, classVar) {
    grid <- seq_len(min(5L, ncol(Train_set)))
    builder <- function(x, y, nt) plsRglm::plsRglm(y, as.data.frame(x), nt = nt, modele = "pls-glm-logistic",
      verbose = FALSE, sparse = TRUE)
    fold <- inner_foldid(Train_set, Train_label[[classVar]])
    losses <- matrix(NA_real_, max(fold), length(grid))
    for (f in seq_len(max(fold))) {
      tr <- fold != f
      va <- !tr
      s <- fit_scaler(Train_set[tr, , drop = FALSE])
      for (j in seq_along(grid)) {
        fit <- builder(apply_scaler(Train_set[tr, , drop = FALSE], s), Train_label[[classVar]][tr], grid[j])
        p <- predict(fit, newdata = as.data.frame(apply_scaler(Train_set[va, , drop = FALSE], s)), type = "response")
        losses[f, j] <- logloss(Train_label[[classVar]][va], p)
      }
    }
    mean_loss <- colSums(losses * tabulate(fold))/length(fold)
    best <- grid[which.min(mean_loss)]
    fit <- builder(Train_set, Train_label[[classVar]], best)
    fit$tuning <- list(nt = best, grid = grid, logloss = mean_loss, folds = max(fold))
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  RunGBM <- function(Train_set, Train_label, mode, classVar) {
    builder <- function(x, y) {
      d <- as.data.frame(x)
      d$.outcome <- y
      gbm::gbm(reformulate(colnames(x), ".outcome"), data = d, distribution = "bernoulli", n.trees = gbm_max_trees,
        interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 0, n.cores = 1, verbose = FALSE)
    }
    fold <- inner_foldid(Train_set, Train_label[[classVar]])
    losses <- matrix(NA_real_, max(fold), length(gbm_tree_grid))
    for (f in seq_len(max(fold))) {
      tr <- fold != f
      va <- !tr
      s <- fit_scaler(Train_set[tr, , drop = FALSE])
      fit <- builder(apply_scaler(Train_set[tr, , drop = FALSE], s), Train_label[[classVar]][tr])
      pp <- predict(fit, newdata = as.data.frame(apply_scaler(Train_set[va, , drop = FALSE], s)), n.trees = gbm_tree_grid,
        type = "response")
      losses[f, ] <- apply(as.matrix(pp), 2, function(p) logloss(Train_label[[classVar]][va], p))
    }
    mean_loss <- colSums(losses * tabulate(fold))/length(fold)
    fit <- builder(Train_set, Train_label[[classVar]])
    fit$best_trees <- gbm_tree_grid[which.min(mean_loss)]
    fit$tuning <- list(n.trees = fit$best_trees, logloss = mean_loss, folds = max(fold))
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  RunXGBoost <- function(Train_set, Train_label, mode, classVar) {
    params <- list(max_depth = 2, eta = 1, nthread = 1, objective = "binary:logistic", eval_metric = "logloss")
    fold <- inner_foldid(Train_set, Train_label[[classVar]])
    losses <- matrix(NA_real_, max(fold), xgb_max_rounds)
    for (f in seq_len(max(fold))) {
      tr <- fold != f
      va <- !tr
      s <- fit_scaler(Train_set[tr, , drop = FALSE])
      dt <- xgboost::xgb.DMatrix(apply_scaler(Train_set[tr, , drop = FALSE], s), label = Train_label[[classVar]][tr])
      dv <- xgboost::xgb.DMatrix(apply_scaler(Train_set[va, , drop = FALSE], s), label = Train_label[[classVar]][va])
      fit <- xgboost::xgb.train(params = params, data = dt, nrounds = xgb_max_rounds, watchlist = list(validation = dv),
        verbose = 0)
      losses[f, ] <- fit$evaluation_log$validation_logloss
    }
    mean_loss <- colSums(losses * tabulate(fold))/length(fold)
    best <- which.min(mean_loss)
    fit <- xgboost::xgb.train(params = params, data = xgboost::xgb.DMatrix(Train_set, label = Train_label[[classVar]]),
      nrounds = best, verbose = 0)
    fit$tuning <- list(nrounds = best, logloss = mean_loss, folds = max(fold))
    fit$subFeature <- colnames(Train_set)
    if (mode == "Variable")
      ExtractVar(fit)
    else fit
  }
  OriginalRunML <- RunML
  RunML <- function(method, Train_set, Train_label, mode = "Model", classVar = "Type") {
    obj <- OriginalRunML(method, Train_set, Train_label, mode, classVar)
    if (mode == "Model")
      obj$selectedFeature <- ExtractVar(obj)
    obj
  }
  fit_candidates <- function(x, y, subjects, methods, stage_seed) {
    current_subjects <<- setNames(subjects, rownames(x))
    labels <- data.frame(Type = y, row.names = rownames(x))
    selectors <- unique(vapply(strsplit(methods, "\\+"), function(z) if (length(z) == 1L)
      "simple"
    else z[1], character(1)))
    selected <- list(simple = colnames(x))
    selector_errors <- list()
    for (s in setdiff(selectors, "simple")) {
      set.seed(stage_seed + match(s, selectors))
      ans <- tryCatch(RunML(s, x, labels, "Variable"), error = function(e) e)
      if (inherits(ans, "error"))
        selector_errors[[s]] <- conditionMessage(ans)
      else selected[[s]] <- ans
    }
    models <- list()
    status <- vector("list", length(methods))
    for (i in seq_along(methods)) {
      parts <- strsplit(methods[i], "\\+")[[1]]
      if (length(parts) == 1L)
        parts <- c("simple", parts)
      state <- "success"
      detail <- ""
      nfeatures <- NA_integer_
      if (length(parts) != 2L) {
        state <- "unsupported_combination"
        detail <- "\u53ea\u652f\u6301\u5355\u6a21\u578b\u6216\u9009\u62e9\u5668+\u5206\u7c7b\u5668"
      }
      else if (!is.null(selector_errors[[parts[1]]])) {
        state <- "selector_error"
        detail <- selector_errors[[parts[1]]]
      }
      else if (!length(selected[[parts[1]]]))
        state <- "empty_selector"
      else {
        set.seed(stage_seed + 1000L + i)
        ans <- tryCatch({
         fit <- RunML(parts[2], x[, selected[[parts[1]]], drop = FALSE], labels)
         nfeatures <- length(ExtractVar(fit))
         if (nfeatures < min_retained_features)
          stop("FEATURE_THRESHOLD: retained=", nfeatures)
         CalPredictScore(fit, x)
         fit
        }, error = function(e) e)
        if (inherits(ans, "error")) {
         detail <- conditionMessage(ans)
         state <- if (startsWith(detail, "FEATURE_THRESHOLD:"))
          "excluded_feature_threshold"
         else "model_error"
        }
        else models[[methods[i]]] <- ans
      }
      status[[i]] <- data.frame(Method = methods[i], status = state, retained_features = nfeatures, detail = detail)
    }
    list(models = models, status = do.call(rbind, status), selectors = selected)
  }
  read_expression <- function(path) {
    z <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"Type" %in% names(z) || anyDuplicated(rownames(z)))
      stop("\u7f3a\u5c11Type\u6216\u6837\u672c\u91cd\u590d")
    y <- suppressWarnings(as.numeric(as.character(z$Type)))
    if (anyNA(y) || !all(y %in% c(0, 1)) || length(unique(y)) != 2L)
      stop("Type\u5fc5\u987b\u5305\u542b0\u548c1")
    x <- as.matrix(z[, setdiff(names(z), "Type"), drop = FALSE])
    storage.mode(x) <- "double"
    if (!ncol(x) || anyDuplicated(colnames(x)) || any(!is.finite(x)))
      stop("\u8868\u8fbe\u77e9\u9635\u6709\u975e\u6570\u503c/NA/\u91cd\u590d\u57fa\u56e0")
    list(x = x, y = y)
  }
  infer_sample_manifest <- function(x, y, dataset) {
    raw <- rownames(x)
    ids <- sub("^Train\\.", "", raw)
    pattern <- "(?:^|_)GSM(S[0-9]+)(DX|day[0-9]+|week[0-9]+)_(NonEndpoint|Endpoint)$"
    matches <- regmatches(ids, regexec(pattern, ids, perl = TRUE, ignore.case = TRUE))
    good <- lengths(matches) == 4L
    if (!all(good))
      stop("\u65e0\u6cd5\u89e3\u6790\u6837\u672cID\uff1a", paste(head(ids[!good], 8), collapse = "; "), "\u3002\u9700\u8981\u7c7b\u4f3c GSE89403part1_GSMS147week4_NonEndpoint \u6216 week24_Endpoint \u7684\u540d\u79f0")
    subject <- toupper(vapply(matches, `[`, character(1), 2L))
    timepoint <- tolower(vapply(matches, `[`, character(1), 3L))
    group <- vapply(matches, `[`, character(1), 4L)
    expected <- ifelse(tolower(group) == "endpoint", 1L, 0L)
    if (any(expected != y))
      stop("\u6837\u672c\u540dNonEndpoint/Endpoint\u4e0eType\u4e0d\u4e00\u81f4\uff1a", paste(head(ids[expected != y], 8), collapse = "; "))
    number <- suppressWarnings(as.numeric(sub("^(day|week)", "", timepoint)))
    week <- ifelse(timepoint == "dx", 0, ifelse(startsWith(timepoint, "day"), number/7, number))
    data.frame(sample_id = ids, subject_id = subject, timepoint = timepoint, sample_week = week, group = ifelse(expected ==
      1, "Endpoint", "NonEndpoint"), dataset = dataset, Type = y, output_sample_id = raw, stringsAsFactors = FALSE)
  }
  run_analysis <- function() {
    inputs <- file.path(input_dir, c("data.train.txt", "data.test.txt"))
    train <- read_expression(inputs[1])
    test <- read_expression(inputs[2])
    tm <- infer_sample_manifest(train$x, train$y, "Training")
    vm <- infer_sample_manifest(test$x, test$y, "Internal_test")
    if (anyDuplicated(paste(tm$subject_id, tm$timepoint)) || anyDuplicated(paste(vm$subject_id, vm$timepoint)))
      stop("\u540c\u4e00\u60a3\u8005\u540c\u4e00\u65f6\u95f4\u70b9\u91cd\u590d\uff0c\u8bf7\u6838\u5b9e\u662f\u5426\u91cd\u590d\u7eb3\u5165")
    if (length(intersect(tm$subject_id, vm$subject_id)) || length(intersect(tm$sample_id, vm$sample_id)))
      stop("\u8bad\u7ec3/\u6d4b\u8bd5\u60a3\u8005\u6216\u6837\u672c\u91cd\u53e0")
    all_manifest <- rbind(tm, vm)
    subject_class_check <- aggregate(Type ~ subject_id, data = all_manifest, FUN = function(z) paste(sort(unique(z)),
      collapse = ""))
    incomplete_subjects <- subject_class_check$subject_id[subject_class_check$Type != "01"]
    if (length(incomplete_subjects))
      stop("\u4ee5\u4e0b\u60a3\u8005\u6ca1\u6709\u540c\u65f6\u63d0\u4f9bNonEndpoint\u548cEndpoint\u6837\u672c\uff1a", paste(incomplete_subjects, collapse = ", "))
    if (any(tm$timepoint == "week24" & tm$Type != 1L) || any(vm$timepoint == "week24" & vm$Type != 1L) || any(tm$timepoint !=
      "week24" & tm$Type != 0L) || any(vm$timepoint != "week24" & vm$Type != 0L)) {
      stop("Type\u4e0e\u6cbb\u7597\u9636\u6bb5\u4e0d\u4e00\u81f4\uff1a\u4ec5week24\u5141\u8bb8Endpoint=1")
    }
    if (!all(colnames(train$x) %in% colnames(test$x)))
      stop("\u6d4b\u8bd5\u96c6\u7f3a\u5c11\u8bad\u7ec3\u57fa\u56e0")
    test$x <- test$x[, colnames(train$x), drop = FALSE]
    methods <- unique(gsub("-| ", "", built_in_methods))
    if (!length(methods) || anyNA(methods) || any(!nzchar(methods)))
      stop("\u7b97\u6cd5\u5217\u8868\u65e0\u6548")
    output_dir <- file.path(output_root, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_seed", ml_seed))
    if (dir.exists(output_dir))
      stop("\u8f93\u51fa\u76ee\u5f55\u5df2\u5b58\u5728\uff0c\u4e0d\u8986\u76d6\u5386\u53f2\u7ed3\u679c")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir))
      stop("\u4e0d\u80fd\u521b\u5efa\u8f93\u51fa\u76ee\u5f55")
    out <- function(name) file.path(output_dir, name)
    write.csv(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))), out("input_checksums.csv"), row.names = FALSE)
    write.csv(rbind(tm, vm), out("sample_manifest.csv"), row.names = FALSE)
    file.copy(inputs, output_dir, overwrite = FALSE)
    write.csv(rbind(tm, vm), out("patient_level_train_test_split_manifest.csv"), row.names = FALSE)
    write.table(data.frame(Model = methods), out("refer.methodLists.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    config <- list(seed = ml_seed, outer_k = outer_k, inner_k = inner_k, min_retained_features = min_retained_features,
      methods = methods, GBM_max_trees = gbm_max_trees, XGB_max_rounds = xgb_max_rounds, positive_class = "Endpoint=week24/M6=1",
      negative_class = "NonEndpoint=DX/day7/week4=0", target = "treatment stage, not sputum conversion or cure outcome")
    saveRDS(config, out("run_config.rds"))
    fold <- make_group_folds(tm$subject_id, train$y, outer_k)
    write.csv(data.frame(sample_id = rownames(train$x), subject_id = tm$subject_id, Type = train$y, outer_fold = fold),
      out("patient_CV_folds.csv"), row.names = FALSE)
    oof <- matrix(NA_real_, nrow(train$x), length(methods), dimnames = list(rownames(train$x), methods))
    statuses <- list()
    for (f in seq_len(max(fold))) {
      message("\u60a3\u8005\u7ea7\u5916\u5c42CV ", f, "/", max(fold), "\uff1b\u7ec4\u5408\u6570 ", length(methods))
      tr <- fold != f
      va <- !tr
      s <- fit_scaler(train$x[tr, , drop = FALSE])
      xx <- apply_scaler(train$x[tr, , drop = FALSE], s)
      vv <- apply_scaler(train$x[va, , drop = FALSE], s)
      fits <- fit_candidates(xx, train$y[tr], tm$subject_id[tr], methods, ml_seed + f * 10000L)
      st <- fits$status
      st$outer_fold <- f
      for (method in names(fits$models)) {
        p <- tryCatch(CalPredictScore(fits$models[[method]], vv), error = function(e) e)
        if (inherits(p, "error")) {
         st$status[st$Method == method] <- "validation_prediction_error"
         st$detail[st$Method == method] <- conditionMessage(p)
        }
        else oof[va, method] <- p
      }
      statuses[[f]] <- st
      write.csv(do.call(rbind, statuses), out("CV_model_status.csv"), row.names = FALSE)
      saveRDS(list(preprocessing = s, selectors = fits$selectors, features = lapply(fits$models, ExtractVar),
        tuning = lapply(fits$models, function(z) z$tuning)), out(paste0("CV_fold", f, "_audit.rds")))
      saveRDS(oof, out("training_OOF_predictions.rds"))
    }
    scaler <- fit_scaler(train$x)
    Train_set <- apply_scaler(train$x, scaler)
    Test_set <- apply_scaler(test$x, scaler)
    message("\u5168\u8bad\u7ec3\u96c6\u62df\u5408\u5019\u9009\u6a21\u578b\uff1b\u6d4b\u8bd5\u96c6\u4e0d\u53c2\u4e0e\u6a21\u578b\u9009\u62e9")
    full <- fit_candidates(Train_set, train$y, tm$subject_id, methods, ml_seed + 90000L)
    saveRDS(full$models, out("model.MLmodel.rds"))
    saveRDS(scaler, out("preprocessing.rds"))
    write.csv(full$status, out("model.status.csv"), row.names = FALSE)
    cvauc <- vapply(seq_along(methods), function(j) patient_balanced_auc(train$y, oof[, j], tm$subject_id),
      numeric(1))
    eligible <- methods %in% names(full$models) & is.finite(cvauc)
    ranking <- data.frame(Method = methods, Train_patient_balanced_OOF_AUC = cvauc, Train_patient_OOF_AUC = cvauc,
      complete_OOF = colSums(is.finite(oof)) == nrow(oof), eligible = eligible)
    ranking <- ranking[order(!ranking$eligible, -ranking$Train_patient_balanced_OOF_AUC, ranking$Method, na.last = TRUE),
      ]
    write.csv(ranking, out("model.training_CV_ranking.csv"), row.names = FALSE)
    write.table(data.frame(id = rownames(oof), oof, check.names = FALSE), out("model.training_OOF_riskMatrix.txt"),
      sep = "\t", row.names = FALSE, quote = FALSE)
    if (!any(eligible))
      stop("\u65e0\u5b8c\u6574OOF\u4e14\u5168\u8bad\u7ec3\u62df\u5408\u6210\u529f\u7684\u6a21\u578b\uff0c\u68c0\u67e5\uff1a", output_dir)
    reported_features <- function(z) if (inherits(z, "rfsrc"))
      z$subFeature
    else ExtractVar(z)
    features <- do.call(rbind, lapply(names(full$models), function(m) data.frame(features = reported_features(full$models[[m]]),
      algorithm = m)))
    write.table(features, out("model.genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    best <- ranking$Method[which(ranking$eligible)[1]]
    fit <- full$models[[best]]
    bundle <- list(method = best, fit = fit, preprocessing = scaler, features = reported_features(fit), required_input_features = fit$subFeature,
      selector_report = ExtractVar(fit), positive_class = 1L, threshold = 0.5, selection = "training patient-balanced OOF AUC",
      config = config)
    saveRDS(bundle, out("selected_model_bundle.rds"))
    writeLines(best, out("selected_method.txt"))
    ptrain <- CalPredictScore(fit, Train_set)
    ptest <- CalPredictScore(fit, Test_set)
    metrics <- data.frame(Method = best, Train_apparent_sample_AUC = fixed_auc(train$y, ptrain), Train_apparent_patient_balanced_AUC = patient_balanced_auc(train$y,
      ptrain, tm$subject_id), Train_patient_balanced_OOF_AUC = ranking$Train_patient_balanced_OOF_AUC[ranking$Method ==
      best], Internal_test_sample_AUC = fixed_auc(test$y, ptest), Internal_test_patient_balanced_AUC = patient_balanced_auc(test$y,
      ptest, vm$subject_id))
    set.seed(ml_seed + 999L)
    subjects <- unique(vm$subject_id)
    boot <- replicate(1000L, {
      drawn <- subjects[sample.int(length(subjects), length(subjects), replace = TRUE)]
      index_list <- lapply(drawn, function(s) which(vm$subject_id == s))
      ii <- unlist(index_list, use.names = FALSE)
      bootstrap_subject <- unlist(lapply(seq_along(index_list), function(i) rep(paste0("bootstrap_patient_",
        i), length(index_list[[i]]))), use.names = FALSE)
      patient_balanced_auc(test$y[ii], ptest[ii], bootstrap_subject)
    })
    good <- boot[is.finite(boot)]
    ci <- if (length(good) >= 900L)
      quantile(good, c(0.025, 0.975), names = FALSE)
    else c(NA_real_, NA_real_)
    metrics$Internal_test_patient_balanced_cluster_bootstrap_CI_low <- ci[1]
    metrics$Internal_test_patient_balanced_cluster_bootstrap_CI_high <- ci[2]
    metrics$valid_bootstrap_runs <- length(good)
    write.csv(metrics, out("selected_model_performance.csv"), row.names = FALSE)
    saveRDS(boot, out("test_patient_bootstrap_auc.rds"))
    write.table(metrics, out("model.AUCmatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    risk <- data.frame(id = c(rownames(Train_set), rownames(Test_set)), probability = c(ptrain, ptest))
    names(risk)[2] <- best
    write.table(risk, out("model.riskMatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    classes <- risk
    classes[[2]] <- as.integer(classes[[2]] >= 0.5)
    write.table(classes, out("model.classMatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    capture.output(sessionInfo(), file = out("sessionInfo.txt"))
    message("\u5b8c\u6210\u3002\u9009\u5b9a\u7ec4\u5408\uff1a", best, "\uff1b\u8f93\u51fa\uff1a", output_dir)
    invisible(list(output_dir = output_dir, metrics = metrics, bundle = bundle))
  }
  built_in_methods <- c("Lasso + Stepglm [both]", "SVM", "glmBoost + SVM", "Ridge", "Lasso + SVM", "glmBoost + Ridge",
    "Enet [alpha=0.1]", "glmBoost + Enet [alpha=0.1]", "Enet [alpha=0.2]", "Enet [alpha=0.3]", "glmBoost + Enet [alpha=0.3]",
    "glmBoost + Enet [alpha=0.2]", "Enet [alpha=0.4]", "glmBoost + Enet [alpha=0.4]", "Lasso + glmBoost", "Enet [alpha=0.5]",
    "glmBoost", "glmBoost + Enet [alpha=0.5]", "Enet [alpha=0.6]", "glmBoost + Enet [alpha=0.6]", "glmBoost + Enet [alpha=0.7]",
    "glmBoost + Enet [alpha=0.8]", "Enet [alpha=0.8]", "Enet [alpha=0.9]", "Lasso", "Enet [alpha=0.7]", "glmBoost + Enet [alpha=0.9]",
    "glmBoost + Lasso", "Lasso + plsRglm", "glmBoost + plsRglm", "glmBoost + Stepglm [forward]", "Lasso + Stepglm [forward]",
    "RF + SVM", "Stepglm [forward]", "plsRglm", "RF + Ridge", "RF + Enet [alpha=0.1]", "RF + plsRglm", "RF + Stepglm [forward]",
    "RF + Enet [alpha=0.2]", "RF + Enet [alpha=0.3]", "RF + Enet [alpha=0.6]", "RF + Lasso", "RF + Enet [alpha=0.7]",
    "RF + Enet [alpha=0.5]", "RF + glmBoost", "RF + Enet [alpha=0.9]", "RF + Enet [alpha=0.4]", "RF + Enet [alpha=0.8]",
    "RF + Stepglm [both]", "RF + Stepglm [backward]", "Stepglm [both] + Ridge", "Stepglm [backward] + Ridge",
    "Stepglm [both] + plsRglm", "Stepglm [backward] + plsRglm", "Stepglm [both] + Enet [alpha=0.9]", "Stepglm [backward] + Enet [alpha=0.9]",
    "Stepglm [both] + Enet [alpha=0.1]", "Stepglm [backward] + Enet [alpha=0.1]", "Stepglm [both] + Enet [alpha=0.8]",
    "Stepglm [backward] + Enet [alpha=0.8]", "Stepglm [both] + Enet [alpha=0.2]", "Stepglm [backward] + Enet [alpha=0.2]",
    "Stepglm [both] + Lasso", "Stepglm [backward] + Lasso", "Stepglm [both] + Enet [alpha=0.6]", "Stepglm [backward] + Enet [alpha=0.6]",
    "glmBoost + GBM", "Stepglm [both] + Enet [alpha=0.7]", "Stepglm [backward] + Enet [alpha=0.7]", "Lasso + Stepglm [backward]",
    "Stepglm [both]", "Stepglm [backward]", "glmBoost + Stepglm [both]", "glmBoost + Stepglm [backward]", "Stepglm [both] + Enet [alpha=0.4]",
    "Stepglm [backward] + Enet [alpha=0.4]", "Stepglm [both] + Enet [alpha=0.3]", "Stepglm [backward] + Enet [alpha=0.3]",
    "Stepglm [both] + glmBoost", "Stepglm [backward] + glmBoost", "Stepglm [both] + Enet [alpha=0.5]", "Stepglm [backward] + Enet [alpha=0.5]",
    "glmBoost + RF", "RF", "Lasso + GBM", "RF + GBM", "GBM", "Stepglm [both] + SVM", "Stepglm [backward] + SVM",
    "Lasso + RF", "Stepglm [both] + GBM", "Stepglm [backward] + GBM", "Stepglm [both] + RF", "LDA", "glmBoost + LDA",
    "RF+LDA", "Stepglm [both]+LDA", "Stepglm [backward]+LDA", "Lasso + LDA", "Stepglm [backward] + RF", "XGBoost",
    "Lasso+XGBoost", "glmBoost+XGBoost", "RF+XGBoost", "Stepglm[both]+XGBoost", "Stepglm[backward]+XGBoost",
    "NaiveBayes", "Lasso+NaiveBayes", "glmBoost+NaiveBayes", "RF+NaiveBayes", "Stepglm[both]+NaiveBayes", "Stepglm[backward]+NaiveBayes")
  if (functions_only)
    return(environment())
  run_analysis()
}

acti_ml_validate <- function(cfg = acti_ml_config()) {
  h <- acti_ml_models(cfg, functions_only = TRUE)
  train <- h$read_expression(file.path(cfg$input_dir, "data.train.txt"))
  test <- h$read_expression(file.path(cfg$input_dir, "data.test.txt"))
  tm <- h$infer_sample_manifest(train$x, train$y, "Training")
  vm <- h$infer_sample_manifest(test$x, test$y, "Internal_test")
  stopifnot(!length(intersect(tm$subject_id, vm$subject_id)), !length(intersect(tm$sample_id, vm$sample_id)),
    all(colnames(train$x) %in% colnames(test$x)))
  meta <- rbind(tm, vm)
  stopifnot(!anyDuplicated(paste(meta$subject_id, meta$timepoint)), all(meta$Type == as.integer(meta$timepoint ==
    "week24")))
  classes <- split(meta$Type, meta$subject_id)
  stopifnot(all(vapply(classes, function(y) setequal(y, c(0, 1)), logical(1))))
  fold <- h$make_group_folds(tm$subject_id, train$y, h$outer_k)
  stopifnot(all(vapply(split(fold, tm$subject_id), function(f) length(unique(f)) == 1L, logical(1))))
  cat("Prepared matrices and patient-grouped CV checks passed.\n")
  cat("Training:", length(unique(tm$subject_id)), "patients /", nrow(tm), "samples\n")
  cat("Internal test:", length(unique(vm$subject_id)), "patients /", nrow(vm), "samples\n")
  cat("Genes:", ncol(train$x), "; algorithms:", length(unique(gsub("-| ", "", h$built_in_methods))), "\n")
  invisible(list(train = tm, test = vm, fold = fold))
}

acti_ml_run <- function(stage, cfg = acti_ml_config()) {
  result <- NULL
  for (s in acti_ml_stages(stage)) {
    message("Running machine-learning stage: ", s)
    result <- switch(s, prepare = acti_ml_prepare(cfg), validate = acti_ml_validate(cfg), train = acti_ml_models(cfg))
  }
  invisible(result)
}

acti_ml_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- "Rscript 4_ML_patient_level_pipeline.R --stage NAME\nStages: all, prepare, train, validate\nUse --check to check inputs; see the configuration function for paths.\n"
  if (!length(args) || "--help" %in% args) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  flags <- c("--source-part-dir", "--pheno-file", "--gene-file", "--output-dir", "--part-dir", "--input-dir",
    "--model-dir", "--split-seed", "--stage")
  values <- list()
  stage <- "all"
  check <- FALSE
  i <- 1L
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
  cfg <- do.call(acti_ml_config, values)
  if (check) {
    if (!acti_ml_check(cfg, stage))
      stop("Preflight incomplete: supply missing inputs/packages.")
    return(invisible(cfg))
  }
  acti_ml_run(stage, cfg)
}

if (sys.nframe() == 0L) acti_ml_main()

