# Internal, external and clinical validation.
# Inputs: Processed expression, phenotype and clinical tables; acti_validation_config().
# Run: Rscript 5_Validation_and_plots.R --stage internal
# source() loads functions without starting an analysis.
# Optional ML heatmap: --ml-run-dir PATH to a completed ML_data run.

acti_validation_config <- function(data_dir = "data", output_dir = file.path("results", "Validation"), bulk_dir = NULL,
  qp_dir = NULL, pre_dir = NULL, pre_source_file = NULL, clinical_output_dir = NULL, ml_run_dir = NULL, save_plots = TRUE,
  overwrite = FALSE) {
  default <- function(x, y) if (is.null(x))
    y
  else x
  pre_dir <- default(pre_dir, file.path(data_dir, "GSE79362"))
  list(data_dir = data_dir, output_dir = output_dir, bulk_dir = default(bulk_dir, file.path(data_dir, "GSE89403")),
    qp_dir = default(qp_dir, file.path(data_dir, "clinical")), pre_dir = pre_dir, pre_source_file = default(pre_source_file,
      file.path(pre_dir, "Figure6_scatter_source_data.csv")), clinical_dir = default(clinical_output_dir,
      file.path(output_dir, "04_clinical")), ml_run_dir = ml_run_dir, save_plots = isTRUE(save_plots))
}

acti_validation_check <- function(cfg = acti_validation_config(), stage = "all") {
  requirements <- list(internal = file.path(cfg$bulk_dir, c("exp2_Definite_Cure.csv", "exp2.csv", "GSE89403_matrix_pheno.csv")),
    external = c(file.path(cfg$data_dir, "GSE40553", c("exp2all.csv", "phenoall.xlsx")), file.path(cfg$data_dir,
      "GSE54992", c("exp2.csv", "pheno.xlsx")), file.path(cfg$data_dir, "GSE31348", c("exp2.csv", "pheno.xlsx"))),
    prediagnosis = c(cfg$pre_source_file, file.path(cfg$pre_dir, "pheno_time.xlsx")), clinical_expression = file.path(cfg$qp_dir,
      c("qpcr_log2_expr_noshifted_merged_allGenes.csv", "qpcr_log2_expr_noShift_20genes.csv", "pheno.xlsx",
        "temp2.xls", "qpcr_20genes_expr_matched_lingchaung2.csv")), clinical_roc = file.path(cfg$qp_dir,
      c("qpcr_log2_expr_noshifted_merged_allGenes.csv", "qpcr_log2_expr_noShift_20genes.csv", "pheno.xlsx")),
    clinical_time = file.path(cfg$clinical_dir, "Scatter_20genes_each_pdf", "Figure7_scatter_source_data.csv"),
    model_matrix = character(), legends = character(), ml_heatmap = if (is.null(cfg$ml_run_dir)) character() else file.path(cfg$ml_run_dir,
      c("data.train.txt", "data.test.txt", "model.training_CV_ranking.csv", "model.training_OOF_riskMatrix.txt",
        "preprocessing.rds", "model.MLmodel.rds", "selected_method.txt")))
  selected <- if (stage == "all")
    names(requirements)
  else stage
  if (any(!selected %in% names(requirements)))
    stop("Unknown stage: ", stage)
  if (stage == "all")
    selected <- setdiff(selected, "clinical_time")
  files <- unique(unlist(requirements[selected], use.names = FALSE))
  print(data.frame(path = files, exists = file.exists(files)), row.names = FALSE)
  if (!all(file.exists(files)))
    stop("Missing required inputs.")
  invisible(TRUE)
}

acti_validation_internal <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(pROC)
    library(patchwork)
  })
  data_dir <- cfg$bulk_dir
  out_dir <- file.path(cfg$output_dir, "01_internal")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(20261013)
  bootstrap_runs <- 1000L
  signature_specs <- list(`ACTI-1` = list(up = c("FCGR1A", "DUSP3", "NPC2", "ANKRD22", "TAP1"), down = character(0),
    origin = "iter_000009; RF+NaiveBayes; main-text model"), `ACTI-2` = list(up = c("DUSP3", "FCGR1A", "TAP1"),
    down = character(0), origin = "iter_000003; RF+Stepglm[backward]"), `ACTI-3` = list(up = c("DUSP3", "FCGR1A",
    "NPC2"), down = character(0), origin = "iter_000007; RF+glmBoost"), `ACTI-4` = list(up = c("DUSP3", "FCGR1A",
    "NPC2", "VAMP5"), down = character(0), origin = "iter_000008; RF+LDA"), `ACTI-5` = list(up = c("FCGR1A",
    "GBP5", "DUSP3"), down = character(0), origin = "iter_000010; Stepglm[backward]+glmBoost"), RISK6 = list(up = c("GBP2",
    "FCGR1B", "SERPING1"), down = c("TUBGCP6", "TRMT2A", "SDR39U1"), origin = "published comparator"), Sweeney3 = list(up = c("GBP5",
    "DUSP3"), down = c("KLF2"), origin = "published comparator"))
  signature_order <- names(signature_specs)
  plot_signature_order <- c("ACTI-1", "ACTI-2", "ACTI-3", "ACTI-4", "ACTI-5", "RISK6", "Sweeney3")
  signature_colors <- c(`ACTI-1` = "#F8766D", `ACTI-2` = "#E69F00", `ACTI-3` = "#0072B2", `ACTI-4` = "#56B4E9",
    `ACTI-5` = "#CC79A7", RISK6 = "#00BA38", Sweeney3 = "#619CFF")
  read_expression <- function(path) {
    x <- read.csv(path, row.names = 1, check.names = FALSE)
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    x
  }
  exp_cure <- read_expression(file.path(data_dir, "exp2_Definite_Cure.csv"))
  exp_all <- read_expression(file.path(data_dir, "exp2.csv"))
  pheno_raw <- read.csv(file.path(data_dir, "GSE89403_matrix_pheno.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  score_signature <- function(exp_mat, sample_ids, spec) {
    up <- intersect(spec$up, rownames(exp_mat))
    down <- intersect(spec$down, rownames(exp_mat))
    if (length(up) == 0L)
      return(rep(NA_real_, length(sample_ids)))
    up_score <- colMeans(exp_mat[up, sample_ids, drop = FALSE], na.rm = TRUE)
    if (length(down) == 0L)
      return(as.numeric(up_score))
    down_score <- colMeans(exp_mat[down, sample_ids, drop = FALSE], na.rm = TRUE)
    as.numeric(up_score - down_score)
  }
  gene_match <- bind_rows(lapply(signature_order, function(nm) {
    spec <- signature_specs[[nm]]
    genes <- c(spec$up, spec$down)
    data.frame(Signature = nm, Direction = c(rep("Up", length(spec$up)), rep("Down", length(spec$down))), Gene = genes,
      Found_in_DefiniteCure = genes %in% rownames(exp_cure), Found_in_AllPatients = genes %in% rownames(exp_all),
      stringsAsFactors = FALSE)
  }))
  write.csv(gene_match, file.path(out_dir, "ACTI_gene_match_audit.csv"), row.names = FALSE)
  if (any(!gene_match$Found_in_DefiniteCure) || any(!gene_match$Found_in_AllPatients)) {
    stop("One or more required genes are missing; see ACTI_gene_match_audit.csv")
  }
  signature_definition <- bind_rows(lapply(signature_order, function(nm) {
    spec <- signature_specs[[nm]]
    data.frame(Signature = nm, N_genes = length(unique(c(spec$up, spec$down))), Up_genes = paste(spec$up, collapse = "; "),
      Down_genes = paste(spec$down, collapse = "; "), Score = "mean(up genes) - mean(down genes)", Origin = spec$origin,
      stringsAsFactors = FALSE)
  }))
  write.csv(signature_definition, file.path(out_dir, "ACTI_signature_definition.csv"), row.names = FALSE)
  sample_manifest <- data.frame(sample_id = colnames(exp_cure), subject_id = sub("_(DX|day_7|week_4|week_24)$",
    "", colnames(exp_cure)), time = case_when(grepl("_DX$", colnames(exp_cure)) ~ "DX", grepl("_day_7$", colnames(exp_cure)) ~
    "day7", grepl("_week_4$", colnames(exp_cure)) ~ "week4", grepl("_week_24$", colnames(exp_cure)) ~ "week24",
    TRUE ~ NA_character_), stringsAsFactors = FALSE) %>% filter(!is.na(time))
  comparison_specs <- list(C = list(early = "DX", later = "day7", title = "Baseline vs Day 7"), D = list(early = "DX",
    later = "week4", title = "Baseline vs Week 4"), E = list(early = "DX", later = "week24", title = "Baseline vs Week 24"),
    F = list(early = "day7", later = "week24", title = "Day 7 vs Week 24"), G = list(early = "week4", later = "week24",
      title = "Week 4 vs Week 24"))
  make_paired_manifest <- function(spec, panel) {
    early_ids <- sample_manifest$subject_id[sample_manifest$time == spec$early]
    later_ids <- sample_manifest$subject_id[sample_manifest$time == spec$later]
    keep_subjects <- intersect(early_ids, later_ids)
    out <- sample_manifest %>% filter(subject_id %in% keep_subjects, time %in% c(spec$early, spec$later)) %>%
      mutate(Panel = panel, Comparison = spec$title, Group = ifelse(time == spec$early, "positive", "negative"))
    stopifnot(all(table(out$subject_id) == 2L))
    out
  }
  paired_manifests <- bind_rows(lapply(names(comparison_specs), function(panel) {
    make_paired_manifest(comparison_specs[[panel]], panel)
  }))
  write.csv(paired_manifests, file.path(out_dir, "ACTI_C_to_G_complete_pair_manifest.csv"), row.names = FALSE)
  calc_roc <- function(group, score) {
    ok <- is.finite(score) & !is.na(group)
    roc(response = factor(group[ok], levels = c("negative", "positive")), predictor = score[ok], levels = c("negative",
      "positive"), direction = "<", quiet = TRUE)
  }
  cluster_bootstrap_auc <- function(dat, B = 1000L, seed = 1L) {
    fast_auc <- function(group, score) {
      ok <- is.finite(score) & !is.na(group)
      group <- group[ok]
      score <- score[ok]
      is_positive <- group == "positive"
      n_positive <- sum(is_positive)
      n_negative <- sum(!is_positive)
      if (n_positive == 0L || n_negative == 0L)
        return(NA_real_)
      ranks <- rank(score, ties.method = "average")
      (sum(ranks[is_positive]) - n_positive * (n_positive + 1)/2)/(n_positive * n_negative)
    }
    set.seed(seed)
    subjects <- unique(dat$subject_id)
    row_indices <- split(seq_len(nrow(dat)), dat$subject_id)
    vals <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      sampled <- sample(subjects, length(subjects), replace = TRUE)
      idx <- unlist(row_indices[sampled], use.names = FALSE)
      vals[b] <- fast_auc(dat$Group[idx], dat$Score[idx])
    }
    valid <- vals[is.finite(vals)]
    c(low = unname(quantile(valid, 0.025, na.rm = TRUE)), high = unname(quantile(valid, 0.975,
      na.rm = TRUE)), valid_runs = length(valid))
  }
  roc_rows <- list()
  auc_rows <- list()
  score_rows <- list()
  row_i <- 0L
  analyze_dataset <- function(dat, exp_mat, panel, comparison, paired, seed_offset = 0L) {
    result_auc <- list()
    result_roc <- list()
    result_scores <- list()
    for (j in seq_along(signature_order)) {
      nm <- signature_order[j]
      z <- dat
      z$Score <- score_signature(exp_mat, z$sample_id, signature_specs[[nm]])
      r <- calc_roc(z$Group, z$Score)
      ci <- cluster_bootstrap_auc(z, B = bootstrap_runs, seed = 20261013L + seed_offset + j)
      paired_p <- NA_real_
      median_delta <- NA_real_
      if (paired) {
        wide <- z %>% select(subject_id, Group, Score) %>% tidyr::pivot_wider(names_from = Group, values_from = Score)
        median_delta <- median(wide$positive - wide$negative, na.rm = TRUE)
        paired_p <- suppressWarnings(wilcox.test(wide$positive, wide$negative, paired = TRUE, exact = FALSE)$p.value)
      }
      result_roc[[j]] <- data.frame(Panel = panel, Comparison = comparison, Signature = nm, FPR = 1 - rev(r$specificities),
        TPR = rev(r$sensitivities), stringsAsFactors = FALSE)
      result_auc[[j]] <- data.frame(Panel = panel, Comparison = comparison, Signature = nm, AUC = as.numeric(auc(r)),
        Cluster_bootstrap_CI_low = ci["low"], Cluster_bootstrap_CI_high = ci["high"], Valid_bootstrap_runs = as.integer(ci["valid_runs"]),
        N_patients = length(unique(z$subject_id)), N_samples = nrow(z), Complete_pair_analysis = paired,
        Median_early_minus_late = median_delta, Paired_Wilcoxon_P = paired_p, stringsAsFactors = FALSE)
      result_scores[[j]] <- data.frame(Panel = panel, Comparison = comparison, Signature = nm, subject_id = z$subject_id,
        sample_id = z$sample_id, Group = z$Group, Score = z$Score, stringsAsFactors = FALSE)
    }
    list(auc = bind_rows(result_auc), roc = bind_rows(result_roc), scores = bind_rows(result_scores))
  }
  time_results <- lapply(seq_along(comparison_specs), function(i) {
    panel <- names(comparison_specs)[i]
    spec <- comparison_specs[[i]]
    dat <- paired_manifests %>% filter(Panel == panel)
    analyze_dataset(dat, exp_cure, panel, spec$title, paired = TRUE, seed_offset = 100L * i)
  })
  parse_time_week <- function(x) {
    x <- as.character(x)
    out <- rep(NA_real_, length(x))
    out[toupper(x) == "DX"] <- 0
    out[grepl("day_7", x, ignore.case = TRUE)] <- 1
    idx <- grepl("week[_ ]?[0-9]+", x, ignore.case = TRUE)
    out[idx] <- as.numeric(sub(".*week[_ ]?([0-9]+).*", "\\1", x[idx], ignore.case = TRUE))
    out
  }
  parse_conversion_week <- function(x) {
    x <- as.character(x)
    out <- rep(NA_real_, length(x))
    idx <- grepl("week[_ ]?[0-9]+", x, ignore.case = TRUE)
    out[idx] <- as.numeric(sub(".*week[_ ]?([0-9]+).*", "\\1", x[idx], ignore.case = TRUE))
    out
  }
  pheno_one <- pheno_raw %>% filter(sample_code %in% colnames(exp_all)) %>% group_by(sample_code) %>% summarise(subject_id = sub("_.*",
    "", first(sample_code)), time = first(time), timetonegativity = first(timetonegativity), treatmentresult = first(treatmentresult),
    .groups = "drop") %>% mutate(time_week = parse_time_week(time), conversion_week = parse_conversion_week(timetonegativity),
    Group = case_when(is.na(time_week) | is.na(conversion_week) ~ NA_character_, time_week < conversion_week ~
      "positive", time_week >= conversion_week ~ "negative"), sample_id = sample_code) %>% filter(!is.na(Group))
  write.csv(pheno_one, file.path(out_dir, "ACTI_H_culture_all_eligible_manifest.csv"), row.names = FALSE)
  culture_primary <- analyze_dataset(pheno_one, exp_all, panel = "H", comparison = "Culture positive vs negative",
    paired = FALSE, seed_offset = 600L)
  culture_cure_dat <- pheno_one %>% filter(grepl("Definite Cure", treatmentresult, ignore.case = TRUE))
  culture_sensitivity <- analyze_dataset(culture_cure_dat, exp_all, panel = "H_sensitivity", comparison = "Culture positive vs negative (Definite Cure only)",
    paired = FALSE, seed_offset = 700L)
  auc_summary <- bind_rows(lapply(time_results, `[[`, "auc"), list(culture_primary$auc)) %>% group_by(Panel) %>%
    mutate(Paired_Wilcoxon_FDR = ifelse(Complete_pair_analysis, p.adjust(Paired_Wilcoxon_P, method = "BH"),
      NA_real_)) %>% ungroup()
  roc_source <- bind_rows(lapply(time_results, `[[`, "roc"), list(culture_primary$roc))
  score_source <- bind_rows(lapply(time_results, `[[`, "scores"), list(culture_primary$scores))
  write.csv(auc_summary, file.path(out_dir, "ACTI_C_to_H_AUC_summary.csv"), row.names = FALSE)
  write.csv(roc_source, file.path(out_dir, "ACTI_C_to_H_ROC_curve_source_data.csv"), row.names = FALSE)
  write.csv(score_source, file.path(out_dir, "ACTI_C_to_H_score_source_data.csv"), row.names = FALSE)
  write.csv(culture_sensitivity$auc, file.path(out_dir, "ACTI_H_DefiniteCure_sensitivity_AUC.csv"), row.names = FALSE)
  make_roc_plot <- function(panel) {
    curve <- roc_source %>% filter(Panel == panel, Signature %in% plot_signature_order) %>% mutate(Signature = factor(Signature,
      levels = plot_signature_order))
    stats <- auc_summary %>% filter(Panel == panel, Signature %in% plot_signature_order) %>% mutate(Signature = factor(Signature,
      levels = plot_signature_order)) %>% arrange(Signature)
    y_offset <- 0.05
    y_step <- 0.07
    labels <- stats %>% mutate(label = paste0(as.character(Signature), ": AUC=", sprintf("%.3f", AUC)), x = 0.95,
      y = y_offset + (seq_len(n()) - 1) * y_step)
    ggplot(curve, aes(x = FPR, y = TPR, colour = Signature)) + geom_line(size = 1.35) + geom_abline(slope = 1,
      intercept = 0, linetype = "dashed", colour = "grey75") + coord_fixed(xlim = c(-0.05,
      1.05), ylim = c(-0.05, 1.05), clip = "off") + scale_x_continuous(breaks = seq(0, 1,
      0.2), limits = c(-0.05, 1.05), expand = c(0, 0)) + scale_y_continuous(breaks = seq(0,
      1, 0.2), limits = c(-0.05, 1.05), expand = c(0, 0)) + scale_colour_manual(values = signature_colors[plot_signature_order],
      breaks = plot_signature_order, name = NULL) + theme_bw(base_size = 18) + labs(x = "False Positive Rate",
      y = "True Positive Rate", colour = NULL) + theme(legend.position = "bottom", panel.border = element_rect(colour = "grey60",
      fill = NA, linewidth = 1), axis.line = element_blank(), axis.ticks = element_line(color = "grey50"),
      plot.title = element_blank(), plot.subtitle = element_blank(), plot.margin = margin(10, 10, 10, 10)) +
      geom_text(data = labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1, vjust = 0,
        size = 5, colour = "black")
  }
  plots <- setNames(lapply(c("C", "D", "E", "F", "G", "H"), make_roc_plot), c("C", "D", "E", "F", "G", "H"))
  for (panel in names(plots)) {
    safe_title <- gsub("[^A-Za-z0-9]+", "_", unique(auc_summary$Comparison[auc_summary$Panel == panel]))
    prefix <- file.path(out_dir, paste0("Panel_", panel, "_", safe_title))
    ggsave(paste0(prefix, ".pdf"), plots[[panel]], width = 6.5, height = 6.5, units = "in")
    ggsave(paste0(prefix, ".png"), plots[[panel]], width = 6.5, height = 6.5, units = "in", dpi = 600, bg = "white")
  }
  combined <- (plots$C | plots$D)/(plots$E | plots$F)/(plots$G | plots$H) + plot_annotation(tag_levels = list(c("C",
    "D", "E", "F", "G", "H"))) & theme(plot.tag = element_text(size = 14, face = "bold"))
  ggsave(file.path(out_dir, "Figure5_C_to_H_ACTI_ROC_combined.pdf"), combined, width = 240, height = 330, units = "mm",
    device = cairo_pdf)
  ggsave(file.path(out_dir, "Figure5_C_to_H_ACTI_ROC_combined.png"), combined, width = 240, height = 330, units = "mm",
    dpi = 400, bg = "white")
  message("Completed: ", out_dir)
}

acti_validation_external <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readxl)
    library(pROC)
  })
  out_dir <- file.path(cfg$output_dir, "02_external")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  signature_order <- c("ACTI-1", "ACTI-2", "ACTI-3", "ACTI-4", "ACTI-5", "RISK6", "Sweeney3")
  signature_colors <- c(`ACTI-1` = "#F8766D", `ACTI-2` = "#E69F00", `ACTI-3` = "#0072B2", `ACTI-4` = "#56B4E9",
    `ACTI-5` = "#CC79A7", RISK6 = "#00BA38", Sweeney3 = "#619CFF")
  signature_specs <- list(`ACTI-1` = list(up = c("FCGR1A", "DUSP3", "NPC2", "ANKRD22", "TAP1"), down = character(0)),
    `ACTI-2` = list(up = c("DUSP3", "FCGR1A", "TAP1"), down = character(0)), `ACTI-3` = list(up = c("DUSP3",
      "FCGR1A", "NPC2"), down = character(0)), `ACTI-4` = list(up = c("DUSP3", "FCGR1A", "NPC2", "VAMP5"),
      down = character(0)), `ACTI-5` = list(up = c("FCGR1A", "GBP5", "DUSP3"), down = character(0)), RISK6 = list(up = c("GBP2",
      "FCGR1B", "SERPING1"), down = c("TUBGCP6", "TRMT2A", "SDR39U1")), Sweeney3 = list(up = c("GBP5", "DUSP3"),
      down = c("KLF2")))
  alias_registry <- list(GBP2 = c("GBP2"), FCGR1B = c("FCGR1B"), SERPING1 = c("SERPING1"), TUBGCP6 = c("TUBGCP6"),
    TRMT2A = c("TRMT2A", "HTF9C"), SDR39U1 = c("SDR39U1", "C14orf124"))
  dataset_cfgs <- list(GSE40553 = list(dir = file.path(cfg$data_dir, "GSE40553"), exp_file = "exp2all.csv", pheno_file = "phenoall.xlsx",
    id_col = "id", group_col = "Sample", comparisons = list(DX_vs_Week_2 = c("active TB pre-treatment", "active TB 2 weeks post treatment"),
      DX_vs_Month_2 = c("active TB pre-treatment", "active TB 2 months post treatment"), DX_vs_Month_6 = c("active TB pre-treatment",
        "active TB 6 months post treatment"), Week_2_vs_Month_2 = c("active TB 2 weeks post treatment",
        "active TB 2 months post treatment"), Week_2_vs_Month_6 = c("active TB 2 weeks post treatment",
        "active TB 6 months post treatment"), Month_2_vs_Month_6 = c("active TB 2 months post treatment",
        "active TB 6 months post treatment"))), GSE54992 = list(dir = file.path(cfg$data_dir, "GSE54992"),
    exp_file = "exp2.csv", pheno_file = "pheno.xlsx", id_col = "id", group_col = "Sample", comparisons = list(DX_vs_Month_3 = c("TB",
      "TB 3m"), DX_vs_Month_6 = c("TB", "TB 6m"))), GSE31348 = list(dir = file.path(cfg$data_dir, "GSE31348"),
    exp_file = "exp2.csv", pheno_file = "pheno.xlsx", id_col = "id", group_col = "Group", recode_group = function(pheno) {
      pheno$Group <- dplyr::case_when(pheno$Sample == "TB_diagnosis" ~ "TB", pheno$Sample == "TB_week1" ~
        "TB_treated_week1", pheno$Sample == "TB_week2" ~ "TB_treated_week2", pheno$Sample == "TB_week4" ~
        "TB_treated_week4", pheno$Sample == "TB_week26" ~ "TB_treated_week26", TRUE ~ NA_character_)
      dplyr::filter(pheno, !is.na(Group))
    }, comparisons = list(DX_vs_Week_1 = c("TB", "TB_treated_week1"), DX_vs_Week_2 = c("TB", "TB_treated_week2"),
      DX_vs_Week_4 = c("TB", "TB_treated_week4"), DX_vs_Week_26 = c("TB", "TB_treated_week26"))))
  display_title <- function(x) gsub("_", " ", x, fixed = TRUE)
  read_exp_numeric <- function(path) {
    dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    gene <- trimws(as.character(dat[[1]]))
    x <- as.data.frame(lapply(dat[-1], function(z) as.numeric(as.character(z))), check.names = FALSE)
    m <- as.matrix(x)
    rownames(m) <- gene
    colnames(m) <- colnames(dat)[-1]
    if (anyDuplicated(rownames(m))) {
      sums <- rowsum(m, group = rownames(m), reorder = FALSE, na.rm = TRUE)
      counts <- table(factor(gene, levels = rownames(sums)))
      m <- sums/as.numeric(counts)
    }
    m
  }
  resolve_gene <- function(gene, available) {
    aliases <- if (gene %in% names(alias_registry))
      alias_registry[[gene]]
    else gene
    idx <- match(toupper(aliases), toupper(available))
    idx <- idx[!is.na(idx)]
    if (!length(idx))
      return(NA_character_)
    available[idx[1]]
  }
  score_signature <- function(exp, ids, spec, mapping) {
    up_rows <- unname(mapping[spec$up])
    down_rows <- unname(mapping[spec$down])
    up_score <- colMeans(exp[up_rows, ids, drop = FALSE], na.rm = TRUE)
    if (length(down_rows)) {
      down_score <- colMeans(exp[down_rows, ids, drop = FALSE], na.rm = TRUE)
    }
    else {
      down_score <- setNames(rep(0, length(ids)), ids)
    }
    as.numeric(up_score - down_score)
  }
  calc_roc <- function(group, score) {
    ok <- !is.na(group) & is.finite(score)
    pROC::roc(response = factor(group[ok], levels = c("negative", "positive")), predictor = score[ok], levels = c("negative",
      "positive"), direction = "<", quiet = TRUE)
  }
  all_auc <- list()
  all_curve <- list()
  all_audit <- list()
  all_counts <- list()
  plot_store <- list()
  for (dataset_name in names(dataset_cfgs)) {
    cfg <- dataset_cfgs[[dataset_name]]
    exp <- read_exp_numeric(file.path(cfg$dir, cfg$exp_file))
    pheno <- as.data.frame(readxl::read_excel(file.path(cfg$dir, cfg$pheno_file)))
    pheno[[cfg$id_col]] <- trimws(as.character(pheno[[cfg$id_col]]))
    if (!is.null(cfg$recode_group))
      pheno <- cfg$recode_group(pheno)
    pheno <- pheno[pheno[[cfg$id_col]] %in% colnames(exp), , drop = FALSE]
    mappings <- list()
    for (sid in signature_order) {
      genes <- unique(c(signature_specs[[sid]]$up, signature_specs[[sid]]$down))
      matched <- setNames(vapply(genes, resolve_gene, character(1), available = rownames(exp)), genes)
      mappings[[sid]] <- matched
      all_audit[[paste(dataset_name, sid, sep = "|")]] <- data.frame(Dataset = dataset_name, Signature = sid,
        CanonicalGene = genes, MatchedRow = unname(matched), Found = !is.na(matched), stringsAsFactors = FALSE)
    }
    for (cmp_name in names(cfg$comparisons)) {
      lv <- cfg$comparisons[[cmp_name]]
      sub <- pheno[pheno[[cfg$group_col]] %in% lv, , drop = FALSE]
      ids <- sub[[cfg$id_col]]
      group <- ifelse(sub[[cfg$group_col]] == lv[1], "positive", "negative")
      counts <- table(factor(group, levels = c("positive", "negative")))
      all_counts[[paste(dataset_name, cmp_name, sep = "|")]] <- data.frame(Dataset = dataset_name, Comparison = cmp_name,
        Early_or_active_group = lv[1], Later_group = lv[2], N_early_or_active = unname(counts[1]), N_later = unname(counts[2]),
        N_total = length(group), stringsAsFactors = FALSE)
      curve_one <- list()
      auc_one <- list()
      for (sid in signature_order) {
        mp <- mappings[[sid]]
        if (anyNA(mp)) {
         missing <- names(mp)[is.na(mp)]
         stop(dataset_name, " / ", sid, " missing: ", paste(missing, collapse = ", "))
        }
        sc <- score_signature(exp, ids, signature_specs[[sid]], mp)
        r <- calc_roc(group, sc)
        ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
        curve_one[[sid]] <- data.frame(Dataset = dataset_name, Comparison = cmp_name, Signature = sid,
         FPR = 1 - rev(r$specificities), TPR = rev(r$sensitivities), stringsAsFactors = FALSE)
        auc_one[[sid]] <- data.frame(Dataset = dataset_name, Comparison = cmp_name, Signature = sid, AUC = as.numeric(pROC::auc(r)),
         CI_low = ci[1], CI_high = ci[3], N_early_or_active = unname(counts[1]), N_later = unname(counts[2]),
         stringsAsFactors = FALSE)
      }
      curve_df <- bind_rows(curve_one) %>% mutate(Signature = factor(Signature, levels = signature_order))
      auc_df <- bind_rows(auc_one) %>% mutate(Signature = factor(Signature, levels = signature_order)) %>%
        arrange(Signature) %>% mutate(label = paste0(as.character(Signature), ": AUC=", sprintf("%.3f",
        AUC)), x = 0.95, y = 0.05 + (seq_len(n()) - 1) * 0.07)
      p <- ggplot(curve_df, aes(FPR, TPR, colour = Signature)) + geom_line(linewidth = 1.8) + geom_abline(slope = 1,
        intercept = 0, linetype = "dashed", colour = "grey75", linewidth = 0.8) + coord_fixed(xlim = c(-0.05,
        1.05), ylim = c(-0.05, 1.05), clip = "off") + scale_x_continuous(breaks = seq(0,
        1, 0.2), limits = c(-0.05, 1.05), expand = c(0, 0)) + scale_y_continuous(breaks = seq(0,
        1, 0.2), limits = c(-0.05, 1.05), expand = c(0, 0)) + scale_colour_manual(values = signature_colors,
        breaks = signature_order, name = NULL) + theme_bw(base_size = 18) + labs(title = NULL, x = "False Positive Rate",
        y = "True Positive Rate", colour = NULL) + theme(legend.position = "none", panel.border = element_rect(colour = "grey30",
        fill = NA, linewidth = 1.2), axis.line = element_blank(), axis.ticks = element_line(colour = "grey50",
        linewidth = 0.8), plot.title = element_blank(), plot.margin = margin(10, 10, 10,
        10)) + geom_text(data = auc_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1,
        vjust = 0, size = 5, colour = "black")
      prefix <- file.path(out_dir, paste0(dataset_name, "_", cmp_name))
      ggsave(paste0(prefix, ".pdf"), p, width = 6.5, height = 6.5, units = "in")
      ggsave(paste0(prefix, ".png"), p, width = 6.5, height = 6.5, units = "in", dpi = 600, bg = "white")
      all_auc[[paste(dataset_name, cmp_name, sep = "|")]] <- auc_df %>% select(Dataset, Comparison, Signature,
        AUC, CI_low, CI_high, N_early_or_active, N_later)
      all_curve[[paste(dataset_name, cmp_name, sep = "|")]] <- curve_df
      plot_store[[paste(dataset_name, cmp_name, sep = "|")]] <- p
    }
  }
  write.csv(bind_rows(all_auc), file.path(out_dir, "Figure6_external_AUC_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(all_curve), file.path(out_dir, "Figure6_external_ROC_source_data.csv"), row.names = FALSE)
  write.csv(bind_rows(all_audit), file.path(out_dir, "Figure6_external_gene_match_audit.csv"), row.names = FALSE)
  write.csv(bind_rows(all_counts), file.path(out_dir, "Figure6_external_sample_counts.csv"), row.names = FALSE)
  message("Completed: ", out_dir)
}

acti_validation_prediagnosis <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(readxl)
  })
  base <- cfg$pre_dir
  out <- file.path(cfg$output_dir, "03_prediagnosis")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  source_file <- cfg$pre_source_file
  pheno <- as.data.frame(read_excel(file.path(base, "pheno_time.xlsx")))
  d <- read.csv(source_file, check.names = FALSE, stringsAsFactors = FALSE)
  genes <- c("DUSP3", "FCGR1A", "NPC2", "TAP1", "ANKRD22", "GBP5", "VAMP5")
  stopifnot(!anyDuplicated(pheno[["Sample Name"]]))
  m <- match(d$Sample_Name, pheno[["Sample Name"]])
  stopifnot(!anyNA(m), all(d[["Pre/Post Diagnosis"]] == pheno[["Pre/Post Diagnosis"]][m]), all(d[["Days to TB Diagnosis"]] ==
    pheno[["Days to TB Diagnosis"]][m]))
  d$Participant_ID <- pheno[["Participant ID"]][m]
  d$Relative_day <- -d[["Days to TB Diagnosis"]]
  pre <- d[d[["Pre/Post Diagnosis"]] == "Pre" & d$Gene %in% genes, ]
  stopifnot(all(pre$Relative_day < 0), !anyDuplicated(pre[c("Sample_Name", "Gene")]))
  stats <- do.call(rbind, lapply(genes, function(g) {
    a <- pre[pre$Gene == g & is.finite(pre$Expression) & is.finite(pre$Relative_day), ]
    sp <- cor.test(a$Relative_day, a$Expression, method = "spearman", exact = FALSE)
    pe <- cor.test(a$Relative_day, a$Expression, method = "pearson")
    orig <- cor.test(-a$Relative_day, a$Expression, method = "spearman", exact = FALSE)
    stopifnot(nrow(a) == 66, abs(unname(sp$estimate) + unname(orig$estimate)) < 1e-12, abs(sp$p.value -
      orig$p.value) < 1e-12)
    data.frame(Gene = g, N_observations = nrow(a), N_participants = length(unique(a$Participant_ID)), Spearman_rho = unname(sp$estimate),
      Spearman_P = sp$p.value, Pearson_r = unname(pe$estimate), Pearson_P = pe$p.value, Original_days_remaining_Spearman_rho = unname(orig$estimate))
  }))
  stats$Spearman_BH_FDR_7 <- p.adjust(stats$Spearman_P, "BH")
  stats$Pearson_BH_FDR_7 <- p.adjust(stats$Pearson_P, "BH")
  write.csv(stats, file.path(out, "Pre_diagnosis_Spearman_Pearson_results.csv"), row.names = FALSE)
  write.csv(pre, file.path(out, "Pre_diagnosis_source_data.csv"), row.names = FALSE)
  plots <- setNames(lapply(genes, function(g) {
    a <- pre[pre$Gene == g, ]
    s <- stats[stats$Gene == g, ]
    lab <- paste0("Pre-diagnosis\nrho = ", sprintf("%.3f", s$Spearman_rho), "\np = ", format(signif(s$Spearman_P,
      3), trim = TRUE), "\nn = ", s$N_observations)
    yr <- range(a$Expression, finite = TRUE)
    dy <- diff(yr)
    if (dy == 0)
      dy <- 1
    ggplot(a, aes(Relative_day, Expression)) + geom_point(colour = "#00BFC4", size = 3, alpha = 0.75) + geom_smooth(method = "loess",
      formula = y ~ x, se = TRUE, colour = "#1F78B4", fill = "#A6CEE3", alpha = 0.16, linewidth = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 1.3) + scale_x_continuous(limits = c(min(pre$Relative_day),
      0), breaks = c(-800, -600, -400, -200, 0)) + labs(x = "Time relative to TB diagnosis (days)", y = "Expression",
      title = NULL) + theme_bw(base_size = 13) + theme(plot.title = element_blank(), legend.position = "none",
      panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1.2),
      axis.ticks = element_line(linewidth = 0.8), aspect.ratio = 1, plot.margin = margin(12,
        12, 12, 12)) + annotate("label", x = min(pre$Relative_day) + 0.025 * abs(min(pre$Relative_day)),
      y = yr[2] + 0.1 * dy, label = lab, hjust = 0, vjust = 1, size = 5.5, fontface = "bold",
      label.size = 0.8, fill = "white", alpha = 0.92, lineheight = 1.05)
  }), genes)
  for (g in genes) {
    stem <- file.path(out, paste0(g, "_pre_diagnosis_time"))
    ggsave(paste0(stem, ".pdf"), plots[[g]], width = 6, height = 6, device = cairo_pdf)
    ggsave(paste0(stem, ".png"), plots[[g]], width = 6, height = 6, dpi = 600, bg = "white")
    if (requireNamespace("svglite", quietly = TRUE))
      ggsave(paste0(stem, ".svg"), plots[[g]], width = 6, height = 6, device = svglite::svglite)
    for (ext in c("pdf", "png", "svg")) {
      f <- paste0(stem, ".", ext)
      if (file.exists(f))
        file.copy(f, file.path(out, paste0(g, "_expression_pre_post_time.", ext)), overwrite = TRUE)
    }
  }
  qa <- wrap_plots(lapply(genes, function(g) plots[[g]] + ggtitle(g) + theme(plot.title = element_text(size = 16,
    face = "bold", hjust = 0.5))), ncol = 4)
  ggsave(file.path(out, "Seven_genes_pre_diagnosis_QA.png"), qa, width = 20, height = 9, dpi = 160, bg = "white")
  ggsave(file.path(out, "Seven_genes_pre_diagnosis_overview.pdf"), wrap_plots(plots, ncol = 4), width = 28, height = 12,
    device = cairo_pdf)
  print(stats, digits = 7, row.names = FALSE)
}

acti_validation_clinical_roc <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(ggplot2)
    library(pROC)
  })
  b <- cfg$qp_dir
  out <- file.path(cfg$clinical_dir, "ROC_ACTI_five_singlePDF")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  expr <- as.matrix(read.csv(file.path(b, "qpcr_log2_expr_noshifted_merged_allGenes.csv"), row.names = 1, check.names = FALSE))
  mode(expr) <- "numeric"
  colnames(expr) <- trimws(colnames(expr))
  rownames(expr) <- trimws(rownames(expr))
  original <- read.csv(file.path(b, "qpcr_log2_expr_noShift_20genes.csv"), row.names = 1, check.names = FALSE)
  ph <- as.data.frame(read_excel(file.path(b, "pheno.xlsx")))
  names(ph) <- trimws(names(ph))
  ph$id <- trimws(as.character(ph$id))
  ph$Group <- trimws(as.character(ph$Group))
  stopifnot(!anyDuplicated(ph$id))
  ph <- ph[ph$id %in% intersect(colnames(original), colnames(expr)) & ph$Group %in% c("A", "B", "C", "D"), ]
  sets <- list(`ACTI-1` = c("FCGR1A", "DUSP3", "NPC2", "ANKRD22", "TAP1"), `ACTI-2` = c("DUSP3", "FCGR1A", "TAP1"),
    `ACTI-3` = c("DUSP3", "FCGR1A", "NPC2"), `ACTI-4` = c("DUSP3", "FCGR1A", "NPC2", "VAMP5"), `ACTI-5` = c("FCGR1A",
      "GBP5", "DUSP3"))
  stopifnot(all(unique(unlist(sets)) %in% rownames(expr)))
  cols <- c(`ACTI-1` = "#F8766D", `ACTI-2` = "#E69F00", `ACTI-3` = "#0072B2", `ACTI-4` = "#56B4E9", `ACTI-5` = "#CC79A7")
  comparisons <- list(Healthy_vs_ATB = c("A", "C"), LTBI_vs_ATB = c("B", "C"), ATB_vs_Other = c("C", "D"))
  labels <- c(Healthy_vs_ATB = "Healthy vs ATB", LTBI_vs_ATB = "LTBI vs ATB", ATB_vs_Other = "ATB vs Other")
  scores <- data.frame(id = ph$id, Group = ph$Group)
  for (sid in names(sets)) scores[[sid]] <- colMeans(expr[sets[[sid]], ph$id, drop = FALSE])
  aucs <- list()
  points <- list()
  for (cmp in names(comparisons)) {
    groups <- comparisons[[cmp]]
    dat <- scores[scores$Group %in% groups, ]
    roc_rows <- list()
    anno <- list()
    for (i in seq_along(sets)) {
      sid <- names(sets)[i]
      good <- is.finite(dat[[sid]])
      dd <- dat[good, ]
      direction <- if (groups[2] == "C")
        "<"
      else ">"
      r <- roc(dd$Group, dd[[sid]], levels = groups, direction = direction, quiet = TRUE)
      ci <- as.numeric(ci.auc(r))
      a <- as.numeric(auc(r))
      aucs[[length(aucs) + 1]] <- data.frame(Comparison = cmp, Signature = sid, AUC = a, CI_low = ci[1],
        CI_high = ci[3], N_group1 = sum(dd$Group == groups[1]), N_group2 = sum(dd$Group == groups[2]),
        Direction = direction)
      rr <- data.frame(FPR = 1 - rev(r$specificities), TPR = rev(r$sensitivities), Signature = sid)
      roc_rows[[sid]] <- rr
      points[[length(points) + 1]] <- cbind(Comparison = cmp, rr)
      anno[[sid]] <- data.frame(label = paste0(sid, ": AUC=", sprintf("%.3f", a)), x = 0.95,
        y = 0.05 + (5 - i) * 0.07)
    }
    rd <- bind_rows(roc_rows)
    rd$Signature <- factor(rd$Signature, levels = names(sets))
    an <- bind_rows(anno)
    p <- ggplot(rd, aes(FPR, TPR, color = Signature)) + geom_path(linewidth = 1.8) + geom_abline(slope = 1,
      intercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.8) + coord_fixed(xlim = c(-0.05,
      1.05), ylim = c(-0.05, 1.05), clip = "off") + scale_x_continuous(breaks = seq(0, 1,
      0.2), limits = c(-0.05, 1.05), expand = c(0, 0)) + scale_y_continuous(breaks = seq(0,
      1, 0.2), limits = c(-0.05, 1.05), expand = c(0, 0)) + scale_color_manual(values = cols,
      drop = FALSE) + theme_bw(base_size = 18) + labs(title = labels[cmp], x = "False Positive Rate", y = "True Positive Rate",
      color = NULL) + theme(legend.position = "none", panel.border = element_rect(colour = "grey30", fill = NA,
      linewidth = 1.2), axis.line = element_blank(), axis.ticks = element_line(color = "grey50", linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, face = "bold"), plot.margin = margin(10, 10, 10, 10)) + geom_text(data = an,
      aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1, vjust = 0, size = 5, color = "black")
    for (ext in c("pdf", "png")) ggsave(file.path(out, paste0("ROC_ACTI_five_", cmp, ".", ext)), p, width = 6.5,
      height = 6.5, dpi = 600)
  }
  write.csv(bind_rows(aucs), file.path(out, "AUC_summary_ACTI_five.csv"), row.names = FALSE)
  write.csv(bind_rows(points), file.path(out, "ROC_source_data.csv"), row.names = FALSE)
  write.csv(scores, file.path(out, "Sample_signature_scores.csv"), row.names = FALSE)
  write.csv(data.frame(Signature = names(sets), Genes = vapply(sets, paste, character(1), collapse = ", ")),
    file.path(out, "ACTI_gene_sets.csv"), row.names = FALSE)
  cat("OUTPUT:", out, "\n")
  print(bind_rows(aucs))
}

acti_validation_clinical_time <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggtext)
  })
  b <- cfg$clinical_dir
  out <- file.path(b, "Scatter_ACTI_five_each_pdf")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- read.csv(file.path(b, "Scatter_20genes_each_pdf", "Figure7_scatter_source_data.csv"), check.names = FALSE)
  sets <- list(`ACTI-1` = c("FCGR1A", "DUSP3", "NPC2", "ANKRD22", "TAP1"), `ACTI-2` = c("DUSP3", "FCGR1A", "TAP1"),
    `ACTI-3` = c("DUSP3", "FCGR1A", "NPC2"), `ACTI-4` = c("DUSP3", "FCGR1A", "NPC2", "VAMP5"), `ACTI-5` = c("FCGR1A",
      "GBP5", "DUSP3"))
  cols <- c(`ACTI-1` = "#F8766D", `ACTI-2` = "#E69F00", `ACTI-3` = "#0072B2", `ACTI-4` = "#56B4E9", `ACTI-5` = "#CC79A7")
  stopifnot(!anyDuplicated(d[c("Gene", "Sample")]), all(unique(unlist(sets)) %in% unique(d$Gene)))
  all_scores <- list()
  all_stats <- list()
  for (sid in names(sets)) {
    dd <- d[d$Gene %in% sets[[sid]], ]
    dd <- dd %>% group_by(Sample) %>% summarise(days = first(days), n_genes = n_distinct(Gene), Score = mean(Expression),
      .groups = "drop")
    stopifnot(all(dd$n_genes == length(sets[[sid]])), all(is.finite(dd$Score)), all(is.finite(dd$days)))
    sp <- suppressWarnings(cor.test(dd$Score, dd$days, method = "spearman"))
    pe <- cor.test(dd$Score, dd$days, method = "pearson")
    sr <- data.frame(Signature = sid, Genes = paste(sets[[sid]], collapse = ", "), N = nrow(dd), Spearman_r = unname(sp$estimate),
      Spearman_p = sp$p.value, Pearson_r = unname(pe$estimate), Pearson_p = pe$p.value, Pearson_CI_low = pe$conf.int[1],
      Pearson_CI_high = pe$conf.int[2])
    all_stats[[sid]] <- sr
    all_scores[[sid]] <- cbind(Signature = sid, dd)
    a <- data.frame(label = sprintf("Spearman r = %.2f (p=%.3g)<br>Pearson r = %.2f (p=%.3g)", sr$Spearman_r,
      sr$Spearman_p, sr$Pearson_r, sr$Pearson_p))
    breaks <- sort(unique(c(pretty(range(dd$days), n = 5), 30)))
    p <- ggplot(dd, aes(days, Score)) + geom_point(size = 2.8, alpha = 0.9,
      colour = cols[sid]) + geom_vline(xintercept = 30, linetype = "dashed", linewidth = 1.3, colour = "grey40") +
      geom_smooth(aes(group = 1), method = "loess", se = FALSE, span = 1.1, linewidth = 1.5,
        colour = "black") + geom_smooth(aes(group = 1), method = "lm", se = FALSE, linewidth = 1.3, colour = cols[sid],
      linetype = "dashed") + ggtext::geom_richtext(data = a, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = 1.04, vjust = 1.04, size = 6, fill = "white", label.color = "black", label.size = 0.8,
      lineheight = 1.15) + scale_x_continuous(breaks = breaks) + theme_bw(base_size = 13) +
      theme(legend.position = "none", plot.title = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey30",
        fill = NA, linewidth = 1.2), axis.ticks = element_line(linewidth = 0.8)) + labs(title = NULL,
      x = "Treatment time (days)", y = "Mean expression (log2, -dCt)")
    for (ext in c("pdf", "png")) ggsave(file.path(out, paste0(sid, "_scatter_time.", ext)), p, width = 6, height = 6,
      dpi = 600, bg = "white")
  }
  stats <- bind_rows(all_stats)
  stats$Spearman_FDR <- p.adjust(stats$Spearman_p, "BH")
  stats$Pearson_FDR <- p.adjust(stats$Pearson_p, "BH")
  write.csv(stats, file.path(out, "ACTI_five_treatment_time_correlation_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(all_scores), file.path(out, "ACTI_five_scatter_source_data.csv"), row.names = FALSE)
  print(stats)
}

acti_validation_model_matrix <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
  })
  out_dir <- file.path(cfg$output_dir, "05_model_matrix")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  model_def <- tibble::tribble(~Model, ~Algorithm, ~Genes, "ACTI-1", "RF + NaiveBayes", "FCGR1A;DUSP3;NPC2;ANKRD22;TAP1",
    "ACTI-2", "RF + Stepglm[backward]", "DUSP3;FCGR1A;TAP1", "ACTI-3", "RF + glmBoost", "DUSP3;FCGR1A;NPC2",
    "ACTI-4", "RF + LDA", "DUSP3;FCGR1A;NPC2;VAMP5", "ACTI-5", "Stepglm[backward] + glmBoost", "FCGR1A;GBP5;DUSP3")
  long_selected <- model_def %>% tidyr::separate_rows(Genes, sep = ";") %>% rename(Gene = Genes) %>% mutate(Selected = TRUE)
  gene_frequency <- long_selected %>% count(Gene, name = "Selection_count") %>% arrange(desc(Selection_count),
    Gene)
  gene_order <- gene_frequency$Gene
  model_order <- model_def$Model
  matrix_df <- tidyr::expand_grid(Model = model_order, Gene = gene_order) %>% left_join(long_selected %>% select(Model,
    Gene, Selected), by = c("Model", "Gene")) %>% mutate(Selected = ifelse(is.na(Selected), FALSE, Selected),
    Model = factor(Model, levels = rev(model_order)), Gene = factor(Gene, levels = gene_order), Membership = factor(ifelse(Selected,
      "Selected", "Not selected"), levels = c("Selected", "Not selected")))
  model_n <- long_selected %>% count(Model, name = "n_genes")
  model_label_df <- model_def %>% left_join(model_n, by = "Model")
  model_labels <- setNames(paste0(model_label_df$Model, "  (", model_label_df$n_genes, " genes)\n", model_label_df$Algorithm),
    model_label_df$Model)
  freq_df <- gene_frequency %>% mutate(Gene = factor(Gene, levels = gene_order))
  common_theme <- theme_bw(base_size = 12) + theme(panel.grid = element_blank(), panel.border = element_rect(colour = "grey55",
    fill = NA, linewidth = 0.8), axis.text = element_text(colour = "black"), plot.margin = margin(3,
    7, 3, 4, unit = "mm"))
  p_frequency <- ggplot(freq_df, aes(x = Gene, y = Selection_count)) + geom_col(width = 0.68,
    fill = "#2F6FA3", colour = "#244F70", linewidth = 0.35) + geom_text(aes(label = paste0(Selection_count,
    "/5")), vjust = -0.35, size = 4.2, fontface = "bold", colour = "black") +
    scale_y_continuous(breaks = 0:5, limits = c(0, 6.8), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(drop = FALSE) + labs(x = NULL, y = "Runs") + common_theme + theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.title.y = element_text(size = 11, face = "bold"), axis.text.y = element_text(size = 10),
    plot.margin = margin(3, 7, 0, 4, unit = "mm"))
  p_matrix <- ggplot(matrix_df, aes(x = Gene, y = Model, fill = Membership)) + geom_tile(colour = "white", linewidth = 1.15,
    width = 0.96, height = 0.92) + scale_fill_manual(values = c(Selected = "#2F6FA3",
    `Not selected` = "#E8EDF1"), drop = FALSE, name = NULL) + scale_x_discrete(drop = FALSE, expand = expansion(add = 0.08)) +
    scale_y_discrete(drop = FALSE, labels = function(x) model_labels[as.character(x)], expand = expansion(add = c(0.12,
      0.4))) + labs(x = NULL, y = NULL) + common_theme + theme(axis.text.x = element_text(size = 10.5,
    face = "italic", angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 8.8,
    face = "bold", hjust = 1, lineheight = 0.95), axis.ticks = element_blank(), legend.position = "bottom",
    legend.direction = "horizontal", legend.text = element_text(size = 10), legend.key.width = grid::unit(7,
      "mm"), legend.key.height = grid::unit(4, "mm"), plot.margin = margin(0, 7, 4, 4, unit = "mm")) + guides(fill = guide_legend(nrow = 1,
    byrow = TRUE))
  final_plot <- p_frequency/p_matrix + plot_layout(heights = c(1.05, 2.15))
  write.csv(model_def, file.path(out_dir, "Figure5C_model_definitions.csv"), row.names = FALSE)
  write.csv(matrix_df, file.path(out_dir, "Figure5C_gene_membership_source_data.csv"), row.names = FALSE)
  write.csv(gene_frequency, file.path(out_dir, "Figure5C_gene_selection_frequency.csv"), row.names = FALSE)
  ggsave(file.path(out_dir, "Figure5C_five_ML_models_gene_matrix.pdf"), final_plot, width = 183, height = 112,
    units = "mm")
  ggsave(file.path(out_dir, "Figure5C_five_ML_models_gene_matrix.png"), final_plot, width = 183, height = 112,
    units = "mm", dpi = 600, bg = "white")
  message("Completed: ", out_dir)
}

acti_validation_clinical_expression <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(ggpubr)
    library(ggtext)
  })
  base_dir <- cfg$qp_dir
  output_dir <- cfg$clinical_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  genes7 <- c("FCGR1A", "DUSP3", "NPC2", "ANKRD22", "TAP1", "VAMP5", "GBP5")
  full <- read.csv(file.path(base_dir, "qpcr_log2_expr_noshifted_merged_allGenes.csv"), row.names = 1, check.names = FALSE)
  small <- read.csv(file.path(base_dir, "qpcr_log2_expr_noShift_20genes.csv"), row.names = 1, check.names = FALSE)
  stopifnot(all(genes7 %in% rownames(full)), all(colnames(small) %in% colnames(full)))
  shared <- intersect(rownames(small), rownames(full))
  stopifnot(isTRUE(all.equal(as.matrix(small[shared, , drop = FALSE]), as.matrix(full[shared, colnames(small),
    drop = FALSE]), check.attributes = FALSE)))
  set.seed(1234)
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(ggpubr)
  })
  expr_file <- file.path(base_dir, "qpcr_log2_expr_noShift_20genes.csv")
  pheno_file <- file.path(base_dir, "pheno.xlsx")
  out_dir <- file.path(output_dir, "Violin_20genes_singlePDF")
  if (!dir.exists(out_dir))
    dir.create(out_dir)
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  colnames(expr) <- trimws(colnames(expr))
  missing_target <- setdiff(genes7, rownames(expr))
  expr <- rbind(expr, as.matrix(full[missing_target, colnames(expr), drop = FALSE]))
  expr <- expr[genes7, , drop = FALSE]
  genes_20 <- genes7
  cat("\u8868\u8fbe\u77e9\u9635\uff1a", nrow(expr), "genes \u00d7", ncol(expr), "samples\n")
  pheno <- read_excel(pheno_file) %>% as.data.frame() %>% rename_with(trimws) %>% mutate(id = trimws(as.character(id)),
    Group = trimws(as.character(Group))) %>% filter(!is.na(id) & id != "")
  common_ids <- intersect(colnames(expr), pheno$id)
  cat("\u4ea4\u96c6\u6837\u672c\u6570\uff1a", length(common_ids), "\n")
  expr_only <- setdiff(colnames(expr), pheno$id)
  if (length(expr_only) > 0) {
    cat("\u26a0\ufe0f \u8868\u8fbe\u91cc\u8fd9\u4e9b\u6837\u672c\u4e0d\u5728 pheno\uff0c\u5c06\u5254\u9664\uff1a\n")
    print(expr_only)
  }
  if (length(common_ids) < 2) {
    stop("\u4ea4\u96c6\u6837\u672c\u592a\u5c11\uff0c\u8bf7\u68c0\u67e5\u8868\u8fbe\u5217\u540d\u4e0e pheno$id \u662f\u5426\u4e00\u81f4\u3002")
  }
  expr_use <- expr[, common_ids, drop = FALSE]
  pheno_use <- pheno %>% filter(id %in% common_ids) %>% arrange(match(id, colnames(expr_use))) %>% mutate(Group = factor(Group,
    levels = c("A", "B", "C", "D")))
  cat("\u4ea4\u96c6\u540e\u5404\u7ec4\u6837\u672c\u6570\uff1a\n")
  print(table(pheno_use$Group))
  df_long <- as.data.frame(expr_use) %>% rownames_to_column("Gene") %>% pivot_longer(-Gene, names_to = "id",
    values_to = "Expression") %>% mutate(id = trimws(as.character(id)), Expression = as.numeric(Expression)) %>%
    left_join(pheno_use %>% dplyr::select(id, Group), by = "id") %>% filter(!is.na(Group), is.finite(Expression)) %>%
    mutate(Group = factor(Group, levels = c("A", "B", "C", "D")))
  fill_cols <- c(A = "#4CAF50", B = "#FFC107", C = "#E53935", D = "#1E88E5")
  group_labels <- c(A = "Healthy", B = "LTBI", C = "ATB", D = "Other")
  plot_one_gene <- function(g) {
    d <- df_long %>% filter(Gene == g)
    if (nrow(d) == 0)
      return(NULL)
    present_groups <- levels(droplevels(d$Group))
    present_groups <- present_groups[present_groups %in% unique(as.character(d$Group))]
    pairwise_list <- NULL
    if (length(present_groups) >= 2) {
      pairwise_list <- combn(present_groups, 2, simplify = FALSE)
    }
    p <- ggplot(d, aes(x = Group, y = Expression, fill = Group)) + geom_violin(trim = FALSE, alpha = 0.85,
      color = NA) + geom_boxplot(linewidth = 0.9, width = 0.14, outlier.shape = NA,
      alpha = 0.35) + geom_jitter(width = 0.12, size = 1.1, alpha = 0.7) +
      scale_fill_manual(values = fill_cols, drop = FALSE) + scale_x_discrete(labels = group_labels, drop = FALSE) +
      labs(title = NULL, x = NULL, y = "Expression") + theme_bw(base_size = 12) + theme(legend.position = "none",
      plot.title = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1, size = 16, face = "bold",
        colour = "black"), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey30",
        fill = NA, linewidth = 1.2), axis.ticks = element_line(linewidth = 0.8))
    if (!is.null(pairwise_list)) {
      p <- p + ggpubr::stat_compare_means(comparisons = pairwise_list, method = "wilcox.test", label = "p.signif",
        hide.ns = FALSE, tip.length = 0.01, bracket.size = 0.8, size = 4.8)
    }
    return(p)
  }
  saved_files <- c()
  for (g in genes_20) {
    p <- plot_one_gene(g)
    if (!is.null(p)) {
      out_pdf <- file.path(out_dir, paste0(g, "_violin_ABCD.pdf"))
      ggsave(out_pdf, p, width = 6, height = 6)
      ggsave(sub("pdf$", "png", out_pdf), p, width = 6, height = 6, dpi = 600)
      saved_files <- c(saved_files, out_pdf)
      cat("\u5df2\u4fdd\u5b58\uff1a", out_pdf, "\n")
    }
    else {
      cat("\u8df3\u8fc7\uff1a", g, "\uff08\u65e0\u53ef\u7528\u6570\u636e\uff09\n")
    }
  }
  summary_df <- data.frame(Gene = genes_20, PDF = file.path(out_dir, paste0(genes_20, "_violin_ABCD.pdf")), Exists_in_plot = genes_20 %in%
    unique(df_long$Gene), stringsAsFactors = FALSE)
  write.csv(summary_df, file.path(out_dir, "Violin_20genes_singlePDF_summary.csv"), row.names = FALSE)
  cat("\u2705 \u5168\u90e8\u5b8c\u6210\n")
  cat("\u8f93\u51fa\u76ee\u5f55\uff1a", normalizePath(out_dir), "\n")
  write.csv(df_long, file.path(out_dir, "Figure7_violin_source_data.csv"), row.names = FALSE)
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggtext)
  })
  expr_file <- file.path(base_dir, "temp2.xls")
  clin_file <- file.path(base_dir, "qpcr_20genes_expr_matched_lingchaung2.csv")
  out_dir <- file.path(output_dir, "Scatter_20genes_each_pdf")
  if (!dir.exists(out_dir))
    dir.create(out_dir)
  genes_20 <- genes7
  expr <- readxl::read_excel(expr_file) %>% as.data.frame()
  colnames(expr)[1] <- "Gene"
  expr$Gene <- trimws(as.character(expr$Gene))
  scatter_ids <- names(expr)[-1]
  stopifnot(all(scatter_ids %in% colnames(full)))
  common_g <- intersect(expr$Gene, rownames(full))
  stopifnot(isTRUE(all.equal(as.matrix(expr[match(common_g, expr$Gene), scatter_ids, drop = FALSE]), as.matrix(full[common_g,
    scatter_ids, drop = FALSE]), check.attributes = FALSE)))
  missing_target <- setdiff(genes7, expr$Gene)
  extra <- data.frame(Gene = missing_target, full[missing_target, scatter_ids, drop = FALSE], check.names = FALSE)
  expr <- rbind(expr, extra)
  expr_long <- expr %>% tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>% dplyr::mutate(Gene = trimws(as.character(Gene)),
    Sample = trimws(as.character(Sample)), Expression = as.numeric(Expression)) %>% dplyr::filter(Gene %in%
    genes_20)
  clin <- read.csv(clin_file, stringsAsFactors = FALSE, check.names = FALSE) %>% dplyr::mutate(Sample = trimws(as.character(Sample)),
    days = suppressWarnings(as.numeric(days))) %>% dplyr::select(Sample, days)
  df <- expr_long %>% dplyr::left_join(clin, by = "Sample") %>% dplyr::filter(is.finite(Expression), is.finite(days)) %>%
    dplyr::mutate(TimeGroup = ifelse(days <= 30, "<=30 days", ">30 days"), TimeGroup = factor(TimeGroup, levels = c("<=30 days",
      ">30 days")))
  missing_days <- expr_long %>% dplyr::left_join(clin, by = "Sample") %>% dplyr::filter(is.na(days)) %>% dplyr::distinct(Sample)
  if (nrow(missing_days) > 0) {
    cat("\u26a0\ufe0f \u4ee5\u4e0b\u6837\u672c\u5728\u4e34\u5e8a\u8868\u4e2d\u6ca1\u6709 days\uff0c\u5df2\u4ece\u5206\u6790\u4e2d\u6392\u9664\uff1a\n")
    print(missing_days)
  }
  missing_genes <- setdiff(genes_20, unique(df$Gene))
  if (length(missing_genes) > 0) {
    cat("\u26a0\ufe0f \u4ee5\u4e0b\u57fa\u56e0\u5728\u8868\u8fbe\u6570\u636e\u4e2d\u672a\u627e\u5230\u6216\u65e0\u6709\u6548 days\uff0c\u5df2\u8df3\u8fc7\uff1a\n")
    print(missing_genes)
  }
  safe_cor <- function(x, y, method = "spearman") {
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    if (n_ok < 3) {
      return(c(r = NA_real_, p = NA_real_, n = n_ok))
    }
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
    c(r = unname(ct$estimate), p = unname(ct$p.value), n = n_ok)
  }
  anno_df <- df %>% dplyr::group_by(Gene) %>% dplyr::summarise(r_s = safe_cor(Expression, days, "spearman")["r"],
    p_s = safe_cor(Expression, days, "spearman")["p"], n_s = safe_cor(Expression, days, "spearman")["n"], r_p = safe_cor(Expression,
      days, "pearson")["r"], p_p = safe_cor(Expression, days, "pearson")["p"], n_p = safe_cor(Expression,
      days, "pearson")["n"], .groups = "drop") %>% dplyr::mutate(label = ifelse(is.na(r_s) | is.na(r_p),
    "Spearman r = NA (p=NA)<br>Pearson r = NA (p=NA)", sprintf("Spearman r = %.2f (p=%.3g)<br>Pearson r = %.2f (p=%.3g)",
      r_s, p_s, r_p, p_p)))
  time_cols <- c(`<=30 days` = "#1E88E5", `>30 days` = "#E53935")
  plot_one_gene <- function(g) {
    d <- df %>% dplyr::filter(Gene == g)
    a <- anno_df %>% dplyr::filter(Gene == g)
    if (nrow(d) == 0)
      return(NULL)
    x_min <- min(d$days, na.rm = TRUE)
    x_max <- max(d$days, na.rm = TRUE)
    x_breaks <- sort(unique(c(pretty(c(x_min, x_max), n = 5), 30)))
    p <- ggplot(d, aes(x = days, y = Expression, color = TimeGroup)) + geom_point(size = 2.8,
      alpha = 0.9) + geom_vline(xintercept = 30, linetype = "dashed", linewidth = 1.3, color = "grey40") +
      geom_smooth(aes(group = 1), method = "loess", se = FALSE, span = 1.1, linewidth = 1.5,
        color = "black") + geom_smooth(aes(group = 1), method = "lm", se = FALSE, linewidth = 1.3, color = "#D62728",
      linetype = "dashed") + ggtext::geom_richtext(data = a, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = 1.04, vjust = 1.04, size = 6, fill = "white", label.color = "black", label.size = 0.8,
      family = "", lineheight = 1.15) + scale_color_manual(values = time_cols, drop = FALSE) +
      scale_x_continuous(breaks = x_breaks) + theme_bw(base_size = 13) + theme(legend.position = "none",
      plot.title = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey30",
        fill = NA, linewidth = 1.2), axis.ticks = element_line(linewidth = 0.8)) + labs(title = NULL,
      x = "Treatment time (days)", y = "Expression (log2, -dCt)")
    return(p)
  }
  saved_files <- c()
  for (g in genes_20) {
    p <- plot_one_gene(g)
    if (!is.null(p)) {
      out_file <- file.path(out_dir, paste0(g, "_scatter_time.pdf"))
      ggsave(out_file, p, width = 6, height = 6)
      ggsave(sub("pdf$", "png", out_file), p, width = 6, height = 6, dpi = 600)
      saved_files <- c(saved_files, out_file)
      cat("\u5df2\u4fdd\u5b58\uff1a", out_file, "\n")
    }
    else {
      cat("\u8df3\u8fc7\uff1a", g, "\uff08\u65e0\u53ef\u7528\u6570\u636e\uff09\n")
    }
  }
  cor_summary <- anno_df %>% dplyr::select(Gene, r_s, p_s, n_s, r_p, p_p, n_p)
  write.csv(cor_summary, file.path(out_dir, "Correlation_summary_20genes.csv"), row.names = FALSE)
  cat("\u2705 \u5168\u90e8\u5b8c\u6210\n")
  cat("\u8f93\u51fa\u76ee\u5f55\uff1a", normalizePath(out_dir), "\n")
  stopifnot(all(genes7 %in% unique(df$Gene)))
  write.csv(df, file.path(out_dir, "Figure7_scatter_source_data.csv"), row.names = FALSE)
}

acti_validation_legends <- function(cfg = acti_validation_config()) {
  if (!cfg$save_plots)
    return(invisible(NULL))
  local({
    suppressPackageStartupMessages({
      library(ggplot2)
      library(grid)
    })
    out_dir <- file.path(cfg$output_dir, "06_legends", "seven_signatures")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    legend_order <- c("ACTI-1", "ACTI-2", "ACTI-3", "ACTI-4", "ACTI-5", "RISK6", "Sweeney3")
    legend_colors <- c(`ACTI-1` = "#F8766D", `ACTI-2` = "#E69F00", `ACTI-3` = "#0072B2", `ACTI-4` = "#56B4E9",
      `ACTI-5` = "#CC79A7", RISK6 = "#00BA38", Sweeney3 = "#619CFF")
    legend_df <- do.call(rbind, lapply(seq_along(legend_order), function(i) {
      data.frame(x = c(1, 2), y = i, Signature = factor(legend_order[i], levels = legend_order))
    }))
    p <- ggplot(legend_df, aes(x = x, y = y, colour = Signature)) + geom_line(size = 1.35) +
      scale_colour_manual(values = legend_colors, breaks = legend_order, limits = legend_order, drop = FALSE,
        name = NULL) + guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1.35))) +
      theme_void() + theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center",
      legend.text = element_text(size = 16, colour = "black"), legend.key.width = unit(13, "mm"), legend.key.height = unit(5,
        "mm"), legend.spacing.x = unit(4.5, "mm"), legend.margin = margin(1.5, 2, 1.5, 2, unit = "mm"),
      plot.margin = margin(0, 0, 0, 0))
    gt <- ggplotGrob(p)
    guide_idx <- which(vapply(gt$grobs, function(x) grepl("^guide-box", x$name), logical(1)))
    guide_idx <- guide_idx[vapply(gt$grobs[guide_idx], function(x) {
      !inherits(x, "zeroGrob")
    }, logical(1))]
    if (!length(guide_idx))
      stop("Legend grob was not found.")
    legend_grob <- gt$grobs[[guide_idx[1]]]
    width_in <- 13.2
    height_in <- 0.72
    pdf(file.path(out_dir, "Figure5D_shared_legend_all7_one_row.pdf"), width = width_in, height = height_in,
      bg = "white")
    grid.newpage()
    grid.draw(legend_grob)
    dev.off()
    png(file.path(out_dir, "Figure5D_shared_legend_all7_one_row_white.png"), width = width_in, height = height_in,
      units = "in", res = 600, bg = "white", type = "cairo")
    grid.newpage()
    grid.draw(legend_grob)
    dev.off()
    png(file.path(out_dir, "Figure5D_shared_legend_all7_one_row_transparent.png"), width = width_in, height = height_in,
      units = "in", res = 600, bg = "transparent", type = "cairo")
    grid.newpage()
    grid.draw(legend_grob)
    dev.off()
    message("Completed: ", out_dir)
  })
  local({
    suppressPackageStartupMessages({
      library(ggplot2)
      library(grid)
    })
    out_dir <- file.path(cfg$output_dir, "06_legends", "five_signatures")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    order <- paste0("ACTI-", 1:5)
    cols <- setNames(c("#F8766D", "#E69F00", "#0072B2", "#56B4E9", "#CC79A7"), order)
    d <- data.frame(x = rep(1:2, 5), y = rep(1:5, each = 2), Signature = factor(rep(order, each = 2), levels = order))
    p <- ggplot(d, aes(x, y, colour = Signature)) + geom_line(linewidth = 1.8) + scale_colour_manual(values = cols,
      breaks = order, name = NULL) + guides(colour = guide_legend(ncol = 1, byrow = TRUE)) + theme_void() +
      theme(legend.position = "right", legend.direction = "vertical", legend.text = element_text(size = 16,
        colour = "black"), legend.key.width = unit(13, "mm"), legend.key.height = unit(10, "mm"), legend.margin = margin(2,
        2, 2, 2, unit = "mm"))
    gt <- ggplotGrob(p)
    idx <- which(vapply(gt$grobs, function(g) grepl("^guide-box", g$name) && !inherits(g, "zeroGrob"), logical(1)))
    stopifnot(length(idx) > 0)
    g <- gt$grobs[[idx[1]]]
    prefix <- file.path(out_dir, "Figure7_shared_legend_ACTI_five_vertical")
    pdf(paste0(prefix, ".pdf"), width = 2.5, height = 2.6, bg = "white")
    grid.newpage()
    grid.draw(g)
    dev.off()
    for (bg in c("white", "transparent")) {
      png(paste0(prefix, "_", bg, ".png"), width = 2.5, height = 2.6, units = "in", res = 600,
        type = "cairo", bg = bg)
      grid.newpage()
      grid.draw(g)
      dev.off()
    }
  })
  local({
    suppressPackageStartupMessages({
      library(ggplot2)
      library(grid)
    })
    out_dir <- file.path(cfg$output_dir, "06_legends", "treatment_time")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    d <- data.frame(TimeGroup = factor(c("<=30 days", ">30 days"), levels = c("<=30 days", ">30 days")), x = 1:2,
      y = 1:2)
    p <- ggplot(d, aes(x, y, colour = TimeGroup)) + geom_point(size = 4, alpha = 0.9) + scale_colour_manual(values = c("#1E88E5",
      "#E53935"), name = "Treatment time", drop = FALSE) + guides(colour = guide_legend(ncol = 1, byrow = TRUE,
      override.aes = list(size = 4.5, alpha = 0.9))) + theme_void() + theme(legend.position = "right",
      legend.direction = "vertical", legend.title = element_text(size = 16, colour = "black"), legend.text = element_text(size = 16,
        colour = "black"), legend.key.height = unit(10, "mm"), legend.key.width = unit(9, "mm"), legend.margin = margin(3,
        3, 3, 3, unit = "mm"))
    gt <- ggplotGrob(p)
    idx <- which(vapply(gt$grobs, function(g) grepl("^guide-box", g$name) && !inherits(g, "zeroGrob"), logical(1)))
    stopifnot(length(idx) > 0)
    g <- gt$grobs[[idx[1]]]
    prefix <- file.path(out_dir, "Figure7_treatment_time_vertical_legend")
    pdf(paste0(prefix, ".pdf"), width = 3.1, height = 1.8, bg = "white")
    grid.newpage()
    grid.draw(g)
    dev.off()
    for (bg in c("white", "transparent")) {
      png(paste0(prefix, "_", bg, ".png"), width = 3.1, height = 1.8, units = "in", res = 600,
        type = "cairo", bg = bg)
      grid.newpage()
      grid.draw(g)
      dev.off()
    }
    cat("OUTPUT:", prefix, "\n")
  })
}

acti_validation_ml_heatmap <- function(cfg = acti_validation_config()) {
  ggsave <- if (cfg$save_plots)
    ggplot2::ggsave
  else function(...) invisible(NULL)
  if (is.null(cfg$ml_run_dir)) {
    message("Optional ML heatmap skipped: supply --ml-run-dir PATH to a completed run.")
    return(invisible(NULL))
  }
  run_dir <- cfg$ml_run_dir
  output_dir <- file.path(cfg$output_dir, "07_ML_heatmap")
  file_prefix <- "Fig5A_five_column_AUC_heatmap"
  prediction_packages <- c("randomForestSRC", "glmnet", "plsRglm", "gbm", "mboost", "e1071", "MASS", "xgboost",
    "pROC")
  available <- vapply(prediction_packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!available))
    stop("Missing prediction packages: ", paste(prediction_packages[!available], collapse = ", "))
  suppressPackageStartupMessages(invisible(lapply(prediction_packages, library, character.only = TRUE)))
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(grid)
    library(pROC)
  })
  apply_scaler <- function(x, s) {
    if (!all(s$features %in% colnames(x)))
      stop("<U+9884><U+6D4B><U+6570><U+636E><U+7F3A><U+5C11><U+57FA><U+56E0>")
    sweep(sweep(x[, s$features, drop = FALSE], 2, s$center, "-"), 2, s$scale, "/")
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
  CalPredictScore <- function(fit, new_data, type = "response") {
    x <- new_data[, fit$subFeature, drop = FALSE]
    p <- switch(class(fit)[1], lognet = predict(fit, newx = as.matrix(x), type = "response"), glm = predict(fit,
      newdata = as.data.frame(x), type = "response"), svm.formula = {
      pp <- attr(predict(fit, as.data.frame(x), probability = TRUE), "probabilities")
      if (is.null(pp) || !"1" %in% colnames(pp)) stop("SVM<U+672A><U+8FD4><U+56DE><U+7C7B><U+522B>1<U+6982><U+7387>")
      pp[, "1"]
    }, lda = predict(fit, as.data.frame(x))$posterior[, "1"], glmboost = predict(fit, newdata = as.data.frame(x),
      type = "response"), plsRglmmodel = predict(fit, newdata = as.data.frame(x), type = "response"), rfsrc = predict(fit,
      as.data.frame(x))$predicted[, "1"], gbm = predict(fit, newdata = as.data.frame(x), n.trees = fit$best_trees,
      type = "response"), xgb.Booster = predict(fit, as.matrix(x)), naiveBayes = predict(fit, newdata = as.data.frame(x),
      type = "raw")[, "1"], stop("<U+672A><U+652F><U+6301><U+6A21><U+578B><U+7C7B><U+522B><U+FF1A>", class(fit)[1]))
    p <- as.numeric(p)
    if (length(p) != nrow(x) || any(!is.finite(p)) || any(p < -1e-8 | p > 1 + 1e-8))
      stop("<U+6982><U+7387><U+65E0><U+6548>")
    setNames(pmin(pmax(p, 0), 1), rownames(x))
  }
  read_expression <- function(path) {
    z <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"Type" %in% names(z) || anyDuplicated(rownames(z)))
      stop("<U+7F3A><U+5C11>Type<U+6216><U+6837><U+672C><U+91CD><U+590D>")
    y <- suppressWarnings(as.numeric(as.character(z$Type)))
    if (anyNA(y) || !all(y %in% c(0, 1)) || length(unique(y)) != 2L)
      stop("Type<U+5FC5><U+987B><U+5305><U+542B>0<U+548C>1")
    x <- as.matrix(z[, setdiff(names(z), "Type"), drop = FALSE])
    storage.mode(x) <- "double"
    if (!ncol(x) || anyDuplicated(colnames(x)) || any(!is.finite(x)))
      stop("<U+8868><U+8FBE><U+77E9><U+9635><U+6709><U+975E><U+6570><U+503C>/NA/<U+91CD><U+590D><U+57FA><U+56E0>")
    list(x = x, y = y)
  }
  infer_sample_manifest <- function(x, y, dataset) {
    raw <- rownames(x)
    ids <- sub("^Train\\.", "", raw)
    pattern <- "(?:^|_)GSM(S[0-9]+)(DX|day[0-9]+|week[0-9]+)_(NonEndpoint|Endpoint)$"
    matches <- regmatches(ids, regexec(pattern, ids, perl = TRUE, ignore.case = TRUE))
    good <- lengths(matches) == 4L
    if (!all(good))
      stop("<U+65E0><U+6CD5><U+89E3><U+6790><U+6837><U+672C>ID<U+FF1A>", paste(head(ids[!good], 8), collapse = "; "),
        "<U+3002><U+9700><U+8981><U+7C7B><U+4F3C> GSE89403part1_GSMS147week4_NonEndpoint <U+6216> week24_Endpoint <U+7684><U+540D><U+79F0>")
    subject <- toupper(vapply(matches, `[`, character(1), 2L))
    timepoint <- tolower(vapply(matches, `[`, character(1), 3L))
    group <- vapply(matches, `[`, character(1), 4L)
    expected <- ifelse(tolower(group) == "endpoint", 1L, 0L)
    if (any(expected != y))
      stop("<U+6837><U+672C><U+540D>NonEndpoint/Endpoint<U+4E0E>Type<U+4E0D><U+4E00><U+81F4><U+FF1A>", paste(head(ids[expected !=
        y], 8), collapse = "; "))
    number <- suppressWarnings(as.numeric(sub("^(day|week)", "", timepoint)))
    week <- ifelse(timepoint == "dx", 0, ifelse(startsWith(timepoint, "day"), number/7, number))
    data.frame(sample_id = ids, subject_id = subject, timepoint = timepoint, sample_week = week, group = ifelse(expected ==
      1, "Endpoint", "NonEndpoint"), dataset = dataset, Type = y, output_sample_id = raw, stringsAsFactors = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  train <- read_expression(file.path(run_dir, "data.train.txt"))
  test <- read_expression(file.path(run_dir, "data.test.txt"))
  train_manifest <- infer_sample_manifest(train$x, train$y, "Training")
  test_manifest <- infer_sample_manifest(test$x, test$y, "Internal_test")
  extract_part <- function(ids) {
    z <- sub(".*GSE89403part([123]).*", "\\1", ids)
    ifelse(z %in% c("1", "2", "3"), paste0("Part ", LETTERS[as.integer(z)]), NA_character_)
  }
  train_manifest$part <- extract_part(train_manifest$output_sample_id)
  if (anyNA(train_manifest$part))
    stop("Some training samples cannot be assigned to Part A/B/C")
  ranking <- read.csv(file.path(run_dir, "model.training_CV_ranking.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  required_ranking_columns <- c("Method", "Train_patient_balanced_OOF_AUC", "complete_OOF", "eligible")
  if (!all(required_ranking_columns %in% names(ranking))) {
    stop("Ranking file is missing required columns")
  }
  method_order <- ranking$Method
  selected_method <- trimws(readLines(file.path(run_dir, "selected_method.txt"), warn = FALSE)[1])
  oof_df <- read.table(file.path(run_dir, "model.training_OOF_riskMatrix.txt"), header = TRUE, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE)
  rownames(oof_df) <- oof_df$id
  oof_df$id <- NULL
  if (!all(rownames(train$x) %in% rownames(oof_df)))
    stop("OOF predictions do not cover all training samples")
  if (!all(method_order %in% colnames(oof_df)))
    stop("OOF predictions do not cover all ranked methods")
  oof <- as.matrix(oof_df[rownames(train$x), method_order, drop = FALSE])
  storage.mode(oof) <- "double"
  part_auc <- sapply(c("Part A", "Part B", "Part C"), function(part_name) {
    ii <- train_manifest$part == part_name
    vapply(method_order, function(method) {
      patient_balanced_auc(train$y[ii], oof[ii, method], train_manifest$subject_id[ii])
    }, numeric(1))
  })
  colnames(part_auc) <- c("GSE89403part1", "GSE89403part2", "GSE89403part3")
  scaler <- readRDS(file.path(run_dir, "preprocessing.rds"))
  models <- readRDS(file.path(run_dir, "model.MLmodel.rds"))
  test_set <- apply_scaler(test$x, scaler)
  test_auc <- setNames(rep(NA_real_, length(method_order)), method_order)
  test_error <- setNames(rep("", length(method_order)), method_order)
  for (method in intersect(method_order, names(models))) {
    pred <- tryCatch(CalPredictScore(models[[method]], test_set), error = function(e) e)
    if (inherits(pred, "error")) {
      test_error[method] <- conditionMessage(pred)
    }
    else {
      test_auc[method] <- patient_balanced_auc(test$y, pred, test_manifest$subject_id)
    }
  }
  overall_train_auc <- setNames(ranking$Train_patient_balanced_OOF_AUC, ranking$Method)[method_order]
  auc_mat <- cbind(part_auc[method_order, , drop = FALSE], `Internal test` = test_auc[method_order], `Training OOF` = overall_train_auc)
  rownames(auc_mat) <- method_order
  eligible <- setNames(as.logical(ranking$eligible), ranking$Method)[method_order]
  complete_oof <- setNames(as.logical(ranking$complete_OOF), ranking$Method)[method_order]
  display_mat <- auc_mat
  display_mat[!(eligible & complete_oof), ] <- NA_real_
  source_data <- data.frame(Method = method_order, Part_A_training_patient_balanced_OOF_AUC = auc_mat[, 1], Part_B_training_patient_balanced_OOF_AUC = auc_mat[,
    2], Part_C_training_patient_balanced_OOF_AUC = auc_mat[, 3], Internal_test_patient_balanced_AUC = auc_mat[,
    4], Overall_training_patient_balanced_OOF_AUC = auc_mat[, 5], complete_OOF = complete_oof, eligible = eligible,
    selected = method_order == selected_method, test_prediction_error = unname(test_error[method_order]), stringsAsFactors = FALSE,
    check.names = FALSE)
  write.csv(source_data, file.path(output_dir, paste0(file_prefix, "_source_data.csv")), row.names = FALSE)
  SimpleHeatmap <- function(Cindex_mat, rank_metric, CohortCol, barCol = "steelblue", cellwidth = 1, cellheight = 0.5,
    selected_method = NULL) {
    column_labels <- colnames(Cindex_mat)
    column_groups <- factor(column_labels, levels = column_labels)
    col_ha <- columnAnnotation(Cohort = column_labels, col = list(Cohort = CohortCol), show_annotation_name = FALSE,
      simple_anno_size = unit(4.5, "mm"), annotation_legend_param = list(Cohort = list(ncol = 1, labels_gp = gpar(fontsize = 9))))
    rank_plot <- rank_metric
    rank_plot[!is.finite(rank_plot)] <- 0
    row_ha <- rowAnnotation(bar = anno_barplot(rank_plot, bar_width = 0.8, border = FALSE,
      ylim = c(0, 1), extend = 0, gp = gpar(fill = barCol, col = NA), add_numbers = FALSE, axis_param = list(at = c(0,
        0.5, 1), labels = c("0", "0.5", "1"), labels_rot = 0), width = unit(3, "cm")), show_annotation_name = FALSE)
    row_face <- ifelse(rownames(Cindex_mat) == selected_method, "bold", "plain")
    Heatmap(as.matrix(Cindex_mat), name = "AUC", right_annotation = row_ha, top_annotation = col_ha, col = c("#4195C1",
      "#FFFFFF", "#FFBC90"), na_col = "#BDBDBD", rect_gp = gpar(col = "black", lwd = 1), cluster_columns = FALSE,
      cluster_rows = FALSE, show_column_names = FALSE, show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 9,
        fontface = row_face), width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"), height = unit(cellheight *
        nrow(Cindex_mat), "cm"), column_split = NULL, cell_fun = function(j, i, x, y, w, h, col) {
        label <- if (is.finite(Cindex_mat[i, j])) {
         format(Cindex_mat[i, j], digits = 3, nsmall = 3)
        }
        else {
         "NA"
        }
        grid.text(label, x, y, gp = gpar(fontsize = 10))
      }, heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9)))
  }
  cohort_col <- brewer.pal(n = 5, name = "Paired")
  names(cohort_col) <- colnames(display_mat)
  hm <- SimpleHeatmap(Cindex_mat = display_mat, rank_metric = overall_train_auc, CohortCol = cohort_col, barCol = "steelblue",
    cellwidth = 1, cellheight = 0.5, selected_method = selected_method)
  draw_heatmap <- function(ht, bar_values) {
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    decorate_annotation("bar", {
      n <- length(bar_values)
      values <- as.numeric(bar_values)
      values[!is.finite(values)] <- 0
      grid.text(sprintf("%.3f", values), x = unit(pmax(values - 0.015, 0.05),
        "native"), y = unit(rev(seq_len(n)), "native"), just = "right", gp = gpar(fontsize = 9, col = "white"))
    })
  }
  pdf_file <- file.path(output_dir, paste0(file_prefix, ".pdf"))
  png_file <- file.path(output_dir, paste0(file_prefix, ".png"))
  if (cfg$save_plots) {
    pdf(pdf_file, width = 1 * ncol(display_mat) + 6, height = 0.5 * nrow(display_mat) * 0.45,
      family = "sans", useDingbats = FALSE)
    draw_heatmap(hm, overall_train_auc)
    dev.off()
    png(png_file, width = 3300, height = 7200, res = 300, bg = "white")
    draw_heatmap(hm, overall_train_auc)
    dev.off()
  }
  message("Completed: ", output_dir)
}

acti_validation_run <- function(stage, cfg = acti_validation_config()) {
  stages <- c("internal", "external", "prediagnosis", "clinical_expression", "clinical_roc", "clinical_time",
    "model_matrix", "legends", "ml_heatmap")
  if (length(stage) != 1L || !(stage %in% c("all", stages)))
    stop("Unknown stage: ", stage)
  chosen <- if (stage == "all")
    stages
  else stage
  result <- NULL
  for (s in chosen) {
    message("Running stage: ", s)
    result <- get(paste0("acti_validation_", s), mode = "function")(cfg)
  }
  invisible(result)
}

acti_validation_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- "Rscript 5_Validation_and_plots.R --stage NAME\nStages: all, internal, external, prediagnosis, clinical_expression, clinical_roc, clinical_time, model_matrix, legends, ml_heatmap\nUse --check to check inputs; see the configuration function for paths.\n"
  if (!length(args) || "--help" %in% args) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  flags <- c("--data-dir", "--bulk-dir", "--qp-dir", "--pre-dir", "--pre-source-file", "--clinical-output-dir",
    "--output-dir", "--ml-run-dir", "--stage")
  values <- list()
  stage <- "all"
  check <- FALSE
  i <- 1L
  while (i <= length(args)) {
    f <- args[[i]]
    if (f == "--check")
      check <- TRUE
    else if (f == "--tables-only")
      values$save_plots <- FALSE
    else if (f == "--overwrite")
      values$overwrite <- TRUE
    else if (f %in% flags) {
      if (i == length(args) || startsWith(args[[i + 1L]], "--"))
        stop("Missing value: ", f)
      i <- i + 1L
      if (f == "--stage")
        stage <- args[[i]]
      else values[[gsub("-", "_", substring(f, 3L), fixed = TRUE)]] <- args[[i]]
    }
    else stop("Unknown argument: ", f)
    i <- i + 1L
  }
  cfg <- do.call(acti_validation_config, values)
  if (check) {
    acti_validation_check(cfg, stage)
    return(invisible(cfg))
  }
  acti_validation_run(stage, cfg)
}

if (sys.nframe() == 0L) acti_validation_main()

