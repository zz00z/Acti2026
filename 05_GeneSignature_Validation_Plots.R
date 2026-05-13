############################################################
## 05_GeneSignature_Validation_Plots.R
##
## Representative validation code for gene signatures.
##
## This script includes three typical analyses:
## 1. Violin plots for gene expression and signature scores
## 2. Correlation analysis between gene expression/signature
##    scores and treatment time
## 3. ROC/AUC analysis for classification performance
##
## Six-gene signature:
##   ANKRD22, LHFPL2, CASP4, GYG1, GBP1, NPC2
##
## Input files should be prepared by users according to the
## templates described below.
############################################################


## =========================================================
## 0. Load required packages
## =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(pROC)
})


## =========================================================
## 1. Set working directory and output folders
## =========================================================
## Please modify this path according to your local environment.

base_dir <- "D:/atraining/blood/analysis/validation"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
setwd(base_dir)

violin_dir <- file.path(base_dir, "validation_violin_plots")
cor_dir    <- file.path(base_dir, "validation_time_correlation")
roc_dir    <- file.path(base_dir, "validation_ROC")

dir.create(violin_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cor_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(roc_dir,    showWarnings = FALSE, recursive = TRUE)


## =========================================================
## 2. Define gene signatures
## =========================================================

six_gene_signature <- c(
  "ANKRD22",
  "LHFPL2",
  "CASP4",
  "GYG1",
  "GBP1",
  "NPC2"
)

ten_candidate_genes <- c(
  "ANKRD22",
  "LHFPL2",
  "CASP4",
  "GYG1",
  "GBP1",
  "NPC2",
  "VAMP5",
  "PSTPIP2",
  "GBP5",
  "FCGR1A"
)

risk6_up_genes <- c(
  "GBP2",
  "FCGR1B",
  "SERPING1"
)

risk6_down_genes <- c(
  "TUBGCP6",
  "TRMT2A",
  "SDR39U1"
)


## =========================================================
## 3. Helper functions
## =========================================================

safe_cor <- function(x, y, method = "spearman") {
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3) {
    return(c(r = NA_real_, p = NA_real_))
  }
  
  ct <- suppressWarnings(
    cor.test(x[ok], y[ok], method = method)
  )
  
  c(
    r = unname(ct$estimate),
    p = unname(ct$p.value)
  )
}


calculate_mean_signature_score <- function(expr_mat, genes) {
  
  genes_found <- intersect(genes, rownames(expr_mat))
  genes_missing <- setdiff(genes, rownames(expr_mat))
  
  cat("Genes found:\n")
  print(genes_found)
  
  if (length(genes_missing) > 0) {
    cat("Genes missing:\n")
    print(genes_missing)
  }
  
  if (length(genes_found) == 0) {
    stop("None of the signature genes were found in the expression matrix.")
  }
  
  score <- colMeans(
    expr_mat[genes_found, , drop = FALSE],
    na.rm = TRUE
  )
  
  return(score)
}


calculate_risk6_score <- function(expr_mat, up_genes, down_genes) {
  
  up_found <- intersect(up_genes, rownames(expr_mat))
  down_found <- intersect(down_genes, rownames(expr_mat))
  
  if (length(up_found) == 0 || length(down_found) == 0) {
    warning("Some RISK6 up/down genes were not found. Returning NA score.")
    return(rep(NA_real_, ncol(expr_mat)))
  }
  
  up_score <- colMeans(
    expr_mat[up_found, , drop = FALSE],
    na.rm = TRUE
  )
  
  down_score <- colMeans(
    expr_mat[down_found, , drop = FALSE],
    na.rm = TRUE
  )
  
  risk6_score <- up_score - down_score
  
  return(risk6_score)
}


read_expression_matrix <- function(expr_file) {
  
  expr <- read.csv(
    expr_file,
    row.names = 1,
    check.names = FALSE
  )
  
  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  
  rownames(expr) <- trimws(rownames(expr))
  colnames(expr) <- trimws(colnames(expr))
  
  return(expr)
}


read_phenotype_table <- function(pheno_file) {
  
  if (grepl("\\.xlsx?$", pheno_file, ignore.case = TRUE)) {
    pheno <- readxl::read_excel(pheno_file) %>% as.data.frame()
  } else {
    pheno <- read.csv(
      pheno_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  pheno <- pheno %>%
    rename_with(trimws)
  
  return(pheno)
}


align_expression_and_pheno <- function(expr, pheno, sample_col = "id") {
  
  pheno[[sample_col]] <- trimws(as.character(pheno[[sample_col]]))
  
  common_ids <- intersect(colnames(expr), pheno[[sample_col]])
  
  if (length(common_ids) < 2) {
    stop("Too few matched samples between expression matrix and phenotype table.")
  }
  
  expr_use <- expr[, common_ids, drop = FALSE]
  
  pheno_use <- pheno %>%
    filter(.data[[sample_col]] %in% common_ids) %>%
    arrange(match(.data[[sample_col]], colnames(expr_use)))
  
  stopifnot(all(pheno_use[[sample_col]] == colnames(expr_use)))
  
  list(
    expr = expr_use,
    pheno = pheno_use
  )
}


## =========================================================
## Part I. Violin plots
## =========================================================
## Purpose:
##   Plot expression differences of candidate genes and
##   signature scores across validation groups.
##
## Required input files:
##   expression_matrix.csv
##     rows    = genes
##     columns = samples
##
##   pheno.xlsx or pheno.csv
##     required columns:
##       id
##       Group
##
## Example Group coding:
##   A = Healthy
##   B = LTBI
##   C = ATB
##   D = Other
##
## You may also use custom group names, such as:
##   HC, LTBI, ATB, OLD
##   TB, Week1, Week2, Month3, Month6
############################################################

run_violin_validation <- function(
    expr_file = "expression_matrix.csv",
    pheno_file = "pheno.xlsx",
    sample_col = "id",
    group_col = "Group",
    group_levels = c("A", "B", "C", "D"),
    group_labels = c(
      A = "Healthy",
      B = "LTBI",
      C = "ATB",
      D = "Other"
    ),
    fill_cols = c(
      A = "#4CAF50",
      B = "#FFC107",
      C = "#E53935",
      D = "#1E88E5"
    ),
    candidate_genes = ten_candidate_genes,
    signature_genes = six_gene_signature,
    output_prefix = "Validation"
) {
  
  cat("\n================ Violin plot validation ================\n")
  
  expr <- read_expression_matrix(expr_file)
  pheno <- read_phenotype_table(pheno_file)
  
  if (!all(c(sample_col, group_col) %in% colnames(pheno))) {
    stop("Phenotype table must contain sample and group columns.")
  }
  
  pheno[[group_col]] <- trimws(as.character(pheno[[group_col]]))
  
  aligned <- align_expression_and_pheno(
    expr = expr,
    pheno = pheno,
    sample_col = sample_col
  )
  
  expr_use <- aligned$expr
  pheno_use <- aligned$pheno
  
  pheno_use[[group_col]] <- factor(
    pheno_use[[group_col]],
    levels = group_levels
  )
  
  cat("Matched group distribution:\n")
  print(table(pheno_use[[group_col]], useNA = "ifany"))
  
  ## -------------------------
  ## A. Candidate gene violin plot
  ## -------------------------
  
  genes_found <- intersect(candidate_genes, rownames(expr_use))
  
  if (length(genes_found) > 0) {
    
    df_gene <- as.data.frame(expr_use[genes_found, , drop = FALSE]) %>%
      rownames_to_column("Gene") %>%
      pivot_longer(
        cols = -Gene,
        names_to = sample_col,
        values_to = "Expression"
      ) %>%
      mutate(
        !!sample_col := trimws(as.character(.data[[sample_col]])),
        Expression = as.numeric(Expression)
      ) %>%
      left_join(
        pheno_use %>% select(all_of(c(sample_col, group_col))),
        by = sample_col
      ) %>%
      filter(!is.na(.data[[group_col]]), is.finite(Expression)) %>%
      mutate(
        GroupForPlot = factor(.data[[group_col]], levels = group_levels)
      )
    
    p_gene <- ggplot(
      df_gene,
      aes(x = GroupForPlot, y = Expression, fill = GroupForPlot)
    ) +
      geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
      geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.35) +
      geom_jitter(width = 0.12, size = 0.8, alpha = 0.65) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
      scale_fill_manual(values = fill_cols, drop = FALSE) +
      scale_x_discrete(labels = group_labels, drop = FALSE) +
      labs(
        x = NULL,
        y = "Expression"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      file.path(violin_dir, paste0(output_prefix, "_Candidate_gene_violin.pdf")),
      p_gene,
      width = 12,
      height = 8
    )
    
    ggsave(
      file.path(violin_dir, paste0(output_prefix, "_Candidate_gene_violin.png")),
      p_gene,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
  
  ## -------------------------
  ## B. Signature score violin plot
  ## -------------------------
  
  signature_score <- calculate_mean_signature_score(
    expr_mat = expr_use,
    genes = signature_genes
  )
  
  score_df <- data.frame(
    id_tmp = names(signature_score),
    Score = as.numeric(signature_score),
    stringsAsFactors = FALSE
  )
  
  colnames(score_df)[1] <- sample_col
  
  score_df <- score_df %>%
    left_join(
      pheno_use %>% select(all_of(c(sample_col, group_col))),
      by = sample_col
    ) %>%
    filter(!is.na(.data[[group_col]]), is.finite(Score)) %>%
    mutate(
      GroupForPlot = factor(.data[[group_col]], levels = group_levels)
    )
  
  write.csv(
    score_df,
    file.path(violin_dir, paste0(output_prefix, "_Signature_score.csv")),
    row.names = FALSE
  )
  
  present_groups <- levels(droplevels(score_df$GroupForPlot))
  present_groups <- present_groups[present_groups %in% unique(as.character(score_df$GroupForPlot))]
  
  pairwise_list <- NULL
  if (length(present_groups) >= 2) {
    pairwise_list <- combn(present_groups, 2, simplify = FALSE)
  }
  
  p_score <- ggplot(
    score_df,
    aes(x = GroupForPlot, y = Score, fill = GroupForPlot)
  ) +
    geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.35) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_x_discrete(labels = group_labels, drop = FALSE) +
    labs(
      title = "Six-gene signature",
      x = NULL,
      y = "Mean expression"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(pairwise_list)) {
    p_score <- p_score +
      ggpubr::stat_compare_means(
        comparisons = pairwise_list,
        method = "wilcox.test",
        method.args = list(exact = FALSE),
        label = "p.signif",
        hide.ns = FALSE,
        tip.length = 0.01,
        size = 5
      )
  }
  
  ggsave(
    file.path(violin_dir, paste0(output_prefix, "_Six_gene_signature_violin.pdf")),
    p_score,
    width = 6,
    height = 6
  )
  
  ggsave(
    file.path(violin_dir, paste0(output_prefix, "_Six_gene_signature_violin.png")),
    p_score,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  cat("Violin plot validation completed.\n")
  
  invisible(
    list(
      expression = expr_use,
      phenotype = pheno_use,
      score = score_df
    )
  )
}


## =========================================================
## Part II. Correlation with treatment time
## =========================================================
## Purpose:
##   Plot correlations between gene expression/signature score
##   and treatment time.
##
## Required input files:
##   longitudinal_expression_matrix.xlsx or .csv
##     first column = Gene
##     other columns = samples
##
##   clinical_time_information.csv
##     required columns:
##       Sample
##       days
############################################################

read_longitudinal_expression <- function(expr_file) {
  
  if (grepl("\\.xlsx?$", expr_file, ignore.case = TRUE)) {
    expr <- readxl::read_excel(expr_file) %>% as.data.frame()
  } else {
    expr <- read.csv(
      expr_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  colnames(expr)[1] <- "Gene"
  
  expr <- expr %>%
    mutate(Gene = trimws(as.character(Gene)))
  
  return(expr)
}


run_time_correlation_validation <- function(
    expr_file = "longitudinal_expression_matrix.xlsx",
    clin_file = "clinical_time_information.csv",
    sample_col = "Sample",
    time_col = "days",
    highlight_genes = c("ANKRD22", "FCGR1A", "GBP5"),
    signature_genes = six_gene_signature,
    method = "spearman",
    output_prefix = "Validation"
) {
  
  cat("\n================ Treatment-time correlation validation ================\n")
  
  expr <- read_longitudinal_expression(expr_file)
  
  expr_long <- expr %>%
    pivot_longer(
      cols = -Gene,
      names_to = sample_col,
      values_to = "Expression"
    ) %>%
    mutate(
      Gene = trimws(as.character(Gene)),
      !!sample_col := trimws(as.character(.data[[sample_col]])),
      Expression = as.numeric(Expression)
    )
  
  clin <- read.csv(
    clin_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    mutate(
      !!sample_col := trimws(as.character(.data[[sample_col]])),
      !!time_col := suppressWarnings(as.numeric(.data[[time_col]]))
    ) %>%
    select(all_of(c(sample_col, time_col)))
  
  df_time <- expr_long %>%
    inner_join(clin, by = sample_col) %>%
    filter(is.finite(Expression), is.finite(.data[[time_col]])) %>%
    mutate(
      TimeGroup = ifelse(.data[[time_col]] <= 30, "<=30 days", ">30 days"),
      TimeGroup = factor(TimeGroup, levels = c("<=30 days", ">30 days"))
    )
  
  ## -------------------------
  ## A. Gene-level correlation
  ## -------------------------
  
  cor_gene <- df_time %>%
    group_by(Gene) %>%
    summarise(
      r = safe_cor(Expression, .data[[time_col]], method = method)["r"],
      p = safe_cor(Expression, .data[[time_col]], method = method)["p"],
      n = sum(is.finite(Expression) & is.finite(.data[[time_col]])),
      .groups = "drop"
    ) %>%
    mutate(
      FDR = p.adjust(p, method = "fdr"),
      label = sprintf(
        "%s r = %.2f\np = %.3g",
        tools::toTitleCase(method),
        r,
        p
      )
    ) %>%
    arrange(p)
  
  write.csv(
    cor_gene,
    file.path(cor_dir, paste0(output_prefix, "_Gene_expression_vs_treatment_time_correlation.csv")),
    row.names = FALSE
  )
  
  plot_df <- df_time %>%
    filter(Gene %in% highlight_genes) %>%
    left_join(cor_gene %>% select(Gene, label), by = "Gene")
  
  anno_df <- cor_gene %>%
    filter(Gene %in% highlight_genes)
  
  if (nrow(plot_df) > 0) {
    
    p_gene_time <- ggplot(
      plot_df,
      aes(x = .data[[time_col]], y = Expression, color = TimeGroup)
    ) +
      geom_point(size = 2.8, alpha = 0.9) +
      geom_smooth(
        aes(group = 1),
        method = "loess",
        se = TRUE,
        span = 0.9,
        linewidth = 0.8,
        color = "black"
      ) +
      geom_smooth(
        aes(group = 1),
        method = "lm",
        se = FALSE,
        linewidth = 0.7,
        color = "#D62728",
        linetype = "dashed"
      ) +
      facet_wrap(~ Gene, scales = "free_y", ncol = length(highlight_genes)) +
      geom_label(
        data = anno_df,
        aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.08,
        size = 4.2,
        label.size = 0.3,
        fill = "white",
        alpha = 0.9
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = "Gene expression vs treatment time",
        x = "Treatment time",
        y = "Expression",
        color = NULL
      )
    
    ggsave(
      file.path(cor_dir, paste0(output_prefix, "_Gene_expression_vs_treatment_time.pdf")),
      p_gene_time,
      width = 10,
      height = 5
    )
    
    ggsave(
      file.path(cor_dir, paste0(output_prefix, "_Gene_expression_vs_treatment_time.png")),
      p_gene_time,
      width = 10,
      height = 5,
      dpi = 300
    )
  }
  
  ## -------------------------
  ## B. Signature score vs treatment time
  ## -------------------------
  
  expr_wide <- expr %>%
    filter(Gene %in% signature_genes) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  mode(expr_wide) <- "numeric"
  
  genes_found <- intersect(signature_genes, rownames(expr_wide))
  
  if (length(genes_found) == 0) {
    stop("No signature genes were found in longitudinal expression matrix.")
  }
  
  score_time <- colMeans(
    expr_wide[genes_found, , drop = FALSE],
    na.rm = TRUE
  )
  
  score_time_df <- data.frame(
    Sample = names(score_time),
    Score = as.numeric(score_time),
    stringsAsFactors = FALSE
  )
  
  colnames(score_time_df)[1] <- sample_col
  
  score_time_df <- score_time_df %>%
    inner_join(clin, by = sample_col) %>%
    filter(is.finite(Score), is.finite(.data[[time_col]])) %>%
    mutate(
      TimeGroup = ifelse(.data[[time_col]] <= 30, "<=30 days", ">30 days"),
      TimeGroup = factor(TimeGroup, levels = c("<=30 days", ">30 days"))
    )
  
  ct_score <- suppressWarnings(
    cor.test(
      score_time_df$Score,
      score_time_df[[time_col]],
      method = method
    )
  )
  
  score_label <- sprintf(
    "%s r = %.2f\np = %.3g",
    tools::toTitleCase(method),
    unname(ct_score$estimate),
    ct_score$p.value
  )
  
  write.csv(
    data.frame(
      method = method,
      r = unname(ct_score$estimate),
      p = ct_score$p.value,
      n = nrow(score_time_df)
    ),
    file.path(cor_dir, paste0(output_prefix, "_Signature_score_vs_treatment_time_correlation.csv")),
    row.names = FALSE
  )
  
  p_score_time <- ggplot(
    score_time_df,
    aes(x = .data[[time_col]], y = Score, color = TimeGroup)
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_smooth(
      aes(group = 1),
      method = "loess",
      se = TRUE,
      span = 0.9,
      linewidth = 0.8,
      color = "black"
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      se = FALSE,
      linewidth = 0.7,
      color = "#D62728",
      linetype = "dashed"
    ) +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = score_label,
      hjust = 1.05,
      vjust = 1.08,
      size = 4.5,
      fill = "white",
      alpha = 0.9
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Gene signature vs treatment time",
      x = "Treatment time",
      y = "Signature score",
      color = NULL
    )
  
  ggsave(
    file.path(cor_dir, paste0(output_prefix, "_Signature_score_vs_treatment_time.pdf")),
    p_score_time,
    width = 6,
    height = 5
  )
  
  ggsave(
    file.path(cor_dir, paste0(output_prefix, "_Signature_score_vs_treatment_time.png")),
    p_score_time,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  cat("Treatment-time correlation validation completed.\n")
  
  invisible(
    list(
      gene_correlation = cor_gene,
      score_correlation = ct_score,
      score_data = score_time_df
    )
  )
}


## =========================================================
## Part III. ROC/AUC analysis
## =========================================================
## Purpose:
##   Evaluate classification performance of gene signatures.
##
## Required input file:
##   signature_score_table.csv
##
## Required columns:
##   Sample
##   Group
##   Score
##
## Example clinical comparisons:
##   LTBI vs ATB
##   Healthy vs ATB
##
## Example external cohort comparisons:
##   TB vs Week1
##   TB vs Month3
##   TB vs Month6
##   Week2 vs Month6
############################################################

run_roc_comparison <- function(
    score_df,
    group_control,
    group_case,
    positive_group = group_case,
    score_col = "Score"
) {
  
  dat <- score_df %>%
    filter(Group %in% c(group_control, group_case)) %>%
    mutate(
      Status = ifelse(Group == positive_group, 1, 0)
    )
  
  if (length(unique(dat$Status)) != 2) {
    stop(paste0("ROC comparison failed: ", group_control, " vs ", group_case))
  }
  
  roc_obj <- pROC::roc(
    response = dat$Status,
    predictor = dat[[score_col]],
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_value <- as.numeric(pROC::auc(roc_obj))
  ci_value <- as.numeric(pROC::ci.auc(roc_obj))
  
  list(
    data = dat,
    roc = roc_obj,
    auc = auc_value,
    ci = ci_value,
    comparison = paste(group_control, "vs", group_case)
  )
}


run_roc_validation <- function(
    score_file = "signature_score_table.csv",
    comparisons = list(
      LTBI_vs_ATB = c("LTBI", "ATB"),
      Healthy_vs_ATB = c("Healthy", "ATB")
    ),
    score_col = "Score",
    sample_col = "Sample",
    group_col = "Group",
    output_prefix = "Validation"
) {
  
  cat("\n================ ROC/AUC validation ================\n")
  
  score_df <- read.csv(
    score_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    mutate(
      !!sample_col := as.character(.data[[sample_col]]),
      !!group_col := as.character(.data[[group_col]]),
      !!score_col := as.numeric(.data[[score_col]])
    ) %>%
    filter(!is.na(.data[[group_col]]), is.finite(.data[[score_col]]))
  
  colnames(score_df)[colnames(score_df) == sample_col] <- "Sample"
  colnames(score_df)[colnames(score_df) == group_col] <- "Group"
  colnames(score_df)[colnames(score_df) == score_col] <- "Score"
  
  roc_list <- list()
  
  for (nm in names(comparisons)) {
    
    g1 <- comparisons[[nm]][1]
    g2 <- comparisons[[nm]][2]
    
    roc_list[[nm]] <- run_roc_comparison(
      score_df = score_df,
      group_control = g1,
      group_case = g2,
      positive_group = g2,
      score_col = "Score"
    )
  }
  
  roc_summary <- bind_rows(lapply(names(roc_list), function(nm) {
    
    x <- roc_list[[nm]]
    
    data.frame(
      Comparison = nm,
      AUC = x$auc,
      CI_lower = x$ci[1],
      CI_median = x$ci[2],
      CI_upper = x$ci[3],
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(
    roc_summary,
    file.path(roc_dir, paste0(output_prefix, "_ROC_AUC_summary.csv")),
    row.names = FALSE
  )
  
  roc_df <- bind_rows(lapply(names(roc_list), function(nm) {
    
    x <- roc_list[[nm]]
    
    data.frame(
      Specificity = x$roc$specificities,
      Sensitivity = x$roc$sensitivities,
      Comparison = sprintf(
        "%s, AUC = %.3f",
        gsub("_", " ", nm),
        x$auc
      )
    )
  }))
  
  p_roc <- ggplot(
    roc_df,
    aes(x = 1 - Specificity, y = Sensitivity, color = Comparison)
  ) +
    geom_line(linewidth = 1.2) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    coord_equal() +
    theme_bw(base_size = 14) +
    theme(
      legend.position = c(0.65, 0.18),
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Gene signature ROC analysis",
      x = "1 - Specificity",
      y = "Sensitivity"
    )
  
  ggsave(
    file.path(roc_dir, paste0(output_prefix, "_ROC_curves.pdf")),
    p_roc,
    width = 6,
    height = 6
  )
  
  ggsave(
    file.path(roc_dir, paste0(output_prefix, "_ROC_curves.png")),
    p_roc,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  cat("ROC/AUC validation completed.\n")
  
  invisible(
    list(
      roc_summary = roc_summary,
      roc_list = roc_list,
      roc_plot = p_roc
    )
  )
}


## =========================================================
## Example usage
## =========================================================
## The following commands are examples. Uncomment and modify
## the file names according to your own validation datasets.


## -------------------------
## Example 1. Violin plots
## -------------------------
## Required:
##   expression_matrix.csv
##   pheno.xlsx
##
## run_violin_validation(
##   expr_file = "expression_matrix.csv",
##   pheno_file = "pheno.xlsx",
##   sample_col = "id",
##   group_col = "Group",
##   group_levels = c("A", "B", "C", "D"),
##   group_labels = c(
##     A = "Healthy",
##     B = "LTBI",
##     C = "ATB",
##     D = "Other"
##   ),
##   fill_cols = c(
##     A = "#4CAF50",
##     B = "#FFC107",
##     C = "#E53935",
##     D = "#1E88E5"
##   ),
##   candidate_genes = ten_candidate_genes,
##   signature_genes = six_gene_signature,
##   output_prefix = "Clinical_validation"
## )


## -------------------------
## Example 2. Treatment-time correlation
## -------------------------
## Required:
##   longitudinal_expression_matrix.xlsx
##   clinical_time_information.csv
##
## run_time_correlation_validation(
##   expr_file = "longitudinal_expression_matrix.xlsx",
##   clin_file = "clinical_time_information.csv",
##   sample_col = "Sample",
##   time_col = "days",
##   highlight_genes = c("ANKRD22", "FCGR1A", "GBP5"),
##   signature_genes = six_gene_signature,
##   method = "spearman",
##   output_prefix = "Clinical_validation"
## )


## -------------------------
## Example 3. ROC/AUC analysis
## -------------------------
## Required:
##   signature_score_table.csv
##
## The score table should include:
##   Sample, Group, Score
##
## Clinical validation example:
##
## run_roc_validation(
##   score_file = "signature_score_table.csv",
##   comparisons = list(
##     LTBI_vs_ATB = c("LTBI", "ATB"),
##     Healthy_vs_ATB = c("Healthy", "ATB")
##   ),
##   score_col = "Score",
##   sample_col = "Sample",
##   group_col = "Group",
##   output_prefix = "Clinical_validation"
## )
##
## External cohort example:
##
## run_roc_validation(
##   score_file = "GSE40553_signature_score_table.csv",
##   comparisons = list(
##     TB_vs_Month6 = c("TB", "Month6"),
##     Week2_vs_Month6 = c("Week2", "Month6")
##   ),
##   score_col = "Score",
##   sample_col = "Sample",
##   group_col = "Group",
##   output_prefix = "GSE40553"
## )
##
## run_roc_validation(
##   score_file = "GSE54992_signature_score_table.csv",
##   comparisons = list(
##     TB_vs_Month3 = c("TB", "Month3"),
##     TB_vs_Month6 = c("TB", "Month6")
##   ),
##   score_col = "Score",
##   sample_col = "Sample",
##   group_col = "Group",
##   output_prefix = "GSE54992"
## )
##
## run_roc_validation(
##   score_file = "GSE31348_signature_score_table.csv",
##   comparisons = list(
##     TB_vs_Week1 = c("TB", "Week1"),
##     TB_vs_Week26 = c("TB", "Week26")
##   ),
##   score_col = "Score",
##   sample_col = "Sample",
##   group_col = "Group",
##   output_prefix = "GSE31348"
## )


cat("Validation plotting script loaded successfully.\n")