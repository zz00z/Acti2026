############################################################
## 06_MR_CASP4_VAMP5_Tuberculosis.R
##
## Two-sample Mendelian randomization analysis for candidate
## genes identified from machine-learning-derived signatures.
##
## Outcome:
##   ukb-b-15622
##
## Main analyses:
## 1. Forward MR analysis for candidate genes
## 2. Forest plot of selected MR results
## 3. Scatter plots for CASP4 and VAMP5
## 4. Single-SNP forest plots for CASP4 and VAMP5
##
## Candidate genes:
##   ANKRD22, FCGR1A, GBP1, GBP5, CASP4, GYG1,
##   LHFPL2, NPC2, PSTPIP2, VAMP5
############################################################


## =========================================================
## 0. Load packages
## =========================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(stringr)
  library(forestploter)
  library(grid)
  library(ggplot2)
})


## =========================================================
## 1. Set directories and outcome
## =========================================================

base_dir  <- "D:/atraining/1medicine/mtb"
cache_dir <- file.path(base_dir, "vcf_cache_txt")

## Outcome used in the final MR analysis
outcome_id  <- "ukb-b-15622"
outcome_txt <- file.path(cache_dir, paste0(outcome_id, ".txt"))

result_root <- file.path(
  base_dir,
  "forward_MR_results_all_genes",
  outcome_id
)

dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

## Local PLINK settings for LD clumping
plink_bin <- "D:/atraining/1medicine/clump/plink/plink.exe"
bfile     <- "D:/atraining/1medicine/clump/1kg.v3/EUR"


## =========================================================
## 2. Candidate gene exposure files
## =========================================================

gene_map <- data.frame(
  gene_txt = c(
    "eqtl-a-ENSG00000152766.txt",
    "eqtl-a-ENSG00000150337.txt",
    "eqtl-a-ENSG00000117228.txt",
    "eqtl-a-ENSG00000154451.txt",
    "eqtl-a-ENSG00000196954.txt",
    "eqtl-a-ENSG00000163754.txt",
    "eqtl-a-ENSG00000145685.txt",
    "eqtl-a-ENSG00000119655.txt",
    "eqtl-a-ENSG00000152229.txt",
    "eqtl-a-ENSG00000168899.txt"
  ),
  gene_symbol = c(
    "ANKRD22",
    "FCGR1A",
    "GBP1",
    "GBP5",
    "CASP4",
    "GYG1",
    "LHFPL2",
    "NPC2",
    "PSTPIP2",
    "VAMP5"
  ),
  stringsAsFactors = FALSE
)

candidate_genes <- gene_map$gene_symbol
priority_genes  <- c("CASP4", "VAMP5")


## =========================================================
## 3. MR parameters
## =========================================================

exp_p_threshold <- 1e-5
F_threshold     <- 10

clump_kb <- 5000
clump_r2 <- 0.01
clump_p  <- 1

method_order <- c(
  "Inverse variance weighted",
  "MR Egger",
  "Weighted median",
  "Weighted mode",
  "Simple mode",
  "Wald ratio"
)


## =========================================================
## 4. Helper functions
## =========================================================

safe_write_csv <- function(x, file) {
  tryCatch({
    write.csv(x, file, row.names = FALSE)
  }, error = function(e) {
    message("Failed to write file: ", file)
  })
}


make_exposure_from_txt <- function(df, exposure_name) {
  
  if ("beta.exposure" %in% colnames(df)) {
    
    out <- df
    
  } else if ("beta.outcome" %in% colnames(df)) {
    
    out <- df %>%
      rename(
        beta.exposure = beta.outcome,
        se.exposure = se.outcome,
        pval.exposure = pval.outcome,
        effect_allele.exposure = effect_allele.outcome,
        other_allele.exposure = other_allele.outcome,
        eaf.exposure = eaf.outcome,
        samplesize.exposure = samplesize.outcome,
        exposure = outcome,
        id.exposure = id.outcome
      )
    
  } else {
    
    stop("No exposure/outcome columns detected in exposure file.")
  }
  
  out$exposure <- exposure_name
  out$id.exposure <- exposure_name
  
  return(out)
}


make_outcome_from_txt <- function(df, outcome_name) {
  
  if ("beta.outcome" %in% colnames(df)) {
    
    out <- df
    
  } else if ("beta.exposure" %in% colnames(df)) {
    
    out <- df %>%
      rename(
        beta.outcome = beta.exposure,
        se.outcome = se.exposure,
        pval.outcome = pval.exposure,
        effect_allele.outcome = effect_allele.exposure,
        other_allele.outcome = other_allele.exposure,
        eaf.outcome = eaf.exposure,
        samplesize.outcome = samplesize.exposure,
        outcome = exposure,
        id.outcome = id.exposure
      )
    
  } else {
    
    stop("No exposure/outcome columns detected in outcome file.")
  }
  
  out$outcome <- outcome_name
  out$id.outcome <- outcome_name
  
  return(out)
}


local_clump_tsmr <- function(
    dat,
    plink_bin,
    bfile,
    clump_kb = 5000,
    clump_r2 = 0.01,
    clump_p = 1
) {
  
  clump_input <- tibble(
    rsid = dat$SNP,
    pval = dat$pval.exposure,
    id   = dat$id.exposure
  )
  
  clumped <- ieugwasr::ld_clump(
    dat = clump_input,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p  = clump_p,
    bfile = bfile,
    plink_bin = plink_bin
  )
  
  dat2 <- dat %>%
    filter(SNP %in% clumped$rsid)
  
  return(dat2)
}


## =========================================================
## 5. Run forward MR for all candidate genes
## =========================================================

run_forward_mr_all_genes <- function() {
  
  ## Read outcome once
  outcome_raw <- fread(
    outcome_txt,
    sep = "\t",
    data.table = FALSE
  )
  
  outcome_raw <- make_outcome_from_txt(
    outcome_raw,
    outcome_id
  )
  
  safe_write_csv(
    outcome_raw,
    file.path(result_root, "00_outcome_full_cached_as_outcome.csv")
  )
  
  all_summary <- list()
  
  for (i in seq_len(nrow(gene_map))) {
    
    gene_symbol <- gene_map$gene_symbol[i]
    exp_txt <- file.path(cache_dir, gene_map$gene_txt[i])
    
    cat("\n========================================\n")
    cat("Processing:", gene_symbol, "->", outcome_id, "\n")
    cat("========================================\n")
    
    result_dir <- file.path(result_root, gene_symbol)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    
    ## -------------------------
    ## 5.1 Read exposure file
    ## -------------------------
    
    exposure_raw <- tryCatch(
      fread(exp_txt, sep = "\t", data.table = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(exposure_raw) || nrow(exposure_raw) == 0) {
      cat("Failed to read exposure file for", gene_symbol, "\n")
      next
    }
    
    exposure_raw <- tryCatch(
      make_exposure_from_txt(exposure_raw, gene_symbol),
      error = function(e) NULL
    )
    
    if (is.null(exposure_raw) || nrow(exposure_raw) == 0) {
      cat("Failed to format exposure for", gene_symbol, "\n")
      next
    }
    
    ## -------------------------
    ## 5.2 Select exposure SNPs
    ## -------------------------
    
    exposure_dat <- exposure_raw %>%
      filter(!is.na(pval.exposure), pval.exposure < exp_p_threshold)
    
    if (nrow(exposure_dat) == 0) {
      cat("No exposure SNPs under p threshold for", gene_symbol, "\n")
      next
    }
    
    exposure_dat$F_stat <- (
      exposure_dat$beta.exposure / exposure_dat$se.exposure
    )^2
    
    exposure_dat <- exposure_dat %>%
      filter(F_stat > F_threshold)
    
    if (nrow(exposure_dat) == 0) {
      cat("No exposure SNPs after F-statistic filtering for", gene_symbol, "\n")
      next
    }
    
    safe_write_csv(
      exposure_dat,
      file.path(result_dir, "01_exposure_pfiltered_Ffiltered.csv")
    )
    
    ## -------------------------
    ## 5.3 LD clumping
    ## -------------------------
    
    exposure_clumped <- tryCatch(
      local_clump_tsmr(
        dat = exposure_dat,
        plink_bin = plink_bin,
        bfile = bfile,
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = clump_p
      ),
      error = function(e) NULL
    )
    
    if (is.null(exposure_clumped) || nrow(exposure_clumped) == 0) {
      cat("No SNPs after clumping for", gene_symbol, "\n")
      next
    }
    
    safe_write_csv(
      exposure_clumped,
      file.path(result_dir, "02_exposure_clumped.csv")
    )
    
    ## -------------------------
    ## 5.4 Match outcome SNPs
    ## -------------------------
    
    outcome_dat <- outcome_raw %>%
      filter(SNP %in% exposure_clumped$SNP)
    
    if (nrow(outcome_dat) == 0) {
      cat("No matched outcome SNPs for", gene_symbol, "\n")
      next
    }
    
    safe_write_csv(
      outcome_dat,
      file.path(result_dir, "03_outcome_matched.csv")
    )
    
    ## -------------------------
    ## 5.5 Harmonise data
    ## -------------------------
    
    harm_dat <- tryCatch(
      harmonise_data(
        exposure_dat = exposure_clumped,
        outcome_dat = outcome_dat,
        action = 2
      ),
      error = function(e) NULL
    )
    
    if (is.null(harm_dat) || nrow(harm_dat) == 0) {
      cat("Harmonisation failed for", gene_symbol, "\n")
      next
    }
    
    harm_dat <- subset(harm_dat, mr_keep == TRUE)
    
    if (nrow(harm_dat) == 0) {
      cat("No MR-eligible SNPs for", gene_symbol, "\n")
      next
    }
    
    safe_write_csv(
      harm_dat,
      file.path(result_dir, "04_harmonised_data.csv")
    )
    
    ## -------------------------
    ## 5.6 Main MR analysis
    ## -------------------------
    
    mr_res <- tryCatch(
      mr(harm_dat),
      error = function(e) NULL
    )
    
    if (is.null(mr_res) || nrow(mr_res) == 0) {
      cat("MR failed for", gene_symbol, "\n")
      next
    }
    
    mr_res$gene <- gene_symbol
    mr_res$outcome_id <- outcome_id
    
    safe_write_csv(
      mr_res,
      file.path(result_dir, "05_mr_beta.csv")
    )
    
    mr_or <- tryCatch(
      generate_odds_ratios(mr_res),
      error = function(e) NULL
    )
    
    if (!is.null(mr_or)) {
      safe_write_csv(
        mr_or,
        file.path(result_dir, "06_mr_or.csv")
      )
    }
    
    ## -------------------------
    ## 5.7 Sensitivity analyses
    ## -------------------------
    
    heter_res <- tryCatch(
      mr_heterogeneity(harm_dat),
      error = function(e) NULL
    )
    
    pleio_res <- tryCatch(
      mr_pleiotropy_test(harm_dat),
      error = function(e) NULL
    )
    
    single_res <- tryCatch(
      mr_singlesnp(harm_dat),
      error = function(e) NULL
    )
    
    loo_res <- tryCatch(
      mr_leaveoneout(harm_dat),
      error = function(e) NULL
    )
    
    if (!is.null(heter_res)) {
      safe_write_csv(
        heter_res,
        file.path(result_dir, "07_heterogeneity.csv")
      )
    }
    
    if (!is.null(pleio_res)) {
      safe_write_csv(
        pleio_res,
        file.path(result_dir, "08_pleiotropy.csv")
      )
    }
    
    if (!is.null(single_res)) {
      safe_write_csv(
        single_res,
        file.path(result_dir, "09_singlesnp.csv")
      )
    }
    
    if (!is.null(loo_res)) {
      safe_write_csv(
        loo_res,
        file.path(result_dir, "10_leaveoneout.csv")
      )
    }
    
    all_summary[[gene_symbol]] <- mr_res
    
    rm(
      exposure_raw,
      exposure_dat,
      exposure_clumped,
      outcome_dat,
      harm_dat,
      mr_res,
      mr_or,
      heter_res,
      pleio_res,
      single_res,
      loo_res
    )
    
    invisible(gc())
  }
  
  if (length(all_summary) > 0) {
    
    final_summary <- data.table::rbindlist(
      all_summary,
      fill = TRUE
    )
    
    safe_write_csv(
      final_summary,
      file.path(result_root, "00_all_genes_summary.csv")
    )
  }
  
  cat("\nForward MR analysis completed.\n")
  cat("Result root:\n")
  cat(result_root, "\n")
}


## Run MR analysis
## Uncomment this line when running for the first time:
## run_forward_mr_all_genes()


## =========================================================
## 6. Forest plot for MR results
## =========================================================
## Figure logic:
##   CASP4 and VAMP5: show all MR methods
##   Other genes: show IVW only

plot_selected_mr_forest <- function() {
  
  mr_summary_file <- file.path(result_root, "00_all_genes_summary.csv")
  
  if (!file.exists(mr_summary_file)) {
    stop("MR summary file not found. Please run run_forward_mr_all_genes() first.")
  }
  
  mr_res <- read.csv(
    mr_summary_file,
    check.names = FALSE
  )
  
  mr_or <- generate_odds_ratios(mr_res)
  
  if (!"gene" %in% colnames(mr_or)) {
    if ("exposure" %in% colnames(mr_or)) {
      mr_or$gene <- mr_or$exposure
    } else {
      stop("No gene or exposure column found in MR result.")
    }
  }
  
  plot_dat <- mr_or %>%
    filter(
      (
        gene %in% priority_genes &
          method %in% method_order
      ) |
        (
          !(gene %in% priority_genes) &
            method == "Inverse variance weighted"
        )
    ) %>%
    mutate(
      method = factor(method, levels = method_order),
      gene = factor(
        gene,
        levels = c(
          "CASP4",
          "VAMP5",
          "ANKRD22",
          "FCGR1A",
          "GBP1",
          "GBP5",
          "GYG1",
          "LHFPL2",
          "NPC2",
          "PSTPIP2"
        )
      )
    ) %>%
    arrange(gene, method) %>%
    mutate(
      OR_CI = sprintf("%.3f (%.3f - %.3f)", or, or_lci95, or_uci95),
      P_text = ifelse(
        pval < 0.001,
        formatC(pval, format = "e", digits = 2),
        sprintf("%.4f", pval)
      ),
      Method = as.character(method),
      SNP_N = nsnp
    )
  
  safe_write_csv(
    plot_dat,
    file.path(result_root, "MR_selected_results_for_forest.csv")
  )
  
  display_table <- plot_dat %>%
    transmute(
      Gene = as.character(gene),
      Method = Method,
      SNP_N = SNP_N,
      ` ` = strrep(" ", 20),
      `OR (95% CI)` = OR_CI,
      P = P_text
    )
  
  safe_write_csv(
    display_table,
    file.path(result_root, "MR_forest_display_table.csv")
  )
  
  xmin <- min(plot_dat$or_lci95, na.rm = TRUE)
  xmax <- max(plot_dat$or_uci95, na.rm = TRUE)
  
  if (!is.finite(xmin)) xmin <- 0.99
  if (!is.finite(xmax)) xmax <- 1.01
  
  ## Since the MR effect sizes in this analysis are small,
  ## a narrow x-axis range is used.
  xlim_use <- c(
    min(0.99, floor(xmin * 1000) / 1000),
    max(1.01, ceiling(xmax * 1000) / 1000)
  )
  
  theme_mr <- forest_theme(
    base_size = 10,
    ci_pch = 20,
    ci_col = "#4575b4",
    ci_lty = 1,
    ci_lwd = 2.3,
    ci_Theight = 0.2,
    refline_lwd = 2,
    refline_lty = "dashed",
    refline_col = "red",
    summary_fill = "#4575b4",
    summary_col = "#4575b4",
    footnote_cex = 1.0,
    footnote_fontface = "italic",
    footnote_col = "blue"
  )
  
  pdf(
    file.path(result_root, "MR_forest_selected_methods.pdf"),
    width = 12,
    height = max(6, 0.45 * nrow(display_table) + 2.5)
  )
  
  p <- forest(
    display_table[, c("Gene", "Method", "SNP_N", " ", "OR (95% CI)", "P")],
    est = plot_dat$or,
    lower = plot_dat$or_lci95,
    upper = plot_dat$or_uci95,
    sizes = 0.6,
    ci_column = 4,
    ref_line = 1,
    xlim = xlim_use,
    ticks_at = pretty(xlim_use, n = 5),
    arrow_lab = c("protective factor", "risk factor"),
    footnote = "CASP4 and VAMP5 are shown using all MR methods; other genes are shown using IVW.",
    theme = theme_mr
  )
  
  print(p)
  dev.off()
  
  cat("MR forest plot completed.\n")
}


## Draw forest plot
## plot_selected_mr_forest()


## =========================================================
## 7. Scatter and single-SNP forest plots for CASP4 and VAMP5
## =========================================================

plot_priority_gene_mr_plots <- function(
    genes_to_plot = priority_genes
) {
  
  plot_dir <- file.path(
    result_root,
    "MR_plots_CASP4_VAMP5"
  )
  
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (g in genes_to_plot) {
    
    cat("Processing MR plots for:", g, "\n")
    
    gene_dir <- file.path(result_root, g)
    
    harm_file   <- file.path(gene_dir, "04_harmonised_data.csv")
    mr_file     <- file.path(gene_dir, "05_mr_beta.csv")
    single_file <- file.path(gene_dir, "09_singlesnp.csv")
    
    if (!file.exists(harm_file)) {
      message("Missing harmonised data for ", g)
      next
    }
    
    if (!file.exists(mr_file)) {
      message("Missing MR result for ", g)
      next
    }
    
    harm_dat <- fread(
      harm_file,
      data.table = FALSE
    )
    
    mr_res <- fread(
      mr_file,
      data.table = FALSE
    )
    
    mr_res$method <- factor(
      mr_res$method,
      levels = method_order
    )
    
    ## -------------------------
    ## MR scatter plot
    ## -------------------------
    
    scatter_list <- tryCatch(
      mr_scatter_plot(mr_res, harm_dat),
      error = function(e) NULL
    )
    
    if (!is.null(scatter_list) && length(scatter_list) > 0) {
      
      pdf(
        file.path(plot_dir, paste0(g, "_MR_scatter_", outcome_id, ".pdf")),
        width = 7,
        height = 7
      )
      
      for (p in scatter_list) {
        print(
          p +
            ggtitle(paste0("MR effect size for '", g, "' on '", outcome_id, "'")) +
            theme_bw(base_size = 12)
        )
      }
      
      dev.off()
    }
    
    ## -------------------------
    ## Single-SNP MR forest plot
    ## -------------------------
    
    if (file.exists(single_file)) {
      
      single_res <- fread(
        single_file,
        data.table = FALSE
      )
      
    } else {
      
      single_res <- tryCatch(
        mr_singlesnp(harm_dat),
        error = function(e) NULL
      )
    }
    
    if (!is.null(single_res) && nrow(single_res) > 0) {
      
      safe_write_csv(
        single_res,
        file.path(plot_dir, paste0(g, "_singleSNP_results_", outcome_id, ".csv"))
      )
      
      single_forest_list <- tryCatch(
        mr_forest_plot(single_res),
        error = function(e) NULL
      )
      
      if (!is.null(single_forest_list) && length(single_forest_list) > 0) {
        
        pdf(
          file.path(plot_dir, paste0(g, "_singleSNP_forest_", outcome_id, ".pdf")),
          width = 7,
          height = 7
        )
        
        for (p in single_forest_list) {
          print(
            p +
              ggtitle(paste0(g, " single-SNP MR on ", outcome_id)) +
              theme_bw(base_size = 12)
          )
        }
        
        dev.off()
      }
    }
  }
  
  cat("CASP4 and VAMP5 MR plots completed.\n")
}


## Draw scatter and single-SNP forest plots
## plot_priority_gene_mr_plots()


## =========================================================
## 8. Optional leave-one-out plots
## =========================================================

plot_priority_gene_leaveoneout <- function(
    genes_to_plot = priority_genes
) {
  
  outdir <- file.path(result_root, "MR_leaveoneout_CASP4_VAMP5")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  for (g in genes_to_plot) {
    
    gene_dir <- file.path(result_root, g)
    harm_file <- file.path(gene_dir, "04_harmonised_data.csv")
    
    if (!file.exists(harm_file)) {
      message("Missing harmonised data for ", g)
      next
    }
    
    harm_dat <- fread(
      harm_file,
      data.table = FALSE
    )
    
    loo <- tryCatch(
      mr_leaveoneout(harm_dat),
      error = function(e) NULL
    )
    
    if (is.null(loo) || nrow(loo) == 0) next
    
    p <- tryCatch(
      mr_leaveoneout_plot(loo),
      error = function(e) NULL
    )
    
    if (!is.null(p) && length(p) > 0) {
      
      pdf(
        file.path(outdir, paste0(g, "_leaveoneout_", outcome_id, ".pdf")),
        width = 7,
        height = 7
      )
      
      print(
        p[[1]] +
          ggtitle(paste0(g, " leave-one-out analysis on ", outcome_id)) +
          theme_bw(base_size = 12)
      )
      
      dev.off()
    }
  }
  
  cat("Leave-one-out plots completed.\n")
}


## Optional:
## plot_priority_gene_leaveoneout()


cat("MR script loaded successfully.\n")