# Seven-gene, two-outcome Mendelian randomization.
# Inputs: eQTL/outcome summary tables, EUR LD reference and PLINK; acti_mr_config().
# Run: Rscript 6_Mendelian_randomization.R --stage instruments
# source() loads functions without starting an analysis.
# Primary P < 1e-5; sensitivity P < 5e-8; F > 10; EUR r2=0.01, 5000 kb.

acti_mr_config <- function(data_dir = file.path("data", "MR"), output_dir = file.path("results", "MR"), result_dir = NULL,
  cache_dir = NULL, plink = "plink", ld_bfile = file.path("data", "LD", "EUR"), save_plots = TRUE, overwrite = FALSE) {
  list(data_dir = data_dir, output_dir = output_dir, result_dir = if (is.null(result_dir)) output_dir else result_dir,
    cache_dir = if (is.null(cache_dir)) file.path(data_dir, "outcome_cache") else cache_dir, plink = plink,
    ld_bfile = ld_bfile, save_plots = isTRUE(save_plots), overwrite = isTRUE(overwrite))
}

acti_mr_initialize <- function(cfg = acti_mr_config()) {
  dir.create(cfg$result_dir, recursive = TRUE, showWarnings = FALSE)
  genes <- data.frame(gene = c("DUSP3", "FCGR1A", "NPC2", "TAP1", "ANKRD22", "GBP5", "VAMP5"), ensembl_id = c("ENSG00000108861",
    "ENSG00000150337", "ENSG00000119655", "ENSG00000168394", "ENSG00000152766", "ENSG00000154451", "ENSG00000168899"),
    opengwas_id = c("eqtl-a-ENSG00000108861", "eqtl-a-ENSG00000150337", "eqtl-a-ENSG00000119655", "eqtl-a-ENSG00000168394",
      "eqtl-a-ENSG00000152766", "eqtl-a-ENSG00000154451", "eqtl-a-ENSG00000168899"), stringsAsFactors = FALSE)
  write.csv(genes, file.path(cfg$result_dir, "gene_manifest.csv"), row.names = FALSE)
  writeLines("[{\"id\": \"ukb-b-15622\", \"trait\": \"Non-cancer illness code, self-reported: tuberculosis (tb)\", \"population\": \"European\", \"sample_size\": 462933, \"ncase\": 2277, \"ncontrol\": 460656, \"year\": 2018, \"unit\": \"SD\", \"nsnp\": 9851867, \"build\": \"HG19/GRCh37\"}, {\"id\": \"finn-b-AB1_RESP_TUBERCU_CONF\", \"trait\": \"Respiratory tuberculosis, bacteriologically and histologically confirmed\", \"build\": \"HG19/GRCh37\", \"group_name\": \"public\", \"category\": \"Binary\", \"subcategory\": \"NA\", \"population\": \"European\", \"sex\": \"Males and Females\", \"author\": \"NA\", \"nsnp\": 16380466, \"year\": 2021, \"ontology\": \"NA\", \"unit\": \"NA\", \"ncase\": 247, \"ncontrol\": 217599, \"consortium\": \"NA\", \"mr\": 1, \"priority\": 0, \"note\": \"AB1_RESP_TUBERCU_CONF\"}]",
    file.path(cfg$result_dir, "outcome_metadata.json"))
  write.csv(data.frame(dataset = c("ukb-b-15622", "finn-b-AB1_RESP_TUBERCU_CONF"), LD_population = "EUR"), file.path(cfg$result_dir,
    "outcome_instrument_status.csv"), row.names = FALSE)
  invisible(genes)
}

acti_mr_instruments <- function(cfg = acti_mr_config()) {
  suppressPackageStartupMessages(library(data.table))
  setDTthreads(2)
  options(stringsAsFactors = FALSE)
  Sys.setenv(OMP_NUM_THREADS = "2")
  base <- cfg$data_dir
  root <- cfg$result_dir
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "cache"), showWarnings = FALSE)
  dir.create(file.path(root, "gene_instruments"), showWarnings = FALSE)
  plink <- cfg$plink
  bfile <- cfg$ld_bfile
  manifest <- fread(file.path(root, "gene_manifest.csv"))
  all_snps <- character()
  qc <- list()
  for (i in seq_len(nrow(manifest))) {
    g <- manifest$gene[i]
    id <- manifest$opengwas_id[i]
    message("PREPARING ", g)
    path <- file.path(base, "exposure", "gene_eQTL_tsv", paste0(g, "__", id, ".tsv.gz"))
    d <- as.data.frame(fread(path))
    d <- d[d$usable_for_mr == 1 & is.finite(d$beta) & is.finite(d$se) & d$se > 0 & is.finite(d$pval) & d$pval <
      1e-5 & (d$beta/d$se)^2 > 10 & grepl("^rs[0-9]+$", d$SNP), , drop = FALSE]
    d$pval <- pmax(d$pval, 1e-200)
    d <- d[!duplicated(d$SNP), , drop = FALSE]
    d$F_stat <- (d$beta/d$se)^2
    gd <- file.path(root, "gene_instruments", g)
    dir.create(gd, showWarnings = FALSE)
    fwrite(d, file.path(gd, "01_p_F_filtered.tsv"), sep = "\t")
    selected <- d[FALSE, , drop = FALSE]
    state <- "no_significant_valid_instrument"
    if (nrow(d)) {
      prefix <- file.path(gd, "plink_EUR")
      inp <- data.frame(SNP = d$SNP, P = pmax(d$pval, 1e-200))
      fwrite(inp, paste0(prefix, "_input.tsv"), sep = "\t")
      a <- c("--bfile", shQuote(bfile), "--clump", shQuote(paste0(prefix, "_input.tsv")), "--clump-p1", "1",
        "--clump-r2", "0.01", "--clump-kb", "5000", "--threads", "2", "--memory", "2048", "--out", shQuote(prefix))
      status <- system2(plink, a, stdout = paste0(prefix, "_stdout.log"), stderr = paste0(prefix, "_stderr.log"))
      cf <- paste0(prefix, ".clumped")
      if (status != 0) {
        state <- paste0("plink_error_exit_", status)
      }
      else if (file.exists(cf)) {
        cdat <- fread(cf, fill = TRUE)
        if ("SNP" %in% names(cdat))
         selected <- d[d$SNP %in% cdat$SNP, , drop = FALSE]
        state <- if (nrow(selected))
         "ready"
        else "no_instrument_in_LD_reference"
      }
      else state <- "no_instrument_in_LD_reference"
    }
    fwrite(selected, file.path(gd, "02_clumped.tsv"), sep = "\t")
    all_snps <- union(all_snps, selected$SNP)
    qc[[g]] <- data.frame(gene = g, opengwas_id = id, p_F_filtered = nrow(d), clumped = nrow(selected), status = state)
    message("GENE ", g, " significant=", nrow(d), " independent=", nrow(selected), " ", state)
  }
  fwrite(rbindlist(qc), file.path(root, "gene_instrument_status.csv"))
  writeLines(sort(all_snps), file.path(root, "cache", "forward_instrument_snps.txt"))
  saveRDS(list(input_base = base, p_threshold = 1e-5, F_threshold = 10, clump_r2 = 0.01, clump_kb = 5000,
    LD_reference = bfile, plink = plink, scope = "all available eQTL loci, following prior 06 script"), file.path(root,
    "analysis_parameters.rds"))
  capture.output(sessionInfo(), file = file.path(root, "R_sessionInfo.txt"))
  message("FORWARD UNION SNPs=", length(all_snps))
  if (any(grepl("^plink_error", rbindlist(qc)$status)))
    stop("PLINK failed. Check the executable, reference files and per-gene logs before running MR.")
}

acti_mr_forward <- function(cfg = acti_mr_config()) {
  suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(data.table)
    library(jsonlite)
  })
  data.table::setDTthreads(2)
  options(stringsAsFactors = FALSE)
  input_root <- cfg$result_dir
  root <- cfg$result_dir
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  base <- cfg$data_dir
  genes <- as.data.frame(fread(file.path(input_root, "gene_manifest.csv")))
  meta <- as.data.frame(fromJSON(file.path(input_root, "outcome_metadata.json")))
  all_dataset_ids <- meta$id
  assigned <- ""
  if (nzchar(assigned))
    meta <- meta[meta$id %in% strsplit(assigned, ",", fixed = TRUE)[[1]], , drop = FALSE]
  gene_qc <- fread(file.path(input_root, "gene_instrument_status.csv"))
  tb_qc <- fread(file.path(input_root, "outcome_instrument_status.csv"))
  dir.create(file.path(input_root, "analyses"), showWarnings = FALSE)
  methods <- c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
  outcome_scale <- setNames(rep("log_odds", nrow(meta)), meta$id)
  outcome_scale["ukb-b-15622"] <- "UKB_SD_scale_no_OR"
  outcome_scale["ebi-a-GCST90038710"] <- "raw_beta_unverified_scale_no_OR"
  write_json(list(min_pval = 1e-200, clump_p1 = 1, clump_p2 = "PLINK default 0.01, following prior ieugwasr::ld_clump",
    clump_input_order = "original VCF order", primary_method = "IVW; single SNP Wald reported separately",
    date = "2026-10-14", genes = genes$gene, outcomes = meta$id, primary_p_threshold = 1e-5,
    strict_sensitivity_p_threshold = 5e-8, F_threshold = 10, clump_r2 = 0.01, clump_kb = 5000,
    harmonisation_action = 2, proxy_lookup = FALSE, primary_instruments = "all available eQTL loci following supplied 06 script; not cis restricted",
    EUR_reference = cfg$ld_bfile, EAS_reference = "OpenGWAS 1000 Genomes EAS API", missing_ancestry_assumption = "GCST90038710: EUR assumed from UK Biobank origin; ancestry metadata NA",
    outcome_effect_scale = as.list(outcome_scale), multiple_testing = "BH within 7 genes per outcome and across all 14 prespecified gene-outcome pairs for each analysis threshold",
    seed = 20261014), file.path(root, "analysis_config.json"), auto_unbox = TRUE, pretty = TRUE)
  clean <- function(d) {
    d <- as.data.frame(d)
    needed <- c("SNP", "beta", "se", "effect_allele", "other_allele", "pval", "eaf", "samplesize")
    if (!all(needed %in% names(d)))
      stop("Missing fields: ", paste(setdiff(needed, names(d)), collapse = ","))
    d <- d[is.finite(d$beta) & is.finite(d$se) & d$se > 0 & is.finite(d$pval) & d$pval >= 0 & grepl("^rs[0-9]+$",
      d$SNP), , drop = FALSE]
    d <- d[order(d$pval, d$se), , drop = FALSE]
    d[!duplicated(d$SNP), , drop = FALSE]
  }
  format_trait <- function(d, type, name, id) {
    d <- clean(d)
    d$phenotype <- name
    d$dataset_id <- id
    d$pval <- pmax(d$pval, 1e-200)
    TwoSampleMR::format_data(d, type = type, phenotype_col = "phenotype", snp_col = "SNP", beta_col = "beta",
      se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
      pval_col = "pval", samplesize_col = "samplesize", id_col = "dataset_id", chr_col = "chr", pos_col = "pos",
      min_pval = 1e-200)
  }
  gene_full <- list()
  gene_ivs <- list()
  tb_full <- list()
  tb_ivs <- list()
  for (i in seq_len(nrow(genes))) {
    g <- genes$gene[i]
    id <- genes$opengwas_id[i]
    message("LOAD_GENE ", g)
    d <- as.data.frame(fread(file.path(base, "exposure", "gene_eQTL_tsv", paste0(g, "__", id, ".tsv.gz"))))
    gene_full[[g]] <- clean(d[d$usable_for_mr == 1, , drop = FALSE])
    gene_ivs[[g]] <- as.data.frame(fread(file.path(input_root, "gene_instruments", g, "02_clumped.tsv")))
  }
  for (id in meta$id) {
    message("LOAD_TB ", id)
    tb_full[[id]] <- clean(fread(file.path(input_root, "cache", paste0(id, "__forward_matches.tsv.gz"))))
  }
  results <- list()
  statuses <- list()
  sensitivities <- list()
  warnings_all <- list()
  counter <- 0L
  run_pair <- function(g, id, direction, analysis) {
    counter <<- ((match(analysis, c("legacy_p1e5", "strict_p5e8")) - 1L) * 320L + (match(direction, c("forward",
      "reverse")) - 1L) * 160L + (match(id, all_dataset_ids) - 1L) * 20L + match(g, genes$gene))
    threshold <- if (analysis == "legacy_p1e5")
      1e-5
    else 5e-8
    ddir <- file.path(input_root, "analyses", analysis, direction, id, g)
    dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
    note <- character()
    st <- data.frame(gene = g, dataset = id, direction = direction, analysis = analysis, p_threshold = threshold,
      exposure_population = if (direction == "forward")
        "European"
      else meta$population[match(id, meta$id)], outcome_population = if (direction == "forward")
        meta$population[match(id, meta$id)]
      else "European", cross_ancestry = id %in% c("bbj-a-149", "ebi-a-GCST90018672"), ancestry_assumed = id ==
        "ebi-a-GCST90038710", LD_population = if (direction == "forward")
        "EUR"
      else tb_qc$LD_population[match(id, tb_qc$dataset)], instruments_before_match = 0L, matched = 0L, harmonised = 0L,
      min_F = NA_real_, median_F = NA_real_, status = "started", detail = "")
    cached_status <- file.path(ddir, "00_status.csv")
    cached_result <- file.path(ddir, "05_MR_results.csv")
    if (file.exists(cached_status)) {
      old <- as.data.frame(fread(cached_status))
      if (nrow(old) == 1 && old$status == "success" && file.exists(cached_result)) {
        r <- as.data.frame(fread(cached_result))
        if (nrow(r) && all(c("b", "se", "pval", "gene", "dataset", "direction", "analysis") %in% names(r))) {
         results[[paste(analysis, direction, id, g, sep = "__")]] <<- r
         for (test in c("heterogeneity", "pleiotropy", "singleSNP", "leaveoneout")) {
          cf <- file.path(ddir, paste0("06_", test, ".csv"))
          if (file.exists(cf)) {
           z <- as.data.frame(fread(cf))
           z$gene <- g
           z$dataset <- id
           z$direction <- direction
           z$analysis <- analysis
           sensitivities[[paste(test, direction, id, g, sep = "__")]] <<- list(test = test, data = z)
          }
         }
         return(old[, names(st), drop = FALSE])
        }
      }
    }
    tryCatch(withCallingHandlers({
      if (direction == "forward") {
        expraw <- gene_ivs[[g]]
        outraw <- tb_full[[id]]
        expname <- paste0(g, " expression")
        expid <- genes$opengwas_id[match(g, genes$gene)]
        outname <- meta$trait[match(id, meta$id)]
        outid <- id
      }
      else {
        expraw <- tb_ivs[[id]]
        outraw <- gene_full[[g]]
        expname <- meta$trait[match(id, meta$id)]
        expid <- id
        outname <- paste0(g, " expression")
        outid <- genes$opengwas_id[match(g, genes$gene)]
      }
      expraw <- expraw[is.finite(expraw$pval) & expraw$pval < threshold, , drop = FALSE]
      st$instruments_before_match <- nrow(expraw)
      if (!nrow(expraw)) {
        st$status <- "no_exposure_instruments"
        return(st)
      }
      outraw <- outraw[outraw$SNP %in% expraw$SNP, , drop = FALSE]
      st$matched <- length(intersect(expraw$SNP, outraw$SNP))
      fwrite(expraw, file.path(ddir, "01_exposure_instruments.tsv"), sep = "\t")
      fwrite(outraw, file.path(ddir, "02_matched_outcome.tsv"), sep = "\t")
      coverage <- data.frame(SNP = expraw$SNP, available_in_outcome = expraw$SNP %in% outraw$SNP)
      fwrite(coverage, file.path(ddir, "00_SNP_coverage.csv"))
      if (!nrow(outraw)) {
        st$status <- "no_matched_outcome_SNP"
        return(st)
      }
      exp <- format_trait(expraw, "exposure", expname, expid)
      out <- format_trait(outraw, "outcome", outname, outid)
      m <- meta[match(id, meta$id), , drop = FALSE]
      if (direction == "forward") {
        out$ncase.outcome <- m$ncase
        out$ncontrol.outcome <- m$ncontrol
      }
      else {
        exp$ncase.exposure <- m$ncase
        exp$ncontrol.exposure <- m$ncontrol
      }
      h <- harmonise_data(exp, out, action = 2)
      fwrite(h, file.path(ddir, "03_harmonised_all.csv"))
      h <- h[h$mr_keep %in% TRUE, , drop = FALSE]
      st$harmonised <- nrow(h)
      if (!nrow(h)) {
        st$status <- "no_SNP_after_harmonisation"
        return(st)
      }
      h$F_stat <- (h$beta.exposure/h$se.exposure)^2
      st$min_F <- min(h$F_stat)
      st$median_F <- median(h$F_stat)
      fwrite(h, file.path(ddir, "04_harmonised_keep.csv"))
      set.seed(20261014 + counter)
      used <- if (nrow(h) == 1)
        "mr_wald_ratio"
      else "mr_ivw"
      res <- mr(h, method_list = used)
      res <- res[is.finite(res$b) & is.finite(res$se) & res$se > 0 & is.finite(res$pval), , drop = FALSE]
      if (!nrow(res)) {
        st$status <- "MR_not_estimable"
        return(st)
      }
      res$gene <- g
      res$dataset <- id
      res$direction <- direction
      res$analysis <- analysis
      res$is_primary <- res$method %in% c("Inverse variance weighted", "Wald ratio")
      res$outcome_effect_scale <- if (direction == "forward")
        outcome_scale[id]
      else "gene_expression_source_scale"
      res$beta_lci95 <- res$b - 1.96 * res$se
      res$beta_uci95 <- res$b + 1.96 * res$se
      res$or <- res$or_lci95 <- res$or_uci95 <- NA_real_
      if (direction == "forward" && outcome_scale[id] == "log_odds") {
        res$or <- exp(res$b)
        res$or_lci95 <- exp(res$beta_lci95)
        res$or_uci95 <- exp(res$beta_uci95)
      }
      res$cross_ancestry <- st$cross_ancestry
      res$ancestry_assumed <- st$ancestry_assumed
      res$min_F <- st$min_F
      res$median_F <- st$median_F
      fwrite(res, file.path(ddir, "05_MR_results.csv"))
      results[[paste(analysis, direction, id, g, sep = "__")]] <<- res
      if (analysis == "legacy_p1e5") {
        calls <- list(heterogeneity = function() if (nrow(h) >= 2) mr_heterogeneity(h, method_list = if (nrow(h) >=
         3) c("mr_ivw", "mr_egger_regression") else "mr_ivw") else NULL, pleiotropy = function() if (nrow(h) >=
         3) mr_pleiotropy_test(h) else NULL, singleSNP = function() mr_singlesnp(h), leaveoneout = function() if (nrow(h) >=
         3) mr_leaveoneout(h) else NULL)
        for (test in names(calls)) {
         z <- tryCatch(calls[[test]](), error = function(e) {
          note <<- c(note, paste(test, conditionMessage(e), sep = ": "))
          NULL
         })
         if (!is.null(z) && nrow(z)) {
          fwrite(z, file.path(ddir, paste0("06_", test, ".csv")))
          z$gene <- g
          z$dataset <- id
          z$direction <- direction
          z$analysis <- analysis
          sensitivities[[paste(test, direction, id, g, sep = "__")]] <<- list(test = test, data = z)
         }
        }
      }
      st$status <- "success"
    }, warning = function(w) {
      note <<- c(note, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      st$status <<- "error"
      st$detail <<- conditionMessage(e)
    })
    if (length(note))
      st$detail <- paste(c(st$detail, unique(note)), collapse = " | ")
    st
  }
  for (analysis in c("legacy_p1e5", "strict_p5e8")) {
    for (direction in c("forward")) {
      for (id in meta$id) {
        for (g in genes$gene) {
         key <- paste(analysis, direction, id, g, sep = "__")
         s <- run_pair(g, id, direction, analysis)
         statuses[[key]] <- s
         fwrite(s, file.path(input_root, "analyses", analysis, direction, id, g, "00_status.csv"))
         message("MR ", analysis, " ", direction, " ", id, " ", g, " n=", s$harmonised, " ", s$status)
        }
        fwrite(rbindlist(statuses, fill = TRUE), file.path(root, "analysis_status.csv"))
      }
    }
  }
  full <- as.data.frame(rbindlist(results, fill = TRUE))
  primary <- full[full$is_primary, , drop = FALSE]
  primary$q_BH_7_genes <- primary$q_BH_14_pairs <- NA_real_
  for (a in unique(primary$analysis)) for (d in unique(primary$direction)) {
    ix <- which(primary$analysis == a & primary$direction == d)
    primary$q_BH_14_pairs[ix] <- p.adjust(primary$pval[ix], "BH", n = 14)
    for (id in meta$id) {
      ii <- ix[primary$dataset[ix] == id]
      if (length(ii))
        primary$q_BH_7_genes[ii] <- p.adjust(primary$pval[ii], "BH", n = 7)
    }
  }
  fwrite(full, file.path(root, "00_all_MR_methods.csv"))
  fwrite(primary, file.path(root, "00_primary_IVW_Wald_results.csv"))
  for (a in c("legacy_p1e5", "strict_p5e8")) for (d in c("forward")) {
    fwrite(primary[primary$analysis == a & primary$direction == d, , drop = FALSE], file.path(root, paste0("00_",
      a, "_", d, "_primary.csv")))
  }
  for (test in c("heterogeneity", "pleiotropy", "singleSNP", "leaveoneout")) {
    ds <- lapply(sensitivities, function(z) if (z$test == test)
      z$data
    else NULL)
    if (any(lengths(ds) > 0))
      fwrite(rbindlist(ds, fill = TRUE), file.path(root, paste0("00_", test, "_all.csv")))
  }
  fwrite(rbindlist(statuses, fill = TRUE), file.path(root, "analysis_status.csv"))
  saveRDS(list(results = full, primary = primary, status = rbindlist(statuses), configuration = read_json(file.path(root,
    "analysis_config.json"))), file.path(root, "MR_analysis_bundle.rds"))
  capture.output(sessionInfo(), file = file.path(root, "R_sessionInfo.txt"))
  message("COMPLETE: ", nrow(primary), " primary results; ", length(statuses), " analysis attempts")
}

acti_mr_forest <- function(cfg = acti_mr_config()) {
  try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)
  suppressPackageStartupMessages({
    library(dplyr)
    library(forestploter)
    library(grid)
  })
  stopifnot(requireNamespace("TwoSampleMR", quietly = TRUE))
  base <- cfg$result_dir
  out <- file.path(cfg$output_dir, "MR_forest_selected_methods")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  primary <- read.csv(file.path(base, "01_main_IVW_14pairs.csv"), check.names = FALSE)
  primary <- primary[primary$analysis == "legacy_p1e5" & primary$method == "Inverse variance weighted", ]
  stopifnot(nrow(primary) == 14, !anyDuplicated(primary[c("gene", "dataset")]))
  priority <- primary[is.finite(primary$pval) & primary$pval < 0.05, ]
  method_order <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode", "Simple mode")
  methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_simple_mode")
  expanded <- list()
  checks <- list()
  set.seed(20261014)
  for (i in seq_len(nrow(priority))) {
    pr <- priority[i, ]
    f <- file.path(base, "analyses", "legacy_p1e5", "forward", pr$dataset, pr$gene, "04_harmonised_keep.csv")
    h <- read.csv(f, check.names = FALSE)
    stopifnot(all(as.logical(h$mr_keep)))
    res <- TwoSampleMR::mr(h, method_list = methods)
    iv <- res[res$method == "Inverse variance weighted", ]
    stopifnot(nrow(iv) == 1, abs(iv$b - pr$b) < 1e-10, abs(iv$se - pr$se) < 1e-10, abs(iv$pval - pr$pval) <
      1e-8)
    checks[[i]] <- data.frame(Gene = pr$gene, Dataset = pr$dataset, IVW_beta_difference = iv$b - pr$b, IVW_se_difference = iv$se -
      pr$se, IVW_p_difference = iv$pval - pr$pval)
    res$gene <- pr$gene
    res$dataset <- pr$dataset
    expanded[[i]] <- res
  }
  extra <- bind_rows(expanded)
  extra$beta <- extra$b
  extra$beta_lci95 <- extra$b - 1.96 * extra$se
  extra$beta_uci95 <- extra$b + 1.96 * extra$se
  write.csv(extra, file.path(out, "Recomputed_all_methods_nominal_IVW_pairs.csv"), row.names = FALSE)
  write.csv(bind_rows(checks), file.path(out, "IVW_reproduction_check.csv"), row.names = FALSE)
  dat <- primary
  for (i in seq_len(nrow(priority))) {
    pr <- priority[i, ]
    ad <- extra[extra$gene == pr$gene & extra$dataset == pr$dataset & extra$method != "Inverse variance weighted",
      ]
    dat <- bind_rows(dat, ad)
  }
  dat$beta_lower <- dat$b - 1.96 * dat$se
  dat$beta_upper <- dat$b + 1.96 * dat$se
  dat$beta <- dat$b
  dat <- dat[, setdiff(names(dat), c("or", "or_lci95", "or_uci95"))]
  dat$expanded_pair <- paste(dat$gene, dat$dataset) %in% paste(priority$gene, priority$dataset)
  dat$method <- factor(dat$method, levels = method_order)
  format_p <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.4f",
    p)))
  for (ds in unique(primary$dataset)) {
    pd <- dat[dat$dataset == ds, ]
    genes <- c(priority$gene[priority$dataset == ds], setdiff(c("DUSP3", "FCGR1A", "NPC2", "TAP1", "ANKRD22",
      "GBP5", "VAMP5"), priority$gene[priority$dataset == ds]))
    pd <- pd[order(match(pd$gene, genes), pd$method), ]
    is_or <- grepl("^finn", ds)
    pd$estimate <- pd$b
    pd$lower <- pd$beta_lower
    pd$upper <- pd$beta_upper
    pd$effect_scale <- if (is_or)
      "Beta (log odds)"
    else "Beta (UKB source scale)"
    digits <- if (is_or)
      3
    else 5
    display <- data.frame(Gene = ifelse(duplicated(pd$gene), "", pd$gene), Method = as.character(pd$method),
      SNP_N = pd$nsnp, space = strrep(" ", 26), effect = sprintf(paste0("%.", digits, "f (%.", digits, "f - %.",
        digits, "f)"), pd$estimate, pd$lower, pd$upper), P = format_p(pd$pval), check.names = FALSE)
    names(display)[4] <- " "
    names(display)[5] <- "Beta (95% CI)"
    names(display)[6] <- "P value"
    ref <- 0
    limits <- range(c(pd$lower, pd$upper, ref), finite = TRUE)
    pad <- diff(limits) * 0.07
    limits <- limits + c(-pad, pad)
    ticks <- if (is_or)
      pretty(limits, n = 3)
    else c(-0.004, 0, 0.004)
    ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
    ft <- forest_theme(base_size = 12, ci_pch = 20, ci_col = "#4575b4", ci_lwd = 2.3, ci_Theight = 0.2,
      refline_lwd = 2, refline_lty = "dashed", refline_col = "red", footnote_gp = gpar(fontsize = 10, col = "black"),
      core = list(padding = unit(c(4, 3), "mm")), colhead = list(fg_params = list(fontface = "bold")))
    foot <- paste0("Primary instruments: P < 1e-5; 95% CI. All methods shown for ", paste(priority$gene[priority$dataset ==
      ds], collapse = ", "), " (nominal IVW P < 0.05); others: IVW only.\nIVW q: BH across 14 gene-outcome pairs; no IVW association survives BH correction.",
      if (!is_or)
        "\nUKB effects are retained on the source beta scale; they are not odds ratios."
      else "")
    p <- forestploter::forest(display, est = pd$estimate, lower = pd$lower, upper = pd$upper, sizes = 0.5,
      ci_column = 4, ref_line = ref, xlim = limits, ticks_at = ticks, xlab = "Beta", theme = ft)
    stem <- file.path(out, paste0("MR_forest_", if (is_or)
      "FinnGen_confirmed_TB_beta"
    else "UKB_self_reported_TB_beta"))
    w <- convertWidth(sum(p$widths), "in", valueOnly = TRUE) + 0.35
    ht <- convertHeight(sum(p$heights), "in", valueOnly = TRUE) + 0.6
    if (cfg$save_plots) {
      cairo_pdf(paste0(stem, ".pdf"), width = w, height = ht, family = "Arial")
      grid.newpage()
      grid.draw(p)
      dev.off()
      ragg::agg_png(paste0(stem, ".png"), width = w, height = ht, units = "in", res = 300, background = "white")
      grid.newpage()
      grid.draw(p)
      dev.off()
    }
    write.csv(pd, paste0(stem, "_source_data.csv"), row.names = FALSE)
    write.csv(display, paste0(stem, "_display_table.csv"), row.names = FALSE)
    writeLines(foot, paste0(stem, "_caption.txt"))
  }
  write.csv(dat, file.path(out, "Selected_forest_all_source_data.csv"), row.names = FALSE)
  write.csv(dat[, c("gene", "dataset", "method", "nsnp", "beta", "se", "beta_lower", "beta_upper", "pval", "q_BH_14_pairs",
    "q_BH_7_genes")], file.path(out, "MR_selected_results_BETA.csv"), row.names = FALSE)
  cat("OUTPUT:", out, "\n")
  print(priority[, c("gene", "dataset", "pval", "q_BH_14_pairs")])
}

acti_mr_cache <- function(cfg = acti_mr_config()) {
  suppressPackageStartupMessages(library(data.table))
  dir.create(file.path(cfg$result_dir, "cache"), recursive = TRUE, showWarnings = FALSE)
  targets <- readLines(file.path(cfg$result_dir, "cache", "forward_instrument_snps.txt"))
  ids <- c("ukb-b-15622", "finn-b-AB1_RESP_TUBERCU_CONF")
  for (id in ids) {
    name <- paste0(id, "__forward_matches.tsv.gz")
    source <- file.path(cfg$cache_dir, name)
    dest <- file.path(cfg$result_dir, "cache", name)
    if (normalizePath(source, winslash = "/", mustWork = FALSE) == normalizePath(dest, winslash = "/", mustWork = FALSE)) {
      stop("Source outcome cache and generated output must use separate directories.")
    }
    if (file.exists(dest) && !cfg$overwrite)
      stop("Cache exists: ", dest, ". Use a new folder or --overwrite.")
    d <- fread(source)
    required <- c("SNP", "chr", "pos", "beta", "se", "effect_allele", "other_allele", "pval", "eaf", "samplesize")
    if (!all(required %in% names(d)))
      stop("Missing summary fields: ", id)
    d <- d[SNP %in% targets]
    fwrite(d, dest, sep = "\t", compress = "gzip")
    message(id, ": ", length(unique(d$SNP)), " matched SNPs; missing SNPs are not proxied.")
  }
  invisible(NULL)
}

acti_mr_summarize <- function(cfg = acti_mr_config()) {
  suppressPackageStartupMessages(library(data.table))
  d <- as.data.frame(fread(file.path(cfg$result_dir, "00_primary_IVW_Wald_results.csv")))
  ids <- c("ukb-b-15622", "finn-b-AB1_RESP_TUBERCU_CONF")
  d <- d[d$dataset %in% ids & d$direction == "forward", , drop = FALSE]
  primary <- d[d$analysis == "legacy_p1e5", , drop = FALSE]
  strict <- d[d$analysis == "strict_p5e8", , drop = FALSE]
  fwrite(primary, file.path(cfg$result_dir, "01_main_IVW_14pairs.csv"))
  fwrite(strict, file.path(cfg$result_dir, "02_sensitivity_IVW_P5e-8.csv"))
  counts <- lapply(ids, function(id) {
    z <- primary[primary$dataset == id, , drop = FALSE]
    data.frame(dataset = id, IVW_n = sum(z$method == "Inverse variance weighted"), nominal_P_lt_0_05 = sum(z$pval <
      0.05), BH7_lt_0_05 = sum(z$q_BH_7_genes < 0.05), BH14_lt_0_05 = sum(z$q_BH_14_pairs <
      0.05))
  })
  fwrite(rbindlist(counts), file.path(cfg$result_dir, "03_two_outcome_overview.csv"))
  s <- fread(file.path(cfg$result_dir, "analysis_status.csv"))
  fwrite(s, file.path(cfg$result_dir, "04_all_28_status.csv"))
  invisible(primary)
}

acti_mr_validate <- function(cfg = acti_mr_config()) {
  d <- read.csv(file.path(cfg$result_dir, "00_primary_IVW_Wald_results.csv"), check.names = FALSE)
  near <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10))
  stopifnot(all(is.finite(d$b)), all(is.finite(d$se) & d$se > 0), near(d$beta_lci95, d$b - 1.96 * d$se), near(d$beta_uci95,
    d$b + 1.96 * d$se))
  for (a in c("legacy_p1e5", "strict_p5e8")) {
    z <- d[d$analysis == a & d$direction == "forward", , drop = FALSE]
    stopifnot(nrow(z) == 14L, !anyDuplicated(z[c("gene", "dataset")]), setequal(z$gene, c("DUSP3", "FCGR1A",
      "NPC2", "TAP1", "ANKRD22", "GBP5", "VAMP5")), setequal(z$dataset, c("ukb-b-15622", "finn-b-AB1_RESP_TUBERCU_CONF")),
      near(z$q_BH_14_pairs, p.adjust(z$pval, "BH", n = 14)), all(z$method == ifelse(z$nsnp == 1L, "Wald ratio",
        "Inverse variance weighted")))
    for (id in unique(z$dataset)) {
      q <- z[z$dataset == id, , drop = FALSE]
      stopifnot(nrow(q) == 7L, near(q$q_BH_7_genes, p.adjust(q$pval, "BH", n = 7)))
    }
  }
  primary <- d[d$analysis == "legacy_p1e5" & d$direction == "forward", , drop = FALSE]
  print(primary[, c("gene", "dataset", "nsnp", "b", "pval", "q_BH_14_pairs")], row.names = FALSE)
  message("Validated 14 primary + 14 strict tests, CIs, methods and BH corrections.")
  invisible(primary)
}

acti_mr_check <- function(cfg = acti_mr_config(), stage = "all") {
  ids <- c("ukb-b-15622", "finn-b-AB1_RESP_TUBERCU_CONF")
  genes <- data.frame(gene = c("DUSP3", "FCGR1A", "NPC2", "TAP1", "ANKRD22", "GBP5", "VAMP5"), ensembl_id = c("ENSG00000108861",
    "ENSG00000150337", "ENSG00000119655", "ENSG00000168394", "ENSG00000152766", "ENSG00000154451", "ENSG00000168899"),
    opengwas_id = c("eqtl-a-ENSG00000108861", "eqtl-a-ENSG00000150337", "eqtl-a-ENSG00000119655", "eqtl-a-ENSG00000168394",
      "eqtl-a-ENSG00000152766", "eqtl-a-ENSG00000154451", "eqtl-a-ENSG00000168899"), stringsAsFactors = FALSE)
  files <- character()
  if (stage %in% c("all", "instruments"))
    files <- c(files, file.path(cfg$data_dir, "exposure", "gene_eQTL_tsv", paste0(genes$gene, "__", genes$opengwas_id,
      ".tsv.gz")), paste0(cfg$ld_bfile, c(".bed", ".bim", ".fam")))
  if (stage %in% c("all", "cache"))
    files <- c(files, file.path(cfg$cache_dir, paste0(ids, "__forward_matches.tsv.gz")))
  if (stage == "forward")
    files <- c(files, file.path(cfg$result_dir, c("gene_manifest.csv", "outcome_metadata.json", "gene_instrument_status.csv")),
      file.path(cfg$result_dir, "gene_instruments", genes$gene, "02_clumped.tsv"), file.path(cfg$result_dir,
        "cache", paste0(ids, "__forward_matches.tsv.gz")))
  if (stage %in% c("summarize", "validate"))
    files <- c(files, file.path(cfg$result_dir, "00_primary_IVW_Wald_results.csv"))
  if (stage == "forest")
    files <- c(files, file.path(cfg$result_dir, "01_main_IVW_14pairs.csv"))
  print(data.frame(path = files, exists = file.exists(files)), row.names = FALSE)
  if (!all(file.exists(files)))
    stop("Missing MR inputs.")
  if (stage %in% c("all", "instruments")) {
    found <- file.exists(cfg$plink) || nzchar(Sys.which(cfg$plink))
    cat("PLINK available:", found, "\n")
    if (!found)
      stop("Supply --plink PATH.")
  }
  invisible(TRUE)
}

acti_mr_run <- function(stage, cfg = acti_mr_config()) {
  stages <- c("initialize", "instruments", "cache", "forward", "summarize", "validate", "forest")
  if (length(stage) != 1L || !(stage %in% c("all", stages)))
    stop("Unknown stage: ", stage)
  chosen <- if (stage == "all")
    stages
  else stage
  result <- NULL
  for (s in chosen) {
    message("Running stage: ", s)
    result <- get(paste0("acti_mr_", s), mode = "function")(cfg)
  }
  invisible(result)
}

acti_mr_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- "Rscript 6_Mendelian_randomization.R --stage NAME\nStages: all, initialize, instruments, cache, forward, summarize, validate, forest\nUse --check to check inputs; see the configuration function for paths.\n"
  if (!length(args) || "--help" %in% args) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  flags <- c("--data-dir", "--output-dir", "--result-dir", "--cache-dir", "--plink", "--ld-bfile", "--stage")
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
  cfg <- do.call(acti_mr_config, values)
  if (check) {
    acti_mr_check(cfg, stage)
    return(invisible(cfg))
  }
  acti_mr_run(stage, cfg)
}

if (sys.nframe() == 0L) acti_mr_main()

