############################################################
## 03_ImagingLinkedHostResponse_Correlation_3DPlot.R
##
## Correlation analysis and three-dimensional visualization
## of imaging-linked host response genes.
##
## Main purpose:
## 1. Use 464 shared genes from scRNA-seq DEGs and WGCNA modules
## 2. Calculate Pearson correlations with TGRV, SUVmax, and treatment time
## 3. Define a composite correlation score:
##    score = -log10(p_TGRV) - log10(p_SUVmax) - log10(p_Time)
## 4. Generate 3D visualization for all 464 genes and top-ranked genes
## 5. Export Supplementary Table 6
############################################################


## =========================================================
## 0. Load packages
## =========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(plotly)
library(viridis)
library(htmlwidgets)


## =========================================================
## 1. Set working directory
## =========================================================

setwd("D:/atraining/blood/GSE89403/PETCTSUVCA")


## =========================================================
## 2. Read input data
## =========================================================
## exp2.csv:
##   Gene expression matrix.
##   Rows = genes, columns = samples.
##
## pet_data.csv:
##   PET-CT TGRV data.
##   Required columns:
##     sample_code
##     TGRV_log10
##
## SUVmax_data.csv:
##   PET-CT SUVmax data.
##   Required columns:
##     sample_code
##     SUVmax
##
## 464.csv:
##   Final 464 shared genes from scRNA-seq DEGs and WGCNA modules.
##   Required column:
##     Gene

expr <- read.csv(
  "exp2.csv",
  row.names = 1,
  check.names = FALSE
)

pet <- read.csv(
  "pet_data.csv",
  check.names = FALSE
)

suv <- read.csv(
  "SUVmax_data.csv",
  check.names = FALSE
)

genes_464 <- read.csv(
  "464.csv",
  check.names = FALSE
)

genes_464$Gene <- as.character(genes_464$Gene)

## Keep only genes present in expression matrix
genes_use <- intersect(genes_464$Gene, rownames(expr))

cat("Number of input shared genes:", length(unique(genes_464$Gene)), "\n")
cat("Number of shared genes found in expression matrix:", length(genes_use), "\n")

expr_464 <- expr[genes_use, , drop = FALSE]


## =========================================================
## 3. Define helper function for gene-wise correlation
## =========================================================

safe_cor_test <- function(x, y, method = "pearson") {
  
  df <- data.frame(
    x = as.numeric(x),
    y = as.numeric(y)
  )
  
  df <- df[complete.cases(df), , drop = FALSE]
  
  if (nrow(df) < 3 || length(unique(df$x)) <= 1 || length(unique(df$y)) <= 1) {
    return(data.frame(
      estimate = NA_real_,
      p.value = NA_real_,
      n = nrow(df)
    ))
  }
  
  ct <- suppressWarnings(
    cor.test(df$x, df$y, method = method)
  )
  
  data.frame(
    estimate = as.numeric(ct$estimate),
    p.value = ct$p.value,
    n = nrow(df)
  )
}


gene_cor <- function(expr_sub, pheno, value_col, method = "pearson") {
  
  res <- lapply(rownames(expr_sub), function(g) {
    
    x <- as.numeric(expr_sub[g, ])
    y <- pheno[[value_col]]
    
    tmp <- safe_cor_test(x, y, method = method)
    
    data.frame(
      Gene = g,
      rho = tmp$estimate,
      p = tmp$p.value,
      n = tmp$n,
      neglog10p = -log10(tmp$p.value),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(res)
}


## =========================================================
## 4. Correlation with PET-CT TGRV
## =========================================================

samples_tgrv <- intersect(colnames(expr_464), pet$sample_code)

expr_tgrv <- expr_464[, samples_tgrv, drop = FALSE]

pet_tgrv <- pet %>%
  filter(sample_code %in% samples_tgrv) %>%
  arrange(match(sample_code, samples_tgrv))

## Ensure sample order is exactly matched
expr_tgrv <- expr_tgrv[, pet_tgrv$sample_code, drop = FALSE]

res_tgrv <- gene_cor(
  expr_sub = expr_tgrv,
  pheno = pet_tgrv,
  value_col = "TGRV_log10",
  method = "pearson"
) %>%
  rename(
    rho_TGRV = rho,
    p_TGRV = p,
    n_TGRV = n,
    neglog10p_TGRV = neglog10p
  )

write.csv(
  res_tgrv,
  "gene_TGRV_correlation_464.csv",
  row.names = FALSE
)


## =========================================================
## 5. Correlation with PET-CT SUVmax
## =========================================================

samples_suv <- intersect(colnames(expr_464), suv$sample_code)

expr_suv <- expr_464[, samples_suv, drop = FALSE]

suv_data <- suv %>%
  filter(sample_code %in% samples_suv) %>%
  arrange(match(sample_code, samples_suv))

## Ensure sample order is exactly matched
expr_suv <- expr_suv[, suv_data$sample_code, drop = FALSE]

res_suv <- gene_cor(
  expr_sub = expr_suv,
  pheno = suv_data,
  value_col = "SUVmax",
  method = "pearson"
) %>%
  rename(
    rho_SUVmax = rho,
    p_SUVmax = p,
    n_SUVmax = n,
    neglog10p_SUVmax = neglog10p
  )

write.csv(
  res_suv,
  "gene_SUVmax_correlation_464.csv",
  row.names = FALSE
)


## =========================================================
## 6. Correlation with treatment time
## =========================================================
## Treatment time was encoded as a continuous variable:
## DX      = 0 weeks
## week_4  = 4 weeks
## week_24 = 24 weeks
##
## day_7 samples were excluded according to the final analysis.

expr_time <- expr_464[, !str_detect(colnames(expr_464), "day_7"), drop = FALSE]

timepoints <- str_extract(
  colnames(expr_time),
  regex("DX|week_\\d+", ignore_case = TRUE)
)

time_continuous <- case_when(
  str_to_upper(timepoints) == "DX" ~ 0,
  timepoints == "week_4" ~ 4,
  timepoints == "week_24" ~ 24,
  TRUE ~ NA_real_
)

cat("Time point distribution:\n")
print(table(timepoints, useNA = "ifany"))

cat("Continuous time coding:\n")
print(table(time_continuous, useNA = "ifany"))

res_time <- lapply(rownames(expr_time), function(g) {
  
  x <- as.numeric(expr_time[g, ])
  y <- time_continuous
  
  tmp <- safe_cor_test(x, y, method = "pearson")
  
  data.frame(
    Gene = g,
    rho_Time = tmp$estimate,
    p_Time = tmp$p.value,
    n_Time = tmp$n,
    neglog10p_Time = -log10(tmp$p.value),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

write.csv(
  res_time,
  "gene_treatment_time_correlation_464.csv",
  row.names = FALSE
)


## =========================================================
## 7. Merge TGRV, SUVmax, and treatment time correlations
## =========================================================

cor_464 <- res_tgrv %>%
  full_join(res_suv, by = "Gene") %>%
  full_join(res_time, by = "Gene") %>%
  mutate(
    neglog10p_TGRV = ifelse(is.infinite(neglog10p_TGRV), NA, neglog10p_TGRV),
    neglog10p_SUVmax = ifelse(is.infinite(neglog10p_SUVmax), NA, neglog10p_SUVmax),
    neglog10p_Time = ifelse(is.infinite(neglog10p_Time), NA, neglog10p_Time),
    
    score = neglog10p_TGRV + neglog10p_SUVmax + neglog10p_Time,
    
    imaging_score = neglog10p_TGRV + neglog10p_SUVmax
  ) %>%
  arrange(desc(score))

write.csv(
  cor_464,
  "Supplementary_Table6_464_genes_correlation_score.csv",
  row.names = FALSE
)

write.csv(
  cor_464,
  "464_deg_correlation.csv",
  row.names = FALSE
)

cat("Number of genes with complete correlation score:\n")
print(sum(!is.na(cor_464$score)))


## =========================================================
## 8. Export top-ranked genes by composite score
## =========================================================
## Choose top_n according to the final figure.
## If Fig. 4F highlights 20 genes, use top_n = 20.
## If Fig. 4F highlights 50 genes, use top_n = 50.

top_n <- 20

top_genes <- cor_464 %>%
  filter(!is.na(score)) %>%
  arrange(desc(score)) %>%
  slice_head(n = top_n)

write.csv(
  top_genes,
  paste0("Top", top_n, "_imaging_linked_host_response_genes.csv"),
  row.names = FALSE
)


## =========================================================
## 9. 3D visualization for all 464 genes
## =========================================================
## Fig. 4B:
## x-axis: correlation with TGRV
## y-axis: correlation with SUVmax
## z-axis: correlation with treatment time
## color: significance of treatment-time correlation
## size: combined imaging-associated significance

plot_all <- cor_464 %>%
  filter(
    !is.na(rho_TGRV),
    !is.na(rho_SUVmax),
    !is.na(rho_Time),
    !is.na(neglog10p_Time),
    !is.na(imaging_score)
  )

fig_all <- plot_ly(
  plot_all,
  x = ~rho_TGRV,
  y = ~rho_SUVmax,
  z = ~rho_Time,
  text = ~paste0(
    "Gene: ", Gene,
    "<br>r_TGRV: ", round(rho_TGRV, 3),
    "<br>p_TGRV: ", signif(p_TGRV, 3),
    "<br>r_SUVmax: ", round(rho_SUVmax, 3),
    "<br>p_SUVmax: ", signif(p_SUVmax, 3),
    "<br>r_Time: ", round(rho_Time, 3),
    "<br>p_Time: ", signif(p_Time, 3),
    "<br>Score: ", round(score, 3)
  ),
  color = ~neglog10p_Time,
  colors = viridis(10),
  size = ~imaging_score,
  marker = list(opacity = 0.8)
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = "Correlation with TGRV"),
      yaxis = list(title = "Correlation with SUVmax"),
      zaxis = list(title = "Correlation with treatment time")
    )
  )

fig_all

htmlwidgets::saveWidget(
  fig_all,
  "Fig4B_3D_correlation_all_464_genes.html",
  selfcontained = TRUE
)


## =========================================================
## 10. 3D visualization for top-ranked genes
## =========================================================
## Fig. 4F:
## Top-ranked genes by composite correlation score.

fig_top <- plot_ly(
  top_genes,
  x = ~rho_TGRV,
  y = ~rho_SUVmax,
  z = ~rho_Time,
  text = ~paste0(
    "Gene: ", Gene,
    "<br>r_TGRV: ", round(rho_TGRV, 3),
    "<br>p_TGRV: ", signif(p_TGRV, 3),
    "<br>r_SUVmax: ", round(rho_SUVmax, 3),
    "<br>p_SUVmax: ", signif(p_SUVmax, 3),
    "<br>r_Time: ", round(rho_Time, 3),
    "<br>p_Time: ", signif(p_Time, 3),
    "<br>Score: ", round(score, 3)
  ),
  color = ~score,
  colors = viridis(10),
  size = ~score,
  marker = list(opacity = 0.9)
) %>%
  add_markers() %>%
  add_text(
    text = ~Gene,
    textposition = "top center",
    showlegend = FALSE,
    textfont = list(size = 8)
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "Correlation with TGRV"),
      yaxis = list(title = "Correlation with SUVmax"),
      zaxis = list(title = "Correlation with treatment time")
    )
  )

fig_top

htmlwidgets::saveWidget(
  fig_top,
  paste0("Fig4F_3D_correlation_top", top_n, "_genes.html"),
  selfcontained = TRUE
)


## =========================================================
## 11. Save R objects
## =========================================================

saveRDS(
  list(
    res_tgrv = res_tgrv,
    res_suv = res_suv,
    res_time = res_time,
    cor_464 = cor_464,
    top_genes = top_genes,
    fig_all = fig_all,
    fig_top = fig_top
  ),
  "Imaging_linked_host_response_correlation_results.rds"
)

cat("Correlation analysis and 3D visualization completed.\n")