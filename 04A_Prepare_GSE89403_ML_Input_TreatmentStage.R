############################################################
## 04A_Prepare_GSE89403_ML_Input_TreatmentStage.R
##
## Prepare expression matrices for machine-learning analysis
## of treatment-stage assessment in GSE89403.
##
## Classification rule:
##   DX / day_7 / week_4  = Control
##   week_24              = Treat
##
## Output:
##   merge.normalize.txt
##   data.train.txt
##   data.test.txt
############################################################


## =========================================================
## 0. Load packages
## =========================================================

library(limma)
library(sva)


## =========================================================
## 1. Set English directory structure
## =========================================================

base_dir <- "D:/atraining/blood/analysis"

raw_dir  <- file.path(base_dir, "data_merge")
prep_dir <- file.path(base_dir, "prepare")
ml_dir   <- file.path(base_dir, "ml_model")

dir.create(base_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(raw_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(prep_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(ml_dir,    showWarnings = FALSE, recursive = TRUE)

set.seed(123)


## =========================================================
## 2. Read three normalized GSE89403 expression matrices
## =========================================================
## Required files in raw_dir:
##   GSE89403part1.normalize.txt
##   GSE89403part2.normalize.txt
##   GSE89403part3.normalize.txt
##
## The first column should contain gene IDs or gene symbols.

setwd(raw_dir)

part1 <- read.table(
  "GSE89403part1.normalize.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

part2 <- read.table(
  "GSE89403part2.normalize.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

part3 <- read.table(
  "GSE89403part3.normalize.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

## Ensure the first column is named id
if (!"id" %in% colnames(part1)) colnames(part1)[1] <- "id"
if (!"id" %in% colnames(part2)) colnames(part2)[1] <- "id"
if (!"id" %in% colnames(part3)) colnames(part3)[1] <- "id"

## Check that gene order is identical across the three files
stopifnot(all(part1$id == part2$id))
stopifnot(all(part1$id == part3$id))


## =========================================================
## 3. Merge all samples
## =========================================================

expr_all <- cbind(
  part1,
  part2[, -1, drop = FALSE],
  part3[, -1, drop = FALSE]
)

cat("Merged expression matrix:\n")
cat(nrow(expr_all), "genes ×", ncol(expr_all) - 1, "samples\n")


## =========================================================
## 4. Randomly split samples into three validation subsets
## =========================================================
## This follows the original analysis strategy:
##   1. Merge all samples
##   2. Randomly shuffle sample columns
##   3. Split samples into three subsets
##   4. Add cohort prefixes for downstream model evaluation

sample_cols <- sample(colnames(expr_all)[-1])

total_cols <- length(sample_cols)
part_size  <- ceiling(total_cols / 3)

cols_A <- sample_cols[1:part_size]
cols_B <- sample_cols[(part_size + 1):min(2 * part_size, total_cols)]
cols_C <- sample_cols[(min(2 * part_size, total_cols) + 1):total_cols]
cols_C <- cols_C[!is.na(cols_C)]

partA <- cbind(expr_all["id"], expr_all[, cols_A, drop = FALSE])
partB <- cbind(expr_all["id"], expr_all[, cols_B, drop = FALSE])
partC <- cbind(expr_all["id"], expr_all[, cols_C, drop = FALSE])

setwd(prep_dir)

write.table(
  partA,
  "GSE89403part1.normalize.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  partB,
  "GSE89403part2.normalize.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  partC,
  "GSE89403part3.normalize.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Three randomly split normalized matrices were saved to prepare directory.\n")


## =========================================================
## 5. Correct Control/Treat suffix
## =========================================================
## Rule:
##   week_24 / week24 = Treat
##   DX / day_7 / week_4 = Control

fix_treatment_stage_label <- function(x) {
  
  base <- sub("_(Control|Treat)$", "", x, ignore.case = TRUE)
  
  is_week24 <- grepl("week_?24|week24", base, ignore.case = TRUE)
  
  suffix <- ifelse(is_week24, "Treat", "Control")
  
  paste0(base, "_", suffix)
}

check_invalid_label <- function(nms) {
  
  any(grepl("week_?24.*_Control$", nms, ignore.case = TRUE)) ||
    any(grepl("(DX|day_?7|week_?4).*_Treat$", nms, ignore.case = TRUE))
}


p1 <- read.table(
  "GSE89403part1.normalize.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

p2 <- read.table(
  "GSE89403part2.normalize.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

p3 <- read.table(
  "GSE89403part3.normalize.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

if (!"id" %in% colnames(p1)) colnames(p1)[1] <- "id"
if (!"id" %in% colnames(p2)) colnames(p2)[1] <- "id"
if (!"id" %in% colnames(p3)) colnames(p3)[1] <- "id"

stopifnot(all(p1$id == p2$id))
stopifnot(all(p1$id == p3$id))

colnames(p1)[-1] <- fix_treatment_stage_label(colnames(p1)[-1])
colnames(p2)[-1] <- fix_treatment_stage_label(colnames(p2)[-1])
colnames(p3)[-1] <- fix_treatment_stage_label(colnames(p3)[-1])

if (
  check_invalid_label(colnames(p1)[-1]) ||
  check_invalid_label(colnames(p2)[-1]) ||
  check_invalid_label(colnames(p3)[-1])
) {
  stop("Invalid Control/Treat labels detected. Please check sample names.")
}


## =========================================================
## 6. Add cohort prefixes
## =========================================================

colnames(p1)[-1] <- paste0("GSE89403part1_", colnames(p1)[-1])
colnames(p2)[-1] <- paste0("GSE89403part2_", colnames(p2)[-1])
colnames(p3)[-1] <- paste0("GSE89403part3_", colnames(p3)[-1])


## =========================================================
## 7. Generate merged matrix for training
## =========================================================
## This file is used as the training cohort after adding
## the prefix "merge_" in the next step.

expr_final <- cbind(
  p1,
  p2[, -1, drop = FALSE],
  p3[, -1, drop = FALSE]
)

colnames(expr_final)[1] <- "geneNames"

write.table(
  expr_final,
  "merge.normalize.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("merge.normalize.txt generated.\n")
cat("Final merged matrix:\n")
cat(nrow(expr_final), "genes ×", ncol(expr_final) - 1, "samples\n")


## =========================================================
## 8. Prepare machine-learning input files
## =========================================================
## Required file in prep_dir:
##   intersect.txt
##
## intersect.txt should contain candidate genes used for
## machine-learning screening.

gene_file <- file.path(prep_dir, "intersect.txt")

if (!file.exists(gene_file)) {
  stop("intersect.txt was not found in prepare directory.")
}

all_files <- dir(
  prep_dir,
  pattern = "normalize.txt$",
  full.names = TRUE
)

cat("Files used for ML input preparation:\n")
print(basename(all_files))


## =========================================================
## 9. Identify common genes across all normalized matrices
## =========================================================

gene_list <- list()

for (file in all_files) {
  
  rt <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  
  gene_names <- as.vector(rt[, 1])
  uniq_gene <- unique(gene_names)
  
  header <- unlist(strsplit(basename(file), "\\.|\\-"))
  gene_list[[header[1]]] <- uniq_gene
}

inter_genes <- Reduce(intersect, gene_list)

cat("Number of common genes across normalized matrices:", length(inter_genes), "\n")


## =========================================================
## 10. Merge expression matrices for ML input
## =========================================================

all_tab <- data.frame()
batch_type <- c()

for (i in seq_along(all_files)) {
  
  input_file <- all_files[i]
  header <- unlist(strsplit(basename(input_file), "\\.|\\-"))
  
  rt <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  
  dimnames <- list(rownames(exp), colnames(exp))
  
  data <- matrix(
    as.numeric(as.matrix(exp)),
    nrow = nrow(exp),
    dimnames = dimnames
  )
  
  ## Collapse duplicated genes by average expression
  rt <- avereps(data)
  
  ## Add dataset prefix
  colnames(rt) <- paste0(header[1], "_", colnames(rt))
  
  if (i == 1) {
    all_tab <- rt[inter_genes, ]
  } else {
    all_tab <- cbind(all_tab, rt[inter_genes, ])
  }
  
  batch_type <- c(batch_type, rep(i, ncol(rt)))
}


## =========================================================
## 11. Optional batch correction
## =========================================================
## Final version:
##   No batch correction was applied.
##
## If needed:
##   sva_tab <- ComBat(all_tab, batch_type, par.prior = TRUE)

sva_tab <- all_tab


## =========================================================
## 12. Keep candidate genes
## =========================================================

gene_rt <- read.table(
  gene_file,
  header = FALSE,
  sep = "\t",
  check.names = FALSE
)

candidate_genes <- as.vector(gene_rt[, 1])

gene_tab <- sva_tab[
  intersect(rownames(sva_tab), candidate_genes),
  ,
  drop = FALSE
]

gene_tab <- t(gene_tab)

cat("Number of candidate genes used for ML:", ncol(gene_tab), "\n")


## =========================================================
## 13. Split training and validation data
## =========================================================
## Original rule:
##   Samples with names beginning with "merge" are used
##   as the training cohort.
##   Other samples are used as validation cohorts.

train_index <- grepl("^merge", rownames(gene_tab), ignore.case = TRUE)

train_exp <- gene_tab[train_index, , drop = FALSE]
test_exp  <- gene_tab[!train_index, , drop = FALSE]

rownames(train_exp) <- gsub("merge_", "Train.", rownames(train_exp))

train_type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(train_exp))
test_type  <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(test_exp))

train_type <- ifelse(train_type == "Control", 0, 1)
test_type  <- ifelse(test_type == "Control", 0, 1)

train_exp <- cbind(train_exp, Type = train_type)
test_exp  <- cbind(test_exp, Type = test_type)

cat("Training samples:", nrow(train_exp), "\n")
cat("Validation samples:", nrow(test_exp), "\n")
cat("Training class distribution:\n")
print(table(train_exp[, "Type"]))
cat("Validation class distribution:\n")
print(table(test_exp[, "Type"]))


## =========================================================
## 14. Export machine-learning input files
## =========================================================

setwd(ml_dir)

train_out <- rbind(id = colnames(train_exp), train_exp)

write.table(
  train_out,
  file = "data.train.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

test_out <- rbind(id = colnames(test_exp), test_exp)

write.table(
  test_out,
  file = "data.test.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

cat("data.train.txt and data.test.txt generated in:\n")
cat(ml_dir, "\n")
cat("ML input preparation completed.\n")