# ACTI2026

R scripts and processed machine-learning data for the study of host gene expression associated with pulmonary tuberculosis activity, PET–CT measurements and treatment stages. The analyses include single-cell RNA sequencing, paired blood co-expression networks, candidate-gene prioritization, machine learning, score validation and exploratory Mendelian randomization.

The six R scripts are placed directly in the repository root. The five archived ML runs are in `ML_data/ACTI1`–`ML_data/ACTI5`.

## Repository contents

```text
Acti2026/
├── README.md
├── 1_ACTI_single_cell_analysis.R
├── 2_WGCNA_analysis.R
├── 3_Scoring_and_3D.R
├── 4_ML_patient_level_pipeline.R
├── 5_Validation_and_plots.R
├── 6_Mendelian_randomization.R
└── ML_data/
    ├── ACTI1/
    ├── ACTI2/
    ├── ACTI3/
    ├── ACTI4/
    └── ACTI5/
```

| Script | Analysis |
|---|---|
| `1_ACTI_single_cell_analysis.R` | Cell quality control, doublet removal, Harmony integration, annotation, patient-paired pseudobulk analysis and single-cell plots. |
| `2_WGCNA_analysis.R` | Paired diagnosis–week 4 and diagnosis–week 24 networks, module selection and shared module genes. |
| `3_Scoring_and_3D.R` | Candidate-gene overlap, imaging/treatment association scores, gene ranking, heatmaps and interactive 3D plots. |
| `4_ML_patient_level_pipeline.R` | Input preparation, patient-level train/test splitting, grouped cross-validation, workflow screening and internal testing. |
| `5_Validation_and_plots.R` | Internal, external, prediagnosis and clinical analyses of fixed ACTI expression scores, with RISK6 and Sweeney3 comparisons. |
| `6_Mendelian_randomization.R` | Exploratory analyses of seven genes and two tuberculosis outcomes, sensitivity analyses and result summaries. |

## Requirements

Run the scripts in R or with `Rscript`. The single-cell workflow targets **R 4.1.3 and Seurat 4.4.0** and rejects Seurat 5. Install the packages needed for the selected analysis:

| Analysis | Main dependencies |
|---|---|
| Single-cell | Seurat, Matrix, SingleCellExperiment, scDblFinder, SingleR, harmony, edgeR, statmod; plotting packages listed in `acti_stage_packages`. |
| WGCNA | WGCNA, dplyr, pheatmap, RColorBrewer, ggplot2. |
| Scoring and 3D plots | data.table, dplyr, tidyr, ggplot2, plotly, htmlwidgets, viridisLite, scales, svglite, ragg. |
| Machine learning | randomForestSRC, glmnet, plsRglm, gbm, mboost, e1071, MASS, xgboost, pROC. |
| Validation | ggplot2, dplyr, tidyr, tibble, pROC, readxl, patchwork, ggpubr, ggtext, ComplexHeatmap, RColorBrewer, svglite. |
| Mendelian randomization | TwoSampleMR, data.table, jsonlite, dplyr, forestploter; PLINK and a European LD reference for instrument clumping. |

Each script supports `--help` and `--check`. Use package versions recorded with the original analyses where available; fitted models and selected genes can vary with package versions and numerical routines.

## Run the supplied ML analyses

Set the working directory to the repository root containing the six R scripts and `ML_data/`. The saved `data.train.txt` and `data.test.txt` files preserve the five original patient allocations and can be used directly:

```bash
# Check the saved input matrices and patient separation.
Rscript 4_ML_patient_level_pipeline.R --stage validate --input-dir ML_data/ACTI1

# Check inputs and required ML packages.
Rscript 4_ML_patient_level_pipeline.R --stage train --check --input-dir ML_data/ACTI1

# Refit workflows using the saved allocation.
Rscript 4_ML_patient_level_pipeline.R --stage train --input-dir ML_data/ACTI1 --model-dir results/ACTI1
```

Replace `ACTI1` with `ACTI2`–`ACTI5` to run another allocation. New model outputs are written to a timestamped subdirectory under the supplied `--model-dir`.

For these saved allocations, start with `validate` and `train`. The `prepare` stage generates an allocation from source expression parts and a phenotype table; it requires additional inputs configured in `acti_ml_config()`. Running `--stage all` includes that preparation step. The model-fitting seed is 1234; the preparation/split seed is a separate setting.

### ML data files

Each `ML_data/ACTI*` directory contains:

| File | Contents |
|---|---|
| `data.train.txt`, `data.test.txt` | Tab-delimited matrices with samples in rows, 20 candidate-gene expression columns, an `id` column and a binary `Type` column. |
| `patient_level_train_test_split_manifest.csv` | Sample IDs, patient IDs, treatment stages and training/test assignments; additional metadata vary between archived runs. |
| `GSE89403part1.normalize.txt`–`GSE89403part3.normalize.txt` | Gene-by-sample processed expression matrices for the three parts. |
| `merge.normalize.txt` | Merged processed expression matrix. |
| `intersect.txt` | The 20 candidate genes supplied to ML. |
| `CV_fold1_audit.rds`–`CV_fold5_audit.rds` | Saved preprocessing, feature-selection and tuning information from the five outer folds. Read with `readRDS()`. |
| `model.genes.txt` | Gene membership by workflow. |
| `five_column_AUC*.csv` | Archived workflow AUC summaries; the filename differs slightly between ACTI1 and ACTI2–ACTI5. |

`Type = 0` denotes diagnosis, day-7 or week-4 samples (**NonEndpoint**); `Type = 1` denotes week-24 samples (**Endpoint**). These labels describe treatment stage. All samples from one patient remain in the same training/test subset. Each run includes 52 training and 21 internal-test patients, with 281 samples in total.

Workflows are ranked by patient-balanced training out-of-fold AUC. Internal-test AUC is used for evaluation. New training outputs include workflow rankings, out-of-fold predictions, fitted models, preprocessing objects, selected features, performance summaries and session information. The archived data folders contain selected inputs and outputs; additional model objects are generated by rerunning the script.

## Other analyses

Public transcriptomic datasets are available through GEO: **GSE192483, GSE135893, GSE89403, GSE79362, GSE40553, GSE54992 and GSE31348**. The repository includes the processed ML inputs above. Other workflows require the corresponding raw or processed data, reference files and intermediate tables. Clinical characteristics, sampling timelines and normalized qPCR values are described in manuscript Supplementary Tables 1, 2 and 8.

Check the configuration function at the beginning of each script and supply paths for the required inputs:

| Script | Configuration function | Additional inputs |
|---|---|---|
| 1 | `acti_config()` | 10x matrices, HPCA reference and a reviewed annotation map or annotated checkpoint. |
| 2 | `acti_wgcna_config()` | GSE89403 expression and phenotype tables. |
| 3 | `acti_fig4_config()` | Single-cell results, WGCNA outputs, blood expression, phenotype and imaging tables. |
| 5 | `acti_validation_config()` | Cohort-specific expression/phenotype tables and clinical qPCR inputs. |
| 6 | `acti_mr_config()` | eQTL and outcome summary data, PLINK executable, LD reference and outcome cache files. |

Use `Rscript <script_name> --help` to list the available stages. Preserve the filename capitalization, including `3_Scoring_and_3D.R`. Script 1 uses `--option=value`; scripts 2–6 use `--option value`. Calling `source()` loads the functions without starting an analysis, so paths can also be set through the configuration functions in an R session. Default input paths are relative to `data/`, and default outputs are under `results/`.

The external and clinical analyses evaluate fixed ACTI expression scores rather than refitting the original ML classifiers. Upstream gene prioritization uses the development cohort, so patient separation during ML fitting does not make the entire discovery pipeline independent. Refer to the manuscript for the study design and interpretation limits.
