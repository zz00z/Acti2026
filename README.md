# PTB-treatment-response-signature
This repository contains representative R scripts used for the analysis of host gene signatures associated with pulmonary tuberculosis (PTB) disease activity, treatment response, PET-CT-linked disease burden, and Mendelian randomization evidence.

The scripts correspond to the main analytical workflow of the manuscript, including single-cell RNA-seq analysis, WGCNA, imaging-linked correlation analysis, machine-learning model screening, external and clinical validation, and Mendelian randomization analysis.

No raw clinical data are included in this repository.

---

## Repository structure

```text
PTB-host-response-signature/
├── README.md
├── .gitignore
└── scripts/
    ├── 01_scRNA_HD_LOW_HIGH_Harmony_Annotation_ImmuneDEG.R
    ├── 02_WGCNA_GSE89403_DX_W4_W24_TreatmentResponse.R
    ├── 03_ImagingLinkedHostResponse_Correlation_3DPlot.R
    ├── 04A_Prepare_GSE89403_ML_Input_TreatmentStage.R
    ├── 04B_ML_TreatmentStage_ModelScreening.R
    ├── 05_GeneSignature_Validation_Plots.R
    └── 06_MR_CASP4_VAMP5_Tuberculosis.R
