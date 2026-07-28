<div align="center">

# BreastCancer-LuminalB-Bioinformatics

### Bioinformatics workflow for inflammatory cytokine profiling in Luminal B breast cancer

[![R](https://img.shields.io/badge/R-4.5+-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-1B9E77)](https://bioconductor.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20Windows-blue)]()
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)]()
[![License](https://img.shields.io/badge/License-MIT-green)]()

</div>

---

# Overview

This repository documents the complete bioinformatics workflow developed during an individual mentoring program focused on **bulk RNA-seq** and **single-cell RNA-seq** analyses in breast cancer.

The project is based on a Master's research investigating the inflammatory cytokine profile in **Luminal B Breast Cancer**, integrating public transcriptomic datasets from **TCGA**, **METABRIC**, and **single-cell RNA sequencing** studies.

The repository is continuously updated throughout the mentoring sessions, emphasizing reproducible research, good programming practices, and biological interpretation of the results.

---

# Scientific Objectives

- Download and process TCGA-BRCA RNA-seq data
- Retrieve and curate clinical metadata
- Filter molecular subtypes (PAM50)
- Perform RNA-seq quality assessment
- Explore inflammatory cytokine expression
- Differential expression analysis
- Survival analysis (Kaplan-Meier & Cox Regression)
- Validation using METABRIC cohort
- Immune deconvolution
- Single-cell RNA-seq analysis
- Cell-cell communication using CellChat

---

# Repository Structure

```text
BreastCancer-LuminalB-Bioinformatics/

├── 01_TCGA/
│
├── 02_Preprocessing/
│
├── 03_Exploratory_Analysis/
│
├── 04_Cytokine_Analysis/
│
├── 05_Survival/
│
├── 06_METABRIC/
│
├── 07_SingleCell_Ozmen2025/
│
├── 08_CellChat/
│
├── Figures/
│
├── Results/
│
├── Articles/
│
└── README.md
```

---

# Public Datasets

| Dataset | Description |
|----------|-------------|
| TCGA-BRCA | RNA-seq and clinical data |
| METABRIC | Validation cohort |
| Ozmen et al. (2025) | Single-cell RNA sequencing |

---

# Main R Packages

### Data acquisition

- TCGAbiolinks
- cBioPortalData

### Data manipulation

- dplyr
- tidyr
- tibble
- janitor

### RNA-seq analysis

- DESeq2
- edgeR
- limma
- SummarizedExperiment

### Visualization

- ggplot2
- ComplexHeatmap
- EnhancedVolcano
- pheatmap

### Functional analysis

- clusterProfiler
- enrichplot
- AnnotationDbi
- org.Hs.eg.db
- biomaRt

### Survival analysis

- survival
- survminer
- timeROC

### Single-cell RNA-seq

- Seurat
- CellChat

---

# Learning Philosophy

Each analysis is developed from both biological and computational perspectives.

Every module includes:

- Biological background
- Statistical concepts
- Reproducible R scripts
- Result interpretation
- Git version control
- Scientific discussion

---

# Planned Workflow

- [ ] TCGA download
- [ ] Clinical metadata
- [ ] PAM50 subtype selection
- [ ] Quality control
- [ ] Exploratory analysis
- [ ] Cytokine panel
- [ ] Differential expression
- [ ] Survival analysis
- [ ] METABRIC validation
- [ ] Immune deconvolution
- [ ] Single-cell RNA-seq
- [ ] CellChat
- [ ] Final figures

---

# Citation

If you use this repository or parts of its workflow, please cite the corresponding public datasets and software packages.

---

# Author

**Jean Resende**

PhD Candidate in Biotechnology

Bioinformatics • Cancer Genomics • Immunogenomics

---

## ⭐ Repository Status

This repository is under continuous development as part of a bioinformatics mentoring program and will be updated as new analyses are completed.
