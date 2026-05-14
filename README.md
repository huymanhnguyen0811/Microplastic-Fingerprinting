# PlastiPrint — Microplastic Computational Fingerprinting Tool

PlastiPrint is an open-source R Shiny application for chemical fingerprinting of microplastics using machine learning. It wraps the complete computational workflow from the publication **"Chemical Fingerprints of New vs. Weathered Microplastics: A Machine Learning Approach"** into an interactive GUI with configurable parameters at every step.

## Overview

The app provides a 9-tab guided workflow that takes users from raw instrumental data (ATD-GC-MS, HPLC-QToF-MS, ICP-MS/MS) through data processing, feature selection, ML classification (Random Forest, SVM, XGBoost), and hierarchical cluster analysis (HCA) — all without writing code.

## License

- **Code**: GNU AGPL v3.0 only (AGPL-3.0-only) — see `LICENSE.txt`
- **Data**: CC BY-SA 4.0 — see `Raw Data - for Github only/LICENSE-DATA.txt`

---

## 💻 System Specifications

Developed and tested on:
- **OS**: Windows 11 Pro 64-bit (Build 26200)
- **Processor**: 11th Gen Intel Core i7-11800H @ 2.30 GHz (16 CPUs)
- **Memory**: 16 GB RAM

---

## Prerequisites

- **R ≥ 4.2** — <https://cran.r-project.org>
- **RStudio** (recommended) — <https://posit.co/downloads/>

---

## 📂 Repository Structure

| File / Folder | Description |
|---|---|
| `PlastiPrint_app-UI.R` | The Shiny application (main entry point) |
| `Helper Functions for Computational Workflow_06Feb2026.R` | Custom functions for compound grouping, imputation, and feature filtering |
| `Complete workflow.Rmd` | The original R Markdown workflow (non-GUI version) |
| `recommended worflow for imputation and normalization and ML algorithms selection_05Feb2026.R` | Diagnostic tests for choosing appropriate ML algorithms |
| `Scripts for Tables and Figures (maintext and SI).Rmd` | Scripts for publication figures and tables |
| `Raw Data - for Github only/` | Raw instrumental data (ATD-GC-MS, HPLC-QToF-MS, ICP-MS/MS) |
| `Supplementary documents for publication/` | Supplementary materials accompanying the publication |
| `R package versions.txt` | Exact package versions for reproducibility |

---

## Installation

### 1. Place files

1. Clone or download this repository.
2. Ensure `PlastiPrint_app-UI.R` and the helper `.R` file are in the project root (same level as `Raw Data - for Github only/`).

### 2. Install R packages

The first launch automatically installs any missing CRAN packages (5–10 min on initial run). To pre-install manually:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "DT",
  "tidyverse", "readxl", "writexl", "data.table",
  "crayon", "grid", "gridExtra", "ggplot2", "pheatmap",
  "viridis", "ggpubr", "scales",
  "stats", "vegan", "cluster", "factoextra", "FactoMineR",
  "ggforce", "MASS", "moments",
  "caret", "randomForest", "ranger", "e1071", "xgboost",
  "missForest", "pROC", "VIM", "softImpute", "RANN",
  "foreach", "doParallel", "parallel"
))
```

<details>
<summary><strong>Full package list by category</strong></summary>

| Category | Packages |
|---|---|
| **GUI** | shiny, shinydashboard, shinyjs, shinyWidgets, DT |
| **Data** | tidyverse, readxl, writexl, data.table |
| **Visualization** | ggplot2, pheatmap, viridis, ggpubr, scales, gridExtra |
| **Statistics** | vegan, cluster, factoextra, FactoMineR, MASS, moments |
| **ML** | caret, randomForest, ranger, e1071, xgboost |
| **Missing data** | missForest, pROC, VIM, softImpute, RANN |
| **Parallel** | foreach, doParallel, parallel |

For exact versions used during development, see `R package versions.txt`.

</details>

---

## Launching the App

**RStudio:** Open `PlastiPrint_app-UI.R` → click **Run App**.

**R console:**
```r
setwd("/path/to/PlastiPrint")
shiny::runApp("PlastiPrint_app-UI.R")
```

**Terminal:**
```bash
cd /path/to/PlastiPrint
Rscript -e 'shiny::runApp("PlastiPrint_app-UI.R", launch.browser = TRUE)'
```

---

## Using PlastiPrint — Step by Step

Work through the sidebar tabs sequentially:

### Welcome Tab
Set your **Project Root Directory** (the folder containing `Raw Data/`) and click **Set Directory**.

### Step 1 — Import & Grouping
Run three sub-steps in order:
1. **Import Data** — reads raw CSV/XLS files from `Raw Data/ATDGCMS` and `Raw Data/HPLCTOFMS`.
2. **Group Compounds** — bins compounds by RT and m/z (adjustable thresholds).
3. **Remove Benchmarks** — strips internal-standard compounds.

### Step 2 — Blank Subtraction
Remove signals below `mean_blank + N × SD_blank` (default SD multiplier: 3×).

### Step 3 — Source Assignment
Automatically labels samples as **Store-Bought** or **Environmental** based on filename patterns.

### Step 4 — Trace Metal Import
Reads ICP-MS Excel files. Configurable: which seawater metals to discard (default: Na, Mg, K, Ca).

### Step 5 — Label Merge
Merges plastic type, subcategory, and polymer labels from the sample-labelling file. Creates all dataset combinations (gc, hplc, icp, and their fusions).

### Step 6 — Feature Reduction & Data Fusion
- **6.1–6.3 Feature Reduction**: choose from regular 80% rule, modified 80% rule, or in-house 4-rule filter.
- **6.4 Data Fusion**: applies in-house filtering per instrument, then fuses multi-instrument datasets.

### Step 7 — Hybrid ML Pipeline
Full ML optimization with configurable parameters:

| Parameter | Default | Description |
|---|---|---|
| Training approach | SB → ENV | Option 1 or Option 2 |
| Datasets | gc, gc_icp | Which cleaned datasets to analyze |
| Imputation | half_min, median, knn | Methods to search over |
| Normalization | none, log, log10, zscore, pareto | Methods to search over |
| Algorithms | ranger, svmRadial, xgbTree | ML classifiers |
| Top % features | 0.20 | Features considered "top" per pipeline |
| Stability threshold | 0.60 | Min stability for robust features |
| CV folds | 5 | Cross-validation folds |
| Top-k pipelines | 5 | Ensemble size |
| Permutation reps | 10 | Permutation importance iterations |

Results appear in sub-tabs: Best Pipeline, Metrics (MCC, accuracy, F1, kappa), Confusion Matrix, Feature Tiers, and Top Features.

### Step 8 — Hierarchical Cluster Analysis
Select Step 7 results, configure distance metric (manhattan/euclidean/bray) and linkage method, then export dendrograms as PDF.

---

## Troubleshooting

| Issue | Solution |
|---|---|
| "Directory not found" | Check the project root path on the Welcome tab |
| "Helper functions file not found" | Ensure the helper `.R` file is in the same directory as the app |
| Step N fails with "Run Step N-1 first" | Steps must be run sequentially |
| Step 7 is very slow | Reduce imputation/normalization/algorithm combinations, or lower `top_k` / `n_permutations` |
| Out-of-memory error | Close other applications or increase R's memory limit on Windows |
| Package install failure | Run `install.packages("package_name")` manually to see the full error |
