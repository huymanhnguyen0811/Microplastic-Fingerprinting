# PlastiPrint — Installation & Usage Guide

## Overview

PlastiPrint is a Shiny-based GUI that wraps the complete computational workflow from *"Chemical Fingerprints of New vs. Weathered Microplastics: A Machine Learning Approach."* Each step of the workflow is exposed as an interactive tab with configurable parameters.

---

## Prerequisites

- **R ≥ 4.2** (download from https://cran.r-project.org)
- **RStudio** (recommended, download from https://posit.co/downloads/)

---

## Installation

### 1. Create the project folder structure

```
PlastiPrint/
├── app.R                                                  ← the GUI file
├── Helper_Functions_for_Computational_Workflow_06Feb2026.R ← your helper file
└── Raw data/                                              ← your existing data folder
    ├── ATDGCMS/
    ├── HPLCTOFMS/
    │   ├── EF_Non-target data_Batch 0/
    │   ├── EF_Non-target data_Batch 1/
    │   └── EF_Non-target data_Batch 2/
    ├── ICPMS_Trace metal/
    └── Sample Labelling_all data_GC+HPLC+ICP_30Sept2025.xlsx
```

### 2. Place files

1. Copy **`app.R`** into the project root folder (same level as `Raw data/`).
2. Copy **`Helper_Functions_for_Computational_Workflow_06Feb2026.R`** into the same folder.
3. Ensure your **`Raw data/`** directory (with all sub-folders and files) is present.

### 3. Install R packages (automatic)

The first time you launch PlastiPrint, it will automatically install any missing CRAN packages. This may take 5–10 minutes on the initial run. The packages include:

| Category | Packages |
|---|---|
| **GUI** | shiny, shinydashboard, shinyjs, shinyWidgets, DT |
| **Data** | tidyverse, readxl, writexl, data.table |
| **Visualization** | ggplot2, pheatmap, viridis, ggpubr, scales, gridExtra |
| **Statistics** | vegan, cluster, factoextra, FactoMineR, MASS, moments |
| **ML** | caret, randomForest, ranger, e1071, xgboost |
| **Missing data** | missForest, pROC, VIM, softImpute, RANN |
| **Parallel** | foreach, doParallel, parallel |

If you prefer to pre-install packages manually, open R and run:

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

---

## Launching the App

### Option A: From RStudio

1. Open `app.R` in RStudio.
2. Click the **"Run App"** button (green triangle at the top right of the source pane).
3. The app will open in your default browser or in RStudio's Viewer.

### Option B: From the R console

```r
setwd("/path/to/PlastiPrint")  # adjust this to your folder
shiny::runApp("app.R")
```

### Option C: From the terminal

```bash
cd /path/to/PlastiPrint
Rscript -e 'shiny::runApp("app.R", launch.browser = TRUE)'
```

---

## Using PlastiPrint — Step-by-Step

When the app launches you will see a sidebar with 9 tabs. Work through them sequentially:

### Welcome Tab
- Set your **Project Root Directory** to the folder containing `Raw data/`.
- Click **"Set Directory"** to confirm.

### Step 1 — Import & Grouping
Run the three sub-steps **in order**:
1. **Import Data** — Reads raw CSV/XLS files from `Raw data/ATDGCMS` and `Raw data/HPLCTOFMS`.
2. **Group Compounds** — Bins compounds by RT and m/z. Adjust thresholds if needed:
   - *GC RT threshold*: Default 0.005 min
   - *HPLC RT thresholds*: Default 0.4 / 0.28 min (split at 6 min)
   - *m/z thresholds*: Default 0.05 (GC) / 0.0005 (HPLC)
3. **Remove Benchmarks** — Strips internal-standard compounds (toluene for GC; 5 known m/z for HPLC).

### Step 2 — Blank Subtraction
- Adjustable parameter: **SD Multiplier** (default 3×). Values below `mean_blank + N × SD_blank` are removed.
- Click **Run Blank Subtraction**.

### Step 3 — Source Assignment
- Automatically labels samples as **Store-Bought** or **Environmental** based on filename patterns (`USE` → Environmental).
- Click **Run Step 3**.

### Step 4 — Trace Metal Import
- Reads ICP-MS Excel files from `Raw data/ICPMS_Trace metal/`.
- Adjustable: **Seawater metals to discard** (default: Na, Mg, K, Ca).
- Click **Import Metals**.

### Step 5 — Label Merge
- Merges `Plastic_type`, `Subcategory`, and `Polymer` labels from the sample-labelling Excel file.
- Creates all dataset combinations: gc, hplc, icp, gc_hplc, gc_icp, hplc_icp, gc_hplc_icp.
- Click **Run Step 5**.

### Step 6 — Feature Reduction & Data Fusion
Two operations:

1. **Feature Reduction** (6.1–6.3):
   - Select a dataset and a filtering method:
     - *Regular 80% rule*: remove features missing in >80% of all samples.
     - *Modified 80% rule*: keep features present ≥80% in *any* plastic type.
     - *In-house*: 4-rule comprehensive filter.
   - Click **Run Feature Reduction**.

2. **Data Fusion** (6.4):
   - Applies in-house filtering to GC, HPLC, ICP individually, then fuses the multi-instrument datasets.
   - Click **Run Data Fusion**.

### Step 7 — Hybrid ML Pipeline (⚠ takes several minutes)
Configure the full ML optimization:

| Parameter | Default | Description |
|---|---|---|
| Training approach | SB → ENV | Option 1 or Option 2 |
| Datasets | gc, gc_icp | Which cleaned datasets to analyze |
| Imputation | half_min, median, knn | Methods to search over |
| Normalization | none, log, log10, zscore, pareto | Methods to search over |
| Algorithms | ranger, svmRadial, xgbTree | ML classifiers |
| Top % features | 0.20 | Features considered "top" per pipeline |
| Stability threshold | 0.60 | Min stability for robust features |
| Mean rank percentile | 0.30 | Max rank percentile for robust features |
| CV folds | 5 | Cross-validation folds |
| Top-k pipelines | 5 | Ensemble size |
| Permutation reps | 10 | Permutation importance iterations |

Click **Run Hybrid Analysis** and wait for the progress bar. Results will appear in the following sub-tabs:
- **Best Pipeline** — imputation / normalization / algorithm
- **Metrics** — MCC, accuracy, balanced accuracy, F1, kappa
- **Confusion Matrix** — publication-ready heatmap
- **Feature Tiers** — stability-ranked bar chart
- **Top Features** — sortable data table

### Step 8 — Hierarchical Cluster Analysis
- Select which Step 7 results to use and configure distance (manhattan / euclidean / bray) and linkage (average / complete / ward.D2 / single).
- Click **Run HCA** to generate the dendrogram.
- Use **Download Dendrogram (PDF)** to export.

---

## Tips & Troubleshooting

| Issue | Solution |
|---|---|
| "Directory not found" | Check the project root path on the Welcome tab |
| "Helper functions file not found" | Ensure the `.R` helper file is in the same directory as `app.R` |
| Step N fails with "Run Step N-1 first" | Steps must be run sequentially — run earlier steps first |
| Step 7 is very slow | Reduce the number of imputation/normalization/algorithm combinations, or lower `top_k` / `n_permutations` |
| Out-of-memory error | Close other applications; reduce dataset size; or increase R's memory limit with `memory.limit(size = ...)` on Windows |
| Package install failure | Run `install.packages("package_name")` manually in R to see the full error |

---

## License

AGPL-3.0-only — same as the original workflow code.
