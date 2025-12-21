# Microplastic Fingerprinting

This repository accompanies the publication: **"Chemical Fingerprints of New vs. Weathered Microplastics: A Machine Learning Approach"**. 

It contains the computational workflow and datasets used to analyze chemical fingerprints of microplastics using various machine learning techniques, including Random Forest and specialized imputation/normalization methods.

## 📂 Repository Structure

* **`Microplastic workflow_20Dec2025_clean4Github.Rmd`**: The main R Markdown workflow file. It handles data import, processing, feature selection, and model training.
* **`Helper Functions for Computational Workflow.R`**: A collection of custom functions used by the main workflow for tasks like compound grouping, imputation, and feature filtering.
* **`Raw Data/`**: Directory containing the raw data files of instrumental analyses: ATD-GC-MS, HPLC-QTOF-MS and ICP-MS/MS.

## 💻 System Specifications
The workflow was developed and tested on the following system configuration:

* **Operating System**: Windows 11 Pro 64-bit (10.0, Build 26200)
* **Processor**: 11th Gen Intel(R) Core(TM) i7-11800H @ 2.30GHz (16 CPUs), ~2.3GHz
* **Memory**: 16 GB RAM (16384MB)

## 📦 Dependencies & Installation
To ensure compatibility and reproducibility, a list of used package versions is generated in `package_versions.txt`.

### Key Libraries
* **Data Manipulation**: `tidyverse`, `dplyr`, `data.table`, `readxl`, `writexl`
* **Machine Learning**: `caret`, `randomForest`, `ranger`, `e1071`, `pROC`
* **Imputation**: `missForest`, `VIM`, `softImpute`
* **Visualization**: `ggplot2`, `pheatmap`, `ggpubr`, `gridExtra`

### Version Tracking
To avoid compatibility issues, please refer to the `R package versions.txt` file (if available) to see the exact package versions used. You can generate this file by running the following command in R:

```r
writeLines(capture.output(sessionInfo()), "session_info.txt")
# OR for a cleaner list of attached packages:
write.table(installed.packages()[,c("Package", "Version")], 
            file = "R package versions.txt", 
            row.names = FALSE, quote = FALSE)
