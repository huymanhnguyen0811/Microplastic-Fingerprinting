# ==============================================================================
# AUTOMATED ANALYSES FOR MULTI-CLASS MODEL SELECTION
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP AND DATA LOADING
# ------------------------------------------------------------------------------
# Install packages if they are not already installed
required_packages <- c("e1071", "FNN", "igraph", "infotheo", "car", "caret")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(e1071)    # For Linear SVM (Linearity Check)
library(FNN)      # For fast k-Nearest Neighbors
library(igraph)   # For Minimum Spanning Tree
library(infotheo) # For Mutual Information
library(car)      # For VIF
library(caret)    # For general data handling
library(dplyr)
library(tidyr)

# --- LOAD DATA (REPLACE THIS WITH YOUR DATA) ---
# NOTE: Ensure your target column is a factor and is the LAST column for this script
data <- hplc_combined_SB_ENV_shared_cols_ntr_clean %>% dplyr::select(-c(Subcategory, technique, Source, Polymer))
target_col_name <- "Plastic_type" 

# 0.000000001 imputation
for (r in 1:nrow(data)) { 
  data[r, which(base::is.na(data[r,]))] <- 0.000000001
}

# Percentage normalization
data_percentage <- as.data.frame(t(apply(data[, 3:ncol(data)],
                                                             MARGIN = 1, 
                                                             function(row) {row/sum(row, na.rm = TRUE)}))) 

data_percentage <- data_percentage %>%
  mutate(Plastic_type = data$Plastic_type) %>%
  relocate(Plastic_type, .before = 1)

# Separate Features (X) and Target (y)
X <- data_percentage[, -which(names(data_percentage ) == target_col_name)]
y <- data_percentage[[target_col_name]]

# Ensure target is a factor
if(!is.factor(y)) { y <- as.factor(y) }

cat("========================================================\n")
cat(" DATASET SUMMARY \n")
cat("========================================================\n")
cat("Features:", ncol(X), "\n")
cat("Samples: ", nrow(X), "\n")
cat("Classes: ", length(unique(y)), "\n\n")


# ------------------------------------------------------------------------------
# 1. DATA COMPLEXITY: LINEAR SEPARABILITY & FISHER'S RATIO
# ------------------------------------------------------------------------------
cat("--- 1. DATA COMPLEXITY MEASURES ---\n")

# A. Linear Separability Index (Using Linear SVM Error)
svm_model <- svm(x = X, y = y, kernel = "linear", cross = 5, cost = 1)
linear_error <- 1 - (svm_model$tot.accuracy / 100)

cat(sprintf("Linearity Defect (Linear SVM Error): %.4f\n", linear_error))
if(linear_error < 0.10) {
  cat("  -> INTERPRETATION: High Linearity. Recommended: Logistic Regression, Linear SVM.\n")
} else {
  cat("  -> INTERPRETATION: Non-Linear. Recommended: Random Forest, XGBoost, Neural Nets.\n")
}

# B. Maximum Fisher's Discriminant Ratio (F1) -- UPDATED LOGIC
calculate_f1 <- function(feature, target) {
  means <- tapply(feature, target, mean)
  overall_mean <- mean(feature)
  # Between-class scatter
  n_classes <- table(target)
  sb <- sum(n_classes * (means - overall_mean)^2)
  # Within-class scatter
  vars <- tapply(feature, target, var)
  vars[is.na(vars)] <- 0 
  sw <- sum((n_classes - 1) * vars)
  
  if (sw == 0) return(Inf)
  return(sb / sw)
}

f1_scores <- apply(X, 2, calculate_f1, target = y)
max_f1 <- max(f1_scores)
best_feature <- names(f1_scores)[which.max(f1_scores)]

cat(sprintf("\nMax Fisher Ratio (F1): %.4f (Best Feature: %s)\n", max_f1, best_feature))

# --- UPDATED F1 INTERPRETATION GUIDE ---
if (max_f1 < 1.0) {
  cat("  -> INTERPRETATION (Low < 1.0): Very Poor Separability.\n")
  cat("     Classes heavily overlap. Single features are weak.\n")
  cat("     RECOMMENDED: Deep Learning, XGBoost (to capture complex interactions).\n")
} else if (max_f1 >= 1.0 && max_f1 <= 3.0) {
  cat("  -> INTERPRETATION (Moderate 1.0-3.0): Ambiguous/Moderate Separability.\n")
  cat("     Classes are distinct but edges touch/mix.\n")
  cat("     RECOMMENDED: Random Forest, SVM (RBF Kernel).\n")
} else {
  cat("  -> INTERPRETATION (High > 3.0): Strong Separability.\n")
  cat("     At least one feature separates classes well.\n")
  cat("     RECOMMENDED: Logistic Regression, Linear SVM, LDA (Avoid Deep Learning).\n")
}


# ------------------------------------------------------------------------------
# 2. CLASS OVERLAP: kNN ERROR & MST
# ------------------------------------------------------------------------------
cat("\n--- 2. CLASS OVERLAP ANALYSIS ---\n")

# A. 1-NN Leave-One-Out Error (Estimating Bayes Error)
knn_res <- knn.cv(train = X, cl = y, k = 1)
knn_error <- mean(as.character(knn_res) != as.character(y))

cat(sprintf("1-NN Cross-Validation Error (Bayes Error Est.): %.4f\n", knn_error))
if(knn_error > 0.30) {
  cat("  -> INTERPRETATION: High Overlap. Data may contain high noise or insufficient features.\n")
}

# B. Minimum Spanning Tree (MST) Class Interleaving
dist_matrix <- as.matrix(dist(X))
g <- graph_from_adjacency_matrix(dist_matrix, mode = "undirected", weighted = TRUE)
mst <- mst(g)

# FIX 1: Use 'names = FALSE' to ensure we get row numbers, not row names
edges <- as_edgelist(mst, names = FALSE) 

# FIX 2: Check for NAs just in case the dataset has missing values
if(any(is.na(edges))) { stop("Edges contain NA values. Please check distance matrix.") }

class_start <- y[edges[,1]]
class_end <- y[edges[,2]]

# FIX 3: Use na.rm = TRUE in sum to be safe
cross_class_edges <- sum(class_start != class_end, na.rm = TRUE)
total_edges <- nrow(edges)
interleaving_ratio <- cross_class_edges / total_edges

cat(sprintf("MST Interleaving Ratio: %.4f (%.0f/%0.f edges cross classes)\n", 
            interleaving_ratio, cross_class_edges, total_edges))

if(is.na(interleaving_ratio)) {
  cat("  -> ERROR: Could not calculate ratio (possible NAs in target variable).\n")
} else if(interleaving_ratio > 0.10) {
  cat("  -> INTERPRETATION: Classes are interleaved (checkerboard). Recommended: Random Forest, kNN.\n")
} else {
  cat("  -> INTERPRETATION: Classes are clustered. Recommended: SVM, LDA.\n")
}


# ------------------------------------------------------------------------------
# 3. FEATURE RELEVANCE & REDUNDANCY (UPDATED WITH NORMALIZED MI)
# ------------------------------------------------------------------------------
cat("\n--- 3. FEATURE ANALYSIS ---\n")

# A. Mutual Information (MI)
# Discretize continuous data for MI calculation
X_disc <- discretize(X, disc = "equalfreq", nbins = 5)

# Calculate raw Mutual Information for each feature
raw_mi_scores <- sapply(X_disc, function(col) mutinformation(col, y))

# Calculate Entropy of the Target (The maximum possible information)
target_entropy <- entropy(y)

# Calculate Normalized MI (Percentage of uncertainty removed)
# Formula: NMI = I(X;Y) / H(Y)
normalized_mi_scores <- raw_mi_scores / target_entropy
avg_nmi <- mean(normalized_mi_scores)

cat(sprintf("Target Entropy: %.4f\n", target_entropy))
cat(sprintf("Average Normalized MI: %.4f (Features explain ~%.1f%% of the target)\n", 
            avg_nmi, avg_nmi * 100))

if(avg_nmi < 0.05) {
  cat("  -> INTERPRETATION (< 5%): Weak Signal. \n")
  cat("     Features individually explain very little.\n")
  cat("     RECOMMENDED: Deep Learning, XGBoost (Needs to aggregate many weak features).\n")
} else if(avg_nmi >= 0.05 && avg_nmi <= 0.20) {
  cat("  -> INTERPRETATION (5-20%): Moderate Signal. \n")
  cat("     Features have decent predictive power.\n")
  cat("     RECOMMENDED: Random Forest, SVM.\n")
} else {
  cat("  -> INTERPRETATION (> 20%): Strong Signal. \n")
  cat("     Features are highly predictive.\n")
  cat("     RECOMMENDED: Decision Trees, Logistic Regression (Simple models may suffice).\n")
}

# B. Multicollinearity (VIF)
fake_num_target <- as.numeric(y) 
vif_values <- tryCatch({
  lm_model <- lm(fake_num_target ~ ., data = as.data.frame(X))
  vif(lm_model)
}, error = function(e) {
  return("Could not compute VIF (singular matrix or too few samples).")
})

if(is.numeric(vif_values)) {
  max_vif <- max(vif_values)
  cat(sprintf("\nMax VIF Score: %.2f\n", max_vif))
  if(max_vif > 10) {
    cat("  -> INTERPRETATION: High Multicollinearity. Recommended: Ridge/Lasso, Random Forest.\n")
  } else {
    cat("  -> INTERPRETATION: Low Redundancy. Linear models are safe.\n")
  }
} else {
  cat("VIF Result:", vif_values, "\n")
}

# ==============================================================================
# DIAGNOSING THE "SINGULAR MATRIX" ERROR
# ==============================================================================

# 1. Check for Perfect Linear Dependencies (The "Alias" Problem)
# We fit a linear model to see which coefficients get "NA" (meaning they are redundant)
# We use a dummy target (1, 2, 3...) just to check the X matrix structure
dummy_y <- 1:nrow(X)
alias_check <- lm(dummy_y ~ ., data = as.data.frame(X))

# Find linearly dependent coefficients (those that are NA)
redundant_vars <- names(which(is.na(coef(alias_check))))

if(length(redundant_vars) > 0) {
  cat("--- CAUSE FOUND: Perfect Multicollinearity ---\n")
  cat("The following variables are exact copies/combinations of others:\n")
  print(redundant_vars)
  cat("\nRECOMMENDATION: Remove these columns or use Random Forest/XGBoost.\n")
} else {
  cat("--- CAUSE FOUND: High Dimensionality or Near-Perfect Correlation ---\n")
  cat("No exact duplicates found, but correlation is extremely high.\n")
}

# 2. Check correlations (just to be sure)
# We look for correlations > 0.99
cor_matrix <- cor(X, use = "pairwise.complete.obs")
# Set diagonal to 0 to ignore self-correlation
diag(cor_matrix) <- 0
high_cor_pairs <- which(abs(cor_matrix) > 0.99, arr.ind = TRUE)

if(nrow(high_cor_pairs) > 0) {
  cat("\n--- HIGH CORRELATION PAIRS (> 0.99) ---\n")
  # Print the first few pairs
  print(head(high_cor_pairs))
}

# ------------------------------------------------------------------------------
# 4. DIMENSIONALITY VS SAMPLE SIZE
# ------------------------------------------------------------------------------
cat("\n--- 4. N/P RATIO ANALYSIS ---\n")

N <- nrow(X)
P <- ncol(X)
ratio <- N / P

cat(sprintf("N (Samples) = %d, P (Features) = %d\n", N, P))
cat(sprintf("Ratio (N/P): %.2f\n", ratio))

if(ratio < 5) {
  cat("  -> INTERPRETATION: High Dimensionality (P > N/5). Risk of overfitting.\n")
  cat("     Recommended: Linear SVM, Regularized Regression (Lasso). Avoid Deep Learning.\n")
} else if (ratio > 50) {
  cat("  -> INTERPRETATION: Data Rich. Deep Learning and Gradient Boosting likely effective.\n")
} else {
  cat("  -> INTERPRETATION: Balanced. Most standard algorithms (RF, SVM) will work.\n")
}

cat("\n========================================================\n")
cat(" ANALYSIS COMPLETE \n")
cat("========================================================\n")