#' ============================================================================
#' UNIFIED DIAGNOSTIC TESTS FOR ML MODEL SELECTION
#' ============================================================================
#' 
#' This function provides automated diagnostic tests to help select appropriate
#' ML algorithms for multi-class classification problems.
#' 
#' Author: Based on Microplastic-Fingerprinting workflow
#' Date: January 2026
#' 
#' ============================================================================
#' PENALTY SCORING RATIONALE
#' ============================================================================
#' 
#' SCORING PHILOSOPHY:
#' - All algorithms start with a base score of 100
#' - Penalties are applied when assumptions are violated
#' - Boosts are applied when conditions favor an algorithm
#' - Final ranking determines recommended algorithm
#' 
#' PENALTY MAGNITUDE GUIDELINES:
#' - SEVERE penalty (-40 to -50): Core mathematical assumption directly violated
#' - MAJOR penalty (-20 to -35): Important assumption violated, significant performance impact
#' - MODERATE penalty (-10 to -19): Assumption partially violated, noticeable impact
#' - MINOR penalty (-3 to -9): Slight assumption concern, minor impact
#' - NEUTRAL (0): No effect on this algorithm
#' - MINOR boost (+3 to +9): Algorithm slightly favored in this condition
#' - MAJOR boost (+10 to +20): Algorithm well-suited for this condition
#' 
#' ALGORITHM CHARACTERISTICS SUMMARY:
#' 
#' 1. LDA (Linear Discriminant Analysis):
#'    - Assumes: Multivariate normality, equal covariance matrices, linear boundaries
#'    - Math: Uses pooled covariance matrix (p(p+1)/2 parameters)
#'    - Sensitive to: Non-normality, unequal covariances, outliers, collinearity
#' 
#' 2. QDA (Quadratic Discriminant Analysis):
#'    - Assumes: Multivariate normality per class (but NOT equal covariances)
#'    - Math: Estimates k separate covariance matrices (k × p(p+1)/2 parameters)
#'    - More parameters than LDA → needs more samples, more sensitive to sample size
#'    - Advantage: Can model quadratic decision boundaries
#' 
#' 3. Logistic Regression:
#'    - Assumes: Linear relationship between log-odds and features
#'    - No normality assumption, but MLE estimation affected by extreme values
#'    - Sensitive to: Multicollinearity (inflates SE), separation (complete/quasi)
#' 
#' 4. Naive Bayes (Gaussian):
#'    - Assumes: Features are INDEPENDENT given class (P(X|Y) = ∏P(Xi|Y))
#'    - Assumes: Each feature follows normal distribution within class
#'    - Severely violated by: Multicollinearity (destroys independence assumption)
#' 
#' 5. k-NN (k-Nearest Neighbors):
#'    - Non-parametric, distance-based
#'    - Sensitive to: Curse of dimensionality, outliers, scale, feature relevance
#'    - No distributional assumptions but distance metrics affected by distribution shape
#' 
#' 6. SVM Linear:
#'    - Finds maximum margin hyperplane
#'    - Works well in high-D, robust to some outliers (soft margin)
#'    - Assumes: Linear separability (or uses slack variables)
#' 
#' 7. SVM RBF (Radial Basis Function):
#'    - Non-linear kernel maps to high-D space
#'    - More flexible than linear, can model complex boundaries
#'    - More robust to non-linearity, less interpretable
#' 
#' 8. Random Forest:
#'    - Ensemble of decision trees with bagging
#'    - Trees use rank-order (splits), invariant to monotonic transformations
#'    - Robust to: Outliers, skewness, collinearity, non-linearity
#'    - Needs: Sufficient samples per class for meaningful splits
#' 
#' 9. XGBoost:
#'    - Gradient boosted trees
#'    - Similar robustness to RF but boosting more prone to overfitting
#'    - Excellent performance but needs more careful tuning
#' 
#' ============================================================================

# Required packages - consolidated list
required_packages <- c(
  "tidyverse", "caret", "ranger", "e1071", "MASS", "class", "xgboost",
  "missForest", "VIM", "mice", "doParallel", "foreach", "mltools",
  "data.table", "yardstick", "pROC", "nortest", "MVN", "heplots",
  "corrplot", "psych", "vegan", "randomForest", "FactoMineR",
  "factoextra", "car", "moments", "cluster"
)

# Check and install missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#' select_ML_diag_tests: Unified diagnostic function for ML model selection
#' 
#' @param data A data frame with features (numeric columns) and a class column
#' @param class_col Name of the class/label column (character)
#' @param alpha Significance level for statistical tests (default = 0.05)
#' @param verbose Whether to print detailed diagnostics (default = TRUE)
#' @param return_recommendations Whether to return algorithm recommendations (default = TRUE)
#' 
#' @return A list containing:
#'   - diagnostics: All diagnostic test results
#'   - assumptions: Summary of assumption violations
#'   - recommendations: Ranked ML algorithm recommendations
#'   - recommended_algorithm: Single best algorithm recommendation
#'   
#' @examples
#' # results <- select_ML_diag_tests(my_data, class_col = "Species")
#' # best_algo <- results$recommended_algorithm

# select_ML_diag_tests <- function(data, 
#                                  class_col, 
#                                  alpha = 0.05,
#                                  verbose = TRUE,
#                                  return_recommendations = TRUE) {
#   
#   # ============================================================================
#   # INITIALIZATION AND DATA VALIDATION
#   # ============================================================================
#   
#   if (verbose) cat("\n", strrep("=", 70), "\n")
#   if (verbose) cat("AUTOMATED ML MODEL SELECTION DIAGNOSTIC TESTS\n")
#   if (verbose) cat(strrep("=", 70), "\n\n")
#   
#   # Validate inputs
#   if (!is.data.frame(data)) {
#     stop("Input 'data' must be a data frame")
#   }
#   
#   if (!class_col %in% names(data)) {
#     stop(paste("Class column '", class_col, "' not found in data"))
#   }
#   
#   # Separate features and class
#   class_vector <- as.factor(data[[class_col]])
#   feature_data <- data %>% 
#     dplyr::select(-all_of(class_col)) %>%
#     dplyr::select(where(is.numeric))
#   
#   n_samples <- nrow(feature_data)
#   n_features <- ncol(feature_data)
#   n_classes <- length(unique(class_vector))
#   class_counts <- table(class_vector)
#   
#   if (verbose) {
#     cat("Data Summary:\n")
#     cat(sprintf("  - Samples: %d\n", n_samples))
#     cat(sprintf("  - Features: %d\n", n_features))
#     cat(sprintf("  - Classes: %d\n", n_classes))
#     cat(sprintf("  - Class distribution: %s\n", 
#                 paste(names(class_counts), "=", class_counts, collapse = ", ")))
#     cat(sprintf("  - Sample-to-feature ratio: %.2f\n", n_samples / n_features))
#     cat("\n")
#   }
#   
#   # Initialize results storage
#   diagnostics <- list()
#   assumptions <- list()
#   scores <- data.frame(
#     Algorithm = c("Random Forest", "SVM (RBF)", "SVM (Linear)", "LDA", 
#                   "QDA", "k-NN", "XGBoost", "Naive Bayes", "Logistic Regression"),
#     Score = rep(100, 9)  # Start with perfect score, penalize violations
#   )
#   
#   # ============================================================================
#   # TEST 1: SAMPLE SIZE ADEQUACY
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Sample size affects algorithms differently based on their parameter estimation needs:
#   #
#   # - QDA (-50 critical, -30 warning): Most affected because it estimates k separate
#   #   covariance matrices, each requiring p(p+1)/2 parameters. With k classes and p 
#   #   features, needs ~10× more samples than LDA. When n/p < 5, covariance matrices
#   #   become singular or poorly conditioned.
#   #
#   # - LDA (-40 critical, -15 warning): Estimates one pooled covariance matrix
#   #   (p(p+1)/2 parameters) plus k class means (k×p parameters). Needs fewer samples
#   #   than QDA but still requires n/p > 5-10 for stable covariance estimation.
#   #
#   # - Logistic Regression (-30 critical, -10 warning): MLE estimation requires 
#   #   sufficient "events per variable" (EPV). Rule of thumb: 10-20 samples per 
#   #   predictor. Low EPV causes biased coefficients and separation problems.
#   #
#   # - k-NN (-25 critical, -15 warning): Requires dense neighborhoods for reliable 
#   #   voting. With few samples, neighbors may be far away and unrepresentative.
#   #   Also affected by curse of dimensionality: distance becomes meaningless
#   #   as dimensions increase relative to samples.
#   #
#   # - Naive Bayes (-20 critical, -10 warning): Estimates only p marginal 
#   #   distributions per class, not full covariance. Less demanding than LDA/QDA
#   #   but still needs enough samples for reliable probability estimates.
#   #
#   # - XGBoost (-15 critical, -10 warning): Boosting iteratively fits residuals,
#   #   prone to overfitting with small samples. More sensitive than RF because
#   #   errors compound across boosting iterations.
#   #
#   # - Random Forest (-10 critical, -5 warning): Bagging provides some protection
#   #   through averaging. Still needs sufficient samples for meaningful tree splits.
#   #   Minimum ~5-10 samples per terminal node for stable predictions.
#   #
#   # - SVM Linear (+10 critical, +5 warning): Maximum margin principle means only
#   #   support vectors matter, not all samples. Regularization prevents overfitting.
#   #   Actually BENEFITS in low-sample scenarios relative to parametric methods.
#   #
#   # - SVM RBF (+5 critical, +0 warning): Similar reasoning to Linear SVM, but
#   #   RBF has more flexibility (gamma parameter) which can cause overfitting
#   #   with very few samples if not tuned properly.
#   #
#   # Per-class sample size additionally penalizes algorithms that need sufficient
#   # samples within each class for reliable class-specific parameter estimation.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 1: Sample Size Adequacy\n", strrep("-", 40), "\n")
#   
#   min_class_size <- min(class_counts)
#   samples_per_feature <- n_samples / n_features
#   samples_per_class_per_feature <- min_class_size / n_features
#   
#   diagnostics$sample_size <- list(
#     n_samples = n_samples,
#     n_features = n_features,
#     n_classes = n_classes,
#     min_class_size = min_class_size,
#     samples_per_feature = samples_per_feature,
#     samples_per_class_per_feature = samples_per_class_per_feature
#   )
#   
#   if (samples_per_feature < 5) {
#     assumptions$sample_size <- "CRITICAL: Very low sample-to-feature ratio"
#     
#     # QDA: -50 (SEVERE) - Must estimate k covariance matrices, each p(p+1)/2 params
#     # With n/p < 5, matrices will be singular or extremely ill-conditioned
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 50
#     
#     # LDA: -40 (SEVERE) - Pooled covariance estimation becomes unstable
#     # Linear discriminant functions unreliable with poorly estimated covariance
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 40
#     
#     # Logistic Regression: -30 (MAJOR) - MLE suffers from separation and bias
#     # Coefficient estimates become unstable, SEs inflate dramatically
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 30
#     
#     # k-NN: -25 (MAJOR) - Sparse neighborhoods, unreliable distance-based voting
#     # Curse of dimensionality severe: distances become nearly equal in high-D
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 25
#     
#     # Naive Bayes: -20 (MODERATE-MAJOR) - Marginal distribution estimates unreliable
#     # Still better than LDA/QDA since doesn't estimate full covariance
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 20
#     
#     # XGBoost: -15 (MODERATE) - Boosting overfits easily with few samples
#     # Sequential error correction amplifies noise in small datasets
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] - 15
#     
#     # Random Forest: -10 (MODERATE) - Bagging helps but still needs samples for splits
#     # Each tree sees ~63% of data; with few samples, trees become similar
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] - 10
#     
#     # SVM Linear: +10 (BOOST) - Margin maximization works well in low-sample regime
#     # Only support vectors define boundary; regularization prevents overfitting
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 10
#     
#     # SVM RBF: +5 (MINOR BOOST) - Still benefits from margin principle
#     # Slightly less boost than linear due to gamma parameter tuning risk
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#     
#   } else if (samples_per_feature < 10) {
#     assumptions$sample_size <- "WARNING: Low sample-to-feature ratio"
#     
#     # Same ordering but reduced penalties - assumptions partially met
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 30
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 15
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 15
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] - 10
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] - 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <-
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 5
#     # SVM (RBF): neutral in this scenario
# 
#   } else {
#     assumptions$sample_size <- "OK: Adequate sample size"
#     
#     # With adequate samples, all methods can perform well
#     # Small boosts to methods that particularly benefit from more data
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 5
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 3
#   }
#   
#   # Additional penalty for small per-class samples
#   # Many algorithms need sufficient samples WITHIN each class
#   if (min_class_size < 10) {
#     # QDA: -25 (additional) - Cannot estimate class-specific covariance with <10 samples
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 25
#     
#     # XGBoost: -20 - Boosting for minority class becomes unstable
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] - 20
#     
#     # Random Forest: -15 - Tree splits for small classes unreliable
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] - 15
#     
#     # k-NN: -15 - Minority class has too few neighbors
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 15
#     
#     # Logistic Regression: -15 - Separation likely with sparse classes
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 15
#     
#     # LDA: -10 - Class mean estimation unreliable
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 10
#     
#     # SVM RBF: -10 - Few support vectors available per class
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] - 10
#     
#     # Naive Bayes: -10 - Class-conditional probability estimates unreliable
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#     
#     # SVM Linear: -5 - Less affected due to simpler boundary
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 5
#     
#   } else if (min_class_size < 20) {
#     # Minor penalties for borderline cases
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 10
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] - 5
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 5
#   }
#   
#   if (verbose) {
#     cat(sprintf("  Samples per feature: %.2f (recommend >= 10)\n", samples_per_feature))
#     cat(sprintf("  Min class samples per feature: %.2f\n", samples_per_class_per_feature))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$sample_size))
#   }
#   
#   # ============================================================================
#   # TEST 2: UNIVARIATE NORMALITY (for each feature)
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Tests whether individual features follow normal distributions. While 
#   # multivariate normality (Test 3) is the true assumption, univariate 
#   # non-normality often indicates multivariate non-normality.
#   #
#   # - QDA (-25 when violated): More parameters to estimate than LDA, making it
#   #   more sensitive to distributional assumptions. Non-normal features lead to
#   #   biased covariance estimates and unreliable class-specific distributions.
#   #
#   # - LDA (-20 when violated): Discriminant functions derived assuming normality.
#   #   The optimal linear boundary (Fisher's criterion) is only optimal under
#   #   normality. Non-normality means the "optimal" boundary isn't actually optimal.
#   #
#   # - Naive Bayes (-15 when violated): Gaussian Naive Bayes explicitly models
#   #   P(Xi|Y) ~ N(μ_k, σ_k). Non-normal features make these probability estimates
#   #   incorrect, though the independence assumption often dominates errors.
#   #
#   # - Logistic Regression (-10 when violated): No normality assumption in theory,
#   #   but extreme skewness can cause numerical instability in MLE optimization
#   #   and affect coefficient interpretation.
#   #
#   # - k-NN (+5 when violated): Non-parametric, no distributional assumptions.
#   #   Actually slightly BENEFITS because it naturally adapts to local density.
#   #
#   # - SVM RBF (+5 when violated): Kernel methods are distribution-free.
#   #   RBF kernel projects to infinite-dimensional space where any distribution
#   #   can be separated.
#   #
#   # - Random Forest (+10 when violated): Trees use rank-order for splits, not
#   #   actual values. Completely invariant to monotonic transformations.
#   #   log(x), x^2, or any non-normal transformation doesn't affect performance.
#   #
#   # - XGBoost (+10 when violated): Same reasoning as Random Forest.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 2: Univariate Normality (Shapiro-Wilk)\n", strrep("-", 40), "\n")
#   
#   normality_results <- data.frame(
#     Feature = character(),
#     W_statistic = numeric(),
#     p_value = numeric(),
#     Normal = logical(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (col in names(feature_data)) {
#     # Shapiro-Wilk test (limited to 5000 samples)
#     test_data <- feature_data[[col]]
#     if (length(test_data) > 5000) {
#       test_data <- sample(test_data, 5000)
#     }
#     
#     # Skip if too few unique values
#     if (length(unique(test_data)) < 3) {
#       normality_results <- rbind(normality_results, data.frame(
#         Feature = col,
#         W_statistic = NA,
#         p_value = NA,
#         Normal = NA
#       ))
#       next
#     }
#     
#     sw_test <- tryCatch(
#       shapiro.test(test_data),
#       error = function(e) list(statistic = NA, p.value = NA)
#     )
#     
#     normality_results <- rbind(normality_results, data.frame(
#       Feature = col,
#       W_statistic = as.numeric(sw_test$statistic),
#       p_value = sw_test$p.value,
#       Normal = ifelse(is.na(sw_test$p.value), NA, sw_test$p.value > alpha)
#     ))
#   }
#   
#   # Bug fix: Handle case where all normality tests return NA (division by zero)
#   non_na_count <- sum(!is.na(normality_results$Normal))
#   pct_normal <- if (non_na_count == 0) {
#     50  # Default to neutral when no tests could be computed
#   } else {
#     sum(normality_results$Normal, na.rm = TRUE) / non_na_count * 100
#   }
# 
#   diagnostics$univariate_normality <- normality_results
#   diagnostics$pct_normal_features <- pct_normal
# 
#   if (pct_normal < 50) {
#     assumptions$univariate_normality <- "VIOLATED: Most features non-normal"
#     
#     # QDA: -25 - Covariance estimation and class distributions unreliable
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 25
#     
#     # LDA: -20 - Discriminant functions no longer optimal
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 20
#     
#     # Naive Bayes: -15 - Gaussian density assumption violated
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 15
#     
#     # Logistic Regression: -10 - Numerical stability and interpretation affected
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     
#     # SVM Linear: 0 - No distributional assumptions, but linear boundary less flexible
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
#     
#     # k-NN: +5 - Non-parametric, adapts to local density
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
#     
#     # SVM RBF: +5 - Kernel methods are distribution-free
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#     
#     # Random Forest: +10 - Only uses rank order, completely robust
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 10
#     
#     # XGBoost: +10 - Same as Random Forest
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 10
#     
#   } else if (pct_normal < 75) {
#     assumptions$univariate_normality <- "PARTIAL: Some features non-normal"
#     
#     # Reduced penalties/boosts for partial violation
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 15
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 10
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 8
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     
#   } else {
#     assumptions$univariate_normality <- "OK: Most features normally distributed"
#     
#     # Parametric methods get boosted when their assumptions are met
#     # LDA: +10 - Discriminant functions are optimal under normality
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 10
#     
#     # Naive Bayes: +8 - Gaussian density assumption validated
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 8
#     
#     # QDA: +8 - Same as LDA but slightly less boost (more parameters = more risk)
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 8
#     
#     # Logistic Regression: +5 - Slight numerical stability benefit
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     
#     # SVM Linear: +3 - Linear boundary often works well with normal data
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     
#     # Non-parametric methods: neutral - they work regardless
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 0
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 0
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 0
#   }
#   
#   if (verbose) {
#     cat(sprintf("  Features tested: %d\n", nrow(normality_results)))
#     cat(sprintf("  Normal features: %.1f%%\n", pct_normal))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$univariate_normality))
#   }
#   
#   # ============================================================================
#   # TEST 3: MULTIVARIATE NORMALITY (Mardia's test)
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Tests whether the joint distribution of all features is multivariate normal.
#   # This is THE key assumption for LDA/QDA, not just univariate normality.
#   #
#   # Mardia's test examines both multivariate skewness and kurtosis:
#   # - Skewness: Tests if data are symmetrically distributed
#   # - Kurtosis: Tests if tails match MVN expectations
#   #
#   # - QDA (-20): Estimates class-specific MVN distributions. If true distribution
#   #   isn't MVN, the quadratic boundaries will be suboptimal.
#   #
#   # - LDA (-15): Uses Mahalanobis distance which is only optimal under MVN.
#   #   Also, Fisher's linear discriminant is derived assuming MVN.
#   #
#   # - Naive Bayes (-10): Not directly MVN assumption, but non-MVN often means
#   #   complex dependencies that violate independence assumption.
#   #
#   # - Logistic Regression (-8): Not strictly required, but MVN often means
#   #   linear log-odds relationship holds better.
#   #
#   # - Tree/SVM methods (+5 to +8): Benefit when parametric assumptions fail.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 3: Multivariate Normality (Mardia's Test)\n", strrep("-", 40), "\n")
#   
#   # Select subset of features if too many (MVN test is computationally expensive)
#   mvn_features <- feature_data
#   if (ncol(mvn_features) > 20) {
#     # Use top 20 features by variance
#     variances <- apply(mvn_features, 2, var, na.rm = TRUE)
#     top_vars <- names(sort(variances, decreasing = TRUE))[1:20]
#     mvn_features <- mvn_features[, top_vars]
#   }
#   
#   # Remove any rows with NA
#   mvn_features <- na.omit(mvn_features)
#   
#   mvn_result <- tryCatch({
#     MVN::mvn(mvn_features, mvnTest = "mardia", univariatePlot = "none",
#              multivariatePlot = "none")
#   }, error = function(e) {
#     list(multivariateNormality = data.frame(
#       Test = c("Mardia Skewness", "Mardia Kurtosis"),
#       Statistic = c(NA, NA),
#       `p value` = c(NA, NA),
#       Result = c("Could not compute", "Could not compute")
#     ))
#   })
#   
#   diagnostics$multivariate_normality <- mvn_result$multivariateNormality
#   
#   mvn_ok <- tryCatch({
#     all(mvn_result$multivariateNormality$Result == "YES")
#   }, error = function(e) FALSE)
#   
#   if (!mvn_ok) {
#     assumptions$multivariate_normality <- "VIOLATED: Data not multivariate normal"
#     
#     # QDA: -20 - MVN is core assumption for class-specific distributions
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 20
#     
#     # LDA: -15 - Mahalanobis-based classification becomes suboptimal
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 15
#     
#     # Naive Bayes: -10 - Non-MVN suggests complex feature interactions
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#     
#     # Logistic Regression: -8 - Linear log-odds may not hold
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 8
#     
#     # SVM Linear: +3 - Margin principle doesn't require MVN
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     
#     # k-NN: +5 - Non-parametric, density-adaptive
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
#     
#     # SVM RBF: +5 - Kernel methods handle arbitrary distributions
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#     
#     # RF/XGBoost: +8 - Trees completely robust to distributional form
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 8
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 8
#     
#   } else {
#     assumptions$multivariate_normality <- "OK: Data is multivariate normal"
#     
#     # LDA: +15 - This is THE assumption LDA is built on
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 15
#     
#     # QDA: +10 - Also benefits but LDA might be simpler (Occam's razor)
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 10
#     
#     # Naive Bayes: +8 - Normality supports the Gaussian version
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 8
#     
#     # Logistic Regression: +5 - Numerical stability benefits
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     
#     # SVM Linear: +5 - Linear boundaries often work well with MVN data
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 5
#     
#     # Non-parametric methods: neutral
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 0
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 0
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 0
#   }
#   
#   if (verbose) {
#     cat("  Mardia's Test Results:\n")
#     print(mvn_result$multivariateNormality)
#     cat(sprintf("  Assessment: %s\n\n", assumptions$multivariate_normality))
#   }
#   
#   # ============================================================================
#   # TEST 4: HOMOGENEITY OF COVARIANCE MATRICES (Box's M test)
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Tests whether all classes have the same covariance matrix (Σ₁ = Σ₂ = ... = Σₖ).
#   # This is the FUNDAMENTAL difference between LDA and QDA.
#   #
#   # - LDA (-20 when violated): DIRECTLY assumes equal covariances. Uses pooled
#   #   Σ = Σ(nₖ-1)Σₖ / (n-k). When violated, the pooled estimate is a poor
#   #   representation of any class, leading to suboptimal boundaries.
#   #
#   # - QDA (+10 when violated): Designed SPECIFICALLY for unequal covariances!
#   #   QDA estimates separate Σₖ per class. When covariances truly differ,
#   #   QDA's quadratic boundaries are more appropriate than LDA's linear ones.
#   #
#   # - Logistic Regression (-10 when violated): Doesn't require equal covariances
#   #   mathematically, but heterogeneous covariances can make the linear
#   #   decision boundary suboptimal.
#   #
#   # - Naive Bayes (-5 when violated): Independence assumption more important,
#   #   but unequal variances within classes can affect probability estimates.
#   #
#   # - SVM/Tree/k-NN: No covariance assumptions at all.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 4: Homogeneity of Covariance (Box's M Test)\n", strrep("-", 40), "\n")
#   
#   # Box's M test - use subset of features if needed
#   boxm_features <- feature_data
#   if (ncol(boxm_features) > 15) {
#     variances <- apply(boxm_features, 2, var, na.rm = TRUE)
#     top_vars <- names(sort(variances, decreasing = TRUE))[1:15]
#     boxm_features <- boxm_features[, top_vars]
#   }
#   
#   boxm_result <- tryCatch({
#     heplots::boxM(boxm_features, class_vector)
#   }, error = function(e) {
#     list(statistic = NA, p.value = NA, error = as.character(e))
#   })
#   
#   diagnostics$boxm_test <- boxm_result
#   
#   boxm_ok <- tryCatch({
#     !is.na(boxm_result$p.value) && boxm_result$p.value > alpha
#   }, error = function(e) NA)
#   
#   if (is.na(boxm_ok)) {
#     assumptions$homogeneity_covariance <- "UNKNOWN: Could not compute Box's M"
#     # No score adjustments when test couldn't be computed
#     
#   } else if (!boxm_ok) {
#     assumptions$homogeneity_covariance <- "VIOLATED: Heterogeneous covariance matrices"
#     
#     # LDA: -20 - Core assumption DIRECTLY violated
#     # Pooled covariance misrepresents all classes
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 20
#     
#     # Logistic Regression: -10 - Linear boundary may be suboptimal
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     
#     # Naive Bayes: -5 - Minor effect through variance estimates
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 5
#     
#     # QDA: +10 (BOOST!) - This is exactly when QDA shines
#     # It's designed for this situation
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 10
#     
#     # Non-parametric methods get minor boost
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 3
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
#     
#   } else {
#     assumptions$homogeneity_covariance <- "OK: Homogeneous covariance matrices"
#     
#     # LDA: +15 - Core assumption validated, LDA is optimal choice
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 15
#     
#     # QDA: +5 - Still works, but LDA is simpler (Occam's razor)
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 5
#     
#     # Logistic Regression: +5 - Linear boundary appropriate
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     
#     # Naive Bayes: +5 - Variance estimates reliable
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 5
#     
#     # SVM Linear: +3 - Linear boundary appropriate
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     
#     # Non-parametric: neutral
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 0
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 0
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 0
#   }
#   
#   if (verbose) {
#     if (!is.null(boxm_result$p.value)) {
#       cat(sprintf("  Box's M statistic: %.2f\n", boxm_result$statistic))
#       cat(sprintf("  p-value: %.4f\n", boxm_result$p.value))
#     }
#     cat(sprintf("  Assessment: %s\n\n", assumptions$homogeneity_covariance))
#   }
#   
#   # ============================================================================
#   # TEST 5: MULTICOLLINEARITY (VIF and Condition Index)
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Multicollinearity occurs when features are highly correlated, meaning they
#   # contain redundant information. This affects algorithms differently:
#   #
#   # - Logistic Regression (-25 severe, -10 moderate): MOST affected. Causes:
#   #   1) Inflated standard errors (coefficients unreliable)
#   #   2) Unstable coefficients (small data changes → large coef changes)
#   #   3) Near-singular Hessian → numerical optimization failures
#   #   The design matrix X'X becomes ill-conditioned.
#   #
#   # - Naive Bayes (-20 severe): DIRECTLY violates the independence assumption!
#   #   If X₁ and X₂ are highly correlated, P(X₁,X₂|Y) ≠ P(X₁|Y)P(X₂|Y).
#   #   This is the "naive" assumption, and collinearity destroys it.
#   #
#   # - QDA (-20 severe): Covariance matrix inversion becomes ill-conditioned.
#   #   Each class needs its own inversion → problems multiply.
#   #
#   # - LDA (-15 severe): Pooled covariance inversion ill-conditioned, but
#   #   only one matrix vs. k matrices for QDA, so slightly less severe.
#   #
#   # - k-NN (-15 severe): Correlated features double-count information in
#   #   distance calculations. If X₁ ≈ X₂, they both contribute to distance,
#   #   effectively doubling the weight of that information.
#   #
#   # - SVM Linear (-10 severe): Regularization helps, but feature weights
#   #   become unstable with high collinearity.
#   #
#   # - SVM RBF (-5 severe): RBF kernel less affected due to non-linear
#   #   projection, but still some impact.
#   #
#   # - Random Forest (+10 severe): BENEFITS from collinearity handling.
#   #   At each split, only √p features considered. Correlated features
#   #   provide backup options; if one is chosen, others are implicitly excluded.
#   #
#   # - XGBoost (+10 severe): Similar feature selection mechanism.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 5: Multicollinearity (Correlation & VIF)\n", strrep("-", 40), "\n")
#   
#   # Correlation matrix analysis
#   cor_matrix <- cor(feature_data, use = "pairwise.complete.obs")
#   cor_matrix[is.na(cor_matrix)] <- 0
#   
#   # Find highly correlated pairs
#   cor_threshold <- 0.9
#   high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & 
#                             upper.tri(cor_matrix), arr.ind = TRUE)
#   
#   n_high_cor <- nrow(high_cor_pairs)
#   pct_high_cor <- n_high_cor / (n_features * (n_features - 1) / 2) * 100
#   
#   # Condition index (eigenvalue-based)
#   # CI = sqrt(λ_max / λ_min), where λ are eigenvalues of correlation matrix
#   # CI > 30 indicates severe multicollinearity
#   eigen_vals <- eigen(cor_matrix)$values
#   condition_index <- sqrt(max(eigen_vals) / min(eigen_vals[eigen_vals > 1e-10]))
#   
#   diagnostics$multicollinearity <- list(
#     n_high_correlations = n_high_cor,
#     pct_high_correlations = pct_high_cor,
#     condition_index = condition_index,
#     correlation_matrix = cor_matrix
#   )
#   
#   if (condition_index > 30 || pct_high_cor > 20) {
#     assumptions$multicollinearity <- "SEVERE: High multicollinearity detected"
#     
#     # Logistic Regression: -25 - Inflated SEs, unstable coefficients
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 25
#     
#     # Naive Bayes: -20 - Independence assumption DIRECTLY violated
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 20
#     
#     # QDA: -20 - Multiple ill-conditioned matrix inversions
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 20
#     
#     # LDA: -15 - Single ill-conditioned inversion
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 15
#     
#     # k-NN: -15 - Distance metric double-counts correlated features
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 15
#     
#     # SVM Linear: -10 - Feature weights unstable
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 10
#     
#     # SVM RBF: -5 - Kernel projection reduces impact
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] - 5
#     
#     # RF/XGBoost: +10 - Feature subsampling handles collinearity naturally
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 10
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 10
#     
#   } else if (condition_index > 15 || pct_high_cor > 10) {
#     assumptions$multicollinearity <- "MODERATE: Some multicollinearity present"
#     
#     # Reduced penalties for moderate collinearity
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 10
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 8
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 8
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 5
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] - 3
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     
#   } else {
#     assumptions$multicollinearity <- "OK: Low multicollinearity"
#     
#     # All methods benefit from clean correlation structure
#     # Naive Bayes gets highest boost since independence is satisfied
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 8
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 5
#     scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 5
#     scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#     # RF/XGBoost: neutral - they handle collinearity either way
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 0
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 0
#   }
#   
#   if (verbose) {
#     cat(sprintf("  High correlation pairs (|r| > %.1f): %d\n", cor_threshold, n_high_cor))
#     cat(sprintf("  Condition index: %.2f (>30 is severe)\n", condition_index))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$multicollinearity))
#   }
#   
#   # ============================================================================
#   # TEST 6: CLASS SEPARABILITY (Linear vs Non-linear)
#   # ============================================================================
#   #
#   # RATIONALE:
#   # This test assesses whether classes can be separated by linear boundaries
#   # or require non-linear decision surfaces.
#   #
#   # The metric used is the ratio of between-class variance to within-class
#   # variance (Fisher's criterion in PCA space):
#   # - High ratio (>1.5): Classes are well-separated, likely linearly separable
#   # - Medium ratio (0.5-1.5): Partial overlap, non-linear methods may help
#   # - Low ratio (<0.5): Significant overlap, need flexible methods
#   #
#   # - LDA (+15 when good): LDA IS Fisher's Linear Discriminant! When classes
#   #   are linearly separable, LDA is theoretically optimal.
#   #
#   # - SVM Linear (+15 when good): Maximum margin linear classifier excels.
#   #
#   # - Logistic Regression (+10 when good): Linear log-odds boundary works well.
#   #
#   # - LDA (-15 when poor): Linear boundary cannot separate overlapping classes.
#   #
#   # - SVM RBF (+15 when poor): RBF kernel can model arbitrary boundaries.
#   #
#   # - RF/XGBoost (+15 when poor): Tree ensembles create step-function
#   #   approximations to any decision boundary.
#   #
#   # - k-NN (+10 when poor): Adapts to local density, naturally non-linear.
#   #
#   # - QDA: Intermediate - can capture quadratic (curved) boundaries.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 6: Class Separability Analysis\n", strrep("-", 40), "\n")
#   
#   # PCA-based separability
#   pca_result <- tryCatch({
#     prcomp(feature_data, scale. = TRUE, center = TRUE)
#   }, error = function(e) NULL)
#   
#   linear_separability <- NA
#   
#   if (!is.null(pca_result)) {
#     # Calculate between-class vs within-class variance in PC space
#     pc_scores <- pca_result$x[, 1:min(5, ncol(pca_result$x))]
#     
#     # Simple separability metric using Fisher's criterion
#     class_centroids <- aggregate(pc_scores, by = list(class_vector), mean)
#     overall_centroid <- colMeans(pc_scores)
#     
#     # Between-class variance: weighted sum of squared distances from grand mean
#     between_var <- sum(sapply(1:n_classes, function(i) {
#       ni <- class_counts[i]
#       sum((as.numeric(class_centroids[i, -1]) - overall_centroid)^2) * ni
#     })) / n_samples
#     
#     # Within-class variance: sum of squared distances from class means
#     within_var <- sum(sapply(levels(class_vector), function(cls) {
#       class_data <- pc_scores[class_vector == cls, , drop = FALSE]
#       class_centroid <- colMeans(class_data)
#       sum(apply(class_data, 1, function(x) sum((x - class_centroid)^2)))
#     })) / n_samples
#     
#     # Separability ratio (Fisher's criterion)
#     linear_separability <- between_var / (within_var + 1e-10)
#   }
#   
#   diagnostics$class_separability <- list(
#     linear_separability_ratio = linear_separability,
#     pca_variance_explained = if(!is.null(pca_result)) 
#       summary(pca_result)$importance[2, 1:min(5, ncol(pca_result$x))] else NA
#   )
#   
#   if (!is.na(linear_separability)) {
#     if (linear_separability > 1.5) {
#       assumptions$class_separability <- "GOOD: Classes appear linearly separable"
#       
#       # Linear methods excel
#       # LDA: +15 - This IS Fisher's discriminant, optimal for linear separation
#       scores$Score[scores$Algorithm == "LDA"] <- 
#         scores$Score[scores$Algorithm == "LDA"] + 15
#       
#       # SVM Linear: +15 - Maximum margin linear classifier
#       scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#         scores$Score[scores$Algorithm == "SVM (Linear)"] + 15
#       
#       # Logistic Regression: +10 - Linear log-odds appropriate
#       scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#         scores$Score[scores$Algorithm == "Logistic Regression"] + 10
#       
#       # QDA: +8 - Works but simpler linear methods might suffice
#       scores$Score[scores$Algorithm == "QDA"] <- 
#         scores$Score[scores$Algorithm == "QDA"] + 8
#       
#       # Non-linear methods: minor boost (they still work)
#       scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#         scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#       scores$Score[scores$Algorithm == "Random Forest"] <- 
#         scores$Score[scores$Algorithm == "Random Forest"] + 5
#       scores$Score[scores$Algorithm == "XGBoost"] <- 
#         scores$Score[scores$Algorithm == "XGBoost"] + 5
#       scores$Score[scores$Algorithm == "k-NN"] <- 
#         scores$Score[scores$Algorithm == "k-NN"] + 5
#       scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#         scores$Score[scores$Algorithm == "Naive Bayes"] + 5
#       
#     } else if (linear_separability > 0.5) {
#       assumptions$class_separability <- "MODERATE: Partial linear separability"
#       
#       # Non-linear methods have advantage
#       scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#         scores$Score[scores$Algorithm == "SVM (RBF)"] + 10
#       scores$Score[scores$Algorithm == "Random Forest"] <- 
#         scores$Score[scores$Algorithm == "Random Forest"] + 10
#       scores$Score[scores$Algorithm == "XGBoost"] <- 
#         scores$Score[scores$Algorithm == "XGBoost"] + 10
#       scores$Score[scores$Algorithm == "k-NN"] <- 
#         scores$Score[scores$Algorithm == "k-NN"] + 8
#       
#       # QDA can capture some curvature
#       scores$Score[scores$Algorithm == "QDA"] <- 
#         scores$Score[scores$Algorithm == "QDA"] + 5
#       
#       # Linear methods: neutral
#       scores$Score[scores$Algorithm == "LDA"] <- 
#         scores$Score[scores$Algorithm == "LDA"] + 0
#       scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#         scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
#       scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#         scores$Score[scores$Algorithm == "Logistic Regression"] + 0
#       scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#         scores$Score[scores$Algorithm == "Naive Bayes"] + 3
#       
#     } else {
#       assumptions$class_separability <- "POOR: Classes likely not linearly separable"
#       
#       # Non-linear methods strongly preferred
#       # SVM RBF: +15 - RBF kernel can model any smooth boundary
#       scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#         scores$Score[scores$Algorithm == "SVM (RBF)"] + 15
#       
#       # RF/XGBoost: +15 - Tree ensembles approximate any boundary
#       scores$Score[scores$Algorithm == "Random Forest"] <- 
#         scores$Score[scores$Algorithm == "Random Forest"] + 15
#       scores$Score[scores$Algorithm == "XGBoost"] <- 
#         scores$Score[scores$Algorithm == "XGBoost"] + 15
#       
#       # k-NN: +10 - Local density-based, inherently non-linear
#       scores$Score[scores$Algorithm == "k-NN"] <- 
#         scores$Score[scores$Algorithm == "k-NN"] + 10
#       
#       # QDA: +5 - Can capture quadratic boundaries
#       scores$Score[scores$Algorithm == "QDA"] <- 
#         scores$Score[scores$Algorithm == "QDA"] + 5
#       
#       # Naive Bayes: neutral
#       scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#         scores$Score[scores$Algorithm == "Naive Bayes"] + 0
#       
#       # Linear methods penalized
#       # LDA: -15 - Cannot capture non-linear boundaries
#       scores$Score[scores$Algorithm == "LDA"] <- 
#         scores$Score[scores$Algorithm == "LDA"] - 15
#       
#       # SVM Linear: -10 - Linear margin insufficient
#       scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#         scores$Score[scores$Algorithm == "SVM (Linear)"] - 10
#       
#       # Logistic Regression: -10 - Linear log-odds inappropriate
#       scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#         scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     }
#   } else {
#     assumptions$class_separability <- "UNKNOWN: Could not compute separability"
#   }
#   
#   if (verbose) {
#     cat(sprintf("  Linear separability ratio: %.3f\n", linear_separability))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$class_separability))
#   }
#   
#   # ============================================================================
#   # TEST 7: CLASS IMBALANCE
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Class imbalance affects algorithms differently based on their learning
#   # mechanism and default optimization criteria.
#   #
#   # Imbalance Ratio = max(class_count) / min(class_count)
#   # - Balanced: ratio < 2
#   # - Moderate: ratio 2-5
#   # - Severe: ratio > 5
#   #
#   # - k-NN (-15 severe): Majority class dominates local neighborhoods.
#   #   Minority class samples are "surrounded" by majority class samples.
#   #   Voting always favors majority class unless k is very small.
#   #
#   # - LDA (-15 severe): Prior probabilities heavily weighted to majority class.
#   #   P(Y=k) estimates become dominated by majority class.
#   #
#   # - QDA (-15 severe): Same prior probability issue as LDA.
#   #
#   # - Logistic Regression (-10 severe): MLE optimizes overall likelihood,
#   #   which is dominated by majority class. Needs class weighting.
#   #
#   # - SVM Linear (-10 severe): Margin optimization biased toward majority.
#   #   Support vectors may not include minority class examples.
#   #
#   # - Naive Bayes (-10 severe): Prior probability P(Y=k) issue similar to LDA.
#   #
#   # - SVM RBF (-8 severe): Slightly more robust due to kernel flexibility,
#   #   but still affected by margin optimization bias.
#   #
#   # - RF/XGBoost (+5 severe): Can use class_weight parameter. Tree splits
#   #   can still find minority class if features are discriminative.
#   #   Stratified sampling in RF helps maintain class representation.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 7: Class Imbalance Assessment\n", strrep("-", 40), "\n")
#   
#   class_proportions <- class_counts / n_samples
#   imbalance_ratio <- max(class_counts) / min(class_counts)
#   
#   # Normalized entropy: 1 = perfect balance, 0 = extreme imbalance
#   entropy <- -sum(class_proportions * log(class_proportions + 1e-10)) / log(n_classes)
#   
#   diagnostics$class_imbalance <- list(
#     class_counts = class_counts,
#     class_proportions = class_proportions,
#     imbalance_ratio = imbalance_ratio,
#     normalized_entropy = entropy
#   )
#   
#   if (imbalance_ratio > 5) {
#     assumptions$class_imbalance <- "SEVERE: Highly imbalanced classes"
#     
#     # k-NN: -15 - Majority class dominates neighborhoods
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] - 15
#     
#     # LDA/QDA: -15 - Prior probability estimation biased
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] - 15
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] - 15
#     
#     # Logistic Regression: -10 - MLE dominated by majority class
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     
#     # Naive Bayes: -10 - Prior probability issue
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#     
#     # SVM Linear: -10 - Margin biased toward majority
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 10
#     
#     # SVM RBF: -8 - Slightly more robust
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] - 8
#     
#     # RF/XGBoost: +5 - Can use class weights, stratified sampling
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     
#   } else if (imbalance_ratio > 2) {
#     assumptions$class_imbalance <- "MODERATE: Some class imbalance"
#     
#     # Reduced penalties
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] - 8
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] - 5
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] - 5
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 5
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] - 3
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 3
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 3
#     
#   } else {
#     assumptions$class_imbalance <- "OK: Classes reasonably balanced"
#     
#     # All methods benefit from balanced classes
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] + 5
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] + 5
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 5
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] + 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 5
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#     # RF/XGBoost: smaller boost (they're robust either way)
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 3
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 3
#   }
#   
#   if (verbose) {
#     cat(sprintf("  Imbalance ratio (max/min): %.2f\n", imbalance_ratio))
#     cat(sprintf("  Normalized entropy: %.3f (1=balanced)\n", entropy))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$class_imbalance))
#   }
#   
#   # ============================================================================
#   # TEST 8: OUTLIERS AND INFLUENTIAL POINTS
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Outliers affect algorithms based on their estimation procedures:
#   #
#   # - LDA (-15 severe): Mean and covariance estimation are NOT robust.
#   #   A single outlier can pull the class mean significantly.
#   #   Covariance matrix is even more sensitive to outliers.
#   #
#   # - QDA (-15 severe): Same sensitivity as LDA, multiplied by having
#   #   separate covariance matrices per class.
#   #
#   # - k-NN (-15 severe): Distance-based. Outliers either become isolated
#   #   (misclassified) or can pull nearby points into wrong class if k is small.
#   #
#   # - Logistic Regression (-10 severe): MLE can be influenced by high-leverage
#   #   points. Outliers in feature space with unusual class labels are problematic.
#   #
#   # - Naive Bayes (-10 severe): Mean/variance estimates affected by outliers.
#   #
#   # - SVM Linear (-8 severe): Soft margin provides some robustness (slack variables),
#   #   but outliers as support vectors can skew the margin.
#   #
#   # - SVM RBF (+5 severe): RBF kernel effectively isolates outliers since they
#   #   are far from other points in kernel space. More robust than linear.
#   #
#   # - Random Forest (+8 severe): Trees are robust due to:
#   #   1) Bagging: outliers only affect ~63% of trees
#   #   2) Averaging: outlier influence is diluted
#   #   3) Splits: outliers often end up in their own leaves
#   #
#   # - XGBoost (+5 severe): Robust but slightly less than RF since errors compound.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 8: Outlier Detection (Mahalanobis Distance)\n", strrep("-", 40), "\n")
#   
#   # Mahalanobis distance for outlier detection
#   outlier_result <- tryCatch({
#     # Use robust covariance estimate if possible
#     center <- colMeans(feature_data, na.rm = TRUE)
#     cov_mat <- cov(feature_data, use = "pairwise.complete.obs")
#     
#     # Regularize if singular
#     if (det(cov_mat) < 1e-10) {
#       cov_mat <- cov_mat + diag(1e-6, ncol(cov_mat))
#     }
#     
#     md <- mahalanobis(feature_data, center, cov_mat)
#     
#     # Chi-square cutoff for outliers (97.5th percentile)
#     cutoff <- qchisq(0.975, df = n_features)
#     n_outliers <- sum(md > cutoff, na.rm = TRUE)
#     pct_outliers <- n_outliers / n_samples * 100
#     
#     list(
#       mahalanobis_distances = md,
#       cutoff = cutoff,
#       n_outliers = n_outliers,
#       pct_outliers = pct_outliers
#     )
#   }, error = function(e) {
#     list(n_outliers = NA, pct_outliers = NA, error = as.character(e))
#   })
#   
#   diagnostics$outliers <- outlier_result
#   
#   if (!is.na(outlier_result$pct_outliers)) {
#     if (outlier_result$pct_outliers > 10) {
#       assumptions$outliers <- "SEVERE: Many potential outliers detected"
#       
#       # LDA/QDA: -15 - Mean and covariance estimation highly sensitive
#       scores$Score[scores$Algorithm == "LDA"] <- 
#         scores$Score[scores$Algorithm == "LDA"] - 15
#       scores$Score[scores$Algorithm == "QDA"] <- 
#         scores$Score[scores$Algorithm == "QDA"] - 15
#       
#       # k-NN: -15 - Distance calculations distorted by outliers
#       scores$Score[scores$Algorithm == "k-NN"] <- 
#         scores$Score[scores$Algorithm == "k-NN"] - 15
#       
#       # Logistic Regression: -10 - MLE influenced by high-leverage points
#       scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#         scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#       
#       # Naive Bayes: -10 - Mean/variance estimates affected
#       scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#         scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#       
#       # SVM Linear: -8 - Soft margin helps but not completely
#       scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#         scores$Score[scores$Algorithm == "SVM (Linear)"] - 8
#       
#       # Random Forest: +8 - Robust due to bagging and averaging
#       scores$Score[scores$Algorithm == "Random Forest"] <- 
#         scores$Score[scores$Algorithm == "Random Forest"] + 8
#       
#       # XGBoost: +5 - Robust but boosting can amplify outlier errors
#       scores$Score[scores$Algorithm == "XGBoost"] <- 
#         scores$Score[scores$Algorithm == "XGBoost"] + 5
#       
#       # SVM RBF: +5 - Kernel isolates outliers in feature space
#       scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#         scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
#       
#     } else if (outlier_result$pct_outliers > 5) {
#       assumptions$outliers <- "MODERATE: Some outliers present"
#       
#       # Reduced penalties/boosts
#       scores$Score[scores$Algorithm == "LDA"] <- 
#         scores$Score[scores$Algorithm == "LDA"] - 5
#       scores$Score[scores$Algorithm == "QDA"] <- 
#         scores$Score[scores$Algorithm == "QDA"] - 5
#       scores$Score[scores$Algorithm == "k-NN"] <- 
#         scores$Score[scores$Algorithm == "k-NN"] - 8
#       scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#         scores$Score[scores$Algorithm == "Logistic Regression"] - 5
#       scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#         scores$Score[scores$Algorithm == "Naive Bayes"] - 5
#       scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#         scores$Score[scores$Algorithm == "SVM (Linear)"] - 3
#       scores$Score[scores$Algorithm == "Random Forest"] <- 
#         scores$Score[scores$Algorithm == "Random Forest"] + 5
#       scores$Score[scores$Algorithm == "XGBoost"] <- 
#         scores$Score[scores$Algorithm == "XGBoost"] + 3
#       scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#         scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#       
#     } else {
#       assumptions$outliers <- "OK: Few outliers detected"
#       
#       # Clean data benefits all methods, especially those sensitive to outliers
#       scores$Score[scores$Algorithm == "k-NN"] <- 
#         scores$Score[scores$Algorithm == "k-NN"] + 5
#       scores$Score[scores$Algorithm == "LDA"] <- 
#         scores$Score[scores$Algorithm == "LDA"] + 3
#       scores$Score[scores$Algorithm == "QDA"] <- 
#         scores$Score[scores$Algorithm == "QDA"] + 3
#       scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#         scores$Score[scores$Algorithm == "Logistic Regression"] + 3
#       scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#         scores$Score[scores$Algorithm == "Naive Bayes"] + 3
#       scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#         scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#       scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#         scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#       # RF/XGBoost: neutral (robust either way)
#       scores$Score[scores$Algorithm == "Random Forest"] <- 
#         scores$Score[scores$Algorithm == "Random Forest"] + 0
#       scores$Score[scores$Algorithm == "XGBoost"] <- 
#         scores$Score[scores$Algorithm == "XGBoost"] + 0
#     }
#   } else {
#     assumptions$outliers <- "UNKNOWN: Could not compute outlier detection"
#   }
#   
#   if (verbose) {
#     if (!is.na(outlier_result$pct_outliers)) {
#       cat(sprintf("  Outliers detected: %d (%.1f%%)\n", 
#                   outlier_result$n_outliers, outlier_result$pct_outliers))
#     }
#     cat(sprintf("  Assessment: %s\n\n", assumptions$outliers))
#   }
#   
#   # ============================================================================
#   # TEST 9: DATA SKEWNESS
#   # ============================================================================
#   #
#   # RATIONALE:
#   # Skewness measures asymmetry of distributions. Many parametric methods assume
#   # symmetric (and often normal) distributions.
#   #
#   # Skewness ranges:
#   # - |skew| < 0.5: Approximately symmetric
#   # - |skew| 0.5-1: Moderately skewed
#   # - |skew| > 1: Highly skewed
#   #
#   # - LDA (-15 severe): Derived assuming normality, which implies symmetry.
#   #   Skewed distributions mean the optimal linear boundary is incorrect.
#   #
#   # - QDA (-15 severe): Same normality assumption as LDA.
#   #
#   # - Naive Bayes (-10 severe): Gaussian assumption violated by skewness.
#   #   Mean and variance don't fully characterize skewed distributions.
#   #
#   # - k-NN (-10 severe): Skewness creates uneven point densities.
#   #   In a skewed distribution, most points cluster on one side, distorting
#   #   the notion of "nearest" neighbors in meaningful ways.
#   #
#   # - Logistic Regression (-8 severe): Numerical optimization may struggle.
#   #   Extreme values from skewed tails can dominate the likelihood.
#   #
#   # - SVM Linear (-5 severe): Less affected, but extreme values can become
#   #   support vectors inappropriately.
#   #
#   # - SVM RBF (+3 severe): Kernel transformation can reduce skewness effect.
#   #
#   # - RF/XGBoost (+10 severe): Trees use RANK ORDER for splits, not actual values.
#   #   Whether data is x, log(x), or x^2, the split point ranking is identical.
#   #   This makes trees completely invariant to monotonic transformations
#   #   including those that would change skewness.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 9: Skewness Analysis\n", strrep("-", 40), "\n")
#   
#   skewness_values <- apply(feature_data, 2, function(x) {
#     moments::skewness(x, na.rm = TRUE)
#   })
#   
#   mean_abs_skewness <- mean(abs(skewness_values), na.rm = TRUE)
#   pct_highly_skewed <- sum(abs(skewness_values) > 1, na.rm = TRUE) / n_features * 100
#   
#   diagnostics$skewness <- list(
#     feature_skewness = skewness_values,
#     mean_abs_skewness = mean_abs_skewness,
#     pct_highly_skewed = pct_highly_skewed
#   )
#   
#   if (mean_abs_skewness > 1 || pct_highly_skewed > 50) {
#     assumptions$skewness <- "SEVERE: Data is highly skewed"
#     
#     # LDA/QDA: -15 - Normality (symmetry) assumption violated
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] - 15
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] - 15
#     
#     # Naive Bayes: -10 - Gaussian assumption violated
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 10
#     
#     # k-NN: -10 - Uneven point densities distort neighborhoods
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] - 10
#     
#     # Logistic Regression: -8 - Numerical issues with extreme values
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 8
#     
#     # SVM Linear: -5 - Extreme values as support vectors
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 5
#     
#     # SVM RBF: +3 - Kernel helps isolate extreme values
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#     
#     # RF/XGBoost: +10 - Completely robust (use rank order only)
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 10
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 10
#     
#   } else if (mean_abs_skewness > 0.5 || pct_highly_skewed > 25) {
#     assumptions$skewness <- "MODERATE: Some features are skewed"
#     
#     # Reduced penalties/boosts
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] - 8
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] - 8
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 5
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] - 5
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 3
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] - 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     
#   } else {
#     assumptions$skewness <- "OK: Data is relatively symmetric"
#     
#     # Parametric methods benefit from symmetric distributions
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] + 5
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] + 5
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 5
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] + 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#     # RF/XGBoost: neutral (robust either way)
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 0
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 0
#   }
#   
#   if (verbose) {
#     cat(sprintf("  Mean absolute skewness: %.3f\n", mean_abs_skewness))
#     cat(sprintf("  Highly skewed features (|skew|>1): %.1f%%\n", pct_highly_skewed))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$skewness))
#   }
#   
#   # ============================================================================
#   # TEST 10: DIMENSIONALITY ASSESSMENT
#   # ============================================================================
#   #
#   # RATIONALE:
#   # High dimensionality (p >> n or p > n) creates unique challenges:
#   #
#   # - QDA (-40 high-D): CANNOT work when p > n/k because class covariance
#   #   matrices become singular. Need at least p samples per class for
#   #   estimable covariance. Most severely affected method.
#   #
#   # - LDA (-30 high-D): Pooled covariance becomes singular when p > n-k.
#   #   Even when estimable, estimates are highly unreliable.
#   #
#   # - Logistic Regression (-25 high-D): More parameters (p) than constraints (n)
#   #   leads to perfect separation and undefined MLE. Regularization required.
#   #
#   # - k-NN (-20 high-D): "Curse of dimensionality" - in high-D, all points
#   #   become equidistant. Distance metrics lose discriminative power.
#   #   The ratio of nearest to farthest neighbor distance → 1 as p → ∞.
#   #
#   # - Naive Bayes (-15 high-D): Estimates p marginal distributions, but
#   #   independence assumption becomes more problematic with more features.
#   #
#   # - RF/XGBoost (-5 high-D): Feature subsampling (√p) helps, but very
#   #   high-D can still cause issues. Each tree sees random subset of features.
#   #
#   # - SVM Linear (+20 high-D): DESIGNED for high-D! Kernel trick works in
#   #   infinite-D space. Maximum margin is well-defined even when p > n.
#   #   This is why SVM excels in text classification (10,000+ features).
#   #
#   # - SVM RBF (+15 high-D): Also handles high-D well, but gamma tuning
#   #   becomes more critical. Slightly less boost than linear.
#   # ============================================================================
#   
#   if (verbose) cat("TEST 10: Dimensionality Assessment\n", strrep("-", 40), "\n")
#   
#   # High dimensionality check
#   is_high_dim <- n_features > n_samples
#   dim_ratio <- n_features / n_samples
#   
#   # Effective dimensionality (based on PCA - how many PCs explain 95% variance)
#   if (!is.null(pca_result)) {
#     var_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
#     effective_dim <- which(var_explained >= 0.95)[1]
#     if (is.na(effective_dim)) effective_dim <- n_features
#   } else {
#     effective_dim <- n_features
#   }
#   
#   diagnostics$dimensionality <- list(
#     n_features = n_features,
#     n_samples = n_samples,
#     dim_ratio = dim_ratio,
#     is_high_dimensional = is_high_dim,
#     effective_dim_95pct = effective_dim
#   )
#   
#   if (is_high_dim) {
#     assumptions$dimensionality <- "HIGH-DIMENSIONAL: More features than samples"
#     
#     # QDA: -40 - Covariance matrices singular, method CANNOT work
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] - 40
#     
#     # LDA: -30 - Pooled covariance singular or near-singular
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] - 30
#     
#     # Logistic Regression: -25 - More parameters than samples, separation issues
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 25
#     
#     # k-NN: -20 - Curse of dimensionality, distances meaningless
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] - 20
#     
#     # Naive Bayes: -15 - Independence assumption strained
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 15
#     
#     # RF/XGBoost: -5 - Feature subsampling helps but still some issues
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] - 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] - 5
#     
#     # SVM Linear: +20 - DESIGNED for high-D, maximum margin works well
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 20
#     
#     # SVM RBF: +15 - Also handles high-D but gamma tuning critical
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 15
#     
#   } else if (dim_ratio > 0.5) {
#     assumptions$dimensionality <- "MODERATE-HIGH: Feature-to-sample ratio >0.5"
#     
#     # Reduced penalties/boosts
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] - 15
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] - 10
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] - 10
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] - 10
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] - 5
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 0
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 0
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 10
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 8
#     
#   } else {
#     assumptions$dimensionality <- "OK: Reasonable dimensionality"
#     
#     # All methods work well with reasonable dimensionality
#     scores$Score[scores$Algorithm == "LDA"] <- 
#       scores$Score[scores$Algorithm == "LDA"] + 5
#     scores$Score[scores$Algorithm == "QDA"] <- 
#       scores$Score[scores$Algorithm == "QDA"] + 5
#     scores$Score[scores$Algorithm == "Logistic Regression"] <- 
#       scores$Score[scores$Algorithm == "Logistic Regression"] + 5
#     scores$Score[scores$Algorithm == "Naive Bayes"] <- 
#       scores$Score[scores$Algorithm == "Naive Bayes"] + 5
#     scores$Score[scores$Algorithm == "k-NN"] <- 
#       scores$Score[scores$Algorithm == "k-NN"] + 5
#     scores$Score[scores$Algorithm == "Random Forest"] <- 
#       scores$Score[scores$Algorithm == "Random Forest"] + 5
#     scores$Score[scores$Algorithm == "XGBoost"] <- 
#       scores$Score[scores$Algorithm == "XGBoost"] + 5
#     scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
#       scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
#     scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
#       scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
#   }
#   
#   if (verbose) {
#     cat(sprintf("  Feature/sample ratio: %.2f\n", dim_ratio))
#     cat(sprintf("  Effective dimensions (95%% var): %d\n", effective_dim))
#     cat(sprintf("  Assessment: %s\n\n", assumptions$dimensionality))
#   }
#   
#   # ============================================================================
#   # FINAL RECOMMENDATIONS
#   # ============================================================================
#   
#   # Ensure no negative scores
#   scores$Score <- pmax(scores$Score, 0)
#   
#   # Rank algorithms
#   scores <- scores %>% arrange(desc(Score))
#   
#   if (verbose) {
#     cat(strrep("=", 70), "\n")
#     cat("ALGORITHM RECOMMENDATIONS (Ranked by Suitability Score)\n")
#     cat(strrep("=", 70), "\n\n")
#     print(scores)
#     cat("\n")
#     
#     cat("SUMMARY OF ASSUMPTION VIOLATIONS:\n")
#     cat(strrep("-", 40), "\n")
#     for (name in names(assumptions)) {
#       status <- assumptions[[name]]
#       if (grepl("VIOLATED|SEVERE|CRITICAL|HIGH-DIM|POOR", status)) {
#         cat(sprintf("  ⚠ %s: %s\n", name, status))
#       } else if (grepl("MODERATE|PARTIAL|WARNING", status)) {
#         cat(sprintf(" ⚡ %s: %s\n", name, status))
#       } else {
#         cat(sprintf(" ✓ %s: %s\n", name, status))
#       }
#     }
#     
#     cat("\n")
#     cat(strrep("=", 70), "\n")
#     cat(sprintf("RECOMMENDED ALGORITHM: %s (Score: %.0f)\n", 
#                 scores$Algorithm[1], scores$Score[1]))
#     cat(strrep("=", 70), "\n")
#   }
#   
#   # ============================================================================
#   # RETURN RESULTS
#   # ============================================================================
#   
#   results <- list(
#     diagnostics = diagnostics,
#     assumptions = assumptions,
#     recommendations = scores,
#     recommended_algorithm = scores$Algorithm[1],
#     recommended_score = scores$Score[1],
#     data_summary = list(
#       n_samples = n_samples,
#       n_features = n_features,
#       n_classes = n_classes,
#       class_distribution = class_counts
#     )
#   )
#   
#   return(results)
# }


#' Helper function: Get algorithm function based on recommendation
#' 
#' @param algorithm_name Name of the recommended algorithm
#' @return A caret-compatible method string

# get_algorithm_method <- function(algorithm_name) {
#   methods <- c(
#     "Random Forest" = "ranger",
#     "SVM (RBF)" = "svmRadial",
#     "SVM (Linear)" = "svmLinear",
#     "LDA" = "lda",
#     "QDA" = "qda",
#     "k-NN" = "knn",
#     "XGBoost" = "xgbTree",
#     "Naive Bayes" = "naive_bayes",
#     "Logistic Regression" = "multinom"
#   )
#   
#   return(methods[algorithm_name])
# }


#' =======================================================================================================================
#' COMPREHENSIVE ANALYSIS: Integreate ML MODEL SELECTION into pipeline for finding best Imputation and normalization combo
#' =======================================================================================================================
#' 
#' This section addresses the key questions regarding:
#' 1. Integration of diagnostic tests into imputation/normalization workflow
#' 2. Alternative evaluation metrics beyond OOB MCC
#' 3. Best practices for multi-class classification with missing/skewed data
#' ============================================================================

# ============================================================================
# SECTION 1: ANALYSIS OF YOUR PROPOSED SOLUTION
# ============================================================================

#' Your Proposed Solution:
#' -----------------------
#' Integrate select_ML_diag_tests() into find_best_impute_normalize() to 
#' recommend the best ML algorithm for each imputation/normalization combo.
#' 
#' VERDICT: This is a GOOD solution, but with some important considerations.
#' 
#' STRENGTHS:
#' ----------
#' 1. Adaptive Algorithm Selection: Different transformations DO change data 
#'    characteristics, making some ML algorithms more suitable than others.
#'    
#' 2. Data-Driven Decisions: Running diagnostics AFTER transformation captures
#'    the actual data structure that the ML model will see.
#'    
#' 3. Handles Edge Cases: Log transforms can linearize relationships, making
#'    LDA more appropriate; Z-score normalization satisfies some parametric
#'    assumptions, etc.
#' 
#' POTENTIAL ISSUES & SOLUTIONS:
#' -----------------------------
#' 
#' Issue 1: Computational Cost
#' ---------------------------
#' Running full diagnostic tests for each of the N×M combinations 
#' (imputation × normalization) is expensive.
#' 
#' Solution: 
#' - Use a "two-stage" approach: 
#'   Stage 1: Quick filter using a consistent algorithm (e.g., Random Forest OOB)
#'   Stage 2: Full diagnostics only on top K candidates
#'   
#' Code example:

#' ============================================================================
#' ENHANCED IMPUTATION & NORMALIZATION SELECTION WITH ADAPTIVE ML EVALUATION
#' ============================================================================
#' 
#' This function evaluates different combinations of data imputation and 
#' normalization techniques, using adaptive ML algorithm selection based on
#' diagnostic tests for each transformed data configuration.
#' 
#' Author: Based on Microplastic-Fingerprinting workflow
#' Date: January 2026
#' ============================================================================

# Note: Packages already loaded at script initialization

#' ============================================================================
#' IMPUTATION METHODS
#' ============================================================================

#' Apply imputation method to data
#' 
#' @param data Data frame with missing values
#' @param method Imputation method: "zero", "min", "mean", "median", "knn", 
#'               "missforest", "mice", "lod" (limit of detection)
#' @param class_col Name of class column to exclude from imputation
#' @param lod_factor Factor for LOD imputation (default 0.5)

apply_imputation <- function(data, method, class_col, lod_factor = 0.5) {
  
  # Separate class and features
  class_data <- data[[class_col]]
  feature_data <- data %>% dplyr::select(-all_of(class_col))
  
  imputed_features <- switch(method,
                             "zero" = {
                               feature_data %>% mutate(across(everything(), ~replace_na(., 0)))
                             },
                             
                             "min" = {
                               # Replace with column minimum
                               feature_data %>% 
                                 mutate(across(everything(), ~replace_na(., min(., na.rm = TRUE))))
                             },
                             
                             "half_min" = {
                               # Replace with half of column minimum (common for LOD)
                               feature_data %>% 
                                 mutate(across(everything(), ~replace_na(., min(., na.rm = TRUE) / 2)))
                             },
                             
                             "mean" = {
                               feature_data %>% 
                                 mutate(across(everything(), ~replace_na(., mean(., na.rm = TRUE))))
                             },
                             
                             "median" = {
                               feature_data %>% 
                                 mutate(across(everything(), ~replace_na(., median(., na.rm = TRUE))))
                             },
                             
                             "knn" = {
                               # k-NN imputation using VIM
                               imputed <- VIM::kNN(feature_data, k = 5, imp_var = FALSE)
                               imputed
                             },
                             
                             "missforest" = {
                               # Random Forest imputation
                               set.seed(123)
                               imputed <- missForest::missForest(as.data.frame(feature_data), 
                                                                 maxiter = 10, ntree = 100, 
                                                                 verbose = FALSE)$ximp
                               imputed
                             },
                             
                             "mice_pmm" = {
                               # Multiple imputation by chained equations (predictive mean matching)
                               set.seed(123)
                               imputed <- mice::mice(feature_data, method = "pmm", m = 1, 
                                                     maxit = 5, printFlag = FALSE)
                               mice::complete(imputed, 1)
                             },
                             
                             "lod" = {
                               # Limit of detection: replace with factor * minimum positive value
                               # Bug fix: Handle columns with no positive values
                               feature_data %>%
                                 mutate(across(everything(), ~{
                                   positive_vals <- .[. > 0 & !is.na(.)]
                                   if (length(positive_vals) == 0) {
                                     # Fallback: use absolute minimum non-NA value or 0
                                     min_val <- min(abs(.[!is.na(.)]))
                                     if (is.infinite(min_val) || is.na(min_val)) min_val <- 1
                                     replace_na(., min_val * lod_factor)
                                   } else {
                                     min_pos <- min(positive_vals)
                                     replace_na(., min_pos * lod_factor)
                                   }
                                 }))
                             },
                             
                             # Default: no imputation (may fail downstream)
                             {
                               warning(paste("Unknown imputation method:", method))
                               feature_data
                             }
  )
  
  # Recombine with class column
  result <- bind_cols(imputed_features, tibble(!!class_col := class_data))
  return(result)
}


#' ============================================================================
#' NORMALIZATION METHODS
#' ============================================================================

#' Apply normalization method to data
#' 
#' @param data Data frame (already imputed)
#' @param method Normalization method: "none", "log", "log2", "log10", 
#'               "sqrt", "minmax", "zscore", "robust", "pareto", "quantile"
#' @param class_col Name of class column to exclude from normalization

apply_normalization <- function(data, method, class_col) {
  
  # Separate class and features
  class_data <- data[[class_col]]
  feature_data <- data %>% dplyr::select(-all_of(class_col))
  
  normalized_features <- switch(method,
                                "none" = {
                                  feature_data
                                },
                                
                                "log" = {
                                  # Natural log (add small constant to avoid log(0))
                                  feature_data %>% 
                                    mutate(across(everything(), ~log(. + 1)))
                                },
                                
                                "log2" = {
                                  feature_data %>% 
                                    mutate(across(everything(), ~log2(. + 1)))
                                },
                                
                                "log10" = {
                                  feature_data %>% 
                                    mutate(across(everything(), ~log10(. + 1)))
                                },
                                
                                "sqrt" = {
                                  feature_data %>% 
                                    mutate(across(everything(), ~sqrt(abs(.))))
                                },
                                
                                "minmax" = {
                                  # Min-max scaling to [0, 1]
                                  feature_data %>% 
                                    mutate(across(everything(), ~{
                                      rng <- range(., na.rm = TRUE)
                                      if (rng[2] - rng[1] == 0) . else (. - rng[1]) / (rng[2] - rng[1])
                                    }))
                                },
                                
                                "zscore" = {
                                  # Z-score standardization
                                  feature_data %>% 
                                    mutate(across(everything(), ~{
                                      s <- sd(., na.rm = TRUE)
                                      if (s == 0) . else (. - mean(., na.rm = TRUE)) / s
                                    }))
                                },
                                
                                "robust" = {
                                  # Robust scaling using median and IQR
                                  feature_data %>% 
                                    mutate(across(everything(), ~{
                                      med <- median(., na.rm = TRUE)
                                      iqr <- IQR(., na.rm = TRUE)
                                      if (iqr == 0) . - med else (. - med) / iqr
                                    }))
                                },
                                
                                "pareto" = {
                                  # Pareto scaling: center, divide by sqrt(sd)
                                  feature_data %>% 
                                    mutate(across(everything(), ~{
                                      s <- sd(., na.rm = TRUE)
                                      if (s == 0) . - mean(., na.rm = TRUE) 
                                      else (. - mean(., na.rm = TRUE)) / sqrt(s)
                                    }))
                                },
                                
                                "quantile" = {
                                  # Quantile normalization (rank-based)
                                  feature_data %>% 
                                    mutate(across(everything(), ~{
                                      ranks <- rank(., na.last = "keep", ties.method = "average")
                                      qnorm((ranks - 0.5) / sum(!is.na(.)))
                                    }))
                                },
                                
                                "box_cox" = {
                                  # Box-Cox transformation (requires positive values)
                                  # Bug fix: Handle constant columns and potential failures
                                  feature_data %>%
                                    mutate(across(everything(), ~{
                                      # Skip constant columns
                                      if (length(unique(.)) <= 1) return(.)
                                      x_pos <- . + abs(min(., na.rm = TRUE)) + 1
                                      tryCatch({
                                        bc <- MASS::boxcox(lm(x_pos ~ 1), plotit = FALSE)
                                        lambda <- bc$x[which.max(bc$y)]
                                        if (abs(lambda) < 0.01) log(x_pos) else (x_pos^lambda - 1) / lambda
                                      }, error = function(e) {
                                        # Fallback to log transform on error
                                        log(x_pos)
                                      })
                                    }))
                                },
                                
                                # Default: no normalization
                                {
                                  warning(paste("Unknown normalization method:", method))
                                  feature_data
                                }
  )
  
  # Recombine with class column
  result <- bind_cols(normalized_features, tibble(!!class_col := class_data))
  return(result)
}

# =============================================================================
# FIXED VERSION - evaluate_combo.R
# All bugs identified and corrected with comments explaining each fix
# =============================================================================

# FIX: Move library calls to the top of the script (not inside functions)
# This is best practice and avoids repeated loading
library(cluster)
library(moments)

#' Multiclass Matthews Correlation Coefficient
#' FIX: This function was called but never defined - adding implementation
multiclass_mcc <- function(conf_mat) {
  # Convert to matrix if table
  conf_mat <- as.matrix(conf_mat)
  n <- sum(conf_mat)
  
  if (n == 0) return(NA)
  
  # Row and column sums
  row_sums <- rowSums(conf_mat)
  col_sums <- colSums(conf_mat)
  
  # Correct predictions (trace)
  correct <- sum(diag(conf_mat))
  
  # MCC formula for multiclass
  numerator <- (correct * n) - sum(row_sums * col_sums)
  
  denominator_left <- sqrt(n^2 - sum(row_sums^2))
  denominator_right <- sqrt(n^2 - sum(col_sums^2))
  
  denominator <- denominator_left * denominator_right
  
  if (denominator == 0) return(0)
  
  return(numerator / denominator)
}


calculate_metrics <- function(predicted, actual, probabilities = NULL) {
  
  # Bug fix: Ensure inputs are factors with consistent levels
  if (!is.factor(actual)) {
    actual <- as.factor(actual)
  }
  if (!is.factor(predicted)) {
    predicted <- factor(predicted, levels = levels(actual))
  }
  
  # FIX: Handle case where both vectors might be empty or all NA
  if (length(actual) == 0 || all(is.na(actual))) {
    return(list(
      accuracy = NA, macro_precision = NA, macro_recall = NA,
      macro_f1 = NA, weighted_f1 = NA, kappa = NA, mcc = NA,
      auc = NA, confusion_matrix = NULL
    ))
  }
  
  # Confusion matrix
  conf_mat <- table(Predicted = predicted, Actual = actual)
  
  # Basic metrics
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  
  # Per-class metrics
  n_classes <- length(levels(actual))
  # FIX: Removed redundant check - after factor conversion, levels always exist
  
  precision_vec <- numeric(n_classes)
  recall_vec <- numeric(n_classes)
  f1_vec <- numeric(n_classes)
  
  # FIX: Use level names for safer indexing in case confusion matrix 
  # doesn't have all classes represented
  class_levels <- levels(actual)
  
  for (i in seq_along(class_levels)) {
    cls <- class_levels[i]
    
    # FIX: Safely access confusion matrix elements using dimnames
    # This handles cases where a class may not appear in predictions
    tp <- if (cls %in% rownames(conf_mat) && cls %in% colnames(conf_mat)) {
      conf_mat[cls, cls]
    } else {
      0
    }
    
    fp <- if (cls %in% rownames(conf_mat)) {
      sum(conf_mat[cls, ]) - tp
    } else {
      0
    }
    
    fn <- if (cls %in% colnames(conf_mat)) {
      sum(conf_mat[, cls]) - tp
    } else {
      0
    }
    
    precision_vec[i] <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall_vec[i] <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1_vec[i] <- ifelse(precision_vec[i] + recall_vec[i] == 0, 0,
                        2 * precision_vec[i] * recall_vec[i] / 
                          (precision_vec[i] + recall_vec[i]))
  }
  
  # Macro-averaged metrics
  macro_precision <- mean(precision_vec)
  macro_recall <- mean(recall_vec)
  macro_f1 <- mean(f1_vec)
  
  # FIX: Ensure class_weights aligns with f1_vec by using factor levels order
  class_counts <- table(factor(actual, levels = class_levels))
  class_weights <- as.numeric(class_counts) / length(actual)
  weighted_f1 <- sum(f1_vec * class_weights)
  
  # Cohen's Kappa
  p_observed <- accuracy
  p_expected <- sum((rowSums(conf_mat) / sum(conf_mat)) * 
                      (colSums(conf_mat) / sum(conf_mat)))
  kappa <- (p_observed - p_expected) / (1 - p_expected)
  
  # MCC (now using the defined function)
  mcc <- multiclass_mcc(conf_mat)
  
  # AUC (if probabilities provided)
  auc <- NA
  # FIX: Check if probabilities is a matrix before calling ncol()
  if (!is.null(probabilities)) {
    if (is.matrix(probabilities) || is.data.frame(probabilities)) {
      if (ncol(probabilities) > 1) {
        auc <- tryCatch({
          pROC::multiclass.roc(actual, probabilities)$auc[1]
        }, error = function(e) NA)
      }
    }
  }
  
  return(list(
    accuracy = accuracy,
    macro_precision = macro_precision,
    macro_recall = macro_recall,
    macro_f1 = macro_f1,
    weighted_f1 = weighted_f1,
    kappa = kappa,
    mcc = mcc,
    auc = auc,
    confusion_matrix = conf_mat
  ))
}


#' ============================================================================
#' B. DATA QUALITY / TRANSFORMATION METRICS
#' ============================================================================

#' 1. CLUSTER SEPARATION QUALITY
calculate_cluster_quality <- function(transformed_data, class_labels) {
  # FIX: Removed library() call from inside function - moved to top of script
  
  # FIX: Handle NA values in class_labels
  valid_idx <- !is.na(class_labels)
  if (sum(valid_idx) < 2) {
    return(list(mean_silhouette = NA, per_class_silhouette = NA))
  }
  
  transformed_data <- transformed_data[valid_idx, , drop = FALSE]
  class_labels <- class_labels[valid_idx]
  
  dist_matrix <- dist(transformed_data)
  sil <- silhouette(as.numeric(as.factor(class_labels)), dist_matrix)
  
  return(list(
    mean_silhouette = mean(sil[, "sil_width"]),
    per_class_silhouette = tapply(sil[, "sil_width"], class_labels, mean)
  ))
}


#' 2. DISCRIMINANT RATIO (Linear Separability)
calculate_discriminant_ratio <- function(transformed_data, class_labels) {
  # FIX: Handle single-observation classes and NA values
  valid_idx <- !is.na(class_labels)
  transformed_data <- transformed_data[valid_idx, , drop = FALSE]
  class_labels <- class_labels[valid_idx]
  
  # Overall mean
  grand_mean <- colMeans(transformed_data)
  
  # Between-class scatter
  classes <- unique(class_labels)
  n_total <- nrow(transformed_data)
  
  between_scatter <- 0
  within_scatter <- 0
  
  for (cls in classes) {
    class_data <- transformed_data[class_labels == cls, , drop = FALSE]
    n_cls <- nrow(class_data)
    
    # FIX: Skip classes with 0 observations
    if (n_cls == 0) next
    
    class_mean <- colMeans(class_data)
    
    # Between-class
    between_scatter <- between_scatter + 
      n_cls * sum((class_mean - grand_mean)^2)
    
    # Within-class
    # FIX: Handle single-observation classes (variance = 0)
    if (n_cls > 1) {
      within_scatter <- within_scatter + 
        sum(apply(class_data, 1, function(x) sum((x - class_mean)^2)))
    }
    # If n_cls == 1, within-class scatter contribution is 0
  }
  
  ratio <- between_scatter / (within_scatter + 1e-10)
  return(ratio)
}


#' 3. CORRELATION STRUCTURE PRESERVATION
calculate_correlation_preservation <- function(original_data, transformed_data) {
  # FIX: Check for dimension compatibility
  # Only compare columns that exist in both datasets
  common_cols <- min(ncol(original_data), ncol(transformed_data))
  
  if (common_cols < 2) {
    return(list(
      frobenius_change = NA,
      correlation_similarity = NA
    ))
  }
  
  # If dimensions differ, we can only compare structure, not exact correspondence
  if (ncol(original_data) != ncol(transformed_data)) {
    warning("Original and transformed data have different number of columns. ",
            "Correlation preservation may not be directly comparable.")
  }
  
  # Handle missing values in original
  original_cor <- cor(original_data[, 1:common_cols, drop = FALSE], 
                      use = "pairwise.complete.obs")
  transformed_cor <- cor(transformed_data[, 1:common_cols, drop = FALSE],
                         use = "pairwise.complete.obs")
  
  # FIX: Handle case where correlation matrices have NAs
  if (all(is.na(original_cor)) || all(is.na(transformed_cor))) {
    return(list(frobenius_change = NA, correlation_similarity = NA))
  }
  
  # Frobenius norm of difference
  cor_change <- sqrt(sum((original_cor - transformed_cor)^2, na.rm = TRUE))
  
  # Correlation between original and transformed correlation matrices
  cor_similarity <- cor(as.vector(original_cor), as.vector(transformed_cor), 
                        use = "complete.obs")
  
  return(list(
    frobenius_change = cor_change,
    correlation_similarity = cor_similarity
  ))
}


#' 4. NORMALITY IMPROVEMENT
calculate_normality_improvement <- function(original_data, transformed_data, alpha = 0.05) {
  
  # Test normality for original (non-missing values)
  original_normal <- sapply(original_data, function(x) {
    x <- na.omit(x)
    if (length(x) < 3 || length(unique(x)) < 3) return(NA)
    if (length(x) > 5000) x <- sample(x, 5000)
    tryCatch(shapiro.test(x)$p.value > alpha, error = function(e) NA)
  })
  
  # Test normality for transformed
  transformed_normal <- sapply(transformed_data, function(x) {
    x <- na.omit(x)  # FIX: Added na.omit for consistency
    if (length(x) < 3 || length(unique(x)) < 3) return(NA)
    if (length(x) > 5000) x <- sample(x, 5000)
    tryCatch(shapiro.test(x)$p.value > alpha, error = function(e) NA)
  })
  
  pct_original <- mean(original_normal, na.rm = TRUE) * 100
  pct_transformed <- mean(transformed_normal, na.rm = TRUE) * 100
  
  return(list(
    original_pct_normal = pct_original,
    transformed_pct_normal = pct_transformed,
    improvement = pct_transformed - pct_original
  ))
}


#' 5. SKEWNESS REDUCTION
calculate_skewness_reduction <- function(original_data, transformed_data) {
  # FIX: Removed library() call from inside function - moved to top of script
  
  original_skew <- sapply(original_data, function(x) skewness(x, na.rm = TRUE))
  # FIX: Use sapply with function wrapper and na.rm = TRUE for consistency
  transformed_skew <- sapply(transformed_data, function(x) skewness(x, na.rm = TRUE))
  
  mean_abs_original <- mean(abs(original_skew), na.rm = TRUE)
  mean_abs_transformed <- mean(abs(transformed_skew), na.rm = TRUE)
  
  # FIX: Handle division by zero when mean_abs_original is 0
  pct_reduction <- if (mean_abs_original == 0) {
    if (mean_abs_transformed == 0) 0 else -Inf
  } else {
    (mean_abs_original - mean_abs_transformed) / mean_abs_original * 100
  }
  
  return(list(
    original_mean_abs_skew = mean_abs_original,
    transformed_mean_abs_skew = mean_abs_transformed,
    skew_reduction = mean_abs_original - mean_abs_transformed,
    pct_reduction = pct_reduction
  ))
}


#' 6. INFORMATION PRESERVATION (for imputation)
calculate_distribution_preservation <- function(original_data, imputed_data, class_col = NULL) {
  
  if (!is.null(class_col)) {
    original_data <- original_data[, !names(original_data) %in% class_col, drop = FALSE]
    imputed_data <- imputed_data[, !names(imputed_data) %in% class_col, drop = FALSE]
  }
  
  ks_pvalues <- sapply(names(original_data), function(col) {
    orig <- na.omit(original_data[[col]])
    imp <- imputed_data[[col]]
    
    if (length(orig) < 2 || length(imp) < 2) return(NA)
    
    tryCatch(
      ks.test(orig, imp)$p.value,
      error = function(e) NA
    )
  })
  
  # Higher p-values = distributions more similar = better imputation
  return(list(
    ks_pvalues = ks_pvalues,
    mean_ks_pvalue = mean(ks_pvalues, na.rm = TRUE),
    pct_preserved = mean(ks_pvalues > 0.05, na.rm = TRUE) * 100
  ))
}

evaluate_combo <- function(data, 
                           class_col,
                           impute_method, 
                           norm_method,
                           ml_diagnostics = NULL,
                           use_adaptive_ml = TRUE,
                           cv_folds = 5,
                           verbose = FALSE) {
  
  set.seed(123)
  
  # Step 1: Apply imputation
  imputed_data <- tryCatch({
    apply_imputation(data, impute_method, class_col)
  }, error = function(e) {
    if (verbose) cat("  Imputation error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(imputed_data)) {
    return(list(
      impute_method = impute_method,
      norm_method = norm_method,
      success = FALSE,
      error = "Imputation failed"
    ))
  }
  
  # Step 2: Apply normalization
  normalized_data <- tryCatch({
    apply_normalization(imputed_data, norm_method, class_col)
  }, error = function(e) {
    if (verbose) cat("  Normalization error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(normalized_data)) {
    return(list(
      impute_method = impute_method,
      norm_method = norm_method,
      success = FALSE,
      error = "Normalization failed"
    ))
  }
  
  # Check for any remaining issues
  feature_data <- normalized_data %>% dplyr::select(-dplyr::all_of(class_col))
  if (any(is.na(feature_data)) || any(is.infinite(as.matrix(feature_data)))) {
    return(list(
      impute_method = impute_method,
      norm_method = norm_method,
      success = FALSE,
      error = "Invalid values after transformation"
    ))
  }
  
  # Step 3: Run ML diagnostics on transformed data (if adaptive)
  if (use_adaptive_ml) {
    diag_results <- tryCatch({
      select_ML_diag_tests(normalized_data, class_col, verbose = FALSE)
    }, error = function(e) {
      if (verbose) cat("  Diagnostic error:", e$message, "\n")
      NULL
    })
    
    recommended_algo <- if (!is.null(diag_results)) {
      diag_results$recommended_algorithm
    } else {
      "Random Forest"  # Default fallback
    }
    
    algo_method <- get_algorithm_method(recommended_algo)
  } else {
    recommended_algo <- "Random Forest"
    algo_method <- "ranger"
    diag_results <- NULL
  }
  
  # Step 4: Evaluate classification metric with cross-validation or OOB MCC
  class_vector <- as.factor(normalized_data[[class_col]])
  
  # Bug fix: Remove NA values from class vector and corresponding features
  na_idx <- is.na(class_vector)
  if (any(na_idx)) {
    class_vector <- class_vector[!na_idx]
    feature_data <- feature_data[!na_idx, , drop = FALSE]
    # Drop unused factor levels after removing NAs
    class_vector <- droplevels(class_vector)
  }
  
  # Determine if we can use full CV or need OOB
  min_class_size <- min(table(class_vector))
  use_oob <- algo_method == "ranger" && min_class_size < cv_folds
  
  if (use_oob && algo_method == "ranger") {
    # Use OOB for Random Forest
    rf_model <- tryCatch({
      ranger::ranger(
        x = feature_data,
        y = class_vector,
        num.trees = 500,
        mtry = floor(sqrt(ncol(feature_data))),
        importance = "impurity",
        oob.error = TRUE,
        probability = TRUE
      )
    }, error = function(e) NULL)
    
    if (is.null(rf_model)) {
      return(list(
        impute_method = impute_method,
        norm_method = norm_method,
        success = FALSE,
        error = "Model training failed"
      ))
    }
    
    # Get OOB predictions
    oob_predictions <- rf_model$predictions
    oob_class <- factor(
      levels(class_vector)[max.col(oob_predictions)],
      levels = levels(class_vector)
    )
    
    metrics <- calculate_metrics(oob_class, class_vector, oob_predictions)
    eval_method <- "OOB"
    
  } else {
    # Use cross-validation
    cv_folds_actual <- min(cv_folds, min_class_size)
    
    train_control <- caret::trainControl(
      method = "cv",
      number = cv_folds_actual,
      classProbs = TRUE,
      savePredictions = "final",
      verboseIter = FALSE
    )
    
    model <- tryCatch({
      caret::train(
        x = feature_data,
        y = class_vector,
        method = algo_method,
        trControl = train_control,
        tuneLength = 3
      )
    }, error = function(e) {
      # Fallback to Random Forest
      if (algo_method != "ranger") {
        caret::train(
          x = feature_data,
          y = class_vector,
          method = "ranger",
          trControl = train_control,
          tuneLength = 3
        )
      } else NULL
    })
    
    if (is.null(model)) {
      return(list(
        impute_method = impute_method,
        norm_method = norm_method,
        success = FALSE,
        error = "Model training failed"
      ))
    }
    
    # Extract CV predictions
    cv_preds <- model$pred
    cv_preds_ordered <- cv_preds[order(cv_preds$rowIndex), ]
    
    predicted <- cv_preds_ordered$pred
    actual <- cv_preds_ordered$obs
    
    # Get probability columns
    prob_cols <- setdiff(names(cv_preds_ordered), 
                         c("pred", "obs", "rowIndex", names(model$bestTune)))
    probabilities <- as.matrix(cv_preds_ordered[, prob_cols, drop = FALSE])
    
    metrics <- calculate_metrics(predicted, actual, probabilities)
    eval_method <- paste0(cv_folds_actual, "-fold CV")
  }
  
  # Step 5: Calculate data quality metrics for the transformation
  original_features <- data %>% dplyr::select(-dplyr::all_of(class_col))
  imputed_features <- imputed_data %>% dplyr::select(-dplyr::all_of(class_col))
  
  # 1. Cluster quality (silhouette score)
  cluster_quality <- tryCatch({
    calculate_cluster_quality(feature_data, class_vector)
  }, error = function(e) list(mean_silhouette = NA, per_class_silhouette = NA))
  
  # 2. Discriminant ratio (between-class / within-class variance)
  discriminant_ratio <- tryCatch({
    calculate_discriminant_ratio(feature_data, class_vector)
  }, error = function(e) NA)
  
  # 3. Correlation structure preservation
  correlation_preservation <- tryCatch({
    calculate_correlation_preservation(original_features, feature_data)
  }, error = function(e) list(frobenius_change = NA, correlation_similarity = NA))
  
  # 4. Normality improvement
  normality_metrics <- tryCatch({
    calculate_normality_improvement(original_features, feature_data)
  }, error = function(e) list(original_pct_normal = NA, transformed_pct_normal = NA, improvement = NA))
  
  # 5. Skewness reduction
  skewness_metrics <- tryCatch({
    calculate_skewness_reduction(original_features, feature_data)
  }, error = function(e) list(original_mean_abs_skew = NA, transformed_mean_abs_skew = NA, 
                              skew_reduction = NA, pct_reduction = NA))
  
  # 6. Distribution preservation (imputation quality via KS test)
  distribution_preservation <- tryCatch({
    calculate_distribution_preservation(original_features, imputed_features)
  }, error = function(e) list(ks_pvalues = NA, mean_ks_pvalue = NA, pct_preserved = NA))
  
  return(list(
    impute_method = impute_method,
    norm_method = norm_method,
    success = TRUE,
    recommended_algorithm = recommended_algo,
    algorithm_used = algo_method,
    evaluation_method = eval_method,
    
    # Classification metrics
    accuracy = metrics$accuracy,
    mcc = metrics$mcc,
    kappa = metrics$kappa,
    macro_precision = metrics$macro_precision,
    macro_recall = metrics$macro_recall,
    balanced_accuracy = metrics$macro_recall,
    macro_f1 = metrics$macro_f1,
    weighted_f1 = metrics$weighted_f1,
    auc = metrics$auc,
    
    # Data quality metrics - Cluster separation
    cluster_mean_silhouette = cluster_quality$mean_silhouette,
    cluster_per_class_silhouette = cluster_quality$per_class_silhouette,
    
    # Data quality metrics - Discriminant ratio
    discriminant_ratio = discriminant_ratio,
    
    # Data quality metrics - Correlation preservation
    correlation_frobenius_change = correlation_preservation$frobenius_change,
    correlation_similarity = correlation_preservation$correlation_similarity,
    
    # Data quality metrics - Normality improvement
    original_pct_normal = normality_metrics$original_pct_normal,
    transformed_pct_normal = normality_metrics$transformed_pct_normal,
    normality_improvement = normality_metrics$improvement,
    
    # Data quality metrics - Skewness reduction
    original_mean_abs_skew = skewness_metrics$original_mean_abs_skew,
    transformed_mean_abs_skew = skewness_metrics$transformed_mean_abs_skew,
    skewness_reduction = skewness_metrics$skew_reduction,
    skewness_pct_reduction = skewness_metrics$pct_reduction,
    
    # Data quality metrics - Distribution preservation (imputation quality)
    imputation_mean_ks_pvalue = distribution_preservation$mean_ks_pvalue,
    imputation_pct_preserved = distribution_preservation$pct_preserved,
    
    # Diagnostics
    ml_diagnostics = diag_results,
    confusion_matrix = metrics$confusion_matrix
  ))
}


#' Calculate Composite Score from evaluate_combo Output
#'
#' @param eval_result The output list from evaluate_combo(), OR separate metrics if using legacy mode
#' @param classification_metrics (Optional) Separate classification metrics list for legacy compatibility
#' @param data_quality_metrics (Optional) Separate data quality metrics list for legacy compatibility
#' @param weights Named list of weights for each metric component
#'
#' @return Numeric composite score between 0 and 1
#'
#' @examples
#' # Direct usage with evaluate_combo output:
#' result <- evaluate_combo(data, "class", "knn", "zscore")
#' score <- calculate_composite_score(result)
#'
#' # Legacy usage with separate metric lists:
#' score <- calculate_composite_score(
#'     classification_metrics = list(mcc = 0.8, balanced_accuracy = 0.85, macro_f1 = 0.82),
#'     data_quality_metrics = list(silhouette = 0.3, discriminant_ratio = 2.5)
#' )

calculate_composite_score <- function(
    eval_result = NULL,
    classification_metrics = NULL,
    data_quality_metrics = NULL,
    weights = list(
      # Classification weights (sum to 0.6)
      mcc = 0.25,
      balanced_accuracy = 0.15,
      macro_f1 = 0.20,

      # Data quality weights (sum to 0.4)
      silhouette = 0.15,
      discriminant_ratio = 0.10,
      normality_improvement = 0.08,
      skew_reduction = 0.07
    )
) {


  # ============================================================================

  # Determine input mode: direct eval_result OR separate metric lists
  # ============================================================================

  if (!is.null(eval_result) && is.list(eval_result)) {
    # Direct mode: extract metrics from evaluate_combo output

    # Check if evaluation was successful
    if (!is.null(eval_result$success) && !eval_result$success) {
      warning("evaluate_combo reported failure: ", eval_result$error)
      return(NA)
    }

    # Map evaluate_combo output to classification_metrics
    classification_metrics <- list(
      accuracy = eval_result$accuracy,
      mcc = eval_result$mcc,
      kappa = eval_result$kappa,
      macro_precision = eval_result$macro_precision,
      macro_recall = eval_result$macro_recall,
      balanced_accuracy = eval_result$balanced_accuracy,
      macro_f1 = eval_result$macro_f1,
      weighted_f1 = eval_result$weighted_f1,
      auc = eval_result$auc
    )

    # Map evaluate_combo output to data_quality_metrics
    data_quality_metrics <- list(
      silhouette = eval_result$cluster_mean_silhouette,
      per_class_silhouette = eval_result$cluster_per_class_silhouette,
      discriminant_ratio = eval_result$discriminant_ratio,
      correlation_frobenius_change = eval_result$correlation_frobenius_change,
      correlation_similarity = eval_result$correlation_similarity,
      original_pct_normal = eval_result$original_pct_normal,
      transformed_pct_normal = eval_result$transformed_pct_normal,
      normality_improvement = eval_result$normality_improvement,
      original_mean_abs_skew = eval_result$original_mean_abs_skew,
      transformed_mean_abs_skew = eval_result$transformed_mean_abs_skew,
      skew_reduction = eval_result$skewness_reduction,
      skew_pct_reduction = eval_result$skewness_pct_reduction,
      imputation_mean_ks_pvalue = eval_result$imputation_mean_ks_pvalue,
      imputation_pct_preserved = eval_result$imputation_pct_preserved
    )
  }

  # Validate that we have metrics to work with
  if (is.null(classification_metrics) || is.null(data_quality_metrics)) {
    stop("Must provide either eval_result OR both classification_metrics and data_quality_metrics")
  }

  # ============================================================================
  # Calculate composite score
  # ============================================================================

  score <- 0
  total_weight_used <- 0

  # ---------------------------------------------------------------------------
  # Classification metrics
  # ---------------------------------------------------------------------------

  # MCC: Matthews Correlation Coefficient [-1, 1] -> [0, 1]
  if (!is.null(classification_metrics$mcc) && !is.na(classification_metrics$mcc)) {
    score <- score + weights$mcc * (classification_metrics$mcc + 1) / 2
    total_weight_used <- total_weight_used + weights$mcc
  }

  # Balanced Accuracy: already in [0, 1]
  if (!is.null(classification_metrics$balanced_accuracy) &&
      !is.na(classification_metrics$balanced_accuracy)) {
    score <- score + weights$balanced_accuracy * classification_metrics$balanced_accuracy
    total_weight_used <- total_weight_used + weights$balanced_accuracy
  }

  # Macro F1: already in [0, 1]
  if (!is.null(classification_metrics$macro_f1) && !is.na(classification_metrics$macro_f1)) {
    score <- score + weights$macro_f1 * classification_metrics$macro_f1
    total_weight_used <- total_weight_used + weights$macro_f1
  }

  # ---------------------------------------------------------------------------
  # Data quality metrics
  # ---------------------------------------------------------------------------

  # Silhouette score: [-1, 1] -> [0, 1]
  if (!is.null(data_quality_metrics$silhouette) && !is.na(data_quality_metrics$silhouette)) {
    score <- score + weights$silhouette * (data_quality_metrics$silhouette + 1) / 2
    total_weight_used <- total_weight_used + weights$silhouette
  }

  # Discriminant ratio: log transform and cap
  if (!is.null(data_quality_metrics$discriminant_ratio) &&
      !is.na(data_quality_metrics$discriminant_ratio)) {
    dr_norm <- min(log1p(data_quality_metrics$discriminant_ratio) / 3, 1)
    score <- score + weights$discriminant_ratio * dr_norm
    total_weight_used <- total_weight_used + weights$discriminant_ratio
  }

  # Normality improvement: percentage difference, normalize to [0, 1]
  # Note: Can be negative if transformation worsened normality
  if (!is.null(data_quality_metrics$normality_improvement) &&
      !is.na(data_quality_metrics$normality_improvement)) {
    # Map from [-100, 100] to [0, 1], where 0% improvement = 0.5
    norm_score <- max(0, min(1, (data_quality_metrics$normality_improvement + 100) / 200))
    score <- score + weights$normality_improvement * norm_score
    total_weight_used <- total_weight_used + weights$normality_improvement
  }

  # Skewness reduction: positive is good
  if (!is.null(data_quality_metrics$skew_reduction) &&
      !is.na(data_quality_metrics$skew_reduction) &&
      !is.infinite(data_quality_metrics$skew_reduction)) {
    # Normalize: assume max reasonable reduction is ~2 units of skewness
    skew_score <- max(0, min(1, data_quality_metrics$skew_reduction / 2))
    score <- score + weights$skew_reduction * skew_score
    total_weight_used <- total_weight_used + weights$skew_reduction
  }

  # ---------------------------------------------------------------------------
  # Normalize score if some metrics were missing
  # ---------------------------------------------------------------------------

  if (total_weight_used > 0 && total_weight_used < 1) {
    # Rescale score to account for missing metrics
    score <- score / total_weight_used
  }

  return(score)
}


#' Wrapper function for easy integration with evaluate_combo
#' 
#' @param ... Arguments passed to evaluate_combo
#' @param weights Custom weights for composite score (optional
#' 
#' @return List containing evaluate_combo results plus composite_score

# evaluate_combo_with_score <- function(..., 
#                                       score_weights = NULL) {
#   
#   # Run evaluate_combo
#   result <- evaluate_combo(...)
#   
#   # Calculate composite score
#   if (!is.null(score_weights)) {
#     result$composite_score <- calculate_composite_score(result, weights = score_weights)
#   } else {
#     result$composite_score <- calculate_composite_score(result)
#   }
#   
#   return(result)
# }

#' ============================================================================
#' MAIN FUNCTION: FIND BEST IMPUTATION & NORMALIZATION
#' ============================================================================

#' Find the best imputation and normalization combination
#' 
#' @param data Data frame with features and class column
#' @param class_col Name of the class column
#' @param imputation_methods Vector of imputation methods to test
#' @param normalization_methods Vector of normalization methods to test
#' @param use_adaptive_ml Use adaptive ML algorithm selection (TRUE/FALSE)
#' @param cv_folds Number of cross-validation folds
#' @param parallel Use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param ranking_metric Primary metric for ranking ("mcc", "accuracy", "macro_f1", etc.)
#' @param verbose Print detailed progress
#' 
#' @return List containing best combination and all results
#' 
#' @examples
#' # results <- find_best_impute_normalize(
#' #   data = my_data,
#' #   class_col = "Product_Type",
#' #   imputation_methods = c("median", "knn", "missforest"),
#' #   normalization_methods = c("none", "log", "zscore"),
#' #   use_adaptive_ml = TRUE
#' # )

find_best_impute_normalize <- function(
    data,
    class_col,
    imputation_methods = c("zero", "half_min", "mean", "median", "knn", "missforest"),
    normalization_methods = c("none", "log", "log10", "sqrt", "minmax", "zscore", "robust", "pareto"),
    use_adaptive_ml = TRUE,
    cv_folds = 5,
    parallel = FALSE,
    n_cores = NULL,
    ranking_metric = "mcc",
    verbose = TRUE
) {
  
  # ============================================================================
  # SETUP
  # ============================================================================
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("IMPUTATION & NORMALIZATION OPTIMIZATION\n")
    cat(strrep("=", 70), "\n\n")
    
    cat("Configuration:\n")
    cat(sprintf("  - Imputation methods: %s\n", paste(imputation_methods, collapse = ", ")))
    cat(sprintf("  - Normalization methods: %s\n", paste(normalization_methods, collapse = ", ")))
    cat(sprintf("  - Adaptive ML: %s\n", use_adaptive_ml))
    cat(sprintf("  - Ranking metric: %s\n", ranking_metric))
    cat(sprintf("  - Total combinations: %d\n", 
                length(imputation_methods) * length(normalization_methods)))
    cat("\n")
  }
  
  # Create all combinations
  combinations <- expand.grid(
    impute = imputation_methods,
    norm = normalization_methods,
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # EVALUATE ALL COMBINATIONS
  # ============================================================================
  
  if (parallel && requireNamespace("doParallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    results_list <- foreach::foreach(
      i = 1:nrow(combinations),
      .packages = c("tidyverse", "caret", "ranger", "missForest", "VIM"),
      .export = c("apply_imputation", "apply_normalization", "evaluate_combo",
                  "select_ML_diag_tests", "get_algorithm_method", "calculate_metrics")
    ) %dopar% {
      evaluate_combo(
        data = data,
        class_col = class_col,
        impute_method = combinations$impute[i],
        norm_method = combinations$norm[i],
        use_adaptive_ml = use_adaptive_ml,
        cv_folds = cv_folds,
        verbose = FALSE
      )
    }
    
    parallel::stopCluster(cl)
    
  } else {
    results_list <- list()
    
    for (i in 1:nrow(combinations)) {
      if (verbose) {
        cat(sprintf("  [%d/%d] Testing %s + %s... ", 
                    i, nrow(combinations),
                    combinations$impute[i], 
                    combinations$norm[i]))
      }
      
      result <- evaluate_combo(
        data = data,
        class_col = class_col,
        impute_method = combinations$impute[i],
        norm_method = combinations$norm[i],
        use_adaptive_ml = use_adaptive_ml,
        cv_folds = cv_folds,
        verbose = FALSE
      )
      
      results_list[[i]] <- result
      
      if (verbose) {
        if (result$success) {
          cat(sprintf("Done. %s=%.3f (Algo: %s)\n", 
                      ranking_metric, 
                      result[[ranking_metric]],
                      result$recommended_algorithm))
        } else {
          cat(sprintf("Failed: %s\n", result$error))
        }
      }
    }
  }
  
  # ============================================================================
  # COMPILE RESULTS
  # ============================================================================
  
  # Filter successful results
  successful_results <- results_list[sapply(results_list, function(x) x$success)]
  
  if (length(successful_results) == 0) {
    warning("No successful combinations found!")
    return(list(
      success = FALSE,
      error = "All combinations failed"
    ))
  }
  
  # Create summary data frame
  summary_df <- data.frame(
    imputation = sapply(successful_results, function(x) x$impute_method),
    normalization = sapply(successful_results, function(x) x$norm_method),
    algorithm = sapply(successful_results, function(x) x$recommended_algorithm),
    eval_method = sapply(successful_results, function(x) x$evaluation_method),
    accuracy = sapply(successful_results, function(x) x$accuracy),
    mcc = sapply(successful_results, function(x) x$mcc),
    kappa = sapply(successful_results, function(x) x$kappa),
    macro_f1 = sapply(successful_results, function(x) x$macro_f1),
    weighted_f1 = sapply(successful_results, function(x) x$weighted_f1),
    auc = sapply(successful_results, function(x) ifelse(is.na(x$auc), NA, x$auc)),
    cluster_quality = sapply(successful_results, function(x) x$cluster_quality),
    correlation_change = sapply(successful_results, function(x) x$correlation_change),
    normality_pct = sapply(successful_results, function(x) x$normality_pct),
    stringsAsFactors = FALSE
  )
  
  # Rank by primary metric
  summary_df <- summary_df %>%
    arrange(desc(!!sym(ranking_metric)))
  
  # Get best combination
  best_idx <- 1
  best_result <- successful_results[[which(
    sapply(successful_results, function(x) x$impute_method) == summary_df$imputation[1] &
      sapply(successful_results, function(x) x$norm_method) == summary_df$normalization[1]
  )[1]]]
  
  # ============================================================================
  # OUTPUT
  # ============================================================================
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("RESULTS SUMMARY\n")
    cat(strrep("=", 70), "\n\n")
    
    cat("Top 5 Combinations (by", ranking_metric, "):\n")
    print(head(summary_df %>% 
                 dplyr::select(imputation, normalization, algorithm, 
                               accuracy, mcc, macro_f1, auc), 5))
    
    cat("\n", strrep("=", 70), "\n")
    cat("BEST COMBINATION:\n")
    cat(strrep("=", 70), "\n")
    cat(sprintf("  Imputation: %s\n", best_result$impute_method))
    cat(sprintf("  Normalization: %s\n", best_result$norm_method))
    cat(sprintf("  Recommended Algorithm: %s\n", best_result$recommended_algorithm))
    cat(sprintf("  Evaluation Method: %s\n", best_result$evaluation_method))
    cat("\n")
    cat("  Performance Metrics:\n")
    cat(sprintf("    - Accuracy: %.4f\n", best_result$accuracy))
    cat(sprintf("    - MCC: %.4f\n", best_result$mcc))
    cat(sprintf("    - Kappa: %.4f\n", best_result$kappa))
    cat(sprintf("    - Macro F1: %.4f\n", best_result$macro_f1))
    if (!is.na(best_result$auc)) {
      cat(sprintf("    - AUC: %.4f\n", best_result$auc))
    }
    cat(strrep("=", 70), "\n")
  }
  
  return(list(
    success = TRUE,
    best_combination = list(
      imputation = best_result$impute_method,
      normalization = best_result$norm_method,
      algorithm = best_result$recommended_algorithm
    ),
    best_metrics = list(
      accuracy = best_result$accuracy,
      mcc = best_result$mcc,
      kappa = best_result$kappa,
      macro_f1 = best_result$macro_f1,
      auc = best_result$auc
    ),
    best_result = best_result,
    all_results = successful_results,
    summary_table = summary_df,
    failed_combinations = results_list[!sapply(results_list, function(x) x$success)]
  ))
}

# =============================================================================
# SECTION 1: COMPOSITE SCORE CALCULATION
# =============================================================================

#' Calculate composite score combining classification and data quality metrics
#' 
#' The composite score balances model performance (60%) with data quality (40%).
#' This ensures we select preprocessing that BOTH performs well AND produces 
#' high-quality transformed data.
#'
#' WEIGHT BREAKDOWN:
#' Classification (60%):
#'   - MCC (0.25): Best overall multiclass metric
#'   - Balanced Accuracy (0.15): Handles class imbalance
#'   - Macro F1 (0.20): Harmonic mean of precision/recall
#' 
#' Data Quality (40%):
#'   - Silhouette (0.15): Cluster separation quality
#'   - Discriminant Ratio (0.10): Between/within class variance
#'   - Normality Improvement (0.08): How much normality improved
#'   - Skewness Reduction (0.07): How much skewness was reduced
#'
#' @param classification_metrics List with mcc, balanced_accuracy, macro_f1
#' @param data_quality_metrics List with silhouette, discriminant_ratio, 
#'        normality_improvement, skew_reduction
#' @param weights Named list of weights (must sum to 1.0)
#' @return Numeric composite score between 0 and 1
calculate_composite_score <- function(
    classification_metrics,
    data_quality_metrics,
    weights = list(
      mcc = 0.25, balanced_accuracy = 0.15, macro_f1 = 0.20,
      silhouette = 0.15, discriminant_ratio = 0.10,
      normality_improvement = 0.08, skew_reduction = 0.07
    )
) {
  
  score <- 0
  total_weight_used <- 0
  
  # ---- Classification metrics ----
  
  # MCC: [-1, 1] -> [0, 1]
  if (!is.null(classification_metrics$mcc) && !is.na(classification_metrics$mcc)) {
    mcc_norm <- (classification_metrics$mcc + 1) / 2
    score <- score + weights$mcc * mcc_norm
    total_weight_used <- total_weight_used + weights$mcc
  }
  
  # Balanced Accuracy: already [0, 1]
  if (!is.null(classification_metrics$balanced_accuracy) && 
      !is.na(classification_metrics$balanced_accuracy)) {
    score <- score + weights$balanced_accuracy * classification_metrics$balanced_accuracy
    total_weight_used <- total_weight_used + weights$balanced_accuracy
  }
  
  # Macro F1: already [0, 1]
  if (!is.null(classification_metrics$macro_f1) && !is.na(classification_metrics$macro_f1)) {
    score <- score + weights$macro_f1 * classification_metrics$macro_f1
    total_weight_used <- total_weight_used + weights$macro_f1
  }
  
  # ---- Data quality metrics ----
  
  # Silhouette: [-1, 1] -> [0, 1]
  if (!is.null(data_quality_metrics$silhouette) && !is.na(data_quality_metrics$silhouette)) {
    sil_norm <- (data_quality_metrics$silhouette + 1) / 2
    score <- score + weights$silhouette * sil_norm
    total_weight_used <- total_weight_used + weights$silhouette
  }
  
  # Discriminant ratio: log transform and cap at 1
  if (!is.null(data_quality_metrics$discriminant_ratio) && 
      !is.na(data_quality_metrics$discriminant_ratio) &&
      data_quality_metrics$discriminant_ratio > 0) {
    dr_norm <- min(log1p(data_quality_metrics$discriminant_ratio) / 3, 1)
    score <- score + weights$discriminant_ratio * dr_norm
    total_weight_used <- total_weight_used + weights$discriminant_ratio
  }
  
  # Normality improvement: [-100, 100] -> [0, 1]
  if (!is.null(data_quality_metrics$normality_improvement) && 
      !is.na(data_quality_metrics$normality_improvement)) {
    norm_score <- max(0, min(1, (data_quality_metrics$normality_improvement + 100) / 200))
    score <- score + weights$normality_improvement * norm_score
    total_weight_used <- total_weight_used + weights$normality_improvement
  }
  
  # Skewness reduction: positive = good
  if (!is.null(data_quality_metrics$skew_reduction) && 
      !is.na(data_quality_metrics$skew_reduction) &&
      !is.infinite(data_quality_metrics$skew_reduction)) {
    skew_score <- max(0, min(1, data_quality_metrics$skew_reduction / 2))
    score <- score + weights$skew_reduction * skew_score
    total_weight_used <- total_weight_used + weights$skew_reduction
  }
  
  # Renormalize if some metrics were missing
  if (total_weight_used > 0 && total_weight_used < 0.99) {
    score <- score / total_weight_used
  }
  
  return(score)
}


# =============================================================================
# SECTION 2: PERMUTATION IMPORTANCE (ALGORITHM-AGNOSTIC)
# =============================================================================

#' Calculate permutation importance for any trained model
#' 
#' Measures feature importance by shuffling each feature and measuring
#' the drop in model performance. Works with ANY model type.
#'
#' @param model Trained caret model
#' @param data Data frame with features (no class column)
#' @param y True class labels (factor)
#' @param metric Performance metric: "mcc", "accuracy", "balanced_accuracy"
#' @param n_permutations Number of permutation repeats for stability
#' @param seed Random seed
#' @return Named vector of importance scores (higher = more important)
calculate_permutation_importance <- function(
    model,
    data,
    y,
    metric = "mcc",
    n_permutations = 10,
    seed = 123
) {
  
  set.seed(seed)
  
  features <- names(data)
  n_features <- length(features)
  
  # Helper for MCC calculation
  calc_mcc <- function(conf_mat) {
    conf_mat <- as.matrix(conf_mat)
    n <- sum(conf_mat)
    if (n == 0) return(NA)
    
    row_sums <- rowSums(conf_mat)
    col_sums <- colSums(conf_mat)
    correct <- sum(diag(conf_mat))
    
    numerator <- (correct * n) - sum(row_sums * col_sums)
    denom_left <- sqrt(n^2 - sum(row_sums^2))
    denom_right <- sqrt(n^2 - sum(col_sums^2))
    denom <- denom_left * denom_right
    
    if (denom == 0) return(0)
    return(numerator / denom)
  }
  
  # Calculate baseline performance
  baseline_preds <- predict(model, data)
  
  baseline_metric <- switch(metric,
                            "mcc" = {
                              conf_mat <- table(Predicted = baseline_preds, Actual = y)
                              calc_mcc(conf_mat)
                            },
                            "accuracy" = mean(baseline_preds == y),
                            "balanced_accuracy" = {
                              conf_mat <- table(Predicted = baseline_preds, Actual = y)
                              recall_per_class <- diag(conf_mat) / colSums(conf_mat)
                              mean(recall_per_class, na.rm = TRUE)
                            },
                            mean(baseline_preds == y)
  )
  
  # Calculate importance for each feature
  importance_scores <- numeric(n_features)
  names(importance_scores) <- features
  
  for (j in seq_along(features)) {
    feature_name <- features[j]
    permuted_metrics <- numeric(n_permutations)
    
    for (perm in 1:n_permutations) {
      data_permuted <- data
      data_permuted[[feature_name]] <- sample(data[[feature_name]])
      
      permuted_preds <- predict(model, data_permuted)
      
      permuted_metrics[perm] <- switch(metric,
                                       "mcc" = {
                                         conf_mat <- table(Predicted = permuted_preds, Actual = y)
                                         calc_mcc(conf_mat)
                                       },
                                       "accuracy" = mean(permuted_preds == y),
                                       "balanced_accuracy" = {
                                         conf_mat <- table(Predicted = permuted_preds, Actual = y)
                                         recall_per_class <- diag(conf_mat) / colSums(conf_mat)
                                         mean(recall_per_class, na.rm = TRUE)
                                       },
                                       mean(permuted_preds == y)
      )
    }
    
    # Importance = drop in performance when feature is shuffled
    importance_scores[j] <- baseline_metric - mean(permuted_metrics, na.rm = TRUE)
  }
  
  return(importance_scores)
}


#' Get feature importance - native or permutation fallback
#' 
#' Tries native importance first (faster), falls back to permutation if unavailable.
#'
#' @param model Trained caret model
#' @param data Feature data
#' @param y Class labels
#' @param n_permutations Permutation repeats if needed
#' @return Named vector of importance scores
get_feature_importance <- function(model, data, y, n_permutations = 10) {
  
  # Try native importance first
  native_importance <- tryCatch({
    
    if (inherits(model, "train")) {
      final_model <- model$finalModel
      
      # Random Forest (ranger)
      if (inherits(final_model, "ranger")) {
        imp <- final_model$variable.importance
        if (!is.null(imp)) return(imp)
      }
      
      # Random Forest (randomForest)
      if (inherits(final_model, "randomForest")) {
        imp <- randomForest::importance(final_model)[, 1]
        if (!is.null(imp)) return(imp)
      }
      
      # XGBoost
      if (inherits(final_model, "xgb.Booster")) {
        imp_df <- xgboost::xgb.importance(model = final_model)
        if (!is.null(imp_df) && nrow(imp_df) > 0) {
          imp_vec <- imp_df$Gain
          names(imp_vec) <- imp_df$Feature
          return(imp_vec)
        }
      }
      
      # Try caret's generic varImp
      var_imp <- caret::varImp(model)
      if (!is.null(var_imp$importance)) {
        imp <- rowMeans(var_imp$importance)
        names(imp) <- rownames(var_imp$importance)
        return(imp)
      }
    }
    
    NULL
    
  }, error = function(e) NULL)
  
  # If native available, return it
  if (!is.null(native_importance)) {
    return(native_importance)
  }
  
  # Fallback to permutation importance
  message("  Using permutation importance (native not available)")
  calculate_permutation_importance(model, data, y, metric = "mcc", 
                                   n_permutations = n_permutations)
}


# =============================================================================
# SECTION 3: ENSEMBLE FEATURE SELECTION WITH STABILITY ANALYSIS
# =============================================================================

#' Ensemble feature selection across top-k pipelines
#' 
#' Aggregates feature importance from multiple pipelines to identify
#' features that are consistently important regardless of preprocessing.
#'
#' @param data Original data frame
#' @param class_col Name of class column
#' @param top_k_pipelines List of top-k pipeline results from find_best_impute_normalize
#' @param top_percent Percentage of features considered "top" per pipeline (0-1)
#' @param stability_threshold Minimum stability for "robust" features (0-1) [CONFIGURABLE]
#' @param mean_rank_percentile Maximum rank percentile for robust features (0-1) [CONFIGURABLE]
#' @param n_permutations Permutation repeats for importance calculation
#' @param verbose Print progress
#' @return List with feature_summary, robust_features, tier classifications
ensemble_feature_selection <- function(
    data,
    class_col,
    top_k_pipelines,
    top_percent = 0.20,
    stability_threshold = 0.60,
    mean_rank_percentile = 0.30,
    n_permutations = 10,
    verbose = TRUE
) {
  
  set.seed(123)
  
  k <- length(top_k_pipelines)
  all_features <- setdiff(names(data), class_col)
  n_features <- length(all_features)
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("ENSEMBLE FEATURE SELECTION WITH STABILITY ANALYSIS\n")
    cat(strrep("=", 70), "\n\n")
    cat("Configuration:\n")
    cat(sprintf("  - Number of pipelines: %d\n", k))
    cat(sprintf("  - Number of features: %d\n", n_features))
    cat(sprintf("  - Top percent per pipeline: %.0f%%\n", top_percent * 100))
    cat(sprintf("  - Stability threshold: %.0f%%\n", stability_threshold * 100))
    cat(sprintf("  - Mean rank percentile: %.0f%%\n", mean_rank_percentile * 100))
    cat("\n")
  }
  
  # Initialize matrices
  importance_matrix <- matrix(NA, nrow = n_features, ncol = k,
                              dimnames = list(all_features, paste0("Pipeline_", 1:k)))
  rank_matrix <- matrix(NA, nrow = n_features, ncol = k)
  selection_matrix <- matrix(0, nrow = n_features, ncol = k)
  
  # ============================================================================
  # Calculate importance for each pipeline
  # ============================================================================
  
  if (verbose) cat("Calculating feature importance for each pipeline...\n")
  
  for (i in 1:k) {
    pipeline <- top_k_pipelines[[i]]
    
    if (verbose) {
      cat(sprintf("  Pipeline %d: %s + %s + %s\n", 
                  i, pipeline$impute_method, pipeline$norm_method, pipeline$best_algorithm))
    }
    
    tryCatch({
      # Get preprocessed data from pipeline
      preprocessed_data <- pipeline$preprocessed_data
      feature_data <- preprocessed_data %>% dplyr::select(-dplyr::all_of(class_col))
      class_vector <- as.factor(preprocessed_data[[class_col]])
      
      # Use stored model
      model <- pipeline$best_model
      
      # Get importance scores
      importance_scores <- get_feature_importance(model, feature_data, class_vector,
                                                  n_permutations)
      
      # Align with all_features
      for (feat in names(importance_scores)) {
        if (feat %in% all_features) {
          importance_matrix[feat, i] <- importance_scores[feat]
        }
      }
      
      # Calculate ranks (1 = most important)
      valid_feats <- !is.na(importance_matrix[, i])
      rank_matrix[valid_feats, i] <- rank(-importance_matrix[valid_feats, i], 
                                          ties.method = "average")
      
      # Determine top features
      n_top <- ceiling(sum(valid_feats) * top_percent)
      sorted_feats <- names(sort(importance_matrix[valid_feats, i], decreasing = TRUE))
      top_features <- sorted_feats[1:min(n_top, length(sorted_feats))]
      selection_matrix[top_features, i] <- 1
      
    }, error = function(e) {
      if (verbose) cat(sprintf("    Error: %s\n", e$message))
    })
  }
  
  # ============================================================================
  # Aggregate results across pipelines
  # ============================================================================
  
  if (verbose) cat("\nAggregating results...\n")
  
  mean_ranks <- rowMeans(rank_matrix, na.rm = TRUE)
  stability_scores <- rowSums(selection_matrix, na.rm = TRUE) / k
  rank_percentile <- mean_ranks / n_features
  
  # Mean normalized importance
  normalized_importance <- apply(importance_matrix, 2, function(x) {
    x_valid <- x[!is.na(x)]
    if (length(x_valid) == 0 || max(abs(x_valid)) == 0) return(rep(NA, length(x)))
    x / max(abs(x_valid), na.rm = TRUE)
  })
  mean_normalized_importance <- rowMeans(normalized_importance, na.rm = TRUE)
  
  # ============================================================================
  # Create feature summary with tier classification
  # ============================================================================
  
  feature_summary <- data.frame(
    feature = all_features,
    mean_rank = mean_ranks,
    rank_percentile = rank_percentile,
    stability_score = stability_scores,
    mean_normalized_importance = mean_normalized_importance,
    times_selected = rowSums(selection_matrix, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Tier classification
  feature_summary$tier <- dplyr::case_when(
    feature_summary$stability_score >= 0.80 ~ "Tier 1: Core Markers",
    feature_summary$stability_score >= 0.60 ~ "Tier 2: Strong Candidates",
    feature_summary$stability_score >= 0.40 ~ "Tier 3: Pipeline-Specific",
    TRUE ~ "Tier 4: Unstable"
  )
  
  # Sort by stability then rank
  feature_summary <- feature_summary %>%
    dplyr::arrange(dplyr::desc(stability_score), mean_rank)
  
  # ============================================================================
  # Identify robust features (meet BOTH criteria)
  # ============================================================================
  
  robust_features <- feature_summary$feature[
    feature_summary$stability_score >= stability_threshold &
      feature_summary$rank_percentile <= mean_rank_percentile
  ]
  
  # Stable features (stability only)
  stable_features <- feature_summary$feature[
    feature_summary$stability_score >= stability_threshold
  ]
  
  # ============================================================================
  # Calculate Nogueira Stability Index
  # ============================================================================
  
  k_bar <- mean(colSums(selection_matrix, na.rm = TRUE))
  selection_freq <- rowSums(selection_matrix, na.rm = TRUE)
  s_squared <- var(selection_freq)
  
  nogueira_index <- if (s_squared > 0) {
    1 - (k_bar * (n_features - k_bar)) / (n_features * s_squared)
  } else { 1 }
  
  # ============================================================================
  # Output
  # ============================================================================
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("FEATURE SELECTION RESULTS\n")
    cat(strrep("=", 70), "\n\n")
    
    cat("Tier Distribution:\n")
    tier_counts <- table(feature_summary$tier)
    for (tier_name in names(tier_counts)) {
      cat(sprintf("  %s: %d features\n", tier_name, tier_counts[tier_name]))
    }
    
    cat(sprintf("\nRobust Features (stability >= %.0f%% AND rank <= %.0f%%): %d\n",
                stability_threshold * 100, mean_rank_percentile * 100, 
                length(robust_features)))
    
    cat(sprintf("Nogueira Stability Index: %.3f ", nogueira_index))
    if (nogueira_index >= 0.7) cat("(Good)\n")
    else if (nogueira_index >= 0.5) cat("(Moderate)\n")
    else cat("(Low)\n")
    
    cat("\nTop 10 Features by Stability:\n")
    print(head(feature_summary[, c("feature", "stability_score", "mean_rank", "tier")], 10))
    
    cat(strrep("=", 70), "\n")
  }
  
  return(list(
    success = TRUE,
    feature_summary = feature_summary,
    robust_features = robust_features,
    stable_features = stable_features,
    importance_matrix = importance_matrix,
    rank_matrix = rank_matrix,
    selection_matrix = selection_matrix,
    nogueira_index = nogueira_index,
    parameters = list(
      top_percent = top_percent,
      stability_threshold = stability_threshold,
      mean_rank_percentile = mean_rank_percentile,
      k_pipelines = k,
      n_features = n_features
    )
  ))
}


# =============================================================================
# SECTION 4: FINAL CLASSIFICATION WITH ROBUST FEATURES
# =============================================================================

#' Train final model with robust features and generate publication outputs
#' 
#' @param data Training data frame
#' @param class_col Name of class column
#' @param robust_features Vector of selected robust feature names
#' @param best_pipeline Best pipeline from find_best_impute_normalize
#' @param test_data Optional separate test data
#' @param cv_folds CV folds if no test_data
#' @param verbose Print progress
#' @return List with final_model, confusion_matrix, metrics, plots
train_final_model <- function(
    data,
    class_col,
    robust_features,
    best_pipeline,
    test_data = NULL,
    cv_folds = 5,
    verbose = TRUE
) {
  
  set.seed(123)
  
  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("FINAL MODEL TRAINING WITH ROBUST FEATURES\n")
    cat(strrep("=", 70), "\n\n")
    cat(sprintf("Robust features: %d\n", length(robust_features)))
    cat(sprintf("Algorithm: %s\n", best_pipeline$best_algorithm))
    cat(sprintf("Imputation: %s\n", best_pipeline$impute_method))
    cat(sprintf("Normalization: %s\n", best_pipeline$norm_method))
    cat("\n")
  }
  
  # Get preprocessed data
  preprocessed_data <- best_pipeline$preprocessed_data
  
  # Select robust features
  available_features <- intersect(robust_features, names(preprocessed_data))
  if (length(available_features) < 3) {
    warning("Too few robust features - using all features")
    available_features <- setdiff(names(preprocessed_data), class_col)
  }
  
  feature_data <- preprocessed_data[, available_features, drop = FALSE]
  class_vector <- as.factor(preprocessed_data[[class_col]])
  
  # Train final model
  min_class <- min(table(class_vector))
  actual_folds <- min(cv_folds, min_class)
  
  train_control <- caret::trainControl(
    method = "cv",
    number = max(2, actual_folds),
    classProbs = TRUE,
    savePredictions = "final",
    verboseIter = FALSE
  )
  
  final_model <- caret::train(
    x = feature_data,
    y = class_vector,
    method = best_pipeline$best_algorithm,
    trControl = train_control,
    tuneLength = 3
  )
  
  # Get CV predictions
  cv_preds <- final_model$pred
  cv_preds_ordered <- cv_preds[order(cv_preds$rowIndex), ]
  
  predicted <- cv_preds_ordered$pred
  actual <- cv_preds_ordered$obs
  
  # Confusion matrix
  conf_mat <- table(Predicted = predicted, Actual = actual)
  
  # Calculate metrics
  prob_cols <- setdiff(names(cv_preds_ordered), 
                       c("pred", "obs", "rowIndex", names(final_model$bestTune)))
  probabilities <- if (length(prob_cols) > 0) {
    as.matrix(cv_preds_ordered[, prob_cols, drop = FALSE])
  } else NULL
  
  metrics <- calculate_metrics(predicted, actual, probabilities)
  
  # Confusion matrix data frame for plotting
  cm_df <- as.data.frame(conf_mat) %>%
    dplyr::group_by(Actual) %>%
    dplyr::mutate(
      Total = sum(Freq),
      Percent = round(Freq / Total * 100, 1)
    ) %>%
    dplyr::ungroup()
  
  if (verbose) {
    cat("CLASSIFICATION RESULTS:\n")
    cat(sprintf("  Accuracy: %.4f\n", metrics$accuracy))
    cat(sprintf("  MCC: %.4f\n", metrics$mcc))
    cat(sprintf("  Balanced Accuracy: %.4f\n", metrics$macro_recall))
    cat(sprintf("  Macro F1: %.4f\n", metrics$macro_f1))
    cat(sprintf("  Kappa: %.4f\n", metrics$kappa))
    
    cat("\nCONFUSION MATRIX:\n")
    print(conf_mat)
    cat(strrep("=", 70), "\n")
  }
  
  return(list(
    success = TRUE,
    final_model = final_model,
    predictions = predicted,
    actuals = actual,
    probabilities = probabilities,
    confusion_matrix = conf_mat,
    confusion_matrix_df = cm_df,
    metrics = metrics,
    features_used = available_features,
    n_features = length(available_features)
  ))
}


# =============================================================================
# SECTION 5: VISUALIZATION FUNCTIONS
# =============================================================================

#' Plot confusion matrix for publication
#' 
#' @param cm_df Confusion matrix data frame
#' @param title Plot title
#' @param mcc MCC value to display
#' @return ggplot object
plot_confusion_matrix <- function(cm_df, title = "Confusion Matrix", mcc = NULL) {
  
  p <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Predicted, y = Actual, fill = Percent)) +
    ggplot2::geom_tile(color = "white", size = 0.5) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Freq, "\n(", Percent, "%)")), 
                       size = 4, color = "black") +
    ggplot2::labs(
      title = title,
      subtitle = if (!is.null(mcc)) sprintf("MCC: %.4f", mcc) else NULL,
      x = "Predicted Class",
      y = "Actual Class",
      fill = "Percentage (%)"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot feature tiers
#' 
#' @param feature_summary Feature summary from ensemble_feature_selection
#' @param top_n Number of features to display
#' @return ggplot object
plot_feature_tiers <- function(feature_summary, top_n = 30) {
  
  plot_data <- head(feature_summary, top_n) %>%
    dplyr::mutate(feature = factor(feature, levels = rev(feature)))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stability_score, y = feature, fill = tier)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c(
        "Tier 1: Core Markers" = "#1a9850",
        "Tier 2: Strong Candidates" = "#91cf60",
        "Tier 3: Pipeline-Specific" = "#fee08b",
        "Tier 4: Unstable" = "#d73027"
      )
    ) +
    ggplot2::geom_vline(xintercept = 0.6, linetype = "dashed", color = "red", size = 0.8) +
    ggplot2::labs(
      title = paste0("Top ", top_n, " Features by Stability Score"),
      subtitle = "Dashed line: 60% stability threshold",
      x = "Stability Score",
      y = "Feature",
      fill = "Tier"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1))
  
  return(p)
}


# =============================================================================
# LOADED MESSAGE
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("HYBRID PIPELINE OPTIMIZATION FUNCTIONS LOADED\n")
cat(strrep("=", 70), "\n")
cat("New Functions:\n")
cat("  - calculate_composite_score(): Composite scoring for pipelines\n")
cat("  - calculate_permutation_importance(): Algorithm-agnostic importance\n")
cat("  - get_feature_importance(): Native or permutation fallback\n")
cat("  - ensemble_feature_selection(): Stability-based feature selection\n")
cat("  - train_final_model(): Final model with robust features\n")
cat("  - plot_confusion_matrix(): Publication-ready confusion matrix\n")
cat("  - plot_feature_tiers(): Feature tier visualization\n")
cat(strrep("=", 70), "\n")
