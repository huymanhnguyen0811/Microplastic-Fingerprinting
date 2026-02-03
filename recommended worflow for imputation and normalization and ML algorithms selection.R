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

# Required packages
required_packages <- c(
  "tidyverse", "caret", "ranger", "e1071", "MASS", "class", "xgboost",
  "missForest", "VIM", "mice", "doParallel", "foreach", "mltools",
  "data.table", "yardstick", "pROC"
)
required_packages <- c(
  "tidyverse", "caret", "e1071", "MASS", "nortest", "MVN", "heplots",
  "corrplot", "psych", "vegan", "randomForest", "xgboost", "ranger",
  "class", "FactoMineR", "factoextra", "car", "moments"
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

select_ML_diag_tests <- function(data, 
                                 class_col, 
                                 alpha = 0.05,
                                 verbose = TRUE,
                                 return_recommendations = TRUE) {
  
  # ============================================================================
  # INITIALIZATION AND DATA VALIDATION
  # ============================================================================
  
  if (verbose) cat("\n", strrep("=", 70), "\n")
  if (verbose) cat("AUTOMATED ML MODEL SELECTION DIAGNOSTIC TESTS\n")
  if (verbose) cat(strrep("=", 70), "\n\n")
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame")
  }
  
  if (!class_col %in% names(data)) {
    stop(paste("Class column '", class_col, "' not found in data"))
  }
  
  # Separate features and class
  class_vector <- as.factor(data[[class_col]])
  feature_data <- data %>% 
    dplyr::select(-all_of(class_col)) %>%
    dplyr::select(where(is.numeric))
  
  n_samples <- nrow(feature_data)
  n_features <- ncol(feature_data)
  n_classes <- length(unique(class_vector))
  class_counts <- table(class_vector)
  
  if (verbose) {
    cat("Data Summary:\n")
    cat(sprintf("  - Samples: %d\n", n_samples))
    cat(sprintf("  - Features: %d\n", n_features))
    cat(sprintf("  - Classes: %d\n", n_classes))
    cat(sprintf("  - Class distribution: %s\n", 
                paste(names(class_counts), "=", class_counts, collapse = ", ")))
    cat(sprintf("  - Sample-to-feature ratio: %.2f\n", n_samples / n_features))
    cat("\n")
  }
  
  # Initialize results storage
  diagnostics <- list()
  assumptions <- list()
  scores <- data.frame(
    Algorithm = c("Random Forest", "SVM (RBF)", "SVM (Linear)", "LDA", 
                  "QDA", "k-NN", "XGBoost", "Naive Bayes", "Logistic Regression"),
    Score = rep(100, 9)  # Start with perfect score, penalize violations
  )
  
  # ============================================================================
  # TEST 1: SAMPLE SIZE ADEQUACY
  # ============================================================================
  #
  # RATIONALE:
  # Sample size affects algorithms differently based on their parameter estimation needs:
  #
  # - QDA (-50 critical, -30 warning): Most affected because it estimates k separate
  #   covariance matrices, each requiring p(p+1)/2 parameters. With k classes and p 
  #   features, needs ~10× more samples than LDA. When n/p < 5, covariance matrices
  #   become singular or poorly conditioned.
  #
  # - LDA (-40 critical, -15 warning): Estimates one pooled covariance matrix
  #   (p(p+1)/2 parameters) plus k class means (k×p parameters). Needs fewer samples
  #   than QDA but still requires n/p > 5-10 for stable covariance estimation.
  #
  # - Logistic Regression (-30 critical, -10 warning): MLE estimation requires 
  #   sufficient "events per variable" (EPV). Rule of thumb: 10-20 samples per 
  #   predictor. Low EPV causes biased coefficients and separation problems.
  #
  # - k-NN (-25 critical, -15 warning): Requires dense neighborhoods for reliable 
  #   voting. With few samples, neighbors may be far away and unrepresentative.
  #   Also affected by curse of dimensionality: distance becomes meaningless
  #   as dimensions increase relative to samples.
  #
  # - Naive Bayes (-20 critical, -10 warning): Estimates only p marginal 
  #   distributions per class, not full covariance. Less demanding than LDA/QDA
  #   but still needs enough samples for reliable probability estimates.
  #
  # - XGBoost (-15 critical, -10 warning): Boosting iteratively fits residuals,
  #   prone to overfitting with small samples. More sensitive than RF because
  #   errors compound across boosting iterations.
  #
  # - Random Forest (-10 critical, -5 warning): Bagging provides some protection
  #   through averaging. Still needs sufficient samples for meaningful tree splits.
  #   Minimum ~5-10 samples per terminal node for stable predictions.
  #
  # - SVM Linear (+10 critical, +5 warning): Maximum margin principle means only
  #   support vectors matter, not all samples. Regularization prevents overfitting.
  #   Actually BENEFITS in low-sample scenarios relative to parametric methods.
  #
  # - SVM RBF (+5 critical, +0 warning): Similar reasoning to Linear SVM, but
  #   RBF has more flexibility (gamma parameter) which can cause overfitting
  #   with very few samples if not tuned properly.
  #
  # Per-class sample size additionally penalizes algorithms that need sufficient
  # samples within each class for reliable class-specific parameter estimation.
  # ============================================================================
  
  if (verbose) cat("TEST 1: Sample Size Adequacy\n", strrep("-", 40), "\n")
  
  min_class_size <- min(class_counts)
  samples_per_feature <- n_samples / n_features
  samples_per_class_per_feature <- min_class_size / n_features
  
  diagnostics$sample_size <- list(
    n_samples = n_samples,
    n_features = n_features,
    n_classes = n_classes,
    min_class_size = min_class_size,
    samples_per_feature = samples_per_feature,
    samples_per_class_per_feature = samples_per_class_per_feature
  )
  
  if (samples_per_feature < 5) {
    assumptions$sample_size <- "CRITICAL: Very low sample-to-feature ratio"
    
    # QDA: -50 (SEVERE) - Must estimate k covariance matrices, each p(p+1)/2 params
    # With n/p < 5, matrices will be singular or extremely ill-conditioned
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 50
    
    # LDA: -40 (SEVERE) - Pooled covariance estimation becomes unstable
    # Linear discriminant functions unreliable with poorly estimated covariance
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 40
    
    # Logistic Regression: -30 (MAJOR) - MLE suffers from separation and bias
    # Coefficient estimates become unstable, SEs inflate dramatically
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 30
    
    # k-NN: -25 (MAJOR) - Sparse neighborhoods, unreliable distance-based voting
    # Curse of dimensionality severe: distances become nearly equal in high-D
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 25
    
    # Naive Bayes: -20 (MODERATE-MAJOR) - Marginal distribution estimates unreliable
    # Still better than LDA/QDA since doesn't estimate full covariance
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 20
    
    # XGBoost: -15 (MODERATE) - Boosting overfits easily with few samples
    # Sequential error correction amplifies noise in small datasets
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] - 15
    
    # Random Forest: -10 (MODERATE) - Bagging helps but still needs samples for splits
    # Each tree sees ~63% of data; with few samples, trees become similar
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] - 10
    
    # SVM Linear: +10 (BOOST) - Margin maximization works well in low-sample regime
    # Only support vectors define boundary; regularization prevents overfitting
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 10
    
    # SVM RBF: +5 (MINOR BOOST) - Still benefits from margin principle
    # Slightly less boost than linear due to gamma parameter tuning risk
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
    
  } else if (samples_per_feature < 10) {
    assumptions$sample_size <- "WARNING: Low sample-to-feature ratio"
    
    # Same ordering but reduced penalties - assumptions partially met
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 30
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 15
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 15
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 10
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] - 10
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] - 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 5
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
    
  } else {
    assumptions$sample_size <- "OK: Adequate sample size"
    
    # With adequate samples, all methods can perform well
    # Small boosts to methods that particularly benefit from more data
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 5
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 3
  }
  
  # Additional penalty for small per-class samples
  # Many algorithms need sufficient samples WITHIN each class
  if (min_class_size < 10) {
    # QDA: -25 (additional) - Cannot estimate class-specific covariance with <10 samples
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 25
    
    # XGBoost: -20 - Boosting for minority class becomes unstable
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] - 20
    
    # Random Forest: -15 - Tree splits for small classes unreliable
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] - 15
    
    # k-NN: -15 - Minority class has too few neighbors
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 15
    
    # Logistic Regression: -15 - Separation likely with sparse classes
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 15
    
    # LDA: -10 - Class mean estimation unreliable
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 10
    
    # SVM RBF: -10 - Few support vectors available per class
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] - 10
    
    # Naive Bayes: -10 - Class-conditional probability estimates unreliable
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 10
    
    # SVM Linear: -5 - Less affected due to simpler boundary
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 5
    
  } else if (min_class_size < 20) {
    # Minor penalties for borderline cases
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 10
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] - 5
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 5
  }
  
  if (verbose) {
    cat(sprintf("  Samples per feature: %.2f (recommend >= 10)\n", samples_per_feature))
    cat(sprintf("  Min class samples per feature: %.2f\n", samples_per_class_per_feature))
    cat(sprintf("  Assessment: %s\n\n", assumptions$sample_size))
  }
  
  # ============================================================================
  # TEST 2: UNIVARIATE NORMALITY (for each feature)
  # ============================================================================
  #
  # RATIONALE:
  # Tests whether individual features follow normal distributions. While 
  # multivariate normality (Test 3) is the true assumption, univariate 
  # non-normality often indicates multivariate non-normality.
  #
  # - QDA (-25 when violated): More parameters to estimate than LDA, making it
  #   more sensitive to distributional assumptions. Non-normal features lead to
  #   biased covariance estimates and unreliable class-specific distributions.
  #
  # - LDA (-20 when violated): Discriminant functions derived assuming normality.
  #   The optimal linear boundary (Fisher's criterion) is only optimal under
  #   normality. Non-normality means the "optimal" boundary isn't actually optimal.
  #
  # - Naive Bayes (-15 when violated): Gaussian Naive Bayes explicitly models
  #   P(Xi|Y) ~ N(μ_k, σ_k). Non-normal features make these probability estimates
  #   incorrect, though the independence assumption often dominates errors.
  #
  # - Logistic Regression (-10 when violated): No normality assumption in theory,
  #   but extreme skewness can cause numerical instability in MLE optimization
  #   and affect coefficient interpretation.
  #
  # - k-NN (+5 when violated): Non-parametric, no distributional assumptions.
  #   Actually slightly BENEFITS because it naturally adapts to local density.
  #
  # - SVM RBF (+5 when violated): Kernel methods are distribution-free.
  #   RBF kernel projects to infinite-dimensional space where any distribution
  #   can be separated.
  #
  # - Random Forest (+10 when violated): Trees use rank-order for splits, not
  #   actual values. Completely invariant to monotonic transformations.
  #   log(x), x^2, or any non-normal transformation doesn't affect performance.
  #
  # - XGBoost (+10 when violated): Same reasoning as Random Forest.
  # ============================================================================
  
  if (verbose) cat("TEST 2: Univariate Normality (Shapiro-Wilk)\n", strrep("-", 40), "\n")
  
  normality_results <- data.frame(
    Feature = character(),
    W_statistic = numeric(),
    p_value = numeric(),
    Normal = logical(),
    stringsAsFactors = FALSE
  )
  
  for (col in names(feature_data)) {
    # Shapiro-Wilk test (limited to 5000 samples)
    test_data <- feature_data[[col]]
    if (length(test_data) > 5000) {
      test_data <- sample(test_data, 5000)
    }
    
    # Skip if too few unique values
    if (length(unique(test_data)) < 3) {
      normality_results <- rbind(normality_results, data.frame(
        Feature = col,
        W_statistic = NA,
        p_value = NA,
        Normal = NA
      ))
      next
    }
    
    sw_test <- tryCatch(
      shapiro.test(test_data),
      error = function(e) list(statistic = NA, p.value = NA)
    )
    
    normality_results <- rbind(normality_results, data.frame(
      Feature = col,
      W_statistic = as.numeric(sw_test$statistic),
      p_value = sw_test$p.value,
      Normal = ifelse(is.na(sw_test$p.value), NA, sw_test$p.value > alpha)
    ))
  }
  
  pct_normal <- sum(normality_results$Normal, na.rm = TRUE) / 
    sum(!is.na(normality_results$Normal)) * 100
  
  diagnostics$univariate_normality <- normality_results
  diagnostics$pct_normal_features <- pct_normal
  
  if (pct_normal < 50) {
    assumptions$univariate_normality <- "VIOLATED: Most features non-normal"
    
    # QDA: -25 - Covariance estimation and class distributions unreliable
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 25
    
    # LDA: -20 - Discriminant functions no longer optimal
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 20
    
    # Naive Bayes: -15 - Gaussian density assumption violated
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 15
    
    # Logistic Regression: -10 - Numerical stability and interpretation affected
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    
    # SVM Linear: 0 - No distributional assumptions, but linear boundary less flexible
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
    
    # k-NN: +5 - Non-parametric, adapts to local density
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
    
    # SVM RBF: +5 - Kernel methods are distribution-free
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
    
    # Random Forest: +10 - Only uses rank order, completely robust
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 10
    
    # XGBoost: +10 - Same as Random Forest
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 10
    
  } else if (pct_normal < 75) {
    assumptions$univariate_normality <- "PARTIAL: Some features non-normal"
    
    # Reduced penalties/boosts for partial violation
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 15
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 10
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 8
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    
  } else {
    assumptions$univariate_normality <- "OK: Most features normally distributed"
    
    # Parametric methods get boosted when their assumptions are met
    # LDA: +10 - Discriminant functions are optimal under normality
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 10
    
    # Naive Bayes: +8 - Gaussian density assumption validated
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 8
    
    # QDA: +8 - Same as LDA but slightly less boost (more parameters = more risk)
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 8
    
    # Logistic Regression: +5 - Slight numerical stability benefit
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    
    # SVM Linear: +3 - Linear boundary often works well with normal data
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    
    # Non-parametric methods: neutral - they work regardless
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 0
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 0
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 0
  }
  
  if (verbose) {
    cat(sprintf("  Features tested: %d\n", nrow(normality_results)))
    cat(sprintf("  Normal features: %.1f%%\n", pct_normal))
    cat(sprintf("  Assessment: %s\n\n", assumptions$univariate_normality))
  }
  
  # ============================================================================
  # TEST 3: MULTIVARIATE NORMALITY (Mardia's test)
  # ============================================================================
  #
  # RATIONALE:
  # Tests whether the joint distribution of all features is multivariate normal.
  # This is THE key assumption for LDA/QDA, not just univariate normality.
  #
  # Mardia's test examines both multivariate skewness and kurtosis:
  # - Skewness: Tests if data are symmetrically distributed
  # - Kurtosis: Tests if tails match MVN expectations
  #
  # - QDA (-20): Estimates class-specific MVN distributions. If true distribution
  #   isn't MVN, the quadratic boundaries will be suboptimal.
  #
  # - LDA (-15): Uses Mahalanobis distance which is only optimal under MVN.
  #   Also, Fisher's linear discriminant is derived assuming MVN.
  #
  # - Naive Bayes (-10): Not directly MVN assumption, but non-MVN often means
  #   complex dependencies that violate independence assumption.
  #
  # - Logistic Regression (-8): Not strictly required, but MVN often means
  #   linear log-odds relationship holds better.
  #
  # - Tree/SVM methods (+5 to +8): Benefit when parametric assumptions fail.
  # ============================================================================
  
  if (verbose) cat("TEST 3: Multivariate Normality (Mardia's Test)\n", strrep("-", 40), "\n")
  
  # Select subset of features if too many (MVN test is computationally expensive)
  mvn_features <- feature_data
  if (ncol(mvn_features) > 20) {
    # Use top 20 features by variance
    variances <- apply(mvn_features, 2, var, na.rm = TRUE)
    top_vars <- names(sort(variances, decreasing = TRUE))[1:20]
    mvn_features <- mvn_features[, top_vars]
  }
  
  # Remove any rows with NA
  mvn_features <- na.omit(mvn_features)
  
  mvn_result <- tryCatch({
    MVN::mvn(mvn_features, mvnTest = "mardia", univariatePlot = "none",
             multivariatePlot = "none")
  }, error = function(e) {
    list(multivariateNormality = data.frame(
      Test = c("Mardia Skewness", "Mardia Kurtosis"),
      Statistic = c(NA, NA),
      `p value` = c(NA, NA),
      Result = c("Could not compute", "Could not compute")
    ))
  })
  
  diagnostics$multivariate_normality <- mvn_result$multivariateNormality
  
  mvn_ok <- tryCatch({
    all(mvn_result$multivariateNormality$Result == "YES")
  }, error = function(e) FALSE)
  
  if (!mvn_ok) {
    assumptions$multivariate_normality <- "VIOLATED: Data not multivariate normal"
    
    # QDA: -20 - MVN is core assumption for class-specific distributions
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 20
    
    # LDA: -15 - Mahalanobis-based classification becomes suboptimal
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 15
    
    # Naive Bayes: -10 - Non-MVN suggests complex feature interactions
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 10
    
    # Logistic Regression: -8 - Linear log-odds may not hold
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 8
    
    # SVM Linear: +3 - Margin principle doesn't require MVN
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    
    # k-NN: +5 - Non-parametric, density-adaptive
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
    
    # SVM RBF: +5 - Kernel methods handle arbitrary distributions
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
    
    # RF/XGBoost: +8 - Trees completely robust to distributional form
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 8
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 8
    
  } else {
    assumptions$multivariate_normality <- "OK: Data is multivariate normal"
    
    # LDA: +15 - This is THE assumption LDA is built on
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 15
    
    # QDA: +10 - Also benefits but LDA might be simpler (Occam's razor)
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 10
    
    # Naive Bayes: +8 - Normality supports the Gaussian version
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 8
    
    # Logistic Regression: +5 - Numerical stability benefits
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    
    # SVM Linear: +5 - Linear boundaries often work well with MVN data
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 5
    
    # Non-parametric methods: neutral
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 0
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 0
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 0
  }
  
  if (verbose) {
    cat("  Mardia's Test Results:\n")
    print(mvn_result$multivariateNormality)
    cat(sprintf("  Assessment: %s\n\n", assumptions$multivariate_normality))
  }
  
  # ============================================================================
  # TEST 4: HOMOGENEITY OF COVARIANCE MATRICES (Box's M test)
  # ============================================================================
  #
  # RATIONALE:
  # Tests whether all classes have the same covariance matrix (Σ₁ = Σ₂ = ... = Σₖ).
  # This is the FUNDAMENTAL difference between LDA and QDA.
  #
  # - LDA (-20 when violated): DIRECTLY assumes equal covariances. Uses pooled
  #   Σ = Σ(nₖ-1)Σₖ / (n-k). When violated, the pooled estimate is a poor
  #   representation of any class, leading to suboptimal boundaries.
  #
  # - QDA (+10 when violated): Designed SPECIFICALLY for unequal covariances!
  #   QDA estimates separate Σₖ per class. When covariances truly differ,
  #   QDA's quadratic boundaries are more appropriate than LDA's linear ones.
  #
  # - Logistic Regression (-10 when violated): Doesn't require equal covariances
  #   mathematically, but heterogeneous covariances can make the linear
  #   decision boundary suboptimal.
  #
  # - Naive Bayes (-5 when violated): Independence assumption more important,
  #   but unequal variances within classes can affect probability estimates.
  #
  # - SVM/Tree/k-NN: No covariance assumptions at all.
  # ============================================================================
  
  if (verbose) cat("TEST 4: Homogeneity of Covariance (Box's M Test)\n", strrep("-", 40), "\n")
  
  # Box's M test - use subset of features if needed
  boxm_features <- feature_data
  if (ncol(boxm_features) > 15) {
    variances <- apply(boxm_features, 2, var, na.rm = TRUE)
    top_vars <- names(sort(variances, decreasing = TRUE))[1:15]
    boxm_features <- boxm_features[, top_vars]
  }
  
  boxm_result <- tryCatch({
    heplots::boxM(boxm_features, class_vector)
  }, error = function(e) {
    list(statistic = NA, p.value = NA, error = as.character(e))
  })
  
  diagnostics$boxm_test <- boxm_result
  
  boxm_ok <- tryCatch({
    !is.na(boxm_result$p.value) && boxm_result$p.value > alpha
  }, error = function(e) NA)
  
  if (is.na(boxm_ok)) {
    assumptions$homogeneity_covariance <- "UNKNOWN: Could not compute Box's M"
    # No score adjustments when test couldn't be computed
    
  } else if (!boxm_ok) {
    assumptions$homogeneity_covariance <- "VIOLATED: Heterogeneous covariance matrices"
    
    # LDA: -20 - Core assumption DIRECTLY violated
    # Pooled covariance misrepresents all classes
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 20
    
    # Logistic Regression: -10 - Linear boundary may be suboptimal
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    
    # Naive Bayes: -5 - Minor effect through variance estimates
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 5
    
    # QDA: +10 (BOOST!) - This is exactly when QDA shines
    # It's designed for this situation
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 10
    
    # Non-parametric methods get minor boost
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 3
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
    
  } else {
    assumptions$homogeneity_covariance <- "OK: Homogeneous covariance matrices"
    
    # LDA: +15 - Core assumption validated, LDA is optimal choice
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 15
    
    # QDA: +5 - Still works, but LDA is simpler (Occam's razor)
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 5
    
    # Logistic Regression: +5 - Linear boundary appropriate
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    
    # Naive Bayes: +5 - Variance estimates reliable
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 5
    
    # SVM Linear: +3 - Linear boundary appropriate
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    
    # Non-parametric: neutral
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 0
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 0
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 0
  }
  
  if (verbose) {
    if (!is.null(boxm_result$p.value)) {
      cat(sprintf("  Box's M statistic: %.2f\n", boxm_result$statistic))
      cat(sprintf("  p-value: %.4f\n", boxm_result$p.value))
    }
    cat(sprintf("  Assessment: %s\n\n", assumptions$homogeneity_covariance))
  }
  
  # ============================================================================
  # TEST 5: MULTICOLLINEARITY (VIF and Condition Index)
  # ============================================================================
  #
  # RATIONALE:
  # Multicollinearity occurs when features are highly correlated, meaning they
  # contain redundant information. This affects algorithms differently:
  #
  # - Logistic Regression (-25 severe, -10 moderate): MOST affected. Causes:
  #   1) Inflated standard errors (coefficients unreliable)
  #   2) Unstable coefficients (small data changes → large coef changes)
  #   3) Near-singular Hessian → numerical optimization failures
  #   The design matrix X'X becomes ill-conditioned.
  #
  # - Naive Bayes (-20 severe): DIRECTLY violates the independence assumption!
  #   If X₁ and X₂ are highly correlated, P(X₁,X₂|Y) ≠ P(X₁|Y)P(X₂|Y).
  #   This is the "naive" assumption, and collinearity destroys it.
  #
  # - QDA (-20 severe): Covariance matrix inversion becomes ill-conditioned.
  #   Each class needs its own inversion → problems multiply.
  #
  # - LDA (-15 severe): Pooled covariance inversion ill-conditioned, but
  #   only one matrix vs. k matrices for QDA, so slightly less severe.
  #
  # - k-NN (-15 severe): Correlated features double-count information in
  #   distance calculations. If X₁ ≈ X₂, they both contribute to distance,
  #   effectively doubling the weight of that information.
  #
  # - SVM Linear (-10 severe): Regularization helps, but feature weights
  #   become unstable with high collinearity.
  #
  # - SVM RBF (-5 severe): RBF kernel less affected due to non-linear
  #   projection, but still some impact.
  #
  # - Random Forest (+10 severe): BENEFITS from collinearity handling.
  #   At each split, only √p features considered. Correlated features
  #   provide backup options; if one is chosen, others are implicitly excluded.
  #
  # - XGBoost (+10 severe): Similar feature selection mechanism.
  # ============================================================================
  
  if (verbose) cat("TEST 5: Multicollinearity (Correlation & VIF)\n", strrep("-", 40), "\n")
  
  # Correlation matrix analysis
  cor_matrix <- cor(feature_data, use = "pairwise.complete.obs")
  cor_matrix[is.na(cor_matrix)] <- 0
  
  # Find highly correlated pairs
  cor_threshold <- 0.9
  high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & 
                            upper.tri(cor_matrix), arr.ind = TRUE)
  
  n_high_cor <- nrow(high_cor_pairs)
  pct_high_cor <- n_high_cor / (n_features * (n_features - 1) / 2) * 100
  
  # Condition index (eigenvalue-based)
  # CI = sqrt(λ_max / λ_min), where λ are eigenvalues of correlation matrix
  # CI > 30 indicates severe multicollinearity
  eigen_vals <- eigen(cor_matrix)$values
  condition_index <- sqrt(max(eigen_vals) / min(eigen_vals[eigen_vals > 1e-10]))
  
  diagnostics$multicollinearity <- list(
    n_high_correlations = n_high_cor,
    pct_high_correlations = pct_high_cor,
    condition_index = condition_index,
    correlation_matrix = cor_matrix
  )
  
  if (condition_index > 30 || pct_high_cor > 20) {
    assumptions$multicollinearity <- "SEVERE: High multicollinearity detected"
    
    # Logistic Regression: -25 - Inflated SEs, unstable coefficients
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 25
    
    # Naive Bayes: -20 - Independence assumption DIRECTLY violated
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 20
    
    # QDA: -20 - Multiple ill-conditioned matrix inversions
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 20
    
    # LDA: -15 - Single ill-conditioned inversion
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 15
    
    # k-NN: -15 - Distance metric double-counts correlated features
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 15
    
    # SVM Linear: -10 - Feature weights unstable
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 10
    
    # SVM RBF: -5 - Kernel projection reduces impact
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] - 5
    
    # RF/XGBoost: +10 - Feature subsampling handles collinearity naturally
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 10
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 10
    
  } else if (condition_index > 15 || pct_high_cor > 10) {
    assumptions$multicollinearity <- "MODERATE: Some multicollinearity present"
    
    # Reduced penalties for moderate collinearity
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 10
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] - 10
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] - 8
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] - 8
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 5
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] - 3
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    
  } else {
    assumptions$multicollinearity <- "OK: Low multicollinearity"
    
    # All methods benefit from clean correlation structure
    # Naive Bayes gets highest boost since independence is satisfied
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 8
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    scores$Score[scores$Algorithm == "LDA"] <- scores$Score[scores$Algorithm == "LDA"] + 5
    scores$Score[scores$Algorithm == "QDA"] <- scores$Score[scores$Algorithm == "QDA"] + 5
    scores$Score[scores$Algorithm == "k-NN"] <- scores$Score[scores$Algorithm == "k-NN"] + 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
    # RF/XGBoost: neutral - they handle collinearity either way
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 0
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 0
  }
  
  if (verbose) {
    cat(sprintf("  High correlation pairs (|r| > %.1f): %d\n", cor_threshold, n_high_cor))
    cat(sprintf("  Condition index: %.2f (>30 is severe)\n", condition_index))
    cat(sprintf("  Assessment: %s\n\n", assumptions$multicollinearity))
  }
  
  # ============================================================================
  # TEST 6: CLASS SEPARABILITY (Linear vs Non-linear)
  # ============================================================================
  #
  # RATIONALE:
  # This test assesses whether classes can be separated by linear boundaries
  # or require non-linear decision surfaces.
  #
  # The metric used is the ratio of between-class variance to within-class
  # variance (Fisher's criterion in PCA space):
  # - High ratio (>1.5): Classes are well-separated, likely linearly separable
  # - Medium ratio (0.5-1.5): Partial overlap, non-linear methods may help
  # - Low ratio (<0.5): Significant overlap, need flexible methods
  #
  # - LDA (+15 when good): LDA IS Fisher's Linear Discriminant! When classes
  #   are linearly separable, LDA is theoretically optimal.
  #
  # - SVM Linear (+15 when good): Maximum margin linear classifier excels.
  #
  # - Logistic Regression (+10 when good): Linear log-odds boundary works well.
  #
  # - LDA (-15 when poor): Linear boundary cannot separate overlapping classes.
  #
  # - SVM RBF (+15 when poor): RBF kernel can model arbitrary boundaries.
  #
  # - RF/XGBoost (+15 when poor): Tree ensembles create step-function
  #   approximations to any decision boundary.
  #
  # - k-NN (+10 when poor): Adapts to local density, naturally non-linear.
  #
  # - QDA: Intermediate - can capture quadratic (curved) boundaries.
  # ============================================================================
  
  if (verbose) cat("TEST 6: Class Separability Analysis\n", strrep("-", 40), "\n")
  
  # PCA-based separability
  pca_result <- tryCatch({
    prcomp(feature_data, scale. = TRUE, center = TRUE)
  }, error = function(e) NULL)
  
  linear_separability <- NA
  
  if (!is.null(pca_result)) {
    # Calculate between-class vs within-class variance in PC space
    pc_scores <- pca_result$x[, 1:min(5, ncol(pca_result$x))]
    
    # Simple separability metric using Fisher's criterion
    class_centroids <- aggregate(pc_scores, by = list(class_vector), mean)
    overall_centroid <- colMeans(pc_scores)
    
    # Between-class variance: weighted sum of squared distances from grand mean
    between_var <- sum(sapply(1:n_classes, function(i) {
      ni <- class_counts[i]
      sum((as.numeric(class_centroids[i, -1]) - overall_centroid)^2) * ni
    })) / n_samples
    
    # Within-class variance: sum of squared distances from class means
    within_var <- sum(sapply(levels(class_vector), function(cls) {
      class_data <- pc_scores[class_vector == cls, , drop = FALSE]
      class_centroid <- colMeans(class_data)
      sum(apply(class_data, 1, function(x) sum((x - class_centroid)^2)))
    })) / n_samples
    
    # Separability ratio (Fisher's criterion)
    linear_separability <- between_var / (within_var + 1e-10)
  }
  
  diagnostics$class_separability <- list(
    linear_separability_ratio = linear_separability,
    pca_variance_explained = if(!is.null(pca_result)) 
      summary(pca_result)$importance[2, 1:min(5, ncol(pca_result$x))] else NA
  )
  
  if (!is.na(linear_separability)) {
    if (linear_separability > 1.5) {
      assumptions$class_separability <- "GOOD: Classes appear linearly separable"
      
      # Linear methods excel
      # LDA: +15 - This IS Fisher's discriminant, optimal for linear separation
      scores$Score[scores$Algorithm == "LDA"] <- 
        scores$Score[scores$Algorithm == "LDA"] + 15
      
      # SVM Linear: +15 - Maximum margin linear classifier
      scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
        scores$Score[scores$Algorithm == "SVM (Linear)"] + 15
      
      # Logistic Regression: +10 - Linear log-odds appropriate
      scores$Score[scores$Algorithm == "Logistic Regression"] <- 
        scores$Score[scores$Algorithm == "Logistic Regression"] + 10
      
      # QDA: +8 - Works but simpler linear methods might suffice
      scores$Score[scores$Algorithm == "QDA"] <- 
        scores$Score[scores$Algorithm == "QDA"] + 8
      
      # Non-linear methods: minor boost (they still work)
      scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
        scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
      scores$Score[scores$Algorithm == "Random Forest"] <- 
        scores$Score[scores$Algorithm == "Random Forest"] + 5
      scores$Score[scores$Algorithm == "XGBoost"] <- 
        scores$Score[scores$Algorithm == "XGBoost"] + 5
      scores$Score[scores$Algorithm == "k-NN"] <- 
        scores$Score[scores$Algorithm == "k-NN"] + 5
      scores$Score[scores$Algorithm == "Naive Bayes"] <- 
        scores$Score[scores$Algorithm == "Naive Bayes"] + 5
      
    } else if (linear_separability > 0.5) {
      assumptions$class_separability <- "MODERATE: Partial linear separability"
      
      # Non-linear methods have advantage
      scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
        scores$Score[scores$Algorithm == "SVM (RBF)"] + 10
      scores$Score[scores$Algorithm == "Random Forest"] <- 
        scores$Score[scores$Algorithm == "Random Forest"] + 10
      scores$Score[scores$Algorithm == "XGBoost"] <- 
        scores$Score[scores$Algorithm == "XGBoost"] + 10
      scores$Score[scores$Algorithm == "k-NN"] <- 
        scores$Score[scores$Algorithm == "k-NN"] + 8
      
      # QDA can capture some curvature
      scores$Score[scores$Algorithm == "QDA"] <- 
        scores$Score[scores$Algorithm == "QDA"] + 5
      
      # Linear methods: neutral
      scores$Score[scores$Algorithm == "LDA"] <- 
        scores$Score[scores$Algorithm == "LDA"] + 0
      scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
        scores$Score[scores$Algorithm == "SVM (Linear)"] + 0
      scores$Score[scores$Algorithm == "Logistic Regression"] <- 
        scores$Score[scores$Algorithm == "Logistic Regression"] + 0
      scores$Score[scores$Algorithm == "Naive Bayes"] <- 
        scores$Score[scores$Algorithm == "Naive Bayes"] + 3
      
    } else {
      assumptions$class_separability <- "POOR: Classes likely not linearly separable"
      
      # Non-linear methods strongly preferred
      # SVM RBF: +15 - RBF kernel can model any smooth boundary
      scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
        scores$Score[scores$Algorithm == "SVM (RBF)"] + 15
      
      # RF/XGBoost: +15 - Tree ensembles approximate any boundary
      scores$Score[scores$Algorithm == "Random Forest"] <- 
        scores$Score[scores$Algorithm == "Random Forest"] + 15
      scores$Score[scores$Algorithm == "XGBoost"] <- 
        scores$Score[scores$Algorithm == "XGBoost"] + 15
      
      # k-NN: +10 - Local density-based, inherently non-linear
      scores$Score[scores$Algorithm == "k-NN"] <- 
        scores$Score[scores$Algorithm == "k-NN"] + 10
      
      # QDA: +5 - Can capture quadratic boundaries
      scores$Score[scores$Algorithm == "QDA"] <- 
        scores$Score[scores$Algorithm == "QDA"] + 5
      
      # Naive Bayes: neutral
      scores$Score[scores$Algorithm == "Naive Bayes"] <- 
        scores$Score[scores$Algorithm == "Naive Bayes"] + 0
      
      # Linear methods penalized
      # LDA: -15 - Cannot capture non-linear boundaries
      scores$Score[scores$Algorithm == "LDA"] <- 
        scores$Score[scores$Algorithm == "LDA"] - 15
      
      # SVM Linear: -10 - Linear margin insufficient
      scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
        scores$Score[scores$Algorithm == "SVM (Linear)"] - 10
      
      # Logistic Regression: -10 - Linear log-odds inappropriate
      scores$Score[scores$Algorithm == "Logistic Regression"] <- 
        scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    }
  } else {
    assumptions$class_separability <- "UNKNOWN: Could not compute separability"
  }
  
  if (verbose) {
    cat(sprintf("  Linear separability ratio: %.3f\n", linear_separability))
    cat(sprintf("  Assessment: %s\n\n", assumptions$class_separability))
  }
  
  # ============================================================================
  # TEST 7: CLASS IMBALANCE
  # ============================================================================
  #
  # RATIONALE:
  # Class imbalance affects algorithms differently based on their learning
  # mechanism and default optimization criteria.
  #
  # Imbalance Ratio = max(class_count) / min(class_count)
  # - Balanced: ratio < 2
  # - Moderate: ratio 2-5
  # - Severe: ratio > 5
  #
  # - k-NN (-15 severe): Majority class dominates local neighborhoods.
  #   Minority class samples are "surrounded" by majority class samples.
  #   Voting always favors majority class unless k is very small.
  #
  # - LDA (-15 severe): Prior probabilities heavily weighted to majority class.
  #   P(Y=k) estimates become dominated by majority class.
  #
  # - QDA (-15 severe): Same prior probability issue as LDA.
  #
  # - Logistic Regression (-10 severe): MLE optimizes overall likelihood,
  #   which is dominated by majority class. Needs class weighting.
  #
  # - SVM Linear (-10 severe): Margin optimization biased toward majority.
  #   Support vectors may not include minority class examples.
  #
  # - Naive Bayes (-10 severe): Prior probability P(Y=k) issue similar to LDA.
  #
  # - SVM RBF (-8 severe): Slightly more robust due to kernel flexibility,
  #   but still affected by margin optimization bias.
  #
  # - RF/XGBoost (+5 severe): Can use class_weight parameter. Tree splits
  #   can still find minority class if features are discriminative.
  #   Stratified sampling in RF helps maintain class representation.
  # ============================================================================
  
  if (verbose) cat("TEST 7: Class Imbalance Assessment\n", strrep("-", 40), "\n")
  
  class_proportions <- class_counts / n_samples
  imbalance_ratio <- max(class_counts) / min(class_counts)
  
  # Normalized entropy: 1 = perfect balance, 0 = extreme imbalance
  entropy <- -sum(class_proportions * log(class_proportions + 1e-10)) / log(n_classes)
  
  diagnostics$class_imbalance <- list(
    class_counts = class_counts,
    class_proportions = class_proportions,
    imbalance_ratio = imbalance_ratio,
    normalized_entropy = entropy
  )
  
  if (imbalance_ratio > 5) {
    assumptions$class_imbalance <- "SEVERE: Highly imbalanced classes"
    
    # k-NN: -15 - Majority class dominates neighborhoods
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] - 15
    
    # LDA/QDA: -15 - Prior probability estimation biased
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] - 15
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] - 15
    
    # Logistic Regression: -10 - MLE dominated by majority class
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    
    # Naive Bayes: -10 - Prior probability issue
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 10
    
    # SVM Linear: -10 - Margin biased toward majority
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 10
    
    # SVM RBF: -8 - Slightly more robust
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] - 8
    
    # RF/XGBoost: +5 - Can use class weights, stratified sampling
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    
  } else if (imbalance_ratio > 2) {
    assumptions$class_imbalance <- "MODERATE: Some class imbalance"
    
    # Reduced penalties
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] - 8
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] - 5
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] - 5
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 5
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] - 3
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 3
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 3
    
  } else {
    assumptions$class_imbalance <- "OK: Classes reasonably balanced"
    
    # All methods benefit from balanced classes
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] + 5
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] + 5
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 5
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] + 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 5
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
    # RF/XGBoost: smaller boost (they're robust either way)
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 3
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 3
  }
  
  if (verbose) {
    cat(sprintf("  Imbalance ratio (max/min): %.2f\n", imbalance_ratio))
    cat(sprintf("  Normalized entropy: %.3f (1=balanced)\n", entropy))
    cat(sprintf("  Assessment: %s\n\n", assumptions$class_imbalance))
  }
  
  # ============================================================================
  # TEST 8: OUTLIERS AND INFLUENTIAL POINTS
  # ============================================================================
  #
  # RATIONALE:
  # Outliers affect algorithms based on their estimation procedures:
  #
  # - LDA (-15 severe): Mean and covariance estimation are NOT robust.
  #   A single outlier can pull the class mean significantly.
  #   Covariance matrix is even more sensitive to outliers.
  #
  # - QDA (-15 severe): Same sensitivity as LDA, multiplied by having
  #   separate covariance matrices per class.
  #
  # - k-NN (-15 severe): Distance-based. Outliers either become isolated
  #   (misclassified) or can pull nearby points into wrong class if k is small.
  #
  # - Logistic Regression (-10 severe): MLE can be influenced by high-leverage
  #   points. Outliers in feature space with unusual class labels are problematic.
  #
  # - Naive Bayes (-10 severe): Mean/variance estimates affected by outliers.
  #
  # - SVM Linear (-8 severe): Soft margin provides some robustness (slack variables),
  #   but outliers as support vectors can skew the margin.
  #
  # - SVM RBF (+5 severe): RBF kernel effectively isolates outliers since they
  #   are far from other points in kernel space. More robust than linear.
  #
  # - Random Forest (+8 severe): Trees are robust due to:
  #   1) Bagging: outliers only affect ~63% of trees
  #   2) Averaging: outlier influence is diluted
  #   3) Splits: outliers often end up in their own leaves
  #
  # - XGBoost (+5 severe): Robust but slightly less than RF since errors compound.
  # ============================================================================
  
  if (verbose) cat("TEST 8: Outlier Detection (Mahalanobis Distance)\n", strrep("-", 40), "\n")
  
  # Mahalanobis distance for outlier detection
  outlier_result <- tryCatch({
    # Use robust covariance estimate if possible
    center <- colMeans(feature_data, na.rm = TRUE)
    cov_mat <- cov(feature_data, use = "pairwise.complete.obs")
    
    # Regularize if singular
    if (det(cov_mat) < 1e-10) {
      cov_mat <- cov_mat + diag(1e-6, ncol(cov_mat))
    }
    
    md <- mahalanobis(feature_data, center, cov_mat)
    
    # Chi-square cutoff for outliers (97.5th percentile)
    cutoff <- qchisq(0.975, df = n_features)
    n_outliers <- sum(md > cutoff, na.rm = TRUE)
    pct_outliers <- n_outliers / n_samples * 100
    
    list(
      mahalanobis_distances = md,
      cutoff = cutoff,
      n_outliers = n_outliers,
      pct_outliers = pct_outliers
    )
  }, error = function(e) {
    list(n_outliers = NA, pct_outliers = NA, error = as.character(e))
  })
  
  diagnostics$outliers <- outlier_result
  
  if (!is.na(outlier_result$pct_outliers)) {
    if (outlier_result$pct_outliers > 10) {
      assumptions$outliers <- "SEVERE: Many potential outliers detected"
      
      # LDA/QDA: -15 - Mean and covariance estimation highly sensitive
      scores$Score[scores$Algorithm == "LDA"] <- 
        scores$Score[scores$Algorithm == "LDA"] - 15
      scores$Score[scores$Algorithm == "QDA"] <- 
        scores$Score[scores$Algorithm == "QDA"] - 15
      
      # k-NN: -15 - Distance calculations distorted by outliers
      scores$Score[scores$Algorithm == "k-NN"] <- 
        scores$Score[scores$Algorithm == "k-NN"] - 15
      
      # Logistic Regression: -10 - MLE influenced by high-leverage points
      scores$Score[scores$Algorithm == "Logistic Regression"] <- 
        scores$Score[scores$Algorithm == "Logistic Regression"] - 10
      
      # Naive Bayes: -10 - Mean/variance estimates affected
      scores$Score[scores$Algorithm == "Naive Bayes"] <- 
        scores$Score[scores$Algorithm == "Naive Bayes"] - 10
      
      # SVM Linear: -8 - Soft margin helps but not completely
      scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
        scores$Score[scores$Algorithm == "SVM (Linear)"] - 8
      
      # Random Forest: +8 - Robust due to bagging and averaging
      scores$Score[scores$Algorithm == "Random Forest"] <- 
        scores$Score[scores$Algorithm == "Random Forest"] + 8
      
      # XGBoost: +5 - Robust but boosting can amplify outlier errors
      scores$Score[scores$Algorithm == "XGBoost"] <- 
        scores$Score[scores$Algorithm == "XGBoost"] + 5
      
      # SVM RBF: +5 - Kernel isolates outliers in feature space
      scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
        scores$Score[scores$Algorithm == "SVM (RBF)"] + 5
      
    } else if (outlier_result$pct_outliers > 5) {
      assumptions$outliers <- "MODERATE: Some outliers present"
      
      # Reduced penalties/boosts
      scores$Score[scores$Algorithm == "LDA"] <- 
        scores$Score[scores$Algorithm == "LDA"] - 5
      scores$Score[scores$Algorithm == "QDA"] <- 
        scores$Score[scores$Algorithm == "QDA"] - 5
      scores$Score[scores$Algorithm == "k-NN"] <- 
        scores$Score[scores$Algorithm == "k-NN"] - 8
      scores$Score[scores$Algorithm == "Logistic Regression"] <- 
        scores$Score[scores$Algorithm == "Logistic Regression"] - 5
      scores$Score[scores$Algorithm == "Naive Bayes"] <- 
        scores$Score[scores$Algorithm == "Naive Bayes"] - 5
      scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
        scores$Score[scores$Algorithm == "SVM (Linear)"] - 3
      scores$Score[scores$Algorithm == "Random Forest"] <- 
        scores$Score[scores$Algorithm == "Random Forest"] + 5
      scores$Score[scores$Algorithm == "XGBoost"] <- 
        scores$Score[scores$Algorithm == "XGBoost"] + 3
      scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
        scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
      
    } else {
      assumptions$outliers <- "OK: Few outliers detected"
      
      # Clean data benefits all methods, especially those sensitive to outliers
      scores$Score[scores$Algorithm == "k-NN"] <- 
        scores$Score[scores$Algorithm == "k-NN"] + 5
      scores$Score[scores$Algorithm == "LDA"] <- 
        scores$Score[scores$Algorithm == "LDA"] + 3
      scores$Score[scores$Algorithm == "QDA"] <- 
        scores$Score[scores$Algorithm == "QDA"] + 3
      scores$Score[scores$Algorithm == "Logistic Regression"] <- 
        scores$Score[scores$Algorithm == "Logistic Regression"] + 3
      scores$Score[scores$Algorithm == "Naive Bayes"] <- 
        scores$Score[scores$Algorithm == "Naive Bayes"] + 3
      scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
        scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
      scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
        scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
      # RF/XGBoost: neutral (robust either way)
      scores$Score[scores$Algorithm == "Random Forest"] <- 
        scores$Score[scores$Algorithm == "Random Forest"] + 0
      scores$Score[scores$Algorithm == "XGBoost"] <- 
        scores$Score[scores$Algorithm == "XGBoost"] + 0
    }
  } else {
    assumptions$outliers <- "UNKNOWN: Could not compute outlier detection"
  }
  
  if (verbose) {
    if (!is.na(outlier_result$pct_outliers)) {
      cat(sprintf("  Outliers detected: %d (%.1f%%)\n", 
                  outlier_result$n_outliers, outlier_result$pct_outliers))
    }
    cat(sprintf("  Assessment: %s\n\n", assumptions$outliers))
  }
  
  # ============================================================================
  # TEST 9: DATA SKEWNESS
  # ============================================================================
  #
  # RATIONALE:
  # Skewness measures asymmetry of distributions. Many parametric methods assume
  # symmetric (and often normal) distributions.
  #
  # Skewness ranges:
  # - |skew| < 0.5: Approximately symmetric
  # - |skew| 0.5-1: Moderately skewed
  # - |skew| > 1: Highly skewed
  #
  # - LDA (-15 severe): Derived assuming normality, which implies symmetry.
  #   Skewed distributions mean the optimal linear boundary is incorrect.
  #
  # - QDA (-15 severe): Same normality assumption as LDA.
  #
  # - Naive Bayes (-10 severe): Gaussian assumption violated by skewness.
  #   Mean and variance don't fully characterize skewed distributions.
  #
  # - k-NN (-10 severe): Skewness creates uneven point densities.
  #   In a skewed distribution, most points cluster on one side, distorting
  #   the notion of "nearest" neighbors in meaningful ways.
  #
  # - Logistic Regression (-8 severe): Numerical optimization may struggle.
  #   Extreme values from skewed tails can dominate the likelihood.
  #
  # - SVM Linear (-5 severe): Less affected, but extreme values can become
  #   support vectors inappropriately.
  #
  # - SVM RBF (+3 severe): Kernel transformation can reduce skewness effect.
  #
  # - RF/XGBoost (+10 severe): Trees use RANK ORDER for splits, not actual values.
  #   Whether data is x, log(x), or x^2, the split point ranking is identical.
  #   This makes trees completely invariant to monotonic transformations
  #   including those that would change skewness.
  # ============================================================================
  
  if (verbose) cat("TEST 9: Skewness Analysis\n", strrep("-", 40), "\n")
  
  skewness_values <- apply(feature_data, 2, function(x) {
    moments::skewness(x, na.rm = TRUE)
  })
  
  mean_abs_skewness <- mean(abs(skewness_values), na.rm = TRUE)
  pct_highly_skewed <- sum(abs(skewness_values) > 1, na.rm = TRUE) / n_features * 100
  
  diagnostics$skewness <- list(
    feature_skewness = skewness_values,
    mean_abs_skewness = mean_abs_skewness,
    pct_highly_skewed = pct_highly_skewed
  )
  
  if (mean_abs_skewness > 1 || pct_highly_skewed > 50) {
    assumptions$skewness <- "SEVERE: Data is highly skewed"
    
    # LDA/QDA: -15 - Normality (symmetry) assumption violated
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] - 15
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] - 15
    
    # Naive Bayes: -10 - Gaussian assumption violated
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 10
    
    # k-NN: -10 - Uneven point densities distort neighborhoods
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] - 10
    
    # Logistic Regression: -8 - Numerical issues with extreme values
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 8
    
    # SVM Linear: -5 - Extreme values as support vectors
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 5
    
    # SVM RBF: +3 - Kernel helps isolate extreme values
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
    
    # RF/XGBoost: +10 - Completely robust (use rank order only)
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 10
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 10
    
  } else if (mean_abs_skewness > 0.5 || pct_highly_skewed > 25) {
    assumptions$skewness <- "MODERATE: Some features are skewed"
    
    # Reduced penalties/boosts
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] - 8
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] - 8
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 5
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] - 5
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 3
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] - 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 0
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    
  } else {
    assumptions$skewness <- "OK: Data is relatively symmetric"
    
    # Parametric methods benefit from symmetric distributions
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] + 5
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] + 5
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 5
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] + 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
    # RF/XGBoost: neutral (robust either way)
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 0
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 0
  }
  
  if (verbose) {
    cat(sprintf("  Mean absolute skewness: %.3f\n", mean_abs_skewness))
    cat(sprintf("  Highly skewed features (|skew|>1): %.1f%%\n", pct_highly_skewed))
    cat(sprintf("  Assessment: %s\n\n", assumptions$skewness))
  }
  
  # ============================================================================
  # TEST 10: DIMENSIONALITY ASSESSMENT
  # ============================================================================
  #
  # RATIONALE:
  # High dimensionality (p >> n or p > n) creates unique challenges:
  #
  # - QDA (-40 high-D): CANNOT work when p > n/k because class covariance
  #   matrices become singular. Need at least p samples per class for
  #   estimable covariance. Most severely affected method.
  #
  # - LDA (-30 high-D): Pooled covariance becomes singular when p > n-k.
  #   Even when estimable, estimates are highly unreliable.
  #
  # - Logistic Regression (-25 high-D): More parameters (p) than constraints (n)
  #   leads to perfect separation and undefined MLE. Regularization required.
  #
  # - k-NN (-20 high-D): "Curse of dimensionality" - in high-D, all points
  #   become equidistant. Distance metrics lose discriminative power.
  #   The ratio of nearest to farthest neighbor distance → 1 as p → ∞.
  #
  # - Naive Bayes (-15 high-D): Estimates p marginal distributions, but
  #   independence assumption becomes more problematic with more features.
  #
  # - RF/XGBoost (-5 high-D): Feature subsampling (√p) helps, but very
  #   high-D can still cause issues. Each tree sees random subset of features.
  #
  # - SVM Linear (+20 high-D): DESIGNED for high-D! Kernel trick works in
  #   infinite-D space. Maximum margin is well-defined even when p > n.
  #   This is why SVM excels in text classification (10,000+ features).
  #
  # - SVM RBF (+15 high-D): Also handles high-D well, but gamma tuning
  #   becomes more critical. Slightly less boost than linear.
  # ============================================================================
  
  if (verbose) cat("TEST 10: Dimensionality Assessment\n", strrep("-", 40), "\n")
  
  # High dimensionality check
  is_high_dim <- n_features > n_samples
  dim_ratio <- n_features / n_samples
  
  # Effective dimensionality (based on PCA - how many PCs explain 95% variance)
  if (!is.null(pca_result)) {
    var_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
    effective_dim <- which(var_explained >= 0.95)[1]
    if (is.na(effective_dim)) effective_dim <- n_features
  } else {
    effective_dim <- n_features
  }
  
  diagnostics$dimensionality <- list(
    n_features = n_features,
    n_samples = n_samples,
    dim_ratio = dim_ratio,
    is_high_dimensional = is_high_dim,
    effective_dim_95pct = effective_dim
  )
  
  if (is_high_dim) {
    assumptions$dimensionality <- "HIGH-DIMENSIONAL: More features than samples"
    
    # QDA: -40 - Covariance matrices singular, method CANNOT work
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] - 40
    
    # LDA: -30 - Pooled covariance singular or near-singular
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] - 30
    
    # Logistic Regression: -25 - More parameters than samples, separation issues
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 25
    
    # k-NN: -20 - Curse of dimensionality, distances meaningless
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] - 20
    
    # Naive Bayes: -15 - Independence assumption strained
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 15
    
    # RF/XGBoost: -5 - Feature subsampling helps but still some issues
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] - 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] - 5
    
    # SVM Linear: +20 - DESIGNED for high-D, maximum margin works well
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 20
    
    # SVM RBF: +15 - Also handles high-D but gamma tuning critical
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 15
    
  } else if (dim_ratio > 0.5) {
    assumptions$dimensionality <- "MODERATE-HIGH: Feature-to-sample ratio >0.5"
    
    # Reduced penalties/boosts
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] - 15
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] - 10
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] - 10
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] - 10
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] - 5
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 0
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 0
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 10
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 8
    
  } else {
    assumptions$dimensionality <- "OK: Reasonable dimensionality"
    
    # All methods work well with reasonable dimensionality
    scores$Score[scores$Algorithm == "LDA"] <- 
      scores$Score[scores$Algorithm == "LDA"] + 5
    scores$Score[scores$Algorithm == "QDA"] <- 
      scores$Score[scores$Algorithm == "QDA"] + 5
    scores$Score[scores$Algorithm == "Logistic Regression"] <- 
      scores$Score[scores$Algorithm == "Logistic Regression"] + 5
    scores$Score[scores$Algorithm == "Naive Bayes"] <- 
      scores$Score[scores$Algorithm == "Naive Bayes"] + 5
    scores$Score[scores$Algorithm == "k-NN"] <- 
      scores$Score[scores$Algorithm == "k-NN"] + 5
    scores$Score[scores$Algorithm == "Random Forest"] <- 
      scores$Score[scores$Algorithm == "Random Forest"] + 5
    scores$Score[scores$Algorithm == "XGBoost"] <- 
      scores$Score[scores$Algorithm == "XGBoost"] + 5
    scores$Score[scores$Algorithm == "SVM (Linear)"] <- 
      scores$Score[scores$Algorithm == "SVM (Linear)"] + 3
    scores$Score[scores$Algorithm == "SVM (RBF)"] <- 
      scores$Score[scores$Algorithm == "SVM (RBF)"] + 3
  }
  
  if (verbose) {
    cat(sprintf("  Feature/sample ratio: %.2f\n", dim_ratio))
    cat(sprintf("  Effective dimensions (95%% var): %d\n", effective_dim))
    cat(sprintf("  Assessment: %s\n\n", assumptions$dimensionality))
  }
  
  # ============================================================================
  # FINAL RECOMMENDATIONS
  # ============================================================================
  
  # Ensure no negative scores
  scores$Score <- pmax(scores$Score, 0)
  
  # Rank algorithms
  scores <- scores %>% arrange(desc(Score))
  
  if (verbose) {
    cat(strrep("=", 70), "\n")
    cat("ALGORITHM RECOMMENDATIONS (Ranked by Suitability Score)\n")
    cat(strrep("=", 70), "\n\n")
    print(scores)
    cat("\n")
    
    cat("SUMMARY OF ASSUMPTION VIOLATIONS:\n")
    cat(strrep("-", 40), "\n")
    for (name in names(assumptions)) {
      status <- assumptions[[name]]
      if (grepl("VIOLATED|SEVERE|CRITICAL|HIGH-DIM|POOR", status)) {
        cat(sprintf("  ⚠ %s: %s\n", name, status))
      } else if (grepl("MODERATE|PARTIAL|WARNING", status)) {
        cat(sprintf(" ⚡ %s: %s\n", name, status))
      } else {
        cat(sprintf(" ✓ %s: %s\n", name, status))
      }
    }
    
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat(sprintf("RECOMMENDED ALGORITHM: %s (Score: %.0f)\n", 
                scores$Algorithm[1], scores$Score[1]))
    cat(strrep("=", 70), "\n")
  }
  
  # ============================================================================
  # RETURN RESULTS
  # ============================================================================
  
  results <- list(
    diagnostics = diagnostics,
    assumptions = assumptions,
    recommendations = scores,
    recommended_algorithm = scores$Algorithm[1],
    recommended_score = scores$Score[1],
    data_summary = list(
      n_samples = n_samples,
      n_features = n_features,
      n_classes = n_classes,
      class_distribution = class_counts
    )
  )
  
  return(results)
}


#' Helper function: Get algorithm function based on recommendation
#' 
#' @param algorithm_name Name of the recommended algorithm
#' @return A caret-compatible method string

get_algorithm_method <- function(algorithm_name) {
  methods <- c(
    "Random Forest" = "ranger",
    "SVM (RBF)" = "svmRadial",
    "SVM (Linear)" = "svmLinear",
    "LDA" = "lda",
    "QDA" = "qda",
    "k-NN" = "knn",
    "XGBoost" = "xgbTree",
    "Naive Bayes" = "naive_bayes",
    "Logistic Regression" = "multinom"
  )
  
  return(methods[algorithm_name])
}


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

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

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
                               feature_data %>% 
                                 mutate(across(everything(), ~{
                                   min_pos <- min(.[. > 0], na.rm = TRUE)
                                   replace_na(., min_pos * lod_factor)
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
                                  feature_data %>% 
                                    mutate(across(everything(), ~{
                                      x_pos <- . + abs(min(., na.rm = TRUE)) + 1
                                      bc <- MASS::boxcox(lm(x_pos ~ 1), plotit = FALSE)
                                      lambda <- bc$x[which.max(bc$y)]
                                      if (abs(lambda) < 0.01) log(x_pos) else (x_pos^lambda - 1) / lambda
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


#' ============================================================================
#' EVALUATION METRICS FOR MULTI-CLASS CLASSIFICATION
#' ============================================================================

#' Calculate Matthews Correlation Coefficient for multi-class
#' 
#' @param conf_matrix Confusion matrix

multiclass_mcc <- function(conf_matrix) {
  # Use mltools implementation or manual calculation
  if (requireNamespace("mltools", quietly = TRUE)) {
    return(mltools::mcc(confusionM = conf_matrix))
  }
  
  # Manual calculation for multi-class MCC
  n <- sum(conf_matrix)
  correct <- sum(diag(conf_matrix))
  
  row_sums <- rowSums(conf_matrix)
  col_sums <- colSums(conf_matrix)
  
  numerator <- n * correct - sum(row_sums * col_sums)
  denominator <- sqrt(
    (n^2 - sum(row_sums^2)) * (n^2 - sum(col_sums^2))
  )
  
  if (denominator == 0) return(0)
  return(numerator / denominator)
}


#' Calculate comprehensive evaluation metrics
#' 
#' @param predicted Predicted class labels
#' @param actual Actual class labels
#' @param probabilities Optional: class probabilities for AUC calculation

calculate_metrics <- function(predicted, actual, probabilities = NULL) {
  
  # Confusion matrix
  conf_mat <- table(Predicted = predicted, Actual = actual)
  
  # Basic metrics
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  
  # Per-class metrics
  n_classes <- nlevels(actual)
  precision_vec <- numeric(n_classes)
  recall_vec <- numeric(n_classes)
  f1_vec <- numeric(n_classes)
  
  for (i in 1:n_classes) {
    tp <- conf_mat[i, i]
    fp <- sum(conf_mat[i, ]) - tp
    fn <- sum(conf_mat[, i]) - tp
    
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
  
  # Weighted F1
  class_weights <- table(actual) / length(actual)
  weighted_f1 <- sum(f1_vec * class_weights)
  
  # Cohen's Kappa
  p_observed <- accuracy
  p_expected <- sum((rowSums(conf_mat) / sum(conf_mat)) * 
                      (colSums(conf_mat) / sum(conf_mat)))
  kappa <- (p_observed - p_expected) / (1 - p_expected)
  
  # MCC
  mcc <- multiclass_mcc(conf_mat)
  
  # AUC (if probabilities provided)
  auc <- NA
  if (!is.null(probabilities) && ncol(probabilities) > 1) {
    auc <- tryCatch({
      pROC::multiclass.roc(actual, probabilities)$auc[1]
    }, error = function(e) NA)
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
#' EVALUATE SINGLE COMBINATION WITH ADAPTIVE ML
#' ============================================================================

#' Evaluate a single imputation-normalization combination
#' 
#' @param data Original data with missing values
#' @param class_col Class column name
#' @param impute_method Imputation method
#' @param norm_method Normalization method
#' @param ml_diagnostics Pre-computed ML diagnostics (optional)
#' @param use_adaptive_ml Whether to use adaptive ML selection
#' @param cv_folds Number of CV folds
#' @param verbose Print progress

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
  feature_data <- normalized_data %>% dplyr::select(-all_of(class_col))
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
  
  # Step 4: Evaluate with cross-validation or OOB
  class_vector <- as.factor(normalized_data[[class_col]])
  
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
    probabilities <- as.matrix(cv_preds_ordered[, prob_cols])
    
    metrics <- calculate_metrics(predicted, actual, probabilities)
    eval_method <- paste0(cv_folds_actual, "-fold CV")
  }
  
  # Step 5: Calculate additional quality metrics for the transformation
  
  # Cluster resolution (silhouette score)
  cluster_quality <- tryCatch({
    dist_mat <- dist(feature_data)
    sil <- cluster::silhouette(as.numeric(class_vector), dist_mat)
    mean(sil[, "sil_width"])
  }, error = function(e) NA)
  
  # Correlation with original data structure
  original_features <- data %>% dplyr::select(-all_of(class_col))
  correlation_preserved <- tryCatch({
    # Compare correlation structures
    orig_cor <- cor(original_features, use = "pairwise.complete.obs")
    new_cor <- cor(feature_data)
    mean(abs(orig_cor - new_cor), na.rm = TRUE)
  }, error = function(e) NA)
  
  # Distributional normality (post-transformation)
  normality_improvement <- tryCatch({
    sw_pvals <- sapply(feature_data, function(x) {
      if (length(unique(x)) > 2 && length(x) <= 5000) {
        shapiro.test(x)$p.value
      } else NA
    })
    mean(sw_pvals > 0.05, na.rm = TRUE)
  }, error = function(e) NA)
  
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
    macro_f1 = metrics$macro_f1,
    weighted_f1 = metrics$weighted_f1,
    auc = metrics$auc,
    
    # Data quality metrics
    cluster_quality = cluster_quality,
    correlation_change = correlation_preserved,
    normality_pct = normality_improvement,
    
    # Diagnostics
    ml_diagnostics = diag_results,
    confusion_matrix = metrics$confusion_matrix
  ))
}


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


#' ============================================================================
#' CONVENIENCE FUNCTIONS
#' ============================================================================

#' Apply the best combination to data
#' 
#' @param data Original data with missing values
#' @param class_col Class column name
#' @param best_combo Result from find_best_impute_normalize

apply_best_combo <- function(data, class_col, best_combo) {
  imputed <- apply_imputation(data, best_combo$imputation, class_col)
  normalized <- apply_normalization(imputed, best_combo$normalization, class_col)
  return(normalized)
}

two_stage_optimization <- function(data, class_col, top_k = 5) {
  # Stage 1: Quick screen with Random Forest OOB
  quick_results <- find_best_impute_normalize(
    data = data,
    class_col = class_col,
    use_adaptive_ml = FALSE,  # Use RF for all
    verbose = FALSE
  )
  
  # Get top K combinations
  top_combos <- head(quick_results$summary_table, top_k)
  
  # Stage 2: Full adaptive evaluation on top candidates
  refined_results <- list()
  for (i in 1:nrow(top_combos)) {
    cat(sprintf("Refining combo %d/%d: %s + %s\n", 
                i, nrow(top_combos),
                top_combos$imputation[i], 
                top_combos$normalization[i]))
    
    result <- evaluate_combo(
      data = data,
      class_col = class_col,
      impute_method = top_combos$imputation[i],
      norm_method = top_combos$normalization[i],
      use_adaptive_ml = TRUE,
      cv_folds = 10  # More thorough CV for final candidates
    )
    refined_results[[i]] <- result
  }
  
  return(refined_results)
}


# ============================================================================
# ALTERNATIVE EVALUATION METRICS for OOB MCC
# ============================================================================

#' Beyond OOB MCC: Alternative Metrics for Imputation/Normalization Selection
#' ==========================================================================
#' 
#' When selecting the best imputation/normalization combo, we need metrics that
#' capture BOTH:
#' A) Predictive performance (classification quality)
#' B) Data quality preservation (statistical properties)
#' 
#' Here are the recommended metrics organized by category:

#' ============================================================================
#' A) CLASSIFICATION-BASED METRICS
#' ============================================================================
#' 
#' 1. Matthews Correlation Coefficient (MCC) - YOUR CURRENT CHOICE
#'    - Pros: Balanced measure, works well with imbalanced classes
#'    - Cons: Can be undefined with certain confusion matrix configurations
#'    - Range: [-1, 1], higher is better
#' 
#' 2. Cohen's Kappa
#'    - Pros: Accounts for chance agreement
#'    - Cons: Can be pessimistic with imbalanced classes
#'    - Range: [-1, 1], >0.8 is excellent
#' 
#' 3. Balanced Accuracy (BA)
#'    - Formula: Mean of per-class recalls
#'    - Pros: Simple, intuitive, handles imbalance
#'    - Range: [0, 1]
#' 
#' 4. Macro-averaged F1 Score
#'    - Pros: Balances precision and recall across classes
#'    - Cons: All classes weighted equally (may not reflect importance)
#' 
#' 5. Weighted F1 Score
#'    - Pros: Accounts for class sizes
#'    - Use when: Some classes are more important/frequent

#' IMPLEMENTATION: Multi-metric scoring function

calculate_all_classification_metrics <- function(predicted, actual, probabilities = NULL) {
  conf_mat <- table(Predicted = predicted, Actual = actual)
  n_classes <- nlevels(actual)
  n_samples <- length(actual)
  
  # Basic accuracy
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  
  # Per-class metrics
  tp <- diag(conf_mat)
  fp <- rowSums(conf_mat) - tp
  fn <- colSums(conf_mat) - tp
  tn <- sum(conf_mat) - tp - fp - fn
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  specificity <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
  f1 <- ifelse(precision + recall == 0, 0, 
               2 * precision * recall / (precision + recall))
  
  # Balanced accuracy
  balanced_accuracy <- mean(recall)
  
  # Macro/Weighted averages
  class_sizes <- table(actual)
  class_weights <- class_sizes / n_samples
  
  macro_precision <- mean(precision)
  macro_recall <- mean(recall)
  macro_f1 <- mean(f1)
  
  weighted_precision <- sum(precision * class_weights)
  weighted_recall <- sum(recall * class_weights)
  weighted_f1 <- sum(f1 * class_weights)
  
  # Cohen's Kappa
  p_observed <- accuracy
  p_expected <- sum((rowSums(conf_mat) / n_samples) * 
                      (colSums(conf_mat) / n_samples))
  kappa <- (p_observed - p_expected) / (1 - p_expected + 1e-10)
  
  # MCC (multiclass)
  mcc <- multiclass_mcc(conf_mat)
  
  # G-mean (geometric mean of recalls)
  g_mean <- exp(mean(log(recall + 1e-10)))
  
  return(list(
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    kappa = kappa,
    mcc = mcc,
    g_mean = g_mean,
    macro_precision = macro_precision,
    macro_recall = macro_recall,
    macro_f1 = macro_f1,
    weighted_precision = weighted_precision,
    weighted_recall = weighted_recall,
    weighted_f1 = weighted_f1,
    per_class = data.frame(
      class = levels(actual),
      precision = precision,
      recall = recall,
      specificity = specificity,
      f1 = f1
    )
  ))
}


#' ============================================================================
#' B. DATA QUALITY / TRANSFORMATION METRICS
#' ============================================================================
#' 
#' These metrics evaluate how well the transformation preserves or improves
#' data quality, independent of a specific ML model.
#' 
#' 1. CLUSTER SEPARATION QUALITY
#'    Silhouette Score: Measures how well samples cluster by class
#'    Range: [-1, 1], higher is better

calculate_cluster_quality <- function(transformed_data, class_labels) {
  library(cluster)
  
  dist_matrix <- dist(transformed_data)
  sil <- silhouette(as.numeric(as.factor(class_labels)), dist_matrix)
  
  return(list(
    mean_silhouette = mean(sil[, "sil_width"]),
    per_class_silhouette = tapply(sil[, "sil_width"], class_labels, mean)
  ))
}


#' 2. DISCRIMINANT RATIO (Linear Separability)
#'    Ratio of between-class to within-class variance

calculate_discriminant_ratio <- function(transformed_data, class_labels) {
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
    class_mean <- colMeans(class_data)
    
    # Between-class
    between_scatter <- between_scatter + 
      n_cls * sum((class_mean - grand_mean)^2)
    
    # Within-class
    within_scatter <- within_scatter + 
      sum(apply(class_data, 1, function(x) sum((x - class_mean)^2)))
  }
  
  ratio <- between_scatter / (within_scatter + 1e-10)
  return(ratio)
}


#' 3. CORRELATION STRUCTURE PRESERVATION
#'    How well does the transformation preserve original variable relationships?

calculate_correlation_preservation <- function(original_data, transformed_data) {
  # Handle missing values in original
  original_cor <- cor(original_data, use = "pairwise.complete.obs")
  transformed_cor <- cor(transformed_data)
  
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
#'    Percentage of features that pass normality tests after transformation

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
    if (length(unique(x)) < 3) return(NA)
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
#'    How well does the transformation reduce skewness?

calculate_skewness_reduction <- function(original_data, transformed_data) {
  library(moments)
  
  original_skew <- sapply(original_data, function(x) skewness(x, na.rm = TRUE))
  transformed_skew <- sapply(transformed_data, skewness)
  
  mean_abs_original <- mean(abs(original_skew), na.rm = TRUE)
  mean_abs_transformed <- mean(abs(transformed_skew), na.rm = TRUE)
  
  return(list(
    original_mean_abs_skew = mean_abs_original,
    transformed_mean_abs_skew = mean_abs_transformed,
    skew_reduction = mean_abs_original - mean_abs_transformed,
    pct_reduction = (mean_abs_original - mean_abs_transformed) / mean_abs_original * 100
  ))
}


#' 6. INFORMATION PRESERVATION (for imputation)
#'    Kolmogorov-Smirnov test comparing distributions

calculate_distribution_preservation <- function(original_data, imputed_data, class_col = NULL) {
  
  if (!is.null(class_col)) {
    original_data <- original_data[, !names(original_data) %in% class_col]
    imputed_data <- imputed_data[, !names(imputed_data) %in% class_col]
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


#' ============================================================================
#' C) COMPOSITE SCORE: COMBINING MULTIPLE METRICS
#' ============================================================================

#' For selecting the best imputation/normalization combo, I recommend using
#' a composite score that balances multiple objectives:

calculate_composite_score <- function(
    classification_metrics,  # List with accuracy, mcc, macro_f1, etc.
    data_quality_metrics,    # List with silhouette, discriminant_ratio, etc.
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
  
  # Normalize metrics to [0, 1] range where needed
  score <- 0
  
  # Classification metrics (already in good ranges)
  score <- score + weights$mcc * (classification_metrics$mcc + 1) / 2  # MCC: [-1,1] -> [0,1]
  score <- score + weights$balanced_accuracy * classification_metrics$balanced_accuracy
  score <- score + weights$macro_f1 * classification_metrics$macro_f1
  
  # Data quality metrics
  if (!is.null(data_quality_metrics$silhouette)) {
    score <- score + weights$silhouette * 
      (data_quality_metrics$silhouette + 1) / 2  # Silhouette: [-1,1] -> [0,1]
  }
  
  if (!is.null(data_quality_metrics$discriminant_ratio)) {
    # Normalize discriminant ratio (log transform, cap at reasonable max)
    dr_norm <- min(log1p(data_quality_metrics$discriminant_ratio) / 3, 1)
    score <- score + weights$discriminant_ratio * dr_norm
  }
  
  if (!is.null(data_quality_metrics$normality_improvement)) {
    # Already in percentage, normalize to [0, 1]
    norm_score <- max(0, min(1, data_quality_metrics$normality_improvement / 100))
    score <- score + weights$normality_improvement * norm_score
  }
  
  if (!is.null(data_quality_metrics$skew_reduction)) {
    # Positive reduction is good, normalize
    skew_score <- max(0, min(1, data_quality_metrics$skew_reduction / 2))
    score <- score + weights$skew_reduction * skew_score
  }
  
  return(score)
}


# ============================================================================
# SECTION 3: RECOMMENDED WORKFLOW FOR Imputation/Normalization Selection
# ============================================================================

#' RECOMMENDED COMPLETE WORKFLOW
#' =============================

recommended_imp_norm_workflow <- function(data, class_col) {
  
  # -------------------------------------------------------------------------
  # STAGE 1: Data Exploration & Quality Check
  # -------------------------------------------------------------------------
  cat("\n=== STAGE 1: Initial Data Assessment ===\n")
  
  # Run initial diagnostics on raw data
  initial_diagnostics <- select_ML_diag_tests(
    data = data,
    class_col = class_col,
    verbose = TRUE
  )
  
  # -------------------------------------------------------------------------
  # STAGE 2: Imputation/Normalization Selection (Two-Stage Approach)
  # -------------------------------------------------------------------------
  cat("\n=== STAGE 2: Imputation & Normalization Optimization ===\n")
  
  # Stage 2a: Quick screen
  cat("Stage 2a: Initial screening with Random Forest OOB...\n")
  
  quick_results <- find_best_impute_normalize(
    data = data,
    class_col = class_col,
    imputation_methods = c("half_min", "median", "knn", "missforest"),
    normalization_methods = c("none", "log", "log10", "sqrt", "zscore", "pareto"),
    use_adaptive_ml = FALSE,  # Use RF for all (faster)
    cv_folds = 5,
    ranking_metric = "mcc",
    verbose = TRUE
  )
  
  # Stage 2b: Refine top candidates with adaptive ML
  cat("\nStage 2b: Refining top 5 candidates with adaptive ML...\n")
  
  top_combos <- head(quick_results$summary_table, 5)
  
  refined_results <- lapply(1:nrow(top_combos), function(i) {
    cat(sprintf("  Evaluating %s + %s...\n", 
                top_combos$imputation[i], 
                top_combos$normalization[i]))
    
    evaluate_combo(
      data = data,
      class_col = class_col,
      impute_method = top_combos$imputation[i],
      norm_method = top_combos$normalization[i],
      use_adaptive_ml = TRUE,
      cv_folds = 10
    )
  })
  
  # -------------------------------------------------------------------------
  # STAGE 3: Final Selection with Composite Scoring
  # -------------------------------------------------------------------------
  cat("\n=== STAGE 3: Final Selection with Composite Score ===\n")
  
  # Calculate comprehensive metrics for each refined candidate
  final_scores <- sapply(refined_results, function(result) {
    if (!result$success) return(-Inf)
    
    # Classification metrics
    class_metrics <- list(
      mcc = result$mcc,
      balanced_accuracy = (result$accuracy + result$macro_recall) / 2,
      macro_f1 = result$macro_f1
    )
    
    # Data quality metrics (from result)
    quality_metrics <- list(
      silhouette = result$cluster_quality,
      normality_improvement = result$normality_pct * 100,
      skew_reduction = 0  # Would need to calculate
    )
    
    calculate_composite_score(class_metrics, quality_metrics)
  })
  
  best_idx <- which.max(final_scores)
  best_result <- refined_results[[best_idx]]
  
  # -------------------------------------------------------------------------
  # STAGE 4: Final Validation
  # -------------------------------------------------------------------------
  cat("\n=== STAGE 4: Final Validation ===\n")
  cat(sprintf("Best combination: %s + %s\n", 
              best_result$impute_method, 
              best_result$norm_method))
  cat(sprintf("Recommended algorithm: %s\n", best_result$recommended_algorithm))
  cat(sprintf("Composite score: %.4f\n", final_scores[best_idx]))
  
  return(list(
    initial_diagnostics = initial_diagnostics,
    quick_screen_results = quick_results,
    refined_results = refined_results,
    final_scores = final_scores,
    best_result = best_result,
    best_combination = list(
      imputation = best_result$impute_method,
      normalization = best_result$norm_method,
      algorithm = best_result$recommended_algorithm
    )
  ))
}


# ============================================================================
# SECTION 4: SUMMARY OF RECOMMENDATIONS
# ============================================================================

#' SUMMARY: KEY RECOMMENDATIONS FOR YOUR WORKFLOW
#' ===============================================
#' 
#' 1. YOUR SOLUTION IS VALID
#'    - Integrating select_ML_diag_tests() into find_best_impute_normalize()
#'      is a sound approach
#'    - Different transformations DO change data characteristics
#'    - Adaptive algorithm selection addresses this well
#' 
#' 2. IMPROVEMENTS TO CONSIDER
#'    a) Two-stage approach to reduce computational cost
#'    b) Use nested CV to avoid overfitting to evaluation metric
#'    c) Consider multiple metrics in a composite score
#' 
#' 3. RECOMMENDED EVALUATION METRICS
#'    Primary: MCC or Balanced Accuracy (for classification)
#'    Secondary: Silhouette score (for cluster quality)
#'    Supporting: Normality improvement, skewness reduction, K-S test
#' 
#' 4. FOR LEFT-SKEWED DATA SPECIFICALLY
#'    - Log transforms (log, log2, log10) often work well
#'    - Square root transform is gentler alternative
#'    - Box-Cox transformation can find optimal power automatically
#'    - Random Forest and XGBoost handle skewness well even without normalization
#' 
#' 5. FOR MISSING DATA
#'    - kNN imputation preserves local structure
#'    - missForest captures non-linear relationships
#'    - Half-minimum or LOD imputation is appropriate for detection limits
#'    - ALWAYS validate that imputation doesn't distort class distributions
#' 
#' 6. FINAL NOTE ON ALGORITHM SELECTION
#'    Given your description (missing values + left-skewed data):
#'    - Random Forest/XGBoost: Most robust, likely top performers
#'    - SVM (RBF): Good if data becomes more separable after transformation
#'    - LDA: Consider if transformations achieve near-normality
#'    - Naive Bayes: Avoid unless independence assumption holds

#' ============================================================================
#' QUICK REFERENCE TABLE: When to Use Each Metric
#' ============================================================================
#'
#' | Situation                      | Primary Metric    | Secondary Metrics    |
#' |-------------------------------|-------------------|---------------------|
#' | Balanced classes              | Accuracy, F1      | AUC, MCC            |
#' | Imbalanced classes            | MCC, Bal. Accuracy| Macro F1, Kappa     |
#' | Care about false positives    | Precision         | Specificity         |
#' | Care about false negatives    | Recall            | Sensitivity         |
#' | Evaluating transformations    | Composite Score   | Silhouette, K-S     |
#' | Final model selection         | MCC + CV variance | Per-class metrics   |
#' ============================================================================

cat("\n=== Functions loaded successfully ===\n")
cat("Available functions:\n")
cat("  - select_ML_diag_tests(): Run diagnostic tests for ML selection\n")
cat("  - find_best_impute_normalize(): Find best imputation/normalization combo\n")
cat("  - recommended_workflow(): Complete recommended workflow\n")
cat("  - calculate_composite_score(): Multi-metric composite scoring\n")
cat("\nExample usage:\n")
cat("  results <- select_ML_diag_tests(my_data, class_col = 'Product_Type')\n")
cat("  best <- find_best_impute_normalize(my_data, class_col = 'Product_Type')\n")

