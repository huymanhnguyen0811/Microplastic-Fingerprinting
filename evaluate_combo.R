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
  
  # Step 4: Evaluate with cross-validation or OOB
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
  
  # Step 5: Calculate additional quality metrics for the transformation
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

evaluate_combo_with_score <- function(..., 
                                      score_weights = NULL) {
  
  # Run evaluate_combo
  result <- evaluate_combo(...)
  
  # Calculate composite score
  if (!is.null(score_weights)) {
    result$composite_score <- calculate_composite_score(result, weights = score_weights)
  } else {
    result$composite_score <- calculate_composite_score(result)
  }
  
  return(result)
}