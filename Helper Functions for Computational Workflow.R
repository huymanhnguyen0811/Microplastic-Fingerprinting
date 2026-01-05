# SPDX-License-Identifier: AGPL-3.0-only

# This code contains helper functions for Complete Workflow.Rmd file.
# Publication: "Chemical Fingerprints of New vs. Weathered Microplastics: A Machine Learning Approach"

# Functions -------------------------------------------------------------------------------------------------------
'%notin%' <- Negate('%in%')

# Groups compounds based on retention time (RT) and mass-to-charge ratio (m/z)
# Uses a fast binning algorithm to efficiently find compounds with similar RT and m/z values
#
# Parameters:
#   data: Data frame containing RT and m.z columns
#   split_time: RT threshold for splitting data into two regions (HPLCTOFMS only)
#   rtthres1: RT tolerance for first region (or entire dataset for ATDGCMS)
#   rtthres2: RT tolerance for second region (HPLCTOFMS only)
#   mzthres: m/z tolerance for grouping
#   type: Either "ATDGCMS" or "HPLCTOFMS"
#
# Returns: Data frame with added Feature column containing group assignments
grouping_comp <- function(data, split_time, rtthres1, rtthres2, mzthres, type) {
  # Initialize data and Feature column
  dat <- data[,, drop = FALSE]
  dat$Feature <- NA_character_
  
  # Validate input type
  if (!(type %in% c("ATDGCMS", "HPLCTOFMS"))) {
    stop("Invalid type specified. Must be either 'ATDGCMS' or 'HPLCTOFMS'")
  }
  
  # Core grouping algorithm: bins by RT, then scans for m/z matches within bins
  assign_groups <- function(dt, rtthres, mzthres, prefix, type_label) {
    n <- nrow(dt)
    if (n == 0) {
      return(dt)
    }
    
    # Extract RT and m/z values as numeric vectors
    rt <- base::as.numeric(dt$RT)
    mz <- base::as.numeric(dt$m.z)
    feat <- rep(NA_character_, n)
    
    # Identify valid (finite) data points
    valid <- is.finite(rt) & is.finite(mz)
    if (!any(valid)) {
      dt$Feature <- feat
      return(dt)
    }
    
    # Validate thresholds
    if (rtthres <= 0 || mzthres <= 0) {
      stop("RT/m.z thresholds must be > 0.")
    }
    
    # BIN CREATION PHASE
    # Divide RT space into bins of size rtthres for faster neighbor searching
    # This avoids O(n²) pairwise comparisons by only checking nearby bins
    min_rt <- min(rt[valid])
    rt_bin <- rep(NA_integer_, n)
    rt_bin[valid] <- floor((rt[valid] - min_rt) / rtthres)
    
    # Group data points by their RT bin
    bin_groups <- split(which(valid), rt_bin[valid], drop = TRUE)
    
    # PRE-SORT PHASE
    # Sort m/z values within each bin to enable binary search later
    bin_idx <- vector("list", length(bin_groups))
    bin_mz <- vector("list", length(bin_groups))
    bin_keys <- names(bin_groups)
    
    for (k in seq_along(bin_groups)) {
      idx <- bin_groups[[k]]
      mz_vals <- mz[idx]
      ord <- order(mz_vals)
      bin_idx[[k]] <- idx[ord]     # Sorted indices
      bin_mz[[k]] <- mz_vals[ord]  # Sorted m/z values
    }
    names(bin_idx) <- bin_keys
    names(bin_mz) <- bin_keys
    
    # GROUP ASSIGNMENT PHASE
    # Scan through all points and assign group IDs to connected components
    counter <- 1L
    
    for (i in seq_len(n)) {
      # Skip if invalid or already assigned
      if (!valid[i] || !is.na(feat[i])) {
        next
      }
      
      rt_i <- rt[i]
      mz_i <- mz[i]
      rtb <- rt_bin[i]
      if (is.na(rtb)) {
        next
      }
      
      # NEIGHBOR SEARCH
      # Check current bin and adjacent bins (offset -1, 0, +1) for RT proximity
      # This ensures we catch compounds near bin boundaries
      cand_parts <- vector("list", 3)
      part_idx <- 0L
      
      for (offset in -1:1) {
        key <- as.character(rtb + offset)
        idx_sorted <- bin_idx[[key]]
        
        if (length(idx_sorted) == 0) {
          next
        }
        
        mz_sorted <- bin_mz[[key]]
        
        # Binary search for m/z range [mz_i - mzthres, mz_i + mzthres]
        # findInterval efficiently locates the boundaries in sorted array
        lo <- findInterval(mz_i - mzthres, mz_sorted) + 1L
        hi <- findInterval(mz_i + mzthres, mz_sorted)
        
        if (hi >= lo) {
          part_idx <- part_idx + 1L
          cand_parts[[part_idx]] <- idx_sorted[lo:hi]
        }
      }
      
      # If no candidates found in any bin, skip
      if (part_idx == 0L) {
        next
      }
      
      # Combine candidates from all checked bins
      cand <- unlist(cand_parts[seq_len(part_idx)], use.names = FALSE)
      if (length(cand) == 0) {
        next
      }
      
      # Filter to unassigned candidates only
      cand <- cand[is.na(feat[cand])]
      if (length(cand) == 0) {
        next
      }
      
      # Apply exact RT and m/z thresholds (binary search gave us a superset)
      cand <- cand[
        abs(rt[cand] - rt_i) <= rtthres & abs(mz[cand] - mz_i) <= mzthres
      ]
      
      # Assign group ID to all matching candidates
      if (length(cand) > 0) {
        feat[cand] <- paste0(prefix, counter, ".", type_label)
        counter <- counter + 1L
      }
    }
    
    dt$Feature <- feat
    dt
  }
  
  # MAIN EXECUTION
  # For HPLC data: split by RT and apply different thresholds to each region
  # For ATDGC data: process entire dataset with single threshold
  if (type %in% "HPLCTOFMS") {
    dat1 <- dat %>% filter(RT <= split_time)
    dat2 <- dat %>% filter(RT > split_time)
    dat1 <- assign_groups(dat1, rtthres1, mzthres, "Compound_RT1_", type)
    dat2 <- assign_groups(dat2, rtthres2, mzthres, "Compound_RT2_", type)
    dat <- rbind(dat1, dat2)
  } else {
    dat <- assign_groups(dat, rtthres1, mzthres, "Compound_", type)
  }
  
  return(dat)
}

# =====================================================================================
# Helper function for feature reduction (regular; modified 80% rule and in-house option)
# =====================================================================================
process_and_filter <- function(
  df_input,
  class_col = "Plastic_type",
  threshold = 0.8,
  filter_mode = "regular80"
) {
  # "regular80", "modified80", "inhouse", or "none"

  # --- Step B: Pivot to Wide Format ---
  df_wide <- df_input %>%
    dplyr::select(
      File,
      Source,
      Polymer,
      Feature,
      technique,
      Values,
      all_of(class_col)
    ) %>%
    dplyr::group_by(
      File,
      Source,
      Polymer,
      technique,
      Feature,
      across(all_of(class_col))
    ) %>%
    dplyr::summarise(Values = mean(Values, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Feature, values_from = Values) %>%
    tibble::column_to_rownames(var = "File") %>%
    dplyr::relocate(
      any_of(c(class_col, "technique", "Source", "Polymer")),
      .before = 1
    )

  # Metadata columns
  meta_cols <- c("Plastic_type", "technique", "Source", "Polymer")

  # --- Step D: Apply Feature Filtering Logic ---

  feature_data <- df_wide %>%
    dplyr::select(all_of(setdiff(names(df_wide), meta_cols)))
  kept_features <- c()

  if (filter_mode == "none") {
    cat("No filtering applied (filter_mode = 'none').\n")
    kept_features <- names(feature_data)
  } else if (filter_mode == "regular80") {
    # Mode 1: Global 80% Rule
    cat("Applying 'Regular' filtering (Global threshold)...\n")
    presence_rate <- colMeans(!is.na(feature_data))
    kept_features <- names(presence_rate[presence_rate >= threshold])
  } else if (filter_mode == "modified80") {
    # Mode 2: Class-Based 80% Rule
    cat(paste0(
      "Applying 'Class' filtering (Threshold met in at least one '",
      class_col,
      "')...\n"
    ))

    class_split <- split(df_wide, df_wide[[class_col]])
    class_kept_list <- list()

    for (cls in names(class_split)) {
      cls_idx <- which(df_wide[[class_col]] == cls)
      cls_dat <- feature_data[cls_idx, ]
      cls_rate <- colMeans(!is.na(cls_dat))
      class_kept_list[[cls]] <- names(cls_rate[cls_rate >= threshold])
    }
    kept_features <- Reduce(union, class_kept_list)
  } else if (filter_mode == "inhouse") {
    # Mode 3: In-House Filtering
    cat("Applying 'In-House' filtering strategy...\n")

    # 1. Calculate Statistics per Feature
    classes <- unique(df_wide[[class_col]])
    n_feats <- ncol(feature_data)

    # Pre-allocate matrices for speed
    class_presence_mat <- matrix(0, nrow = length(classes), ncol = n_feats)
    class_counts_mat <- matrix(0, nrow = length(classes), ncol = n_feats)

    for (i in seq_along(classes)) {
      cls <- classes[i]
      cls_rows <- which(df_wide[[class_col]] == cls)
      cls_dat <- feature_data[cls_rows, , drop = FALSE]

      # A. % Present in this class (for Rule 1)
      class_presence_mat[i, ] <- colMeans(!is.na(cls_dat))

      # B. Exact Count of data points in this class (for Rule 4)
      class_counts_mat[i, ] <- colSums(!is.na(cls_dat))
    }

    # Feature Stats
    global_presence <- colMeans(!is.na(feature_data))
    max_class_presence <- apply(class_presence_mat, 2, max)

    total_data_points <- colSums(class_counts_mat) # Sum of all counts
    n_classes_with_data <- colSums(class_counts_mat > 0) # How many classes have at least 1 point

    # Key Metric for Rule 4: How many classes have >= 3 data points?
    n_classes_with_3plus <- colSums(class_counts_mat >= 3)

    # 2. Apply The 4 Rules Logic

    # Rule 1: Missingness (Keep if Global > 10% OR Max Class > 90%)
    keep_r1 <- (global_presence > 0.10) | (max_class_presence >= 0.90)

    # Rule 2: Singleton Sample (Keep if > 1 data point total)
    keep_r2 <- total_data_points > 1

    # Rule 3: Singleton Category (Keep if > 1 category... UNLESS diagnostic)
    is_diagnostic <- max_class_presence >= 0.90
    keep_r3 <- (n_classes_with_data > 1) | is_diagnostic

    # Rule 4: Sparse Distribution
    # Logic: If present in >= 2 categories, we REQUIRE at least 2 categories to have >= 3 points.
    # If present in < 2 categories (i.e. 1), this rule is skipped (True).
    keep_r4 <- (n_classes_with_data < 2) | (n_classes_with_3plus >= 2)

    # FINAL DECISION: Feature must pass ALL checks
    final_mask <- keep_r1 & keep_r2 & keep_r3 & keep_r4
    kept_features <- colnames(feature_data)[final_mask]
  } else {
    stop(
      "Invalid filter_mode. Choose 'regular80', 'modified80', 'inhouse', or 'none'."
    )
  }

  # --- Step E: Build Final Dataframe ---
  final_cols <- unique(c(meta_cols, kept_features))
  df_filtered <- df_wide %>% dplyr::select(any_of(final_cols))

  return(list(
    original_data = df_wide,
    filtered_data = df_filtered,
    kept_features = kept_features
  ))
}

# Find shared features between SB and ENV
shared_feature_filtering <- function(dataset, feat_reduc_mode) {
  SB <- dataset %>% dplyr::filter(Source == "Store-Bought")
  ENV <- dataset %>% dplyr::filter(Source == "Environmental")

  shared_cols_SB_ENV <- Reduce(
    intersect,
    list(unique(SB$Feature), unique(ENV$Feature))
  )

  combined_SB_ENV_shared_cols_ntr <- rbind(
    SB %>%
      dplyr::filter(Feature %in% shared_cols_SB_ENV),
    ENV %>%
      dplyr::filter(Feature %in% shared_cols_SB_ENV)
  )
  output <- process_and_filter(
    combined_SB_ENV_shared_cols_ntr,
    class_col = "Plastic_type",
    threshold = 0.8,
    filter_mode = feat_reduc_mode
  )
  return(output)
}

# ==========================================================
# Helper function for # Imputation and Normalization methods
# ==========================================================

impute_train_apply <- function(X_train, X_test, imp_name, k = 5) {
  if (imp_name == "mean_imputation") {
    mu <- sapply(X_train, function(x) mean(x, na.rm = TRUE))
    fill <- function(X) {
      X2 <- X
      for (nm in names(X2)) {
        idx <- is.na(X2[[nm]])
        if (any(idx)) X2[[nm]][idx] <- mu[[nm]]
      }
      X2
    }
    return(list(train = fill(X_train), test = fill(X_test)))
  }

  if (imp_name == "median_imputation") {
    md <- sapply(X_train, function(x) median(x, na.rm = TRUE))
    fill <- function(X) {
      X2 <- X
      for (nm in names(X2)) {
        idx <- is.na(X2[[nm]])
        if (any(idx)) X2[[nm]][idx] <- md[[nm]]
      }
      X2
    }
    return(list(train = fill(X_train), test = fill(X_test)))
  }

  if (imp_name == "zero_imputation") {
    fill <- function(X) {
      X2 <- X
      for (nm in names(X2)) {
        idx <- is.na(X2[[nm]])
        if (any(idx)) X2[[nm]][idx] <- 0
      }
      X2
    }
    return(list(train = fill(X_train), test = fill(X_test)))
  }

  if (imp_name == "half_min_imputation") {
    halfmins <- sapply(X_train, function(x) {
      nz <- x[!is.na(x) & x > 0]
      if (length(nz) > 0) min(nz) / 2 else 0
    })
    fill <- function(X) {
      X2 <- X
      for (nm in names(X2)) {
        idx <- is.na(X2[[nm]])
        if (any(idx)) X2[[nm]][idx] <- halfmins[[nm]]
      }
      X2
    }
    return(list(train = fill(X_train), test = fill(X_test)))
  }

  # --- OPTIMIZED Complex Imputations ---

  if (imp_name == "random_forest_imputation") {
    # Identify columns with missing values in TRAIN
    na_cols <- names(which(colSums(is.na(X_train)) > 0))

    # We copy the datasets to avoid modifying originals in place
    train_imp <- X_train
    test_imp <- X_test

    # We must start with a rough initialization (e.g. mean) for the PREDICTORS
    train_filled <- X_train
    for (j in 1:ncol(train_filled)) {
      # Safety check: if column is all NA, fill with 0 to prevent crash
      if (all(is.na(train_filled[, j]))) {
        train_filled[, j] <- 0
      } else {
        train_filled[is.na(train_filled[, j]), j] <- mean(
          train_filled[, j],
          na.rm = TRUE
        )
      }
    }

    for (target in na_cols) {
      # 1. Train model only on observed data in X_train
      valid_rows <- !is.na(X_train[[target]])

      # Skip if column is entirely NA
      if (sum(valid_rows) == 0) {
        next
      }

      # [FIX] Use dependent.variable.name instead of formula
      # This prevents crashes if 'target' name has spaces or symbols (e.g., "m/z")
      rf_model <- ranger::ranger(
        dependent.variable.name = target,
        data = train_filled[valid_rows, ],
        num.trees = 50,
        min.node.size = 5,
        verbose = FALSE,
        num.threads = 1 # Avoid nested parallelism if called inside parallel loop
      )

      # 2. Predict missing slots in Train
      idx_train <- is.na(X_train[[target]])
      if (any(idx_train)) {
        preds_train <- predict(
          rf_model,
          data = train_filled[idx_train, ]
        )$predictions
        train_imp[idx_train, target] <- preds_train
      }

      # 3. Predict missing slots in Test
      test_filled <- X_test
      for (j in names(test_filled)) {
        if (j != target) {
          # If test col has NA, fill with Train Mean (to prevent leakage)
          mu_val <- mean(X_train[[j]], na.rm = TRUE)
          # Handle edge case where train col was all NA
          if (is.na(mu_val)) {
            mu_val <- 0
          }
          test_filled[is.na(test_filled[, j]), j] <- mu_val
        }
      }

      idx_test <- is.na(X_test[[target]])
      if (any(idx_test)) {
        preds_test <- predict(
          rf_model,
          data = test_filled[idx_test, ]
        )$predictions
        test_imp[idx_test, target] <- preds_test
      }
    }

    return(list(train = train_imp, test = test_imp))
  }

  if (imp_name == "knn_imputation") {
    # Optimization: Check for parallel backend.
    if (!getDoParRegistered()) {
      message(
        "Tip: Register a parallel backend (doParallel) to speed up KNN imputation."
      )
    }

    # [FIX] CRITICAL: Remove Zero-Variance Columns Before KNN
    # knnImpute requires scaling. Scaling fails if StdDev is 0.
    nzv_metrics <- caret::nearZeroVar(X_train, saveMetrics = TRUE)
    zero_var_cols <- rownames(nzv_metrics[nzv_metrics$zeroVar == TRUE, ])

    if (length(zero_var_cols) > 0) {
      cat(sprintf(
        "\n[Auto-Fix] Dropping %d zero-variance cols before KNN (e.g., %s...)\n",
        length(zero_var_cols),
        zero_var_cols[1]
      ))
      # Remove from both Train and Test so columns match
      X_train <- X_train[, !(names(X_train) %in% zero_var_cols)]
      X_test <- X_test[, !(names(X_test) %in% zero_var_cols)]
    }

    # Check if we have any columns left
    if (ncol(X_train) == 0) {
      warning(
        "All columns had zero variance and were dropped. Returning empty datasets."
      )
      return(list(train = X_train, test = X_test))
    }

    # Now it is safe to run preProcess
    pp <- caret::preProcess(
      X_train,
      method = c("center", "scale", "knnImpute"),
      k = k
    )

    return(list(train = predict(pp, X_train), test = predict(pp, X_test)))
  }

  stop("Unsupported imputer in train-apply mode: ", imp_name)

  # ==============================================================================
  # [FINAL SAFETY CHECK] Force Test Columns to Match Train Columns Exactly
  # ==============================================================================
  # This prevents crashes in Random Forest due to column mismatch

  # 1. Identify columns in Train
  train_cols <- names(X_train)

  # 2. Subset Test to only keep columns that exist in Train
  #    (Drops any extra "junk" columns in Test that aren't in Train)
  X_test_safe <- X_test[, names(X_test) %in% train_cols, drop = FALSE]

  # 3. Check if Test is missing any columns that Train has
  missing_cols <- setdiff(train_cols, names(X_test_safe))

  # 4. If Test is missing columns, add them and fill with 0 (or Mean)
  #    (This happens if a compound was in Train but completely absent in Test)
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      # Fill with 0 (safe fallback) or the mean of the training column
      # Here we use 0 to be safe (assuming missing = not detected)
      X_test_safe[[col]] <- 0
    }
  }

  # 5. Re-order Test columns to match Train order exactly
  X_test_safe <- X_test_safe[, train_cols, drop = FALSE]

  # ==============================================================================

  # Return the aligned datasets
  return(list(train = X_train, test = X_test_safe))
}

no_normalization <- function(x) x

normalize_train_apply <- function(X_train, X_test, norm_name) {
  if (norm_name == "No_normalization") {
    return(list(train = X_train, test = X_test))
  }

  if (norm_name == "TSN_Normalization") {
    tsn <- function(X) {
      rs <- rowSums(X, na.rm = TRUE)
      rs[rs == 0] <- 1
      as.data.frame(sweep(X, 1, rs, "/"))
    }
    return(list(train = tsn(X_train), test = tsn(X_test)))
  }

  if (norm_name == "Min_Max_Scaler") {
    mins <- apply(X_train, 2, min, na.rm = TRUE)
    maxs <- apply(X_train, 2, max, na.rm = TRUE)
    rng <- maxs - mins
    rng[rng == 0] <- 1

    scale_mm <- function(X) {
      Xmat <- as.matrix(X)
      Xmat <- sweep(Xmat, 2, mins, "-")
      Xmat <- sweep(Xmat, 2, rng, "/")
      as.data.frame(Xmat)
    }
    return(list(train = scale_mm(X_train), test = scale_mm(X_test)))
  }

  if (norm_name == "Z_Score_Normalization") {
    mu <- colMeans(X_train, na.rm = TRUE)
    sdv <- apply(X_train, 2, sd, na.rm = TRUE)
    sdv[sdv == 0] <- 1

    zfun <- function(X) {
      Xmat <- as.matrix(X)
      Xmat <- sweep(Xmat, 2, mu, "-")
      Xmat <- sweep(Xmat, 2, sdv, "/")
      as.data.frame(Xmat)
    }
    return(list(train = zfun(X_train), test = zfun(X_test)))
  }

  stop("Unsupported normalizer in train-apply mode: ", norm_name)
}

impute_norm_train_apply <- function(X_train, X_test, imp_name, norm_name) {
  imp <- impute_train_apply(X_train, X_test, imp_name)
  normalize_train_apply(imp$train, imp$test, norm_name)
}

no_normalization <- function(x) x

###############################
##  4. Cluster resolution  ##
###############################

## Calculate cluster resolution using the first two PCA components
calculate_cluster_resolution <- function(data1, data2) {
  centroid1 <- colMeans(data1)
  centroid2 <- colMeans(data2)

  # Computes the Euclidean distance With Dim.1 as x-axis and Dim.2 as y-axis
  centroid_distance <- sqrt(
    (centroid1[1][[1]] - centroid2[1][[1]])^2 +
      (centroid1[2][[1]] - centroid2[2][[1]])^2
  )

  avg_distance1 <- mean(sqrt(rowSums(
    (as.matrix(data1) -
      matrix(centroid1, nrow = nrow(data1), ncol = ncol(data1), byrow = TRUE))^2
  )))
  avg_distance2 <- mean(sqrt(rowSums(
    (as.matrix(data2) -
      matrix(centroid2, nrow = nrow(data2), ncol = ncol(data2), byrow = TRUE))^2
  )))

  resolution <- centroid_distance / (avg_distance1 + avg_distance2)
  return(resolution)
}

## Compute average cluster resolution over all pairs of groups
calculate_average_cluster_resolution <- function(X, group_vector) {
  pca_model <- res.pca <- FactoMineR::PCA(
    X,
    scale.unit = FALSE,
    graph = FALSE
  )
  scores <- as.data.frame(pca_model$ind$coord[, 1:2])

  group <- as.factor(group_vector)
  group_levels <- levels(group)

  pair_resolutions <- c()

  for (i in 1:(length(group_levels) - 1)) {
    for (j in (i + 1):length(group_levels)) {
      data1 <- scores[group == group_levels[i], , drop = FALSE]
      data2 <- scores[group == group_levels[j], , drop = FALSE]

      if (nrow(data1) > 1 && nrow(data2) > 1) {
        res <- calculate_cluster_resolution(data1, data2)
        pair_resolutions <- c(pair_resolutions, res)
      }
    }
  }

  avg_resolution <- mean(pair_resolutions, na.rm = TRUE)
  return(avg_resolution)
}

mcc_multiclass <- function(true, pred) {
  true <- factor(true)
  pred <- factor(pred, levels = levels(true))
  cm <- table(true, pred)

  # Gorodkin multi-class MCC
  t_sum <- rowSums(cm)
  p_sum <- colSums(cm)
  n <- sum(cm)
  c <- sum(diag(cm))

  sum_pk_tk <- sum(p_sum * t_sum)
  sum_pk2 <- sum(p_sum^2)
  sum_tk2 <- sum(t_sum^2)

  numerator <- (c * n) - sum_pk_tk
  denom <- sqrt((n^2 - sum_pk2) * (n^2 - sum_tk2))

  if (denom == 0) {
    return(NA_real_)
  }
  numerator / denom
}

evaluate_combo_oob <- function(
  df_train,
  type_col,
  remove_cols,
  imp_name,
  norm_name,
  ntree = 500,
  seed = 123,
  num_threads = 1
) {
  X <- df_train %>%
    dplyr::select(-all_of(c(remove_cols, type_col))) %>%
    as.data.frame()
  y <- df_train[[type_col]]

  keep_cols <- names(X)[!sapply(X, function(v) all(is.na(v)))]
  X <- X[, keep_cols, drop = FALSE]

  if (ncol(X) < 1) {
    return(list(OOB_MCC = NA_real_))
  }
  if (length(unique(y)) < 2) {
    return(list(OOB_MCC = NA_real_))
  }

  pp <- tryCatch(
    impute_norm_train_apply(X, X, imp_name, norm_name),
    error = function(e) NULL
  )
  if (is.null(pp)) {
    return(list(OOB_MCC = NA_real_))
  }

  Xp <- pp$train
  if (any(is.na(Xp)) || any(!is.finite(as.matrix(Xp)))) {
    return(list(OOB_MCC = NA_real_))
  }

  y_fac <- factor(y)
  p <- ncol(Xp)
  mtry <- max(1, floor(sqrt(p)))

  train_df <- data.frame(.Class = y_fac, Xp, check.names = FALSE)

  set.seed(seed)
  rf <- ranger::ranger(
    dependent.variable.name = ".Class",
    data = train_df,
    num.trees = ntree,
    mtry = mtry,
    splitrule = "gini",
    importance = "impurity",
    oob.error = TRUE,
    probability = FALSE,
    num.threads = num_threads
  )

  oob_pred <- rf$predictions
  if (is.null(oob_pred)) {
    oob_pred <- predict(rf, data = train_df)$predictions
  }

  oob_pred <- factor(oob_pred, levels = levels(y_fac))
  mcc <- mcc_multiclass(y_fac, oob_pred)

  list(OOB_MCC = mcc)
}

##################################################################
##   6. Iterate Over Datasets, Imputation, and Normalization    ##
##################################################################
imputation_methods <- c(
  "mean_imputation",
  "median_imputation",
  "zero_imputation",
  "half_min_imputation",
  "random_forest_imputation",
  "knn_imputation"
)

normalization_methods <- c(
  "No_normalization",
  "TSN_Normalization",
  "Min_Max_Scaler",
  "Z_Score_Normalization"
)

find_best_impute_normalize <- function(
  df_train,
  type_col,
  remove_cols,
  ntree_screen = 500
) {
  X <- df_train %>%
    dplyr::select(-all_of(c(remove_cols, type_col))) %>%
    as.data.frame()
  y <- df_train[[type_col]]

  results_list <- list()

  for (imp_name in imputation_methods) {
    for (norm_name in normalization_methods) {
      oob_res <- evaluate_combo_oob(
        df_train,
        type_col,
        remove_cols,
        imp_name,
        norm_name,
        ntree = ntree_screen
      )

      # Unsupervised metrics on train only
      pp_self <- tryCatch(
        impute_norm_train_apply(X, X, imp_name, norm_name),
        error = function(e) NULL
      )
      if (is.null(pp_self)) {
        next
      }

      X_norm <- pp_self$train
      if (any(is.na(X_norm)) || any(!is.finite(as.matrix(X_norm)))) {
        next
      }

      cluster_res <- calculate_average_cluster_resolution(X_norm, y)

      orig_corr <- tryCatch(
        cor(X, use = "pairwise.complete.obs"),
        error = function(e) NULL
      )
      imp_corr <- tryCatch(
        cor(X_norm, use = "pairwise.complete.obs"),
        error = function(e) NULL
      )

      corr_score <- NA_real_
      if (
        !is.null(orig_corr) &&
          !is.null(imp_corr) &&
          all(dim(orig_corr) == dim(imp_corr)) &&
          ncol(orig_corr) > 1
      ) {
        corr_score <- suppressWarnings(
          cor(
            orig_corr[upper.tri(orig_corr)],
            imp_corr[upper.tri(imp_corr)],
            use = "complete.obs"
          )
        )
      }

      # KS on train only
      observed_all <- c()
      imputed_all <- c()
      for (col in names(X)) {
        if (any(is.na(X[[col]])) && any(!is.na(X[[col]]))) {
          observed_all <- c(observed_all, X[[col]][!is.na(X[[col]])])
          imputed_all <- c(imputed_all, X_norm[[col]])
        }
      }

      results_list[[length(results_list) + 1]] <- data.frame(
        Imputation = imp_name,
        Normalization = norm_name,
        OOB_MCC = oob_res$OOB_MCC,
        ClusterRes = cluster_res,
        CorrScore = corr_score,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results_list) == 0) {
    return(list(results_table = data.frame(), best_combo = data.frame()))
  }

  df_results <- do.call(rbind, results_list)

  ord <- order(
    -df_results$OOB_MCC,
    -df_results$ClusterRes,
    -df_results$CorrScore,
  )
  df_results$combined_rank <- NA_integer_
  df_results$combined_rank[ord] <- seq_along(ord)

  best_row <- df_results[which.min(df_results$combined_rank), ]
  list(results_table = df_results, best_combo = best_row)
}

#####################################################
##  Data Fusion to create Multi-instrumental datasets
#####################################################
data_fusion <- function(df_input, source_feats, class_col = "Plastic_type") {
  tech_map <- df_input %>%
    dplyr::distinct(Feature, technique)

  all_techs <- unique(tech_map$technique)

  # --- Step B: Pivot to Wide Format ---
  df_wide <- df_input %>%
    dplyr::select(
      File,
      Source,
      Polymer,
      Feature,
      technique,
      Values,
      all_of(class_col)
    ) %>%
    dplyr::group_by(
      File,
      Source,
      Polymer,
      technique,
      Feature,
      across(all_of(class_col))
    ) %>%
    dplyr::summarise(Values = mean(Values, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Feature, values_from = Values) %>%
    tibble::column_to_rownames(var = "File") %>%
    dplyr::relocate(
      any_of(c(class_col, "technique", "Source", "Polymer")),
      .before = 1
    )

  # Metadata columns
  meta_cols <- c("Plastic_type", "technique", "Source", "Polymer")

  # --- Step C: Fast Zero Replacement (Matrix Method) ---
  if (length(all_techs) > 1) {
    cat(
      "Detected multiple techniques (",
      paste(all_techs, collapse = ", "),
      "). Applying Fast Zero Replacement...\n"
    )

    feat_cols <- setdiff(names(df_wide), meta_cols)
    mat <- as.matrix(df_wide[, feat_cols])
    row_techs <- df_wide$technique

    for (tech in all_techs) {
      tech_feats <- tech_map$Feature[tech_map$technique == tech]
      target_cols <- intersect(colnames(mat), tech_feats)

      if (length(target_cols) > 0) {
        target_rows <- which(row_techs != tech)
        if (length(target_rows) > 0) {
          sub_mat <- mat[target_rows, target_cols, drop = FALSE]
          sub_mat[is.na(sub_mat)] <- 0
          mat[target_rows, target_cols] <- sub_mat
        }
      }
    }
    df_wide[, feat_cols] <- mat
  } else {
    cat("Detected single technique. Skipping Zero Replacement.\n")
  }

  final_feat_cols <- c()
  for (source in source_feats) {
    f <- setdiff(names(source), meta_cols)
    final_feat_cols <- c(final_feat_cols, f)
  }
  output <- df_wide[, c(meta_cols, unique(final_feat_cols))]
  return(output)
}

###########################
### Helper functions for RF
###########################
pairwise_significance_tests <- function(input_df, group_col, start_col_index) {
  # For each feature from start_col_index onward,
  #   for each pair of group categories,
  #     check if both groups have >=3 data points
  #     do Shapiro test for normality
  #     if not normal => Mann-Whitney (wilcox.test)
  #     else => t.test
  # Collect p-values and do multiple corrections

  df_results <- data.frame(
    Feature = character(),
    comparison_pair = character(),
    pval = numeric()
  )

  # get unique group combos
  group_pairs <- utils::combn(unique(as.character(input_df[[group_col]])), 2)

  for (feature_col in start_col_index:ncol(input_df)) {
    for (col in 1:ncol(group_pairs)) {
      p_1 <- group_pairs[1, col]
      p_2 <- group_pairs[2, col]

      vec1 <- as.numeric(unlist(input_df[
        which(input_df[[group_col]] == p_1),
        feature_col
      ]))
      vec2 <- as.numeric(unlist(input_df[
        which(input_df[[group_col]] == p_2),
        feature_col
      ]))

      # ensure at least 3 non-NA data points
      if (sum(!is.na(vec1)) >= 3 && sum(!is.na(vec2)) >= 3) {
        # skip if all values are the same
        # if(all(vec1 == vec1[1]) || all(vec2 == vec2[1])){
        #   next
        # }

        # normality tests
        p_shapiro1 <- tryCatch(shapiro.test(vec1)$p.value, error = function(e) {
          NA
        })
        p_shapiro2 <- tryCatch(shapiro.test(vec2)$p.value, error = function(e) {
          NA
        })

        # decide parametric vs non-parametric
        if (is.na(p_shapiro1) || is.na(p_shapiro2)) {
          # if for some reason test failed, skip
          next
        }

        if (p_shapiro1 < 0.05 || p_shapiro2 < 0.05) {
          # mann-whitney
          pval_test <- tryCatch(
            wilcox.test(vec1, vec2)$p.value,
            error = function(e) NA
          )
        } else {
          # t-test
          pval_test <- tryCatch(
            t.test(vec1, vec2, var.equal = FALSE)$p.value,
            error = function(e) NA
          )
        }

        if (!is.na(pval_test)) {
          # append to df
          df_results <- rbind(
            df_results,
            data.frame(
              Feature = colnames(input_df)[feature_col],
              Category_1 = p_1,
              Category_2 = p_2,
              pval = as.numeric(pval_test),
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  # multiple corrections
  pvalue_correction_methods <- "bonferroni" #, c("holm", "hochberg", "hommel", "BH", "BY")

  # For convenience, store all results in a list
  # results_list <- list()
  # j <- 1
  # for(method in pvalue_correction_methods){
  temp <- df_results
  temp$adjusted_pvalue <- p.adjust(
    temp$pval,
    method = pvalue_correction_methods
  )
  significant <- temp %>%
    dplyr::filter(adjusted_pvalue < 0.05) %>%
    arrange(adjusted_pvalue)
  # results_list[[j]] <- significant

  if (nrow(significant) == 0) {
    cat(paste0(
      pvalue_correction_methods,
      " did not result in any significant compound\n"
    ))
  } else {
    # cat(paste0("Significant results for method: ", method, "\n"))
    # print(significant)
    cat(paste0(
      pvalue_correction_methods,
      " has ",
      length(unique(significant$Feature)),
      " significant features\n"
    ))
    # j <- j + 1
  }
  # }

  return(list(
    uncorrected = df_results,
    corrected = unique(significant$Feature) # results_list
  ))
}

#─────────────────────────────────────────────────────────────────────────────
# Helper: retune mtry (based on subset size) + ntree for any train/test split
#─────────────────────────────────────────────────────────────────────────────
tune_rf_subset <- function(
  X_train,
  y_train,
  X_test,
  y_test,
  ntree_candidates,
  metric = "Accuracy",
  parallel = TRUE,
  n_cores = NULL,
  num_threads = 1,
  min_node_size = 5 # <--- NEW PARAMETER: Default 5 for Regularization
) {
  y_train <- as.factor(y_train)
  y_test <- as.factor(y_test)

  k_inner <- min(table(y_train))

  p <- ncol(X_train)
  mtry_vals <- unique(pmax(1, c(floor(sqrt(p)), floor(log2(p)))))

  ranger_predict_probs <- function(rf, newdata) {
    pr <- predict(rf, data = newdata)$predictions
    if (is.matrix(pr)) {
      prob_mat <- pr
      pred_chr <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
      pred_fac <- factor(pred_chr, levels = colnames(prob_mat))
      return(list(preds = pred_fac, probs = prob_mat))
    }
    pred_fac <- factor(pr, levels = levels(y_train))
    return(list(preds = pred_fac, probs = NULL))
  }

  # Choose cores
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  # ---- Branch A: enough data for CV ----
  if (k_inner > 1) {
    ctrl <- caret::trainControl(
      method = "cv",
      number = k_inner,
      classProbs = TRUE,
      summaryFunction = defaultSummary,
      savePredictions = "final",
      allowParallel = parallel
    )

    rf_grid <- expand.grid(
      mtry = mtry_vals,
      splitrule = "gini",
      min.node.size = 1
    )

    best_mod <- NULL
    best_val <- -Inf
    best_trees <- NA

    # If you parallelize caret resamples, keep ranger internal threads modest
    ranger_threads_for_caret <- if (parallel) 1 else num_threads

    for (nt in ntree_candidates) {
      set.seed(123)

      # ERROR FIX: 'probability = TRUE' removed below because caret handles it automatically
      tmp <- caret::train(
        x = X_train,
        y = y_train,
        method = "ranger",
        trControl = ctrl,
        tuneGrid = rf_grid,
        metric = metric,
        num.trees = nt,
        importance = "impurity",
        # probability = TRUE,  <-- REMOVED THIS LINE
        num.threads = ranger_threads_for_caret
      )

      this_val <- suppressWarnings(max(tmp$results[[metric]], na.rm = TRUE))
      if (is.finite(this_val) && this_val > best_val) {
        best_val <- this_val
        best_mod <- tmp
        best_trees <- nt
      }
    }

    caret_model <- best_mod
    final_rf <- caret_model$finalModel
  } else {
    # ---- Branch B: too few per class -> OOB tuning grid ----

    train_df <- data.frame(.Class = y_train, X_train, check.names = FALSE)
    grid <- expand.grid(ntree = ntree_candidates, 
                        mtry = mtry_vals,
                        min.node.size = min_node_size) # <--- Added to OOB grid

    # Parallel backend
    created_cluster <- FALSE
    if (parallel) {
      if (!foreach::getDoParRegistered() || foreach::getDoParWorkers() < 2) {
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
        created_cluster <- TRUE
      }
    }

    on.exit(
      {
        if (parallel && created_cluster) {
          try(parallel::stopCluster(cl), silent = TRUE)
          foreach::registerDoSEQ()
        }
      },
      add = TRUE
    )

    eval_one <- function(nt, m, ns) {
      set.seed(123)
      rf_mod <- ranger::ranger(
        dependent.variable.name = ".Class",
        data = train_df,
        num.trees = nt,
        mtry = m,
        min.node.size = ns,      # <--- Applies Regularization
        splitrule = "gini",
        importance = "impurity",
        probability = FALSE,
        oob.error = TRUE,
        num.threads = 1
      )
      data.frame(
        ntree = nt, 
        mtry = m, 
        min.node.size = ns,
        oob_err = rf_mod$prediction.error
      )
    }
    
    if (parallel) {
      res_grid <- foreach::foreach(
        i = seq_len(nrow(grid)),
        .combine = rbind,
        .packages = "ranger"
      ) %dopar% {
        eval_one(grid$ntree[i], grid$mtry[i], grid$min.node.size[i])
      }
    } else {
      res_grid <- do.call(
        rbind,
        lapply(seq_len(nrow(grid)), function(i) {
          eval_one(grid$ntree[i], grid$mtry[i], grid$min.node.size[i])
        })
      )
    }
    
    res_grid <- res_grid[is.finite(res_grid$oob_err), , drop = FALSE]
    
    if (nrow(res_grid) == 0) {
      return(list(
        model = NULL,
        ntree = NA,
        preds = factor(NA),
        p_m = matrix(NA, nrow = 0, ncol = 0),
        eval_metrics = list(),
        Accuracy = NA,
        Kappa = NA,
        AUC = NA,
        F1_Macro = NA,
        F1_Weighted = NA,
        MCC = NA
      ))
    }

    best_row <- res_grid[which.min(res_grid$oob_err), ]
    best_trees <- best_row$ntree
    best_mtry <- best_row$mtry
    best_node_size <- best_row$min.node.size # Capture best regularizer

    set.seed(123)
    best_rf <- ranger::ranger(
      dependent.variable.name = ".Class",
      data = train_df,
      num.trees = best_trees,
      mtry = best_mtry,
      min.node.size = best_node_size, # <--- Use best node size
      splitrule = "gini",
      importance = "impurity",
      probability = TRUE,
      oob.error = TRUE,
      num.threads = num_threads
    )

    caret_model <- list(
      finalModel = best_rf,
      bestTune = data.frame(
        mtry = best_mtry,
        splitrule = "gini",
        min.node.size = best_node_size
      )
    )
    final_rf <- best_rf
  }

  # ---- Predict & score ----
  pred_out <- ranger_predict_probs(final_rf, X_test)
  preds <- pred_out$preds
  prob_mat <- pred_out$probs

  if (is.null(prob_mat)) {
    cls <- levels(y_train)
    prob_mat <- matrix(0, nrow = length(preds), ncol = length(cls))
    colnames(prob_mat) <- cls
    rownames(prob_mat) <- rownames(X_test)
    prob_mat[cbind(seq_along(preds), match(as.character(preds), cls))] <- 1
  }

  ev <- eval_metrics(
    true_labels = y_test,
    prob_matrix = prob_mat,
    pred_labels = preds
  )

  list(
    model = caret_model,
    ntree = best_trees,
    preds = preds,
    p_m = prob_mat,
    eval_metrics = ev,
    Accuracy = ev$Accuracy,
    Kappa = ev$Kappa,
    AUC = ev$AUC,
    F1_Macro = ev$F1_Macro,
    F1_Weighted = ev$F1_Weighted,
    MCC = ev$MCC_Multiclass
  )
}

# Recursive Feature Addition -------------
recursive_feature_addition <- function(
  X_train,
  y_train,
  X_test,
  y_test,
  ntree_candidates,
  metric = "Accuracy",
  parallel = TRUE,
  n_cores = NULL,
  num_threads_model = 1
) {
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  # 1) Fit full model to get initial importance
  full_res <- tune_rf_subset(
    X_train,
    y_train,
    X_test,
    y_test,
    ntree_candidates,
    metric,
    parallel = parallel,
    n_cores = n_cores,
    num_threads = num_threads_model
  )

  imp_vec <- full_res$model$finalModel$variable.importance
  init_feat <- names(imp_vec)[which.max(imp_vec)]

  best_feats <- init_feat
  remaining <- setdiff(colnames(X_train), best_feats)

  # 2) Baseline performance
  base_res <- tune_rf_subset(
    X_train[, best_feats, drop = FALSE],
    y_train,
    X_test[, best_feats, drop = FALSE],
    y_test,
    ntree_candidates,
    metric,
    parallel = FALSE, # avoid nested parallel
    num_threads = num_threads_model
  )
  best_mcc <- base_res$MCC

  # Parallel backend for feature-try loop
  created_cluster <- FALSE
  if (parallel) {
    if (!foreach::getDoParRegistered() || foreach::getDoParWorkers() < 2) {
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      created_cluster <- TRUE
    }
  }

  on.exit(
    {
      if (parallel && created_cluster) {
        try(parallel::stopCluster(cl), silent = TRUE)
        foreach::registerDoSEQ()
      }
    },
    add = TRUE
  )

  # 3) Iterative greedy addition
  while (length(remaining) > 0) {
    # Evaluate all candidates in parallel processing to optimize computational time
    if (parallel) {
      cand_res <- foreach::foreach(
        feat = remaining,
        .combine = rbind,
        # FIX IS HERE: Added .export
        .export = c("tune_rf_subset", "eval_metrics"),
        .packages = c("ranger", "caret", "pROC", "dplyr")
      ) %dopar%
        {
          trial_feats <- c(best_feats, feat)

          # Now the worker node can find this function:
          res <- tune_rf_subset(
            X_train[, trial_feats, drop = FALSE],
            y_train,
            X_test[, trial_feats, drop = FALSE],
            y_test,
            ntree_candidates,
            metric,
            parallel = FALSE, # critical
            num_threads = 1 # critical inside workers
          )
          data.frame(feat = feat, MCC = res$MCC, stringsAsFactors = FALSE)
        }
    } else {
      cand_res <- do.call(
        rbind,
        lapply(remaining, function(feat) {
          trial_feats <- c(best_feats, feat)
          res <- tune_rf_subset(
            X_train[, trial_feats, drop = FALSE],
            y_train,
            X_test[, trial_feats, drop = FALSE],
            y_test,
            ntree_candidates,
            metric,
            parallel = FALSE,
            num_threads = num_threads_model
          )
          data.frame(feat = feat, MCC = res$MCC, stringsAsFactors = FALSE)
        })
      )
    }

    cand_res <- cand_res[is.finite(cand_res$MCC), , drop = FALSE]
    if (nrow(cand_res) == 0) {
      break
    }

    # Pick best improvement
    best_cand <- cand_res[which.max(cand_res$MCC), , drop = FALSE]

    if (best_cand$MCC > best_mcc) {
      best_mcc <- best_cand$MCC
      best_feats <- c(best_feats, best_cand$feat)
      remaining <- setdiff(remaining, best_cand$feat)
    } else {
      break
    }
  }

  # 4) Final model
  final_res <- tune_rf_subset(
    X_train[, best_feats, drop = FALSE],
    y_train,
    X_test[, best_feats, drop = FALSE],
    y_test,
    ntree_candidates,
    metric,
    parallel = FALSE,
    num_threads = num_threads_model
  )

  list(
    best_mcc = final_res$MCC,
    best_acc = final_res$Accuracy,
    best_features_rf = best_feats,
    final_model = final_res$model,
    prob_matrix = final_res$p_m,
    predictions = final_res$preds,
    eval_metrics = final_res$eval_metrics
  )
}

#─────────────────────────────────────────────────────────────────────────────
library(caret)
library(pROC)

eval_metrics <- function(true_labels, prob_matrix, pred_labels) {
  # 0) Ensure factor levels align
  classes <- union(levels(true_labels), levels(pred_labels))
  y_true <- factor(true_labels, levels = classes)
  y_pred <- factor(pred_labels, levels = classes)

  # 1) Build the confusion‐matrix table
  cm <- table(Actual = y_true, Predicted = y_pred)
  N <- sum(cm) # total samples
  diag_cm <- diag(cm)

  # 2) Accuracy
  accuracy <- sum(diag_cm) / N

  # 3) Cohen’s Kappa (via caret)
  kappa <- tryCatch(
    {
      cm_obj <- confusionMatrix(y_pred, y_true)
      as.numeric(cm_obj$overall["Kappa"])
    },
    error = function(e) NA
  )

  # 4) One‐vs‐all multiclass AUC
  auc <- tryCatch(
    {
      roc_obj <- pROC::multiclass.roc(
        response = true_labels,
        predictor = prob_matrix
      )
      as.numeric(roc_obj$auc)
    },
    error = function(e) NA
  )

  # 5) Per‐class precision & recall
  col_sums <- colSums(cm) # predicted totals for each class
  row_sums <- rowSums(cm) # actual   totals for each class

  precision_per_class <- ifelse(col_sums > 0, diag_cm / col_sums, 0) # no preds → p=0
  recall_per_class <- ifelse(row_sums > 0, diag_cm / row_sums, 0) # no actual → r=0

  # 6) Per‐class F1 & weighted F1
  f1_per_class <- ifelse(
    (precision_per_class + recall_per_class) > 0,
    2 *
      precision_per_class *
      recall_per_class /
      (precision_per_class + recall_per_class),
    0
  )
  support <- row_sums # how many true examples of each class
  f1_weighted <- sum(f1_per_class * support) / sum(support)

  # 7) Multiclass MCC (Gorodkin’s formula)
  t_i <- row_sums
  p_j <- col_sums
  c_val <- sum(diag_cm) # sum of true positives
  num <- N * c_val - sum(t_i * p_j)
  den <- sqrt((N^2 - sum(t_i^2)) * (N^2 - sum(p_j^2)))
  mcc <- ifelse(den > 0, num / den, 0)

  return(list(
    Accuracy = accuracy,
    Kappa = kappa,
    AUC = auc,
    F1_Weighted = f1_weighted,
    MCC_Multiclass = mcc
  ))
}


####################################################
##########        CONFUSION MATRIX       ###########
####################################################
conf_mat_plot <- function(y_test, preds, conf_mat_title, accuracy) {
  # 1) Compute raw counts table: rows = treated labels, cols = predicted pristine labels
  all_classes <- unique(c(levels(preds), levels(y_test)))
  y_f <- factor(y_test, levels = all_classes)
  p_f <- factor(preds, levels = all_classes)

  cm_tab <- table(Actual = y_f, Predicted = p_f)

  # 2) Turn it into a data.frame for ggplot
  cm_df <- as.data.frame(cm_tab, stringsAsFactors = FALSE) %>%
    dplyr::rename(Freq = Freq)

  # 3) Add total-per-row and percent-per-cell
  cm_df$Total <- NA
  cm_df$Percent <- NA
  for (label in unique(cm_df$Actual)) {
    cm_df[which(cm_df$Actual == label), ]$Total <- sum(
      cm_df[which(cm_df$Actual == label), ]$Freq
    )
  }

  cm_df <- cm_df %>%
    dplyr::mutate(
      Percent = round(Freq / Total * 100, 2)
    ) %>%
    dplyr::mutate(Percent = ifelse(is.nan(Percent), 0, Percent))

  cm_df <- cm_df %>%
    dplyr::mutate(
      Label = ifelse(Predicted == Actual, paste0(round(Percent, 1), "%"), "")
    )

  # 5) Plot
  conf_mat <- ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Percent)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") + #, limits = c(0,1)) +
    geom_text(aes(label = Percent), size = 6, show.legend = FALSE) +
    # explicitly drop any size guide, keep only the fill (Recall) legend
    guides(size = "none") +
    labs(
      title = paste0("Confusion Matrix of ENV → SB Mapping | ", conf_mat_title),
      subtitle = paste0(sprintf(
        "Overall linking accuracy: %.2f%%\n",
        accuracy
      )),
      x = "Predicted Class (Store-Bought)",
      y = "Actual Class (Environmental)",
      fill = "Percentage of\nClassification (%)" # "Recall %"
    ) +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  return(list(conf_mat = conf_mat, cm_df = cm_df))
}

#################
### RF classifier
#################
run_rf_analysis_manuscript1 <- function(
  data,
  type_col,
  remove_cols,
  train_proportion,
  use_source_split = FALSE,
  use_store_vs_environmental_split = FALSE,
  use_polymer = FALSE,
  do_pairwise_test = FALSE,
  do_top_importance_selection = FALSE,
  top_importance_counts = c(100, 50, 25, 10),
  data_name,
  seed = 123,
  ntree_candidates = c(100, 500, 1000, 2500),
  metric = "Accuracy",
  do_rfe = FALSE,
  do_rfa = FALSE,
  min_sample_number,
  do_impute_norm_screen = TRUE,
  group_for_significance,
  K_outer_max = 5,
  ntree_screen = 500
) {
  set.seed(seed)
  use_crayon <- requireNamespace("crayon", quietly = TRUE)
  style_info <- function(text) {
    if (use_crayon) {
      return(crayon::cyan(text))
    }
    text
  }
  style_warn <- function(text) {
    if (use_crayon) {
      return(crayon::yellow(text))
    }
    text
  }
  style_step <- function(text) {
    if (use_crayon) {
      return(crayon::bold(crayon::blue(text)))
    }
    text
  }

  # ---- 1) Clean classes ----
  data[[type_col]] <- as.factor(data[[type_col]])
  levels(data[[type_col]]) <- make.names(
    levels(data[[type_col]]),
    unique = TRUE
  )

  class_counts <- table(data[[type_col]])
  keep_classes <- names(class_counts[class_counts >= min_sample_number])
  data <- data[data[[type_col]] %in% keep_classes, ]
  data[[type_col]] <- droplevels(data[[type_col]])

  rownames(data) <- NULL

  data <- data %>%
    mutate(sample_id = paste0(Plastic_type, "_", technique)) %>%
    mutate(
      sample_id = ave(sample_id, sample_id, FUN = function(x) {
        if (length(x) > 1) paste0(x, "_rep", seq_along(x)) else x
      })
    ) %>%
    tibble::column_to_rownames(var = "sample_id")

  # ---- 2) Build OUTER folds first ----
  outer_folds <- list()

  if (use_store_vs_environmental_split) {
    store_bought_set <- data[data$Source == "Store-Bought", ]
    environmental_set <- data[data$Source == "Environmental", ]

    # Outer folds on STORE-Bought only (training stability)
    sb_y <- store_bought_set[[type_col]]
    min_sb <- min(table(as.character(sb_y)))
    K <- max(2, min(K_outer_max, min_sb))

    sb_test_folds <- caret::createFolds(sb_y, k = K, returnTrain = FALSE)

    for (f in seq_len(K)) {
      test_idx_sb <- sb_test_folds[[f]]
      train_idx_sb <- setdiff(seq_len(nrow(store_bought_set)), test_idx_sb)

      outer_folds[[f]] <- list(
        train_data = store_bought_set[train_idx_sb, , drop = FALSE],
        test_data = environmental_set # fixed test
      )
    }
  } else if (use_source_split) {
    store_bought_set <- data[data$Source == "Store-Bought", ]
    environmental_set <- data[data$Source == "Environmental", ]

    env_y <- environmental_set[[type_col]]
    min_env <- min(table(as.character(env_y)))
    K <- max(2, min(K_outer_max, min_env))

    env_test_folds <- caret::createFolds(env_y, k = K, returnTrain = FALSE)

    for (f in seq_len(K)) {
      test_idx_env <- env_test_folds[[f]]
      train_idx_env <- setdiff(seq_len(nrow(environmental_set)), test_idx_env)

      outer_folds[[f]] <- list(
        train_data = rbind(
          store_bought_set,
          environmental_set[train_idx_env, , drop = FALSE]
        ),
        test_data = environmental_set[test_idx_env, , drop = FALSE]
      )
    }
  } else {
    y <- data[[type_col]]
    min_all <- min(table(y))
    K <- max(2, min(K_outer_max, min_all))
    cat(
      style_info(paste0("min_all is: ", min_all, "_Number of fold: ", K)),
      "\n"
    )

    test_folds <- caret::createFolds(y, k = K, returnTrain = FALSE)

    for (f in seq_len(K)) {
      test_idx <- test_folds[[f]]
      train_idx <- setdiff(seq_len(nrow(data)), test_idx)

      outer_folds[[f]] <- list(
        train_data = data[train_idx, , drop = FALSE],
        test_data = data[test_idx, , drop = FALSE]
      )
    }
  }

  # ---- 3) Fit/evaluate per OUTER fold ----
  final_final_res_list <- vector("list", length = length(outer_folds))

  for (f in seq_along(outer_folds)) {
    train_data <- outer_folds[[f]]$train_data
    test_data <- outer_folds[[f]]$test_data

    train_data[[type_col]] <- droplevels(train_data[[type_col]])
    test_data[[type_col]] <- droplevels(test_data[[type_col]])

    # --- 3a) Impute+Norm screening on TRAIN only (OOB MCC) ---
    if (do_impute_norm_screen) {
      t0 <- Sys.time()

      best_res_imp_norm <- find_best_impute_normalize(
        df_train = train_data,
        type_col = type_col,
        remove_cols = remove_cols,
        ntree_screen = ntree_screen
      )

      best_imp <- best_res_imp_norm$best_combo$Imputation
      best_norm <- best_res_imp_norm$best_combo$Normalization

      numeric_cols <- names(which(sapply(train_data, is.numeric)))
      numeric_cols <- setdiff(numeric_cols, c(remove_cols, type_col))

      X_train_orig <- train_data[, numeric_cols, drop = FALSE]
      X_test_orig <- test_data[, numeric_cols, drop = FALSE]

      pp_final <- impute_norm_train_apply(
        X_train_orig,
        X_test_orig,
        best_imp,
        best_norm
      )

      train_data[, numeric_cols] <- pp_final$train
      test_data[, numeric_cols] <- pp_final$test

      t1 <- Sys.time()
      time_imp_norm <- as.numeric(difftime(t1, t0, units = "secs"))
    } else {
      best_res_imp_norm <- NULL
      best_imp <- NA_character_
      best_norm <- NA_character_
      time_imp_norm <- NA_real_
    }

    # --- 3c) Build X/y ---
    X_train <- train_data %>%
      dplyr::select(-c(all_of(remove_cols), all_of(type_col)))
    y_train <- train_data[[type_col]]
    X_test <- test_data %>%
      dplyr::select(-c(all_of(remove_cols), all_of(type_col)))
    y_test <- test_data[[type_col]]

    # --- 3d) Full features ---
    cat(style_step("### Train/test with Full original features ###"), "\n")
    start_time_full_features <- Sys.time()

    all_feats <- tune_rf_subset(
      X_train,
      y_train,
      X_test,
      y_test,
      ntree_candidates,
      metric
    )

    best_model_all_feats <- all_feats$model
    mcc_full <- all_feats$MCC
    eval_metrics_all_feats <- all_feats$eval_metrics
    # --- Generate Heatmap (Full Features) ---
    prob_matrix_all_feats <- all_feats$p_m
    conf_mat <- conf_mat_plot(
      y_test,
      all_feats$preds,
      conf_mat_title = paste0("All Features - ", data_name),
      accuracy = mcc_full
    )
    conf_mat_all_feats <- conf_mat$conf_mat
    cm_df_all_feats <- conf_mat$cm_df
    end_time_full_features <- Sys.time()
    time_all_feats <- as.numeric(difftime(
      end_time_full_features,
      start_time_full_features,
      units = "secs"
    ))

    # --- 3e) Pairwise features (if any) ---
    if (do_pairwise_test) {
      cat(style_step("### Train/test with Pair-wise features ###"), "\n")
      start_time_sig <- Sys.time()
      pairwise_res <- pairwise_significance_tests(
        input_df = train_data,
        group_col = group_for_significance,
        start_col_index = 6
      )
      sig_feats <- pairwise_res$corrected

      if (length(sig_feats) > 0) {
        # Case A: Features found - Run RF
        tmp_sig <- tune_rf_subset(
          X_train[, sig_feats, drop = FALSE],
          y_train,
          X_test[, sig_feats, drop = FALSE],
          y_test,
          ntree_candidates,
          metric
        )
        mcc_sig <- tmp_sig$MCC
        eval_metrics_sig <- tmp_sig$eval_metrics

        prob_matrix_sig <- tmp_sig$p_m
        conf_mat <- conf_mat_plot(
          y_test,
          tmp_sig$preds,
          conf_mat_title = paste0(
            "Pair-wise Significance-Based Features - ",
            data_name
          ),
          accuracy = mcc_sig
        )
        conf_mat_sig <- conf_mat$conf_mat
        cm_df_sig <- conf_mat$cm_df

        # Save predictions to a common variable
        preds_sig <- tmp_sig$preds
      } else {
        # Case B: No features found - Initialize empty placeholders
        mcc_sig <- NA_real_
        eval_metrics_sig <- NULL
        prob_matrix_sig <- NULL
        conf_mat_sig <- NULL
        cm_df_sig <- NULL
        preds_sig <- NULL # Crucial: prevents "object 'tmp_sig' not found" later
      }
      end_time_sig <- Sys.time()
      time_sig <- as.numeric(difftime(
        end_time_sig,
        start_time_sig,
        units = "secs"
      ))
    } else {
      # If do_pairwise_test is FALSE entirely
      pairwise_sig_features <- character(0)
      sig_feats <- intersect(pairwise_sig_features, colnames(X_train))
    }

    # --- 3f) RFA / Top-N, etc. (your existing blocks can follow) ---
    # Recursive Feature Addition
    if (do_rfa) {
      cat(style_step("### Train/test with RFA features ###"), "\n")
      start_time_rfa <- Sys.time()

      rfa_res <- recursive_feature_addition(
        X_train,
        y_train,
        X_test,
        y_test,
        ntree_candidates = ntree_candidates,
        metric = metric
      )

      rfa_selected_feats <- rfa_res$best_features_rf
      mcc_rfa <- rfa_res$best_mcc
      final_rf_rfa <- rfa_res$final_model$finalModel

      prob_matrix_rfa <- rfa_res$prob_matrix

      conf_mat <- conf_mat_plot(
        y_test,
        rfa_res$predictions,
        conf_mat_title = paste0("RFA Features - ", type_col),
        accuracy = mcc_rfa
      )
      conf_mat_rfa <- conf_mat$conf_mat
      cm_df_rfa <- conf_mat$cm_df

      eval_metrics_rfa <- rfa_res$eval_metrics

      end_time_rfa <- Sys.time()
      time_rfa <- as.numeric(difftime(
        end_time_rfa,
        start_time_rfa,
        units = "secs"
      ))
    }

    # 10) Classification with Top-N Features by Importance (if requested)
    if (do_top_importance_selection) {
      cat(
        style_step(
          "### Running Classification on Top-N Features by Importance ###"
        ),
        "\n"
      )

      imp_vec <- best_model_all_feats$finalModel$variable.importance

      rf_importances <- data.frame(
        Feature = names(imp_vec),
        Importance = as.numeric(imp_vec),
        stringsAsFactors = FALSE
      )

      rf_importances <- rf_importances[
        order(rf_importances$Importance, decreasing = TRUE),
      ]

      top_importance_results <- list()
      total_features <- nrow(rf_importances)

      for (top_n in top_importance_counts) {
        start_time_top_100_50_25_10 <- base::Sys.time()

        if (total_features < top_n) {
          selected_features <- rf_importances$Feature
          cat(
            style_warn(sprintf(
              "Only %d features available; using all.",
              total_features
            )),
            "\n"
          )
        } else {
          selected_features <- head(rf_importances$Feature, top_n)
        }

        tmp_top <- tune_rf_subset(
          X_train[, selected_features, drop = FALSE],
          y_train,
          X_test[, selected_features, drop = FALSE],
          y_test,
          ntree_candidates,
          metric
        )

        mcc_top <- tmp_top$MCC
        final_rf_top <- tmp_top$model$finalModel

        prob_matrix_top <- tmp_top$p_m
        conf_mat <- conf_mat_plot(
          y_test,
          tmp_top$preds,
          conf_mat_title = paste0(
            "Top ",
            length(selected_features),
            " Features - ",
            data_name
          ),
          accuracy = mcc_top
        )
        conf_mat_top <- conf_mat$conf_mat
        cm_df_top <- conf_mat$cm_df

        eval_metrics_top <- tmp_top$eval_metrics

        end_time_top_100_50_25_10 <- Sys.time()
        time_top_100_50_25_10 <- as.numeric(difftime(
          end_time_top_100_50_25_10,
          start_time_top_100_50_25_10,
          units = "secs"
        ))

        top_importance_results[[paste0("top_", top_n)]] <- list(
          best_model = final_rf_top,
          selected_features = selected_features,
          preds_for_buidling_confmat_top = tmp_top$preds,
          prob_matrix_top = prob_matrix_top,
          conf_mat_top = conf_mat_top,
          eval_metrics_top = eval_metrics_top,
          cm_df_top = cm_df_top,
          time_top_100_50_25_10 = time_top_100_50_25_10
        )
      }
    }

    # --- 3g) Package fold output ---
    result_list <- list(
      fold_id = f,
      final_imp_norm_train = train_data,
      final_imp_norm_test = test_data,
      best_model_all_feats = best_model_all_feats,
      final_rf_all_feats = best_model_all_feats$finalModel,
      all_features_mcc = mcc_full,
      eval_metrics_all_feats = eval_metrics_all_feats,
      prob_matrix = prob_matrix_all_feats,
      time_all_features = time_all_feats,
      conf_mat_all_feats = conf_mat_all_feats,
      cm_df_all_feats = cm_df_all_feats,
      preds_for_buidling_confmat_all = all_feats$preds
    )

    if (do_impute_norm_screen) {
      result_list$imp_norm_res_table <- best_res_imp_norm$results_table
      result_list$time_imp_norm <- time_imp_norm
      result_list$best_imputation <- best_imp
      result_list$best_normalization <- best_norm
    }

    if (do_pairwise_test) {
      result_list$mcc_sig <- mcc_sig
      result_list$time_sig <- time_sig
      result_list$prob_matrix_sig <- prob_matrix_sig
      result_list$eval_metrics_sig <- eval_metrics_sig
      result_list$sig_selected_feats <- sig_feats
      result_list$conf_mat_sig <- conf_mat_sig
      result_list$cm_df_sig <- cm_df_sig
      # CHANGED: Use the common variable 'preds_sig' defined above
      result_list$preds_for_buidling_confmat_sig <- preds_sig
    }

    if (do_rfa) {
      result_list$rfa_selected_feats <- rfa_selected_feats
      result_list$prob_matrix_rfa <- prob_matrix_rfa
      result_list$mcc_rfa <- mcc_rfa
      result_list$time_rfa <- time_rfa
      result_list$eval_metrics_rfa <- eval_metrics_rfa
      result_list$conf_mat_rfa <- conf_mat_rfa
      result_list$cm_df_rfa <- cm_df_rfa
      result_list$preds_for_buidling_confmat_rfa <- rfa_res$predictions
    }

    if (do_top_importance_selection) {
      result_list$top_importance_results <- top_importance_results
    }

    # Save the result for this fold into the master list
    final_final_res_list[[f]] <- result_list
  }

  # -----------------------------------------------------------
  # NEW: Name the list elements as fold1, fold2, etc.
  # -----------------------------------------------------------
  names(final_final_res_list) <- paste0("fold", seq_along(final_final_res_list))

  return(final_final_res_list)
}
