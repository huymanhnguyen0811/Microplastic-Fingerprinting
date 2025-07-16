# Recursive Feature Addition -------------
recursive_feature_addition <- function(X_train, y_train,
                                       X_test,  y_test,
                                       ntree_candidates,
                                       metric) {
  # 1) Fit full model to get initial feature importance
  full_res <- tune_rf_subset(
    X_train, y_train,
    X_test,  y_test, ntree_candidates, metric
  )
  # Extract importance from the randomForest final model
  imp_mat <- importance(full_res$model$finalModel)
  # Choose first top feature by MeanDecreaseAccuracy (or use MeanDecreaseGini) to use as base feature to initiate RFA
  init_feat <- rownames(imp_mat)[which.max(imp_mat[, "MeanDecreaseAccuracy"])]
  
  # Initialize feature sets
  best_feats <- init_feat
  remaining  <- setdiff(colnames(X_train), best_feats)
  
  # 2) Baseline performance with the single top feature
  base_res      <- tune_rf_subset(
    X_train[, best_feats, drop = FALSE], y_train,
    X_test [, best_feats, drop = FALSE], y_test,
    ntree_candidates, metric
  )
  best_mcc <- base_res$MCC_Multiclass
  
  # 3) Iterative addition of remaining features
  while(length(remaining) > 0) {
    improved <- FALSE
    for(feat in remaining) {
      trial_feats <- c(best_feats, feat)
      res <- tune_rf_subset(
        X_train[, trial_feats, drop = FALSE], y_train,
        X_test [, trial_feats, drop = FALSE], y_test, ntree_candidates, metric
      )
      # Add feature only if all three metrics improve
      if(res$MCC_Multiclass > best_mcc) {
        best_mcc <- res$MCC_Multiclass
        best_feats    <- trial_feats
        improved      <- TRUE
        break
      }
    }
    # Drop the last-tested feature from remaining
    remaining <- setdiff(remaining, feat)
    if(!improved) break
  }
  
  # 4) Final model on the chosen feature set
  final_res <- tune_rf_subset(
    X_train[, best_feats, drop = FALSE], y_train,
    X_test [, best_feats, drop = FALSE], y_test,
    ntree_candidates, metric
  )
  
  # 5) Return performance and selected features
  list(
    best_mcc    = final_res$MCC_Multiclass,
    best_acc = final_res$Accuracy,
    best_features_rf = best_feats,
    final_model      = final_res$model,
    prob_matrix = final_res$p_m,
    predictions = final_res$preds,
    eval_metrics = final_res$eval_metrics
  )
}

#─────────────────────────────────────────────────────────────────────────────
# Helper: retune mtry (based on subset size) + ntree for any train/test split
tune_rf_subset <- function(X_train, y_train, X_test, y_test,
                           ntree_candidates, metric) {
  # 1) Check class counts
  k_inner <- min(table(y_train))
  
  # 2) Build the mtry grid
  p        <- ncol(X_train)
  mtry_vals <- unique(pmax(1, c(floor(sqrt(p)), floor(log2(p)))))
  
  # 3) Branch: enough data for CV?
  if (k_inner > 1) {
    # —> use caret CV
    ctrl <- trainControl(
      method         = "cv",
      number         = k_inner,
      classProbs     = TRUE,
      summaryFunction= defaultSummary
    )
    rf_grid <- expand.grid(mtry = mtry_vals)
    
    best_mod   <- NULL
    best_val   <- -Inf
    best_trees <- NA
    
    for (nt in ntree_candidates) {
      set.seed(123)
      tmp <- caret::train(
        x         = X_train,
        y         = as.factor(y_train),
        method    = "rf",
        trControl = ctrl,
        tuneGrid  = rf_grid,
        metric    = metric,
        ntree     = nt,
        importance= TRUE
      )
      this_val <- max(tmp$results[[metric]])
      if (this_val > best_val) {
        best_val   <- this_val
        best_mod   <- tmp
        best_trees <- nt
      }
    }
    
    caret_model <- best_mod
    final_rf    <- caret_model$finalModel
    
  } else {
    # —> too few per class: do OOB-based tuning
    best_oob   <- Inf
    best_mtry  <- NA
    best_trees <- NA
    best_rf    <- NULL
    
    for (nt in ntree_candidates) {
      for (m in mtry_vals) {
        set.seed(123)
        rf_mod <- randomForest(
          x          = X_train,
          y          = as.factor(y_train),
          mtry       = m,
          ntree      = nt,
          importance = TRUE
        )
        # final OOB error rate:
        oob_err <- tail(rf_mod$err.rate[ , "OOB"], 1)
        if (oob_err < best_oob) {
          best_oob   <- oob_err
          best_mtry  <- m
          best_trees <- nt
          best_rf    <- rf_mod
        }
      }
    }
    # wrap into a dummy caret-like object so that downstream calls to predict(<train>,…) and predict(<train>, type="prob") continue to work.
    caret_model <- list(
      finalModel = best_rf,
      # so that later code that uses predict(<train>,…) still works:
      methods = list(predict = "randomForest")
    )
    final_rf <- best_rf
  }
  
  # 4) Predict & score
  preds       <- predict(final_rf, newdata = X_test, type = "response")
  prob_matrix <- as.matrix(predict(final_rf, newdata = X_test, type = "prob"))
  ev          <- eval_metrics(true_labels = y_test, 
                              prob_matrix = prob_matrix, 
                              pred_labels = preds)
  
  return(list( model        = caret_model,
               ntree        = best_trees,
               preds        = preds,
               p_m          = prob_matrix,
               eval_metrics = ev,
               Accuracy = ev$Accuracy,
               Kappa = ev$Kappa,
               AUC = ev$AUC,
               F1_Macro     = ev$F1_Macro,
               F1_Weighted = ev$F1_Weighted,
               MCC          = ev$MCC_Multiclass
  ))
}

#─────────────────────────────────────────────────────────────────────────────
library(caret)
library(pROC)

eval_metrics <- function(true_labels, prob_matrix, pred_labels){
  # 0) Ensure factor levels align
  classes <- union(levels(true_labels), levels(pred_labels))
  y_true <- factor(true_labels, levels = classes)
  y_pred <- factor(pred_labels,  levels = classes)
  
  # 1) Build the confusion‐matrix table
  cm <- table(Actual = y_true, Predicted = y_pred)
  N  <- sum(cm)                       # total samples
  diag_cm <- diag(cm)
  
  # 2) Accuracy
  accuracy <- sum(diag_cm) / N
  
  # 3) Cohen’s Kappa (via caret)
  kappa <- tryCatch({
    cm_obj <- confusionMatrix(y_pred, y_true)
    as.numeric(cm_obj$overall["Kappa"])
  }, error = function(e) NA)
  
  # 4) One‐vs‐all multiclass AUC
  auc <- tryCatch({
    roc_obj <- multiclass.roc(response  = true_labels,
                              predictor = prob_matrix)
    as.numeric(roc_obj$auc)
  }, error = function(e) NA)
  
  # 5) Per‐class precision & recall
  col_sums <- colSums(cm)   # predicted totals for each class
  row_sums <- rowSums(cm)   # actual   totals for each class
  
  precision_per_class <- ifelse(col_sums > 0,
                                diag_cm / col_sums,
                                0)        # no preds → p=0
  recall_per_class    <- ifelse(row_sums > 0,
                                diag_cm / row_sums,
                                0)        # no actual → r=0
  
  # 6) Per‐class F1 & weighted F1
  f1_per_class <- ifelse(
    (precision_per_class + recall_per_class) > 0,
    2 * precision_per_class * recall_per_class /
      (precision_per_class + recall_per_class),
    0
  )
  support <- row_sums   # how many true examples of each class
  f1_weighted <- sum(f1_per_class * support) / sum(support)
  
  # 7) Multiclass MCC (Gorodkin’s formula)
  t_i <- row_sums
  p_j <- col_sums
  c_val <- sum(diag_cm)         # sum of true positives
  num   <- N * c_val - sum(t_i * p_j)
  den   <- sqrt((N^2 - sum(t_i^2)) * (N^2 - sum(p_j^2)))
  mcc   <- ifelse(den > 0, num / den, 0)
  
  return(list(
    Accuracy         = accuracy,
    Kappa            = kappa,
    AUC              = auc,
    F1_Weighted      = f1_weighted,
    MCC_Multiclass   = mcc
  ))
}


####################################################
##########        CONFUSION MATRIX       ###########
####################################################
conf_mat_plot <- function(y_test, preds, conf_mat_title, accuracy) {
  # 1) Compute raw counts table: rows = treated labels, cols = predicted pristine labels
  cm_tab <- table(
    Actual    = as.character(y_test),
    Predicted = as.character(preds)
  )
  
  # 2) Turn it into a data.frame for ggplot
  cm_df <- as.data.frame(cm_tab, stringsAsFactors = FALSE) %>%
    dplyr::rename(Freq = Freq)
  
  # 3) Add total-per-row and percent-per-cell
  cm_df$Total <- NA
  cm_df$Percent <- NA
  for (label in unique(cm_df$Actual)){
    cm_df[which(cm_df$Actual == label),]$Total <- sum(cm_df[which(cm_df$Actual == label),]$Freq)
  }
  
  cm_df <- cm_df %>%
    mutate(
      Percent = round(Freq / Total * 100, 2)
    )
  
  
  cm_df <- cm_df %>%
    mutate(
      Label = ifelse(Predicted == Actual,
                     paste0(round(Percent,1),"%"),
                     ""
      )
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
      subtitle = paste0(
      sprintf("Overall linking accuracy: %.2f%%\n", accuracy)
      # # " Recall per treated class (N=", cm_df$Total[1], 
      # " replicates each)"
      ),
      x = "Predicted Class (Store-Bought)",
      y = "Actual Class (Environmental)",
      fill = "Percentage of\nClassification (%)" # "Recall %"
    ) +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title   = element_text(face = "bold", hjust = 0.5),
      plot.subtitle= element_text(hjust = 0.5)
    )
  return(list(conf_mat = conf_mat, cm_df = cm_df))
}

#─────────────────────────────────────────────────────────────────────────────

run_rf_analysis_manuscript1 <- function(data,
                                        type_col,
                                        remove_cols,
                                        train_proportion        = 0.5,
                                        use_source_split        = FALSE,               # Option 1: environmental 50:50 split and combine with store-bought
                                        use_store_vs_environmental_split = FALSE,       # Option 2: use store-bought for training and environmental for testing
                                        use_polymer             = FALSE,               # Whether to include "Polymer" as a predictor
                                        do_pairwise_test        = FALSE,               # Run pairwise/ANOVA tests on imputed/normalized data
                                        do_top_importance_selection = FALSE,           # Run additional classification using top-N features by importance
                                        top_importance_counts   = c(100, 50, 25, 10),   # Which top features to try
                                        data_name,              # Name of the input data combination for figure titles
                                        seed                    = 123,
                                        ntree_candidates        = c(100, 500, 1000, 2500),
                                        metric                  = "Accuracy",
                                        do_rfe                  = FALSE,
                                        do_rfa                  = FALSE,
                                        min_sample_number,
                                        do_impute_norm_screen   = TRUE,
                                        group_for_significance) {
  
  
  # 1) Prepare Target Variable and Remove Classes with a Single Sample
  data[[type_col]] <- as.factor(data[[type_col]])
  levels(data[[type_col]]) <- make.names(levels(data[[type_col]]), unique = TRUE)
  
  # Set minimum number of observation for plastic category
  class_counts  <- table(data[[type_col]])
  keep_classes  <- names(class_counts[class_counts >= min_sample_number])
  data          <- data[data[[type_col]] %in% keep_classes, ]
  data[[type_col]] <- droplevels(data[[type_col]])
  rownames(data) <- NULL
  data <- data %>% 
    mutate(sample_id = paste0(Plastic_type, "_", technique)) %>% 
    mutate(sample_id = ave(sample_id, sample_id, FUN = function(x) {
      if (length(x) > 1) paste0(x, "_rep", seq_along(x)) else x
    })) %>% 
    tibble::column_to_rownames(., var = "sample_id")
  
  # 3) Imputation and Normalization on Clean Data
  best_imp  <- "No_imputation"
  best_norm <- "No_normalization"
  if(do_impute_norm_screen){
    start_time_imp_norm <- Sys.time()
    cat("\n### Screening Best Imputation + Normalization Combo ###\n")
    
    # Identify numeric columns (excluding remove_cols and target)
    numeric_cols <- names(which(sapply(data, is.numeric)))
    numeric_cols <- setdiff(numeric_cols, c(remove_cols, type_col))
    
    best_res_imp_norm <- find_best_impute_normalize(
      df = data,
      type_col = type_col,
      group_for_significance = group_for_significance,
      remove_cols = remove_cols,
      use_store_vs_environmental_split = use_store_vs_environmental_split,
      use_source_split = use_source_split
    )
    cat("\n*** Summary of Imputation/Normalization Combos ***\n")
    print(best_res_imp_norm$results_table)
    cat("\n*** Best Combo ***\n")
    print(best_res_imp_norm$best_combo)
    
    best_imp  <- best_res_imp_norm$best_combo$Imputation
    best_norm <- best_res_imp_norm$best_combo$Normalization
    
    X_original <- data[, numeric_cols, drop = FALSE]
    X_imputed <- imputation_methods[[best_imp]](as.data.frame(X_original))
    X_final <- normalization_methods[[best_norm]](X_imputed)
    data[, numeric_cols] <- X_final
    
    end_time_imp_norm   <- Sys.time()
    time_imp_norm      <- as.numeric(difftime(end_time_imp_norm,
                                              start_time_imp_norm, 
                                              units="secs"))
  } else {
    cat("\n### Skipping Imputation + Normalization; Using Raw Data ###\n")
  }
  
  cat("NAs in original data:", sum(is.na(data)), "\n")
  cat("NAs after imputation+normalization:", sum(is.na(data)), "\n")
  
  # 4) Run Pairwise/ANOVA Tests on Cleaned Data (if selected)
  if(do_pairwise_test) {
    cat("\n### Running Pairwise Significance Tests on Cleaned Data ###\n")
    pairwise_res <- pairwise_significance_tests(input_df = data, 
                                                group_col = group_for_significance, 
                                                start_col_index = 6)
    pairwise_sig_features <- pairwise_res$corrected
    if(length(pairwise_sig_features) == 0) {
      cat("No significant features found from pairwise tests. Skipping pairwise feature selection.\n")
      do_pairwise_test <- FALSE
    }
  }
  
  # 5) Train/Test Split
  if(use_store_vs_environmental_split){
    cat("\n### Using Store-vs-Environmental Splitting ###\n")
    
    store_bought_set <- data[data$Source == "Store-Bought", ]
    environmental_set <- data[data$Source == "Environmental", ]
    
    # Setup k-fold
    set.seed(123)
    min_class_count_ENV <- min(table(as.character(environmental_set[[type_col]])))
    min_class_count_SB <- min(table(as.character(environmental_set[[type_col]])))
    K <- ifelse(min(min_class_count_ENV, min_class_count_SB) > 2, min_class_count_ENV, 2)
    # stratified folds on the full target vector
    K_folds_SB <- createFolds(store_bought_set[[type_col]], k = K, returnTrain = TRUE)
    K_folds_ENV <- createFolds(environmental_set[[type_col]], k = K, returnTrain = TRUE)
    
  } else if(use_source_split){
    cat("\n### Using Source-based Splitting with Environmental 50:50 Split ###\n")
    
    store_bought_set <- data[data$Source == "Store-Bought", ]
    environmental_set <- data[data$Source == "Environmental", ]
    
    # Setup k-fold
    set.seed(123)
    min_class_count_ENV <- min(table(as.character(environmental_set[[type_col]])))
    min_class_count_SB <- min(table(as.character(environmental_set[[type_col]])))
    K <- ifelse(min(min_class_count_ENV, min_class_count_SB) > 2, min_class_count_ENV, 2)
    # stratified folds on the full target vector
    K_folds_SB <- createFolds(store_bought_set[[type_col]], k = K, returnTrain = TRUE)
    K_folds_ENV <- createFolds(environmental_set[[type_col]], k = K, returnTrain = TRUE)
  }
  
  # Setup K-folds CV 
  final_final_res_list <- vector("list", length = K)
  i <- 1
  
  for(f in 1:K) {
    if(use_store_vs_environmental_split){
      train_idx <- K_folds_SB[[f]]
      test_idx  <- K_folds_ENV[[f]]
      
      train_data <- store_bought_set[train_idx, ]
      test_data  <- environmental_set[test_idx, ]
    } else if(use_source_split){
      train_idx_SB <- K_folds_SB[[f]]
      train_idx_ENV  <- K_folds_ENV[[f]]
      test_idx_ENV  <- setdiff(seq_len(nrow(environmental_set)), K_folds_ENV[[f]])
      
      train_data <- rbind(store_bought_set[train_idx_SB, ], environmental_set[train_idx_ENV, ])
      test_data  <- environmental_set[test_idx_ENV, ]
    }
    # Drop unused factor levels after splitting
    train_data[[type_col]] <- droplevels(train_data[[type_col]])
    test_data[[type_col]]  <- droplevels(test_data[[type_col]])
    
    # 6) Prepare Predictors and Response for Full-Feature Classification
    X_train <- train_data %>% dplyr::select(-c(all_of(remove_cols), all_of(type_col)))
    y_train <- train_data[[type_col]]
    X_test  <- test_data %>% dplyr::select(-c(all_of(remove_cols), all_of(type_col)))
    y_test  <- test_data[[type_col]]
    
    # 8) Full-Feature Classification: Grid Search over (mtry, ntree)
    cat("\n### Train/test with Full original features ###\n")
    start_time_full_features <- Sys.time()
    
    all_feats <- tune_rf_subset(
      X_train, y_train, X_test, y_test,
      ntree_candidates, metric
    )
    
    best_model_all_feats      <- all_feats$model
    best_trees_all_feats      <- all_feats$ntree
    acc_full        <- all_feats$Accuracy
    
    # --- Generate Heatmap (Full Features) ---
    {prob_matrix_all_feats <- all_feats$p_m
      conf_mat <- conf_mat_plot(y_test, all_feats$preds, 
                                          conf_mat_title = paste0("All Features - ", data_name), 
                                          accuracy = acc_full)
      conf_mat_all_feats <- conf_mat$conf_mat
      cm_df_all_feats <- conf_mat$cm_df
    }
    
    eval_metrics_all_feats <- all_feats$eval_metrics
    
    end_time_full_features   <- Sys.time()
    time_all_feats       <- as.numeric(difftime(end_time_full_features,
                                                start_time_full_features, 
                                                units="secs"))
    
    # 9) Classification with Pairwise/ANOVA Selected Features (if requested)
    cat("\n### Train/test with Pair-wise features ###\n")
    
    if(do_pairwise_test){
      start_time_sig <- Sys.time()
      pairwise_res <- pairwise_significance_tests(input_df = data, 
                                                  group_col = group_for_significance, 
                                                  start_col_index = 6)
      pairwise_sig_features <- pairwise_res$corrected
      sig_feats <- intersect(pairwise_sig_features, colnames(X_train))
      
      if(length(sig_feats)>0) {
        
        tmp_sig <- tune_rf_subset(
          X_train[, sig_feats, drop=FALSE], y_train,
          X_test [, sig_feats, drop=FALSE], y_test,
          ntree_candidates, metric
        )
        acc_sig         <- tmp_sig$Accuracy
        final_sig_model <- tmp_sig$model$finalModel
      }
      
      {prob_matrix_sig <- tmp_sig$p_m
        conf_mat <- conf_mat_plot(y_test, tmp_sig$preds, 
                                      conf_mat_title = paste0("Pair-wise Significance-Based Features - ", data_name), 
                                      accuracy = acc_sig)
        conf_mat_sig <- conf_mat$conf_mat
        cm_df_sig <- conf_mat$cm_df 
      }
      
      eval_metrics_sig <- tmp_sig$eval_metrics
      
      end_time_sig   <- Sys.time()
      time_sig       <- as.numeric(difftime(end_time_sig, start_time_sig, units="secs"))
    }
    
    # Recursive Feature Addition
    if(do_rfa){
      cat("\n### Train/test with RFA features ###\n")
      start_time_rfa <- Sys.time()
      
      rfa_res   <- recursive_feature_addition(
        X_train, y_train,
        X_test,  y_test,
        ntree_candidates = ntree_candidates,
        metric           = metric
      )
      
      rfa_selected_feats <- rfa_res$best_features_rf
      acc_rfa                 <- rfa_res$best_accuracy
      final_rf_rfa            <- rfa_res$final_model$finalModel
      
      prob_matrix_rfa <- rfa_res$prob_matrix
      
      conf_mat <- conf_mat_plot(y_test, rfa_res$predictions, 
                                    conf_mat_title = paste0("RFA Features - ", type_col), 
                                    accuracy = acc_rfa)
      conf_mat_rfa <- conf_mat$conf_mat
      cm_df_rfa <- conf_mat$cm_df
      
      eval_metrics_rfa <- rfa_res$eval_metrics
      
      end_time_rfa   <- Sys.time()
      time_rfa       <- as.numeric(difftime(end_time_rfa, start_time_rfa, units="secs"))
    }
    
    # 10) Classification with Top-N Features by Importance (if requested)
    if(do_top_importance_selection){
      
      cat("\n### Running Classification on Top-N Features by Importance ###\n")
      
      rf_importances <- as.data.frame(best_model_all_feats$finalModel$importance)
      rf_importances$Feature <- rownames(rf_importances)
      rf_importances <- rf_importances[order(rf_importances$MeanDecreaseGini, decreasing = TRUE), ]
      
      top_importance_results <- list()
      total_features <- nrow(rf_importances)
      
      for(top_n in top_importance_counts){
        start_time_top_100_50_25_10 <- base::Sys.time()
        
        if(total_features < top_n){
          selected_features <- rf_importances$Feature
          cat(sprintf("Only %d features available; using all.\n", total_features))
        } else {
          selected_features <- head(rf_importances$Feature, top_n)
        }
        
        tmp_top <- tune_rf_subset(
          X_train[, selected_features, drop = FALSE], y_train,
          X_test[, selected_features, drop = FALSE], y_test,
          ntree_candidates, metric
        )
        
        acc_top         <- tmp_top$Accuracy
        final_rf_top <- tmp_top$model$finalModel
        
        prob_matrix_top <- tmp_top$p_m
        conf_mat <- conf_mat_plot(y_test, tmp_top$preds, 
                                      conf_mat_title = paste0("Top ", length(selected_features), " Features - ", data_name),
                                      accuracy = acc_top)
        conf_mat_top <- conf_mat$conf_mat
        cm_df_top <- conf_mat$cm_df
        
        eval_metrics_top <- tmp_top$eval_metrics
        
        end_time_top_100_50_25_10   <- Sys.time()
        time_top_100_50_25_10   <- as.numeric(difftime(end_time_top_100_50_25_10,
                                                       start_time_top_100_50_25_10, 
                                                       units="secs"))
        
        top_importance_results[[paste0("top_", top_n)]] <- list(
          best_model = final_rf_top,
          selected_features = selected_features,
          prob_matrix_top = prob_matrix_top,
          conf_mat_top=conf_mat_top,
          eval_metrics_top = eval_metrics_top,
          cm_df_top = cm_df_top,
          time_top_100_50_25_10 = time_top_100_50_25_10
        )
      }
      
    }
    
    # 11) (Optional) RFE with Fixed Hyperparameters --------------
    if(do_rfe){
      cat("\n### Train/test with RFE features ###\n")
      start_time_rfe <- Sys.time()
      
      myFuncs <- rfFuncs
      myFuncs$fit <- function(x, y, first, last, ...) {
        res <- tune_rf_subset(
          X_train = x, y_train = y,
          X_test  = x, y_test  = y,
          ntree_candidates, metric
        )
        res$model$finalModel
      }
      
      # compute the smallest class count
      rfe_k <- min(table(train_data[[type_col]]))
      rfe_k <- if (rfe_k > 2) rfe_k else 2
      
      rfe_ctl <- rfeControl(
        functions   = myFuncs,
        method      = "cv",
        number      = rfe_k,
        saveDetails = TRUE,
        returnResamp= "final"
      )
      rfe_out <- rfe(
        x          = X_train,
        y          = y_train,
        sizes      = seq_len(ncol(X_train)),
        rfeControl = rfe_ctl
      )
      rfe_feats <- rfe_out$optVariables
      
      # final retune+test on RFE set
      tmp_rfe <- tune_rf_subset(
        X_train[, rfe_feats, drop=FALSE], y_train,
        X_test[, rfe_feats, drop=FALSE], y_test,
        ntree_candidates, metric
      )
      acc_rfe          <- tmp_rfe$Accuracy
      final_rfe_model  <- tmp_rfe$model$finalModel
      
      # --- Generate Heatmap (RFE) ---
      {prob_matrix_rfe <- tmp_rfe$p_m
        conf_mat <-  conf_mat_plot(y_test, tmp_rfe$preds, 
                                       conf_mat_title = paste0("RFE Features - ", data_name), 
                                       accuracy = acc_rfe)
        conf_mat_rfe <- conf_mat$conf_mat
        cm_df_rfe <- conf_mat$cm_df
      }
      
      eval_metrics_rfe <- tmp_rfe$eval_metrics
      
      end_time_rfe   <- Sys.time()
      time_rfe       <- as.numeric(difftime(end_time_rfe, start_time_rfe, units="secs"))
    }
    
    # Return a list of final results
    result_list <- list(
      final_imp_norm_dat    = data,
      best_model_all_feats      = best_model_all_feats,
      final_rf_all_feats    = best_model_all_feats$finalModel,
      prob_matrix = prob_matrix_all_feats,
      all_features_acc     = acc_full,
      time_all_features = time_all_feats,
      eval_metrics_all_feats = eval_metrics_all_feats,
      conf_mat_all_feats = conf_mat_all_feats, 
      cm_df_all_feats = cm_df_all_feats
    )
    
    if(do_impute_norm_screen) {
      result_list$imp_norm_res_table <- best_res_imp_norm$results_table
      result_list$time_imp_norm <- time_imp_norm
      result_list$best_imputatation     <- best_imp
      result_list$best_normalization    <- best_norm
    }
    
    if(do_pairwise_test){
      result_list$sig_model    <- final_sig_model
      result_list$sig_selected_feats <- sig_feats
      result_list$prob_matrix_sig <- prob_matrix_sig
      result_list$acc_sig <- acc_sig
      result_list$time_sig          <- time_sig
      result_list$eval_metrics_sig <- eval_metrics_sig
      result_list$conf_mat_sig <- conf_mat_sig
      result_list$cm_df_sig <- cm_df_sig
    }
    
    if(do_rfa) {
      result_list$rfa_model <- final_rf_rfa
      result_list$rfa_selected_feats <- rfa_selected_feats
      result_list$prob_matrix_rfa    <- prob_matrix_rfa
      result_list$acc_rfa  <- acc_rfa
      result_list$time_rfa                <- time_rfa
      result_list$eval_metrics_rfa <- eval_metrics_rfa
      result_list$conf_mat_rfa <- conf_mat_rfa
      result_list$cm_df_rfa <- cm_df_rfa
    }
    
    if(do_rfe) {
      result_list$rfe_model             <- rfe_out
      result_list$final_rf_rfe          <- final_rfe_model
      result_list$rfe_selected_features <- rfe_feats
      result_list$prob_matrix_rfe    <- prob_matrix_rfe
      result_list$acc_rfe  <- acc_rfe
      result_list$time_rfe              <- time_rfe
      result_list$eval_metrics_rfe <- eval_metrics_rfe
      result_list$conf_mat_rfe <- conf_mat_rfe
      result_list$cm_df_rfe <- cm_df_rfe
    }
    
    if(do_top_importance_selection) {
      result_list$top_importance_results <- top_importance_results
    }
    
    final_final_res_list[[i]] <- result_list
    i <- i + 1
  } 
  return(final_final_res_list)
}


# Option 1: Use SB to predict ENV ##########################################################
data_combinations <- list(
  gc                 = gc_combined_SB_ENV_shared_cols_ntr_clean,
  hplc               = hplc_combined_SB_ENV_shared_cols_ntr_clean,
  icp                = icp_combined_SB_ENV_shared_cols_ntr_clean ,
  gc_hplc            = gc_hplc_combined_SB_ENV_shared_cols_ntr_clean,
  gc_icp             = gc_icp_combined_SB_ENV_shared_cols_ntr_clean,
  hplc_icp           = hplc_icp_combined_SB_ENV_shared_cols_ntr_clean,
  gc_hplc_icp        = gc_hplc_icp_combined_SB_ENV_shared_cols_ntr_clean
)

# Predefine the result list with names from data_combinations
run_rf_analysis_list_SB_train <- vector("list", length(data_combinations))
names(run_rf_analysis_list_SB_train) <- c("gc_USSB_train",
                                          "hplc_USSB_train",
                                          "icp_USSB_train",
                                          "gc_hplc_USSB_train",
                                          "gc_icp_USSB_train",
                                          "hplc_icp_USSB_train",
                                          "gc_hplc_icp_USSB_train")

# Iterate through each data combination by index
for(i in seq_along(data_combinations)) {
  cat("Processing dataset:", names(data_combinations)[i], "\n")
  
  # Run the analysis using the corresponding data and pass the custom name for figure titles
  run_rf_analysis_list_SB_train[[names(run_rf_analysis_list_SB_train)[i]]] <- run_rf_analysis_manuscript1(
    data = data_combinations[[i]],
    data_name = names(data_combinations)[i], 
    use_store_vs_environmental_split = TRUE,
    type_col                = "Plastic_type", 
    remove_cols             = c("Source", "Polymer", "technique", "Subcategory"),
    group_for_significance = "Plastic_type", 
    do_pairwise_test = TRUE,             # Adjust parameters as needed
    do_top_importance_selection = TRUE,
    do_rfe = TRUE,
    do_rfa = TRUE,
    train_proportion = 0.8,
    do_impute_norm_screen = TRUE,
    seed = 123,
    min_sample_number = 3)
}

## Option 2: Use Both SB and ENV to predict ENV ##########################################################
data_combinations <- list(
  gc                 = gc_combined_SB_ENV_shared_cols_ntr_clean,
  hplc               = hplc_combined_SB_ENV_shared_cols_ntr_clean,
  icp                = icp_combined_SB_ENV_shared_cols_ntr_clean,
  gc_hplc            = gc_hplc_combined_SB_ENV_shared_cols_ntr_clean,
  gc_icp             = gc_icp_combined_SB_ENV_shared_cols_ntr_clean,
  hplc_icp           = hplc_icp_combined_SB_ENV_shared_cols_ntr_clean,
  gc_hplc_icp        = gc_hplc_icp_combined_SB_ENV_shared_cols_ntr_clean
)

# Predefine the result list with names from data_combinations
run_rf_analysis_list_SB_and_ENV_train <- vector("list", length(data_combinations))
names(run_rf_analysis_list_SB_and_ENV_train) <- c("gc_SB_and_ENV_train",
                                                  "hplc_SB_and_ENV_train",
                                                  "icp_SB_and_ENV_train",
                                                  "gc_hplc_SB_and_ENV_train",
                                                  "gc_icp_SB_and_ENV_train",
                                                  "hplc_icp_SB_and_ENV_train",
                                                  "gc_hplc_icp_SB_and_ENV_train")

# Iterate through each data combination by index
for(i in seq_along(data_combinations)) {
  cat("Processing dataset:", names(data_combinations)[i], "\n")
  
  # Run the analysis using the corresponding data and pass the custom name for figure titles
  run_rf_analysis_list_SB_and_ENV_train[[names(run_rf_analysis_list_SB_and_ENV_train)[i]]] <- run_rf_analysis_manuscript1(
    data = data_combinations[[i]],
    data_name = names(data_combinations)[i], 
    use_source_split = TRUE,
    type_col                = "Plastic_type", 
    remove_cols             = c("Source", "Polymer", "technique", "Subcategory"),
    group_for_significance = "Plastic_type",
    do_pairwise_test = TRUE,             # Adjust parameters as needed
    do_top_importance_selection = TRUE,
    do_rfe = TRUE,
    do_rfa = TRUE,
    train_proportion = 0.8,
    do_impute_norm_screen = TRUE,
    seed = 123,
    min_sample_number = 3)
}

