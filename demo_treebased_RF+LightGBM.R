# RandomForestSRC ----------------
library(randomForestSRC)
# https://www.randomforestsrc.org/articles/getstarted.html

# WARNING!!! -> Values from `Percent_Area` are not uniquely identified -> ASK ROXANA!!!!
View(shared_comp_plastic_type %>%
       dplyr::group_by(File, collapsed_compound) %>%
       dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
       dplyr::filter(n > 1L))

# PCA-based feature data ==========================
# add plastic_type as factor
ptype_n_samp <- function(data, use) {
  newdata <- data %>%
    tibble::rownames_to_column(., var = "File") %>%
    dplyr::mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                                        ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                               ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                                      ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                             ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                                    ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
    dplyr::mutate(Type = ifelse(str_detect(File, "USSB"), "Store-Bought", "Environmental")) %>%
    dplyr::mutate(plastic_type = factor(plastic_type, levels = unique(plastic_type))) %>%
    dplyr::mutate(Type = factor(Type, levels = unique(Type))) %>%
    dplyr::relocate(c(Type, plastic_type), .before = 1) %>%
    tibble::column_to_rownames(., var = "File")
  
  if (use == "min") {
    # identify plastic type with minimum number of samples
    min_samp_size <- as.numeric(newdata %>% count(plastic_type) %>% summarise(min(n)))
    
    # select random number of rows for each plastic type == minimum number of samples
    newdf <- list()
    i <- 1
    for (plt in unique(newdata$plastic_type)) {
      tmpdf <- newdata[which(newdata$plastic_type == plt),]
      # random sampling within the levels of y when y is a factor to balance the class distributions within the splits.
      idx <- caret::createDataPartition(tmpdf$Type, 
                                        p = base::round(min_samp_size/nrow(tmpdf), 1), 
                                        list = F)
      newdf[[i]] <- tmpdf[idx, ]
      i <- i + 1
    }
    
    df <- bind_rows(newdf)
  }
  else {
    df <- newdata
  }
  return(df)
}

PCAtools_alldf <- ptype_n_samp(PCAtools_mergePC, use = "all")
e1071_alldf <- ptype_n_samp(e1071_merge_PC, use = "all")
PCAtools_mindf <- ptype_n_samp(PCAtools_mergePC, use = "min")
e1071_mindf <- ptype_n_samp(e1071_merge_PC, use = "min")

rfsrc.result <- function(dat, split.ratio){
  set.seed(1234)
  # Partitioning train&test sets / training / predict on test set
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  # Train
  mergePC.rf <- rfsrc(plastic_type ~ ., 
                      ntree=20000, 
                      splitrule = "auc", 
                      nodesize = 1, #Minumum size of terminal node for classification (1)
                      # mtry = 21,
                      importance = "permute", 
                      samptype = "swr",
                      membership = TRUE,
                      perf.type="misclass",
                      block.size = 1, # cumulative error rate on every tree
                      data = plastic_trn
  )
  print(mergePC.rf)

  # Prediction results
  pred_res <- predict(mergePC.rf, newdata = plastic_tst, type = "prob")$predicted
  rownames(pred_res) <- rownames(plastic_tst)
  
  imp_var <- vimp(mergePC.rf, importance = "permute")$importance
  # selection of the best feature candidates
  # md.obj <- max.subtree(mergePC.rf)
  # best.feature <- md.obj$topvars # extracts the names of the variables in the object md.obj
  
  return(list(mergePC.rf, pred_res, imp_var))
}

PCAtools_alldf.RFresult <- rfsrc.result(PCAtools_alldf, split.ratio = 0.6)
e1071_alld.RFresult <- rfsrc.result(PCAtools_alldf, split.ratio = 0.6)
PCAtools_mindf.RFresult <- rfsrc.result(PCAtools_alldf, split.ratio = 0.6)
e1071_mindf.RFresult <- rfsrc.result(e1071_mindf, split.ratio = 0.6)

# my.data.RFresult <- rfsrc.result(my.data, split.ratio = 0.6)

rfsrc.plots <- function(preddat, rf.mod) {
  newdat <- as.data.frame(preddat) %>%
    tibble::rownames_to_column(., var = "File") %>%
    tidyr::pivot_longer(., cols = 2:ncol(.), names_to = "plastic_type", values_to = "prob")
  
  plots <- list()
  # Make bar graph of prediction probability -----
  plots[[1]] <- ggplot(data = newdat) + 
    geom_col(aes(x = File, y = prob, fill = plastic_type), 
             position = "dodge" # separating stacking prob cols
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(n.breaks = 10) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))
  
  # plot OOB error rate against the number of trees -------
  plots[[2]] <-plot(ggRandomForests::gg_error(rf.mod))
  
  # Estimate the variables importance --------
  my.rf.vimp <- ggRandomForests::gg_vimp(rf.mod, nvar = 100) # provides the predictor's importance of top 100 predictors
  plots[[3]] <- plot(my.rf.vimp, theme(axis.text.x = element_text(angle = 90))) # visualises the predictor’s importance
  
  grid.arrange(arrangeGrob(
    plots[[3]], 
    arrangeGrob(plots[[1]], plots[[2]]),
    ncol=2 ,widths=c(2,1)))
}

rfsrc.plots(preddat = PCAtools_alldf.RFresult[[2]],
            rf.mod = PCAtools_alldf.RFresult[[1]])

# Partial dependence plots - aka. response ~ each predictor --------------
my.rf.part.plot <- plot.variable(mergePC.rf, partial=TRUE, sorted=FALSE,
                                 show.plots=FALSE, 
                                 nvar = 10 # look at top 10 predictor's importance
                                 )
gg.part <- ggRandomForests::gg_partial(my.rf.part.plot)
plot(gg.part, xvar=names(plastic_trn[,-1]), panel=TRUE, se=TRUE)


# LightGBM for Multiclass Classification ==================
library(lightgbm)

lightgbm.result <- function(dat, split.ratio) {
  # Must convert plastic_type from factors to numeric
  test <- copy(dat)
  test$plastic_type <- as.numeric(as.factor(test$plastic_type)) - 1L
  # Split train and test sets
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(test$plastic_type, p = 0.6, list = F)
  plastic_trn <- as.matrix(test[plastic_idx, ])
  plastic_tst <- as.matrix(test[-plastic_idx, ])
  
  dtrain <- lgb.Dataset(data = plastic_trn[, 3:ncol(plastic_trn)], 
                        label = plastic_trn[, 2])
  dtest <- lgb.Dataset.create.valid(dtrain, 
                                    data = plastic_tst[, 3:ncol(plastic_tst)], 
                                    label = plastic_tst[, 2])
  
  # Validation set to be used during training, not testing!!!
  valids <-  list(test = dtest)
  
  # Setup parameters
  params <- list(
    min_data = 1L
    , learning_rate = 1
    , objective = "multiclass"
    , metric = "multi_error"
    , num_class = 7L
  ) 
  
  lgb.model <- lgb.train(params,
                         dtrain,
                         nrounds = 1000,
                         valids, 
                         early_stopping_rounds = 25L
  )
  
  # prediction
  my_preds <- round(predict(object = lgb.model, 
                            data = plastic_tst[, 3:ncol(plastic_tst)],
                            reshape = TRUE), 3)

  colnames(my_preds) <- levels(dat$plastic_type)
  rownames(my_preds) <- rownames(plastic_tst)
  
  plots <- list()
  # Important features
  tree_imp <- lgb.importance(lgb.model, percentage = T)
  varimpmat <- lgb.plot.importance(tree_imp, measure = "Gain")
  plots[[1]] <- ggplot(data = varimpmat, aes(x = Gain, y = Feature)) +
    geom_col() +
    title("Feature importance")

  newdat <- as.data.frame(my_preds) %>%
    tibble::rownames_to_column(., var = "File") %>%
    tidyr::pivot_longer(., cols = 2:ncol(.), names_to = "plastic_type", values_to = "prob")
  
  # Make bar graph of prediction probability
  plots[[2]] <- ggplot(data = newdat) + 
    geom_col(aes(x = File, y = prob, fill = plastic_type), 
             position = "dodge" # separating stacking prob cols
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(n.breaks = 10) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))
  
  grid.arrange(grobs = plots,
               ncol = 2)
  
  return(my_preds)
}

PCAtools_mergePC.LGBMresult <- lightgbm.result(PCAtools_alldf, split.ratio = 0.6)
e1071_merge_PC.LGBMresult <- lightgbm.result(e1071_alldf, split.ratio = 0.6)
