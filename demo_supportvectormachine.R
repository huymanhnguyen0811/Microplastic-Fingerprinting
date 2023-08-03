library(e1071)
library(caret)
library(ROCR)
library(MLmetrics)


# e1071 package: ===========================
# PCA-based feature data
# add plastic_type as factor
add_p_type <- function(data) {
  newdata <- data %>%
    rownames_to_column(., var = "File") %>%
    mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                                 ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                        ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                               ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                      ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                             ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
    mutate(plastic_type = factor(plastic_type, levels = unique(plastic_type))) %>%
    relocate(plastic_type, .before = 1) %>%
    column_to_rownames(., var = "File")
  
  return(newdata)
}

PCAtools_mergePC <- add_p_type(PCAtools_mergePC)
e1071_merge_PC <- add_p_type(e1071_merge_PC)

# PArtitioning train&test sets / training / predict on test set
e1071.SVM.result <- function(dat, split.ratio){
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  tune.out <- e1071::tune(e1071::svm, plastic_type ~ ., data = plastic_trn, kernel = "radial",
                          ranges = list(cost = seq(from = 0.1, to = 1, by = 0.01),
                                        gamma = 0.01,
                                        epsilon = seq(from = 0.0001, to = 0.001, by = 0.0001)),
                          decision.values = TRUE, 
                          probability = TRUE)
  bestmod <- tune.out$best.model
  
  # Prediction results
  pred_prob <- predict(bestmod, newdata = plastic_tst,
                       decision.values = TRUE, probability = TRUE)
  pred_res <- attr(pred_prob, "probabilities")
  return(pred_res)
}

PCAtools_mergePC.SVMresult <- e1071.SVM.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.SVMresult <- e1071.SVM.result(e1071_merge_PC, split.ratio = 0.6)
my.data.SVMresult <- e1071.SVM.result(my.data, split.ratio = 0.6)

plot.dat <- function(dat) {
  newdat <- as.data.frame(dat) %>%
    tibble::rownames_to_column(., var = "File") %>%
    tidyr::pivot_longer(., cols = 2:ncol(.), names_to = "plastic_type", values_to = "prob")
  return(newdat)
}   

PCAtoolsdat <- plot.dat(PCAtools_mergePC.SVMresult)
e1071dat <- plot.dat(e1071_merge_PC.SVMresult)
mydat <- plot.dat(my.data.SVMresult)

# Make bar graph of prediction probability
ggplot(data = mydat) + 
  geom_col(aes(x = File, y = prob, fill = plastic_type), 
           position = "dodge" # separating stacking prob cols
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(n.breaks = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))


# Caret package: ===================================
caret.SVM.result <- function(dat, split.ratio) {
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = 0.6, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]
  
  # Control params for SVM
  ctrl <- caret::trainControl(
    method = "cv", 
    number = 10, # Performing 10-fold CrossValidation
    classProbs = TRUE,                 
    summaryFunction = multiClassSummary  # also needed for AUC/ROC
  )
  
  # Tune an SVM with caret
  caret_svm_auc <- caret::train(
    plastic_type ~ ., 
    data = plastic_trn,
    method = "svmRadial",               
    preProcess = c("center", "scale"),
    metric = "ROC",  # area under ROC curve (AUC)       
    trControl = ctrl,
    tuneLength = 10
  )
  
  
  # Prediction
  pred_prob <-predict(caret_svm_auc, newdata = plastic_tst, type = "prob")
  rownames(pred_prob) <- rownames(plastic_tst)
  
  # Confusion Matrix
  # caret::confusionMatrix(caret_svm_auc)
  
  # Feature importance
  
  vip_all <- c()
  for (pt in unique(plastic_trn$plastic_type)) {
    set.seed(282)  # for reproducibility
    
    # Add prediction wrapper function to return the predicted class probabilities for the reference class of interest.
    prob <- function(object, newdata) {
      predict(object, newdata = newdata, type = "prob")[, pt] # Food_Packaging_Waste
    }
    
    # Variable importance plot
    vip_res <- vip::vip(caret_svm_auc, method = "permute", nsim = 5, train = plastic_trn, 
                        target = "plastic_type", metric = "auc", reference_class = pt,
                        pred_wrapper = prob)
    vip_all <- c(vip_all, vip_res$data[1:3,]$Variable)
  }
  
  # construct PDPs for the top four features of reference_class
  # features <- vip_all
  # 
  # pdps <- lapply(features, function(x) {
  #   partial(caret_svm_auc, pred.var = x, which.class = 2,  
  #           prob = TRUE, plot = TRUE, plot.engine = "ggplot2") +
  #     coord_flip()
  # })
  # 
  # grid.arrange(grobs = pdps,  ncol = 5)
  # 
  
  return(list(pred_prob, vip_all))
}

PCAtools_mergePC.SVMresult <- caret.SVM.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.SVMresult <- caret.SVM.result(e1071_merge_PC, split.ratio = 0.6)
my.data.SVMresult <- caret.SVM.result(my.data, split.ratio = 0.6)

View(PCAtools_mergePC.SVMresult[[1]])
View(e1071_merge_PC.SVMresult[[1]])
View(my.data.SVMresult[[1]])

PCAtoolsdat <- plot.dat(PCAtools_mergePC.SVMresult[[1]])
e1071dat <- plot.dat(e1071_merge_PC.SVMresult[[1]])
mydat <- plot.dat(my.data.SVMresult[[1]])

# Make bar graph of prediction probability
ggplot(data = mydat) + 
  geom_col(aes(x = File, y = prob, fill = plastic_type), 
           position = "dodge" # separating stacking prob cols
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(n.breaks = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))
