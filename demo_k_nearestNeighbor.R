# With caret package ---------------------------
library(caret)

#k-NN cannot handle missing values
for (c in 2:ncol(my.data)) { 
  my.data[which(base::is.na(my.data[,c])), c] <- runif(length(which(base::is.na(my.data[,c]))),
                                                       min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                       max = sort(shared_comp_plastic_type$Percent_Area)[2])
}

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
caret.kNN.result <- function(dat, split.ratio){
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  default_knn_mod <- caret::train(
    plastic_type ~ .,
    data = plastic_trn,
    method = "knn",
    trControl = trainControl(method = "cv", 
                             number = 5, # using 5-fold cross-validation
                             classProbs = TRUE), 
    # preProcess = c("center", "scale"),
    tuneGrid = expand.grid(k = seq(1, 9, by = 1))
  )
  
  
  # Prediction results
  knn_pred_res <- predict(default_knn_mod, newdata = plastic_tst, type = "prob")
  rownames(knn_pred_res) <- rownames(plastic_tst)
  
  pROC.res <- pROC::multiclass.roc(response = plastic_tst$plastic_type, predictor = knn_pred_res)
  
  return(list(knn_pred_res, pROC.res))
}

PCAtools_mergePC.kNNresult <- caret.kNN.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.kNNresult <- caret.kNN.result(e1071_merge_PC, split.ratio = 0.6)

View(PCAtools_mergePC.kNNresult[[2]])
View(e1071_merge_PC.kNNresult[[2]])


# With mlr3verse ---------------------
library(mlr3verse)
library(stats)
library(nnet)

ml3verse.kNN.result <- function(dat, split.ratio) {
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  task_plastic_type <- TaskClassif$new(id="plastic_type",
                                       backend = plastic_trn,
                                       target = "plastic_type")
  
  model_knn <- lrn("classif.kknn",kernel="rectangular",predict_type="prob")
  
  #model KNN
  model_knn$train(task = task_plastic_type)
  pre_knn_new <- model_knn$predict_newdata(newdata = plastic_tst)
  
  return(as.data.table(pre_knn_new))
} 

PCAtools_mergePC.kNNresult <- ml3verse.kNN.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.kNNresult <- ml3verse.kNN.result(e1071_merge_PC, split.ratio = 0.6)
my.data.kNN.result <- ml3verse.kNN.result(my.data, split.ratio = 0.6)

PCAtools_mergePC.kNNresult$row_ids <- rownames(plastic_tst)
e1071_merge_PC.kNNresult$row_ids <- rownames(plastic_tst)
my.data.kNN.result$row_ids <- rownames(plastic_tst)

plot.dat <- function(dat) {
  newdat <- as.data.frame(dat) %>%
    tidyr::pivot_longer(., cols = 4:ncol(.), names_to = "plastic_type", values_to = "prob")
  return(newdat)
}   

PCAtoolsdat <- plot.dat(PCAtools_mergePC.kNNresult)
PCAtoolsdat$plastic_type <- gsub("prob.", "",PCAtoolsdat$plastic_type)
e1071dat <- plot.dat(e1071_merge_PC.kNNresult)
e1071dat$plastic_type <- gsub("prob.", "",e1071dat$plastic_type)
mydat <- plot.dat(my.data.kNN.result)
mydat$plastic_type <- gsub("prob.", "",mydat$plastic_type)

# Make bar graph of prediction probability
ggplot(data = mydat) + 
  geom_col(aes(x = row_ids, y = prob, fill = plastic_type), 
           position = "dodge" # separating stacking prob cols
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(n.breaks = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))
