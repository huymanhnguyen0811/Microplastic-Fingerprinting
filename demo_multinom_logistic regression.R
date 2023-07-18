library(class)
library(caret)
library(nnet)


my.data <- shared_comp_plastic_type  %>%
  dplyr::select(File, collapsed_compound, Percent_Area) %>%
  group_by(File, collapsed_compound) %>%
  # , we summarize these values by taking the mean of them
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                               ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                      ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                             ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                    ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                           ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
  mutate(plastic_type = factor(plastic_type, levels = unique(plastic_type))) %>%
  relocate(plastic_type, .before = 1) %>%
  column_to_rownames(., var = "File")

#k-NN cannot handle missing values
for (c in 2:ncol(my.data)) { 
  my.data[which(base::is.na(my.data[,c])), c] <- runif(length(which(base::is.na(my.data[,c]))),
                                                       min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                       max = sort(shared_comp_plastic_type$Percent_Area)[2])
}

# With PCA-based feature data
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
caret.multinomlog.result <- function(dat, split.ratio){
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  plastic_multinomlog_mod <- train(
    plastic_type ~ .,
    data = plastic_trn,
    method = "multinom", # Penalized Multinomial Regression
    trControl = trainControl(method = "cv", number = 5),
    trace = FALSE,
    MaxNWts = 8 * (ncol(plastic_trn)) # maximum allowable number of weights
  )
  
  # Prediction results
  pred_res <- round(predict(plastic_multinomlog_mod, newdata = plastic_tst, type = "prob"), 3)

  return(pred_res)
}

PCAtools_mergePC.SVMresult <- caret.multinomlog.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.SVMresult <- caret.multinomlog.result(e1071_merge_PC, split.ratio = 0.6)


