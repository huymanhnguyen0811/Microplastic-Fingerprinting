library(class)
library(caret)
library(nnet)

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
  
  plastic_multinomlog_mod <- caret::train(
    plastic_type ~ .,
    data = plastic_trn,
    method = "multinom", # Penalized Multinomial Regression -> https://remiller1450.github.io/s230f19/caret3.html
    trControl = trainControl(method = "repeatedcv", number = 5, 
                             repeats = 10, verboseIter = TRUE),
    trace = FALSE,
    MaxNWts = 8 * (ncol(plastic_trn)) # maximum allowable number of weights
  )
  
  # Prediction results
  pred_res <- round(predict(plastic_multinomlog_mod, newdata = plastic_tst, type = "prob"), 5)

  return(pred_res)
}

PCAtools_mergePC.multinomlogresult <- caret.multinomlog.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.multinomlogresult <- caret.multinomlog.result(e1071_merge_PC, split.ratio = 0.6)
my.data.multinomlogresult <- caret.multinomlog.result(my.data, split.ratio = 0.6)

plot.dat <- function(dat) {
  newdat <- as.data.frame(dat) %>%
    tibble::rownames_to_column(., var = "File") %>%
    tidyr::pivot_longer(., cols = 2:ncol(.), names_to = "plastic_type", values_to = "prob")
  return(newdat)
}   

PCAtoolsdat <- plot.dat(PCAtools_mergePC.multinomlogresult)
e1071dat <- plot.dat(e1071_merge_PC.multinomlogresult)
mydat <- plot.dat(my.data.multinomlogresult)

# Make bar graph of prediction probability
ggplot(data = e1071dat) + 
  geom_col(aes(x = File, y = prob, fill = plastic_type), 
           position = "dodge" # separating stacking prob cols
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(n.breaks = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))
