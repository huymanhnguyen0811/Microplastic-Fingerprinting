library(e1071)
library(caret)
library(ROCR)
library(MLmetrics)

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

# Fill NA with LOD ------------------------------------------------
filled.data <- copy(my.data)
for (c in 2:ncol(filled.data)) { 
  filled.data[which(base::is.na(filled.data[,c])), c] <- runif(length(which(base::is.na(filled.data[,c]))),
                                                               min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                               max = sort(shared_comp_plastic_type$Percent_Area)[2])
}

# split train and test data 60:40 , random sampling is done within the levels of plastic_type
samp <- caret::createDataPartition(filled.data$plastic_type, p = 0.5, list = F)
train <- filled.data[samp,]
test <- filled.data[-samp,]

# Train and Select best cost parameter -------------------
tune.out <- tune(svm, plastic_type ~ ., data = train, kernel = "polynomial",
                 ranges = list(cost = c(0.001 , 0.01 , 0.1 , 1, 5, 10, 100)),
                               # gamma = c(0.5 , 1, 2, 3, 4)),
                 decision.values = TRUE)
# summary(tune.out)
bestmod <- tune.out$best.model
View(table(bestmod$fitted, train$plastic_type))

# Test
ypred <- predict(bestmod, test)

View(table(ypred, test$plastic_type))

# Caret: Control params for SVM
ctrl <- caret::trainControl(
  method = "cv", 
  number = 10, 
  classProbs = TRUE,                 
  summaryFunction = multiClassSummary  # also needed for AUC/ROC
)

# Tune an SVM with caret
svm_auc <- caret::train(
  plastic_type ~ ., 
  data = train,
  method = "svmPoly",               
  preProcess = c("center", "scale"),
  metric = "ROC",
  trControl = ctrl,
  tuneLength = 10
)

# Print results
svm_auc$results

caret::confusionMatrix(svm_auc)



# Feature importance
prob_yes <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "Yes"]
}

# Variable importance plot
set.seed(2827)  # for reproducibility
vip(bestmod, method = "permute", nsim = 5, train = train, 
    target = "plastic_type", metric = "auc", reference_class = "Yes", 
    pred_wrapper = prob_yes)
