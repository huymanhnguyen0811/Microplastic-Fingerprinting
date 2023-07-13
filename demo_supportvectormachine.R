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

# e1071 package: ===========================

samp <- caret::createDataPartition(filled.data$plastic_type, p = 0.7, list = F)
train <- filled.data[samp,]
test <- filled.data[-samp,]


# Train and Select best cost parameter
tune.out <- e1071::tune(e1071::svm, plastic_type ~ ., data = train, kernel = "radial",
                 ranges = list(cost = c(0.001 , 0.01 , 0.1 , 1, 5, 10, 100),
                               gamma = seq(from=0, to=0.02, by= 0.005),
                               epsilon = c(0.001 , 0.01 , 0.1)),
                 decision.values = TRUE)

# summary(tune.out)
bestmod <- tune.out$best.model
View(table(bestmod$fitted, train$plastic_type))

# Test
ypred <- predict(bestmod, test)

View(table(ypred, test$plastic_type))

# Caret package: ===================================
# split train and test data 60:40 with stratification on plastic_type
set.seed(123)  # for reproducibility
# Conduct stratified sampling
churn_split <- rsample::initial_split(filled.data, prop = 0.7, strata = "plastic_type")
# Assign train & test set to variables
churn_train <- rsample::training(churn_split)
churn_test  <- rsample::testing(churn_split)

# Control params for SVM
ctrl <- caret::trainControl(
  method = "cv", 
  number = 10, # Performing 10-fold CrossValidation
  classProbs = TRUE,                 
  summaryFunction = multiClassSummary  # also needed for AUC/ROC
)

# Tune an SVM with caret
churn_svm_auc <- caret::train(
  plastic_type ~ ., 
  data = filled.data,
  method = "svmRadial",               
  # preProcess = c("center", "scale"),  
  metric = "ROC",  # area under ROC curve (AUC)       
  trControl = ctrl,
  tuneLength = 10
)

churn_svm_auc$results

# Confusion Matrix
caret::confusionMatrix(churn_svm_auc)

# Feature importance


vip_all <- c()
for (pt in unique(filled.data$plastic_type)) {
  set.seed(282)  # for reproducibility
  
  # Add prediction wrapper function to return the predicted class probabilities for the reference class of interest.
  prob <- function(object, newdata) {
    predict(object, newdata = newdata, type = "prob")[, pt] # Food_Packaging_Waste
  }
  
  # Variable importance plot
  vip_res <- vip::vip(churn_svm_auc, method = "permute", nsim = 5, train = filled.data, 
                      target = "plastic_type", metric = "auc", reference_class = pt,
                      pred_wrapper = prob)
  vip_all <- c(vip_all, vip_res$data[1:3,]$Variable)
}

# construct PDPs for the top four features of reference_class
features <- vip_all

pdps <- lapply(features, function(x) {
  partial(churn_svm_auc, pred.var = x, which.class = 2,  
          prob = TRUE, plot = TRUE, plot.engine = "ggplot2") +
    coord_flip()
})

grid.arrange(grobs = pdps,  ncol = 5)
