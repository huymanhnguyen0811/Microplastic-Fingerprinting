# With class package ---------------------------
library(class)
library(caret)


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
class.kNN.result <- function(dat, split.ratio){
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  default_knn_mod <- caret::train(
    plastic_type ~ .,
    data = plastic_trn,
    method = "knn",
    trControl = trainControl(method = "cv", number = 5), # using 5-fold cross-validation
    # preProcess = c("center", "scale"),
    tuneGrid = expand.grid(k = seq(1, 9, by = 2))
  )
  
  
  # Prediction results
  knn_pred_res <- predict(default_knn_mod, newdata = plastic_tst, type = "prob")
  rownames(knn_pred_res) <- rownames(plastic_tst)
  
  return(knn_pred_res)
}

PCAtools_mergePC.SVMresult <- class.kNN.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.SVMresult <- class.kNN.result(e1071_merge_PC, split.ratio = 0.6)


# With mlr3verse ---------------------
library(mlr3verse)
library(stats)
library(nnet)


inputdata <- shared_comp_plastic_type %>%
  dplyr::select(plastic_type, Percent_Area, Percent_Height, Ions, m.z, RT) %>%
  mutate(plastic_type = factor(plastic_type, levels = unique(plastic_type)))

task_plastic_type <- TaskClassif$new(id="plastic_type",
                                    backend = inputdata,
                                    target = "plastic_type")


model_knn <- lrn("classif.kknn",kernel="rectangular",predict_type="prob")
model_knn

model_multinomreg <- lrn("classif.multinom",predict_type="prob")
model_multinomreg

model_multinomreg_onevsrest <- mlr3pipelines::pipeline_ovr(model_multinomreg)

resample_plastic_cv = rsmp("cv",folds=10) 
set.seed(21193)
resample_plastic_cv$instantiate(task = task_plastic_type)

model_multinomreg$train(task = task_plastic_type) 
model_multinomreg$model


learner_plastic <- list(model_multinomreg_onevsrest,
                        model_knn,
                        model_multinomreg
)

design <- benchmark_grid(tasks = task_plastic_type,
                         learners = learner_plastic,
                         resamplings = resample_plastic_cv 
)


bmr = benchmark(design,store_models = TRUE)


result = bmr$aggregate(list(msr("classif.acc"),
                            msr("classif.bacc"),
                            msr("classif.mbrier")
)
)
result


#model KNN
model_knn$train(task = task_iris)
prediksi_knn_new <- model_knn$predict_newdata(newdata = iris_baru)
as.data.table(prediksi_knn_new)

