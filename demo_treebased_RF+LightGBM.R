# RandomForestSRC ----------------
library(randomForestSRC)

set.seed (1234) # sets a numerical starting point; will be set randomly if not set by the user

# WARNING!!! -> Values from `Percent_Area` are not uniquely identified -> ASK ROXANA!!!!
View(shared_comp_plastic_type %>%
       dplyr::group_by(File, collapsed_compound) %>%
       dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
       dplyr::filter(n > 1L))

# since we have duplicates with different values of the same compound in some samples
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

# Not fill NA with LOD --------------------------------------------
my.rf <- rfsrc(plastic_type ~ ., ntree=2000, splitrule = "auc", 
               na.action = "na.impute", 
               # nimpute = 2,
               nodesize = 1, importance = "random", data=my.data) # function to run RF

my.rf # output of model details, e.g. goodness-of-fit, out-of-bag (OBB) error

# Fill NA with LOD ------------------------------------------------
filled.data <- copy(my.data)
for (c in 2:ncol(filled.data)) { 
  filled.data[which(base::is.na(filled.data[,c])), c] <- runif(length(which(base::is.na(filled.data[,c]))),
                                                               min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                               max = sort(shared_comp_plastic_type$Percent_Area)[2])
}

filled.rf <- rfsrc(plastic_type ~ ., ntree=2000, splitrule = "auc", 
                   nodesize = 1,
                   importance = "random", data = filled.data) # function to run RF

filled.rf

# PCA-based feature data ==========================
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
rfsrc.result <- function(dat, split.ratio){
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$plastic_type, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  mergePC.rf <- rfsrc(plastic_type ~ ., ntree=2000, splitrule = "auc", 
                      nodesize = 1,
                      importance = "random", data = plastic_trn) # function to run RF
  
  mergePC.rf
  
  # Prediction results
  pred_res <- predict(mergePC.rf, newdata = plastic_tst, type = "prob")$predicted
  rownames(pred_res) <- rownames(plastic_tst)
  
  # selection of the best feature candidates ------------------
  md.obj <- max.subtree(mergePC.rf)
  best.feature <- md.obj$topvars # extracts the names of the variables in the object md.obj
  
  return(list(pred_res, best.feature))
}

PCAtools_mergePC.RFresult <- rfsrc.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.RFresult <- rfsrc.result(e1071_merge_PC, split.ratio = 0.6)


# plot OOB error rate against the number of trees -------
plot(ggRandomForests::gg_error(my.rf)) 
plot(ggRandomForests::gg_error(filled.rf)) 

# Estimate the variables importance --------
my.rf.vimp <- ggRandomForests::gg_vimp(mergePC.rf, nvar = 100) # provides the predictor's importance of top 100 predictors
plot(my.rf.vimp) # visualises the predictor’s importance
 
# Plot the response variable against each predictor variable, we can generate partial dependance plots --------------
my.rf.part.plot <- plot.variable(mergePC.rf, partial=TRUE, sorted=FALSE,
                                 show.plots=FALSE, 
                                 nvar = 10 # look at top 10 predictor's importance
                                 )
gg.part <- ggRandomForests::gg_partial(my.rf.part.plot)
plot(gg.part, xvar=names(plastic_trn[,-1]), panel=TRUE, se=TRUE)



# my.rf.interaction <- find.interaction(my.rf, xvar.names=md.obj$topvars,
                                      # importance="random", method="vimp", nrep=3) 


# LightGBM for Multiclass Classification ==================
library(lightgbm)

lightgbm.result <- function(dat, split.ratio) {
  
  # Must convert plastic_type from factors to numeric
  test <- copy(dat)
  test$plastic_type <- as.numeric(as.factor(test$plastic_type)) - 1L
  # Split train and test sets
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(test$plastic_type, p = split.ratio, list = F)
  plastic_trn <- as.matrix(test[plastic_idx, ])
  plastic_tst <- as.matrix(test[-plastic_idx, ])
  
  dtrain <- lgb.Dataset(data = plastic_trn[, 2:ncol(plastic_trn)], label = plastic_trn[, 1])
  dtest <- lgb.Dataset.create.valid(dtrain, data = plastic_tst[, 2:ncol(plastic_tst)], label = plastic_tst[, 1])
  valids <- list(test = dtest)
  
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
                    nrounds = 100,
                    valids,
                    early_stopping_rounds = 25L
                    )
  
  # prediction
  my_preds <- round(predict(lgb.model, plastic_tst[, 2:ncol(plastic_tst)]), 3)
  smaller_my_preds <- split(my_preds, rep(1:29, each = 7))
  # create prediction result df
  testing <- data.frame(matrix(nrow=0, ncol = 7))
  for (i in 1:length(smaller_my_preds)) {
    testing <- rbind(testing, smaller_my_preds[[i]])
  }

  colnames(testing) <- levels(dat$plastic_type)
  rownames(testing) <- rownames(plastic_tst)

  # IMportant features
  tree_imp = lgb.importance(lgb.model, percentage = T)
  imp_var <- lgb.plot.importance(tree_imp, measure = "Gain")
  
  return(list(testing, imp_var))
}


PCAtools_mergePC.LGBMresult <- lightgbm.result(PCAtools_mergePC, split.ratio = 0.6)
e1071_merge_PC.LGBMresult <- lightgbm.result(e1071_merge_PC, split.ratio = 0.6)

View(PCAtools_mergePC.LGBMresult[[1]])
View(e1071_merge_PC.LGBMresult[[1]])

# Reserve Code -------------------
library(tree)

# Fit classification tree to predict plastic_type 
tree.model <- tree::tree(plastic_type ~ .,shared_comp_normalized)
summary(tree.model)

plot(tree.model)
# pretty = 0 instructs R to include the category names for any qualitative predictors, rather than simply displaying a letter for each category.
text(tree.model, pretty = 0)


# Does it make sense to use ONLY collapsed compounds as predictor variables, 
# if yes then what is the data frame needed for randomforest 

set.seed(7)
train <- sample(1: nrow(shared_comp_normalized), 1500)
dat.test <- shared_comp_normalized[-train, ]
plastic.true <- shared_comp_normalized$plastic_type[-train]
tree.plastic <- tree(plastic_type ~ . ,shared_comp_normalized,
                     subset = train)
tree.pred <- predict(tree.plastic, dat.test,
                     type = "class")
table(tree.pred, plastic.true)

# Improve with tree prunning
set.seed (7)
cv.plastic <- cv.tree(tree.plastic, FUN = prune.misclass)
plot(cv.plastic$size , cv.plastic$dev , type = "b")
plot(cv.plastic$k, cv.plastic$dev , type = "b")
prune <- prune.misclass(tree.plastic, best = 9)
plot(prune)
text(prune, pretty = 0)

# How many predictor variables in total? 

# Calculate MSE 
mean((yhat - boston.test)^2)

# Bagging 
bag.plas <- randomForest(plastic_type ~ .,data = shared_comp_normalized, subset = train, mtry = 12, importance = TRUE)
yhat.bag <- predict(bag.plas, newdata = shared_comp_normalized[-train , ])
plas.test <- shared_comp_normalized[-train , "plastic_type"]
plot(yhat.bag, plas.test)
mean((yhat.bag - plas.test)^2)

bag.plas <- randomForest(plastic_type ~ ., data = shared_comp_normalized,
                         subset = train, mtry = 12, ntree = 25)
yhat.bag <- predict(bag.plas, newdata = shared_comp_normalized[-train , ])
mean((yhat.bag - plas.test)^2)

# Random Forests 
rf.plastic <- randomForest(plastic_type ~ Percent_Area+m.z+Ions+RT+collapsed_compound, data = shared_comp_normalized,
                           subset = train, importance = TRUE)
yhat.rf <- predict(rf.plastic, newdata = shared_comp_normalized[-train , ])
mean((yhat.rf - plas.test)^2)
importance(rf.plastic) 
# %IncMSE: mean decrease of accuracy in predictions on the out of bag samples when a given variable is permuted
# IncNodePurity: measure of the total decrease in node impurity that results from splits over that variable, averaged over all trees
varImpPlot(rf.plastic)

# Boosting
library(gbm)
boost.plastic <- gbm(medv ~ ., data = Boston[train, ],
                     distribution = "gaussian", n.trees = 5000,
                     interaction.depth = 4)
summary(boost.boston)

# PArtial dependence plot: illustrate the marginal effect of the selected variables on the response after integrating out the other variables
# median house prices are increasing with rm and decreasing with lstat.
plot(boost.boston, i = "rm")
plot(boost.boston, i = "lstat")

yhat.boost <- predict(boost.boston ,
                      newdata = Boston[-train , ], n.trees = 5000)

# Tuning Shrinkage parameter λ 
boost.boston <- gbm(medv ~ ., data = Boston[train, ],
                    distribution = "gaussian", n.trees = 5000,
                    interaction.depth = 4, shrinkage = 0.2, verbose = F)
yhat.boost <- predict(boost.boston ,
                      newdata = Boston[-train , ], n.trees = 5000)
mean((yhat.boost - boston.test)^2)

# Bayesian Additive Regression Trees
library(BART)
x <- Boston[, 1:12]
y <- Boston[, "medv"]
xtrain <- x[train , ]
ytrain <- y[train]
xtest <- x[-train , ]
ytest <- y[- train]
set.seed(1)
bartfit <- gbart(xtrain , ytrain , x.test = xtest)

yhat.bart <- bartfit$yhat.test.mean
mean((ytest - yhat.bart)^2)
# check how many times each variable appeared in the collection of trees.
ord <- order(bartfit$varcount.mean , decreasing = T)
bartfit$varcount.mean[ord]

