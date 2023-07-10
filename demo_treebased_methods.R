# RandomForestSRC ----------------
library(randomForestSRC)

my.data <- my_data %>%
  rownames_to_column(., var = "File") %>%
  mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                               ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                      ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                             ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                    ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                           ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
  relocate(plastic_type, .before = 1) %>%
  column_to_rownames(., var = "File") %>%
  mutate(plastic_type = factor(plastic_type, levels = unique(plastic_type)))

set.seed (1234) # sets a numerical starting point; will be set randomly if not set by the user

my.rf <- rfsrc(plastic_type ~ ., mtry=5, ntree=2000, splitrule = "gini", na.action = "na.impute",
               importance = "random", data=my.data) # function to run RF

my.rf # output of model details, e.g. goodness-of-fit, out-of-bag (OBB) error

plot(ggRandomForests::gg_error(my.rf)) # plot OOB error rate against the number of trees

# Estimate the variables importance --------
my.rf.vimp <- ggRandomForests::gg_vimp(my.rf, nvar = 100) # provides the predictor's importance of top 100 predictors
my.rf.vimp
plot(my.rf.vimp) # visualises the predictor’s importance
 
# Plot the response variable against each predictor variable, we can generate partial dependance plots --------------
my.rf.part.plot <- plot.variable(my.rf, partial=TRUE, sorted=FALSE,
                                  show.plots=FALSE, nvar = 10)
gg.part <- gg_partial(my.rf.part.plot)
plot(gg.part, xvar=names(my.data[,-1]), panel=TRUE, se=TRUE)

# selection of the best candidates ------------------
md.obj <- max.subtree(my.rf)
md.obj$topvars # extracts the names of the variables in the object md.obj

my.rf.interaction <- find.interaction(my.rf, xvar.names=md.obj$topvars,
                                      importance="random", method="vimp", nrep=3) 

# Boosted Regression Trees
library(gbm)
library(dismo)

my.brt <- gbm.step(data=my.data, gbm.x=c(1:3, 5), gbm.y=6,
                    family="gaussian", tree.complexity=5, learning.rate=0.005,
                    bag.fraction=0.7)



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

# LightGBM for Multiclass Classification
library(lightgbm)
# split into train and test  
indexes = caret::createDataPartition(shared_comp_normalized$plastic_type, p = .8, list = F)

# split dataset into train and test parts that both contain feature and label parts.
train = as.matrix(shared_comp_normalized[-indexes, ])
test = as.matrix(shared_comp_normalized[indexes, ])

train_x = train[, -9]
train_y = train[, 9]

test_x = test[, -9]
test_y = test[, 9]

# load the train and test data into the LightGBM dataset object
dtrain = lightgbm::lgb.Dataset(train_x, label = train_y)
dtest = lightgbm::lgb.Dataset.create.valid(dtrain, data = test_x, label = test_y)

# Setup parameters
params = list(
  objective= 'multiclass',
  metric = "multi_error", 
  num_class = 7) 

valids = list(test = dtest)
model = lgb.train(params,
                  dtrain,
                  nrounds = 100,
                  valids,
                  min_data=1,
                  learning_rate = 1,
                  early_stopping_rounds = 10)

# prediction
pred = predict(model, test_x, reshape=T)
pred_y = max.col(pred)-1

# Confusion matrix
confusionMatrix(as.factor(test_y), as.factor(pred_y))

tree_imp = lgb.importance(model, percentage = T)
lgb.plot.importance(tree_imp, measure = "Gain")

