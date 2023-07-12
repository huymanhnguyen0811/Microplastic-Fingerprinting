library(e1071)
library(caret)
library(ROCR)

# split train and test data 60:40 , random sampling is done within the levels of plastic_type
samp <- caret::createDataPartition(as.factor(filled.data$plastic_type), p = 0.5, list = F)
train <- filled.data[samp,]
test <- filled.data[-samp,]

# Train and Select best cost parameter
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



