library(bnlearn)
library(klaR)
library(e1071)
library(naivebayes)


samp <- caret::createDataPartition(filled.data$plastic_type, p = 0.7, list = F)
train <- filled.data[samp,]
test <- filled.data[-samp,]

# with e1071 package ------------------
NB_e1071 <- e1071::naiveBayes(plastic_type ~ .,
                  data = train,
                  laplace = 0 # without Laplace smoothing
                  )

predict(NB, test,
        type = "raw")

# with naivebayes package -----------------
NB_naivebayes <- naivebayes::multinomial_naive_bayes(x = as.matrix(train[, -1]),
                                                     y = train[, 1])

predict(NB_naivebayes, newdata = as.matrix(test[, -1]), type = "class")
