library(tree)

library ( ISLR2 )
attach ( Carseats )
High <- factor(ifelse(Sales <= 8, "No", " Yes "))
Carseats <- data.frame (Carseats , High)

tree.model <- tree::tree(High ~ . -Sales,  Carseats)
summary(tree.model)
tree.model

plot(tree.model)
# pretty = 0 instructs R to include the category names for any qualitative predictors, rather than simply displaying a letter for each category.
text(tree.model, pretty = 0)


# Does it make sense to use ONLY collapsed compounds as predictor variables, 
# if yes then what is the data frame needed for randomforest 

set.seed(7)
train <- sample(1: nrow(Carseats), 200)
Carseats.test <- Carseats[-train, ]
High.test <- High[- train]
tree.carseats <- tree(High ~ . - Sales ,Carseats,
                      subset = train )
tree.pred <- predict(tree.carseats, Carseats.test ,
                     type = "class")
table(tree.pred, High.test)

# Improve with tree prunning
set.seed (7)
cv.carseats <- cv.tree (tree.carseats, FUN = prune.misclass)
names(cv.carseats )