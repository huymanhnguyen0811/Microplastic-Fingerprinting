library(e1071)

library(ISLR2)
# Khan Gene Data: consists of tissue samples (as rows of xtrain and xtest) corresponding to four distinct types of small round blue cell tumors.
# For each tissue sample, 2308 gene expression measurements are available
names(Khan)
table(Khan$ytrain)
table(Khan$ytest)
