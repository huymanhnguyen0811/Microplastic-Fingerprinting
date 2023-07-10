# Check which compounds are positively correlated and combine them into a compound group
# https://stats.stackexchange.com/questions/116853/combining-merging-correlated-variables - Factor analysis
# https://socialsciences.mcmaster.ca/jfox/Misc/sem/SEM-paper.pdf - Structural Equation Modelling


library(sem)


Klein$P.lag <- c(NA, Klein$P[-22])
Klein$X.lag <- c(NA, Klein$X[-22])
