# Check which compounds are positively correlated and combine them into a compound group
# https://stats.stackexchange.com/questions/116853/combining-merging-correlated-variables - Factor analysis
# https://socialsciences.mcmaster.ca/jfox/Misc/sem/SEM-paper.pdf - Structural Equation Modelling


library(lavaan)

HS.model <- ' visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 '

fit <- lavaan::cfa(model = HS.model, data = HolzingerSwineford1939)
summary(fit, fit.measures = TRUE)
