# library(R.utils)
`%notin%` <- Negate(`%in%`)

# Metadata of compounds that appear in at least 2 samples ===================
stats_rq1 <- shared_comp_sample %>%
  select(File, collapsed_compound, Percent_Area) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = File, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

# Fill in missing value with LOD
for (r in 1:nrow(stats_rq1)) { 
  stats_rq1[r, which(base::is.na(stats_rq1[r,]))] <- runif(length(which(base::is.na(stats_rq1[r,]))),
                                                           min = sort(shared_comp_sample$Percent_Area)[1],
                                                           max = sort(shared_comp_sample$Percent_Area)[2])
}

#  Metadata of compounds that appear in at least 2 plastic types ===============
stats_rq1b <- shared_comp_plastic_type %>%
  select(File, collapsed_compound, Percent_Area) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = File, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

# Fill in missing value with LOD
for (r in 1:nrow(stats_rq1b)) { 
  stats_rq1b[r, which(base::is.na(stats_rq1b[r,]))] <- runif(length(which(base::is.na(stats_rq1b[r,]))),
                                                           min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                           max = sort(shared_comp_plastic_type$Percent_Area)[2])
}


# Reconstruct data frame
transpose_df <- data.table::transpose(stats_rq1b)
rownames(transpose_df) <- colnames(stats_rq1b)
colnames(transpose_df) <- rownames(stats_rq1b)

transpose_df <- transpose_df %>%
  tibble::rownames_to_column(., var = "File") %>%
  dplyr::mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                               ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                      ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                             ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                    ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                           ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
  relocate(plastic_type, .after = File)

# Extract balloons from transpose data frame and make them class numeric 
balloons <- as.vector(t(transpose_df[1:18, 3:ncol(transpose_df)]))
hist(balloons, col='steelblue', main='Balloons')

### Examining for normal distribution ----------------
## Histogram 
# Create empty list 
histlist <- list()
# For loop creating histograms of each plastic type
i <- 1
for (plastic in unique(transpose_df$plastic_type)) {
  grouped_plastic <- which(transpose_df$plastic_type == plastic)
  histlist[[i]] <- ggplot() +
    geom_histogram(mapping = aes(as.vector(t(transpose_df[grouped_plastic, 3:ncol(transpose_df)])))) + 
    labs(x = plastic, title = plastic)
  i <- i + 1
}
# Appending each histogram to the empty list in order to view all of them
gridExtra::grid.arrange(grobs = histlist)

## QQ Plot 
# Create empty list 
qqlist <- list()
# For loop creating QQ plots of each plastic type
i <- 1 
# Appending QQ plots together 
par(mfrow=c(3,3))
for (plastic in unique(transpose_df$plastic_type)) {
  grouped_plastic <- which(transpose_df$plastic_type == plastic)
  qqlist[[i]] <- qqnorm(as.vector(t(transpose_df[grouped_plastic, 3:ncol(transpose_df)])), main = plastic, xlab = plastic, col = 'steelblue')
  qqline(as.vector(t(transpose_df[grouped_plastic, 3:ncol(transpose_df)])), col = 'red')
  i <- i + 1
}

### Conclusion = they are not normally distributed 



### Examining for equality of variance -----------------------------
# Levene’s test-non-normally distributed data (car::leveneTest()), significant if p-value < 0.05 
library(car)

results <- list()
for (col in 1:ncol(utils::combn(unique(transpose_df$plastic_type), 2))) {
  # extract the combinations of plastic type pairs
  plastic_type_1 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][1]
  plastic_type_2 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][2]
  idx1 <- which(transpose_df$plastic_type == plastic_type_1)
  idx2 <- which(transpose_df$plastic_type == plastic_type_2)
  V1 <- as.vector(t(transpose_df[idx1, 3:ncol(transpose_df)]))
  V2 <- as.vector(t(transpose_df[idx2, 3:ncol(transpose_df)]))
  datat <- c(V1, V2)
  grouped <- as.factor(c(rep(plastic_type_1, times = length(V1)), rep(plastic_type_2, times = length(V2))))
  non_norm <- data.frame(datat, grouped)
  results[[paste0(plastic_type_1, plastic_type_2)]] <- car::leveneTest(datat ~ grouped, data = non_norm)
}
# -> p-value is 0.4096 which is greater the significance level of 0.05 therefore can conclude there is no significant difference between the variances


# Fligner-Killeen's test (fligner.test()), significant if p-value < 0.05 
fligner.test(datat ~ grouped, data = non_norm)
# p-value less than significance level therefore there is significant difference between variances.
# Test to determine the homogeneity of group variances. 


# K-S Test non-parametric data ---------
# The KS test compared the CDF of these samples. The null hypothesis is that the two distributions are the same. 

stats::ks.test(V1, V2)
### With alternative "two-sided" which is default, p-value less than significance value. We can reject the null hypothesis and claim the two 
#distributions are not the same. D = 0.05 which is close to 0, however, means the two samples' distribution are pretty similar even if they 
#are not the same. The D statistic is a calculation of the maximum difference between the two samples' empirical distribution 
#functions.

# KS Test asks if x and y come from the same distribution. KS test measures the spread and shape of the data sets, 
#i.e. if the data is symmetric or skewed. The Wilcoxon test is interested in the location differences of the data sets, not so much 
#the spread of the data but where individually are their medians and means comparatively. 

# Wilcoxon test to verify non normality of samples -----------
# Examine a difference in median value between distributions. 
wilcox.test(datat ~ grouped, data = non_norm, paired = FALSE)
# -> p-value smaller than significance level meaning two samples differ in distribution and median


summary_table_of_tests <- data.frame(comp=character(), pair=character(), wilcox_test=integer(), ks_test=integer())

for (comp in colnames(transpose_df[,3:ncol(transpose_df)])) {
  # looping through the combinations of plastic types 
  for (col in 1:ncol(utils::combn(unique(transpose_df$plastic_type), 2))) {
    # extract the combinations of plastic type pairs
    plastic_type_1 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][1]
    plastic_type_2 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][2]
    # calculating the p-value between each plastic type pair 
    wilcox_test <- wilcox.test(transpose_df[which(transpose_df$plastic_type == plastic_type_1), comp], 
                               transpose_df[which(transpose_df$plastic_type == plastic_type_2), comp])$p.value
    ks_test <- ks.test(transpose_df[which(transpose_df$plastic_type == plastic_type_1), comp],
                       transpose_df[which(transpose_df$plastic_type == plastic_type_2), comp])$p.value
    # assigning row information
    summary_table_of_tests[nrow(summary_table_of_tests) + 1,] <- c(comp, paste0(plastic_type_1, "&", plastic_type_2), wilcox_test, ks_test)
  }
}

summary_table_of_tests$adjusted_pvalue_holm <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "holm")
summary_table_of_tests$adjusted_pvalue_hochberg <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "hochberg")
summary_table_of_tests$adjusted_pvalue_hommel <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "hommel")
summary_table_of_tests$adjusted_pvalue_bonferroni <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "bonferroni")
summary_table_of_tests$adjusted_pvalue_BH <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "BH")
summary_table_of_tests$adjusted_pvalue_BY <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "BY")
summary_table_of_tests$adjusted_pvalue_fdr <- stats::p.adjust(summary_table_of_tests$wilcox_test, method = "fdr")

alpha <- summary_table_of_tests %>%
  filter(., adjusted_pvalue_holm < 0.05) %>%
  arrange(adjusted_pvalue_holm)

summary_p0.05 <- summary_table_of_tests %>% filter(wilcox_test < 0.05 & ks_test < 0.05)

unique(summary_p0.05$comp)

top <- summary_p0.05 %>% arrange(wilcox_test, ks_test)

top_wil <- summary_p0.05 %>% arrange(wilcox_test)
top_ks <- summary_p0.05 %>% arrange(ks_test)


#get first 10 rows
top_10 <- top[1:10,]
#then make ggplot with unique of each compound 

C16 <- transpose_df %>% select(., c("plastic_type", "Compound_16.")) %>% filter(plastic_type %in% c("Plastic_Cups", 
                                                                                                    "Plastic_Drinking_Straws"))
ggplot(data = C16) + geom_boxplot(aes(x = plastic_type, y = Compound_16., fill = plastic_type))

C1164 <- transpose_df %>% select(., c("plastic_type", "Compound_1164.")) %>% filter(plastic_type %in% c("Food_Packaging_Waste", 
                                                                                                    "Plastic_Cups"))
ggplot(data = C1164) + geom_boxplot(aes(x = plastic_type, y = Compound_1164., fill = plastic_type))



result <- c()
for (i in unique(summary_p0.05$comp)) {
  if (length(which(summary_p0.05$comp == i)) > 4) {
    result <- c(result, i)
  }
  else {
    next
  }
}


ggplot(data = transpose_df %>% 
         dplyr::select(., c("plastic_type", result)) %>%
         pivot_longer(., cols = 2:ncol(.), names_to = "compound", values_to = "percent_area"),
       aes(x = plastic_type, y = percent_area, fill = plastic_type)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ compound) +
  geom_boxplot()
  


# 5 assumptions of multiple linear regression --------------------
must_convert<-sapply(transpose_df, is.factor)
new_pt <- sapply(transpose_df[, must_convert], unclass)
out<-cbind(transpose_df[,!must_convert],new_pt) 

out <- out %>% 
  relocate(new_pt, .after = File) %>%
  column_to_rownames(., var = "File") 

# Linear relationship (predicting the plastic type)
linearity <- lm(new_pt ~ ., data = out)
summary(linearity)

 # No multicollinearity 

library(car)
vif(linearity)

  # check correlation between predictor variables
cor_matrix <- out[, 2:ncol(out)]

cor(cor_matrix)
  # create heat map to see correlation 
library(reshape2)
cor_mat <- round(cor(cor_matrix), 3)
melt_cor_mat <- melt(cor_mat)

#ggplot(data = melt_cor_mat, aes(x=Var1, y=Var2, fill=value)) +
 # geom_tile()

install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot::ggcorrplot(cor(cor_matrix))

 # Independence 
  # the Durbin-Watson Test
library(car)
durbinWatsonTest(linearity)
  # reject null hypothesis because autocorrelation is not 0,
  # alternative hypothesis is true (true autocorrelation is greater than 0).

 # Homoscedasticity 
plot(fitted(linearity), resid(linearity), xlab='Fitted Values', ylab='Residuals')
abline(0,0)
   # does not work because lm model is filled with NA's, must fix this. 
   # all compounds that had a number for st deviation are important. must remove all NA's (because of collinearity). 
  # slice into coefficients in linearity and take out all compounds where is.na = FALSE because we want all the values. 
  # go back to out matrix and subset all the columns that are important 

 # Multivariate Normality 



# Multivariate parametric test -> anova (assumptions must be met, earlier tests were non-parametric)




# Examine whether 2 samples have equal variance and normally distributed.-----------------------------------------------
## Normally distributed by histogram, Q-Q plots and Shapiro-Wilk test, if p-value > 0.05 
#-> normally distributed, not for each compound but for the whole
GasData <- as.vector(t(cat_5[,c(1:21)])) # 100149 data points
DieselData <- as.vector(t(cat_5[,c(22:25)])) # 19076 data point



# Histogram
hist(GasData, col='steelblue', main='Gas')
hist(DieselData, col='steelblue', main='Diesel')

# Q-Q plots aka. Normal Probability plots
stats::qqnorm(GasData, main='Gas')
stats::qqline(GasData)

stats::qqnorm(DieselData, main='Diesel')
stats::qqline(DieselData)


# CONCLUSION: both Gas and Diesel populations are NOT normally distrubuted

# Equality of variance between Gas and Diesel populations by (https://www.r-bloggers.com/2021/06/equality-of-variances-in-r-homogeneity-test-quick-guide/)

# F-test-normally distributed (var.test()); - unusable b/c data non-normally distributed
# Bartlett’s test-normally distributed (bartlett.test()) - unusable b/c data non-normally distributed

# Levene’s test-non-normally distributed data (car::leveneTest()); p-val < 0.05 
# -> there is significant difference between the tested sample variances => Equality of variance is NOT satisfied
library(car)
data <- c(GasData, DieselData)
group <- as.factor(c(rep("Gas", times = length(GasData)), rep("Diesel", times = length(DieselData))))
non_norm_dist_data <- data.frame(data, group)
car::leveneTest(data ~ group, data = non_norm_dist_data)

# Fligner-Killeen test (fligner.test()) - p-value < 0.05 
# -> there is significant difference between tested sample variance => Equality of variance is NOT satisfied
stats::fligner.test(data ~ group, data = non_norm_dist_data)

  
# Wilcoxon Test ========================================================================================
# !! vary sample size gene/compound dataset, three solutions: exact wilcoxon, permutation test, SAM (stanford)

# create empty data.frame
summary_table <- data.frame(collapsed_compound = character(), pair=character(), pvalue=integer())

# for loop through each compound
for (comp in colnames(transpose_df[,3:ncol(transpose_df)])) {
  # looping through the combinations of plastic types 
  for (col in 1:ncol(utils::combn(unique(transpose_df$plastic_type), 2))) {
    # extract the combinations of plastic type pairs
    plastic_type_1 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][1]
    plastic_type_2 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][2]
    # calculating the p-value between each plastic type pair 
    p.value.res <- wilcox.test(transpose_df[which(transpose_df$plastic_type == plastic_type_1), comp], 
                               transpose_df[which(transpose_df$plastic_type == plastic_type_2), comp])$p.value
    # assigning row information
    summary_table[nrow(summary_table) + 1,] <- c(comp, paste0(plastic_type_1, "_", plastic_type_2), p.value.res)
  }
}

summary_table$adjusted_pvalue_holm <- stats::p.adjust(summary_table$pvalue, method = "holm")
summary_table$adjusted_pvalue_hochberg <- stats::p.adjust(summary_table$pvalue, method = "hochberg")
summary_table$adjusted_pvalue_hommel <- stats::p.adjust(summary_table$pvalue, method = "hommel")
summary_table$adjusted_pvalue_bonferroni <- stats::p.adjust(summary_table$pvalue, method = "bonferroni")
summary_table$adjusted_pvalue_BH <- stats::p.adjust(summary_table$pvalue, method = "BH")
summary_table$adjusted_pvalue_BY <- stats::p.adjust(summary_table$pvalue, method = "BY")
summary_table$adjusted_pvalue_fdr <- stats::p.adjust(summary_table$pvalue, method = "fdr")

alpha <- summary_table %>%
  filter(., adjusted_pvalue_holm < 0.1) %>%
  arrange(adjusted_pvalue_holm)



# Examine whether groups have equal variance and normally distributed.-----------------------------------------------
## Normally distributed by histogram, Q-Q plots and Shapiro-Wilk test, if p-value > 0.05 
#-> normally distributed, not for each compound but for the whole
gas138_Data <- as.vector(t(df_stats_rq2_rt10.1 %>% 
                             filter(., gas_station %in% c("Station_1", "Station_3", "Station_8")) %>% 
                             select(., 3:ncol(.))))
gas579_Data <- as.vector(t(df_stats_rq2_rt10.1 %>%
                             filter(., gas_station %in% c("Station_5", "Station_7", "Station_9")) %>% 
                             select(., 3:ncol(.))))

# Histogram
hist(gas138_Data, col='steelblue', main='Gas stations 1, 3, 8')
hist(gas579_Data, col='steelblue', main='Gas stations 5, 7, 9')

# Q-Q plots aka. Normal Probability plots
stats::qqnorm(gas138_Data, main='Gas stations 1, 3, 8')
stats::qqline(gas138_Data)

stats::qqnorm(gas579_Data, main='Gas stations 5, 7, 9')
stats::qqline(gas579_Data)

# Shapiro-Wilk tests - Error in shapiro.test(GasData) : sample size must be between 3 and 5000 -> do Normal Probability plot instead=
# shapiro.test(GasData)
# shapiro.test(DieselData)

# CONCLUSION: both Gas station group 1 (1, 3, 8) and Gas station group 2 (5, 7, 9) populations are NOT normally distrubuted

# Equality of variance between Gas and Diesel populations by (https://www.r-bloggers.com/2021/06/equality-of-variances-in-r-homogeneity-test-quick-guide/)

# F-test-normally distributed (var.test()); - unusable b/c data non-normally distributed
# Bartlett’s test-normally distributed (bartlett.test()) - unusable b/c data non-normally distributed

# Levene’s test-non-normally distributed data (car::leveneTest()); p-val > 0.05 
# -> there is no significant difference between the tested sample variances => Equality of variance is satisfied
library(car)
data <- c(gas138_Data, gas579_Data)
group <- as.factor(c(rep("Gas138", times = length(gas138_Data)), rep("Gas579", times = length(gas579_Data))))
non_norm_dist_data <- data.frame(data, group)
car::leveneTest(data ~ group, data = non_norm_dist_data)

# Fligner-Killeen test (fligner.test()) - p-value < 0.05 
# -> there is significant difference between tested sample variance => Equality of variance is NOT satisfied
stats::fligner.test(data ~ group, data = non_norm_dist_data)