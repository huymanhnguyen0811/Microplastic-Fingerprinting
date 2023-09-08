# library(R.utils)
`%notin%` <- Negate(`%in%`)

# Metadata of compounds that appear in at least 2 samples HPLCTOFMS ===================
stats_rq1 <- gc_hplc %>%
  select(File, collapsed_compound, product_cat, Values) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, product_cat, collapsed_compound) %>%
  summarise(across(Values, mean)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Values) %>%
  column_to_rownames(., var = "File") 

select_col <- c()
for (c in 2:ncol(stats_rq1)) {
  non_na <- c()
  for (prod in unique(stats_rq1$product_cat)) {
    non_na <- c(non_na, sum(!is.na(stats_rq1[which(stats_rq1$product_cat == prod), c])))
  }
  if (sum(non_na > 2) > 2) {
    # print(paste0("non_na", non_na))
    next
  } else {
    # print(paste0("c", c))
    select_col <- c(select_col, c)
  }
}

new_stats <- stats_rq1[, select_col]
new_stats$product_cat <- stats_rq1$product_cat

# Fill in missing value with LOD
for (r in 1:nrow(new_stats)) { 
  new_stats[r, which(base::is.na(new_stats[r,]))] <- runif(length(which(base::is.na(new_stats[r,]))),
                                                           min = sort(gc_hplc$Values)[1],
                                                           max = sort(gc_hplc$Values)[2])
}

# stats_rq1 <- stats_rq1 %>%
#   rownames_to_column(., var = "product_cat")
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
for (product in unique(transpose_df$product_cat)) {
  grouped_product <- which(transpose_df$product_cat == product)
  qqlist[[i]] <- qqnorm(as.vector(t(transpose_df[grouped_product, 3:ncol(transpose_df)])), main = product, xlab = product, col = 'steelblue')
  qqline(as.vector(t(transpose_df[grouped_product, 3:ncol(transpose_df)])), col = 'red')
  i <- i + 1
}

### Conclusion = they are not normally distributed 



### Examining for equality of variance -----------------------------
# Levene’s test-non-normally distributed data (car::leveneTest()), significant if p-value < 0.05 
library(car)

results <- list()
for (col in 1:ncol(utils::combn(unique(transpose_df$product_cat), 2))) {
  # extract the combinations of plastic type pairs
  p_1 <- utils::combn(unique(transpose_df$product_cat), 2)[,col][1]
  p_2 <- utils::combn(unique(transpose_df$product_cat), 2)[,col][2]
  idx1 <- which(transpose_df$product_cat == p_1)
  idx2 <- which(transpose_df$product_cat == p_2)
  V1 <- as.vector(t(transpose_df[idx1, 3:ncol(transpose_df)]))
  V2 <- as.vector(t(transpose_df[idx2, 3:ncol(transpose_df)]))
  datat <- c(V1, V2)
  grouped <- as.factor(c(rep(p_1, times = length(V1)), rep(p_2, times = length(V2))))
  non_norm <- data.frame(datat, grouped)
  results[[paste0(p_1, p_2)]] <- car::leveneTest(datat ~ grouped, data = non_norm)
}
# -> p-value is 0.4096 which is greater the significance level of 0.05 therefore can conclude there is no significant difference between the variances


# Fligner-Killeen's test (fligner.test()), significant if p-value < 0.05 
fligner.test(datat ~ grouped, data = non_norm)
# p-value less than significance level therefore there is significant difference between variances.
# Test to determine the homogeneity of group variances. 



# Wilcoxon test to verify non normality of samples -----------
# Examine a difference in median value between distributions. 

# -> p-value smaller than significance level meaning two samples differ in distribution and median
summary_table_of_tests <- data.frame(comp=character(), pair=character(), wilcox_test=integer(), ks_test=integer())

for (comp in colnames(new_stats[,3:100])) {
  # looping through the combinations of product categories 
  for (col in 1:ncol(utils::combn(unique(new_stats$product_cat), 2))) {
    # extract the combinations of product category pairs
    p_1 <- utils::combn(unique(new_stats$product_cat), 2)[,col][1]
    p_2 <- utils::combn(unique(new_stats$product_cat), 2)[,col][2]
    # calculating the p-value between each product category pair 
    wilcox_test <- wilcox.test(new_stats[which(new_stats$product_cat == p_1), comp], 
                               new_stats[which(new_stats$product_cat == p_2), comp])$p.value
    ks_test <- ks.test(new_stats[which(new_stats$product_cat == p_1), comp],
                       new_stats[which(new_stats$product_cat == p_2), comp])$p.value
    # assigning row information
    summary_table_of_tests[nrow(summary_table_of_tests) + 1,] <- c(comp, paste0(p_1, "&", p_2), wilcox_test, ks_test)
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
  filter(., adjusted_pvalue_holm < 0.5) %>%
  arrange(adjusted_pvalue_holm)

summary_p0.05 <- summary_table_of_tests %>% filter(wilcox_test < 0.05 & ks_test < 0.05)

unique(summary_p0.05$comp)

top <- summary_p0.05 %>% arrange(wilcox_test, ks_test)

top_wil <- summary_p0.05 %>% arrange(wilcox_test)
top_ks <- summary_p0.05 %>% arrange(ks_test)

# 5 assumptions of multiple linear regression --------------------
transpose_df <- transpose_df %>% 
  mutate(plastic_type = factor(x = plastic_type, levels=unique(plastic_type))) %>%
  mutate(File = factor(x = File, levels=unique(File)))

must_convert<-sapply(transpose_df, is.factor)
new_pt <- sapply(transpose_df[,must_convert], unclass)
  # Using FactoMineR
out <- cbind(new_pt,var$contrib)
rownames(out) <- NULL

out <- as.data.frame(out) %>% 
  column_to_rownames(., var = "File") 

 # Linear relationship (predicting the plastic type)
linearity <- stats::lm(plastic_type ~ ., data = out)
summary(linearity)
#No linearity between predictor and response variables


 # Multicollinearity
library(car)
vif(linearity)
# all VIFs are between 1 and 5 which indicates moderate correlation between variables
# but not severe enough. 
  # check correlation between predictor variables
cor_matrix <- out[, 2:ncol(out)]
cor(cor_matrix)
  # create heat map to see correlation 
library(reshape2)
cor_mat <- round(cor(cor_matrix), 3)
melt_cor_mat <- melt(cor_mat)
#ggplot(data = melt_cor_mat, aes(x=Var1, y=Var2, fill=value)) +
#  geom_tile()
install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot::ggcorrplot(cor(cor_matrix))


 # Independence 
  # the Durbin-Watson Test
    # Used to detect autocorrelation in residuals. 
    # Ranges from 0 to 4, 2 is no autocorrelation, below 2 is positive autocorrelation 
    # and above 2 is negative correlation.
library(car)
durbinWatsonTest(linearity)
  # The statistic is 0.46 and p-value is 0 so 
  # reject null hypothesis and conclude residuals are auto correlated.
  # alternative hypothesis is true (true autocorrelation is greater than 0).

 # Homoscedasticity 
plot(fitted(linearity), resid(linearity), xlab='Fitted Values', ylab='Residuals')
abline(0,0)
  # residuals seem to follow an equal variance throughout the plot.

  # Formally tests for heteroscedasticity using Breusch-Pagan test 
install.packages("lmtest")
library(zoo)
library(lmtest)
bptest(linearity)
   # since the p-value is greater than 0.05 we don't reject the null hypothesis 
   # and conclude homoscedasticity is present and not violated.


 # Multivariate Normality 
   # Use a QQ plot to determine whether or not residuals of a model follow normal
   # distribution, if points form straight line then assumption is met. 

# qq norm takes numerical data like our data set out, but not in data frames. 

out %>% mutate_all(out)
