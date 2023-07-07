library(Hmisc)
library(corrplot)

# Correlation Analysis:
# Calculate the correlation coefficients between each collasped_compound and the plastic_type

# Examine correlation between plastic types -----------------
my_data <- shared_comp_normalized %>%
  select(plastic_type, collapsed_compound, Percent_Area) %>%
  mutate(plastic_type = factor(plastic_type, levels = c(unique(plastic_type)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(plastic_type, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = plastic_type, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

# With Pearson correlation coefficient
res <- stats::cor(my_data, 
                  method = "pearson", 
                  use = "pairwise.complete.obs") # FYI: ?stats::cor


View(res)

# Correlation < 0.5 is considered weak/no linear relationship

res2 <- Hmisc::rcorr(as.matrix(my_data), 
                     type = "pearson")
flattenCorrMatrix(res2$r, res2$P)


corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01)


# Dimensionality Reduction Technique