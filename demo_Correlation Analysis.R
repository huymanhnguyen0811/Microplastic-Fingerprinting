library(Hmisc)
library(corrplot)

# Correlation Analysis:
# Calculate the correlation coefficients between each collasped_compound

# Transform data frame for cor() function; file names as rows & collasped_compound as columns
my_data <- shared_comp_plastic_type %>%
  select(File, collapsed_compound, Percent_Area) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "File")
             
# Merge all correlated variables ===============================================

res2 <- Hmisc::rcorr(as.matrix(my_data), 
                     type = "spearman")

# NOTE:: Select all combination of compounds that have positive correlation > 0.7
# Create matrix with all combinations 
corr_mat <- flattenCorrMatrix(res2$r, res2$P) %>% filter(., (cor > 0.8 & p < 0.05))


