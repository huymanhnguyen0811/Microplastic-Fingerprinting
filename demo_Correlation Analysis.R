library(Hmisc)
library(corrplot)

# Correlation Analysis:
# Calculate the correlation coefficients between each collasped_compound

# Transform data frame for cor() function; file names as rows & collasped_compound as columns
my_data <- shared_comp_plastic_type %>%
  dplyr::select(File, collapsed_compound, Percent_Area) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  tibble::column_to_rownames(., var = "File")
             
# Fill NA
for (c in 1:ncol(my_data)) { 
  my_data[which(base::is.na(my_data[,c])),c] <- runif(length(which(base::is.na(my_data[,c]))),
                                                      min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                      max = sort(shared_comp_plastic_type$Percent_Area)[2])
}

# correlated variables ===============================================

res2 <- Hmisc::rcorr(as.matrix(my_data), 
                     type = "spearman")

res3 <- Hmisc::rcorr(as.matrix(my_data), 
                     type = "spearman")

# NOTE:: Select all combination of compounds that have positive correlation --------
 
corr_mat <- flattenCorrMatrix(res2$r, res2$P) %>% filter(., (cor > 0.8 & p < 0.05))

corr_mat_filled <- flattenCorrMatrix(res3$r, res3$P) %>% filter(., (cor > 0.8 & p < 0.05))


corrplot(cor(my_data), order = "hclust", tl.col = "black", tl.srt = 45)


View(as.matrix(stats::dist(as.matrix(my_data), method = "minkowski", upper = TRUE)))
