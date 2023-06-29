df_hc <- shared_comp_normalized %>%
  select(File, collapsed_compound, Percent_Area) %>%
  mutate(File = factor(File, levels = c(unique(File)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have multiple different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = File, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

for (r in 1:nrow(df_hc)) { 
  df_hc[r, which(base::is.na(df_hc[r,]))] <- runif(length(which(base::is.na(df_hc[r,]))),
                                                   min = sort(shared_comp_normalized$Percent_Area)[1],
                                                   max = sort(shared_comp_normalized$Percent_Area)[2])
}

plastic.data <- data.table::transpose(df_hc)
rownames(plastic.data) <- colnames(df_hc)
colnames(plastic.data) <- rownames(df_hc)
plastic.labs <- rownames(plastic.data)

# Distance matrix
data.dist_euclidean <- stats::dist(plastic.data)
data.dist_maximum <- stats::dist(plastic.data, method = "maximum")
data.dist_manhattan <- stats::dist(plastic.data, method = "manhattan")
data.dist_canberra <- stats::dist(plastic.data, method = "canberra")
data.dist_minkowski <- stats::dist(plastic.data, method = "minkowski")


plot(stats::hclust(data.dist_euclidean), xlab = "", sub = "", ylab = "",
     labels = plastic.labs, main = " Complete Linkage ")

plot(stats::hclust(data.dist_euclidean , method = "average"),
       labels = plastic.labs , main = " Average Linkage ",
       xlab = "", sub = "", ylab = "")

plot(stats::hclust(data.dist_euclidean , method = "single"),
       labels = plastic.labs , main = " Single Linkage ",
       xlab = "", sub = "", ylab = "")

# Extract important compounds that contribute to the split in dendrogram
