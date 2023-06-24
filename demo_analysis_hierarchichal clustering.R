plastic.data <- data.table::transpose(df_pca[[1]])
rownames(plastic.data) <- colnames(df_pca[[1]])
colnames(plastic.data) <- rownames(df_pca[[1]])
plastic.labs <- rownames(plastic.data)


data.dist_euclidean <- stats::dist(plastic.data)
data.dist_maximum <- stats::dist(plastic.data, method = "maximum")
data.dist_manhattan <- stats::dist(plastic.data, method = "manhattan")
data.dist_canberra <- stats::dist(plastic.data, method = "canberra")
data.dist_minkowski <- stats::dist(plastic.data, method = "minkowski")


plot(stats::hclust(data.dist_minkowski), xlab = "", sub = "", ylab = "",
     labels = plastic.labs, main = " Complete Linkage ")

plot ( hclust(data.dist_minkowski , method = "average"),
       labels = plastic.labs , main = " Average Linkage ",
       xlab = "", sub = "", ylab = "")

plot (hclust(data.dist_minkowski , method = "single"),
       labels = plastic.labs , main = " Single Linkage ",
       xlab = "", sub = "", ylab = "")

# each row of .data is a sample name
# 