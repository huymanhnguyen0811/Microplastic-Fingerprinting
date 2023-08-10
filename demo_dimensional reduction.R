library(stats)
library(vegan)

input_df <- function(data) {
  # create sample df
  df_X_rq1 <- data %>%
    dplyr::select(File, collapsed_compound, Percent_Area) %>%
    # filter(., !str_detect(File, "PC_Sample")) %>% # Exclude Plastic_cups to examine further clustering of others types
    # filter(., !str_detect(File, "Balloons_Sample")) %>%
    mutate(File = factor(File, levels = unique(File))) %>%
    # since we have multiple different values of the same compound in some samples, we summarize these values by taking the mean of them
    group_by(File, collapsed_compound) %>%
    summarise(across(Percent_Area, mean)) %>%
    pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
    column_to_rownames(., var = "File")
  
  
  for (r in 1:nrow(df_X_rq1)) {
    df_X_rq1[r, which(base::is.na(df_X_rq1[r,]))] <- runif(length(which(base::is.na(df_X_rq1[r,]))),
                                                           min = sort(data$Percent_Area)[1],
                                                           max = sort(data$Percent_Area)[2])
  }
  
  # table for information (rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq1 <- data.frame(unique(data$File)) 
  colnames(metadata_X_rq1) <- c('File')
  metadata_X_rq1 <- metadata_X_rq1 %>% 
    # filter(., !str_detect(File, "PC_Sample")) %>% # Exclude Plastic_cups to examine further clustering of others types
    # filter(., !str_detect(File, "Balloons_Sample")) %>%
    mutate(plastic_type = ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                 # ifelse(str_detect(File, "Balloons"), "Balloons",
                                 ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                        ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                               # ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                               ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))) %>%
    column_to_rownames(., var = "File")
  
  return(list(df_X_rq1 ,metadata_X_rq1))
}

inputdf <- input_df(shared_comp_plastic_type)[[1]]
# Multidimensional scaling =====================================================


#calculate distance matrix
d <- stats::dist(inputdf)

# Testing multiple  maximum dimension 'k' for multidimensional scaling
for (testk in seq(from = 3, to = 75, length.out = 5)) {
  fit <- stats::cmdscale(d,
                         eig = TRUE,
                         k = testk)
  
  #extract first and second coordinates of multidimensional scaling
  x <- fit$points[,1]
  y <- fit$points[,2]
  
  #create scatter plot
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="Multidimensional Scaling Results", type="n")
  
  #add row names of data frame as labels
  text(x, y, labels=row.names(inputdf))
}

# Isometric Feature mapping ordination =========================================
# dis <- vegan::vegdist(inputdf)
# ord <- isomap(dis, ndim = 5, k = 4
#               # fragmentedOK = TRUE
#               )
# 
# # extract the ordination scores.
# ord_score <- vegan::scores(ord)
# 
# pl <- plot(ord, main="isomap k=3")

library(RDRToolbox)
isomap.res <- RDRToolbox::Isomap(data = as.matrix(inputdf),
                                 dims = 1:30,
                                 k = 5,
                                 plotResiduals = TRUE)

# extract the ordination scores.
ord_score <- isomap.res$dim30
rownames(ord_score) <- rownames(inputdf)


# t-Distributed Stochastic Neighbor Embedding (t-SNE) ==========================
library(Rtsne)

rtsne_out <- Rtsne(as.matrix(inputdf), dims = 2,
                   pca = FALSE, verbose = TRUE,
                   theta = 0.0, max_iter = 2000, 
                   Y_init = NULL, perplexity = 20)

tsne_coord <- as.data.frame(rtsne_out$Y) %>%
  cbind(File = rownames(inputdf)) %>%
  left_join(inputdf %>%
              rownames_to_column(., var = "File"), 
            by = 'File') %>% 
  mutate(plastic_type = ifelse(str_detect(File, "FPW_"), 0,
                               ifelse(str_detect(File, "Balloons"), 1,
                                      ifelse(str_detect(File, "MPW_"), 2, 
                                             ifelse(str_detect(File, "PBBC_"), 3,
                                                    ifelse(str_detect(File, "PC_Sample"),4,
                                                           ifelse(str_detect(File, "PDS_Sample"),5, 6)))))))
          
library(ggrepel)

gg <- ggplot() +
  # labs(title = "") +
  scale_fill_gradient(low = 'darkblue', high = 'red') +
  coord_fixed(ratio = 1) +
  theme_void() +
  stat_summary_hex(data = tsne_coord, aes(x=V1, y=V2, z = plastic_type), bins=10, fun = mean, alpha = 0.9) +
  geom_point(data = filter(tsne_coord, plastic_type == 0), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'black') +
  geom_point(data = filter(tsne_coord, plastic_type == 1), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'white') +
  geom_point(data = filter(tsne_coord, plastic_type == 2), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'green') +
  geom_point(data = filter(tsne_coord, plastic_type == 3), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'blue') +
  geom_point(data = filter(tsne_coord, plastic_type == 4), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'yellow') +
  geom_point(data = filter(tsne_coord, plastic_type == 5), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'orange') +
  geom_point(data = filter(tsne_coord, plastic_type == 6), aes(x = V1, y = V2), alpha = 0.7, size = 1, col = 'purple') +
  theme(plot.title = element_text(hjust = 0.5, family = 'Calibri'),
        legend.title.align=0.5) +
  geom_text_repel(data = tsne_coord, aes(x = V1, y = V2, label = plastic_type), size = 5, col = 'black')

gg


# K-means clustering for gene expression =======================================
# Goal is to cluster X collapsed_compounds into n (n << X) clusters
# perform PCA
pca = prcomp(t(test), center=TRUE, scale=TRUE)
rotation = data.frame(pca$x)
plot(rotation[1:3], pch=16, cex=0.6, cex.main=0.9)

# perform Kmeans in PCA space
wss = (nrow(rotation)-1)*sum(apply(rotation,2,var))
for (i in 2:9) wss[i] <- sum(kmeans(rotation, centers=i)$withinss)
plot(1:9, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
km = kmeans(rotation, 2)
cluster_assignment = as.factor(km$cluster)

# plot k-means cluster assignment in PCA space
library(ggplot2)
ggplot(rotation, aes(rotation$PC1, rotation$PC2, color=cluster_assignment, label=rownames(rotation))) + geom_point() + geom_text()

threshold <- 1
nms <- as.character()
for(n in 1:n_clusters) {
  
  for(i in 1:ncol(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),])){
    if(mean(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),i]) > threshold){
      nms <- c(nms, colnames(t(test)[row.names(t(test)) %in% names(cluster_assignment[cluster_assignment == n]),])[i])
    }
  }
  if(n == 1) result <- data.frame(GENE = nms, CLUSTER = rep(n,length(nms))); rm(nms); nms <- as.character()
  if(n >  1) result <- rbind(result, data.frame(GENE = nms, CLUSTER = rep(n,length(nms))))
}

result


# Hierarchichal Clustering of collapsed_compounds ==============================
# Do not scale the colaapsed_compounds bofre compute dissimilarity, because we already normalized data 
# and need to preserve chemical patterns as much as possible . Also because all the unit are already on the same scale
# https://stats.stackexchange.com/questions/110428/clustering-a-noisy-data-or-with-outliers

library(rsample)

clus.cv.res <- rsample::clustering_cv()