library(PCAtools)


# PCA -----------------
# prep input
df_pca <- function(data) {
  # create sample df
  df_X_rq1 <- data %>%
    dplyr::select(File, collapsed_compound, Percent_Area) %>%
    filter(., !str_detect(File, "PC_Sample")) %>% # Exclude Plastic_cups to examine further clustering of others types
    filter(., !str_detect(File, "Balloons_Sample")) %>%
    mutate(File = factor(File, levels = unique(File))) %>%
    # since we have multiple different values of the same compound in some samples, we summarize these values by taking the mean of them
    group_by(File, collapsed_compound) %>%
    summarise(across(Percent_Area, mean)) %>%
    pivot_wider(names_from = File, values_from = Percent_Area) %>%
    column_to_rownames(., var = "collapsed_compound")
  

  for (r in 1:nrow(df_X_rq1)) { 
    df_X_rq1[r, which(base::is.na(df_X_rq1[r,]))] <- runif(length(which(base::is.na(df_X_rq1[r,]))),
                                                           min = sort(data$Percent_Area)[1],
                                                           max = sort(data$Percent_Area)[2])
  }
  
  # table for information (rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq1 <- data.frame(unique(data$File)) 
  colnames(metadata_X_rq1) <- c('File')
  metadata_X_rq1 <- metadata_X_rq1 %>% 
    filter(., !str_detect(File, "PC_Sample")) %>% # Exclude Plastic_cups to examine further clustering of others types
    filter(., !str_detect(File, "Balloons_Sample")) %>%
    mutate(plastic_type = ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                 # ifelse(str_detect(File, "Balloons"), "Balloons",
                                        ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                               ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                      # ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                             ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))) %>%
    column_to_rownames(., var = "File")
  
  return(list(df_X_rq1 ,metadata_X_rq1))
}

df_pca <- df_pca(shared_comp_plastic_type)

# PCA with PCAtools::pca ===========
colnames(df_pca[[2]]) <- c("Plastic type")

# PCAtools::pca requires mat input (columns as sample name, rows as collapsed_compound)
p <- PCAtools::pca(mat = df_pca[[1]], metadata = df_pca[[2]])

# Retrieve PC and add as new variables to data frame 
PCAtools_mergePC <- p$rotated

# PCA with stats::prcomp ===========
# stats::prcomp requires input x as df (columns as collapsed_compound, rows as sample name) -> change function df_pca pivot_wider(names_from=..)

prcomp_res <- stats::prcomp(df_pca[[1]])

# stats::biplot(prcomp_res)

# Retrieve PC and add as new variables to data frame 
e1071_merge_PC <- as.data.frame(prcomp_res$x)


# PCA further visualizations ----------------------------------------------------------------
# Scree plot
screeplot(p, components = getComponents(p, 1:30),
          hline = 80, vline = 27, axisLabSize = 14, titleLabSize = 20,
          returnPlot = FALSE) +
  geom_label(aes(20, 80, label = '80% explained variation', vjust = -1, size = 8))

# A bi-plot
PCAtools::biplot(p,
                 lab = NULL, 
                 colby = "Plastic type",
                 hline = 0, vline = 0,
                 legendPosition = 'right', labSize = 5,
                 sizeLoadingsNames = 5,
                 showLoadings = TRUE,
                 # showLoadingsNames = FALSE,
                 ntopLoadings = 10,
                 pointSize = 4, 
                 legendLabSize = 15,
                 legendTitleSize = 16,
                 legendIconSize = 6)

# Retrieve compound name of top 100 loading
loadingS_sorted <- p$loadings %>% arrange(PC1, PC2)
toploadings <- rownames(loadingS_sorted[c(1:50, nrow(loadingS_sorted):(nrow(loadingS_sorted) - 50)),1:2])


# Pairs plot
pairsplot(p,
          components = getComponents(p, c(1:5)),
          triangle = FALSE,
          trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 1.5,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Plastic type',
          title = 'Pairs plot',
          axisLabSize = 14, plotaxes = TRUE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


# explore further the collapsed_compounds that are driving these differences along each PC.
plotloadings(p,
             rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
             caption = 'Top 10% variables',
             labSize = 4)


eigencorplot(p, metavars = c('fuel_type', 'gas_station'))

p.prcomp <- list(sdev = p$sdev,
                 rotation = data.matrix(p$loadings),
                 x = data.matrix(p$rotated),
                 center = TRUE, scale = TRUE)
predict(p.prcomp, newdata = newdata)[,1:5]
