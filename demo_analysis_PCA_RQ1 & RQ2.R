library(PCAtools)

# PCA -----------------
df_pca <- function(data) {
  # create sample df
  df_X_rq1 <- data %>%
    select(File, collapsed_compound, Percent_Area) %>%
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
  
  # table for information ((rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq1 <- data.frame(unique(data$File)) 
  colnames(metadata_X_rq1) <- c('File')
  metadata_X_rq1 <- metadata_X_rq1 %>%
    mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                                 ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                        ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                               ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                      ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                             ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
    column_to_rownames(., var = "File")
  
  return(list(df_X_rq1 ,metadata_X_rq1))
}

df_pca <- df_pca(shared_comp_sample)

# Conduct principal component analysis (PCA):
colnames(df_pca[[2]]) <- c("Plastic type")

p <- pca(mat = df_pca[[1]], metadata = df_pca[[2]])

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
loadingS_rq1_sorted <- p_rq1$loadings %>% arrange(PC1, PC2)
toploadings_rq1 <- rownames(loadingS_rq1_sorted[c(1:50, nrow(loadingS_rq1_sorted):(nrow(loadingS_rq1_sorted) - 50)),1:2])


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
plotloadings(p_rq1,
             rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
             caption = 'Top 10% variables',
             labSize = 4)


eigencorplot(p, metavars = c('fuel_type', 'gas_station'))