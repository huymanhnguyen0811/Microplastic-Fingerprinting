library(PCAtools)
library(stats)

# Bins values of collapsed_compounds into bins
?Hmisc::cut2

# PCA -----------------
# prep input
input_df <- function(data, pkg) {
  # create sample df
  if (pkg == "PCAtools") {
    df_X_rq1 <- data  %>%
      # dplyr::select(File, collapsed_compound, Values, product_cat) %>% # the options for select can be: File, plastic_type, product_cat, be careful since it will compressed obs. differently in dplyr::summarise() function subsequently
      # since we have multiple different values of the same compound in some samples, we summarize these values by taking the mean of them
      group_by(product_cat, collapsed_compound) %>%
      summarise(across(Values, mean)) %>%
      pivot_wider(names_from = product_cat, values_from = Values) %>%
      column_to_rownames(., var = "collapsed_compound")
  }
  
  for (r in 1:nrow(df_X_rq1)) {
    df_X_rq1[r, which(base::is.na(df_X_rq1[r,]))] <- runif(length(which(base::is.na(df_X_rq1[r,]))),
                                                           min = sort(data$Values)[1],
                                                           max = sort(data$Values)[2])
  }
  
  # table for information (rows are sample IDs, columns are sample information) -----------------------
  # metadata_X_rq1 <- data.frame(File = colnames(df_X_rq1))
  # #   data.frame(unique((data 
  # #                                      # filter(., !str_detect(File, "_USE"))
  # # )$File))
  # product_cat <- c()
  # for (row in 1:nrow(metadata_X_rq1)) {
  #   product_cat<- c(product_cat, unique(data[which(data$File == metadata_X_rq1[row, 'File']),]$product_cat))
  # }
  # metadata_X_rq1$product_cat <- product_cat
  # 
  # metadata_X_rq1 <- metadata_X_rq1 %>%
  # mutate(product_cat = ifelse(str_detect(File, "USE-01"), "Food contact materials",
  #                              ifelse(str_detect(File, "USE-02"), "Mixed_Plastic_Waste",
  #                                     ifelse(str_detect(File, "USE-03"), "Food contact materials",
  #                                            ifelse(str_detect(File, "USE-05"), "Cigarettes",
  #                                                   ifelse(str_detect(File, "USE-07"),"Food contact materials",
  #                                                          ifelse(str_detect(File, "USE-09"),"Food contact materials",
  #                                                                 ifelse(str_detect(File, "USE-11"),"Toys",
  #                                                                        ifelse(str_detect(File, "USE-13"),"Food contact materials",
  #                                                                               ifelse(str_detect(File, "USE-14"),"Food contact materials","Other")))))))))) %>%
  # mutate(product_cat = ifelse(str_detect(File, "Balloons"), "Toys",
  #                             ifelse(str_detect(File, "FPW_"), "Food contact materials",
  #                                    ifelse(str_detect(File, "Pbal_Sample"), "Toys",
  #                                           ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste",
  #                                                  ifelse(str_detect(File, "PBBC_"), "Food contact materials",
  #                                                         ifelse(str_detect(File, "Pbag_"),"Ziploc bags",
  #                                                                ifelse(str_detect(File, "PDS_Sample"),"Food contact materials",
  #                                                                       ifelse(str_detect(File, "Pcut_Sample"), "Food contact materials",
  #                                                                              ifelse(str_detect(File, "PC_Sample"), "Food contact materials",
  #                                                                                     ifelse(str_detect(File, "Cigs_"), "Cigarettes",
  #                                                                                            ifelse(str_detect(File, "Cmat"), "Construction materials",
  # ifelse(str_detect(File, "Mask_Sample"), "Clothes", "Misc"))))))))))))) %>%
  # column_to_rownames(., var = "File")
  
  # IF data$plastic_type then we look at mat = 7 rows of each plastic type; if data$product_cat then mat = ~6-7 rows of each product_cat
  # ATTENTION!!!  If group data like this, the observations are compressed and somewhat lost their variation because we calculated mean of all obs. for df_X_rq1 dataframe
  metadata_X_rq1 <- data.frame(product_cat = colnames(df_X_rq1)) # product_cat
  metadata_X_rq1$product_cat2 <- metadata_X_rq1$product_cat
  metadata_X_rq1 <- metadata_X_rq1 %>%
    column_to_rownames(., var = "product_cat")
  
  return(list(df_X_rq1 ,metadata_X_rq1))
}

df_pca <- input_df(gc_hplc, pkg = "PCAtools")

# PCA with PCAtools::pca ===========
# colnames(df_pca[[2]])[2] <- c("Plastic type")

# PCAtools::pca requires mat input (columns as sample name, rows as collapsed_compound)
p <- PCAtools::pca(mat = df_pca[[1]], 
                   metadata = df_pca[[2]], 
                   # center = FALSE,
                   scale = FALSE 
)

# Retrieve PC and add as new variables to data frame 
PCAtools_mergePC <- p$rotated

# PCA with stats::prcomp ===========
# stats::prcomp requires input df (columns as collapsed_compound, rows as sample name) -> change function df_pca pivot_wider(names_from=..)
df_pca <- input_df(merge_df, pkg = "stats")
prcomp_res <- stats::prcomp(df_pca[[1]], center = FALSE)
stats::biplot(x = prcomp_res)

# Retrieve PC and add as new variables to data frame 
e1071_merge_PC <- as.data.frame(prcomp_res$x)


# PCA further visualizations ----------------------------------------------------------------
# Scree plot
screeplot(p, components = getComponents(p),
          hline = 80, vline = 27, axisLabSize = 14, titleLabSize = 20,
          returnPlot = FALSE) 

# A bi-plot
PCAtools::biplot(p,
                 lab = rownames(p$metadata), # NULL 
                 colby = "product_cat2", 
                 hline = 0, vline = 0,
                 legendPosition = 'right', labSize = 5,
                 sizeLoadingsNames = 5,
                 showLoadings = TRUE,
                 showLoadingsNames = FALSE,
                 ntopLoadings = 10,
                 pointSize = 4, 
                 legendLabSize = 15,
                 legendTitleSize = 16,
                 legendIconSize = 6)


# Pairs plot
# pairsplot(p,
#           components = getComponents(p, c(1:5)),
#           triangle = FALSE,
#           trianglelabSize = 12,
#           hline = 0, vline = 0,
#           pointSize = 1.5,
#           gridlines.major = FALSE, gridlines.minor = FALSE,
#           colby = 'Plastic type',
#           title = 'Pairs plot',
#           axisLabSize = 14, plotaxes = TRUE,
#           margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


# explore further the collapsed_compounds that are driving these differences along each PC.
# plotloadings(p,
#              rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
#              caption = 'Top 10% variables',
#              labSize = 4)


# p.prcomp <- list(sdev = p$sdev,
#                  rotation = data.matrix(p$loadings),
#                  x = data.matrix(p$rotated),
#                  center = TRUE, scale = TRUE)
# predict(p.prcomp, newdata = newdata)[,1:5]
