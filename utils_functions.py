#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def filtering(data, filter_list):
    clean_data = data.copy()
    for ele in filter_list:
        clean_data = clean_data[-clean_data['Compound'].str.contains(ele, flags = re.IGNORECASE)]
    return clean_data

# # Grouping compounds based on RT1, RT2, and Ion1
# grouping_comp <- function(data) {

#   # create empty list, each sub-list is a compound group with following criteria:
#   # RT1 threshold +- 0.4
#   # RT2 threshold +- 0.1
#   # Ion1 threshold +-0.1
  
#   # Initialize the compound_group column filled with NA values
#   data$compound_group <- NA
#   i <- 1
#   for (row in 1:nrow(data)) {
#     # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
#     rt1 <- data[row,]$RT1
#     rt2 <- data[row,]$RT2
#     ion1 <- data[row,]$`Ion 1`
#     idx <- which(data$RT1 < (rt1 + 0.4) & data$RT1 > (rt1 - 0.4) & 
#                  data$RT2 < (rt2 + 0.125) & data$RT2 > (rt2 - 0.125) & 
#                  data$`Ion 1` < (ion1 + 0.1) & data$`Ion 1` > (ion1 - 0.1) & 
#                  is.na(data$compound_group))
#     if (identical(idx, integer(0))) {
#       next
#     }
#     else {
#       data[idx, "compound_group"] <- paste0("Compound_", i, ".")
#       i <- i + 1
#     }
#     rm(rt1)
#     rm(rt2)
#     rm(ion1)
#   }
#   return(data)
# }

# # Filtering similar and unique compound
# comp_filter <- function(data, file_list) {
#   all_similar_compounds_idx <- c()
#   all_other_compounds_idx <- c()
#   all_unique_compounds_idx <- c()
  
#   for (comp_grp in unique(data$compound_group)) {
#     # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    
#     idx <- which(grepl(paste0("^", comp_grp, "$"), data$compound_group))
    
#     if (length(unique(data[idx,]$sample_name)) > (length(file_list) - 1)) {
#       all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
#     }
#     else if (length(unique(data[idx,]$sample_name)) < 2) {
#       all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
#     }
#     else {
#       all_other_compounds_idx <- c(all_other_compounds_idx, idx)
#     }
#   }
#   return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
# }

# # Probabilistic Quotient Normalization
# pqn <- function(X, n = "median", QC = NULL) {
#   X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
#   colnames(X.norm) <- colnames(X)
#   rownames(X.norm) <- rownames(X)
  
#   if (!is.null(QC)) {
#     # if QC vector exists, use this as reference spectrum
#     if (length(QC) == 1) {
#       # only 1 reference sample given
#       mX <- as.numeric(X[QC, ])
#     } else {
#       if (n == "mean") {
#         mX <- as.numeric(colMeans(X[QC, ]))
#       }
#       if (n == "median") {
#         mX <- as.numeric(apply(X[QC, ], 2, median))
#       }
#     }
#   } else {
#     # otherwise use the mean or median of all samples as reference sample
#     if (n == "mean") {
#       mX <- as.numeric(colMeans(X))
#     }
#     if (n == "median") {
#       mX <- as.numeric(apply(X, 2, median))
#     }
#   }
  
#   # do the actual normalisation
#   for (a in 1:nrow(X)) {
#     X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
#   }
  
#   return(X.norm)
# }

# # RLA plots
# RlaPlots <- function(inputdata, type=c("ag", "wg"), cols=NULL,
#                      cex.axis=0.8, las=2, ylim=c(-2, 2), oma=c(7, 4, 4, 2) + 0.1, ...) {
#   type <- match.arg(type)
#   groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
#   unique.groups <- levels(groups)
#   if (is.null(cols)) 
#     cols <- ColList(length(unique.groups))
#   box_cols <- c(rep(NA, length(rownames(inputdata))))
#   for (ii in 1:length(inputdata[, 1])) 
#     box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
  
#   # Within groups
#   if(type == "wg") {
#     out_data<-data.frame()
#     for (grp in unique.groups) {
#       submat <- inputdata[which(inputdata[, 1] == grp), -1]
#       med_vals <- apply(submat, 2, median)
#       swept_mat <- sweep(submat, 2, med_vals, "-")
#       out_data <- rbind(out_data, swept_mat)
#     }
#     # Across groups (i.e. type == "ag")
#   } else  {
#     med_vals <- apply(inputdata[, -1], 2, median)
#     out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
#   }
  
#   boxplot(t(out_data),
#           cex.axis=cex.axis,                 # font size
#           las=las,                           # label orientation
#           col=box_cols,                      # colours
#           ylim=ylim,                         # y-axis range
#           oma=oma,                           # outer margin size
#           ...
#   )
  
#   abline(h=0)
# }

