# Loading Packages --------------------------------------------------------
library(ggplot2)
library(readxl)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(lubridate)
library(dplyr)
library(data.table)
library(purrr)
library(stringr)
library(stringi)
library(plotly)
library(writexl)

# Functions -------------------------------------------------------------------------------------------------------
# Filtering matched compound names
filtering <- function(df, filter_list) {
  clean_data <-  copy(df)
  for (ele in filter_list) {
    clean_data <- clean_data %>%
      filter(!grepl(ele, Compound))
  }
  return(clean_data)
}

# Filtering limit of observations
limit_obser <- function(df_list, file_list, cap) {
  df_list_filter_area <- list()
  df_list_removed_area <- list()
  for (i in 1:length(df_list)) {
    df_list_filter_area[[i]] <- df_list[[i]] %>%
      filter(., Area > cap) %>%
      mutate(sample_name = file_list[[i]])
      
    
    df_list_removed_area[[i]] <- df_list[[i]] %>%
      filter(., Area <= cap) %>%
      mutate(sample_name = file_list[[i]])
  }
  return(list(df_list_filter_area, df_list_removed_area))
}

# Notin function
`%notin%` <- Negate(`%in%`)

# Grouping compounds based on RT1, RT2, and Ion1 - Version 1
grouping_comp_ver1 <- function(data, rt1thres, rt2thres, ion1thres, ion2thres) {
  
  # create empty list, each sub-list is a compound group with following criteria:
  # rt1thres: RT1 threshold window
  # rt2thres: RT2 threshold window
  # ion1thres: Ion1 threshold window
  # ion2thres: Ion2 threshold window
  # region_applied: list of x-y coordinate for different regions to applied different threshold window, x-axis is RT1, y-axis is RT2
  dat <- copy(data) %>% 
    arrange(RT1, RT2)
  
  # Initialize the compound column filled with NA values
  dat$collapsed_compound <- NA
  i <- 1
  
  for (row in 1:nrow(dat)) {
    # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    rt1 <- dat[row,]$RT1
    rt2 <- dat[row,]$RT2
    ion1 <- dat[row,]$Ion1
    ion2 <- dat[row,]$Ion2
    
    idx_thres <- which(dat$RT1 <= (rt1 + rt1thres) & dat$RT1 >= (rt1 - rt1thres) & 
                         dat$RT2 <= (rt2 + rt2thres) & dat$RT2 >= (rt2 - rt2thres) & 
                         dat$Ion1 <= (ion1 + ion1thres) & dat$Ion1 >= (ion1 - ion1thres) & 
                         dat$Ion2 <= (ion2 + ion2thres) & dat$Ion2 >= (ion2 - ion2thres) &
                         is.na(dat$collapsed_compound))
    
    if (identical(idx_thres, integer(0))) {
      next
    }
    else {
      dat[idx_thres, "collapsed_compound"] <- paste0("Compound_", i, ".")
      i <- i + 1
    }  
  }
  
  return(dat)
}

# Grouping compounds based on RT1, RT2, and Ion1 - Version 2
grouping_comp_ver2 <- function(data) {
  
  # create empty list, each sub-list is a compound group with following criteria:
  # rt1thres: RT1 threshold window
  # rt2thres: RT2 threshold window
  # ion1thres: Ion1 threshold window
  # ion2thres: Ion2 threshold window
  # region_applied: list of x-y coordinate for different regions to applied different threshold window, x-axis is RT1, y-axis is RT2
  data <- copy(all_data_pre_norm_filter_area)
  
  region_num <- as.numeric(base::readline("Please input the number of dividing region: "))
  rt1thres <- as.numeric(base::readline("Please input the RT1 window threshold for the region applied: "))
  rt2thres <- as.numeric(base::readline("Please input the RT2 window threshold for the region applied: "))
  ion1thres <- as.numeric(base::readline("Please input the Ion1 window threshold for the region applied: "))
  ion2thres <- as.numeric(base::readline("Please input the Ion2 window threshold for the region applied: "))
  
  # Initialize the compound column filled with NA values
  data$compound <- NA
  i <- 1
  
  for (reg_num in 1:region_num) {
    # User will input the coordinate of region that they want to applied a specific threshold of RT1
    region_x1 <- as.numeric(base::readline("Please input the bottom-left coordinate of the region applied: "))
    region_x2 <- as.numeric(base::readline("Please input the bottom-right coordinate of the region applied: "))
    region_y1 <- as.numeric(base::readline("Please input the top-left coordinate of the region applied: "))
    region_y2 <- as.numeric(base::readline("Please input the top-right coordinate of the region applied: "))
    
    region_applied <- c(region_x1, region_x2, region_y1, region_y2)
    
    # all_data_pre_norm_filter_area <- bind_rows(list_remaining_area) %>%
    #   arrange(RT1, RT2)
    # 
    # rt1thres <- 0.2
    # rt2thres <- 0.2
    # ion1thres <- 0.05
    # region_applied <- c(2, 72, 1, 6)
    
    # Filter data frame so that it only contain data in the region applied 
    idx_region <- which(data$RT1 >= region_applied[1] & data$RT1 <= region_applied[2] & 
                          data$RT2 >= region_applied[3] & data$RT2 <= region_applied[4])
    
    for (row in idx_region) {
      # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
      rt1 <- data[idx_region,][row,]$RT1
      rt2 <- data[idx_region,][row,]$RT2
      ion1 <- data[idx_region,][row,]$Ion1
      ion2 <- data[idx_region,][row,]$Ion2
      
      idx_thres <- which(data[idx_region,]$RT1 <= (rt1 + rt1thres) & data[idx_region,]$RT1 >= (rt1 - rt1thres) & 
                           data[idx_region,]$RT2 <= (rt2 + rt2thres) & data[idx_region,]$RT2 >= (rt2 - rt2thres) & 
                           data[idx_region,]$Ion1 <= (ion1 + ion1thres) & data[idx_region,]$Ion1 >= (ion1 - ion1thres) & 
                           data[idx_region,]$Ion2 <= (ion2 + ion2thres) & data[idx_region,]$Ion2 >= (ion2 - ion2thres) &
                           is.na(data[idx_region,]$compound))
      
      if (identical(idx_thres, integer(0))) {
        next
      }
      else {
        data[idx_region,][idx_thres, "compound"] <- paste0("Compound_", i, ".")
        i <- i + 1
      }  
    }
  }
  
  return(data)
}

# Filtering similar and unique compound
comp_filter_ver1 <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$collapsed_compound, fixed = TRUE))
    
    if (length(unique(data[idx,]$sample_name)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$sample_name)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  similar_compounds <- data[all_similar_compounds_idx,] 
  other_compounds <- data[all_other_compounds_idx,] 
  unique_compounds <- data[all_unique_compounds_idx,]
  
  return(list(similar_compounds, other_compounds, unique_compounds))
}


# Probabilistic Quotient Normalization
pqn <- function(X, n = "median", QC = NULL) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }
  
  # do the actual normalisation
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
  }
  
  return(X.norm)
}

# Relative log abundance plots
RlaPlots <- function(inputdata, type=c("ag", "wg"), cols=NULL,
                     cex.axis=0.8, las=2, ylim=c(-2, 2), oma=c(7, 4, 4, 2) + 0.1, ...) {
  type <- match.arg(type)
  groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
  unique.groups <- levels(groups)
  if (is.null(cols)) 
    cols <- ColList(length(unique.groups))
  box_cols <- c(rep(NA, length(rownames(inputdata))))
  for (ii in 1:length(inputdata[, 1])) 
    box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
  
  # Within groups
  if(type == "wg") {
    out_data<-data.frame()
    for (grp in unique.groups) {
      submat <- inputdata[which(inputdata[, 1] == grp), -1]
      med_vals <- apply(submat, 2, median)
      swept_mat <- sweep(submat, 2, med_vals, "-")
      out_data <- rbind(out_data, swept_mat)
    }
    # Across groups (i.e. type == "ag")
  } else  {
    med_vals <- apply(inputdata[, -1], 2, median)
    out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
  }
  
  boxplot(t(out_data),
          cex.axis=cex.axis,                 # font size
          las=las,                           # label orientation
          col=box_cols,                      # colours
          ylim=ylim,                         # y-axis range
          oma=oma,                           # outer margin size
          ...
  )
  
  abline(h=0)
}


# STEP 1.1: Data import --------------------------------------------
setwd("C:/Users/huyng/Desktop/Microplastic/Exported_CSV")

file_list <- list.files(pattern = '*.csv') %>%
  .[!str_detect(., "Blank")]

# Import samples to list
df_list_step1.1 <- purrr::map(file_list, read.csv, skip = 2)

# df_step1.1 <- dplyr::bind_rows(df_list_step1.1)

# summary(df_step1.1)

# Plotting data distribution pre-removal------------------
data_plot_pre_removal <- list() # NOTE: Data sets are all heavy left-skewed
for (i in 1:30) { # length(df_list_step1.1)
  filter_area <- df_list_step1.1[[i]]
  data_plot_pre_removal[[i]] <- ggplot(data = filter_area,
                                       aes(x = Area)) +
    geom_histogram(bins = 100) +
    ggtitle(file_list[[i]]) +
    # scale_x_continuous(limits = c(0, 2000000)) +
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "hidden", 
          axis.title = element_text(size = 5),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 12))
}

y <- grid::textGrob("Count", rot = 90, gp = gpar(fontsize = 20))
x <- grid::textGrob("Peak Area", gp = gpar(fontsize = 20))
grid.arrange(grobs = data_plot_pre_removal, ncol = 5, left = y, bottom = x)


# QUALITY CONTROL A OF STEP 1.2B: Plot Percentage coverage after removal of limit observation----------------------------
plot_a <- list()
j <- 1
for (i in 31:60) { # length(df_list_step1.1)
  coverage <- c()
  for (threshold in c(seq(from = 0, to = 200000, by = 50000))) {
    df_filter_area <- df_list_step1.1[[i]] %>%
      filter(., Area > threshold)
    coverage <- c(coverage, sum(df_filter_area$Area)*100/sum(df_list_step1.1[[i]]$Area))
  }
  df <- data.frame(thres = seq(from = 0, to = 200000, by = 50000), cover = coverage)
  plot_a[[j]] <- ggplot(data = df,
                        aes(x = thres, y = cover)) +
    geom_col() +
    theme(text = element_text(size = 15)) +
    geom_text(aes(label = round(cover, digits = 3)), color = "green", angle = 90, hjust = 1, size = 6) +
    scale_x_continuous(breaks = seq(from = 0, to = 200000, by = 50000),
                       # remove space between plotted data and xy-axes
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, by = 25), 
                       # remove space between plotted data and xy-axes
                       expand = c(0,0)) +
    ggtitle(file_list[[i]]) +
    labs(x = NULL, y = NULL) 
  j <- j + 1
}

y <- textGrob("Percentage coverage of remaining peaks after removal", rot = 90, gp = gpar(fontsize = 20))
x <- textGrob("Threshold of removal for limit observations", gp = gpar(fontsize = 20))

grid.arrange(grobs = plot_a, ncol = 5, 
             left = y,
             bottom = x)

# QUALITY CONTROL B OF STEP 1.2B: Plot number of peak remains after removal of limit observation----------------------------
plot_b <- list()
j <- 1
for (i in 31:60) { # length(df_list_step1.1)
  peak_remain <- c()
  for (threshold in c(seq(from = 0, to = 200000, by = 50000))) {
    df_filter_area <- df_list_step1.1[[i]] %>%
      filter(., Area > threshold)
    peak_remain <- c(peak_remain, dim(df_filter_area)[1])
  }
  df <- data.frame(thres = seq(from = 0, to = 200000, by = 50000), remain = peak_remain)
  plot_b[[j]] <- ggplot(data = df,
                        aes(x = thres, y = remain)) +
    geom_col() +
    geom_text(aes(label = remain), color = "green", vjust = 1.2, size = 5) +
    scale_x_continuous(breaks = seq(from = 0, to = 200000, by = 50000), 
                       # remove space between plotted data and xy-axes
                       expand = c(0,0)) +
    ggtitle(file_list[[i]]) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)) +
    labs(x = NULL, y = NULL) +
    theme_classic()
  j <- j + 1
}

y <- textGrob("Number of peak remains after removal of limit observation", rot = 90, gp = gpar(fontsize = 15))
x <- textGrob("Threshold of removal for limit observations", gp = gpar(fontsize = 15))

grid.arrange(grobs = plot_b, ncol = 5, 
             left = y,
             bottom = x)


# STEP 1.2B Filtering out limit of observations---------------------------- 
# Area  = 100000 seems to be the right cutoff after inspection with quality control A and B
list_remaining_area <- limit_obser(df_list_step1.1, file_list, cap = 100000)[[1]]
list_removed_area <- limit_obser(df_list_step1.1, file_list, cap = 100000)[[2]]

# STEP 1.3: Grouping compounds based on RT1, RT2, Ion1 -----------------------------------------------------------------------
# STEP 1.3A: Generate 1 grand data frame of all 31 IL samples

df_step1.3 <- bind_rows(list_remaining_area)


# STEP 1.3B: Grouping compounds based on RT1, RT2, Ion1 threshold 
# rt10.2 <- grouping_comp_ver1(df_step1.3,
#                              rt1thres = 0.2,
#                              rt2thres = 0.15,
#                              ion1thres = 0.05, # Ion 1 and 2 indicates molecular structure (2 most prevalent mass-to-charge)
#                              ion2thres = 0.05)

rt10.1 <- grouping_comp_ver1(df_step1.3,
                             rt1thres = 0.1,
                             rt2thres = 0.15,
                             ion1thres = 0.05, # Ion 1 and 2 indicates molecular structure (2 most prevalent mass-to-charge)
                             ion2thres = 0.05)


# rt10.3 <- grouping_comp_ver1(df_step1.3,
#                              rt1thres = 0.3,
#                              rt2thres = 0.15,
#                              ion1thres = 0.05, # Ion 1 and 2 indicates molecular structure (2 most prevalent mass-to-charge)
#                              ion2thres = 0.05)

# STEP 2: Identify shared and unique compound groups across samples ------------------------------------------------
# RT1 0.2
# filter_rt10.2 <- comp_filter_ver1(rt10.2, 
#                                   length(ILR_file_list))
# 
# shared_comp_rt10.2 <- bind_rows(filter_rt10.2[[1]], filter_rt10.2[[2]])

# RT1 0.1
filter_rt10.1 <- comp_filter_ver1(rt10.1, 
                                  length(ILR_file_list))

shared_comp_rt10.1 <- bind_rows(filter_rt10.1[[1]], filter_rt10.1[[2]])

# RT1 0.3
# filter_rt10.3 <- comp_filter_ver1(rt10.3, 
#                                   length(ILR_file_list))
# 
# shared_comp_rt10.3 <- bind_rows(filter_rt10.3[[1]], filter_rt10.3[[2]])

# STEP 3: Data Normalization ================================================================================
add_data_normalization <- function(data) {
  temp_list <- list()
  i <- 1
  for (sample in unique(data$sample_name)) {
    df <- data[which(data$sample_name == sample),] %>%
      # Log-based normalization
      mutate(Log_Area = log10(Area)) %>%
      mutate(Log_Height = log10(Height)) %>%
      # TSN - Percent-based normalization
      mutate(Percent_Area = Area/sum(.$Area)) %>%
      mutate(Percent_Height = Height/sum(.$Height))
    # Box-Cox normalization
    df$boxcox_area <- bestNormalize::boxcox(df$Area)[[1]]
    df$boxcox_height <- bestNormalize::boxcox(df$Height)[[1]]
    temp_list[[i]] <- df
    i <- i + 1
  }
  # Then combine data again to 1 grand data frame
  newdata <- dplyr::bind_rows(temp_list)
  return(newdata)
}

shared_comp_normalized_rt10.1 <- add_data_normalization(shared_comp_rt10.1)
# shared_comp_normalized_rt10.2 <- add_data_normalization(shared_comp_rt10.2)
# shared_comp_normalized_rt10.3 <- add_data_normalization(shared_comp_rt10.3)

