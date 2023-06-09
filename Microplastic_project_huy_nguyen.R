# Loading Packages --------------------------------------------------------
library(ggplot2)
library(readxl)
library(ggpubr)
library(gtable)
library(gridExtra)
library(viridis)
library(wesanderson)
library(tidyverse)
library(lubridate)
library(dplyr)
library(data.table)
library(corrplot)
library(purrr)
library(stringr)
library(stringi)
library(grid)
library(plotly)
library(umap)
library(mdatools)
library(writexl)
library(ggsignif)


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
      filter(., Area > cap)
      # mutate(sample_name = file_list[[i]])
      
    df_list_removed_area[[i]] <- df_list[[i]] %>%
      filter(., Area <= cap)
      # mutate(sample_name = file_list[[i]])
  }
  return(list(df_list_filter_area, df_list_removed_area))
}
  
# Notin function
`%notin%` <- Negate(`%in%`)

# Grouping compounds based on RT1, RT2, and Ion1 - Version 1
grouping_comp_ver1 <- function(data, rtthres, mzthres) {
  
  # create empty list, each sub-list is a compound group with following criteria:
  # rtthres: RT threshold window
  # mzthres: mz threshold window
  dat <- copy(data)
  
  # Initialize the compound column filled with NA values
  dat$collapsed_compound <- NA
  i <- 1
  
  for (row in 1:nrow(dat)) {
    # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    rt <- dat[row,]$RT
    mz <- dat[row,]$m.z
    
    idx <- which(dat$RT <= (rt + rtthres) & dat$RT >= (rt - rtthres) &
                   dat$m.z <= (mz + mzthres) & dat$m.z >= (mz - mzthres) &
                   is.na(dat$collapsed_compound))
    
    if (identical(idx, integer(0))) {
      next
    }
    else {
      dat[idx, "collapsed_compound"] <- paste0("Compound_", i, ".")
      i <- i + 1
    }  
  }
  
  return(dat)
}

# Filtering similar and unique compound 
# compound appear in at least 2 samples

comp_filter_ver1 <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$collapsed_compound, fixed = TRUE))
    
    if (length(unique(data[idx,]$File)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$File)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
}

# compound appear in at least 2 plastic types

comp_filter_ver2 <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$collapsed_compound, fixed = TRUE))
    
    if (length(unique(data[idx,]$plastic_type)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$plastic_type)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
}

# TSN - Percent-based normalization 
data_normalization <- function(data) {
  temp_list <- list()
  i <- 1
  # Normalize Peak Area for each sample 
  for (sample in unique(data$File)) {
    df <- data[which(data$File == sample),] %>%
      mutate(Percent_Area = Area/sum(.$Area)) %>%
      mutate(Percent_Height = Height/sum(.$Height))
    temp_list[[i]] <- df
    i <- i + 1
  }
  # Then combine data again to 1 grand data frame
  newdata <- dplyr::bind_rows(temp_list)
  return(newdata)
}


# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# STEP 1.1: Data import --------------------------------------------
setwd("C:/Users/huyng/OneDrive - Toronto Metropolitan University/Microplastic/Microplastic-Fingerprinting/CSV_Export_2023-05-23_Cleaned")

file_list <- list.files(pattern = '*.csv') %>%
  .[!str_detect(., "Blank")]

# Blank samples 
blank_list <- list.files(pattern = '*.csv') %>%
  .[str_detect(., "Blank")]

# Import samples to list
df_list_step1.1 <- purrr::map(file_list, read.csv)

df_list_blank <- purrr::map(blank_list, read.csv)

# df_step1.1 <- dplyr::bind_rows(df_list_step1.1)
df_blank <- dplyr::bind_rows(df_list_blank)
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


# QUALITY CONTROL C (CONFIRMATION) OF STEP 1.2B: -----------------------
## Plotting post-removal data distribution
data_plot_post_removal <- list() 
j <- 1
for (i in 31:60) { # length(list_remaining_area)
  filter_area <- list_remaining_area[[i]]
  data_plot_post_removal[[j]] <- ggplot(data = filter_area,
                                    aes(x = Area)) +
    geom_histogram(bins = 100) +
    ggtitle(file_list[[i]]) +
    # scale_x_continuous(limits = c(0, 2000000)) +
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "hidden", axis.title = element_text(size = 5))
  j <- j + 1
}

y <- grid::textGrob("Count", rot = 90, gp = gpar(fontsize = 15))
x <- grid::textGrob("Peak Area", gp = gpar(fontsize = 15))
grid.arrange(grobs = data_plot_post_removal, ncol = 5, left = y, bottom = x)


# STEP 1.3: Grouping compounds based on Retention time and molecular ions  -----------------------------------------------------------------------
# STEP 1.3A: Generate 1 grand data frame of all 31 IL samples

df_step1.3 <- bind_rows(list_remaining_area) %>%
  select(-c("Start", "End", "Width", "Base.Peak")) %>%
  mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                            ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                   ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                          ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                 ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                        ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other")))))))

df_blank <- df_blank %>%
  select(-c("Start", "End", "Width", "Base.Peak")) %>%
  mutate(plastic_type = "Blanks")

combined_df <- rbind(df_step1.3, df_blank) %>% arrange(RT)


# STEP 1.3B: Collapsing compounds based on RT1, RT2, Ion1 threshold

combined_df_grouped <- grouping_comp_ver1(combined_df,
                                          rtthres = 0.05,
                                          mzthres = 0.05)



# Plotting data distribution pre-removal -----------------------------------
data_plot_pre_removal <- list() 
i <- 1
for (sample in unique(other_compounds_filter_area_samples$sample_name)) {
  data_plot_pre_removal[[i]] <- ggplot(data = other_compounds_filter_area_samples %>% 
                                         filter(., sample_name %in% sample),
                                       aes(x = Area)) +
    geom_histogram(bins = 50) +
    ggtitle(sample) +
    # scale_x_continuous() + # limits = c(0, 2000000)
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "hidden", axis.title = element_text(size = 5))
  i <- i + 1
}

y <- grid::textGrob("Count", rot = 90, gp = gpar(fontsize = 15))
x <- grid::textGrob("Peak Area", gp = gpar(fontsize = 15))
grid.arrange(grobs = data_plot_pre_removal, ncol = 5, left = y, bottom = x)


# STEP 2: Normalizing data accordingly to different data frames of interest --------------------------------------------------

comp_normalized <- data_normalization(combined_df_grouped)

# Step 3: Readjust compound RA (sample) by average blank RA ======================
# Create list to store temp dfs
temp_list <- list()
i <- 1
# Iterate through each collapsed_compound
for (comp in unique(comp_normalized$collapsed_compound)) {
  temp <- comp_normalized[which(comp_normalized$collapsed_compound == comp),]
  # if compound does not exist in blanks then skip the compounds
  if (identical(which(temp$plastic_type == "Blanks"), integer(0))) {
    temp_list[[i]] <- temp
    i <- i + 1
    next
  }
  else {
    # Calculate avg_blank for that compound across all blanks
    avg_blank <- mean(temp[which(temp$plastic_type == "Blanks"),]$Percent_Area)
    temp <- temp[which(temp$plastic_type != "Blanks"),]
    # iterate through each sample
    for (sample in unique(temp$File)) {
      # Adjust RA for each compound of each sample = RA (sample) - avg_blank
      temp[which(temp$File == sample),]$Percent_Area <- temp[which(temp$File == sample),]$Percent_Area - avg_blank
    } 
  }
  # Append current temp df to temp_list
  temp_list[[i]] <- temp
  i <- i + 1
}

adjusted_df <- bind_rows(temp_list)

# Step 4: Replace negative adjusted RA values with LOD =============================
adjusted_df$Percent_Area[adjusted_df$Percent_Area < 0] <- runif(length(adjusted_df$Percent_Area[adjusted_df$Percent_Area < 0]),
                                                                min = sort(adjusted_df$Percent_Area[adjusted_df$Percent_Area > 0])[1],
                                                                max = sort(adjusted_df$Percent_Area[adjusted_df$Percent_Area > 0])[2])

# STEP 5: Identify shared and unique compound groups across samples ------------------------------------------------
# at least in 2 samples
idx_list_filter_samples <- comp_filter_ver1(adjusted_df, 
                                            length(file_list))

# at least in 2 plastic types
idx_list_filter_plastic_types <- comp_filter_ver2(adjusted_df, 
                                                  length(unique(adjusted_df$plastic_type)))

# Combine compounds that occur in at least 2 samples
shared_comp_sample <- adjusted_df[c(idx_list_filter_samples[[1]], idx_list_filter_samples[[2]]),]

# Combine compounds that occur in at least 2 plastic types 
shared_comp_plastic_type <- adjusted_df[c(idx_list_filter_plastic_types[[1]], idx_list_filter_plastic_types[[2]]),]


# Plotting data distribution post-removal ------------------------------------------------
data_plot_post_removal <- list() 
i <- 1
for (sample in unique(other_compounds_filter_area_samples_normalized$sample_name)) {
  data_plot_post_removal[[i]] <- ggplot(data = other_compounds_filter_area_samples_normalized %>% 
                                         filter(., sample_name %in% sample),
                                       aes(x = Percent_Area)) +
    geom_histogram(bins = 50) +
    ggtitle(sample) +
    # scale_x_continuous() + # limits = c(0, 2000000)
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "hidden", axis.title = element_text(size = 5))
  i <- i + 1
}

y <- grid::textGrob("Count", rot = 90, gp = gpar(fontsize = 15))
x <- grid::textGrob("Peak Area", gp = gpar(fontsize = 15))
grid.arrange(grobs = data_plot_post_removal, ncol = 5, left = y, bottom = x)

# scipy.stats.norm.ppf function from Python 

# STEP 3.1: Examine data distribution post-normalization ===================================================

# View(whole_df %>% group_by(compound, fuel_type) %>% summarize(var(Percent_Area)))
# View(whole_df %>% group_by(compound, fuel_type) %>% summarize(var(Log_Area)))

