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
library(amap)
library(plotly)
library(umap)
library(sqldf)
library(multiway)
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

data_normalization <- function(data) {
  temp_list <- list()
  i <- 1
  # Normalize Peak Area for each sample 
  for (sample in unique(data$File)) {
    df <- shared_comp_sample[which(data$File == sample),] %>%
      # Log-based normalization
      mutate(Log_Area = log10(Area)) %>%
      mutate(Log_Height = log10(Height)) %>%
      # TSN - Percent-based normalization
      mutate(Percent_Area = Area/sum(.$Area)) %>%
      mutate(Percent_Height = Height/sum(.$Height))
    temp_list[[i]] <- df
    i <- i + 1
  }
  # Then combine data again to 1 grand data frame
  newdata <- dplyr::bind_rows(temp_list)
  return(newdata)
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
  arrange(RT) %>%
  mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                            ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                   ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                          ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                 ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                        ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other")))))))


# QUALITY CONTROL PRIOR TO STEP 1.3B ==========================================================================
# Identify linear/non-linear relationship in retention time 
# -> predict the effect of RT1 and RT2 on count of compound
# For instance, low RT1 (&RT2) "bins" will include more compounds.
# Should do this for individual fuel type sample ??????

# Supplementary Materials: Check distribution of spacing of retention time of remaining peaks --------------------
rt_dif_list <- list()
plot_dif_RT <- list()
for (n in 1:length(list_remaining_area)) {
  dif_rt <- c()
  df <- list_remaining_area[[n]] %>%
    arrange(RT)
  
  # Calculate difference in RT of each peak in a sample
  for (i in 1:(dim(df)[[1]] - 1)) {
    dif <- base::abs(base::diff(c(df[i,]$RT, df[i + 1,]$RT)))
    dif_rt <- c(dif, dif_rt)
  }
  
  rt_dif_list[[n]] <- dif_rt
  names(rt_dif_list)[n] <- paste0(unique(list_remaining_area[[n]]$File), "_diff_RT")

  plot_dif_RT[[n]] <- hist(dif_rt, breaks = 200, main = names(rt_dif_list)[n])
}

grid.arrange(grobs = plot_dif_RT, ncol = 5)



# Plot peak width against retention time ----------------------
ggplot(data = list_remaining_area[[1]]  %>%
         mutate(peak_width = .$Area/.$Height)) + 
  geom_point(aes(x = Center.X, y = peak_width)) + 
  labs(title = unique(list_remaining_area[[1]]$sample_name)) + 
  theme_bw(base_size = 15)

# Data Visualization of RT1 and RT2: -----------------------------------------
# Plot Peak width (=Peak Area/Peak Height) against RT1 and RT2 for each sample
plot_peak_width <- list()
for (i in 1:length(list_remaining_area)) {
  temp <- list_remaining_area[[i]] %>%
    mutate(peak_width = .$Area/.$Height)
  plot_peak_width[[i]] <- ggscatter(temp, x = c("RT1", "RT2"), y = "peak_width", color = "sample_name", palette = "jco")  %>% 
    invoke(ggarrange, .)
}

grid.arrange(grobs = plot_peak_width, ncol = 4)

# Plot histogram of count of RT1 and RT2 (same graph) for each sample
ggplot(data = df_remaining_area %>% 
         pivot_longer(cols = c(RT1, RT2), names_to = "RT", values_to = "value"), 
       aes(x = value, fill = RT)) +
  geom_histogram(color = '#e9ecef', alpha = 0.6, position = 'identity', bins = 120) +
  facet_wrap(~sample_name) + 
  scale_x_continuous(limits = c(3, 6))

# Pseudo-Heatmap scatter plot of RT1 RT2 as xy-axis and Peak Width as heat color intensity
plot_heatmap <- list()
i <- 1
for (item in 1:length(list_remaining_area)) {
  if (unique(list_remaining_area[[item]]$sample_name) %in% gas_only_file_list) {
    data <- list_remaining_area[[item]] %>%
      mutate(peak_width = .$Area/.$Height) # %>% filter(., peak_width <= quantile(peak_width)[2])
    plot_heatmap[[i]] <- ggplot(data, aes(x = RT1, y = RT2, color = peak_width)) +
      geom_point(size = 1.5) +
      theme_minimal() +
      scale_colour_gradientn(colors=c("red", "yellow", "green")) + 
      labs(title = unique(data$sample_name))
    i <- i + 1
  }
  else {
    next
  }
}

# grid.arrange(grobs = plot_heatmap, ncol = 3)
plot_heatmap[[1]]
# => random distribution of peak-width (not any widening pattern of peak width from low to high retention time) -> GOOD!!!


## K-means clustering on Center.X ----------------------------------------------------
x <- matrix(c(df_step1.3$Center.X), ncol = 1)
set.seed(runif(1, min = 1, max = 5))
km.out <- amap::Kmeans(x, centers = 20, nstart = 100, method = "manhattan")
plot(x = x, y = km.out$cluster, col = (km.out$cluster + 1),
      pch = 20, cex = 2)

# STEP 1.3B: Collapsing compounds based on RT1, RT2, Ion1 threshold ----------------------------------------
# Test integrity of function grouping_comp_ver1 by scrambling data frame in multiple ways
number_collapsedcomp <- c()
for (i in 1:100) {
  shuffled_df <- df_step1.3[base::sample(1:nrow(df_step1.3)), ]
  
  collapsed_shuffled_df <- grouping_comp_ver1(shuffled_df,
                                    rt1thres = 0.2,
                                    rt2thres = 0.125,
                                    ion1thres = 0.05,
                                    ion2thres = 0.05)
  
  number_collapsedcomp <- c(length(unique(collapsed_shuffled_df$collapsed_compound)), number_collapsedcomp)
}

df_step1.3_grouped <- grouping_comp_ver1(df_step1.3,
                                         rtthres = 0.05,
                                         mzthres = 0.05)

# STEP 2: Identify shared and unique compound groups across samples ------------------------------------------------
idx_list_filter_area_samples <- comp_filter_ver1(df_step1.3_grouped, 
                                                 length(file_list))

similar_compounds_filter_area_samples <- df_step1.3_grouped[idx_list_filter_area_samples[[1]],] 
other_compounds_filter_area_samples <- df_step1.3_grouped[idx_list_filter_area_samples[[2]],] 
unique_compounds_filter_area_samples <- df_step1.3_grouped[idx_list_filter_area_samples[[3]],]

# Combine similar_compounds_filter_area and other_compounds_filter_area to one data frame 
shared_comp_sample <- bind_rows(similar_compounds_filter_area_samples, other_compounds_filter_area_samples)


# STEP 3: Data Normalization ================================================================================
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


# Normalizing data accordingly to different data frames of interest --------------------------------------------------

shared_comp_normalized <- data_normalization(shared_comp_sample)

# QUALITY CONTROL STEP 3A: For each collapsed compounds we need at least 2 values of that compound for each plastic_type =======
# What is the min number of observation of collapsed compounds ?
# Create list to store number of occurence of each collasped compounds in each plastic_type
count_obs_list <- list()

for (comp in unique(shared_comp_normalized$collapsed_compound)) {
  temp <- list()
  i <- 1
  idx <- which(shared_comp_normalized$collapsed_compound == comp)
  subset_df <- shared_comp_normalized[idx,]
  for (type in unique(subset_df$plastic_type)) {
    # summary <- c()
    # summary <- c(summary, type)
    count_area <- sum(shared_comp_normalized[idx,]$plastic_type == type)
    # summary <- c(summary, count)
    temp[type] <- count_area
    i <- i + 1
  }
  count_obs_list[comp] <- temp
}


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

# Histogram Percentage-normalized Area distribution of each sample/ gas_station/ fuel_type

ggplot(data = shared_comp_normalized,
       aes(x = Percent_Area)) +
  facet_wrap(~File, scales = "free_y") + 
  geom_histogram(bins = 120) + 
  labs(x = "Percentage Area") +
  ggtitle("Percentage-normalized Area distribution") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(vjust = 0.5))
