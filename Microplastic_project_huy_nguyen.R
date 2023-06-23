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
setwd("C:/Users/huyng/OneDrive - Toronto Metropolitan University/Microplastic-Fingerprinting/CSV_Export_2023-05-23_Cleaned")

file_list <- list.files(pattern = '*.csv') %>%
  .[!str_detect(., "Blank")]

# Blank samples 
blank <- list.files(pattern = '*.csv') %>%
  .[str_detect(., "Blank")]

# Import samples to list
df_list_step1.1 <- purrr::map(file_list, read.csv)

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
add_data_normalization <- function(data) {
  temp_list <- list()
  i <- 1
  for (sample in unique(data$File)) {
    df <- data[which(data$File == sample),] %>%
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

shared_comp_normalized <- add_data_normalization(shared_comp_sample)

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

# STEP 4: Identifying potential biomarkers via PCA-HCPC and supervised classifier -----------------------------------

#@@# See how oneway-ANOVA prediction match with the data.=============================================

onewayanovamodel <- stats::lm(Percent_Area ~ fuel_type, data = df)

groups <- data.frame(fuel_type = unique(df$fuel_type))

yCI <- predict(onewayanovamodel, newdata = groups, interval = "confidence")

pred_fuelCI <- data.frame(groups, yCI)

ggplot(df, aes(x = fuel_type, y = Percent_Area, color = fuel_type)) +
  geom_violin(aes(fill = fuel_type), alpha = 0.05) +
  geom_jitter(size = 1.5, width = 0.05) +
  geom_point(mapping = aes(x = fuel_type, y = fit), 
             data = pred_fuelCI,
             size = 4) +
  geom_segment(aes(x = fuel_type, xend = fuel_type, y = lwr, yend = upr),
               data = pred_fuelCI, 
               lwd = 1.5, 
               alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Treatment group") +
  ylab("Weight (g)") +
  labs(fill = "Group", color = "Group")

# The ANOVA prediction is wrong since it underrated the Gas variation and overrated the other fuels variation
# The means were predicted correctly but the variation are wrong!!

# Two-way ANOVA for each compound group using Percent_Area AND Percent_Height ------------------------------------

df <- copy(whole_df) 

# Boxplot distribution of fuel type in each compound group
df %>% 
  ggplot(aes(x = compound, y = Percent_Area, color = fuel_type, fill = fuel_type)) +
  geom_boxplot(alpha = 0.10) +
  geom_point(position = position_jitterdodge(.1)) +
  xlab('compound') +
  ylab('Percent Area') +
  theme_bw() +
  labs(fill = "Fuel type", color = "Fuel type") +
  theme(axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 3),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank()
  )

for (row in 1:nrow(df)) {
  df[which(is.na(df[row,])), "Percent_Area"] <- runif(length(which(is.na(df[row,]))),
                                                      min = 0,
                                                      max = min(whole_df$Percent_Area))
}

for (row in 1:nrow(df)) {
  df[which(is.na(df[row,])), "Percent_Height"] <- runif(length(which(is.na(df[row,]))),
                                                        min = 0,
                                                        max = min(whole_df$Percent_Height))
}

# Two way-ANOVA with main effects only (fuel_type + compound)
twowayanovamodel1 <- stats::lm(Percent_Area ~ fuel_type + compound, data = df)

summary(twowayanovamodel1)
anova(twowayanovamodel1)

X <- model.matrix(twowayanovamodel1)

groups <- data.frame(
  with(copy(whole_df), unique(data.frame(fuel_type, compound)))
)

fuel_y_pred <- predict(
  twowayanovamodel1, newdata = groups, interval = "confidence"
)

pred_fuel <- data.frame(groups, fuel_y_pred)

# Plot model prediction without interaction term between fuel_type and compound against raw data
ggplot(copy(whole_df), aes(x = fuel_type, y = Percent_Area, color = fuel_type)) +
  geom_violin(aes(fill=fuel_type), alpha = 0.1) +
  geom_jitter(size = 1, width = 0.005) +
  geom_point(mapping = aes(x = fuel_type, y = fit), 
             data = pred_fuel,
             color = 'black',
             size = 2) +
  geom_segment(
    aes(x = fuel_type, xend = fuel_type, y = lwr, yend = upr),
    data = pred_fuel,
    color = 'black') +  
  facet_wrap(~compound) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Fuel_type") +
  ylab("Percent_Area") +
  labs(fill = "Group", color = "Group")
# RESULT: Some of the prediction of the mean of fuel type for several compound groups are not correct

# Examine interaction terms between fuel_type and compound
twowayanovamodel2 <- stats::lm(Percent_Area ~ fuel_type * compound, data = df)
summary(twowayanovamodel2)
anova(twowayanovamodel2) # -> RESULT:  there is significant in the interaction between fuel_type and 
# compound, we would say that percent_area changes differently between fuel_type within compound.

# Making prediction from model with interaction between fuel_type and compound
int_y_pred <- predict(
  twowayanovamodel2, newdata = groups, interval = "confidence"
)

int_pred <- data.frame(groups, int_y_pred)

# Plot model prediction with interaction between fuel_type and compound against raw data
ggplot(copy(whole_df), aes(x = fuel_type, y = Percent_Area, color = fuel_type)) +
  geom_violin(aes(fill=fuel_type), alpha = 0.1) +
  geom_jitter(size = 1, width = 0.005) +
  geom_point(mapping = aes(x = fuel_type, y = fit), 
             data = int_pred,
             color = 'black',
             size = 2) +
  geom_segment(
    aes(x = fuel_type, xend = fuel_type, y = lwr, yend = upr),
    data = int_pred,
    color = 'black') +  
  facet_wrap(~compound) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Fuel_type") +
  ylab("Percent_Area") +
  labs(fill = "Group", color = "Group")
# RESULT: All of our predicted means match up better with the observed data in the violins.

# Histogram and Boxplots of residuals of each fuel type in each compound group
histplots <- list()
boxplots <- list()
i <- 1
for (comp_grp in unique(whole_df$compound)) {
  means <- whole_df %>%
    filter(., compound %in% comp_grp) %>%
    group_by(fuel_type) %>%
    summarise(mean = mean(Percent_Area))
  
  resid_df <- merge(whole_df %>%
                      filter(., compound %in% comp_grp), means)
  
  resid_df$epsilon <- resid_df$Percent_Area - resid_df$mean
  temphistplots <- list()
  j <- 1
  for (fueltype in unique(resid_df$fuel_type)) {
    temphistplots[[j]] <- ggplot(resid_df %>%
                                   filter(., fuel_type %in% fueltype), aes(x = epsilon)) +
      geom_histogram() +
      ggtitle(paste0(comp_grp, "_", fueltype))
    j <- j + 1
  }
  
  histplots[[i]] <- grid.arrange(grobs = temphistplots, ncol = 2)
  boxplots[[i]] <- ggplot(resid_df, aes(x = fuel_type, y = epsilon)) +
    geom_boxplot() +
    ggtitle(comp_grp) +
    theme(axis.text.x = element_text(size = 1))
  i <- i + 1
}

grid.arrange(grobs = histplots, ncol = 5)
grid.arrange(grobs = boxplots, ncol = 5)

# Generalized Linear Models ----------------------------------
mods <- list()
# Example model call
mods[[1]] <- glm(path ~ year + hatchery + length + flow, family = binomial, data = choice)

for (i in 1:length(mods)) {
  names(mods)[i] <- as.character(mods[[i]]$call$formula)[3]
}

library(AICcmodavg)

# Make the model selection table
mod_table <- aictab(cand.set = mods, modnames = names(mods))


# PCA-HCPC test on compounds that result from pairwise_test_similar_other_compounds-----------------------

# Input data Research Question 1: Can we distinguish fuel types? ----------------------------------------------------
# Compounds occuring in Gas only OR Diesel only

# Gas only
idx_gas <- c()
for (comp in unique((shared_comp_samples %>% filter(., fuel_type %in% c("Gas", "Diesel")))$collapsed_compound)) {
  idx <- which(grepl(comp, shared_comp_samples$collapsed_compound, fixed = TRUE))
  
  if (colSums(shared_comp_samples[idx, 'fuel_type'] == "Diesel") < 1 & 
      colSums(shared_comp_samples[idx, 'fuel_type'] == "GasComp") < 1 &
      colSums(shared_comp_samples[idx, 'fuel_type'] == "DieselComp") < 1) {
  
    idx_gas <- c(idx_gas, idx)
  } else {
    next
  }
}

# Diesel only
idx_diesel <- c()
for (comp in unique(other_compounds_filter_area_samples$collapsed_compound)) {
  idx <- which(grepl(comp, other_compounds_filter_area_samples$collapsed_compound, fixed = TRUE))
  
  if (sum(str_detect(other_compounds_filter_area_samples[idx,]$fuel_type, "Gas")) < 1 & 
      sum(str_detect(other_compounds_filter_area_samples[idx,]$fuel_type, "GasComp")) < 1 &
      sum(str_detect(other_compounds_filter_area_samples[idx,]$fuel_type, "DieselComp")) < 1) {

    idx_diesel <- c(idx_diesel, idx)
  } else {
    next
  }
}

# DieselComp only
# idx_dieselcomp <- c()
# for (comp in unique(other_compounds_filter_area_samples_normalized$collapsed_compound)) {
#   idx <- which(grepl(comp, other_compounds_filter_area_samples_normalized$collapsed_compound, fixed = TRUE))
#   
#   if (sum(str_detect(other_compounds_filter_area_samples_normalized[idx,]$fuel_type, "Gas")) < 1 & 
#       sum(str_detect(other_compounds_filter_area_samples_normalized[idx,]$fuel_type, "Diesel")) < 1 &
#       sum(str_detect(other_compounds_filter_area_samples_normalized[idx,]$fuel_type, "GasComp")) < 1) {
#     
#     idx_dieselcomp <- c(idx_dieselcomp, idx)
#   } else {
#     next
#   }
# }
# RESULT: no compound found in only just DieselComp, this makes sense b/c DieselComp is a mixture of Diesel samples 
# (compound in DieselComp must/should be found in both DieselComp and Diesel samples)

# GasComp only
# idx_gascomp <- c()
# for (comp in unique(other_compounds_filter_area_samples_normalized$collapsed_compound)) {
#   idx <- which(grepl(comp, other_compounds_filter_area_samples_normalized$collapsed_compound, fixed = TRUE))
#   
#   if (sum(str_detect(other_compounds_filter_area_samples_normalized[idx,]$fuel_type, "Gas")) < 1 & 
#       sum(str_detect(other_compounds_filter_area_samples_normalized[idx,]$fuel_type, "Diesel")) < 1 &
#       sum(str_detect(other_compounds_filter_area_samples_normalized[idx,]$fuel_type, "DieselComp")) < 1) {
#     
#     idx_gascomp <- c(idx_gascomp, idx)
#   } else {
#     next
#   }
# }

# RESULT: no compound found in only just GasComp, this makes sense b/c GasComp is a mixture of Gas samples
# (compound in GasComp must/should be found in both GasComp and Gas samples)

# Create df of Only_Gas and Only_Diesel

only_gas <- other_compounds_filter_area_samples[idx_gas,] # , other_compounds_filter_area_samples[idx_diesel,]
length(unique(only_gas$collapsed_compound))

# Input data Research Question 2: Can we distinguish gas stations? ----------------------------------------------------
# Remove all GasComp, Diesel, DieselComp
gas_only_normalized <- other_compounds_filter_area_samples_normalized[idx_gas,]

# PCA with different Research Question 1 & 2 input data -------------------------------------------------------
# group_by sample_name because we want to preserve the detail variation between samples rather than combine them into bigger group of fuel_type or gas_stations

# Research question 1: Distinguish fuel type?
pcadfTSN_1 <- other_compounds_filter_area_samples_normalized %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  group_by(sample_name, collapsed_compound) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "sample_name")

# Research question 2: Distinguish gas stations?
pcadfTSN_2 <- gas_only_normalized %>%
  mutate(gas_station = factor(gas_station, levels = c(unique(gas_station)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  group_by(gas_station, collapsed_compound) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "gas_station")


### Fill in missing values with LOD
for (col in 1:ncol(pcadfTSN_1)) { 
  pcadfTSN_1[which(is.na(pcadfTSN_1[,col])), col] <- runif(length(which(is.na(pcadfTSN_1[,col]))),
                                                       min = sort(other_compounds_filter_area_samples_normalized$Percent_Area)[1],
                                                       max = sort(other_compounds_filter_area_samples_normalized$Percent_Area)[2]) 
}

for (col in 1:ncol(pcadfTSN_2)) { 
  pcadfTSN_2[which(is.na(pcadfTSN_2[,col])), col] <- runif(length(which(is.na(pcadfTSN_2[,col]))),
                                                       min = sort(other_compounds_filter_area_samples_normalized$Percent_Area)[1],
                                                       max = sort(other_compounds_filter_area_samples_normalized$Percent_Area)[2]) 
}

### PCA
res.pcaTSN_1 <- PCA(pcadfTSN_1, 
               scale.unit = TRUE, 
               graph = FALSE)

res.pcaTSN_2 <- PCA(pcadfTSN_2, 
                  scale.unit = TRUE, 
                  graph = FALSE)

# subset <- pcadfLog %>%
#   rownames_to_column(., "sample_name") %>%
#   mutate(sample_name = factor(sample_name, levels = c(unique(sample_name))))
# 
# fviz_pca_biplot(res.pcaLog,
#                 select.var = list(contrib = 20),
#                 repel = TRUE,
#                 label = "ind",
#                 labelsize = 4,
#                 habillage = subset$sample_name,
#                 dpi = 900)

### HCPC
hcpcTSN <- HCPC(res.pcaTSN_1, nb.clust = -1,
             graph = TRUE, 
             consol = TRUE)

hcpcTSN <- HCPC(res.pcaTSN_2, nb.clust = -1,
                graph = TRUE, 
                consol = TRUE)

# Which ASTM compounds are found in data frame of Research question 1: Distinguish fuel type? -----------------------------------
# Looking for ASTM match by RT1, RT2, Ion1, Ion2 from ASTM_list
View(only_gas_only_diesel_normalized)
View(ASTM_list)

# Which ASTM compounds are found in data frame of Research question 1: Distinguish fuel type? -----------------------------------
View(gas_only_normalized)
View(ASTM_list)

# Research Question 3: Beyond ASTM -> t-test and One-way ANOVA followed by "post-hoc" Tukey HSD test ----------------------

# Description: We used a one-way analysis of variance (ANOVA) on each compound group in whole_df
# to estimate the effects of fuel type on the peak area assuming a Type-I error rate of alpha = 0.05. 
# Our null hypothesis was that all fuel type means were equal:
# Rmarkdown: \mu_{Gas} = \mu_{Diesel} = \mu_{DieselComp} = \mu_{GasComp} 
# If we want to do further comparisons between groups (other than just comparing each fuel type to the Intercept by 
# themselves), then we need to add on a little “post-hoc” testing to find out which groups differ.
# One tool that allows us to make multiple comparisons between groups while adjusting for elevated Type-I error is 
# the Tukey HSD (honest significant difference) test. This test makes comparisons between each pairing of groups 
# while controlling for Type-I error

# TSN data
# Research Question 1: Beyond ASTM to distinguish different fuel_type? --------------------------------------------------------------
# Pair wise comparison ------------------------------
stat_test_fuel_type <- function(data, p_val_threshold, test_choice = "t.test", min_pairwise_samplesize = 2) {
  list_pairwise_fuel_type_TSN <- list()
  postANOVA_Tukey_fuel_type_names_list_TSN <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    df <- data %>%
      filter(., collapsed_compound %in% comp_grp)
    idx <- list()
    j <- 1
    for (fueltype in unique(df$fuel_type)) {
      idx[[j]] <- which(df$fuel_type == fueltype)
      j <- j + 1
    }
    if (length(unique(df$fuel_type)) == 1) {
      next
    } else if (length(unique(df$fuel_type)) == 2) {
      # iterates through every fuel_type
      if (any(lengths(idx) < min_pairwise_samplesize)) {
        next
      } else {
        if (test_choice == "ks") {
          # 37 comp passed
          pvalue <- ks.test(x = df[idx[[1]],]$Percent_Area,
                            y = df[idx[[2]],]$Percent_Area,
                            alternative = "two.sided")$p.value
        } else if (test_choice == "wilcox") {
          # 46 comp passed
          pvalue <- wilcox.test(x = df[idx[[1]],]$Percent_Area,
                                y = df[idx[[2]],]$Percent_Area,
                                alternative = "two.sided")$p.value
        # } else {
        #   # 215 compounds passed
        #   pvalue <- t.test(x = df[idx[[1]],]$Log_Area, # Percent_Area
        #                    y = df[idx[[2]],]$Log_Area)$p.value
        }
      }
      if (pvalue > 0.05) {
        next
      }
      else {
        # list_ttest_fuel_type_TSN[[j]] <- templist3 # names(list_ttest_fuel_type_TSN)[j] <- comp_grp
        list_pairwise_fuel_type_TSN[comp_grp] <- pvalue
      }
    } else if (length(unique(df$fuel_type)) > 2) {
      # if none of the idx has >1 observation -> skip to next compound
      if (all(lengths(idx) < min_pairwise_samplesize)) {
        next 
      } else if (sum(lengths(idx) > 1) == 2) { # if only two (out of x) sets of idx has >1 observation -> do t-test for these pairs
        if (test_choice == "ks") {
          # 37 comp passed
          pvalue <- ks.test(x = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                            y = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                            alternative = "two.sided")$p.value
        } else if (test_choice == "wilcox") {
          # 46 comp passed
          pvalue <- wilcox.test(x = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                                y = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                                alternative = "two.sided")$p.value
        } else {
          # 215 compounds passed
          pvalue <- t.test(x = df[idx[which(lengths(idx) > 1)][[1]],]$Log_Area, # Percent_Area
                           y = df[idx[which(lengths(idx) > 1)][[2]],]$Log_Area)$p.value
        }
        if (pvalue > 0.05) {
          next
        } else {
          list_pairwise_fuel_type_TSN[comp_grp] <- pvalue
        }
      } else if (sum(lengths(idx) > 1) > 2) { # if more than 2 sets of idx has >1 observation -> do one-way ANOVA followed by posthoc Tukey-HSD on them
        combined_idx <- unlist(idx[which(lengths(idx) > 1)])
        anova_df <- df[combined_idx,]
        tukey <- as.data.frame(stats::TukeyHSD(aov(Area ~ fuel_type, data = anova_df))[[1]])
        if (sum(tukey$`p adj` < 0.05) > (base::choose(length(idx), 2) - 1)) { # base::choose(length(idx), 2)
          postANOVA_Tukey_fuel_type_names_list_TSN <- c(postANOVA_Tukey_fuel_type_names_list_TSN, comp_grp)
        } else {
          next
        }
      } else {
        next
      }
    } else {
      next
    }
  }
  return(list(list_pairwise_fuel_type_TSN, postANOVA_Tukey_fuel_type_names_list_TSN))
}

# other_compounds_filter_area_samples_normalized
q1_ttest <- stat_test_fuel_type(other_compounds_filter_area_samples_normalized, p_val_threshold = 0.05, test_choice = "t.test", min_pairwise_samplesize = 2)
q1_ks <- stat_test_fuel_type(other_compounds_filter_area_samples_normalized, p_val_threshold = 0.05, test_choice = "ks", min_pairwise_samplesize = 2)
q1_wilcox <- stat_test_fuel_type(other_compounds_filter_area_samples_normalized, p_val_threshold = 0.05, test_choice = "wilcox", min_pairwise_samplesize = 2)

significant_compounds_pairwise_fuel_type_TSN <- unique(c(names(q1_ttest[[1]]), names(q1_ks[[1]]), names(q1_wilcox[[1]]),
                                                         q1_ttest[[2]], q1_ks[[2]], q1_wilcox[[2]]))
                                        
# Create data frame of significant compounds
significant_idx <- c()

for (comp in significant_compounds_pairwise_fuel_type_TSN) {
  idx <- which(other_compounds_filter_area_samples$collapsed_compound %in% comp)
  significant_idx <- c(significant_idx, idx)
}

significant_df <- other_compounds_filter_area_samples[significant_idx,]

# Is there any compound that only have 1 observation? ==> NONE
for (comp in unique(significant_df$collapsed_compound)) {
  if (length(significant_df[which(significant_df$collapsed_compound == comp),]) < 2) {
    print(comp)
  }
}

#### Which are the compound beyond ASTM in significant_df
significant_ASTM_fuel_type <- c("Compound_2038.", "Compound_2501.", "Compound_2712.", "Compound_3009.",
                                "Compound_3042.", "Compound_3395.", "Compound_4493.", "Compound_4996.",
                                "Compound_5842.", "Compound_7237.", "Compound_3136.", "Compound_3398.",
                                "Compound_3943.", "Compound_13807.", "Compound_14025.", "Compound_14215.")

beyond_astm_fuel_type_df <- significant_df %>%
  filter(., collapsed_compound %notin% significant_ASTM_fuel_type)

# Which significant ASTM passed one-way ANOVA-TukeyHSD ? What are they? 
# => Benzene, 1-ethyl-2,4-dimethyl- / Benzene, 1-methyl-4-(1-methylpropyl)- / Benzene, (1-methylpropyl)-
q1_ttest[[2]][which(q1_ttest[[2]] %in% significant_ASTM_fuel_type)]

# Research Question 1: Separated Wilcoxon and KS test for only Gas and Diesel ==================================================
# Iterating through each compound and only compare between Gas and Diesel, even if GasComp and/or DieselComp found in the compound

# comp_with_1_observation_per_sample <- c()
# comp_in_ks_and_wilcox <- c()
# comp_fail_ks_test <- c()
# comp_pass_ks_test <- c()
# comp_fail_wilcox_test <- c()
# comp_pass_wilcox_test <- c()
# p_value_list_ks <- list()
# p_value_list_wilcox <- list()
# 
# for (comp_grp in unique(shared_comp_normalized$collapsed_compound)) {
#   # Select data frame with each compound and selecting only Gas and Diesel observations
#   df <- shared_comp_normalized %>%
#     filter(., collapsed_compound %in% comp_grp) %>%
#     # remove all observation o DieselComp and GasComp -> dont care about DieselComp-GasComp or Diesel-DieselComp or Gas-GasComp
#     filter(., fuel_type %in% c("Gas", "Diesel"))
# 
#   if (length(which(df$fuel_type == "Gas")) < 2 | length(which(df$fuel_type == "Diesel")) < 2) {
#     comp_with_1_observation_per_sample <- c(comp_with_1_observation_per_sample, comp_grp)
#     next
#   } else {
#     comp_in_ks_and_wilcox <- c(comp_in_ks_and_wilcox, comp_grp)
#     pvalue_ks <- ks.test(x = df[which(df$fuel_type == "Gas"),]$Percent_Area,
#                          y = df[which(df$fuel_type == "Diesel"),]$Percent_Area,
#                          alternative = "two.sided")$p.value
# 
#     pvalue_wilcox <- wilcox.test(x = df[which(df$fuel_type == "Gas"),]$Percent_Area,
#                                  y = df[which(df$fuel_type == "Diesel"),]$Percent_Area,
#                                  alternative = "two.sided")$p.value
# 
#     p_value_list_ks[comp_grp] <- pvalue_ks
#     p_value_list_wilcox[comp_grp] <- pvalue_wilcox
#   }
#   if (pvalue_ks > 0.05 & pvalue_wilcox < 0.05) {
#     comp_fail_ks_test <- c(comp_fail_ks_test, comp_grp)
#     comp_pass_wilcox_test <- c(comp_pass_wilcox_test, comp_grp)
#   } else if (pvalue_ks < 0.05 & pvalue_wilcox > 0.05) {
#     comp_fail_wilcox_test <- c(comp_fail_wilcox_test, comp_grp)
#     comp_pass_ks_test <- c(comp_pass_ks_test, comp_grp)
#   } else if (pvalue_ks < 0.05 & pvalue_wilcox < 0.05) {
#     comp_pass_ks_test <- c(comp_pass_ks_test, comp_grp)
#     comp_pass_wilcox_test <- c(comp_pass_wilcox_test, comp_grp)
#   } else {
#     comp_fail_ks_test <- c(comp_fail_ks_test, comp_grp)
#     comp_fail_wilcox_test <- c(comp_fail_wilcox_test, comp_grp)
#     next
#   }
# }

# Create data frame for permutation - pivot wider - Research Question1 PCA
# table for the areas (rows are compound IDs, columns are sample IDs)) ------------------------------
df_X <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  select(sample_name, collapsed_compound, Percent_Area) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(sample_name, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = sample_name, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

for (col in 1:ncol(df_X)) { 
  df_X[which(is.na(df_X[,col])), col] <- runif(length(which(is.na(df_X[,col]))), # <- 0
                                               min = sort(shared_comp_normalized$Percent_Area)[1],
                                               max = sort(shared_comp_normalized$Percent_Area)[2])
}


# table for information ((rows are sample IDs, columns are sample information)-----------------------
metadata_X <- data.frame(unique((shared_comp_normalized %>%
                                filter(., fuel_type %in% c("Gas", "Diesel")))$sample_name)) 
colnames(metadata_X) <- c('sample_name')
metadata_X <- metadata_X %>%
  mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                            ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                   ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))) %>%
  # # Grouping samples into respective Gas stations
  mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                              ifelse(str_detect(sample_name, "F001"), "Station_1",
                                     ifelse(str_detect(sample_name, "F007"), "Station_7",
                                            ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                   ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                          ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
  column_to_rownames(., var = "sample_name") 

# check that sample names match exactly between pdata and expression data 
all(colnames(df_X) == rownames(metadata_X))

# Conduct principal component analysis (PCA): ---------------------------------
library(PCAtools)

p <- pca(df_X, metadata = metadata_X, removeVar = 0.1)

# A bi-plot
biplot(p,
       lab = row.names(p$metadata),
       colby = 'fuel_type',
       hline = 0, vline = 0,
       legendPosition = 'right', labSize = 5,
       showLoadings = TRUE)

# Pairs plot
pairsplot(p,
          components = getComponents(p, c(1:5)),
          triangle = FALSE,
          trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 1.5,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'fuel_type',
          title = 'Pairs plot',
          axisLabSize = 14, plotaxes = TRUE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


# explore further the collapsed_compounds that are driving these differences along each PC.
plotloadings(p,
             rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
             caption = 'Top 10% variables',
             labSize = 3)


eigencorplot(p, metavars = c('fuel_type', 'gas_station'))

# Testing significance of clusters on PCA plot -> Hotellings T-test


# Do Wilcoxon and KS test like Prof. You Liang suggested -----------------------------------------
df_X_stat <- df_X %>%
  relocate(`0220F001D.xlsx`, `0220F009D.xlsx`, `0220F009-2D.xlsx`, `0220F005D.xlsx`, .after = last_col())

p_value_list_ks <- list()
p_value_list_wilcox <- list()

for (row in 1:nrow(df_X_stat)) {
  pvalue_ks <- ks.test(x = as.numeric(df_X_stat[row,1:21]),
                       y = as.numeric(df_X_stat[row,22:25]),
                       alternative = "two.sided")$p.value

  pvalue_wilcox <- wilcox.test(x = as.numeric(df_X_stat[row,1:21]),
                               y = as.numeric(df_X_stat[row,22:25]),
                               alternative = "two.sided")$p.value

  p_value_list_ks[rownames(df_X_stat[row,])] <- pvalue_ks
  p_value_list_wilcox[rownames(df_X_stat[row,])] <- pvalue_wilcox
}

# Send data frame to Prof You Liang ---------------------------

df_X_new <- df_X %>% rownames_to_column(., var = "collapsed_compound") # df_X derived from shared_comp_normalized
write_xlsx(df_X_new, path = paste0(getwd(), "/data_table_24Nov2022.xlsx"))

# NO p-value CORRECTION: Check ASTM compounds present in compounds pass Wilcoxon test ---------
View(shared_comp_normalized %>% 
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% comp_pass_wilcox_test))

# Create data frame of p-value to do p-value adjustment
p_value_ks_df <- bind_rows(p_value_list_ks) %>% 
  pivot_longer(., cols = c(1:ncol(.)), names_to = "collapsed_compound", values_to = "p_value") %>%
  arrange(p_value)
p_value_wilcox_df <- bind_rows(p_value_list_wilcox) %>% 
  pivot_longer(., cols = c(1:ncol(.)), names_to = "collapsed_compound", values_to = "p_value") %>%
  arrange(p_value)

# P-value correction with Bonferroni method ---------------------------------------------
adjusted_pvalue_Bonferroni <- 0.05/nrow(p_value_ks_df)

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique((p_value_ks_df %>% filter(.,p_value_ks_df$p_value < adjusted_pvalue_Bonferroni))$collapsed_compound)))

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique((p_value_wilcox_df %>% filter(.,p_value_wilcox_df$p_value < adjusted_pvalue_Bonferroni))$collapsed_compound)))


# P-value correction by controlling False Discovery Rate ---------------------------------------------
fdrtest_ks <- NULL
for (n in 1:nrow(p_value_ks_df)) {
  fdrtest_ks <- c(fdrtest_ks, p_value_ks_df[n,]$p_value < (n*0.05/nrow(p_value_ks_df)))
}

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique(p_value_ks_df[fdrtest_ks,]$collapsed_compound))) # 128 significant compounds

fdrtest_wilcox <- NULL
for (n in 1:nrow(p_value_wilcox_df)) {
  fdrtest_wilcox <- c(fdrtest_wilcox, p_value_wilcox_df[n,]$p_value < (n*0.05/nrow(p_value_wilcox_df)))
}

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique(p_value_wilcox_df[fdrtest_wilcox,]$collapsed_compound))) # 119 significant compounds

# P-value correction by controlling positive False Discovery Rate ---------------------------------------------
library(qvalue)
pfdrtest_ks <- p_value_ks_df[which(qvalue(p_value_ks_df$p_value)$qvalues < 0.05),] 
pfdrtest_wilcox <- p_value_wilcox_df[which(qvalue(p_value_wilcox_df$p_value)$qvalues < 0.05),]

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique(pfdrtest_ks$collapsed_compound))) # 230 significant compounds

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique(pfdrtest_wilcox$collapsed_compound))) # 200 significant compounds

# P-value correction with ŠIDÁK-HOLM method ---------------------------------------------
sidak_holm_ks <- NULL
for (n in 0:(nrow(p_value_ks_df) - 1)) {
  sidak_holm_ks <- c(sidak_holm_ks, (1-(1-p_value_ks_df[n+1,]$p_value)^(nrow(p_value_ks_df) - n)) < 0.05)
  }

p_value_ks_df[sidak_holm_ks,] # 11 compounds

sidak_holm_wilcox <- NULL
for (n in 0:(nrow(p_value_wilcox_df) - 1)) {
  sidak_holm_wilcox <- c(sidak_holm_wilcox, (1-(1-p_value_wilcox_df[n+1,]$p_value)^(nrow(p_value_wilcox_df) - n)) < 0.05)
}

View(shared_comp_normalized %>%
       filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
       filter(., collapsed_compound %in% unique(p_value_wilcox_df[sidak_holm_wilcox,]$collapsed_compound))) # 11 compounds

# P-value correction with Hochberg-Step-up method ---------------------------------------------

# P-value correction with Westfall/Young resampling method ---------------------------------------------
library(multtest)
  
westfall_young_wilcoxon <- multtest::mt.minP(df_X, 
                                             classlabel = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1),
                                             test = "wilcoxon", 
                                             nonpara = "y") # 1.43% of all cases that the minimal p-value is smaller than 0.0001 

length(which(westfall_young_wilcoxon$adjp == 0.0143)) # 393 compounds
length(which(westfall_young_wilcoxon$adjp == 0.0528)) 



# Plot the distribution of the compounds passed t-test and Tukey-HSD test ----------------------------------------
# Box plot of compound passed t-test significant p-value < 0.05
plot_ttest_fuel_type_TSN <- list()
i <- 1

for (comp in unique(c(names(q1_ttest[[1]]), names(q1_ks[[1]])))) {
  plotdata <- shared_comp_normalized %>%
    filter(., collapsed_compound %in% comp)
  
  plot_ttest_fuel_type_TSN[[i]] <- ggplot(data = plotdata, 
                                          aes(x = fuel_type, y = Percent_Area, 
                                              color = fuel_type, fill = fuel_type)) +
    geom_boxplot(alpha = 0.15, width = .25) +
    geom_jitter(width = 0.05) +
    theme_light() +
    ggtitle(comp) +
    labs(x = NULL, y = NULL) +
    theme(axis.title.x = element_text(vjust = -1),
          axis.title.y = element_text(vjust = 3),
          legend.position = "hidden")
  
  # Iterative addition of significant codes to the graph
  # y_position <- max(plotdata$Percent_Area) + 0.0025
  # for (pair in signif_pvalue_pair) {
  #   signif_layers <- ggsignif::geom_signif(comparisons = strsplit(pair, "-"),
  #                                          annotations = "***", y_position = y_position)
  #   plot_1wayanova_compare[[element]] <- plot_1wayanova_compare[[element]] + signif_layers
  #   y_position <- y_position + 0.0025
  # }
  i <- i + 1 
}

grid.arrange(grobs = plot_ttest_fuel_type_TSN, ncol = 9, top = "T-test-TSN-normalized", left = "Percent_Area", bottom = "Fuel Type")

# Box plot of compound passed one-way ANOVA-TukeyHSD significant p-value < 0.05
# plot_ANOVA_Tukey_fuel_type_TSN <- list()
# i <- 1
# 
# for (comp in q1_ks[[2]]){
#   plotdata <- other_compounds_filter_area_samples_normalized %>%
#     filter(., collapsed_compound %in% comp)
#   
#   plot_ANOVA_Tukey_fuel_type_TSN[[i]] <- ggplot(data = plotdata, 
#                                                 aes(x = fuel_type, y = Percent_Area, 
#                                                     color = fuel_type, fill = fuel_type)) +
#     geom_boxplot(alpha = 0.15, width = .25) +
#     geom_jitter(width = 0.05) +
#     theme_light() +
#     ggtitle(comp) +
#     labs(x = NULL, y = NULL) +
#     theme(axis.title.x = element_text(vjust = -1),
#           axis.title.y = element_text(vjust = 3),
#           legend.position = "hidden")
#   i <- i + 1 
# }
# 
# grid.arrange(grobs = plot_ANOVA_Tukey_fuel_type_TSN, ncol = 4, top = "ANOVA-Tukey-HSD-TSN-normalized", left = "Percent_Area", bottom = "Fuel Type")

# Case 2: Beyond ASTM to distinguish different gas_station? ------------------------------------------------------------

stat_test_gas_station <- function(data, p_val_threshold, test_choice = "t.test") {
  list_ttest_gas_station_TSN <- list()
  postANOVA_Tukey_gas_station_names_list_TSN <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    df <- data %>%
      filter(., collapsed_compound %in% comp_grp)
    idx <- list()
    j <- 1
    for (station in unique(df$gas_station)) {
      idx[[j]] <- which(df$gas_station == station)
      j <- j + 1
    }
    if (length(unique(df$gas_station)) == 1) {
      next
    } else if (length(unique(df$gas_station)) == 2) { # if compound only found in only 2 gas stations
      # count elements in each nested  using lengths function to see if any element has only 1 value
      if (any(lengths(idx) < 2)) {
        next
      } else {
        if (test_choice == "ks") {
          # 37 comp passed
          pvalue <- ks.test(x = df[idx[[1]],]$Percent_Area,
                            y = df[idx[[2]],]$Percent_Area,
                            alternative = "two.sided")$p.value
        } else if (test_choice == "wilcox") {
          # 46 comp passed
          pvalue <- wilcox.test(x = df[idx[[1]],]$Percent_Area,
                                y = df[idx[[2]],]$Percent_Area,
                                alternative = "two.sided")$p.value
        } else {
          # 215 compounds passed
          pvalue <- t.test(x = df[idx[[1]],]$Log_Area, # Percent_Area
                           y = df[idx[[2]],]$Log_Area)$p.value
        }
      }
      if (pvalue > p_val_threshold) { # any(templist1 > 0.05) | any(templist2 > 0.05) | 
        next
      } else {
        list_ttest_gas_station_TSN[comp_grp] <- pvalue
      }
    } else if (length(unique(df$gas_station)) > 2) { # if compound only found in > 2 gas stations
      # if none of the idx has >1 observation -> skip to next compound
      if (all(lengths(idx) < 2)) {
        next
      } else if (sum(lengths(idx) > 1) == 2) { # if only two (out of x) sets of idx has >1 observation -> do t-test for these pairs
        if (test_choice == "ks") {
          # 37 comp passed
          pvalue <- ks.test(x = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                            y = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                            alternative = "two.sided")$p.value
        } else if (test_choice == "wilcox") {
          # 46 comp passed
          pvalue <- wilcox.test(x = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                                y = df[idx[which(lengths(idx) > 1)][[1]],]$Percent_Area,
                                alternative = "two.sided")$p.value
        } else {
          # 215 compounds passed
          pvalue <- t.test(x = df[idx[which(lengths(idx) > 1)][[1]],]$Log_Area, # Percent_Area
                           y = df[idx[which(lengths(idx) > 1)][[2]],]$Log_Area)$p.value
        }
        if (pvalue > p_val_threshold) {
          next
        } else {
          list_ttest_gas_station_TSN[comp_grp] <- pvalue
        }
      } else if (sum(lengths(idx) > 1) > 2) { # if more than 2 sets of idx has >1 observation -> do one-way ANOVA followed by posthoc Tukey-HSD on them
        combined_idx <- unlist(idx[which(lengths(idx) > 1)])
        anova_df <- df[combined_idx,]
        # onewayanovamodel <- stats::lm(Percent_Area ~ gas_station, data = anova_df)
        tukey <- as.data.frame(stats::TukeyHSD(aov(Percent_Area ~ gas_station, data = anova_df))[[1]])
        # Kruskal_wallis test
        pvalue <- stats::kruskal.test(Percent_Area ~ gas_station, data = anova_df)$p.value
        
        if (sum(tukey$`p adj` < p_val_threshold) > (base::choose(length(idx), 2) - 1)) {
          postANOVA_Tukey_gas_station_names_list_TSN <- c(postANOVA_Tukey_gas_station_names_list_TSN, comp_grp)
        } else {
          next
        }
      } else {
        next
      }
    }
  }
}

q2_ttest <- stat_test_gas_station(gas_only_normalized, p_val_threshold = 0.05, test_choice = "t.test")
q2_ks <- stat_test_gas_station(gas_only_normalized, p_val_threshold = 0.05, test_choice = "ks")
q2_wilcox <- stat_test_gas_station(gas_only_normalized, p_val_threshold = 0.05, test_choice = "wilcox")

significant_compounds_pairwise_gas_station_TSN <- unique(c(names(q1_ttest[[1]]), names(q1_ks[[1]]), names(q1_wilcox[[1]]),
                                                           q1_ttest[[2]], q1_ks[[2]], q1_wilcox[[2]]))

# UPDATE 16-11-2022, no compounds in gas_only_normalized pass the t-test and one-way ANOVA-TukeyHSD -> the individual compound alone cannot help us distinguish between gas stations,
# due to either missing data points or the peak area is just similar

# Create data frame of significant compounds
significant_idx <- c()

for (comp in significant_compounds_pairwise_gas_station_TSN) {
  idx <- which(other_compounds_filter_area_samples$collapsed_compound %in% comp)
  significant_idx <- c(significant_idx, idx)
}

significant_df <- other_compounds_filter_area_samples[significant_idx,]

# Box plot the Compounds that pass t-test and ANOVA-TukeyHSD
plot_ANOVA_Tukey_gas_station_TSN <- list()
i <- 1
significant_comp <- c(names(list_ttest_gas_station_TSN), postANOVA_Tukey_gas_station_names_list_TSN)

for (comp in significant_comp){
  plotdata <- gas_only_normalized %>%
    filter(., collapsed_compound %in% comp)
  
  plot_ANOVA_Tukey_gas_station_TSN[[i]] <- ggplot(data = plotdata, 
                                                  aes(x = gas_station, y = Percent_Area, 
                                                      color = gas_station, fill = gas_station)) +
    geom_boxplot(alpha = 0.15, width = .25) +
    geom_jitter(width = 0.05) +
    theme_light() +
    ggtitle(comp) +
    labs(x = NULL, y = NULL) +
    theme(axis.title.x = element_text(vjust = -1),
          axis.title.y = element_text(vjust = 3),
          legend.position = "hidden")
  
  # Iterative addition of significant codes to the graph
  # y_position <- max(plotdata$Percent_Area) + 0.0025
  # for (pair in signif_pvalue_pair) {
  #   signif_layers <- ggsignif::geom_signif(comparisons = strsplit(pair, "-"),
  #                                          annotations = "***", y_position = y_position)
  #   plot_1wayanova_compare[[element]] <- plot_1wayanova_compare[[element]] + signif_layers
  #   y_position <- y_position + 0.0025
  # }
  i <- i + 1 
}

grid.arrange(grobs = plot_ANOVA_Tukey_gas_station_TSN, ncol = 4, top = "ANOVA-Tukey-HSD-TSN-normalized", left = "Percent_Area", bottom = "Gas Station")



# Visualize distribution of compound (named by GCMS software) in each compound group that were used in oneway-ANOVA
plottingdata <- whole_df %>%
  filter(., collapsed_compound %in% names(list_1wayanovaresult)) %>% # "Compound_279."
  group_by(collapsed_compound, Compound) %>%
  tally(., name = "frequency_of_occurence")

plot <- list()
i <- 1

for (comp in unique(plottingdata$collapsed_compound)) {
  plot[[i]] <- ggplot(data = plottingdata %>%
                        filter(., collapsed_compound %in% comp)
                      
                      , aes(Compound, frequency_of_occurence)) +
    geom_col() +
    ggtitle(comp) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 2, "cm"))
  plot(plot[[i]])
  i <- i + 1
}

# Testing post one-way ANOVA + Tukey-HSD data with PCA-HCPC
# PCA with different input data -------------------------------------------------------

pcadfTSN <- postANOVA_Tukeydata_TSN %>%
  mutate(sample_name = factor(gas_station, levels = c(unique(gas_station)))) %>%
  mutate(compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  group_by(gas_station, collapsed_compound) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "gas_station")

# pcadfLog <- postANOVA_Tukeydata_Log10 %>%
#   mutate(sample_name = factor(gas_station, levels = c(unique(gas_station)))) %>%
#   mutate(compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
#   group_by(gas_station, collapsed_compound) %>%
#   summarise(across(Log_Area, median)) %>%
#   pivot_wider(names_from = collapsed_compound, values_from = Log_Area) %>%
#   column_to_rownames(., var = "gas_station")
# 
# pcadfBoxcox <- postANOVA_Tukeydata_boxcox %>%
#   mutate(sample_name = factor(gas_station, levels = c(unique(gas_station)))) %>%
#   mutate(compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
#   group_by(gas_station, collapsed_compound) %>%
#   summarise(across(boxcox_area, median)) %>%
#   pivot_wider(names_from = collapsed_compound, values_from = boxcox_area) %>%
#   column_to_rownames(., var = "gas_station")

### PCA
res.pcaTSN <- PCA(pcadfTSN, 
                  scale.unit = TRUE, 
                  graph = FALSE)

# res.pcaLog <- PCA(pcadfLog, 
#                   scale.unit = TRUE, 
#                   graph = FALSE)
# 
# res.pcaBoxcox <- PCA(pcadfBoxcox, 
#                      scale.unit = TRUE, 
#                      graph = FALSE)
### HCPC
hcpcTSN <- HCPC(res.pcaTSN, nb.clust = -1,
                # metric = "euclidean",
                # method = "complete",
                graph = TRUE, 
                consol = TRUE)

# hcpcLog <- HCPC(res.pcaLog, nb.clust = -1,
#                 # metric = "euclidean",
#                 # method = "complete",
#                 graph = TRUE, 
#                 consol = TRUE)
# 
# hcpcBoxcox <- HCPC(res.pcaBoxcox, nb.clust = -1,
#                    # metric = "euclidean",
#                    # method = "complete",
#                    graph = TRUE, 
#                    consol = TRUE)

# Random Forest --------------------------
dat <- pcadfTSN %>%
  rownames_to_column(., var = "sample_name")

dat <- dat %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name))))

set.seed(123)
model.rf <- randomForest(sample_name~., dat, ntree=1000,
                         importance=TRUE, nodesize=5)
model.rf
varImpPlot(model.rf,
           n.var = 20)


# Transpose pcadf
transpose_pca <- base::t(pcadfTSN)
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------


# Pivot wider for other type of Data Normalization ---------------------------------------------------------------------------------
# all_data_pre_norm_grouped_wider <- all_data_pre_norm_grouped %>%
#   mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
#   mutate(compound = factor(compound, levels = c(unique(compound)))) %>%
#   # for a sample, if there are multiple occurences of a compound, then impute with mean of %Area and %Height 
#   group_by(sample_name, compound) %>% # fuel_type 
#   # Here we collapse the duplicates compound by calculate the mean of Percent Area,
#   # assuming that duplicates of similar compounds has the normal distribution
#   summarise(across(Percent_Area, median)) %>%
#   pivot_wider(names_from = compound, values_from = Percent_Area)

# View(all_data_pre_norm_grouped_wider)
  
# MAYBE IRRELEVANT___Differentiation of isomers from the same compound groups in GasComp and DieselComp samples ----------------------------------------------------
# Selecting compound groups from similar_compounds df that has higher mean Percent_Area in GasComp and/or DieselComp

sim1 <- similar_compounds %>%
  group_by(compound, Compound, fuel_type, sample_name, Percent_Area) %>%
  tally(., name = "frequency_of_occurence")

for (comp in unique(sim1$compound)) {
  gasmean <- mean(sim1[which(sim1$compound == comp & 
                                sim1$fuel_type == "Gas"),]$Percent_Area)
  gascompmean <- mean(sim1[which(sim1$compound == comp & 
                                   sim1$fuel_type == "GasComp"),]$Percent_Area)
  dieselmean <- mean(sim1[which(sim1$compound == comp & 
                                  sim1$fuel_type == "Diesel"),]$Percent_Area)
  dieselcompmean <- mean(sim1[which(sim1$compound == comp & 
                                      sim1$fuel_type == "DieselComp"),]$Percent_Area)
  if (dieselcompmean > dieselmean || gascompmean > gasmean) {
    # if there are more than 1 compound in the compound group that satisfied above condition, then differentiate isomer with suffixes
    if (length(unique(sim1[which(sim1$compound == comp),]$Compound)) < 2) {
      next
    }
    else {
      i <- 1
      for (compound in unique(sim1[which(sim1$compound == comp),]$Compound)) {
        idx <- which(sim1[which(sim1$compound == comp),]$Compound == compound)
        sim1[idx,]$compound <- paste0(compound, i)
      }
    }
  }
  else {
    next
  }
}

# Selecting compound groups from other_compounds that has higher Percent_Area and occurrence in GasComp and/or DieselComp 

# has higher Percent_Area in GasComp and/or DieselComp

# Data Summary after Compound Grouping ----------------------------------------------------------------------------------------------------

# Visualize distribution of compound (named by GCMS software) in each compound group pf similar_compounds
plottingdata <- similar_compounds %>%
  # filter(., compound %in% colnames(pcainput1)) %>%
  group_by(compound, Compound) %>%
  tally(., name = "frequency_of_occurence")

plot <- list()
i <- 1
for (comp in unique(plottingdata$compound)) {
  plot[[i]] <- ggplot(data = plottingdata %>%
                        filter(., compound %in% comp)
                      , aes(Compound, frequency_of_occurence)) +
    geom_col() +
    ggtitle(comp) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 2, "cm"))
  plot(plot[[i]])
  i <- i + 1
}

# Distribution of %Area of each compound in composite sample only --> histogram plot
### Checkpoint For GasComp and DieselComp
plotdat2 <- similar_compounds %>%
  filter(., fuel_type %in% "GasComp") %>%
  mutate(., combined_names = paste0(Compound, "/", sample_name))

plotdat3 <- similar_compounds %>%
  filter(., fuel_type %in% "DieselComp") %>%
  mutate(., combined_names = paste0(Compound, "/", sample_name))

plot2 <- list()
plot3 <- list()
i <- 1

for (comp in unique(plotdat2$compound)) {
  plot2[[i]] <- ggplot(data = plotdat2 %>%
                         filter(., compound %in% comp)
                       , aes(combined_names, Percent_Area)) +
    geom_col() +
    ggtitle(paste0(comp, "GasComp")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 3, "cm"))
  
  plot3[[i]] <- ggplot(data = plotdat3 %>%
                         filter(., compound %in% comp)
                       , aes(combined_names, Percent_Area)) +
    geom_col() +
    ggtitle(paste0(comp, "DieselComp")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 3, "cm"))

  grid.arrange(plot2[[i]], plot3[[i]], ncol = 2)
  i <- i + 1
}

# 4) Distribution of RT1, RT2, Ion1
# summarydata1 <- similar_compounds %>%
#   group_by(fuel_type, compound, sample_name) %>%
#   summarise(
#     RT1 = as.double(unlist(across(RT1, mean))),
#     RT2 = as.double(unlist(across(RT2, mean))),
#     `Ion 1` = as.double(unlist(across(`Ion 1`, mean))),
#     Percent_Area = as.double(unlist(across(Percent_Area, mean)))
#     # n = n(compound)
#   )

# # Distribution of compound in different compound
# summarydata2 <- similar_compounds %>%
#   group_by(compound, Compound, fuel_type) %>%
#   summarise(count = n(Compound)) 


# Pair wise test all data
# pairwise_all_data_ttest <- pairwise_test(all_data_pre_norm_grouped, 
#                                          p_val_threshold = 0.15,
#                                          test_choice = "t.test")
# 
# pairwise_all_data_ks <- pairwise_test(all_data_pre_norm_grouped, 
#                                          p_val_threshold = 0.15,
#                                          test_choice = "ks")
# pairwise_all_data_mn <- pairwise_test(all_data_pre_norm_grouped, 
#                                          p_val_threshold = 0.15,
#                                          test_choice = "mn")
# 
# pairwise_all_data_clean_mn <- c()
# for (len in 1:length(pairwise_all_data_mn)) {
#   if (length(pairwise_all_data_mn[[len]]) > 0) { # > 5
#     pairwise_all_data_clean_mn <- c(pairwise_all_data_clean_mn, pairwise_all_data_mn[len])
#   }
# }
# 
# pairwise_all_data <- unique(names(c(pairwise_all_data_clean_ttest, 
#                                     pairwise_all_data_clean_ks,
#                                     pairwise_all_data_clean_mn)))

# Pair wise test similar compounds =====================================================
pairwise_similar_compounds_ks <- pairwise_test(similar_compounds_filter_area, 
                                               p_val_threshold = 0.05,
                                               test_choice = "ks")
pairwise_similar_compounds_mn <- pairwise_test(similar_compounds_filter_area, 
                                               p_val_threshold = 0.05,
                                               test_choice = "mn")
pairwise_similar_compounds_ttest <- pairwise_test(similar_compounds_filter_area, 
                                                  p_val_threshold = 0.05,
                                                  test_choice = "t.test")

# Combine all resulting compound groups from KS, MN and t-test
all_similar_compounds_pairwise <- unique(names(c(pairwise_similar_compounds_ks, 
                                                 pairwise_similar_compounds_mn,
                                                 pairwise_similar_compounds_ttest)))

# Violin + Boxplot of data distribution of compounds in each compound group-=-=-=-=-=-=-=-=-=-=-=-=-=
# List of Compounds in distinctive compound groups from similar compounds
plottingdata <- similar_compounds_filter_area %>%
  filter(., compound %in% all_similar_compounds_pairwise) %>%
  group_by(compound, Compound) %>%
  tally(., name = "frequency_of_occurence")

plot <- list()
i <- 1
for (comp in unique(plottingdata$compound)) {
  plot[[i]] <- ggplot(data = plottingdata %>%
                        filter(., compound %in% comp)
                      , aes(Compound, frequency_of_occurence)) +
    geom_col() +
    ggtitle(comp) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 2, "cm"))
  plot(plot[[i]])
  i <- i + 1
}
# Visualize distribution of data in each fuel_type of each compound group selected by pairwise test
ggplot(data = similar_compounds_filter_area
       %>% filter(., compound %in% all_similar_compounds_pairwise)
       , aes(fuel_type, Percent_Area, fill = fuel_type)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.25, color = "red", alpha = 0.05) +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~compound, scales = "free_y") +
  ggtitle("Data distribution of each compound group selected by pairwise test of similar compounds") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Filter with resulting compound groups from pair-wise comparisons and pivot wider for PCA
similar_compounds_wider <- similar_compounds_filter_area %>%
  filter(., compound %in% all_similar_compounds_pairwise) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound = factor(compound, levels = c(unique(compound)))) %>%
  group_by(sample_name, compound) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "sample_name")

res.pca <- PCA(similar_compounds_wider, 
               scale.unit = TRUE, 
               graph = FALSE)

subset <- similar_compounds_wider %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

fviz_pca_biplot(res.pca,
                # select.var = list(cos2 = 5),
                repel = TRUE,
                # axes = c(1,2),
                label = "ind",
                labelsize = 6,
                # col.ind = "cos2",
                # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                # generate color palette with ggpubr: https://rpkgs.datanovia.com/ggpubr/reference/get_palette.html ;
                # https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
                # palette = ggpubr::get_palette(c("#ae017e", "#dd1c77", "#2c7fb8", "#fc8d59", 
                #                                 "#c51b8a", "#2ca25f", "#e6550d", "#5ab4ac",
                #                                 "#d7191c", "#bae4bc", "#980043", "#7a0177",
                #                                 "#a50f15", "#253494", "#810f7c", "#a6d854"), 31),
                habillage = subset$sample_name,
                # title = "PCA_Biplot Similar compound groups filtering via Pair-wise test at p-value 0.05",
                dpi = 900
                )

# Hierarchical Clustering on Principle Components - Similar Compounds only
hcpc <- HCPC(res.pca, nb.clust = -1,
             metric = "euclidean",
             method = "complete",
             graph = TRUE)

# Pair wise test Other compounds =======================================================
# Filter only compounds that found in all 4 fuel_type
other_compounds_4_fueltype <- other_compounds_filter_area %>%
  group_by(compound) %>%
  filter(length(unique(fuel_type)) > 3)

# BEfore imputing missing values--> pick out compound groups that most likely to has some differences between 4 fuel types
pairwise_other_compounds_ks <- pairwise_test(other_compounds_4_fueltype, 
                                             p_val_threshold = 0.05,
                                             test_choice = "ks")
pairwise_other_compounds_mn <- pairwise_test(other_compounds_4_fueltype, 
                                             p_val_threshold = 0.05,
                                             test_choice = "mn")
pairwise_other_compounds_ttest <- pairwise_test(other_compounds_4_fueltype, 
                                                p_val_threshold = 0.05,
                                                test_choice = "ttest")

# Remove compound that have less than 6 cross comparison, since some of the compound groups do not have enough data point to run pairwise test
pairwise_other_compounds_clean_ks <- c()
for (len in 1:length(pairwise_other_compounds_ks)) {
  if (length(pairwise_other_compounds_ks[[len]]) > 5) { 
    pairwise_other_compounds_clean_ks <- c(pairwise_other_compounds_clean_ks, pairwise_other_compounds_ks[len])
  }
}
pairwise_other_compounds_clean_mn <- c()
for (len in 1:length(pairwise_other_compounds_mn)) {
  if (length(pairwise_other_compounds_mn[[len]]) > 5) { 
    pairwise_other_compounds_clean_mn <- c(pairwise_other_compounds_clean_mn, pairwise_other_compounds_mn[len])
  }
}
pairwise_other_compounds_clean_ttest <- c()
for (len in 1:length(pairwise_other_compounds_ttest)) {
  if (length(pairwise_other_compounds_ttest[[len]]) > 5) { 
    pairwise_other_compounds_clean_ttest <- c(pairwise_other_compounds_clean_ttest, pairwise_other_compounds_ttest[len])
  }
}

all_other_compounds_pairwise <- unique(names(c(pairwise_other_compounds_clean_ks, 
                                               pairwise_other_compounds_clean_mn,
                                               pairwise_other_compounds_clean_ttest)))

# Violin + Boxplot of data distribution of compounds in each compound group-=-=-=-=-=-=-=-=-=-=-=-=-=
# Visualize distribution of data in each fuel_type of each compound group before and after imputation
ggplot(data = all_data_pre_norm_grouped_filter_area
       %>% filter(., compound %in% all_other_compounds_pairwise)
       , aes(fuel_type, Percent_Area, fill = fuel_type)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.25, color = "red", alpha = 0.05) +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~compound, scales = "free_y") +
  ggtitle("other compounds before imputation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# OPTION 1: replace the missing values with LOD range ----------------------
# Create wide df contain only compounds that has some differences between 4 fuel types, but still have missing values
other_compounds_wider1 <- other_compounds_4_fueltype %>%
  filter(., compound %in% all_other_compounds_pairwise) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound = factor(compound, levels = c(unique(compound)))) %>%
  group_by(sample_name, compound) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "sample_name")

# Impute missing value with LOD
for (col in 1:ncol(other_compounds_wider1)) {
  other_compounds_wider1[which(is.na(other_compounds_wider1[,col])), col] <- runif(length(which(is.na(other_compounds_wider1[,col]))),
                                                                                   min = 0,
                                                                                   max = min(all_data_pre_norm_grouped$Percent_Area))
}

### Visualize distribution of data in each fuel_type of each compound group after LOD imputation
ggplot(data = other_compounds_wider1 
       %>% rownames_to_column(., var = "sample_name") 
       %>% pivot_longer(cols = 2:23, names_to = "compound", values_to = "Percent_Area")
       %>% mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                     ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                            ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
       , aes(fuel_type, Percent_Area, fill = fuel_type)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.25, color = "red", alpha = 0.05) +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~compound, scales = "free_y") +
  ggtitle("other compounds After LOD imputation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Initiate pca input by combining imputed LOD of other compounds with filtered similar compounds
pcainput1 <- full_join(other_compounds_wider1 %>%
                         rownames_to_column(., var = "sample_name"), 
                       similar_compounds_wider %>%
                         rownames_to_column(., var = "sample_name")) %>%
  column_to_rownames(., var = "sample_name")

# TRy PCA with these selected & imputed compound groups
res.pca1 <- PCA(pcainput1, 
               scale.unit = TRUE, 
               graph = FALSE)

# Scree plot
fviz_eig(res.pca1,
         addlabels = TRUE)
# Biplot
subset1 <- pcainput1 %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

fviz_pca_biplot(res.pca1,
                # select.var = list(cos2 = 5),
                repel = TRUE,
                # axes = c(1,2),
                label = "ind",
                labelsize = 6,
                # palette = ggpubr::get_palette(c("#ae017e", "#dd1c77", "#2c7fb8", "#fc8d59", 
                #                                 "#c51b8a", "#2ca25f", "#e6550d", "#5ab4ac",
                #                                 "#d7191c", "#bae4bc", "#980043", "#7a0177",
                #                                 "#a50f15", "#253494", "#810f7c", "#a6d854"), 31),
                habillage = subset1$sample_name,
                # title = "PCA_Biplot_LOD_Compound groups filtering via Pair-wise test at p-value 0.05",
                dpi = 900
                )

# Hierarchical Clustering on Principle Components--Other compounds_LOD + Similar compounds
hcpc <- HCPC(res.pca1,
             # metric = "euclidean",
             # method = "complete",
             # graph = TRUE,
             nb.clust = -1) 

# OPTION 2: REplace missing value with imputePCA ---------------------------
# Create wide df contain only compounds that has some differences between 4 fuel types, but still have missing values
other_compounds_wider2 <- other_compounds_4_fueltype %>%
  filter(., compound %in% all_other_compounds_pairwise) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(compound = factor(compound, levels = c(unique(compound)))) %>%
  group_by(sample_name, compound) %>%
  summarise(across(Percent_Area, median)) %>%
  pivot_wider(names_from = compound, values_from = Percent_Area) %>%
  column_to_rownames(., var = "sample_name")

# Impute missing values with imputePCA
imputed_other_compounds <- imputePCA(other_compounds_wider2,
                                     scale = TRUE,
                                     maxiter = 2000,
                                     method = "Regularized",
                                     seed = 123)

### Visualize distribution of data in each fuel_type of each compound group after imputePCA() imputation
ggplot(data = data.frame(imputed_other_compounds$completeObs) 
       %>% rownames_to_column(., var = "sample_name") 
       %>% pivot_longer(cols = 2:23, names_to = "compound", values_to = "Percent_Area")
       %>% mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                     ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                            ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
       , aes(fuel_type, Percent_Area, fill = fuel_type)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.25, color = "red", alpha = 0.05) +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~compound, scales = "free_y") +
  ggtitle("other compounds After imputePCA() imputation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# COmbine compound groups selected after pairwise test of similar_compounds and other_compounds
pcainput2 <- full_join(data.frame(imputed_other_compounds$completeObs) %>%
                        rownames_to_column(., var = "sample_name"), 
                       similar_compounds_wider %>%
                         rownames_to_column(., var = "sample_name")) %>%
  column_to_rownames(., var = "sample_name")

# TRy PCA with these selected & imputed compound groups
res.pca2 <- PCA(pcainput2, 
               scale.unit = TRUE, 
               graph = FALSE)

# Scree plot
fviz_eig(res.pca2,
         addlabels = TRUE)
# Biplot
subset2 <- pcainput2 %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

fviz_pca_biplot(res.pca2,
                # select.var = list(cos2 = 5),
                repel = TRUE,
                # axes = c(1,2),
                label = "ind",
                labelsize = 6,
                habillage = subset2$sample_name,
                # title = "PCA_Biplot_imputePCA_Compound groups filtering via Pair-wise test at p-value 0.01",
                dpi = 900
                )

# Hierarchical Clustering on Principle Components--Other compounds imputePCA + Similar conpounds
hcpc <- HCPC(res.pca2, 
             nb.clust = -1,
             metric = "euclidean",
             method = "complete",
             graph = TRUE)



# t-SNE clustering ------------------------------------------------------------------------------------------------
# REFERENCES VISUALIZATION: https://plotly.com/r/t-sne-and-umap-projections/
# https://distill.pub/2016/misread-tsne/
features <- subset(subset2, select = -c(sample_name)) # pcasubset
# subset_filterquantile_similar - produced dissimilar result to PCA on the same dataset
tsne <- tsne(features,
             initial_dims = 3, 
             k = 3, 
             perplexity = 15, # Hyperparameter: perplexity < number of data points
             max_iter = 2000
             )
             # pca = FALSE, perplexity=10, theta=0.5, dims=2,
             # check_duplicates = FALSE)

pdb <- cbind(data.frame(tsne),subset2$sample_name)
options(warn = -1)
tsne_plot <- plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, 
               color = ~subset2$sample_name) %>% 
  add_markers(size = 8) %>%
  layout( 
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'), 
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    scene =list(bgcolor = "#e5ecf6"))
tsne_plot


# UMAP clustering -------------------------------------------------------------------------------------------------
umap <- umap(features, n_components = 3, random_state = 15)

layout <- cbind(data.frame(umap[["layout"]]), subset2$sample_name)
umap_plot <- plot_ly(layout, x = ~X1, y = ~X2, z = ~X3, 
                color = ~subset2$sample_name) %>% 
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'x-axis'), 
                                   yaxis = list(title = 'y-axis'), 
                                   zaxis = list(title = 'z-axis'))) 
umap_plot

# REDUNDANT CODEs
# Regression Classification PCR/PLS-DA, etc. ----------------------------------------------------------------------
# NEED MORE DATA TO TRAIN CLASSIFICATION
library(caret)
subset2$lab_enc <- ifelse(str_detect(subset2$sample_name, "DieselCom"), 1, 
                          ifelse(str_detect(subset2$sample_name, "GasComp"), 2,
                                 ifelse(str_detect(subset2$sample_name, "D"), 3, 4)))
# move label to the front of df
subset2 <- subset2 %>%
  relocate(lab_enc, .before = 2)

# Partitioning data set 
classification_data <- subset(subset_filterquantile2, select = -c(sample_name))
inTraining <- createDataPartition(classification_data$lab_enc,
                                  p = .60, list = FALSE)
training <- classification_data[inTraining,]
testing  <- classification_data[-inTraining,]

ctrl <- trainControl(
  method = "cv",
  number = 10,
)

model <- train(lab_enc~.,
               data = training,
               method = "pls", # "lm", pls", "lasso", rf"
               preProcess = c("center", "scale", "pca"),
               trControl = ctrl)

plot(model)

test.features <- subset(testing, select = -c(lab_enc))
test.target <- subset(testing, select = lab_enc)[,1]

predictions <- predict(model, newdata = test.features)
cor(test.target, predictions) ^ 2


var <- get_pca_var(res.pca)
var_coord <- var$coord
var_contrib <- var$contrib
# var_cos2 <- var$cos2
# corrplot(var$cos2, is.corr = FALSE)


# Selecting representative Diesel compounds for each diesel sample ---------------------------------------------------------------------------
# Get coordinates of variables
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

# Quadrant 1  - dim1 [0 to 1], dim2 [0 to -0.25] - Influencer for sample 0220F001D
sample_0220F001D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.7 & Dim.2 > -0.4 & Dim.2 < -0.15)

# Quadrant 2 - dim1 [0 to 1], dim2 [0 to 0.15] - Influencer for sample 0220F009-2D 
sample_0220F009_2D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.5 & Dim.2 > 0.05 & Dim.2 < 0.1)

# Quadrant 3 - dim1 [0 to 1], dim2 [0.15 to 0.25] - Influencer for sample 0220F009D
sample_0220F009D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0.5 & Dim.1 > 0.25 & Dim.2 > 0.3 & Dim.2 < 0.4)

# Quadrant 4  - dim1 [0 to -1], dim2 [0.6 to 0.9]- Influencer for sample 0220F005D
sample_0220F005D_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.3 & Dim.1 > -0.5 & Dim.2 > 0.6 & Dim.2 < 0.9)

# Quadrant 5  - dim1 [0 to -1], dim2 [-0.5 to -0.75]- Influencer for sample 0220FDieselComp1,2
sample_0220FDieselComp1_2_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.7 & Dim.2 < -0.5 & Dim.2 > -0.75)

# Quadrant 6  - dim1 [0 to -1], dim2 [-0.75 to -1]- Influencer for sample 0220FDieselComp3
sample_0220FDieselComp3_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.15 & Dim.2 < -0.7)

# Biplot
name <- list(name = sample_0220FDieselComp3_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <- fviz_pca_biplot(subset_filterquantilePCA, 
                              select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                              repel = TRUE,
                              axes = c(1,2),
                              label = "var",
                              habillage = subset_filterquantile$sample_name,
                              # addEllipses=TRUE,
                              dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220FDieselComp3_influencer.png"),
              biplot_pca,
              height = 8,
              width = 15)

# Selecting representative Gasoline Composite compounds for each Gasoline sample ---------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord
sample_0220GasComp1_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.025 & Dim.1 < 0.06 & Dim.2 > 0)

sample_0220GasComp2_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.85 & Dim.2 < -0.425 & Dim.2 > -0.47)

sample_0220GasComp3_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.7 & Dim.2 < -0.5 & Dim.2 > -0.55)

# Biplot
name <- list(name = sample_0220GasComp2_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                              select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                              repel = TRUE,
                              axes = c(1,2),
                              label = "var",
                              habillage = subset_filterquantile$sample_name,
                              # addEllipses=TRUE,
                              dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220GasComp2_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)


# Selecting representative DieselComp compounds for each dieselcomp sample ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord
# Quadrant 1  - dim1 [0 to 1], dim2 [0 to -0.25] - Influencer for sample 0220F001D
sample_0220FDieselComp1_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.75 & Dim.2 > 0.55 & Dim.2 < 0.6)

# Quadrant 2 - dim1 [0 to 1], dim2 [0 to 0.15] - Influencer for sample 0220F009-2D 
sample_0220FDieselComp2_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.15 & Dim.1 > -0.19 & Dim.2 < -0.75)

# Quadrant 3 - dim1 [0 to 1], dim2 [0.15 to 0.25] - Influencer for sample 0220F009D
sample_0220FDieselComp3_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.75 & Dim.2 > 0.35 & Dim.2 < 0.4)

# Biplot
name <- list(name = sample_0220FDieselComp3_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                              select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                              repel = TRUE,
                              axes = c(1,2),
                              label = "var",
                              habillage = subset_filterquantile$sample_name,
                              # addEllipses=TRUE,
                              dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220FDieselComp3_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)




# Selecting representative Gasoline station 1, 3, 8 by gas 91 ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

sample_0220F00191_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0 & Dim.2 < -0.4 & Dim.2 > -0.55)

sample_0220F00391_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0 & Dim.2 < -0.4 & Dim.2 > -0.6)

sample_0220F00891_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0 & Dim.1 > -0.1 & Dim.2 > 0.75)

# Biplot
name <- list(name = sample_0220F00891_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                  select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                  repel = TRUE,
                  axes = c(1,2),
                  label = "var",
                  habillage = subset_filterquantile$sample_name,
                  # addEllipses=TRUE,
                  dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220F00891_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)

# Selecting representative Gasoline station 1, 3, 8 by gas 87 ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

sample_0220F00187_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0 & Dim.2 > 0.15 & Dim.2 < 0.25)

sample_0220F00387_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.5 & Dim.2 <  -0.6)

sample_0220F00887_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < -0.5 & Dim.2 > 0.66)

# Biplot
name <- list(name = sample_0220F00887_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                  select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                  repel = TRUE,
                  axes = c(1,2),
                  label = "var",
                  habillage = subset_filterquantile$sample_name,
                  # addEllipses=TRUE,
                  dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220F00887_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)




# Selecting representative Gasoline station 1, 3, 8 by gas 89 ---------------------------------------------------------------------------
var_coor <- get_pca_var(subset_filterquantilePCA)$coord

sample_0220F00189_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0.15 & Dim.1 < 0.35 & Dim.2 < -0.8)

sample_0220F00389_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 < 0 & Dim.2 > 0.2 & Dim.2 < 0.4)

sample_0220F00889_influencer <- var_coor %>% 
  as_tibble(rownames = "compound") %>%
  filter(Dim.1 > 0 & Dim.2 > 0.56 & Dim.2 < 0.66)

# Biplot
name <- list(name = sample_0220F00889_influencer$compound)
# fviz_cos2(subset_filterquantilePCA, choice = "var", axes = 1, top = 20)
biplot_pca <-
  fviz_pca_biplot(subset_filterquantilePCA, 
                  select.var = name, # list(cos2 = 140), # Top x active variables with the highest cos2
                  repel = TRUE,
                  axes = c(1,2),
                  label = "var",
                  habillage = subset_filterquantile$sample_name,
                  # addEllipses=TRUE,
                  dpi = 480)

ggsave(paste0(getwd(), "/PCA graphs/sample_0220F00889_influencer.png"),
       biplot_pca,
       height = 8,
       width = 15)




# PCA reserve code--------------------------------------------------------------------
# Using Percent_Area and Percent_height value on compounds
# Diesel and DieselComp cluster
# diesel_sample <- c("0220F001D.xlsx", "0220F005D.xlsx", "0220F009-2D.xlsx", "0220F009D.xlsx",
#                    "0220FDieselComp1.xlsx", "0220FDieselComp2.xlsx", "0220FDieselComp3.xlsx")
# dieselcomp_sample <- c("0220FDieselComp1.xlsx", "0220FDieselComp2.xlsx", "0220FDieselComp3.xlsx")
# 
# # All gasoline
# gas_clusall <- c(
#   "0220F00187-2.xlsx","0220F00187-3.xlsx","0220F00187.xlsx","0220F00387.xlsx", "0220F00887.xlsx"
#   ,"0220F00189.xlsx","0220F00389.xlsx", "0220F00889.xlsx"
#   ,"0220F00191.xlsx", "0220F00391.xlsx","0220F00891.xlsx"
#   ,"0220F00894.xlsx", "0220F00587.xlsx","0220F00589.xlsx","0220F00591.xlsx",
#   "0220F00787.xlsx","0220F00789.xlsx","0220F00791.xlsx",
#   "0220F00987.xlsx", "0220F00989.xlsx", "0220F00991.xlsx"
# )
# 
# # Gasoline station 1, 3, 8 cluster - cluster 1
# gas_clus1_sample <- c(
# "0220F00187-2.xlsx","0220F00187-3.xlsx","0220F00187.xlsx","0220F00387.xlsx", "0220F00887.xlsx"
# ,"0220F00189.xlsx","0220F00389.xlsx", "0220F00889.xlsx"
# ,"0220F00191.xlsx", "0220F00391.xlsx","0220F00891.xlsx"
# ,"0220F00894.xlsx", 
# )
# 
# # Gasoline station 5, 7, 9 cluster - cluster 2
# gas_clus2_sample <- c("0220F00587.xlsx","0220F00589.xlsx","0220F00591.xlsx",
#                        "0220F00787.xlsx","0220F00789.xlsx","0220F00791.xlsx",
#                        "0220F00987.xlsx", "0220F00989.xlsx", "0220F00991.xlsx")
# 
# # Checkpoint for %Area and %Height distribution of unique compounds in all_similar_compounds
# ggplot(all_similar_compounds, aes(x = Percent_Height)) + 
#   facet_wrap(~Compound, scales = "free_y") + 
#   geom_histogram(binwidth = 0.0005) 
# Checkpoint Result: OK, it seems like both %Area and %Height have centralized distribution->it"s safe to calculate mean
# Individual IL files
# Data frame for sorting percent_area & percent_height from highest to lowest
for (i in 1:length(df_list_clean)) {
  testdf <- df_list_clean[[i]] %>%
    mutate(Percent_Area = Area/sum(Area)) %>%
    mutate(Percent_Height = Height/sum(Height)) %>%
    arrange(desc(Percent_Height), desc(Percent_Area)) %>%
    # Compound column convert to rownames
    mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", Area, "-", Height))
  
  
  # # subset data based on the largest number of iteration
  for (row_num in 1:nrow(testdf)) {
    # slice data based on condition of cumulative sum of percent_height, limit ~ 80%
    if (sum(testdf[1:row_num,]$Percent_Height) > 0.99) {
      testdf <- slice_head(testdf, n = row_num)
      break
    }
  }
  
  # cov_df_scaled <- cov(scale(testdf %>%
  #                      column_to_rownames(., var = "Compound"),
  #                    center = TRUE, scale = TRUE))
  # view(cov_df_scaled)
  # Grouping compound types
  testdf <- all_subset_clean %>%                                    # ignore.case -> case insensitive
    # mutate(compound_type = ifelse(grepl("bromo", Compound, ignore.case = TRUE),"bromo",
    #                               ifelse(grepl("cyclo", Compound, ignore.case = TRUE),"cyclo",
    #                                      ifelse(grepl("chlor", Compound, ignore.case = TRUE),"chloro",
    #                                             ifelse(grepl("phosph", Compound, ignore.case = TRUE),"phospho",
    #                                                    ifelse(grepl("sulf", Compound, ignore.case = TRUE),"sulfur",
    #                                                           ifelse(grepl("amin", Compound, ignore.case = TRUE),"amine",
    #                                                                  ifelse(grepl("naphthal", Compound, ignore.case = TRUE),"naphthalene",
    #                                                                         ifelse(grepl("Benze", Compound, ignore.case = TRUE), "benzene","others"))))))))) %>%
    # filter(!grepl("others", compound_type)) %>%
    column_to_rownames(., var = "Compound")
  
  # MUST CONVERT GROUPING COLUMN to factor type, otherwise error "undefined columns selected" will happen for 'habillage'
  testdf$sample_name <- factor(testdf$sample_name, levels = c(unique(testdf$sample_name)))
  
  filter_quantile_pca <- PCA(testdf[c(3,4,11,12)], scale.unit = TRUE, graph = FALSE)
  
  # Investigate grouping of compounds
  # REFERENCE: https://pca4ds.github.io/biplot-and-pca.html
  ILR_pca <- fviz_pca_biplot(filter_quantile_pca, repel = TRUE, label = "var",
                  habillage = testdf$sample_name,
                  # palette = "Dark2",
                  addEllipses=TRUE,
                  dpi = 480)
  # Scree plot
  fviz_eig(filter_quantile_pca, addlabels = TRUE)
  
  # Cos2 aka. Quality of representation
  # Source: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials
  var <- get_pca_var(filter_quantile_pca)
  corrplot(var$cos2, is.corr = FALSE)
  
  fviz_cos2(filter_quantile_pca, choice = "var", axes = 1, top = 20)
  
  # Top variables (RT1, RT2,etc.) and compounds with highest contribution
  fviz_contrib(filter_quantile_pca, choice = "var", axes = 1) # contrib of var to PC1
  fviz_contrib(filter_quantile_pca, choice = "var", axes = 2) # contrib of var to PC2
  fviz_contrib(filter_quantile_pca, choice = "ind", axes = 2, top = 20) + # contrib of indi to PC2
    theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 3.5, "cm"))
  # ggsave(paste0(getwd(), "/PCA graphs/ILR_pca.png"),
  #        ILR_pca,
  #        height = 8,
  #        width = 15)
}

# PCA on indi_IL_type ---------------------------------------------------------------------------------------------
# Using only Percent_height value on compounds that exist in only one sample
unique_subsetdf <- all_unique_compounds %>%                                    # ignore.case -> case insensitive
  # mutate(compound_type = ifelse(grepl("bromo", Compound, ignore.case = TRUE),"bromo",
  #                               ifelse(grepl("cyclo", Compound, ignore.case = TRUE),"cyclo",
  #                                      ifelse(grepl("chlor", Compound, ignore.case = TRUE),"chloro",
  #                                             ifelse(grepl("phosph", Compound, ignore.case = TRUE),"phospho",
  #                                                    ifelse(grepl("sulf", Compound, ignore.case = TRUE),"sulfur",
  #                                                           ifelse(grepl("amin", Compound, ignore.case = TRUE),"amine",
  #                                                                  ifelse(grepl("naphthal", Compound, ignore.case = TRUE),"naphthalene",
  #                                                                         ifelse(grepl("Benze", Compound, ignore.case = TRUE), "benzene","others"))))))))) %>%
  # filter(!grepl("others", compound_type)) %>%
  mutate(Compound = paste(Compound, "-", MF, "-", RMF, "-", Area, "-", Height)) %>%
  column_to_rownames(., var = "Compound")

# MUST CONVERT GROUPING COLUMN to factor type, otherwise error "undefined columns selected" will happen for 'habillage'
unique_subsetdf$sample_name <- factor(unique_subsetdf$sample_name, levels = c(unique(unique_subsetdf$sample_name)))

# Histogram plot for Percent Height 

hist(all_similar_compounds$Percent_Height,
     # xlim = c(0, 0.002), 
     breaks = 3000)
quantile(all_similar_compounds$Percent_Height)
filter_quantile <- subset(all_similar_compounds, 
                          Percent_Height >  0.0005 &
                            Percent_Height < 0.002)
length(unique(filter_quantile$sample_name))

# Manual checkpoint for correlation of different variables (RT1,Rt2, %Area, %height, etc.) via covariance matrix 
# CAUTION!!! - If the variables are not strongly correlated (abs. covariance value must > 0.75), then there is no point to use PCA
view(cov(scale(filter_quantile[c(1:12)], center = TRUE, scale = TRUE)))
view(cor(scale(filter_quantile[c(1:12)], center = TRUE, scale = TRUE), method = "spearman"))

# Compare the correlation matrix our data to iris data
view(cor(scale(iris[c(1:4)], center = TRUE, scale = TRUE), method = "spearman"))

filter_quantile_pca <- PCA(filter_quantile[c(3,4,11,12)]
                           , scale.unit = TRUE, graph = FALSE)
# Investigate grouping of compounds- REFERENCE: https://pca4ds.github.io/biplot-and-pca.html
fviz_pca_biplot(filter_quantile_pca, repel = TRUE, label = "var",
                habillage = filter_quantile$sample_name,
                # addEllipses=TRUE,
                dpi = 480)

# Compare the PCA of our data to iris data
iris_pca <- PCA(iris[c(1:4)], scale.unit = TRUE)
fviz_pca_biplot(iris_pca, repel = TRUE, label = "var",
                habillage = iris$Species,
                dpi = 480)
# ggsave(paste0(getwd(), "/PCA graphs/dieselcomp_pca.png"),
#        pca,
#        height = 8,
#        width = 15)


# PCA with all_similar_compounds - Worse PC percentage explained variability than imputePCA -----------------------
subset_filterquantile_similar <- all_similar_compounds %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  # for a sample, if there are multiple occurences of a compound, then impute with mean of %Area and %Height
  group_by(sample_name, Compound) %>%
  summarise(across(Percent_Area, mean)) %>% # c(Percent_Area, Percent_Height) assuming that duplicates of similar compounds has the normal distribution
  # filter(sample_name %in% gas_clusall) %>%
  pivot_wider(names_from = Compound, values_from = Percent_Area)# c(Percent_Area, Percent_Height)

all_similar_compounds_PCA <- PCA(subset_filterquantile_similar[c(2:dim(subset_filterquantile_similar)[2])],
                                 scale.unit = TRUE,
                                 graph = FALSE)

# Scree plot
fviz_eig(subset_filterquantile_similar,
         addlabels = TRUE)

fviz_pca_var(subset_filterquantile_similar,
             col.var = "black",
             select.var = list(cos2 = 20))

fviz_pca_biplot(all_similar_compounds_PCA,
                select.var = list(cos2 = 50),# name, # list(cos2 = 140), # Top x active variables with the highest cos2
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset_filterquantile_similar$sample_name,
                # addEllipses=TRUE,
                dpi = 480)

# ggsave(paste0(getwd(), "/PCA graphs/Gasoline_station5_7_9.png"),
#        Gasoline_station5_7_9,
#        height = 8,
#        width = 15)



# Without imputePCA function, PCA input full of NA values - REALLY BAD CLUSTERING than imputePCA ------------------
# subset_filterquantile1 <- rownames_to_column(subset_filterquantile_removecol, "sample_name") 
# subset_filterquantile1$sample_name <- factor(subset_filterquantile1$sample_name, 
#                                              levels = c(unique(subset_filterquantile1$sample_name)))
# subset_filterquantilePCA <- PCA(subset_filterquantile1[c(2:dim(subset_filterquantile1)[2])], 
#                                 scale.unit = TRUE, 
#                                 graph = FALSE)
# 
# fviz_pca_biplot(subset_filterquantilePCA,
#                 select.var = list(cos2 = 50),# name, # list(cos2 = 140), # Top x active variables with the highest cos2
#                 repel = TRUE,
#                 axes = c(1,2),
#                 label = "ind",
#                 habillage = subset_filterquantile1$sample_name,
#                 # addEllipses=TRUE,
#                 dpi = 480)



# Reserve code ----------------------------------------------------------------------------------------------------
# all_subset_clean / all_subset_clean_height
all_similar_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0))
colnames(all_similar_compounds) <- colnames(all_subset_clean)
all_similar_compounds$Compound <- as.character(all_similar_compounds$Compound)
all_similar_compounds$sample_name <- as.character(all_similar_compounds$sample_name)
all_similar_compounds$fuel_type <- as.character(all_similar_compounds$fuel_type)

all_different_compounds <- data.frame(matrix(ncol = ncol(all_subset_clean), nrow = 0))
colnames(all_different_compounds) <- colnames(all_subset_clean)
all_different_compounds$Compound <- as.character(all_different_compounds$Compound)
all_different_compounds$sample_name <- as.character(all_different_compounds$sample_name)
all_different_compounds$fuel_type <- as.character(all_different_compounds$fuel_type)

system.time({for (compound_name in unique(all_subset_clean$Compound)) { #all_subset_clean
  # https://ashleytinsleyaileen.blogspot.com/2020/05/syntax-error-in-regexp-pattern.html?msclkid=7e2f2593d15b11ecbf9464b31d04ea64
  
  # Filter 1: baseR::grepl 
  slice_df1 <- all_subset_clean[which(grepl(compound_name, all_subset_clean$Compound, fixed = TRUE)),]
  
  # !!! --> If try() produces error, then next iteration in for loop
  if (class(try(which(grepl(paste0("^", compound_name, "$"), all_subset_clean$Compound)), silent = TRUE)) %in% "try-error") {
    # put weird compound names in unique subset
    all_different_compounds <- bind_rows(all_different_compounds, slice_df1, .id = NULL)
    next
  } 
  # Filter 2: stringr::str_which
  else {
    slice_df2 <- slice_df1[which(grepl(paste0("^", compound_name, "$"), slice_df1$Compound)),]
  }
  rm(slice_df1)
  rm(compound_name)
  # if compound was found in another sample, then append it to new dataframe: all_subset_clean_similar
  if (length(unique(slice_df2$sample_name)) > (length(indi_IL_file_list) - 1)) {
    # append the remove slice_df to all_subset_clean_similar
    all_similar_compounds <- bind_rows(all_similar_compounds, slice_df2, .id = NULL)
  }
  else {
    all_different_compounds <- bind_rows(all_different_compounds, slice_df2, .id = NULL)
  }
  rm(slice_df2)
}})

all_unique_compounds <- data.frame(matrix(ncol = ncol(all_different_compounds), nrow = 0))
colnames(all_unique_compounds) <- colnames(all_different_compounds)
all_unique_compounds$Compound <- as.character(all_unique_compounds$Compound)
all_unique_compounds$sample_name <- as.character(all_unique_compounds$sample_name)
all_unique_compounds$fuel_type <- as.character(all_unique_compounds$fuel_type)

system.time({
  for (compound_name in unique(all_different_compounds$Compound)) {
    slice_df <- all_different_compounds[which(grepl(compound_name, all_different_compounds$Compound, fixed = TRUE)),]
    if (length(unique(slice_df$sample_name)) < 2) {
      # append the remove slice_df to all_subset_clean_similar
      all_unique_compounds <- bind_rows(all_unique_compounds, slice_df, .id = NULL)
    } 
  }
}) 

length(unique(all_similar_compounds$Compound))
length(unique(all_unique_compounds$Compound))


# system.time({for (comp_grp in unique(all_subset_clean_grouped$compound)) {
#   # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
#   
#   idx <- which(grepl(paste0("^", comp_grp, "$"), all_subset_clean_grouped$compound))
#   
#   if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) > (length(indi_IL_file_list) - 1)) {
#     all_similar_compounds_idx1 <- c(all_similar_compounds_idx1, idx)
#   }
#   else if (length(unique(all_subset_clean_grouped[idx,]$sample_name)) < 2) {
#     all_unique_compounds_idx1 <- c(all_unique_compounds_idx1, idx)
#   }
#   else {
#     all_different_compounds_idx1 <- c(all_different_compounds_idx1, idx)
#   }
# }})
# Label compound type based on chemical structure - MAYBE IRRELEVANT -----------------------------------------------------------------
# REGEX REFERENCES:
# https://regenerativetoday.com/a-complete-beginners-guide-to-regular-expressions-in-r/
# https://en.wikibooks.org/wiki/R_Programming/Text_Processing#Regular_Expressions
# https://towardsdatascience.com/a-gentle-introduction-to-regular-expressions-with-r-df5e897ca432

# test1 <- df_list_clean[[1]] %>%
#   select(-ends_with(c("Area %", "Ion 1", "Ion 2", "Ion 3"))) %>%
#   mutate(compound_type = ifelse(grepl("cyclo", Compound, ignore.case = TRUE),"cyclo", # ignore.case -> case insensitive
#                                 ifelse(grepl("bromo", Compound),"bromo",
#                                        ifelse(grepl("chloro", Compound, ignore.case = TRUE),"chloro",
#                                               ifelse(grepl("phospho", Compound),"phospho",
#                                                      ifelse(grepl("sulf", Compound),"sulfur",
#                                                             ifelse(grepl("amin", Compound),"amine",
#                                                                    ifelse(grepl("naphthal", Compound),"naphthalene",
#                                                                           ifelse(grepl("Benze", Compound), "benzene", "others")))))))))


#--------------------------------------------------------------------------------------
# 28 out of 28 gasoline samples share 25 common compounds
# more than 15 out of 28 gasoline samples share 96 common compounds
# 5 out of 5 diesel samples share 50 common compounds
# more than 15 out of 39 IL samples share 150 common compounds
# 963 compounds unique for 5 diesel compounds

# examine the cumulative peak height and peak area per sample of compounds found across 39 samples
# imputed_other_compounds_long <- data.frame(imputed_other_compounds$completeObs) %>%
#   rownames_to_column(., "sample_name") %>%
#   pivot_longer(cols = c(2:length(.)),
#                names_to = "compound",
#                values_to = "Percent_Area") %>%
#   mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp",
#                             ifelse(str_detect(sample_name, "GasComp"), "GasComp",
#                                    ifelse(str_detect(sample_name, "D"), "Diesel", "Gas"))))
# 
# system.time({pairwise_imputed_other_compounds <- pairwise_test(imputed_other_compounds_long, 
#                                                                p_val_threshold = 0.1,
#                                                                test_choice = "t.test")})

# REserve code of PCA -------------------------------------------------------------------------------------------------------------
# PCA with data normalized by TSN (Percent_Area)
pcaTSN <- all_data_pre_norm_grouped_wider %>%
  column_to_rownames(., var = "sample_name") # must do before select_col and imputePCA()

# PCA with data normalized by Median Normalization
# pcaMN <- MN_wider_data %>%
#   column_to_rownames(., var = "sample_name")

# remove columns that has less than 5 unique values, including NA as a unique value
# Aka. we remove compounds that exist in less than x samples ("lower bound compound filter")
# Since Regularized approach of imputePCA drawn initial value from Gaussian distribution, we need at
# least 
# REF: https://marketing.astm.org/acton/attachment/9652/f-f77f2c0b-9bdd-43c4-b29e-a5dc68c3a4b1/1/-/-/-/-/ja17dp.pdf#:~:text=What%20is%20the%20minimum%20number%20of%20data%20points,common%20answer%20from%20most%20statistical%20professionals%20is%20%E2%80%9C30.%E2%80%9D

select_col <- c()
for (col in 1:ncol(pcaTSN)) {
  if (sum(!is.na(pcaTSN[,col])) > 30) { # select compounds that exist in every 31 samples
    select_col <- c(select_col, col)
  }
}

pcasubset_removecol <- subset(pcaTSN, select = select_col)
# selectcolcomp <- colnames(pcasubset_removecol)
# selectcolcomp %in% unique(similar_compounds$compound)

# Boxplot distribution
ggplot(data = pcasubset_removecol %>%
         rownames_to_column(., "sample_name") %>%
         pivot_longer(cols = c(2:length(.)),
                      names_to = "compound",
                      values_to = "Percent_Area") %>%
         mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                   ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                          ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))), 
       aes(fuel_type, Percent_Area)) +
  geom_boxplot() +
  facet_wrap(~compound, scales = "free_y")

# Apply imputePCA function, since PCA input cannot have NA values
PCA_impute <- imputePCA(pcasubset_removecol,
                        scale = TRUE,
                        maxiter = 2000, # need to optimise for best max iteration
                        method = "Regularized", #iterative approach-less overfitting
                        seed = 123)

# For plotting biplot later
subset2 <- pcasubset_removecol %>%
  rownames_to_column(., "sample_name") %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) 

# PCA section
pca_input <- data.frame(PCA_impute$completeObs)

res.pca <- PCA(pcasubset_removecol,  #  pca_input / pcasubset
               scale.unit = TRUE, 
               graph = FALSE)

# Scree plot
fviz_eig(res.pca,
         addlabels = TRUE)
# Biplot
fviz_pca_biplot(res.pca,
                select.var = list(cos2 = 5),# name, # Top x active variables with the highest cos2
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset2$sample_name, #  subset2 / pcasubset
                # addEllipses=TRUE,
                dpi = 900)

# Hierarchical Clustering on Principle Components
hcpc <- HCPC(res.pca, nb.clust = -1)

# Top variables (RT1, RT2,etc.) and compounds with highest contribution
fviz_contrib(res.pca, choice = "var",
             top = 1500,
             axes = 1:2) + # contrib of var to PC1 and 2
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 3.5, "cm"))

##### Extract the top 1100 compounds contribute the most to PC1:PC2
var_contrib_sorted <- data.frame(var_contrib) %>%
  rownames_to_column(., var = "Compound") %>%
  mutate_at("Dim.1", funs(sort(., decreasing = TRUE))) %>% # sort descending percent_area
  mutate_at("Dim.2",funs(sort(., decreasing = TRUE)))

var_contrib_sorted <- slice_head(df, n = 1500)

# Iterative loop removing variables Method 1: Remove_col that have less than x number of values -----------------------------------------------------
summary_list <- list()
i <- 1

system.time({for (id in 1:length(indi_IL_file_list)) {
  
  templist <- list()
  templist <- append(templist, id)
  
  select_col <- c()
  for (col in 1:ncol(pcaTSN)) {
    if (sum(!is.na(pcaTSN[,col])) > id) { # the amount of non-NA values of compounds must > x 
      select_col <- c(select_col, col)
    }
  }
  
  pcasubset_removecol <- subset(pcaTSN, select = select_col)
  
  # Apply imputePCA function, since PCA input cannot have NA values
  PCA_impute <- imputePCA(pcasubset_removecol,
                          scale = TRUE,
                          maxiter = 2000,
                          method = "Regularized", # iterative approach-less overfitting
                          seed = 123)
  
  # PCA section
  pca_input <- data.frame(PCA_impute$completeObs)
  
  res.pca <- PCA(pca_input,
                 scale.unit = TRUE, 
                 graph = FALSE)
  
  templist <- append(templist, get_eigenvalue(res.pca)[2,3]) # get the sum of PC1 & PC2
  summary_list[[i]] <- templist
  i <- i + 1
}})

# summary_list_df <- bind_cols(summary_list)
# count(n < 7)

# Iterative loop removing variables Method 2: remove one column at a time ----------------------------------
summary_list2 <- c()
i <- 1

select_col <- c()
n <- c()
for (col in 1:ncol(pcaTSN)) {
  n <- c(n, sum(!is.na(pcaTSN[,col])))
  if (sum(!is.na(pcaTSN[,col])) > 1) { # the amount of non-NA values of compounds must > x 
    select_col <- c(select_col, col)
  }
}

pcasubset_removecol <- subset(pcaTSN, select = select_col)

# Apply imputePCA function, since PCA input cannot have NA values
PCA_impute <- imputePCA(pcasubset_removecol,
                        scale = TRUE,
                        maxiter = 2000, # need to optimise for best max iteration
                        method = "Regularized", #iterative approach-less overfitting
                        seed = 123)

# system.time({for (colnum in 1:length(pcasubset_removecol)) {
# Remove one column at a time

# Check in the number of observation of compounds in 31 sample
# n <- c()
# for (col in 4771:length(pcasubset_removecol)) {
#   n <- c(n, sum(!is.na(pcasubset_removecol[,col])))}
# min(n)
# max(n)

pca_input <- data.frame(PCA_impute$completeObs)[, -c(1:4770)] 


res.pca <- PCA(pca_input,
               scale.unit = TRUE, 
               graph = FALSE)

subset2 <- rownames_to_column(pca_input,
                              "sample_name")
subset2 <- subset2 %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name))))

fviz_pca_biplot(res.pca,
                select.var = list(cos2 = 5),
                repel = TRUE,
                axes = c(1,2),
                label = "ind",
                habillage = subset2$sample_name,
                dpi = 900)

hcpc <- HCPC(res.pca, nb.clust = -1)

# Boxplot of top 101 compounds values before imputePCA
# Extract compounds from pca_input 
top101compounds <- colnames(pca_input)
top101compoundsdf <- all_data_pre_norm_grouped[which(all_data_pre_norm_grouped$compound %in% top100compounds),]
ggplot(data = top101compoundsdf[, -c(2,3,6:8)], 
       aes(fuel_type, Percent_Area)) +
  geom_boxplot() +
  facet_wrap(~compound, scales = "free_y")

# Boxplot of top 101 compounds values after imputePCA
ggplot(data = pca_input %>%
         rownames_to_column(., "sample_name") %>%
         pivot_longer(cols = c(2:length(.)),
                      names_to = "compound",
                      values_to = "Percent_Area") %>%
         mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                   ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                          ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))), 
       aes(fuel_type, Percent_Area)) +
  geom_boxplot() +
  facet_wrap(~compound, scales = "free_y")

# summary_list2 <- c(summary_list2, get_eigenvalue(res.pca)[2,3])

# }})


# When include 99% of cumulative peak height, all diesel samples share 304 compounds in common
# When include 99% of cumulative peak height, all gasoline samples share 39 compounds in common
# When include 99% of cumulative peak height, all diesel composite samples share 357 compounds in common
# When include 99% of cumulative peak height, all gasoline composite samples share 248 compounds in common
# When include 99% of cumulative peak height, all IL samples share 29 compounds in common (b4 compound grouping)
# After grouping compounds based on RT1, RT2, Ion1 threshold, all 31 IL samples share 13 compound "groups" in common 

# STEP 4 OPTION 1: Select compound groups that have widest/largest interquartile range -----------------------------------------
comp_gas_rank <- list()
for (comp in unique(whole_df$compound)) {
  comp_gas_rank[comp] <- stats::IQR(whole_df[which(whole_df$compound == comp),]$Percent_Area)
}

# STEP 4 OPTION 1: Inspect the IQR ranking of compound groups from largest:smallest IQR
View(pivot_longer(as.data.frame(comp_gas_rank), names_to = "compound", values_to = "IQR", 
                  cols = 1:length(as.data.frame(comp_gas_rank))) %>%
       arrange(desc(IQR)))

# STEP 4 OPTION 2: Univariate statistical test (pairwise test) for Gas vs. Diesel on similar_other_compounds_filter_area_all_gas_diesel-------------------------
pairwise_test_similar_other_compounds <- list()
i <- 1
for (com_grp in unique(whole_df$compound)) {
  templist <- list()
  idx1 <- which(whole_df$fuel_type == "Gas" & 
                  whole_df$compound == com_grp)
  if (length(idx1) < 2) {
    next
  }
  else {
    idx2 <- which(whole_df$fuel_type == "Diesel" & 
                    whole_df$compound == com_grp) 
    if (length(idx2) < 2) {
      next
    }
    else {
      templist[paste0("Gas", "-", "Diesel")] <- ks.test(x = whole_df[idx1,]$Percent_Area,
                                                        y = whole_df[idx2,]$Percent_Area,
                                                        alternative = "two.sided")$p.value
      
      # templist[paste0("Gas", "-", "Diesel")] <- wilcox.test(x = similar_other_compounds_filter_area_all_gas_diesel[idx1,]$Percent_Area,
      #                                                    y = similar_other_compounds_filter_area_all_gas_diesel[idx2,]$Percent_Area,
      #                                                    alternative = "two.sided")$p.value
      # 
      # 
      # templist[paste0("Gas", "-", "Diesel")] <- t.test(x = similar_other_compounds_filter_area_all_gas_diesel[idx1,]$Percent_Area,
      #                                               y = similar_other_compounds_filter_area_all_gas_diesel[idx2,]$Percent_Area)$p.value
      
    }}
  
  if (any(templist > 0.05)) {
    next
  }
  else {
    pairwise_test_similar_other_compounds[com_grp] <- templist
  }
  i <- i + 1
}

# STEP 4 OPTION 2: Confirmation Box plot of compound groups with distinct distribution of Percent_area between Gas & Diesel ==================
ggplot(data = reshape2::melt(whole_df %>%
                               filter(., compound %in% names(pairwise_test_similar_other_compounds)), id = 1:14)
       , aes(x = fuel_type, y = value)) +
  geom_violin(aes(color = variable), width = 1) +
  geom_boxplot(aes(fill = variable),
               width = 0.25) +
  facet_wrap(~compound, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# STEP 4 OPTION 2: For cross-check with ASTM target compound (compounds with wide distribution in Gas) ++++++++++++++++++++++++++++  
unique((whole_df %>%
          filter(., compound %in% names(pairwise_test_similar_other_compounds)))$Compound)

# Arson work flow reserve code ------------------------------------------------------------------------------------
library(DiagrammeR)
library(DiagrammeRsvg)
library(htmltools)

# m1 <- mermaid('
#   gantt
#     title Gantt With Custom Config
#     section Objective 1
#       check with Rich         :a1, 2015-03-09, 2d
#       add custom config param :a2, after a1, 20d
#       get bresler feedback    :a3, after a2,  2d
# ',height = 200)
# 
# # make a copy so we can compare in a tag list later
# m2 <- m1
# 
# m2$x$config = list(ganttConfig = list(
#   # a little tricky setup in what is already a hack
#   #  treat this like a filter function with filter as second component in array
#   #  and the time formatter in the first
#   #  more than likely you will want to know your scale
#   axisFormatter = list(list(
#     "%b %d, %Y" # date format to return; see d3 time formatting
#     ,htmlwidgets::JS(
#       'function(d){ return d.getDay() == 1 }' # filter for Mondays or day of week = 1
#     )
#   ))
# ))
# 
# 
# html_print(tagList(
#   tags$h1("Default Behavior")
#   ,m1
#   ,tags$h1("Hacked Behavior")
#   ,m2
# ))

DiagrammeR::grViz("digraph my_flowchart  {
      graph[splines = ortho] # to get 90 degree angles and straight lines.
      node [layout = dot, fontname = Helvetica, shape = box, fixedsize = false, width = 4, height = 1]

      node1[label = <<FONT COLOR='blue' POINT-SIZE='25'><b>Data import </b></FONT>>];
      node2[label = <<b><font color='blue' POINT-SIZE='25'>Filtering out column bleed / solvent / BTEX and <br/>minimum area limit observations (associated with GCxGC system)     </font></b>>]
      
      node4[label = <<b><font color='blue' POINT-SIZE='25'>Collapsing compounds <br/>based on RT1, RT2, Ion1 threshold      </font></b>>]
      node5[label = <<b><font color='blue' POINT-SIZE='25'>Identify shared and unique compound    <br/>groups across all 31 IL samples        </font></b>>]
      node6[label = <<b><font color='blue' POINT-SIZE='25'>Data Normalization  </font></b>>]
      node6a[label = <<b><font color='red' POINT-SIZE='25'>Total Sum Normalization   </font></b>>]
      node6b[label = <<b><font color='red' POINT-SIZE='25'>Log Normalization   </font></b>>]
      node6c[label = <<b><font color='red' POINT-SIZE='25'>Box-Cox Normalization   </font></b>>]
      node7[label = <<b><font color='blue' POINT-SIZE='25'>Statistical tests to identify potential markers      </font></b>>]
      node8[label = <<b><font color='blue' POINT-SIZE='25'>Clustering analysis      </font></b>>]
      node8a[label = <<b><font color='red' POINT-SIZE='25'>Principle Component Analysis    </font></b>>]
      node8b[label = <<b><font color='red' POINT-SIZE='25'>Hierarchical Clustering on   <br/>Principle Components     </font></b>>]
      node8c[label = <<b><font color='red' POINT-SIZE='25'>T-distributed stochastic <br/>neighbor embedding (t-SNE)     </font></b>>]
      node8d[label = <<b><font color='red' POINT-SIZE='25'>Uniform Manifold Approximation    <br/>and Projection (UMAP)</font></b>>]
      node9[label = <<b><font color='blue' POINT-SIZE='25'>Regression Classification   </font></b>>]
      node10[label = <<b><font color='blue' POINT-SIZE='25'>Quality Control <br/>(cross-validation with ASTM compound list) </font></b>>]

      blank1[label = '', width = 0.01, height = 0.01]
      blank1a[label = '', width = 0.01, height = 0.01]
      blank1b[label = '', width = 0.01, height = 0.01]
      blank1c[label = '', width = 0.01, height = 0.01]
      m1 [label = <<b><font color='red' POINT-SIZE='25'>Pair-wise comparison</font></b>>]
      m2 [label = <<b><font color='red' POINT-SIZE='25'>Multiple comparison</font></b>>]

      node7 -> blank1[dir = none];
      blank1 -> blank1a[dir = none, minlen = 12];
      {rank = same; blank1 blank1a}
      blank1a -> blank1b[dir = none, maxlen = 1];
      blank1a -> blank1c[dir = none, maxlen = 1];
      {rank = same; blank1b blank1c}
      blank1b -> m1[maxlen = 1];
      {rank = same; blank1b m1}
      blank1c -> m2[maxlen = 1];
      {rank = same; blank1c m2}
      

      m1a [label = <<b><font color='darkgreen' POINT-SIZE='25'>Kolmogorov-Smirnov test    </font></b>>]
      m1b [label = <<b><font color='darkgreen' POINT-SIZE='25'>Mann-Whitney <br/>Rank Sum test</font></b>>]
      m1c [label = <<b><font color='darkgreen' POINT-SIZE='25'>Student t-Test</font></b>>]

      blank2[label = '', width = 0.01, height = 0.01]
      m1 -> blank2[dir = none, maxlen = 1];
      blank2 -> m1a[maxlen = 1];
      blank2 -> m1b[maxlen = 1];
      blank2 -> m1c[maxlen = 1];
      {rank = same; m1a m1b m1c}

      m2a [label = <<b><font color='darkgreen' POINT-SIZE='25'>ANOVA</font></b>>]
      m2b [label = <<b><font color='darkgreen' POINT-SIZE='25'>ANCOVA</font></b>>]

      blank3[label = '', width = 0.01, height = 0.01]
      m2 -> blank3[dir = none, maxlen = 1];
      blank3 -> m2a[maxlen = 1];
      blank3 -> m2b[maxlen = 1];
      {rank = same; m2a m2b}

      blank4[label = '', width = 0.01, height = 0.01]
      node6 -> blank4[dir = none];
      blank4 -> node6a[minlen = 6];
      {rank = same; blank4 node6a};
      blank4 -> node6b[minlen = 6];
      {rank = same; blank4 node6b};
      blank4 -> node6c[minlen = 6];
      {rank = same; blank4 node6c}
      blank4 -> node7
      blank4 -> node8


      node8 -> node8a;

      node8a -> node8b;

      node8 -> node8c;

      node8 -> node8d;


     # create undirected edge from source to dummy node
      edge [dir=normal]
      node1 -> node2 -> node4 -> node5 -> node10 -> node6;
      node8 -> node9[minlen = 6];
    }", height = '100%', width = '100%')

# Grouping with k-means clustering

# factoextra::fviz_nbclust(scale(all_data_pre_norm_asc[c(4,5,9)]), kmeans, method = "wss")
# 
# kmeans_test <- stats::kmeans(x = scale(all_data_pre_norm_asc[c(4,5,9)]),
#        centers = 500,
#        iter.max = 1000,
#        nstart = 1,
#        algorithm = "Hartigan-Wong",
#        trace = FALSE)
# 
# # visualizing the clusters
# 
# fviz_cluster(kmeans_test,
#              geom = "point",
#              data = scale(all_data_pre_norm_asc[c(4,5,9)])) + 
#   labs(title = "K-means Cluster Plot",
#        fill = "Cluster",
#        shape = "Cluster",
#        color = "Cluster") + 
#   theme_bw()
# 
# # Grouping with PCA
# pca_test <- FactoMineR::PCA(all_data_pre_norm_asc[c(4,5,9)],
#                             scale.unit = TRUE, 
#                             graph = FALSE)
# 
# # Scree plot
# fviz_eig(pca_test,
#          addlabels = TRUE)
# 
# # Biplot
# fviz_pca_biplot(pca_test,
#                 repel = TRUE,
#                 dpi = 900)


# Log10 normalized data
# list_1wayanovaresult_log10 <- list()
# i <- 1
# for (comp_grp in unique(input5$collapsed_compound)) {
# df <- input5 %>%
#   filter(., collapsed_compound %in% comp_grp)
# 
# onewayanovamodel <- stats::lm(Log_Area ~ gas_station, data = df)
# list_1wayanovaresult_log10[[i]] <- onewayanovamodel
# names(list_1wayanovaresult_log10)[i] <- comp_grp
# i <- i + 1


# Box-Cox normalized data
# list_1wayanovaresult_boxcox <- list()
# i <- 1
# for (comp_grp in unique(input5$collapsed_compound)) {
#   df <- input5 %>%
#     filter(., collapsed_compound %in% comp_grp)
#   
#   onewayanovamodel <- stats::lm(boxcox_area ~ gas_station, data = df)
#   list_1wayanovaresult_boxcox[[i]] <- onewayanovamodel
#   names(list_1wayanovaresult_boxcox)[i] <- comp_grp
#   i <- i + 1
# }

# Sub-setting the compound that satisfy p-value 0.05 of one-way ANOVA & Tukey-HSD test-----------------
# TSN data
# postANOVA_Tukeynames_list_TSN <- c()
# 
# for (element in 1:length(list_1wayanovaresult_TSN)){
#   tukey <- as.data.frame(stats::TukeyHSD(aov(list_1wayanovaresult_TSN[[element]]))[[1]])
#   
#   # At least 1 pair of gas station has significant data distribution
#   if (is.na(sum(tukey$`p adj` < 0.05))) {
#     next
#   }
#   else if (sum(tukey$`p adj` < 0.05) < 2) {
#     next
#   }
#   else {
#     postANOVA_Tukeynames_list_TSN <- c(postANOVA_Tukeynames_list_TSN, names(list_1wayanovaresult_TSN[element]))
#   }
# }
# 
# postANOVA_Tukeydata_TSN <- input5 %>%
#   filter(., collapsed_compound %in% postANOVA_Tukeynames_list_TSN)

# Log10 data
# postANOVA_Tukeynames_list_Log10 <- c()
# 
# for (element in 1:length(list_1wayanovaresult_log10)){
#   tukey <- as.data.frame(stats::TukeyHSD(aov(list_1wayanovaresult_log10[[element]]))[[1]])
#   
#   # At least 1 pair of gas station has significant data distribution
#   if (is.na(sum(tukey$`p adj` < 0.05))) {
#     next
#   }
#   else if (sum(tukey$`p adj` < 0.05) < 2) {
#     next
#   }
#   else {
#     postANOVA_Tukeynames_list_Log10 <- c(postANOVA_Tukeynames_list_Log10, names(list_1wayanovaresult_log10[element]))
#   }
# }
# 
# postANOVA_Tukeydata_Log10 <- input5 %>%
#   filter(., collapsed_compound %in% postANOVA_Tukeynames_list_Log10)
# 
# # Box-Cox data
# postANOVA_Tukeynames_list_boxcox <- c()
# 
# for (element in 1:length(list_1wayanovaresult_boxcox)){
#   tukey <- as.data.frame(stats::TukeyHSD(aov(list_1wayanovaresult_boxcox[[element]]))[[1]])
#  
#   # At least 1 pair of gas station has significant data distribution
#   if (is.na(sum(tukey$`p adj` < 0.05))) {
#     next
#   }
#   else if (sum(tukey$`p adj` < 0.05) < 2) {
#     next
#   }
#   else {
#     postANOVA_Tukeynames_list_boxcox <- c(postANOVA_Tukeynames_list_boxcox, names(list_1wayanovaresult_boxcox[element]))
#   }
# }
# 
# postANOVA_Tukeydata_boxcox <- input5 %>%
#   filter(., collapsed_compound %in% postANOVA_Tukeynames_list_boxcox)

# Plot difference between each pair of gas station per compound that pass one-way ANOVA and Tukey-HSD ================
plot_1wayanova <- function(names_list, normalized_variable, left_label, title) {
  plot_1wayanova_compare_TSN <- list()
  
  for (comp in names_list){
    plotdata <- input5 %>%
      filter(., collapsed_compound %in% comp)
    
    plot_1wayanova_compare_TSN[[element]] <- ggplot(data = plotdata, 
                                                    aes(x = gas_station, y = normalized_variable, 
                                                        color = gas_station, fill = gas_station)) +
      geom_boxplot(alpha = 0.15, width = .25) +
      geom_jitter(width = 0.05) +
      # labs(fill = "Fuel Type", color = "Fuel Type") + 
      theme_light() +
      ggtitle(comp) +
      labs(x = NULL, y = NULL) +
      theme(axis.title.x = element_text(vjust = -1),
            axis.title.y = element_text(vjust = 3),
            legend.position = "hidden")
    
    # Iterative addition of significant codes to the graph
    # y_position <- max(plotdata$Percent_Area) + 0.0025
    # for (pair in signif_pvalue_pair) {
    #   signif_layers <- ggsignif::geom_signif(comparisons = strsplit(pair, "-"),
    #                                          annotations = "***", y_position = y_position)
    #   plot_1wayanova_compare[[element]] <- plot_1wayanova_compare[[element]] + signif_layers
    #   y_position <- y_position + 0.0025
    # }
    
  }
  legend <- cowplot::get_legend(ggplot(whole_df %>%
                                         filter(., collapsed_compound %in% names(list_1wayanovaresult[1])), 
                                       aes(x = fuel_type, y = Percent_Area, 
                                           color = fuel_type, fill = fuel_type)) +
                                  geom_boxplot(alpha = 0.15, width = .25) +
                                  geom_jitter(width = 0.1) +
                                  # labs(fill = "Fuel Type", color = "Fuel Type") + 
                                  theme_bw() +
                                  labs(x = NULL, y = NULL) +
                                  theme(axis.title.x = element_text(vjust = -1),
                                        axis.title.y = element_text(vjust = 3),
                                        legend.position = "right"))
  
  grid.arrange(grobs = plot_1wayanova_compare_TSN, ncol = 6, top = title,
               right = legend, left = left_label, bottom = "Gas Station")
}

plot_1wayanova(list_1wayanovaresult, Percent_Area, "Percent_Area", title = "TSN-normalized")
plot_1wayanova(list_1wayanovaresult, Log_Area , "Log_Area", title = "Log-normalized")
plot_1wayanova(list_1wayanovaresult, boxcox_area , "boxcox_area", title = "BoxCox-normalized")
# Display  the topX (top100, top50, top25, top10) quantitative variables that associated the most with each cluster-----------

# 1) TSN normalization
hcpcTSN$desc.var$quanti[[2]][1:10,] # DieselComp1:  Compound_5155. Compound_9523. Compound_10667. Compound_8162. Compound_10693.
hcpcTSN$desc.var$quanti[[3]][1:10,] # DieselComp2: Compound_1469. Compound_4108. Compound_13583. Compound_5244. Compound_7823.
hcpcTSN$desc.var$quanti[[4]][1:10,] # F001D: Compound_8278. Compound_4629. Compound_8817. Compound_7143. Compound_6207.
hcpcTSN$desc.var$quanti[[5]][1:10,] # DieselComp3: Compound_5834. Compound_9923. Compound_10396. Compound_10018. Compound_10024.
hcpcTSN$desc.var$quanti[[6]][1:10,] # F005D: Compound_9467. Compound_11781. Compound_8976. Compound_9786. Compound_7439.
hcpcTSN$desc.var$quanti[[7]][1:10,] # F009D: Compound_7980. Compound_8153. Compound_10456. Compound_6541. Compound_5896.
hcpcTSN$desc.var$quanti[[8]][1:10,] # F009D-2D: Compound_11508. Compound_10789. Compound_8282. Compound_10291. Compound_7342.

important_compounds <- row.names(hcpcTSN$desc.var$quanti[[2]][1:10,])

# Does ASTM compounds found in the top 10?
important_compounds_GCMS <- unique((whole_df %>%
                                      filter(collapsed_compound %in% important_compounds))$Compound)
# Compound name
print(important_compounds_GCMS)

# Retention time summary of those top 10 compounds
important_compounds_df <- whole_df %>%
  filter(., collapsed_compound %in% "Compound_5155.") # important_compounds

summary(important_compounds_df$RT1)
summary(important_compounds_df$RT2)

# Violin+Box plot of each compound groups in cluster ------------------------------------
plotdata_TSN <-  pcadfTSN %>%
  dplyr::select(important_compounds) %>% 
  rownames_to_column(., var = "sample_name") %>%
  pivot_longer(., cols = 2:ncol(.), 
               names_to = "collapsed_compound", 
               values_to = "Percent_Area")

# Melt data to pre-plotting if use both Area and Height
# plotdata_TSN <- reshape2::melt(pcadfTSN %>%
#                                  dplyr::select(important_compounds) %>%
#                                  rownames_to_column(., var = "sample_name") %>%
#                                  pivot_longer(., cols = 2:11, 
#                                               names_to = "collapsed_compound", 
#                                               values_to = "Percent_Area")
#                                , id = 1:2)

ggplot(data = plotdata_TSN,
       aes(x = sample_name, y = Percent_Area)) +
  geom_violin(position = position_dodge(width = 1.05)) +
  geom_boxplot(width = 0.25) +
  facet_wrap(~collapsed_compound, scales = "free_y") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# What are the GCMS name of these top10 compounds? ------------------
plottingdata <- whole_df %>%
  filter(collapsed_compound %in% important_compounds) %>%
  group_by(collapsed_compound, Compound) %>%
  tally(., name = "frequency_of_occurence")

plot <- list()
i <- 1
for (comp in unique(plottingdata$collapsed_compound)) {
  plot[[i]] <- ggplot(data = plottingdata %>%
                        filter(., collapsed_compound %in% comp)
                      , aes(Compound, frequency_of_occurence)) +
    geom_col() +
    ggtitle(comp) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 2, "cm"))
  plot(plot[[i]])
  i <- i + 1
}


