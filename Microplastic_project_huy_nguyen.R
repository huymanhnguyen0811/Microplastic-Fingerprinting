# Loading Packages --------------------------------------------------------
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(data.table)
library(Hmisc)

# Functions -------------------------------------------------------------------------------------------------------
# Grouping compounds based on RT1, RT2, and Ion1 - Version 1
grouping_comp_ver1 <- function(data, rtthres, mzthres, type) {
  
  # create empty list, each sub-list is a compound group with following criteria:
  # rtthres: RT threshold window
  # mzthres: mz threshold window
  dat <- copy(data)
  
  # Initialize the compound column filled with NA values
  dat$collapsed_compound <- NA
  i <- 1
  if (type == "ATDGCMS") {
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
        dat[idx, "collapsed_compound"] <- paste0("Compound_", i, ".ATDGCMS")
        i <- i + 1
      }  
    }
  }
  if (type == "HPLCTOFMS") {
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
        dat[idx, "collapsed_compound"] <- paste0("Compound_", i, ".HPLCTOFMS")
        i <- i + 1
      }  
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
    
    if (length(unique(data[idx,]$product_cat)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$product_cat)) < 2) {
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
# STEP 1.1: Data import --------------------------------------------

# ATDGCMS
setwd("C:/Users/huyng/OneDrive - Toronto Metropolitan University/Microplastic/Microplastic-Fingerprinting/data/ATDGCMS")

file_list1 <- list.files(pattern = '*.csv') %>%
  .[!str_detect(., "Blank")]
  # .[!str_detect(., "_USE")] # exclude environmental samples

# Blank samples 
blank_list1 <- list.files(pattern = '*.csv') %>%
  .[str_detect(., "Blank")]

# Import samples to list
df_list1_step1.1 <- purrr::map(file_list1, read.csv)

df_list1_blank <- purrr::map(blank_list1, read.csv)

df_blank1 <- dplyr::bind_rows(df_list1_blank)

# Sample information ATDGCMS
sampleinfo1 <- readxl::read_excel(paste0(getwd(), '/SampleInfo.xlsx'))
colnames(sampleinfo1)[1] <- 'File'

sampleinfo1$`Collection Date (YYYY-MM-DD)` <- as.Date(as.numeric(sampleinfo1$`Collection Date (YYYY-MM-DD)`), origin = "1899-12-30")

# HPLCTOFMS
setwd("C:/Users/huyng/OneDrive - Toronto Metropolitan University/Microplastic/Microplastic-Fingerprinting/data/HPLCTOFMS")

file_list2 <- list.files(pattern = '*.xls') %>%
  .[!str_detect(., "Blank")] %>%
  .[!str_detect(., "Info")]
# .[!str_detect(., "_USSB")] # exclude environmental samples

# Blank samples 
blank_list2 <- list.files(pattern = '*.xls') %>%
  .[str_detect(., "Blank")]

# Import samples to list
df_list2_step1.1 <- purrr::map(file_list2, read_xls, skip = 1)

df_list2_blank <- purrr::map(blank_list2, read_xls, skip = 1)

df_blank2 <- dplyr::bind_rows(df_list2_blank)

# Sample information HPLCMS
sampleinfo2 <- readxl::read_excel(paste0(getwd(), '/Plastic Product Info EF.xlsx'))
colnames(sampleinfo2)[1] <- 'Sample_ID'
sampleinfo2 <- sampleinfo2 %>%
  filter(str_detect(Sample_ID, "USE"))

# STEP 1.2B Filtering out limit of observations---------------------------- 

# ATDGCMS

list1_remaining_area <- list()
list1_removed_area <- list()
for (i in 1:length(df_list1_step1.1)) {
  list1_remaining_area[[i]] <- df_list1_step1.1[[i]] %>%
    filter(., Area > 100000)
  
  list1_removed_area[[i]] <- df_list1_step1.1[[i]] %>%
    filter(., Area <= 100000)
}

# HPLCTOFMS
list2_remaining_area <- list()
list2_removed_area <- list()
for (i in 1:length(df_list2_step1.1)) {
  list2_remaining_area[[i]] <- df_list2_step1.1[[i]] %>%
    filter(., Height > 5000)
  
  list2_removed_area[[i]] <- df_list2_step1.1[[i]] %>%
    filter(., Height <= 5000)
}

# STEP 1.3: Grouping compounds based on Retention time and molecular ions  -----------------------------------------------------------------------
# STEP 1.3A: Generate 1 grand data frame

# ATDGCMS
df1_step1.3 <- bind_rows(list1_remaining_area) %>%
  select(-c("Start", "End", "Width", "Base.Peak")) %>%
  mutate(product_cat = ifelse(str_detect(File, "Balloons"), "Toys", 
                              ifelse(str_detect(File, "FPW_"), "Food contact materials",
                                     ifelse(str_detect(File, "Pbal_Sample"), "Toys",
                                            ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                                   ifelse(str_detect(File, "PBBC_"), "Food contact materials",
                                                          ifelse(str_detect(File, "Pbag_"),"Food contact materials",
                                                                 ifelse(str_detect(File, "PDS_Sample"),"Food contact materials", 
                                                                        ifelse(str_detect(File, "Pcut_Sample"), "Food contact materials",
                                                                               ifelse(str_detect(File, "PC_Sample"), "Food contact materials",
                                                                                      ifelse(str_detect(File, "Cigs_"), "Cigarettes", 
                                                                                             ifelse(str_detect(File, "Cmat"), "Construction materials",
                                                                                                    ifelse(str_detect(File, "Mask_Sample"), "Clothes", "Misc")))))))))))))

df_blank1 <- df_blank1 %>%
  select(-c("Start", "End", "Width", "Base.Peak")) %>%
  mutate(product_cat = "Blanks")

combined_df1 <- rbind(df1_step1.3, df_blank1) %>% arrange(RT)

# HPLCTOFMS
df2_step1.3 <- bind_rows(list2_remaining_area) %>%
  select(c("m/z", "RT", "Height", "File")) %>%
  mutate(type = "Sample")


df_blank2 <- df_blank2 %>%
  select(c("m/z", "RT", "Height", "File")) %>%
  mutate(type = "Blanks")

combined_df2 <- rbind(df2_step1.3, df_blank2) %>% 
  arrange(RT)

colnames(combined_df2)[[1]] <- "m.z"

# STEP 1.3B: Collapsing compounds based on RT1, RT2, Ion1 threshold

# Statistical Description for selecting rtthres and mzthres
Hmisc::describe(combined_df1)

# ATDGCMS
combined_df1_grouped <- grouping_comp_ver1(combined_df1,
                                           rtthres = 0.05,
                                           mzthres = 0.05,
                                           type = "ATDGCMS")

# Statistical Description for selecting rtthres and mzthres
Hmisc::describe(combined_df2)

# HPLCTOFMS
combined_df2_grouped <- grouping_comp_ver1(combined_df2,
                                           rtthres = 0.1,
                                           mzthres = 0.00003,
                                           type = "HPLCTOFMS")

# Step 2: Readjust compound RA (sample) by average blank RA ======================

# ATDGCMS

# Create list to store temp dfs
temp_list <- list()
i <- 1
# Iterate through each collapsed_compound
for (comp in unique(comp_normalized1$collapsed_compound)) {
  temp <- comp_normalized1[which(comp_normalized1$collapsed_compound == comp),]
  # if compound does not exist in blanks then skip the compounds
  if (identical(which(temp$product_cat == "Blanks"), integer(0))) {
    temp_list[[i]] <- temp
    i <- i + 1
    next
  }
  else {
    # Calculate avg_blank for that compound across all blanks
    avg_blank <- mean(temp[which(temp$product_cat == "Blanks"),]$Area)
    temp <- temp[which(temp$product_cat != "Blanks"),]
    # iterate through each sample
    for (sample in unique(temp$File)) {
      # Adjust RA for each compound of each sample = RA (sample) - avg_blank
      temp[which(temp$File == sample),]$Area <- temp[which(temp$File == sample),]$Area - avg_blank
    } 
  }
  # Append current temp df to temp_list
  temp_list[[i]] <- temp
  i <- i + 1
}

adjusted_df1 <- bind_rows(temp_list)

# HPLCTOFMS

# Create list to store temp dfs
temp_list <- list()
i <- 1
# Iterate through each collapsed_compound
for (comp in unique(comp_normalized2$collapsed_compound)) {
  temp <- comp_normalized2[which(comp_normalized2$collapsed_compound == comp),]
  # if compound does not exist in blanks then skip the compounds
  if (identical(which(temp$type == "Blanks"), integer(0))) {
    temp_list[[i]] <- temp
    i <- i + 1
    next
  }
  else {
    # Calculate avg_blank for that compound across all blanks
    avg_blank <- mean(temp[which(temp$type == "Blanks"),]$Height)
    temp <- temp[which(temp$type != "Blanks"),]
    # iterate through each sample
    for (sample in unique(temp$File)) {
      # Adjust RA for each compound of each sample = RA (sample) - avg_blank
      temp[which(temp$File == sample),]$Height <- temp[which(temp$File == sample),]$Height - avg_blank
    } 
  }
  # Append current temp df to temp_list
  temp_list[[i]] <- temp
  i <- i + 1
}

adjusted_df2 <- bind_rows(temp_list) %>%
  mutate(product_cat = ifelse(str_detect(File, "USE-01"), "Food contact materials",
                              ifelse(str_detect(File, "USE-02"), "Mixed_Plastic_Waste",
                                     ifelse(str_detect(File, "USE-03"), "Food contact materials",
                                            ifelse(str_detect(File, "USE-05"), "Cigarettes",
                                                   ifelse(str_detect(File, "USE-07"),"Food contact materials",
                                                          ifelse(str_detect(File, "USE-09"),"Food contact materials", 
                                                                 ifelse(str_detect(File, "USE-11"),"Toys", 
                                                                        ifelse(str_detect(File, "USE-13"),"Food contact materials", 
                                                                               ifelse(str_detect(File, "USE-14"),"Food contact materials","Misc")))))))))) %>%
  mutate(Sample_ID = ifelse(str_detect(File, "USE-01"), "USE-01",
                            ifelse(str_detect(File, "USE-02"), "USE-02",
                                   ifelse(str_detect(File, "USE-03"), "USE-03",
                                          ifelse(str_detect(File, "USE-05"), "USE-05",
                                                 ifelse(str_detect(File, "USE-06"), "USE-06",
                                                        ifelse(str_detect(File, "USE-07"), "USE-07",
                                                               ifelse(str_detect(File, "USE-09"), "USE-09",
                                                                      ifelse(str_detect(File, "USE-11"), "USE-11",
                                                                             ifelse(str_detect(File, "USE-12"), "USE-12",
                                                                                    ifelse(str_detect(File, "USE-13"), "USE-13",
                                                                                           ifelse(str_detect(File, "USE-14"), "USE-14", "USE-15"))))))))))))


# STEP 3: Normalizing data accordingly to different data frames of interest on with the one with positive Area values -----------------------------------

# ATDGCMS
temp_list <- list()
i <- 1
# Normalize Peak Area for each sample 
for (sample in unique(adjusted_df1$File)) {
  df <- adjusted_df1[which(adjusted_df1$File == sample),] %>%
    filter(Area > 0) %>%
    mutate(Percent_Area = Area/sum(.$Area))
  temp_list[[i]] <- df
  i <- i + 1
}
# Then combine data again to 1 grand data frame
comp_normalized1 <- dplyr::bind_rows(temp_list)

# HPLCTOFMS
temp_list <- list()
i <- 1
# Normalize Peak Area for each sample 
for (sample in unique(adjusted_df2$File)) {
  df <- adjusted_df2[which(adjusted_df2$File == sample),] %>%
    filter(Height > 0) %>%
    mutate(Percent_Height = Height/sum(.$Height))
  temp_list[[i]] <- df
  i <- i + 1
}
# Then combine data again to 1 grand data frame
comp_normalized2 <- dplyr::bind_rows(temp_list)


# STEP 4: Identify shared and unique compound groups across samples ------------------------------------------------
# at least in 2 samples
# # ATDGCMS
# idx_list_filter_samples1 <- comp_filter_ver1(adjusted_df1, 
#                                             length(file_list1))
# 
# # HPLCTOFMS
# idx_list_filter_samples2 <- comp_filter_ver1(adjusted_df2, 
#                                             length(file_list2))

# Combine compounds that occur in at least 2 samples
# # ATDGCMS
# shared_comp_sample1 <- adjusted_df1[c(idx_list_filter_samples1[[1]], idx_list_filter_samples1[[2]]),]
# 
# # HPLCTOFMS
# shared_comp_sample2 <- adjusted_df2[c(idx_list_filter_samples2[[1]], idx_list_filter_samples2[[2]]),]

# at least in 2 product cats
# ATDGCMS
idx_list_filter_product_cat1 <- comp_filter_ver2(comp_normalized1, 
                                                   length(unique(comp_normalized1$product_cat)))

# HPLCTOFMS
idx_list_filter_product_cat2 <- comp_filter_ver2(comp_normalized2, 
                                                   length(unique(comp_normalized2$product_cat)))

# Combine compounds that occur in at least 2 plastic types
# ATDGCMS
shared_comp_product_cat1 <- comp_normalized1[c(idx_list_filter_product_cat1[[1]],
                                           idx_list_filter_product_cat1[[2]]
)
,]

# HPLCTOFMS
shared_comp_product_cat2 <- comp_normalized2[c(idx_list_filter_product_cat2[[1]],
                                           idx_list_filter_product_cat2[[2]]
)
,]

# Step 5: Merging Sample info with shared df ===========================================
# ATDGCMS
merge_df1 <- dplyr::full_join(x = sampleinfo1,  y = shared_comp_product_cat1, by = 'File') %>%
  filter(., !is.na(collapsed_compound))

merge_df1$File <- gsub("_", "-", merge_df1$File)

new_merge_df1 <- merge_df1 %>%
  select('File', 'product_cat', 'collapsed_compound', 'Percent_Area')
colnames(new_merge_df1)[4] <- 'Values'

# HPLCTOFMS
merge_df2 <- dplyr::full_join(x = sampleinfo2, y = shared_comp_product_cat2, by = 'Sample_ID') %>%
  filter(., !is.na(collapsed_compound))

# grep(pattern = paste(unique(merge_df2$Sample_ID), collapse = "|"), x = merge_df1$File, value = TRUE)

new_merge_df2 <- merge_df2 %>%
  select('File', 'product_cat', 'collapsed_compound', 'Percent_Height')
colnames(new_merge_df2)[4] <- 'Values'

gc_hplc <- rbind(new_merge_df1, new_merge_df2)









# Plotting data distribution pre-removal------------------
data_plot_pre_removal <- list() # NOTE: Data sets are all heavy left-skewed
for (i in 1:30) { # length(df_list)
  filter_area <- df_list2_step1.1[[i]]
  data_plot_pre_removal[[i]] <- ggplot(data = filter_area,
                                       aes(x = Height)) +
    geom_histogram(bins = 100) +
    ggtitle(file_list2[[i]]) +
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
library(grid)
library(gridExtra)

plot_a <- list()
j <- 1
for (i in 1:25) { # length(df_list)
  coverage <- c()
  for (threshold in c(seq(from = 0, to = 200000, by = 50000))) {
    df_filter_area <- df_list1_step1.1[[i]] %>%
      filter(., Area > threshold)
    coverage <- c(coverage, sum(df_filter_area$Area)*100/sum(df_list1_step1.1[[i]]$Area))
  }
  df <- data.frame(thres = seq(from = 0,to = 200000, by = 50000), cover = coverage)
  plot_a[[j]] <- ggplot(data = df,
                        aes(x = thres, y = cover)) +
    geom_col() +
    theme(text = element_text(size = 15)) +
    geom_text(aes(label = round(cover, digits = 3)), color = "green", angle = 90, hjust = 1, size = 4) +
    scale_x_continuous(breaks = seq(from = 0, to = 200000, by = 50000),
                       # remove space between plotted data and xy-axes
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, by = 25), 
                       # remove space between plotted data and xy-axes
                       expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45)) +
    ggtitle(file_list1[[i]]) +
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
for (i in 1:25) { # length(df_list_step1.1)
  peak_remain <- c()
  for (threshold in c(seq(from = 0, to = 200000, by = 50000))) {
    df_filter_area <- df_list1_step1.1[[i]] %>%
      filter(., Area > threshold)
    peak_remain <- c(peak_remain, dim(df_filter_area)[1])
  }
  df <- data.frame(thres = seq(from = 0, to = 200000, by = 50000), remain = peak_remain)
  plot_b[[j]] <- ggplot(data = df,
                        aes(x = thres, y = remain)) +
    geom_col() +
    geom_text(aes(label = remain), color = "green", vjust = 1.2, size = 3) +
    scale_x_continuous(breaks = seq(from = 0, to = 200000, by = 50000), 
                       # remove space between plotted data and xy-axes
                       expand = c(0,0)) +
    ggtitle(file_list1[[i]]) +
    theme(axis.text.x = element_text(size = 20, angle = 90),
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

# RESERVE CODE ===========
# p <- "C:/Users/huyng/OneDrive - Toronto Metropolitan University/Microplastic/Microplastic-Fingerprinting/data/Table of product categorization.xlsx"
# 
# # Create namedf to use as reference for changing name of file name in Table of product categorization.xlsx
# namedf <- rbind(data.frame(File = unique(comp_normalized1$File), NewFile = unique(comp_normalized1$NewFile)), 
#                 data.frame(File = unique(comp_normalized2$File), NewFile = unique(comp_normalized2$NewFile)))
# 
# new_excel <- list()
# for (i in 1:length(excel_sheets(path = p))) {
#   # Table of categorization
#   sampinfo <- readxl::read_excel(path = p, 
#                                  sheet = excel_sheets(path = p)[i])
#   
#   newfile <- c()
#   for (row in 1:nrow(sampinfo)) {
#     newfile <- c(newfile, namedf[match(sampinfo[row,]$File, namedf$File),]$NewFile)
#   }
#   sampinfo$NewFile <- newfile
#   new_excel[[paste0("Grouping", i)]] <- sampinfo
# }
# 
# writexl::write_xlsx(x = new_excel, path = "Table of product categorization_NewFileName.xlsx")

