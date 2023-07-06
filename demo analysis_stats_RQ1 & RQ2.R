# library(R.utils)
`%notin%` <- Negate(`%in%`)

# Research Question 1
stats_rq1 <- shared_comp_sample %>%
  select(File, collapsed_compound, Percent_Area) %>%
  # mutate(File = factor(File, levels = c(unique(File)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = File, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

# Fill in missing value with LOD
for (r in 1:nrow(stats_rq1)) { 
  stats_rq1[r, which(base::is.na(stats_rq1[r,]))] <- runif(length(which(base::is.na(stats_rq1[r,]))),
                                                           min = sort(shared_comp_sample$Percent_Area)[1],
                                                           max = sort(shared_comp_sample$Percent_Area)[2])
}

# Reconstruct data frame
transpose_df <- data.table::transpose(stats_rq1)
rownames(transpose_df) <- colnames(stats_rq1)
colnames(transpose_df) <- rownames(stats_rq1)

transpose_df <- transpose_df %>%
  rownames_to_column(., var = "File") %>%
  mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                               ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                      ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                             ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                    ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                           ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other"))))))) %>%
  relocate(plastic_type, .after = File)


# Examine whether 2 samples have equal variance and normally distributed.-----------------------------------------------
## Normally distributed by histogram, Q-Q plots and Shapiro-Wilk test, if p-value > 0.05 
#-> normally distributed, not for each compound but for the whole
GasData <- as.vector(t(cat_5[,c(1:21)])) # 100149 data points
DieselData <- as.vector(t(cat_5[,c(22:25)])) # 19076 data point

# Histogram
hist(GasData, col='steelblue', main='Gas')
hist(DieselData, col='steelblue', main='Diesel')

# Q-Q plots aka. Normal Probability plots
stats::qqnorm(GasData, main='Gas')
stats::qqline(GasData)

stats::qqnorm(DieselData, main='Diesel')
stats::qqline(DieselData)

# Shapiro-Wilk tests - Error in shapiro.test(GasData) : sample size must be between 3 and 5000 -> do Normal Probability plot instead=
# shapiro.test(GasData)
# shapiro.test(DieselData)

# CONCLUSION: both Gas and Diesel populations are NOT normally distrubuted

# Equality of variance between Gas and Diesel populations by (https://www.r-bloggers.com/2021/06/equality-of-variances-in-r-homogeneity-test-quick-guide/)

# F-test-normally distributed (var.test()); - unusable b/c data non-normally distributed
# Bartlett’s test-normally distributed (bartlett.test()) - unusable b/c data non-normally distributed

# Levene’s test-non-normally distributed data (car::leveneTest()); p-val < 0.05 
# -> there is significant difference between the tested sample variances => Equality of variance is NOT satisfied
library(car)
data <- c(GasData, DieselData)
group <- as.factor(c(rep("Gas", times = length(GasData)), rep("Diesel", times = length(DieselData))))
non_norm_dist_data <- data.frame(data, group)
car::leveneTest(data ~ group, data = non_norm_dist_data)

# Fligner-Killeen test (fligner.test()) - p-value < 0.05 
# -> there is significant difference between tested sample variance => Equality of variance is NOT satisfied
stats::fligner.test(data ~ group, data = non_norm_dist_data)

  
# Wilcoxon Test ========================================================================================
# !! vary sample size gene/compound dataset, three solutions: exact wilcoxon, permutation test, SAM (stanford)

# create empty data.frame
summary_table <- data.frame(collapsed_compound = character(), pair=character(), pvalue=integer())

# for loop through each compound
for (comp in colnames(transpose_df[,3:ncol(transpose_df)])) {
  # looping through the combinations of plastic types 
  for (col in 1:ncol(utils::combn(unique(transpose_df$plastic_type), 2))) {
    # extract the combinations of plastic type pairs
    plastic_type_1 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][1]
    plastic_type_2 <- utils::combn(unique(transpose_df$plastic_type), 2)[,col][2]
    # calculating the p-value between each plastic type pair 
    p.value.res <- wilcox.test(transpose_df[which(transpose_df$plastic_type == plastic_type_1), comp], 
                               transpose_df[which(transpose_df$plastic_type == plastic_type_2), comp])$p.value
    # assigning row information
    summary_table[nrow(summary_table) + 1,] <- c(comp, paste0(plastic_type_1, "_", plastic_type_2), p.value.res)
  }
}

summary_table$adjusted_pvalue_holm <- stats::p.adjust(summary_table$pvalue, method = "holm")
summary_table$adjusted_pvalue_hochberg <- stats::p.adjust(summary_table$pvalue, method = "hochberg")
summary_table$adjusted_pvalue_hommel <- stats::p.adjust(summary_table$pvalue, method = "hommel")
summary_table$adjusted_pvalue_bonferroni <- stats::p.adjust(summary_table$pvalue, method = "bonferroni")
summary_table$adjusted_pvalue_BH <- stats::p.adjust(summary_table$pvalue, method = "BH")
summary_table$adjusted_pvalue_BY <- stats::p.adjust(summary_table$pvalue, method = "BY")
summary_table$adjusted_pvalue_fdr <- stats::p.adjust(summary_table$pvalue, method = "fdr")

alpha <- summary_table %>%
  filter(., adjusted_pvalue_holm < 0.1) %>%
  arrange(adjusted_pvalue_holm)

## Beyond ASTM compounds RQ1 ----------------------
beyondASTMrq1 <- rq1_alpha0.1 %>% filter(., collapsed_compound %notin% c("Compound_481.", "Compound_2269.", "Compound_4588.",
                                                                         "Compound_1557.", "Compound_1893.", "Compound_1981.",
                                                                         "Compound_2005.", "Compound_2065.", "Compound_2166.", 
                                                                         "Compound_2377.", "Compound_2758.", "Compound_3011.",
                                                                         "Compound_3045.", "Compound_3454.", "Compound_4259.",
                                                                         "Compound_4664.", "Compound_4823.", "Compound_4781.",
                                                                         "Compound_2969.", "Compound_4862.", "Compound_7863.", 
                                                                         "Compound_8048.", "Compound_8261.", "Compound_8327.",
                                                                         "Compound_8581."))
beyondASTMrq1_df <- shared_comp_normalized_rt10.1 %>% 
  filter(., collapsed_compound %in% unique(beyondASTMrq1$collapsed_compound)) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = median(RT1), 
            rt2 = median(RT2),
            ion1 = median(Ion1),
            ion2 = median(Ion2)) %>%
  arrange(rt1, rt2)

colnames(beyondASTMrq1_df) <- c("Collapsed_compound", "Retention time 1", "Retention time 2", 
                                "Molecular Ion1", "Molecular Ion2")

writexl::write_xlsx(beyondASTMrq1_df, path = paste0(getwd(), "/beyondASTMrq1_df.xlsx"))


# Research Question 2 ---------------
df_rq2 <- function(data) {
  mydata2 <- data %>%
    filter(., fuel_type %in% "Gas") %>%
    select(sample_name, collapsed_compound, Percent_Area) %>%
    mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
    mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
    # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
    group_by(sample_name, collapsed_compound) %>%
    summarise(across(Percent_Area, mean)) %>%
    pivot_wider(names_from = sample_name, values_from = Percent_Area) %>%
    column_to_rownames(., var = "collapsed_compound")
  
  # transpose the rows and columns
  transpose_mydata2 <- data.table::transpose(mydata2)
  rownames(transpose_mydata2) <- colnames(mydata2)
  colnames(transpose_mydata2) <- rownames(mydata2)
  transpose_mydata2 <- transpose_mydata2 %>%
    rownames_to_column(., var = "sample_name")
  
  transpose_mydata2_new <- transpose_mydata2 %>%
    # # Grouping samples into respective Gas stations
    mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                                ifelse(str_detect(sample_name, "F001"), "Station_1",
                                       ifelse(str_detect(sample_name, "F007"), "Station_7",
                                              ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                     ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                            ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
    relocate(gas_station, .after = sample_name)
  
  # Category 1 (rq2_cat1) : Compound found in only 1 gas station and not in any other
  rq2_cat1 <- data.frame(matrix(nrow = nrow(transpose_mydata2_new)))
  # Category 2 (rq2_cat2): Compound found in >=2 gas stations
  rq2_cat2 <- data.frame(matrix(nrow = nrow(transpose_mydata2_new)))
  
  # if compounds has only 1 record
  rq2_cat1 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 1]
  # if compounds has 2 record
  temp_1 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 2]
  rq2_cat2_col_id <- c()
  rq2_cat1_col_id <- c()
  for (col in 1:ncol(temp_1)) {
    idx1 <- which(!is.na(temp_1[,col]))
    # if 2 records from 2 different gas stations
    if (length(unique(transpose_mydata2_new[idx1]$gas_station)) > 1) {
      rq2_cat2_col_id <- c(rq2_cat2_col_id, col)
    } else { # if 2 records from same gas stations
      rq2_cat1_col_id <- c(rq2_cat1_col_id, col)
    }
  }
  rq2_cat1 <- base::cbind(rq2_cat1, temp_1[, rq2_cat1_col_id])
  rq2_cat2 <- temp_1[, rq2_cat2_col_id]
  # if compounds have 3 records
  temp_2 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 3]
  rq2_cat2_col_id <- c()
  rq2_cat1_col_id <- c()
  for (col in 1:ncol(temp_2)) {
    idx2 <- which(!is.na(temp_2[,col]))
    # if 3 records from different gas stations
    if (length(unique(transpose_mydata2_new[idx2,]$gas_station)) > 1) {
      rq2_cat2_col_id <- c(rq2_cat2_col_id, col)
    } else {  # if 3 records from same gas stations
      rq2_cat1_col_id <- c(rq2_cat1_col_id, col)
    }
  }
  rq2_cat1 <- base::cbind(rq2_cat1, temp_2[, rq2_cat1_col_id])
  rq2_cat2 <- base::cbind(rq2_cat2, temp_2[, rq2_cat2_col_id])
  
  # if compounds have >=4 records -> it definitely appear in >=2 Gas stations (712 compounds)
  temp_3 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) > 3]
  rq2_cat2 <- base::cbind(rq2_cat2, temp_3)
  
  # Insert sample name and gas station codes to Cat1 and Cat2 data frames ====================================
  rq2_cat1 <- base::cbind(rq2_cat1, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  rq2_cat2 <- base::cbind(rq2_cat2, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  # Subset compounds that have >= 2 obs per gas station for Wilcoxon test ------------------------------------------------
  temp_4 <- base::cbind(temp_3, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  col_id <- c()
  
  for (col in 3:ncol(temp_4)) {
    record_count <- c()
    idx <- which(!is.na(temp_4[,col]))
    
    for (station in unique(temp_4[idx, "gas_station"])) {
      count <- sum(temp_4[idx,]$gas_station == station)
      record_count <- c(record_count, count) 
    }
    
    if (sum(record_count >= 2) >= 2) {
      col_id <- c(col_id, col)
    } else {
      next
    }
  }
  
  # Data frame of compounds that have >= 2 obs per gas station
  rq2_cat2_stats <- temp_4[, col_id] # 506 compounds
  
  # Impute NA with LOD  for RQ1 Category 2: compounds in >=2 Gas stations
  for (col in 1:ncol(rq2_cat2_stats)) {
    rq2_cat2_stats[which(is.na(rq2_cat2_stats[,col])), col] <- runif(length(which(is.na(rq2_cat2_stats[,col]))),
                                                                     min = sort(data$Percent_Area)[1],
                                                                     max = sort(data$Percent_Area)[2])
  }


  rq2_cat2_stats_imputed <- base::cbind(rq2_cat2_stats, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  return(rq2_cat2_stats_imputed)
}  

df_stats_rq2_rt10.1 <- df_rq2(shared_comp_normalized_rt10.1)
# df_stats_rq2_rt10.2 <- df_rq2(shared_comp_normalized_rt10.2)
# df_stats_rq2_rt10.3 <- df_rq2(shared_comp_normalized_rt10.3)

# Examine whether 2 gas station groups have equal variance and normally distributed.-----------------------------------------------
## Normally distributed by histogram, Q-Q plots and Shapiro-Wilk test, if p-value > 0.05 
#-> normally distributed, not for each compound but for the whole
gas138_Data <- as.vector(t(df_stats_rq2_rt10.1 %>% 
                             filter(., gas_station %in% c("Station_1", "Station_3", "Station_8")) %>% 
                             select(., 3:ncol(.))))
gas579_Data <- as.vector(t(df_stats_rq2_rt10.1 %>%
                             filter(., gas_station %in% c("Station_5", "Station_7", "Station_9")) %>% 
                             select(., 3:ncol(.))))

# Histogram
hist(gas138_Data, col='steelblue', main='Gas stations 1, 3, 8')
hist(gas579_Data, col='steelblue', main='Gas stations 5, 7, 9')

# Q-Q plots aka. Normal Probability plots
stats::qqnorm(gas138_Data, main='Gas stations 1, 3, 8')
stats::qqline(gas138_Data)

stats::qqnorm(gas579_Data, main='Gas stations 5, 7, 9')
stats::qqline(gas579_Data)

# Shapiro-Wilk tests - Error in shapiro.test(GasData) : sample size must be between 3 and 5000 -> do Normal Probability plot instead=
# shapiro.test(GasData)
# shapiro.test(DieselData)

# CONCLUSION: both Gas station group 1 (1, 3, 8) and Gas station group 2 (5, 7, 9) populations are NOT normally distrubuted

# Equality of variance between Gas and Diesel populations by (https://www.r-bloggers.com/2021/06/equality-of-variances-in-r-homogeneity-test-quick-guide/)

# F-test-normally distributed (var.test()); - unusable b/c data non-normally distributed
# Bartlett’s test-normally distributed (bartlett.test()) - unusable b/c data non-normally distributed

# Levene’s test-non-normally distributed data (car::leveneTest()); p-val > 0.05 
# -> there is no significant difference between the tested sample variances => Equality of variance is satisfied
library(car)
data <- c(gas138_Data, gas579_Data)
group <- as.factor(c(rep("Gas138", times = length(gas138_Data)), rep("Gas579", times = length(gas579_Data))))
non_norm_dist_data <- data.frame(data, group)
car::leveneTest(data ~ group, data = non_norm_dist_data)

# Fligner-Killeen test (fligner.test()) - p-value < 0.05 
# -> there is significant difference between tested sample variance => Equality of variance is NOT satisfied
stats::fligner.test(data ~ group, data = non_norm_dist_data)

# Wilcoxon exact test -----------------------
# Compare between station 1,3,8 and 5,7,9
wilcox_test_rq2 <- function(data) {
  pvalue.w_rq2 <- c()
  i <- 1
  for (col in 3:ncol(data)){
    pvalue.w_rq2[i] <- stats::wilcox.test(as.numeric(data[which(data$gas_station %in% c("Station_5", "Station_7", "Station_9")), col]),
                                   as.numeric(data[which(data$gas_station %in% c("Station_1", "Station_3", "Station_8")), col]))$p.value
    i <- i + 1
  }
  
  rq2_cat2_summary_table <- cbind(colnames(data[,3:ncol(data)]), 
                                  as.data.frame(stats::p.adjust(pvalue.w_rq2,
                                                                       method = "holm"))) %>% # default for p.adjust() method is Holm
    arrange(p.adjust(pvalue.w_rq2))# arrange p value by ascending order
  
  colnames(rq2_cat2_summary_table) <- c("collapsed_compound", "adjusted_pvalue")
  
  rq2_cat2_alpha0.05 <- rq2_cat2_summary_table %>% filter(., adjusted_pvalue < 0.05) # 61 compounds statistically significant
  rq2_cat2_alpha0.1 <- rq2_cat2_summary_table %>% filter(., adjusted_pvalue < 0.1) # 71 compounds statistically significant0
  
  return(list(rq2_cat2_alpha0.05, rq2_cat2_alpha0.1))
}

wilcox_result_rt10.1 <- wilcox_test_rq2(df_stats_rq2_rt10.1)
wilcox_result_rt10.2 <- wilcox_test_rq2(df_stats_rq2_rt10.2)
wilcox_result_rt10.3 <- wilcox_test_rq2(df_stats_rq2_rt10.3)

# sum(is.na(df_stats_rq2_rt10.3[,3:ncol(df_stats_rq2_rt10.3)])) # no NA found so no need to impute data

# Compounds includes extra into rq2_cat2_alpha0.1 from rq2_cat2_alpha0.05
setdiff(rq2_cat2_alpha0.1, rq2_cat2_alpha0.05)


## Beyond ASTM compounds RQ2 ----------------
beyondASTMrq2 <- wilcox_result_rt10.1[[2]] %>% filter(., collapsed_compound %notin% c("Compound_2269.", "Compound_1557.", "Compound_1893.",
                                                                                      "Compound_1981.", "Compound_2005.", "Compound_2065.",
                                                                                      "Compound_2166.", "Compound_2377.", "Compound_2508.",
                                                                                      "Compound_2680.", "Compound_2754.", "Compound_2758.",
                                                                                      "Compound_3011.", "Compound_3045.", "Compound_3129.",
                                                                                      "Compound_3233.", "Compound_3334.", "Compound_3373.",
                                                                                      "Compound_3504.", "Compound_3454.", "Compound_3691.",
                                                                                      "Compound_3803.", "Compound_3863.", "Compound_4259.", 
                                                                                      "Compound_4245.", "Compound_4722.", "Compound_4823.", 
                                                                                      "Compound_2969.", "Compound_4862.", "Compound_8261.", 
                                                                                      "Compound_8327."))

beyondASTMrq2_df <- shared_comp_normalized_rt10.1 %>% 
  filter(., collapsed_compound %in% unique(beyondASTMrq2$collapsed_compound)) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = median(RT1), 
            rt2 = median(RT2),
            ion1 = median(Ion1),
            ion2 = median(Ion2)) %>%
  arrange(rt1, rt2)

colnames(beyondASTMrq2_df) <- c("Collapsed_compound", "Retention time 1", "Retention time 2", 
                                "Molecular Ion1", "Molecular Ion2")

writexl::write_xlsx(beyondASTMrq2_df, path = paste0(getwd(), "/beyondASTMrq2_df.xlsx"))

# RESERVED CODE ---------------------------
# RQ1 Beyond_ASTM data frame that significant Wilcoxon alpha 0.1 
idx <- which(str_detect(ASTM_list$`RQ1: Category 5: Compounds with >=2 Gas record and >=2 Diesel record, PASS WILCOXON TEST (data imputed with LOD), alpha threshold < 0.1 (rt1thres = 0.2)`, 
                        "x"))

name <- unlist(flatten(ASTM_list[idx, 16]))

ASTM_rq1_cat5_wilcoxon_stats_alpha0.1 <- gsub("x_", "", name)

beyond_ASTM_rq1_cat5_wilcox_alpha0.1_names <- unique((alpha0.1 %>%
                                                        filter(., collapsed_compound %notin% ASTM_rq1_cat5_wilcoxon_stats_alpha0.1))$collapsed_compound)

beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% beyond_ASTM_rq1_cat5_wilcox_alpha0.1_names)

# beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df <- beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df[, -c(2,3, 6:8, 16:19)]
write_xlsx(beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df, path = paste0(getwd(), "/RQ1_Beyond_ASTM_passed_Wilcoxon_alpha0.1.xlsx"))

# Which compounds among 118 compounds are ASTM compounds? - 21 out of 80 ASTM found 
# Insert ASTM compounds names to alpha0.1 table
alpha0.1$ASTM <- NA

alpha0.1[alpha0.1$collapsed_compound == "Compound_4720.",]$ASTM <- "Dimethylindane (24.1321)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_4663.",]$ASTM <- "Dimethylindane (24)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_4607.",]$ASTM <- "Dimethylindane OR 4,7-Dimethylindane"
alpha0.1[alpha0.1$collapsed_compound == "Compound_4204.",]$ASTM <- "Pentylbenzene"
alpha0.1[alpha0.1$collapsed_compound == "Compound_3379.",]$ASTM <- "1-Ethyl-2,4-dimethylbenzene (2-1)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_2985.",]$ASTM <- "3- & 4-Propyltoluene (Gang of Four)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_2719.",]$ASTM <- "1,2,3-Trimethylbenzene"
alpha0.1[alpha0.1$collapsed_compound == "Compound_1523.",]$ASTM <- "Isopropylbenzene"
alpha0.1[alpha0.1$collapsed_compound == "Compound_7918.",]$ASTM <- "2,6- & 1,3- & 1,7-Dimethylnaphthalene (Five Fingers)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_8432.",]$ASTM <- "1,4-Dimethylnaphthalene (Five Fingers)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_8122.",]$ASTM <- "1,6-Dimethylnaphthalene (Five Fingers)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_7736.",]$ASTM <- "1-Ethylnaphthalene (Five Fingers)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_7741.",]$ASTM <- "2-Ethylnaphthalene (Five Fingers)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_7237.",]$ASTM <- "Tetradecane (n-alkanes)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_4793.",]$ASTM <- "Naphthalene (PAHs)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_3743.",]$ASTM <- "Isopentylbenzene"
alpha0.1[alpha0.1$collapsed_compound == "Compound_2950.",]$ASTM <- "1,3-Diethylbenzene (Gang of Four)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_2914.",]$ASTM <- "Indane"
alpha0.1[alpha0.1$collapsed_compound == "Compound_2346.",]$ASTM <- "1,2,4-Trimethylbenzene"
# alpha0.1[alpha0.1$collapsed_compound == "Compound_2349.",]$ASTM <- "Decane (n-alkanes)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_2135.",]$ASTM <- "2-Ethyltoluene (Castle Group)"
alpha0.1[alpha0.1$collapsed_compound == "Compound_1862.",]$ASTM <- "Propylbenzene (Castle Group)"

# Insert min-max RT1 to alpha0.1 table
alpha0.1$RT1min <- NA
alpha0.1$RT1max <- NA
alpha0.1$RT2min <- NA
alpha0.1$RT2max <- NA
ASTM_alpha0.1 <- c("Compound_4720.", "Compound_4663.","Compound_4204.", "Compound_3379.", "Compound_2719.", 
                   "Compound_1523.", "Compound_7918.", "Compound_8432.", "Compound_8122.","Compound_7736.",
                   "Compound_7741.", "Compound_7237.","Compound_4793.", "Compound_4607.", "Compound_3743.",
                   "Compound_2985.", "Compound_2950.","Compound_2914.", "Compound_2346.","Compound_2349.",
                   "Compound_2135.","Compound_1862.") 

for (comp in ASTM_alpha0.1) {
  temp <- shared_comp_normalized %>% 
    filter(., collapsed_compound == comp)
  alpha0.1[alpha0.1$collapsed_compound == comp,]$RT1min <- base::sort(temp$RT1)[1]
  alpha0.1[alpha0.1$collapsed_compound == comp,]$RT1max <- base::sort(temp$RT1, decreasing = TRUE)[1]
  alpha0.1[alpha0.1$collapsed_compound == comp,]$RT2min <- base::sort(temp$RT2)[1]
  alpha0.1[alpha0.1$collapsed_compound == comp,]$RT2max <- base::sort(temp$RT2, decreasing = TRUE)[1]
  
}

# rq2_cat2_alpha0.05 -> Which compounds among 61 compounds are ASTM compounds?
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% unique(rq2_cat2_alpha0.05$collapsed_compound)))

# rq2_cat2_alpha0.1 -> Which compounds among 71 compounds are ASTM compounds?
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% unique(rq2_cat2_alpha0.1$collapsed_compound)))

# RQ2 Beyond_ASTM data frame that pass Wilcoxon alpha 0.1 
idx <- which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 Gas record FOR EACH GAS STATION, AFTER Wilcoxon test (data imputed with LOD), alpha threshold < 0.1 (rt1_thres = 0.2)`, 
                        "x"))

name <- unlist(flatten(ASTM_list[idx, 22]))

ASTM_rq2_cat2_wilcoxon_stats_alpha0.1 <- gsub("x_", "", name)

beyond_ASTM_rq2_cat2_wilcox_alpha0.1_names <- unique((rq2_cat2_alpha0.1 %>%
                                                        filter(., collapsed_compound %notin% unique(ASTM_rq2_cat2_wilcoxon_stats_alpha0.1)))$collapsed_compound)

beyond_ASTM_rq2_cat2_wilcox_alpha0.1_df <- shared_comp_normalized %>%
  filter(., collapsed_compound %in% beyond_ASTM_rq2_cat2_wilcox_alpha0.1_names)

beyond_ASTM_rq2_cat2_wilcox_alpha0.1_df <- beyond_ASTM_rq2_cat2_wilcox_alpha0.1_df[, -c(2,3, 6:8, 16:19)]

write_xlsx(beyond_ASTM_rq2_cat2_wilcox_alpha0.1_df, path = paste0(getwd(), "/RQ2_Beyond_ASTM_passed_Wilcoxon_alpha0.1.xlsx"))

# All Beyond ASTM compounds from RQ1 and RQ2
all_beyond_ASTM <- unique(c(setdiff(beyond_ASTM_rq1_cat5_wilcox_alpha0.1_names, 
                                    beyond_ASTM_rq2_cat2_wilcox_alpha0.1_names), 
                            beyond_ASTM_rq2_cat2_wilcox_alpha0.1_names))

all_beyond_ASTM_df <- shared_comp_normalized %>%
  filter(., collapsed_compound %in% all_beyond_ASTM)

all_beyond_ASTM_df <- all_beyond_ASTM_df[, -c(2,3, 6:8, 16:19)]

write_xlsx(all_beyond_ASTM_df, path = paste0(getwd(), "/All_Beyond_ASTM_passed_Wilcoxon_alpha0.1.xlsx"))
