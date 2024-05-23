library(ggplot2)
library(vegan)
library(readxl)
library(tidyverse)
library(dplyr)
library(data.table)
library(writexl)
library(tidyr)
library(grid)
library(gridExtra)

library(stats)
library(FactoMineR)
library(factoextra)
library(compositions)
library(ggforce)
library(latticeExtra)
library(cluster)

polymer <- readxl::read_xlsx(paste0(getwd(), "/data/TMU_Polymer_trainingsetA_Thomas_crookes - after remove columns with unknown label.xlsx"), 
                             .name_repair = "minimal" # No name repair or checks, beyond basic existence,
                             )


group_columns <- function(df) {
  # Get the unique column names
  unique_names <- unique(names(df))
  
  # Initialize an empty list to store grouped column names
  grouped_columns <- list()
  
  # Loop through unique column names
  for (name in unique_names) {
    # Find columns with the same name
    matching_columns <- which(names(df) == name)
    
    # If there are duplicates, group them
    if (length(matching_columns) > 1) {
      grouped_columns[[name]] <- matching_columns
    }
  }
  
  return(grouped_columns)
}

# Create a list of columns with the same name (aka. replicates of the same measurement)
column_groups <- group_columns(polymer)

df_list <- list()
count <- 1
# For each element in the list, make a dataframe that contains wavenumber, sample_name and values
for (i in 1:length(column_groups)) {
  temp <- polymer[, c(1, column_groups[i][[1]])]
  for (j in 2:ncol(temp)) {
    colnames(temp)[j] <- paste0(colnames(temp)[j], "_rep", j-1)
  }
  
  temp <- temp %>% pivot_longer(cols = 2:ncol(.), names_to = "sample", values_to = "values")
  df_list[[count]] <- temp
  count <- count + 1
}

# Combine them all into big df
grand_df <- do.call(rbind, df_list) %>%
  mutate(polymer = ifelse(grepl("ULDPE", sample, ignore.case = TRUE), "ULDPE",
                          ifelse(grepl("LLDPE1", sample, ignore.case = TRUE), "LLDPE1",
                                 ifelse(grepl("LLDPE2", sample, ignore.case = TRUE), "LLDPE2",
                                        ifelse(grepl("LDPE1", sample,, ignore.case = TRUE), "LDPE1",
                                               ifelse(grepl("LDPE2", sample, ignore.case = TRUE), "LDPE2",
                                                      ifelse(grepl("MDPE", sample, ignore.case = TRUE), "MDPE",
                                                             ifelse(grepl("HDPE1", sample, ignore.case = TRUE), "HDPE1",
                                                                    ifelse(grepl("HDPE2", sample,ignore.case = TRUE), "HDPE2", 
                                                                           ifelse(grepl("PP", sample,ignore.case = TRUE), "PP", 
                                                                                  ifelse(grepl("PEST", sample,ignore.case = TRUE), "PEST", 
                                                                                         ifelse(grepl("PET1", sample,ignore.case = TRUE), "PET1", 
                                                                                                ifelse(grepl("PET2", sample,ignore.case = TRUE), "PET2", 
                                                                                                       ifelse(grepl("EVA", sample,ignore.case = TRUE), "EVA",
                                                                                                              ifelse(grepl("ABS", sample,ignore.case = TRUE), "ABS",
                                                                                                                     ifelse(grepl("EPS", sample,ignore.case = TRUE), "EPS",
                                                                                                                            ifelse(grepl("PS", sample,ignore.case = TRUE), "PS",
                                                                                                                                   ifelse(grepl("PA6", sample,ignore.case = TRUE), "PA6",
                                                                                                                                          ifelse(grepl("PA66", sample,ignore.case = TRUE), "PA66",
                                                                                                                                                 ifelse(grepl("PVC1", sample,ignore.case = TRUE), "PVC1",
                                                                                                                                                        ifelse(grepl("PVC2", sample, ignore.case = TRUE), "PVC2",
                                                                                                                                                               ifelse(grepl("CR", sample, ignore.case = TRUE), "CR", "CA")))))))))))))))))))))) %>%
  pivot_wider(names_from = `Wavenumber (cm{¹)`, values_from = values)

# Export grand_df to excel
writexl::write_xlsx(grand_df, path = paste0(getwd(), "/Thomas_Polymer_grand_df.xlsx"))


# PCA
res.pca <- FactoMineR::PCA(
  data %>% select(-c("Subcategory")),
  # impute_list[[4]][[1]],
  # test2[, -1],
  # clr_transformed_pca,
  scale.unit = FALSE,
  graph = FALSE)

# Scree plot
fviz_screeplot(res.pca, ncp=10)

# Biplot
factoextra::fviz_pca_biplot(res.pca,  
                            geom = c("point", "text"),
                            label = "none", 
                            invisible = "var", 
                            repel = TRUE,
                            labelsize = 10, 
                            habillage = factor(data$Subcategory),
                            # addEllipses = TRUE,
                            # ellipse.level=0.95,
                            ggtheme = ggplot2::theme_minimal(base_size = 30),
                            title = "",
                            xlim = c(-0.5, 0.5),
                            ylim = c(-0.25, 0.2)
) + 
  guides(fill=guide_legend(ncol=5)) +
  # if has error "Too few points to calculate an ellipse"
  ggforce::geom_mark_ellipse(aes(fill = Groups,
                                 color = Groups),
                             label.buffer = unit(40, 'mm')) +
  theme(legend.position = 'bottom',
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
  )
                                                                    