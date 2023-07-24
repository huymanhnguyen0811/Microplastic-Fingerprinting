my.data <- shared_comp_plastic_type  %>%
  dplyr::select(File, collapsed_compound, Percent_Area) %>%
  filter(., !str_detect(File, "PC_Sample")) %>%
  group_by(File, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>% # , we summarize these values by taking the mean of them
  pivot_wider(names_from = collapsed_compound, values_from = Percent_Area) %>%
  mutate(plastic_type = ifelse(str_detect(File, "Balloons"), "Balloons", 
                               ifelse(str_detect(File, "FPW_"), "Food_Packaging_Waste",
                                      ifelse(str_detect(File, "MPW_"), "Mixed_Plastic_Waste", 
                                             ifelse(str_detect(File, "PBBC_"), "Plastic_Bottles_and_Bottle_Caps",
                                                    # ifelse(str_detect(File, "PC_Sample"),"Plastic_Cups",
                                                           ifelse(str_detect(File, "PDS_Sample"),"Plastic_Drinking_Straws", "Other")))))) %>%
  mutate(plastic_type = factor(plastic_type, levels = unique(plastic_type))) %>%
  relocate(plastic_type, .before = 1) %>%
  column_to_rownames(., var = "File")

# Fill NA with LOD
for (c in 2:ncol(my.data)) { 
  my.data[which(base::is.na(my.data[,c])), c] <- runif(length(which(base::is.na(my.data[,c]))),
                                                       min = sort(shared_comp_plastic_type$Percent_Area)[1],
                                                       max = sort(shared_comp_plastic_type$Percent_Area)[2])
}
