library(e1071)
library(ROCR)
library(MLmetrics)


# e1071 package: ===========================

# PCA-based feature data
df_pca <- input_df(atdgcms)

# PCAtools::pca requires mat input (columns as "sample name", rows as "collapsed_compound")
p <- PCAtools::pca(mat = df_pca[[1]], 
                   metadata = df_pca[[2]], 
                   # center = FALSE,
                   scale = FALSE 
)

PCAtools_atdgcms <- p$rotated %>%
  tibble::rownames_to_column(., var = "File") %>%
  dplyr:: mutate(product_cat = ifelse(str_detect(File, "Balloons"), "Toys", 
                                      ifelse(str_detect(File, "FPW"), "Food contact materials",
                                             ifelse(str_detect(File, "Pbal"), "Toys",
                                                    ifelse(str_detect(File, "MPW"), "Mixed_Plastic_Waste", 
                                                           ifelse(str_detect(File, "PBBC"), "Food contact materials",
                                                                  ifelse(str_detect(File, "Pbag"),"Food contact materials",
                                                                         ifelse(str_detect(File, "PDS"),"Food contact materials", 
                                                                                ifelse(str_detect(File, "Pcut"), "Food contact materials",
                                                                                       ifelse(str_detect(File, "PC"), "Food contact materials",
                                                                                              ifelse(str_detect(File, "Cigs"), "Cigarettes", 
                                                                                                     ifelse(str_detect(File, "Cmat"), "Construction materials",
                                                                                                            ifelse(str_detect(File, "Mask"), "Clothes", "Misc"))))))))))))) %>%
  dplyr::mutate(product_cat = factor(product_cat, levels = unique(product_cat))) %>%
  dplyr::relocate(product_cat, .before = 1) %>%
  tibble::column_to_rownames(., var = "File")


# PArtitioning train&test sets / training / predict on test set
e1071.SVM.result <- function(dat, split.ratio){
  set.seed(1234)
  plastic_idx <- caret::createDataPartition(dat$product_cat, p = split.ratio, list = F)
  plastic_trn <- dat[plastic_idx, ]
  plastic_tst <- dat[-plastic_idx, ]  
  
  # Create 1 million evenly spaced values on a log scale between 10^-6 and 10^2
  lseq <- function(from = 0.000001, to = 100, length.out=1000000) {
    # logarithmic spaced sequence
    exp(seq(log(from), log(to), length.out = length.out))
  }
  
  # Model tuning
  tune.out <- e1071::tune(e1071::svm, product_cat ~ ., 
                          data = plastic_trn,
                          kernel = "radial",
                          scale = TRUE,
                          # ranges = list(cost = lseq()),
                          decision.values = TRUE, 
                          probability = TRUE)
  bestmod <- tune.out$best.model
  
  # Prediction results
  pred_prob <- predict(bestmod, newdata = plastic_tst,
                       decision.values = TRUE, probability = TRUE)
  pred_res <- attr(pred_prob, "probabilities")
  
  vip_all <- c()
  for (pt in unique(plastic_trn$product_cat)) {
    set.seed(282)  # for reproducibility

    # Add prediction wrapper function to return the predicted class probabilities for the reference class of interest.
    prob <- function(object, newdata) {
      predict(object, newdata = newdata, type = "prob")
    }

    # Variable importance plot
    vip_res <- vip::vip(bestmod, method = "permute", nsim = 10, train = plastic_trn,
                        target = "product_cat", metric = "auc", reference_class = pt,
                        pred_wrapper = prob)
    
    # Get all the top variable for each product_cat as reference_class
    vip_all <- c(vip_all, vip_res$data[1:3,]$Variable)
  }
  
  return(list(pred_res, unique(vip_all)))
}

PCAtools_atdgcms.SVMresult <- e1071.SVM.result(PCAtools_atdgcms, split.ratio = 0.6)


## Original data
df <- atdgcms %>%
  dplyr::select(File, product_cat, collapsed_compound, Values) %>%
  mutate(product_cat = factor(product_cat, levels = unique(product_cat))) %>%
  # since we have multiple different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(File, product_cat, collapsed_compound) %>%
  summarise(across(Values, mean)) %>%
  pivot_wider(names_from = collapsed_compound, values_from = Values) %>%
  column_to_rownames(., var = "File")

for (r in 1:nrow(df)) {
  df[r, which(base::is.na(df[r,]))] <- runif(length(which(base::is.na(df[r,]))),
                                             min = sort(atdgcms$Values)[1],
                                             max = sort(atdgcms$Values)[2])
}

atdgcms.SVMresult <- e1071.SVM.result(df, split.ratio = 0.6)


plot.dat <- function(dat) {
  newdat <- as.data.frame(dat) %>%
    tibble::rownames_to_column(., var = "File") %>%
    tidyr::pivot_longer(., cols = 2:ncol(.), names_to = "product_cat", values_to = "prob")
  return(newdat)
}   

PCAtoolsdat <- plot.dat(PCAtools_atdgcms.SVMresult[[1]])
mydat <- plot.dat(atdgcms.SVMresult[[1]])

# Make bar graph of prediction probability
ggplot(data = mydat) + 
  geom_col(aes(x = File, y = prob, fill = product_cat), 
           position = "dodge" # separating stacking prob cols
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(n.breaks = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))


