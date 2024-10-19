##########################################
#                Setup                   #
##########################################

# setting working directory to current file location
wd <- getActiveDocumentContext()$path
setwd(dirname(dirname(wd)))

# Specify the path to the folder
figures_folder <- paste0(dirname(dirname(wd)), "/Figures")

# Check if the folder exists
if (!dir.exists(figures_folder)) {
  # If the folder does not exist, create it
  dir.create(figures_folder)
  message("Folder created: ", figures_folder)
} else {
  message("Folder already exists: ", figures_folder)
}

##########################################
#            Data Importing              #
##########################################

# importing metabolomics and meta data
metabolomics <- as.data.frame(read_excel("Data/EDRN_MetabolomicsData_2023_08_27.xlsx", col_names = T, na = "", sheet = 2))
meta_data <- as.data.frame(read_excel("Data/EDRN_ClinicalData.xlsx", col_names = T, na = "", sheet = 2))

# finding which rows do not have any metabolomics measurements
unmeasured_sample_id <- meta_data$SAMPLE_ID[!meta_data$SAMPLE_ID %in% metabolomics$CLIENT_SAMPLE_ID]

# removing meta_data samples without metabolomics measurements
meta_data <- meta_data[!meta_data$SAMPLE_ID %in% unmeasured_sample_id,]

# adding patient id to rowname
rownames(metabolomics) <- metabolomics$CLIENT_SAMPLE_ID
metabolomics <- metabolomics %>% select(-"CLIENT_SAMPLE_ID")


##########################################
#            Missing Features            #
##########################################

# number of missing values per feature
na_feature <- colSums(is.na(metabolomics))
na_feature_ratio <- data.frame(ratios = na_feature/nrow(metabolomics))

# remove metabolites with more than 50% missing
metabolomics_filtered <- metabolomics[,!na_feature_ratio > 0.5]

##########################################
#              Transform                 #
##########################################

# check if any 0 is present in the metabolomics data
if (any(metabolomics_filtered == 0, na.rm = T)) {
  print("zeroes need to be resolved before log transformation")
} else {
  metabolomics_log <- log(metabolomics_filtered)
}

set.seed(13)
# quantile regression imputation of left-censored data (QRILC)
sigma_grid <- seq(0.3, 1, 0.1)
# iterate over all sigma values, this may take a while. Around 5 minutes per iteration
for (i in 1:length(sigma_grid)){
  metabolomics_imputed <- impute.QRILC(metabolomics_log, tune.sigma = sigma_grid[i])[[1]] # which sigma to choose?
  metabolomics_normalized <- pqn(metabolomics_imputed, n = "median", QC = NULL)
  # autoscaling the data for PCA
  PCA_data <- apply(metabolomics_normalized, 2, function(col){
    mean_col <- mean(col, na.rm = T)
    sd_col <- sd(col, na.rm = T)
    (col-mean_col)/sqrt(sd_col)
  })
  
  pca_result<- robpca(PCA_data, alpha = 0.75)
  PCA_scores <- as.data.frame(pca_result$scores)
  
  pca_plot <- ggplot(PCA_scores, aes(x = PC1, y=PC2)) + 
    geom_point(aes(color=meta_data$Diagnosis), alpha = 0.7) +
    theme(legend.position = "none")
  
  # Define the filename for the plot
  file_name <- paste0("Figures/PCA_plot_sigma_", sigma_grid[i], ".png")
  
  # Save the plot as a PNG file
  ggsave(filename = file_name, plot = pca_plot, width = 8, height = 6)
}

# all of these PCA score plots do not demonstrate markedly different distributions. Therefore it was concluded
# that the value of sigma will not affect downstream 