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

# checking proportions of patients and control after removal (313 control 267 patients in: https://www.nature.com/articles/s41597-023-02750-7)
table(meta_data$Diagnosis)

# adding patient id to rownames
rownames(metabolomics) <- metabolomics$CLIENT_SAMPLE_ID

# removing patient IDs from the metabolomics measurements
metabolomics <- metabolomics %>% select(-"CLIENT_SAMPLE_ID")

##########################################
#               Duplicates               #
##########################################

# check if there are any duplicate patient IDs
any(duplicated(rownames(metabolomics)))
any(duplicated(meta_data$SAMPLE_ID))

# check if any rows are exact duplicates of each other
any(duplicated(metabolomics))
any(duplicated(meta_data[, 2:420]))

##########################################
#           Inconsistencies              #
##########################################

# check each value in the dataframe if commas are present
commas_present <- any(as.data.frame(lapply(metabolomics, function(x) grepl(",", x))))
if (commas_present){
  # substitute a decimal point for any comma
  metabolomics_dec <- apply(metabolomics, 2, function(x) gsub(",", ".", x))
  # change values to numeric
  metabolomics_dec <- as.data.frame(apply(metabolomics_dec, 2, function(x) as.numeric(x)))
  # add patient ID as rownames
  rownames(metabolomics_dec) <- rownames(metabolomics) # the apply function should not change the order of samples
  # replace metabolomics object with corrected dataframe
  metabolomics <- metabolomics_dec
  # remove corrected dataframe to free up memory
  rm(metabolomics_dec)
}

##########################################
#            Missing Features            #
##########################################

# number of missing values per feature
na_feature <- colSums(is.na(metabolomics))
na_feature_ratio <- data.frame(ratios = na_feature/nrow(metabolomics))

# plot missingness ratio vs frequency
ggplot(na_feature_ratio, aes(ratios)) +
  stat_ecdf(geom = "step", color = "blue", linewidth = 1) +
  labs(title = "Cumulative frequency of missingness", x = "Missingness ratio", y = "Cumulative Frequency") +
  theme_minimal()

# save the plot to the figures folder
ggsave(
  filename = "Figures/Missingness_per_variable.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

# the plot shows that a good cutoff value can be picked at around 0.5 missingness, as it retains more than 75% of features

# investigating missingness per group #

# split metabolomics measurements into groups
if (all(meta_data$SAMPLE_ID == rownames(metabolomics))){ # ensure that the order of samples is the same across measurements and metadata
  cases <- metabolomics[meta_data$Diagnosis == 1, ]
  controls <- metabolomics[meta_data$Diagnosis == 0, ]
} else {
  stop("samples are not in the same order")
}

# ratio of missingness per feature per group
na_per_group <- bind_rows(
  data.frame(ratios = colSums(is.na(cases)) / nrow(cases), group = "Cases"),
  data.frame(ratios = colSums(is.na(controls)) / nrow(controls), group = "Controls")
)

# plotting missingness vs frequency between groups
ggplot(na_per_group, aes(x = ratios, color = group)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(title = "Cumulative Frequency Plot", x = "Ratio", y = "Cumulative Frequency") +
  theme_minimal()

# save the plot to the figures folder
ggsave(
  filename = "Figures/Missingness_between_groups.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

# calculating difference in missingness ratio per feature
diff <- abs(na_per_group[na_per_group$group == "Cases", "ratios"] - na_per_group[na_per_group$group == "Controls", "ratios"])

# plotting missingness per metabolite between groups
plot(diff) # the maximum difference in missingness is less than 12%, therefore deemed not very different between groups

# remove metabolites with more than 50% missing
metabolomics_filtered <- metabolomics[ ,!na_feature_ratio > 0.5]

# investigating missingness per sample #

# calculate sum of missing values per row/sample
na_samples <- rowSums(is.na(metabolomics_filtered))
# calculating ratio of missingness per sample
na_samples_ratio <- data.frame(ratios = na_samples/ncol(metabolomics_filtered))

# plotting sample missingness ratio and cumulative frequency
ggplot(na_samples_ratio, aes(ratios)) +
  stat_ecdf(geom = "step", color = "blue", linewidth = 1) +
  labs(title = "Cumulative Frequency Plot", x = "Ratio", y = "Cumulative Frequency") +
  theme_minimal()

# save the plot to the figures folder
ggsave(
  filename = "Figures/Missingness_per_sample.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)
# the plot shows no extreme values of missingness, maximally 0.25 in very few samples. Therefore no samples are removed

##########################################
#            HM imputation               #
##########################################

# find the minimum value of each metabolite
half_min <- apply(metabolomics_filtered, 2, function(x) min(x, na.rm = T))

# compute half of the minimum value
half_min <- half_min/2

# create a new dataframe to store imputed values
metabolomics_halfmin <- metabolomics_filtered

# replace each NA with half of the minimal value of that metabolite
for (i in 1:ncol(metabolomics_halfmin)) {
  metabolomics_halfmin[is.na(metabolomics_halfmin[, i]), i] <- half_min[i]
}

##########################################
#              Transform                 #
##########################################

# check if any 0 is present in the metabolomics data
if (any(metabolomics_halfmin == 0, na.rm = T)) {
  print("zeroes need to be resolved before log transformation")
} else { # if not, log-transform the dataframe
  metabolomics_log <- log(metabolomics_halfmin)
}

##########################################
#            Train-test split            #
##########################################

# split the data into both classes of Diagnosis
Diagnosis_class_1 <- metabolomics_log[meta_data$Diagnosis == 1, ]
Diagnosis_class_0 <- metabolomics_log[meta_data$Diagnosis == 0, ]

# create a 70/30 train/test split using the Kennard-Stone algorithm
index_list_Diagnosis_1 <- kenStone(Diagnosis_class_1, metric = "mahal", k = nrow(Diagnosis_class_1) *0.7, pc = 0.8, .center = T, .scale = T)
index_list_Diagnosis_0 <- kenStone(Diagnosis_class_0, metric = "mahal", k = nrow(Diagnosis_class_0) *0.7, pc = 0.8, .center = T, .scale = T)

# use indices to subset data to train and test sets for each class
train_Diagnosis_1 <- Diagnosis_class_1[index_list_Diagnosis_1$model, ]
test_Diagnosis_1 <- Diagnosis_class_1[index_list_Diagnosis_1$test, ]

train_Diagnosis_0 <- Diagnosis_class_0[index_list_Diagnosis_0$model, ]
test_Diagnosis_0 <- Diagnosis_class_0[index_list_Diagnosis_0$test, ]

# Combine train and test sets from both classes
metabolomics_train <- rbind(train_Diagnosis_1, train_Diagnosis_0)
metabolomics_test <- rbind(test_Diagnosis_1, test_Diagnosis_0)

# convert train and test sets back to not imputed data (because another imputation method will be used after)
metabolomics_train <- log(metabolomics_filtered[which(rownames(metabolomics_filtered) %in% rownames(metabolomics_train)), ])
metabolomics_test <- log(metabolomics_filtered[which(rownames(metabolomics_filtered) %in% rownames(metabolomics_test)), ])

# extract the appropriate meta data for both train and test sets
meta_data_train <- meta_data[which(meta_data$SAMPLE_ID %in% rownames(metabolomics_train)), ]
meta_data_test <- meta_data[which(meta_data$SAMPLE_ID %in% rownames(metabolomics_test)), ]

# ensure that original sample ID and class label combinations are preserved after splitting and re-ordering (very important!!)
check_discpancies <- function(meta_data_split, meta_data_original) {
  # extract sample ID and diagnosis class
  selected_split <- meta_data_split[, c("SAMPLE_ID", "Diagnosis")]
  selected_original <- meta_data_original[, c("SAMPLE_ID", "Diagnosis")]
  # merge the meta_data_train with the original meta_data based on SAMPLE_ID to compare
  merged_train_data <- merge(selected_split, selected_original, by = "SAMPLE_ID", suffixes = c("_train", "_original")) # merging assures correct matches of sample ID regardless of order
  
  # Check if there are any discrepancies between the Diagnosis in the meta data before and after the split
  discrepancies <- merged_train_data[merged_train_data$Diagnosis_train != merged_train_data$Diagnosis_original, ]
  if (nrow(discrepancies) > 0){
    print("STOP! there are discrepancies between the meta_data before and after the split")
  } else {
    print("no discrepancies detected")
    # removing objects created by the function to free up memory
    rm(selected_split, selected_original, merged_train_data, discrepancies)
  }
}

check_discpancies(meta_data_train, meta_data)
check_discpancies(meta_data_test, meta_data)

##########################################
#           QRILC imputation             # QRILC will be used due to literature recommendations, (see: R script "comparing_QRILC_and_HM_imputation")
##########################################

# set random seed for reproducible results 
set.seed(13)

# re-define the QRILC function to be compatible with train-test data
QRILC_split <- function (train_data, test_data, tune.sigma = 1) {
  dataSet.mvs <- train_data
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2]
  nFeatures.test = dim(test_data)[1]
  dataSet.imputed = dataSet.mvs
  test_data.imputed <- test_data # define test data to be imputed
  QR.obj = list()
  for (i in 1:nSamples) {
    curr.sample = dataSet.mvs[, i]
    pNAs = length(which(is.na(curr.sample)))/length(curr.sample)
    upper.q = 0.99
    q.normal = qnorm(seq((pNAs + 0.001), (upper.q + 0.001), 
                         (upper.q - pNAs)/(upper.q * 100)), mean = 0, sd = 1)
    q.curr.sample = quantile(curr.sample, probs = seq(0.001, 
                                                      (upper.q + 0.001), 0.01), na.rm = T)
    temp.QR = lm(q.curr.sample ~ q.normal)
    QR.obj[[i]] = temp.QR
    mean.CDD = temp.QR$coefficients[1]
    sd.CDD = as.numeric(temp.QR$coefficients[2])
    data.to.imp = rtmvnorm(n = nFeatures, mean = mean.CDD, 
                           sigma = sd.CDD * tune.sigma, upper = qnorm((pNAs + 
                                                                         0.001), mean = mean.CDD, sd = sd.CDD), algorithm = c("gibbs"))
    
    curr.sample.imputed = curr.sample
    curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(curr.sample))]
    dataSet.imputed[, i] = curr.sample.imputed
    
    # repeat sampling the distribution for the test sets based on the distribution of train set
    test.data.to.imp = rtmvnorm(n = nFeatures.test, mean = mean.CDD, 
                                sigma = sd.CDD * tune.sigma, upper = qnorm((pNAs + 
                                                                              0.001), mean = mean.CDD, sd = sd.CDD), algorithm = c("gibbs"))
    
    test.sample.imputed = test_data[, i]
    test.sample.imputed[which(is.na(test_data[, i]))] = test.data.to.imp[which(is.na(test_data[, i]))]
    test_data.imputed[, i] = test.sample.imputed
  }
  results = list(dataSet.imputed, test_data.imputed, QR.obj)
  return(results)
}

# perform quantile regression imputation of left-censored data (QRILC) imputation on the train and test data
QRILC_results <- QRILC_split(metabolomics_train, metabolomics_test, tune.sigma = 1) # value of sigma does not affect distribution markedly (see: R script "comparing_sigma_values")
metabolomics_imputed_train <- QRILC_results[[1]]
metabolomics_imputed_test <- QRILC_results[[2]]

##########################################
#           PQN normalization            #
##########################################

# re-define the pqn function to be compatible with train-test data
pqn_train <- function (X, n = "median", QC = NULL) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  if (!is.null(QC)) {
    if (length(QC) == 1) {
      mX <- as.numeric(X[QC, ])
    }
    else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  }
  else {
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ]/median(as.numeric(X[a, 
    ]/mX)))
  }
  return(list(X.norm, mX)) # return the normalized train data and median values to be used later for test-data normalization
}

# define a function based on the pqn function to impute the test data using medians of the train data
pqn_test <- function (X, mX) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ]/median(as.numeric(X[a, 
    ]/mX)))
  }
  return(X.norm)
}

# Probabilistic quotient normalization (this may take a short while)
pqn_res <- pqn_train(metabolomics_imputed_train, n = "median", QC = NULL)
# take the first returned element to be the normalized train data
metabolomics_normalized_train <- pqn_res[[1]]
# take the second returned element to normalize the test data 
metabolomics_normalized_test <- pqn_test(metabolomics_imputed_test, pqn_res[[2]])

##########################################
#             Robust PCA                 #
##########################################

# pareto-scaling the data for PCA
PCA_data <- apply(metabolomics_normalized_train, 2, function(col){
  mean_col <- mean(col, na.rm = T)
  sd_col <- sd(col, na.rm = T)
  (col-mean_col)/sqrt(sd_col)
})

# perform robust PCA with 75% of the data at each iteration
pca_result<- robpca(PCA_data, alpha = 0.75)
PCA_scores <- as.data.frame(pca_result$scores)

# Calculate the proportion of variance explained for each principal component
variance_explained <- pca_result$eigenvalues / sum(pca_result$eigenvalues) * 100

# Create dynamic labels with variance explained for PC1 and PC2 (change this to the PCs of interest)
x_label <- paste0("PC1 (", sprintf("%.2f", variance_explained[1]), "%)")
y_label <- paste0("PC2 (", sprintf("%.2f", variance_explained[2]), "%)")

# Plot the PCA scores with axis labels
ggplot(PCA_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = as.factor(meta_data_train$Diagnosis)), alpha = 0.7) +
  coord_fixed(1) +
  xlab(x_label) + 
  ylab(y_label) 

# save the plot to the figures folder
ggsave(
  filename = "Figures/robPCA_PC1_PC2_diagnosis.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

# projecting the test data onto the PCA #
# the PCA scores plot shows that the test data is similarly distributed as the training data, ensuring realistic validation outcomes

# pareto-scaling the test data for PCA
PCA_test <- apply(metabolomics_normalized_test, 2, function(col){
  mean_col <- mean(col, na.rm = T)
  sd_col <- sd(col, na.rm = T)
  (col-mean_col)/sqrt(sd_col)
})

# calculating scores by multiplying with previously obtained loadings
PCA_test <- as.data.frame(as.matrix(PCA_test) %*% as.matrix(pca_result$loadings))
PCA_test$Type <- "Test"
PCA_scores$Type <- "Train"
# combining previous PCA scores and test_data scores to one object for plotting
combined_PCA <- rbind(PCA_scores, PCA_test)
combined_PCA$Diagnosis <- c(meta_data_train$Diagnosis, meta_data_test$Diagnosis)
# changing the labels for cancer and non-cancer for a nicer legend in the PCA
combined_PCA[combined_PCA$Diagnosis == 0, "Diagnosis"] <- "Control"
combined_PCA[combined_PCA$Diagnosis == 1, "Diagnosis"] <- "Case"

# plotting the PCA scores plot
ggplot(combined_PCA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Diagnosis, shape = Type), alpha = 0.7) +
  coord_fixed(1) +
  xlab(x_label) + 
  ylab(y_label)

# save the plot to the figures folder
ggsave(
  filename = "Figures/projected_test_robPCA_PC1_PC2_diagnosis.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

##########################################
#          Isolation Forest              # isolation forest was used due to its robust performance and literature recommendations see: https://www.sciencedirect.com/science/article/pii/S0031320317303916
##########################################

# define a range of trees
n_trees <- 2^(7:12)
# initialize a list to store outliers per iteration
outlier_samples <- list()

# iterate over all values of n_trees to build isolation forests of varying sizes
for (i in 1:length(n_trees)){
  # build an extended isolation forest (extended due to ndim=2)
  ext_iforest <- isolation.forest(
    metabolomics_normalized_train, ndim=2, ntrees=n_trees[i],
    missing_action="fail", sample_size = nrow(metabolomics_normalized_train),
    sample_with_replacement = T, 
  )
  # extract the outlier score for each sample
  pred <- predict.isolation_forest(ext_iforest, metabolomics_normalized_train)
  # store the samples with an outlier score above a threshold 0.5 was chosen as this yielded a reasonable number of outliers
  outlier_samples[[i]] <- rownames(metabolomics_normalized_train[which(pred > 0.5),])
}

# name the list components for plotting later
names(outlier_samples) <- paste0("ntrees_", 2^(7:12))

# create a binary dataframe indicating presence/absence of each sample in the outlier sets
binary_df <- data.frame(
  row.names = unique(unlist(outlier_samples)), # set rownames to all unique samples across all isolation forests
  # apply a function over each outlier set created by each isolation forest
  sapply(outlier_samples, function(x) {
    # extract all unique outliers of all outlier sets
    unique_samples <- unique(unlist(outlier_samples))
    # find which outliers out of all outliers are in the current outlier set
    unique_samples %in% x # this creates a logical vector for each sample indicating whether it is present in this set
  })
)

# convert logical values to numeric to use in an upset plot
binary_df <- data.frame(apply(binary_df, MARGIN = 2, FUN = as.numeric))

# create the upset plot
upset_plot <- upset(
  binary_df, 
  sets = names(outlier_samples), 
  mainbar.y.label = "Intersection size", 
  sets.x.label = "Number of outliers", 
  keep.order = T, 
  nsets = 6
)


# based on the previous exploration ntrees does not matter much and 128 were chosen for computational efficiency (which is in line with the authors of isolation forest)
ext_iforest <- isolation.forest(
  metabolomics_normalized_train, ndim=2, ntrees=128,
  missing_action="fail", sample_size = nrow(metabolomics_normalized_train), # sample size does not matter according to the authors, taking all samples seems therefore appropriate and is the default
  sample_with_replacement = T, 
)

pred <- predict.isolation_forest(ext_iforest, metabolomics_normalized_train)

# visualize outlier scores on previous PCA (look similar on PC1 and PC2 but PC3 shows the difference)
ggplot(PCA_scores, aes(x = PC1, y=PC3)) + 
  geom_point(aes(color=as.factor(pred > 0.5)), alpha = 0.7) + # 0.5 as a cutoff for outlier score
  scale_color_discrete(name = "anomaly score above 0.5") # Change the legend title  

# save the plot to the figures folder
ggsave(
  filename = "Figures/anomaly_score_PCA.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

# remove the samples which have an outlier score above 0.5
metabolomics_outliers_removed <- metabolomics_normalized_train[-which(pred > 0.5), ]
meta_data_outliers_removed <- meta_data_train[-which(pred > 0.5), ]

# checking proportion patients and controls after outlier removal
table(meta_data_outliers_removed$Diagnosis)

# scaling data for PCA after outlier removal (this may impact the mean somewhat)
PCA_data <- apply(metabolomics_outliers_removed, 2, function(col){
  mean_col <- mean(col, na.rm = T)
  sd_col <- sd(col, na.rm = T)
  (col-mean_col)/sqrt(sd_col)
})

# re-computing PCA for data without outliers
pca_result<- robpca(PCA_data, alpha = 1)
PCA_scores <- as.data.frame(pca_result$scores)

# Calculate the proportion of variance explained for each principal component
variance_explained <- pca_result$eigenvalues / sum(pca_result$eigenvalues) * 100

# Create dynamic labels with variance explained for PC1 and PC2 (change this to the PCs of interest)
x_label <- paste0("PC1 (", sprintf("%.2f", variance_explained[1]), "%)")
y_label <- paste0("PC2 (", sprintf("%.2f", variance_explained[2]), "%)")

# plot the PCA scores
ggplot(PCA_scores, aes(x = PC1, y=PC2)) + 
  geom_point(aes(color=as.factor(meta_data_outliers_removed$Diagnosis)), alpha = 0.7) +
  coord_fixed(1) +
  xlab(x_label) + 
  ylab(y_label) +
  scale_color_discrete(name = "Diagnosis")

# save the plot to the figures folder
ggsave(
  filename = "Figures/PCA_outliers_removed.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

##########################################
#           Writing Data                 #
##########################################

# make sure the sample ID in the metabolomics and meta data are in the same order
all(rownames(metabolomics_outliers_removed) == meta_data_outliers_removed$SAMPLE_ID)
# ensure that class labels are preserved after removing outliers
check_discpancies(meta_data_outliers_removed, meta_data)

# saving train metabolomics and metadata for downstream analysis in python
write.csv(metabolomics_outliers_removed, "Data/train_preprocessed.csv", col.names = T, row.names = T)
write.csv(meta_data_outliers_removed, "Data/train_metadata.csv", col.names = T, row.names = T)

# make sure the sample ID in the metabolomics and meta data are in the same order
all(rownames(metabolomics_normalized_test) == meta_data_test$SAMPLE_ID)
# ensure that class labels are preserved in the test set after normalizing
check_discpancies(meta_data_test, meta_data)

# saving test metabolomics and metadata for downstream analysis in python
write.csv(metabolomics_normalized_test, "Data/test_preprocessed.csv", col.names = T, row.names = T)
write.csv(meta_data_test, "Data/test_metadata.csv", col.names = T, row.names = T)