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
#            Missing Features            #
##########################################

# number of missing values per feature
na_feature <- colSums(is.na(metabolomics))
na_feature_ratio <- data.frame(ratios = na_feature/nrow(metabolomics))

# remove metabolites with more than 50% missing
metabolomics_filtered <- metabolomics[ ,!na_feature_ratio > 0.5]

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

# extract the appropriate meta data for both train and test sets
meta_data_train <- meta_data[which(meta_data$SAMPLE_ID %in% rownames(metabolomics_train)), ]
meta_data_test <- meta_data[which(meta_data$SAMPLE_ID %in% rownames(metabolomics_test)), ]

# convert train and test sets back to not imputed data (because another imputation method will be used after)
metabolomics_train_QRILC <- log(metabolomics_filtered[which(rownames(metabolomics_filtered) %in% rownames(metabolomics_train)), ])
metabolomics_test_QRILC <- log(metabolomics_filtered[which(rownames(metabolomics_filtered) %in% rownames(metabolomics_test)), ])

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
#           QRILC imputation             #
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
QRILC_results <- QRILC_split(metabolomics_train_QRILC, metabolomics_test_QRILC, tune.sigma = 1) # value of sigma does not affect distribution markedly (see: R script "finding sigma")
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

# Probabilistic quotient normalization for both imputation methods
pqn_res_QRILC <- pqn_train(metabolomics_imputed_train, n = "median", QC = NULL)
pqn_res_HM <- pqn_train(metabolomics_train, n = "median", QC = NULL)

# extract the results of PQN for QRILC imputation
QRILC_train <- pqn_res_QRILC[[1]]
QRILC_test <- pqn_test(metabolomics_imputed_test, pqn_res_QRILC[[2]])

# extract the results of PQN for HM imputation
HM_train <- pqn_res_HM[[1]]
HM_test <- pqn_test(metabolomics_imputed_test, pqn_res_HM[[2]])

##########################################
#         Robust PCA QRILC               #
##########################################

# pareto-scaling the data for PCA
PCA_data <- apply(QRILC_train, 2, function(col){
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

# projecting the test data onto the PCA #

# pareto-scaling the test data for PCA
PCA_test <- apply(QRILC_test, 2, function(col){
  mean_col <- mean(col, na.rm = T)
  sd_col <- sd(col, na.rm = T)
  (col-mean_col)/sqrt(sd_col)
})

# calculating scores by multiplying with previously obtained loadings
PCA_test <- as.data.frame(as.matrix(PCA_test) %*% as.matrix(pca_result$loadings))
PCA_test$Type <- "Test"
PCA_scores$Type <- "Train"
combined_PCA <- rbind(PCA_scores, PCA_test)
combined_PCA$Diagnosis <- c(meta_data_train$Diagnosis, meta_data_test$Diagnosis)
combined_PCA[combined_PCA$Diagnosis == 0, "Diagnosis"] <- "Control"
combined_PCA[combined_PCA$Diagnosis == 1, "Diagnosis"] <- "Case"

ggplot(combined_PCA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Diagnosis, shape = Type), alpha = 0.7) +
  coord_fixed(1) +
  xlab(x_label) + 
  ylab(y_label) +
  ggtitle("PCA on QRILC imputation")

##########################################
#          Robust PCA HM                 #
##########################################

# pareto-scaling the data for PCA
PCA_data <- apply(HM_train, 2, function(col){
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

# projecting the test data onto the PCA #

# pareto-scaling the test data for PCA
PCA_test <- apply(HM_test, 2, function(col){
  mean_col <- mean(col, na.rm = T)
  sd_col <- sd(col, na.rm = T)
  (col-mean_col)/sqrt(sd_col)
})

# calculating scores by multiplying with previously obtained loadings
PCA_test <- as.data.frame(as.matrix(PCA_test) %*% as.matrix(pca_result$loadings))
PCA_test$Type <- "Test"
PCA_scores$Type <- "Train"
combined_PCA <- rbind(PCA_scores, PCA_test)
combined_PCA$Diagnosis <- c(meta_data_train$Diagnosis, meta_data_test$Diagnosis)
combined_PCA[combined_PCA$Diagnosis == 0, "Diagnosis"] <- "Control"
combined_PCA[combined_PCA$Diagnosis == 1, "Diagnosis"] <- "Case"

ggplot(combined_PCA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Diagnosis, shape = Type), alpha = 0.7) +
  coord_fixed(1) +
  xlab(x_label) + 
  ylab(y_label) +
  ggtitle("PCA on HM imputation")

# There is no clear difference visible between both imputation methods on the PCA. However due to its recommendations in literature (https://www.nature.com/articles/s41598-017-19120-0)
# QRILC will be used for the final imputation as it preserves the original distribution better, which may be difficult to see on such
# high dimensional data.