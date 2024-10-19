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
#         Data Exploration               #
##########################################
# note that many more characteristics are collected in the meta data, feel free to explore other variables
# suggestions are: Arthritis, Diabetes, IBD, Ulcer, ChronicLungDisease, HeartFailure, Stroke, Hypertension, 
#                 HeartAttack, Angina, Cirrhosis and Depression

# creating separate plotting data in case manipulations are necessary that may intervene with downstream analysis
plotting_data <- meta_data

# changing diagnosis from 0 and 1 to no cancer and cancer
plotting_data$Diagnosis <-ifelse(meta_data$Diagnosis == 0, "no cancer", "cancer")

# changing SmokingCurrently from 0 and 1 to non-smoking and smoking
plotting_data$SmokingCurrently <-ifelse(meta_data$SmokingCurrently == 0, "non-smoker", "smoker")

# changing SmokingCurrently from 0 and 1 to non-smoking and smoking
plotting_data$SmokingEver <-ifelse(meta_data$SmokingEver == 0, "never smoked", "smoked")

# plot the age distribution of both groups
ggplot(plotting_data, aes(x = as.factor(Diagnosis), y = Age, fill = as.factor(Diagnosis))) +
  geom_violin(trim = FALSE) +  # create the violin plot 
  geom_boxplot(width = 0.03, position = position_dodge(0.9), outlier.shape = "*") + # create a boxplot in the middle
  labs(title = "Age Grouped by Diagnosis",
       x = "Diagnosis",
       y = "Age") +
  theme_minimal() +
  theme(legend.position = "none")

# save the plot to the figures folder
ggsave(
  filename = "Figures/violin_plot_age_diagnosis.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

# Plotting previous smoking status between groups
ggplot(plotting_data, aes(x = Diagnosis, fill = SmokingEver)) + 
  geom_bar(position = "fill") +
  labs(title = "Previous Smoking Status by Diagnosis",
       x = "Diagnosis",
       y = "Proportion") +
  scale_fill_manual(values = c("smoked" = "lightcoral", "never smoked" = "lightblue")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# save the plot to the figures folder
ggsave(
  filename = "Figures/Bar_plot_previous_smoking_diagnosis.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)


# Plotting current smoking status between groups
ggplot(plotting_data, aes(x = Diagnosis, fill = SmokingCurrently)) + 
  geom_bar(position = "fill") + 
  labs(title = "Current Smoking Status by Diagnosis",
       x = "Diagnosis",
       y = "Proportion") +
  scale_fill_manual(values = c("smoker" = "lightcoral", "non-smoker" = "lightblue")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# save the plot to the figures folder
ggsave(
  filename = "Figures/Bar_plot_current_smoking_diagnosis.png",  
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)
