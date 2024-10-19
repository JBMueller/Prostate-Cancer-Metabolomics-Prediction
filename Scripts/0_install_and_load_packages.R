##########################################
#            Load libraries              #
##########################################

# Packages necessary to run this script
package_list <- c("rstudioapi",   # version: 0.16.0
                  "readxl",       # version: 1.4.3
                  "dplyr",        # version: 1.1.4
                  "dbplyr",       # version: 2.5.0
                  "ggplot2",      # version: 3.5.1
                  "imputeLCMD",   # version: 2.1
                  "remotes",      # version: 2.5.0
                  "rospca",       # version: 1.1.0
                  "isotree",      # version: 0.6.1.1
                  "UpSetR",       # version: 1.4.0
                  "prospectr")    # version: 0.2.7


# Set repositories R can access packages from
setRepositories(ind = c(1, 2)) # Sets CRAN and Bioconductor as repositorieslib

# Load each package in the list
for(package in package_list){
  if (!require(package,character.only = TRUE)) {
    # If the package is not yet installed, install it
    install.packages(package)
  }
  library(package,character.only = TRUE)
} 

# this package needs to be installed from GitHub (uncomment if installation is required)
#remotes::install_github("ricoderks/Rcpm") # note that both git and rtools42 need to be installed on the machine prior
library(Rcpm) # version: 1.0.4
