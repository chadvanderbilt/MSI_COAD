#!/opt/common/CentOS_7/R/R-4.2.0/bin/Rscript
# List of packages to install and load
packages <- c(
  "varhandle",
  "argparse",
  "dplyr",
  "ggplot2",
  "tibble",
  "tidyr",
  "stringr",
  "stringi",
  "ggpubr",
  "fmsb"
)

# Install packages if not already installed
install_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(install_packages) > 0) {
  install.packages(install_packages, dependencies = TRUE)
}

