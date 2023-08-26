# This script combines and automates the entire data analysis for LD50/LT50 experiments

# Files must be named in the format Dataxx.csv where xx is a number (01, 02, 03... 10)
# All data used in the analysis should be in one folder with no other extra files
# Run code using ctrl + shift + s (Source) to see uncluttered statistics in the console

# Load or install required packages
if (!require("Cairo")) {
  install.packages("Cairo")
  library(Cairo)
}

source("non_linear_regression.R")

CairoWin()

main_title <- "2h tubeworm oxic pressure-nopressure"

# Directories where to scan for csv data
dir_paths <- list(
  "data/cruise/t_test/tubeworms_mussels/oxic"
)

########### ATTENTION!!!!!!! ##########
# params start values for nls model (must be the same element numbers as dir_paths)
nls_param_list <- list(
  list(100, 30, 4)
)

k <- 0
# Iterate over directories
for (dir_path in dir_paths) {
  k <- k + 1
  chart_subtitle <- chart_subtitle_func(dir_path)
  non_linear_regression_result <- non_linear_regression_dir_func(
    dir_path = dir_path,
    nls_param_list = nls_param_list,
    k = k,
    chart_subtitle = chart_subtitle,
    main_title = main_title
  )
}

# ANOVA model analysis test
# p1 Estimate first model
beta1 <- non_linear_regression_result$nls_coefficients[[1]]
# p1 Std. Error first model
se1 <- non_linear_regression_result$nls_coefficients[[4]]
# p1 Estimate second model
beta2 <- non_linear_regression_result$nls_coefficients[[13]]
# p1 Std. Error second model
se2 <- non_linear_regression_result$nls_coefficients[[16]]
t_statistics <- (beta1 - beta2) / sqrt(se1^2 + se2^2)
p_value <- 2 * (1 - pt(abs(t_statistics), non_linear_regression_result$min_df))
cat("########## ANOVA ANALYSIS ##########\n\n")
cat("P Value: ", p_value, "\n")
