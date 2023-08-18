# This script combines and automates the entire data analysis for LD50/LT50 experiments

# Files must be named in the format Dataxx.csv where xx is a number (01, 02, 03... 10)
# All data used in the analysis should be in one folder with no other extra files
# Run code using ctrl + shift + s (Source) to see uncluttered statistics in the console

# Load or install required packages
if (!require("MASS")) {
  install.packages("MASS")
  library(MASS)
}

if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

source("chart_subtitle.R")
source("non_linear_regression.R")
source("linear_regression.R")

# Directories where to scan for csv data
dir_paths <- list(
  "data/spedizione/tubeworms_mussels/Oxic Pressure",
  #"data/spedizione/tubeworms_mussels/Anoxic Pressure - No10"
  "data/spedizione/chimney/Oxic Pressure - No10",
  #"data/spedizione/chimney/Oxic Pressure",
  "data/nioz"
)

########### ATTENTION!!!!!!! ##########
# params start values for nls model (must be the same element numbers as dir_paths)
nls_param_list <- list(
  list(100, 30, 4),
  list(100, 30, 4),
  list(0.5, 33, 77)
)

anova_data <- data.frame()
anova_slopes <- data.frame()
t_test_data <- data.frame()
k <- 0
# Iterate over directories
for (dir_path in dir_paths) {
  k <- k + 1
  chart_subtitle <- chart_subtitle_func(dir_path)
  non_linear_regression_result <- non_linear_regression_dir_func(
    dir_path = dir_path,
    nls_param_list = nls_param_list,
    k = k,
    chart_subtitle = chart_subtitle
  )
  linear_regression_result <- linear_regression_func(
    dir_path = dir_path,
    ld50 = non_linear_regression_result$ld50,
    time_list = non_linear_regression_result$time_list,
    chart_subtitle = chart_subtitle
  )
  anova_data <- rbind(anova_data, linear_regression_result$anova_data)
  anova_slopes <- rbind(anova_slopes, linear_regression_result$anova_slopes)
  t_test_data <- rbind(t_test_data, linear_regression_result$t_test_data)
}

anova_analysis(
  anova_data = anova_data,
  anova_slopes = anova_slopes
)

t_test_func(
  t_test_data = t_test_data
)
