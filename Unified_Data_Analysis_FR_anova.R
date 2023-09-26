# This script combines and automates the entire data analysis for LD50/LT50 experiments

# Files must be named in the format Dataxx.csv where xx is a number (01, 02, 03... 10)
# All data used in the analysis should be in one folder with no other extra files
# Run code using ctrl + shift + s (Source) to see uncluttered statistics in the console

source("plot_utilities.R")
source("non_linear_regression.R")
source("linear_regression.R")

check_os_func()

main_title <- "tdt pompeii worm-tubeworm-intertidal oxic"

# Directories where to scan for csv data
dir_paths <- list(
  "data/cruise/tubeworms_mussels/Oxic Pressure",
  #"data/cruise/tubeworms_mussels/Anoxic Pressure"
  "data/cruise/pompeii_worms/Oxic Pressure",
  "data/nioz/intertidal/Oxic Pressure"
)

########### ATTENTION!!!!!!! ##########
# params start values for nls model (must be the same element numbers as dir_paths)
nls_param_list <- list(
  list(100, 30, 4),
  list(100, 30, 4),
  list(0.5, 33, 77)
)

cat("########## EXECUTION FOR: ", main_title, " ##########\n\n")

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
    chart_subtitle = chart_subtitle,
    main_title = main_title
  )
  linear_regression_result <- linear_regression_func(
    dir_path = dir_path,
    ld50 = non_linear_regression_result$ld50,
    time_list = non_linear_regression_result$time_list,
    main_title = main_title,
    chart_subtitle = chart_subtitle
  )
  anova_data <- rbind(anova_data, linear_regression_result$anova_data)
  anova_slopes <- rbind(anova_slopes, linear_regression_result$anova_slopes)
  t_test_data <- rbind(t_test_data, linear_regression_result$t_test_data)
}

anova_analysis_func(
  anova_data = anova_data,
  anova_slopes = anova_slopes,
  main_title = main_title
)

t_test_func(
  t_test_data = t_test_data,
  main_title = main_title
)
