source("plot_utilities.R")
source("non_linear_regression.R")
source("linear_regression.R")

check_os_func()

dir_path <- "data/cruise/pompeii_worms/Oxic Pressure"
csv_file <- c("Data08.csv")

non_linear_regression_dir_func(
  dir_path = dir_path,
  nls_param_list = list(list(100, 30, 4)),
  k = 1,
  chart_subtitle = chart_subtitle_func(dir_path),
  csv_file_list = c(csv_file)
)