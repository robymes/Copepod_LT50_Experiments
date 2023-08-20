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

source("plot_utilities.R")
source("non_linear_regression.R")
source("linear_regression.R")

survival_param_1 <- 0
survival_param_2 <- 0
survival_param_3 <- 0
survival_param_4 <- 0
survival_param_z <- 0

# Directories where to scan for csv ata
dir_path <- "data/spedizione/tubeworms_mussels/Anoxic Pressure - No10"

########### ATTENTION!!!!!!! ##########
# params start values for nls model (must be the same element numbers as dir_paths)
nls_param_list <- list(
  list(100, 30, 4)
  #list(0.5, 33, 77)
)

survival_func <- function(x, y) {
  return(survival_param_1 /
           (1 + (((x - (survival_param_z * log10(y / survival_param_4))) / survival_param_2)^survival_param_3)))
}

chart_subtitle <- chart_subtitle_func(dir_path)
non_linear_regression_result <- non_linear_regression_dir_func(
  dir_path = dir_path,
  nls_param_list = nls_param_list,
  k = 1,
  chart_subtitle = chart_subtitle
)

linear_regression_result <- linear_regression_func(
  dir_path = dir_path,
  ld50 = non_linear_regression_result$ld50,
  time_list = non_linear_regression_result$time_list,
  chart_subtitle = chart_subtitle
)

survival_param_1 <- non_linear_regression_result$survival_params[1]
survival_param_2 <- non_linear_regression_result$survival_params[2]
survival_param_3 <- non_linear_regression_result$survival_params[3]
survival_param_4 <- non_linear_regression_result$survival_params[4]
survival_param_z <- linear_regression_result$survival_param_z

# Create the 3D Thermal tolerance landscape plot

# Create a sequence of temperature values from 25.0 to 45.0 with an increment of 0.5.
temperature <- seq(0.0, 60.0, by = 0.5)

# Create a sequence of time values from 3.0 to 102.0 with an increment of 3.0.
# Add two more values 0.1 and 1 at the start of the time sequence.
time <- seq(survival_param_4, 100.0, by = 1)
time <- c(1, time)

# Calculate the survival rate for each combination of temperature and time values.
survival_rate <- outer(temperature, time, survival_func)

# Define a custom color scale that varies from red to yellow to green.
color_scale <- list(c(0, 0.5, 1), c("red", "yellow", "green"))

# Create a 3D plot of the survival rate as a function of temperature and time.
fig <- plot_ly(
  type = "surface",
  opacity = 0.5, # Set surface opacity
  colorscale = color_scale, # Use custom color scale
  # Show contour lines on the surface
  contours = list(
    x = list(
      show = TRUE,
      start = time[1],
      end = tail(time, n = 1),
      size = time[4] - time[3],
      color = "white"
    ),
    y = list(
      show = TRUE,
      start = temperature[1],
      end = tail(temperature, n = 1),
      size = temperature[4] - temperature[3],
      color = "white"
    ),
    z = list(
      show = TRUE,
      start = 0,
      end = 1,
      size = 0.05
    )
  ),
  x = ~time, # Assign time values to the x-axis
  y = ~temperature, # Assign temperature values to the y-axis
  z = ~survival_rate # Assign survival rate values to the z-axis
)
fig <- fig %>% layout(
  scene = list(
    xaxis = list(nticks = 5),
    zaxis = list(nticks = 5),
    camera = list(eye = list(x = 0, y = -1, z = 0.5)),
    aspectratio = list(x = .8, y = .9, z = 0.5)
  )
)

fig
