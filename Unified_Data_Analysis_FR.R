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

# Specify the directory
dir_path <- "data/spedizione/chimney/Oxic Pressure - No10"

# Get a list of all files in the directory
csv_files <- list.files(path = dir_path)

# Exposure time is in the file name: it must be extracted and converted to a numeric value
time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
time_list <- as.numeric(time_character)

# These will be replaced with their true values in a following loop
ld50 <- c()
param_1 <- 0
param_2 <- 0
param_3 <- 0
param_4 <- time_list[1]
param_z <- 0

tdt_line <- function(z) {
  return(tdt_intercept + tdt_slope * z)
}

survival_func <- function(x, y) {
  return(param_1 / (1 + (((x - (param_z * log10(y / param_4))) / param_2)^param_3)))
}

#######################################################################################################################

# Create the LT50 Survival curve plot and calculate relevant parameter values

# This is used to find the correct axis values for the plot
min_temperature <- .Machine$double.xmax
max_temperature <- .Machine$double.xmin

# Iterate over the files to create the plot
for (csv_file in csv_files) {
  # Find correct axis values to assign the plot
  csv_path <- sprintf("%s\\%s", dir_path, csv_file)
  data <- read.csv(csv_path, row.names = NULL)

  if (min(data$Temperature) < min_temperature) {
    min_temperature <- min(data$Temperature)
  }
  if (max(data$Temperature) > max_temperature) {
    max_temperature <- max(data$Temperature)
  }
}

# Create empty plot with correct axis values and labels
plot(1,
  type = "n",
  xlim = c(min_temperature, max_temperature),
  ylim = c(0.0, 1.0),
  xlab = "Temperature (°C)",
  ylab = "Proportional Survival",
  main = paste("LT50 Survival Curve\n", dir_path)
)

# Assign a color to each data set
lines_color <- c("green", "red", "blue", "orange", "purple")

# This is used to number things in the loop
i <- 0

# Iterate data analysis and plotting over the files
cat("########## NON-LINEAR REGRESSION ##########\n\n")
for (csv_file in csv_files) {
  i <- i + 1

  # Fit binomial distribution
  csv_path <- sprintf("%s/%s", dir_path, csv_file)
  cat("########## FILE: ", csv_path, " ##########\n\n")
  data <- read.csv(csv_path, row.names = NULL)
  data$Survival <- data$Alive / (data$Alive + data$Dead)

  nls_p1_start <- 100
  nls_p2_start <- 30
  nls_p3_start <- 4
  xdata <- data[, 1]
  ydata <- data[, 4]
  fit <- nls(ydata ~ p1 / (1 + (xdata / p2)^p3),
    start = list(p1 = nls_p1_start, p2 = nls_p2_start, p3 = nls_p3_start),
    control = nls.control(maxiter = 500)
  )
  new <- data.frame(xdata = seq(min(xdata), max(xdata), len = 200))

  # The summary() function gives important info for statistical analysis

  fit_params <- summary(fit)
  cat("Summary statistics for", time_list[i], "h exposure\n")
  print(summary(fit))

  # Add the model_dose value to the LD50 list
  ld50 <- c(ld50, fit_params$parameters[2])

  # Plot non-linear chart
  model_dose <- fit_params$parameters[2]
  rounded_dose <- round(model_dose, digits = 2)
  model_dose_y <- predict(fit, newdata = data.frame(xdata = model_dose))
  model_dose_stderr <- fit_params$parameters[2, 2]
  lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])
  points(jitter(xdata, factor = 0.5), ydata, col = lines_color[i])
  text(
    x = rounded_dose,
    y = model_dose_y,
    labels = sprintf("%.2f°C", model_dose),
    pos = cos(pi * i) + 2,
    cex = 0.8
  )
  # Adding flex points
  points(
    x = model_dose,
    y = model_dose_y,
    type = "p",
    pch = 19,
    col = lines_color[i]
  )
  # Adding std error ranges
  arrows(
    x0 = model_dose - model_dose_stderr,
    y0 = model_dose_y,
    x1 = model_dose + model_dose_stderr,
    y1 = model_dose_y,
    code = 3,
    angle = 90,
    length = 0.05,
    col = lines_color[i]
  )

  # In the 3D plot the parameters of the shortest exposure are used in the function
  if (i < 2) {
    param_1 <- fit_params$parameters[1]
    param_2 <- fit_params$parameters[2]
    param_3 <- fit_params$parameters[3]
  }
}

# Add a legend to the plot
j <- 0
exposures <- c()

for (n in time_list) {
  j <- j + 1
  legend_line <- c(paste(as.character(time_list[j]), "h exposure"))
  exposures <- c(exposures, legend_line)
}

legend("topright",
  inset = 0.02,
  legend = exposures,
  col = lines_color,
  lty = 1,
  lwd = 2,
  cex = 0.7
)

# Create the Thermal death time curve (TDT) plot

# Create a data frame with LD50 and LT50 to construct TDT curve
lt50 <- data.frame(ld50, time_list)
cat("########## LINEAR REGRESSION ##########\n\n")
cat("LD50 values and corresponding LT50 exposure values\n")
print(lt50)

# Run a linear regression
tdt_results <- lm(lt50$ld50 ~ log10(lt50$time_list), data = lt50)
tdt_slope <- tdt_results$coefficients["log10(lt50$time_list)"]
tdt_intercept <- tdt_results$coefficients["(Intercept)"]

# Parameter used in 3D plot
param_z <- tdt_slope

# Create plot with LD50/LT50 values in log-scale
plot(log10(lt50$time_list), lt50$ld50,
  xlim = c(0, 2), ylim = c(0, 55),
  xlab = expression("log"[10] * "Time (LT50, hours)"), ylab = "Temperature (LD50, °C)",
  main = "Thermal death time (TDT) Curve"
)

# Create a vector of LD50 values and use it as input to plot the curve
z <- seq(0, 2, 0.01)
lines(z, tdt_line(z), new = TRUE)

# Calculate R squared goodness of fit for the TDT curve
cat("Summary statistics for TDT curve\n")
print(summary(tdt_results))
cat("R squared value for TDT curve: ", summary(tdt_results)$r.squared, "\n\n")

# Create the 3D Thermal tolerance landscape plot

# Create a sequence of temperature values from 25.0 to 45.0 with an increment of 0.5.
temperature <- seq(0.0, 70.0, by = 0.5)

# Create a sequence of time values from 3.0 to 102.0 with an increment of 3.0.
# Add two more values 0.1 and 1 at the start of the time sequence.
time <- seq(param_4, 100.0, by = 1)
time <- c(0.1, 1, time)

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
  z = ~survival_rate
) # Assign survival rate values to the z-axis
fig <- fig %>% layout(
  scene = list(
    xaxis = list(nticks = 5),
    zaxis = list(nticks = 5),
    camera = list(eye = list(x = 0, y = -1, z = 0.5)),
    aspectratio = list(x = .8, y = .9, z = 0.5)
  )
)

fig
