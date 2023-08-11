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

# Define the functions necessary to fit the data to the desired distribution
logistic_line <- function(x, model_intercept, model_slope) {
  eta <- model_intercept + model_slope * x;
  return (1 / (1 + exp(-eta)))
}

TDT_line = function(z) {
  eta = TDT_intercept + TDT_slope * z
}

survival_func <- function(x, y) {
  param1
  param2
  param3
  paramZ
  param4
  return (param1 / (1 + (((x - (paramZ * log10(y / param4))) / param2) ^ param3)))
}

# Specify the directory
dir_path <- "data/spedizione/chimney/Oxic Pressure"
#dir_path <- "data/spedizione/chimney/Oxic Pressure"

# Get a list of all files in the directory
csv_files <- list.files(path = dir_path)

# Exposure time is in the file name: it must be extracted and converted to a numeric value
time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
time_list <- as.numeric(time_character)

# These will be replaced with their true values in a following loop
LD50 <- c()
param1 <- 0
param2 <- 0
parma3 <- 0
paramZ <- 0
param4 <- time_list[1]

#######################################################################################################################

# Create the LT50 Survival curve plot and calculate relevant parameter values

# This is used to find the correct axis values for the plot
min_temperature <- .Machine$double.xmax
max_temperature <- .Machine$double.xmin

# Iterate over the files to create the plot
for (csv_file in csv_files){
  
  # Find correct axis values to assign the plot
  csv_path <- sprintf("%s\\%s",dir_path, csv_file)
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
     main = "LT50 Survival Curve"
)

# Assign a color to each data set
lines_color <- c("green", "red", "blue", "orange", "purple")

# This is used to number things in the loop
i <- 0

# Iterate data analysis and plotting over the files
for(csv_file in csv_files){
  i <- i + 1
  
  # Fit binomial distribution
  csv_path <- sprintf("%s/%s",dir_path, csv_file)
  data <- read.csv(csv_path, row.names = NULL)
  data$Survival = data$Alive / (data$Alive + data$Dead)
  
  p1=100
  p2=30
  p3=4
  xdata<-data[,1]
  ydata<-data[,4]
  fit = nls(ydata ~ p1 / (1 + (xdata / p2) ^ p3), start=list(p1 = p1, p2 = p2, p3 = p3),
    control = nls.control(maxiter = 500))
  new = data.frame(xdata = seq(min(xdata), max(xdata), len=200))
  
  lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])
  points(xdata, ydata, col = lines_color[i])
  
  #lines(jitter(data$Temperature, factor = 0.5), (data$Alive / (data$Alive + data$Dead)), 
  #      text(rounded_dose + 0.1, 0.5, sprintf("%.2f°C", model_dose), pos = cos(pi*i) + 2, cex = 0.8), 
  #      las = 1, col = lines_color[i])
  
  #points(x = model_dose, y = 0.5, type = "p", pch = 19, col = lines_color[i])
  
  #arrows(x0 = model_dose - model_stderr, y0 = 0.5, 
  #       x1 = model_dose + model_stderr, y1 = 0.5, 
  #       code = 3, angle = 90, length = 0.05)
  
  # The summary() function gives important info for statistical analysis
  cat("\n")
  cat("Summary statistics for", time_list[i], "h exposure")
  #print(summary(model_results))
  fit_params = summary(fit)
  print(summary(fit))
  #confint(fit, level = 0.90)
  
  # Calculate R squared goodness of fit 
  #cat("R squared value for", time_list[i], "h exposure:", with(summary(model_results), 1 - deviance/null.deviance), "\n")
  #cat("-------------------------------------------", "\n")
  
  # Add the model_dose value to the LD50 list
  LD50 = c(LD50, fit_params$parameters[2])
  
  # In the 3D plot the parameters of the shortest exposure are used in the function
  if(i < 2){
    param1 = fit_params$parameters[1]
    param2 = fit_params$parameters[2]
    param3 = fit_params$parameters[3]
  }
}

#Add a legend to the plot
j <- 0
exposures <- c()

for(n in time_list){
  j <- j + 1
  legend_line = c(paste(as.character(time_list[j]), "h exposure"))
  exposures = c(exposures, legend_line)
}

legend("topright", inset = 0.02, 
       legend = exposures, 
       col = lines_color, 
       lty = 1, lwd = 3, cex = 0.9)

#######################################################################################################################

# Create the Thermal death time curve (TDT) plot

# Create a data frame with LD50 and LT50 to construct TDT curve
LT50 <- data.frame(LD50, time_list)
cat("\n")
cat("LD50 values and corresponding LT50 exposure values", "\n")
print(LT50)
cat("\n")

# Run a linear regression
TDT_results <- lm(LT50$LD50 ~ log10(LT50$time_list), data = LT50)
TDT_slope <- TDT_results$coefficients["log10(LT50$time_list)"]
TDT_intercept <- TDT_results$coefficients["(Intercept)"]

# Parameter used in 3D plot
paramZ = TDT_slope

# Create plot with LD50/LT50 values in log-scale
plot(log10(LT50$time_list), LT50$LD50, 
     xlim = c(0, 2), ylim = c(0, 55), 
     xlab = expression("log"[10]*"Time (LT50, hours)"), ylab = "Temperature (LD50, °C)",
     main = "Thermal death time (TDT) Curve")

# Create a vector of LD50 values and use it as input to plot the curve
z = seq(0, 2, 0.01)
lines(z, TDT_line(z), new = TRUE)

# Calculate R squared goodness of fit for the TDT curve
cat("Summary statistics for TDT curve")
print(summary(TDT_results))
cat("R squared value for TDT curve:", summary(TDT_results)$r.squared)

#######################################################################################################################

# Create the 3D Thermal tolerance landscape plot

# Create a sequence of temperature values from 25.0 to 45.0 with an increment of 0.5.
temperature = seq(0.0, 70.0, by =0.5)

# Create a sequence of time values from 3.0 to 102.0 with an increment of 3.0.
# Add two more values 0.1 and 1 at the start of the time sequence.
time = seq(param4, 100.0, by = 1)
time = c(0.1, 1, time)

# Calculate the survival rate for each combination of temperature and time values.
survival_rate <- outer(temperature, time, survival_func)

# Define a custom color scale that varies from red to yellow to green.
color_scale <- list(c(0, 0.5, 1), c('red', 'yellow', 'green'))

# Create a 3D plot of the survival rate as a function of temperature and time.
fig <- plot_ly(
  type = 'surface',
  opacity = 0.5,  # Set surface opacity
  colorscale = color_scale,  # Use custom color scale
  # Show contour lines on the surface
  contours = list(
    x = list(show = TRUE, start = time[1], end = tail(time, n = 1), size = time[4] - time[3], color = 'white'),
    y = list(show = TRUE, start = temperature[1], end = tail(temperature, n = 1), size = temperature[4] - temperature[3], color = 'white'),
    z = list(show = TRUE, start = 0, end = 1, size = 0.05)),
  x = ~time,  # Assign time values to the x-axis
  y = ~temperature,  # Assign temperature values to the y-axis
  z = ~survival_rate)  # Assign survival rate values to the z-axis
fig <- fig %>% layout(
  scene = list(
    xaxis = list(nticks = 5),
    zaxis = list(nticks = 20),
    camera = list(eye = list(x = 0, y = -1, z = 0.5)),
    aspectratio = list(x = .8, y = .9, z = 0.5)))

fig

