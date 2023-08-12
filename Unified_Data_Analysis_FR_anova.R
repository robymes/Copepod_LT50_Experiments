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

TDT_line = function(slope, intercept, z) {
  eta = intercept + slope * z
}

survival_func <- function(x, y) {
  param1
  param2
  param3
  paramZ
  param4
  return (param1 / (1 + (((x - (paramZ * log10(y / param4))) / param2) ^ param3)))
}

# Directories where to scan for csv data
dir_paths = list(
  "data/spedizione/tubeworms _mussels/Oxic Pressure", 
  "data/spedizione/chimney/Oxic Pressure"
  #"data/nioz"
)

# params start values for nls model (must be the same element numbers as dir_paths)
nls_param_list <- list(
  list(100, 30, 4),
  list(100, 30, 4)
  #list(0.5, 33, 77)
)

# Initialize dataframes for ANOVA analysis 
# Dataframe for LD50
anova_data = data.frame(LD50 = numeric(0), time_list = numeric(0), dir = character((0)))
# Dataframe for linear regression params
anova_slopes <- data.frame(slope = numeric(0), intercept = numeric(0), dir = character(0))
# Dataframe for T test
t_test_data <- data.frame(slope = numeric(0), slope_std_err = numeric(0),
  intercept = numeric(0), intercept_std_err = numeric(0), data_points_length = numeric(0),
  dir = character(0))

k = 0
# Iterate over directories
k = 0
for (dir_path in dir_paths){
  k = k + 1
  # Get a list of all files in the directory
  csv_files <- list.files(path = dir_path)
  
  # Exposure time is in the file name: it must be extracted and converted to a numeric value
  time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
  time_list <- as.numeric(time_character)
  
  # Params initialization
  # These will be replaced with their true values in a following loop
  LD50 <- c()
  
  #######################################################################################################################
  
  # Create the LT50 Survival curve plot and calculate relevant parameter values
  
  # This is used to find the correct axis values for the plot
  min_temperature <- .Machine$double.xmax
  max_temperature <- .Machine$double.xmin
  
  # Iterate over the csv files to find plot x-axis extension
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
  
  # Iterate over csv files and create plot
  for(csv_file in csv_files){
    i <- i + 1
    
    # Fit binomial distribution
    csv_path <- sprintf("%s/%s",dir_path, csv_file)
    data <- read.csv(csv_path, row.names = NULL)
    # Calculate survival rate
    data$Survival = data$Alive / (data$Alive + data$Dead)
    
    # Initialize regression params
    p1_start = nls_param_list[[k]][[1]]
    p2_start = nls_param_list[[k]][[2]]
    p3_start = nls_param_list[[k]][[3]]
    xdata<-data[,1]
    ydata<-data[,4]
    # Execute regression
    fit = nls(ydata ~ p1 / (1 + (xdata / p2) ^ p3), start=list(p1 = p1_start, p2 = p2_start, p3 = p3_start),
      control = nls.control(maxiter = 500))
    
    # Syntetic dataframe initialization to plot distrubution line
    new = data.frame(xdata = seq(min(xdata), max(xdata), len = 200))
    lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])
    # Adding real points
    points(xdata, ydata, col = lines_color[i])
    
    # The summary() function gives important info for statistical analysis
    cat("\n")
    cat("Summary statistics for", time_list[i], "h exposure")
    fit_params = summary(fit)
    print(summary(fit))
    
    # Add the model_dose value to the LD50 list
    LD50 = c(LD50, fit_params$parameters[2])
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
  LT50_temp <- LT50
  LT50_temp$dir = dir_path 
  anova_data = rbind(anova_data, LT50_temp)
  cat("\n")
  cat("LD50 values and corresponding LT50 exposure values", "\n")
  print(LT50)
  cat("\n")
  
  # Run a linear regression
  TDT_results <- lm(LT50$LD50 ~ log10(LT50$time_list), data = LT50)
  TDT_slope <- TDT_results$coefficients["log10(LT50$time_list)"]
  TDT_intercept <- TDT_results$coefficients["(Intercept)"]
  
  new_anova_element <- data.frame(slope = TDT_slope, intercept = TDT_intercept, dir = dir_path)
  anova_slopes <- rbind(anova_slopes, new_anova_element)
  
  new_t_test_data <- data.frame(
    slope = summary(TDT_results)$coefficients[2, 1],
    slope_std_err = summary(TDT_results)$coefficients[2, 2],
    intercept = summary(TDT_results)$coefficients[1, 1], 
    intercept_std_err = summary(TDT_results)$coefficients[1, 2], 
    data_points_length = length(LT50_temp$time_list),
    dir = dir_path)
  t_test_data <- rbind(t_test_data, new_t_test_data)
  
  # Parameter used in 3D plot
  paramZ = TDT_slope
  
  # Create plot with LD50/LT50 values in log-scale
  plot(log10(LT50$time_list), LT50$LD50, 
       xlim = c(0, 2), ylim = c(0, 55), 
       xlab = expression("log"[10]*"Time (LT50, hours)"), ylab = "Temperature (LD50, °C)",
       main = "Thermal death time (TDT) Curve")
  
  # Create a vector of LD50 values and use it as input to plot the curve
  z = seq(0, 2, 0.01)
  lines(z, TDT_line(TDT_slope, TDT_intercept, z), new = TRUE)
  
  # Calculate R squared goodness of fit for the TDT curve
  cat("Summary statistics for TDT curve")
  print(summary(TDT_results))
  cat("R squared value for TDT curve:", summary(TDT_results)$r.squared) 
}

########################################################################################
# Plot ANOVA analysis
anova_data$time_list = log10(anova_data$time_list)
anova_plot <- ggplot(
  anova_data, 
  aes(
    x = time_list, 
    y = LD50, 
    color = factor(dir)
    )
  ) +
  geom_point() +
  labs(
    x = "LOG10 Time (LT50, hours)", 
    y = "Temperature (LD50, °C)", 
    title = "Thermal death time (TDT) Curve") +
  scale_color_discrete(name = "dir")
for (i in 1:nrow(anova_slopes)){
  z = seq(0, 2, 0.01)
  anova_slope_data <- TDT_line(anova_slopes[i, "slope"], anova_slopes[i, "intercept"], z)
  anova_slope_df <- data.frame(LD50 = anova_slope_data, time_list = z, dir = anova_slopes[i, "dir"])
  anova_plot <- anova_plot +
    geom_line(data = anova_slope_df, aes(group = dir))
}
print(anova_plot)
print("ANOVA verification")
anova(lm(LD50~log(time_list)*dir, anova_data))

######################################################################################
# T test slope
t_value <- (t_test_data$slope[1] - t_test_data$slope[2]) / sqrt(t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)
df_value <- (t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)^2 / (t_test_data$slope_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) + t_test_data$slope_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
p_value <- 2 * (1 - pt(abs(t_value), df_value))
print(paste("P Value slopes: ", p_value))

t_value <- (t_test_data$intercept[1] - t_test_data$intercept[2]) / sqrt(t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)
df_value <- (t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)^2 / (t_test_data$intercept_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) + t_test_data$intercept_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
p_value <- 2 * (1 - pt(abs(t_value), df_value))
print(paste("P Value intercept: ", p_value))