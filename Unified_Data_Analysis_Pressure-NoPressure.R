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
  "data/spedizione/anova_analysis/chimney")

# Initialize dataframes for ANOVA analysis 
# Dataframe for LD50
anova_data = data.frame(LD50 = numeric(0), time_list = numeric(0), dir = character((0)))
# Dataframe for linear regression params
anova_slopes <- data.frame(slope = numeric(0), intercept = numeric(0), dir = character(0))
nls_coefficients = list()
min_df = .Machine$double.xmax

# Iterate over directories
for (dir_path in dir_paths){
  # Get a list of all files in the directory
  csv_files <- list.files(path = dir_path)
  
  # Exposure time is in the file name: it must be extracted and converted to a numeric value
  time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
  file_prefixes <- gsub("^(.{3}).*", "\\1", csv_files)
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
    p1_start = 100
    p2_start = 30
    p3_start = 4
    xdata<-data[,1]
    ydata<-data[,4]
    # Execute regression
    fit = nls(ydata ~ p1 / (1 + (xdata / p2) ^ p3), start=list(p1 = p1_start, p2 = p2_start, p3 = p3_start),
      control = nls.control(maxiter = 500))
    current_df = df.residual(fit)
    if (current_df < min_df){
      min_df <- current_df
    }
    
    # Storing model coefficients
    nls_coefficients <- append(nls_coefficients, coef(summary(fit)))
    
    # Syntetic dataframe initialization to plot distrubution line
    new = data.frame(xdata = seq(min(xdata), max(xdata), len = 200))
    lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])
    # Adding real points
    points(xdata, ydata, col = lines_color[i])
    
    # The summary() function gives important info for statistical analysis
    cat("\n")
    cat("Summary statistics for", substr(csv_file, 1, 3), " ", time_list[i], "h exposure")
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
    legend_line = c(paste(file_prefixes[j], " ", as.character(time_list[j]), "h exposure"))
    exposures = c(exposures, legend_line)
  }
  
  legend(
    "topright", 
    inset = 0.02, 
    legend = exposures, 
    col = lines_color, 
    lty = 1, 
    lwd = 3, 
    cex = 0.9)
}

###########################################################################################
# ANOVA model analysis test
# p1 Estimate first model
beta1 <- nls_coefficients[[1]]
# p1 Std. Error first model
se1 <- nls_coefficients[[4]]
# p1 Estimate second model
beta2 <- nls_coefficients[[13]]
# p1 Std. Error second model
se2 <- nls_coefficients[[16]]
t_statistics <- (beta1 - beta2) / sqrt(se1^2 + se2^2)
p_value <- 2 * (1 - pt(abs(t_statistics), min_df))