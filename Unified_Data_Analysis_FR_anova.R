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
  eta <- model_intercept + model_slope * x
  return(1 / (1 + exp(-eta)))
}

tdt_line <- function(slope, intercept, z) {
  return(intercept + slope * z)
}

# Directories where to scan for csv data
dir_paths <- list(
  "data/spedizione/tubeworms_mussels/Oxic Pressure",
  "data/spedizione/tubeworms_mussels/Anoxic Pressure - No10"
  #"data/spedizione/chimney/Oxic Pressure - No10"
  #"data/spedizione/chimney/Oxic Pressure"
  #"data/nioz"
)

# params start values for nls model (must be the same element numbers as dir_paths)
nls_param_list <- list(
  list(100, 30, 4),
  list(100, 30, 4)
  # list(0.5, 33, 77)
)

# Check if dir_paths and nls_param_list are the same length
stopifnot(length(dir_paths) == length(nls_param_list), "dir_paths and nls_param_list are not the same length!")

# Initialize dataframes for ANOVA analysis
# Dataframe for LD50
anova_data <- data.frame(
  ld50 = numeric(0),
  time_list = numeric(0),
  dir = character((0))
)
# Dataframe for linear regression params
anova_slopes <- data.frame(
  slope = numeric(0),
  intercept = numeric(0),
  dir = character(0)
)
# Dataframe for T test
t_test_data <- data.frame(
  slope = numeric(0),
  slope_std_err = numeric(0),
  intercept = numeric(0),
  intercept_std_err = numeric(0),
  data_points_length = numeric(0),
  dir = character(0)
)

k <- 0
# Iterate over directories
k <- 0
for (dir_path in dir_paths) {
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")
  k <- k + 1
  # Get a list of all files in the directory
  csv_files <- list.files(path = dir_path)

  # Exposure time is in the file name: it must be extracted and converted to a numeric value
  time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
  time_list <- as.numeric(time_character)

  # Params initialization
  # These will be replaced with their true values in a following loop
  ld50 <- c()

  cat("########## NON LINEAR REGRESSION ##########\n\n")

  # Create the LT50 Survival curve plot and calculate relevant parameter values

  # This is used to find the correct axis values for the plot
  min_temperature <- .Machine$double.xmax
  max_temperature <- .Machine$double.xmin

  # Iterate over the csv files to find plot x-axis extension
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
    main = "LT50 Survival Curve"
  )

  # Assign a color to each data set
  lines_color <- c("green", "red", "blue", "orange", "purple")

  # This is used to number things in the loop
  i <- 0

  # Iterate over csv files and create plot
  for (csv_file in csv_files) {
    i <- i + 1

    # Fit binomial distribution
    csv_path <- sprintf("%s/%s", dir_path, csv_file)
    cat("########## FILE: ", csv_path, " ###########\n\n")
    data <- read.csv(csv_path, row.names = NULL)
    # Calculate survival rate
    data$Survival <- data$Alive / (data$Alive + data$Dead)

    # Initialize regression params
    p1_start <- nls_param_list[[k]][[1]]
    p2_start <- nls_param_list[[k]][[2]]
    p3_start <- nls_param_list[[k]][[3]]
    xdata <- data[, 1]
    ydata <- data[, 4]
    # Execute non-linear regression
    fit <- nls(ydata ~ p1 / (1 + (xdata / p2)^p3),
      start = list(p1 = p1_start, p2 = p2_start, p3 = p3_start),
      control = nls.control(maxiter = 500)
    )

    # Syntetic dataframe initialization to plot distrubution line
    new <- data.frame(xdata = seq(min(xdata), max(xdata), len = 200))

    # The summary() function gives important info for statistical analysis
    cat("Summary statistics for ", time_list[i], "h exposure\n")
    fit_params <- summary(fit)
    print(summary(fit))
    cat("\n")

    # Plot non-linear chart
    model_dose <- fit_params$parameters[2]
    rounded_dose <- round(model_dose, digits = 2)
    model_dose_y <- predict(fit, newdata = data.frame(xdata = model_dose))
    model_dose_stderr <- fit_params$parameters[2, 2]
    lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])
    points(xdata, ydata, col = lines_color[i])
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

    # Add the model_dose value to the LD50 list
    ld50 <- c(ld50, fit_params$parameters[2])

    # R squared calculation for fit model
    rss <- sum(residuals(fit)^2)
    tss <- sum((ydata - mean(ydata))^2)
    r_squared <- 1 - (rss / tss)
    cat("R squared for fit model: ", r_squared, "\n\n")
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
    lty = 1, lwd = 3, cex = 0.9
  )

  # Create the Thermal death time curve (TDT) plot

  cat("########## LINEAR REGRESSION ##########\n\n")
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")

  # Create a data frame with LD50 and LT50 to construct TDT curve
  lt50 <- data.frame(ld50, time_list)
  lt50_temp <- lt50
  lt50_temp$dir <- dir_path
  anova_data <- rbind(anova_data, lt50_temp)
  cat("LD50 values and corresponding LT50 exposure values\n")
  print(lt50)
  cat("\n")

  # Run a linear regression with time in hours
  tdt_results_hours <- lm(lt50$ld50 ~ log10(lt50$time_list), data = lt50)
  tdt_slope <- tdt_results_hours$coefficients["log10(lt50$time_list)"]
  tdt_intercept_hours <- tdt_results_hours$coefficients["(Intercept)"]

  new_anova_element <- data.frame(
    slope = tdt_slope,
    intercept = tdt_intercept_hours,
    dir = dir_path
  )
  anova_slopes <- rbind(anova_slopes, new_anova_element)

  new_t_test_data <- data.frame(
    slope = summary(tdt_results_hours)$coefficients[2, 1],
    slope_std_err = summary(tdt_results_hours)$coefficients[2, 2],
    intercept = summary(tdt_results_hours)$coefficients[1, 1],
    intercept_std_err = summary(tdt_results_hours)$coefficients[1, 2],
    data_points_length = length(lt50_temp$time_list),
    dir = dir_path
  )
  t_test_data <- rbind(t_test_data, new_t_test_data)

  # Create plot with LD50/LT50 values in log-scale
  plot(log10(lt50$time_list), lt50$ld50,
    xlim = c(0, 2), ylim = c(0, 55),
    xlab = expression("log"[10] * "Time (LT50, hours)"),
    ylab = "Temperature (LD50, °C)",
    main = "Thermal death time (TDT) Curve"
  )

  # Create a vector of LD50 values and use it as input to plot the curve
  z <- seq(0, 2, 0.01)
  lines(z, tdt_line(tdt_slope, tdt_intercept_hours, z), new = TRUE)

  # Calculate R squared goodness of fit for the TDT curve
  cat("Summary statistics for TDT curve\n")
  print(summary(tdt_results_hours))
  cat("R squared value for TDT curve: ", summary(tdt_results_hours)$r.squared, "\n\n")

  # Run a linear regression with time in minutes
  lt50_minutes <- lt50
  lt50_minutes$time_list <- lt50_minutes$time_list * 60
  tdt_results_minutes <- lm(lt50_minutes$ld50 ~ log10(lt50_minutes$time_list), data = lt50_minutes)
  tdt_intercept_minutes <- tdt_results_minutes$coefficients["(Intercept)"]
  cat("TDT intercept in minutes: ", tdt_intercept_minutes, "\n\n")
}

# Plot ANOVA analysis
anova_data$time_list <- log10(anova_data$time_list)
anova_plot <- ggplot(
  anova_data,
  aes(
    x = time_list,
    y = ld50,
    color = factor(dir)
  )
) +
  geom_point() +
  labs(
    x = "LOG10 Time (LT50, hours)",
    y = "Temperature (LD50, °C)",
    title = "Thermal death time (TDT) Curve"
  ) +
  scale_color_discrete(name = "dir")
for (i in 1:nrow(anova_slopes)) { # nolint: seq_linter.
  z <- seq(0, 2, 0.01)
  anova_slope_data <- tdt_line(anova_slopes[i, "slope"], anova_slopes[i, "intercept"], z)
  anova_slope_df <- data.frame(ld50 = anova_slope_data, time_list = z, dir = anova_slopes[i, "dir"])
  anova_plot <- anova_plot +
    geom_line(data = anova_slope_df, aes(group = dir))
}
print(anova_plot)
anova_analysis <- anova(lm(ld50 ~ log(time_list) * dir, anova_data))
cat("########## ANOVA ANALYSIS ##########\n\n")
print(anova_analysis)
cat("\n")

# T test slope
t_value_slope <-
  (t_test_data$slope[1] - t_test_data$slope[2]) /
  sqrt(t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)
df_value_slope <-
  (t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)^2 /
  (t_test_data$slope_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) +
   t_test_data$slope_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
p_value_slope <- 2 * (1 - pt(abs(t_value_slope), df_value_slope))

t_value_intercept <-
  (t_test_data$intercept[1] - t_test_data$intercept[2]) /
  sqrt(t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)
df_value_intercept <-
  (t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)^2 /
  (t_test_data$intercept_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) +
   t_test_data$intercept_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
p_value_intercept <- 2 * (1 - pt(abs(t_value_intercept), df_value_intercept))
cat("########## T TEST ##########\n\n")
cat("P Value slopes: ", p_value_slope, "\n")
cat("P Value intercept: ", p_value_intercept, "\n\n")
