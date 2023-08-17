if (!require("MASS")) {
  install.packages("MASS")
  library(MASS)
}

if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

non_linear_regression_file_func <- function(dir_path, csv_file, nls_param_list, time_list, lines_color, k, i) {
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

  # Add the model_dose value to the LD50 list
  ld50 <- c(fit_params$parameters[2])

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

  # R squared calculation for fit model
  rss <- sum(residuals(fit)^2)
  tss <- sum((ydata - mean(ydata))^2)
  r_squared <- 1 - (rss / tss)
  cat("R squared for fit model: ", r_squared, "\n\n")

  return(list(
    ld50 = ld50,
    fit_params = fit_params,
    df = df.residual(fit)
  ))
}

non_linear_regression_dir_func <- function(dir_path, nls_param_list, k, chart_subtitle, csv_file_list = NULL) {
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")

  if (is.null(csv_file_list)) {
    # Get a list of all files in the directory
    csv_files <- list.files(path = dir_path)
  } else {
    csv_files <- csv_file_list
  }

  # Exposure time is in the file name: it must be extracted and converted to a numeric value
  time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
  time_list <- as.numeric(time_character)
  ld50 <- c()
  nls_coefficients <- list()
  min_df <- .Machine$double.xmax
  survival_params <- c(0, 0, 0, 0)

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
    main = paste("LT50 Survival Curve\n", chart_subtitle)
  )

  # Assign a color to each data set
  lines_color <- c("green", "red", "blue", "orange", "purple")

  # This is used to number things in the loop
  i <- 0

  # Iterate over csv files and create plot
  for (csv_file in csv_files) {
    i <- i + 1
    non_linear_regression_file_result <- non_linear_regression_file_func(
      dir_path = dir_path,
      csv_file = csv_file,
      nls_param_list = nls_param_list,
      time_list = time_list,
      lines_color = lines_color,
      k = k,
      i = i
    )
    ld50 <- c(ld50, non_linear_regression_file_result$ld50)
    nls_coefficients <- append(nls_coefficients, coef(non_linear_regression_file_result$fit_params))
    if (non_linear_regression_file_result$df < min_df) {
      min_df <- non_linear_regression_file_result$df
    }
    # In the 3D plot the parameters of the shortest exposure are used in the function
    if (i == 1) {
      survival_params <- c(
        non_linear_regression_file_result$fit_params$parameters[1],
        non_linear_regression_file_result$fit_params$parameters[2],
        non_linear_regression_file_result$fit_params$parameters[3],
        time_list[1]
      )
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

  return(list(
    ld50 = ld50,
    time_list = time_list,
    nls_coefficients = nls_coefficients,
    min_df = min_df,
    survival_params = survival_params
  ))
}