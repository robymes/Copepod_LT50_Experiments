if (!require("MASS")) {
  install.packages("MASS")
  library(MASS)
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("Cairo")) {
  install.packages("Cairo")
  library(Cairo)
}

source("plot_utilities.R")

non_linear_regression_file_func <- function(dir_path, csv_file, nls_param_list, time_list, lines_color, k, i) {
  csv_path <- sprintf("%s/%s", dir_path, csv_file)
  cat("########## FILE: ", csv_path, " ###########\n\n")
  data <- read.csv(csv_path, row.names = NULL)
  data$Survival <- data$Alive / (data$Alive + data$Dead)

  p1_start <- nls_param_list[[k]][[1]]
  p2_start <- nls_param_list[[k]][[2]]
  p3_start <- nls_param_list[[k]][[3]]
  xdata <- data[, 1]
  ydata <- data[, 4]

  fit <- nls(ydata ~ p1 / (1 + (xdata / p2)^p3),
             start = list(p1 = p1_start, p2 = p2_start, p3 = p3_start),
             control = nls.control(maxiter = 500))

  new <- data.frame(xdata = seq(min(xdata), max(xdata), len = 200))
  fit_params <- summary(fit)
  print(summary(fit))
  cat("\n")

  ld50 <- c(fit_params$parameters[2])

  rss <- sum(residuals(fit)^2)
  tss <- sum((ydata - mean(ydata))^2)
  r_squared <- 1 - (rss / tss)
  cat("R squared for fit model: ", r_squared, "\n\n")

  model_dose <- fit_params$parameters[2]
  model_dose_y <- predict(fit, newdata = data.frame(xdata = model_dose))
  model_dose_stderr <- fit_params$parameters[2, 2]

  # Extract the left part of the csv file name in order to verify if data has to be splitted by this part of the name
  left_part_csv_file <- ""
  if (grepl("_Data", csv_file)) {
    left_part_csv_file <- sub("(.*)_Data.*", "\\1", csv_file)
    left_part_csv_file <- gsub("_", " ", left_part_csv_file)
  }

  return(list(
    ld50 = ld50,
    fit_params = fit_params,
    df = df.residual(fit),
    new_data = new,
    fit_pred = predict(fit, newdata = new),
    time = ifelse(
      left_part_csv_file == "",
      paste(sprintf("%02d", time_list[i]), "h exposure"),
      paste(left_part_csv_file, "-", sprintf("%02d", time_list[i]), "h exposure")
    ),
    jitter_xdata = jitter(xdata, factor = 0.5),
    ydata = ydata,
    model_dose = model_dose,
    model_dose_y = model_dose_y,
    model_dose_stderr = model_dose_stderr
  ))
}

non_linear_regression_dir_func <- function(
  dir_path,
  nls_param_list,
  k,
  main_title,
  chart_subtitle,
  csv_file_list = NULL
) {
  cat("########## NON LINEAR REGRESSION ##########\n\n")
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")

  if (is.null(csv_file_list)) {
    csv_files <- list.files(path = dir_path)
  } else {
    csv_files <- csv_file_list
  }

  time_character <- c(regmatches(csv_files, regexpr("[0-9]*[0-9]", csv_files)))
  time_list <- as.numeric(time_character)
  ld50 <- c()
  nls_coefficients <- list()
  min_df <- .Machine$double.xmax
  survival_params <- c(0, 0, 0, 0)

  lines_color <- c("green", "red", "blue", "orange", "purple")
  i <- 0

  data_list <- list()

  for (csv_file in csv_files) {
    i <- i + 1
    non_linear_regression_file_result <- non_linear_regression_file_func(
      dir_path = dir_path,
      csv_file = csv_file,
      nls_param_list = nls_param_list,
      time_list = time_list,
      k = k,
      i = i
    )
    ld50 <- c(ld50, non_linear_regression_file_result$ld50)
    nls_coefficients <- append(nls_coefficients, coef(non_linear_regression_file_result$fit_params))
    if (non_linear_regression_file_result$df < min_df) {
      min_df <- non_linear_regression_file_result$df
    }
    if (i == 1) {
      survival_params <- c(
        non_linear_regression_file_result$fit_params$parameters[1],
        non_linear_regression_file_result$fit_params$parameters[2],
        non_linear_regression_file_result$fit_params$parameters[3],
        time_list[1]
      )
    }

    data_list[[i]] <- data.frame(
      xdata = non_linear_regression_file_result$new_data$xdata,
      ydata = non_linear_regression_file_result$fit_pred,
      time = non_linear_regression_file_result$time
    )
  }

  plot_data <- do.call(rbind, data_list)

  # Create a plot adding the non-linear regression lines
  non_linear_regression_plot <- ggplot(plot_data, aes(x = xdata, y = ydata, color = as.factor(time))) +
    geom_line() +
    labs(title = sprintf("LT50 Survival Curve\n%s", chart_subtitle),
         x = "Temperature (°C)", y = "Proportional Survival") +
    scale_color_manual(values = lines_color, name = "Exposure Time (h)")

  for (i in 1:length(data_list)) {
    non_linear_regression_file_result <- non_linear_regression_file_func(
      dir_path = dir_path,
      csv_file = csv_files[i],
      nls_param_list = nls_param_list,
      time_list = time_list,
      k = k,
      i = i
    )

    # Add to plot original sample data points
    non_linear_regression_plot <- non_linear_regression_plot + geom_point(data = data.frame(
      x = non_linear_regression_file_result$jitter_xdata,
      y = non_linear_regression_file_result$ydata,
      time = non_linear_regression_file_result$time
    ), aes(x = x, y = y, color = as.factor(time)), size = 0.5)

    # Add to plot model dose points
    non_linear_regression_plot <- non_linear_regression_plot + geom_point(data = data.frame(
      x = non_linear_regression_file_result$model_dose,
      y = non_linear_regression_file_result$model_dose_y,
      time = non_linear_regression_file_result$time
    ), aes(x = x, y = y, color = as.factor(time)), shape = 18, size = 4,
    inherit.aes = FALSE)

    # Add to plot std. error bars
    non_linear_regression_plot <- non_linear_regression_plot + geom_errorbarh(data = data.frame(
      xmin = non_linear_regression_file_result$model_dose - non_linear_regression_file_result$model_dose_stderr,
      xmax = non_linear_regression_file_result$model_dose + non_linear_regression_file_result$model_dose_stderr,
      y = non_linear_regression_file_result$model_dose_y
    ), aes(xmin = xmin, xmax = xmax, y = y), height = 0.02,
    inherit.aes = FALSE)

    # Add to plot model dose temperature values
    non_linear_regression_plot <- non_linear_regression_plot + geom_text(data = data.frame(
      x = non_linear_regression_file_result$model_dose,
      y = non_linear_regression_file_result$model_dose_y,
      label = sprintf("%.1f", non_linear_regression_file_result$model_dose)
    ), aes(x = x, y = y, label = label),
    position = position_nudge(y = -0.05),
    size = 2,
    inherit.aes = FALSE)
  }

  print(non_linear_regression_plot)

  save_plot_func(
    plot = non_linear_regression_plot,
    path = paste("charts/", gsub("\\\\", "/", gsub(" ", "", tools::toTitleCase(trimws(main_title)))), sep = ""),
    filename = paste(gsub(" ", "", trimws(chart_subtitle)), "_lt50", sep = ""),
    width = 1920,
    height = 1080
  )

  return(list(
    ld50 = ld50,
    time_list = time_list,
    nls_coefficients = nls_coefficients,
    min_df = min_df,
    survival_params = survival_params
  ))
}