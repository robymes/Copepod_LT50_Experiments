# Script to verify assumptions for t-test on non-linear model parameters
# Generic version for comparing any two experimental conditions
# Checks model residuals, convergence, and parameter validity for thermal tolerance analysis

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("gridExtra")) {
  install.packages("gridExtra")
  library(gridExtra)
}

if (!require("car")) {
  install.packages("car")
  library(car)
}

if (!require("boot")) {
  install.packages("boot")
  library(boot)
}

# Function to verify assumptions for non-linear model parameter t-test
# Now generic for any two experimental conditions
verify_nonlinear_ttest_assumptions <- function(data_path, 
                                               file_group1, 
                                               file_group2, 
                                               group1_name = "Group 1", 
                                               group2_name = "Group 2",
                                               nls_param_list = list(list(100, 30, 4))) {
  
  cat("################################################\n")
  cat("### NON-LINEAR MODEL T-TEST ASSUMPTIONS ###\n")
  cat("### Folder:", data_path, "###\n")
  cat("### Comparing:", group1_name, "vs", group2_name, "###\n")
  cat("################################################\n\n")
  
  # Complete file paths
  path_group1 <- file.path(data_path, file_group1)
  path_group2 <- file.path(data_path, file_group2)
  
  # Check file existence
  if (!file.exists(path_group1)) {
    stop("File not found: ", path_group1)
  }
  if (!file.exists(path_group2)) {
    stop("File not found: ", path_group2)
  }
  
  cat("Files to analyze:\n")
  cat("- ", group1_name, ":", file_group1, "\n")
  cat("- ", group2_name, ":", file_group2, "\n\n")
  
  cat("### 1. FITTING NON-LINEAR MODELS ###\n")
  
  # Fit models for both conditions
  model_results_group1 <- fit_survival_model(path_group1, nls_param_list, group1_name)
  model_results_group2 <- fit_survival_model(path_group2, nls_param_list, group2_name)
  
  cat("\n### 2. MODEL CONVERGENCE AND FIT QUALITY ###\n")
  
  # Check convergence
  convergence_group1 <- model_results_group1$convergence
  convergence_group2 <- model_results_group2$convergence
  
  cat(group1_name, "model convergence:", ifelse(convergence_group1, "✓ SUCCESS", "✗ FAILED"), "\n")
  cat(group2_name, "model convergence:", ifelse(convergence_group2, "✓ SUCCESS", "✗ FAILED"), "\n")
  
  if (!convergence_group1 || !convergence_group2) {
    cat("⚠ WARNING: Model convergence failed. Results may be unreliable.\n")
  }
  
  # R-squared values
  cat(group1_name, "R²:", round(model_results_group1$r_squared, 4), "\n")
  cat(group2_name, "R²:", round(model_results_group2$r_squared, 4), "\n\n")
  
  cat("### 3. RESIDUAL ANALYSIS ###\n")
  
  # Analyze residuals for both models
  analyze_residuals(model_results_group1, group1_name)
  analyze_residuals(model_results_group2, group2_name)
  
  cat("\n### 4. PARAMETER ESTIMATES AND UNCERTAINTY ###\n")
  
  # Extract parameter estimates and standard errors for p2 (LD50)
  p2_group1 <- model_results_group1$parameters[2, "Estimate"]
  p2_se_group1 <- model_results_group1$parameters[2, "Std. Error"]
  
  p2_group2 <- model_results_group2$parameters[2, "Estimate"]  
  p2_se_group2 <- model_results_group2$parameters[2, "Std. Error"]
  
  cat("Parameter p2 (LD50 - Temperature at 50% Survival):\n")
  cat("  ", group1_name, ": ", round(p2_group1, 4), "±", round(p2_se_group1, 4), "°C\n")
  cat("  ", group2_name, ": ", round(p2_group2, 4), "±", round(p2_se_group2, 4), "°C\n")
  
  # Calculate confidence intervals
  df_group1 <- model_results_group1$df
  df_group2 <- model_results_group2$df
  
  ci_group1 <- p2_group1 + c(-1, 1) * qt(0.975, df_group1) * p2_se_group1
  ci_group2 <- p2_group2 + c(-1, 1) * qt(0.975, df_group2) * p2_se_group2
  
  cat("  95% CI", group1_name, ": [", round(ci_group1[1], 4), ",", round(ci_group1[2], 4), "] °C\n")
  cat("  95% CI", group2_name, ": [", round(ci_group2[1], 4), ",", round(ci_group2[2], 4), "] °C\n\n")
  
  cat("### 5. T-TEST ON PARAMETERS (ORIGINAL METHOD) ###\n")
  
  # Perform the original t-test as in the source code (now for p2 - LD50)
  t_statistic <- (p2_group1 - p2_group2) / sqrt(p2_se_group1^2 + p2_se_group2^2)
  min_df <- min(df_group1, df_group2)
  p_value_original <- 2 * (1 - pt(abs(t_statistic), min_df))
  
  cat("Original t-test on LD50 (Welch-Satterthwaite approximation):\n")
  cat("  t =", round(t_statistic, 4), "\n")
  cat("  df =", min_df, "\n")
  cat("  p-value =", format(p_value_original, scientific = TRUE), "\n")
  
  # Interpretation considering that higher LD50 = greater heat tolerance
  if (p_value_original < 0.05) {
    if (t_statistic > 0) {
      interpretation <- paste(group1_name, "shows significantly higher thermal tolerance than", group2_name)
    } else {
      interpretation <- paste(group2_name, "shows significantly higher thermal tolerance than", group1_name)
    }
  } else {
    interpretation <- "No statistically significant difference in thermal tolerance"
  }
  
  cat("  Interpretation:", interpretation, "\n\n")
  
  cat("### 6. BOOTSTRAP CONFIDENCE INTERVALS ###\n")
  
  # Perform bootstrap analysis
  bootstrap_results <- perform_bootstrap_analysis(
    path_group1, path_group2, nls_param_list, group1_name, group2_name
  )
  
  cat("Bootstrap results (1000 iterations):\n")
  cat("  Difference in LD50 (", group1_name, "-", group2_name, "): ", round(bootstrap_results$difference_mean, 4), "°C\n")
  cat("  95% Bootstrap CI: [", round(bootstrap_results$ci[1], 4), ",", 
      round(bootstrap_results$ci[2], 4), "] °C\n")
  cat("  Bootstrap p-value: ", format(bootstrap_results$p_value, scientific = TRUE), "\n")
  
  # Interpretation considering that positive difference means group1 has higher thermal tolerance
  if (bootstrap_results$p_value < 0.05) {
    if (bootstrap_results$difference_mean > 0) {
      interpretation <- paste("Statistically significant difference:", group1_name, "shows higher thermal tolerance")
    } else {
      interpretation <- paste("Statistically significant difference:", group2_name, "shows higher thermal tolerance")
    }
  } else {
    interpretation <- "No statistically significant difference in thermal tolerance (Bootstrap)"
  }
  
  cat("  Interpretation:", interpretation, "\n\n")
  
  cat("### 7. MODEL COMPARISON AND RECOMMENDATIONS ###\n")
  
  # Overall assessment
  residuals_ok_group1 <- model_results_group1$residuals_normal && model_results_group1$residuals_homoscedastic
  residuals_ok_group2 <- model_results_group2$residuals_normal && model_results_group2$residuals_homoscedastic
  
  cat("Model validity assessment:\n")
  cat("  ", group1_name, "model: ", ifelse(residuals_ok_group1, "✓ VALID", "⚠ QUESTIONABLE"), "\n")
  cat("  ", group2_name, "model: ", ifelse(residuals_ok_group2, "✓ VALID", "⚠ QUESTIONABLE"), "\n")
  
  if (residuals_ok_group1 && residuals_ok_group2) {
    cat("\n✓ RECOMMENDATION: Both parametric and bootstrap methods are valid.\n")
    cat("  The original t-test approach is statistically sound.\n")
  } else {
    cat("\n⚠ RECOMMENDATION: Use bootstrap confidence intervals.\n")
    cat("  Model assumptions may be violated - bootstrap provides more robust inference.\n")
  }
  
  # Create diagnostic plots
  cat("\n### 8. CREATING DIAGNOSTIC PLOTS ###\n")
  create_model_diagnostic_plots(model_results_group1, model_results_group2, 
                                bootstrap_results, data_path, group1_name, group2_name)
  
  # Return comprehensive results
  return(list(
    model_group1 = model_results_group1,
    model_group2 = model_results_group2,
    original_ttest = list(
      t_statistic = t_statistic,
      p_value = p_value_original,
      df = min_df
    ),
    bootstrap_results = bootstrap_results,
    recommendation = ifelse(residuals_ok_group1 && residuals_ok_group2, 
                           "parametric", "bootstrap"),
    comparison = paste(group1_name, "vs", group2_name)
  ))
}

# Function to fit survival model to data
fit_survival_model <- function(file_path, nls_param_list, condition_name) {
  
  cat("Fitting model for", condition_name, "...\n")
  
  # Read data
  data <- read.csv(file_path, row.names = NULL)
  data$Survival <- data$Alive / (data$Alive + data$Dead)
  
  # Extract initial parameters
  p1_start <- nls_param_list[[1]][[1]]
  p2_start <- nls_param_list[[1]][[2]]
  p3_start <- nls_param_list[[1]][[3]]
  
  xdata <- data$Temperature
  ydata <- data$Survival
  
  # Fit non-linear model
  convergence <- TRUE
  model <- NULL
  
  tryCatch({
    model <- nls(ydata ~ p1 / (1 + (xdata / p2)^p3),
                 start = list(p1 = p1_start, p2 = p2_start, p3 = p3_start),
                 control = nls.control(maxiter = 500))
  }, error = function(e) {
    convergence <<- FALSE
    cat("  Error in model fitting:", e$message, "\n")
  })
  
  if (!convergence) {
    return(list(convergence = FALSE))
  }
  
  # Calculate R-squared
  residuals_model <- residuals(model)
  rss <- sum(residuals_model^2)
  tss <- sum((ydata - mean(ydata))^2)
  r_squared <- 1 - (rss / tss)
  
  # Analyze residuals
  fitted_values <- fitted(model)
  standardized_residuals <- residuals_model / sqrt(sum(residuals_model^2) / df.residual(model))
  
  # Test residual normality
  shapiro_test <- shapiro.test(standardized_residuals)
  residuals_normal <- shapiro_test$p.value > 0.05
  
  # Test homoscedasticity (Breusch-Pagan test approximation)
  bp_test <- lm(residuals_model^2 ~ fitted_values)
  bp_p_value <- summary(bp_test)$fstatistic
  if (!is.null(bp_p_value)) {
    bp_p_value <- pf(bp_p_value[1], bp_p_value[2], bp_p_value[3], lower.tail = FALSE)
    residuals_homoscedastic <- bp_p_value > 0.05
  } else {
    residuals_homoscedastic <- TRUE
  }
  
  cat("  Convergence: ✓\n")
  cat("  R²:", round(r_squared, 4), "\n")
  cat("  Residuals normal:", ifelse(residuals_normal, "✓", "✗"), 
      "(p =", round(shapiro_test$p.value, 4), ")\n")
  cat("  Residuals homoscedastic:", ifelse(residuals_homoscedastic, "✓", "✗"), "\n")
  
  return(list(
    model = model,
    convergence = TRUE,
    parameters = summary(model)$parameters,
    r_squared = r_squared,
    residuals = residuals_model,
    fitted = fitted_values,
    standardized_residuals = standardized_residuals,
    df = df.residual(model),
    data = data,
    residuals_normal = residuals_normal,
    residuals_homoscedastic = residuals_homoscedastic,
    shapiro_p = shapiro_test$p.value,
    condition = condition_name
  ))
}

# Function to analyze residuals
analyze_residuals <- function(model_results, condition_name) {
  
  if (!model_results$convergence) {
    cat(condition_name, "model failed to converge - skipping residual analysis\n")
    return()
  }
  
  cat(condition_name, "residual analysis:\n")
  cat("  Normality test (Shapiro-Wilk): p =", round(model_results$shapiro_p, 4), "\n")
  cat("  Interpretation:", ifelse(model_results$residuals_normal,
                                  "Residuals are normally distributed",
                                  "Residuals are NOT normally distributed"), "\n")
  cat("  Homoscedasticity:", ifelse(model_results$residuals_homoscedastic,
                                   "✓ Constant variance",
                                   "✗ Non-constant variance"), "\n")
}

# Function to perform bootstrap analysis
perform_bootstrap_analysis <- function(path_group1, path_group2, nls_param_list, 
                                       group1_name, group2_name, n_bootstrap = 1000) {
  
  cat("Performing bootstrap analysis (", n_bootstrap, "iterations)...\n")
  
  # Read original data
  data_group1 <- read.csv(path_group1, row.names = NULL)
  data_group2 <- read.csv(path_group2, row.names = NULL)
  
  # Bootstrap function
  bootstrap_difference <- function(data1, data2, indices1, indices2) {
    
    # Resample data
    boot_data1 <- data1[indices1, ]
    boot_data2 <- data2[indices2, ]
    
    # Calculate survival proportions
    boot_data1$Survival <- boot_data1$Alive / (boot_data1$Alive + boot_data1$Dead)
    boot_data2$Survival <- boot_data2$Alive / (boot_data2$Alive + boot_data2$Dead)
    
    # Fit models
    tryCatch({
      model1 <- nls(Survival ~ p1 / (1 + (Temperature / p2)^p3),
                    data = boot_data1,
                    start = list(p1 = nls_param_list[[1]][[1]], 
                                p2 = nls_param_list[[1]][[2]], 
                                p3 = nls_param_list[[1]][[3]]),
                    control = nls.control(maxiter = 200))
      
      model2 <- nls(Survival ~ p1 / (1 + (Temperature / p2)^p3),
                    data = boot_data2,
                    start = list(p1 = nls_param_list[[1]][[1]], 
                                p2 = nls_param_list[[1]][[2]], 
                                p3 = nls_param_list[[1]][[3]]),
                    control = nls.control(maxiter = 200))
      
      # Return difference in p1 parameters
      return(coef(model1)[1] - coef(model2)[1])
      
    }, error = function(e) {
      return(NA)
    })
  }
  
  # Perform bootstrap
  n_obs1 <- nrow(data_group1)
  n_obs2 <- nrow(data_group2)
  
  bootstrap_diffs <- replicate(n_bootstrap, {
    indices1 <- sample(1:n_obs1, n_obs1, replace = TRUE)
    indices2 <- sample(1:n_obs2, n_obs2, replace = TRUE)
    bootstrap_difference(data_group1, data_group2, indices1, indices2)
  })
  
  # Remove failed iterations
  bootstrap_diffs <- bootstrap_diffs[!is.na(bootstrap_diffs)]
  
  if (length(bootstrap_diffs) < n_bootstrap * 0.5) {
    cat("  Warning: Many bootstrap iterations failed (", 
        length(bootstrap_diffs), "successful out of", n_bootstrap, ")\n")
  }
  
  # Calculate statistics
  difference_mean <- mean(bootstrap_diffs)
  ci <- quantile(bootstrap_diffs, c(0.025, 0.975))
  
  # Bootstrap p-value (two-tailed)
  p_value_bootstrap <- 2 * min(
    mean(bootstrap_diffs >= 0),
    mean(bootstrap_diffs <= 0)
  )
  
  cat("  Successful iterations:", length(bootstrap_diffs), "\n")
  
  return(list(
    differences = bootstrap_diffs,
    difference_mean = difference_mean,
    ci = ci,
    p_value = p_value_bootstrap,
    comparison = paste(group1_name, "vs", group2_name)
  ))
}

# Function to create diagnostic plots
create_model_diagnostic_plots <- function(model_group1, model_group2, bootstrap_results, 
                                          data_path, group1_name, group2_name) {
  
  if (!model_group1$convergence || !model_group2$convergence) {
    cat("Cannot create plots - model convergence failed\n")
    return()
  }
  
  # 1. Residual plots
  residual_plot1 <- ggplot(data = data.frame(
    fitted = model_group1$fitted,
    residuals = model_group1$standardized_residuals
  ), aes(x = fitted, y = residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(se = FALSE, color = "blue") +
    labs(title = paste("Residuals vs Fitted -", group1_name),
         x = "Fitted Values", y = "Standardized Residuals") +
    theme_minimal()
  
  residual_plot2 <- ggplot(data = data.frame(
    fitted = model_group2$fitted,
    residuals = model_group2$standardized_residuals
  ), aes(x = fitted, y = residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(se = FALSE, color = "blue") +
    labs(title = paste("Residuals vs Fitted -", group2_name),
         x = "Fitted Values", y = "Standardized Residuals") +
    theme_minimal()
  
  # 2. Q-Q plots for residuals
  qq_plot1 <- ggplot(data = data.frame(residuals = model_group1$standardized_residuals),
                     aes(sample = residuals)) +
    stat_qq() + stat_qq_line(color = "red") +
    labs(title = paste("Q-Q Plot Residuals -", group1_name)) +
    theme_minimal()
  
  qq_plot2 <- ggplot(data = data.frame(residuals = model_group2$standardized_residuals),
                     aes(sample = residuals)) +
    stat_qq() + stat_qq_line(color = "red") +
    labs(title = paste("Q-Q Plot Residuals -", group2_name)) +
    theme_minimal()
  
  # 3. Model fit plots
  fit_plot1 <- ggplot(model_group1$data, aes(x = Temperature, y = Survival)) +
    geom_point() +
    geom_line(aes(y = model_group1$fitted), color = "red", size = 1) +
    labs(title = paste("Model Fit -", group1_name),
         x = "Temperature (°C)", y = "Survival Proportion") +
    theme_minimal()
  
  fit_plot2 <- ggplot(model_group2$data, aes(x = Temperature, y = Survival)) +
    geom_point() +
    geom_line(aes(y = model_group2$fitted), color = "red", size = 1) +
    labs(title = paste("Model Fit -", group2_name), 
         x = "Temperature (°C)", y = "Survival Proportion") +
    theme_minimal()
  
  # 4. Bootstrap distribution
  bootstrap_plot <- ggplot(data = data.frame(diff = bootstrap_results$differences),
                          aes(x = diff)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7, fill = "skyblue") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = bootstrap_results$ci[1], color = "blue", linetype = "dotted") +
    geom_vline(xintercept = bootstrap_results$ci[2], color = "blue", linetype = "dotted") +
    labs(title = paste("Bootstrap Distribution:", group1_name, "-", group2_name),
         x = paste("Difference in p1 (", group1_name, "-", group2_name, ")"), y = "Density") +
    theme_minimal()
  
  # Combine plots
  combined_plots <- grid.arrange(
    residual_plot1, residual_plot2,
    qq_plot1, qq_plot2,
    fit_plot1, fit_plot2,
    bootstrap_plot,
    ncol = 2, nrow = 4,
    heights = c(1, 1, 1, 1)
  )
  
  # Save plots
  condition_name <- basename(data_path)
  comparison_name <- paste(gsub(" ", "", group1_name), "vs", gsub(" ", "", group2_name), sep = "_")
  ggsave(
    filename = paste0("model_diagnostics_", gsub(" ", "_", condition_name), "_", comparison_name, ".png"),
    plot = combined_plots,
    width = 12, height = 16, dpi = 300
  )
  
  cat("Model diagnostic plots saved as 'model_diagnostics_", 
      gsub(" ", "_", condition_name), "_", comparison_name, ".png'\n\n")
}

# USAGE EXAMPLES:

# Example 2: Pressure comparison (original case)
data_path <- "data/cruise/t_test/pompeii_worms/oxic"
results_pressure <- verify_nonlinear_ttest_assumptions(
  data_path = data_path,
  file_group1 = "Pressure_Data02.csv",
  file_group2 = "No_Pressure_Data02.csv",
  group1_name = "Pressure", 
  group2_name = "No Pressure"
)