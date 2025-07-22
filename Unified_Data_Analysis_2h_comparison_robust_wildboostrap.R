# WILD BOOTSTRAP ANALYSIS WITH COMPREHENSIVE DIAGNOSTICS
# Optimized for thermal tolerance NLS models with assumption testing

if (!require("car")) {
  install.packages("car")
  library(car)
}

if (!require("lmtest")) {
  install.packages("lmtest")
  library(lmtest)
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("gridExtra")) {
  install.packages("gridExtra")
  library(gridExtra)
}

# =====================================================================
# MAIN FUNCTION: WILD BOOTSTRAP WITH COMPREHENSIVE DIAGNOSTICS
# =====================================================================

wild_bootstrap_analysis <- function(data_path, file_group1, file_group2,
                                   group1_name = "Group 1", group2_name = "Group 2",
                                   n_bootstrap = 2000,
                                   nls_start_params = list(100, 30, 4),
                                   alpha = 0.05) {
  
  cat("##################################################\n")
  cat("### WILD BOOTSTRAP ANALYSIS WITH DIAGNOSTICS ###\n")
  cat("### FOR THERMAL TOLERANCE NLS MODELS ###\n")
  cat("##################################################\n\n")
  
  cat("Folder:", data_path, "\n")
  cat("Comparing:", group1_name, "vs", group2_name, "\n")
  cat("Bootstrap iterations:", n_bootstrap, "\n")
  cat("Significance level:", alpha, "\n\n")
  
  # Read and prepare data
  cat("### 1. LOADING AND PREPARING DATA ###\n")
  
  path_group1 <- file.path(data_path, file_group1)
  path_group2 <- file.path(data_path, file_group2)
  
  if (!file.exists(path_group1) || !file.exists(path_group2)) {
    stop("One or both data files not found. Check file paths.")
  }
  
  data1 <- read.csv(path_group1, row.names = NULL)
  data2 <- read.csv(path_group2, row.names = NULL)
  
  # Calculate survival proportions
  data1$Survival <- data1$Alive / (data1$Alive + data1$Dead)
  data2$Survival <- data2$Alive / (data2$Alive + data2$Dead)
  
  cat("Data loaded successfully:\n")
  cat(" ", group1_name, ":", nrow(data1), "observations\n")
  cat(" ", group2_name, ":", nrow(data2), "observations\n\n")
  
  # Fit original models
  cat("### 2. FITTING ORIGINAL NLS MODELS ###\n")
  
  model1 <- fit_nls_safe(data1, nls_start_params, group1_name)
  model2 <- fit_nls_safe(data2, nls_start_params, group2_name)
  
  if (!model1$converged || !model2$converged) {
    stop("One or both original models failed to converge. Check data and starting parameters.")
  }
  
  # Extract LD50 values
  ld50_group1 <- coef(model1$model)[2]
  ld50_group2 <- coef(model2$model)[2]
  original_difference <- ld50_group1 - ld50_group2
  
  cat("Model fitting successful:\n")
  cat(" ", group1_name, "LD50:", round(ld50_group1, 2), "°C\n")
  cat(" ", group2_name, "LD50:", round(ld50_group2, 2), "°C\n")
  cat("  Original difference:", round(original_difference, 3), "°C\n\n")
  
  # Model quality assessment
  r2_group1 <- get_r_squared(model1$model, data1$Survival)
  r2_group2 <- get_r_squared(model2$model, data2$Survival)
  
  cat("Model fit quality:\n")
  cat(" ", group1_name, "R²:", round(r2_group1, 4), "\n")
  cat(" ", group2_name, "R²:", round(r2_group2, 4), "\n\n")
  
  # Comprehensive residual diagnostics
  cat("### 3. RESIDUAL DIAGNOSTICS AND ASSUMPTION TESTING ###\n")
  
  diagnostics1 <- comprehensive_residual_diagnostics(model1$model, data1, group1_name)
  diagnostics2 <- comprehensive_residual_diagnostics(model2$model, data2, group2_name)
  
  # Method recommendation based on diagnostics
  wild_bootstrap_recommended <- recommend_wild_bootstrap(diagnostics1, diagnostics2)
  
  cat("\n### 4. METHOD RECOMMENDATION ###\n")
  cat("====================================\n")
  
  if (wild_bootstrap_recommended$use_wild_bootstrap) {
    cat("✅ RECOMMENDATION: Wild Bootstrap is APPROPRIATE\n")
    cat("Reasons:\n")
    for (reason in wild_bootstrap_recommended$reasons) {
      cat(" • ", reason, "\n")
    }
  } else {
    cat("⚠️  RECOMMENDATION: Wild Bootstrap may not be optimal\n")
    cat("Issues detected:\n")
    for (issue in wild_bootstrap_recommended$issues) {
      cat(" • ", issue, "\n")
    }
    cat("\nProceeding with Wild Bootstrap but interpret results cautiously.\n")
  }
  
  cat("\n### 5. WILD BOOTSTRAP ANALYSIS ###\n")
  cat("====================================\n")
  
  # Perform wild bootstrap
  cat("Running", n_bootstrap, "wild bootstrap iterations...\n")
  
  bootstrap_results <- wild_bootstrap_core(
    model1$model, model2$model, 
    data1, data2, 
    n_bootstrap, nls_start_params
  )
  
  cat("Bootstrap completed:\n")
  cat("  Successful iterations:", bootstrap_results$successful, "/", n_bootstrap,
      "(", round(100 * bootstrap_results$successful / n_bootstrap, 1), "%)\n")
  
  if (bootstrap_results$successful < n_bootstrap * 0.8) {
    cat("⚠️  WARNING: Low bootstrap success rate (<80%). Results may be unstable.\n")
  }
  
  # Calculate statistics
  boot_differences <- bootstrap_results$differences
  boot_mean <- mean(boot_differences, na.rm = TRUE)
  boot_ci <- quantile(boot_differences, c(alpha/2, 1-alpha/2), na.rm = TRUE)
  
  # Two-tailed p-value
  p_value <- 2 * min(
    mean(boot_differences >= 0, na.rm = TRUE),
    mean(boot_differences <= 0, na.rm = TRUE)
  )
  
  cat("\n### 6. BOOTSTRAP RESULTS ###\n")
  cat("=============================\n")
  cat("Original difference:", round(original_difference, 3), "°C\n")
  cat("Bootstrap mean:", round(boot_mean, 3), "°C\n")
  cat("Bootstrap SD:", round(sd(boot_differences, na.rm = TRUE), 3), "°C\n")
  cat("95% Confidence Interval: [", round(boot_ci[1], 3), ",", round(boot_ci[2], 3), "] °C\n")
  cat("P-value (two-tailed):", format(p_value, scientific = TRUE), "\n")
  
  # Statistical interpretation
  cat("\n### 7. STATISTICAL INTERPRETATION ###\n")
  cat("======================================\n")
  
  if (p_value < alpha) {
    if (original_difference > 0) {
      interpretation <- paste("Statistically significant difference:", group1_name, "shows higher thermal tolerance")
    } else {
      interpretation <- paste("Statistically significant difference:", group2_name, "shows higher thermal tolerance")
    }
  } else {
    interpretation <- "No statistically significant difference in thermal tolerance"
  }
  
  cat("Interpretation:", interpretation, "\n")
  
  # Effect size assessment (using raw residuals for NLS)
  residuals_1 <- residuals(model1$model)
  residuals_2 <- residuals(model2$model)
  pooled_sd <- sqrt(((nrow(data1)-1)*var(residuals_1) + 
                     (nrow(data2)-1)*var(residuals_2)) / 
                    (nrow(data1) + nrow(data2) - 2))
  effect_size <- abs(original_difference) / pooled_sd
  
  cat("Effect size (standardized):", round(effect_size, 3), "\n")
  cat("Effect magnitude:", 
      ifelse(effect_size < 0.2, "Negligible",
             ifelse(effect_size < 0.5, "Small",
                    ifelse(effect_size < 0.8, "Medium", "Large"))), "\n")
  
  # Biological interpretation
  cat("\n### 8. BIOLOGICAL INTERPRETATION ###\n")
  cat("=====================================\n")
  
  if (abs(original_difference) > 1) {
    biological_relevance <- "Potentially biologically relevant (>1°C difference)"
  } else {
    biological_relevance <- "Limited biological relevance (<1°C difference)"
  }
  
  cat("Magnitude assessment:", biological_relevance, "\n")
  
  if (p_value < alpha && abs(original_difference) > 1) {
    cat("Conclusion: Both statistically significant and biologically relevant effect\n")
  } else if (p_value < alpha) {
    cat("Conclusion: Statistically significant but small biological effect\n")
  } else if (abs(original_difference) > 1) {
    cat("Conclusion: Potentially relevant effect but not statistically significant\n")
  } else {
    cat("Conclusion: No significant or relevant effect detected\n")
  }
  
  # Generate diagnostic plots
  cat("\n### 9. GENERATING DIAGNOSTIC PLOTS ###\n")
  cat("=======================================\n")
  
  plots <- generate_diagnostic_plots(
    model1$model, model2$model, 
    data1, data2, 
    boot_differences, 
    group1_name, group2_name
  )
  
  cat("Diagnostic plots generated successfully.\n")
  
  # Compile final results
  final_results <- list(
    # Data and models
    data1 = data1,
    data2 = data2,
    model1 = model1$model,
    model2 = model2$model,
    
    # Original estimates
    ld50_group1 = ld50_group1,
    ld50_group2 = ld50_group2,
    original_difference = original_difference,
    
    # Model quality
    r2_group1 = r2_group1,
    r2_group2 = r2_group2,
    
    # Diagnostics
    diagnostics_group1 = diagnostics1,
    diagnostics_group2 = diagnostics2,
    recommendation = wild_bootstrap_recommended,
    
    # Bootstrap results
    bootstrap_differences = boot_differences,
    bootstrap_mean = boot_mean,
    bootstrap_ci = boot_ci,
    p_value = p_value,
    successful_iterations = bootstrap_results$successful,
    
    # Interpretation
    interpretation = interpretation,
    effect_size = effect_size,
    biological_relevance = biological_relevance,
    
    # Plots
    diagnostic_plots = plots,
    
    # Metadata
    method = "Wild Bootstrap",
    n_bootstrap = n_bootstrap,
    alpha = alpha,
    group_names = c(group1_name, group2_name)
  )
  
  return(final_results)
}

# =====================================================================
# CORE WILD BOOTSTRAP FUNCTION
# =====================================================================

wild_bootstrap_core <- function(model1, model2, data1, data2, n_bootstrap, start_params) {
  
  # Extract fitted values and residuals
  fitted1 <- fitted(model1)
  fitted2 <- fitted(model2)
  residuals1 <- residuals(model1)
  residuals2 <- residuals(model2)
  
  # Initialize storage
  differences <- numeric(n_bootstrap)
  successful <- 0
  
  # Progress reporting
  progress_interval <- max(1, floor(n_bootstrap / 10))
  
  for (i in 1:n_bootstrap) {
    
    # Generate Rademacher weights {-1, +1}
    weights1 <- sample(c(-1, 1), length(residuals1), replace = TRUE)
    weights2 <- sample(c(-1, 1), length(residuals2), replace = TRUE)
    
    # Create wild bootstrap responses
    y1_boot <- fitted1 + weights1 * residuals1
    y2_boot <- fitted2 + weights2 * residuals2
    
    # Constrain to [0,1] for survival proportions
    y1_boot <- pmax(0, pmin(1, y1_boot))
    y2_boot <- pmax(0, pmin(1, y2_boot))
    
    # Create bootstrap datasets
    data1_boot <- data.frame(Temperature = data1$Temperature, Survival = y1_boot)
    data2_boot <- data.frame(Temperature = data2$Temperature, Survival = y2_boot)
    
    # Fit bootstrap models
    boot_model1 <- fit_nls_safe(data1_boot, start_params, paste("Boot", i, "G1"), silent = TRUE)
    boot_model2 <- fit_nls_safe(data2_boot, start_params, paste("Boot", i, "G2"), silent = TRUE)
    
    # Store result if both converged
    if (boot_model1$converged && boot_model2$converged) {
      differences[i] <- coef(boot_model1$model)[2] - coef(boot_model2$model)[2]
      successful <- successful + 1
    } else {
      differences[i] <- NA
    }
    
    # Progress report
    if (i %% progress_interval == 0) {
      cat("Completed", i, "/", n_bootstrap, "iterations\n")
    }
  }
  
  return(list(
    differences = differences[!is.na(differences)],
    successful = successful
  ))
}

# =====================================================================
# COMPREHENSIVE RESIDUAL DIAGNOSTICS
# =====================================================================

comprehensive_residual_diagnostics <- function(model, data, group_name) {
  
  cat("\n--- DIAGNOSTICS FOR", group_name, "---\n")
  
  # Get residuals and fitted values for NLS models
  residuals_raw <- residuals(model)
  fitted_vals <- fitted(model)
  
  # Calculate standardized residuals manually for NLS
  # (rstandard doesn't work with nls objects)
  mse <- sum(residuals_raw^2) / df.residual(model)
  residuals_std <- residuals_raw / sqrt(mse)
  
  # Test 1: Normality (Shapiro-Wilk)
  shapiro_test <- shapiro.test(residuals_std)
  normality_ok <- shapiro_test$p.value > 0.05
  
  cat("1. Normality Test (Shapiro-Wilk):\n")
  cat("   W =", round(shapiro_test$statistic, 4), 
      ", p-value =", format(shapiro_test$p.value, scientific = TRUE), "\n")
  cat("   Result:", ifelse(normality_ok, "✅ Normal", "❌ Non-normal"), "\n")
  
  # Test 2: Heteroscedasticity (Breusch-Pagan)
  # Create auxiliary regression: residuals² ~ fitted
  aux_data <- data.frame(
    residuals_squared = residuals_std^2,
    fitted = fitted_vals
  )
  
  aux_model <- lm(residuals_squared ~ fitted, data = aux_data)
  bp_test <- lmtest::bptest(aux_model)
  homoscedastic_ok <- bp_test$p.value > 0.05
  
  cat("2. Heteroscedasticity Test (Breusch-Pagan):\n")
  cat("   LM =", round(bp_test$statistic, 4), 
      ", p-value =", format(bp_test$p.value, scientific = TRUE), "\n")
  cat("   Result:", ifelse(homoscedastic_ok, "✅ Homoscedastic", "❌ Heteroscedastic"), "\n")
  
  # Test 3: Outliers (based on standardized residuals)
  outliers <- abs(residuals_std) > 2.5  # Conservative threshold
  n_outliers <- sum(outliers)
  outliers_ok <- n_outliers <= max(1, round(0.05 * length(residuals_std)))
  
  cat("3. Outlier Assessment:\n")
  cat("   Extreme residuals (|z| > 2.5):", n_outliers, "/", length(residuals_std), "\n")
  cat("   Result:", ifelse(outliers_ok, "✅ No excessive outliers", "⚠️  Potential outliers"), "\n")
  
  # Test 4: Model fit assessment
  r_squared <- get_r_squared(model, data$Survival)
  fit_ok <- r_squared > 0.7  # Threshold for good fit
  
  cat("4. Model Fit Quality:\n")
  cat("   R² =", round(r_squared, 4), "\n")
  cat("   Result:", ifelse(fit_ok, "✅ Good fit", "⚠️  Moderate fit"), "\n")
  
  # Additional diagnostic: Residuals vs Temperature pattern
  temp_pattern <- cor(data$Temperature, abs(residuals_std))
  pattern_ok <- abs(temp_pattern) < 0.3  # Low correlation with temperature
  
  cat("5. Temperature-Residual Pattern:\n")
  cat("   |Correlation| =", round(abs(temp_pattern), 3), "\n")
  cat("   Result:", ifelse(pattern_ok, "✅ No systematic pattern", "⚠️  Temperature-dependent variance"), "\n")
  
  return(list(
    shapiro_test = shapiro_test,
    bp_test = bp_test,
    normality_ok = normality_ok,
    homoscedastic_ok = homoscedastic_ok,
    outliers_ok = outliers_ok,
    fit_ok = fit_ok,
    pattern_ok = pattern_ok,
    residuals_std = residuals_std,
    fitted_vals = fitted_vals,
    r_squared = r_squared,
    n_outliers = n_outliers,
    temp_correlation = temp_pattern
  ))
}

# =====================================================================
# WILD BOOTSTRAP RECOMMENDATION SYSTEM
# =====================================================================

recommend_wild_bootstrap <- function(diag1, diag2) {
  
  reasons <- character(0)
  issues <- character(0)
  use_wild_bootstrap <- TRUE
  
  # Check normality
  if (!diag1$normality_ok || !diag2$normality_ok) {
    reasons <- c(reasons, "Non-normal residuals detected")
  } else {
    issues <- c(issues, "Residuals appear normal (parametric methods might be sufficient)")
  }
  
  # Check heteroscedasticity
  if (!diag1$homoscedastic_ok || !diag2$homoscedastic_ok) {
    reasons <- c(reasons, "Heteroscedasticity detected")
  }
  
  # Check model fit
  if (!diag1$fit_ok || !diag2$fit_ok) {
    issues <- c(issues, "Poor model fit detected (R² < 0.7)")
  }
  
  # Check outliers
  if (!diag1$outliers_ok || !diag2$outliers_ok) {
    reasons <- c(reasons, "Potential outliers present")
  }
  
  # Check temperature patterns
  if (!diag1$pattern_ok || !diag2$pattern_ok) {
    reasons <- c(reasons, "Temperature-dependent residual patterns")
  }
  
  # Overall recommendation
  n_issues <- length(issues)
  n_reasons <- length(reasons)
  
  if (n_reasons >= 2) {
    reasons <- c(reasons, "Multiple assumption violations favor robust methods")
  }
  
  if (n_issues > n_reasons) {
    use_wild_bootstrap <- FALSE
    issues <- c(issues, "Consider alternative approaches or model re-specification")
  }
  
  if (length(reasons) == 0) {
    reasons <- c("Precautionary robustness for non-linear models")
  }
  
  return(list(
    use_wild_bootstrap = use_wild_bootstrap,
    reasons = reasons,
    issues = issues,
    confidence = ifelse(n_reasons >= 2, "High", ifelse(n_reasons >= 1, "Medium", "Low"))
  ))
}

# =====================================================================
# DIAGNOSTIC PLOT GENERATION
# =====================================================================

generate_diagnostic_plots <- function(model1, model2, data1, data2, 
                                     bootstrap_diffs, group1_name, group2_name) {
  
  # Extract residuals and fitted values for NLS models
  # Calculate standardized residuals manually
  residuals_raw1 <- residuals(model1)
  residuals_raw2 <- residuals(model2)
  fitted1 <- fitted(model1)
  fitted2 <- fitted(model2)
  
  # Standardize residuals manually for NLS
  mse1 <- sum(residuals_raw1^2) / df.residual(model1)
  mse2 <- sum(residuals_raw2^2) / df.residual(model2)
  residuals1 <- residuals_raw1 / sqrt(mse1)
  residuals2 <- residuals_raw2 / sqrt(mse2)
  
  # 1. Residuals vs Fitted
  plot_data_resid <- rbind(
    data.frame(fitted = fitted1, residuals = residuals1, group = group1_name),
    data.frame(fitted = fitted2, residuals = residuals2, group = group2_name)
  )
  
  p1 <- ggplot(plot_data_resid, aes(x = fitted, y = residuals, color = group)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_smooth(method = "loess", se = FALSE, linetype = "solid") +
    labs(title = "Residuals vs Fitted Values",
         x = "Fitted Values", y = "Standardized Residuals",
         color = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 2. Q-Q Plots for normality
  qq_data1 <- data.frame(
    theoretical = qnorm(ppoints(length(residuals1))),
    sample = sort(residuals1),
    group = group1_name
  )
  qq_data2 <- data.frame(
    theoretical = qnorm(ppoints(length(residuals2))),
    sample = sort(residuals2),
    group = group2_name
  )
  qq_data <- rbind(qq_data1, qq_data2)
  
  p2 <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~group) +
    labs(title = "Q-Q Plots for Normality",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # 3. Bootstrap distribution
  p3 <- ggplot(data.frame(diff = bootstrap_diffs), aes(x = diff)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7, fill = "skyblue") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = quantile(bootstrap_diffs, c(0.025, 0.975)), 
               color = "blue", linetype = "dotted", size = 1) +
    labs(title = paste("Bootstrap Distribution:", group1_name, "-", group2_name),
         x = paste("Difference in LD50 (°C)"), y = "Density") +
    theme_minimal()
  
  # 4. Model fits
  temp_range <- seq(min(c(data1$Temperature, data2$Temperature)), 
                    max(c(data1$Temperature, data2$Temperature)), 
                    length.out = 100)
  
  # Extract parameters for prediction
  params1 <- coef(model1)
  params2 <- coef(model2)
  
  pred1 <- params1[1] / (1 + (temp_range / params1[2])^params1[3])
  pred2 <- params2[1] / (1 + (temp_range / params2[2])^params2[3])
  
  fit_data <- rbind(
    data.frame(Temperature = data1$Temperature, Survival = data1$Survival, group = group1_name, type = "Observed"),
    data.frame(Temperature = data2$Temperature, Survival = data2$Survival, group = group2_name, type = "Observed"),
    data.frame(Temperature = temp_range, Survival = pred1, group = group1_name, type = "Fitted"),
    data.frame(Temperature = temp_range, Survival = pred2, group = group2_name, type = "Fitted")
  )
  
  p4 <- ggplot(fit_data, aes(x = Temperature, y = Survival, color = group)) +
    geom_point(data = subset(fit_data, type == "Observed"), alpha = 0.7, size = 2) +
    geom_line(data = subset(fit_data, type == "Fitted"), size = 1) +
    labs(title = "Model Fits to Data",
         x = "Temperature (°C)", y = "Survival Proportion",
         color = "Group") +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(
    residuals_vs_fitted = p1,
    qq_plots = p2,
    bootstrap_distribution = p3,
    model_fits = p4
  ))
}

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

fit_nls_safe <- function(data, start_params, name = "Model", silent = FALSE) {
  
  convergence <- TRUE
  model <- NULL
  
  tryCatch({
    model <- nls(Survival ~ p1 / (1 + (Temperature / p2)^p3),
                data = data,
                start = list(p1 = start_params[1], 
                            p2 = start_params[2], 
                            p3 = start_params[3]),
                control = nls.control(maxiter = 200, warnOnly = TRUE))
  }, error = function(e) {
    convergence <<- FALSE
    if (!silent) cat("Error fitting", name, ":", e$message, "\n")
  })
  
  return(list(
    model = model,
    converged = convergence,
    name = name
  ))
}

get_r_squared <- function(model, observed) {
  fitted_vals <- fitted(model)
  ss_res <- sum((observed - fitted_vals)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  return(1 - ss_res/ss_tot)
}

# =====================================================================
# USAGE EXAMPLES
# =====================================================================

# Example 1: Pompeii worms under oxic conditions
results <- wild_bootstrap_analysis(
  data_path = "data/cruise/t_test/tubeworms_mussels/oxic",
  file_group1 = "No_Pressure_Data02.csv",
  file_group2 = "Pressure_Data02.csv",
  group1_name = "No Pressure",
  group2_name = "Pressure",
  n_bootstrap = 2000
)

# Example 2: Pompeii worms under anoxic conditions  
# results_anoxic <- wild_bootstrap_analysis(
#   data_path = "data/cruise/t_test/pompeii_worms/anoxic",
#   file_group1 = "No_Pressure_Data02.csv",
#   file_group2 = "Pressure_Data02.csv",
#   group1_name = "No Pressure",
#   group2_name = "Pressure",
#   n_bootstrap = 2000
# )

# View diagnostic plots
# grid.arrange(results_oxic$diagnostic_plots$residuals_vs_fitted,
#              results_oxic$diagnostic_plots$qq_plots,
#              results_oxic$diagnostic_plots$bootstrap_distribution,
#              results_oxic$diagnostic_plots$model_fits,
#              ncol = 2)

# =====================================================================
# BATCH ANALYSIS FOR MULTIPLE COMPARISONS
# =====================================================================

batch_wild_bootstrap <- function(analysis_list) {
  
  cat("### BATCH WILD BOOTSTRAP ANALYSIS ###\n")
  cat("======================================\n\n")
  
  results_summary <- data.frame()
  all_results <- list()
  
  for (i in 1:length(analysis_list)) {
    
    analysis <- analysis_list[[i]]
    cat("Analysis", i, ":", analysis$description, "\n")
    cat("-------------------------------------------\n")
    
    # Run analysis
    result <- wild_bootstrap_analysis(
      data_path = analysis$data_path,
      file_group1 = analysis$file_group1,
      file_group2 = analysis$file_group2,
      group1_name = analysis$group1_name,
      group2_name = analysis$group2_name,
      n_bootstrap = analysis$n_bootstrap %||% 2000
    )
    
    # Store results
    all_results[[analysis$description]] <- result
    
    # Add to summary
    summary_row <- data.frame(
      Analysis = analysis$description,
      Group1_LD50 = round(result$ld50_group1, 2),
      Group2_LD50 = round(result$ld50_group2, 2),
      Difference = round(result$original_difference, 3),
      P_value = result$p_value,
      Significant = result$p_value < 0.05,
      Effect_Size = round(result$effect_size, 3),
      Wild_Bootstrap_Recommended = result$recommendation$use_wild_bootstrap
    )
    
    results_summary <- rbind(results_summary, summary_row)
    
    cat("\n")
  }
  
  cat("### BATCH ANALYSIS SUMMARY ###\n")
  print(results_summary)
  
  return(list(
    summary = results_summary,
    detailed_results = all_results
  ))
}

# Helper function for default values
`%||%` <- function(x, y) if (is.null(x)) y else x