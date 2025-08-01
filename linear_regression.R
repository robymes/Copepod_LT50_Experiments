if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("lmtest")) {
  install.packages("lmtest")
  library(lmtest)
}

source("plot_utilities.R")

# Initialize dataframes to store results for ANOVA statistical analysis
# These global dataframes accumulate data across multiple function calls
# for comparative analysis of different experimental conditions

# Dataframe to store LD50 values and their corresponding exposure times
# This is the primary data for time-to-death regression analysis
anova_data <- data.frame(
  ld50 = numeric(0),
  time_list = numeric(0),
  dir = character((0))
)

# Dataframe to store parameters from linear regressions
# Used to compare slopes and intercepts between different experimental conditions
anova_slopes <- data.frame(
  slope = numeric(0),
  intercept = numeric(0),
  dir = character(0)
)

# Dataframe to store regression statistics for t-tests
# Includes standard errors and degrees of freedom for comparing regression coefficients
t_test_data <- data.frame(
  slope = numeric(0),
  slope_std_err = numeric(0),
  intercept = numeric(0),
  intercept_std_err = numeric(0),
  data_points_length = numeric(0),
  dir = character(0)
)

# Helper function to calculate points along the Thermal Death Time (TDT) line
# Parameters:
# - slope: The z-value or temperature coefficient from regression
# - intercept: The temperature at reference time
# - z: The log-transformed time value
# Returns: The predicted LD50 temperature at given log-time
tdt_line_func <- function(slope, intercept, z) {
  return(intercept + slope * z)
}

# Main function for performing linear regression on thermal death time data
# Analyzes the relationship between LD50 temperatures and log-transformed exposure times
#
# Parameters:
# - dir_path: Path to the data directory (used for identification in output)
# - ld50: Vector of LD50 values from non-linear regression
# - time_list: Vector of exposure times corresponding to the LD50 values
# - main_title: Main title for generated charts
# - chart_subtitle: Subtitle for generated charts identifying the condition
#
# Returns a list containing:
# - anova_data: Data formatted for ANOVA analysis
# - anova_slopes: Regression parameters for slope comparisons
# - t_test_data: Statistical information for t-tests
# - survival_param_z: The z-value (temperature coefficient) for survival modeling
linear_regression_func <- function(dir_path, ld50, time_list, main_title, chart_subtitle) {
  # Create the Thermal death time curve (TDT) plot
  
  cat("########## LINEAR REGRESSION ##########\n\n")
  cat("########## DIRECTORY: ", dir_path, " ##########\n\n")
  
  # Create a data frame with LD50 and exposure times to construct TDT curve
  lt50 <- data.frame(ld50, time_list)
  lt50_temp <- lt50
  lt50_temp$dir <- chart_subtitle
  anova_data <- rbind(anova_data, lt50_temp)
  cat("LD50 values and corresponding LT50 exposure values\n")
  print(lt50)
  cat("\n")
  
  # Run linear regression with time in hours on log10 scale
  # This fits the classic TDT model: LD50 = intercept + slope * log10(time)
  # The slope represents the z-value (temperature coefficient)
  tdt_results_hours <- lm(lt50$ld50 ~ log10(lt50$time_list), data = lt50)
  tdt_slope <- tdt_results_hours$coefficients["log10(lt50$time_list)"]
  tdt_intercept_hours <- tdt_results_hours$coefficients["(Intercept)"]
  
  # Store regression parameters for comparative analysis
  new_anova_element <- data.frame(
    slope = tdt_slope,
    intercept = tdt_intercept_hours,
    dir = chart_subtitle
  )
  
  # Store detailed regression statistics for t-test comparisons
  # Includes standard errors which are needed for significance testing
  new_t_test_data <- data.frame(
    slope = summary(tdt_results_hours)$coefficients[2, 1],
    slope_std_err = summary(tdt_results_hours)$coefficients[2, 2],
    intercept = summary(tdt_results_hours)$coefficients[1, 1],
    intercept_std_err = summary(tdt_results_hours)$coefficients[1, 2],
    data_points_length = length(lt50_temp$time_list),
    dir = dir_path
  )
  
  # Prepare data for plotting
  df <- data.frame(x = log10(lt50$time_list), y = lt50$ld50)
  
  # Create plot with LD50/LT50 values on log-scale
  # Generate sequence of log10 time values for the fitted line
  z <- seq(0, 2, 0.01)
  df_lines <- data.frame(z = z, y = tdt_line_func(tdt_slope, tdt_intercept_hours, z))
  
  # Build the TDT plot using ggplot2
  linear_regression_plot <- ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    geom_line(data = df_lines, aes(x = z, y = y)) +
    xlim(0, 2) +
    ylim(0, 55) +
    labs(
      x = expression("log"[10] * "Time (LT50, hours)"),
      y = "Temperature (LD50, °C)",
      title = paste("Thermal death time (TDT) Curve\n", chart_subtitle)
    )
  
  # Display the plot
  print(linear_regression_plot)
  
  # Save the plot to disk with standardized naming and formatting
  save_plot_func(
    plot = linear_regression_plot,
    path = paste("charts/", gsub("\\\\", "/", gsub(" ", "", tools::toTitleCase(trimws(main_title)))), sep = ""),
    filename = paste(gsub(" ", "", trimws(chart_subtitle)), "_tdt", sep = ""),
    width = 1920,
    height = 1080
  )
  
  # Output detailed statistical information about the regression
  cat("Summary statistics for TDT curve\n")
  print(summary(tdt_results_hours))
  cat("R squared value for TDT curve: ", summary(tdt_results_hours)$r.squared, "\n\n")
  
  # Run a second linear regression with time in minutes 
  # This provides an alternative parameterization that might be more intuitive
  # for some applications (e.g., comparing to literature values)
  lt50_minutes <- lt50
  lt50_minutes$time_list <- lt50_minutes$time_list * 60
  tdt_results_minutes <- lm(lt50_minutes$ld50 ~ log10(lt50_minutes$time_list), data = lt50_minutes)
  tdt_intercept_minutes <- tdt_results_minutes$coefficients["(Intercept)"]
  cat("TDT intercept in minutes: ", tdt_intercept_minutes, "\n\n")
  
  # NEW: Verify normality
  normality_p <- verify_normality(ld50, chart_subtitle)
  
  # Save results for later comparisons
  normality_results <- list(
    p_value = normality_p,
    ld50_values = ld50,
    condition = chart_subtitle
  )
  
  # Return comprehensive results for further analysis
  return(list(
    anova_data = anova_data,
    anova_slopes = new_anova_element,
    t_test_data = new_t_test_data,
    survival_param_z = tdt_slope,
    normality_results = normality_results
  ))
}

# Function for comparative analysis of TDT curves using ANOVA
# Creates a combined plot of multiple TDT curves and performs statistical testing
# to determine if slopes differ significantly between experimental conditions
#
# Parameters:
# - anova_data: Combined dataframe of LD50 values from multiple conditions
# - anova_slopes: Dataframe containing regression parameters from multiple conditions
# - main_title: Title for the generated chart
anova_analysis_func <- function(anova_data, anova_slopes, main_title) {
  # Convert exposure time to log10 scale for analysis
  anova_data$time_list <- log10(anova_data$time_list)
  
  # Create a combined plot with all TDT curves for visual comparison
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
    scale_color_discrete(name = "Sample")
  
  # Add regression lines for each experimental condition
  for (i in 1:nrow(anova_slopes)) { # nolint: seq_linter.
    z <- seq(0, 2, 0.01)
    anova_slope_data <- tdt_line_func(anova_slopes[i, "slope"], anova_slopes[i, "intercept"], z)
    anova_slope_df <- data.frame(ld50 = anova_slope_data, time_list = z, dir = anova_slopes[i, "dir"])
    anova_plot <- anova_plot +
      geom_line(data = anova_slope_df, aes(group = dir))
  }
  
  # Display the combined plot
  print(anova_plot)
  
  # Save the plot to disk with standardized formatting
  save_plot_func(
    plot = anova_plot,
    path = paste("charts/", gsub("\\\\", "/", gsub(" ", "", tools::toTitleCase(trimws(main_title)))), sep = ""),
    filename = "tdt",
    width = 1920,
    height = 1080
  )
  
  # Perform ANOVA to test for significant differences in the temperature-time relationship
  # between different experimental conditions
  # Tests both main effects and interaction effects (different slopes)
  anova_analysis <- anova(lm(ld50 ~ log(time_list) * dir, anova_data))
  cat("########## ANOVA ANALYSIS ##########\n\n")
  print(anova_analysis)
  cat("\n")
}

# Function to perform t-tests comparing regression parameters between conditions
# Used when comparing exactly two experimental conditions to determine if 
# their TDT curves differ significantly
#
# Parameters:
# - t_test_data: Dataframe containing regression statistics from different conditions
# - main_title: Title for the analysis (used for logging)
t_test_func <- function(t_test_data, normality_results_list, main_title) {
  # Display header for t-test results
  cat("########## T TEST ##########\n\n")
  
  # Only perform test if exactly two conditions are being compared
  if (length(t_test_data$slope) == 2) {
    
    # NEW: Verify assumptions before t-test
    cat("############################################\n")
    cat("### Statistical Assumptions Verification ###\n")
    cat("############################################\n")
    
    # Extract LD50 values from both groups
    ld50_group1 <- normality_results_list[[1]]$ld50_values
    ld50_group2 <- normality_results_list[[2]]$ld50_values
    
    # Verify variance homogeneity
    homogeneity_results <- verify_homogeneity(
      ld50_group1, ld50_group2,
      normality_results_list[[1]]$condition,
      normality_results_list[[2]]$condition
    )
    
    # Decision on test type
    use_welch <- homogeneity_results$levene_test < 0.05
    
    if (use_welch) {
      cat("\n>>> Using Welch's t-test (unequal variances) <<<\n\n")
    } else {
      cat("\n>>> Using Student's t-test (equal variances) <<<\n\n")
    }
    
    # If data are not normal, consider alternatives
    if (any(c(normality_results_list[[1]]$p_value, 
              normality_results_list[[2]]$p_value) < 0.05)) {
      cat("\n!!! WARNING: Data not normally distributed !!!\n")
      cat("Consider using Mann-Whitney U test as alternative:\n\n")
      
      # Mann-Whitney U test
      mw_test <- wilcox.test(ld50_group1, ld50_group2)
      cat("Mann-Whitney U test p-value:", mw_test$p.value, "\n\n")
    }
    
    # Calculate t-statistic for slope comparison using the formula:
    # t = (β₁ - β₂) / sqrt(SE₁² + SE₂²)
    t_value_slope <-
      (t_test_data$slope[1] - t_test_data$slope[2]) /
      sqrt(t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)
    
    # Calculate degrees of freedom using Welch-Satterthwaite approximation
    # This accounts for potentially unequal variances and sample sizes
    df_value_slope <-
      (t_test_data$slope_std_err[1]^2 + t_test_data$slope_std_err[2]^2)^2 /
      (t_test_data$slope_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) +
         t_test_data$slope_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
    
    # Calculate two-tailed p-value
    p_value_slope <- 2 * (1 - pt(abs(t_value_slope), df_value_slope))
    
    # Repeat the process for intercepts
    t_value_intercept <-
      (t_test_data$intercept[1] - t_test_data$intercept[2]) /
      sqrt(t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)
    df_value_intercept <-
      (t_test_data$intercept_std_err[1]^2 + t_test_data$intercept_std_err[2]^2)^2 /
      (t_test_data$intercept_std_err[1]^4 / (t_test_data$data_points_length[1] - 2) +
         t_test_data$intercept_std_err[2]^4 / (t_test_data$data_points_length[2] - 2))
    p_value_intercept <- 2 * (1 - pt(abs(t_value_intercept), df_value_intercept))
    
    # Output p-values for both tests
    t_test_threshold <- 0.05
    cat("P Value slopes: ", p_value_slope, "\n")
    if (p_value_slope < t_test_threshold) {
      cat("P Value slopes: is statistically significant\n")
    } else {
      cat("P Value slopes: is NOT statistically significant\n")
    }
    cat("P Value intercept: ", p_value_intercept, "\n")
    if (p_value_intercept < t_test_threshold) {
      cat("P Value intercept: is statistically significant\n")
    } else {
      cat("P Value intercept: is NOT statistically significant\n")
    }
    cat("\n")
  } else {
    # Warn if the t-test cannot be performed due to incorrect number of conditions
    cat("T test cannot be executed with ", length(t_test_data$slope), " samples\n")
  }
}

# Function to verify normality of LD50 data
verify_normality <- function(ld50_values, condition_name) {
  cat("###############################################\n")
  cat("\n### Normality Test for", condition_name, "###\n")
  cat("###############################################\n")
  
  # Shapiro-Wilk test (for small samples, n < 50)
  shapiro_test <- shapiro.test(ld50_values)
  cat("Shapiro-Wilk test:\n")
  cat("  W =", shapiro_test$statistic, "\n")
  cat("  p-value =", shapiro_test$p.value, "\n")
  cat("  Interpretation:", ifelse(shapiro_test$p.value > 0.05, 
                                  "Data compatible with normal distribution", 
                                  "Data NOT normally distributed"), "\n\n")
  
  # Q-Q plot for visual assessment
  qqnorm(ld50_values, main = paste("Q-Q Plot -", condition_name))
  qqline(ld50_values, col = "red")
  
  # Histogram with overlaid normal curve
  hist(ld50_values, probability = TRUE, 
       main = paste("LD50 Histogram -", condition_name),
       xlab = "LD50 (°C)")
  curve(dnorm(x, mean = mean(ld50_values), sd = sd(ld50_values)), 
        add = TRUE, col = "blue", lwd = 2)
  
  return(shapiro_test$p.value)
}

# Function to verify homogeneity of variance using Breusch-Pagan test
# Updated to maintain consistency with wild bootstrap analysis methodology
verify_homogeneity <- function(ld50_group1, ld50_group2, 
                               name_group1, name_group2) {
  cat("\n### Homogeneity of Variance Test ###\n")
  
  # F-test to compare variances (kept for compatibility)
  f_test <- var.test(ld50_group1, ld50_group2)
  cat("F-test:\n")
  cat("  F =", f_test$statistic, "\n")
  cat("  p-value =", f_test$p.value, "\n")
  cat("  Interpretation:", ifelse(f_test$p.value > 0.05, 
                                  "Homogeneous variances", 
                                  "Non-homogeneous variances"), "\n\n")
  
  # Breusch-Pagan test (replaces Levene's test for consistency)
  # Create combined dataset for regression-based heteroscedasticity test
  combined_data <- data.frame(
    ld50 = c(ld50_group1, ld50_group2),
    group = factor(c(rep(name_group1, length(ld50_group1)),
                     rep(name_group2, length(ld50_group2))))
  )
  
  # Fit linear model: LD50 ~ group
  linear_model <- lm(ld50 ~ group, data = combined_data)
  
  # Apply Breusch-Pagan test to residuals
  if (!require("lmtest", quietly = TRUE)) {
    install.packages("lmtest")
    library(lmtest)
  }
  
  bp_test <- lmtest::bptest(linear_model)
  
  cat("Breusch-Pagan test:\n")
  cat("  LM statistic:", round(bp_test$statistic, 4), "\n")
  cat("  Degrees of freedom:", bp_test$parameter, "\n")
  cat("  P-value:", format(bp_test$p.value, scientific = TRUE), "\n")
  cat("  Interpretation:", ifelse(bp_test$p.value > 0.05,
                                  "Homoscedastic (constant variance)",
                                  "Heteroscedastic (non-constant variance)"), "\n\n")
  
  return(list(f_test = f_test$p.value, 
              levene_test = bp_test$p.value))  # Keep same return structure
}