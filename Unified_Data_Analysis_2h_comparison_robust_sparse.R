# METODI ROBUSTI PER CONFRONTO PARAMETRI NLS CON RESIDUI NON NORMALI
# Quando t-test e bootstrap case resampling falliscono

if (!require("permute")) {
  install.packages("permute")
  library(permute)
}

if (!require("coin")) {
  install.packages("coin")
  library(coin)
}

if (!require("boot")) {
  install.packages("boot")
  library(boot)
}

# =====================================================================
# 1. PERMUTATION TEST (METODO PIÙ ROBUSTO)
# =====================================================================

permutation_test_nls <- function(data_path, file_group1, file_group2, 
                                 group1_name = "Group 1", group2_name = "Group 2",
                                 n_permutations = 5000,
                                 nls_param_list = list(list(100, 30, 4))) {
  
  cat("=== PERMUTATION TEST FOR NLS PARAMETERS ===\n")
  cat("Folder:", data_path, "\n")
  cat("Comparing:", group1_name, "vs", group2_name, "\n")
  cat("Permutations:", n_permutations, "\n\n")
  
  # Read data
  path_group1 <- file.path(data_path, file_group1)
  path_group2 <- file.path(data_path, file_group2)
  
  data1 <- read.csv(path_group1, row.names = NULL)
  data2 <- read.csv(path_group2, row.names = NULL)
  
  data1$Survival <- data1$Alive / (data1$Alive + data1$Dead)
  data2$Survival <- data2$Alive / (data2$Alive + data2$Dead)
  
  # Original models and difference
  model1 <- fit_nls_safe(data1, nls_param_list[[1]], group1_name)
  model2 <- fit_nls_safe(data2, nls_param_list[[1]], group2_name)
  
  if (!model1$converged || !model2$converged) {
    stop("Original models failed to converge")
  }
  
  original_diff <- coef(model1$model)[2] - coef(model2$model)[2]
  
  cat("Original LD50 difference:", round(original_diff, 4), "°C\n")
  cat("(", group1_name, ":", round(coef(model1$model)[2], 2), "°C -", 
      group2_name, ":", round(coef(model2$model)[2], 2), "°C )\n\n")
  
  # Combined dataset for permutation
  combined_data <- rbind(
    transform(data1, group = group1_name),
    transform(data2, group = group2_name)
  )
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  # Permutation test
  cat("Running", n_permutations, "permutations...\n")
  
  perm_differences <- numeric(n_permutations)
  successful_perms <- 0
  
  for (i in 1:n_permutations) {
    # Random permutation of group labels
    perm_indices <- sample(nrow(combined_data))
    perm_data <- combined_data[perm_indices, ]
    
    # Split back into two groups
    perm_data1 <- perm_data[1:n1, ]
    perm_data2 <- perm_data[(n1+1):(n1+n2), ]
    
    # Fit models
    perm_model1 <- fit_nls_safe(perm_data1, nls_param_list[[1]], paste("Perm", i, "G1"))
    perm_model2 <- fit_nls_safe(perm_data2, nls_param_list[[1]], paste("Perm", i, "G2"))
    
    if (perm_model1$converged && perm_model2$converged) {
      perm_differences[i] <- coef(perm_model1$model)[2] - coef(perm_model2$model)[2]
      successful_perms <- successful_perms + 1
    } else {
      perm_differences[i] <- NA
    }
    
    if (i %% 1000 == 0) cat("Completed", i, "permutations...\n")
  }
  
  # Remove failed permutations
  perm_differences <- perm_differences[!is.na(perm_differences)]
  
  cat("\nSuccessful permutations:", length(perm_differences), "/", n_permutations, 
      "(", round(100 * length(perm_differences) / n_permutations, 1), "%)\n")
  
  # Calculate p-value
  p_value <- mean(abs(perm_differences) >= abs(original_diff))
  
  # Results
  cat("\n=== PERMUTATION TEST RESULTS ===\n")
  cat("Original difference:", round(original_diff, 4), "°C\n")
  cat("P-value (two-tailed):", format(p_value, scientific = TRUE), "\n")
  cat("Interpretation:", ifelse(p_value < 0.05,
                               paste("Statistically significant difference:", 
                                     ifelse(original_diff > 0, group1_name, group2_name), 
                                     "shows higher thermal tolerance"),
                               "No statistically significant difference"), "\n")
  
  # Distribution statistics
  cat("\nPermutation distribution:\n")
  cat("Mean:", round(mean(perm_differences), 4), "°C\n")
  cat("SD:", round(sd(perm_differences), 4), "°C\n")
  cat("95% CI:", round(quantile(perm_differences, c(0.025, 0.975)), 4), "°C\n")
  
  return(list(
    original_difference = original_diff,
    p_value = p_value,
    permutation_differences = perm_differences,
    successful_permutations = length(perm_differences),
    method = "Permutation Test"
  ))
}

# =====================================================================
# 2. ROBUST BOOTSTRAP (WILD BOOTSTRAP)
# =====================================================================

wild_bootstrap_nls <- function(data_path, file_group1, file_group2,
                               group1_name = "Group 1", group2_name = "Group 2",
                               n_bootstrap = 2000,
                               nls_param_list = list(list(100, 30, 4))) {
  
  cat("=== WILD BOOTSTRAP FOR NLS PARAMETERS ===\n")
  
  # Read and fit original models
  path_group1 <- file.path(data_path, file_group1)
  path_group2 <- file.path(data_path, file_group2)
  
  data1 <- read.csv(path_group1, row.names = NULL)
  data2 <- read.csv(path_group2, row.names = NULL)
  
  data1$Survival <- data1$Alive / (data1$Alive + data1$Dead)
  data2$Survival <- data2$Alive / (data2$Alive + data2$Dead)
  
  model1 <- fit_nls_safe(data1, nls_param_list[[1]], group1_name)
  model2 <- fit_nls_safe(data2, nls_param_list[[1]], group2_name)
  
  if (!model1$converged || !model2$converged) {
    stop("Original models failed to converge")
  }
  
  # Get fitted values and residuals
  fitted1 <- fitted(model1$model)
  fitted2 <- fitted(model2$model)
  residuals1 <- residuals(model1$model)
  residuals2 <- residuals(model2$model)
  
  original_diff <- coef(model1$model)[2] - coef(model2$model)[2]
  
  # Wild bootstrap
  boot_differences <- numeric(n_bootstrap)
  successful_boots <- 0
  
  cat("Running", n_bootstrap, "wild bootstrap iterations...\n")
  
  for (i in 1:n_bootstrap) {
    # Generate wild bootstrap weights (Rademacher distribution)
    weights1 <- sample(c(-1, 1), length(residuals1), replace = TRUE)
    weights2 <- sample(c(-1, 1), length(residuals2), replace = TRUE)
    
    # Create bootstrap samples
    y1_boot <- fitted1 + weights1 * residuals1
    y2_boot <- fitted2 + weights2 * residuals2
    
    # Constrain to [0,1] for survival proportions
    y1_boot <- pmax(0, pmin(1, y1_boot))
    y2_boot <- pmax(0, pmin(1, y2_boot))
    
    data1_boot <- data.frame(Temperature = data1$Temperature, Survival = y1_boot)
    data2_boot <- data.frame(Temperature = data2$Temperature, Survival = y2_boot)
    
    # Fit bootstrap models
    boot_model1 <- fit_nls_safe(data1_boot, nls_param_list[[1]], paste("Boot", i, "G1"))
    boot_model2 <- fit_nls_safe(data2_boot, nls_param_list[[1]], paste("Boot", i, "G2"))
    
    if (boot_model1$converged && boot_model2$converged) {
      boot_differences[i] <- coef(boot_model1$model)[2] - coef(boot_model2$model)[2]
      successful_boots <- successful_boots + 1
    } else {
      boot_differences[i] <- NA
    }
    
    if (i %% 500 == 0) cat("Completed", i, "bootstrap iterations...\n")
  }
  
  # Remove failed iterations
  boot_differences <- boot_differences[!is.na(boot_differences)]
  
  cat("\nSuccessful bootstrap iterations:", length(boot_differences), "/", n_bootstrap,
      "(", round(100 * length(boot_differences) / n_bootstrap, 1), "%)\n")
  
  # Calculate statistics
  p_value <- 2 * min(mean(boot_differences >= 0), mean(boot_differences <= 0))
  ci <- quantile(boot_differences, c(0.025, 0.975))
  
  cat("\n=== WILD BOOTSTRAP RESULTS ===\n")
  cat("Original difference:", round(original_diff, 4), "°C\n")
  cat("Bootstrap mean:", round(mean(boot_differences), 4), "°C\n")
  cat("95% CI:", round(ci, 4), "°C\n")
  cat("P-value:", format(p_value, scientific = TRUE), "\n")
  cat("Interpretation:", ifelse(p_value < 0.05,
                               "Statistically significant difference",
                               "No statistically significant difference"), "\n")
  
  return(list(
    original_difference = original_diff,
    bootstrap_differences = boot_differences,
    p_value = p_value,
    ci = ci,
    method = "Wild Bootstrap"
  ))
}

# =====================================================================
# 3. SIMULATION-BASED INFERENCE
# =====================================================================

simulation_based_test <- function(data_path, file_group1, file_group2,
                                  group1_name = "Group 1", group2_name = "Group 2",
                                  n_simulations = 3000,
                                  nls_param_list = list(list(100, 30, 4))) {
  
  cat("=== SIMULATION-BASED INFERENCE ===\n")
  
  # Fit original models
  path_group1 <- file.path(data_path, file_group1)
  path_group2 <- file.path(data_path, file_group2)
  
  data1 <- read.csv(path_group1, row.names = NULL)
  data2 <- read.csv(path_group2, row.names = NULL)
  
  data1$Survival <- data1$Alive / (data1$Alive + data1$Dead)
  data2$Survival <- data2$Alive / (data2$Alive + data2$Dead)
  
  model1 <- fit_nls_safe(data1, nls_param_list[[1]], group1_name)
  model2 <- fit_nls_safe(data2, nls_param_list[[1]], group2_name)
  
  if (!model1$converged || !model2$converged) {
    stop("Original models failed to converge")
  }
  
  original_diff <- coef(model1$model)[2] - coef(model2$model)[2]
  
  # Estimate error distribution from residuals (using empirical distribution)
  residuals1 <- residuals(model1$model)
  residuals2 <- residuals(model2$model)
  
  # Simulate under null hypothesis (no difference in LD50)
  # Use average LD50 for both groups
  avg_ld50 <- mean(c(coef(model1$model)[2], coef(model2$model)[2]))
  
  sim_differences <- numeric(n_simulations)
  successful_sims <- 0
  
  cat("Running", n_simulations, "simulations under null hypothesis...\n")
  
  for (i in 1:n_simulations) {
    # Generate data under null (same LD50 for both groups)
    # Use original p1 and p3, but common p2 (LD50)
    
    p1_1 <- coef(model1$model)[1]
    p1_2 <- coef(model2$model)[1]
    p3_1 <- coef(model1$model)[3]
    p3_2 <- coef(model2$model)[3]
    
    # Simulate survival under null hypothesis
    sim_fitted1 <- p1_1 / (1 + (data1$Temperature / avg_ld50)^p3_1)
    sim_fitted2 <- p1_2 / (1 + (data2$Temperature / avg_ld50)^p3_2)
    
    # Add noise from empirical residual distribution
    sim_y1 <- sim_fitted1 + sample(residuals1, length(sim_fitted1), replace = TRUE)
    sim_y2 <- sim_fitted2 + sample(residuals2, length(sim_fitted2), replace = TRUE)
    
    # Constrain to [0,1]
    sim_y1 <- pmax(0, pmin(1, sim_y1))
    sim_y2 <- pmax(0, pmin(1, sim_y2))
    
    sim_data1 <- data.frame(Temperature = data1$Temperature, Survival = sim_y1)
    sim_data2 <- data.frame(Temperature = data2$Temperature, Survival = sim_y2)
    
    # Fit models to simulated data
    sim_model1 <- fit_nls_safe(sim_data1, nls_param_list[[1]], paste("Sim", i, "G1"))
    sim_model2 <- fit_nls_safe(sim_data2, nls_param_list[[1]], paste("Sim", i, "G2"))
    
    if (sim_model1$converged && sim_model2$converged) {
      sim_differences[i] <- coef(sim_model1$model)[2] - coef(sim_model2$model)[2]
      successful_sims <- successful_sims + 1
    } else {
      sim_differences[i] <- NA
    }
    
    if (i %% 500 == 0) cat("Completed", i, "simulations...\n")
  }
  
  # Remove failed simulations
  sim_differences <- sim_differences[!is.na(sim_differences)]
  
  cat("\nSuccessful simulations:", length(sim_differences), "/", n_simulations,
      "(", round(100 * length(sim_differences) / n_simulations, 1), "%)\n")
  
  # Calculate p-value
  p_value <- mean(abs(sim_differences) >= abs(original_diff))
  
  cat("\n=== SIMULATION-BASED TEST RESULTS ===\n")
  cat("Original difference:", round(original_diff, 4), "°C\n")
  cat("Null distribution mean:", round(mean(sim_differences), 4), "°C\n")
  cat("Null distribution SD:", round(sd(sim_differences), 4), "°C\n")
  cat("P-value:", format(p_value, scientific = TRUE), "\n")
  cat("Interpretation:", ifelse(p_value < 0.05,
                               "Statistically significant difference",
                               "No statistically significant difference"), "\n")
  
  return(list(
    original_difference = original_diff,
    null_differences = sim_differences,
    p_value = p_value,
    method = "Simulation-based Test"
  ))
}

# =====================================================================
# HELPER FUNCTION: SAFE NLS FITTING
# =====================================================================

fit_nls_safe <- function(data, start_params, name = "Model", silent = TRUE) {
  
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

# =====================================================================
# MAIN FUNCTION: COMPREHENSIVE ROBUST ANALYSIS
# =====================================================================

robust_nls_analysis <- function(data_path, file_group1, file_group2,
                                group1_name = "Group 1", group2_name = "Group 2",
                                nls_param_list = list(list(100, 30, 4))) {
  
  cat("##################################################\n")
  cat("### COMPREHENSIVE ROBUST NLS PARAMETER ANALYSIS ###\n")
  cat("### FOR NON-NORMAL RESIDUALS ###\n")
  cat("##################################################\n\n")
  
  # Run all three methods
  cat("1. PERMUTATION TEST (Most Robust)\n")
  cat("====================================\n")
  perm_results <- permutation_test_nls(data_path, file_group1, file_group2,
                                       group1_name, group2_name, 
                                       n_permutations = 3000, nls_param_list)
  
  cat("\n\n2. WILD BOOTSTRAP (Robust to Heteroscedasticity)\n")
  cat("=================================================\n")
  wild_results <- wild_bootstrap_nls(data_path, file_group1, file_group2,
                                     group1_name, group2_name,
                                     n_bootstrap = 1500, nls_param_list)
  
  cat("\n\n3. SIMULATION-BASED TEST (Model-Based)\n")
  cat("=======================================\n")
  sim_results <- simulation_based_test(data_path, file_group1, file_group2,
                                       group1_name, group2_name,
                                       n_simulations = 2000, nls_param_list)
  
  # Summary comparison
  cat("\n\n### SUMMARY OF ALL METHODS ###\n")
  cat("===============================\n")
  
  methods_summary <- data.frame(
    Method = c("Permutation Test", "Wild Bootstrap", "Simulation-based"),
    P_value = c(perm_results$p_value, wild_results$p_value, sim_results$p_value),
    Significant = c(perm_results$p_value < 0.05, 
                   wild_results$p_value < 0.05, 
                   sim_results$p_value < 0.05),
    Recommended = c("Primary", "Secondary", "Supplementary")
  )
  
  print(methods_summary)
  
  cat("\n### FINAL RECOMMENDATION ###\n")
  if (sum(methods_summary$Significant) >= 2) {
    cat("✅ CONCLUSION: Statistically significant difference\n")
    cat("Multiple robust methods agree on significance.\n")
  } else if (sum(methods_summary$Significant) == 1) {
    cat("⚠️  CONCLUSION: Borderline significance\n") 
    cat("Methods disagree - consider larger sample size.\n")
  } else {
    cat("❌ CONCLUSION: No statistically significant difference\n")
    cat("All robust methods agree on non-significance.\n")
  }
  
  return(list(
    permutation = perm_results,
    wild_bootstrap = wild_results,
    simulation = sim_results,
    summary = methods_summary
  ))
}

data_path <- "data/cruise/t_test/tubeworms_mussels/oxic"
results_pressure <- robust_nls_analysis(
  data_path = data_path,
  file_group1 = "No_Pressure_Data02.csv",
  file_group2 = "Pressure_Data02.csv",
  group1_name = "No Pressure", 
  group2_name = "Pressure"
)

# =====================================================================
# USAGE EXAMPLES
# =====================================================================

# Example: Tubeworms with non-normal residuals
# data_path <- "data/cruise/t_test/tubeworms_mussels/oxic"
# results <- robust_nls_analysis(
#   data_path = data_path,
#   file_group1 = "No_Pressure_Data02.csv",
#   file_group2 = "Pressure_Data02.csv",
#   group1_name = "No Pressure",
#   group2_name = "Pressure"
# )

# Example: Pompeii worms oxygen comparison
# data_path <- "data/cruise/t_test/pompeii_worms/anoxic"
# results <- robust_nls_analysis(
#   data_path = data_path,
#   file_group1 = "No_Pressure_Data02.csv", 
#   file_group2 = "Pressure_Data02.csv",
#   group1_name = "No Pressure",
#   group2_name = "Pressure"
# )