# Install necessary packages (run only once if not already installed)
# install.packages("lavaan")

library(lavaan)

# --- Simulation Parameters ---
num_replications <- 100 # Number of Monte Carlo replications per scenario
sample_sizes <- c(100, 300) # Sample sizes to investigate
num_items <- 10 # Number of items in the CFA model

# True parameters for the CFA model (1 factor, 10 items)
true_factor_loadings_cfa <- rep(0.7, num_items) # True factor loadings for all items
true_factor_variance <- 1 # Latent factor variance (fixed to 1 for model identification)

# --- Data Distribution Scenarios ---
# Theta is always Normal (latent factor) to isolate error impact.
# Error distributions are designed to be more extreme and are *not* scaled to a fixed variance.
distribution_scenarios <- list(
  "Normal_Errors" = list(
    theta_dist = list(type = "norm", mean = 0, sd = 1),
    error_dist = list(type = "norm_scaled", sd = sqrt(1 - 0.7^2)) # Explicitly scaled to true error variance
  ),
  "Skewed_Errors_LogNormal" = list( # More severely skewed errors
    theta_dist = list(type = "norm", mean = 0, sd = 1),
    error_dist = list(type = "lognorm_raw", meanlog = 0, sdlog = 1.5) # Increased sdlog for more skew
  ),
  "Heavy_Tailed_Errors_t3df" = list( # Very heavy-tailed errors (df=3)
    theta_dist = list(type = "norm", mean = 0, sd = 1),
    error_dist = list(type = "t_raw", df = 3) # Extremely heavy tails, variance defined but large
  ),
  "Heavy_Tailed_Errors_t2df" = list( # Extremely heavy-tailed errors (df=2), infinite variance
    theta_dist = list(type = "norm", mean = 0, sd = 1),
    error_dist = list(type = "t_raw", df = 2) # Variance undefined/infinite, should challenge ML
  )
)

# --- List to store final summarized CFA results ---
all_cfa_results_summary <- list()

# --- Main Simulation Loop ---
scenario_counter <- 1 # Counter for storing results in the lists

for (dist_name in names(distribution_scenarios)) {
  current_scenario <- distribution_scenarios[[dist_name]]

  for (N in sample_sizes) {
    cat(paste0("\nRunning Scenario: Distribution = ", dist_name, ", N = ", N, "\n"))

    # --- Matrices to store results for each replication in this scenario (pre-filled with NA) ---
    # CFA (for both ML and MLR)
    cfa_ml_bias_loadings_repl <- matrix(NA, nrow = num_replications, ncol = num_items)
    cfa_ml_rmse_loadings_repl_sq <- matrix(NA, nrow = num_replications, ncol = num_items)
    cfa_ml_fit_indices_repl <- matrix(NA, nrow = num_replications, ncol = 3, dimnames = list(NULL, c("CFI", "RMSEA", "TLI")))
    cfa_ml_converged_repl <- rep(FALSE, num_replications)

    cfa_mlr_bias_loadings_repl <- matrix(NA, nrow = num_replications, ncol = num_items)
    cfa_mlr_rmse_loadings_repl_sq <- matrix(NA, nrow = num_replications, ncol = num_items)
    cfa_mlr_fit_indices_repl <- matrix(NA, nrow = num_replications, ncol = 3, dimnames = list(NULL, c("CFI", "RMSEA", "TLI")))
    cfa_mlr_converged_repl <- rep(FALSE, num_replications)

    for (i in 1:num_replications) {
      # Display simulation progress
      if (i %% 10 == 0) cat(paste0("  Replication: ", i, "/", num_replications, "\n"))

      # --- 1. Data Generation ---
      # Generate latent factor (theta) for individuals (always normal in this setup)
      theta <- rnorm(N, mean = current_scenario$theta_dist$mean, sd = current_scenario$theta_dist$sd)

      # Generate observed data (item scores)
      data_simulated <- matrix(NA, nrow = N, ncol = num_items)
      for (j in 1:num_items) {
        errors <- NULL # Initialize errors

        # Generate errors (residuals) based on the current distribution scenario
        if (current_scenario$error_dist$type == "norm_scaled") {
          errors <- rnorm(N, 0, current_scenario$error_dist$sd)
        } else if (current_scenario$error_dist$type == "lognorm_raw") {
          error_raw <- rlnorm(N, meanlog = current_scenario$error_dist$meanlog, sdlog = current_scenario$error_dist$sdlog)
          errors <- error_raw - mean(error_raw) # Center to zero mean, but let variance be its natural log-normal variance
          # Add a small amount of noise if variance is near zero to prevent issues
          if (sd(errors) < 1e-6) errors <- rnorm(N, 0, 0.1)
        } else if (current_scenario$error_dist$type == "t_raw") {
          error_raw <- rt(N, df = current_scenario$error_dist$df)
          errors <- error_raw - mean(error_raw) # Center to zero mean, let variance be its natural t-distribution variance
          # Add a small amount of noise if variance is near zero to prevent issues
          if (sd(errors) < 1e-6) errors <- rnorm(N, 0, 0.1)
        }
        
        # CFA: Item = Factor_Loading * Latent_Factor + Error
        data_simulated[,j] <- true_factor_loadings_cfa[j] * theta + errors
      }

      # Convert simulated data matrix to a data.frame for lavaan
      data_simulated_df <- as.data.frame(data_simulated)
      colnames(data_simulated_df) <- paste0("item", 1:num_items)

      # --- 2. Run CFA (ML estimator) ---
      cfa_model <- '
        F1 =~ item1 + item2 + item3 + item4 + item5 + item6 + item7 + item8 + item9 + item10
      '
      fit_cfa_ml <- tryCatch({
        cfa(cfa_model, data = data_simulated_df, std.lv = TRUE, estimator = "ML", se = "standard")
        }, error = function(e) { NULL }, warning = function(w) { NULL })

      if (!is.null(fit_cfa_ml) && lavInspect(fit_cfa_ml, "converged")) {
        params_cfa_ml <- parameterEstimates(fit_cfa_ml)
        loadings_cfa_ml <- params_cfa_ml[params_cfa_ml$op == "=~", "est"]
        
        if (length(loadings_cfa_ml) == num_items) {
            cfa_ml_bias_loadings_repl[i,] <- loadings_cfa_ml - true_factor_loadings_cfa
            cfa_ml_rmse_loadings_repl_sq[i,] <- (loadings_cfa_ml - true_factor_loadings_cfa)^2
            fit_measures_cfa_ml <- fitMeasures(fit_cfa_ml, c("cfi", "rmsea", "tli"))
            cfa_ml_fit_indices_repl[i,] <- fit_measures_cfa_ml
            cfa_ml_converged_repl[i] <- TRUE
        }
      }

      # --- 2b. Run CFA (MLR estimator - Robust) ---
      fit_cfa_mlr <- tryCatch({
        cfa(cfa_model, data = data_simulated_df, std.lv = TRUE, estimator = "MLR") # MLR implies robust SE and scaled test stats
        }, error = function(e) { NULL }, warning = function(w) { NULL })

      if (!is.null(fit_cfa_mlr) && lavInspect(fit_cfa_mlr, "converged")) {
        params_cfa_mlr <- parameterEstimates(fit_cfa_mlr)
        loadings_cfa_mlr <- params_cfa_mlr[params_cfa_mlr$op == "=~", "est"]
        
        if (length(loadings_cfa_mlr) == num_items) {
            cfa_mlr_bias_loadings_repl[i,] <- loadings_cfa_mlr - true_factor_loadings_cfa
            cfa_mlr_rmse_loadings_repl_sq[i,] <- (loadings_cfa_mlr - true_factor_loadings_cfa)^2
            fit_measures_cfa_mlr <- fitMeasures(fit_cfa_mlr, c("cfi", "rmsea", "tli")) # MLR CFI/RMSEA/TLI are also robust
            cfa_mlr_fit_indices_repl[i,] <- fit_measures_cfa_mlr
            cfa_mlr_converged_repl[i] <- TRUE
        }
      }
    } # End of replication loop

    # --- 4. Calculate average metrics for this scenario and store in final list ---
    
    # CFA ML Results
    avg_cfa_ml_bias_loadings <- colMeans(cfa_ml_bias_loadings_repl, na.rm = TRUE)
    avg_cfa_ml_rmse_loadings <- sqrt(colMeans(cfa_ml_rmse_loadings_repl_sq, na.rm = TRUE))
    avg_cfa_ml_fit_indices <- colMeans(cfa_ml_fit_indices_repl, na.rm = TRUE)
    num_cfa_ml_converged <- sum(cfa_ml_converged_repl)

    # CFA MLR Results
    avg_cfa_mlr_bias_loadings <- colMeans(cfa_mlr_bias_loadings_repl, na.rm = TRUE)
    avg_cfa_mlr_rmse_loadings <- sqrt(colMeans(cfa_mlr_rmse_loadings_repl_sq, na.rm = TRUE))
    avg_cfa_mlr_fit_indices <- colMeans(cfa_mlr_fit_indices_repl, na.rm = TRUE)
    num_cfa_mlr_converged <- sum(cfa_mlr_converged_repl)

    # Store results for this scenario
    all_cfa_results_summary[[scenario_counter]] <- data.frame(
      Distribution = dist_name,
      N = N,
      Estimator = c("ML", "MLR"), # Indicate which estimator the results belong to
      Avg_Bias_Loadings_Item1 = c(avg_cfa_ml_bias_loadings[1], avg_cfa_mlr_bias_loadings[1]),
      Avg_RMSE_Loadings_Item1 = c(avg_cfa_ml_rmse_loadings[1], avg_cfa_mlr_rmse_loadings[1]),
      CFI = c(avg_cfa_ml_fit_indices["CFI"], avg_cfa_mlr_fit_indices["CFI"]),
      RMSEA = c(avg_cfa_ml_fit_indices["RMSEA"], avg_cfa_mlr_fit_indices["RMSEA"]),
      TLI = c(avg_cfa_ml_fit_indices["TLI"], avg_cfa_mlr_fit_indices["TLI"]),
      Converged_Count = c(num_cfa_ml_converged, num_cfa_mlr_converged),
      stringsAsFactors = FALSE
    )

    scenario_counter <- scenario_counter + 1
  } # End of N (sample size) loop
} # End of distribution scenario loop

# --- Display and Analyze Final Simulation Results ---
# Convert the list of summarized results into final data.frame
df_cfa_results_final <- do.call(rbind, all_cfa_results_summary)

cat("\n--- Final Summary of Confirmatory Factor Analysis (CFA) Results ---\n")
print(df_cfa_results_final)
