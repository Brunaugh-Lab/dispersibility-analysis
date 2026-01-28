# ==============================================================================
# 02_wasserstein_core.R
# Core Wasserstein-1 Distance Calculations for Dispersibility Analysis
#
# Purpose: Calculate W1 distance between test and reference particle size
#          distributions following proper methodology (pool replicates first)
#          Automatically reads processed data and saves results.
#
# Auto-execution: Script automatically runs when sourced
#   - Reads data_v2/tidy_v2/standardized_data_with_conditions.csv (from script 01)
#   - Calculates W1 distances for all formulations
#   - Saves to results_v2/wasserstein_results.csv
#
# Input: data_v2/tidy_v2/standardized_data_with_conditions.csv (from 01_data_import.R)
# Output: results_v2/wasserstein_results.csv
#
# Methodology: Wasserstein-1 (Earth Mover's) distance quantifies the minimum
#              redistribution work needed to transform one distribution into another
#
# Usage:
#   source("scripts/02_wasserstein_core.R")  # That's it!
#
# Reference: Brunaugh et al. (2025) Pharmaceutics - "A Distribution-Based Metric
#            for Quantifying Dispersibility in Dry Powder Inhalers"
# ==============================================================================

library(tidyverse)

# ==============================================================================
# HELPER FUNCTION: Pool replicate CDFs before W1 calculation
# ==============================================================================

#' Pool Technical Replicates into Single Representative CDF
#'
#' CRITICAL: Wasserstein distance should be calculated between POOLED replicate
#' distributions, not between individual replicates. This function averages
#' replicate CDFs at each particle size to create a single representative
#' distribution per condition.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param formulation Character string specifying formulation to pool
#' @param module Character string specifying module ("INHALER" or "RODOS")
#'
#' @return Tibble with pooled CDF:
#'   - particle_size_um: Common size grid
#'   - q3_cdf_mean: Average CDF across replicates (0-1)
#'   - q3_cdf_sd: Standard deviation across replicates
#'   - n_replicates: Number of replicates pooled
#'
pool_replicate_cdfs <- function(data, formulation, module, device_resistance = NULL, pressure_drop = NULL) {

  # Start with base filter
  pooled_cdf <- data %>%
    filter(
      formulation == !!formulation,
      module == !!module
    )

  # Add device/pressure filters ONLY for INHALER data
  if (!is.null(device_resistance)) {
    pooled_cdf <- pooled_cdf %>%
      filter(device_resistance == !!device_resistance)
  }

  if (!is.null(pressure_drop)) {
    pooled_cdf <- pooled_cdf %>%
      filter(pressure_drop_clean == !!pressure_drop)
  }

  # Pool replicates
  pooled_cdf <- pooled_cdf %>%
    group_by(particle_size_um) %>%
    summarise(
      q3_cdf_mean = mean(q3_cdf, na.rm = TRUE),
      q3_cdf_sd = sd(q3_cdf, na.rm = TRUE),
      n_replicates = n(),
      .groups = 'drop'
    ) %>%
    arrange(particle_size_um)

  return(pooled_cdf)
}


# ==============================================================================
# CORE FUNCTION: Calculate Wasserstein-1 distance
# ==============================================================================

#' Calculate 1D Wasserstein Distance Between Two CDFs
#'
#' Computes the Wasserstein-1 (Earth Mover's) distance between two cumulative
#' distribution functions. This is the integral of the absolute difference
#' between CDFs across the particle size axis.
#'
#' Mathematical formulation:
#'   W1 = ∫ |F_test(x) - F_ref(x)| dx
#'
#' where:
#'   - F_test(x) = test condition CDF (e.g., INHALER)
#'   - F_ref(x) = reference condition CDF (e.g., RODOS)
#'   - x = particle size in micrometers
#'
#' The result has units of micrometers and represents the average particle
#' size shift needed to transform the test distribution into the reference.
#'
#' @param size_grid Numeric vector of particle sizes (µm) - must be same for both CDFs
#' @param cdf_test Numeric vector of test CDF values (0-1) at each size
#' @param cdf_ref Numeric vector of reference CDF values (0-1) at each size
#'
#' @return Numeric scalar: W1 distance in micrometers
#'
#' @examples
#' # Simple example with artificial CDFs
#' sizes <- c(1, 2, 3, 4, 5)
#' cdf_inhaler <- c(0.1, 0.3, 0.5, 0.8, 1.0)
#' cdf_rodos <- c(0.2, 0.5, 0.7, 0.9, 1.0)
#'
#' w1 <- calculate_wasserstein_1d(sizes, cdf_inhaler, cdf_rodos)
#' print(paste("W1 =", round(w1, 3), "µm"))
#'
calculate_wasserstein_1d <- function(size_grid, cdf_test, cdf_ref) {

  # Validation checks
  if (length(size_grid) != length(cdf_test) || length(size_grid) != length(cdf_ref)) {
    stop("All inputs must have the same length")
  }

  if (length(size_grid) < 2) {
    stop("Need at least 2 data points to calculate W1")
  }

  if (any(is.na(size_grid)) || any(is.na(cdf_test)) || any(is.na(cdf_ref))) {
    stop("NA values detected in inputs")
  }

  if (any(cdf_test < 0) || any(cdf_test > 1) || any(cdf_ref < 0) || any(cdf_ref > 1)) {
    stop("CDF values must be between 0 and 1")
  }

  # Ensure size grid is sorted
  sort_order <- order(size_grid)
  size_grid <- size_grid[sort_order]
  cdf_test <- cdf_test[sort_order]
  cdf_ref <- cdf_ref[sort_order]

  # Calculate absolute difference between CDFs at each size point
  cdf_diff <- abs(cdf_test - cdf_ref)

  # Trapezoidal integration: sum of |CDF difference| × size increment
  n <- length(size_grid)
  w1 <- 0

  for (i in 1:(n - 1)) {
    dx <- size_grid[i + 1] - size_grid[i]  # Size increment

    if (dx > 0) {  # Only integrate over positive intervals
      # Trapezoidal rule: average height × width
      w1 <- w1 + dx * (cdf_diff[i] + cdf_diff[i + 1]) / 2
    }
  }

  return(w1)
}


# ==============================================================================
# MAIN FUNCTION: Calculate W1 for all formulations
# ==============================================================================

#' Calculate Pairwise Wasserstein Distances Between Test and Reference
#'
#' Compares test condition (e.g., INHALER) to reference condition (e.g., RODOS)
#' for all formulations in the dataset. Properly pools replicates before
#' calculating W1 distance.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param reference_module Character string for reference condition (default: "RODOS")
#' @param test_module Character string for test condition (default: "INHALER")
#' @param output_dir Directory to save results. Default: "results/"
#' @param save_output Should results be saved to CSV? Default: TRUE
#' @param output_filename Name of output file. Default: "wasserstein_results.csv"
#' @param verbose Print progress messages? (default: TRUE)
#'
#' @return Tibble with W1 results:
#'   - formulation: Formulation identifier
#'   - W1_micrometers: Absolute W1 distance in µm
#'   - d50_reference_um: Median particle size of reference (µm)
#'   - d50_test_um: Median particle size of test (µm)
#'   - W1_normalized: W1 divided by reference d50 (dimensionless)
#'   - d50_shift_um: Test d50 minus reference d50 (µm)
#'
#' @examples
#' # After running 01_data_import.R
#' source("scripts/02_wasserstein_core.R")
#' w1_results <- calculate_pairwise_wasserstein(data)
#'
calculate_pairwise_wasserstein <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    output_dir = "results_v2",
    save_output = TRUE,
    output_filename = "wasserstein_results.csv",
    verbose = TRUE
) {

  # Create output directory if needed
  if (save_output && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) {
      cat("Created output directory:", output_dir, "\n")
    }
  }

  if (verbose) {
    cat("\n========================================================================\n")
    cat("CALCULATING WASSERSTEIN-1 DISTANCES\n")
    cat("========================================================================\n")
    cat("Reference condition:", reference_module, "\n")
    cat("Test condition:", test_module, "\n")
    cat("Methodology: Pool replicates → Calculate W1\n")
    if (save_output) {
      cat("Output file:", file.path(output_dir, output_filename), "\n")
    }
    cat("------------------------------------------------------------------------\n\n")
  }

  # Get unique formulations
  formulations <- unique(data$formulation)

  if (verbose) {
    cat("Formulations to process:", length(formulations), "\n")
    cat("Formulation IDs:", paste(formulations, collapse = ", "), "\n\n")
  }

  # Initialize results list
  results_list <- list()

  # Loop through each formulation
  for (form in formulations) {

    if (verbose) {
      cat("Processing:", form, "... ")
    }

    # Pool replicates for reference condition
    ref_pooled <- pool_replicate_cdfs(data, form, reference_module)

    # Pool replicates for test condition
    test_pooled <- pool_replicate_cdfs(data, form, test_module)

    # Check that we have data
    if (nrow(ref_pooled) == 0 || nrow(test_pooled) == 0) {
      if (verbose) {
        cat("FAILED (no data)\n")
      }
      warning(paste("No data found for", form))
      next
    }

    # Calculate W1 distance
    w1 <- calculate_wasserstein_1d(
      size_grid = ref_pooled$particle_size_um,
      cdf_test = test_pooled$q3_cdf_mean,
      cdf_ref = ref_pooled$q3_cdf_mean
    )

    # Calculate d50 values for normalization
    d50_ref <- calculate_d50(ref_pooled$particle_size_um, ref_pooled$q3_cdf_mean)
    d50_test <- calculate_d50(test_pooled$particle_size_um, test_pooled$q3_cdf_mean)

    # Store results
    results_list[[length(results_list) + 1]] <- tibble(
      formulation = form,
      W1_micrometers = w1,
      d50_reference_um = d50_ref,
      d50_test_um = d50_test,
      W1_normalized = w1 / d50_ref,
      d50_shift_um = d50_test - d50_ref
    )

    if (verbose) {
      cat(sprintf("W1 = %.4f µm, W1/d50 = %.4f\n", w1, w1 / d50_ref))
    }
  }

  # Combine all results
  if (length(results_list) == 0) {
    stop("No successful W1 calculations. Check your data.")
  }

  results <- bind_rows(results_list) %>%
    arrange(formulation)

  if (verbose) {
    cat("\n========================================================================\n")
    cat("WASSERSTEIN CALCULATION COMPLETE\n")
    cat("========================================================================\n")
    cat("Successfully calculated W1 for", nrow(results), "formulations\n\n")

    cat("Summary statistics:\n")
    cat(sprintf("  W1 range: %.4f - %.4f µm\n",
                min(results$W1_micrometers), max(results$W1_micrometers)))
    cat(sprintf("  W1/d50 range: %.4f - %.4f\n",
                min(results$W1_normalized), max(results$W1_normalized)))
    cat(sprintf("  Mean W1: %.4f µm (SD = %.4f)\n",
                mean(results$W1_micrometers), sd(results$W1_micrometers)))
    cat(sprintf("  Mean W1/d50: %.4f (SD = %.4f)\n",
                mean(results$W1_normalized), sd(results$W1_normalized)))
    cat("========================================================================\n\n")
  }

  # Save output if requested
  if (save_output) {
    output_path <- file.path(output_dir, output_filename)
    write_csv(results, output_path)

    if (verbose) {
      cat("✓ Results saved to:", output_path, "\n")
      cat("  Use W1_micrometers for DoE modeling\n\n")
    }
  }

  return(results)
}


# ==============================================================================
# HELPER FUNCTION: Calculate d50 (median particle diameter)
# ==============================================================================

#' Calculate d50 from CDF
#'
#' Finds the particle diameter where cumulative distribution equals 0.5 (50%)
#' using linear interpolation.
#'
#' @param size_grid Numeric vector of particle sizes (µm)
#' @param cdf Numeric vector of CDF values (0-1)
#'
#' @return Numeric scalar: d50 in micrometers
#'
calculate_d50 <- function(size_grid, cdf) {

  # Find where CDF crosses 0.5
  idx <- which(cdf >= 0.5)[1]

  if (is.na(idx)) {
    warning("CDF never reaches 0.5 - check data")
    return(NA)
  }

  # If exactly 0.5, return that size
  if (cdf[idx] == 0.5) {
    return(size_grid[idx])
  }

  # If first point is already > 0.5, return first size
  if (idx == 1) {
    return(size_grid[1])
  }

  # Linear interpolation between bracketing points
  x1 <- size_grid[idx - 1]
  x2 <- size_grid[idx]
  y1 <- cdf[idx - 1]
  y2 <- cdf[idx]

  # Interpolate to find exact d50
  d50 <- x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)

  return(d50)
}


# ==============================================================================
# VALIDATION FUNCTION: Check W1 results quality
# ==============================================================================

#' Validate Wasserstein Distance Results
#'
#' Checks W1 results for common issues and provides diagnostic information
#'
#' @param w1_results Tibble from calculate_pairwise_wasserstein()
#' @param max_w1_normalized Maximum reasonable W1/d50 value (default: 2.0)
#'   Values above this threshold are flagged as potentially problematic
#'
#' @return Invisibly returns TRUE if validation passes
#'
validate_wasserstein_results <- function(w1_results, max_w1_normalized = 2.0) {

  cat("\n========================================================================\n")
  cat("VALIDATING WASSERSTEIN RESULTS\n")
  cat("========================================================================\n")

  all_valid <- TRUE

  # Check for NA values
  na_counts <- w1_results %>%
    summarise(across(everything(), ~sum(is.na(.))))

  if (any(na_counts > 0)) {
    cat("\nWARNING: NA values detected:\n")
    print(na_counts %>% select(where(~. > 0)))
    all_valid <- FALSE
  } else {
    cat("✓ No NA values in results\n")
  }

  # Check for negative W1 (impossible)
  if (any(w1_results$W1_micrometers < 0, na.rm = TRUE)) {
    cat("\nERROR: Negative W1 values detected (this is impossible!)\n")
    all_valid <- FALSE
  } else {
    cat("✓ All W1 values are positive\n")
  }

  # Check for unreasonably large W1/d50
  large_w1 <- w1_results %>%
    filter(W1_normalized > max_w1_normalized)

  if (nrow(large_w1) > 0) {
    cat(sprintf("\nWARNING: Some W1/d50 values exceed %.2f:\n", max_w1_normalized))
    print(large_w1 %>% select(formulation, W1_normalized, d50_reference_um))
    cat("  → This suggests very poor dispersibility or data quality issues\n")
    all_valid <- FALSE
  } else {
    cat(sprintf("✓ All W1/d50 values are reasonable (< %.2f)\n", max_w1_normalized))
  }

  # Check for negative d50 shifts (test finer than reference - unusual)
  negative_shifts <- w1_results %>%
    filter(d50_shift_um < 0)

  if (nrow(negative_shifts) > 0) {
    cat("\nNOTE: Some formulations show negative d50 shifts (test finer than reference):\n")
    print(negative_shifts %>% select(formulation, d50_shift_um))
    cat("  → This is unusual but possible if test dispersion is more efficient\n")
  }

  # Summary statistics
  cat("\n------------------------------------------------------------------------\n")
  cat("DISTRIBUTION SUMMARY:\n")
  cat(sprintf("  W1: %.4f ± %.4f µm (range: %.4f - %.4f)\n",
              mean(w1_results$W1_micrometers, na.rm = TRUE),
              sd(w1_results$W1_micrometers, na.rm = TRUE),
              min(w1_results$W1_micrometers, na.rm = TRUE),
              max(w1_results$W1_micrometers, na.rm = TRUE)))

  cat(sprintf("  W1/d50: %.4f ± %.4f (range: %.4f - %.4f)\n",
              mean(w1_results$W1_normalized, na.rm = TRUE),
              sd(w1_results$W1_normalized, na.rm = TRUE),
              min(w1_results$W1_normalized, na.rm = TRUE),
              max(w1_results$W1_normalized, na.rm = TRUE)))

  cat("------------------------------------------------------------------------\n")
  if (all_valid) {
    cat("VALIDATION PASSED: Results look reasonable\n")
  } else {
    cat("VALIDATION WARNINGS: Review issues above\n")
  }
  cat("========================================================================\n\n")

  return(invisible(all_valid))
}


# ==============================================================================
# CONVENIENCE FUNCTIONS: Load and run with defaults
# ==============================================================================

#' Load Previously Calculated Wasserstein Results
#'
#' Quickly loads the W1 results saved by calculate_pairwise_wasserstein()
#'
#' @param results_dir Directory containing results. Default: "results/"
#' @param filename Name of results file. Default: "wasserstein_results.csv"
#' @param verbose Print loading message? Default: TRUE
#'
#' @return Tibble with W1 results
#'
#' @examples
#' # Load previously calculated W1 distances
#' w1_results <- load_wasserstein_results()
#'
load_wasserstein_results <- function(
    results_dir = "results",
    filename = "wasserstein_results.csv",
    verbose = TRUE
) {

  file_path <- file.path(results_dir, filename)

  if (!file.exists(file_path)) {
    stop("Results file not found: ", file_path,
         "\nRun calculate_pairwise_wasserstein() first to create this file.")
  }

  if (verbose) {
    cat("Loading W1 results from:", file_path, "\n")
  }

  results <- read_csv(file_path, show_col_types = FALSE)

  if (verbose) {
    cat("✓ Loaded results for", nrow(results), "formulations\n")
    cat("  W1 range:", sprintf("%.4f - %.4f µm\n",
                              min(results$W1_micrometers),
                              max(results$W1_micrometers)))
  }

  return(results)
}


#' Run Complete Wasserstein Analysis with Auto-Load
#'
#' Convenience wrapper that:
#' 1. Auto-loads data/tidy/standardized_data.csv
#' 2. Calculates W1 distances
#' 3. Auto-saves to results/wasserstein_results.csv
#' 4. Validates results
#'
#' @param processed_dir Directory with tidy data. Default: "data/tidy/"
#' @param results_dir Directory to save results. Default: "results/"
#' @param reference_module Reference condition name. Default: "RODOS"
#' @param test_module Test condition name. Default: "INHALER"
#' @param verbose Print progress? Default: TRUE
#'
#' @return Tibble with W1 results
#'
#' @examples
#' # Simple one-liner after running 01_data_import.R
#' source("scripts/02_wasserstein_core.R")
#' w1_results <- run_wasserstein_analysis()
#'
run_wasserstein_analysis <- function(
    processed_dir = "data_v2/tidy",
    results_dir = "results_v2",
    reference_module = "RODOS",
    test_module = "INHALER",
    verbose = TRUE
) {

  # Load processed data
  data_file <- file.path(processed_dir, "standardized_data.csv")

  if (!file.exists(data_file)) {
    stop("Processed data not found: ", data_file,
         "\nRun 01_data_import.R first to create this file.")
  }

  if (verbose) {
    cat("Loading processed data from:", data_file, "\n")
  }

  data <- read_csv(data_file, show_col_types = FALSE)

  if (verbose) {
    cat("✓ Loaded", nrow(data), "rows\n")
    cat("  Formulations:", n_distinct(data$formulation), "\n")
    cat("  Modules:", paste(unique(data$module), collapse = ", "), "\n")
  }

  # Calculate W1 distances
  w1_results <- calculate_pairwise_wasserstein(
    data,
    reference_module = reference_module,
    test_module = test_module,
    output_dir = results_dir,
    save_output = TRUE,
    verbose = verbose
  )

  # Validate results
  validate_wasserstein_results(w1_results)

  return(w1_results)
}


# ==============================================================================
# AUTO-EXECUTION: Run analysis when script is sourced
# ==============================================================================

# Check if processed data exists
if (file.exists("data_v2/tidy/standardized_data_with_conditions.csv")) {

  cat("\n========================================================================\n")
  cat("AUTO-RUNNING WASSERSTEIN ANALYSIS\n")
  cat("========================================================================\n")
  cat("Reading: data_v2/tidy/standardized_data_with_conditions.csv\n")
  cat("Saving to: results_v2/wasserstein_results.csv\n")
  cat("------------------------------------------------------------------------\n")

  # Run the complete analysis
  .w1_results <- run_wasserstein_analysis(verbose = TRUE)

  cat("\n========================================================================\n")
  cat("ANALYSIS COMPLETE - Results saved to results/wasserstein_results.csv\n")
  cat("========================================================================\n")
  cat("Use W1_micrometers column for DoE modeling\n")
  cat("------------------------------------------------------------------------\n")
  cat("To reload results later:\n")
  cat("  source('scripts/02_wasserstein_core.R')\n")
  cat("  w1_results <- load_wasserstein_results()\n")
  cat("========================================================================\n\n")

  # Clean up the auto-generated variable (optional - keeps workspace clean)
  # Uncomment if you don't want .w1_results in the environment
  # rm(.w1_results)

} else {
  cat("\n========================================================================\n")
  cat("WASSERSTEIN ANALYSIS - WAITING FOR INPUT DATA\n")
  cat("========================================================================\n")
  cat("Tidy data not found: data_v2/tidy/standardized_data_with_conditions.csv\n")
  cat("\nPlease run 01_data_import.R first:\n")
  cat("  source('scripts/01_data_import.R')\n")
  cat("\nThen run this script again:\n")
  cat("  source('scripts/02_wasserstein_core.R')\n")
  cat("========================================================================\n\n")
}
