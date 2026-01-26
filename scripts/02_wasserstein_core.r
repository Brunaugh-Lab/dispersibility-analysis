# ==============================================================================
# 02_wasserstein_core.R
# Core Wasserstein-1 Distance Calculations for Dispersibility Analysis
#
# Purpose: Calculate W1 distance between test and reference particle size
#          distributions following proper methodology (pool replicates first)
#
# Methodology: Wasserstein-1 (Earth Mover's) distance quantifies the minimum
#              redistribution work needed to transform one distribution into another
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
pool_replicate_cdfs <- function(data, formulation, module) {

  pooled_cdf <- data %>%
    filter(
      formulation == !!formulation,
      module == !!module
    ) %>%
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
#' # Basic usage after importing data
#' source("scripts/01_data_import.R")
#' data <- read_ld_data_from_structure("data/")
#'
#' # Calculate W1 distances
#' w1_results <- calculate_pairwise_wasserstein(
#'   data,
#'   reference_module = "RODOS",
#'   test_module = "INHALER"
#' )
#'
calculate_pairwise_wasserstein <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    verbose = TRUE
) {

  if (verbose) {
    cat("\n========================================================================\n")
    cat("CALCULATING WASSERSTEIN-1 DISTANCES\n")
    cat("========================================================================\n")
    cat("Reference condition:", reference_module, "\n")
    cat("Test condition:", test_module, "\n")
    cat("Methodology: Pool replicates → Calculate W1\n")
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
# EXAMPLE USAGE
# ==============================================================================

# Uncomment to test with your data:
#
# # Step 1: Load data import script
# source("scripts/01_data_import.R")
#
# # Step 2: Import data
# data <- read_ld_data_from_structure(
#   data_directory = "data/",
#   formulation_pattern = ".*",
#   replicate_pattern = "[Rr]ep_?\\d+",
#   verbose = TRUE
# )
#
# # Step 3: Standardize replicate names
# data <- data %>% mutate(replicate = tolower(replicate))
#
# # Step 4: Calculate Wasserstein distances
# w1_results <- calculate_pairwise_wasserstein(
#   data,
#   reference_module = "RODOS",
#   test_module = "INHALER",
#   verbose = TRUE
# )
#
# # Step 5: Validate results
# validate_wasserstein_results(w1_results)
#
# # Step 6: View results
# print(w1_results, n = Inf)
#
# # Step 7: Save results
# write_csv(w1_results, "wasserstein_results.csv")
