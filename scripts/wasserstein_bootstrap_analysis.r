# ==============================================================================
# 04_bootstrap_analysis.R
# Bootstrap Resampling for Wasserstein Distance Confidence Intervals
#
# Purpose: Generate empirical sampling distributions for W1 distances through
#          bootstrap resampling of replicate reference and inhaler distributions with replacement
#          Pool or average across resampled triplicate Q3 values at each particle size to get
#          RODOS + INHALER empirical distributions.
#          For each RODOS + INHALER distribution pair, the W1 distance is calculated.
#          Calculate confidence intervals and effect-to-noise ratios to
#          quantify measurement uncertainty vs true condition effects.
#
# Auto-execution: Script automatically runs when sourced
#   - Reads data/tidy/standardized_data.csv (from script 01)
#   - Performs bootstrap resampling (default: 2000 iterations)
#   - Calculates 95% CIs and standard errors for each formulation
#   - Computes effect-to-noise ratios for device resistance/pressure effects
#   - Saves to results/bootstrap_results.csv and results/effect_noise_ratios.csv
#
# Input: data/tidy/standardized_data.csv (from 01_data_import.R)
# Output:
#   - results/bootstrap_results.csv (W1 distributions with CIs)
#   - results/effect_noise_ratios.csv (signal vs noise quantification)
#
# Methodology:
#   1. Pool replicates → create empirical particle size distributions
#   2. Bootstrap resample with replacement from pooled distributions
#   3. Calculate W1 for each bootstrap iteration
#   4. Estimate sampling distribution parameters (mean, SD, 95% CI)
#   5. Compare bootstrap variability to between-condition effects
#
# Usage:
#   source("scripts/wassserstein_bootstrap_analysis.R")
#
# ==============================================================================
library(tidyverse)
library(broom)

# ==============================================================================
# HELPER FUNCTION: Get individual replicate CDFs for a condition
# ==============================================================================

#' Extract Individual Replicate CDFs for Bootstrap Resampling
#'
#' Retrieves the individual replicate CDFs (before averaging) for a specific
#' formulation and module combination. These will be resampled during bootstrap.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param formulation Character string specifying formulation
#' @param module Character string specifying module ("INHALER" or "RODOS")
#' @param common_size_grid Optional common size grid for interpolation
#'
#' @return List of CDF vectors (one per replicate)
#'
get_replicate_cdfs <- function(data, formulation, module, common_size_grid = NULL) {

  # Filter data for this condition
  condition_data <- data %>%
    filter(
      formulation == !!formulation,
      module == !!module
    )

  if (nrow(condition_data) == 0) {
    warning(sprintf("No data found for formulation %s, module %s", formulation, module))
    return(list())
  }

  # Get unique replicates
  replicates <- unique(condition_data$replicate)

  # If no common grid provided, use union of all size points
  if (is.null(common_size_grid)) {
    common_size_grid <- sort(unique(condition_data$particle_size_um))
  }

  # Extract CDF for each replicate, interpolated to common grid
  replicate_cdfs <- map(replicates, function(rep) {
    rep_data <- condition_data %>%
      filter(replicate == !!rep) %>%
      arrange(particle_size_um)

    # Interpolate to common grid
    if (length(rep_data$particle_size_um) > 1) {
      approx(
        rep_data$particle_size_um,
        rep_data$q3_cdf,
        xout = common_size_grid,
        rule = 2  # Use endpoint values for extrapolation
      )$y
    } else {
      # Single point - not enough for interpolation
      rep(rep_data$q3_cdf[1], length(common_size_grid))
    }
  })

  names(replicate_cdfs) <- replicates

  return(list(
    cdfs = replicate_cdfs,
    size_grid = common_size_grid,
    replicates = replicates
  ))
}
# ==============================================================================
# HELPER FUNCTION: Bootstrap resample replicate CDFs and pool
# ==============================================================================

#' Bootstrap Resample Replicate CDFs and Create Pooled CDF
#'
#' Samples replicate CDFs with replacement and averages them to create
#' a bootstrap version of the pooled CDF, following the same methodology
#' as the original analysis.
#'
#' @param replicate_cdfs_info List from get_replicate_cdfs()
#' @param n_replicates Number of replicates to sample (default: original count)
#'
#' @return Vector of pooled CDF values
#'
bootstrap_pool_replicates <- function(replicate_cdfs_info, n_replicates = NULL) {

  cdfs <- replicate_cdfs_info$cdfs

  if (length(cdfs) == 0) {
    return(rep(NA_real_, length(replicate_cdfs_info$size_grid)))
  }

  # Default to original number of replicates
  if (is.null(n_replicates)) {
    n_replicates <- length(cdfs)
  }

  # Sample replicates with replacement
  sampled_indices <- sample(seq_along(cdfs), size = n_replicates, replace = TRUE)

  # Average the selected CDFs
  sampled_cdfs <- cdfs[sampled_indices]
  pooled_cdf <- Reduce("+", sampled_cdfs) / length(sampled_cdfs)

  return(pooled_cdf)
}


# ==============================================================================
# CORE FUNCTION: Bootstrap W1 calculation for single formulation
# ==============================================================================

#' Bootstrap Wasserstein Distance for Single Formulation
#'
#' Performs bootstrap resampling by selecting replicate CDFs with replacement,
#' pooling them (averaging), and calculating W1. This captures uncertainty
#' in the pooling process due to which specific replicates were measured.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param formulation Character string specifying formulation
#' @param reference_module Reference condition (default: "RODOS")
#' @param test_module Test condition (default: "INHALER")
#' @param n_bootstrap Number of bootstrap iterations (default: 2000)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages
#'
#' @return Tibble with bootstrap results:
#'   - formulation: Formulation identifier
#'   - w1_mean: Mean W1 across bootstrap iterations
#'   - w1_sd: Standard deviation (standard error)
#'   - w1_ci_lower: Lower 95% confidence interval
#'   - w1_ci_upper: Upper 95% confidence interval
#'   - w1_observed: Original W1 from pooled data
#'   - bootstrap_samples: List column with all W1 values
#'
bootstrap_w1_single <- function(data, formulation, reference_module = "RODOS",
                               test_module = "INHALER", n_bootstrap = 2000,
                               seed = NULL, verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  if (verbose) {
    cat(sprintf("Bootstrap resampling: %s (%d iterations)\n", formulation, n_bootstrap))
  }

  # Load calculate_wasserstein_1d function from script 02
  if (!exists("calculate_wasserstein_1d")) {
    source("/mnt/user-data/uploads/02_wasserstein_core.r", local = TRUE)
  }

  # Get individual replicate CDFs for both conditions
  ref_replicates <- get_replicate_cdfs(data, formulation, reference_module)
  test_replicates <- get_replicate_cdfs(data, formulation, test_module)

  # Check that we have data for both conditions
  if (length(ref_replicates$cdfs) == 0 || length(test_replicates$cdfs) == 0) {
    warning(sprintf("Missing replicate data for formulation %s", formulation))
    return(tibble(
      formulation = formulation,
      w1_mean = NA_real_, w1_sd = NA_real_,
      w1_ci_lower = NA_real_, w1_ci_upper = NA_real_,
      w1_observed = NA_real_, n_bootstrap = n_bootstrap,
      n_ref_replicates = length(ref_replicates$cdfs),
      n_test_replicates = length(test_replicates$cdfs),
      bootstrap_samples = list(numeric(0))
    ))
  }

  # Create common size grid (union of both conditions)
  common_sizes <- sort(unique(c(ref_replicates$size_grid, test_replicates$size_grid)))

  # Re-extract with common grid
  ref_replicates <- get_replicate_cdfs(data, formulation, reference_module, common_sizes)
  test_replicates <- get_replicate_cdfs(data, formulation, test_module, common_sizes)

  # Calculate observed W1 from normally pooled data (average of all replicates)
  ref_pooled_observed <- Reduce("+", ref_replicates$cdfs) / length(ref_replicates$cdfs)
  test_pooled_observed <- Reduce("+", test_replicates$cdfs) / length(test_replicates$cdfs)
  w1_observed <- calculate_wasserstein_1d(common_sizes, test_pooled_observed, ref_pooled_observed)

  # Bootstrap resampling
  w1_bootstrap <- numeric(n_bootstrap)

  for (i in seq_len(n_bootstrap)) {
    # Bootstrap sample and pool reference condition
    ref_pooled_boot <- bootstrap_pool_replicates(ref_replicates)

    # Bootstrap sample and pool test condition
    test_pooled_boot <- bootstrap_pool_replicates(test_replicates)

    # Calculate W1 for this bootstrap iteration
    w1_bootstrap[i] <- calculate_wasserstein_1d(common_sizes, test_pooled_boot, ref_pooled_boot)
  }

  # Calculate summary statistics
  w1_mean <- mean(w1_bootstrap, na.rm = TRUE)
  w1_sd <- sd(w1_bootstrap, na.rm = TRUE)
  w1_ci <- quantile(w1_bootstrap, c(0.025, 0.975), na.rm = TRUE)

  # Return results
  tibble(
    formulation = formulation,
    w1_mean = w1_mean,
    w1_sd = w1_sd,
    w1_ci_lower = w1_ci[[1]],
    w1_ci_upper = w1_ci[[2]],
    w1_observed = w1_observed,
    n_bootstrap = n_bootstrap,
    n_ref_replicates = length(ref_replicates$cdfs),
    n_test_replicates = length(test_replicates$cdfs),
    bootstrap_samples = list(w1_bootstrap)
  )
}
 ==============================================================================
# MAIN FUNCTION: Bootstrap analysis for all formulations
# ==============================================================================

#' Bootstrap Wasserstein Analysis for All Formulations
#'
#' Performs bootstrap resampling analysis for all formulations in dataset by
#' resampling replicate CDFs with replacement, pooling, and calculating W1.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param reference_module Reference condition (default: "RODOS")
#' @param test_module Test condition (default: "INHALER")
#' @param n_bootstrap Number of bootstrap iterations per formulation (default: 2000)
#' @param output_dir Directory to save results (default: "results")
#' @param save_output Should results be saved? (default: TRUE)
#' @param output_filename Output filename (default: "bootstrap_results.csv")
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Tibble with bootstrap results for all formulations
#'
bootstrap_w1_analysis <- function(data, reference_module = "RODOS",
                                 test_module = "INHALER", n_bootstrap = 2000,
                                 output_dir = "results", save_output = TRUE,
                                 output_filename = "bootstrap_results.csv",
                                 seed = 42, verbose = TRUE) {

  # Create output directory if needed
  if (save_output && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", output_dir, "\n")
  }

  if (verbose) {
    cat("\n========================================================================\n")
    cat("BOOTSTRAP WASSERSTEIN ANALYSIS\n")
    cat("========================================================================\n")
    cat("Bootstrap iterations:", n_bootstrap, "\n")
    cat("Method: Resample replicate CDFs → Pool → Calculate W1\n")
    cat("Reference condition:", reference_module, "\n")
    cat("Test condition:", test_module, "\n")
    if (save_output) {
      cat("Output file:", file.path(output_dir, output_filename), "\n")
    }
    cat("------------------------------------------------------------------------\n\n")
  }

  # Get unique formulations
  formulations <- unique(data$formulation)
  n_formulations <- length(formulations)

  if (verbose) {
    cat("Formulations to process:", n_formulations, "\n")
    cat("Total bootstrap samples:", n_formulations * n_bootstrap, "\n\n")
  }

  # Run bootstrap analysis for each formulation
  bootstrap_results <- map_dfr(formulations, function(form) {
    bootstrap_w1_single(
      data = data,
      formulation = form,
      reference_module = reference_module,
      test_module = test_module,
      n_bootstrap = n_bootstrap,
      seed = seed,
      verbose = verbose
    )
  })

  # Add coefficient of variation and other metrics
  bootstrap_results <- bootstrap_results %>%
    mutate(
      w1_cv = w1_sd / w1_mean,  # Coefficient of variation
      w1_ci_width = w1_ci_upper - w1_ci_lower,  # CI width
      w1_relative_se = w1_sd / w1_observed  # Relative standard error
    )

  if (verbose) {
    cat("\n------------------------------------------------------------------------\n")
    cat("BOOTSTRAP SUMMARY:\n")
    cat(sprintf("  Mean W1: %.4f ± %.4f µm\n",
                mean(bootstrap_results$w1_mean, na.rm = TRUE),
                sd(bootstrap_results$w1_mean, na.rm = TRUE)))
    cat(sprintf("  Mean SE: %.4f µm (%.1f%% CV)\n",
                mean(bootstrap_results$w1_sd, na.rm = TRUE),
                100 * mean(bootstrap_results$w1_cv, na.rm = TRUE)))
    cat(sprintf("  Mean CI width: %.4f µm\n",
                mean(bootstrap_results$w1_ci_width, na.rm = TRUE)))
    cat("------------------------------------------------------------------------\n")
  }

  # Save results
  if (save_output) {
    # Save main results (without bootstrap samples to keep file size reasonable)
    bootstrap_summary <- bootstrap_results %>%
      select(-bootstrap_samples)

    output_path <- file.path(output_dir, output_filename)
    write_csv(bootstrap_summary, output_path)

    if (verbose) {
      cat("✓ Bootstrap results saved to:", output_path, "\n")
    }
  }

  return(bootstrap_results)
}
