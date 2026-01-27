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
