# ==============================================================================
# wasserstein_bootstrap_analysis.R
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
#   - Saves to results_v2/bootstrap_results.csv and results_v2/effect_noise_ratios.csv
#
# Input: data_v2/tidy/standardized_data.csv (from 01_data_import.R)
# Output:
#   - results_v2/bootstrap_results.csv (W1 distributions with CIs)
#   - results_v2/effect_noise_ratios.csv (signal vs noise quantification)
#
# Methodology:
#   1. Pool replicates → create empirical particle size distributions
#   2. Bootstrap resample with replacement from pooled distributions
#   3. Calculate W1 for each bootstrap iteration
#   4. Estimate sampling distribution parameters (mean, SD, 95% CI)
#   5. Compare bootstrap variability to between-condition effects
#
# Usage:
#   source("scripts/wasserstein_bootstrap_analysis.R")
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
get_replicate_cdfs <- function(data, formulation, module,
                               device_resistance = NULL,
                               pressure_drop = NULL,
                               common_size_grid = NULL) {

  # Filter data for this condition
  condition_data <- data %>%
    filter(
      formulation == !!formulation,
      module == !!module
    )

  # Filter by device resistance if specified (for INHALER conditions)
  if (!is.null(device_resistance)) {
    condition_data <- condition_data %>%
      filter(device_resistance == !!device_resistance)
  }

  # Filter by pressure drop if specified (for INHALER conditions)
  if (!is.null(pressure_drop)) {
    condition_data <- condition_data %>%
      filter(pressure_drop_clean == !!pressure_drop)
  }

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
                               test_module = "INHALER", device_resistance = NULL,
                              pressure_drop = NULL, n_bootstrap = 2000,
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
  test_replicates <- get_replicate_cdfs(data, formulation, test_module, device_resistance = device_resistance, pressure_drop = pressure_drop)

  # Check that we have data for both conditions
  if (length(ref_replicates$cdfs) == 0 || length(test_replicates$cdfs) == 0) {
    warning(sprintf("Missing replicate data for formulation %s", formulation))
    return(tibble(
      formulation = formulation,
      device_resistance = device_resistance %||% NA_character_,
      pressure_drop = pressure_drop %||% NA_character_,
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

# ==============================================================================
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
                                 output_dir = "results_v2", save_output = TRUE,
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
# ==============================================================================
# FUNCTION: Calculate effect-to-noise ratios
# ==============================================================================

#' Calculate Overall Effect-to-Noise Ratio
#'
#' Compares the overall variability of W1 values across all formulations
#' to the average bootstrap-estimated measurement uncertainty. Higher ratios
#' indicate that formulation differences exceed measurement variability.
#'
#' @param bootstrap_results Results from bootstrap_w1_analysis()
#' @param output_dir Directory to save results (default: "results")
#' @param save_output Should results be saved? (default: TRUE)
#' @param output_filename Output filename (default: "effect_noise_ratios.csv")
#' @param verbose Print progress (default: TRUE)
#'
#' @return Tibble with effect-to-noise ratio
#'
calculate_effect_noise_ratios <- function(bootstrap_results,
                                         output_dir = "results",
                                         save_output = TRUE,
                                         output_filename = "effect_noise_ratios.csv",
                                         verbose = TRUE) {

  if (verbose) {
    cat("\n========================================================================\n")
    cat("EFFECT-TO-NOISE RATIO ANALYSIS\n")
    cat("========================================================================\n")
    cat("Analysis: Overall formulation variability vs measurement uncertainty\n")
  }

  # Calculate effect magnitude (between-formulation variability)
  effect_magnitude <- sd(bootstrap_results$w1_observed, na.rm = TRUE)

  # Calculate noise level (average measurement uncertainty)
  noise_level <- mean(bootstrap_results$w1_sd, na.rm = TRUE)

  # Calculate ratio
  effect_to_noise_ratio <- effect_magnitude / noise_level

  # Create results
  ratio_results <- tibble(
    analysis_type = "overall_formulation_variability",
    effect_magnitude_um = effect_magnitude,
    noise_level_um = noise_level,
    effect_to_noise_ratio = effect_to_noise_ratio,
    n_formulations = nrow(bootstrap_results),
    interpretation = case_when(
      effect_to_noise_ratio >= 3 ~ "Strong signal: Formulation differences >> measurement noise",
      effect_to_noise_ratio >= 2 ~ "Moderate signal: Formulation differences > measurement noise",
      effect_to_noise_ratio >= 1 ~ "Weak signal: Formulation differences ~ measurement noise",
      TRUE ~ "Poor signal: Formulation differences < measurement noise"
    ),
    # Additional metrics
    mean_w1_um = mean(bootstrap_results$w1_observed, na.rm = TRUE),
    min_w1_um = min(bootstrap_results$w1_observed, na.rm = TRUE),
    max_w1_um = max(bootstrap_results$w1_observed, na.rm = TRUE),
    cv_between_formulations = effect_magnitude / mean(bootstrap_results$w1_observed, na.rm = TRUE),
    mean_relative_uncertainty = mean(bootstrap_results$w1_sd / bootstrap_results$w1_observed, na.rm = TRUE)
  )

  if (verbose) {
    cat("\nRESULTS:\n")
    cat(sprintf("  Effect magnitude (between-formulation SD): %.4f µm\n", effect_magnitude))
    cat(sprintf("  Noise level (mean bootstrap SE): %.4f µm\n", noise_level))
    cat(sprintf("  Effect-to-noise ratio: %.2f\n", effect_to_noise_ratio))
    cat(sprintf("  Number of formulations: %d\n", nrow(bootstrap_results)))
    cat("\nINTERPRETATION:\n")
    cat(sprintf("  %s\n", ratio_results$interpretation))
    cat("\nADDITIONAL METRICS:\n")
    cat(sprintf("  W1 range: %.4f - %.4f µm (mean: %.4f µm)\n",
                ratio_results$min_w1_um, ratio_results$max_w1_um, ratio_results$mean_w1_um))
    cat(sprintf("  CV between formulations: %.1f%%\n", 100 * ratio_results$cv_between_formulations))
    cat(sprintf("  Mean relative uncertainty: %.1f%%\n", 100 * ratio_results$mean_relative_uncertainty))
    cat("------------------------------------------------------------------------\n")
  }

  # Save results
  if (save_output) {
    output_path <- file.path(output_dir, output_filename)
    write_csv(ratio_results, output_path)

    if (verbose) {
      cat("✓ Effect-to-noise ratio saved to:", output_path, "\n")
    }
  }

  return(ratio_results)
}

# ==============================================================================
# FUNCTION: Visualization of bootstrap results
# ==============================================================================

#' Plot Bootstrap Distributions and Confidence Intervals
#'
#' Creates diagnostic plots for bootstrap analysis results.
#'
#' @param bootstrap_results Results from bootstrap_w1_analysis()
#' @param output_dir Directory to save plots (default: "figures")
#' @param save_plots Should plots be saved? (default: TRUE)
#' @param verbose Print progress (default: TRUE)
#'
#' @return List of ggplot objects
#'
plot_bootstrap_results <- function(bootstrap_results, output_dir = "figures",
                                  save_plots = TRUE, verbose = TRUE) {

  library(ggplot2)
  library(patchwork)

  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Plot 1: Confidence intervals
  p1 <- bootstrap_results %>%
    mutate(formulation = fct_reorder(formulation, w1_observed)) %>%
    ggplot(aes(x = formulation)) +
    geom_errorbar(aes(ymin = w1_ci_lower, ymax = w1_ci_upper),
                  width = 0.3, alpha = 0.7) +
    geom_point(aes(y = w1_observed), size = 3, color = "red") +
    geom_point(aes(y = w1_mean), size = 2, color = "blue", alpha = 0.7) +
    labs(
      x = "Formulation",
      y = "Wasserstein Distance (µm)",
      title = "Bootstrap Confidence Intervals for W1",
      subtitle = "Red: Observed W1, Blue: Bootstrap Mean, Bars: 95% CI"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  # Plot 2: Standard errors
  p2 <- bootstrap_results %>%
    mutate(formulation = fct_reorder(formulation, w1_sd)) %>%
    ggplot(aes(x = formulation, y = w1_sd)) +
    geom_col(fill = "skyblue", alpha = 0.7) +
    labs(
      x = "Formulation",
      y = "Bootstrap Standard Error (µm)",
      title = "Measurement Uncertainty by Formulation"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  # Plot 3: Relative standard error
  p3 <- bootstrap_results %>%
    mutate(formulation = fct_reorder(formulation, w1_relative_se)) %>%
    ggplot(aes(x = formulation, y = 100 * w1_relative_se)) +
    geom_col(fill = "lightcoral", alpha = 0.7) +
    geom_hline(yintercept = c(5, 10, 20), linetype = "dashed", alpha = 0.5) +
    labs(
      x = "Formulation",
      y = "Relative Standard Error (%)",
      title = "Measurement Precision by Formulation",
      subtitle = "Lower values indicate more precise measurements"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  # Combine plots
  combined <- (p1 / p2 / p3) +
    plot_annotation(
      title = "Bootstrap Analysis Summary",
      tag_levels = 'A'
    )

  plots <- list(confidence_intervals = p1, standard_errors = p2,
                relative_errors = p3, combined = combined)

  # Save plots
  if (save_plots) {
    output_path <- file.path(output_dir, "bootstrap_analysis.pdf")
    ggsave(output_path, combined, width = 12, height = 10, device = "pdf")

    if (verbose) {
      cat("✓ Bootstrap plots saved to:", output_path, "\n")
    }
  }

  return(plots)
}

# ==============================================================================
# FUNCTION: Effect-to-noise ratio visualization
# ==============================================================================

#' Plot Effect-to-Noise Ratios
#'
#' Creates visualizations showing the relationship between formulation
#' variability and bootstrap measurement uncertainty.
#'
#' @param bootstrap_results Results from bootstrap_w1_analysis()
#' @param effect_noise_results Results from calculate_effect_noise_ratios()
#' @param output_dir Directory to save plots (default: "figures")
#' @param save_plots Should plots be saved? (default: TRUE)
#' @param verbose Print progress (default: TRUE)
#'
#' @return List of ggplot objects
#'
plot_effect_noise_analysis <- function(bootstrap_results, effect_noise_results,
                                      output_dir = "figures", save_plots = TRUE,
                                      verbose = TRUE) {

  library(ggplot2)
  library(patchwork)

  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Plot 1: W1 values with uncertainty bars
  p1 <- bootstrap_results %>%
    mutate(formulation = fct_reorder(formulation, w1_observed)) %>%
    ggplot(aes(x = formulation)) +
    geom_errorbar(aes(ymin = w1_ci_lower, ymax = w1_ci_upper),
                  width = 0.3, alpha = 0.7, color = "gray60") +
    geom_point(aes(y = w1_observed), size = 3, color = "steelblue") +
    labs(
      x = "Formulation",
      y = "Wasserstein Distance (µm)",
      title = "W1 Values Across Formulations",
      subtitle = "Points: Observed W1, Bars: 95% Bootstrap CI"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  # Plot 2: Effect vs Noise comparison
  effect_mag <- effect_noise_results$effect_magnitude_um
  noise_level <- effect_noise_results$noise_level_um
  ratio <- effect_noise_results$effect_to_noise_ratio

  p2 <- tibble(
    metric = c("Formulation\nVariability\n(Effect)", "Measurement\nUncertainty\n(Noise)"),
    value = c(effect_mag, noise_level),
    color = c("Effect", "Noise")
  ) %>%
    ggplot(aes(x = metric, y = value, fill = color)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.4f µm", value)),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(
      values = c("Effect" = "darkgreen", "Noise" = "orange"),
      guide = "none"
    ) +
    labs(
      x = "",
      y = "Standard Deviation (µm)",
      title = sprintf("Effect-to-Noise Ratio = %.2f", ratio),
      subtitle = effect_noise_results$interpretation
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.x = element_text(size = 11)
    )

  # Plot 3: Relative uncertainty by formulation
  p3 <- bootstrap_results %>%
    mutate(
      formulation = fct_reorder(formulation, w1_sd / w1_observed),
      relative_se_pct = 100 * w1_sd / w1_observed
    ) %>%
    ggplot(aes(x = formulation, y = relative_se_pct)) +
    geom_col(fill = "lightcoral", alpha = 0.7) +
    geom_hline(yintercept = c(5, 10, 20), linetype = "dashed", alpha = 0.5) +
    labs(
      x = "Formulation",
      y = "Relative Standard Error (%)",
      title = "Measurement Precision by Formulation",
      subtitle = "Lower values indicate more precise W1 estimates"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  # Plot 4: Signal strength interpretation
  interpretation_color <- case_when(
    ratio >= 3 ~ "darkgreen",
    ratio >= 2 ~ "orange",
    ratio >= 1 ~ "gold",
    TRUE ~ "red"
  )

  p4 <- tibble(
    x = 1, y = 1,
    ratio = ratio,
    interpretation = effect_noise_results$interpretation
  ) %>%
    ggplot(aes(x, y)) +
    geom_point(size = 50, color = interpretation_color, alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", ratio)),
              size = 8, fontface = "bold", color = "white") +
    labs(
      title = "Signal Strength Assessment",
      subtitle = str_wrap(effect_noise_results$interpretation, 40)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    xlim(0.5, 1.5) + ylim(0.5, 1.5)

  # Combine plots
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Effect-to-Noise Analysis Summary",
      tag_levels = 'A'
    )

  plots <- list(w1_values = p1, effect_vs_noise = p2,
                precision = p3, signal_strength = p4, combined = combined)

  # Save plots
  if (save_plots) {
    output_path <- file.path(output_dir, "effect_noise_analysis.pdf")
    ggsave(output_path, combined, width = 16, height = 12, device = "pdf")

    if (verbose) {
      cat("✓ Effect-noise analysis plots saved to:", output_path, "\n")
    }
  }

  return(plots)
}


# ==============================================================================
# CONVENIENCE FUNCTION: Run complete bootstrap analysis
# ==============================================================================

#' Run Complete Bootstrap Analysis Pipeline
#'
#' Convenience wrapper that runs the full bootstrap analysis:
#' 1. Load data
#' 2. Perform bootstrap resampling of replicate CDFs
#' 3. Calculate effect-to-noise ratios
#' 4. Generate diagnostic plots
#' 5. Save all results
#'
#' @param data_file Path to standardized data (default: auto-detect)
#' @param n_bootstrap Number of bootstrap iterations (default: 2000)
#' @param reference_module Reference condition (default: "RODOS")
#' @param test_module Test condition (default: "INHALER")
#' @param seed Random seed (default: 42)
#' @param verbose Print progress (default: TRUE)
#'
#' @return List with all analysis results
#'
run_bootstrap_analysis <- function(data_file = NULL, n_bootstrap = 2000,
                                   reference_module = "RODOS",
                                   test_module = "INHALER",
                                   seed = 42, verbose = TRUE) {

  # Load data
  if (is.null(data_file)) {
    data_file <- "data/tidy/standardized_data.csv"
  }

  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file,
         "\nRun 01_data_import.R first to create this file.")
  }

  if (verbose) {
    cat("Loading data from:", data_file, "\n")
  }

  data <- read_csv(data_file, show_col_types = FALSE)

  # Bootstrap analysis
  bootstrap_results <- bootstrap_w1_analysis(
    data = data,
    reference_module = reference_module,
    test_module = test_module,
    n_bootstrap = n_bootstrap,
    seed = seed,
    verbose = verbose
  )

  # Effect-to-noise analysis
  effect_noise_results <- calculate_effect_noise_ratios(
    bootstrap_results = bootstrap_results,
    verbose = verbose
  )

  # Diagnostic plots
  bootstrap_plots <- plot_bootstrap_results(
    bootstrap_results = bootstrap_results,
    verbose = verbose
  )

  # Effect-to-noise plots
  effect_noise_plots <- plot_effect_noise_analysis(
    bootstrap_results = bootstrap_results,
    effect_noise_results = effect_noise_results,
    verbose = verbose
  )

  if (verbose) {
    cat("\n========================================================================\n")
    cat("BOOTSTRAP ANALYSIS COMPLETE\n")
    cat("========================================================================\n")
    cat("Files created:\n")
    cat("  - results/bootstrap_results.csv\n")
    cat("  - results/effect_noise_ratios.csv\n")
    cat("  - figures/bootstrap_analysis.pdf\n")
    cat("  - figures/effect_noise_analysis.pdf\n")
    cat("------------------------------------------------------------------------\n")
    cat("Bootstrap summary:\n")
    cat(sprintf("  %d formulations analyzed\n", nrow(bootstrap_results)))
    cat(sprintf("  %d bootstrap samples per formulation\n", n_bootstrap))
    cat(sprintf("  Mean measurement uncertainty: %.4f µm\n",
                mean(bootstrap_results$w1_sd, na.rm = TRUE)))
    cat("------------------------------------------------------------------------\n")
    cat("Effect-to-noise summary:\n")
    if (nrow(effect_noise_results) > 0) {
      max_ratio <- max(effect_noise_results$effect_to_noise_ratio, na.rm = TRUE)
      cat(sprintf("  Maximum effect-to-noise ratio: %.2f\n", max_ratio))
      strong_effects <- sum(effect_noise_results$effect_to_noise_ratio >= 3, na.rm = TRUE)
      cat(sprintf("  Strong effects (ratio ≥ 3): %d/%d\n",
                  strong_effects, nrow(effect_noise_results)))
    }
    cat("========================================================================\n\n")
  }

  return(list(
    bootstrap_results = bootstrap_results,
    effect_noise_ratios = effect_noise_results,
    bootstrap_plots = bootstrap_plots,
    effect_noise_plots = effect_noise_plots,
    data = data
  ))
}


# ==============================================================================
# AUTO-EXECUTION: Run analysis when script is sourced
# ==============================================================================

# Check if processed data exists
if (file.exists("data/tidy/standardized_data.csv")) {

  cat("\n========================================================================\n")
  cat("AUTO-RUNNING BOOTSTRAP ANALYSIS\n")
  cat("========================================================================\n")
  cat("Reading: data/tidy/standardized_data.csv\n")
  cat("Bootstrap iterations: 2000 per formulation\n")
  cat("Saving to: results/bootstrap_results.csv\n")
  cat("------------------------------------------------------------------------\n")

  # Run the complete analysis
  .bootstrap_analysis <- run_bootstrap_analysis(verbose = TRUE)

  cat("\n========================================================================\n")
  cat("BOOTSTRAP ANALYSIS COMPLETE\n")
  cat("========================================================================\n")
  cat("Next steps:\n")
  cat("  - Review confidence intervals in results/bootstrap_results.csv\n")
  cat("  - Check diagnostic plots in figures/bootstrap_analysis.pdf\n")
  cat("  - Compare effect-to-noise ratios in results/effect_noise_ratios.csv\n")
  cat("  - Examine device condition effects in figures/effect_noise_analysis.pdf\n")
  cat("------------------------------------------------------------------------\n")
  cat("To reload results later:\n")
  cat("  source('scripts/04_bootstrap_analysis.R')\n")
  cat("  bootstrap_results <- read_csv('results/bootstrap_results.csv')\n")
  cat("  effect_noise_ratios <- read_csv('results/effect_noise_ratios.csv')\n")
  cat("========================================================================\n\n")

} else {
  cat("\n========================================================================\n")
  cat("BOOTSTRAP ANALYSIS - WAITING FOR INPUT DATA\n")
  cat("========================================================================\n")
  cat("Standardized data not found: data/tidy/standardized_data.csv\n")
  cat("\nPlease run the data processing pipeline first:\n")
  cat("  source('scripts/01_data_import.R')\n")
  cat("  source('scripts/02_wasserstein_core.R')\n")
  cat("  source('scripts/04_bootstrap_analysis.R')\n")
  cat("========================================================================\n\n")
}
