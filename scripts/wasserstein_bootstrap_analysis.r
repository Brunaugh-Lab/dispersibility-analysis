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
                               test_module = "INHALER",
                               device_resistance = NULL,
                               pressure_drop = NULL,
                               n_bootstrap = 2000,
                               seed = NULL, verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  if (verbose) {
    condition_label <- paste(
      formulation,
      if (!is.null(device_resistance)) paste0("Device:", device_resistance) else "",
      if (!is.null(pressure_drop)) paste0("Pressure:", pressure_drop) else ""
    )
    cat(sprintf("Bootstrap resampling: %s (%d iterations)\n", condition_label, n_bootstrap))
  }

  # Load calculate_wasserstein_1d function from script 02
  if (!exists("calculate_wasserstein_1d")) {
    source("scripts/02_wasserstein_core.R", local = TRUE)
  }

  # Get individual replicate CDFs for both conditions
  # RODOS has no device/pressure factors - pool all replicates
  ref_replicates <- get_replicate_cdfs(data, formulation, reference_module)

  # INHALER filtered by device and pressure
  test_replicates <- get_replicate_cdfs(data, formulation, test_module,
                                        device_resistance = device_resistance,
                                        pressure_drop = pressure_drop)

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
  ref_replicates <- get_replicate_cdfs(data, formulation, reference_module,
                                       common_size_grid = common_sizes)
  test_replicates <- get_replicate_cdfs(data, formulation, test_module,
                                        device_resistance = device_resistance,
                                        pressure_drop = pressure_drop,
                                        common_size_grid = common_sizes)

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

  # Return bootstrap results
  return(tibble(
    formulation = formulation,
    device_resistance = device_resistance %||% NA_character_,
    pressure_drop = pressure_drop %||% NA_character_,
    w1_mean = mean(w1_bootstrap, na.rm = TRUE),
    w1_sd = sd(w1_bootstrap, na.rm = TRUE),
    w1_ci_lower = quantile(w1_bootstrap, 0.025, na.rm = TRUE),
    w1_ci_upper = quantile(w1_bootstrap, 0.975, na.rm = TRUE),
    w1_observed = w1_observed,
    n_bootstrap = n_bootstrap,
    bootstrap_samples = list(w1_bootstrap)
  ))
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

  # Get all unique combinations of formulation × device × pressure
  test_conditions <- data %>%
    filter(module == test_module) %>%
    distinct(formulation, device_resistance, pressure_drop_clean)

  n_formulations <- n_distinct(test_conditions$formulation)
  n_conditions <- nrow(test_conditions)

  if (verbose) {
    cat("Formulations to process:", n_formulations, "\n")
    cat("Total conditions (formulation × device × pressure):", n_conditions, "\n")
    cat("Total bootstrap samples:", n_conditions * n_bootstrap, "\n\n")
  }

  # Run bootstrap analysis for each condition combination
  bootstrap_results <- pmap_dfr(test_conditions, function(formulation, device_resistance, pressure_drop_clean) {
    bootstrap_w1_single(
      data = data,
      formulation = formulation,
      reference_module = reference_module,
      test_module = test_module,
      device_resistance = device_resistance,
      pressure_drop = pressure_drop_clean,
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
    cat("Analyzing: Formulation, Device, and Pressure effects\n")
  }

    # ============================================================================
  # 1. FORMULATION EFFECT: Variability across formulations
  # ============================================================================
  formulation_effect <- bootstrap_results %>%
    group_by(device_resistance, pressure_drop) %>%
    summarise(
      effect_size = sd(w1_mean, na.rm = TRUE),
      avg_noise = mean(w1_sd, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    summarise(
      effect_magnitude_um = mean(effect_size, na.rm = TRUE),
      noise_level_um = mean(avg_noise, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      factor_type = "formulation",
      effect_to_noise_ratio = effect_magnitude_um / noise_level_um,
      n_levels = n_distinct(bootstrap_results$formulation),
      interpretation = case_when(
        effect_to_noise_ratio >= 3 ~ "Strong formulation effect: Differences >> measurement noise",
        effect_to_noise_ratio >= 2 ~ "Moderate formulation effect: Differences > measurement noise",
        effect_to_noise_ratio >= 1 ~ "Weak formulation effect: Differences ~ measurement noise",
        TRUE ~ "Poor formulation signal: Differences < measurement noise"
      )
    )

  # ============================================================================
  # 2. DEVICE RESISTANCE EFFECT: Variability across device levels
  # ============================================================================
  # Only calculate if device_resistance column exists AND has >1 unique level
  if ("device_resistance" %in% names(bootstrap_results) &&
      n_distinct(bootstrap_results$device_resistance, na.rm = TRUE) > 1) {

    device_effect <- bootstrap_results %>%
      group_by(formulation, pressure_drop) %>%
      summarise(
        effect_size = max(w1_mean, na.rm = TRUE) - min(w1_mean, na.rm = TRUE),
        avg_noise = mean(w1_sd, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      summarise(
        effect_magnitude_um = mean(effect_size, na.rm = TRUE),
        noise_level_um = mean(avg_noise, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        factor_type = "device_resistance",
        effect_to_noise_ratio = effect_magnitude_um / noise_level_um,
        n_levels = n_distinct(bootstrap_results$device_resistance, na.rm = TRUE),
        interpretation = case_when(
          effect_to_noise_ratio >= 3 ~ "Strong device effect: Differences >> measurement noise",
          effect_to_noise_ratio >= 2 ~ "Moderate device effect: Differences > measurement noise",
          effect_to_noise_ratio >= 1 ~ "Weak device effect: Differences ~ measurement noise",
          TRUE ~ "Poor device signal: Differences < measurement noise"
        )
      )
  } else {
    device_effect <- NULL
    if (verbose) {
      cat("  Note: Device resistance effect not calculated (only 1 level or column missing)\n")
    }
  }

  # ============================================================================
  # 3. PRESSURE DROP EFFECT: Variability across pressure levels
  # ============================================================================
  # Only calculate if pressure_drop column exists AND has >1 unique level
  if ("pressure_drop" %in% names(bootstrap_results) &&
      n_distinct(bootstrap_results$pressure_drop, na.rm = TRUE) > 1) {

    pressure_effect <- bootstrap_results %>%
      group_by(formulation, device_resistance) %>%
      summarise(
        effect_size = max(w1_mean, na.rm = TRUE) - min(w1_mean, na.rm = TRUE),
        avg_noise = mean(w1_sd, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      summarise(
        effect_magnitude_um = mean(effect_size, na.rm = TRUE),
        noise_level_um = mean(avg_noise, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        factor_type = "pressure_drop",
        effect_to_noise_ratio = effect_magnitude_um / noise_level_um,
        n_levels = n_distinct(bootstrap_results$pressure_drop, na.rm = TRUE),
        interpretation = case_when(
          effect_to_noise_ratio >= 3 ~ "Strong pressure effect: Differences >> measurement noise",
          effect_to_noise_ratio >= 2 ~ "Moderate pressure effect: Differences > measurement noise",
          effect_to_noise_ratio >= 1 ~ "Weak pressure effect: Differences ~ measurement noise",
          TRUE ~ "Poor pressure signal: Differences < measurement noise"
        )
      )
  } else {
    pressure_effect <- NULL
    if (verbose) {
      cat("  Note: Pressure drop effect not calculated (only 1 level or column missing)\n")
    }
  }

  # Combine all effect-to-noise ratios (only include non-NULL)
  ratio_results <- bind_rows(
    formulation_effect,
    device_effect,
    pressure_effect
  )

  if (verbose) {
    cat("------------------------------------------------------------------------\n")
    cat("Effect-to-Noise Ratios:\n")
    for (i in 1:nrow(ratio_results)) {
      cat(sprintf("  %s: %.2f (%s)\n",
                  str_to_title(ratio_results$factor_type[i]),
                  ratio_results$effect_to_noise_ratio[i],
                  ratio_results$interpretation[i]))
    }
    cat("------------------------------------------------------------------------\n")
  }

  # Save results
  if (save_output) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    output_path <- file.path(output_dir, output_filename)
    write_csv(ratio_results, output_path)

    if (verbose) {
      cat("✓ Effect-to-noise ratios saved to:", output_path, "\n")
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

  # Detect number of factor levels
  n_devices <- n_distinct(bootstrap_results$device_resistance, na.rm = TRUE)
  n_pressures <- n_distinct(bootstrap_results$pressure_drop, na.rm = TRUE)
  has_device <- "device_resistance" %in% names(bootstrap_results) && n_devices > 1
  has_pressure <- "pressure_drop" %in% names(bootstrap_results) && n_pressures > 1

  # Create condition label for x-axis (combines all non-varying factors)
  bootstrap_results <- bootstrap_results %>%
    mutate(
      condition_label = case_when(
        has_device & has_pressure ~ paste0(device_resistance, "\n", pressure_drop),
        has_device ~ as.character(device_resistance),
        has_pressure ~ as.character(pressure_drop),
        TRUE ~ "All"
      )
    )

  # Plot 1: Confidence intervals with smart faceting
  p1 <- bootstrap_results %>%
    ggplot(aes(x = formulation, color = formulation)) +
    geom_errorbar(aes(ymin = w1_ci_lower, ymax = w1_ci_upper),
                  width = 0.3, alpha = 0.7) +
    geom_point(aes(y = w1_observed), size = 3) +
    geom_point(aes(y = w1_mean), size = 2, alpha = 0.5, shape = 1) +
    labs(
      x = "Formulation",
      y = "Wasserstein Distance (µm)",
      title = "Bootstrap Confidence Intervals for W1",
      subtitle = "Filled: Observed W1, Open: Bootstrap Mean, Bars: 95% CI"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )

  # Add faceting if we have varying factors
  if (has_device & has_pressure) {
    p1 <- p1 + facet_grid(device_resistance ~ pressure_drop,
                          labeller = labeller(
                            device_resistance = ~paste("Resistance:", .),
                            pressure_drop = ~paste("Pressure:", .)
                          )) +
      theme(strip.text = element_text(size = 9))
  } else if (has_device) {
    p1 <- p1 + facet_wrap(~device_resistance,
                          labeller = labeller(device_resistance = ~paste("Resistance:", .))) +
      theme(strip.text = element_text(size = 9))
  } else if (has_pressure) {
    p1 <- p1 + facet_wrap(~pressure_drop,
                          labeller = labeller(pressure_drop = ~paste("Pressure:", .))) +
      theme(strip.text = element_text(size = 9))
  }

  # Plot 2: Standard errors - show by condition
  p2 <- bootstrap_results %>%
    ggplot(aes(x = formulation, y = w1_sd, fill = formulation)) +
    geom_col(alpha = 0.7) +
    labs(
      x = "Formulation",
      y = "Bootstrap Standard Error (µm)",
      title = "Measurement Uncertainty by Formulation"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )

  if (has_device & has_pressure) {
    p2 <- p2 + facet_grid(device_resistance ~ pressure_drop,
                          labeller = labeller(
                            device_resistance = ~paste("Resistance:", .),
                            pressure_drop = ~paste("Pressure:", .)
                          )) +
      theme(strip.text = element_text(size = 9))
  } else if (has_device) {
    p2 <- p2 + facet_wrap(~device_resistance,
                          labeller = labeller(device_resistance = ~paste("Resistance:", .))) +
      theme(strip.text = element_text(size = 9))
  } else if (has_pressure) {
    p2 <- p2 + facet_wrap(~pressure_drop,
                          labeller = labeller(pressure_drop = ~paste("Pressure:", .))) +
      theme(strip.text = element_text(size = 9))
  }

  # Plot 3: Relative standard error
  p3 <- bootstrap_results %>%
    ggplot(aes(x = formulation, y = 100 * w1_relative_se, fill = formulation)) +
    geom_col(alpha = 0.7) +
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
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )

  if (has_device & has_pressure) {
    p3 <- p3 + facet_grid(device_resistance ~ pressure_drop,
                          labeller = labeller(
                            device_resistance = ~paste("Resistance:", .),
                            pressure_drop = ~paste("Pressure:", .)
                          )) +
      theme(strip.text = element_text(size = 9))
  } else if (has_device) {
    p3 <- p3 + facet_wrap(~device_resistance,
                          labeller = labeller(device_resistance = ~paste("Resistance:", .))) +
      theme(strip.text = element_text(size = 9))
  } else if (has_pressure) {
    p3 <- p3 + facet_wrap(~pressure_drop,
                          labeller = labeller(pressure_drop = ~paste("Pressure:", .))) +
      theme(strip.text = element_text(size = 9))
  }

  # Combine plots with dynamic sizing
  combined <- (p1 / p2 / p3) +
    plot_annotation(
      title = "Bootstrap Analysis Summary",
      tag_levels = 'A'
    )

  # Calculate dynamic plot dimensions
  plot_width <- if (has_device & has_pressure) {
    max(14, 5 + n_pressures * 3)
  } else if (has_device | has_pressure) {
    max(12, 8 + max(n_devices, n_pressures) * 2)
  } else {
    12
  }

  plot_height <- if (has_device & has_pressure) {
    max(15, 12 + n_devices * 2)
  } else {
    12
  }

  plots <- list(confidence_intervals = p1, standard_errors = p2,
                relative_errors = p3, combined = combined)

  # Save plots
  if (save_plots) {
    output_path <- file.path(output_dir, "bootstrap_analysis.pdf")
    ggsave(output_path, combined, width = plot_width, height = plot_height, device = "pdf")

    if (verbose) {
      cat("✓ Bootstrap plots saved to:", output_path, "\n")
    }
  }  # <-- ADD THIS CLOSING BRACE

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

  # Plot 2: Effect vs Noise comparison - improved bar chart
  p2 <- effect_noise_results %>%
    mutate(
      factor_label = case_when(
        factor_type == "formulation" ~ "Formulation",
        factor_type == "device_resistance" ~ "Device Resistance",
        factor_type == "pressure_drop" ~ "Pressure Drop",
        TRUE ~ str_to_title(factor_type)
      ),
      factor_label = fct_reorder(factor_label, effect_to_noise_ratio)
    ) %>%
    pivot_longer(
      cols = c(effect_magnitude_um, noise_level_um),
      names_to = "metric_type",
      values_to = "value"
    ) %>%
    mutate(
      metric_label = if_else(metric_type == "effect_magnitude_um", "Effect Size", "Noise Level")
    ) %>%
    ggplot(aes(x = factor_label, y = value, fill = metric_label)) +
    geom_col(alpha = 0.8, position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", value)),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3) +
    scale_fill_manual(
      values = c("Effect Size" = "darkgreen", "Noise Level" = "orange"),
      name = ""
    ) +
    labs(
      x = "Factor",
      y = "Magnitude (µm)",
      title = "Effect Size vs Measurement Noise"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "top"
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

  # Plot 4: Effect-to-Noise Ratio bar chart (clearer than circles)
  p4 <- effect_noise_results %>%
    mutate(
      factor_label = case_when(
        factor_type == "formulation" ~ "Formulation",
        factor_type == "device_resistance" ~ "Device Resistance",
        factor_type == "pressure_drop" ~ "Pressure Drop",
        TRUE ~ str_to_title(factor_type)
      ),
      factor_label = fct_reorder(factor_label, effect_to_noise_ratio),
      interpretation_color = case_when(
        effect_to_noise_ratio >= 3 ~ "darkgreen",
        effect_to_noise_ratio >= 2 ~ "orange",
        effect_to_noise_ratio >= 1 ~ "gold",
        TRUE ~ "red"
      ),
      signal_category = case_when(
        effect_to_noise_ratio >= 3 ~ "Strong (≥3)",
        effect_to_noise_ratio >= 2 ~ "Moderate (2-3)",
        effect_to_noise_ratio >= 1 ~ "Weak (1-2)",
        TRUE ~ "Poor (<1)"
      )
    ) %>%
    ggplot(aes(x = factor_label, y = effect_to_noise_ratio)) +
    geom_col(aes(fill = interpretation_color), alpha = 0.8, width = 0.6) +
    geom_hline(yintercept = c(1, 2, 3), linetype = "dashed", alpha = 0.4) +
    geom_text(aes(label = sprintf("%.2f", effect_to_noise_ratio)),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_identity() +
    labs(
      x = "Factor",
      y = "Effect-to-Noise Ratio",
      title = "Signal Strength by Factor",
      subtitle = "Higher ratios indicate factor effects exceed measurement noise"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

  # Combine plots - improved layout
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
    data_file <- "data/tidy/standardized_data_with_conditions.csv"
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
if (file.exists("data/tidy/standardized_data_with_conditions.csv")) {

  cat("\n========================================================================\n")
  cat("AUTO-RUNNING BOOTSTRAP ANALYSIS\n")
  cat("========================================================================\n")
  cat("Reading: data/tidy/standardized_data_with_conditions.csv\n")
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
  cat("Standardized data not found: data/tidy/standardized_data_with_conditions.csv\n")
  cat("\nPlease run the data processing pipeline first:\n")
  cat("  source('scripts/01_data_import.R')\n")
  cat("  source('scripts/02_wasserstein_core.R')\n")
  cat("  source('scripts/04_bootstrap_analysis.R')\n")
  cat("========================================================================\n\n")
}
