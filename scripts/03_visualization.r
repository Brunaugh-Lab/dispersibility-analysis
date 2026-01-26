# ==============================================================================
# 03_visualization.R
# Visualization Functions for Dispersibility Analysis
#
# Purpose: Publication-ready plots for comparing particle size distributions
#          and dispersibility metrics across formulations
#
# Designed for: Single test condition vs reference (e.g., INHALER vs RODOS)
#               Focus on formulation-level comparisons
# ==============================================================================

library(tidyverse)
library(viridis)
library(patchwork)  # For combining plots

# ==============================================================================
# FUNCTION 1: Plot CDF Comparison
# ==============================================================================

#' Plot Cumulative Distribution Function Comparison
#'
#' Overlays test and reference CDFs for visual comparison of dispersibility.
#' Shows mean ± SD ribbons across technical replicates.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param reference_module Character string for reference (default: "RODOS")
#' @param test_module Character string for test (default: "INHALER")
#' @param formulations Optional character vector to subset formulations
#'   If NULL (default), plots all formulations
#' @param color_palette Named vector of colors for modules
#'   Default: c("RODOS" = "#E31A1C", "INHALER" = "#1F78B4")
#' @param facet_by Facet plots by formulation? (default: TRUE)
#' @param ncol Number of columns for faceting (default: 3)
#'
#' @return ggplot object
#'
#' @examples
#' # Basic usage
#' p <- plot_cdf_comparison(data)
#' print(p)
#' ggsave("cdf_comparison.png", p, width = 12, height = 8, dpi = 300)
#'
#' # Subset specific formulations
#' p <- plot_cdf_comparison(data, formulations = c("FormA", "FormB", "FormC"))
#'
plot_cdf_comparison <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    formulations = NULL,
    color_palette = c("RODOS" = "#E31A1C", "INHALER" = "#1F78B4"),
    facet_by = TRUE,
    ncol = 3
) {

  # Subset formulations if specified
  if (!is.null(formulations)) {
    data <- data %>% filter(formulation %in% formulations)
  }

  # Filter to only reference and test modules
  plot_data <- data %>%
    filter(module %in% c(reference_module, test_module))

  # Calculate mean and SD for each formulation-module-size combination
  cdf_summary <- plot_data %>%
    group_by(formulation, module, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  # Ensure module is a factor with correct ordering (reference first)
  cdf_summary$module <- factor(
    cdf_summary$module,
    levels = c(reference_module, test_module)
  )

  # Create plot
  p <- ggplot(cdf_summary, aes(x = particle_size_um, y = q3_percent_mean,
                                color = module, fill = module)) +
    geom_ribbon(
      aes(ymin = q3_percent_mean - q3_percent_sd,
          ymax = q3_percent_mean + q3_percent_sd),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.0) +
    scale_x_log10(
      limits = c(0.5, 100),
      breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
      labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20)
    ) +
    scale_color_manual(
      values = color_palette,
      name = "Module"
    ) +
    scale_fill_manual(
      values = color_palette,
      guide = "none"
    ) +
    labs(
      x = "Particle Size (µm)",
      y = "Cumulative Distribution Q₃ (%)",
      title = "Cumulative Particle Size Distributions",
      subtitle = paste0(test_module, " vs ", reference_module, " Comparison")
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )

  # Add faceting if requested
  if (facet_by) {
    p <- p + facet_wrap(~ formulation, ncol = ncol)
  }

  return(p)
}


# ==============================================================================
# FUNCTION 2: Plot Particle Size Density (PSD)
# ==============================================================================

#' Plot Particle Size Density Distribution
#'
#' Shows the derivative of the CDF (density) on log scale, which better
#' reveals multimodal distributions and fine particle populations.
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param reference_module Character string for reference (default: "RODOS")
#' @param test_module Character string for test (default: "INHALER")
#' @param formulations Optional character vector to subset formulations
#' @param color_palette Named vector of colors for modules
#' @param facet_by Facet plots by formulation? (default: TRUE)
#' @param ncol Number of columns for faceting (default: 3)
#'
#' @return ggplot object
#'
plot_psd_density <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    formulations = NULL,
    color_palette = c("RODOS" = "#E31A1C", "INHALER" = "#1F78B4"),
    facet_by = TRUE,
    ncol = 3
) {

  # Subset formulations if specified
  if (!is.null(formulations)) {
    data <- data %>% filter(formulation %in% formulations)
  }

  # Filter to only reference and test modules
  plot_data <- data %>%
    filter(module %in% c(reference_module, test_module))

  # Calculate density (dQ3/d(log x)) by numerical differentiation
  psd_data <- plot_data %>%
    arrange(formulation, module, replicate, particle_size_um) %>%
    group_by(formulation, module, replicate) %>%
    mutate(
      log_size = log10(particle_size_um),
      dQ3_dlogx = c(0, diff(q3_percent) / diff(log_size))  # Numerical derivative
    ) %>%
    ungroup()

  # Calculate mean and SD
  psd_summary <- psd_data %>%
    group_by(formulation, module, particle_size_um) %>%
    summarise(
      dQ3_dlogx_mean = mean(dQ3_dlogx, na.rm = TRUE),
      dQ3_dlogx_sd = sd(dQ3_dlogx, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(dQ3_dlogx_sd = ifelse(is.na(dQ3_dlogx_sd), 0, dQ3_dlogx_sd))

  # Ensure module is a factor
  psd_summary$module <- factor(
    psd_summary$module,
    levels = c(reference_module, test_module)
  )

  # Create plot
  p <- ggplot(psd_summary, aes(x = particle_size_um, y = dQ3_dlogx_mean,
                                color = module, fill = module)) +
    geom_ribbon(
      aes(ymin = dQ3_dlogx_mean - dQ3_dlogx_sd,
          ymax = dQ3_dlogx_mean + dQ3_dlogx_sd),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 0.8) +
    scale_x_log10(
      limits = c(0.5, 100),
      breaks = c(1, 10, 100),
      minor_breaks = c(0.5, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90),
      labels = c("1", "10", "100")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = color_palette, name = "Module") +
    scale_fill_manual(values = color_palette, guide = "none") +
    labs(
      x = "Particle Size (µm)",
      y = "Particle Size Density (dQ₃/d log x)",
      title = "Particle Size Density Distributions"
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (facet_by) {
    p <- p + facet_wrap(~ formulation, ncol = ncol, scales = "free_y")
  }

  return(p)
}


# ==============================================================================
# FUNCTION 3: Plot Wasserstein Distance Bar Chart
# ==============================================================================

#' Plot Wasserstein Distance Bar Chart
#'
#' Visualizes W1 distances across formulations, ranked by dispersibility.
#' Lower W1 = better dispersibility.
#'
#' @param w1_results Results tibble from calculate_pairwise_wasserstein()
#' @param metric Character string: "W1_micrometers" or "W1_normalized"
#' @param sort_by Should bars be sorted by W1 value? (default: TRUE)
#' @param bar_color Color for bars (default: "#1F78B4")
#' @param show_values Show W1 values on bars? (default: TRUE)
#'
#' @return ggplot object
#'
#' @examples
#' # Absolute W1, sorted by dispersibility
#' p1 <- plot_w1_bars(w1_results, metric = "W1_micrometers")
#'
#' # Normalized W1
#' p2 <- plot_w1_bars(w1_results, metric = "W1_normalized")
#'
plot_w1_bars <- function(
    w1_results,
    metric = "W1_micrometers",
    sort_by = TRUE,
    bar_color = "#1F78B4",
    show_values = TRUE
) {

  # Validate metric
  if (!metric %in% c("W1_micrometers", "W1_normalized")) {
    stop("metric must be 'W1_micrometers' or 'W1_normalized'")
  }

  # Prepare data
  plot_data <- w1_results %>%
    select(formulation, all_of(metric))

  # Sort if requested
  if (sort_by) {
    plot_data <- plot_data %>%
      arrange(!!sym(metric)) %>%
      mutate(formulation = factor(formulation, levels = formulation))
  }

  # Determine y-axis label
  y_label <- if (metric == "W1_micrometers") {
    "Wasserstein Distance W₁ (µm)"
  } else {
    "Normalized Wasserstein Distance (W₁/d₅₀)"
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = formulation, y = !!sym(metric))) +
    geom_col(fill = bar_color, color = "black", linewidth = 0.3) +
    labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility Ranking by Wasserstein Distance",
      subtitle = "Lower W₁ = Better Dispersibility"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )

  # Add value labels if requested
  if (show_values) {
    p <- p + geom_text(
      aes(label = sprintf("%.3f", !!sym(metric))),
      vjust = -0.5,
      size = 3.5,
      fontface = "bold"
    )
  }

  return(p)
}


# ==============================================================================
# FUNCTION 4: Plot d50 Comparison
# ==============================================================================

#' Plot d50 Comparison Bar Chart
#'
#' Compares median particle diameters between reference and test conditions.
#' Shows the d50 shift (test - reference) visually.
#'
#' @param w1_results Results tibble from calculate_pairwise_wasserstein()
#' @param sort_by Sort by d50_shift? (default: TRUE)
#' @param reference_color Color for reference bars (default: "#E31A1C")
#' @param test_color Color for test bars (default: "#1F78B4")
#'
#' @return ggplot object
#'
plot_d50_comparison <- function(
    w1_results,
    sort_by = TRUE,
    reference_color = "#E31A1C",
    test_color = "#1F78B4"
) {

  # Reshape data for plotting
  plot_data <- w1_results %>%
    select(formulation, d50_reference_um, d50_test_um) %>%
    pivot_longer(
      cols = c(d50_reference_um, d50_test_um),
      names_to = "condition",
      values_to = "d50"
    ) %>%
    mutate(
      condition = recode(condition,
                        d50_reference_um = "Reference",
                        d50_test_um = "Test")
    )

  # Sort by d50 shift if requested
  if (sort_by) {
    order_levels <- w1_results %>%
      arrange(d50_shift_um) %>%
      pull(formulation)

    plot_data <- plot_data %>%
      mutate(formulation = factor(formulation, levels = order_levels))
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = formulation, y = d50, fill = condition)) +
    geom_col(position = position_dodge(width = 0.8),
             color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Reference" = reference_color, "Test" = test_color),
      name = "Condition"
    ) +
    labs(
      x = "Formulation",
      y = "Median Diameter d₅₀ (µm)",
      title = "Median Particle Size Comparison",
      subtitle = "Reference vs Test Conditions"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )

  return(p)
}


# ==============================================================================
# FUNCTION 5: Create Publication Figure Panel
# ==============================================================================

#' Create Multi-Panel Publication Figure
#'
#' Combines multiple plots into a publication-ready figure using patchwork.
#' Default layout: CDF comparison + W1 bars + d50 comparison
#'
#' @param data Standardized data from read_ld_data_from_structure()
#' @param w1_results Results tibble from calculate_pairwise_wasserstein()
#' @param layout Character string specifying layout:
#'   "horizontal" (default), "vertical", or "grid"
#' @param reference_module Reference module name (default: "RODOS")
#' @param test_module Test module name (default: "INHALER")
#'
#' @return patchwork object (combined ggplot)
#'
#' @examples
#' # Create full panel
#' fig <- create_publication_panel(data, w1_results)
#' ggsave("Figure_Dispersibility.png", fig, width = 16, height = 12, dpi = 300)
#'
create_publication_panel <- function(
    data,
    w1_results,
    layout = "horizontal",
    reference_module = "RODOS",
    test_module = "INHALER"
) {

  # Create individual plots
  p1 <- plot_cdf_comparison(data, reference_module, test_module, facet_by = TRUE)
  p2 <- plot_w1_bars(w1_results, metric = "W1_micrometers")
  p3 <- plot_d50_comparison(w1_results)

  # Combine based on layout
  if (layout == "horizontal") {
    combined <- p1 / (p2 | p3)
  } else if (layout == "vertical") {
    combined <- p1 / p2 / p3
  } else if (layout == "grid") {
    combined <- (p1 | p2) / p3
  } else {
    stop("layout must be 'horizontal', 'vertical', or 'grid'")
  }

  # Add panel labels
  combined <- combined +
    plot_annotation(
      tag_levels = 'A',
      theme = theme(plot.tag = element_text(face = "bold", size = 16))
    )

  return(combined)
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# Uncomment to test with your data:
#
# # Load previous scripts
# source("scripts/01_data_import.R")
# source("scripts/02_wasserstein_core.R")
#
# # Import data and calculate W1
# data <- read_ld_data_from_structure("data/") %>%
#   mutate(replicate = tolower(replicate))
#
# w1_results <- calculate_pairwise_wasserstein(data)
#
# # Load visualization script
# source("scripts/03_visualization.R")
#
# # Create individual plots
# p_cdf <- plot_cdf_comparison(data)
# p_w1 <- plot_w1_bars(w1_results)
# p_d50 <- plot_d50_comparison(w1_results)
# p_psd <- plot_psd_density(data)
#
# # Save plots
# ggsave("CDF_comparison.png", p_cdf, width = 12, height = 8, dpi = 300)
# ggsave("W1_ranking.png", p_w1, width = 10, height = 6, dpi = 300)
# ggsave("d50_comparison.png", p_d50, width = 10, height = 6, dpi = 300)
# ggsave("PSD_density.png", p_psd, width = 12, height = 8, dpi = 300)
#
# # Create combined publication figure
# fig <- create_publication_panel(data, w1_results, layout = "horizontal")
# ggsave("Figure_Dispersibility_Panel.png", fig, width = 16, height = 12, dpi = 300)
