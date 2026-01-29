# ==============================================================================
# 03_visualization.R
# Visualization Functions for Dispersibility Analysis
#
# Purpose: Publication-ready plots for comparing particle size distributions
#          and dispersibility metrics across formulations. Automatically loads
#          data and generates figures.
#
# Auto-execution: Script automatically runs when sourced
#   - Reads data_v2/tidy/standardized_data.csv (from script 01)
#   - Reads results_v2/wasserstein_results.csv (from script 02)
#   - Generates figures and saves to figures/
#     * One PDF per formulation (RODOS vs INHALER comparison)
#     * One overlay PDF (all INHALER distributions)
#     * One PDF (W1 ranking bars)
#
# Additional plots available as functions (call manually if needed):
#   - plot_d50_comparison() - Median particle size comparison
#   - create_publication_panel() - Combined multi-panel figure
#   - plot_psd_density() - Density distributions
#
# - **NEW**: 4 factor-faceted plotting functions for simplified comparisons
# - Creates publication-ready PDFs with dynamic sizing
# - Auto-detects experimental design (devices × pressures × formulations)

# **New Plot Types:**
# - `CDF_by_device.pdf` - All formulations compared across device resistance levels
# - `CDF_by_pressure.pdf` - All formulations compared across pressure drop levels
# - `W1_by_device.pdf` - Dispersibility metrics faceted by device resistance
# - `W1_by_pressure.pdf` - Dispersibility metrics faceted by pressure drop

# **Features:**
# - Flexible design: Works with any n×m experimental design
# - Smart sorting: Natural ordering of factors (low→med→high, numeric pressures)
# - Dynamic sizing: Plot dimensions scale automatically
# - Publication formatting: Consistent themes, proper axis labels
# ```

# #### **2b. Update figures_v2/ Directory**
# The script will create these new files when run:
# ```
# figures_v2/
# ├── CDF_by_device.pdf        ← NEW
# ├── CDF_by_pressure.pdf      ← NEW
# ├── W1_by_device.pdf         ← NEW
# └── W1_by_pressure.pdf       ← NEW
#
# Input:
#   - data_v2/tidy/standardized_data.csv
#   - results_v2/wasserstein_results.csv
# Output: figures/*.pdf and figures/*.png
#
# Designed for: Single test condition vs reference (e.g., INHALER vs RODOS)
#               Focus on formulation-level comparisons
#
# Usage:
#   source("scripts/03_visualization.R")  # That's it!
#
# ==============================================================================

library(tidyverse)
library(viridis)
library(patchwork)  # For combining plots

# ==============================================================================
# FUNCTION 1: Plot Individual Formulation Comparison (PDF per formulation)
# ==============================================================================

plot_individual_formulation_pdfs <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    formulations = NULL,
    color_palette = c("RODOS" = "#E31A1C", "INHALER" = "#1F78B4"),
    output_dir = "figures_v2",
    width = 8,
    height = 6,
    verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (!is.null(formulations)) {
    data <- data %>% filter(formulation %in% formulations)
  }

  plot_data <- data %>%
    filter(module %in% c(reference_module, test_module))

  formulation_list <- unique(plot_data$formulation)

  plots <- list()

  for (form in formulation_list) {

    form_data <- plot_data %>%
      filter(formulation == form)

    # Pool RODOS reference (no device/pressure filtering)
    rodos_summary <- form_data %>%
      filter(module == reference_module) %>%
      group_by(particle_size_um) %>%
      summarise(
        q3_percent_mean = mean(q3_percent, na.rm = TRUE),
        q3_percent_sd = sd(q3_percent, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd),
        module = reference_module
      )

    # Pool INHALER by device×pressure combination
    inhaler_summary <- form_data %>%
      filter(module == test_module) %>%
      group_by(device_resistance, pressure_drop_clean, particle_size_um) %>%
      summarise(
        q3_percent_mean = mean(q3_percent, na.rm = TRUE),
        q3_percent_sd = sd(q3_percent, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd),
        module = test_module
      )

    # Get unique device-pressure combinations
    device_pressure_combos <- inhaler_summary %>%
      distinct(device_resistance, pressure_drop_clean)

    # Expand RODOS to match all device-pressure combos (for faceting)
    rodos_expanded <- device_pressure_combos %>%
      cross_join(rodos_summary)

    # Combine for plotting
    summary_data <- bind_rows(rodos_expanded, inhaler_summary) %>%
      mutate(module = factor(module, levels = c(reference_module, test_module)))

    # FLEXIBLE: Auto-detect factor levels and create natural ordering
    # For device_resistance: natural order (low, medium, high) or alphabetical
    device_levels <- summary_data %>%
      distinct(device_resistance) %>%
      arrange(device_resistance) %>%
      pull(device_resistance)

    # Try to sort as: low < medium < high if those levels exist
    if (all(c("low", "medium", "high") %in% device_levels)) {
      device_levels <- c("low", "medium", "high")
    }

    # For pressure_drop: extract numeric values and sort
    pressure_levels <- summary_data %>%
      distinct(pressure_drop_clean) %>%
      mutate(
        numeric_pressure = as.numeric(str_extract(pressure_drop_clean, "\\d+"))
      ) %>%
      arrange(numeric_pressure) %>%
      pull(pressure_drop_clean)

    summary_data <- summary_data %>%
      mutate(
        device_resistance = factor(device_resistance, levels = device_levels),
        pressure_drop_clean = factor(pressure_drop_clean, levels = pressure_levels)
      )

    # FLEXIBLE: Create human-readable labels
    device_labels <- setNames(
      str_to_title(str_replace_all(device_levels, "_", " ")),
      device_levels
    )

    pressure_labels <- setNames(
      str_replace(pressure_levels, "_", " "),
      pressure_levels
    )

    # FLEXIBLE: Calculate optimal plot dimensions based on number of facets
    n_devices <- length(device_levels)
    n_pressures <- length(pressure_levels)
    plot_width <- max(10, 4 + n_pressures * 3)  # Min 10", scales with columns
    plot_height <- max(8, 3 + n_devices * 2.5)  # Min 8", scales with rows

    p <- ggplot(summary_data, aes(x = particle_size_um, y = q3_percent_mean,
                                   color = module, fill = module)) +
      geom_ribbon(
        aes(ymin = q3_percent_mean - q3_percent_sd,
            ymax = q3_percent_mean + q3_percent_sd),
        alpha = 0.15, color = NA
      ) +
      geom_line(linewidth = 1.0) +
      facet_grid(device_resistance ~ pressure_drop_clean,
                 labeller = labeller(
                   device_resistance = device_labels,
                   pressure_drop_clean = pressure_labels
                 )) +
      scale_x_log10(
        limits = c(0.5, 100),
        breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
        labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
      ) +
      scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, 20)
      ) +
      scale_color_manual(values = color_palette, name = "Module") +
      scale_fill_manual(values = color_palette, guide = "none") +
      labs(
        x = "Particle Size (µm)",
        y = expression("Cumulative Distribution " * Q[3] * " (%)"),
        title = paste("Formulation:", form),
        subtitle = sprintf("INHALER (%d devices × %d pressures) vs %s Reference",
                          n_devices, n_pressures, reference_module)
      ) +
      theme_classic(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 10)
      )

    filename <- paste0(form, "_comparison_faceted.pdf")
    output_path <- file.path(output_dir, filename)
    ggsave(output_path, p, width = plot_width, height = plot_height, device = "pdf")

    if (verbose) {
      cat(sprintf("✓ Saved: %s (%d×%d grid)\n", output_path, n_devices, n_pressures))
    }

    plots[[form]] <- p
  }

  return(invisible(plots))
}

# ==============================================================================
# FUNCTION 2: Plot All INHALER Distributions Overlay
# ==============================================================================

plot_all_inhaler_overlay <- function(
    data,
    test_module = "INHALER",
    formulations = NULL,
    output_dir = "figures_v2",
    filename = "all_inhaler_overlay.pdf",
    width = 10,
    height = 6,
    verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (!is.null(formulations)) {
    data <- data %>% filter(formulation %in% formulations)
  }

  plot_data <- data %>%
    filter(module == test_module)

  summary_data <- plot_data %>%
    group_by(formulation, device_resistance, pressure_drop_clean, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  # FLEXIBLE: Auto-detect factor levels
  device_levels <- summary_data %>%
    distinct(device_resistance) %>%
    arrange(device_resistance) %>%
    pull(device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- summary_data %>%
    distinct(pressure_drop_clean) %>%
    mutate(numeric_pressure = as.numeric(str_extract(pressure_drop_clean, "\\d+"))) %>%
    arrange(numeric_pressure) %>%
    pull(pressure_drop_clean)

  summary_data <- summary_data %>%
    mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop_clean = factor(pressure_drop_clean, levels = pressure_levels)
    )

  # Create labels
  device_labels <- setNames(
    str_to_title(str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- setNames(
    str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  n_formulations <- n_distinct(summary_data$formulation)
  n_devices <- length(device_levels)
  n_pressures <- length(pressure_levels)

  # Dynamic plot sizing
  plot_width <- max(12, 5 + n_pressures * 3)
  plot_height <- max(8, 3 + n_devices * 2.5)

  p <- ggplot(summary_data, aes(x = particle_size_um, y = q3_percent_mean,
                                 color = formulation, fill = formulation)) +
    geom_ribbon(
      aes(ymin = q3_percent_mean - q3_percent_sd,
          ymax = q3_percent_mean + q3_percent_sd),
      alpha = 0.1, color = NA
    ) +
    geom_line(linewidth = 1.0) +
    facet_grid(device_resistance ~ pressure_drop_clean,
               labeller = labeller(
                 device_resistance = device_labels,
                 pressure_drop_clean = pressure_labels
               )) +
    scale_x_log10(
      limits = c(0.5, 100),
      breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
      labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20)
    ) +
    scale_color_viridis_d(name = "Formulation", option = "turbo") +
    scale_fill_viridis_d(guide = "none", option = "turbo") +
    labs(
      x = "Particle Size (µm)",
      y = expression("Cumulative Distribution " * Q[3] * " (%)"),
      title = sprintf("All %s Distributions by Device Resistance × Pressure Drop", test_module),
      subtitle = sprintf("%d formulations across %d conditions (%d×%d grid)",
                        n_formulations, n_devices * n_pressures, n_devices, n_pressures)
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 10)
    )

  output_path <- file.path(output_dir, filename)
  ggsave(output_path, p, width = plot_width, height = plot_height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d×%d grid)\n", output_path, n_devices, n_pressures))
  }

  return(p)
}

# ==============================================================================
# FUNCTION 2: Plot Particle Size Density (PSD)
# ==============================================================================

plot_psd_density <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    formulations = NULL,
    color_palette = c("RODOS" = "#E31A1C", "INHALER" = "#1F78B4"),
    facet_by = TRUE,
    ncol = 3
) {

  if (!is.null(formulations)) {
    data <- data %>% filter(formulation %in% formulations)
  }

  plot_data <- data %>%
    filter(module %in% c(reference_module, test_module))

  psd_data <- plot_data %>%
    arrange(formulation, module, replicate, particle_size_um) %>%
    group_by(formulation, module, replicate) %>%
    mutate(
      log_size = log10(particle_size_um),
      dQ3_dlogx = c(0, diff(q3_percent) / diff(log_size))
    ) %>%
    ungroup()

  psd_summary <- psd_data %>%
    group_by(formulation, module, particle_size_um) %>%
    summarise(
      dQ3_dlogx_mean = mean(dQ3_dlogx, na.rm = TRUE),
      dQ3_dlogx_sd = sd(dQ3_dlogx, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(dQ3_dlogx_sd = ifelse(is.na(dQ3_dlogx_sd), 0, dQ3_dlogx_sd))

  psd_summary$module <- factor(
    psd_summary$module,
    levels = c(reference_module, test_module)
  )

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

plot_w1_bars <- function(
    w1_results,
    metric = "W1_micrometers",
    sort_by = TRUE,
    save_plot = FALSE,
    output_dir = "figures_v2",
    filename = "w1_ranking.pdf",
    width = NULL,  # Now auto-calculated if NULL
    height = NULL, # Now auto-calculated if NULL
    dpi = 300
) {

  metric_labels <- list(
    W1_micrometers = "W₁ Distance (µm)",
    W1_normalized = "W₁/d₅₀ (Normalized)",
    d50_shift_um = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]

  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um")
  }

  # FLEXIBLE: Auto-detect factor levels
  device_levels <- w1_results %>%
    distinct(device_resistance) %>%
    arrange(device_resistance) %>%
    pull(device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- w1_results %>%
    distinct(pressure_drop) %>%
    mutate(numeric_pressure = as.numeric(str_extract(pressure_drop, "\\d+"))) %>%
    arrange(numeric_pressure) %>%
    pull(pressure_drop)

  plot_data <- w1_results %>%
    mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop = factor(pressure_drop, levels = pressure_levels)
    )

  if (sort_by) {
    # Sort formulations by mean metric across all conditions
    order_levels <- w1_results %>%
      group_by(formulation) %>%
      summarise(mean_metric = mean(.data[[metric]], na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(mean_metric)) %>%
      pull(formulation)

    plot_data <- plot_data %>%
      mutate(formulation = factor(formulation, levels = order_levels))
  }

  # Create labels
  device_labels <- setNames(
    str_to_title(str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- setNames(
    str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  n_devices <- length(device_levels)
  n_pressures <- length(pressure_levels)
  n_formulations <- n_distinct(plot_data$formulation)

  # FLEXIBLE: Auto-calculate dimensions if not provided
  if (is.null(width)) {
    width <- max(10, 5 + n_pressures * 2 + n_formulations * 0.5)
  }
  if (is.null(height)) {
    height <- max(8, 3 + n_devices * 2.5)
  }

  p <- ggplot(plot_data, aes(x = formulation, y = .data[[metric]])) +
    geom_col(fill = "#1F78B4", color = "black", linewidth = 0.3) +
    facet_grid(device_resistance ~ pressure_drop,
               labeller = labeller(
                 device_resistance = device_labels,
                 pressure_drop = pressure_labels
               )) +
    labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility Ranking by Device Resistance × Pressure Drop",
      subtitle = sprintf("Lower W₁ = Better dispersibility (%d×%d conditions)",
                        n_devices, n_pressures)
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 10)
    )

  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_path <- file.path(output_dir, filename)
    ggsave(output_path, p, width = width, height = height, dpi = dpi)
    cat(sprintf("✓ Saved: %s (%d×%d grid)\n", output_path, n_devices, n_pressures))
  }

  return(p)
}

# ==============================================================================
# FUNCTION 4: Plot d50 Comparison
# ==============================================================================

plot_d50_comparison <- function(
    w1_results,
    sort_by = TRUE,
    reference_color = "#E31A1C",
    test_color = "#1F78B4",
    save_plot = FALSE,
    output_dir = "figures_v2",
    filename = "d50_comparison.png",
    width = 10,
    height = 6,
    dpi = 300
) {

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

  if (sort_by) {
    order_levels <- w1_results %>%
      arrange(d50_shift_um) %>%
      pull(formulation)

    plot_data <- plot_data %>%
      mutate(formulation = factor(formulation, levels = order_levels))
  }

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

  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_path <- file.path(output_dir, filename)
    ggsave(output_path, p, width = width, height = height, dpi = dpi)
    cat("✓ Saved:", output_path, "\n")
  }

  return(p)
}

# ==============================================================================
# NEW FUNCTION: Plot CDFs Faceted by Device Resistance (all formulations)
# ==============================================================================

plot_cdf_by_device <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    formulation_colors = NULL,
    pressure_linetypes = NULL,
    output_dir = "figures_v2",
    filename = "CDF_by_device.pdf",
    width = NULL,
    height = NULL,
    verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Prepare data
  plot_data <- data %>%
    filter(module %in% c(reference_module, test_module))

  # Pool RODOS (add to all device panels)
  rodos_summary <- plot_data %>%
    filter(module == reference_module) %>%
    group_by(formulation, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  # Pool INHALER by device×pressure
  inhaler_summary <- plot_data %>%
    filter(module == test_module) %>%
    group_by(formulation, device_resistance, pressure_drop_clean, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  # Auto-detect device levels
  device_levels <- inhaler_summary %>%
    distinct(device_resistance) %>%
    arrange(device_resistance) %>%
    pull(device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  # Auto-detect pressure levels and sort numerically
  pressure_levels <- inhaler_summary %>%
    distinct(pressure_drop_clean) %>%
    mutate(numeric_pressure = as.numeric(str_extract(pressure_drop_clean, "\\d+"))) %>%
    arrange(numeric_pressure) %>%
    pull(pressure_drop_clean)

    # Auto-generate linetypes if not provided
  if (is.null(pressure_linetypes)) {
    linetype_options <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    n_pressures_detected <- length(pressure_levels)
    pressure_linetypes <- setNames(
      linetype_options[1:n_pressures_detected],
      pressure_levels
    )
  }

  # Factor levels
  inhaler_summary <- inhaler_summary %>%
    mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop_clean = factor(pressure_drop_clean, levels = pressure_levels)
    )

  # Create device labels
  device_labels <- setNames(
    str_to_title(str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  # Auto-generate formulation colors if not provided
  n_formulations <- n_distinct(inhaler_summary$formulation)
  if (is.null(formulation_colors)) {
    formulation_colors <- viridis_pal(option = "turbo")(n_formulations)
    names(formulation_colors) <- sort(unique(inhaler_summary$formulation))
  }

  # Calculate dimensions
  n_devices <- length(device_levels)
  if (is.null(width)) width <- max(12, 4 * n_devices)
  if (is.null(height)) height <- 6

  p <- ggplot() +
    # INHALER data only (varies by device and pressure)
    geom_line(data = inhaler_summary,
              aes(x = particle_size_um, y = q3_percent_mean,
                  color = formulation, linetype = pressure_drop_clean),
              linewidth = 0.8) +
    facet_wrap(~ device_resistance, nrow = 1,
               labeller = labeller(device_resistance = device_labels)) +
    scale_x_log10(
      limits = c(0.5, 100),
      breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
      labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25)
    ) +
    scale_color_manual(values = formulation_colors, name = "Formulation") +
    scale_linetype_manual(
      values = pressure_linetypes,
      name = "Device Pressure Drop",
      labels = str_replace(names(pressure_linetypes), "_", " "),
      guide = guide_legend(
        override.aes = list(linewidth = 1.2),
        keywidth = unit(3, "cm")
      )
    ) +
    labs(
      x = "Particle Size (µm)",
      y = expression("Cumulative Distribution " * Q[3] * " (%)"),
      title = paste0(test_module, " Dispersibility by Device Resistance")
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 12)
    )

  output_path <- file.path(output_dir, filename)
  ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d devices)\n", output_path, n_devices))
  }

  return(p)
}

# ==============================================================================
# NEW FUNCTION: Plot CDFs Faceted by Pressure Drop (all formulations)
# ==============================================================================

plot_cdf_by_pressure <- function(
    data,
    reference_module = "RODOS",
    test_module = "INHALER",
    formulation_colors = NULL,
    device_linetypes = NULL,
    output_dir = "figures_v2",
    filename = "CDF_by_pressure.pdf",
    width = NULL,
    height = NULL,
    verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Prepare data (same as plot_cdf_by_device)
  plot_data <- data %>%
    filter(module %in% c(reference_module, test_module))

  rodos_summary <- plot_data %>%
    filter(module == reference_module) %>%
    group_by(formulation, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  inhaler_summary <- plot_data %>%
    filter(module == test_module) %>%
    group_by(formulation, device_resistance, pressure_drop_clean, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  # Auto-detect levels
  device_levels <- inhaler_summary %>%
    distinct(device_resistance) %>%
    arrange(device_resistance) %>%
    pull(device_resistance)

# Auto-generate linetypes if not provided
  if (is.null(device_linetypes)) {
    linetype_options <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    n_devices_detected <- length(device_levels)
    device_linetypes <- setNames(
      linetype_options[1:n_devices_detected],
      device_levels
    )
  }

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- inhaler_summary %>%
    distinct(pressure_drop_clean) %>%
    mutate(numeric_pressure = as.numeric(str_extract(pressure_drop_clean, "\\d+"))) %>%
    arrange(numeric_pressure) %>%
    pull(pressure_drop_clean)

  inhaler_summary <- inhaler_summary %>%
    mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop_clean = factor(pressure_drop_clean, levels = pressure_levels)
    )

  # Create pressure labels
  pressure_labels <- setNames(
    str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  # Auto-generate colors
  n_formulations <- n_distinct(inhaler_summary$formulation)
  if (is.null(formulation_colors)) {
    formulation_colors <- viridis_pal(option = "turbo")(n_formulations)
    names(formulation_colors) <- sort(unique(inhaler_summary$formulation))
  }

  # Calculate dimensions
  n_pressures <- length(pressure_levels)
  if (is.null(width)) width <- max(12, 4 * n_pressures)
  if (is.null(height)) height <- 6

  p <- ggplot() +
    # INHALER data only
    geom_line(data = inhaler_summary,
              aes(x = particle_size_um, y = q3_percent_mean,
                  color = formulation, linetype = device_resistance),
              linewidth = 0.8) +
    facet_wrap(~ pressure_drop_clean, nrow = 1,
               labeller = labeller(pressure_drop_clean = pressure_labels)) +
    scale_x_log10(
      limits = c(0.5, 100),
      breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
      labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25)
    ) +
    scale_color_manual(values = formulation_colors, name = "Formulation") +
    scale_linetype_manual(
      values = device_linetypes,
      name = "Device Resistance",
      labels = str_to_title(str_replace_all(names(device_linetypes), "_", " ")),
      guide = guide_legend(
        override.aes = list(linewidth = 1.2),
        keywidth = unit(3, "cm")
      )
    ) +
    labs(
      x = "Particle Size (µm)",
      y = expression("Cumulative Distribution " * Q[3] * " (%)"),
      title = paste0(test_module, " Dispersibility by Pressure Drop")
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 12)
    )

  output_path <- file.path(output_dir, filename)
  ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d pressures)\n", output_path, n_pressures))
  }

  return(p)
}

# ==============================================================================
# NEW FUNCTION: Plot W1 Bars Faceted by Device Resistance
# ==============================================================================

plot_w1_by_device <- function(
    w1_results,
    metric = "W1_micrometers",
    output_dir = "figures_v2",
    filename = "W1_by_device.pdf",
    width = NULL,
    height = NULL,
    verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  metric_labels <- list(
    W1_micrometers = "W₁ Distance (µm)",
    W1_normalized = "W₁/d₅₀ (Normalized)",
    d50_shift_um = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]
  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um")
  }

  # Auto-detect levels
  device_levels <- w1_results %>%
    distinct(device_resistance) %>%
    arrange(device_resistance) %>%
    pull(device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- w1_results %>%
    distinct(pressure_drop) %>%
    mutate(numeric_pressure = as.numeric(str_extract(pressure_drop, "\\d+"))) %>%
    arrange(numeric_pressure) %>%
    pull(pressure_drop)

  plot_data <- w1_results %>%
    mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop = factor(pressure_drop, levels = pressure_levels)
    )

  # Create labels
  device_labels <- setNames(
    str_to_title(str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- setNames(
    str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  # Calculate dimensions
  n_devices <- length(device_levels)
  n_pressures <- length(pressure_levels)
  n_formulations <- n_distinct(plot_data$formulation)

  if (is.null(width)) width <- max(12, 3 * n_devices + n_formulations * 0.3)
  if (is.null(height)) height <- 6

  p <- ggplot(plot_data, aes(x = formulation, y = .data[[metric]], fill = pressure_drop)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", linewidth = 0.3) +
    facet_wrap(~ device_resistance, nrow = 1,
               labeller = labeller(device_resistance = device_labels)) +
    scale_fill_viridis_d(option = "plasma", name = "Device Pressure Drop",
                         labels = pressure_labels) +
    labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility by Device Resistance"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 16),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )

  output_path <- file.path(output_dir, filename)
  ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d devices)\n", output_path, n_devices))
  }

  return(p)
}

# ==============================================================================
# NEW FUNCTION: Plot W1 Bars Faceted by Pressure Drop
# ==============================================================================

plot_w1_by_pressure <- function(
    w1_results,
    metric = "W1_micrometers",
    output_dir = "figures_v2",
    filename = "W1_by_pressure.pdf",
    width = NULL,
    height = NULL,
    verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  metric_labels <- list(
    W1_micrometers = "W₁ Distance (µm)",
    W1_normalized = "W₁/d₅₀ (Normalized)",
    d50_shift_um = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]
  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um")
  }

  # Auto-detect levels
  device_levels <- w1_results %>%
    distinct(device_resistance) %>%
    arrange(device_resistance) %>%
    pull(device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- w1_results %>%
    distinct(pressure_drop) %>%
    mutate(numeric_pressure = as.numeric(str_extract(pressure_drop, "\\d+"))) %>%
    arrange(numeric_pressure) %>%
    pull(pressure_drop)

  plot_data <- w1_results %>%
    mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop = factor(pressure_drop, levels = pressure_levels)
    )

  # Create labels
  device_labels <- setNames(
    str_to_title(str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- setNames(
    str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  # Calculate dimensions
  n_devices <- length(device_levels)
  n_pressures <- length(pressure_levels)
  n_formulations <- n_distinct(plot_data$formulation)

  if (is.null(width)) width <- max(12, 3 * n_pressures + n_formulations * 0.3)
  if (is.null(height)) height <- 6

  p <- ggplot(plot_data, aes(x = formulation, y = .data[[metric]], fill = device_resistance)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", linewidth = 0.3) +
    facet_wrap(~ pressure_drop, nrow = 1,
               labeller = labeller(pressure_drop = pressure_labels)) +
    scale_fill_viridis_d(option = "plasma", name = "Device Resistance",
                         labels = device_labels) +
    labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility by Pressure Drop"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 16),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )

  output_path <- file.path(output_dir, filename)
  ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d pressures)\n", output_path, n_pressures))
  }

  return(p)
}

# ==============================================================================
# FUNCTION 5: Create Publication Figure Panel
# ==============================================================================

create_publication_panel <- function(
    data,
    w1_results,
    layout = "horizontal",
    reference_module = "RODOS",
    test_module = "INHALER",
    save_plot = FALSE,
    output_dir = "figures_v2",
    filename = "dispersibility_panel.png",
    width = 16,
    height = 12,
    dpi = 300
) {

  p1 <- plot_cdf_comparison(data, reference_module, test_module, facet_by = TRUE)
  p2 <- plot_w1_bars(w1_results, metric = "W1_micrometers")
  p3 <- plot_d50_comparison(w1_results)

  if (layout == "horizontal") {
    combined <- p1 / (p2 | p3)
  } else if (layout == "vertical") {
    combined <- p1 / p2 / p3
  } else if (layout == "grid") {
    combined <- (p1 | p2) / p3
  } else {
    stop("layout must be 'horizontal', 'vertical', or 'grid'")
  }

  combined <- combined +
    plot_annotation(
      tag_levels = 'A',
      theme = theme(plot.tag = element_text(face = "bold", size = 16))
    )

  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_path <- file.path(output_dir, filename)
    ggsave(output_path, combined, width = width, height = height, dpi = dpi)
    cat("✓ Saved:", output_path, "\n")
  }

  return(combined)
}



# ==============================================================================
# CONVENIENCE FUNCTION: Generate all standard plots
# ==============================================================================

generate_all_plots <- function(
    data,
    w1_results,
    output_dir = "figures_v2",
    reference_module = "RODOS",
    test_module = "INHALER",
    verbose = TRUE
) {

  if (verbose) {
    cat("\n========================================================================\n")
    cat("GENERATING DISPERSIBILITY PLOTS\n")
    cat("========================================================================\n")
    cat("Output directory:", output_dir, "\n")
    cat("------------------------------------------------------------------------\n")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", output_dir, "\n")
  }

  if (verbose) cat("Creating individual formulation comparison PDFs...\n")
  p_individual <- plot_individual_formulation_pdfs(
    data,
    reference_module = reference_module,
    test_module = test_module,
    output_dir = output_dir,
    verbose = verbose
  )

  if (verbose) cat("Creating INHALER overlay PDF...\n")
  p_overlay <- plot_all_inhaler_overlay(
    data,
    test_module = test_module,
    output_dir = output_dir,
    filename = "all_inhaler_overlay.pdf",
    verbose = verbose
  )

  if (verbose) cat("Creating W1 ranking plot...\n")
  p_w1 <- plot_w1_bars(
    w1_results,
    metric = "W1_micrometers",
    save_plot = TRUE,
    output_dir = output_dir,
    filename = "w1_ranking.pdf",  # <-- CHANGED to PDF
    width = 10,
    height = 6
  )
  # NEW: Factor-faceted plots
  if (verbose) cat("Creating CDF plot faceted by device...\n")
  p_cdf_device <- plot_cdf_by_device(
    data,
    reference_module = reference_module,
    test_module = test_module,
    output_dir = output_dir,
    verbose = verbose
  )

  if (verbose) cat("Creating CDF plot faceted by pressure...\n")
  p_cdf_pressure <- plot_cdf_by_pressure(
    data,
    reference_module = reference_module,
    test_module = test_module,
    output_dir = output_dir,
    verbose = verbose
  )

  if (verbose) cat("Creating W1 plot faceted by device...\n")
  p_w1_device <- plot_w1_by_device(
    w1_results,
    metric = "W1_micrometers",
    output_dir = output_dir,
    verbose = verbose
  )

  if (verbose) cat("Creating W1 plot faceted by pressure...\n")
  p_w1_pressure <- plot_w1_by_pressure(
    w1_results,
    metric = "W1_micrometers",
    output_dir = output_dir,
    verbose = verbose
  )

  if (verbose) {
    cat("------------------------------------------------------------------------\n")
    cat("PLOTS COMPLETE - Figures saved to", output_dir, "\n")
    cat("========================================================================\n\n")
  }

  return(invisible(list(
    individual = p_individual,
    overlay = p_overlay,
    w1 = p_w1,
    cdf_device = p_cdf_device,
    cdf_pressure = p_cdf_pressure,
    w1_device = p_w1_device,
    w1_pressure = p_w1_pressure
  )))
}

# ==============================================================================
# AUTO-EXECUTION: Generate plots when script is sourced
# ==============================================================================

tidy_data_exists <- file.exists("data_v2/tidy/standardized_data_with_conditions.csv")
results_exist <- file.exists("results_v2/wasserstein_results.csv")

if (tidy_data_exists && results_exist) {

  cat("\n========================================================================\n")
  cat("AUTO-RUNNING VISUALIZATION\n")
  cat("========================================================================\n")
  cat("Reading: data_v2/tidy/standardized_data_with_conditions.csv\n")
  cat("Reading: results_v2/wasserstein_results.csv\n")
  cat("Saving to: figures_v2/\n")
  cat("------------------------------------------------------------------------\n")

  .viz_data <- read_csv("data_v2/tidy/standardized_data_with_conditions.csv", show_col_types = FALSE)
  .viz_results <- read_csv("results_v2/wasserstein_results.csv", show_col_types = FALSE)

  .viz_plots <- generate_all_plots(
    data = .viz_data,
    w1_results = .viz_results,
    verbose = TRUE
  )

  cat("\n========================================================================\n")
  cat("VISUALIZATION COMPLETE\n")
  cat("========================================================================\n")
  cat("Generated files:\n")
  n_formulations <- n_distinct(.viz_data$formulation)
  n_devices <- n_distinct(.viz_data$device_resistance)
  n_pressures <- n_distinct(.viz_data$pressure_drop_clean)
  cat(sprintf("  - figures_v2/*_comparison_faceted.pdf (%d formulations, %d×%d grids)\n",
              n_formulations, n_devices, n_pressures))
  cat(sprintf("  - figures_v2/all_inhaler_overlay.pdf (%d×%d facets)\n", n_devices, n_pressures))
  cat(sprintf("  - figures_v2/w1_ranking.pdf (%d×%d facets)\n", n_devices, n_pressures))
  cat("\n  NEW Factor-Faceted Plots:\n")
  cat(sprintf("  - figures_v2/CDF_by_device.pdf (%d device panels)\n", n_devices))
  cat(sprintf("  - figures_v2/CDF_by_pressure.pdf (%d pressure panels)\n", n_pressures))
  cat(sprintf("  - figures_v2/W1_by_device.pdf (%d device panels)\n", n_devices))
  cat(sprintf("  - figures_v2/W1_by_pressure.pdf (%d pressure panels)\n", n_pressures))
  cat("------------------------------------------------------------------------\n")
  cat("Additional plots available via functions:\n")
  cat("  - plot_d50_comparison() for median size comparison\n")
  cat("  - create_publication_panel() for combined figures\n")
  cat("  - plot_cdf_by_device() for device-specific CDF comparison\n")
  cat("  - plot_cdf_by_pressure() for pressure-specific CDF comparison\n")
  cat("  - plot_w1_by_device() for device-specific W1 comparison\n")
  cat("  - plot_w1_by_pressure() for pressure-specific W1 comparison\n")
  cat("========================================================================\n\n")

} else {
  cat("\n========================================================================\n")
  cat("VISUALIZATION - WAITING FOR INPUT DATA\n")
  cat("========================================================================\n")

  if (!tidy_data_exists) {
    cat("✗ Tidy data not found: data_v2/tidy/standardized_data_with_conditions.csv\n")
    cat("  Run: source('scripts/01_data_import.R')\n\n")
  }

  if (!results_exist) {
    cat("✗ Results not found: results_v2/wasserstein_results.csv\n")
    cat("  Run: source('scripts/02_wasserstein_core.R')\n\n")
  }

  cat("Complete pipeline:\n")
  cat("  source('scripts/01_data_import.R')\n")
  cat("  source('scripts/02_wasserstein_core.R')\n")
  cat("  source('scripts/03_visualization.R')\n")
  cat("========================================================================\n\n")
}
