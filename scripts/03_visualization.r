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
    group_by(formulation, particle_size_um) %>%
    summarise(
      q3_percent_mean = mean(q3_percent, na.rm = TRUE),
      q3_percent_sd = sd(q3_percent, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(q3_percent_sd = ifelse(is.na(q3_percent_sd), 0, q3_percent_sd))

  n_formulations <- n_distinct(summary_data$formulation)

  p <- ggplot(summary_data, aes(x = particle_size_um, y = q3_percent_mean,
                                 color = formulation, fill = formulation)) +
    geom_ribbon(
      aes(ymin = q3_percent_mean - q3_percent_sd,
          ymax = q3_percent_mean + q3_percent_sd),
      alpha = 0.1, color = NA
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
    scale_color_viridis_d(name = "Formulation", option = "turbo") +
    scale_fill_viridis_d(guide = "none", option = "turbo") +
    labs(
      x = "Particle Size (µm)",
      y = expression("Cumulative Distribution " * Q[3] * " (%)"),
      title = paste0("All ", test_module, " Distributions Overlay"),
      subtitle = paste("Comparing", n_formulations, "formulations")
    ) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )

  output_path <- file.path(output_dir, filename)
  ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat("✓ Saved:", output_path, "\n")
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
    bar_color = "#1F78B4",
    show_values = TRUE,
    save_plot = FALSE,
    output_dir = "figures_v2",
    filename = "w1_ranking.pdf",  # <-- CHANGED default to PDF
    width = 10,
    height = 6,
    dpi = 300
) {

  if (!metric %in% c("W1_micrometers", "W1_normalized")) {
    stop("metric must be 'W1_micrometers' or 'W1_normalized'")
  }

  plot_data <- w1_results %>%
    select(formulation, all_of(metric))

  if (sort_by) {
    plot_data <- plot_data %>%
      arrange(!!sym(metric)) %>%
      mutate(formulation = factor(formulation, levels = formulation))
  }

  y_label <- if (metric == "W1_micrometers") {
    "Wasserstein Distance W₁ (µm)"
  } else {
    "Normalized Wasserstein Distance (W₁/d₅₀)"
  }

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

  if (show_values) {
    p <- p + geom_text(
      aes(label = sprintf("%.3f", !!sym(metric))),
      vjust = -0.5,
      size = 3.5,
      fontface = "bold"
    )
  }

  # Save if requested
  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    output_path <- file.path(output_dir, filename)

    # Force base PDF device for reliability
    ggsave(
      filename = output_path,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = "pdf"
    )

    # Assert it actually wrote
    if (!file.exists(output_path)) {
      stop("ggsave completed but file not found at: ", normalizePath(output_path, winslash = "/"))
    }

    cat("✓ Saved:", normalizePath(output_path, winslash = "/"), "\n")
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

  if (verbose) {
    cat("------------------------------------------------------------------------\n")
    cat("PLOTS COMPLETE - Figures saved to", output_dir, "\n")
    cat("========================================================================\n\n")
  }

  return(invisible(list(
    individual = p_individual,
    overlay = p_overlay,
    w1 = p_w1
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
  cat("  - figures/*_comparison.pdf (", n_formulations, " individual formulation PDFs)\n", sep = "")
  cat("  - figures/all_inhaler_overlay.pdf (all INHALER distributions)\n")
  cat("  - figures/w1_ranking.pdf (dispersibility ranking)\n")  # <-- CHANGED to PDF
  cat("------------------------------------------------------------------------\n")
  cat("Additional plots available via functions:\n")
  cat("  - plot_d50_comparison() for median size comparison\n")
  cat("  - create_publication_panel() for combined figures\n")
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
