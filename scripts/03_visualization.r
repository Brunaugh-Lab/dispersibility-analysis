# ==============================================================================
# 03_visualization.R
# Visualization Functions for Dispersibility Analysis
#
# PURPOSE
#   Publication-ready plotting utilities for particle size distributions
#   and dispersibility metrics. All CDF-based plots follow the same
#   replicate-pooling logic used in Wasserstein-1 calculations (script 02).
#
#   This script defines visualization functions only and does NOT
#   auto-execute when sourced.
#
# ------------------------------------------------------------------
# HOW TO USE (MANUAL EXECUTION)
#
#   source("scripts/03_visualization.R")
#
#   data <- readr::read_csv(
#     file.path(data_dir, "tidy", "standardized_data_with_conditions.csv"),
#     show_col_types = FALSE
#   )
#
#   w1_results <- readr::read_csv(
#     file.path(results_dir, "wasserstein_results.csv"),
#     show_col_types = FALSE
#   )
#
#   generate_all_plots(
#     data = data,
#     w1_results = w1_results,
#     output_dir = figures_dir
#   )
#
# ------------------------------------------------------------------
# INPUTS (from upstream pipeline)
#   - data/tidy/standardized_data_with_conditions.csv   (script 01)
#   - results/wasserstein_results.csv                   (script 02)
#
# OUTPUTS
#   - figures/*.pdf
#
#   Examples:
#     * One PDF per formulation showing pooled reference (RODOS) vs
#       all available INHALER test conditions
#     * Per-formulation reference vs single-condition comparison PDFs
#     * W₁ ranking and factor-faceted bar plots
#
# ------------------------------------------------------------------
# AVAILABLE PLOTTING FUNCTIONS
#
#   Core CDF plots:
#     - plot_reference_vs_test()
#     - export_pairwise_condition_pdfs()
#     - export_formulation_overlay_reference_plus_all_tests()
#
#   Summary / metric plots:
#     - plot_w1_bars()
#
#   Factor-faceted plots:
#     - plot_cdf_by_device()
#     - plot_cdf_by_pressure()
#     - plot_w1_by_device()
#     - plot_w1_by_pressure()
#
# ------------------------------------------------------------------
# DESIGN PRINCIPLES
#   - Consistent with W1 methodology: replicates are pooled before plotting
#   - Explicit reference vs test comparisons (no implicit overlays)
#   - Modular functions with no side effects on source()
#   - Publication-oriented defaults with defensive input checks
#
# ==============================================================================

.require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      "\nInstall with:\n  install.packages(c(",
      paste0('"', missing, '"', collapse = ", "),
      "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.require_pkgs(c(
  "ggplot2", "dplyr", "readr", "stringr", "tidyr", "purrr", "forcats", "tibble",
  "viridis", "patchwork", "grid", "rlang"
))

# ==============================================================================
# PATH CONFIG (canonical pipeline directories)
# ==============================================================================
data_dir    <- "data"
results_dir <- "results"
figures_dir <- "figures"

# Canonical input file paths
tidy_data_path  <- file.path(data_dir, "tidy", "standardized_data_with_conditions.csv")
w1_results_path <- file.path(results_dir, "wasserstein_results.csv")

# ==============================================================================
# CORE PLOTTER: Reference vs Single Test Condition (returns ggplot object)
# ==============================================================================

plot_reference_vs_test <- function(
  ref_summary,
  test_summary,
  reference_label = "RODOS",
  test_label = "Test",
  color_palette = c("reference" = "#E31A1C", "test" = "#1F78B4")
) {

  required <- c("particle_size_um", "q3_percent_mean", "q3_percent_sd")
  missing_ref  <- setdiff(required, names(ref_summary))
  missing_test <- setdiff(required, names(test_summary))

  if (length(missing_ref) > 0) {
    stop("ref_summary is missing: ", paste(missing_ref, collapse = ", "), call. = FALSE)
  }
  if (length(missing_test) > 0) {
    stop("test_summary is missing: ", paste(missing_test, collapse = ", "), call. = FALSE)
  }

  summary_data <- dplyr::bind_rows(
    dplyr::mutate(ref_summary,  curve = "reference"),
    dplyr::mutate(test_summary, curve = "test")
  ) |>
    dplyr::mutate(
      curve = factor(curve, levels = c("reference", "test")),
      curve = dplyr::recode(curve, reference = reference_label, test = test_label)
    )

  # IMPORTANT: names must match the *actual* curve labels after recode()
  curve_cols <- stats::setNames(
    unname(color_palette[c("reference", "test")]),
    c(reference_label, test_label)
  )

  ggplot2::ggplot(
    summary_data,
    ggplot2::aes(
      x = .data$particle_size_um,
      y = .data$q3_percent_mean,
      color = .data$curve,
      fill  = .data$curve
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$q3_percent_mean - .data$q3_percent_sd,
        ymax = .data$q3_percent_mean + .data$q3_percent_sd
      ),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::scale_x_log10(
      limits = c(0.5, 100),
      breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::scale_color_manual(values = curve_cols, name = NULL) +
    ggplot2::scale_fill_manual(values = curve_cols, guide = "none") +
    ggplot2::labs(
      x = "Particle Size (µm)",
      y = expression("Cumulative Distribution " * Q[3] * " (%)")
    ) +
    ggplot2::theme_classic(base_size = 12)
}

# ==============================================================================
# EXPORTER: Save Reference vs Test PDFs for All Formulations × Conditions
#   (Pooling matches 02_wasserstein_core.R: average within replicate first,
#    then average across replicates)
# ==============================================================================

export_pairwise_condition_pdfs <- function(
  data,
  reference_module = "RODOS",
  test_module = "INHALER",
  condition_cols = c("device_resistance", "pressure_drop_clean"),
  output_dir = figures_dir,
  verbose = TRUE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Defensive check: replicate column must exist for replicate-first pooling
  if (!"replicate" %in% names(data)) {
    stop("Expected column 'replicate' not found in data.", call. = FALSE)
  }

  formulations <- unique(data$formulation)

  for (form in formulations) {

    form_data <- dplyr::filter(data, .data$formulation == form)

    # ---- reference (pooled once per formulation; replicate-first) ----
    ref_summary <- form_data |>
      dplyr::filter(.data$module == reference_module) |>
      # stage 1: one value per (replicate, size)
      dplyr::group_by(.data$particle_size_um, .data$replicate) |>
      dplyr::summarise(
        q3_percent_rep = mean(.data$q3_percent, na.rm = TRUE),
        .groups = "drop"
      ) |>
      # stage 2: pool across replicates at each size
      dplyr::group_by(.data$particle_size_um) |>
      dplyr::summarise(
        q3_percent_mean = mean(.data$q3_percent_rep, na.rm = TRUE),
        q3_percent_sd   = stats::sd(.data$q3_percent_rep, na.rm = TRUE),
        n_replicates    = dplyr::n(),
        .groups = "drop"
      )

    if (nrow(ref_summary) == 0) {
      if (isTRUE(verbose)) cat("Skipping ", form, ": no reference rows\n", sep = "")
      next
    }

    # ---- discover test conditions ----
    test_conditions <- form_data |>
      dplyr::filter(.data$module == test_module) |>
      dplyr::distinct(dplyr::across(dplyr::all_of(condition_cols)))

    if (nrow(test_conditions) == 0) {
      if (isTRUE(verbose)) cat("Skipping ", form, ": no test rows\n", sep = "")
      next
    }

    for (i in seq_len(nrow(test_conditions))) {

      condition <- test_conditions[i, , drop = FALSE]
      condition_vals <- unlist(condition, use.names = FALSE)

      # subset to this test condition (your explicit loop is fine + readable)
      test_subset <- form_data |>
        dplyr::filter(.data$module == test_module)

      for (j in seq_along(condition_cols)) {
        col <- condition_cols[j]
        val <- condition[[col]][[1]]
        test_subset <- dplyr::filter(test_subset, .data[[col]] == val)
      }

      # ---- test summary (replicate-first) ----
      test_summary <- test_subset |>
        # stage 1: one value per (replicate, size)
        dplyr::group_by(.data$particle_size_um, .data$replicate) |>
        dplyr::summarise(
          q3_percent_rep = mean(.data$q3_percent, na.rm = TRUE),
          .groups = "drop"
        ) |>
        # stage 2: pool across replicates at each size
        dplyr::group_by(.data$particle_size_um) |>
        dplyr::summarise(
          q3_percent_mean = mean(.data$q3_percent_rep, na.rm = TRUE),
          q3_percent_sd   = stats::sd(.data$q3_percent_rep, na.rm = TRUE),
          n_replicates    = dplyr::n(),
          .groups = "drop"
        )

      if (nrow(test_summary) == 0) next

      # ---- build plot ----
      p <- plot_reference_vs_test(
        ref_summary,
        test_summary,
        reference_label = reference_module,
        test_label = paste(condition_vals, collapse = ", ")
      )

      # ---- filename ----
      condition_slug <- paste(
        paste(condition_cols, condition_vals, sep = "_"),
        collapse = "__"
      )

      filename <- paste0(
        as.character(form), "__", condition_slug, "__",
        reference_module, "_vs_", test_module, ".pdf"
      )

      ggplot2::ggsave(
        file.path(output_dir, filename),
        p,
        width = 8,
        height = 6,
        device = "pdf"
      )

      if (isTRUE(verbose)) {
        cat("✓ Saved:", filename, "\n")
      }
    }
  }
}

# ==============================================================================
# EXPORTER 2: Per-formulation overlay
#   One PDF per formulation showing:
#     - Reference (e.g., RODOS) pooled once per formulation
#     - All test conditions (e.g., INHALER device×pressure) overlaid
#   Pooling matches 02_wasserstein_core.R: average within replicate first,
#   then average across replicates.
# ==============================================================================

export_formulation_overlay_reference_plus_all_tests <- function(
  data,
  reference_module = "RODOS",
  test_module = "INHALER",
  condition_cols = c("device_resistance", "pressure_drop_clean"),
  output_dir = figures_dir,
  filename_suffix = "reference_plus_all_tests.pdf",
  width = 10,
  height = 7,
  verbose = TRUE
) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (!"replicate" %in% names(data)) {
    stop("Expected column 'replicate' not found in data.", call. = FALSE)
  }

  # Defensive check: required columns exist
  needed_cols <- c("formulation", "module", "particle_size_um", "q3_percent", "replicate", condition_cols)
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Helper: replicate-first pooling for q3_percent
  pool_q3_percent <- function(df) {
    df |>
      dplyr::group_by(.data$particle_size_um, .data$replicate) |>
      dplyr::summarise(
        q3_percent_rep = mean(.data$q3_percent, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::group_by(.data$particle_size_um) |>
      dplyr::summarise(
        q3_percent_mean = mean(.data$q3_percent_rep, na.rm = TRUE),
        q3_percent_sd   = stats::sd(.data$q3_percent_rep, na.rm = TRUE),
        n_replicates    = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$particle_size_um)
  }

  formulations <- unique(data$formulation)

  for (form in formulations) {

    form_data <- dplyr::filter(data, .data$formulation == form)

    # ---- reference pooled once per formulation ----
    ref_df <- dplyr::filter(form_data, .data$module == reference_module)

    if (nrow(ref_df) == 0) {
      if (isTRUE(verbose)) cat("Skipping ", form, ": no reference rows\n", sep = "")
      next
    }

    ref_summary <- pool_q3_percent(ref_df) |>
      dplyr::mutate(curve = reference_module)

    # ---- test conditions present for this formulation ----
    test_df <- dplyr::filter(form_data, .data$module == test_module)

    if (nrow(test_df) == 0) {
      if (isTRUE(verbose)) cat("Skipping ", form, ": no test rows\n", sep = "")
      next
    }

    test_conditions <- test_df |>
      dplyr::distinct(dplyr::across(dplyr::all_of(condition_cols)))

    if (nrow(test_conditions) == 0) {
      if (isTRUE(verbose)) cat("Skipping ", form, ": no test conditions\n", sep = "")
      next
    }

    # ---- build pooled test curves for each condition ----
    test_summaries <- vector("list", nrow(test_conditions))

    for (i in seq_len(nrow(test_conditions))) {

      cond <- test_conditions[i, , drop = FALSE]

      cond_subset <- test_df
      for (j in seq_along(condition_cols)) {
        col <- condition_cols[j]
        val <- cond[[col]][[1]]
        cond_subset <- dplyr::filter(cond_subset, .data[[col]] == val)
      }

      if (nrow(cond_subset) == 0) next

      # Robust, unambiguous label (includes column names)
      cond_label <- paste(
        paste0(
          condition_cols, "=",
          vapply(condition_cols, \(cc) as.character(cond[[cc]][[1]]), character(1))
        ),
        collapse = " | "
      )

      test_summaries[[i]] <- pool_q3_percent(cond_subset) |>
        dplyr::mutate(curve = cond_label)
    }

    # Drop NULL entries defensively
    test_summaries <- purrr::compact(test_summaries)
    test_summary_all <- dplyr::bind_rows(test_summaries)

    if (nrow(test_summary_all) == 0) next

    # ---- combine + deterministic factor levels ----
    plot_data <- dplyr::bind_rows(ref_summary, test_summary_all)

    test_levels  <- sort(unique(test_summary_all$curve))
    curve_levels <- c(reference_module, test_levels)

    plot_data <- plot_data |>
      dplyr::mutate(curve = factor(.data$curve, levels = curve_levels))

    # ---- deterministic colors: reference fixed + tests viridis (generated ONCE) ----
    test_cols <- viridis::viridis(length(test_levels), option = "turbo")
    names(test_cols) <- test_levels

    curve_cols <- c(setNames("#E31A1C", reference_module), test_cols)

    # ---- plot ----
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$particle_size_um,
        y = .data$q3_percent_mean,
        color = .data$curve,
        fill  = .data$curve
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$q3_percent_mean - dplyr::coalesce(.data$q3_percent_sd, 0),
          ymax = .data$q3_percent_mean + dplyr::coalesce(.data$q3_percent_sd, 0)
        ),
        alpha = 0.12,
        color = NA
      ) +
      ggplot2::geom_line(linewidth = 1.0) +
      ggplot2::scale_x_log10(
        limits = c(0.5, 100),
        breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
        labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      ggplot2::scale_color_manual(values = curve_cols, name = NULL) +
      ggplot2::scale_fill_manual(values = curve_cols, guide = "none") +
      ggplot2::labs(
        x = "Particle Size (µm)",
        y = expression("Cumulative Distribution " * Q[3] * " (%)"),
        title = paste0(form, ": ", reference_module, " reference vs all ", test_module, " conditions"),
        subtitle = paste0(
          "Overlay of ", nrow(test_conditions), " test conditions (",
          paste(condition_cols, collapse = " × "), ")"
        )
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = "right",
        legend.text = ggplot2::element_text(size = 9),
        plot.title = ggplot2::element_text(face = "bold", size = 14)
      )

    out_file <- paste0(as.character(form), "__", filename_suffix)
    out_path <- file.path(output_dir, out_file)
    ggplot2::ggsave(out_path, p, width = width, height = height, device = "pdf")

    if (isTRUE(verbose)) cat("✓ Saved:", out_file, "\n")
  }

  invisible(TRUE)
}

# ==============================================================================
# FUNCTION 3: Plot Wasserstein Distance Bar Chart
# ==============================================================================

plot_w1_bars <- function(
    w1_results,
    metric = "W1_micrometers",
    sort_by = TRUE,
    save_plot = FALSE,
    output_dir = figures_dir,
    filename = "w1_ranking.pdf",
    width = NULL,
    height = NULL,
    dpi = 300
) {

  metric_labels <- list(
    W1_micrometers = "W₁ Distance (µm)",
    W1_normalized  = "W₁/d₅₀ (Normalized)",
    d50_shift_um   = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]
  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um", call. = FALSE)
  }

  device_levels <- w1_results |>
    dplyr::distinct(device_resistance) |>
    dplyr::arrange(device_resistance) |>
    dplyr::pull(device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- w1_results |>
    dplyr::distinct(pressure_drop) |>
    dplyr::mutate(numeric_pressure = as.numeric(stringr::str_extract(pressure_drop, "\\d+"))) |>
    dplyr::arrange(numeric_pressure) |>
    dplyr::pull(pressure_drop)

  plot_data <- w1_results |>
    dplyr::mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop     = factor(pressure_drop, levels = pressure_levels)
    )

  if (isTRUE(sort_by)) {
    order_levels <- w1_results |>
      dplyr::group_by(formulation) |>
      dplyr::summarise(mean_metric = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(mean_metric)) |>
      dplyr::pull(formulation)

    plot_data <- plot_data |>
      dplyr::mutate(formulation = factor(formulation, levels = order_levels))
  }

  device_labels <- stats::setNames(
    stringr::str_to_title(stringr::str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- stats::setNames(
    stringr::str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  n_devices      <- length(device_levels)
  n_pressures    <- length(pressure_levels)
  n_formulations <- dplyr::n_distinct(plot_data$formulation)

  if (is.null(width)) {
    width <- max(10, 5 + n_pressures * 2 + n_formulations * 0.5)
  }
  if (is.null(height)) {
    height <- max(8, 3 + n_devices * 2.5)
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = formulation, y = .data[[metric]])) +
    ggplot2::geom_col(fill = "#1F78B4", color = "black", linewidth = 0.3) +
    ggplot2::facet_grid(
      device_resistance ~ pressure_drop,
      labeller = ggplot2::labeller(
        device_resistance = device_labels,
        pressure_drop     = pressure_labels
      )
    ) +
    ggplot2::labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility Ranking by Device Resistance × Pressure Drop",
      subtitle = sprintf(
        "Lower W₁ = Better dispersibility (%d×%d conditions)",
        n_devices, n_pressures
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 9),
      axis.title  = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.title    = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 11),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text       = ggplot2::element_text(face = "bold", size = 10)
    )

  if (isTRUE(save_plot)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    output_path <- file.path(output_dir, filename)
    ggplot2::ggsave(output_path, p, width = width, height = height, dpi = dpi)
    cat(sprintf("✓ Saved: %s (%d×%d grid)\n", output_path, n_devices, n_pressures))
  }

  return(p)
}




# ==============================================================================
# FUNCTION 6: Plot W1 Bars Faceted by Device Resistance
# ==============================================================================

plot_w1_by_device <- function(
    w1_results,
    metric = "W1_micrometers",
    output_dir = figures_dir,
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
    W1_normalized  = "W₁/d₅₀ (Normalized)",
    d50_shift_um   = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]
  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um", call. = FALSE)
  }

  device_levels <- w1_results |>
    dplyr::distinct(.data$device_resistance) |>
    dplyr::arrange(.data$device_resistance) |>
    dplyr::pull(.data$device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- w1_results |>
    dplyr::distinct(.data$pressure_drop) |>
    dplyr::mutate(
      numeric_pressure = as.numeric(stringr::str_extract(.data$pressure_drop, "\\d+"))
    ) |>
    dplyr::arrange(.data$numeric_pressure) |>
    dplyr::pull(.data$pressure_drop)

  plot_data <- w1_results |>
    dplyr::mutate(
      device_resistance = factor(.data$device_resistance, levels = device_levels),
      pressure_drop     = factor(.data$pressure_drop,     levels = pressure_levels)
    )

  device_labels <- setNames(
    stringr::str_to_title(stringr::str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- setNames(
    stringr::str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  n_devices      <- length(device_levels)
  n_formulations <- dplyr::n_distinct(plot_data$formulation)

  if (is.null(width))  width  <- max(12, 3 * n_devices + n_formulations * 0.3)
  if (is.null(height)) height <- 6

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$formulation, y = .data[[metric]], fill = .data$pressure_drop)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.9),
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::facet_wrap(
      ~ device_resistance,
      nrow = 1,
      labeller = ggplot2::labeller(device_resistance = device_labels)
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "plasma",
      name = "Device Pressure Drop",
      breaks = pressure_levels,
      labels = unname(pressure_labels[pressure_levels])
    ) +
    ggplot2::labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility by Device Resistance"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.title         = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.title         = ggplot2::element_text(face = "bold", size = 16),
      strip.background   = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text         = ggplot2::element_text(face = "bold", size = 12),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_text(face = "bold")
    )

  output_path <- file.path(output_dir, filename)
  ggplot2::ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d devices)\n", output_path, n_devices))
  }

  return(p)
}

# ==============================================================================
# FUNCTION 7: Plot W1 Bars Faceted by Pressure Drop
# ==============================================================================

plot_w1_by_pressure <- function(
    w1_results,
    metric = "W1_micrometers",
    output_dir = figures_dir,
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
    W1_normalized  = "W₁/d₅₀ (Normalized)",
    d50_shift_um   = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]
  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um", call. = FALSE)
  }

  device_levels <- w1_results |>
    dplyr::distinct(.data$device_resistance) |>
    dplyr::arrange(.data$device_resistance) |>
    dplyr::pull(.data$device_resistance)

  if (all(c("low", "medium", "high") %in% device_levels)) {
    device_levels <- c("low", "medium", "high")
  }

  pressure_levels <- w1_results |>
    dplyr::distinct(.data$pressure_drop) |>
    dplyr::mutate(
      numeric_pressure = as.numeric(stringr::str_extract(.data$pressure_drop, "\\d+"))
    ) |>
    dplyr::arrange(.data$numeric_pressure) |>
    dplyr::pull(.data$pressure_drop)

  plot_data <- w1_results |>
    dplyr::mutate(
      device_resistance = factor(.data$device_resistance, levels = device_levels),
      pressure_drop     = factor(.data$pressure_drop,     levels = pressure_levels)
    )

  device_labels <- setNames(
    stringr::str_to_title(stringr::str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- setNames(
    stringr::str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  n_pressures    <- length(pressure_levels)
  n_formulations <- dplyr::n_distinct(plot_data$formulation)

  if (is.null(width))  width  <- max(12, 3 * n_pressures + n_formulations * 0.3)
  if (is.null(height)) height <- 6

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$formulation, y = .data[[metric]], fill = .data$device_resistance)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.9),
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::facet_wrap(
      ~ pressure_drop,
      nrow = 1,
      labeller = ggplot2::labeller(pressure_drop = pressure_labels)
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "plasma",
      name = "Device Resistance",
      breaks = device_levels,
      labels = unname(device_labels[device_levels])
    ) +
    ggplot2::labs(
      x = "Formulation",
      y = y_label,
      title = "Dispersibility by Pressure Drop"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.title         = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.title         = ggplot2::element_text(face = "bold", size = 16),
      strip.background   = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text         = ggplot2::element_text(face = "bold", size = 12),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_text(face = "bold")
    )

  output_path <- file.path(output_dir, filename)
  ggplot2::ggsave(output_path, p, width = width, height = height, device = "pdf")

  if (verbose) {
    cat(sprintf("✓ Saved: %s (%d pressures)\n", output_path, n_pressures))
  }

  return(p)
}

# ==============================================================================
# CONVENIENCE FUNCTION: Generate all standard plots
# ==============================================================================
generate_all_plots <- function(
    data,
    w1_results,
    output_dir = figures_dir,
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

  # ---- Exporter 1: pairwise PDFs ----
  if (verbose) cat("Creating pairwise reference vs test PDFs...\n")
  export_pairwise_condition_pdfs(
    data,
    reference_module = reference_module,
    test_module = test_module,
    output_dir = output_dir,
    verbose = verbose
  )
  p_individual <- NULL

  # ---- Exporter 2: per-formulation overlay ----
  if (verbose) cat("Creating per-formulation reference + test overlay PDFs...\n")
  export_formulation_overlay_reference_plus_all_tests(
    data,
    reference_module = reference_module,
    test_module = test_module,
    output_dir = output_dir,
    verbose = verbose
  )
  p_overlay <- NULL

  # ---- W1 ranking figure ----
  if (verbose) cat("Creating W1 ranking plot...\n")
  p_w1 <- plot_w1_bars(
    w1_results,
    metric = "W1_micrometers",
    save_plot = TRUE,
    output_dir = output_dir,
    filename = "w1_ranking.pdf",
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
