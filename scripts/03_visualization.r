# ==============================================================================
# 03_visualization.R
# Visualization Functions for Dispersibility Analysis
#
# PURPOSE
#   Publication-ready plotting utilities for particle size distributions
#   and dispersibility metrics. All CDF-based plots follow the same
#   replicate-pooling logic used in Wasserstein-1 calculations (script 02).
#
#
# ------------------------------------------------------------------
# HOW TO USE (MANUAL EXECUTION)
#
#   source("scripts/03_visualization.R")
#
#   data <- readr::read_csv(file.path(data_dir, "tidy", "standardized_data_with_conditions.csv"), show_col_types = FALSE)
#
#   w1_results <- readr::read_csv(file.path(results_dir, "wasserstein_results.csv"),show_col_types = FALSE)
#
#   generate_all_plots(data = data, w1_results = w1_results, output_dir = figures_dir)
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
#       all available test conditions (e.g., INHALER device × pressure)
#     * One PDF per formulation × test condition (pairwise reference vs test)
#     * W₁ dispersibility ranking plot
#
# ------------------------------------------------------------------
# AVAILABLE PLOTTING FUNCTIONS
#
#   Core CDF exporters (publication-facing):
#     - plot_reference_vs_test()
#     - export_pairwise_condition_pdfs()
#     - export_formulation_overlay_reference_plus_all_tests()
#
#   Summary / metric plots:
#     - plot_w1_bars()
#
#   (Optional / internal QC utilities may be added separately, but are
#    not part of the default plotting workflow.)
#
# ------------------------------------------------------------------
# DESIGN PRINCIPLES
#   - Methodological consistency with W₁ calculations:
#       replicates are pooled before visualization
#   - Explicit reference vs test comparisons (no ambiguous global overlays)
#   - Exporter-based workflow: figures are written to disk by default
#   - Modular, side-effect-free functions (no execution on source())
#   - Publication-oriented defaults with defensive input validation
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
# FUNCTION: Plot Wasserstein Distance Bars (auto-adapts to dataset shape)
# ==============================================================================

plot_w1_bars <- function(
  w1_results,
  metric = "W1_micrometers",
  sort_by = TRUE,
  facet_mode = c("auto", "grid", "device", "pressure", "none"),
  save_plot = FALSE,
  output_dir = figures_dir,
  filename = "w1_ranking.pdf",
  width = NULL,
  height = NULL,
  dpi = 300
) {

  facet_mode <- match.arg(facet_mode)

  metric_labels <- list(
    W1_micrometers = "W₁ Distance (µm)",
    W1_normalized  = "W₁/d₅₀ (Normalized)",
    d50_shift_um   = "d₅₀ Shift (µm)"
  )

  y_label <- metric_labels[[metric]]
  if (is.null(y_label)) {
    stop("metric must be one of: W1_micrometers, W1_normalized, d50_shift_um", call. = FALSE)
  }

  required_cols <- c("formulation", "device_resistance", "pressure_drop", metric)
  missing_cols <- setdiff(required_cols, names(w1_results))
  if (length(missing_cols) > 0) {
    stop("w1_results is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # ---- levels + ordering ----
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

  device_labels <- stats::setNames(
    stringr::str_to_title(stringr::str_replace_all(device_levels, "_", " ")),
    device_levels
  )

  pressure_labels <- stats::setNames(
    stringr::str_replace(pressure_levels, "_", " "),
    pressure_levels
  )

  plot_data <- w1_results |>
    dplyr::mutate(
      device_resistance = factor(device_resistance, levels = device_levels),
      pressure_drop     = factor(pressure_drop,     levels = pressure_levels),
      # A human-readable condition key for “single formulation / many conditions” cases
      condition = paste(as.character(device_resistance), as.character(pressure_drop), sep = ", ")
    )

  n_formulations <- dplyr::n_distinct(plot_data$formulation)
  n_devices      <- length(device_levels)
  n_pressures    <- length(pressure_levels)

  # ---- sorting ----
  if (isTRUE(sort_by)) {
    if (n_formulations > 1) {
      order_levels <- plot_data |>
        dplyr::group_by(formulation) |>
        dplyr::summarise(mean_metric = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(mean_metric)) |>
        dplyr::pull(formulation)

      plot_data <- plot_data |>
        dplyr::mutate(formulation = factor(formulation, levels = order_levels))
    } else {
      # one formulation -> sort conditions instead
      order_levels <- plot_data |>
        dplyr::group_by(condition) |>
        dplyr::summarise(mean_metric = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(mean_metric)) |>
        dplyr::pull(condition)

      plot_data <- plot_data |>
        dplyr::mutate(condition = factor(condition, levels = order_levels))
    }
  }

  # ---- AUTO facet decision ----
  if (facet_mode == "auto") {
    facet_mode <- if (n_devices > 1 && n_pressures > 1) {
      "grid"
    } else if (n_devices > 1) {
      "device"
    } else if (n_pressures > 1) {
      "pressure"
    } else {
      "none"
    }
  }

  # ---- choose x-axis based on dataset shape ----
  x_var <- if (n_formulations > 1) "formulation" else "condition"
  x_lab <- if (n_formulations > 1) "Formulation" else "Condition (device, pressure)"

  # ---- build plot ----
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_var]], y = .data[[metric]])
  ) +
    ggplot2::geom_col(fill = "#1F78B4", color = "black", linewidth = 0.3) +
    ggplot2::labs(
      x = x_lab,
      y = y_label,
      title = "Dispersibility Ranking (W₁)",
      subtitle = sprintf(
        "Data shape: %d formulation(s), %d device level(s), %d pressure level(s)",
        n_formulations, n_devices, n_pressures
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

  # ---- facets (optional) ----
  if (facet_mode == "grid") {
    p <- p + ggplot2::facet_grid(
      device_resistance ~ pressure_drop,
      labeller = ggplot2::labeller(
        device_resistance = device_labels,
        pressure_drop     = pressure_labels
      )
    )
  } else if (facet_mode == "device") {
    p <- p + ggplot2::facet_wrap(
      ~ device_resistance,
      nrow = 1,
      labeller = ggplot2::labeller(device_resistance = device_labels)
    )
  } else if (facet_mode == "pressure") {
    p <- p + ggplot2::facet_wrap(
      ~ pressure_drop,
      nrow = 1,
      labeller = ggplot2::labeller(pressure_drop = pressure_labels)
    )
  }

  # ---- sizing defaults ----
  if (is.null(width)) {
    width <- if (x_var == "formulation") max(10, 4 + n_formulations * 1.2) else max(10, 4 + dplyr::n_distinct(plot_data$condition) * 0.8)
  }
  if (is.null(height)) {
    height <- if (facet_mode == "grid") max(6, 2 + n_devices * 2.0) else 6
  }

  if (isTRUE(save_plot)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    output_path <- file.path(output_dir, filename)
    ggplot2::ggsave(output_path, p, width = width, height = height, dpi = dpi)
    cat(sprintf("✓ Saved: %s (facet_mode = %s)\n", output_path, facet_mode))
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
    facet_mode = "auto",   # <--- key change
    save_plot = TRUE,
    output_dir = output_dir,
    filename = "w1_ranking.pdf"
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
