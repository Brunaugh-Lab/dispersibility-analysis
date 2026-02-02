# ==============================================================================
# 01_data_import.R — Laser Diffraction Data Import + Standardization
#
# What this does
#   - Reads Sympatec PAQXOS CSV exports for RODOS (reference) and INHALER (test)
#   - Standardizes columns and extracts metadata needed for downstream W₁ analysis
#   - Writes a single tidy CSV for the rest of the pipeline
#
# Inputs expected
#   data_dir/
#     RODOS/<formulation>/*.csv          (replicates; no condition metadata in file)
#     INHALER/*.csv                      (each file contains condition metadata)
#
# Output
#   <output_dir>/<output_filename>
#   Required columns include:
#     particle_size_um, q3_percent, q3_cdf,
#     formulation, module, replicate,
#     device_resistance, pressure_drop_clean,
#     measurement_time, source_file
#
# Key design choices (do not change without intent)
#   - q3_cdf = q3_percent / 100  (W₁ operates on CDFs in [0,1])
#   - RODOS is treated as formulation-level reference (condition-independent)
#   - INHALER replicates are file-level; replicate labels are assigned deterministically
#
# How to run
#   source("scripts/01_data_import.R")   # loads functions
#   data <- run_data_import("data")      # executes import + validation + write
#
# ==============================================================================


# ==============================================================================
# Dependencies (no library(); use pkg::fn everywhere)
# ==============================================================================

.required_packages <- c("readr", "dplyr", "tidyr", "purrr", "stringr", "janitor")

.missing_packages <- .required_packages[
  !vapply(.required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(.missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(.missing_packages, collapse = ", "),
    "\nInstall with:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', .missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}


# ==============================================================================
# CORE FUNCTION: Read and standardize laser diffraction data
# ==============================================================================

read_ld_data_from_structure <- function(
  data_directory,
  formulation_pattern = ".*",            # Default: use entire folder name
  replicate_pattern   = "[Rr]ep_?\\d+",  # rep1, Rep1, rep_1, Rep_1
  skip_rows           = 2,
  module_folders      = c("inhaler", "INHALER", "rodos", "RODOS"),  # kept for API compatibility
  output_dir          = file.path(data_directory, "tidy"),
  save_output         = TRUE,
  output_filename     = "standardized_data_with_conditions.csv",
  verbose             = TRUE
) {

  # Validate inputs
  if (!dir.exists(data_directory)) {
    stop("Data directory does not exist: ", data_directory, call. = FALSE)
  }

  # Create output directory if it doesn't exist
  if (save_output && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory: ", output_dir, "\n", sep = "")
  }

  # Find all CSV files recursively
  file_paths <- list.files(
    path = data_directory,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(file_paths) == 0) {
    stop("No CSV files found in ", data_directory, call. = FALSE)
  }

  # Normalize paths so detection works on Windows too
  file_paths <- normalizePath(file_paths, winslash = "/", mustWork = FALSE)

  if (verbose) {
    cat("\n========================================================================\n")
    cat("READING LASER DIFFRACTION DATA\n")
    cat("========================================================================\n")
    cat("Data directory: ", data_directory, "\n", sep = "")
    cat("CSV files found: ", length(file_paths), "\n", sep = "")
    cat("Skip rows: ", skip_rows, "\n", sep = "")
    if (save_output) {
      cat("Output directory: ", output_dir, "\n", sep = "")
      cat("Output file: ", file.path(output_dir, output_filename), "\n", sep = "")
    }
    cat("------------------------------------------------------------------------\n\n")
  }

  # ---------------------------------------------------------------------------
  # Split INHALER vs RODOS by path
  # ---------------------------------------------------------------------------
  inhaler_files <- file_paths[stringr::str_detect(file_paths, "(?i)(^|/)inhaler(/|$)")]
  rodos_files   <- file_paths[stringr::str_detect(file_paths, "(?i)(^|/)rodos(/|$)")]

  # If a file matches neither, keep it out (public-safe behavior)
  other_files <- setdiff(file_paths, c(inhaler_files, rodos_files))
  if (verbose && length(other_files) > 0) {
    cat("NOTE: Skipping CSVs not under INHALER/ or RODOS/:\n")
    cat(paste0("  - ", other_files, collapse = "\n"), "\n\n")
  }

  # ---------------------------------------------------------------------------
  # INHALER: read metadata rows (first 2 rows) + data (skip 2 rows)
  # ---------------------------------------------------------------------------
  if (length(inhaler_files) > 0) {

    inhaler_metadata <- purrr::map_dfr(inhaler_files, function(file) {
      headers <- readr::read_csv(file, n_max = 1, col_names = FALSE, show_col_types = FALSE)
      values  <- readr::read_csv(file, skip = 1, n_max = 1, col_names = FALSE, show_col_types = FALSE)

      min_cols <- min(ncol(headers), ncol(values))
      headers <- headers[, 1:min_cols, drop = FALSE]
      values  <- values[, 1:min_cols, drop = FALSE]

      # Name the value row using header row
      metadata_row <- as.list(values)
      names(metadata_row) <- as.character(headers[1, ])

      metadata_row$source_file <- file
      dplyr::as_tibble(metadata_row)
    })

    inhaler_data <- readr::read_csv(
      inhaler_files,
      id = "source_file",
      skip = 2,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    ) |>
      janitor::clean_names()

    inhaler_combined <- inhaler_data |>
      dplyr::left_join(inhaler_metadata, by = "source_file")

    if (verbose && nrow(inhaler_combined) > 0) {
      cat("DEBUG: INHALER columns after joining:\n")
      cat(paste(names(inhaler_combined), collapse = ", "), "\n\n")
    }

  } else {
    inhaler_combined <- dplyr::tibble()
  }

  # ---------------------------------------------------------------------------
  # RODOS: read normally (no metadata in file; folder encodes formulation)
  # ---------------------------------------------------------------------------
  if (length(rodos_files) > 0) {
    rodos_combined <- readr::read_csv(
      rodos_files,
      id = "source_file",
      skip = skip_rows,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    ) |>
      janitor::clean_names()
  } else {
    rodos_combined <- dplyr::tibble()
  }

  # ---------------------------------------------------------------------------
  # Combine + standardize numeric columns
  # ---------------------------------------------------------------------------
  combined_data <- dplyr::bind_rows(inhaler_combined, rodos_combined)

  # Basic contract checks (public-facing helpful failures)
  required_raw <- c("xo_mm", "q3_percent")
  missing_raw <- setdiff(required_raw, names(combined_data))
  if (length(missing_raw) > 0) {
    stop(
      "Missing expected PAQXOS columns: ", paste(missing_raw, collapse = ", "),
      "\nCheck your PAQXOS export format and/or skip_rows setting.",
      call. = FALSE
    )
  }

  combined_data <- combined_data |>
    dplyr::mutate(
      particle_size_um = suppressWarnings(as.numeric(xo_mm)),
      q3_percent       = suppressWarnings(as.numeric(q3_percent)),
      q3_cdf           = q3_percent / 100
    ) |>
    dplyr::filter(!is.na(particle_size_um))

  # ---------------------------------------------------------------------------
  # Extract metadata (INHALER from embedded columns; RODOS from folder structure)
  # ---------------------------------------------------------------------------
  # Normalize source_file path separator already done above.

  combined_data <- combined_data |>
    dplyr::mutate(
      is_inhaler = stringr::str_detect(source_file, "(?i)(^|/)inhaler(/|$)"),
      is_rodos   = stringr::str_detect(source_file, "(?i)(^|/)rodos(/|$)"),

      formulation = dplyr::case_when(
        is_inhaler ~ formulation_id,
        TRUE ~ stringr::str_extract(basename(dirname(source_file)), formulation_pattern)
      ),

      module = dplyr::case_when(
        is_inhaler ~ "INHALER",
        is_rodos   ~ "RODOS",
        TRUE       ~ "UNKNOWN"
      ),

      device_resistance = dplyr::case_when(
        is_inhaler & !is.na(Device) & stringr::str_detect(Device, "(?i)low")    ~ "low",
        is_inhaler & !is.na(Device) & stringr::str_detect(Device, "(?i)medium") ~ "medium",
        is_inhaler & !is.na(Device) & stringr::str_detect(Device, "(?i)high")   ~ "high",
        !is_inhaler ~ "reference",
        TRUE        ~ "unknown"
      ),

      pressure_drop_clean = dplyr::case_when(
        is_inhaler ~ stringr::str_extract(
          dplyr::coalesce(pressure_drop, `pressure-drop`),
          "\\d+"
        ),
        !is_inhaler ~ "reference",
        TRUE        ~ "unknown"
      ),

      measurement_time = dplyr::case_when(
        is_inhaler ~ as.character(Time),
        TRUE       ~ NA_character_
      ),

      filename = basename(source_file),

      replicate = dplyr::case_when(
        !is_inhaler ~ stringr::str_extract(filename, replicate_pattern),
        TRUE        ~ NA_character_
      )
    )

  # ---------------------------------------------------------------------------
  # Assign INHALER replicate labels deterministically (file-level)
  # ---------------------------------------------------------------------------
  combined_data <- combined_data |>
    dplyr::group_by(formulation, module, device_resistance, pressure_drop_clean) |>
    dplyr::mutate(
      # rank unique files within each condition; stable + deterministic
      .file_rank = dplyr::dense_rank(source_file),
      replicate = dplyr::case_when(
        is_inhaler ~ paste0("rep", .file_rank),
        TRUE       ~ replicate
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      particle_size_um,
      q3_percent,
      q3_cdf,
      formulation,
      module,
      device_resistance,
      pressure_drop_clean,
      replicate,
      measurement_time,
      source_file
    ) |>
    dplyr::mutate(replicate = tolower(replicate))

  # ---------------------------------------------------------------------------
  # Validate extraction (warnings, not hard stops)
  # ---------------------------------------------------------------------------
  if (any(is.na(combined_data$formulation))) {
    warning("Some files have NA formulation - check formulation_pattern or INHALER formulation_id column")
  }
  if (any(is.na(combined_data$module))) {
    warning("Some files have NA module - check directory structure")
  }
  if (any(is.na(combined_data$replicate))) {
    warning("Some files have NA replicate - check replicate_pattern (RODOS) or INHALER file parsing")
  }
  if (any(is.na(combined_data$device_resistance))) {
    warning("Some INHALER files have NA device_resistance - check Device column")
  }
  if (any(is.na(combined_data$pressure_drop_clean))) {
    warning("Some INHALER files have NA pressure_drop - check pressure_drop column")
  }

  # ---------------------------------------------------------------------------
  # Print summary
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("Data extraction summary:\n")
    cat("\nFormulations found: ", dplyr::n_distinct(combined_data$formulation), "\n", sep = "")
    print(unique(combined_data$formulation))

    cat("\nModules found: ", dplyr::n_distinct(combined_data$module), "\n", sep = "")
    print(unique(combined_data$module))

    cat("\nReplicates found: ", dplyr::n_distinct(combined_data$replicate), "\n", sep = "")
    print(unique(combined_data$replicate))

    cat("\nDevice resistances found: ", dplyr::n_distinct(combined_data$device_resistance), "\n", sep = "")
    print(unique(combined_data$device_resistance))

    cat("\nPressure drops found: ", dplyr::n_distinct(combined_data$pressure_drop_clean), "\n", sep = "")
    print(unique(combined_data$pressure_drop_clean))

    cat("\nFiles per formulation-module combination:\n")
    print(
      combined_data |>
        dplyr::distinct(source_file, formulation, module) |>
        dplyr::count(formulation, module) |>
        tidyr::pivot_wider(names_from = module, values_from = n, values_fill = 0)
    )

    cat("\n========================================================================\n")
    cat("DATA IMPORT COMPLETE\n")
    cat("Total rows: ", nrow(combined_data), "\n", sep = "")
    cat("Formulations: ", dplyr::n_distinct(combined_data$formulation), "\n", sep = "")
    cat("Modules: ", dplyr::n_distinct(combined_data$module), "\n", sep = "")
    cat("Files processed: ", dplyr::n_distinct(combined_data$source_file), "\n", sep = "")
    cat("========================================================================\n\n")
  }

  # ---------------------------------------------------------------------------
  # Save output
  # ---------------------------------------------------------------------------
  if (save_output) {
    output_path <- file.path(output_dir, output_filename)
    readr::write_csv(combined_data, output_path)

    if (verbose) {
      cat("✓ Standardized data saved to: ", output_path, "\n", sep = "")
      cat("  File size: ", format(utils::object.size(combined_data), units = "MB"), "\n", sep = "")
      cat("  This file can be loaded by subsequent analysis scripts\n\n")
    }
  }

  combined_data
}


# ==============================================================================
# HELPER FUNCTION: Validate data structure
# ==============================================================================

validate_ld_data <- function(data, check_replicates = TRUE, min_replicates = 3) {

  cat("\n========================================================================\n")
  cat("VALIDATING DATA STRUCTURE\n")
  cat("========================================================================\n")

  all_valid <- TRUE

  required_cols <- c(
    "particle_size_um", "q3_percent", "q3_cdf",
    "formulation", "module", "device_resistance",
    "pressure_drop_clean", "replicate", "source_file"
  )

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    cat("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "), "\n", sep = "")
    all_valid <- FALSE
  } else {
    cat("✓ All required columns present\n")
  }

  na_counts <- data |>
    dplyr::summarise(
      dplyr::across(
        c(particle_size_um, q3_percent, formulation, module,
          device_resistance, pressure_drop_clean, replicate),
        ~ sum(is.na(.x))
      )
    )

  if (any(na_counts > 0)) {
    cat("\nWARNING: NA values detected:\n")
    print(na_counts)
    all_valid <- FALSE
  } else {
    cat("✓ No NA values in key columns\n")
  }

  if (any(data$q3_cdf < 0, na.rm = TRUE) || any(data$q3_cdf > 1, na.rm = TRUE)) {
    cat("\nWARNING: q3_cdf values outside [0,1] range\n")
    all_valid <- FALSE
  } else {
    cat("✓ CDF values within [0,1] range\n")
  }

  if (check_replicates) {
    replicate_counts <- data |>
      dplyr::distinct(source_file, formulation, module, replicate) |>
      dplyr::count(formulation, module, name = "n_replicates")

    if (any(replicate_counts$n_replicates < min_replicates)) {
      cat("\nWARNING: Some conditions have fewer than ", min_replicates, " replicates:\n", sep = "")
      print(replicate_counts |>
        dplyr::filter(n_replicates < min_replicates)
      )
      all_valid <- FALSE
    } else {
      cat("✓ All conditions have ≥ ", min_replicates, " replicates\n", sep = "")
    }
  }

  duplicate_files <- data |>
    dplyr::group_by(source_file, particle_size_um) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::distinct(source_file) |>
    dplyr::ungroup()

  if (nrow(duplicate_files) > 0) {
    cat("\nWARNING: Duplicate entries detected for some files\n")
    all_valid <- FALSE
  } else {
    cat("✓ No duplicate file entries\n")
  }

  cat("------------------------------------------------------------------------\n")
  if (all_valid) {
    cat("VALIDATION PASSED: Data structure is ready for analysis\n")
  } else {
    cat("VALIDATION FAILED: Please review warnings above\n")
  }
  cat("========================================================================\n\n")

  invisible(all_valid)
}


# ==============================================================================
# CONVENIENCE FUNCTION: Load previously saved standardized data
# ==============================================================================

load_standardized_data <- function(
  processed_dir = "data/tidy",
  filename = "standardized_data_with_conditions.csv",
  verbose = TRUE
) {

  file_path <- file.path(processed_dir, filename)

  if (!file.exists(file_path)) {
    stop(
      "Standardized data file not found: ", file_path,
      "\nRun read_ld_data_from_structure() first to create this file.",
      call. = FALSE
    )
  }

  if (verbose) cat("Loading standardized data from: ", file_path, "\n", sep = "")

  data <- readr::read_csv(file_path, show_col_types = FALSE)

  if (verbose) {
    cat("✓ Loaded ", nrow(data), " rows\n", sep = "")
    cat("  Formulations: ", dplyr::n_distinct(data$formulation), "\n", sep = "")
    cat("  Modules: ", paste(unique(data$module), collapse = ", "), "\n\n", sep = "")
  }

  data
}


# ==============================================================================
# CONVENIENCE FUNCTION: Run import with project defaults
# ==============================================================================

run_data_import <- function(
  data_directory = "data",
  formulation_pattern = ".*",
  replicate_pattern = "[Rr]ep_?\\d+",
  verbose = TRUE
) {

  data <- read_ld_data_from_structure(
    data_directory = data_directory,
    formulation_pattern = formulation_pattern,
    replicate_pattern = replicate_pattern,
    output_dir = file.path(data_directory, "tidy"),
    save_output = TRUE,
    verbose = verbose
  )

  validate_ld_data(data, check_replicates = TRUE, min_replicates = 3)

  data
}


# ==============================================================================
# Example usage (no auto-execution)
# ==============================================================================
#
# source("scripts/01_data_import.R")
# data <- run_data_import("data")
#
# # Fast reload:
# # data <- load_standardized_data(
# #   processed_dir = "data/tidy",
# #   filename = "standardized_data_with_conditions.csv"
# # )
#
