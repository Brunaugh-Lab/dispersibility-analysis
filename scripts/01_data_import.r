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
#   - INHALER replicates are identified at the file level; replicate labels are derived
#     deterministically (see replicate assignment section below)
#
# How to run
#   source("scripts/01_data_import.R")
#   # or call run_data_import("data")
#
# Configuration
#   - Default data directory: data/
#   - Default output: data/tidy/standardized_data_with_conditions.csv
#   - Edit defaults inside run_data_import() if your project layout differs
#
# Maintainers
#   - Brunaugh Lab, University of Michigan
#   - Code: Grace Xia, Ashlee Brunaugh
#   - Last updated: 2026-02-02
# ==============================================================================

library(tidyverse)
library(janitor)

AUTO_RUN <- FALSE  # set TRUE for lab convenience; keep FALSE for public use
DEFAULT_DATA_DIR <- "data"

# ==============================================================================
# CORE FUNCTION: Read and standardize laser diffraction data
# ==============================================================================

read_ld_data_from_structure <- function(
    data_directory,
    formulation_pattern = ".*",  # Default: use entire folder name
    replicate_pattern = "[Rr]ep_?\\d+",  # Flexible: rep1, Rep1, rep_1, Rep_1
    skip_rows = 2,
    module_folders = c("inhaler", "INHALER", "rodos", "RODOS"),
    output_dir = file.path(data_directory, "tidy"),
    save_output = TRUE,
    output_filename = "standardized_data_with_conditions.csv",
    verbose = TRUE
) {

  # Validate inputs
  if (!dir.exists(data_directory)) {
    stop("Data directory does not exist: ", data_directory)
  }

  # Create output directory if it doesn't exist
  if (save_output && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) {
      cat("Created output directory:", output_dir, "\n")
    }
  }

  # Find all CSV files recursively
  file_paths <- list.files(
    path = data_directory,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(file_paths) == 0) {
    stop("No CSV files found in ", data_directory)
  }

  if (verbose) {
    cat("\n========================================================================\n")
    cat("READING LASER DIFFRACTION DATA\n")
    cat("========================================================================\n")
    cat("Data directory:", data_directory, "\n")
    cat("CSV files found:", length(file_paths), "\n")
    cat("Skip rows:", skip_rows, "\n")
    if (save_output) {
      cat("Output directory:", output_dir, "\n")
      cat("Output file:", file.path(output_dir, output_filename), "\n")
    }
    cat("------------------------------------------------------------------------\n\n")
  }

  # ... everything else unchanged ...
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
    stop("Standardized data file not found: ", file_path,
         "\nRun read_ld_data_from_structure() first to create this file.")
  }

  if (verbose) {
    cat("Loading standardized data from:", file_path, "\n")
  }

  data <- read_csv(file_path, show_col_types = FALSE)

  if (verbose) {
    cat("✓ Loaded", nrow(data), "rows\n")
    cat("  Formulations:", n_distinct(data$formulation), "\n")
    cat("  Modules:", paste(unique(data$module), collapse = ", "), "\n\n")
  }

  return(data)
}

# ==============================================================================
# CONVENIENCE FUNCTION: Run import with project defaults
# ==============================================================================

run_data_import <- function(
    data_directory = "data",
    formulation_pattern = ".*",  # Use entire folder name
    replicate_pattern = "[Rr]ep_?\\d+",  # Flexible replicate matching
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

  return(data)
}

# ==============================================================================
# AUTO-EXECUTION (optional)
# ==============================================================================

if (AUTO_RUN) {
  if (dir.exists(DEFAULT_DATA_DIR)) {

    cat("\n========================================================================\n")
    cat("AUTO-RUNNING DATA IMPORT\n")
    cat("========================================================================\n")
    cat("Reading from: ", DEFAULT_DATA_DIR, "/\n", sep = "")
    cat("Saving to: ", file.path(DEFAULT_DATA_DIR, "tidy", "standardized_data_with_conditions.csv"), "\n", sep = "")
    cat("------------------------------------------------------------------------\n")

    .standardized_data <- run_data_import(DEFAULT_DATA_DIR, verbose = TRUE)

    cat("\n========================================================================\n")
    cat("IMPORT COMPLETE\n")
    cat("========================================================================\n")
    cat("Next step: Run Wasserstein analysis\n")
    cat("  source('scripts/02_wasserstein_core.R')\n")
    cat("------------------------------------------------------------------------\n")
    cat("To reload data later without re-importing:\n")
    cat("  data <- load_standardized_data(processed_dir = 'data/tidy',\n")
    cat("                                filename = 'standardized_data_with_conditions.csv')\n")
    cat("========================================================================\n\n")

  } else {
    cat("\n========================================================================\n")
    cat("DATA IMPORT - WAITING FOR DATA FOLDER\n")
    cat("========================================================================\n")
    cat("Data directory not found: ", DEFAULT_DATA_DIR, "/\n", sep = "")
    cat("Create it and re-run.\n")
    cat("========================================================================\n\n")
  }
}
