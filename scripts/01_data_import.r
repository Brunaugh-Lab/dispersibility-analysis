# ==============================================================================
# 01_data_import.R
# Data Import and Standardization for Laser Diffraction Dispersibility Analysis
#
# Purpose: Flexible reading of Sympatec PAQXOS CSV exports with metadata
#          extraction from directory structure
#
# Expected Directory Structure:
#   data/
#   ├── FormulationA/
#   │   ├── inhaler/  (or INHALER)
#   │   │   ├── rep1.csv
#   │   │   ├── rep2.csv
#   │   │   └── rep3.csv
#   │   └── rodos/    (or RODOS)
#   │       ├── rep1.csv
#   │       └── ...
#   └── FormulationB/
#       └── ...
#
# ==============================================================================

library(tidyverse)
library(janitor)

# ==============================================================================
# CORE FUNCTION: Read and standardize laser diffraction data
# ==============================================================================

#' Read Laser Diffraction CSV Files with Metadata from Directory Structure
#'
#' @param data_directory Path to directory containing subdirectories organized
#'   by formulation and dispersion module (e.g., "data/", "./raw_data/")
#' @param formulation_pattern Regex pattern to extract formulation ID from
#'   folder name. Default extracts entire folder name. Examples:
#'   - ".*" (default): Uses entire folder name as formulation ID
#'   - "F\\d+" : Extracts F2, F3, F4, etc.
#'   - "\\d+_IMT" : Extracts "132067_IMT", "231067_IMT", etc.
#' @param replicate_pattern Regex pattern to extract replicate ID from filename.
#'   Default: "rep\\d+" extracts rep1, rep2, rep3, etc.
#'   Alternative: "Rep_\\d+" for Rep_1, Rep_2, etc.
#' @param skip_rows Number of header rows to skip in CSV files.
#'   Default: 2 (standard for Sympatec PAQXOS exports)
#' @param module_folders Character vector of folder names that indicate dispersion
#'   modules. Default: c("inhaler", "INHALER", "rodos", "RODOS")
#'   Function will standardize these to uppercase for consistency.
#' @param verbose Print progress messages? Default: TRUE
#'
#' @return Tibble with standardized columns:
#'   - particle_size_um: Particle diameter in micrometers (from xo column)
#'   - q3_percent: Cumulative volume distribution, 0-100%
#'   - q3_cdf: Cumulative distribution function, 0-1 (for Wasserstein calculation)
#'   - formulation: Formulation identifier extracted from folder structure
#'   - module: Dispersion module (INHALER or RODOS, standardized to uppercase)
#'   - replicate: Replicate identifier extracted from filename
#'   - source_file: Full path to original CSV file for traceability
#'
#' @examples
#' # Basic usage with default patterns
#' data <- read_ld_data_from_structure("data/")
#'
#' # Custom formulation pattern (extract only numeric portion)
#' data <- read_ld_data_from_structure(
#'   "data/",
#'   formulation_pattern = "\\d+_IMT"
#' )
#'
#' # Custom replicate pattern for uppercase Rep_N format
#' data <- read_ld_data_from_structure(
#'   "data/",
#'   replicate_pattern = "Rep_\\d+"
#' )
#'
read_ld_data_from_structure <- function(
    data_directory,
    formulation_pattern = ".*",  # Default: use entire folder name
    replicate_pattern = "rep\\d+",
    skip_rows = 2,
    module_folders = c("inhaler", "INHALER", "rodos", "RODOS"),
    verbose = TRUE
) {

  # Validate inputs
  if (!dir.exists(data_directory)) {
    stop("Data directory does not exist: ", data_directory)
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
    cat("------------------------------------------------------------------------\n\n")
  }

  # Read all CSV files
  combined_data <- read_csv(
    file_paths,
    id = "source_file",
    skip = skip_rows,
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  ) %>%
    clean_names() %>%
    mutate(
      # Convert size and cumulative distribution to numeric
      particle_size_um = as.numeric(xo_mm),  # xo_mm is particle size in µm
      q3_percent = as.numeric(q3_percent),
      # CRITICAL: Convert Q3 from percent (0-100) to probability (0-1) for CDF
      q3_cdf = q3_percent / 100
    ) %>%
    filter(!is.na(particle_size_um))

  # Extract metadata from directory structure
  combined_data <- combined_data %>%
    mutate(
      # Extract formulation from parent directory name
      # Path structure: .../FormulationFolder/module/file.csv
      formulation_folder = basename(dirname(dirname(source_file))),
      formulation = str_extract(formulation_folder, formulation_pattern),

      # Extract module from immediate parent directory
      module_folder = basename(dirname(source_file)),
      module = case_when(
        tolower(module_folder) == "inhaler" ~ "INHALER",
        tolower(module_folder) == "rodos" ~ "RODOS",
        TRUE ~ toupper(module_folder)  # Standardize to uppercase
      ),

      # Extract replicate from filename
      filename = basename(source_file),
      replicate = str_extract(filename, replicate_pattern)
    ) %>%
    select(
      particle_size_um,
      q3_percent,
      q3_cdf,
      formulation,
      module,
      replicate,
      source_file
    )

  # Validate extraction
  if (any(is.na(combined_data$formulation))) {
    warning("Some files have NA formulation - check formulation_pattern")
  }
  if (any(is.na(combined_data$module))) {
    warning("Some files have NA module - check directory structure")
  }
  if (any(is.na(combined_data$replicate))) {
    warning("Some files have NA replicate - check replicate_pattern")
  }

  # Print summary
  if (verbose) {
    cat("Data extraction summary:\n")
    cat("\nFormulations found:", n_distinct(combined_data$formulation), "\n")
    print(unique(combined_data$formulation))

    cat("\nModules found:", n_distinct(combined_data$module), "\n")
    print(unique(combined_data$module))

    cat("\nReplicates found:", n_distinct(combined_data$replicate), "\n")
    print(unique(combined_data$replicate))

    cat("\nFiles per formulation-module combination:\n")
    print(
      combined_data %>%
        distinct(source_file, formulation, module) %>%
        count(formulation, module) %>%
        pivot_wider(names_from = module, values_from = n, values_fill = 0)
    )

    cat("\n========================================================================\n")
    cat("DATA IMPORT COMPLETE\n")
    cat("Total rows:", nrow(combined_data), "\n")
    cat("Formulations:", n_distinct(combined_data$formulation), "\n")
    cat("Modules:", n_distinct(combined_data$module), "\n")
    cat("Files processed:", n_distinct(combined_data$source_file), "\n")
    cat("========================================================================\n\n")
  }

  return(combined_data)
}


# ==============================================================================
# HELPER FUNCTION: Validate data structure
# ==============================================================================

#' Validate Laser Diffraction Data Structure
#'
#' Checks that imported data has required columns and reasonable values
#'
#' @param data Tibble from read_ld_data_from_structure()
#' @param check_replicates Should function check for balanced replicates? Default: TRUE
#' @param min_replicates Minimum expected replicates per condition. Default: 3
#'
#' @return Invisibly returns TRUE if validation passes, otherwise prints warnings
#'
validate_ld_data <- function(data, check_replicates = TRUE, min_replicates = 3) {

  cat("\n========================================================================\n")
  cat("VALIDATING DATA STRUCTURE\n")
  cat("========================================================================\n")

  all_valid <- TRUE

  # Check required columns
  required_cols <- c("particle_size_um", "q3_percent", "q3_cdf",
                     "formulation", "module", "replicate", "source_file")
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    all_valid <- FALSE
  } else {
    cat("✓ All required columns present\n")
  }

  # Check for NA values in key columns
  na_counts <- data %>%
    summarise(across(c(particle_size_um, q3_percent, formulation, module, replicate),
                     ~sum(is.na(.))))

  if (any(na_counts > 0)) {
    cat("\nWARNING: NA values detected:\n")
    print(na_counts)
    all_valid <- FALSE
  } else {
    cat("✓ No NA values in key columns\n")
  }

  # Check CDF bounds (should be 0-1)
  if (any(data$q3_cdf < 0, na.rm = TRUE) || any(data$q3_cdf > 1, na.rm = TRUE)) {
    cat("\nWARNING: q3_cdf values outside [0,1] range\n")
    all_valid <- FALSE
  } else {
    cat("✓ CDF values within [0,1] range\n")
  }

  # Check for balanced replicates
  if (check_replicates) {
    replicate_counts <- data %>%
      distinct(source_file, formulation, module, replicate) %>%
      count(formulation, module) %>%
      rename(n_replicates = n)

    if (any(replicate_counts$n_replicates < min_replicates)) {
      cat("\nWARNING: Some conditions have fewer than", min_replicates, "replicates:\n")
      print(replicate_counts %>% filter(n_replicates < min_replicates))
      all_valid <- FALSE
    } else {
      cat("✓ All conditions have ≥", min_replicates, "replicates\n")
    }
  }

  # Check for duplicate files
  duplicate_files <- data %>%
    group_by(source_file, particle_size_um) %>%
    filter(n() > 1) %>%
    distinct(source_file)

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

  return(invisible(all_valid))
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# Uncomment to test with your data structure:
# data <- read_ld_data_from_structure(
#   data_directory = "data/",
#   formulation_pattern = "\\d+_IMT",  # Extracts "132067_IMT", etc.
#   replicate_pattern = "rep\\d+",      # Extracts "rep1", "rep2", etc.
#   verbose = TRUE
# )
#
# validate_ld_data(data, check_replicates = TRUE, min_replicates = 3)
