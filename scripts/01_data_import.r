# ==============================================================================
# 01_data_import.R
# Data Import and Standardization for Laser Diffraction Dispersibility Analysis
#
# Purpose: Enhanced reading of Sympatec PAQXOS CSV exports with support for
#          device operating conditions (resistance/pressure) and metadata
#          extraction from both directory structure and CSV headers.
#
# Auto-execution: Script automatically runs when sourced
#   - Reads all CSV files from data_v2/ directory
#   - Extracts device conditions from INHALER CSV metadata
#   - Assigns replicates by timestamp for INHALER files
#   - Saves to data_v2/tidy/standardized_data_with_conditions.csv
#
# Enhanced Directory Structure:
#   data_v2/
#   ├── RODOS/                 ← Reference measurements
#   │   ├── FormulationA/      ← Folder name = Formulation ID
#   │   │   ├── rep1.csv
#   │   │   ├── rep2.csv
#   │   │   └── rep3.csv
#   │   ├── FormulationB/
#   │   │   ├── rep1.csv
#   │   │   ├── rep2.csv
#   │   │   └── rep3.csv
#   │   └── FormulationC/
#   │       └── ...
#   └── INHALER/               ← Test measurements with device conditions
#       ├── file1.csv          ← Contains metadata: formulation_id, Device, pressure_drop, Time
#       ├── file2.csv          ← Device: RS01-M7-low/medium/high (flexible)
#       ├── file3.csv          ← Pressure: 1_kPa, 2_kPa, 4_kPa (flexible)
#       └── ...fileN.csv       ← Any number of files with various condition combinations
#
# Output Structure (auto-created):
#   data_v2/
#   └── tidy/
#       └── standardized_data_with_conditions.csv
#
# Usage:
#   source("scripts/01_data_import.R")  # That's it!
#
# ==============================================================================

library(tidyverse)
library(janitor)

# ==============================================================================
# CORE FUNCTION: Read and standardize laser diffraction data
# ==============================================================================

#' Read Laser Diffraction CSV Files with Enhanced Metadata Support
#'
#' Supports two data structures:
#' 1. RODOS: Formulation folders with replicate CSVs (existing structure)
#' 2. INHALER: Flat folder with metadata inside each CSV file
#'
#' @param data_directory Path to directory containing RODOS and INHALER subdirectories
#'   (e.g., "data_v2/")
#' @param formulation_pattern Regex pattern to extract formulation ID from
#'   folder name. Default: ".*" (uses entire folder name - RECOMMENDED)
#'   Only customize if you need to extract a portion:
#'   - "\\d+_IMT" : Extracts "132067_IMT" from "132067_IMT_batch1"
#'   - "Form[A-Z]" : Extracts "FormA" from "FormA_replicate_set"
#' @param replicate_pattern Regex pattern to extract replicate ID from filename.
#'   Default: "[Rr]ep_?\\d+" matches rep1, Rep1, rep_1, Rep_1
#' @param skip_rows Number of header rows to skip in CSV files.
#'   Default: 2 (standard for Sympatec PAQXOS exports)
#' @param module_folders Character vector of folder names that indicate dispersion
#'   modules. Default: c("inhaler", "INHALER", "rodos", "RODOS")
#'   Function will standardize these to uppercase for consistency.
#' @param output_dir Directory to save tidy data. Default: "data/tidy/"
#' @param save_output Should standardized data be saved to CSV? Default: TRUE
#' @param output_filename Name of output file. Default: "standardized_data.csv"
#' @param verbose Print progress messages? Default: TRUE
#'
#' @return Tibble with standardized columns:
#'   - particle_size_um: Particle diameter in micrometers (from xo column)
#'   - q3_percent: Cumulative volume distribution, 0-100%
#'   - q3_cdf: Cumulative distribution function, 0-1 (for Wasserstein calculation)
#'   - formulation: Formulation identifier (from folder name)
#'   - module: Dispersion module (INHALER or RODOS, standardized to uppercase)
#'   - replicate: Replicate identifier (auto-standardized to lowercase)
#'   - source_file: Full path to original CSV file for traceability
#'
#' @details
#' This function automatically:
#' - Creates data/tidy/ directory if needed
#' - Saves standardized_data.csv for downstream scripts
#' - Standardizes replicate names to lowercase
#' - Extracts formulation ID from folder name (entire name by default)
#'
#' @examples
#' # Recommended - uses entire folder name as formulation ID
#' data <- read_ld_data_from_structure("data/")
#'
#' # Only needed if extracting portion of folder name
#' data <- read_ld_data_from_structure(
#'   "data/",
#'   formulation_pattern = "\\d+_IMT"
#' )
#'
read_ld_data_from_structure <- function(
    data_directory,
    formulation_pattern = ".*",  # Default: use entire folder name
    replicate_pattern = "[Rr]ep_?\\d+",  # Flexible: rep1, Rep1, rep_1, Rep_1
    skip_rows = 2,
    module_folders = c("inhaler", "INHALER", "rodos", "RODOS"),
    output_dir = "data_v2/tidy",
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

  # Read metadata and data separately for INHALER files, normal reading for RODOS
inhaler_files <- file_paths[str_detect(file_paths, "(?i)inhaler")]
rodos_files <- file_paths[!str_detect(file_paths, "(?i)inhaler")]

# Read INHALER files with metadata extraction
if (length(inhaler_files) > 0) {
  # Read metadata (first 2 rows) for INHALER files
  inhaler_metadata <- map_dfr(inhaler_files, function(file) {
    headers <- read_csv(file, n_max = 1, col_names = FALSE, show_col_types = FALSE)
    values <- read_csv(file, skip = 1, n_max = 1, col_names = FALSE, show_col_types = FALSE)

    # Ensure both have same number of columns
    min_cols <- min(ncol(headers), ncol(values))
  headers <- headers[, 1:min_cols]
  values <- values[, 1:min_cols]

  # Combine headers and values
  metadata_row <- as.list(values)
  names(metadata_row) <- as.character(headers[1,])

  metadata_row$source_file <- file
  return(as_tibble(metadata_row))
  })

  # Read data portion (skip metadata rows) for INHALER files
  inhaler_data <- read_csv(
    inhaler_files,
    id = "source_file",
    skip = 2,  # Skip metadata rows
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  ) %>%
    clean_names()

  # Join metadata with data
  inhaler_combined <- inhaler_data %>%
    left_join(inhaler_metadata, by = "source_file")
} else {
  inhaler_combined <- tibble()
}
# Join metadata with data
inhaler_combined <- inhaler_data %>%
  left_join(inhaler_metadata, by = "source_file")

# DEBUG: Check column names after joining
if (verbose && nrow(inhaler_combined) > 0) {
  cat("DEBUG: INHALER columns after joining:\n")
  cat(paste(names(inhaler_combined), collapse = ", "), "\n\n")
}

# Read RODOS files normally (no metadata in CSV)
if (length(rodos_files) > 0) {
  rodos_combined <- read_csv(
    rodos_files,
    id = "source_file",
    skip = skip_rows,
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  ) %>%
    clean_names()
} else {
  rodos_combined <- tibble()
}

# Combine INHALER and RODOS data
combined_data <- bind_rows(inhaler_combined, rodos_combined) %>%
  mutate(
    # Convert size and cumulative distribution to numeric
    particle_size_um = as.numeric(xo_mm),  # xo_mm is particle size in µm
    q3_percent = as.numeric(q3_percent),
    # CRITICAL: Convert Q3 from percent (0-100) to probability (0-1) for CDF
    q3_cdf = q3_percent / 100
  ) %>%
  filter(!is.na(particle_size_um))

  # Extract metadata with enhanced support for INHALER files
combined_data <- combined_data %>%
  mutate(
    # Determine if this is INHALER or RODOS based on path
    # Check if INHALER or RODOS is anywhere in the path
    is_inhaler = str_detect(source_file, "(?i)/inhaler/"),
    is_rodos = str_detect(source_file, "(?i)/rodos/"),

    # Extract metadata differently for INHALER vs RODOS
    formulation = case_when(
    # INHALER: Extract from formulation_id column in CSV
    is_inhaler ~ formulation_id,
    # RODOS: Extract from immediate parent folder (one level up)
    TRUE ~ str_extract(basename(dirname(source_file)), formulation_pattern)
    ),

    module = case_when(
    is_inhaler ~ "INHALER",
    is_rodos ~ "RODOS",
    TRUE ~ "UNKNOWN"
  ),

    # Extract device information for INHALER files
    device_resistance = case_when(
    is_inhaler & str_detect(Device, "(?i)low") ~ "low",
    is_inhaler & str_detect(Device, "(?i)medium") ~ "medium",
    is_inhaler & str_detect(Device, "(?i)high") ~ "high",
      !is_inhaler ~ "reference",  # RODOS is reference
      TRUE ~ "unknown"
    ),

    pressure_drop_clean = case_when(
    is_inhaler ~ str_extract(coalesce(pressure_drop, `pressure-drop`), "\\d+"),  # Handle both naming variants
    !is_inhaler ~ "reference",  # RODOS is reference
    TRUE ~ "unknown"
  ),

    # Extract timestamp for INHALER files
    measurement_time = case_when(
    is_inhaler ~ as.character(Time),  # Convert datetime to character
    TRUE ~ NA_character_
  ),

    # Extract replicate from filename for RODOS only
    filename = basename(source_file),
    replicate = case_when(
      # RODOS: Extract from filename using pattern (existing logic)
      !is_inhaler ~ str_extract(filename, replicate_pattern),
      # INHALER: Will be assigned by timestamp in next step
      TRUE ~ NA_character_
    )
  ) %>%
  # Auto-assign replicate numbers for INHALER files based on timestamps
  group_by(formulation, module, device_resistance, pressure_drop_clean) %>%
  arrange(measurement_time) %>%
  mutate(
    replicate = case_when(
    is_inhaler ~ paste0("rep", dense_rank(source_file)),  # This counts unique files!
    TRUE ~ replicate
  )
) %>%
  ungroup() %>%
  select(
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
  )

# Standardize replicate names to lowercase automatically
combined_data <- combined_data %>%
  mutate(replicate = tolower(replicate))

# Validate extraction
if (any(is.na(combined_data$formulation))) {
  warning("Some files have NA formulation - check formulation_pattern or formulation_id column")
}
if (any(is.na(combined_data$module))) {
  warning("Some files have NA module - check directory structure")
}
if (any(is.na(combined_data$replicate))) {
  warning("Some files have NA replicate - check replicate_pattern or timestamp assignment")
}
if (any(is.na(combined_data$device_resistance))) {
  warning("Some INHALER files have NA device_resistance - check Device column")
}
if (any(is.na(combined_data$pressure_drop_clean))) {
  warning("Some INHALER files have NA pressure_drop - check pressure_drop column")
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

    cat("\nDevice resistances found:", n_distinct(combined_data$device_resistance), "\n")
    print(unique(combined_data$device_resistance))

    cat("\nPressure drops found:", n_distinct(combined_data$pressure_drop_clean), "\n")
    print(unique(combined_data$pressure_drop_clean))

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

  # Save output if requested
  if (save_output) {
    output_path <- file.path(output_dir, output_filename)
    write_csv(combined_data, output_path)

    if (verbose) {
      cat("✓ Standardized data saved to:", output_path, "\n")
      cat("  File size:", format(object.size(combined_data), units = "MB"), "\n")
      cat("  This file can be loaded by subsequent analysis scripts\n\n")
    }
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
                     "formulation", "module", "device_resistance",
                     "pressure_drop_clean", "replicate", "source_file")
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    all_valid <- FALSE
  } else {
    cat("✓ All required columns present\n")
  }

  # Check for NA values in key columns
  na_counts <- data %>%
  summarise(across(c(particle_size_um, q3_percent, formulation, module,
                     device_resistance, pressure_drop_clean, replicate),
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
# CONVENIENCE FUNCTION: Load previously saved standardized data
# ==============================================================================

#' Load Standardized Data from Previous Run
#'
#' Quickly loads the standardized data saved by read_ld_data_from_structure()
#' without re-reading all raw CSV files.
#'
#' @param processed_dir Directory containing tidy data. Default: "data/tidy/"
#' @param filename Name of standardized data file. Default: "standardized_data.csv"
#' @param verbose Print loading message? Default: TRUE
#'
#' @return Tibble with standardized data
#'
#' @examples
#' # Load previously processed data (much faster than re-importing)
#' data <- load_standardized_data()
#'
load_standardized_data <- function(
    processed_dir = "data/tidy",
    filename = "standardized_data.csv",
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

#' Run Data Import with Sensible Defaults
#'
#' Convenience wrapper that uses standard settings:
#' - Entire folder name becomes formulation ID
#' - Flexible replicate matching (rep1, Rep1, rep_1, etc.)
#' - Auto-saves to processed/standardized_data.csv
#' - Auto-validates
#'
#' @param data_directory Path to data folder. Default: "data/"
#' @param formulation_pattern Regex for formulation ID. Default: ".*" (entire folder name)
#' @param replicate_pattern Regex for replicate ID. Default: "[Rr]ep_?\\d+"
#' @param verbose Print progress? Default: TRUE
#'
#' @return Tibble with standardized data
#'
#' @examples
#' # Simple usage with all defaults
#' data <- run_data_import()
#'
#' # Custom data directory
#' data <- run_data_import("raw_data/")
#'
run_data_import <- function(
    data_directory = "data",
    formulation_pattern = ".*",  # Use entire folder name
    replicate_pattern = "[Rr]ep_?\\d+",  # Flexible replicate matching
    verbose = TRUE
) {

  # Run the full import
  data <- read_ld_data_from_structure(
    data_directory = data_directory,
    formulation_pattern = formulation_pattern,
    replicate_pattern = replicate_pattern,
    save_output = TRUE,
    verbose = verbose
  )

  # Validate
  validate_ld_data(data, check_replicates = TRUE, min_replicates = 3)

  return(data)
}


# ==============================================================================
# AUTO-EXECUTION: Run import when script is sourced
# ==============================================================================

# Check if data directory exists
if (dir.exists("data_v2")) {

  cat("\n========================================================================\n")
  cat("AUTO-RUNNING DATA IMPORT\n")
  cat("========================================================================\n")
  cat("Reading from: data_v2/\n")
  cat("Saving to: data_v2/tidy/standardized_data_with_conditions.csv\n")
  cat("------------------------------------------------------------------------\n")

  # Run the import with defaults
  .standardized_data <- run_data_import("data_v2", verbose = TRUE)

  cat("\n========================================================================\n")
  cat("IMPORT COMPLETE - Data saved to data_v2/tidy/standardized_data_with_conditions.csv\n")
  cat("========================================================================\n")
  cat("Next step: Run Wasserstein analysis\n")
  cat("  source('scripts/02_wasserstein_core.R')\n")
  cat("------------------------------------------------------------------------\n")
  cat("To reload data later without re-importing:\n")
  cat("  source('scripts/01_data_import.R')\n")
  cat("  data <- load_standardized_data()\n")
  cat("========================================================================\n\n")

  # Clean up the auto-generated variable (optional)
  # Uncomment if you don't want .standardized_data in the environment
  # rm(.standardized_data)

} else {
  cat("\n========================================================================\n")
  cat("DATA IMPORT - WAITING FOR DATA FOLDER\n")
  cat("========================================================================\n")
  cat("Data directory not found: data/\n")
  cat("\nPlease create a data/ folder with your laser diffraction files:\n")
  cat("  data/\n")
  cat("  ├── FormulationA/\n")
  cat("  │   ├── inhaler/\n")
  cat("  │   │   ├── rep1.csv\n")
  cat("  │   │   ├── rep2.csv\n")
  cat("  │   │   └── rep3.csv\n")
  cat("  │   └── rodos/\n")
  cat("  │       └── ...\n")
  cat("  └── FormulationB/\n")
  cat("      └── ...\n")
  cat("\nFolder names will become formulation IDs.\n")
  cat("Then run this script again:\n")
  cat("  source('scripts/01_data_import.R')\n")
  cat("========================================================================\n\n")
}
