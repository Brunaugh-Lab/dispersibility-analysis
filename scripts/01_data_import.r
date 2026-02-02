# ==============================================================================
# 01_data_import.R — Laser Diffraction Data Import + Standardization
#
# What this does
#   - Reads Sympatec PAQXOS CSV exports for RODOS (reference) and INHALER (test)
#   - Reads the 2-row PAQXOS metadata block for BOTH dispersers (consistent ingestion)
#   - Determines disperser type primarily from CSV metadata ("Dispersing system")
#   - AUTO-DETECTS the data-block header line per file (skip_rows becomes per-file)
#   - Standardizes columns and extracts metadata needed for downstream W₁ analysis
#   - Writes a single tidy CSV for the rest of the pipeline
#
# Inputs expected
#   data_dir/  (any structure; subfolders optional)
#     *.csv
#
# Output
#   <data_dir>/tidy/standardized_data_with_conditions.csv   (default)
#
# How to run
#   source("scripts/01_data_import.R")   # loads functions
#   data <- run_data_import("data")      # executes import + validation + write
# ==============================================================================


# ==============================================================================
# Dependencies
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
# INTERNAL: Detect PAQXOS data-table header row for a single file
#   - Returns skip value for readr::read_csv() so that the distribution header
#     (e.g., "xo / µm,Q₃ / %,...") becomes the column header line.
# ==============================================================================
.detect_paqxos_skip <- function(file, max_lines = 40, default_skip = 2, verbose = FALSE) {

  lines <- tryCatch(
    readLines(file, n = max_lines, warn = FALSE),
    error = function(e) character(0)
  )

  if (length(lines) == 0) {
    if (verbose) cat("WARN: Could not read lines for skip detection: ", file, "\n", sep = "")
    return(default_skip)
  }

  # Normalize: lowercase, trim whitespace
  l <- stringr::str_trim(tolower(lines))

  # Stronger heuristic: true PAQXOS header must contain BOTH xo and q3 on the same line,
  # plus commas. This avoids false positives from metadata rows.
  has_xo     <- stringr::str_detect(l, "\\bxo\\b")
  has_q3     <- stringr::str_detect(l, "\\bq\\s*3\\b|q₃|q3")
  has_comma <- stringr::str_detect(l, ",")

  is_header <- stringr::str_detect(l, "^xo\\s*([,/]|\\s)") &
               has_xo & has_q3 & has_comma

  idx <- which(is_header)[1]

  if (is.na(idx)) {
    if (verbose) {
      cat(
        "WARN: Could not auto-detect PAQXOS header row; using default skip=",
        default_skip, " for ", file, "\n", sep = ""
      )
    }
    return(default_skip)
  }

  # readr::read_csv(skip = k-1) will read line k as header row
  skip <- max(idx - 1, 0)
  skip
}


# ==============================================================================
# INTERNAL: Read 2-row PAQXOS metadata block (row 1 headers, row 2 values)
#   - Returns one row per file
#   - All metadata kept as character
#   - Column names cleaned for stable downstream use
# ==============================================================================
.read_paqxos_metadata <- function(files) {
  purrr::map_dfr(files, function(file) {

    headers <- readr::read_csv(
      file,
      n_max = 1,
      col_names = FALSE,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    )

    values <- readr::read_csv(
      file,
      skip = 1,
      n_max = 1,
      col_names = FALSE,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    )

    min_cols <- min(ncol(headers), ncol(values))
    headers <- headers[, 1:min_cols, drop = FALSE]
    values  <- values[, 1:min_cols, drop = FALSE]

    row <- as.list(values)
    names(row) <- as.character(headers[1, ])

    out <- dplyr::as_tibble(row)
    out$source_file <- file

    janitor::clean_names(out, replace = c("µ" = "u", "μ" = "u", "\u00b5" = "u"))
  })
}

# ==============================================================================
# INTERNAL: Standardize PAQXOS column name variants after janitor::clean_names()
# ==============================================================================
.standardize_ld_columns <- function(df) {

  # size bin column (xo / µm) can normalize differently across systems
  size_candidates <- c("xo_mm", "xo_um", "xo_m", "xo")
  size_found <- size_candidates[size_candidates %in% names(df)][1]
  if (!is.na(size_found) && size_found != "xo_mm") {
    df <- dplyr::rename(df, xo_mm = dplyr::all_of(size_found))
  }

  # Q3 percent column (Q₃ / %) sometimes becomes q_3_percent
  q3_candidates <- c("q3_percent", "q_3_percent", "q3_pct", "q_3_pct")
  q3_found <- q3_candidates[q3_candidates %in% names(df)][1]
  if (!is.na(q3_found) && q3_found != "q3_percent") {
    df <- dplyr::rename(df, q3_percent = dplyr::all_of(q3_found))
  }

  df
}

# ==============================================================================
# INTERNAL: Read PAQXOS distribution data block for each file (auto-skip)
# ==============================================================================
.read_paqxos_data_block <- function(files, default_skip = 2, verbose = FALSE) {

  bad_files <- character(0)

  out <- purrr::map_dfr(files, function(file) {

    skip <- .detect_paqxos_skip(file, default_skip = default_skip, verbose = verbose)

    df <- suppressMessages(
      readr::read_csv(
        file,
        skip = skip,
        col_types = readr::cols(.default = "c"),
        show_col_types = FALSE,
        name_repair = "minimal"
      )
    ) |>
      janitor::clean_names(replace = c("µ" = "u", "μ" = "u", "\u00b5" = "u")) |>
      .standardize_ld_columns()

    # Basic per-file contract check
    required <- c("xo_mm", "q3_percent")
    if (!all(required %in% names(df))) {
      bad_files <<- c(bad_files, file)
      return(dplyr::tibble())  # skip this file cleanly
    }

    df$source_file <- file
    df
  })

  # Emit a single useful warning with the file list
  if (length(bad_files) > 0) {
    warning(
      "Skipped ", length(bad_files), " CSV(s) that did not contain expected PAQXOS columns ",
      "(xo_mm + q3_percent) after auto-detect.\n",
      "First few:\n  - ",
      paste(utils::head(bad_files, 10), collapse = "\n  - "),
      if (length(bad_files) > 10) "\n  ... (see full list in warnings())" else "",
      call. = FALSE
    )
  }

  out
}


# ==============================================================================
# CORE FUNCTION: Read and standardize laser diffraction data
# ==============================================================================
read_ld_data_from_structure <- function(
  data_directory,
  formulation_pattern = ".*",            # Default: use entire folder name (RODOS fallback)
  replicate_pattern   = "[Rr]ep_?\\d+",  # rep1, Rep1, rep_1, Rep_1 (RODOS fallback)
  skip_rows           = 2,               # used as default fallback if auto-detect fails
  module_folders      = c("inhaler", "INHALER", "rodos", "RODOS"),  # API compatibility
  output_dir          = file.path(data_directory, "tidy"),
  save_output         = TRUE,
  output_filename     = "standardized_data_with_conditions.csv",
  verbose             = TRUE
) {

  if (!dir.exists(data_directory)) {
    stop("Data directory does not exist: ", data_directory, call. = FALSE)
  }

  if (save_output && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory: ", output_dir, "\n", sep = "")
  }

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

  # Exclude generated outputs from being re-read as inputs
  output_dir_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  file_paths <- file_paths[
    !stringr::str_detect(file_paths, paste0("^", stringr::fixed(output_dir_norm), "(/|$)"))
  ]

  # Extra safety: also exclude the output filename anywhere it appears
  file_paths <- file_paths[basename(file_paths) != output_filename]

  if (length(file_paths) == 0) {
    stop(
      "No input PAQXOS CSV files found after excluding output_dir/output file.\n",
      "Check your data_directory and folder structure.",
      call. = FALSE
    )
  }


  if (verbose) {
    cat("\n========================================================================\n")
    cat("READING LASER DIFFRACTION DATA\n")
    cat("========================================================================\n")
    cat("Data directory: ", data_directory, "\n", sep = "")
    cat("CSV files found: ", length(file_paths), "\n", sep = "")
    cat("Default skip_rows fallback: ", skip_rows, "\n", sep = "")
    if (save_output) {
      cat("Output directory: ", output_dir, "\n", sep = "")
      cat("Output file: ", file.path(output_dir, output_filename), "\n", sep = "")
    }
    cat("------------------------------------------------------------------------\n\n")
  }

  # Read metadata for ALL files
  metadata <- .read_paqxos_metadata(file_paths)

  # Read distribution data block for ALL files (auto-detect skip per file)
  data_block <- .read_paqxos_data_block(
    file_paths,
    default_skip = skip_rows,
    verbose = FALSE
  )

  if (nrow(data_block) == 0) {
    stop(
      "No valid PAQXOS distribution tables were read.\n",
      "Auto-detect likely failed for all files, or files are not PAQXOS exports.",
      call. = FALSE
    )
  }

  combined_data <- data_block |>
    dplyr::left_join(metadata, by = "source_file")

  # Determine module primarily from metadata "dispersing_system"
  combined_data <- combined_data |>
    dplyr::mutate(
      .ds = tolower(dplyr::coalesce(dispersing_system, NA_character_)),
      module_from_metadata = dplyr::case_when(
        !is.na(.ds) & stringr::str_detect(.ds, "inhaler") ~ "INHALER",
        !is.na(.ds) & stringr::str_detect(.ds, "rodos")   ~ "RODOS",
        TRUE ~ NA_character_
      ),
      module_from_path = dplyr::case_when(
        stringr::str_detect(source_file, "(?i)(^|/)inhaler(/|$)") ~ "INHALER",
        stringr::str_detect(source_file, "(?i)(^|/)rodos(/|$)")   ~ "RODOS",
        TRUE ~ NA_character_
      ),
      module = dplyr::coalesce(module_from_metadata, module_from_path)
    )

  # Keep only PAQXOS-like RODOS/INHALER exports
  kept <- combined_data |>
    dplyr::distinct(source_file, module) |>
    dplyr::filter(!is.na(module)) |>
    dplyr::pull(source_file)

  skipped <- setdiff(unique(combined_data$source_file), kept)

  if (verbose && length(skipped) > 0) {
    cat("NOTE: Skipping CSVs that do not identify as RODOS/INHALER via metadata or folder:\n")
    cat(paste0("  - ", skipped, collapse = "\n"), "\n\n")
  }

  combined_data <- combined_data |>
    dplyr::filter(source_file %in% kept)

  # Contract checks
  required_raw <- c("xo_mm", "q3_percent")
  missing_raw <- setdiff(required_raw, names(combined_data))
  if (length(missing_raw) > 0) {
    stop(
      "Missing expected PAQXOS columns: ", paste(missing_raw, collapse = ", "),
      "\nAuto-detect may have failed to find the distribution header row for some files.",
      "\nTry increasing max_lines in .detect_paqxos_skip() or inspect one failing export.",
      call. = FALSE
    )
  }

  # Standardize numeric columns
  combined_data <- combined_data |>
    dplyr::mutate(
      particle_size_um = suppressWarnings(as.numeric(xo_mm)),
      q3_percent       = suppressWarnings(as.numeric(q3_percent)),
      q3_cdf           = q3_percent / 100
    ) |>
    dplyr::filter(!is.na(particle_size_um))


  # Ensure expected metadata columns exist (avoids hard-fail on missing columns)
  if (!("pressure_drop" %in% names(combined_data))) combined_data$pressure_drop <- NA_character_
  if (!("device"        %in% names(combined_data))) combined_data$device        <- NA_character_
  if (!("formulation_id" %in% names(combined_data))) combined_data$formulation_id <- NA_character_

  # Extract metadata used downstream
  combined_data <- combined_data |>
    dplyr::mutate(
      is_inhaler = module == "INHALER",
      is_rodos   = module == "RODOS",

      # formulation_id is now the primary source for BOTH dispersers
      .form_id = dplyr::na_if(formulation_id, ""),
      .form_id = dplyr::na_if(.form_id, "NA"),
      .form_id = dplyr::na_if(.form_id, "na"),

      formulation = dplyr::coalesce(
        .form_id,
        # legacy/public fallback: folder name
        stringr::str_extract(basename(dirname(source_file)), formulation_pattern)
      ),

      device_resistance = dplyr::case_when(
        is_inhaler & !is.na(device) & stringr::str_detect(device, "(?i)low")    ~ "low",
        is_inhaler & !is.na(device) & stringr::str_detect(device, "(?i)medium") ~ "medium",
        is_inhaler & !is.na(device) & stringr::str_detect(device, "(?i)high")   ~ "high",
        is_rodos   ~ "reference",
        TRUE       ~ "unknown"
      ),

      pressure_drop_clean = dplyr::case_when(
        is_inhaler ~ stringr::str_extract(pressure_drop, "\\d+"),
        is_rodos   ~ "reference",
        TRUE       ~ "unknown"
      ),

      measurement_time = dplyr::coalesce(
        as.character(time),
        as.character(identifier)
      ),

      filename = basename(source_file),

      replicate = dplyr::case_when(
        is_rodos ~ stringr::str_extract(filename, replicate_pattern),
        TRUE     ~ NA_character_
      )
    ) |>
    dplyr::select(-.form_id)

  # Deterministic INHALER replicate labels
  combined_data <- combined_data |>
    dplyr::group_by(formulation, module, device_resistance, pressure_drop_clean) |>
    dplyr::mutate(
      .time_key = dplyr::if_else(
        !is.na(measurement_time),
        measurement_time,
        source_file
      )
    ) |>
    dplyr::arrange(.time_key, .by_group = TRUE) |>
    dplyr::mutate(
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

  # Warnings
  if (any(is.na(combined_data$formulation))) {
    warning(
      "Some files have NA formulation. Expected: formulation_id in metadata. ",
      "Check PAQXOS export settings or folder fallback pattern."
    )
  }
  if (any(is.na(combined_data$replicate))) {
    warning("Some files have NA replicate - check replicate_pattern (RODOS) or INHALER file parsing.")
  }

  # Summary
  if (verbose) {
    cat("Data extraction summary:\n")
    cat("\nFormulations found: ", dplyr::n_distinct(combined_data$formulation), "\n", sep = "")
    print(unique(combined_data$formulation))

    cat("\nModules found: ", dplyr::n_distinct(combined_data$module), "\n", sep = "")
    print(unique(combined_data$module))

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

  # Save output
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

  na_counts_vec <- unlist(na_counts, use.names = FALSE)

  if (any(na_counts_vec > 0)) {
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
