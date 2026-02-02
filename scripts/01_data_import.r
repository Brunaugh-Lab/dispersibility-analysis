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

  # Robust heuristic:
  #  - header line starts with "xo"
  #  - contains commas (CSV header)
  #  - contains a q-column marker (q, q3, q₃, etc.)
  is_header <- stringr::str_detect(l, "^\\s*xo\\s*[/,]") &
    stringr::str_detect(l, ",") &
    stringr::str_detect(l, "q")

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

    janitor::clean_names(out, replace = c(
      "µ" = "u", "μ" = "u", "\u00b5" = "u",
      "₃" = "3", "³" = "3"
    ))
  })
}

# ==============================================================================
# INTERNAL: Standardize PAQXOS column name variants after janitor::clean_names()
#   FIXED VERSION - more flexible pattern matching
# ==============================================================================
.standardize_ld_columns <- function(df, verbose = FALSE) {

  nms <- names(df)

  if (verbose) {
    cat("DEBUG: Column names before standardization:\n")
    cat("  ", paste(nms, collapse = ", "), "\n", sep = "")
  }

  # --- size bin column (xo / µm) ---
  # janitor::clean_names() converts "xo / µm" to something like "xo_um" or "xo_u_m"
  # Expanded candidates to catch more variants
  size_candidates <- c("xo_mm", "xo_um", "xo_u_m", "xo_m", "xo")
  size_found <- size_candidates[size_candidates %in% nms][1]

  # If not found, fall back to pattern (starts with xo)
  if (is.na(size_found)) {
    size_found <- nms[stringr::str_detect(nms, "^xo($|_)")][1]
  }

  if (!is.na(size_found) && size_found != "xo_mm") {
    df <- dplyr::rename(df, xo_mm = dplyr::all_of(size_found))
  }

  # --- Q3 percent column (Q₃ / %) ---
  # janitor::clean_names() converts "Q₃ / %" to something like "q_3_percent" or "q3_percent"
  # Expanded candidates and patterns
  q3_candidates <- c("q3_percent", "q_3_percent", "q3_pct", "q_3_pct", "q_3")
  q3_found <- q3_candidates[q3_candidates %in% nms][1]

  # Pattern fallback: more flexible - starts with q and optionally contains 3, percent, or pct
  if (is.na(q3_found)) {
    # First try: q followed by optional underscore/3 and contains "percent"
    q3_found <- nms[
      stringr::str_detect(nms, "^q(_?3)?(_|$)") &
        stringr::str_detect(nms, "percent|pct")
    ][1]
  }

  # Second fallback: just starts with q and has 3 somewhere
  if (is.na(q3_found)) {
    q3_found <- nms[
      stringr::str_detect(nms, "^q") &
        stringr::str_detect(nms, "3")
    ][1]
  }

  if (is.na(q3_found)) {
    if (verbose) {
      cat("DEBUG: Could not find Q3 column\n")
      cat("       Available columns: ", paste(nms, collapse = ", "), "\n", sep = "")
    }
  } else {
    if (verbose) cat("DEBUG: Found Q3 column: ", q3_found, "\n", sep = "")
    if (q3_found != "q3_percent") {
      df <- dplyr::rename(df, q3_percent = dplyr::all_of(q3_found))
    }
  }

  df
}

# ==============================================================================
# INTERNAL: Read PAQXOS distribution data block for each file (auto-skip)
#   FIXED VERSION - includes diagnostic output
# ==============================================================================
.read_paqxos_data_block <- function(files, default_skip = 2, verbose = FALSE) {

  bad_files <- character(0)
  diagnostic_info <- list()

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
    )

    # Store original column names for diagnostics
    original_names <- names(df)

    # Clean names
    df <- df |>
      janitor::clean_names(replace = c(
        "µ" = "u", "μ" = "u", "\u00b5" = "u",
        "₃" = "3", "³" = "3"
      ))

    cleaned_names <- names(df)

    # Standardize
    df <- .standardize_ld_columns(df, verbose = verbose)

    # Basic per-file contract check
    required <- c("xo_mm", "q3_percent")
    if (!all(required %in% names(df))) {
      bad_files <<- c(bad_files, file)
      diagnostic_info[[file]] <<- list(
        original = original_names,
        cleaned = cleaned_names,
        final = names(df),
        missing = setdiff(required, names(df))
      )
      return(dplyr::tibble())  # skip this file cleanly
    }

    df$source_file <- file
    df
  })

  # Emit detailed diagnostic warning
  if (length(bad_files) > 0) {
    warning(
      "Skipped ", length(bad_files), " CSV(s) that did not contain expected PAQXOS columns ",
      "(xo_mm + q3_percent) after auto-detect.\n\n",
      "DIAGNOSTIC INFO for first file:\n",
      if (length(diagnostic_info) > 0) {
        first_file <- names(diagnostic_info)[1]
        info <- diagnostic_info[[first_file]]
        paste0(
          "  File: ", basename(first_file), "\n",
          "  Original columns: ", paste(info$original[1:min(5, length(info$original))], collapse = ", "), "\n",
          "  After clean_names: ", paste(info$cleaned[1:min(5, length(info$cleaned))], collapse = ", "), "\n",
          "  After standardize: ", paste(info$final[1:min(5, length(info$final))], collapse = ", "), "\n",
          "  Missing: ", paste(info$missing, collapse = ", "), "\n"
        )
      } else "",
      "\nFirst few file paths:\n  - ",
      paste(utils::head(basename(bad_files), 5), collapse = "\n  - "),
      if (length(bad_files) > 5) "\n  ... (and more)" else "",
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
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) cat("Created output directory: ", output_dir, "\n", sep = "")
  }

  if (verbose) {
    cat("\n========================================================================\n")
    cat("READING LASER DIFFRACTION DATA\n")
    cat("========================================================================\n")
    cat("Data directory: ", data_directory, "\n", sep = "")
  }

  # --- File discovery (exclude tidy/ outputs to avoid re-ingestion) ---
  all_csv_files <- list.files(
    data_directory,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  tidy_pattern <- paste0(
    "[\\/\\\\]", basename(output_dir), "[\\/\\\\]"
  )

  csv_files <- all_csv_files[
    !stringr::str_detect(all_csv_files, tidy_pattern)
  ]

  if (length(csv_files) == 0) {
    stop("No CSV files found in ", data_directory, call. = FALSE)
  }

  if (verbose) {
    cat("Total CSV files found: ", length(csv_files), "\n", sep = "")
    if (length(all_csv_files) > length(csv_files)) {
      cat("Excluded outputs in /tidy/: ",
          length(all_csv_files) - length(csv_files), "\n", sep = "")
    }
  }

  # --- Read metadata block (rows 1–2) ---
  if (verbose) cat("\nReading PAQXOS metadata (rows 1-2)...\n")
  metadata <- .read_paqxos_metadata(csv_files)

  if (verbose) {
    cat("Files with metadata: ", nrow(metadata), "\n", sep = "")
  }

  # --- Read distribution data block (auto-detect header row) ---
  if (verbose) cat("\nReading distribution data blocks (auto-detecting skip rows)...\n")
  data_block <- .read_paqxos_data_block(
    csv_files,
    default_skip = skip_rows,
    verbose = verbose
  )

  if (nrow(data_block) == 0) {
    stop(
      "No valid PAQXOS distribution data found. Check CSV structure and column names.",
      call. = FALSE
    )
  }

  if (verbose) {
    cat("Files with valid data blocks: ",
        dplyr::n_distinct(data_block$source_file), "\n", sep = "")
  }

  # --- Join metadata + data ---
  combined_data <- dplyr::inner_join(
    data_block,
    metadata,
    by = "source_file",
    suffix = c("", "_meta")
  )

  if (nrow(combined_data) == 0) {
    stop(
      "No data after joining metadata and distribution blocks. ",
      "Check that source_file paths match exactly.",
      call. = FALSE
    )
  }

  if (verbose) {
    cat("Rows after joining metadata + data: ", nrow(combined_data), "\n", sep = "")
  }

  # --- Numeric conversion + CDF ---
  combined_data <- combined_data |>
    dplyr::mutate(
      particle_size_um = suppressWarnings(as.numeric(xo_mm)),
      q3_percent       = suppressWarnings(as.numeric(q3_percent)),
      q3_cdf           = q3_percent / 100
    ) |>
    dplyr::filter(
      !is.na(particle_size_um),
      !is.na(q3_percent),
      particle_size_um > 0
    )

  # --- Module classification (disperser type) ---
  combined_data <- combined_data |>
    dplyr::mutate(
      dispersing_system_clean = stringr::str_trim(
        tolower(
          dplyr::coalesce(dispersing_system, "")
        )
      ),
      is_rodos   = stringr::str_detect(dispersing_system_clean, "rodos"),
      is_inhaler = stringr::str_detect(dispersing_system_clean, "inhaler"),

      # Folder-based fallback
      folder = dirname(source_file),
      is_rodos = dplyr::case_when(
        is_rodos ~ TRUE,
        stringr::str_detect(tolower(folder), "rodos") ~ TRUE,
        TRUE ~ is_rodos
      ),
      is_inhaler = dplyr::case_when(
        is_inhaler ~ TRUE,
        stringr::str_detect(tolower(folder), "inhaler") ~ TRUE,
        TRUE ~ is_inhaler
      ),

      # Require exactly one disperser type
      .is_ambiguous = is_rodos & is_inhaler,
      .is_unknown   = !(is_rodos | is_inhaler),

      module = dplyr::case_when(
        .is_ambiguous ~ NA_character_,
        .is_unknown   ~ NA_character_,
        is_rodos      ~ "RODOS",
        is_inhaler    ~ "INHALER",
        TRUE          ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(module))

  if (verbose) {
    cat("Files classified as RODOS or INHALER: ",
        dplyr::n_distinct(combined_data$source_file), "\n", sep = "")
  }

  # --- Extract analysis-relevant metadata ---
  combined_data <- combined_data |>
    dplyr::mutate(
      .form_id = dplyr::coalesce(formulation_id, NA_character_),
      formulation = dplyr::case_when(
        !is.na(.form_id) ~ .form_id,
        TRUE ~ stringr::str_extract(
          basename(dirname(source_file)),
          formulation_pattern
        )
      ),

      device_resistance = dplyr::case_when(
        is_rodos ~ "reference",
        TRUE     ~ dplyr::case_when(
          stringr::str_detect(tolower(device), "low")    ~ "low",
          stringr::str_detect(tolower(device), "medium") ~ "medium",
          stringr::str_detect(tolower(device), "high")   ~ "high",
          TRUE ~ NA_character_
        )
      ),

      pressure_drop_clean = dplyr::case_when(
        is_rodos ~ NA_real_,
        TRUE     ~ suppressWarnings(
          as.numeric(
            stringr::str_extract(pressure_drop, "\\d+(\\.\\d+)?")
          )
        )
      ),

      measurement_time = dplyr::coalesce(
        time,
        stringr::str_extract(identifier, "\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}")
      ),

      filename = basename(source_file),

      replicate = NA_character_
    ) |>
    dplyr::select(-.form_id)

  # Deterministic  replicate labels
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
      replicate = paste0("rep", .file_rank)  # Apply to all modules
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
