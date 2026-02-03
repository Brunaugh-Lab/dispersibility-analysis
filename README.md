# Dispersibility Analysis Toolkit

A reproducible R-based toolkit for quantifying dry powder inhaler dispersibility using the Wasserstein distance metric.

## What This Does

This toolkit calculates **dispersibility metrics** for dry powder inhalers by:
1. Importing Sympatec PAQXOS laser diffraction exports with automatic header detection
2. Comparing particle size distributions from inhaler-based dispersion (INHALER) to maximally dispersed reference conditions (RODOS)
3. Computing the Wasserstein-1 (W₁) distance - a measure of how much particle redistribution is needed to achieve full dispersion
4. Generating publication-ready visualizations of cumulative distribution functions and dispersibility rankings

**Physical interpretation:** Lower W₁ = better dispersibility (closer to fully dispersed state)

**Pipeline workflow:**
- **Script 01** (`01_data_import.R`): Reads raw PAQXOS CSVs → Saves `data/tidy/standardized_data_with_conditions.csv`
  - Auto-detects data header row per file (handles variations in PAQXOS export format)
  - Assigns replicate numbers deterministically based on measurement timestamps
  - Works with any folder structure - no rigid naming requirements
  - Validates data quality and flags unexpected issues

- **Script 02** (`02_wasserstein_core.R`): Reads tidy data → Saves `results/wasserstein_results.csv`
  - Pools technical replicates before W1 calculation (proper methodology)
  - Handles multiple test conditions (device resistance × pressure drop combinations)
  - Calculates both absolute W₁ (µm) and normalized W₁/d₅₀ metrics

- **Script 03** (`03_visualization.R`): Reads both files → Generates PDFs in `figures/`
  - One PDF per formulation showing reference vs all test conditions
  - Pairwise comparison PDFs for each formulation-condition combination
  - W₁ dispersibility ranking plots
  - Uses same replicate pooling as W₁ calculations for methodological consistency

Each script is **standalone** and auto-creates needed folders/files.

**Complete pipeline (3 commands):**
```r
source("scripts/01_data_import.R")
data <- run_data_import("data")                      # Auto-imports → data/tidy/

source("scripts/02_wasserstein_core.R")
w1_results <- run_wasserstein_analysis()             # Auto-calculates → results/

source("scripts/03_visualization.R")
generate_all_plots(data, w1_results)                 # Auto-plots → figures/
```

---

## Required Folder Structure

**Flexible structure** - The pipeline works with any organization as long as:
1. CSV files are somewhere under `data/`
2. Filenames or folder paths contain "RODOS" or "INHALER" (case-insensitive)
3. Files have `formulation_id` in their PAQXOS metadata (preferred) OR use folder name as fallback

**Recommended structure:**
```
Wasserstein_DPI/
├── data/
│   ├── RODOS/                 # Reference disperser data
│   │   ├── file1.csv          # Timestamp determines replicate order
│   │   ├── file2.csv          # (earliest = rep1, next = rep2, etc.)
│   │   └── file3.csv
│   ├── INHALER/               # Test disperser data
│   │   ├── high_1kPa_file1.csv
│   │   ├── high_1kPa_file2.csv
│   │   ├── high_1kPa_file3.csv
│   │   ├── high_2kPa_file1.csv
│   │   └── ...                # Multiple conditions supported
│   └── tidy/                  # Auto-created by 01_data_import.R
│       └── standardized_data_with_conditions.csv
├── results/                   # Auto-created by 02_wasserstein_core.R
│   └── wasserstein_results.csv
├── figures/                   # Auto-created by 03_visualization.R
│   ├── FormA_reference_vs_high_1kPa.pdf
│   ├── FormA_reference_vs_high_2kPa.pdf
│   └── ...
└── scripts/
    ├── 01_data_import.R
    ├── 02_wasserstein_core.R
    └── 03_visualization.R
```

**Alternative valid structures:**
```
data/
├── Formulation_A/
│   ├── RODOS/
│   │   └── *.csv
│   └── INHALER/
│       └── *.csv
├── Formulation_B/
│   ├── rodos/              # Lowercase works too
│   │   └── *.csv
│   └── inhaler/
│       └── *.csv
```

**Key flexibility features:**
- ✅ **No rigid naming requirements** - any folder/file names work
- ✅ **Replicate assignment is automatic** - based on measurement timestamps in PAQXOS metadata
- ✅ **Case-insensitive** - "RODOS", "rodos", "Rodos" all work
- ✅ **Multiple test conditions** - handles device resistance and pressure drop variations
- ✅ **Formulation ID detection** - uses PAQXOS metadata field or folder name as fallback

**Important notes:**
- **CSV files** must be Sympatec PAQXOS exports (standard format: 2 metadata rows + distribution data)
- **Module detection** happens via:
  1. PAQXOS metadata field "Dispersing system" (primary)
  2. Folder path containing "rodos" or "inhaler" (fallback)
- **Replicate assignment** is deterministic:
  - Files grouped by (formulation, module, device, pressure_drop)
  - Sorted by measurement timestamp
  - Assigned rep1, rep2, rep3... sequentially
  - NO manual replicate labeling in filenames required!

---

## Quick Start Guide

### Step 1: Install Required R Packages

Open R or RStudio and run:

```r
# Core packages for data import and Wasserstein calculation
install.packages(c("tidyverse", "janitor"))

# Additional packages for visualization
install.packages(c("viridis", "patchwork"))
```

### Step 2: Set Your Working Directory

```r
# Set working directory to your repository
setwd("~/Documents/GitHub/Wasserstein_DPI")
# Or on Windows: setwd("C:/Users/YourName/Documents/GitHub/Wasserstein_DPI")
```

### Step 3: Import Your Data

```r
# Load the data import script
source("scripts/01_data_import.R")

# Run the complete import + validation pipeline
data <- run_data_import(
  data_directory = "data",
  verbose = TRUE
)
```

**What happens automatically:**
- ✓ Recursively finds all CSV files under `data/`
- ✓ Auto-detects the distribution data header row for each file
- ✓ Extracts formulation ID from PAQXOS metadata (or uses folder name)
- ✓ Determines module type (RODOS vs INHALER) from metadata/path
- ✓ Assigns replicates deterministically based on timestamps
- ✓ Creates `data/tidy/` folder
- ✓ Saves standardized data to `data/tidy/standardized_data_with_conditions.csv`
- ✓ Validates data structure (checks for NAs, CDF ranges, replicate counts)

**Expected console output:**
```
========================================================================
READING LASER DIFFRACTION DATA
========================================================================
Data directory: data
Total CSV files found: 90
Excluded outputs in /tidy/: 0

Reading PAQXOS metadata (rows 1-2)...
Files with metadata: 90

Reading distribution data blocks (auto-detecting skip rows)...
Files with valid data blocks: 90
Rows after joining metadata + data: 2880
Files classified as RODOS or INHALER: 90

Data extraction summary:
Formulations found: 3
[1] "INU_SULFB_20-1" "MAN_SULFB_20-1" "TRE_SULFB_20-1"

Modules found: 2
[1] "INHALER" "RODOS"

Files per formulation-module combination:
# A tibble: 3 × 3
  formulation    INHALER RODOS
1 INU_SULFB_20-1      27     3
2 MAN_SULFB_20-1      27     3
3 TRE_SULFB_20-1      27     3

========================================================================
VALIDATING DATA STRUCTURE
========================================================================
✓ All required columns present
✓ No unexpected NA values in key columns
✓ CDF values within [0,1] range
✓ All conditions have ≥ 3 replicates
✓ No duplicate file entries
------------------------------------------------------------------------
VALIDATION PASSED: Data structure is ready for analysis
========================================================================
```

**Loading previously processed data (much faster!):**

If you've already run the import once, you can quickly reload:

```r
# Fast reload without re-reading raw CSV files
data <- load_standardized_data(
  processed_dir = "data/tidy",
  verbose = TRUE
)
```

### Step 4: Calculate Wasserstein Distances

```r
# Load the Wasserstein calculation script
source("scripts/02_wasserstein_core.R")

# Run the complete W1 analysis pipeline
w1_results <- run_wasserstein_analysis(
  tidy_data_path = "data/tidy/standardized_data_with_conditions.csv",
  output_dir = "results",
  verbose = TRUE
)
```

**What happens automatically:**
- ✓ Loads standardized data from script 01
- ✓ Groups by formulation and test condition (device + pressure)
- ✓ Pools technical replicates (averages CDFs before W1 calculation)
- ✓ Calculates W1 distance between each test condition and its reference
- ✓ Computes both absolute W₁ (µm) and normalized W₁/d₅₀
- ✓ Creates `results/` folder
- ✓ Saves results to `results/wasserstein_results.csv`
- ✓ Validates results (checks for negative W1, extreme values, etc.)

**Expected console output:**
```
========================================================================
CALCULATING WASSERSTEIN-1 DISTANCES
========================================================================
Reference condition: RODOS
Test conditions: INHALER (multiple device/pressure combinations)
Methodology: Pool replicates → Calculate W1
Output file: results/wasserstein_results.csv
------------------------------------------------------------------------

Formulations to process: 3
Formulation IDs: INU_SULFB_20-1, MAN_SULFB_20-1, TRE_SULFB_20-1

Processing formulation: INU_SULFB_20-1
  Reference: 3 replicates pooled
  Test condition: high_1kPa → W1 = 2.34 µm
  Test condition: high_2kPa → W1 = 1.89 µm
  Test condition: high_4kPa → W1 = 1.45 µm
  ...

========================================================================
WASSERSTEIN CALCULATION COMPLETE
========================================================================
Total comparisons: 81 (3 formulations × 27 test conditions)
Results saved to: results/wasserstein_results.csv
```

**Results table structure:**
```
formulation      device_resistance  pressure_drop  W1_micrometers  W1_normalized  d50_reference_um  d50_test_um
INU_SULFB_20-1   high              1              2.34            0.68           3.45              5.12
INU_SULFB_20-1   high              2              1.89            0.55           3.45              4.68
INU_SULFB_20-1   high              4              1.45            0.42           3.45              4.23
...
```

### Step 5: Generate Visualizations

```r
# Load the visualization script
source("scripts/03_visualization.R")

# Generate all publication-ready plots
generate_all_plots(
  data = data,
  w1_results = w1_results,
  output_dir = "figures"
)
```

**What happens automatically:**
- ✓ Loads standardized data and W1 results
- ✓ Pools replicates (same methodology as W1 calculations)
- ✓ Creates `figures/` folder
- ✓ Generates individual PDFs for each formulation-condition pair
- ✓ Creates overlay PDFs showing reference vs all test conditions
- ✓ Produces W1 ranking plots
- ✓ All plots use publication-ready formatting (high resolution, proper labels, color schemes)

**Generated files:**
```
figures/
├── INU_SULFB_20-1_reference_vs_high_1kPa.pdf
├── INU_SULFB_20-1_reference_vs_high_2kPa.pdf
├── INU_SULFB_20-1_reference_vs_high_4kPa.pdf
├── INU_SULFB_20-1_all_test_conditions.pdf    # Overlay of all conditions
├── w1_ranking_by_condition.pdf               # Dispersibility comparison
└── ...
```

---

## Advanced Usage

### Custom Formulation ID Extraction

If your folder names contain extra information you don't need:

```r
# Example: Folders named "Batch1_FormA_2025" → Extract "FormA"
data <- read_ld_data_from_structure(
  data_directory = "data",
  formulation_pattern = "Form[A-Z]",  # Extracts FormA, FormB, etc.
  verbose = TRUE
)

# Example: Numeric codes only from "132067_IMT_highTemp"
formulation_pattern = "\\d+"  # Extracts 132067

# Default: Use entire folder name (most common)
formulation_pattern = ".*"
```

### Manual Function Calls (More Control)

Instead of `run_data_import()`, you can call functions individually:

```r
# Step 1: Import with custom settings
data <- read_ld_data_from_structure(
  data_directory = "data",
  output_dir = "data/processed",           # Custom output location
  output_filename = "my_data.csv",         # Custom filename
  save_output = TRUE,
  verbose = TRUE
)

# Step 2: Validate separately
validate_ld_data(data, check_replicates = TRUE, min_replicates = 3)

# Step 3: Calculate W1 with custom parameters
w1_results <- calculate_pairwise_wasserstein(
  data = data,
  reference_module = "RODOS",
  test_module = "INHALER",
  output_dir = "custom_results",
  output_filename = "w1_analysis.csv",
  save_output = TRUE,
  verbose = TRUE
)
```

### Subset Analysis

Analyze only specific formulations or conditions:

```r
# Filter to specific formulations
data_subset <- data %>%
  filter(formulation %in% c("FormA", "FormB"))

# Calculate W1 for subset
w1_subset <- calculate_pairwise_wasserstein(
  data = data_subset,
  verbose = TRUE
)

# Visualize subset
generate_all_plots(
  data = data_subset,
  w1_results = w1_subset,
  output_dir = "figures/subset"
)
```

---

## Understanding the Output

### Standardized Data CSV
**Location:** `data/tidy/standardized_data_with_conditions.csv`

**Key columns:**
- `particle_size_um` - Particle diameter (µm)
- `q3_percent` - Cumulative volume % (0-100)
- `q3_cdf` - Cumulative distribution function (0-1)
- `formulation` - Formulation identifier
- `module` - RODOS or INHALER
- `device_resistance` - low/medium/high/reference
- `pressure_drop_clean` - Pressure drop (kPa) - NA for RODOS
- `replicate` - rep1, rep2, rep3... (assigned by timestamp)
- `measurement_time` - Timestamp from PAQXOS
- `source_file` - Original CSV path

**Expected NAs:**
- `pressure_drop_clean` should be NA for all RODOS rows (reference has no pressure drop)
- This is correct behavior and will pass validation

### Wasserstein Results CSV
**Location:** `results/wasserstein_results.csv`

**Key columns:**
- `W1_micrometers` - **Primary metric for DoE analysis** - Absolute dispersibility distance
- `W1_normalized` - W1/d50 ratio for cross-formulation comparison
- `d50_reference_um` - Median diameter of reference (RODOS) distribution
- `d50_test_um` - Median diameter of test (INHALER) distribution
- `d50_shift_um` - Difference (test - reference)
- `device_resistance` - Test device resistance level
- `pressure_drop` - Test pressure drop condition

**Interpretation:**
- **Lower W1 = Better dispersibility** (closer to fully dispersed state)
- W1 ≈ 0 µm would mean perfect dispersion (INHALER = RODOS)
- W1 > 5 µm typically indicates poor dispersibility
- Use `W1_micrometers` for statistical modeling (DoE, regression)
- Use `W1_normalized` for comparing formulations with different particle sizes

---

## Troubleshooting

### Import Issues

**Problem:** "No CSV files found"
```r
# Solution: Check your working directory
getwd()
list.files("data", recursive = TRUE, pattern = "\\.csv$")
```

**Problem:** "Could not auto-detect PAQXOS header row"
- **Cause:** CSV file doesn't have expected PAQXOS structure
- **Solution:** Verify files are exported from Sympatec PAQXOS (not manually edited)

**Problem:** "Some files have NA formulation"
- **Cause:** PAQXOS metadata missing `formulation_id` field AND folder path doesn't match pattern
- **Solution:** Either add `formulation_id` to PAQXOS exports or adjust `formulation_pattern`

### Validation Issues

**Problem:** "VALIDATION FAILED: NA values detected"
- **Check the warning message** - it shows which columns have unexpected NAs
- `pressure_drop_clean` NAs for RODOS rows are **expected and correct**
- Other NAs indicate data quality issues

**Problem:** "Some conditions have fewer than 3 replicates"
```r
# Solution: Check which combinations are incomplete
data %>%
  distinct(source_file, formulation, module, device_resistance, pressure_drop_clean) %>%
  count(formulation, module, device_resistance, pressure_drop_clean) %>%
  filter(n < 3)
```

### Wasserstein Calculation Issues

**Problem:** "Negative W1 values detected"
- **Cause:** Mathematically impossible - indicates data error
- **Solution:** Check CDF values are properly formed (0 to 1, monotonically increasing)

```r
# Diagnostic check
data %>%
  group_by(source_file) %>%
  summarise(
    min_cdf = min(q3_cdf),
    max_cdf = max(q3_cdf),
    is_monotonic = all(diff(q3_cdf) >= 0)
  ) %>%
  filter(min_cdf < 0 | max_cdf > 1 | !is_monotonic)
```

**Problem:** "W1 calculation fails for specific formulation"
```r
# Check data availability
data %>%
  filter(formulation == "ProblematicFormulation") %>%
  count(module, device_resistance, pressure_drop_clean, replicate)

# Check if pooling works
pool_replicate_cdfs(
  data,
  formulation_value = "ProblematicFormulation",
  module_value = "RODOS"
)
```

---

## Statistical Analysis

Once you have W1 results, you can use them for design-of-experiments analysis:

```r
# Load W1 results
w1_results <- read_csv("results/wasserstein_results.csv")

# Use W1_micrometers as response variable in your DoE model
# Example: If you have compositional predictors (x1, x2, x3)
model <- lm(W1_micrometers ~ x1 + x2 + x3 + x1:x2 + x1:x3 + x2:x3,
            data = design_matrix)

# Or for mixture designs:
library(mixexp)
model <- MixtureLM(W1_micrometers ~ x1 + x2 + x3,
                   data = design_matrix)
```

---

## Key Methodological Notes

**Replicate Pooling:**
The toolkit follows best practices by pooling technical replicates before calculating Wasserstein distances. This means:
1. For each condition, replicate CDFs are averaged at each particle size
2. W1 is calculated once between the pooled distributions
3. This approach treats replicates as measurement uncertainty, not biological variation

**Why timestamp-based replicate assignment?**
- Eliminates manual filename labeling errors
- Deterministic and reproducible
- Works across different naming conventions
- Handles multiple test conditions automatically

**Validation philosophy:**
- Expected NAs (RODOS pressure drops) are allowed
- Unexpected NAs trigger warnings but not errors
- Data quality checks happen at every step
- Clear messages help diagnose issues quickly

**Reproducibility and visualization scope:**
This repository provides transparent, reusable implementations of the dispersibility metric and core analysis workflow. Visualization functions prioritize sensible defaults across diverse experimental designs. Manuscript figures may differ in presentation details (faceting, ordering, aesthetics) but use identical underlying calculations.

---

## Computational Performance

**Typical runtime (90 files, 3 formulations):**
- Script 01 (import): ~5-10 seconds
- Script 02 (W1 calculation): ~1-2 seconds
- Script 03 (visualization): ~10-20 seconds
- Total: ~20-30 seconds

**Memory requirements:**
- Minimal - tested with datasets up to 1000 files
- Peak memory usage typically < 500 MB

**Scaling considerations:**
- Runtime scales linearly with number of files
- Visualization time scales with number of conditions
- Consider subset analysis for very large datasets (>1000 files)

---

## Citation

If you use this toolkit in your research, please cite the relevant publications and the software archive as appropriate.

**Software (this toolkit):**  
Brunaugh AD, Xia G. *Dispersibility Analysis Toolkit*. Zenodo.  
DOI: https://doi.org/10.5281/zenodo.18475101

![DOI](https://zenodo.org/badge/1096698618.svg)

**Dispersibility methodology:**  
Xia G, Dechayont B, Che L, Comfort I, Brunaugh AD.  
“A Distribution-Based Metric for Quantifying Dispersibility in Dry Powder Inhalers.”  
*Pharmaceutics* (submitted 2025).

**Application example:**  
Xia G, Bennett N, Watts A, Brunaugh AD.  
“Mapping a Ternary Carbohydrate Design Space for Stable and Dispersible Protein Dry Powders.”  
*Molecular Pharmaceutics* (submitted 2025).

---

## License

MIT License - Free to use for academic and commercial purposes.

---

## Authors

- **Grace Xia** - Method development, data analysis, and software implementation
- **Ashlee D. Brunaugh** - Principal investigator, project design, and scientific oversight

University of Michigan, College of Pharmacy, Department of Pharmaceutical Sciences

---

## Contact

Questions or feedback? Email: brunaugh@umich.edu

---

