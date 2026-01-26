# Dispersibility Analysis Toolkit

A reproducible R-based toolkit for quantifying dry powder inhaler dispersibility using the Wasserstein distance metric.

## 📋 What This Does

This toolkit calculates **dispersibility metrics** for dry powder inhalers by:
1. Comparing particle size distributions from inhaler-based dispersion to maximally dispersed reference conditions (RODOS)
2. Computing the Wasserstein-1 (W₁) distance - a measure of how much particle redistribution is needed to achieve full dispersion
3. Normalizing W₁ by reference particle size (d₅₀) to enable fair cross-formulation comparisons

**Physical interpretation:** Lower W₁ = better dispersibility (closer to fully dispersed state)

**Pipeline workflow:**
- Script 01: Reads raw data → Saves `data/tidy/standardized_data.csv`
- Script 02: Reads tidy data → Saves `results/wasserstein_results.csv`
- Script 03: Reads both → Generates `figures/*.png`

Each script is **standalone** and auto-creates needed folders/files.

**Complete pipeline (3 lines):**
```r
source("scripts/01_data_import.R")       # Auto-imports → data/tidy/
source("scripts/02_wasserstein_core.R")  # Auto-calculates → results/
source("scripts/03_visualization.R")     # Auto-plots → figures/
```

---

## 🗂️ Required Folder Structure

Organize your laser diffraction data like this:

```
Wasserstein_DPI/
├── data/                      # Raw CSV files and tidy data
│   ├── FormulationA/          # Any formulation naming scheme works!
│   │   ├── inhaler/            # Inhaler dispersion data (lowercase or UPPERCASE)
│   │   │   ├── rep1.csv
│   │   │   ├── rep2.csv
│   │   │   └── rep3.csv
│   │   └── rodos/              # RODOS reference data (lowercase or UPPERCASE)
│   │       ├── rep1.csv
│   │       ├── rep2.csv
│   │       └── rep3.csv
│   ├── FormulationB/
│   │   ├── inhaler/
│   │   └── rodos/
│   ├── ...
│   └── tidy/                   # Auto-created by 01_data_import.R
│       └── standardized_data.csv  # Cleaned data for analysis
├── results/                   # Auto-created by 02_wasserstein_core.R
│   └── wasserstein_results.csv
├── figures/                   # Auto-created by 03_visualization.R
│   ├── cdf_comparison.png
│   ├── w1_ranking.png
│   ├── d50_comparison.png
│   └── dispersibility_panel.png
└── scripts/
    ├── 01_data_import.R
    ├── 02_wasserstein_core.R
    └── 03_visualization.R
```

**Formulation naming is flexible!** Examples of valid folder names:
- Simple: `FormA`, `FormB`, `FormC`
- Numbered: `Run2`, `Run3`, `Run4`
- Descriptive: `Trehalose_High`, `Mannitol_Low`
- Coded: `132067_IMT`, `231067_IMT` (current project example)
- Complex: `F2_40-20-40_IMT`, `F3_33-00-67_IMT`

By default, the entire folder name becomes the formulation ID.

**Important notes:**
- **Formulation folder names** = your formulation IDs (by default, entire name is used)
- **Subfolder names** must contain "inhaler" or "rodos" (case-insensitive)
- **Replicate filenames** must contain "rep1", "rep2", "rep3" (or "Rep1", "Rep_1", etc.)
- **CSV files** must be Sympatec PAQXOS exports (standard format with 2 header rows)
- **Flexibility:** You can extract only part of folder names using custom patterns (see Step 4)

---

## 🚀 Quick Start Guide

### Step 1: Install Required R Packages

Open R or RStudio and run:

```r
# Install packages if you don't have them
install.packages("tidyverse")
install.packages("janitor")
```

### Step 2: Set Your Working Directory

```r
# Set working directory to your repository
setwd("~/Documents/GitHub/Wasserstein_DPI")
# Or on Windows: setwd("C:/Users/YourName/Documents/GitHub/Wasserstein_DPI")
```

### Step 3: Load the Data Import Script

```r
# Load the data import functions
source("scripts/01_data_import.R")
```

### Step 4: Import Your Data

**RECOMMENDED - Basic usage (auto-saves to data/tidy/standardized_data.csv):**

```r
data <- read_ld_data_from_structure(
  data_directory = "data/",
  formulation_pattern = ".*",            # Uses entire folder name (default)
  replicate_pattern = "[Rr]ep_?\\d+",    # Matches rep1, Rep1, rep_1, Rep_1
  save_output = TRUE,                    # Default: saves cleaned data
  verbose = TRUE
)
```

**What this does automatically:**
- ✓ Creates `data/tidy/` folder if it doesn't exist
- ✓ Saves cleaned data to `data/tidy/standardized_data.csv`
- ✓ Standardizes replicate names to lowercase (rep1, rep2, rep3)
- ✓ Validates data structure

This works for ANY folder naming scheme - the folder names become your formulation IDs.

**Loading previously processed data (much faster!):**

If you've already run the import once, you can quickly reload:

```r
# Fast reload without re-reading 40+ raw CSV files
data <- load_standardized_data()
```

**Example project-specific patterns:**

If you need to extract only PART of the folder name, customize `formulation_pattern`:

```r
# Example 1: Folders named "132067_IMT", "231067_IMT" → Extract "132067_IMT"
formulation_pattern = "\\d+_IMT"

# Example 2: Folders named "FormA_batch1_data" → Extract "FormA"
formulation_pattern = "Form[A-Z]"

# Example 3: Folders named "Run2", "Run3", "Run4" → Extract "Run2", "Run3", "Run4"
formulation_pattern = "Run\\d+"

# Example 4: Folders named "F2_40-20-40_IMT" → Extract entire name
formulation_pattern = ".*"  # (default - use this!)
```

**When to customize patterns:**
- Your folder names have extra text you don't want in formulation IDs
- You want to extract specific portions (like numeric codes)
- You need to match a specific naming convention from another lab

**When to use default (`".*"`):**
- Your folder names ARE your formulation IDs (most common!)
- You want to keep the entire folder name
- You're not sure - start with this!

### Step 5: Validate Your Data

```r
# Check that everything imported correctly
validate_ld_data(data, check_replicates = TRUE, min_replicates = 3)

# Quick summary: Should show 3 replicates per formulation-module combination
data %>%
  distinct(formulation, module, replicate) %>%
  count(formulation, module) %>%
  pivot_wider(names_from = module, values_from = n, values_fill = 0)
```

**Expected output:**
```
  formulation  INHALER RODOS
1 FormA              3     3
2 FormB              3     3
3 FormC              3     3
...
```

(Formulation names will match your folder names)

### Step 6: Inspect Your Data

```r
# View first few rows
head(data)

# Open in RStudio viewer
View(data)

# Check how many data points per file
data %>%
  count(formulation, module, replicate) %>%
  print(n = Inf)
```

### Step 7: Calculate Wasserstein Distances

```r
# Load the Wasserstein calculation functions
source("scripts/02_wasserstein_core.R")

# Calculate W1 distances (compares INHALER to RODOS for each formulation)
w1_results <- calculate_pairwise_wasserstein(
  data,
  reference_module = "RODOS",
  test_module = "INHALER",
  verbose = TRUE
)

# Validate results (checks for errors and unreasonable values)
validate_wasserstein_results(w1_results)

# View results
print(w1_results, n = Inf)

# Save results to CSV
write_csv(w1_results, "wasserstein_results.csv")
```

**Expected output:**
```
========================================================================
CALCULATING WASSERSTEIN-1 DISTANCES
========================================================================
Reference condition: RODOS
Test condition: INHALER
Methodology: Pool replicates → Calculate W1
------------------------------------------------------------------------

Formulations to process: 7
Processing: FormA ... W1 = 0.3245 µm, W1/d50 = 0.0891
Processing: FormB ... W1 = 0.4123 µm, W1/d50 = 0.1156
...

========================================================================
WASSERSTEIN CALCULATION COMPLETE
========================================================================
Successfully calculated W1 for 7 formulations
...
```

**Results table columns:**
- `W1_micrometers` - **Use this for DoE analysis** (absolute dispersibility)
- `d50_reference_um` - Reference median diameter
- `d50_test_um` - Test median diameter
- `W1_normalized` - W1/d50 ratio (for cross-formulation comparison)
- `d50_shift_um` - Difference in median diameters

---

## 📊 What You Get

### After Data Import (Step 4)

**In-memory R object** - The imported `data` tibble with these columns:

| Column | Description | Example Values |
|--------|-------------|----------------|
| `particle_size_um` | Particle diameter in micrometers | 0.5, 1.0, 2.0, ... 100 |
| `q3_percent` | Cumulative volume distribution (0-100%) | 0, 10.5, 45.2, ... 100 |
| `q3_cdf` | Cumulative distribution function (0-1) | 0, 0.105, 0.452, ... 1.0 |
| `formulation` | Formulation identifier (from folder name) | "FormA", "Run2", "132067_IMT" |
| `module` | Dispersion module (standardized) | "INHALER", "RODOS" |
| `replicate` | Replicate identifier (auto-standardized) | "rep1", "rep2", "rep3" |
| `source_file` | Full path to original CSV | "data/FormA/inhaler/rep1.csv" |

**Saved file** - `data/tidy/standardized_data.csv`
- Cleaned and validated data ready for analysis
- Can be quickly reloaded with `load_standardized_data()`
- Used by downstream scripts (02, 03)

### After Wasserstein Calculation (Step 7)

The `w1_results` object is a tibble with dispersibility metrics:

| Column | Description | Use For |
|--------|-------------|---------|
| `formulation` | Formulation identifier | Grouping |
| `W1_micrometers` | **Absolute W1 distance in µm** | **DoE analysis** |
| `d50_reference_um` | Reference median diameter (µm) | Context |
| `d50_test_um` | Test median diameter (µm) | Context |
| `W1_normalized` | W1/d50 ratio (dimensionless) | Cross-study comparison |
| `d50_shift_um` | Test d50 - Reference d50 (µm) | Understanding shift |

**For Design of Experiments (DoE):** Use `W1_micrometers` as your response variable. Lower values indicate better dispersibility (less redistribution needed to match fully dispersed state).

---

## 🔧 Troubleshooting

### Problem: "Some files have NA formulation"

**Cause:** The `formulation_pattern` doesn't match your folder names.

**Solution:** Try using the entire folder name:
```r
data <- read_ld_data_from_structure(
  data_directory = "data/",
  formulation_pattern = ".*",  # Use entire folder name
  verbose = TRUE
)
```

### Problem: "Some files have NA module"

**Cause:** Subfolders are not named "inhaler" or "rodos".

**Solution:** Check your folder structure:
```r
# List all subdirectories to see what they're called
list.dirs("data/", recursive = TRUE)
```

Rename folders to contain "inhaler" or "rodos" (case doesn't matter).

### Problem: "Some files have NA replicate"

**Cause:** Filenames don't contain "rep1", "rep2", "rep3".

**Solution:**

**Option A - Adjust the pattern** to match your naming:
```r
# For uppercase Rep_1 format
replicate_pattern = "[Rr]ep_?\\d+"

# For "replicate1" format
replicate_pattern = "replicate\\d+"

# For "r1" format
replicate_pattern = "r\\d+"
```

**Option B - Find which files are problematic:**
```r
data %>%
  filter(is.na(replicate)) %>%
  distinct(source_file)
```

Then either rename those files or adjust the pattern.

### Problem: "Some conditions have fewer than 3 replicates"

**Cause:** Missing CSV files for some formulations.

**Solution:** Check which formulation-module combinations are incomplete:
```r
data %>%
  distinct(formulation, module, replicate) %>%
  count(formulation, module) %>%
  filter(n < 3)
```

Find the missing files and add them to the correct folders.

### Problem: "No CSV files found"

**Cause:** Wrong directory path or files are in wrong location.

**Solution:**
```r
# Check your current working directory
getwd()

# List files to verify structure
list.files("data/", recursive = TRUE, pattern = "\\.csv$")
```

Make sure you're in the repository root and CSV files are in `data/` subdirectories.

---

## 🔧 Troubleshooting Wasserstein Calculations

### Problem: "Negative W1 values detected"

**Cause:** This is mathematically impossible and indicates a calculation error.

**Solution:** Check your data:
```r
# Verify CDFs are properly formed (0 to 1, monotonic increasing)
data %>%
  group_by(source_file) %>%
  summarise(
    min_cdf = min(q3_cdf),
    max_cdf = max(q3_cdf),
    is_monotonic = all(diff(q3_cdf) >= 0)
  ) %>%
  filter(min_cdf < 0 | max_cdf > 1 | !is_monotonic)
```

### Problem: "W1/d50 values exceed 2.0"

**Cause:** Very poor dispersibility or potential data quality issues.

**Solution:** This is a warning, not necessarily an error. Large W1/d50 values can be legitimate for highly cohesive powders, but you should:
1. Visually inspect the CDFs for the flagged formulation
2. Verify the raw CSV files are correct
3. Check if RODOS reference shows proper dispersion

### Problem: "Some formulations show negative d50 shifts"

**Cause:** Test condition (INHALER) produced finer aerosol than reference (RODOS).

**Solution:** This is unusual but possible. It might indicate:
- Inhaler is more efficient than expected (good!)
- RODOS didn't fully disperse the powder (check pressure)
- Data quality issue (verify raw files)

### Problem: W1 calculation fails for specific formulation

**Cause:** Usually missing data or file read errors.

**Solution:**
```r
# Check which formulation failed
formulations <- unique(data$formulation)

# Verify data exists for that formulation
data %>%
  filter(formulation == "ProblemFormulation") %>%
  count(module, replicate)

# Check if pooling worked
pool_replicate_cdfs(data, "ProblemFormulation", "RODOS")
pool_replicate_cdfs(data, "ProblemFormulation", "INHALER")
```

---

## 📋 Complete Workflow Example

Here's the full pipeline from raw data to W1 results:

```r
# ============================================================================
# COMPLETE DISPERSIBILITY ANALYSIS WORKFLOW
# ============================================================================

# Set working directory
setwd("~/Documents/GitHub/Wasserstein_DPI")

# Install packages (only needed once)
# install.packages("tidyverse")
# install.packages("janitor")

library(tidyverse)

# ----------------------------------------------------------------------------
# STEP 1: IMPORT DATA
# ----------------------------------------------------------------------------
source("scripts/01_data_import.R")

# Option A: Full import from raw data (first time or when data changes)
data <- read_ld_data_from_structure(
  data_directory = "data/",
  formulation_pattern = ".*",           # Use entire folder name
  replicate_pattern = "[Rr]ep_?\\d+",   # Flexible replicate matching
  save_output = TRUE,                   # Auto-saves to data/tidy/
  verbose = TRUE
)

# Option B: Load previously processed data (subsequent runs - much faster!)
# data <- load_standardized_data()

# Validate data structure
validate_ld_data(data, check_replicates = TRUE, min_replicates = 3)

# Quick check
data %>%
  distinct(formulation, module, replicate) %>%
  count(formulation, module) %>%
  pivot_wider(names_from = module, values_from = n, values_fill = 0)

# ----------------------------------------------------------------------------
# STEP 2: CALCULATE WASSERSTEIN DISTANCES
# ----------------------------------------------------------------------------
source("scripts/02_wasserstein_core.R")

w1_results <- calculate_pairwise_wasserstein(
  data,
  reference_module = "RODOS",
  test_module = "INHALER",
  verbose = TRUE
)

# Validate W1 results
validate_wasserstein_results(w1_results)

# View results
print(w1_results, n = Inf)

# ----------------------------------------------------------------------------
# STEP 3: SAVE RESULTS
# ----------------------------------------------------------------------------
# Save W1 results for DoE analysis
write_csv(w1_results, "wasserstein_results.csv")

# Save a summary for quick reference
w1_summary <- w1_results %>%
  select(formulation, W1_micrometers, W1_normalized) %>%
  arrange(W1_micrometers)  # Sort by dispersibility (best first)

write_csv(w1_summary, "dispersibility_ranking.csv")

cat("\n========================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================================================\n")
cat("Files saved:\n")
cat("  - wasserstein_results.csv (full results)\n")
cat("  - dispersibility_ranking.csv (summary)\n")
cat("\nNext steps:\n")
cat("  1. Use W1_micrometers for DoE modeling\n")
cat("  2. Generate visualizations (scripts coming soon)\n")
cat("  3. Run statistical analysis (scripts coming soon)\n")
cat("========================================================================\n")
```

---

## 🎯 Next Steps

Once your data is successfully imported and W1 distances calculated:

1. ✅ **Import data** (`01_data_import.R`) - Complete!
2. ✅ **Calculate Wasserstein distances** (`02_wasserstein_core.R`) - Complete!
3. ✅ **Generate visualizations** (`03_visualization.R`) - Complete!
   - CDF comparison plots
   - W1 bar charts
   - d50 comparisons
   - Publication-ready panels
4. **Run statistical analysis** (script coming soon)
   - Mixture model fitting
   - Component effects
   - Composition-response relationships

---

## 📚 Citation

If you use this toolkit in your research, please cite:

**Paper in preparation:**
Xia G, Dechayont B, Che L, Comfort I, Brunaugh AD. "A Distribution-Based Metric for Quantifying Dispersibility in Dry Powder Inhalers." *Pharmaceutics* (submitted 2025).

**Methodology reference:**
Xia G, Bennett N, Watts A, Brunaugh AD. "Mapping a Ternary Carbohydrate Design Space for Stable and Dispersible Protein Dry Powders." *Molecular Pharmaceutics* (submitted 2025).

---

## 🤝 Contributing

Found a bug or have a suggestion? Please open an issue on GitHub!

---

## 📝 License

[Add your license here - MIT, GPL, etc.]

---

## 👥 Authors

- **Grace Xia** - Data analysis and method development
- **Ashlee D. Brunaugh** - Principal investigator and toolkit design

University of Michigan, College of Pharmacy, Department of Pharmaceutical Sciences

---

## 📧 Contact

Questions? Email: brunaugh@umich.edu
