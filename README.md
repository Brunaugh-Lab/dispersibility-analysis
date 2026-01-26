# Dispersibility Analysis Toolkit

A reproducible R-based toolkit for quantifying dry powder inhaler dispersibility using the Wasserstein distance metric.

## 📋 What This Does

This toolkit calculates **dispersibility metrics** for dry powder inhalers by:
1. Comparing particle size distributions from inhaler-based dispersion to maximally dispersed reference conditions (RODOS)
2. Computing the Wasserstein-1 (W₁) distance - a measure of how much particle redistribution is needed to achieve full dispersion
3. Normalizing W₁ by reference particle size (d₅₀) to enable fair cross-formulation comparisons

**Physical interpretation:** Lower W₁ = better dispersibility (closer to fully dispersed state)

---

## 🗂️ Required Folder Structure

Organize your laser diffraction data like this:

```
Wasserstein_DPI/
├── data/
│   ├── FormulationA/          # Each formulation gets its own folder
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
│   └── ...
└── scripts/
    ├── 01_data_import.R
    ├── 02_wasserstein_core.R  (coming soon)
    └── ...
```

**Important notes:**
- Formulation folder names will be used as formulation IDs
- Subfolder names must contain "inhaler" or "rodos" (case-insensitive)
- Replicate filenames must contain "rep1", "rep2", "rep3" (or "Rep1", "Rep_1", etc.)
- CSV files must be Sympatec PAQXOS exports (standard format with 2 header rows)

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

**Basic usage** (if your formulation folders are simple names):

```r
data <- read_ld_data_from_structure(
  data_directory = "data/",
  verbose = TRUE
)
```

**For formulations named like "132067_IMT", "231067_IMT", etc:**

```r
data <- read_ld_data_from_structure(
  data_directory = "data/",
  formulation_pattern = "\\d+_IMT",      # Extracts the numeric_IMT part
  replicate_pattern = "[Rr]ep_?\\d+",    # Matches rep1, Rep1, rep_1, Rep_1
  verbose = TRUE
)

# Standardize replicate names to lowercase for consistency
data <- data %>%
  mutate(replicate = tolower(replicate))
```

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
  formulation INHALER RODOS
1 132067_IMT        3     3
2 231067_IMT        3     3
3 262054_IMT        3     3
...
```

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

---

## 📊 What You Get

The imported `data` object is a tibble (data frame) with these columns:

| Column | Description | Example Values |
|--------|-------------|----------------|
| `particle_size_um` | Particle diameter in micrometers | 0.5, 1.0, 2.0, ... 100 |
| `q3_percent` | Cumulative volume distribution (0-100%) | 0, 10.5, 45.2, ... 100 |
| `q3_cdf` | Cumulative distribution function (0-1) | 0, 0.105, 0.452, ... 1.0 |
| `formulation` | Formulation identifier | "132067_IMT", "231067_IMT" |
| `module` | Dispersion module (standardized) | "INHALER", "RODOS" |
| `replicate` | Replicate identifier | "rep1", "rep2", "rep3" |
| `source_file` | Full path to original CSV | "data/132067_IMT/inhaler/rep1.csv" |

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

## 🎯 Next Steps

Once your data is successfully imported and validated:

1. **Calculate Wasserstein distances** (script coming soon: `02_wasserstein_core.R`)
2. **Compute dispersibility metrics** (script coming soon: `03_dispersibility_metrics.R`)
3. **Generate visualizations** (script coming soon: `04_visualization.R`)
4. **Run statistical analysis** (script coming soon: `05_statistical_analysis.R`)

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
