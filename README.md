# Wasserstein Dispersibility Toolkit

A reproducible R-based toolkit for quantifying dry powder inhaler dispersibility using distribution-level metrics based on the Wasserstein distance.

---

## 01. Data Import & Standardization

Utilities for importing and standardizing Sympatec PAQXOS laser diffraction data.

`01_data_import.R` reads raw RODOS (reference) and INHALER (test) CSV exports, extracts and standardizes metadata, converts cumulative volume distributions to proper CDFs, and saves a single tidy dataset for downstream analysis.

RODOS measurements are treated as formulation-level reference states; INHALER measurements may vary by device resistance and pressure. Replicate identity is preserved and pooling is performed downstream.

**Usage**
```r
source("scripts/01_data_import.R")
data <- run_data_import("data")
```
Creates `data/tidy/standardized_data_with_conditions.csv`.

---

## 02. Wasserstein Distance Calculation *(coming soon)*

Tools for computing first-order Wasserstein (W₁) distances between inhaler-generated and reference particle size distributions.

This module will:
- pool replicates prior to distance calculation
- support multi-factor experimental designs (formulation × device × pressure)
- report absolute and normalized W₁ metrics suitable for DoE analysis

---

## 03. Visualization *(coming soon)*

Functions for visualizing particle size distributions and dispersibility metrics, including:
- reference vs inhaler CDF comparisons
- formulation-level distribution overlays
- W₁-based ranking plots

---

## 04. Statistical Analysis *(planned)*

Extensions for uncertainty quantification and statistical interrogation of dispersibility metrics, including bootstrap-based analyses and effect-to-noise comparisons.

---

## Scope

This toolkit is intended to support distribution-level reasoning about aerosol dispersibility and to complement, not replace, pharmacopeial aerodynamic endpoints.

---

## Authors

Grace Xia · Ashlee D. Brunaugh
University of Michigan, College of Pharmacy
