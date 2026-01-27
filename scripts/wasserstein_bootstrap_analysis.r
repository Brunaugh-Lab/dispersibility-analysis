# ==============================================================================
# 04_bootstrap_analysis.R
# Bootstrap Resampling for Wasserstein Distance Confidence Intervals
#
# Purpose: Generate empirical sampling distributions for W1 distances through
#          bootstrap resampling of replicate reference and inhaler distributions with replacement
#          Pool or average across resampled triplicate Q3 values at each particle size to get
#          RODOS + INHALER empirical distributions.
#          For each RODOS + INHALER distribution pair, the W1 distance is calculated.
#          Calculate confidence intervals and effect-to-noise ratios to
#          quantify measurement uncertainty vs true condition effects.
#
# Auto-execution: Script automatically runs when sourced
#   - Reads data/tidy/standardized_data.csv (from script 01)
#   - Performs bootstrap resampling (default: 2000 iterations)
#   - Calculates 95% CIs and standard errors for each formulation
#   - Computes effect-to-noise ratios for device resistance/pressure effects
#   - Saves to results/bootstrap_results.csv and results/effect_noise_ratios.csv
#
# Input: data/tidy/standardized_data.csv (from 01_data_import.R)
# Output:
#   - results/bootstrap_results.csv (W1 distributions with CIs)
#   - results/effect_noise_ratios.csv (signal vs noise quantification)
#
# Methodology:
#   1. Pool replicates → create empirical particle size distributions
#   2. Bootstrap resample with replacement from pooled distributions
#   3. Calculate W1 for each bootstrap iteration
#   4. Estimate sampling distribution parameters (mean, SD, 95% CI)
#   5. Compare bootstrap variability to between-condition effects
#
# Usage:
#   source("scripts/wassserstein_bootstrap_analysis.R")
#
# ==============================================================================
