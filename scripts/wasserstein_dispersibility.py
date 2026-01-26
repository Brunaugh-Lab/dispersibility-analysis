from pathlib import Path
import re
import numpy as np
import pandas as pd
from collections import defaultdict

def load_sympatec_cdf(csv_path):
    """
    Load a Sympatec CSV and return a dataframe with:
        - size_um: particle size grid (µm)
        - cdf: cumulative fraction (0–1)
    """
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path, skiprows=2)

    # Rename the two columns we care about
    df = df.rename(columns={
        "xo / µm": "size_um",
        "Q₃ / %": "q3_percent",
    })

    # Keep only relevant columns and drop trailing NA row
    df = df[["size_um", "q3_percent"]].dropna(subset=["size_um"])

    # Convert to numeric
    df["size_um"] = pd.to_numeric(df["size_um"], errors="raise")
    df["q3_percent"] = pd.to_numeric(df["q3_percent"], errors="raise")

    # Sort just in case
    df = df.sort_values("size_um").reset_index(drop=True)

    # Convert percent → fraction
    df["cdf"] = df["q3_percent"] / 100.0

    return df[["size_um", "cdf"]]


def wasserstein_1d_from_cdfs(df_inhaler, df_rodos):
    """
    Compute 1-Wasserstein distance between INHALER and RODOS
    CDFs given on (ideally) the same size grid.

    Returns W1 in units of µm.
    """

    x_I = df_inhaler["size_um"].to_numpy()
    x_R = df_rodos["size_um"].to_numpy()

    # If grids are not identical, you *can* interpolate one onto the other.
    # For now, assert they match (which they should if recorded identically).
    if not np.allclose(x_I, x_R, rtol=1e-6, atol=1e-9):
        raise ValueError("Size grids do not match between INHALER and RODOS.")

    x = x_I
    F_I = df_inhaler["cdf"].to_numpy()
    F_R = df_rodos["cdf"].to_numpy()

    # Stepwise integral of |F_I - F_R| over x
    dx = np.diff(x)
    diff_F = np.abs(F_I[:-1] - F_R[:-1])

    W1 = np.sum(diff_F * dx)
    return W1


def find_inhaler_rodos_pairs(base_dir):
    """
    Search recursively under base_dir for CSV files, identify INHALER and RODOS
    measurements, and pair them by (formulation_id, run_id, replicate).

    Assumes a folder layout like:
        base_dir / <formulation_id> / Rep_X / <filename>.csv

    Filenames must contain:
        - 'INHALER' or 'RODOS'
        - a run ID like 'Run4'
        - a replicate like 'rep1', 'Rep1', 'rep_1', etc.

    Returns a list of dicts with keys:
        - formulation_id
        - run_id
        - replicate   (canonical form: 'rep1', 'rep2', ...)
        - inhaler_path
        - rodos_path
    """
    base_dir = Path(base_dir)
    entries = []

    for csv_path in base_dir.rglob("*.csv"):
        name = csv_path.name

        # Identify module type (INHALER vs RODOS)
        module_match = re.search(r"(INHALER|RODOS)", name)
        if not module_match:
            continue
        module = module_match.group(1)

        # Extract run ID, e.g. "Run4"
        run_match = re.search(r"Run\d+", name)
        run_id = run_match.group(0) if run_match else None

        # Extract replicate, e.g. "rep1", "Rep1", "rep_1"
        rep_match = re.search(r"[Rr]ep[_-]?\d+", name)
        if rep_match:
            token = rep_match.group(0)          # e.g. "Rep3" or "rep_3"
            num_match = re.search(r"\d+", token)
            rep_num = num_match.group(0) if num_match else None
            replicate = f"rep{rep_num}" if rep_num is not None else None
        else:
            replicate = None

        # Formulation folder: .../<formulation_id>/Rep_X/file.csv
        # parent = Rep_X, parent.parent = formulation_id
        formulation_id = csv_path.parents[1].name if len(csv_path.parents) >= 2 else None

        entries.append({
            "formulation_id": formulation_id,
            "run_id": run_id,
            "replicate": replicate,  # canonical
            "module": module,
            "path": csv_path,
        })

    # Group by (formulation_id, run_id, replicate) and pair INHALER/RODOS
    grouped = defaultdict(dict)
    for e in entries:
        key = (e["formulation_id"], e["run_id"], e["replicate"])
        grouped[key][e["module"]] = e["path"]

    pairs = []
    for (formulation_id, run_id, replicate), modules in grouped.items():
        inh = modules.get("INHALER")
        rod = modules.get("RODOS")
        if inh and rod:
            pairs.append({
                "formulation_id": formulation_id,
                "run_id": run_id,
                "replicate": replicate,
                "inhaler_path": inh,
                "rodos_path": rod,
            })

    return pairs


def compute_w1_table(base_dir):
    """
    For all INHALER/RODOS CSV pairs under base_dir, compute the 1D 1-Wasserstein
    distance (in µm) and return a pandas DataFrame summarizing the results.

    Columns:
        - formulation_id
        - run_id
        - replicate
        - inhaler_file
        - rodos_file
        - W1_um
    """
    pairs = find_inhaler_rodos_pairs(base_dir)
    rows = []

    for pair in pairs:
        df_inh = load_sympatec_cdf(pair["inhaler_path"])
        df_rod = load_sympatec_cdf(pair["rodos_path"])
        w1 = wasserstein_1d_from_cdfs(df_inh, df_rod)

        rows.append({
            "formulation_id": pair["formulation_id"],
            "run_id": pair["run_id"],
            "replicate": pair["replicate"],
            "inhaler_file": str(pair["inhaler_path"]),
            "rodos_file": str(pair["rodos_path"]),
            "W1_um": w1,
        })

    return pd.DataFrame(rows)

def compute_d50_from_cdf(df):
    """
    Given a dataframe with columns 'size_um' and 'cdf' (0–1),
    compute the volume-based D50 by linear interpolation.

    Assumes 'cdf' is monotone increasing from 0 to ~1.
    Returns D50 in µm.
    """
    x = df["size_um"].to_numpy()
    F = df["cdf"].to_numpy()

    # Find first index where CDF >= 0.5
    idx = np.searchsorted(F, 0.5)

    if idx == 0:
        return x[0]
    if idx >= len(x):
        return x[-1]

    x0, x1 = x[idx - 1], x[idx]
    F0, F1 = F[idx - 1], F[idx]

    # Linear interpolation between (x0, F0) and (x1, F1)
    if F1 == F0:
        return x0
    frac = (0.5 - F0) / (F1 - F0)
    return x0 + frac * (x1 - x0)


def compute_rodos_d50_table(base_dir):
    """
    For all INHALER/RODOS CSV pairs under base_dir, compute the RODOS D50
    (in µm) and return a pandas DataFrame summarizing the results.

    One row per replicate-level condition:
        - formulation_id
        - run_id
        - replicate
        - rodos_d50_um
    """
    pairs = find_inhaler_rodos_pairs(base_dir)
    rows = []

    for pair in pairs:
        df_rod = load_sympatec_cdf(pair["rodos_path"])
        d50 = compute_d50_from_cdf(df_rod)

        rows.append({
            "formulation_id": pair["formulation_id"],
            "run_id": pair["run_id"],
            "replicate": pair["replicate"],
            "rodos_d50_um": d50,
        })

    return pd.DataFrame(rows)


def compute_condition_level_summary(base_dir):
    """
    Compute condition-level (formulation × run) means for W1 and D50,R.

    Returns a DataFrame with columns:
        - formulation_id
        - run_id
        - W1_mean_um
        - W1_sd_um
        - n_reps
        - D50R_mean_um
        - D50R_sd_um
    """
    df_w1 = compute_w1_table(base_dir)
    df_d50 = compute_rodos_d50_table(base_dir)

    # Merge replicate-level W1 and D50 tables
    df = df_w1.merge(
        df_d50,
        on=["formulation_id", "run_id", "replicate"],
        how="inner",
    )

    # Condition-level means
    df_cond = (
        df
        .groupby(["formulation_id", "run_id"], as_index=False)
        .agg(
            W1_mean_um=("W1_um", "mean"),
            W1_sd_um=("W1_um", "std"),
            n_reps=("W1_um", "size"),
            D50R_mean_um=("rodos_d50_um", "mean"),
            D50R_sd_um=("rodos_d50_um", "std"),
        )
    )

    return df_cond


def normalization_diagnostics(df_cond, cv_threshold=0.15, corr_threshold=0.5):
    """
    Given a condition-level summary dataframe (from compute_condition_level_summary),
    compute:
        - CV_D50 across conditions
        - Pearson correlation between W1_mean_um and D50R_mean_um

    Returns a dict:
        {
            "CV_D50": float,
            "corr_W1_D50": float,
            "normalize": bool
        }

    Normalization is recommended if:
        CV_D50 >= cv_threshold AND abs(corr_W1_D50) >= corr_threshold.
    """
    d50 = df_cond["D50R_mean_um"]
    w1 = df_cond["W1_mean_um"]

    cv_d50 = d50.std(ddof=1) / d50.mean()
    corr = w1.corr(d50)

    normalize = (cv_d50 >= cv_threshold) and (abs(corr) >= corr_threshold)

    return {
        "CV_D50": cv_d50,
        "corr_W1_D50": corr,
        "normalize": normalize,
    }


def run_normalization_analysis(base_dir, cv_threshold=0.15, corr_threshold=0.5):
    """
    High-level helper to:
        1) Compute condition-level means for W1 and D50,R
        2) Evaluate whether normalization is recommended
        3) If recommended, add a 'W1_norm_by_D50R' column

    Returns:
        df_cond : DataFrame with condition-level summary
                  (and W1_norm_by_D50R if normalization is recommended)
        diagnostics : dict from normalization_diagnostics(...)
    """
    df_cond = compute_condition_level_summary(base_dir)
    diagnostics = normalization_diagnostics(
        df_cond,
        cv_threshold=cv_threshold,
        corr_threshold=corr_threshold,
    )

    if diagnostics["normalize"]:
        df_cond["W1_norm_by_D50R"] = df_cond["W1_mean_um"] / df_cond["D50R_mean_um"]

    return df_cond, diagnostics