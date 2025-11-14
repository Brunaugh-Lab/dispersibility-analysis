from pathlib import Path
import numpy as np
import pandas as pd

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