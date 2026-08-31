"""
Universal parser for experimental and computational adsorption isotherm CSV/TSV/TXT tables.
"""

from typing import Dict, Any, Tuple
import os
import pandas as pd
import numpy as np


def parse_isotherm_csv(filepath: str) -> Dict[str, Any]:
    """
    Parses a CSV/TSV file containing adsorption isotherm data.

    Parameters
    ----------
    filepath : str
        Path to tabular file.

    Returns
    -------
    data : dict
        Pressure array, loading array, detected units, and column mappings.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    sep = ","
    if filepath.endswith(".tsv") or filepath.endswith(".tab"):
        sep = "\t"
    elif filepath.endswith(".txt"):
        with open(filepath, "r", encoding="utf-8") as f:
            first_line = f.readline()
            if "\t" in first_line:
                sep = "\t"
            elif ";" in first_line:
                sep = ";"
            elif "," in first_line:
                sep = ","
            else:
                sep = r"\s+"

    df = pd.read_csv(filepath, sep=sep, engine="python")
    
    # Identify pressure column
    p_col = None
    q_col = None
    
    for c in df.columns:
        c_l = str(c).lower()
        if any(term in c_l for term in ["press", "p (", "p_", "p[", "pressure"]):
            p_col = c
        elif any(term in c_l for term in ["load", "q (", "q_", "q[", "uptake", "amount", "capacity"]):
            q_col = c

    # Fallback to first two numeric columns
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if p_col is None and len(numeric_cols) >= 1:
        p_col = numeric_cols[0]
    if q_col is None and len(numeric_cols) >= 2:
        q_col = numeric_cols[1]

    if p_col is None or q_col is None:
        raise ValueError(f"Could not automatically identify pressure and loading columns in {filepath}.")

    pressures = df[p_col].dropna().values.astype(float)
    loadings = df[q_col].dropna().values.astype(float)
    
    # Ensure equal length
    min_len = min(len(pressures), len(loadings))
    pressures = pressures[:min_len]
    loadings = loadings[:min_len]

    return {
        "pressure": pressures,
        "loading": loadings,
        "pressure_col": str(p_col),
        "loading_col": str(q_col),
        "n_points": min_len
    }
