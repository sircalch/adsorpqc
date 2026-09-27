"""
Henry coefficient (K_H) and isosteric heat of adsorption (q_st) calculations.
"""

from typing import Dict, Optional
import numpy as np

R_KJ = 8.314462618e-3  # kJ/(mol K)


def calculate_henry_constant(
    pressure: np.ndarray,
    loading: np.ndarray,
    n_points_linear: int = 4,
    max_fraction_of_max_loading: float = 0.15,
) -> Dict[str, float]:
    """
    Henry adsorption constant K_H = lim_{P->0} q/P from the low-pressure part of an isotherm.

    q/P is regressed linearly on q (first-order virial / Langmuir expansion,
    q/P = K_H - b q + ...) using the lowest-pressure points with q below
    max_fraction_of_max_loading of the largest loading (at least 2, at most n_points_linear
    unless more points satisfy the loading criterion); K_H is the intercept. This removes the
    downward bias of a straight line q = K_H P through points that are already curving.

    Parameters
    ----------
    pressure, loading : np.ndarray
        Isotherm (any consistent units; K_H is returned in loading/pressure units).
    n_points_linear : int, default 4
        Minimum number of low-pressure points to use when available.
    max_fraction_of_max_loading : float, default 0.15
        Points with q above this fraction of max(q) are excluded from the extrapolation.

    Returns
    -------
    result : dict
        henry_constant, r_squared of the q/P vs q regression, n_points_used, slope_b.
    """
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    keep = p > 0
    p, q = p[keep], q[keep]
    order = np.argsort(p)
    p, q = p[order], q[order]
    if len(p) < 2:
        return {"henry_constant": float("nan"), "r_squared": float("nan"), "n_points_used": len(p), "slope_b": float("nan")}

    n_low = int(np.sum(q <= max_fraction_of_max_loading * np.max(q)))
    k = n_low if n_low >= 3 else min(len(p), max(2, n_points_linear))
    p_sub, q_sub = p[:k], q[:k]
    y = q_sub / p_sub
    if np.ptp(q_sub) > 0 and k >= 3:
        slope, intercept = np.polyfit(q_sub, y, 1)
        y_pred = intercept + slope * q_sub
    else:
        slope, intercept = 0.0, float(np.mean(y))
        y_pred = np.full_like(y, intercept)
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return {
        "henry_constant": float(intercept),
        "r_squared": float(r_sq),
        "n_points_used": int(k),
        "slope_b": float(-slope),
    }


def calculate_isosteric_heat_fluctuations(
    energy_trajectory: np.ndarray,
    loading_trajectory: np.ndarray,
    temperature_k: float = 298.15,
    energy_unit: str = "kJ/mol",
) -> float:
    """
    Isosteric heat of adsorption from particle-energy fluctuations in the grand-canonical ensemble
    (Vlugt et al., J. Chem. Theory Comput. 4, 1107 (2008)):

        q_st = R T - (<U N> - <U><N>) / (<N^2> - <N>^2)

    Parameters
    ----------
    energy_trajectory : np.ndarray
        Total adsorbate potential energy U (adsorbate-framework + adsorbate-adsorbate) of the
        simulation box, per sample, in `energy_unit` ("kJ/mol" or "K", i.e. U/k_B as RASPA writes it).
    loading_trajectory : np.ndarray
        NUMBER OF MOLECULES N in the same box, per sample (not mol/kg).
    temperature_k : float
    energy_unit : {"kJ/mol", "K"}

    Returns
    -------
    q_st : float
        Isosteric heat in kJ/mol (positive for exothermic adsorption); NaN if it cannot be computed.
    """
    u = np.asarray(energy_trajectory, dtype=float)
    n = np.asarray(loading_trajectory, dtype=float)
    if energy_unit == "K":
        u = u * R_KJ
    elif energy_unit != "kJ/mol":
        raise ValueError("energy_unit must be 'kJ/mol' or 'K'")
    if len(n) < 10 or len(u) != len(n):
        return float("nan")
    var_n = float(np.var(n))
    if var_n < 1e-12:
        return float("nan")
    cov_un = float(np.mean(u * n) - np.mean(u) * np.mean(n))
    return float(R_KJ * temperature_k - cov_un / var_n)


def calculate_isosteric_heat_clausius_clapeyron(
    p1: float,
    p2: float,
    t1: float,
    t2: float
) -> float:
    """
    Isosteric heat between two temperatures at EQUAL loading (Clausius-Clapeyron):
        q_st = -R ln(P2 / P1) / (1/T2 - 1/T1)      [kJ/mol]
    Returns NaN for invalid input; the sign is kept (negative values signal inconsistent data).
    """
    if p1 <= 0 or p2 <= 0 or t1 == t2:
        return float("nan")
    denom = (1.0 / t2) - (1.0 / t1)
    return float(-R_KJ * np.log(p2 / p1) / denom)


def calculate_isosteric_heat(
    energy_series: Optional[np.ndarray] = None,
    loading_series: Optional[np.ndarray] = None,
    temperature_k: float = 298.15,
    energy_unit: str = "kJ/mol",
) -> float:
    """q_st from a GCMC (U, N) series; NaN when no series is available."""
    if energy_series is not None and loading_series is not None:
        return calculate_isosteric_heat_fluctuations(energy_series, loading_series, temperature_k, energy_unit)
    return float("nan")
