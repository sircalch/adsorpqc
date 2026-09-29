"""
Henry coefficient (K_H) and isosteric heat of adsorption (q_st) calculations.
"""

from typing import Dict, Any, Optional, Tuple, List
import numpy as np


def calculate_henry_constant(
    pressure: np.ndarray,
    loading: np.ndarray,
    n_points_linear: int = 4
) -> Dict[str, float]:
    """
    Calculates the Henry adsorption constant K_H by linear regression in the low-pressure regime.

    Parameters
    ----------
    pressure : np.ndarray
        Pressure values.
    loading : np.ndarray
        Loading values.
    n_points_linear : int, default 4
        Number of initial points to fit.

    Returns
    -------
    result : dict
        K_H slope, R^2, and intercept.
    """
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    
    # Sort by pressure
    sort_idx = np.argsort(p)
    p_s = p[sort_idx]
    q_s = q[sort_idx]
    
    k_pts = min(len(p_s), max(2, n_points_linear))
    p_sub = p_s[:k_pts]
    q_sub = q_s[:k_pts]
    
    # Linear fit through origin: q = K_H * P
    slope, residuals, _, _ = np.linalg.lstsq(p_sub[:, np.newaxis], q_sub, rcond=None)
    k_h = float(slope[0])
    
    # Unconstrained linear fit for R^2
    poly = np.polyfit(p_sub, q_sub, 1)
    y_pred = np.polyval(poly, p_sub)
    ss_res = np.sum((q_sub - y_pred) ** 2)
    ss_tot = np.sum((q_sub - np.mean(q_sub)) ** 2)
    r_sq = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    
    return {
        "henry_constant": max(0.0, k_h),
        "r_squared": float(r_sq),
        "n_points_used": k_pts
    }


def calculate_isosteric_heat_fluctuations(
    energy_trajectory: np.ndarray,
    loading_trajectory: np.ndarray,
    temperature_k: float = 298.15
) -> float:
    """
    Calculates the isosteric heat of adsorption q_st via particle-energy fluctuations
    in the Grand Canonical Monte Carlo (mu, V, T) ensemble:
        q_st = R * T - (<U * N> - <U><N>) / (<N^2> - <N>^2)

    Parameters
    ----------
    energy_trajectory : np.ndarray
        Total host-adsorbate interaction energy (kJ/mol or K).
    loading_trajectory : np.ndarray
        Number of adsorbed molecules N.
    temperature_k : float, default 298.15
        Temperature in Kelvin.

    Returns
    -------
    q_st : float
        Isosteric heat of adsorption in kJ/mol.
    """
    u = np.asarray(energy_trajectory, dtype=float)
    n = np.asarray(loading_trajectory, dtype=float)
    
    if len(n) < 10:
        return 0.0
        
    var_n = np.var(n)
    if var_n < 1e-8:
        return 0.0
        
    cov_un = np.mean(u * n) - np.mean(u) * np.mean(n)
    r_gas_kj = 8.314462618e-3  # kJ/(mol*K)
    
    # q_st = R*T - Cov(U, N)/Var(N)
    q_st = r_gas_kj * temperature_k - (cov_un / var_n)
    return float(abs(q_st))


def calculate_isosteric_heat_clausius_clapeyron(
    p1: float,
    p2: float,
    t1: float,
    t2: float
) -> float:
    """
    Calculates q_st between two temperatures at equal loading via the Clausius-Clapeyron equation:
        q_st = -R * ln(P2 / P1) / (1/T2 - 1/T1)
    """
    r_gas_kj = 8.314462618e-3
    if p1 <= 0 or p2 <= 0 or t1 == t2:
        return 0.0
    denom = (1.0 / t2) - (1.0 / t1)
    q_st = -r_gas_kj * np.log(p2 / p1) / denom
    return float(abs(q_st))


def calculate_isosteric_heat(
    energy_series: Optional[np.ndarray] = None,
    loading_series: Optional[np.ndarray] = None,
    temperature_k: float = 298.15
) -> float:
    """Calculates isosteric heat of adsorption with fallback."""
    if energy_series is not None and loading_series is not None:
        return calculate_isosteric_heat_fluctuations(energy_series, loading_series, temperature_k)
    return 0.0
