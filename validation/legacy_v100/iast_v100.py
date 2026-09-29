"""
Ideal Adsorbed Solution Theory (IAST) and binary selectivity with bootstrap uncertainty.
"""

from typing import Dict, Any, Optional, Tuple, List
import numpy as np
from adsorpqc.core.isotherms import fit_langmuir, IsothermFitResult


def calculate_iast_selectivity(
    isotherm_a_p: np.ndarray,
    isotherm_a_q: np.ndarray,
    isotherm_b_p: np.ndarray,
    isotherm_b_q: np.ndarray,
    gas_mole_fraction_a: float = 0.5,
    gas_mole_fraction_b: float = 0.5,
    total_pressure: float = 1.0,
    n_bootstrap: int = 500
) -> Dict[str, Any]:
    """
    Computes binary mixture selectivity S_{A/B} = (x_A / x_B) / (y_A / y_B) and
    IAST selectivity with 95% bootstrap confidence intervals.

    Parameters
    ----------
    isotherm_a_p, isotherm_a_q : np.ndarray
        Single-component isotherm data for adsorbate A.
    isotherm_b_p, isotherm_b_q : np.ndarray
        Single-component isotherm data for adsorbate B.
    gas_mole_fraction_a : float, default 0.5
        Gas phase mole fraction y_A.
    gas_mole_fraction_b : float, default 0.5
        Gas phase mole fraction y_B.
    total_pressure : float, default 1.0
        Total system pressure.
    n_bootstrap : int, default 500
        Number of bootstrap replicates.

    Returns
    -------
    result : dict
        Selectivity value, 95% CI, fitted affinity ratios.
    """
    y_a = float(gas_mole_fraction_a)
    y_b = float(gas_mole_fraction_b)
    
    # Fit single components
    fit_a = fit_langmuir(isotherm_a_p, isotherm_a_q)
    fit_b = fit_langmuir(isotherm_b_p, isotherm_b_q)
    
    if fit_a is None or fit_b is None:
        # Fallback to empirical ratio at closest pressure
        p_eval = total_pressure * y_a
        q_a = float(np.interp(p_eval, isotherm_a_p, isotherm_a_q))
        q_b = float(np.interp(total_pressure * y_b, isotherm_b_p, isotherm_b_q))
        if q_b > 0 and y_b > 0:
            sel = (q_a / q_b) / (y_a / y_b)
        else:
            sel = 1.0
        return {
            "selectivity_a_b": float(sel),
            "ci_lower_95": float(sel),
            "ci_upper_95": float(sel),
            "method": "Empirical Loading Ratio",
            "k_a": 0.0,
            "k_b": 0.0
        }
        
    k_a = fit_a.parameters["K"]
    k_b = fit_b.parameters["K"]
    q_sat_a = fit_a.parameters["q_sat"]
    q_sat_b = fit_b.parameters["q_sat"]
    
    # Langmuirian IAST selectivity S_AB = (q_sat_A * K_A) / (q_sat_B * K_B) or K_A / K_B
    selectivity_nominal = (k_a / (k_b + 1e-12)) * (q_sat_a / (q_sat_b + 1e-12))
    
    # Bootstrap uncertainty
    n_a = len(isotherm_a_p)
    n_b = len(isotherm_b_p)
    rng = np.random.default_rng(42)
    
    boot_sels = []
    for _ in range(n_bootstrap):
        idx_a = rng.choice(n_a, size=n_a, replace=True)
        idx_b = rng.choice(n_b, size=n_b, replace=True)
        
        b_fit_a = fit_langmuir(isotherm_a_p[idx_a], isotherm_a_q[idx_a])
        b_fit_b = fit_langmuir(isotherm_b_p[idx_b], isotherm_b_q[idx_b])
        
        if b_fit_a and b_fit_b:
            ka = b_fit_a.parameters["K"]
            kb = b_fit_b.parameters["K"]
            qsa = b_fit_a.parameters["q_sat"]
            qsb = b_fit_b.parameters["q_sat"]
            s_val = (ka / (kb + 1e-12)) * (qsa / (qsb + 1e-12))
            if 0 < s_val < 1e6:
                boot_sels.append(s_val)

    if boot_sels:
        ci_low = float(np.percentile(boot_sels, 2.5))
        ci_high = float(np.percentile(boot_sels, 97.5))
    else:
        ci_low = selectivity_nominal
        ci_high = selectivity_nominal

    return {
        "selectivity_a_b": float(selectivity_nominal),
        "ci_lower_95": ci_low,
        "ci_upper_95": ci_high,
        "method": "IAST (Langmuirian)",
        "k_a": float(k_a),
        "k_b": float(k_b),
        "q_sat_a": float(q_sat_a),
        "q_sat_b": float(q_sat_b)
    }
