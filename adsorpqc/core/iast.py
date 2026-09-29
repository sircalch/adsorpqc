"""
Ideal Adsorbed Solution Theory (IAST; Myers & Prausnitz, AIChE J. 11, 121, 1965) for binary
mixtures, from fitted single-component isotherms, with bootstrap uncertainty.

For gas-phase mole fractions y_i at total pressure P, IAST finds adsorbed-phase mole fractions x_i
and pure-component pressures P_i^0 such that

    P y_i = P_i^0 x_i          and      pi_1(P_1^0) = pi_2(P_2^0),
    pi_i(P^0) = int_0^{P^0} q_i(p) / p dp       (reduced spreading pressure),

and the total loading is 1 / q_T = sum_i x_i / q_i(P_i^0). The selectivity is
S_12 = (x_1 / x_2) / (y_1 / y_2). Versions before 1.1.0 returned (q_sat,1 K_1) / (q_sat,2 K_2) under
this name, which is the Henry-limit selectivity and ignores pressure and composition.
"""

from typing import Dict, Any, Callable, Optional
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from adsorpqc.core.isotherms import (fit_langmuir, fit_dual_site_langmuir, fit_sips, fit_toth,
                                     fit_all_isotherm_models, IsothermFitResult)

_FITTERS = {"Langmuir": fit_langmuir, "Dual-Site Langmuir": fit_dual_site_langmuir,
            "Sips": fit_sips, "Toth": fit_toth}


def reduced_spreading_pressure(q: Callable[[float], float], p0: float) -> float:
    """pi(P^0) = int_0^{P^0} q(p)/p dp, integrated in ln p so that the low-pressure end is resolved."""
    if p0 <= 0:
        return 0.0
    lo = np.log(p0) - 40.0                       # q(p)/p -> K_H: the tail below e^-40 P^0 is negligible
    val, _ = quad(lambda u: float(q(np.exp(u))), lo, np.log(p0), limit=400)
    return float(val)


def iast_binary(q1: Callable[[float], float], q2: Callable[[float], float], y1: float, total_pressure: float) -> Dict[str, float]:
    """
    Solves binary IAST for pure-component isotherms q1(p), q2(p) (same pressure and loading units).

    Returns x1, x2, the component loadings q1_mix, q2_mix, the total loading, the pure-component
    pressures P1^0, P2^0 and the selectivity S_12.
    """
    y2 = 1.0 - y1
    if not (0.0 < y1 < 1.0) or total_pressure <= 0:
        raise ValueError("need 0 < y1 < 1 and total_pressure > 0")
    P = float(total_pressure)

    def f(x1):
        return reduced_spreading_pressure(q1, P * y1 / x1) - reduced_spreading_pressure(q2, P * y2 / (1.0 - x1))

    eps = 1e-10
    x1 = brentq(f, eps, 1.0 - eps, xtol=1e-12, maxiter=500)
    x2 = 1.0 - x1
    p1, p2 = P * y1 / x1, P * y2 / x2
    q_total = 1.0 / (x1 / float(q1(p1)) + x2 / float(q2(p2)))
    return {"x1": x1, "x2": x2, "q1_mix": x1 * q_total, "q2_mix": x2 * q_total, "q_total": q_total,
            "p1_pure": p1, "p2_pure": p2, "selectivity_12": (x1 / x2) / (y1 / y2)}


def _fit(p, q, model: str) -> Optional[IsothermFitResult]:
    if model == "best":
        return fit_all_isotherm_models(p, q).get("best_model")
    return _FITTERS[model](p, q)


def calculate_iast_selectivity(
    isotherm_a_p: np.ndarray,
    isotherm_a_q: np.ndarray,
    isotherm_b_p: np.ndarray,
    isotherm_b_q: np.ndarray,
    gas_mole_fraction_a: float = 0.5,
    gas_mole_fraction_b: float = 0.5,
    total_pressure: float = 1.0,
    n_bootstrap: int = 500,
    model: str = "best",
    random_state: int = 42,
) -> Dict[str, Any]:
    """
    Binary IAST selectivity S_A/B = (x_A/x_B)/(y_A/y_B) at the given total pressure and composition.

    Each single-component isotherm is fitted (model: "best" = lowest AICc among Langmuir, dual-site
    Langmuir, Sips and Toth, or one of those names) and IAST is solved with the fitted functions.
    The 95% interval comes from resampling the isotherm points of each component and repeating the
    fit and the IAST solution (pairs bootstrap); it reflects the fit uncertainty only.

    Returns
    -------
    dict with selectivity_a_b, ci_lower_95, ci_upper_95, x_a, loadings of both components in the
    mixture, the models used, and the Henry-limit selectivity for comparison.
    """
    pa, qa = np.asarray(isotherm_a_p, float), np.asarray(isotherm_a_q, float)
    pb, qb = np.asarray(isotherm_b_p, float), np.asarray(isotherm_b_q, float)
    y_a = float(gas_mole_fraction_a)
    if gas_mole_fraction_b is not None and abs(y_a + float(gas_mole_fraction_b) - 1.0) > 1e-9:
        raise ValueError("gas-phase mole fractions must add up to 1 for a binary mixture")

    fa, fb = _fit(pa, qa, model), _fit(pb, qb, model)
    if fa is None or fb is None:
        raise ValueError("single-component isotherm fit failed; IAST needs a fitted model for both components")
    sol = iast_binary(fa.fitted_function, fb.fitted_function, y_a, total_pressure)

    # Henry-limit selectivity (the quantity returned by versions < 1.1.0)
    h = 1e-9 * max(pa.max(), pb.max())
    henry_sel = float(fa.fitted_function(h) / h) / float(fb.fitted_function(h) / h)

    rng = np.random.default_rng(random_state)
    boot = []
    for _ in range(n_bootstrap):
        ia = rng.choice(len(pa), size=len(pa), replace=True)
        ib = rng.choice(len(pb), size=len(pb), replace=True)
        if len(np.unique(ia)) < 3 or len(np.unique(ib)) < 3:
            continue
        ba, bb = _fit(pa[ia], qa[ia], fa.model_name), _fit(pb[ib], qb[ib], fb.model_name)
        if ba is None or bb is None:
            continue
        try:
            boot.append(iast_binary(ba.fitted_function, bb.fitted_function, y_a, total_pressure)["selectivity_12"])
        except (ValueError, RuntimeError, ZeroDivisionError):
            continue
    lo, hi = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))) if len(boot) >= 20 else (float("nan"), float("nan"))

    return {
        "selectivity_a_b": float(sol["selectivity_12"]),
        "ci_lower_95": lo,
        "ci_upper_95": hi,
        "n_bootstrap_ok": len(boot),
        "method": "IAST",
        "model_a": fa.model_name,
        "model_b": fb.model_name,
        "x_a": sol["x1"],
        "q_a_mix": sol["q1_mix"],
        "q_b_mix": sol["q2_mix"],
        "q_total": sol["q_total"],
        "henry_limit_selectivity": henry_sel,
        "total_pressure": float(total_pressure),
        "y_a": y_a,
    }
