"""
Non-linear adsorption isotherm model fitting and Akaike Information Criterion (AIC) selection.
"""

from typing import Dict, Any, Optional, Tuple, List, Callable
from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import curve_fit


@dataclass
class IsothermFitResult:
    model_name: str
    parameters: Dict[str, float]
    r_squared: float
    rmse: float
    aic: float
    bic: float
    fitted_function: Callable[[np.ndarray], np.ndarray]

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d.pop("fitted_function", None)
        return d


# 1. Model equations
def langmuir_eq(p: np.ndarray, q_sat: float, k: float) -> np.ndarray:
    return q_sat * (k * p) / (1.0 + k * p)


def dual_site_langmuir_eq(p: np.ndarray, q1: float, k1: float, q2: float, k2: float) -> np.ndarray:
    return q1 * (k1 * p) / (1.0 + k1 * p) + q2 * (k2 * p) / (1.0 + k2 * p)


def sips_eq(p: np.ndarray, q_sat: float, k: float, n: float) -> np.ndarray:
    kp_n = (np.maximum(1e-12, k * p)) ** n
    return q_sat * kp_n / (1.0 + kp_n)


def toth_eq(p: np.ndarray, q_sat: float, k: float, t: float) -> np.ndarray:
    kp = np.maximum(1e-12, k * p)
    return q_sat * kp / ((1.0 + kp ** t) ** (1.0 / t))


def freundlich_eq(p: np.ndarray, k: float, n_inv: float) -> np.ndarray:
    return k * (np.maximum(1e-12, p) ** n_inv)


def _compute_fit_metrics(
    y_obs: np.ndarray,
    y_pred: np.ndarray,
    k_params: int
) -> Tuple[float, float, float, float]:
    """
    Computes R^2, RMSE, AIC and BIC. The AIC is the small-sample corrected AICc
    (Hurvich & Tsai 1989), AIC + 2k(k+1)/(n-k-1), since isotherms rarely have n/k > 40;
    it is infinite when n <= k + 1.
    """
    n = len(y_obs)
    residuals = y_obs - y_pred
    rss = float(np.sum(residuals ** 2))
    
    ss_tot = float(np.sum((y_obs - np.mean(y_obs)) ** 2))
    r_sq = 1.0 - (rss / ss_tot) if ss_tot > 0 else 0.0
    rmse = np.sqrt(rss / n)
    
    # AICc = 2k + n ln(RSS/n) + 2k(k+1)/(n-k-1)
    if rss > 0 and n > k_params + 1:
        aic = 2.0 * k_params + n * np.log(rss / n) + 2.0 * k_params * (k_params + 1) / (n - k_params - 1)
        bic = k_params * np.log(n) + n * np.log(rss / n)
    else:
        aic = float("inf")
        bic = float("inf")
        
    return float(r_sq), float(rmse), float(aic), float(bic)


def fit_langmuir(pressure: np.ndarray, loading: np.ndarray) -> Optional[IsothermFitResult]:
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    if len(p) < 3:
        return None
    try:
        q_max = float(np.max(q))
        p0 = [q_max * 1.1, 1.0 / (np.median(p) + 1e-6)]
        bounds = ([0.0, 0.0], [q_max * 10.0, 1e6])
        popt, _ = curve_fit(langmuir_eq, p, q, p0=p0, bounds=bounds, maxfev=5000)
        q_sat, k = popt
        y_pred = langmuir_eq(p, q_sat, k)
        r_sq, rmse, aic, bic = _compute_fit_metrics(q, y_pred, k_params=2)
        return IsothermFitResult(
            model_name="Langmuir",
            parameters={"q_sat": float(q_sat), "K": float(k)},
            r_squared=r_sq,
            rmse=rmse,
            aic=aic,
            bic=bic,
            fitted_function=lambda p_eval: langmuir_eq(p_eval, q_sat, k)
        )
    except Exception:
        return None


def fit_dual_site_langmuir(pressure: np.ndarray, loading: np.ndarray) -> Optional[IsothermFitResult]:
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    if len(p) < 5:
        return None
    try:
        q_max = float(np.max(q))
        p0 = [q_max * 0.5, 10.0 / (np.mean(p) + 1e-6), q_max * 0.5, 0.1 / (np.mean(p) + 1e-6)]
        bounds = ([0.0, 0.0, 0.0, 0.0], [q_max * 5.0, 1e6, q_max * 5.0, 1e6])
        popt, _ = curve_fit(dual_site_langmuir_eq, p, q, p0=p0, bounds=bounds, maxfev=10000)
        q1, k1, q2, k2 = popt
        y_pred = dual_site_langmuir_eq(p, q1, k1, q2, k2)
        r_sq, rmse, aic, bic = _compute_fit_metrics(q, y_pred, k_params=4)
        return IsothermFitResult(
            model_name="Dual-Site Langmuir",
            parameters={"q1": float(q1), "K1": float(k1), "q2": float(q2), "K2": float(k2), "q_sat_total": float(q1 + q2)},
            r_squared=r_sq,
            rmse=rmse,
            aic=aic,
            bic=bic,
            fitted_function=lambda p_eval: dual_site_langmuir_eq(p_eval, q1, k1, q2, k2)
        )
    except Exception:
        return None


def fit_sips(pressure: np.ndarray, loading: np.ndarray) -> Optional[IsothermFitResult]:
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    if len(p) < 4:
        return None
    try:
        q_max = float(np.max(q))
        p0 = [q_max * 1.1, 1.0 / (np.median(p) + 1e-6), 1.0]
        bounds = ([0.0, 0.0, 0.1], [q_max * 10.0, 1e6, 5.0])
        popt, _ = curve_fit(sips_eq, p, q, p0=p0, bounds=bounds, maxfev=8000)
        q_sat, k, n = popt
        y_pred = sips_eq(p, q_sat, k, n)
        r_sq, rmse, aic, bic = _compute_fit_metrics(q, y_pred, k_params=3)
        return IsothermFitResult(
            model_name="Sips",
            parameters={"q_sat": float(q_sat), "K": float(k), "n": float(n)},
            r_squared=r_sq,
            rmse=rmse,
            aic=aic,
            bic=bic,
            fitted_function=lambda p_eval: sips_eq(p_eval, q_sat, k, n)
        )
    except Exception:
        return None


def fit_toth(pressure: np.ndarray, loading: np.ndarray) -> Optional[IsothermFitResult]:
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    if len(p) < 4:
        return None
    try:
        q_max = float(np.max(q))
        p0 = [q_max * 1.1, 1.0 / (np.median(p) + 1e-6), 1.0]
        bounds = ([0.0, 0.0, 0.1], [q_max * 10.0, 1e6, 3.0])
        popt, _ = curve_fit(toth_eq, p, q, p0=p0, bounds=bounds, maxfev=8000)
        q_sat, k, t = popt
        y_pred = toth_eq(p, q_sat, k, t)
        r_sq, rmse, aic, bic = _compute_fit_metrics(q, y_pred, k_params=3)
        return IsothermFitResult(
            model_name="Toth",
            parameters={"q_sat": float(q_sat), "K": float(k), "t": float(t)},
            r_squared=r_sq,
            rmse=rmse,
            aic=aic,
            bic=bic,
            fitted_function=lambda p_eval: toth_eq(p_eval, q_sat, k, t)
        )
    except Exception:
        return None


def fit_all_isotherm_models(
    pressure: np.ndarray,
    loading: np.ndarray
) -> Dict[str, Any]:
    """
    Fits the Langmuir, dual-site Langmuir, Sips and Toth models and selects the one with the lowest AICc.

    Parameters
    ----------
    pressure : np.ndarray
        Equilibrium pressure points (bar, Pa, kPa).
    loading : np.ndarray
        Adsorbate loading values (mol/kg, mg/g, molec/uc).

    Returns
    -------
    result : dict
        All fitted models, metrics, and identified optimal model.
    """
    p = np.asarray(pressure, dtype=float)
    q = np.asarray(loading, dtype=float)
    
    models = {}
    
    m_lang = fit_langmuir(p, q)
    if m_lang:
        models["Langmuir"] = m_lang
        
    m_dslg = fit_dual_site_langmuir(p, q)
    if m_dslg:
        models["Dual-Site Langmuir"] = m_dslg
        
    m_sips = fit_sips(p, q)
    if m_sips:
        models["Sips"] = m_sips
        
    m_toth = fit_toth(p, q)
    if m_toth:
        models["Toth"] = m_toth

    if not models:
        return {"models": {}, "best_model": None, "best_model_name": "None"}
        
    # Pick lowest AIC
    best_name = min(models.keys(), key=lambda k: models[k].aic)
    best_fit = models[best_name]
    
    return {
        "models": models,
        "best_model": best_fit,
        "best_model_name": best_name,
        "best_r_squared": best_fit.r_squared,
        "best_rmse": best_fit.rmse,
        "best_aic": best_fit.aic
    }
