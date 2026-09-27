"""
GCMC cycle burn-in detection, production stationarity, and loading drift evaluation.
"""

from typing import Dict, Any, Optional, Tuple, List
from dataclasses import dataclass
import numpy as np


@dataclass
class GCMCBurninResult:
    n_total_cycles: int
    n_burnin_cycles: int
    n_production_cycles: int
    burnin_fraction: float
    production_mean_loading: float
    production_std_loading: float
    production_ci_lower_95: float
    production_ci_upper_95: float
    loading_drift_pct: float
    initial_20pct_mean: float
    final_20pct_mean: float
    tau_int_cycles: float
    n_eff: float
    status: str  # 'PASS', 'WARNING', 'FAIL'
    diagnostic_message: str


def compute_autocorr_1d(series: np.ndarray, max_lag: int = 200) -> np.ndarray:
    """Computes normalized autocorrelation for GCMC timeseries."""
    x = np.asarray(series, dtype=float)
    n = len(x)
    if n < 4:
        return np.array([1.0])
    x_zm = x - np.mean(x)
    var = np.var(x)
    if var == 0:
        return np.ones(min(n, max_lag))
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    fx = np.fft.fft(x_zm, n=n_fft)
    px = fx * np.conjugate(fx)
    autocov = np.fft.ifft(px).real[:min(n, max_lag)]
    denom = (n - np.arange(len(autocov))) * var
    acf = autocov / denom
    if acf[0] != 0:
        acf = acf / acf[0]
    return acf


def estimate_tau_int(series: np.ndarray, c_window: float = 6.0) -> Tuple[float, float]:
    """Calculates integrated autocorrelation time and statistical inefficiency."""
    x = np.asarray(series, dtype=float)
    n = len(x)
    if n < 10:
        return 0.5, 1.0
    acf = compute_autocorr_1d(x, max_lag=min(n // 2, 500))
    tau = 0.5
    for k in range(1, len(acf)):
        if acf[k] < 0:
            break
        tau += acf[k]
        if k >= c_window * tau:
            break
    tau = max(0.5, float(tau))
    g = max(1.0, 1.0 + 2.0 * (tau - 0.5))
    return tau, g


def assess_loading_drift(
    production_series: np.ndarray,
    warn_threshold_pct: float = 3.0,
    fail_threshold_pct: float = 8.0,
    z_warn: float = 2.0,
    z_fail: float = 3.0,
) -> Tuple[float, float, float, str, str]:
    """
    Drift between the first and last 20% of the production phase, judged both statistically and
    in magnitude.

    The difference of the two segment means is compared with its standard error. The segment
    standard errors use the scatter and the statistical inefficiency of the linearly detrended
    series, so that a genuine trend is not absorbed into the noise estimate and pure noise in a
    short or low-loading run is not mistaken for drift. A drift is reported only when it is both
    significant (z above z_warn / z_fail) and large relative to the mean (above the percentage
    thresholds).

    Returns
    -------
    drift_pct, mean_first_20, mean_last_20, status, message
    """
    y = np.asarray(production_series, dtype=float)
    n = len(y)
    if n < 10:
        mean_all = float(np.mean(y)) if n > 0 else 0.0
        return 0.0, mean_all, mean_all, "PASS", "Too few samples to evaluate drift."

    n_seg = max(2, int(0.20 * n))
    mean_init = float(np.mean(y[:n_seg]))
    mean_final = float(np.mean(y[-n_seg:]))
    mean_prod = float(np.mean(y))

    t = np.arange(n, dtype=float)
    resid = y - np.polyval(np.polyfit(t, y, 1), t)
    sd = float(np.std(resid, ddof=2))
    _, g = estimate_tau_int(resid)
    se_seg = sd * np.sqrt(g / n_seg)
    diff = abs(mean_final - mean_init)
    z = diff / (np.sqrt(2.0) * se_seg) if se_seg > 0 else (np.inf if diff > 0 else 0.0)
    drift_pct = diff / max(abs(mean_prod), 1e-12) * 100.0 if mean_prod != 0 else (0.0 if diff == 0 else np.inf)

    if z > z_fail and drift_pct > fail_threshold_pct:
        status = "FAIL"
        msg = (f"CRITICAL: significant loading drift ({drift_pct:.2f}% of the mean, z = {z:.1f}). "
               "DO NOT REPORT PRODUCTION AVERAGE: GCMC remained in a non-equilibrium state.")
    elif z > z_warn and drift_pct > warn_threshold_pct:
        status = "WARNING"
        msg = f"Loading drift detected ({drift_pct:.2f}% of the mean, z = {z:.1f}). Consider extending equilibration cycles."
    else:
        status = "PASS"
        msg = (f"Stationary production phase (first/last 20% differ by {drift_pct:.2f}% of the mean, "
               f"z = {z:.1f}, within statistical noise or below {warn_threshold_pct}%).")

    return float(drift_pct), mean_init, mean_final, status, msg


def detect_gcmc_burnin(
    loading_series: np.ndarray,
    cycle_indices: Optional[np.ndarray] = None,
    step_search: int = 10,
    warn_drift_pct: float = 3.0,
    fail_drift_pct: float = 8.0
) -> GCMCBurninResult:
    """
    Automatically detects the optimal GCMC burn-in / equilibration cutoff and verifies production stationarity.

    Parameters
    ----------
    loading_series : np.ndarray
        Array of adsorbate loadings per cycle (e.g. molecules/unit cell, mol/kg).
    cycle_indices : np.ndarray, optional
        Cycle numbers.
    step_search : int, default 10
        Step increment for candidate cutoff search.

    Returns
    -------
    result : GCMCBurninResult
        Detailed burn-in, effective sampling, drift, and quality certification.
    """
    y = np.asarray(loading_series, dtype=float)
    n = len(y)
    
    if n < 20:
        tau, g = estimate_tau_int(y)
        m = float(np.mean(y)) if n > 0 else 0.0
        return GCMCBurninResult(
            n_total_cycles=n,
            n_burnin_cycles=0,
            n_production_cycles=n,
            burnin_fraction=0.0,
            production_mean_loading=m,
            production_std_loading=float(np.std(y, ddof=1)) if n > 1 else 0.0,
            production_ci_lower_95=m,
            production_ci_upper_95=m,
            loading_drift_pct=0.0,
            initial_20pct_mean=m,
            final_20pct_mean=m,
            tau_int_cycles=tau,
            n_eff=float(n / g),
            status="PASS",
            diagnostic_message="Short simulation run."
        )
        
    # Search cutoff t_0 up to 80% of data
    max_idx = int(0.80 * n)
    indices = np.arange(0, max_idx, max(1, step_search))
    
    best_t0 = 0
    best_neff = -1.0
    best_tau = 0.5
    best_g = 1.0
    
    for idx in indices:
        y_sub = y[idx:]
        if len(y_sub) < 10:
            break
        tau, g = estimate_tau_int(y_sub)
        neff = float(len(y_sub)) / g
        if neff > best_neff:
            best_neff = neff
            best_t0 = idx
            best_tau = tau
            best_g = g

    y_prod = y[best_t0:]
    n_prod = len(y_prod)
    burnin_frac = float(best_t0) / float(n)
    
    mean_prod = float(np.mean(y_prod))
    std_prod = float(np.std(y_prod, ddof=1)) if n_prod > 1 else 0.0
    
    # 95% Confidence Interval for mean using standard error accounting for g
    se_mean = (std_prod * np.sqrt(best_g)) / np.sqrt(n_prod) if n_prod > 1 else 0.0
    ci_low = mean_prod - 1.96 * se_mean
    ci_high = mean_prod + 1.96 * se_mean
    
    # Drift evaluation
    drift_pct, m_init, m_final, drift_st, drift_msg = assess_loading_drift(
        y_prod, warn_threshold_pct=warn_drift_pct, fail_threshold_pct=fail_drift_pct
    )
    
    status = drift_st
    if best_neff < 50:
        if status != "FAIL":
            status = "WARNING"
        diag_msg = f"{drift_msg} Marginal statistical sampling (N_eff = {best_neff:.0f} < 50 independent cycles)."
    else:
        diag_msg = f"{drift_msg} Burn-in discarded: {best_t0} cycles ({burnin_frac*100:.1f}%). N_eff = {best_neff:.0f} independent samples."

    return GCMCBurninResult(
        n_total_cycles=n,
        n_burnin_cycles=int(best_t0),
        n_production_cycles=int(n_prod),
        burnin_fraction=burnin_frac,
        production_mean_loading=mean_prod,
        production_std_loading=std_prod,
        production_ci_lower_95=ci_low,
        production_ci_upper_95=ci_high,
        loading_drift_pct=drift_pct,
        initial_20pct_mean=m_init,
        final_20pct_mean=m_final,
        tau_int_cycles=best_tau,
        n_eff=best_neff,
        status=status,
        diagnostic_message=diag_msg
    )
