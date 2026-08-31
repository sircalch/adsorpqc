"""
Tests for GCMC burn-in detection and loading drift evaluation.
"""

import numpy as np
import pytest
from adsorpqc.core.burnin import detect_gcmc_burnin, assess_loading_drift


def test_detect_gcmc_burnin_converged():
    # Synthetic trajectory: 1000 steps initial approach + 3000 steps stationary production around mean = 10.0
    rng = np.random.default_rng(42)
    n = 3000
    t = np.arange(n)
    loading = 10.0 * (1.0 - np.exp(-t / 200.0)) + rng.normal(0, 0.2, size=n)
    
    res = detect_gcmc_burnin(loading)
    assert res.status == "PASS"
    assert res.n_burnin_cycles > 100
    assert np.isclose(res.production_mean_loading, 10.0, atol=0.2)
    assert res.loading_drift_pct < 3.0


def test_assess_loading_drift_severe_drift():
    # Continual linear rise from 5.0 to 15.0 -> severe drift -> FAIL
    n = 2000
    loading = np.linspace(5.0, 15.0, n)
    
    drift_pct, m_init, m_final, status, msg = assess_loading_drift(loading)
    assert status == "FAIL"
    assert drift_pct > 8.0
    assert "DO NOT REPORT" in msg
