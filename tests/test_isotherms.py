"""
Tests for non-linear isotherm model fitting (Langmuir, Dual-Site, Sips, Toth).
"""

import numpy as np
import pytest
from adsorpqc.core.isotherms import (
    fit_langmuir,
    fit_dual_site_langmuir,
    fit_sips,
    fit_all_isotherm_models
)


def test_fit_langmuir_synthetic():
    p = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    # Exact Langmuir: q_sat = 5.0, K = 2.0
    q_exact = 5.0 * (2.0 * p) / (1.0 + 2.0 * p)
    
    res = fit_langmuir(p, q_exact)
    assert res is not None
    assert np.isclose(res.parameters["q_sat"], 5.0, atol=0.01)
    assert np.isclose(res.parameters["K"], 2.0, atol=0.01)
    assert res.r_squared > 0.999


def test_fit_all_models_selection():
    p = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    # Dual-site Langmuir synthetic data
    q = 4.0 * (10.0 * p) / (1.0 + 10.0 * p) + 2.0 * (0.5 * p) / (1.0 + 0.5 * p)
    
    res_all = fit_all_isotherm_models(p, q)
    assert "models" in res_all
    assert len(res_all["models"]) >= 2
    assert res_all["best_r_squared"] > 0.98
    assert res_all["best_model_name"] in ["Dual-Site Langmuir", "Sips", "Toth", "Langmuir"]
