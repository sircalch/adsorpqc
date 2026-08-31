"""
Tests for Henry constant K_H and isosteric heat q_st.
"""

import numpy as np
import pytest
from adsorpqc.core.energetics import (
    calculate_henry_constant,
    calculate_isosteric_heat_fluctuations,
    calculate_isosteric_heat_clausius_clapeyron
)


def test_henry_constant():
    # Low pressure linear regime: q = 0.5 * P
    p = np.array([0.01, 0.02, 0.05, 0.1, 1.0, 5.0])
    q = 0.5 * p
    
    res = calculate_henry_constant(p, q, n_points_linear=4)
    assert np.isclose(res["henry_constant"], 0.5, atol=1e-3)
    assert res["r_squared"] > 0.999


def test_isosteric_heat_fluctuations():
    # Synthetic GCMC energy and particle count with covariance
    rng = np.random.default_rng(42)
    n_pts = 1000
    n_particles = rng.poisson(lam=20, size=n_pts).astype(float)
    # Energy proportional to particles: U = -30 * N + noise
    u_energy = -30.0 * n_particles + rng.normal(0, 2.0, size=n_pts)
    
    qst = calculate_isosteric_heat_fluctuations(u_energy, n_particles, temperature_k=298.15)
    # q_st should be approx 30 kJ/mol + RT
    assert 25.0 < qst < 35.0


def test_isosteric_heat_clausius_clapeyron():
    # T1 = 273.15 K, P1 = 0.5 bar; T2 = 298.15 K, P2 = 1.2 bar
    qst = calculate_isosteric_heat_clausius_clapeyron(0.5, 1.2, 273.15, 298.15)
    assert qst > 0.0
    assert 20.0 < qst < 35.0
