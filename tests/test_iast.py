"""
Tests for binary IAST.

References: for two Langmuir components with equal saturation capacity, IAST reduces exactly to
extended Langmuir (S = K_A / K_B at every pressure and composition). For other cases, the mixture
loadings computed by pyIAST 1.4 (Simon et al., Comput. Phys. Commun. 2016) with the exact model
parameters (validation/iast_vs_pyiast.py) are used.
"""

import numpy as np
import pytest
from adsorpqc.core.iast import calculate_iast_selectivity, iast_binary, reduced_spreading_pressure


def langmuir(qs, k):
    return lambda p: qs * k * p / (1 + k * p)


def test_iast_selectivity_bootstrap():
    pressures = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    # equal q_sat: IAST = extended Langmuir, S = K_A / K_B = 10 at any pressure
    q_a = 5.0 * (5.0 * pressures) / (1.0 + 5.0 * pressures)
    q_b = 5.0 * (0.5 * pressures) / (1.0 + 0.5 * pressures)
    res = calculate_iast_selectivity(pressures, q_a, pressures, q_b, gas_mole_fraction_a=0.5,
                                     gas_mole_fraction_b=0.5, total_pressure=2.0, n_bootstrap=100,
                                     model="Langmuir")
    assert res["selectivity_a_b"] == pytest.approx(10.0, rel=1e-4)
    assert res["ci_lower_95"] <= res["selectivity_a_b"] <= res["ci_upper_95"]


def test_spreading_pressure_langmuir_analytic():
    assert reduced_spreading_pressure(langmuir(3.0, 2.0), 1.5) == pytest.approx(3.0 * np.log(1 + 2.0 * 1.5), rel=1e-10)


@pytest.mark.parametrize("f1,f2,q1_ref,q2_ref", [
    (langmuir(2.0, 8.0), langmuir(6.0, 0.3), 1.4560578873852428, 0.3318493287896544),
    (lambda p: 1.5 * 20 * p / (1 + 20 * p) + 2.0 * 0.5 * p / (1 + 0.5 * p),
     lambda p: 3.0 * 0.2 * p / (1 + 0.2 * p) + 1.0 * 2.0 * p / (1 + 2.0 * p),
     1.6564196804054263, 0.2047101888398731),
])
def test_iast_matches_pyiast(f1, f2, q1_ref, q2_ref):
    sol = iast_binary(f1, f2, 0.5, 1.0)
    assert sol["q1_mix"] == pytest.approx(q1_ref, rel=1e-9)
    assert sol["q2_mix"] == pytest.approx(q2_ref, rel=1e-9)


def test_selectivity_depends_on_pressure():
    """The Henry-limit ratio returned by versions < 1.1.0 is not the IAST selectivity."""
    p = np.logspace(-2, 1, 12)
    res_lo = calculate_iast_selectivity(p, langmuir(2.0, 8.0)(p), p, langmuir(6.0, 0.3)(p),
                                        total_pressure=0.01, n_bootstrap=0, model="Langmuir")
    res_hi = calculate_iast_selectivity(p, langmuir(2.0, 8.0)(p), p, langmuir(6.0, 0.3)(p),
                                        total_pressure=10.0, n_bootstrap=0, model="Langmuir")
    assert res_lo["selectivity_a_b"] == pytest.approx(8.76, abs=0.01)
    assert res_hi["selectivity_a_b"] == pytest.approx(1.15, abs=0.01)
    assert res_hi["henry_limit_selectivity"] == pytest.approx(8.89, abs=0.01)
