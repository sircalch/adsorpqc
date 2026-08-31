"""
Tests for IAST mixture selectivity and bootstrap confidence intervals.
"""

import numpy as np
import pytest
from adsorpqc.core.iast import calculate_iast_selectivity


def test_iast_selectivity_bootstrap():
    pressures = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    
    # Gas A: Strong adsorption (q_sat=5.0, K=5.0)
    # Gas B: Weak adsorption (q_sat=5.0, K=0.5)
    # Expected selectivity S_AB = (5.0 * 5.0) / (5.0 * 0.5) = 10.0
    q_a = 5.0 * (5.0 * pressures) / (1.0 + 5.0 * pressures)
    q_b = 5.0 * (0.5 * pressures) / (1.0 + 0.5 * pressures)
    
    res = calculate_iast_selectivity(
        pressures, q_a, pressures, q_b,
        gas_mole_fraction_a=0.5, gas_mole_fraction_b=0.5,
        n_bootstrap=100
    )
    
    assert np.isclose(res["selectivity_a_b"], 10.0, atol=0.5)
    assert res["ci_lower_95"] <= res["selectivity_a_b"] <= res["ci_upper_95"]
