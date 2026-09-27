"""
Real RASPA3 3.0.29 runs: methane in silicalite-1 (MFI) at 300 K, GCMC at 1 kPa - 1 MPa and a Widom
run (validation/raspa3_mfi_methane.py). References are RASPA3's own block averages and its
Widom Henry coefficient.
"""
import glob
import os
import numpy as np
import pytest
from adsorpqc.parsers.raspa import parse_raspa_output, collect_raspa_isotherm
from adsorpqc.core.energetics import calculate_henry_constant, calculate_isosteric_heat_fluctuations
from adsorpqc.core.burnin import detect_gcmc_burnin

D = os.path.join(os.path.dirname(__file__), "data", "raspa3_mfi_ch4")


def _one(name):
    return parse_raspa_output(glob.glob(os.path.join(D, name, "*.txt"))[0])


def test_parse_gcmc_output():
    d = _one("p_100000")
    m = d["metadata"]
    assert (m["engine"], m["framework"], m["adsorbate"]) == ("RASPA3", "MFI_SI", "methane")
    assert m["temperature_k"] == 300.0 and m["pressure_pa"] == 1e5
    assert (m["init_cycles"], m["prod_cycles"]) == (5000, 40000)
    assert d["loading_mol_kg"] == pytest.approx(0.50371, abs=1e-5)
    assert d["isosteric_heat_kj_mol"] == pytest.approx(19.29, abs=0.01)
    n = len(d["cycle_loadings"])
    assert n == len(d["cycle_molecules"]) == len(d["cycle_energies"]) == 45
    assert d["n_init_snapshots"] == 5


def test_isotherm_and_widom():
    c = collect_raspa_isotherm([D])
    assert list(c["pressure_pa"]) == [1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]    # Widom run excluded
    assert np.all(np.diff(c["loading_mol_kg"]) > 0)
    kh_widom = c["widom"][0]["henry_coefficient"]
    err = c["widom"][0]["henry_coefficient_err"]
    assert kh_widom == pytest.approx(6.10841e-06)
    kh = calculate_henry_constant(c["pressure_pa"], c["loading_mol_kg"])["henry_constant"]
    assert abs(kh - kh_widom) < 3 * err                       # v1.0.0 method: 4.5% (6 SE) low
    old = np.linalg.lstsq(c["pressure_pa"][:4, None], c["loading_mol_kg"][:4], rcond=None)[0][0]
    assert abs(old - kh_widom) > 3 * err


@pytest.mark.parametrize("name", ["p_300000", "p_1000000"])
def test_fluctuation_qst_matches_raspa(name):
    d = _one(name)
    k = d["n_init_snapshots"]
    q = calculate_isosteric_heat_fluctuations(d["cycle_energies"][k:], d["cycle_molecules"][k:], 300.0, "K")
    # 40 snapshots vs RASPA's average over all 40 000 cycles
    assert q == pytest.approx(d["isosteric_heat_kj_mol"], abs=0.5)


def test_no_false_drift_on_stationary_runs():
    for name in ["p_1000", "p_3000", "p_10000", "p_30000", "p_100000", "p_300000", "p_1000000"]:
        b = detect_gcmc_burnin(_one(name)["cycle_loadings"])
        assert b.status != "FAIL", (name, b.diagnostic_message)   # v1.0.0: 5 of 7 FAIL
