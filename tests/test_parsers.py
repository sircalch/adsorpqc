"""
Tests for RASPA output parser and generic CSV isotherm parser.
"""

import os
import tempfile
import pytest
import numpy as np
from adsorpqc.parsers.raspa import parse_raspa_output
from adsorpqc.parsers.generic_isotherm import parse_isotherm_csv


def test_raspa_parser():
    raspa_content = """
================================================================================
                                R A S P A - 2 . 0
================================================================================

Framework name: Mg-MOF-74
Component 0 [CO2]
Temperature: 298.15 [K]
Pressure: 100000.000000 [Pa]
Number of initialization cycles: 1000
Number of production cycles: 5000

Average loading absolute [mol/kg framework]   8.45120 +/- 0.05120
Average loading absolute [molecules/unit cell]   12.35000 +/- 0.25000
Heat of adsorption: 38.54 +/- 0.45 [kJ/mol]
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".data", delete=False) as f:
        f.write(raspa_content)
        f_path = f.name
        
    try:
        data = parse_raspa_output(f_path)
        assert data["metadata"]["framework"] == "Mg-MOF-74"
        assert data["metadata"]["adsorbate"] == "CO2"
        assert data["metadata"]["temperature_k"] == 298.15
        assert np.isclose(data["metadata"]["pressure_bar"], 1.0)
        assert np.isclose(data["loading_mol_kg"], 8.4512)
        assert np.isclose(data["isosteric_heat_kj_mol"], 38.54)
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)


def test_isotherm_csv_parser():
    csv_content = """Pressure_bar,Loading_mol_kg
0.01,0.52
0.05,2.15
0.10,3.80
0.50,6.90
1.00,8.20
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(csv_content)
        f_path = f.name
        
    try:
        data = parse_isotherm_csv(f_path)
        assert len(data["pressure"]) == 5
        assert len(data["loading"]) == 5
        assert data["pressure"][0] == 0.01
        assert data["loading"][-1] == 8.20
    finally:
        if os.path.exists(f_path):
            os.remove(f_path)
