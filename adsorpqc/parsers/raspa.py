"""
Parser for RASPA Grand Canonical Monte Carlo (GCMC) and Henry coefficient output files.
"""

from typing import Dict, Any, List, Optional
import os
import re
import numpy as np


def parse_raspa_output(filepath: str) -> Dict[str, Any]:
    """
    Parses a RASPA output file (output_*.data).

    Parameters
    ----------
    filepath : str
        Path to the RASPA output file.

    Returns
    -------
    data : dict
        Parsed simulation metadata, loadings, cycle counts, energy, and Henry constants.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    metadata: Dict[str, Any] = {
        "engine": "RASPA",
        "framework": "Unknown",
        "adsorbate": "Unknown",
        "temperature_k": 298.15,
        "pressure_bar": 1.0,
        "init_cycles": 0,
        "prod_cycles": 0
    }

    # 1. Framework & Component
    fw_match = re.search(r"Framework name:\s+([A-Za-z0-9\-_\.]+)", content)
    if fw_match:
        metadata["framework"] = fw_match.group(1)
        
    comp_match = re.search(r"Component\s+\d+\s+\[([A-Za-z0-9\-_\.]+)\]", content)
    if comp_match:
        metadata["adsorbate"] = comp_match.group(1)

    # 2. Temperature & Pressure
    t_match = re.search(r"Temperature:\s+([\d\.]+)\s+\[K\]", content)
    if t_match:
        metadata["temperature_k"] = float(t_match.group(1))
        
    p_match = re.search(r"Pressure:\s+([\d\.]+)\s+\[Pa\]", content)
    if p_match:
        metadata["pressure_bar"] = float(p_match.group(1)) / 1e5
    else:
        p_match_bar = re.search(r"Pressure:\s+([\d\.]+)\s+\[bar\]", content)
        if p_match_bar:
            metadata["pressure_bar"] = float(p_match_bar.group(1))

    # 3. Cycles
    init_match = re.search(r"Number of initialization cycles:\s+(\d+)", content)
    if init_match:
        metadata["init_cycles"] = int(init_match.group(1))
        
    prod_match = re.search(r"Number of production cycles:\s+(\d+)", content)
    if prod_match:
        metadata["prod_cycles"] = int(prod_match.group(1))

    # 4. Average Loading
    # Example: "Average loading absolute [mol/kg framework]   2.45120 +/- 0.05120"
    # Example: "Average loading absolute [milligram/gram framework]   107.85 +/- 2.25"
    # Example: "Average loading absolute [molecules/unit cell]   12.35 +/- 0.25"
    loading_mol_kg = None
    loading_molec_uc = None
    
    mol_kg_match = re.search(r"Average loading absolute \[mol/kg framework\]\s+([\d\.]+)\s+\+/-", content)
    if mol_kg_match:
        loading_mol_kg = float(mol_kg_match.group(1))
        
    molec_match = re.search(r"Average loading absolute \[molecules/unit cell\]\s+([\d\.]+)\s+\+/-", content)
    if molec_match:
        loading_molec_uc = float(molec_match.group(1))

    # 5. Henry coefficient (if Henry simulation)
    # "Henry coefficient: 0.00123456 [mol/kg/Pa]"
    henry_val = None
    henry_match = re.search(r"Henry coefficient:\s+([\d\.eE\-+]+)\s+\[mol/kg/Pa\]", content)
    if henry_match:
        henry_val = float(henry_match.group(1))

    # 6. Heat of adsorption
    # "Heat of adsorption: 28.54 +/- 0.45 [kJ/mol]"
    qst_val = None
    qst_match = re.search(r"Heat of adsorption:\s+([\d\.]+)\s+\+/-", content)
    if qst_match:
        qst_val = float(qst_match.group(1))

    # 7. Check if cycle trajectory is present
    # "Cycle 1: Loading = 1.25, Energy = -45.2"
    cycle_loadings = []
    cycle_energies = []
    lines = content.splitlines()
    for l in lines:
        if "Current loading" in l or "Loading per cycle" in l:
            parts = re.findall(r"[\-\d\.]+", l)
            if parts:
                cycle_loadings.append(float(parts[-1]))

    return {
        "metadata": metadata,
        "loading_mol_kg": loading_mol_kg,
        "loading_molec_uc": loading_molec_uc,
        "henry_coefficient": henry_val,
        "isosteric_heat_kj_mol": qst_val,
        "cycle_loadings": np.asarray(cycle_loadings) if cycle_loadings else None,
        "cycle_energies": np.asarray(cycle_energies) if cycle_energies else None
    }
