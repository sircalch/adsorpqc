"""
Parsers for RASPA3 and RASPA2 Grand Canonical Monte Carlo (GCMC) and Widom (Henry coefficient) output.

RASPA3 (output/output_<T>_<P>.s0.txt) is validated against real RASPA3 3.0.29 runs. The RASPA2
patterns (Output/System_0/output_*.data) follow the documented RASPA2 output format but have not
been checked against real RASPA2 output.

Only the first component is read (single-component isotherms).
"""

from typing import Dict, Any, List, Optional
import glob
import os
import re
import numpy as np

_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
KB_KJ_PER_MOL_K = 0.00831446261815324


def _search(pattern: str, text: str, flags: int = 0) -> Optional[re.Match]:
    return re.search(pattern, text, flags)


def _parse_raspa3(content: str, meta: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    m = _search(r"^Framework \d+ \[([^\]]+)\]", content, re.M)
    if m:
        meta["framework"] = m.group(1)
    m = _search(r"^Component 0 \[([^\]]+)\]", content, re.M)
    if m:
        meta["adsorbate"] = m.group(1)
    m = _search(rf"^Temperature:\s+({_NUM}) \[K\]", content, re.M)
    if m:
        meta["temperature_k"] = float(m.group(1))
    m = _search(rf"^Pressure:\s+({_NUM}) \[Pa\]", content, re.M)
    if m:
        meta["pressure_pa"] = float(m.group(1))
        meta["pressure_bar"] = float(m.group(1)) / 1e5
    m = _search(r"Initialization: Current cycle: \d+ out of (\d+)", content)
    if m:
        meta["init_cycles"] = int(m.group(1))
    m = _search(r"^Current cycle: \d+ out of (\d+)", content, re.M)
    if m:
        meta["prod_cycles"] = int(m.group(1))

    # Block-averaged loadings (final summary)
    for kind, key in (("Abs.", "loading"), ("Excess", "excess_loading")):
        for unit, suffix in (("mol/kg-framework", "mol_kg"), ("molecules/uc", "molec_uc"), ("mg/g-framework", "mg_g")):
            m = _search(rf"{re.escape(kind)} loading average\s+({_NUM}) \+/-\s+({_NUM}) \[{re.escape(unit)}\]", content)
            if m:
                out[f"{key}_{suffix}"] = float(m.group(1))
                out[f"{key}_{suffix}_err"] = float(m.group(2))

    # Enthalpy of adsorption (GCMC fluctuation formula). RASPA3 prints Delta H (negative);
    # the isosteric heat is Q_st = -Delta H.
    m = _search(rf"Enthalpy of adsorption:\s+{_NUM} \+/-\s+{_NUM} \[K\]\s*\n\s*({_NUM}) \+/-\s+({_NUM}) \[kJ/mol\]", content)
    if m:
        out["enthalpy_kj_mol"] = float(m.group(1))
        out["enthalpy_kj_mol_err"] = float(m.group(2))

    # Widom insertion: Henry coefficient
    m = _search(rf"Average Henry coefficient:\s+({_NUM}) \+/-\s+({_NUM}) \[mol/kg/Pa\]", content)
    if m:
        out["henry_mol_kg_pa"] = float(m.group(1))
        out["henry_mol_kg_pa_err"] = float(m.group(2))

    # Snapshot series written every PrintEvery cycles (initialization + production). Each block
    # starts at a "Current cycle:" line; the three quantities are read inside the same block so the
    # series stay aligned (the initialization and production blocks are formatted differently).
    marks = [(mm.start(), mm.group(1) is not None) for mm in
             re.finditer(r"^(Initialization: )?Current cycle: \d+ out of \d+", content, re.M)]
    final = content.find("Final state after")
    loads, mols, ens, n_init = [], [], [], 0
    for k, (start, is_init) in enumerate(marks):
        end = marks[k + 1][0] if k + 1 < len(marks) else (final if final > start else len(content))
        blk = content[start:end]
        m_n = re.search(rf"absolute adsorption:\s+({_NUM}) molecules", blk)
        m_q = re.search(rf"absolute adsorption:(?:.*\n){{2}}\s+({_NUM}) mol/kg", blk)
        m_u = re.search(rf"^Total potential energy/k\S*\s+({_NUM})", blk, re.M)
        if not (m_n and m_q and m_u):
            continue
        mols.append(float(m_n.group(1)))
        loads.append(float(m_q.group(1)))
        ens.append(float(m_u.group(1)))
        n_init += int(is_init)
    out["cycle_loadings"] = np.array(loads) if loads else None          # mol/kg
    out["cycle_molecules"] = np.array(mols) if mols else None           # molecules in the box
    out["cycle_energies"] = np.array(ens) if ens else None              # U/k_B [K]
    out["n_init_snapshots"] = n_init
    out["n_prod_snapshots"] = len(loads) - n_init
    return out


def _parse_raspa2(content: str, meta: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    m = _search(r"Framework name:\s+(\S+)", content)
    if m:
        meta["framework"] = m.group(1)
    m = _search(r"Component\s+0\s+\[([^\]]+)\]", content)
    if m:
        meta["adsorbate"] = m.group(1)
    m = _search(rf"External temperature:\s+({_NUM}) \[K\]", content) or _search(rf"Temperature:\s+({_NUM}) \[K\]", content)
    if m:
        meta["temperature_k"] = float(m.group(1))
    m = _search(rf"External Pressure:\s+({_NUM}) \[Pa\]", content) or _search(rf"Pressure:\s+({_NUM}) \[Pa\]", content)
    if m:
        meta["pressure_pa"] = float(m.group(1))
        meta["pressure_bar"] = float(m.group(1)) / 1e5
    m = _search(r"Number of initialization cycles:\s+(\d+)", content)
    if m:
        meta["init_cycles"] = int(m.group(1))
    m = _search(r"Number of cycles:\s+(\d+)", content) or _search(r"Number of production cycles:\s+(\d+)", content)
    if m:
        meta["prod_cycles"] = int(m.group(1))
    for kind, key in (("absolute", "loading"), ("excess", "excess_loading")):
        for unit, suffix in (("mol/kg framework", "mol_kg"), ("molecules/unit cell", "molec_uc"), ("milligram/gram framework", "mg_g")):
            m = _search(rf"Average loading {kind} \[{re.escape(unit)}\]\s+({_NUM}) \+/-\s+({_NUM})", content)
            if m:
                out[f"{key}_{suffix}"] = float(m.group(1))
                out[f"{key}_{suffix}_err"] = float(m.group(2))
    m = _search(rf"\[[^\]]+\] Average Henry coefficient:\s+({_NUM}) \+/-\s+({_NUM}) \[mol/kg/Pa\]", content)
    if m:
        out["henry_mol_kg_pa"] = float(m.group(1))
        out["henry_mol_kg_pa_err"] = float(m.group(2))
    m = _search(rf"Enthalpy of adsorption:.*?\n.*?({_NUM}) \+/-\s+({_NUM}) \[KJ/MOL\]", content, re.S | re.I)
    if m:
        out["enthalpy_kj_mol"] = float(m.group(1))
        out["enthalpy_kj_mol_err"] = float(m.group(2))
    else:
        # legacy one-line summary "Heat of adsorption: Q +/- dQ [kJ/mol]" (Q_st, positive)
        m = _search(rf"Heat of adsorption:\s+({_NUM}) \+/-\s+({_NUM}) \[kJ/mol\]", content)
        if m:
            out["enthalpy_kj_mol"] = -float(m.group(1))
            out["enthalpy_kj_mol_err"] = float(m.group(2))
    out["cycle_loadings"] = None
    out["cycle_energies"] = None
    return out


def parse_raspa_output(filepath: str) -> Dict[str, Any]:
    """
    Parses a RASPA3 (output/*.s0.txt) or RASPA2 (output_*.data) output file.

    Returns
    -------
    data : dict
        metadata (framework, adsorbate, temperature_k, pressure_pa, pressure_bar, cycles, engine),
        loading_mol_kg (+ _err), loading_molec_uc, excess loadings, henry_coefficient [mol/kg/Pa],
        enthalpy_kj_mol (Delta H, negative), isosteric_heat_kj_mol (Q_st = -Delta H),
        cycle_loadings [mol/kg], cycle_molecules [molecules in the box] and cycle_energies [U/k_B, K]:
        RASPA3 snapshot series written every PrintEvery cycles (initialization + production).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    meta: Dict[str, Any] = {"framework": "Unknown", "adsorbate": "Unknown", "temperature_k": None,
                            "pressure_pa": None, "pressure_bar": None, "init_cycles": None, "prod_cycles": None}
    is_raspa3 = "LoadingData" in content or "Abs. loading average" in content or "Widom insertion" in content
    if is_raspa3:
        meta["engine"] = "RASPA3"
        d = _parse_raspa3(content, meta)
    else:
        meta["engine"] = "RASPA2"
        d = _parse_raspa2(content, meta)

    enthalpy = d.get("enthalpy_kj_mol")
    return {
        "metadata": meta,
        "loading_mol_kg": d.get("loading_mol_kg"),
        "loading_mol_kg_err": d.get("loading_mol_kg_err"),
        "loading_molec_uc": d.get("loading_molec_uc"),
        "excess_loading_mol_kg": d.get("excess_loading_mol_kg"),
        "henry_coefficient": d.get("henry_mol_kg_pa"),
        "henry_coefficient_err": d.get("henry_mol_kg_pa_err"),
        "enthalpy_kj_mol": enthalpy,
        "enthalpy_kj_mol_err": d.get("enthalpy_kj_mol_err"),
        "isosteric_heat_kj_mol": -enthalpy if enthalpy is not None else None,
        "cycle_loadings": d.get("cycle_loadings"),
        "cycle_energies": d.get("cycle_energies"),
        "cycle_molecules": d.get("cycle_molecules"),
        "n_init_snapshots": d.get("n_init_snapshots"),
    }


def collect_raspa_isotherm(paths: List[str]) -> Dict[str, Any]:
    """
    Builds an isotherm from several RASPA outputs (one pressure each). `paths` may contain files,
    directories (searched recursively for RASPA3/RASPA2 outputs) or glob patterns.

    Returns pressures [Pa], absolute loadings [mol/kg] and their block errors, sorted by pressure,
    plus the list of parsed outputs.
    """
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files += glob.glob(os.path.join(p, "**", "output_*.s0.txt"), recursive=True)
            files += glob.glob(os.path.join(p, "**", "output_*.data"), recursive=True)
        else:
            files += glob.glob(p) or [p]
    parsed = [parse_raspa_output(f) for f in sorted(set(files))]
    points = [(d["metadata"]["pressure_pa"], d["loading_mol_kg"], d["loading_mol_kg_err"] or 0.0, d)
              for d in parsed if d["metadata"]["pressure_pa"] and d["loading_mol_kg"] is not None
              and d["henry_coefficient"] is None]          # Widom-only runs (P = 0) are not isotherm points
    points.sort(key=lambda t: t[0])
    return {
        "pressure_pa": np.array([p[0] for p in points]),
        "loading_mol_kg": np.array([p[1] for p in points]),
        "loading_err": np.array([p[2] for p in points]),
        "outputs": [p[3] for p in points],
        "widom": [d for d in parsed if d["henry_coefficient"] is not None],
    }
