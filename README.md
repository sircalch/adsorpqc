# AdsorpQC

[![CI](https://github.com/sircalch/adsorpqc/actions/workflows/test.yml/badge.svg)](https://github.com/sircalch/adsorpqc/actions)
[![PyPI version](https://img.shields.io/pypi/v/adsorpqc.svg?color=blue)](https://pypi.org/project/adsorpqc/)
[![Python versions](https://img.shields.io/pypi/pyversions/adsorpqc.svg)](https://pypi.org/project/adsorpqc/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22217574.svg)](https://doi.org/10.5281/zenodo.22217574)

> **Automated Quality-Control, GCMC Burn-in Detection, Isotherm Model Fitting, and Reproducibility Toolkit for Adsorption Simulations in Nanoporous Materials.**

---

## Overview

**AdsorpQC** is an open-source scientific software toolkit tailored for the rigorous validation and certification of molecular adsorption simulations (**Grand Canonical Monte Carlo, GCMC**, RASPA, LAMMPS) and experimental gas isotherms in **MOFs, COFs, Zeolites, and Porous Carbons**.

In GCMC simulations, reporting average loadings before reaching true thermodynamic equilibrium (*under-equilibrated burn-in*) or during unresolved systematic drift leads to invalid capacities. `adsorpqc` addresses this with an automated quality-control pipeline in a single command:

- ⏳ **Automated GCMC Burn-in Detection**:
  - Maximizes the number of statistically independent production samples ($N_{\text{eff}}$).
  - Automatically identifies and truncates un-equilibrated initialization cycles.
- 📉 **Loading Drift & Stationarity Audit**:
  - Compares the first and last 20% of the production phase. The drift is flagged only when it is statistically significant against the autocorrelation-corrected noise ($z > 2$ / $z > 3$) and large relative to the mean.
  - Generates clear safety alerts: `"DO NOT REPORT PRODUCTION AVERAGE"` if the drift is significant ($z > 3$) and larger than 8% of the mean.
- 🧪 **Non-Linear Isotherm Fitting & Model Selection**:
  - Fits **Langmuir**, **Dual-Site Langmuir (DSLG)**, **Sips (Langmuir-Freundlich)** and **Toth** models.
  - Selects the model with the lowest small-sample corrected Akaike Information Criterion (**AICc**), and reports $R^2$ and RMSE.
- ⚡ **Thermodynamic Energetics & Selectivity**:
  - **Henry coefficient ($K_H$)**, obtained by extrapolating $q/P$ to zero loading. On RASPA3 methane/MFI data it agrees with Widom insertion within 0.1%.
  - **Isosteric heat of adsorption ($q_{\text{st}}$)**, from the $\mu VT$ particle–energy fluctuation formula. It needs the number of molecules and the energy with explicit units; RASPA3's own value is read when available.
  - **Ideal Adsorbed Solution Theory (IAST)** binary mixture selectivity ($S_{A/B}$) with **95% Bootstrap Confidence Intervals**. This part has not yet been validated against an independent IAST implementation.
- 📑 **Publication Deliverables**:
  - Interactive self-contained `report.html` dashboard.
  - Vector publication plots (Isotherm fit curves, GCMC burn-in shaded trajectory) in SVG, PDF, PNG (300 DPI).
  - Ready-to-compile LaTeX summary tables (`.tex`).
  - Draft **Methods** text paragraph and BibTeX citation (`citation.bib`).

```
    GCMC Simulation (RASPA) / Isotherm CSV
                      │
                      ▼
  ┌───────────────────────────────────────────────────────────┐
  │                         AdsorpQC                          │
  │  ├── Automated GCMC Burn-in / Equilibration Cutoff        │
  │  ├── Production Loading Drift & Stationarity Audit        │
  │  ├── Non-Linear Isotherm Fitting (Langmuir, Sips, Toth)   │
  │  ├── Henry Constant (K_H) & Isosteric Heat (q_st)        │
  │  └── IAST Selectivity (S_A/B) with Bootstrap 95% CI       │
  └───────────────────────────────────────────────────────────┘
                      │
                      ▼
  ┌───────────────────────────────────────────────────────────┐
  │                   Publication Deliverables                │
  │  ├── report.html (Interactive Dashboard & Badges)         │
  │  ├── adsorpqc_isotherm_fit.pdf/svg/png                    │
  │  ├── adsorpqc_gcmc_burnin_trajectory.pdf/svg/png         │
  │  ├── adsorpqc_summary_table.tex / .csv                    │
  │  ├── methods_snippet.txt (Ready for Manuscript)           │
  │  └── citation.bib (BibTeX Reference)                      │
  └───────────────────────────────────────────────────────────┘
```

---

## Installation

### From PyPI
```bash
pip install adsorpqc
```

### From Source
```bash
git clone https://github.com/sircalch/adsorpqc.git
cd adsorpqc
pip install -e .[dev]
```

---

## Quickstart (CLI)

### 1. Run Benchmark Demo (Instant GCMC + Isotherm Simulation)
```bash
adsorpqc demo -o my_adsorption_audit/
```
Open `my_adsorption_audit/report.html` in your browser to explore the interactive dashboard!

### 2. Assess Isotherm CSV File
```bash
adsorpqc assess -i my_isotherm.csv --framework "Mg-MOF-74" --adsorbate "CO2" --temp 298.15 -o output_report/
```

### 3. Assess RASPA Simulation Output
```bash
# one RASPA3 run (output/output_<T>_<P>.s0.txt) or RASPA2 run (output_*.data)
adsorpqc assess -i output/output_300_1e+05.s0.txt -o raspa_audit/
# a folder of RASPA runs, one pressure each; Widom runs are compared with the isotherm K_H
adsorpqc assess -i runs/ -o raspa_isotherm/
```

---

## Python API Usage

```python
import numpy as np
from adsorpqc import assess_adsorption_quality
from adsorpqc.parsers import parse_isotherm_csv
from adsorpqc.reporters import generate_adsorption_figures, generate_adsorption_manuscript_assets, generate_adsorption_html_report

# 1. Parse isotherm data
data = parse_isotherm_csv("co2_isotherm.csv")

# 2. Assess quality and fit models
report = assess_adsorption_quality(
    metadata={"framework": "HKUST-1", "adsorbate": "CO2", "temperature_k": 298.15},
    pressure_isotherm=data["pressure"],
    loading_isotherm=data["loading"],
    temperature_k=298.15
)

print(f"Overall Certification: {report.overall_status}")
print(f"Optimal Model: {report.isotherm_fits['best_model_name']} (R^2 = {report.isotherm_fits['best_r_squared']:.4f})")

# 3. Export all publication assets
generate_adsorption_figures(report, "output_dir/", pressure_raw=data["pressure"], loading_raw=data["loading"])
generate_adsorption_manuscript_assets(report, "output_dir/")
generate_adsorption_html_report(report, "output_dir/report.html")
```

---

## Citation

If you use AdsorpQC in your publications, please cite:

```bibtex
@software{monreal2026adsorpqc,
  author = {Monreal-Hern{\'a}ndez, Andre},
  title = {{AdsorpQC: An Open-Source Toolkit for Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and Reproducibility Assessment of Adsorption Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/adsorpqc}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

