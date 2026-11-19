# Validation

The `validation/` folder of the repository contains the scripts and the result files behind these comparisons, so they can be repeated.

| Quantity | Reference | Result | Files |
|---|---|---|---|
| Henry coefficient | Widom insertion, RASPA3, MFI | Agreement within 0.11 % (CH\(_4\)) and 0.16 % (CO\(_2\)), inside the Widom uncertainty of about 0.7 % | `raspa3_mfi_methane.py`, `raspa3_co2_ch4_mfi.py`, `make_figures.py`, `results/henry_summary.csv` |
| IAST | pyIAST | Relative difference below \(10^{-11}\) across the tested model combinations | `iast_vs_pyiast.py`, `results/iast_vs_pyiast.csv` |
| IAST | Explicit binary GCMC, CO\(_2\)/CH\(_4\) in MFI | Adsorbed amounts within 3.8 % | `iast_vs_gcmc.py`, `results/iast_vs_gcmc.csv` |
| Drift test | Synthetic series with known drift | False-flag and detection rates by regime and run length | `known_answer_studies.py`, `results/drift_study.csv` |
| Model selection | Synthetic isotherms from known models | Frequency with which each criterion recovers the generating model | `known_answer_studies.py`, `results/model_selection_study.csv` |

A detailed comparison with RASPA3 is available as a preprint (not peer reviewed):

> Monreal Hernández, A., Franco Amaya, S. L., & Martínez Osorio, C. I. (2026). *Checking the post-processing of adsorption simulations: Henry coefficients, isosteric heats, IAST and stationarity tests validated against RASPA3.* Research Square. <https://doi.org/10.21203/rs.3.rs-11626447/v1>
