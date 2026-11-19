"""
Worked example: from a folder of RASPA3 runs to a quality-control report.

Uses the RASPA3 results for methane in MFI bundled with the repository
(tests/data/raspa3_mfi_ch4: one run per pressure plus a Widom insertion run in `henry/`).

Run from the repository root:

    python examples/raspa3_to_report.py

The same analysis is available from the command line:

    adsorpqc assess -i tests/data/raspa3_mfi_ch4 -o raspa3_example_output
"""

import os

from adsorpqc import assess_adsorption_quality, detect_gcmc_burnin
from adsorpqc.parsers.raspa import collect_raspa_isotherm
from adsorpqc.reporters import (
    generate_adsorption_figures,
    generate_adsorption_html_report,
    generate_adsorption_manuscript_assets,
)

DATA_DIR = os.path.join("tests", "data", "raspa3_mfi_ch4")
OUTPUT_DIR = "raspa3_example_output"


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Read every RASPA run in the folder and build the isotherm.
    iso = collect_raspa_isotherm([DATA_DIR])
    pressure_bar = iso["pressure_pa"] / 1e5
    loading = iso["loading_mol_kg"]
    print(f"Runs read: {len(iso['outputs'])} GCMC, {len(iso['widom'])} Widom")

    # 2. Burn-in and drift check for each pressure.
    print(f"\n{'P (Pa)':>12} {'loading (mol/kg)':>18} {'burn-in':>8} {'N_eff':>8}  status")
    for run in iso["outputs"]:
        series = run["cycle_loadings"]
        if series is None:
            continue
        burn = detect_gcmc_burnin(series)
        print(f"{run['metadata']['pressure_pa']:>12.0f} {run['loading_mol_kg']:>18.4f} "
              f"{burn.n_burnin_cycles:>8d} {burn.n_eff:>8.0f}  {burn.status}")

    # 3. Isotherm model selection, Henry coefficient and overall status.
    meta = dict(iso["outputs"][0]["metadata"])
    meta["henry_unit"] = "mol/kg/bar"
    report = assess_adsorption_quality(
        metadata=meta,
        pressure_isotherm=pressure_bar,
        loading_isotherm=loading,
        temperature_k=meta.get("temperature_k") or 298.15,
    )
    best = report.isotherm_fits["best_model"]
    print(f"\nSelected model: {best.model_name} (R^2 = {best.r_squared:.4f}, RMSE = {best.rmse:.4f})")

    # 4. Compare the Henry coefficient with Widom insertion (both in mol/kg/Pa).
    if iso["widom"] and report.henry_data:
        k_iso = report.henry_data["henry_constant"] / 1e5   # mol/kg/bar -> mol/kg/Pa
        k_widom = iso["widom"][0]["henry_coefficient"]
        k_err = iso["widom"][0]["henry_coefficient_err"]
        print(f"Henry coefficient: isotherm {k_iso:.4e} vs Widom {k_widom:.4e} +/- {k_err:.1e} mol/kg/Pa "
              f"({100 * (k_iso - k_widom) / k_widom:+.2f} %)")

    print(f"Overall status: {report.overall_status}")

    # 5. Write the report, figures, tables and methods paragraph.
    generate_adsorption_figures(report, OUTPUT_DIR, pressure_raw=pressure_bar, loading_raw=loading)
    generate_adsorption_manuscript_assets(report, OUTPUT_DIR)
    generate_adsorption_html_report(report, os.path.join(OUTPUT_DIR, "report.html"))
    print(f"\nReport written to {os.path.join(OUTPUT_DIR, 'report.html')}")


if __name__ == "__main__":
    main()
