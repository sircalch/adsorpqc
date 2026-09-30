"""
IAST predictions from single-component GCMC isotherms against binary-mixture GCMC (RASPA3), for
CO2/CH4 in silicalite-1 at 300 K, and the Henry-limit ratio of AdsorpQC 1.0.0 for comparison.

Usage: python iast_vs_gcmc.py CH4_RUNS_DIR CO2_MIX_RUNS_DIR OUT_DIR
  CH4_RUNS_DIR      output of raspa3_mfi_methane.py (p_<Pa>/)
  CO2_MIX_RUNS_DIR  output of raspa3_co2_ch4_mfi.py (co2_p_<Pa>/, co2_henry/, mix_*/)
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from adsorpqc.parsers.raspa import parse_raspa_output, collect_raspa_isotherm       # noqa: E402
from adsorpqc.core.isotherms import fit_all_isotherm_models                         # noqa: E402
from adsorpqc.core.iast import iast_binary, calculate_iast_selectivity              # noqa: E402
from adsorpqc.core.energetics import calculate_henry_constant                       # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from legacy_v100 import iast_v100                                                    # noqa: E402


def isotherm(dirs):
    iso = collect_raspa_isotherm(dirs)
    return iso["pressure_pa"], iso["loading_mol_kg"], iso["loading_err"], iso


def main(ch4_dir, co2_dir, out):
    os.makedirs(out, exist_ok=True)
    p_ch4, q_ch4, e_ch4, _ = isotherm(glob.glob(os.path.join(ch4_dir, "p_*")))
    p_co2, q_co2, e_co2, iso_co2 = isotherm(glob.glob(os.path.join(co2_dir, "co2_*")))
    fits = {}
    rows_iso = []
    for name, p, q, e in (("CH4", p_ch4, q_ch4, e_ch4), ("CO2", p_co2, q_co2, e_co2)):
        r = fit_all_isotherm_models(p / 1e5, q)                   # pressures in bar for the fits
        fits[name] = r["best_model"]
        kh = calculate_henry_constant(p, q)["henry_constant"]
        for pp, qq, ee in zip(p, q, e):
            rows_iso.append({"component": name, "p_pa": pp, "q_mol_kg": qq, "err": ee,
                             "fit": float(r["best_model"].fitted_function(pp / 1e5))})
        print(name, "best model", r["best_model"].model_name, {k: round(v, 4) for k, v in r["best_model"].parameters.items()},
              "R2", round(r["best_model"].r_squared, 5), "K_H isotherm", kh)
    widom = [d for d in iso_co2["widom"]]
    if widom:
        print("CO2 Widom K_H", widom[0]["henry_coefficient"], "+/-", widom[0]["henry_coefficient_err"])
    pd.DataFrame(rows_iso).to_csv(os.path.join(out, "pure_isotherms.csv"), index=False)

    rows = []
    for d in sorted(glob.glob(os.path.join(co2_dir, "mix_*"))):
        m = re.match(r"mix_y([\d.]+)_p_(\d+)", os.path.basename(d))
        y, P = float(m.group(1)), float(m.group(2))
        out_txt = glob.glob(os.path.join(d, "output", "*.txt"))[0]
        g = {c["name"]: c for c in parse_raspa_output(out_txt)["components"]}
        sol = iast_binary(fits["CO2"].fitted_function, fits["CH4"].fitted_function, y, P / 1e5)
        full = calculate_iast_selectivity(p_co2 / 1e5, q_co2, p_ch4 / 1e5, q_ch4, gas_mole_fraction_a=y,
                                          gas_mole_fraction_b=1 - y, total_pressure=P / 1e5, n_bootstrap=300)
        v100 = iast_v100.calculate_iast_selectivity(p_co2 / 1e5, q_co2, p_ch4 / 1e5, q_ch4, gas_mole_fraction_a=y,
                                                    gas_mole_fraction_b=1 - y, total_pressure=P / 1e5, n_bootstrap=50)
        q1, q2 = g["CO2"]["loading_mol_kg"], g["CH4"]["loading_mol_kg"]
        e1, e2 = g["CO2"]["loading_mol_kg_err"], g["CH4"]["loading_mol_kg_err"]
        s_gcmc = (q1 / q2) / (y / (1 - y))
        s_err = s_gcmc * np.sqrt((e1 / q1) ** 2 + (e2 / q2) ** 2)
        rows.append({"y_co2": y, "p_pa": P, "gcmc_q_co2": q1, "gcmc_q_co2_err": e1, "gcmc_q_ch4": q2,
                     "gcmc_q_ch4_err": e2, "iast_q_co2": sol["q1_mix"], "iast_q_ch4": sol["q2_mix"],
                     "gcmc_selectivity": s_gcmc, "gcmc_selectivity_err": s_err,
                     "iast_selectivity": sol["selectivity_12"], "iast_ci_low": full["ci_lower_95"],
                     "iast_ci_high": full["ci_upper_95"], "v100_henry_ratio": full["henry_limit_selectivity"], "v100_output": v100["selectivity_a_b"],
                     "dev_q_co2_pct": 100 * (sol["q1_mix"] / q1 - 1), "dev_q_ch4_pct": 100 * (sol["q2_mix"] / q2 - 1)})
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(out, "iast_vs_gcmc.csv"), index=False)
    print(d.round(4).to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:4])
