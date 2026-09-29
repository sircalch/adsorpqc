"""
Cross-check of AdsorpQC's IAST solver against pyIAST (Simon, Smit & Haranczyk, Comput. Phys.
Commun. 200, 364, 2016) for model isotherms with known parameters, where IAST differs from extended
Langmuir: unequal saturation capacities, and dual-site Langmuir components.

Usage: python iast_vs_pyiast.py OUT_CSV
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from adsorpqc.core.iast import iast_binary      # noqa: E402

import pyiast                                    # noqa: E402

P_GRID = np.logspace(-3, 2, 60)                   # bar


def langmuir(qs, k):
    return lambda p: qs * k * p / (1 + k * p)


def dsl(q1, k1, q2, k2):
    return lambda p: q1 * k1 * p / (1 + k1 * p) + q2 * k2 * p / (1 + k2 * p)


CASES = {
    "Langmuir, equal q_sat": ((langmuir(5.0, 5.0), {"M": 5.0, "K": 5.0}, "Langmuir"),
                              (langmuir(5.0, 0.5), {"M": 5.0, "K": 0.5}, "Langmuir")),
    "Langmuir, q_sat 2 vs 6": ((langmuir(2.0, 8.0), {"M": 2.0, "K": 8.0}, "Langmuir"),
                               (langmuir(6.0, 0.3), {"M": 6.0, "K": 0.3}, "Langmuir")),
    "dual-site Langmuir": ((dsl(1.5, 20.0, 2.0, 0.5), {"M1": 1.5, "K1": 20.0, "M2": 2.0, "K2": 0.5}, "DSLangmuir"),
                           (dsl(3.0, 0.2, 1.0, 2.0), {"M1": 3.0, "K1": 0.2, "M2": 1.0, "K2": 2.0}, "DSLangmuir")),
}


def pyiast_model(fn, params, name):
    df = pd.DataFrame({"P": P_GRID, "q": fn(P_GRID)})
    iso = pyiast.ModelIsotherm(df, loading_key="q", pressure_key="P", model=name, param_guess=params,
                               optimization_method="Nelder-Mead")
    iso.params = dict(params)                      # use the exact parameters, not a refit
    return iso


def main(out):
    rows = []
    for case, ((f1, p1, m1), (f2, p2, m2)) in CASES.items():
        i1, i2 = pyiast_model(f1, p1, m1), pyiast_model(f2, p2, m2)
        for y1 in (0.1, 0.5, 0.9):
            for P in (0.01, 0.1, 1.0, 10.0):
                ours = iast_binary(f1, f2, y1, P)
                ref = pyiast.iast(np.array([P * y1, P * (1 - y1)]), [i1, i2], verboseflag=False)
                rows.append({"case": case, "y1": y1, "P_bar": P, "q1_ours": ours["q1_mix"], "q2_ours": ours["q2_mix"],
                             "q1_pyiast": ref[0], "q2_pyiast": ref[1],
                             "rel_diff_q1": abs(ours["q1_mix"] - ref[0]) / ref[0],
                             "rel_diff_q2": abs(ours["q2_mix"] - ref[1]) / ref[1],
                             "selectivity": ours["selectivity_12"],
                             "henry_limit_selectivity": (f1(1e-9) / 1e-9) / (f2(1e-9) / 1e-9)})
    d = pd.DataFrame(rows)
    d.to_csv(out, index=False)
    print(d.groupby("case")[["rel_diff_q1", "rel_diff_q2"]].max())
    print(d[["case", "y1", "P_bar", "selectivity", "henry_limit_selectivity"]].to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
