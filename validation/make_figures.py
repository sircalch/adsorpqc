"""
Figures and tables for the AdsorpQC v2 manuscript, built only from validation/results/ and the
archived RASPA3 outputs (tests/data/raspa3_mfi_ch4, validation/results/raspa3_co2_ch4_outputs.zip).

    python validation/make_figures.py

Style: Springer double-column width 174 mm, 8 pt sans-serif text; one fixed colour and marker per
series in every figure (validated categorical palette; identity never relies on colour alone).
"""
import glob
import os
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from adsorpqc.parsers.raspa import collect_raspa_isotherm                    # noqa: E402
from adsorpqc.core.isotherms import fit_all_isotherm_models                  # noqa: E402
from adsorpqc.core.energetics import calculate_henry_constant                # noqa: E402
from adsorpqc.core.iast import iast_binary                                   # noqa: E402
from legacy_v100 import energetics_v100                                      # noqa: E402

RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
MM = 1 / 25.4
DOUBLE = 174 * MM
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C = {"CO2": ("#2a78d6", "o"), "CH4": ("#eb6834", "s"), "new": "#1baf7a", "old": "#e87ba4", "gray": "#8a8984"}


def setup():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2, "lines.markersize": 4,
        "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42})


def panel(ax, letter):
    ax.text(-0.15, 1.03, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def load_isotherms():
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(os.path.join(RES, "raspa3_co2_ch4_outputs.zip")) as z:
        z.extractall(tmp)
    ch4 = collect_raspa_isotherm([os.path.join(ROOT, "tests", "data", "raspa3_mfi_ch4")])
    co2 = collect_raspa_isotherm(glob.glob(os.path.join(tmp, "co2_*")))
    return {"CH4": ch4, "CO2": co2}


def fig_pure(iso):
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 64 * MM), gridspec_kw={"width_ratios": [1.25, 1]})
    rows = []
    for gas in ("CO2", "CH4"):
        col, mk = C[gas]
        d = iso[gas]
        P, q, e = d["pressure_pa"], d["loading_mol_kg"], d["loading_err"]
        fit = fit_all_isotherm_models(P / 1e5, q)["best_model"]
        kw = d["widom"][0]["henry_coefficient"]
        ke = d["widom"][0]["henry_coefficient_err"]
        grid = np.logspace(2.5, 6.05, 200)
        a.errorbar(P, q, yerr=1.96 * e, fmt=mk, color=col, ms=4.5, mec="white", mew=0.5, capsize=1.5, elinewidth=0.6, zorder=3)
        a.plot(grid, fit.fitted_function(grid / 1e5), color=col, lw=1.0, zorder=2)
        a.plot(grid[grid < 3e4], kw * grid[grid < 3e4], color=col, lw=0.8, ls=":", zorder=2)
        lab = {"CO2": r"CO$_2$", "CH4": r"CH$_4$"}[gas]
        a.text(P[-1] * 1.35, q[-1], lab, color=col, va="center", ha="left", fontsize=7)
        new = calculate_henry_constant(P, q)["henry_constant"]
        old = energetics_v100.calculate_henry_constant(P, q)["henry_constant"]
        rows.append({"gas": gas, "widom": kw, "widom_err": ke, "new": new, "old": old, "model": fit.model_name})
    a.set(xscale="log", yscale="log", xlabel="pressure (Pa)", ylabel="absolute loading (mol kg$^{-1}$)",
          xlim=(5e2, 4e6), ylim=(3e-3, 5))
    a.legend(handles=[Line2D([], [], color=INK2, marker="o", ls="none", ms=4, label="GCMC (RASPA3), 95% CI"),
                      Line2D([], [], color=INK2, lw=1.0, label="best fit (AICc)"),
                      Line2D([], [], color=INK2, lw=0.8, ls=":", label="Widom Henry law")], loc="upper left")
    panel(a, "a")

    h = pd.DataFrame(rows)
    x = np.arange(len(h))
    w = 0.34
    dev_new = 100 * (h.new / h.widom - 1)
    dev_old = 100 * (h.old / h.widom - 1)
    b.bar(x - w / 2, dev_new, w * 0.92, color=C["new"], label="AdsorpQC 1.1 (virial extrapolation)")
    b.bar(x + w / 2, dev_old, w * 0.92, color=C["old"], label="AdsorpQC 1.0.0 (line through origin)")
    for xi, v in zip(x - w / 2, dev_new):
        b.text(xi, v - 1.2 if v < 0 else v + 0.4, f"{v:+.2f}%", ha="center", va="top" if v < 0 else "bottom", fontsize=6.5)
    for xi, v in zip(x + w / 2, dev_old):
        if v < -10:
            b.text(xi, v / 2, f"{v:+.1f}%", ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")
        else:
            b.text(xi, v - 0.6, f"{v:+.1f}%", ha="center", va="top", fontsize=6.5)
    b.axhline(0, color=INK2, lw=0.6)
    b.set_xticks(x, [r"CO$_2$" if g == "CO2" else r"CH$_4$" for g in h.gas])
    b.set(ylabel="Henry coefficient vs Widom (%)", ylim=(-27, 4))
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig2_pure")
    h.to_csv(os.path.join(RES, "henry_summary.csv"), index=False)
    return h


def fig_iast(iso):
    mix = pd.read_csv(os.path.join(RES, "iast_vs_gcmc.csv"))
    fits = {g: fit_all_isotherm_models(iso[g]["pressure_pa"] / 1e5, iso[g]["loading_mol_kg"])["best_model"] for g in iso}
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 66 * MM))
    lim = (0.02, 3.5)
    a.plot(lim, lim, color=INK2, lw=0.7, ls="--", zorder=1)
    for gas, gcol, icol, ecol in (("CO2", "gcmc_q_co2", "iast_q_co2", "gcmc_q_co2_err"),
                                  ("CH4", "gcmc_q_ch4", "iast_q_ch4", "gcmc_q_ch4_err")):
        col, mk = C[gas]
        a.errorbar(mix[gcol], mix[icol], xerr=1.96 * mix[ecol], fmt=mk, color=col, ms=5, mec="white", mew=0.5,
                   capsize=1.5, elinewidth=0.6, label=r"CO$_2$" if gas == "CO2" else r"CH$_4$", zorder=3)
    a.set(xscale="log", yscale="log", xlim=lim, ylim=lim, xlabel="mixture GCMC loading (mol kg$^{-1}$)",
          ylabel="IAST prediction (mol kg$^{-1}$)")
    a.set_aspect("equal")
    a.legend(loc="upper left", title="component", title_fontsize=6.5)
    panel(a, "a")

    grid = np.logspace(3.5, 6.2, 40)
    s_iast = [iast_binary(fits["CO2"].fitted_function, fits["CH4"].fitted_function, 0.5, P / 1e5)["selectivity_12"] for P in grid]
    b.plot(grid, s_iast, color=C["new"], lw=1.2, label="IAST, AdsorpQC 1.1")
    b.axhline(mix.v100_output.iloc[0], color=C["old"], lw=1.2, ls="--", label="AdsorpQC 1.0.0 (Henry ratio)")
    e = mix[mix.y_co2 == 0.5]
    b.errorbar(e.p_pa, e.gcmc_selectivity, yerr=1.96 * e.gcmc_selectivity_err, fmt="D", color=INK, ms=4.5,
               mec="white", mew=0.5, capsize=1.5, elinewidth=0.6, label="mixture GCMC, $y_{\\mathrm{CO_2}}$ = 0.5", zorder=4)
    for _, r in e.iterrows():
        b.plot([r.p_pa, r.p_pa], [r.iast_ci_low, r.iast_ci_high], color=C["new"], lw=3, alpha=0.3, solid_capstyle="butt")
    b.set(xscale="log", xlabel="total pressure (Pa)", ylabel=r"CO$_2$/CH$_4$ selectivity", ylim=(5.0, 8.0))
    h, l = b.get_legend_handles_labels()
    h.insert(1, Line2D([], [], color=C["new"], lw=3, alpha=0.3))
    l.insert(1, "IAST 95% CI (bootstrap of the fits)")
    b.legend(h, l, loc="upper left")
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig3_iast")


def fig_known():
    dr = pd.read_csv(os.path.join(RES, "drift_study.csv"))
    ms = pd.read_csv(os.path.join(RES, "model_selection_study.csv"))
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 66 * MM), gridspec_kw={"width_ratios": [1.3, 1]})
    cases = [("high loading (AR(1), CV 5%)", 40, "high, 40"), ("high loading (AR(1), CV 5%)", 400, "high, 400"),
             ("low loading (Poisson, mean 2)", 40, "low, 40"), ("low loading (Poisson, mean 2)", 400, "low, 400")]
    x = np.arange(len(cases))
    w = 0.2
    for i, (ver, drift, col, lab) in enumerate((("1.0.0", 0.0, C["old"], "1.0.0, no drift (false FAIL)"),
                                                ("1.1", 0.0, C["new"], "1.1, no drift (false FAIL)"),
                                                ("1.0.0", 0.2, "#f3b5cf", "1.0.0, 20% drift (FAIL)"),
                                                ("1.1", 0.2, "#8fd9bb", "1.1, 20% drift (FAIL)"))):
        vals = [float(dr[(dr.regime == r) & (dr.n == n) & (dr.drift == drift) & (dr.version == ver)].fail_rate.iloc[0]) for r, n, _ in cases]
        a.bar(x + (i - 1.5) * w, vals, w * 0.92, color=col, label=lab)
    a.set_xticks(x, [c[2] for c in cases])
    a.set(xlabel="regime, samples per run", ylabel="fraction of runs failed", ylim=(0, 1.05))
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=2, columnspacing=0.8)
    panel(a, "a")

    lab = []
    vals = {"AIC": [], "AICc": []}
    for truth in ("Langmuir", "Dual-Site Langmuir"):
        for pts in (7, 12):
            lab.append(f"{'Langmuir' if truth == 'Langmuir' else 'dual-site'}\n{pts} points")
            for crit in vals:
                vals[crit].append(float(ms[(ms.truth == truth) & (ms.points == pts) & (ms.criterion == crit) & (ms.selected == truth)].fraction.iloc[0]))
    x = np.arange(len(lab))
    b.bar(x - 0.18, vals["AIC"], 0.34, color=C["gray"], label="AIC")
    b.bar(x + 0.18, vals["AICc"], 0.34, color=C["CO2"][0], label="AICc")
    b.set_xticks(x, lab)
    b.set(ylabel="true model selected (fraction)", ylim=(0, 1.05))
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig4_known_answer")


def table_iast():
    m = pd.read_csv(os.path.join(RES, "iast_vs_gcmc.csv")).sort_values(["y_co2", "p_pa"])
    lines = [r"\begin{tabular}{rrrrrrrrr}", r"\toprule",
             r"$y_{\mathrm{CO_2}}$ & $P$ (kPa) & \multicolumn{2}{c}{$q_{\mathrm{CO_2}}$ (mol kg$^{-1}$)} & \multicolumn{2}{c}{$q_{\mathrm{CH_4}}$ (mol kg$^{-1}$)} & \multicolumn{3}{c}{selectivity CO$_2$/CH$_4$} \\",
             r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-9}",
             r" & & GCMC & IAST & GCMC & IAST & GCMC & IAST [95\% CI] & 1.0.0 \\", r"\midrule"]
    for _, r in m.iterrows():
        lines.append(f"{r.y_co2:g} & {r.p_pa / 1e3:g} & {r.gcmc_q_co2:.3f}({round(r.gcmc_q_co2_err * 1e3):d}) & {r.iast_q_co2:.3f} & "
                     f"{r.gcmc_q_ch4:.4f}({round(r.gcmc_q_ch4_err * 1e4):d}) & {r.iast_q_ch4:.4f} & "
                     f"{r.gcmc_selectivity:.2f}({round(r.gcmc_selectivity_err * 100):d}) & "
                     f"{r.iast_selectivity:.2f} [{r.iast_ci_low:.2f}, {r.iast_ci_high:.2f}] & {r.v100_output:.2f} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_iast.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def table_pure(iso):
    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
             r"$P$ (kPa) & \multicolumn{2}{c}{CO$_2$} & \multicolumn{2}{c}{CH$_4$} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
             r" & $q$ (mol kg$^{-1}$) & $q_{\mathrm{st}}$ (kJ mol$^{-1}$) & $q$ (mol kg$^{-1}$) & $q_{\mathrm{st}}$ (kJ mol$^{-1}$) \\",
             r"\midrule"]
    by = {g: {o["metadata"]["pressure_pa"]: o for o in iso[g]["outputs"]} for g in iso}
    for P in sorted(by["CO2"]):
        cells = []
        for g in ("CO2", "CH4"):
            o = by[g][P]
            cells.append(f"{o['loading_mol_kg']:.4f} $\\pm$ {o['loading_mol_kg_err']:.4f} & {o['isosteric_heat_kj_mol']:.2f} $\\pm$ {o['enthalpy_kj_mol_err']:.2f}")
        lines.append(f"{P / 1e3:g} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_pure.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def si_tables():
    """Supplementary tables S2-S4 from the known-answer results."""
    p = pd.read_csv(os.path.join(RES, "iast_vs_pyiast.csv"))
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"Isotherm pair & $y_1$ & $P$ (bar) & $S_{12}$ (IAST) & $S_{12}$ (1.0.0) & max.\ rel.\ difference from pyIAST \\",
             r"\midrule"]
    for _, r in p.iterrows():
        case = r.case.replace("q_sat", r"$q_{\mathrm{sat}}$")
        lines.append(f"{case} & {r.y1:g} & {r.P_bar:g} & {r.selectivity:.3f} & {r.henry_limit_selectivity:.3f} & "
                     f"{max(r.rel_diff_q1, r.rel_diff_q2):.1e} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_s_pyiast.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    d = pd.read_csv(os.path.join(RES, "drift_study.csv"))
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"Regime & samples & drift & FAIL (1.0.0) & FAIL (1.1) & flagged (1.0.0) & flagged (1.1) \\",
             r"\midrule"]
    for (reg, n, dr), g in d.groupby(["regime", "n", "drift"], sort=False):
        o, w = g[g.version == "1.0.0"].iloc[0], g[g.version == "1.1"].iloc[0]
        lines.append(f"{reg} & {n} & {dr:.0%} & {o.fail_rate:.3f} & {w.fail_rate:.3f} & {o.flag_rate:.3f} & {w.flag_rate:.3f} ".replace("%", r"\%") + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_s_drift.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    m = pd.read_csv(os.path.join(RES, "model_selection_study.csv"))
    piv = m.pivot_table(index=["truth", "points", "criterion"], columns="selected", values="fraction").reset_index()
    cols = ["Langmuir", "Sips", "Toth", "Dual-Site Langmuir"]
    lines = [r"\begin{tabular}{lrlrrrr}", r"\toprule",
             r"Generating model & points & criterion & Langmuir & Sips & Toth & dual-site Langmuir \\", r"\midrule"]
    for _, r in piv.iterrows():
        lines.append(f"{r.truth} & {r.points} & {r.criterion.replace('AICc', 'AIC$_c$')} & " + " & ".join(f"{r[c]:.3f}" for c in cols) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_s_models.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    for d in (FIG, TAB):
        os.makedirs(d, exist_ok=True)
    setup()
    iso = load_isotherms()
    print(fig_pure(iso))
    fig_iast(iso)
    fig_known()
    table_iast()
    table_pure(iso)
    si_tables()
    print("figures and tables written")


if __name__ == "__main__":
    main()
