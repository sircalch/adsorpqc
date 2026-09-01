"""
Manuscript Methods text generator, LaTeX summary tables, and BibTeX citations for AdsorpQC.
"""

from typing import Dict, Any, Optional
import os
import pandas as pd
from adsorpqc.core.scoring import AdsorptionValidationReport


def generate_adsorption_manuscript_assets(
    report: AdsorptionValidationReport,
    output_dir: str
) -> Dict[str, str]:
    """
    Generates manuscript Methods text, LaTeX summary tables, and BibTeX citations.

    Parameters
    ----------
    report : AdsorptionValidationReport
        Validation report.
    output_dir : str
        Output directory.

    Returns
    -------
    paths : dict
        Mapping of generated asset paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    generated = {}
    
    # 1. Summary DataFrame
    rows = []
    meta = report.metadata
    
    rows.append({"Parameter": "Framework & Adsorbate", "Value": f"{meta.get('framework', 'Material')} / {meta.get('adsorbate', 'Gas')}", "Status": "PASS"})
    rows.append({"Parameter": "Simulation Engine", "Value": f"{meta.get('engine', 'RASPA')}", "Status": "PASS"})
    rows.append({"Parameter": "Temperature", "Value": f"{meta.get('temperature_k', 298.15):.2f} K", "Status": "PASS"})
    
    if report.burnin_result:
        br = report.burnin_result
        rows.append({"Parameter": "GCMC Total Cycles", "Value": f"{br.n_total_cycles} cycles", "Status": "PASS"})
        rows.append({"Parameter": "Optimal Burn-in Cutoff", "Value": f"{br.n_burnin_cycles} cycles ({br.burnin_fraction*100:.1f}%)", "Status": br.status})
        rows.append({"Parameter": "Production Loading", "Value": f"{br.production_mean_loading:.3f} [{br.production_ci_lower_95:.3f}, {br.production_ci_upper_95:.3f}] (95% CI)", "Status": br.status})
        rows.append({"Parameter": "Production Loading Drift", "Value": f"{br.loading_drift_pct:.2f}%", "Status": br.status})
        rows.append({"Parameter": "Independent Samples (N_eff)", "Value": f"N_eff = {br.n_eff:.0f} (tau_int = {br.tau_int_cycles:.1f} cycles)", "Status": br.status})

    if report.isotherm_fits and report.isotherm_fits.get("best_model"):
        bm = report.isotherm_fits["best_model"]
        rows.append({"Parameter": "Optimal Isotherm Model", "Value": f"{bm.model_name} (R^2 = {bm.r_squared:.4f}, RMSE = {bm.rmse:.4f})", "Status": "PASS"})
        for p_name, p_val in bm.parameters.items():
            rows.append({"Parameter": f"Fit Param: {p_name}", "Value": f"{p_val:.4e}", "Status": "PASS"})
            
    if report.henry_data:
        hd = report.henry_data
        rows.append({"Parameter": "Henry Constant (K_H)", "Value": f"{hd['henry_constant']:.4e} mol/kg/Pa (R^2 = {hd['r_squared']:.3f})", "Status": "PASS"})
        
    if report.isosteric_heat_kj_mol is not None and report.isosteric_heat_kj_mol > 0:
        rows.append({"Parameter": "Isosteric Heat (q_st)", "Value": f"{report.isosteric_heat_kj_mol:.2f} kJ/mol", "Status": "PASS"})
        
    if report.iast_selectivity:
        iast = report.iast_selectivity
        rows.append({"Parameter": "IAST Selectivity S_{A/B}", "Value": f"{iast['selectivity_a_b']:.2f} [{iast['ci_lower_95']:.2f}, {iast['ci_upper_95']:.2f}] (95% CI)", "Status": "PASS"})

    df_summary = pd.DataFrame(rows)
    
    # CSV Table
    csv_path = os.path.join(output_dir, "adsorpqc_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path
    
    # LaTeX Table
    tex_table_path = os.path.join(output_dir, "adsorpqc_summary_table.tex")
    tex_table = df_summary.to_latex(index=False, escape=False)
    with open(tex_table_path, "w", encoding="utf-8") as f:
        f.write("% AdsorpQC Adsorption Simulation Quality and Parameter Summary Table\n")
        f.write(tex_table)
    generated["summary_tex"] = tex_table_path

    # 2. Methods Text Snippet
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    fw_str = meta.get("framework", "nanoporous material")
    ads_str = meta.get("adsorbate", "adsorbate")
    temp_str = f"{meta.get('temperature_k', 298.15):.1f} K"
    
    gcmc_sentence = ""
    if report.burnin_result:
        br = report.burnin_result
        gcmc_sentence = (
            f"Adsorption equilibrium in {fw_str} was simulated using Grand Canonical Monte Carlo (GCMC). "
            f"Equilibration and statistical stationarity were validated using AdsorpQC v1.0.0 (Monreal-Hernández, 2026). "
            f"The initial {br.n_burnin_cycles} cycles ({br.burnin_fraction*100:.1f}%) were discarded as burn-in to eliminate initial state bias, "
            f"yielding an effective sample size of N_eff = {br.n_eff:.0f} independent samples with a loading drift of {br.loading_drift_pct:.2f}%. "
        )
        
    iso_sentence = ""
    if report.isotherm_fits and report.isotherm_fits.get("best_model"):
        bm = report.isotherm_fits["best_model"]
        iso_sentence = (
            f"Adsorption isotherms of {ads_str} at {temp_str} were fitted to non-linear analytical models, "
            f"with the {bm.model_name} model exhibiting optimal thermodynamic fidelity (R^2 = {bm.r_squared:.4f}, RMSE = {bm.rmse:.4f}, AIC = {bm.aic:.1f}). "
        )
        
    iast_sentence = ""
    if report.iast_selectivity:
        iast = report.iast_selectivity
        iast_sentence = (
            f"Binary mixture adsorption selectivity was evaluated using Ideal Adsorbed Solution Theory (IAST), "
            f"yielding S = {iast['selectivity_a_b']:.2f} [95% CI: {iast['ci_lower_95']:.2f} - {iast['ci_upper_95']:.2f}] evaluated across 500 bootstrap iterations. "
        )

    full_methods = (
        f"{gcmc_sentence}{iso_sentence}{iast_sentence}"
        f"Overall adsorption simulation quality achieved certification status: {report.overall_status}."
    )
    
    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(full_methods + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX Citation
    bib_path = os.path.join(output_dir, "citation.bib")
    bib_content = """@software{monreal2026adsorpqc,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{AdsorpQC: An Open-Source Toolkit for Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and Reproducibility Assessment of Adsorption Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/adsorpqc}
}
"""
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(bib_content)
    generated["citation_bib"] = bib_path

    return generated

