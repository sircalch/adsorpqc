"""
Command Line Interface (CLI) for AdsorpQC.
"""

import sys
import os
import argparse
import numpy as np

from adsorpqc import __version__
from adsorpqc.parsers.raspa import parse_raspa_output
from adsorpqc.parsers.generic_isotherm import parse_isotherm_csv
from adsorpqc.core.scoring import assess_adsorption_quality
from adsorpqc.reporters.plot_generator import generate_adsorption_figures
from adsorpqc.reporters.manuscript_prep import generate_adsorption_manuscript_assets
from adsorpqc.reporters.html_report import generate_adsorption_html_report


def print_banner():
    banner = rf"""
     _       _                       ____   _____ 
    / \   __| |___  ___  _ __ _ __  / __ \ / ____|
   / _ \ / _` / __|/ _ \| '__| '_ \| |  | | |     
  / ___ \ (_| \__ \ (_) | |  | |_) | |__| | |____ 
 /_/   \_\__,_|___/\___/|_|  | .__/ \___\_\\_____| v{__version__}
                             |_|                  
 Adsorption Simulation & Isotherm Quality-Control Toolkit
 Monreal-Hernández et al., 2026
"""
    print(banner)


def run_demo(output_dir: str = "adsorpqc_demo_output"):
    """
    Simulates a benchmark GCMC simulation of CO2 adsorption in Mg-MOF-74 at 298.15 K
    with raw cycles, burn-in equilibration, non-linear isotherm fitting, K_H, q_st, and IAST selectivity.
    """
    print(f"\n[AdsorpQC] Running demonstration benchmark (CO2 in Mg-MOF-74 at 298 K)...")
    os.makedirs(output_dir, exist_ok=True)
    
    metadata = {
        "engine": "SYNTHETIC DEMO DATA (RASPA-like GCMC output; not a real simulation)",
        "framework": "Mg-MOF-74",
        "adsorbate": "CO2",
        "temperature_k": 298.15,
        "pressure_bar": 1.0,
        "init_cycles": 5000,
        "prod_cycles": 20000
    }
    
    # 1. Generate realistic GCMC equilibration trajectory (5000 cycles)
    # Reaching saturation mean = 8.45 mol/kg from initial 0
    rng = np.random.default_rng(42)
    n_cycles = 4000
    t_arr = np.arange(n_cycles)
    
    # Approach to equilibrium: loading(t) = 8.45 * (1 - exp(-t / 400)) + noise
    ideal_traj = 8.45 * (1.0 - np.exp(-t_arr / 350.0))
    noise = rng.normal(0, 0.35, size=n_cycles)
    loading_traj = np.maximum(0.0, ideal_traj + noise)
    
    # Energy trajectory: ~ -38.5 kJ/mol host-guest interaction
    energy_traj = -38.5 * loading_traj + rng.normal(0, 5.0, size=n_cycles)
    
    # 2. Generate multi-pressure isotherm data (CO2 in Mg-MOF-74)
    pressures = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])  # bar
    # Dual-site Langmuir: q1 = 6.2, K1 = 15.0, q2 = 3.5, K2 = 0.2
    q_co2 = 6.2 * (15.0 * pressures) / (1.0 + 15.0 * pressures) + 3.5 * (0.2 * pressures) / (1.0 + 0.2 * pressures)
    q_co2 += rng.normal(0, 0.08, size=len(pressures))
    
    # Secondary gas for IAST (N2: weaker adsorption, single site q_sat=3.0, K=0.08)
    q_n2 = 3.0 * (0.08 * pressures) / (1.0 + 0.08 * pressures)
    
    print("  -> Performing automated burn-in detection, production drift assessment, and isotherm fitting...")
    report = assess_adsorption_quality(
        metadata=metadata,
        gcmc_loading_series=loading_traj,
        gcmc_energy_series=energy_traj,
        pressure_isotherm=pressures,
        loading_isotherm=q_co2,
        gas_mixture={"y_a": 0.15, "y_b": 0.85},  # Flue gas: 15% CO2, 85% N2
        isotherm_b_pressure=pressures,
        isotherm_b_loading=q_n2,
        temperature_k=298.15
    )
    
    print("  -> Generating publication-ready vector figures (isotherm fits and GCMC trajectory)...")
    generate_adsorption_figures(
        report, output_dir,
        pressure_raw=pressures,
        loading_raw=q_co2,
        gcmc_trajectory=loading_traj
    )
    
    print("  -> Drafting manuscript Methods text snippet, summary LaTeX tables, and BibTeX citations...")
    assets = generate_adsorption_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing interactive report to {html_p}...")
    generate_adsorption_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print("\n" + "="*70)
    print(f" [RESULT] Overall Adsorption Simulation Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    print(f" * Framework & Gas   : {report.metadata['framework']} / {report.metadata['adsorbate']} at {report.metadata['temperature_k']} K")
    if report.burnin_result:
        br = report.burnin_result
        print(f" * GCMC Burn-in      : Cutoff at {br.n_burnin_cycles} cycles ({br.burnin_fraction*100:.1f}% discarded) | Status: {br.status}")
        print(f" * Production Loading: {br.production_mean_loading:.3f} [{br.production_ci_lower_95:.3f}, {br.production_ci_upper_95:.3f}] mol/kg (95% CI)")
        print(f" * Loading Drift     : {br.loading_drift_pct:.2f}% | N_eff = {br.n_eff:.0f} independent samples")
    if report.isotherm_fits and report.isotherm_fits.get("best_model"):
        bm = report.isotherm_fits["best_model"]
        print(f" * Optimal Model     : {bm.model_name} (R^2 = {bm.r_squared:.4f}, RMSE = {bm.rmse:.4f})")
    if report.henry_data:
        print(f" * Henry Const (K_H) : {report.henry_data['henry_constant']:.4e} mol/kg/Pa (R^2 = {report.henry_data['r_squared']:.3f})")
    if report.isosteric_heat_kj_mol:
        print(f" * Isosteric Heat qst: {report.isosteric_heat_kj_mol:.2f} kJ/mol")
    if report.iast_selectivity:
        iast = report.iast_selectivity
        print(f" * IAST Selectivity  : S = {iast['selectivity_a_b']:.2f} [{iast['ci_lower_95']:.2f}, {iast['ci_upper_95']:.2f}] (95% CI)")
    print("="*70)
    print(f"\nAll outputs successfully saved to: {os.path.abspath(output_dir)}/")
    print(f"Open {os.path.abspath(html_p)} in your browser to inspect the full report.\n")


def run_assess(args):
    """
    Evaluates user-provided isotherm CSV or RASPA calculation files.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    input_file = args.input
    if not input_file:
        print("[Error] Please specify an isotherm CSV or RASPA output file with --input.", file=sys.stderr)
        sys.exit(1)
        
    print(f"\n[AdsorpQC] Parsing input from {input_file}...")
    
    pressures = None
    loadings = None
    gcmc_loadings = None
    gcmc_energies = None
    meta = {
        "framework": args.framework or "Nanoporous Material",
        "adsorbate": args.adsorbate or "Adsorbate",
        "temperature_k": args.temp or 298.15,
        "engine": "User Input"
    }
    
    if input_file.endswith(".data") or "output_" in input_file:
        parsed_raspa = parse_raspa_output(input_file)
        meta.update(parsed_raspa["metadata"])
        if parsed_raspa["cycle_loadings"] is not None:
            gcmc_loadings = parsed_raspa["cycle_loadings"]
            gcmc_energies = parsed_raspa["cycle_energies"]
    else:
        parsed_csv = parse_isotherm_csv(input_file)
        pressures = parsed_csv["pressure"]
        loadings = parsed_csv["loading"]
        
    print("  -> Performing quality certification, burn-in validation, and model fitting...")
    report = assess_adsorption_quality(
        metadata=meta,
        gcmc_loading_series=gcmc_loadings,
        gcmc_energy_series=gcmc_energies,
        pressure_isotherm=pressures,
        loading_isotherm=loadings,
        temperature_k=meta.get("temperature_k", 298.15)
    )
    
    print("  -> Generating publication figures...")
    generate_adsorption_figures(
        report, output_dir,
        pressure_raw=pressures,
        loading_raw=loadings,
        gcmc_trajectory=gcmc_loadings
    )
    
    print("  -> Generating manuscript text, LaTeX summary table, and BibTeX citations...")
    assets = generate_adsorption_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    print(f"  -> Writing HTML quality report to {html_p}...")
    generate_adsorption_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print("\n" + "="*70)
    print(f" [RESULT] Overall Adsorption Certification: {report.overall_status}")
    print(f" [SCORE]  {report.validation_score}")
    print("="*70)
    if report.burnin_result:
        print(f" * Burn-in Status    : {report.burnin_result.status} (Drift = {report.burnin_result.loading_drift_pct:.2f}%)")
    if report.isotherm_fits and report.isotherm_fits.get("best_model"):
        bm = report.isotherm_fits["best_model"]
        print(f" * Optimal Model     : {bm.model_name} (R^2 = {bm.r_squared:.4f})")
    print("="*70)
    print(f"\nReport ready at: {os.path.abspath(html_p)}\n")


def print_citation():
    bib = """@software{monreal2026adsorpqc,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{AdsorpQC: An Open-Source Toolkit for Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and Reproducibility Assessment of Adsorption Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/adsorpqc}
}"""
    print("\nIf you use AdsorpQC in your publications, please cite:\n")
    print("APA Style:")
    print("Monreal-Hernández, A. (2026). AdsorpQC: An Open-Source Toolkit for Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and Reproducibility Assessment of Adsorption Simulations (v1.0.0). Zenodo. https://github.com/sircalch/adsorpqc\n")
    print("BibTeX:")
    print(bib)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="adsorpqc",
        description="AdsorpQC: Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and Reproducibility Toolkit."
    )
    parser.add_argument("-v", "--version", action="version", version=f"adsorpqc {__version__}")
    
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")
    
    # Assess command
    assess_parser = subparsers.add_parser("assess", help="Assess isotherm CSV or RASPA simulation file")
    assess_parser.add_argument("-i", "--input", required=True, help="Path to isotherm CSV or RASPA output_*.data")
    assess_parser.add_argument("-o", "--output", default="adsorpqc_output", help="Directory for output report and assets (default: adsorpqc_output)")
    assess_parser.add_argument("--framework", default=None, help="Framework name (e.g. MOF-5, HKUST-1)")
    assess_parser.add_argument("--adsorbate", default=None, help="Adsorbate gas name (e.g. CO2, CH4, N2)")
    assess_parser.add_argument("--temp", type=float, default=298.15, help="Simulation / experimental temperature in K")
    
    # Demo command
    demo_parser = subparsers.add_parser("demo", help="Run AdsorpQC on a benchmark GCMC / isotherm dataset")
    demo_parser.add_argument("-o", "--output", default="adsorpqc_demo_output", help="Output directory (default: adsorpqc_demo_output)")
    
    # Cite command
    subparsers.add_parser("cite", help="Display BibTeX and APA citation details")
    
    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)
        
    args = parser.parse_args()
    
    if args.command == "assess":
        print_banner()
        run_assess(args)
    elif args.command == "demo":
        print_banner()
        run_demo(args.output)
    elif args.command == "cite":
        print_banner()
        print_citation()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

