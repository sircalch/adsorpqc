"""
Quickstart API tutorial for AdsorpQC.
"""

import os
from adsorpqc import assess_adsorption_quality
from adsorpqc.parsers import parse_isotherm_csv
from adsorpqc.reporters import (
    generate_adsorption_figures,
    generate_adsorption_manuscript_assets,
    generate_adsorption_html_report
)
from examples.generate_sample_raspa_data import generate_sample_isotherm_csv


def main():
    print("Running AdsorpQC Python API quickstart tutorial...")
    output_dir = "quickstart_adsorpqc_output"
    os.makedirs(output_dir, exist_ok=True)
    
    csv_file = "sample_co2_isotherm.csv"
    generate_sample_isotherm_csv(csv_file)
    
    # 1. Parse isotherm
    data = parse_isotherm_csv(csv_file)
    
    # 2. Assess calculation
    report = assess_adsorption_quality(
        metadata={"framework": "Mg-MOF-74", "adsorbate": "CO2", "temperature_k": 298.15},
        pressure_isotherm=data["pressure"],
        loading_isotherm=data["loading"],
        temperature_k=298.15
    )
    
    print(f"\nOverall Certification: {report.overall_status}")
    print(f"Validation Score: {report.validation_score}")
    if report.isotherm_fits and report.isotherm_fits.get("best_model"):
        bm = report.isotherm_fits["best_model"]
        print(f"Optimal Model: {bm.model_name} (R^2 = {bm.r_squared:.4f}, RMSE = {bm.rmse:.4f})")
    if report.henry_data:
        print(f"Henry Constant K_H: {report.henry_data['henry_constant']:.4e} mol/kg/Pa")
        
    # 3. Export all publication assets
    generate_adsorption_figures(report, output_dir, pressure_raw=data["pressure"], loading_raw=data["loading"])
    assets = generate_adsorption_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    generate_adsorption_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print(f"\nCompleted! HTML report available at: {os.path.abspath(html_p)}")


if __name__ == "__main__":
    main()
