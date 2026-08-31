"""
Tests for quality scoring, reporting, and CLI demo execution in AdsorpQC.
"""

import os
import tempfile
import numpy as np
import pytest
from adsorpqc.core.scoring import assess_adsorption_quality
from adsorpqc.reporters.plot_generator import generate_adsorption_figures
from adsorpqc.reporters.manuscript_prep import generate_adsorption_manuscript_assets
from adsorpqc.reporters.html_report import generate_adsorption_html_report
from adsorpqc.cli import run_demo


def test_full_adsorption_validation_pipeline():
    meta = {
        "engine": "RASPA",
        "framework": "HKUST-1",
        "adsorbate": "CH4",
        "temperature_k": 298.15
    }
    
    pressures = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    loadings = 12.0 * (0.3 * pressures) / (1.0 + 0.3 * pressures)
    
    report = assess_adsorption_quality(
        metadata=meta,
        pressure_isotherm=pressures,
        loading_isotherm=loadings,
        temperature_k=298.15
    )
    
    assert report.overall_status == "PASS"
    assert report.isotherm_fits is not None
    assert report.isotherm_fits["best_r_squared"] > 0.99
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Plot generation
        plots = generate_adsorption_figures(report, tmpdir, pressure_raw=pressures, loading_raw=loadings, formats=["png", "svg"])
        assert len(plots) > 0
        for p in plots:
            assert os.path.exists(p)
            
        # Manuscript assets
        assets = generate_adsorption_manuscript_assets(report, tmpdir)
        assert os.path.exists(assets["summary_csv"])
        assert os.path.exists(assets["summary_tex"])
        assert os.path.exists(assets["methods_text"])
        assert os.path.exists(assets["citation_bib"])
        
        # HTML report
        html_p = os.path.join(tmpdir, "report.html")
        generate_adsorption_html_report(report, html_p, methods_text="Sample methods", citation_bib="@software{}")
        assert os.path.exists(html_p)
        assert os.path.getsize(html_p) > 500


def test_cli_demo_execution():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_demo(output_dir=tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "report.html"))
        assert os.path.exists(os.path.join(tmpdir, "adsorpqc_summary_table.csv"))
        assert os.path.exists(os.path.join(tmpdir, "citation.bib"))
