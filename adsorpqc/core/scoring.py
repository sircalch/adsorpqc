"""
Adsorption Simulation and Isotherm Quality Scoring Engine.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
import numpy as np

from adsorpqc.core.burnin import GCMCBurninResult, detect_gcmc_burnin
from adsorpqc.core.isotherms import fit_all_isotherm_models, IsothermFitResult
from adsorpqc.core.energetics import calculate_henry_constant, calculate_isosteric_heat
from adsorpqc.core.iast import calculate_iast_selectivity


@dataclass
class AdsorptionValidationReport:
    overall_status: str  # 'PASS', 'WARNING', 'FAIL'
    validation_score: str
    metadata: Dict[str, Any]
    burnin_result: Optional[GCMCBurninResult]
    isotherm_fits: Optional[Dict[str, Any]]
    henry_data: Optional[Dict[str, float]]
    isosteric_heat_kj_mol: Optional[float]
    iast_selectivity: Optional[Dict[str, Any]]
    recommendations: List[str]
    provenance: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_adsorption_quality(
    metadata: Dict[str, Any],
    gcmc_loading_series: Optional[np.ndarray] = None,
    gcmc_energy_series: Optional[np.ndarray] = None,
    pressure_isotherm: Optional[np.ndarray] = None,
    loading_isotherm: Optional[np.ndarray] = None,
    gas_mixture: Optional[Dict[str, float]] = None,
    isotherm_b_pressure: Optional[np.ndarray] = None,
    isotherm_b_loading: Optional[np.ndarray] = None,
    temperature_k: float = 298.15
) -> AdsorptionValidationReport:
    """
    Assesses the quality, statistical stationarity, and thermodynamic consistency of adsorption simulations and isotherms.

    Parameters
    ----------
    metadata : dict
        Material name, framework, adsorbate, simulation engine.
    gcmc_loading_series : np.ndarray, optional
        Raw timeseries of loading per cycle.
    gcmc_energy_series : np.ndarray, optional
        Raw timeseries of energy per cycle.
    pressure_isotherm : np.ndarray, optional
        Pressure grid.
    loading_isotherm : np.ndarray, optional
        Equilibrium loading values across pressures.
    temperature_k : float, default 298.15
        Temperature in Kelvin.

    Returns
    -------
    report : AdsorptionValidationReport
        Complete certification report.
    """
    statuses = []
    recommendations = []
    
    # 1. GCMC Burn-in and Stationarity Check
    burn_res = None
    q_st = None
    if gcmc_loading_series is not None and len(gcmc_loading_series) > 0:
        burn_res = detect_gcmc_burnin(gcmc_loading_series)
        statuses.append(burn_res.status)
        if burn_res.status != "PASS":
            recommendations.append(burn_res.diagnostic_message)
            
        # Energetics (q_st)
        if gcmc_energy_series is not None:
            # Use production slice for q_st
            prod_u = gcmc_energy_series[burn_res.n_burnin_cycles:]
            prod_n = gcmc_loading_series[burn_res.n_burnin_cycles:]
            q_st = calculate_isosteric_heat(prod_u, prod_n, temperature_k=temperature_k)

    # 2. Isotherm Fitting & Henry's Law
    iso_res = None
    henry_res = None
    if pressure_isotherm is not None and loading_isotherm is not None and len(pressure_isotherm) >= 3:
        iso_res = fit_all_isotherm_models(pressure_isotherm, loading_isotherm)
        henry_res = calculate_henry_constant(pressure_isotherm, loading_isotherm)
        
        if iso_res["best_r_squared"] < 0.95:
            statuses.append("WARNING")
            recommendations.append(f"Low isotherm fitting quality (Best model {iso_res['best_model_name']} R^2 = {iso_res['best_r_squared']:.3f} < 0.95). Check for stepped adsorption / gate opening.")
        else:
            statuses.append("PASS")

    # 3. IAST Selectivity Check
    iast_res = None
    if (pressure_isotherm is not None and loading_isotherm is not None and
        isotherm_b_pressure is not None and isotherm_b_loading is not None):
        ya = gas_mixture.get("y_a", 0.5) if gas_mixture else 0.5
        yb = gas_mixture.get("y_b", 0.5) if gas_mixture else 0.5
        iast_res = calculate_iast_selectivity(
            pressure_isotherm, loading_isotherm,
            isotherm_b_pressure, isotherm_b_loading,
            gas_mole_fraction_a=ya, gas_mole_fraction_b=yb
        )

    # Overall scoring
    if "FAIL" in statuses:
        overall_status = "FAIL"
        validation_score = "ADSORPTION SIMULATION QUALITY = REJECTED (NON-EQUILIBRIUM / SEVERE DRIFT)"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        validation_score = "ADSORPTION SIMULATION QUALITY = ACCEPTABLE WITH WARNINGS"
    else:
        overall_status = "PASS"
        validation_score = "ADSORPTION SIMULATION QUALITY = FULLY CERTIFIED"

    return AdsorptionValidationReport(
        overall_status=overall_status,
        validation_score=validation_score,
        metadata=metadata,
        burnin_result=burn_res,
        isotherm_fits=iso_res,
        henry_data=henry_res,
        isosteric_heat_kj_mol=q_st,
        iast_selectivity=iast_res,
        recommendations=recommendations,
        provenance={
            "tool": "AdsorpQC",
            "version": "1.0.0",
            "citation": "Monreal-Hernández, A. (2026). AdsorpQC: An Open-Source Toolkit for Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and Reproducibility Assessment of Adsorption Simulations."
        }
    )
