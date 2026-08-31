"""
AdsorpQC: Automated Quality-Control, GCMC Burn-in Detection, Isotherm Fitting, and
Reproducibility Assessment for Adsorption Simulations.
"""

__version__ = "1.0.0"
__author__ = "Andre Monreal-Hernández"
__license__ = "MIT"

from adsorpqc.core.burnin import (
    detect_gcmc_burnin,
    assess_loading_drift,
    GCMCBurninResult
)
from adsorpqc.core.isotherms import (
    fit_all_isotherm_models,
    fit_langmuir,
    fit_dual_site_langmuir,
    fit_sips,
    fit_toth,
    IsothermFitResult
)
from adsorpqc.core.energetics import (
    calculate_henry_constant,
    calculate_isosteric_heat
)
from adsorpqc.core.iast import calculate_iast_selectivity
from adsorpqc.core.scoring import assess_adsorption_quality, AdsorptionValidationReport

__all__ = [
    "__version__",
    "detect_gcmc_burnin",
    "assess_loading_drift",
    "GCMCBurninResult",
    "fit_all_isotherm_models",
    "fit_langmuir",
    "fit_dual_site_langmuir",
    "fit_sips",
    "fit_toth",
    "IsothermFitResult",
    "calculate_henry_constant",
    "calculate_isosteric_heat",
    "calculate_iast_selectivity",
    "assess_adsorption_quality",
    "AdsorptionValidationReport"
]
