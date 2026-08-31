"""
Core algorithms and scientific calculations for AdsorpQC.
"""

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
