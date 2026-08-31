"""
Vector figure generation for adsorption isotherms and GCMC trajectory analysis.
"""

from typing import List, Optional
import os
import numpy as np
import matplotlib.pyplot as plt
from adsorpqc.core.scoring import AdsorptionValidationReport

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'lines.linewidth': 2.0,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})


def generate_adsorption_figures(
    report: AdsorptionValidationReport,
    output_dir: str,
    pressure_raw: Optional[np.ndarray] = None,
    loading_raw: Optional[np.ndarray] = None,
    gcmc_trajectory: Optional[np.ndarray] = None,
    formats: List[str] = ("png", "svg", "pdf")
) -> List[str]:
    """
    Generates publication figures for adsorption isotherms and GCMC convergence.

    Parameters
    ----------
    report : AdsorptionValidationReport
        Validation report.
    output_dir : str
        Directory to save figures.
    pressure_raw : np.ndarray, optional
        Raw experimental/simulation pressures.
    loading_raw : np.ndarray, optional
        Raw loadings.
    gcmc_trajectory : np.ndarray, optional
        GCMC cycles loading series.
    formats : list of str
        Output file formats.

    Returns
    -------
    saved_paths : list of str
        List of generated file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = []

    # 1. Adsorption Isotherm Fits Figure
    if report.isotherm_fits and report.isotherm_fits.get("models") and pressure_raw is not None and loading_raw is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        
        # Plot data points
        ax.scatter(pressure_raw, loading_raw, color="#0f172a", s=45, zorder=5, label="Simulation / Exp. Data")
        
        p_dense = np.linspace(0, float(np.max(pressure_raw)) * 1.05, 300)
        models = report.isotherm_fits["models"]
        
        colors = {"Langmuir": "#2563eb", "Dual-Site Langmuir": "#16a34a", "Sips": "#d97706", "Toth": "#9333ea"}
        
        for name, fit_obj in models.items():
            col = colors.get(name, "#64748b")
            lw = 2.5 if name == report.isotherm_fits["best_model_name"] else 1.5
            ls = "-" if name == report.isotherm_fits["best_model_name"] else "--"
            
            y_pred = fit_obj.fitted_function(p_dense)
            ax.plot(p_dense, y_pred, color=col, linewidth=lw, linestyle=ls,
                    label=f"{name} ($R^2$ = {fit_obj.r_squared:.3f})")

        ax.set_xlabel("Pressure (bar)")
        ax.set_ylabel("Adsorbate Loading (mol/kg)")
        fw = report.metadata.get("framework", "Material")
        ads = report.metadata.get("adsorbate", "Adsorbate")
        ax.set_title(f"Adsorption Isotherm & Non-Linear Model Fits ({ads} in {fw})")
        ax.grid(True)
        ax.legend(frameon=True, loc="lower right")
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"adsorpqc_isotherm_fit.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    # 2. GCMC Equilibration Trajectory Figure
    if gcmc_trajectory is not None and len(gcmc_trajectory) > 10 and report.burnin_result:
        fig, ax = plt.subplots(figsize=(9, 4.5))
        
        y = np.asarray(gcmc_trajectory, dtype=float)
        n = len(y)
        cycles = np.arange(1, n + 1)
        
        t_burn = report.burnin_result.n_burnin_cycles
        
        # Plot burn-in phase
        if t_burn > 0:
            ax.plot(cycles[:t_burn], y[:t_burn], color="#ef4444", alpha=0.6, label="Burn-in Phase (Discarded)")
            ax.axvspan(1, t_burn, color="#fee2e2", alpha=0.5)
            
        # Plot production phase
        ax.plot(cycles[t_burn:], y[t_burn:], color="#10b981", alpha=0.85, label="Production Phase (Stationary)")
        
        # Production mean line
        m_prod = report.burnin_result.production_mean_loading
        ax.axhline(m_prod, color="#047857", linestyle="--", linewidth=2.0,
                   label=f"Production Mean ({m_prod:.3f} ± {report.burnin_result.production_std_loading:.3f})")
        
        ax.set_xlabel("GCMC Simulation Cycle")
        ax.set_ylabel("Adsorbate Loading (molec/uc or mol/kg)")
        ax.set_title(f"GCMC Equilibration Trajectory & Burn-in Quality ({report.burnin_result.status})")
        ax.grid(True)
        ax.legend(frameon=True, loc="lower right")
        
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"adsorpqc_gcmc_burnin_trajectory.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    return saved_files
