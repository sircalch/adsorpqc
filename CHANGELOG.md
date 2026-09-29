# Changelog

## 1.1.0 (2026-09-27)

Validated against real RASPA3 3.0.29 simulations of methane in silicalite-1 (MFI) at 300 K
(`validation/raspa3_mfi_methane.py`):
- a Widom run giving the Henry coefficient;
- GCMC at 7 pressures from 1 kPa to 1 MPa, each with 5 000 initialization and 40 000 production cycles.

The references are RASPA3's own block averages, its Widom Henry coefficient and its enthalpy of
adsorption.

### Fixed
- **RASPA3 output was not read at all.** The parser only knew RASPA2-style lines, and some of those
  (`Current loading`, `Heat of adsorption`) do not exist in either program. With a RASPA3 file it
  returned no loading and no cycle counts, a default pressure, and the framework "Unknown". It now
  reads from the RASPA3 text output: metadata, absolute and excess loadings with their block errors,
  the enthalpy of adsorption, the Widom Henry coefficient, and aligned snapshot series (loading,
  molecule count, energy) for the initialization and production phases.
- **The CLI discarded what it read.** With a RASPA file it used only a per-cycle series that was
  never found, so the loading, K_H and q_st were ignored. `assess -i` now accepts:
  - a RASPA output, with burn-in check and q_st;
  - a folder of RASPA runs, from which it builds the isotherm, fits it, checks the burn-in of every
    run and compares the isotherm K_H with the Widom K_H. The burn-in of every run enters the overall verdict.
- **Henry coefficient.** A line through the origin fitted to the four lowest-pressure points, whether
  or not they were still linear, was 4.5 % (6 standard errors) below the Widom value. K_H is now the
  zero-loading intercept of q/P against q, using only points below 15 % of the maximum loading. It
  agrees with Widom within 0.1 %.
- **Isosteric heat.**
  - The fluctuation formula needs the number of molecules N. It was given the loading in mol/kg,
    which scales the result.
  - Energies in K (RASPA writes U/k_B) were combined with RT in kJ/mol.
  - `abs()` hid sign errors, and 0 was returned when the value could not be computed.
  
  Now `energy_unit` is explicit, the molecule series is required, and NaN is returned when the value
  cannot be computed. On the RASPA3 runs, 40 snapshots give 19.62 and 20.16 kJ/mol, against RASPA3's
  19.67 and 20.18 kJ/mol. The Clausius–Clapeyron function keeps its sign.
- **Drift test.** Drift was the difference between the first and last 20 % of production relative to
  the mean, with no account of noise, so stationary low-loading runs failed. Five of the seven RASPA3
  runs gave FAIL, with "drifts" of 10–76 %. A drift is now flagged only if it is significant against
  the autocorrelation-corrected standard error of the detrended series (z > 2 / 3) and above the
  percentage thresholds. No stationary run fails now; a linear ramp still fails.
- **Model selection** uses the small-sample corrected AICc; the plain AIC favours over-parameterised
  models on the typical 5–10 isotherm points.
- **README.** It no longer claims Freundlich and BET fits, which are not performed. It states that
  IAST has not yet been validated against an independent implementation, and that the RASPA2 patterns
  have not been checked against real RASPA2 output.

- **IAST was not IAST.** `calculate_iast_selectivity` returned (q_sat,A K_A)/(q_sat,B K_B), the
  Henry-limit selectivity, and ignored pressure and composition. It now solves binary IAST with
  equal reduced spreading pressures for any fitted model (`iast_binary`). It agrees with pyIAST 1.4
  to within 3·10⁻¹² in relative terms, and exactly with extended Langmuir when both saturation
  capacities are equal. On model isotherms, the old value was up to 7 times too high at high pressure.
- The RASPA3 parser reads the loading of every component of a mixture (`components`).

### Added
- `collect_raspa_isotherm()`.
- The `energy_unit` and `gcmc_molecule_series` arguments.
- Tests on the real RASPA3 outputs (`tests/test_raspa3_real.py`, data in `tests/data/raspa3_mfi_ch4`).
