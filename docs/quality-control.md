# How the checks work

Each check returns a status (`PASS`, `WARNING` or `FAIL`) and a short diagnostic. The overall status is the worst of the individual statuses. A `WARNING` asks you to look at a result; a `FAIL` means it should not be reported as it is.

## Burn-in cut-off

GCMC runs start far from equilibrium, and the first cycles must be discarded. AdsorpQC tries candidate cut-offs (every `step_search` cycles, 10 by default) and, for each one, estimates the statistical inefficiency \(g\) of the remaining production series from its integrated autocorrelation time. It keeps the cut-off that **maximises the number of statistically independent samples**, \(N_{\text{eff}} = n / g\).

Production with fewer than 50 independent samples is reported as marginal sampling.

Function: [`detect_gcmc_burnin`](api.md).

## Drift in the production phase

After the cut-off, the mean loading of the first 20 % of the production phase is compared with that of the last 20 %. The difference is judged in two ways at once:

- **Statistically**, against its standard error. The segment errors use the scatter and the statistical inefficiency of the linearly detrended series, so a genuine trend is not absorbed into the noise estimate.
- **In size**, relative to the mean loading.

| Status | Condition (defaults) |
|---|---|
| `WARNING` | \(z > 2\) and drift larger than 3 % of the mean |
| `FAIL` | \(z > 3\) and drift larger than 8 % of the mean |

Requiring both conditions avoids false alarms on short or low-loading runs that are noisy but stationary. The thresholds are arguments of [`assess_loading_drift`](api.md) (`z_warn`, `z_fail`, `warn_threshold_pct`, `fail_threshold_pct`) and of `detect_gcmc_burnin`.

## Isotherm model selection

The Langmuir, dual-site Langmuir, Sips and Toth models are fitted by non-linear least squares. The model with the lowest **small-sample corrected AIC** (AICc) is selected, and \(R^2\) and RMSE are reported for each model. AICc penalises extra parameters more strongly than AIC when there are few isotherm points, which is the usual case.

Function: [`fit_all_isotherm_models`](api.md).

## Henry coefficient

\(K_H = \lim_{P \to 0} q/P\) is obtained by regressing \(q/P\) on \(q\) for the lowest-pressure points (first-order virial expansion, \(q/P = K_H - b\,q + \dots\)) and taking the intercept. This removes the bias of simply dividing the lowest loading by its pressure. When a Widom insertion run is available, the two estimates are compared.

Function: [`calculate_henry_constant`](api.md).

## Isosteric heat

From a GCMC series of the adsorbate energy \(U\) and number of molecules \(N\), the isosteric heat is computed with the particle–energy fluctuation formula (Vlugt et al., 2008):

\[
q_{\text{st}} = RT - \frac{\langle U N\rangle - \langle U\rangle\langle N\rangle}{\langle N^2\rangle - \langle N\rangle^2}
\]

The energy unit must be stated explicitly. When RASPA3 reports its own value, it is read directly.

Function: [`calculate_isosteric_heat`](api.md).

## IAST selectivity

Each single-component isotherm is fitted (by default with the best model by AICc), and binary IAST is solved with the fitted functions. The selectivity is \(S_{A/B} = (x_A/x_B)/(y_A/y_B)\). The 95 % confidence interval comes from resampling the isotherm points of each component and repeating the fit and the IAST solution.

Function: [`calculate_iast_selectivity`](api.md).
