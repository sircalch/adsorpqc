# AdsorpQC

**Quality control for the post-processing of adsorption simulations and gas isotherms.**

AdsorpQC reads the output of grand canonical Monte Carlo (GCMC) simulations, such as RASPA2 or RASPA3 runs, or a simple isotherm table, and checks the steps that usually happen in private scripts before results are reported:

- **Burn-in (equilibration) cut-off**, chosen to maximise the number of statistically independent production samples.
- **Drift in the production phase**, judged against autocorrelation-corrected noise and against its size relative to the mean.
- **Isotherm model selection** among Langmuir, dual-site Langmuir, Sips and Toth by the small-sample corrected Akaike information criterion (AICc).
- **Henry coefficient** from the low-pressure limit of the isotherm and, when available, comparison with Widom insertion.
- **Isosteric heat of adsorption** from particle–energy fluctuations.
- **Binary mixture selectivity** with Ideal Adsorbed Solution Theory (IAST) and bootstrap confidence intervals.

Every decision is recorded, and the results are written as a self-contained HTML report, publication figures, LaTeX and CSV tables and a draft methods paragraph.

```bash
pip install adsorpqc
adsorpqc demo -o my_first_report/
```

Open `my_first_report/report.html` in a browser to see a full report built from bundled benchmark data.

## Where to go next

- [Installation](installation.md)
- [From RASPA output to a report](tutorials/raspa.md)
- [Isotherm table, model selection and IAST](tutorials/isotherm-csv.md)
- [How the checks work](quality-control.md): what each status means and how to change the thresholds
- [Validation](validation.md): comparisons with Widom insertion, pyIAST and binary GCMC
- [API reference](api.md)

## Scope

AdsorpQC checks the **post-processing** of adsorption data. It does not run simulations and does not replace inspection by an expert: a `PASS` means that the checks implemented here found no problem, not that the simulation set-up (force field, framework model, charges) is correct.

## Getting help

Questions and bug reports are welcome in the [issue tracker](https://github.com/sircalch/adsorpqc/issues). See [`CONTRIBUTING.md`](https://github.com/sircalch/adsorpqc/blob/master/CONTRIBUTING.md) to propose changes.
