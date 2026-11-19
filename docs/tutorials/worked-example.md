# Worked example: RASPA3 methane in MFI

This example goes step by step through a real dataset: RASPA3 simulations of methane in silicalite-1 (MFI) at 300 K, one GCMC run per pressure from 1 kPa to 1 MPa plus a Widom insertion run for the Henry coefficient. The files are bundled with the repository in `tests/data/raspa3_mfi_ch4`.

Run it from the repository root:

```bash
python examples/raspa3_to_report.py
```

or, with the command line only:

```bash
adsorpqc assess -i tests/data/raspa3_mfi_ch4 -o raspa3_example_output
```

## 1. Reading the runs

```text
Runs read: 7 GCMC, 1 Widom
```

AdsorpQC finds every RASPA output in the folder, reads the pressure, the loading with its block error and the per-sample loading series of each GCMC run, and recognises the Widom run.

## 2. Burn-in and drift for each pressure

```text
      P (Pa)   loading (mol/kg)  burn-in    N_eff  status
        1000             0.0061        0       45  WARNING
        3000             0.0183        0       45  WARNING
       10000             0.0602        0       38  WARNING
       30000             0.1743        0       45  WARNING
      100000             0.5037        0       45  WARNING
      300000             1.1178        0       39  WARNING
     1000000             1.8959        0       45  WARNING
```

No burn-in had to be discarded, and the drift test finds every production phase stationary. For example, at 100 kPa the first and last 20 % differ by 4 % of the mean, but with \(z = 0.3\), well within the statistical noise.

Every run still gets a `WARNING`, for a different reason: the bundled files keep only about 45 samples per run, so there are fewer than 50 statistically independent samples. AdsorpQC reports this as **marginal sampling**. The averages may well be right, but the series is too short to show it. With the full production series from your own simulations this warning disappears. This is the kind of situation the checks are meant to make visible rather than hide.

## 3. Isotherm model and Henry coefficient

```text
Selected model: Toth (R^2 = 1.0000, RMSE = 0.0012)
Henry coefficient: isotherm 6.1154e-06 vs Widom 6.1084e-06 +/- 4.4e-08 mol/kg/Pa (+0.11 %)
Overall status: PASS
```

Among Langmuir, dual-site Langmuir, Sips and Toth, the Toth model has the lowest AICc. The Henry coefficient extrapolated from the isotherm agrees with the Widom insertion value within 0.11 %, inside the Widom uncertainty.

The overall status printed by the script covers the isotherm analysis. The command-line tool also includes the burn-in check of every run in its verdict, so for this dataset `adsorpqc assess` reports `WARNING` because of the marginal sampling above.

## 4. The report

The script writes the same files as the command line (see [Outputs](../outputs.md)). In `report.html`:

- **Statistical & Thermodynamic Validation** lists each check with its status and diagnostic message: the burn-in and drift results, the fitted models with their \(R^2\), RMSE and AICc, and the Henry coefficient.
- **Computational Details (Publication Ready Methods)** contains a draft methods paragraph with the decisions taken. Edit it so that it describes your study.
- **BibTeX Citation** gives the entry to cite AdsorpQC.
