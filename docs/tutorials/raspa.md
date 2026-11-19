# From RASPA output to a report

This tutorial runs AdsorpQC on RASPA output and explains the summary it prints.

## A single run

Point `adsorpqc assess` to one RASPA3 output file (`output/output_<T>_<P>.s0.txt`) or one RASPA2 output file (`output_*.data`):

```bash
adsorpqc assess -i output/output_300_1e+05.s0.txt -o raspa_audit/
```

For a single run, AdsorpQC reads the loading per cycle, detects the burn-in cut-off, tests the production phase for drift and, when the energy and number of molecules are reported, computes the isosteric heat.

## A folder of runs (an isotherm)

When each pressure is simulated in its own run, put the runs in one folder and pass the folder:

```bash
adsorpqc assess -i runs/ -o raspa_isotherm/
```

AdsorpQC assesses every run, builds the isotherm from their production averages, fits and ranks the isotherm models and computes the Henry coefficient. If the folder also contains a Widom insertion run, its Henry coefficient is compared with the one extrapolated from the isotherm.

The repository includes such a folder with RASPA3 results for methane in MFI (`tests/data/raspa3_mfi_ch4`, one sub-folder per pressure plus `henry` with the Widom run). From a clone of the repository:

```bash
adsorpqc assess -i tests/data/raspa3_mfi_ch4 -o raspa_isotherm/
```

Among other lines, the summary reports the selected model and the comparison with Widom insertion:

```text
 * Optimal Model     : Toth (R^2 = 1.0000)
 * Henry Const (K_H) : 6.1154e-01 mol/kg/bar (virial extrapolation)
 * Henry coefficient : isotherm 6.1154e-06 vs Widom 6.1084e-06 +/- 4.4e-08 mol/kg/Pa (+0.1%)
```

The overall status for this example is `WARNING`; open `raspa_isotherm/report.html` to see which check raised it and why. A warning asks you to look at a result, not to discard it.

## Options

| Option | Meaning |
|---|---|
| `-i`, `--input` | Isotherm table, one RASPA3/RASPA2 output file, or a folder of RASPA runs |
| `-o`, `--output` | Output directory (default `adsorpqc_output`) |
| `--framework` | Framework name for the report, e.g. `HKUST-1` |
| `--adsorbate` | Adsorbate name for the report, e.g. `CO2` |
| `--temp` | Temperature in K (default 298.15) |

See [Outputs](../outputs.md) for the files written to the output directory and [How the checks work](../quality-control.md) for the meaning of each status.
