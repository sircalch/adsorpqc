# Isotherm table, model selection and IAST

AdsorpQC also works with an adsorption isotherm stored as a table, whether it comes from simulations or from experiments.

## The input table

A CSV file with one pressure column and one loading column is enough. The pressure column is recognised by its name (for example `Pressure_bar`, `pressure`, `P (kPa)`); otherwise the first two numeric columns are used. Tab-separated files (`.tsv`, `.tab`) and `.txt` files separated by tabs, semicolons, commas or spaces are also read.

```text
Pressure_bar,Loading_mol_kg
0.01,0.8156816801180248
0.05,2.691796322489392
0.1,3.7886274509803926
```

The repository includes this example as `sample_co2_isotherm.csv`.

## Command line

```bash
adsorpqc assess -i sample_co2_isotherm.csv --framework "HKUST-1" --adsorbate CO2 --temp 298.15 -o csv_report/
```

The summary reports the selected model and the Henry coefficient, which is given in the loading and pressure units of the input:

```text
 [RESULT] Overall Adsorption Certification: PASS
 * Optimal Model     : Dual-Site Langmuir (R^2 = 1.0000)
 * Henry Const (K_H) : 9.3248e+01 loading/pressure units of the input (virial extrapolation)
```

## Python API

The same analysis is available from Python. This example fits and ranks the isotherm models and builds the full report object:

```python
from adsorpqc import assess_adsorption_quality, fit_all_isotherm_models
from adsorpqc.parsers import parse_isotherm_csv

data = parse_isotherm_csv("sample_co2_isotherm.csv")

fits = fit_all_isotherm_models(data["pressure"], data["loading"])
print(fits["best_model_name"])

report = assess_adsorption_quality(
    metadata={"framework": "HKUST-1", "adsorbate": "CO2"},
    pressure_isotherm=data["pressure"],
    loading_isotherm=data["loading"],
    temperature_k=298.15,
)
print(report.overall_status)
```

## Binary selectivity with IAST

Given the single-component isotherms of two gases measured at the same temperature, `calculate_iast_selectivity` fits each isotherm, solves IAST at the requested composition and total pressure, and estimates a 95 % confidence interval by bootstrap:

```python
from adsorpqc import calculate_iast_selectivity

result = calculate_iast_selectivity(
    isotherm_a_p=p_co2, isotherm_a_q=q_co2,
    isotherm_b_p=p_ch4, isotherm_b_q=q_ch4,
    gas_mole_fraction_a=0.5, gas_mole_fraction_b=0.5,
    total_pressure=1.0,      # same pressure units as the isotherms
)
print(result["selectivity_a_b"], result["ci_lower_95"], result["ci_upper_95"])
```

Both isotherms must use the same pressure and loading units. See the [API reference](../api.md) for all arguments and the returned values, and [Validation](../validation.md) for how the IAST implementation compares with pyIAST and binary GCMC.
