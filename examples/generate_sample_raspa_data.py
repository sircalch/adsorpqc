"""
Generates sample RASPA output files and CSV isotherms for testing and tutorials.
"""

import os
import pandas as pd
import numpy as np


def generate_sample_raspa_file(filepath: str = "sample_raspa_output.data"):
    content = """
================================================================================
                                R A S P A - 2 . 0
================================================================================

Framework name: Mg-MOF-74
Component 0 [CO2]
Temperature: 298.15 [K]
Pressure: 100000.000000 [Pa]
Number of initialization cycles: 2000
Number of production cycles: 10000

Average loading absolute [mol/kg framework]   8.45120 +/- 0.05120
Average loading absolute [milligram/gram framework]   371.93000 +/- 2.25000
Average loading absolute [molecules/unit cell]   12.35000 +/- 0.25000

Heat of adsorption: 38.54 +/- 0.45 [kJ/mol]
"""
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Sample RASPA output written to: {os.path.abspath(filepath)}")


def generate_sample_isotherm_csv(filepath: str = "sample_co2_isotherm.csv"):
    pressures = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    # Dual-site Langmuir
    q_co2 = 6.2 * (15.0 * pressures) / (1.0 + 15.0 * pressures) + 3.5 * (0.2 * pressures) / (1.0 + 0.2 * pressures)
    df = pd.DataFrame({"Pressure_bar": pressures, "Loading_mol_kg": q_co2})
    df.to_csv(filepath, index=False)
    print(f"Sample Isotherm CSV written to: {os.path.abspath(filepath)}")


if __name__ == "__main__":
    generate_sample_raspa_file()
    generate_sample_isotherm_csv()
