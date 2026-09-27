"""RASPA3 3.0.29 validation campaign for AdsorpQC: methane in MFI (TraPPE / RASPA example force field), 300 K.

  henry/          Widom insertion -> Henry coefficient
  p_<Pa>/         GCMC isotherm points
Usage: python raspa3_mfi_methane.py [init_cycles prod_cycles]   (defaults 5000 20000; the tests used 5000 40000)
Requires raspa3 3.0.29 (conda-forge); set RASPA below to its executable. Outputs are copied to
tests/data/raspa3_mfi_ch4/.
"""
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = HERE   # MFI_SI.cif, force_field.json and methane.json from the RASPA3 example 8_mc_adsorption_of_methane_in_mfi
RASPA = r"C:\Users\Andre\AppData\Local\r3\Library\bin\raspa3.exe"
if not os.path.exists(RASPA):
    RASPA = r"C:\Users\Andre\AppData\Local\r3\bin\raspa3.exe"
PRESSURES = [1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]
INIT, PROD = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (5000, 20000)


def setup(name, sim):
    d = os.path.join(HERE, "runs", name)
    os.makedirs(d, exist_ok=True)
    for f in ("MFI_SI.cif", "force_field.json", "methane.json"):
        shutil.copy(os.path.join(SRC, f), d)
    json.dump(sim, open(os.path.join(d, "simulation.json"), "w"), indent=2)
    return d


def gcmc(p):
    return {
        "SimulationType": "MonteCarlo", "NumberOfCycles": PROD, "NumberOfInitializationCycles": INIT,
        "PrintEvery": 1000, "RandomSeed": 12345,
        "Systems": [{"Type": "Framework", "Name": "MFI_SI", "NumberOfUnitCells": [2, 2, 2],
                     "ExternalTemperature": 300.0, "ExternalPressure": p, "ChargeMethod": "None"}],
        "Components": [{"Name": "methane", "FugacityCoefficient": 1.0, "IdealGasRosenbluthWeight": 1.0,
                        "TranslationProbability": 0.5, "ReinsertionProbability": 0.5, "SwapProbability": 1.0,
                        "CreateNumberOfMolecules": 0}],
    }


def widom():
    return {
        "SimulationType": "MonteCarlo", "NumberOfCycles": PROD, "NumberOfInitializationCycles": 0,
        "PrintEvery": 1000, "RandomSeed": 12345,
        "Systems": [{"Type": "Framework", "Name": "MFI_SI", "NumberOfUnitCells": [2, 2, 2],
                     "ExternalTemperature": 300.0, "ChargeMethod": "None"}],
        "Components": [{"Name": "methane", "IdealGasRosenbluthWeight": 1.0, "WidomProbability": 1.0,
                        "CreateNumberOfMolecules": 0}],
    }


ENV = dict(os.environ)
_R3 = r"C:\Users\Andre\AppData\Local\r3"
ENV["PATH"] = os.pathsep.join([_R3, os.path.join(_R3, "bin"), os.path.join(_R3, "Library", "bin"), ENV["PATH"]])
ENV["RASPA_DIR"] = os.path.join(_R3, "share", "raspa3")


def run(d):
    with open(os.path.join(d, "stdout.log"), "w") as fh:
        r = subprocess.run([RASPA], cwd=d, stdout=fh, stderr=subprocess.STDOUT, env=ENV)
    return os.path.basename(d), r.returncode


if __name__ == "__main__":
    dirs = [setup("henry", widom())] + [setup(f"p_{int(p)}", gcmc(p)) for p in PRESSURES]
    with ThreadPoolExecutor(max_workers=len(dirs)) as ex:
        for name, rc in ex.map(run, dirs):
            print(name, "exit", rc, flush=True)
