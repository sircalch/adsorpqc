"""
RASPA3 3.0.29 campaign for the AdsorpQC IAST validation: CO2, and CO2/CH4 mixtures, in silicalite-1
(MFI, 2x2x2 unit cells) at 300 K. Force field: the RASPA3 examples 8 (CH4) and 10 (CO2, TraPPE-like,
framework charges O -1.025, Si +2.05, Ewald), which share the framework parameters.

  co2_henry/      Widom insertion of CO2 (Henry coefficient)
  co2_p_<Pa>/     pure-CO2 GCMC isotherm points
  mix_y<yCO2>_p_<Pa>/  CO2/CH4 mixture GCMC (gas-phase mole fraction yCO2)
The pure-CH4 isotherm is raspa3_mfi_methane.py (same framework and CH4 parameters).

Usage: python raspa3_co2_ch4_mfi.py OUT_DIR [init prod] [--workers 8]
Requires raspa3 (conda-forge); RASPA3_PREFIX points to the environment (default %LOCALAPPDATA%\\r3).
"""
import argparse
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
PREFIX = os.environ.get("RASPA3_PREFIX", os.path.join(os.environ.get("LOCALAPPDATA", ""), "r3"))
EXAMPLES = os.path.join(PREFIX, "share", "raspa3", "examples", "basic")
RASPA = os.path.join(PREFIX, "bin", "raspa3.exe") if os.name == "nt" else "raspa3"
PRESSURES = [1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]
MIXTURES = [(0.5, 1e4), (0.5, 1e5), (0.5, 1e6), (0.1, 1e5)]
T = 300.0


def keep_awake():
    if os.name == "nt":
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)


def force_field():
    co2 = json.load(open(os.path.join(EXAMPLES, "10_mc_adorption_co2_in_mfi", "force_field.json")))
    ch4 = json.load(open(os.path.join(EXAMPLES, "8_mc_adsorption_of_methane_in_mfi", "force_field.json")))
    names = {p["name"] for p in co2["PseudoAtoms"]}
    co2["PseudoAtoms"] += [p for p in ch4["PseudoAtoms"] if p["name"] not in names]
    names = {s["name"] for s in co2["SelfInteractions"]}
    co2["SelfInteractions"] += [s for s in ch4["SelfInteractions"] if s["name"] not in names]
    return co2


def component(name, mol_fraction=None, widom=False):
    c = {"Name": name, "FugacityCoefficient": 1.0, "CreateNumberOfMolecules": 0}
    if widom:
        c.update({"WidomProbability": 1.0, "IdealGasRosenbluthWeight": 1.0})
        return c
    c.update({"TranslationProbability": 0.5, "ReinsertionProbability": 0.5, "SwapProbability": 1.0})
    if name == "CO2":
        c["RotationProbability"] = 0.5
    if mol_fraction is not None:            # as in the RASPA3 example non_basic/1 (CO2/CH4 in IRMOF-1)
        c["MolFraction"] = mol_fraction
    return c


def simulation(components, pressure, init, prod):
    system = {"Type": "Framework", "Name": "MFI_SI", "NumberOfUnitCells": [2, 2, 2],
              "ExternalTemperature": T, "ChargeMethod": "Ewald", "HeliumVoidFraction": 0.3}
    if pressure is not None:
        system["ExternalPressure"] = pressure
    return {"SimulationType": "MonteCarlo", "NumberOfCycles": prod, "NumberOfInitializationCycles": init,
            "PrintEvery": 1000, "RandomSeed": 20260929, "Systems": [system], "Components": components}


def setup(out, name, sim, ff):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    src = os.path.join(EXAMPLES, "10_mc_adorption_co2_in_mfi")
    for f in ("MFI_SI.cif", "CO2.json"):
        shutil.copy(os.path.join(src, f), d)
    shutil.copy(os.path.join(EXAMPLES, "8_mc_adsorption_of_methane_in_mfi", "methane.json"), os.path.join(d, "CH4.json"))
    json.dump(ff, open(os.path.join(d, "force_field.json"), "w"), indent=2)
    json.dump(sim, open(os.path.join(d, "simulation.json"), "w"), indent=2)
    return d


def run(d):
    if any(f.endswith(".txt") and "Final state" in open(os.path.join(d, "output", f), errors="ignore").read()
           for f in (os.listdir(os.path.join(d, "output")) if os.path.isdir(os.path.join(d, "output")) else [])):
        return os.path.basename(d), "cached"
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([PREFIX, os.path.join(PREFIX, "bin"), os.path.join(PREFIX, "Library", "bin"), env["PATH"]])
    env["RASPA_DIR"] = os.path.join(PREFIX, "share", "raspa3")
    with open(os.path.join(d, "stdout.log"), "w") as fh:
        r = subprocess.run([RASPA], cwd=d, stdout=fh, stderr=subprocess.STDOUT, env=env)
    return os.path.basename(d), r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("init", nargs="?", type=int, default=5000)
    ap.add_argument("prod", nargs="?", type=int, default=30000)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    keep_awake()
    ff = force_field()
    dirs = [setup(a.out, "co2_henry", simulation([component("CO2", widom=True)], None, 0, a.prod), ff)]
    for p in PRESSURES:
        dirs.append(setup(a.out, f"co2_p_{int(p)}", simulation([component("CO2")], p, a.init, a.prod), ff))
    for y, p in MIXTURES:
        comps = [component("CO2", y), component("CH4", round(1 - y, 6))]
        dirs.append(setup(a.out, f"mix_y{y:g}_p_{int(p)}", simulation(comps, p, a.init, a.prod), ff))
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for name, rc in ex.map(run, dirs):
            print(name, "exit", rc, flush=True)


if __name__ == "__main__":
    main()
