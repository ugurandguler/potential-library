#!/usr/bin/env python3
"""
The equation of state: how the energy behaves away from the fitted volume.

Everything the fit constrains lives at one volume.  The cohesive energy and
the pressure are imposed there, the bulk modulus is the curvature there, and
the elastic constants are curvature there.  The first derivative of the bulk
modulus with pressure, B', is not fitted anywhere and is measurable, so it is
the cheapest honest question to put to a record that was fitted at a point:
does the well have the right shape a few per cent away from its bottom?

    E(V) on a grid of isotropic scalings, Birch-Murnaghan third order fitted
    to it, and B0 and B' read off.  The fitted B0 doubles as a check on the
    fit: it must come back as the bulk modulus the record was constrained to.

Our arms and the published potentials go through the SAME door - the same
cells, the same LAMMPS build, the same static evaluation - because a B'
computed one way for ours and another way for theirs compares two methods,
not two potentials.

The cell is not relaxed: the shape is held and only the scale changes, which
is what an isotropic equation of state is.  For hcp that holds c/a at its
fitted value, so B' here is the hydrostatic response at fixed shape.

    python eos.py                      # every element, tap and tap_ug + baselines
    python eos.py Cu W --sets tap      # a subset
    python eos.py --span 0.06 --n 13   # +/- 6 per cent in 13 points
"""
import io
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
os.environ.setdefault("LMP", f"{HOME}/lammps/src/lmp_serial")
os.environ.setdefault("BASEFILE", os.path.join(ROOT, "standalone", "baselines.json"))

import numpy as np         # noqa: E402
import latdyn as L         # noqa: E402
import cellfile           # noqa: E402
import refdata            # noqa: E402
import elastic_T as E     # noqa: E402
import struct_rank as SR   # noqa: E402  (potential_for and the wsl path helper)

#  Our .ugur files are written per run from the parameter pack, but a baseline
#  is a file that has to be found.  The baselines live in the LAMMPS
#  distribution's own potentials directory, which on this machine sits inside
#  WSL; POTDIR overrides the search.  The first directory that actually holds
#  one of them wins, so a tree without the baselines still runs our arms.
def _potdir():
    unc = "\\\\wsl.localhost\\Ubuntu" + HOME.replace("/", "\\") + "\\lammps\\potentials"
    names = [f for v in E.BASE.values() for f, _ in v]
    for d in (os.environ.get("POTDIR"), os.path.join(HERE, "potentials"), unc, E.POTDIR):
        if d and os.path.isdir(d) and any(os.path.exists(os.path.join(d, n)) for n in names):
            return d
    return os.path.join(HERE, "potentials")


E.POTDIR = _potdir()

SKIN = 2.0
OUT = os.path.join(HERE, "eos.json")

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       cell.data
pair_style      {style}
{coeff}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
run             0
variable        e equal pe/atoms
variable        v equal vol/atoms
print           "E ${{e}}"
print           "V ${{v}}"
"""


def one(el, tag, scale):
    """static energy per atom and volume per atom at one isotropic scaling"""
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "eosruns", f"{el}_{safe}_{scale:.4f}")
    os.makedirs(d, exist_ok=True)
    p = E.PACK[el]
    e = refdata.ELEMENTS[el]
    try:
        cry = L.Crystal(p["struct"], scale * float(p["a0"]), e.get("c_over_a"),
                        mass=p["mass"])
        style, pot = SR.potential_for(d, el, tag)
        box, _ = cellfile.orthogonal_cell(cry)
        rc = max(p[tag]["rcut2"], p[tag]["rcut3"]) if tag in E.OURS else 8.6
        rep = tuple(max(3, int(np.ceil(1.8 * (rc + SKIN) / b))) for b in box)
        cellfile.write_data(os.path.join(d, "cell.data"), cry, p["mass"], rep)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    open(os.path.join(d, "in.eos"), "w").write(
        IN.format(style=style, coeff=E.coeff_line(style, pot, el), skin=SKIN))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {SR.wsl(d)} && {E.LMP} -in in.eos > out.txt 2>&1"],
                   capture_output=True, text=True)
    lg = os.path.join(d, "log.lammps")
    if not os.path.exists(lg):
        return {"error": "no log"}
    t = io.open(lg, errors="ignore").read()
    me = re.search(r"^E\s+([-\d.eE+]+)", t, re.M)
    mv = re.search(r"^V\s+([-\d.eE+]+)", t, re.M)
    if not (me and mv):
        err = [l for l in t.splitlines() if "ERROR" in l]
        return {"error": err[0][:70] if err else "no energy"}
    return {"E": float(me.group(1)), "V": float(mv.group(1))}


def birch_murnaghan(V, V0, E0, B0, Bp):
    """third order, energy form; B0 in eV/A^3"""
    x = (V0 / V) ** (2.0 / 3.0) - 1.0
    return E0 + 9.0 * V0 * B0 / 16.0 * (x ** 3 * Bp + x ** 2 * (6.0 - 4.0 * (V0 / V) ** (2.0 / 3.0)))


def fit_bm(V, Ea, B_guess_GPa):
    """least squares on the four parameters, seeded from the fitted values"""
    from scipy.optimize import least_squares
    V = np.asarray(V, float)
    Ea = np.asarray(Ea, float)
    i = int(np.argmin(Ea))
    p0 = [V[i], Ea[i], B_guess_GPa / 160.21766208, 4.0]
    r = least_squares(lambda p: birch_murnaghan(V, *p) - Ea, p0,
                      bounds=([0.5 * V[i], -np.inf, 1e-4, -40.0],
                              [2.0 * V[i], np.inf, 10.0, 60.0]))
    V0, E0, B0, Bp = r.x
    res = birch_murnaghan(V, V0, E0, B0, Bp) - Ea
    return dict(V0=float(V0), E0=float(E0), B0_GPa=float(B0 * 160.21766208), Bp=float(Bp),
                rms_meV=float(1000 * np.sqrt(np.mean(res ** 2))))


def main():
    argv = sys.argv[1:]
    sets = ["tap", "tap_ug"]
    span, npts = 0.05, 11
    for name, cast in (("--sets", None), ("--span", float), ("--n", int)):
        if name in argv:
            i = argv.index(name)
            v = argv[i + 1]
            del argv[i:i + 2]
            if name == "--sets":
                sets = v.split(",")
            elif name == "--span":
                span = cast(v)
            else:
                npts = cast(v)
    with_base = "--nobase" not in argv
    els = [a for a in argv if not a.startswith("--")] or sorted(E.PACK)
    SR.add_hard_arms(sets)

    #  equal steps in the linear scale, so the volume grid is symmetric in a
    scales = np.linspace(1.0 - span, 1.0 + span, npts)
    jobs = []
    for el in els:
        for tag in sets:
            if tag in E.PACK.get(el, {}):
                jobs += [(el, tag, s) for s in scales]
        if with_base:
            for fn, _ in E.BASE.get(el, []):
                if all(os.path.exists(os.path.join(E.POTDIR, g)) for g in fn.split("+")):
                    jobs += [(el, "base|" + fn, s) for s in scales]
    if not jobs:
        print("0 runs - element or set name not recognised; nothing was written")
        return
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"{len(jobs)} runs, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda j: one(*j), jobs))

    got = {}
    for (el, tag, s), r in zip(jobs, res):
        got.setdefault((el, tag), []).append(r)

    out = {}
    if os.path.exists(OUT):
        out = json.load(io.open(OUT, encoding="utf-8"))
    print(f"\n{'el':4s}{'set':>26s}{'B0 fit':>9s}{'B0 target':>11s}{'B prime':>9s}{'rms meV':>9s}")
    print("-" * 68)
    for (el, tag), rows in sorted(got.items()):
        good = [r for r in rows if "E" in r]
        if len(good) < 5:
            print(f"{el:4s}{tag[:26]:>26s}   only {len(good)} of {len(rows)} runs returned")
            continue
        B_target = float(E.PACK[el]["Cij_exp"]["B"]) if "B" in E.PACK[el].get("Cij_exp", {}) \
            else float(refdata.ELEMENTS[el]["B"])
        f = fit_bm([r["V"] for r in good], [r["E"] for r in good], B_target)
        f["B_target_GPa"] = B_target
        f["span"] = span
        f["n"] = len(good)
        out[f"{el}|{tag}"] = f
        print(f"{el:4s}{tag[:26]:>26s}{f['B0_GPa']:9.1f}{B_target:11.1f}"
              f"{f['Bp']:9.2f}{f['rms_meV']:9.3f}")
    json.dump(out, io.open(OUT, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print(f"\n{len(got)} records written; file now holds {len(out)}  -> {OUT}")


if __name__ == "__main__":
    main()
