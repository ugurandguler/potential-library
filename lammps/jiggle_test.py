#!/usr/bin/env python3
"""
Does the fitted lattice survive being nudged?

Every stability test in this project so far asks a question at zero amplitude.
The fit's screen diagonalises the force constants and asks whether the
reference structure is a harmonic minimum.  The phonon spectrum asks the same
thing on a mesh.  Both are answers about an infinitesimal displacement, and
three parameter sets pass them while failing something much weaker.

Niobium's switched set is the clearest.  Its lowest phonon is +23.5 cm^-1,
every mode real, no imaginary frequency anywhere on an 8^3 mesh.  Displace its
atoms by 1e-5 A - a ten-thousandth of a thermal amplitude - relax, and the
crystal does not come back.  It settles into a different configuration,
degenerate in energy to twelve figures, whose elastic constants are 161/101/41
against the 245/132/28 the fit was scored on.  At 1e-7 A it does come back, and
at 1e-3 A it reaches the same place as at 1e-5 to ten significant figures, so
the thing on the other side is a genuine basin and not numerical noise.  The
barrier is of order 1e-6 A and no harmonic analysis can see it.

That is why this test exists and why it is not redundant with the phonon
screen: the defect is ANHARMONIC.  Stable at zero amplitude, unstable at any
amplitude that matters.

The measure is the displacement that survives relaxation.  Nudge every atom by
`amp`, minimise with the cell fixed, and compare the result with where they
started.  A sound crystal returns essentially to the lattice; an unsound one
keeps a finite distortion.  Reported as the RMS displacement per atom relative
to the nudge, so the number is a ratio and comparable between elements.

Cheap on purpose - one minimisation, not the thirteen an elastic-constant run
costs - because the point is to be affordable enough to run on every record.

    python jiggle_test.py                 # tapered sets, every element
    python jiggle_test.py Nb Cr --set tap
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
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import cellfile         # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"
SKIN = 2.0
AMP = 1.0e-5
#  A lattice that returns keeps 0.2-0.3 of the nudge; one that does not keeps
#  8000-14000 of it and drops 0.2-1.3 meV/atom.  Four orders of magnitude
#  separate the two, so the threshold is not a judgement call - it was set at
#  0.20 first, which put the alkalis on the wrong side of a line drawn through
#  the middle of the "returned" cluster.  Both conditions are required: a large
#  residual displacement with no energy gain would be a flat direction, which
#  is a different (and harmless) thing.
KEEP_MAX = 10.0
DE_MAX = -5.0e-5      # eV/atom

SETS = {"hard": (None, "ugur"), "tap": ("tap", "ugur"),
        "ug": ("ug", "ugur/ang"), "tap_ug": ("tap_ug", "ugur/ang")}

POTFILE = """# {el}, written by jiggle_test.py
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

IN = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       cell.data
pair_style      {style}
pair_coeff      * * pot.ugur {el}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
min_modify      dmax 1.0e-3 line quadratic

#  compute displace/atom fixes its reference at the moment it is defined,
#  which is exactly what is wanted here and also handles periodic images -
#  an atom-style variable on xu is not allowed, and one on x would measure a
#  wrap as a displacement.
compute         disp all displace/atom

displace_atoms  all random {amp} {amp} {amp} 87287 units box
variable        e1 equal pe/atoms
run             0
print           "E_NUDGED ${{e1}}"

#  the cell is held fixed: this asks whether the STRUCTURE returns, and
#  letting the box move would confound that with an equation-of-state question
minimize        0 1e-12 100000 1000000
variable        e2 equal pe/atoms
print           "E_RELAXED ${{e2}}"

variable        d2 atom c_disp[4]*c_disp[4]
compute         sd all reduce sum v_d2
variable        rms equal sqrt(c_sd/atoms)
run             0
print           "RMS ${{rms}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def one(job):
    el, name = job
    key, style = SETS[name]
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    rec = lib[el] if key is None else lib[el].get(key)
    if not rec or "m" not in rec:
        return None
    d = os.path.join(HERE, "jigruns", f"{el}_{name}")
    os.makedirs(d, exist_ok=True)
    p = dict(rec)
    p.setdefault("lam2", 0.0)
    p.setdefault("lam4", 0.0)
    p["taper"] = rec.get("taper") or -1.0
    open(os.path.join(d, "pot.ugur"), "w").write(POTFILE.format(
        el=el, **{k: p[k] for k in ("m", "D", "alpha", "r0", "gamma", "C",
                                    "alpha3", "rcut2", "rcut3", "taper",
                                    "lam2", "lam4")}))
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    box, _ = cellfile.orthogonal_cell(cry)
    rc = max(rec["rcut2"], rec["rcut3"])
    rep = tuple(max(3, int(np.ceil(1.8 * (rc + SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, "cell.data"), cry,
                        refdata.MASSES[el], rep)
    open(os.path.join(d, "in.jig"), "w").write(IN.format(
        style=style, el=el, skin=SKIN, amp=AMP))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in in.jig > o.txt 2>&1"],
                   capture_output=True, text=True)
    lg = os.path.join(d, "log.lammps")
    if not os.path.exists(lg):
        return {"err": "kosmadi"}
    t = io.open(lg, errors="ignore").read()
    g = {k: float(m.group(1)) for k, m in
         ((k, re.search(rf"{k}\s+([-\d.eE+]+)", t))
          for k in ("E_NUDGED", "E_RELAXED", "RMS")) if m}
    if "RMS" not in g:
        err = [l for l in t.splitlines() if "ERROR" in l]
        return {"err": err[0][:70] if err else "RMS unreadable"}
    #  the nudge is uniform in [-amp, amp] per axis, so its own RMS is
    #  amp*sqrt(3/3) = amp/sqrt(3)*sqrt(3) -> amp/sqrt(3) per axis, amp overall
    keep = g["RMS"] / (AMP / np.sqrt(3) * np.sqrt(3))
    return {"rms": g["RMS"], "keep": keep,
            "dE": g["E_RELAXED"] - g["E_NUDGED"],
            "ok": bool(keep < KEEP_MAX
                       or g["E_RELAXED"] - g["E_NUDGED"] > DE_MAX)}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    keys = ["tap", "tap_ug"]
    for a in sys.argv[1:]:
        if a.startswith("--set"):
            keys = [a.split("=")[1]] if "=" in a else keys
    if "--set" in sys.argv:
        keys = [sys.argv[sys.argv.index("--set") + 1]]
        args = [a for a in args if a not in keys]
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = args or sorted(lib)
    jobs = [(el, k) for k in keys for el in els]
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"Nudge test: every atom is displaced by {AMP:g} A, then minimised")
    print(f"with the cell held fixed, to see whether it returns. ({len(jobs)} records)\n")
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(one, jobs))
    out, bad = {}, []
    print(f"{'el':4s}{'set':9s}{'kalan/sarsma':>14s}{'dE meV/atom':>13s}   verdict")
    print("-" * 52)
    for (el, k), r in zip(jobs, res):
        if r is None:
            continue
        if "err" in r:
            print(f"{el:4s}{k:9s}   {r['err']}")
            continue
        out[f"{el}|{k}"] = r
        if not r["ok"]:
            bad.append(f"{el}/{k}")
        print(f"{el:4s}{k:9s}{r['keep']:14.3f}{r['dE'] * 1000:13.4f}   "
              f"{'saglam' if r['ok'] else 'DONMUYOR'}")
    p = os.path.join(HERE, "jiggle_test.json")
    all_ = {}
    if os.path.exists(p):
        try:
            all_ = json.load(open(p))
        except Exception:
            all_ = {}
    all_.update(out)
    json.dump(all_, open(p, "w"), indent=1, sort_keys=True)
    print(f"\ndid not return {len(bad)}: {' '.join(bad) or 'none'}")
    print(f"-> jiggle_test.json ({len(all_)} records)")


if __name__ == "__main__":
    main()
