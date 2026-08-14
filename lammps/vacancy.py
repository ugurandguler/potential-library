#!/usr/bin/env python3
"""
Vacancy formation energy: the first test of anything this was not fitted to.

Everything the library has been checked against so far is either a fitted
target - the elastic constants, the lattice constant, the cohesive energy - or
a vibrational property that follows from the same second derivatives.  The
phonons are genuinely out of sample and they agree to 6-32 per cent, but they
are still the curvature of the same energy surface at the same geometry.

A vacancy is not.  Removing an atom asks what the potential does when the
coordination changes, which is the question every defect, surface and
diffusion property reduces to, and it is the question a pair potential plus an
angle-blind three-body term has the least reason to get right.  It is also the
cheapest such test there is.

    E_vac = E(N-1 atoms, one site empty, relaxed) - (N-1)/N * E(N atoms)

The relaxation matters and is easy to get wrong in both directions.  Positions
must relax - an unrelaxed vacancy overestimates by a few tenths of an eV - and
the cell volume must NOT, because a single vacancy in a real crystal sits in an
infinite matrix that does not contract around it.  `fix box/relax` here would
measure something else.

Compared against 420 published-potential values from JARVIS-FF, the same
database used for the elastic comparison, for the same reason: a number with no
distribution around it cannot be judged.

    python vacancy.py               # every element, tapered set
    python vacancy.py Cu Al --set hard
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

SETS = {"hard": (None, "ugur"), "tap": ("tap", "ugur"),
        "ug": ("ug", "ugur/ang"), "tap_ug": ("tap_ug", "ugur/ang")}

POTFILE = """# {el}, written by vacancy.py
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

IN = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       {data}
pair_style      {style}
pair_coeff      * * {pot} {el}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
{remove}
#  positions relax, the cell does not: one vacancy sits in an infinite matrix
minimize        0 1e-10 5000 50000
variable        etot equal pe
variable        nat  equal count(all)
print           "OUT ${{etot}} ${{nat}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def one(el, rec, style, tag, remove):
    d = os.path.join(HERE, "vacruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    p = dict(rec)
    p.setdefault("lam2", 0.0)
    p.setdefault("lam4", 0.0)
    p["taper"] = rec.get("taper") or -1.0
    open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(
        el=el, **{k: p[k] for k in ("m", "D", "alpha", "r0", "gamma", "C",
                                    "alpha3", "rcut2", "rcut3", "taper",
                                    "lam2", "lam4")}))
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    box, _ = cellfile.orthogonal_cell(cry)
    rc = max(rec["rcut2"], rec["rcut3"])
    #  the vacancy must not see its own periodic image; 2.2x the cutoff on
    #  every axis is what makes that true rather than nearly true
    rep = tuple(max(3, int(np.ceil(2.2 * (rc + SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep)
    open(os.path.join(d, "in.vac"), "w").write(IN.format(
        data=f"{el}.data", style=style, pot=f"{el}.ugur", el=el, skin=SKIN,
        remove=remove))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in in.vac 2>&1"],
                   capture_output=True, text=True)
    lg = os.path.join(d, "log.lammps")
    if not os.path.exists(lg):
        return None
    m = re.search(r"OUT\s+([-\d.eE+]+)\s+([-\d.eE+]+)",
                  io.open(lg, errors="ignore").read())
    if not m:
        return None
    return float(m.group(1)), int(float(m.group(2)))


def evac(el, rec, style):
    perfect = one(el, rec, style, "perfect", "")
    if perfect is None:
        return None
    withvac = one(el, rec, style, "vac",
                  "group v id 1\ndelete_atoms group v compress no")
    if withvac is None:
        return None
    E0, N = perfect
    E1, M = withvac
    if M != N - 1:
        return None
    return E1 - (M / N) * E0


def main():
    args = sys.argv[1:]
    which = "tap"
    if "--set" in args:
        i = args.index("--set")
        which = args[i + 1]
        del args[i:i + 2]
    key, style = SETS[which]
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    ref = {}
    rp = os.path.join(ROOT, "standalone", "jarvis_vac.json")
    if os.path.exists(rp):
        ref = json.load(open(rp))
    els = args or sorted(lib)

    jobs = []
    for el in els:
        rec = lib[el] if key is None else lib[el].get(key)
        if rec and "m" in rec:
            jobs.append((el, dict(rec)))
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"Bosluk olusum enerjisi, {which} seti  ({len(jobs)} element, "
          f"{nw} in parallel)")
    print("Positions relax, cell volume fixed.\n")
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda j: evac(j[0], j[1], style), jobs))

    print(f"{'el':4s}{'ours':>9s}{'pub. lowest':>15s}{'median':>10s}"
          f"{'en yuksek':>11s}{'n':>4s}   konum")
    print("-" * 60)
    out = {}
    for (el, rec), v in zip(jobs, res):
        if v is None:
            print(f"{el:4s}   hesaplanamadi")
            continue
        out[el] = v
        pub = sorted(x["Ev"] for x in ref.get(el, []))
        if not pub:
            print(f"{el:4s}{v:9.2f}{'':>15s}{'':>10s}{'':>11s}{0:4d}"
                  f"   none published")
            continue
        import statistics as st
        below = sum(1 for x in pub if x < v)
        where = ("dagilimin ICINDE" if pub[0] <= v <= pub[-1]
                 else "ALTINDA" if v < pub[0] else "USTUNDE")
        print(f"{el:4s}{v:9.2f}{pub[0]:15.2f}{st.median(pub):10.2f}"
              f"{pub[-1]:11.2f}{len(pub):4d}   {where}")
    json.dump(out, open(os.path.join(HERE, f"vacancy_{which}.json"), "w"),
              indent=1, sort_keys=True)
    print(f"\n-> vacancy_{which}.json")


if __name__ == "__main__":
    main()
