#!/usr/bin/env python3
"""
The self-interstitial formation energy: the quantity beside the vacancy.

    E_I = E(N+1) - (N+1)/N * E(N)

in the same cell, with the same relaxation, as vacancy.py - positions free,
box fixed.  The vacancy is the defect this form is known to get wrong; the
interstitial is the other end of the same question, an atom crammed in rather
than taken out, and it has never been computed for the shipped library.

WHICH SPLIT.  Which dumbbell is stable depends on the potential: bcc metals
mostly prefer <111> and chromium takes <110>, fcc metals mostly <100>.
Choosing one would answer the question in advance, so all three are tried from
the usual starting guesses and the lowest is kept, with its orientation.
Hexagonal is not attempted: its interstitial sites are a longer list (basal
and non-basal split, octahedral, tetrahedral, crowdion) and getting them wrong
silently is easier than getting them right.

The comparison is not to experiment.  Formation energies of self-interstitials
are not measured directly; the reference here is density-functional theory,
which puts them at several eV for the transition metals, and the useful test
at this level is the SIGN and the order of magnitude.

    python interstitial.py                 # every cubic element, switched arm
    python interstitial.py Nb W --set tap
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
LMP = os.environ.get("LMP", f"{HOME}/lammps/src/lmp_serial")
SKIN = 2.0
SETS = {"hard": (None, "ugur"), "tap": ("tap", "ugur"),
        "ug": ("ug", "ugur/ang"), "tap_ug": ("tap_ug", "ugur/ang")}
SITES = {"<111>": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
         "<110>": np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
         "<100>": np.array([1.0, 0.0, 0.0])}

POTFILE = """# {el}, written by interstitial.py
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
#  positions relax, the cell does not: one defect in an infinite matrix
minimize        0 1e-10 5000 50000
variable        etot equal pe
variable        nat  equal count(all)
print           "OUT ${{etot}} ${{nat}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def run_cell(el, rec, style, tag, extra_atom=None):
    """one static relaxation; extra_atom is (moved_first_atom, added_atom)"""
    d = os.path.join(HERE, "siaruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    q = dict(rec)
    q.setdefault("taper", -1.0)
    q.setdefault("lam2", 0.0)
    q.setdefault("lam4", 0.0)
    open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(el=el, **q))
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"), mass=refdata.MASSES[el])
    box, _ = cellfile.orthogonal_cell(cry)
    rc = max(rec["rcut2"], rec["rcut3"])
    rep = tuple(max(3, int(np.ceil(2.2 * (rc + SKIN) / b))) for b in box)
    path = os.path.join(d, f"{el}.data")
    cellfile.write_data(path, cry, refdata.MASSES[el], rep)
    if extra_atom is not None:
        add_atom(path, *extra_atom)
    open(os.path.join(d, "in.sia"), "w").write(IN.format(
        data=f"{el}.data", style=style, pot=f"{el}.ugur", el=el, skin=SKIN))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in in.sia > out.txt 2>&1"],
                   capture_output=True, text=True)
    lg = os.path.join(d, "log.lammps")
    if not os.path.exists(lg):
        return None
    m = re.search(r"OUT\s+([-\d.eE+]+)\s+([-\d.eE+]+)", io.open(lg, errors="ignore").read())
    return (float(m.group(1)), int(float(m.group(2)))) if m else None


def add_atom(path, shift, offset):
    """rewrite a LAMMPS data file: move atom 1 by -shift and append one at +shift

    The file is edited rather than rebuilt so that the perfect cell and the
    defective one differ by exactly one atom and nothing else."""
    lines = io.open(path, encoding="utf-8", errors="ignore").read().splitlines()
    ia = next(i for i, l in enumerate(lines) if l.strip().startswith("Atoms"))
    head = lines[:ia]
    body = [l for l in lines[ia + 1:] if l.strip()]
    first = body[0].split()
    x, y, z = (float(first[2]), float(first[3]), float(first[4]))
    body[0] = f"{first[0]} {first[1]} {x - shift[0]:.10f} {y - shift[1]:.10f} {z - shift[2]:.10f}"
    n = len(body)
    body.append(f"{n + 1} {first[1]} {x + offset[0]:.10f} {y + offset[1]:.10f} "
                f"{z + offset[2]:.10f}")
    head = [re.sub(r"^\s*\d+\s+atoms", f"{n + 1} atoms", l) for l in head]
    io.open(path, "w", encoding="utf-8", newline="\n").write(
        "\n".join(head) + "\nAtoms # atomic\n\n" + "\n".join(body) + "\n")


def one(el, lib, which):
    key, style = SETS[which]
    rec = lib[el] if key is None else lib[el].get(key)
    if not isinstance(rec, dict) or rec.get("D") is None:
        return None
    if refdata.ELEMENTS[el]["struct"] == "hcp":
        return {"skipped": "hexagonal"}
    perfect = run_cell(el, rec, style, "perfect")
    if perfect is None:
        return {"error": "perfect cell did not run"}
    E0, N = perfect
    a0 = refdata.ELEMENTS[el]["a0"]
    out = {}
    for tag, u in SITES.items():
        off = 0.30 * a0 * u
        r = run_cell(el, rec, style, tag.strip("<>"), extra_atom=(off, off))
        if r is None:
            continue
        E1, M = r
        if M != N + 1:
            continue
        out[tag] = E1 - (M / N) * E0
    if not out:
        return {"error": "no dumbbell ran"}
    best = min(out, key=out.get)
    return dict(E_I=out[best], site=best, per_site=out, natoms=N, set=which)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    which = "tap"
    if "--set" in sys.argv:
        which = sys.argv[sys.argv.index("--set") + 1]
    lib = json.load(io.open(os.path.join(ROOT, "standalone", "library.json"), encoding="utf-8"))
    els = args or sorted(e for e in lib if refdata.ELEMENTS[e]["struct"] != "hcp")
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"{len(els)} elements, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda e: one(e, lib, which), els))
    out, path = {}, os.path.join(HERE, f"interstitial_{which}.json")
    if os.path.exists(path):
        out = json.load(io.open(path, encoding="utf-8"))
    print(f"\n{'el':4}{'E_I (eV)':>10}{'site':>8}  per site")
    for el, r in zip(els, res):
        if not r or "E_I" not in r:
            print(f"{el:4}  {r.get('skipped') or r.get('error') if r else 'no record'}")
            continue
        out[el] = r
        print(f"{el:4}{r['E_I']:10.3f}{r['site']:>8}  "
              + ", ".join(f"{k} {v:.3f}" for k, v in sorted(r["per_site"].items())))
    out["_note"] = ("E_I = E(N+1) - (N+1)/N E(N), positions relaxed at a fixed cell, the lowest "
                    "of the <111>, <110> and <100> dumbbells. Hexagonal elements are not "
                    "attempted. Same cell rule as vacancy.py: at least 2.2 interaction ranges "
                    "on every axis.")
    json.dump(out, io.open(path, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print(f"\n{len([k for k in out if not k.startswith('_')])} records -> {path}")


if __name__ == "__main__":
    main()
