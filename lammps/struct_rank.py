#!/usr/bin/env python3
"""
Is the fitted structure this potential's ground state?

The Bain scan says that for vanadium it is not: bcc sits in a well 24 meV/atom
deep with a structure 49 meV/atom lower on the far side, and thermal motion at
109 K carries 9 meV.  That is the mechanism behind the finite-temperature
elastic collapse, and it is worth asking the question directly rather than
along one particular path - so each candidate structure is built, its lattice
constants relaxed, and the energies compared.

Every published potential tested alongside puts the right structure lowest,
by 135 to 232 meV/atom.  That is not a coincidence: a potential meant for
dynamics is checked against the competing phases as a matter of course, and
this library's objective never looked at any structure but the fitted one.

The comparison is per atom at the relaxed lattice constant of each structure,
so it is independent of the cell and of the fitted volume.  hcp relaxes c/a as
well, since holding it at the ideal value would report an hcp energy that
belongs to no minimum.

    python struct_rank.py                # every element in the library
    python struct_rank.py V Nb Fe --sets tap,tap_ug
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
os.environ.setdefault("BASEFILE",
                      os.path.join(ROOT, "standalone", "baselines.json"))

import numpy as np         # noqa: E402
import latdyn as L         # noqa: E402
import cellfile           # noqa: E402
import refdata            # noqa: E402
import elastic_T as E     # noqa: E402

if os.path.isdir(os.path.join(HERE, "potentials")):
    E.POTDIR = os.path.join(HERE, "potentials")

CANDIDATES = ("bcc", "fcc", "hcp")
SKIN = 2.0

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       cell.data
pair_style      {style}
{coeff}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
#  aniso, not iso: hcp has two independent lattice parameters and holding c/a
#  at the ideal value reports an energy that belongs to no minimum
fix             1 all box/relax aniso 0.0 vmax 0.001
#  etol must not be zero here.  With etol=0 the only stopping condition is a
#  force tolerance that a box/relax minimisation of a 432-atom cell never
#  reaches, so every run walked the full 100000 iterations and the sweep did
#  six of fifty-one in thirty-five minutes.  1e-12 eV on a per-cell energy is
#  four orders below the meV/atom differences being compared.
minimize        1e-12 1e-10 20000 200000
variable        e equal pe/atoms
variable        lx equal lx
variable        ly equal ly
variable        lz equal lz
print           "E ${{e}}"
print           "LX ${{lx}}"
print           "LY ${{ly}}"
print           "LZ ${{lz}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def potential_for(d, el, tag):
    """write the potential file into d and return (pair_style, file name)"""
    if tag in E.OURS:
        style, ext = E.OURS[tag]
        name = f"{el}.{ext}"
        p = E.PACK[el]
        open(os.path.join(d, name), "w").write(
            E.POTFILE.format(el=el, **{k: p[tag][k] for k in
                                       ("m", "D", "alpha", "r0", "gamma", "C",
                                        "alpha3", "rcut2", "rcut3", "taper",
                                        "lam2", "lam4")}))
        return style, name
    fn = tag.split("|", 1)[1]
    style = next(s for f, s in E.BASE[el] if f == fn)
    for one_f in fn.split("+"):
        open(os.path.join(d, one_f), "wb").write(
            open(os.path.join(E.POTDIR, one_f), "rb").read())
    return style, fn


def one(el, tag, struct):
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "structrank", f"{el}_{safe}_{struct}")
    os.makedirs(d, exist_ok=True)
    p = E.PACK[el]
    #  every candidate starts from a lattice constant that gives the same
    #  atomic volume as the fitted structure, so none of them begins the
    #  minimisation at an absurd density
    nper = {"bcc": 2, "fcc": 4, "hcp": 2}
    a0 = float(p["a0"])
    vat = (a0 ** 3 / nper[p["struct"]] if p["struct"] != "hcp"
           else (np.sqrt(3) / 2) * a0 ** 3 * float(p.get("c_over_a", 1.633)) / 2)
    if struct == "hcp":
        coa = 1.633
        a = (2 * vat / ((np.sqrt(3) / 2) * coa)) ** (1 / 3.0)
    else:
        coa = None
        a = (nper[struct] * vat) ** (1 / 3.0)
    try:
        cry = L.Crystal(struct, a, coa, mass=p["mass"])
        style, pot = potential_for(d, el, tag)
        box, _ = cellfile.orthogonal_cell(cry)
        rc = max(p[tag]["rcut2"], p[tag]["rcut3"]) if tag in E.OURS else 6.0
        rep = tuple(max(3, int(np.ceil(1.8 * (rc + SKIN) / b))) for b in box)
        cellfile.write_data(os.path.join(d, "cell.data"), cry, p["mass"], rep)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    open(os.path.join(d, "in.rank"), "w").write(
        IN.format(style=style, coeff=E.coeff_line(style, pot, el), skin=SKIN))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {E.LMP} -in in.rank > out.txt 2>&1"],
                   capture_output=True, text=True)
    lg = os.path.join(d, "log.lammps")
    if not os.path.exists(lg):
        return {"error": "no log"}
    t = io.open(lg, errors="ignore").read()
    m = re.search(r"^E\s+([-\d.eE+]+)", t, re.M)
    if not m:
        err = [l for l in t.splitlines() if "ERROR" in l]
        return {"error": err[0][:70] if err else "enerji cikmadi"}
    box = {}
    for k in ("LX", "LY", "LZ"):
        mm = re.search(rf"^{k}\s+([-\d.eE+]+)", t, re.M)
        if mm:
            box[k.lower()] = float(mm.group(1))
    #  box/relax aniso moves the three axes independently, so the cell is free
    #  to leave the symmetry it was built in.  If it does, the energy in the
    #  "hcp" column is not an hcp energy and the comparison says nothing.  The
    #  shape is therefore checked, not assumed: cubic cells must come back with
    #  all three axes in the ratio they started, and the orthogonal hcp cell
    #  must keep ly/lx = sqrt(3) * (ny/nx).
    nx, ny, nz = rep
    if box:
        want = {"bcc": 1.0, "fcc": 1.0, "hcp": np.sqrt(3)}[struct] * ny / nx
        got = box["ly"] / box["lx"]
        box["shape_ok"] = bool(abs(got / want - 1) < 2e-3)
        box["ly_lx"] = got
        box["ly_lx_want"] = want
        if struct != "hcp":
            zw = float(nz) / nx
            box["shape_ok"] = bool(box["shape_ok"]
                                   and abs((box["lz"] / box["lx"]) / zw - 1)
                                   < 2e-3)
    return {"E": float(m.group(1)), **box}


def main():
    #  An option's VALUE is not a positional argument.  Filtering on the "--"
    #  prefix alone leaves "tap,tap_ug" in the element list, so the sweep looks
    #  for an element by that name, finds none, runs zero jobs and reports
    #  success - and then overwrites the previous, good, struct_rank.json with
    #  an empty one.  The same mistake cost hours in elastic_T.py.
    argv = sys.argv[1:]
    sets = ["tap", "tap_ug"]
    if "--sets" in argv:
        i = argv.index("--sets")
        sets = argv[i + 1].split(",")
        del argv[i:i + 2]
    with_base = "--nobase" not in argv
    args = [a for a in argv if not a.startswith("--")]
    els = args or sorted(E.PACK)

    jobs = []
    for el in els:
        for tag in sets:
            if tag in E.PACK.get(el, {}):
                jobs += [(el, tag, s) for s in CANDIDATES]
        if with_base:
            for fn, _ in E.BASE.get(el, []):
                if all(os.path.exists(os.path.join(E.POTDIR, g))
                       for g in fn.split("+")):
                    jobs += [(el, "base|" + fn, s) for s in CANDIDATES]

    nw = max(1, (os.cpu_count() or 4) - 2)
    #  and refuse to write an empty result over a good one
    if not jobs:
        print("0 runs - element or set name not recognised; nothing was written")
        print(f"  istenen elementler: {els}")
        print(f"  istenen setler:     {sets}")
        return
    print(f"{len(jobs)} runs, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda j: one(*j), jobs))
    got = {}
    for (el, tag, s), r in zip(jobs, res):
        got.setdefault((el, tag), {})[s] = r

    print()
    print(f"{'el':4s}{'should be':>15s}{'source':>28s}"
          f"{'bcc':>9s}{'fcc':>9s}{'hcp':>9s}   en dusuk")
    print("-" * 92)
    wrong = []
    for (el, tag), d in got.items():
        want = refdata.ELEMENTS[el]["struct"]
        vals = {s: d[s].get("E") for s in CANDIDATES}
        if any(v is None for v in vals.values()):
            bad = [s for s in CANDIDATES if vals[s] is None]
            print(f"{el:4s}{want:>15s}{tag[:28]:>28s}   eksik: "
                  f"{', '.join(f'{s}={d[s].get(chr(101)+chr(114)+chr(114)+chr(111)+chr(114))}' for s in bad)[:50]}")
            continue
        low = min(vals, key=vals.get)
        rel = {s: 1000 * (vals[s] - vals[want]) for s in CANDIDATES}
        lab = {"tap": "MAU", "tap_ug": "UG"}.get(tag, tag.replace("base|", ""))
        mark = "" if low == want else f"  <-- {low.upper()} DAHA DUSUK"
        #  a cell that left its symmetry is not the structure in the heading
        bent = [s for s in CANDIDATES if d[s].get("shape_ok") is False]
        if bent:
            mark += f"   [simetri bozuldu: {','.join(bent)}]"
        print(f"{el:4s}{want:>15s}{lab[:28]:>28s}"
              f"{rel['bcc']:9.1f}{rel['fcc']:9.1f}{rel['hcp']:9.1f}"
              f"   {low}{mark}")
        if low != want and not bent:
            wrong.append((el, lab, low, rel[low]))
    print("\nsayilar meV/atom, olmasi gereken yapiya gore.  "
          "Negative = that structure is lower.")
    print(f"\nwrong ground state: {len(wrong)}/{len(got)}")
    for el, lab, low, dv in sorted(wrong, key=lambda x: x[3]):
        print(f"  {el:3s} {lab:26s} {low} {dv:+8.1f} meV/atom daha dusuk")
    #  Merge, never overwrite.  A five-element diagnostic run must not erase a
    #  ninety-five-record sweep, and in this project that has now happened to
    #  md_screen.json, compression.json, library.json and this file.  Writing
    #  through a temp file as well, so a kill mid-write leaves the old one.
    out = os.path.join(HERE, "struct_rank.json")
    old = {}
    if os.path.exists(out):
        try:
            old = json.load(open(out))
        except Exception:
            old = {}
    old.update({f"{el}|{tag}": {s: d[s] for s in CANDIDATES}
                for (el, tag), d in got.items()})
    tmp = out + ".tmp"
    json.dump(old, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, out)
    print(f"\n{len(got)} records written; file now holds {len(old)}  -> {out}")


if __name__ == "__main__":
    main()
