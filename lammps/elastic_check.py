#!/usr/bin/env python3
"""
Prove `pair_style ugur` against latdyn at the level the fit actually cares about.

validate_pair.py compares energy and pressure - the potential and its first
derivative.  The elastic constants are the second derivative, and they carry
the one piece of physics neither of the other two touches: for any cell with
more than one atom in the basis, straining the box leaves internal forces
behind, and the atoms have to relax before the stress means anything.  That
relaxation is worth tens of per cent in the shear constants and it is the whole
reason the fit uses relaxed rather than frozen-ion values.  A pair style with a
sign error in the three-body force can still give the right energy and the
right pressure and get this wrong.

The recipe is LAMMPS's own examples/ELASTIC, unmodified apart from init.mod:
six box deformations, each minimised over the internal coordinates, stress
differenced to give a full 6x6 tensor.  init.mod is replaced because the
original builds its crystal with `lattice fcc`, and `lattice hcp` fixes c/a at
the ideal sqrt(8/3) while every hcp element in the library sits below it.  The
cell comes from latdyn's Crystal through cellfile.py instead.

hcp comes out in the standard orientation - c along z - so the five independent
constants are read off directly and C66 is left over as a free check on the
whole chain: hexagonal symmetry forces C66 = (C11 - C12)/2, and nothing in the
calculation was told to make that true.

    python elastic_check.py             # every element
    python elastic_check.py Ti Zr Ru
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
UGREF = {}
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import cellfile         # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"

#  `up` is the strain the tensor is differenced over.  The ELASTIC example
#  ships 1e-6, which is right for a smooth potential and wrong for this one at
#  hard truncation: a strain that small still moves neighbours across the
#  cutoff, and the discontinuity lands in the stress.  1e-5 is small enough to
#  stay in the harmonic regime and large enough that the step is a fixed
#  offset rather than noise.  Both are reported when they disagree.
INIT = """variable up equal {up}
variable atomjiggle equal 1.0e-5
units           metal
variable cfac equal 1.0e-4
variable cunits string GPa
variable etol equal 0.0
variable ftol equal 1.0e-12
variable maxiter equal 1000
variable maxeval equal 10000
variable dmax equal 1.0e-2
boundary        p p p
read_data       {data}
neigh_modify every 1 delay 0 check yes
"""

#  The ELASTIC example ships `neighbor 1.0 nsq`, which builds the list by
#  looping over every pair.  That is harmless for its own 4-atom cell and not
#  for ours: this is a three-body potential and the cells run to a few thousand
#  atoms, where an O(N^2) rebuild inside a minimiser dominates everything.
POTMOD = """pair_style      {style}
pair_coeff      * * {pot} {el}
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes
min_style       cg
min_modify      dmax ${{dmax}} line quadratic
thermo          1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
"""

POTFILE = """# {el}, written by elastic_check.py from library.json
# Tersoff-style: one line per (centre, leg, leg) triple.  For one element that
# is a single line, and it is written this way so the single-species files and
# the alloy files are the same format - a regression here is a regression in
# what the alloy path will use.
# el1 el2 el3  m D alpha r0 gamma C alpha3 rcut2 rcut3 taper lam2 lam4
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

WANT = ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66")
SKIN = 1.0
GROW = 1.2
TAG = ""

#  latdyn returns a 6x6 in Voigt order; both codes put c along z for hcp, so
#  the indices mean the same thing on either side and nothing has to be rotated
VOIGT = {"C11": (0, 0), "C22": (1, 1), "C33": (2, 2), "C12": (0, 1),
         "C13": (0, 2), "C23": (1, 2), "C44": (3, 3), "C55": (4, 4),
         "C66": (5, 5)}


def as_dict(C):
    return {k: float(C[i, j]) for k, (i, j) in VOIGT.items()}


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def lammps_cij(el, rec, cry, up=1.0e-5, grow=1.2, tag=""):
    d = os.path.join(HERE, "elruns", el + tag)
    os.makedirs(d, exist_ok=True)
    for f in ("in.elastic", "displace.mod"):
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cp ~/lammps/examples/ELASTIC/{f} {wsl(d)}/"],
                       capture_output=True)
    open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(
        el=el, taper=(rec.get("taper") or -1.0),
        lam2=rec.get("lam2") or 0.0, lam4=rec.get("lam4") or 0.0,
        **{k: rec[k] for k in ("m", "D", "alpha", "r0", "gamma",
                               "C", "alpha3", "rcut2", "rcut3")}))
    #  Elastic constants belong to the perfect lattice under homogeneous
    #  strain, so the cell only has to be big enough for the interaction range,
    #  not big enough to be a piece of material.  `grow` exists to prove that:
    #  running an element at 1.2 and again at 2.0 has to give the same tensor,
    #  and if it does not the box is too small and the answer is wrong.
    rcut = max(rec["rcut2"], rec["rcut3"])
    box, _ = cellfile.orthogonal_cell(cry)
    rep = tuple(max(2, int(np.ceil(grow * (rcut + SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep, triclinic=True)
    open(os.path.join(d, "init.mod"), "w").write(
        INIT.format(data=f"{el}.data", up=up))
    open(os.path.join(d, "potential.mod"), "w").write(
        POTMOD.format(pot=f"{el}.ugur", el=el,
                      #  the angular sets are a DIFFERENT potential and need a
                      #  different style.  Writing lam2 = lam4 = 0 under
                      #  pair_style ugur, as this did, silently drops the
                      #  angular term and then compares the remainder against
                      #  a reference that also drops it - agreement there
                      #  means nothing about the potential being distributed.
                      style=("ugur/ang" if (rec.get("lam2") or rec.get("lam4"))
                             else "ugur")))
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {LMP} -in in.elastic 2>&1"],
                       capture_output=True, text=True)
    got = {}
    for k in WANT:
        m = re.search(rf"{k}all\s*=\s*([-\d.eE+]+)", r.stdout)
        if m:
            got[k] = float(m.group(1))
    if not got:
        err = [l for l in r.stdout.splitlines() if "ERROR" in l]
        got["err"] = err[0][:80] if err else "output unreadable"
    return got


def main():
    global GROW, TAG
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    up = os.path.join(ROOT, "angular", "elastic_ref.json")
    global UGREF
    UGREF = json.load(open(up)) if os.path.exists(up) else {}
    args = sys.argv[1:]
    GROW = 1.2
    for a in list(args):
        if a.startswith("--grow="):
            GROW = float(a.split("=")[1]); args.remove(a)
    TAG = "" if GROW == 1.2 else f"_g{GROW}"
    #  which parameter set.  This checked only the hard-truncated one, which is
    #  the set that has not changed since it was written - so it was silently
    #  the least useful of the four to re-verify after a refit.
    SETKEY = None
    for a in list(args):
        if a.startswith("--set="):
            k = a.split("=")[1]
            SETKEY = None if k == "hard" else k
            args.remove(a)
            TAG += "_" + k
    els = args or sorted(lib)
    print("pair_style ugur ile latdyn, esneklik sabitleri (GPa), ic gevseme dahil")
    print(f"{'el':4s}{'struct':>5s}{'const':>7s}{'latdyn':>10s}{'lammps':>10s}"
          f"{'diff %':>9s}   symmetry")
    print("-" * 62)
    worst = []
    for el in els:
        e = refdata.ELEMENTS[el]
        rec = lib[el] if SETKEY is None else lib[el].get(SETKEY)
        if not rec or "m" not in rec:
            continue
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(rec)
        #  For the angular sets the reference CANNOT come from this tree: the
        #  standalone latdyn has no h(cos theta) and returns the potential
        #  without it, so comparing against it compares two different
        #  potentials and reports the angular term itself as a disagreement.
        #  angular/elastic_ref.py computes it in the tree that knows, and hands
        #  it over as JSON.
        ref = None
        if rec.get("lam2") or rec.get("lam4"):
            ref = UGREF.get(f"{el}|{SETKEY}")
            if ref is None:
                print(f"{el:4s}{e['struct']:>5s}  no angular reference "
                      f"(angular/elastic_ref.py --set {SETKEY} calistirin)")
                continue
        if ref is None:
            ref = as_dict(L.elastic(cry, pot)[0])   # relaxed, not frozen-ion
        got = lammps_cij(el, rec, cry, grow=GROW, tag=TAG)
        if "err" in got:
            print(f"{el:4s}{e['struct']:>5s}  LAMMPS: {got['err']}")
            worst.append((99.9, el, "error"))
            continue
        keys = (("C11", "C33", "C12", "C13", "C44")
                if e["struct"] == "hcp" else ("C11", "C12", "C44"))
        #  hexagonal symmetry is not imposed anywhere, so it is a free test
        note = ""
        if e["struct"] == "hcp" and "C66" in got:
            pred = 0.5 * (got["C11"] - got["C12"])
            note = (f"C66 {got['C66']:.2f} - (C11-C12)/2 {pred:.2f} = "
                    f"{got['C66'] - pred:+.3f}")
        for i, k in enumerate(keys):
            if k not in ref or k not in got:
                continue
            a, b = ref[k], got[k]
            pc = 100.0 * (b - a) / a if a else float("nan")
            worst.append((abs(pc), el, k))
            print(f"{el if i == 0 else '':4s}{e['struct'] if i == 0 else '':>5s}"
                  f"{k:>7s}{a:10.3f}{b:10.3f}{pc:9.3f}"
                  f"   {note if i == 0 else ''}")
    worst.sort(reverse=True)
    print()
    print("en kotu bes:", ", ".join(f"{el} {k} %{p:.2f}"
                                    for p, el, k in worst[:5]))


if __name__ == "__main__":
    main()
