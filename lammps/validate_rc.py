#!/usr/bin/env python3
"""
Do LAMMPS and latdyn agree about the RE-CUT potentials?

validate_pair.py asks this of the published arm - it reads lib[el] and
nothing else.  The candidates live in lib[el]["rc"] and lib[el]["rc_ug"] and
have never been put through it, which matters now because the two sides
disagree about what those potentials do.

The disagreement: npt_expansion, which runs in LAMMPS, returns a thermal
expansion of about -37e-6/K for the re-cut alkalis.  latdyn, on the same
records, says every one of them is dynamically stable, has a Grueneisen
parameter near +1.5 with not one negative mode out of 213, and a single
energy minimum.  Positive Grueneisen means positive expansion.  Both cannot
be describing the same potential.

Thermostat damping was the first suspect and has been ruled out: rerunning
caesium with the damping matched to its own phonon period, 1.0 ps instead of
0.1 ps, gives -36.8 against -36.7.

What is left is the possibility that the two implementations of the potential
differ, and the re-cut fits are where that would first show.  They pushed the
pair exponent m from about 1.2 to between 5 and 14 - caesium's is 13.68 -
which is a region the published parameters never visit and the pair style has
therefore never been tested in.

    python validate_rc.py [El ...]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))

import validate_pair as V          # noqa: E402
import latdyn as L                 # noqa: E402
import refdata                     # noqa: E402


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or ["Cs", "Rb", "K", "Na", "Li", "Cu", "Mo"]
    print(f"{'el':4s}{'arm':7s}{'m':>8s}{'E latdyn':>11s}{'E lammps':>11s}"
          f"{'dE':>11s}{'P latdyn':>10s}{'P lammps':>10s}   durum")
    print("-" * 79)
    for el in els:
        e = refdata.ELEMENTS[el]
        for arm in ("tap", "rc"):
            rec = lib[el].get(arm)
            if not rec or "m" not in rec:
                continue
            pot = L.Potential.from_record(rec)
            cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                            mass=refdata.MASSES[el])
            E_ref = L.energy(cry, pot)
            g = V.lammps(el, rec, e, cry)
            if "E" not in g:
                print(f"{el:4s}{arm:7s}{rec['m']:8.2f}   LAMMPS hata: "
                      f"{g.get('err','?')[:40]}")
                continue
            h = 1e-5
            nat = len(cry.frac)
            ep = L.energy(cry.strained(np.eye(3) * h), pot)
            em = L.energy(cry.strained(np.eye(3) * (-h)), pot)
            P_ref = ((ep - em) / (2 * h)
                     / (3 * cry.vol / nat)) * L.EV_A3_TO_GPA
            P_lmp = -g["P"] / V.BAR_PER_GPA
            dE = g["E"] - E_ref
            ok = abs(dE) < 1e-5 * max(abs(E_ref), 1.0) and abs(P_lmp - P_ref) < 0.05
            print(f"{el:4s}{arm:7s}{rec['m']:8.2f}{E_ref:11.5f}{g['E']:11.5f}"
                  f"{dE:11.2e}{P_ref:10.4f}{P_lmp:10.4f}   "
                  + ("tamam" if ok else "AYRISIYOR"))


if __name__ == "__main__":
    main()
