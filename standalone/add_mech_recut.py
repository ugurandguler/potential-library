#!/usr/bin/env python3
"""
Voigt-Reuss-Hill averages for the re-cut candidate arms.

The mechanical table compares the top-level MAU record against `ug`, because
those are the two arms `mech` was ever computed for.  The candidates carry
their full elastic tensor `C` and nothing derived from it, so the table has
been showing two columns where it could show four.

This runs the SAME mechanics.analyse that build_library runs, on the same
tensor, so the numbers are comparable rather than merely similar.  Nothing is
refitted and `C` is not touched.

Note what this is NOT: the re-cut arms do not belong in the "MAU against UG"
table, which exists to isolate the ANGULAR term and already refuses to draw
UG when the two fits used different three-body cutoffs.  A re-cut arm differs
by cutoff by construction, so putting it there would put back exactly the
confound that table was built to exclude.  A Voigt-Reuss-Hill average carries
no such confound: it is what one tensor implies, whatever cutoff produced it.

    python add_mech_recut.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import latdyn as L                       # noqa: E402
import mechanics                         # noqa: E402
import refdata                           # noqa: E402


def tensor(struct, c):
    """the 6x6 from the independent constants of that symmetry

    The candidate records keep `Cij` - three constants for a cubic cell, five
    for hcp - and not the `Cfull` build_library stores, so the tensor is
    rebuilt here.  Nothing is approximated: those are exactly the independent
    components, and for hcp C66 = (C11 - C12) / 2 is a symmetry identity, not
    a fit.  `C` in these records is the three-body coefficient of the
    potential and NOT an elastic constant; reading it as one is how the first
    version of this file made every tensor singular.
    """
    C = np.zeros((6, 6))
    if struct == "hcp":
        c11, c12, c13, c33, c44 = (c["C11"], c["C12"], c["C13"],
                                   c["C33"], c["C44"])
        C[0, 0] = C[1, 1] = c11
        C[2, 2] = c33
        C[0, 1] = C[1, 0] = c12
        C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = c13
        C[3, 3] = C[4, 4] = c44
        C[5, 5] = 0.5 * (c11 - c12)
    else:
        c11, c12, c44 = c["C11"], c["C12"], c["C44"]
        for i in range(3):
            C[i, i] = c11
            for j in range(3):
                if i != j:
                    C[i, j] = c12
        C[3, 3] = C[4, 4] = C[5, 5] = c44
    return C


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    n, bad = 0, []
    print(f"{'el':4s}{'arm':8s}{'B_H':>8s}{'G_H':>8s}{'E_H':>8s}{'nu_H':>8s}"
          f"{'A_U':>8s}")
    print("-" * 48)
    for el in sorted(lib):
        v = lib[el]
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        nat = len(cry.frac)
        for arm in ("rc", "rc_ug"):
            rec = v.get(arm)
            if not isinstance(rec, dict) or not isinstance(rec.get("Cij"), dict):
                continue
            need = (("C11", "C12", "C13", "C33", "C44") if e["struct"] == "hcp"
                    else ("C11", "C12", "C44"))
            if any(k not in rec["Cij"] for k in need):
                continue
            C = tensor(e["struct"], rec["Cij"])
            try:
                m = mechanics.analyse(C, mass_amu=refdata.MASSES[el] * nat,
                                      volume_A3=cry.vol, natoms=nat)
                pl = {p: mechanics.plane_curves(C, p)
                      for p in ("xy", "xz", "yz")}
            except np.linalg.LinAlgError:
                bad.append(f"{el}/{arm}")     # singular C, as build_library has it
                rec["mech"] = None
                continue
            rec["mech"] = m
            rec["mech_planes"] = pl
            n += 1
            print(f"{el:4s}{arm:8s}{m['B_H']:8.1f}{m['G_H']:8.1f}"
                  f"{m['E_H']:8.1f}{m['nu_H']:8.3f}{m['A_U']:8.3f}")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print()
    print(f"{n} candidate arms now carry a mechanical analysis"
          + (f"; singular tensor for {' '.join(bad)}" if bad else ""))


if __name__ == "__main__":
    main()
