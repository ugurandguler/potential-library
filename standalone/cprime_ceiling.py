#!/usr/bin/env python3
"""
Which part of the (C', C44) plane can this functional form actually reach?

Every bad fit in the library fails the same way: C44 lands on experiment and
C' = (C11 - C12)/2 comes out far too small.  The first version of this probe
maximised C' on its own and found it can be driven to 550 GPa - but C44 went to
1245 GPa with it.  So the form has no trouble making C' large; what it cannot do
is make C' large *relative to* C44.

The bulk modulus is a hard constraint, so C11 + 2 C12 is pinned and the two
remaining degrees of freedom in a cubic tensor are exactly C' and C44.  The
question is therefore the reachable region in that plane, and the scale-free
statement of it is the range of

    R = C44 / C'                    (half the Zener anisotropy is C44/C')

Elements fit well when their experimental R sits inside the reachable range and
badly when it does not.  Normalising by B makes the clouds comparable between
elements.

    python cprime_ceiling.py             # the interesting cubic elements
    python cprime_ceiling.py Nb V Cu

Writes cprime_region.json (the point clouds) for plotting.
"""
import json
import os
import sys

import numpy as np

import fit as F
import refdata

#  The question is what the functional form can reach, so the stability screen
#  is switched off - an unstable point still proves reachability, and leaving it
#  on would confuse "the form cannot do it" with "the form can, but only
#  unstably".  Both are reported separately below.
F.REQUIRE_DYNAMICAL_STABILITY = False

HERE = os.path.dirname(os.path.abspath(__file__))

#  coarse but wide: the point is the boundary of the region, not a converged
#  optimum anywhere inside it
M_GRID = [1.2, 1.6, 2.2, 3.0, 4.0, 5.5, 7.5, 10.0, 13.0, 16.0, 20.0]
G_GRID = [0.0, 0.25, 0.5, 1.0, 1.6, 2.4, 3.5, 5.0]
S_GRID = [0.4, 0.7, 1.0, 1.4, 2.0, 3.0]
C_GRID = [-60, -40, -25, -15, -8, -4, -1.5, 0.0, 1.5, 4, 8, 15, 25, 40, 60]


def scan(el):
    """every valid (C', C44) the grid can reach, in units of B"""
    e = refdata.ELEMENTS[el]
    B = e["B"]
    pts = []
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for m in M_GRID:
            for g in G_GRID:
                for s in S_GRID:
                    for C in C_GRID:
                        try:
                            r = F.evaluate(el, e, m, g, C, s)
                        except (OverflowError, ValueError,
                                FloatingPointError, np.linalg.LinAlgError):
                            continue
                        if not r:
                            continue
                        c = r["Cij"]
                        cp = 0.5 * (c["C11"] - c["C12"])
                        if cp <= 1e-6 or c["C44"] <= 1e-6:
                            continue
                        pts.append((cp / B, c["C44"] / B))
    return np.array(pts), B


def main():
    els = sys.argv[1:] or ["Nb", "Cr", "Mo", "V", "W", "Ta", "Al",
                           "Cu", "Ni", "Ag", "Pd", "Pt", "Fe", "K"]
    print(f"{'el':4s}{'R_exp':>8s}{'R_min':>8s}{'R_max':>8s}"
          f"{'inside?':>9s}{'C/B exp':>9s}{'C/B max':>9s}{'pts':>7s}   RMS")
    print("-" * 78)

    out, verdicts = {}, []
    for el in els:
        e = refdata.ELEMENTS[el]
        if e["struct"] == "hcp":
            print(f"{el:4s}  (hcp, skipped - four independent constants)")
            continue
        c = e["Cij"]
        cp_exp = 0.5 * (c["C11"] - c["C12"])
        R_exp = c["C44"] / cp_exp
        pts, B = scan(el)
        if not len(pts):
            print(f"{el:4s}  no valid solution anywhere on the grid")
            continue
        R = pts[:, 1] / pts[:, 0]
        inside = R.min() <= R_exp <= R.max()
        verdicts.append((el, R_exp, inside))
        out[el] = dict(B=B, R_exp=R_exp, R_min=float(R.min()),
                       R_max=float(R.max()), pts=pts.tolist())
        print(f"{el:4s}{R_exp:8.2f}{R.min():8.2f}{R.max():8.1f}"
              f"{('yes' if inside else 'NO'):>9s}"
              f"{cp_exp / B:9.3f}{pts[:, 0].max():9.3f}{len(pts):7d}")

    #  merge, never replace: this is normally run on a few elements at a time
    p = os.path.join(HERE, "cprime_region.json")
    old = json.load(open(p)) if os.path.exists(p) else {}
    old.update(out)
    tmp = p + ".tmp"
    json.dump(old, open(tmp, "w"))
    os.replace(tmp, p)

    if verdicts:
        bad = [el for el, _, ins in verdicts if not ins]
        print(f"\nexperimental C44/C' unreachable for: {bad or 'none'}")
        print("reachable for: "
              f"{[el for el, _, ins in verdicts if ins]}")
        print("\nB is a hard constraint, so C11 + 2 C12 is pinned and C' and C44"
              "\nare the only two free cubic constants - this is the whole "
              "reachable set.")


if __name__ == "__main__":
    main()
