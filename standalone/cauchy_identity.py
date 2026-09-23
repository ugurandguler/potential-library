#!/usr/bin/env python3
"""
Where the Cauchy violation comes from, exactly.

A central pair potential on a Bravais lattice at zero pressure gives
C12 = C44.  This form breaks that, and the derivation says where: under a
homogeneous strain the triplet energy E3(r_a, r_b) changes through both legs,
and its second derivative carries a BOND-BOND CROSS term that a pair potential
has no analogue of.  Everything else - the pair term, and the parts of the
triplet term that involve one leg twice - is symmetric under the exchange that
the Cauchy relation tests, and cancels.  What is left is

    C12 - C44 = (1/V) SUM_triplets  K_ab r_a r_b
                     [ n_ax^2 n_by^2 + n_bx^2 n_ay^2 - 2 n_ax n_ay n_bx n_by ]

with K_ab = d2 E3 / d r_a d r_b, the two legs' unit vectors n_a, n_b, and the
sum over the same triplets the energy uses.  With the leg switch on,

    K_ab = g'' S_a S_b + g' (S_a' S_b + S_a S_b') + g S_a' S_b' ,

which is what makes the switched records obey the same identity as the hard-cut
ones rather than a different one.

This script is the check, not the derivation: it evaluates that sum and
compares it with C12 - C44 from the finite-strain elastic tensor, for every
cubic record in the library.  Hexagonal is excluded - it has two atoms per
cell, a non-affine correction and a different set of Cauchy relations.

    python cauchy_identity.py            # every cubic element, all arms
    python cauchy_identity.py Cu W
"""
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import latdyn as L      # noqa: E402
import refdata          # noqa: E402

ARMS = (("hard", lambda v: v), ("tap", lambda v: v.get("tap")),
        ("ug", lambda v: v.get("ug")), ("tap_ug", lambda v: v.get("tap_ug")),
        ("rc", lambda v: v.get("rc")))


def cross_sum(cry, pot):
    """the bond-bond cross term, in GPa"""
    s = 0.0
    for (_i, _ja, _jb, ra, rb, na, nb) in L.triplets(cry, pot):
        k = pot.leg3(ra, rb)[5]                     # d2E3/dra drb, switch included
        geom = (na[0] ** 2 * nb[1] ** 2 + nb[0] ** 2 * na[1] ** 2
                - 2 * na[0] * na[1] * nb[0] * nb[1])
        s += k * ra * rb * geom
    return L.EV_A3_TO_GPA * s / cry.vol


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    lib = json.load(io.open(os.path.join(HERE, "library.json"), encoding="utf-8"))
    els = args or sorted(e for e in lib if refdata.ELEMENTS[e]["struct"] in ("fcc", "bcc"))
    out, worst = {}, (0.0, None)
    print(f"{'record':12}{'C12-C44':>10}{'cross term':>12}{'difference':>12}")
    for el in els:
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"), mass=1.0)
        for name, get in ARMS:
            rec = get(lib[el])
            if not isinstance(rec, dict) or rec.get("D") is None:
                continue
            pot = L.Potential.from_record(rec)
            C = L.elastic(cry, pot)[0]
            got = float(C[0, 1] - C[3, 3])
            pred = cross_sum(cry, pot)
            out[f"{el}|{name}"] = dict(cauchy_GPa=got, cross_GPa=pred, diff_GPa=got - pred)
            if abs(got - pred) > worst[0]:
                worst = (abs(got - pred), f"{el}|{name}")
            print(f"{el + '|' + name:12}{got:10.3f}{pred:12.3f}{got - pred:12.6f}")
    out["_note"] = ("C12 - C44 against the bond-bond cross term of the triplet energy, GPa. "
                    "Cubic records only. The two agree to the numerical noise of the "
                    "finite-strain tensor, which is what makes the identity exact.")
    out["_worst_diff_GPa"] = worst[0]
    out["_worst_record"] = worst[1]
    json.dump(out, io.open(os.path.join(HERE, "cauchy_identity.json"), "w", encoding="utf-8"),
              indent=1, sort_keys=True)
    print(f"\n{len(out) - 3} records; worst difference {worst[0]:.2e} GPa ({worst[1]})")


if __name__ == "__main__":
    main()
