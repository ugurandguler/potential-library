#!/usr/bin/env python3
"""
Add dynamical-stability data to library.json.

Born stability and dynamical stability are different criteria and, in a crystal,
largely independent of each other: Born tests the elastic matrix, which is the
q -> 0 limit only, while dynamical stability requires every mode at every q to
be real.  The fit only ever sees elastic constants, so nothing in it prevents a
potential that is elastically fine and dynamically unstable - which is exactly
what Li turned out to be: 24.7 % of modes imaginary, down to -120 cm^-1, over
most of the Brillouin zone.

Stores, per element, under "dyn":
    imag_frac   fraction of modes with omega^2 < 0 on the mesh
    most_neg    the most negative frequency, cm^-1 (0 if stable)
    stable      imag_frac == 0

    python add_dynstab.py [nq ...]
"""
import json, os, sys
import numpy as np
import latdyn as L
import refdata
#  the SAME path build_library.py draws the dispersion panel on, so the badge
#  and the picture cannot disagree
from build_library import sc_segments, NQ_PATH

HERE = os.path.dirname(os.path.abspath(__file__))

#  a mode counts as imaginary below this, in cm^-1
IMAG_TOL_CM1 = -1.0


def main(meshes=(8, 9)):
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print(f"{'el':4s}{'imag modes':>12s}{'most neg':>11s}  status")
    print("-" * 44)
    bad = []
    for el in sorted(lib):
        v = lib[el]
        cry = L.Crystal(v["struct"], v["a0"], v.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(v)
        #  One mesh is not enough.  L.spectrum uses a half-shifted
        #  Monkhorst-Pack grid, so 4^3, 8^3 and 16^3 are nearly nested and step
        #  over the same places - an 8^3 mesh alone calls Ta and V stable when
        #  9^3 finds modes at -37 and -67 cm^-1, near the bcc N point where the
        #  known anomalies of Nb, V, W and the alkalis live.  Take the union of
        #  an even and an odd mesh.
        #  A Monkhorst-Pack mesh cannot see the symmetry lines.  L.mesh is
        #  half-shifted - (arange(n)+0.5)/n - 0.5 - so it never lands on Gamma,
        #  on a zone-boundary point, or on the line between two of them, and
        #  refining it does not help because a shifted grid steps over the same
        #  places however fine it gets.  An instability confined to a line is
        #  invisible to it.  Iridium is the case that showed this: -8.1 cm^-1
        #  half way along K-Gamma, called stable here for as long as only the
        #  mesh was sampled - while the dispersion panel on the same record
        #  drew the negative branch for anyone who looked.  Sample the path as
        #  well, and the same one the panel uses.
        qs = []
        for (_, ka, __, kb) in sc_segments(refdata.ELEMENTS[el]["struct"]):
            qs.extend(ka + np.linspace(0, 1, NQ_PATH)[:, None] * (kb - ka))
        f = np.concatenate(
            [L.spectrum(cry, pot, nq=n).ravel() for n in meshes]
            + [L.frequencies_many(cry, pot, np.array(qs)).ravel()]) * L.CM1
        #  The tolerance is PHYSICAL, not a floating-point epsilon.  The three
        #  acoustic branches are exactly zero at Gamma and come back from the
        #  diagonalisation a shade under it, so a cut at -1e-6 THz flags modes
        #  at -1e-4 cm^-1 - four orders below anything real, on the 0.08 % of
        #  modes that ARE the Gamma point.  refit/dynscreen.py hit exactly that
        #  and called beryllium and chromium unstable.  Real instabilities in
        #  this library run -6 to -63 cm^-1, so the range between is empty.
        neg = int((f < IMAG_TOL_CM1).sum())
        frac = neg / f.size
        most = float(f.min()) if neg else 0.0
        v["dyn"] = {"imag_frac": round(frac, 4),
                    "most_neg_cm1": round(most, 1),
                    "stable": neg == 0,
                    "nq": "+".join(map(str, meshes)),
                    #  the viewer renders nq as "8^3 and 9^3", so the path
                    #  travels as its own field rather than inside that string
                    "path": "Setyawan-Curtarolo"}
        tag = "stable" if neg == 0 else "DYNAMICALLY UNSTABLE"
        if neg:
            bad.append((frac, el, most))
        print(f"{el:4s}{100*frac:11.1f}%{most:11.1f}  {tag}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\n{len(lib)-len(bad)}/{len(lib)} dynamically stable")
    for frac, el, most in sorted(bad, reverse=True):
        print(f"   {el}: {100*frac:.1f}% imaginary, down to {most:.0f} cm-1")
    print(f"\nmerged into {path}")


if __name__ == "__main__":
    main(tuple(int(a) for a in sys.argv[1:]) or (8, 9))
