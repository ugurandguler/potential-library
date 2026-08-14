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

HERE = os.path.dirname(os.path.abspath(__file__))


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
        f = np.concatenate([L.spectrum(cry, pot, nq=n).ravel()
                            for n in meshes])
        neg = int((f < -1e-6).sum())
        frac = neg / f.size
        most = float(f.min()*L.CM1) if neg else 0.0
        v["dyn"] = {"imag_frac": round(frac, 4),
                    "most_neg_cm1": round(most, 1),
                    "stable": neg == 0, "nq": "+".join(map(str, meshes))}
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
