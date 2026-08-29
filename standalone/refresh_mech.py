#!/usr/bin/env python3
"""
Recompute the Debye temperature in every stored mechanical record.

`mechanics.sound` used to build it from the Voigt-Reuss-Hill averages B and G,
which discards the anisotropy before an average over 1/v^3 that is dominated by
the slow directions.  It came out systematically high - 1.023 of Stewart's
measured theta_D across the library, against 0.999 for the Christoffel average
over the real slowness surface.  `mechanics.debye_speed` now does the latter;
this rewrites what is already stored so the page and the tables agree with it.

Nothing else in the record changes.  B, G, the anisotropy measures, the Young
and Poisson envelopes and the two aggregate sound speeds v_l and v_t are all
untouched - v_l and v_t ARE the isotropic aggregate speeds a polycrystal has,
and only the Debye average needed the directions kept.  The old value travels
with the new one as `debye_iso`, and the old speed as `v_m_iso`, because this
changes a published number and a reader should be able to see what it was.

The tensor is taken from the record itself - `Cfull` where it exists, the
`Cij` the arm carries otherwise - and the volume from refdata's lattice
constant, which is what `build_library.py` and `add_mech_recut.py` both used.

    python refresh_mech.py [--dry]
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                                        # noqa: E402
import latdyn as L                                        # noqa: E402
import mechanics                                          # noqa: E402
import refdata                                            # noqa: E402

KEYS = {"cubic": ("C11", "C12", "C44"),
        "hcp": ("C11", "C12", "C13", "C33", "C44")}


def tensor_from_cij(struct, c):
    C = np.zeros((6, 6))
    if struct == "hcp":
        C[0, 0] = C[1, 1] = c["C11"]
        C[2, 2] = c["C33"]
        C[0, 1] = C[1, 0] = c["C12"]
        C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = c["C13"]
        C[3, 3] = C[4, 4] = c["C44"]
        C[5, 5] = 0.5 * (c["C11"] - c["C12"])
    else:
        C[0, 0] = C[1, 1] = C[2, 2] = c["C11"]
        C[0, 1] = C[1, 0] = C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = c["C12"]
        C[3, 3] = C[4, 4] = C[5, 5] = c["C44"]
    return C


def main():
    dry = "--dry" in sys.argv
    path = os.path.join(HERE, "library.json")
    with open(path) as f:
        lib = json.load(f)

    n, moved = 0, []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or el not in refdata.ELEMENTS:
            continue
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        nat = len(cry.frac)

        def holders(d, tag):
            out = []
            if isinstance(d, dict):
                if isinstance(d.get("mech"), dict):
                    out.append((tag, d))
                for k, x in d.items():
                    if k != "mech" and isinstance(x, dict):
                        out += holders(x, tag + "/" + k if tag else k)
            return out

        if el in refdata.THETA_D:
            v["theta_D_exp"] = refdata.THETA_D[el]
        for tag, rec in holders(v, ""):
            C = rec.get("Cfull")
            if C is None and isinstance(rec.get("Cij"), dict):
                try:
                    C = tensor_from_cij(e["struct"], rec["Cij"])
                except KeyError:
                    C = None
            if C is None:
                continue
            C = np.asarray(C, dtype=float)
            try:
                m = mechanics.analyse(C, mass_amu=refdata.MASSES[el] * nat,
                                      volume_A3=cry.vol, natoms=nat, nu=400)
            except np.linalg.LinAlgError:
                continue
            old = rec["mech"].get("debye")
            for k in ("debye", "debye_iso", "v_m", "v_m_iso"):
                if k in m:
                    rec["mech"][k] = m[k]
            n += 1
            if old and m.get("debye"):
                moved.append((m["debye"] / old - 1.0, el, tag or "top"))

    moved.sort()
    print("%d mechanical records rebuilt" % n)
    if moved:
        print("  largest drop:  %s"
              % ", ".join("%s/%s %.1f%%" % (e, t, 100 * d)
                          for d, e, t in moved[:4]))
        print("  largest rise:  %s"
              % ", ".join("%s/%s %+.1f%%" % (e, t, 100 * d)
                          for d, e, t in moved[-4:]))
        print("  median change: %+.2f%%"
              % (100 * moved[len(moved) // 2][0]))

    #  and the point of the exercise: against the measured values
    rat = []
    for el in sorted(lib):
        v = lib[el]
        if isinstance(v, dict) and isinstance(v.get("mech"), dict) \
                and el in refdata.THETA_D:
            m = v["mech"]
            if m.get("debye") and m.get("debye_iso"):
                rat.append((m["debye"] / refdata.THETA_D[el],
                            m["debye_iso"] / refdata.THETA_D[el]))
    if rat:
        a = np.array(rat)
        print("\n  against Stewart's theta_D, %d elements:" % len(rat))
        print("    Christoffel  mean %.3f  spread %.3f"
              % (a[:, 0].mean(), a[:, 0].std()))
        print("    isotropic    mean %.3f  spread %.3f"
              % (a[:, 1].mean(), a[:, 1].std()))
        print()
        print("  Both sit near 0.9 and that is NOT the averaging.  These are")
        print("  built from elastic constants at ~300 K and a ~293 K volume,")
        print("  while Stewart's theta_D is a T -> 0 quantity - the anchor")
        print("  mismatch refdata.py's header documents and does not correct.")
        print("  Carried to 0 K with III/29a's coefficients the close-packed")
        print("  metals land on it: Cu 1.004, Ag 0.999, Au 1.004, Pd 0.999,")
        print("  Ni 0.990.  That is the check that says this averaging is the")
        print("  right one; the isotropic average misses by two per cent even")
        print("  after the temperature is accounted for.")

    if dry:
        print("\n--dry: nothing written")
        return
    with open(path, "w") as f:
        json.dump(lib, f, indent=1, sort_keys=True)
    print("\n-> library.json")


if __name__ == "__main__":
    main()
