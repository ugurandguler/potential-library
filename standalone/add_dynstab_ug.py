#!/usr/bin/env python3
"""
Dynamical stability of the UG (angular) records, on the path as well as the
mesh.

`angular/export_ug.py` decides `ug.dyn.stable` from the mesh minimum and a
sample near Gamma.  Neither can see an instability confined to a symmetry line
- a half-shifted Monkhorst-Pack grid never lands on one, and refining it steps
over the same places - so aluminium and sodium are recorded stable while the
dispersion drawn for them on the same page dips to -3.4 and -1.5 cm^-1.  The
picture was right and the badge was wrong, which is the same defect iridium
had in the MAU arm.

Force constants come from angfc explicitly, as export_ug does: latdyn
delegates there when the Legendre coefficients are non-zero, but a silent
fallback to the angle-free constants would call an unstable fit stable.

    python add_dynstab_ug.py [nq ...]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "angular")))

import latdyn as L          # noqa: E402  - the ANGULAR one; asserted below
import angfc                # noqa: E402
import refdata              # noqa: E402
from build_library import sc_segments, NQ_PATH   # noqa: E402

#  the two latdyn modules differ in exactly the field that matters here
assert "lam2" in open(L.__file__, encoding="utf-8").read(), \
    f"wrong latdyn on the path: {L.__file__} has no angular terms"

IMAG_TOL_CM1 = -1.0


def main(meshes=(8, 9)):
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print(f"{'el':4s}{'mesh':>9s}{'path':>9s}   was    now")
    print("-" * 40)
    changed, bad = [], []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or not isinstance(v.get("ug"), dict):
            continue
        u = v["ug"]
        if "m" not in u or el not in refdata.ELEMENTS:
            continue
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(u)
        Phi = angfc.force_constants(cry, pot)
        qs = []
        for (_, ka, __, kb) in sc_segments(e["struct"]):
            qs.extend(ka + np.linspace(0, 1, NQ_PATH)[:, None] * (kb - ka))
        fm = np.concatenate([L.spectrum(cry, pot, nq=n).ravel()
                             for n in meshes]) * L.CM1
        fp = L.frequencies_many(cry, pot, np.array(qs), Phi).ravel() * L.CM1
        f = np.concatenate([fm, fp])
        neg = int((f < IMAG_TOL_CM1).sum())
        was = u.get("dyn", {}).get("stable")
        d = dict(u.get("dyn") or {})
        d.update({"min_mesh_cm1": round(float(fm.min()), 3),
                  "min_path_cm1": round(float(fp.min()), 3),
                  "imag_frac": round(neg / f.size, 4),
                  "most_neg_cm1": round(float(f.min()), 1) if neg else 0.0,
                  "stable": neg == 0,
                  "nq": "+".join(map(str, meshes)),
                  "path": "Setyawan-Curtarolo"})
        u["dyn"] = d
        #  `stable` means the same thing on every record: mesh AND path.
        #  refresh_measured.py used to leave the mesh-only verdict here, which
        #  is how Al and Na read stable beside a dispersion dipping below zero.
        u["stable"] = neg == 0
        if bool(was) != (neg == 0):
            changed.append((el, was, neg == 0))
        if neg:
            bad.append((el, float(fp.min())))
        print(f"{el:4s}{fm.min():9.2f}{fp.min():9.2f}   {str(was):5s}  "
              f"{str(neg == 0):5s}{'  <- CHANGED' if bool(was) != (neg == 0) else ''}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    n = sum(1 for el in lib if isinstance(lib.get(el), dict)
            and isinstance(lib[el].get("ug"), dict)
            and "dyn" in lib[el]["ug"])
    print(f"\n{n - len(bad)}/{n} dynamically stable on mesh and path")
    for el, m in sorted(bad, key=lambda r: r[1]):
        print(f"   {el}: down to {m:.1f} cm-1 on the path")
    print(f"verdict changed for {len(changed)}: "
          + ", ".join(f"{e} {a}->{b}" for e, a, b in changed))


if __name__ == "__main__":
    main(tuple(int(a) for a in sys.argv[1:]) or (8, 9))
