#!/usr/bin/env python3
"""
The symmetry-path stability screen, for the arms that never had one.

add_dynstab does this for the top-level MAU record and add_dynstab_ug for the
UG one.  Both of those are HARD-CUT arms.  The tapered arms - `tap` and
`tap_ug`, the ones the library actually ships for molecular dynamics - carry
only the check the fit ran on itself, and iridium is the standing proof that
that check is not enough: it reports +6.1 cm^-1 for a record whose symmetry
path dips to -8.1.  A mode confined to a line is invisible to a half-shifted
Monkhorst-Pack mesh however fine, and refining steps over the same places.

So the shipped arms were the least screened ones on the page.

    python add_dynstab_arms.py tap        # pair latdyn
    python add_dynstab_arms.py tap_ug     # ANGULAR latdyn
    python add_dynstab_arms.py rc         # and the candidates,
    python add_dynstab_arms.py rc_ug      # so every arm is screened alike

The two cannot run in one process - both trees define latdyn and only one
carries the Legendre terms - so the arm is chosen before the import and the
module that gets loaded is asserted against it.  Writes `dyn` into the arm's
own record, in the same shape add_dynstab writes.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ARM = (sys.argv[1] if len(sys.argv) > 1 else "tap")
if ARM not in ("tap", "tap_ug", "rc", "rc_ug"):
    raise SystemExit("arm must be tap, tap_ug, rc or rc_ug")
ANG = ARM.endswith("_ug")
if ANG:
    sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "angular")))
else:
    sys.path.insert(0, HERE)

import latdyn as L                                   # noqa: E402
import refdata                                       # noqa: E402
from build_library import sc_segments, NQ_PATH       # noqa: E402

src = open(L.__file__, encoding="utf-8").read()
assert ("lam2" in src) == ANG, \
    f"wrong latdyn for arm {ARM}: {L.__file__}"

IMAG_TOL_CM1 = -1.0


def main(meshes=(8, 9)):
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print(f"screening arm {ARM} with {L.__file__}")
    print(f"{'el':4s}{'mesh':>9s}{'path':>9s}   was     now")
    print("-" * 42)
    changed, bad = [], []
    for el in sorted(lib):
        v = lib[el]
        rec = v.get(ARM) if isinstance(v, dict) else None
        if not isinstance(rec, dict) or "m" not in rec:
            continue
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(rec)
        f_mesh = np.concatenate([L.spectrum(cry, pot, nq=n).ravel()
                                 for n in meshes]) * L.CM1
        qs = []
        for (_, ka, __, kb) in sc_segments(e["struct"]):
            qs.extend(ka + np.linspace(0, 1, NQ_PATH)[:, None] * (kb - ka))
        f_path = L.frequencies_many(cry, pot, np.array(qs)).ravel() * L.CM1
        f = np.concatenate([f_mesh, f_path])
        mn_m, mn_p = float(f_mesh.min()), float(f_path.min())
        ok = bool(f.min() >= IMAG_TOL_CM1)
        was = rec.get("stable")
        rec["dyn"] = {"stable": ok, "min_mesh_cm1": round(mn_m, 3),
                      "min_path_cm1": round(mn_p, 3),
                      "most_neg_cm1": round(float(f.min()), 3),
                      "imag_frac": round(float((f < IMAG_TOL_CM1).mean()), 5),
                      "nq": "+".join(map(str, meshes)),
                      "path": "Setyawan-Curtarolo"}
        print(f"{el:4s}{mn_m:9.2f}{mn_p:9.2f}{str(was):>7s}{str(ok):>8s}")
        if not ok:
            bad.append((el, round(mn_p, 1)))
        if was is not None and bool(was) != ok:
            changed.append(f"{el} {was}->{ok}")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    n = sum(1 for v in lib.values()
            if isinstance(v, dict) and isinstance(v.get(ARM), dict)
            and "dyn" in v[ARM])
    print()
    print(f"{n - len(bad)}/{n} stable on mesh and path")
    for el, x in bad:
        print(f"   {el}: down to {x} cm-1 on the path")
    if changed:
        print("verdict changed for " + ", ".join(changed))


if __name__ == "__main__":
    main()
