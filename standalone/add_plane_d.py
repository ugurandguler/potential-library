#!/usr/bin/env python3
"""
Add the (1 -1 0) section to every arm that already has mechanical planes.

The panel drew x-y, x-z and y-z.  Between them those contain [100], [110] and
[101] and no member of <111>, and <111> is where a cubic crystal's Young's
modulus is extremal - so the panel was missing the stiffest direction of every
cubic metal in the library.  Measured, not asserted: copper's three sections
span 66.7 to 130.3 GPa while E along [111] is 191.1, iridium's 329.6 to 526.6
against 657.8, and sodium's 1.9 to 5.0 against 10.6, less than half.

The fit is not re-run.  Each arm's own tensor is used: `Cfull` where the
record carries it, and for the candidate arms - which keep only the
independent constants - the tensor is rebuilt from those, which reproduces the
stored analyses to 0.0000 GPa.

A hexagonal cell does not need this plane: its surface is a body of revolution
about c and the x-z and y-z sections already agree to 0.0000 GPa.  It is
computed anyway, so every record carries the same four and the viewer needs no
per-symmetry case.

    python add_plane_d.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import mechanics                              # noqa: E402
import refdata                                # noqa: E402
from add_mech_recut import tensor             # noqa: E402


def tensor_of(el, rec, struct):
    if rec.get("Cfull"):
        return np.array(rec["Cfull"], dtype=float)
    c = rec.get("Cij")
    if not isinstance(c, dict):
        return None
    need = (("C11", "C12", "C13", "C33", "C44") if struct == "hcp"
            else ("C11", "C12", "C44"))
    if any(k not in c for k in need):
        return None
    return tensor(struct, c)


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    n, skipped = 0, []
    print(f"{'el':4s}{'arm':10s}{'3 kesitin en buyugu':>21s}{'yeni kesit max':>16s}"
          f"{'kazanc':>9s}")
    print("-" * 60)
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict):
            continue
        struct = refdata.ELEMENTS[el]["struct"]
        for name, rec in (("MAU", v), ("UG", v.get("ug")),
                          ("MAU re-cut", v.get("rc")),
                          ("UG re-cut", v.get("rc_ug"))):
            if not isinstance(rec, dict) or not rec.get("mech_planes"):
                continue
            if "d" in rec["mech_planes"]:
                continue
            C = tensor_of(el, rec, struct)
            if C is None:
                skipped.append(f"{el}/{name}")
                continue
            try:
                rec["mech_planes"]["d"] = mechanics.plane_curves(C, "d")
            except np.linalg.LinAlgError:
                skipped.append(f"{el}/{name}")
                continue
            n += 1
            if name == "MAU" and struct != "hcp":
                old = max(x for p in ("xy", "xz", "yz")
                          for x in rec["mech_planes"][p]["E"])
                new = max(rec["mech_planes"]["d"]["E"])
                print(f"{el:4s}{name:10s}{old:21.1f}{new:16.1f}"
                      f"{100 * (new - old) / old:8.0f}%")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print()
    print(f"{n} arm records gained the (1 -1 0) section"
          + (f"; no tensor for {' '.join(skipped)}" if skipped else ""))


if __name__ == "__main__":
    main()
