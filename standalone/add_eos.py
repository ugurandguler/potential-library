#!/usr/bin/env python3
"""
Fold the equation-of-state scan into the library.

Writes lib[el]["eos"] from lammps/eos.json and standalone/eos_ref.json: for
each parameter set, the Birch-Murnaghan fit of energy against volume over the
scanned span - V0, B0, and the pressure derivative B0' - beside the measured
values.

**B0' is the point of this.** The bulk modulus is a fitted target and agreeing
with it says nothing; its pressure derivative is not fitted anywhere, so it is
the first number on the page that the fit was never shown.  A published
potential is scanned in the same loop where one exists, so the comparison is
against something and not against nothing.

The fit residual travels with every set and matters more here than usual: the
hard-truncated sets step at the cut-off, the step lands inside the scanned
volume range, and a fit through it returns a B0' that is a property of the
discontinuity rather than of the metal.  Those records show a residual two to
three orders of magnitude above the switched ones, and the page says so rather
than drawing them as though they were comparable.

    python add_eos.py
    python add_eos.py --dry
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCAN = os.path.join(ROOT, "lammps", "eos.json")
REF = os.path.join(HERE, "eos_ref.json")
KEEP = ("B0_GPa", "Bp", "V0", "E0", "rms_meV", "n", "span", "B_target_GPa")


def main():
    dry = "--dry" in sys.argv
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    scan = json.load(open(SCAN, encoding="utf-8"))
    ref = json.load(open(REF, encoding="utf-8"))

    #  "El|arm" for our own sets, "El|base|file" for a published potential
    by_el = {}
    for key, rec in scan.items():
        if key.startswith("_"):
            continue
        parts = key.split("|")
        el, arm = parts[0], parts[1]
        slot = by_el.setdefault(el, dict(sets={}, published={}))
        kept = {k: rec[k] for k in KEEP if k in rec}
        if arm == "base":
            slot["published"][parts[2] if len(parts) > 2 else "published"] = kept
        else:
            slot["sets"][arm] = kept

    have, lack, refd = [], [], []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or "a0" not in v:
            continue
        slot = by_el.get(el)
        if not slot or not slot["sets"]:
            v.pop("eos", None)
            lack.append(el)
            continue
        out = dict(sets=slot["sets"], published=slot["published"])
        r = ref.get(el)
        if isinstance(r, list) and len(r) == 4:
            out["ref"] = dict(V0=r[0], Ecoh=r[1], B0_GPa=r[2], Bp=r[3])
            out["ref_source"] = ref.get("_source", "")
            refd.append(el)
        v["eos"] = out
        have.append(el)

    if not dry:
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(lib, fh, indent=1, sort_keys=True, default=str)
        os.replace(tmp, path)
    print("eos written for %d elements, %d of them with a measured B0' to "
          "compare against; without a scan: %s%s"
          % (len(have), len(refd), " ".join(lack) or "none",
             "  [dry]" if dry else ""))


if __name__ == "__main__":
    main()
