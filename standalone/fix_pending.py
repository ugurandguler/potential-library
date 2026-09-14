#!/usr/bin/env python3
"""
Recompute `lammps_pending` from what the record actually has.

add_candidates writes it as a fixed list of eight panel names when a candidate
record is first created, and nothing has updated it since - so every candidate
arm still claims that elasticT, Bain, the gamma surface and the expansion
"have not been run for this set", which is what the viewer then prints.  Those
runs finished and were merged, and those panels now draw the candidate arms,
so the sentence contradicts the plots directly above it.

The names are not all fields.  `gamma` is the generalized stacking fault, and
its result lives in `stacking` alongside it - `rec["gamma"]` is the three-body
parameter of the potential and would always look present.  elasticT lives at
lib[el]["elasticT"][arm] rather than on the arm.  Both are handled here rather
than by matching names to keys.

NOT APPLICABLE is not outstanding.  The generalized stacking fault and the
intrinsic stacking fault are fcc/hcp quantities: no bcc record in the library
has a `stacking` block, the published tapered arm included.  So for the eight
bcc candidates (Ba Cs K Li Mo Na Rb W) `gamma` and `stacking` could never be
found, and the page kept saying those runs had not been done - for runs that
are not defined for that lattice.  They are dropped for bcc, not marked done.

The key is REMOVED when nothing is outstanding, so the viewer's own
`P.lammps_pending ? ... : ""` says nothing at all.

    python fix_pending.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

#  panel name -> how to find its result for arm `a` of element record `v`
WHERE = {
    "elasticT":  lambda v, a, r: a in (v.get("elasticT") or {}),
    "bain":      lambda v, a, r: bool(r.get("bain")),
    "gamma":     lambda v, a, r: bool(r.get("stacking")),
    "stacking":  lambda v, a, r: bool(r.get("stacking")),
    "expansion": lambda v, a, r: bool(r.get("expansion")),
    "surface":   lambda v, a, r: bool(r.get("surface")),
    "md_screen": lambda v, a, r: bool(r.get("md_screen")),
    "jiggle":    lambda v, a, r: bool(r.get("jiggle")),
}


NOT_FOR_BCC = ("gamma", "stacking")


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    cleared, left = 0, {}
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict):
            continue
        for arm in ("rc", "rc_ug"):
            r = v.get(arm)
            if not isinstance(r, dict) or not r.get("lammps_pending"):
                continue
            na = NOT_FOR_BCC if v.get("struct") == "bcc" else ()
            still = [k for k in r["lammps_pending"]
                     if k in WHERE and k not in na and not WHERE[k](v, arm, r)]
            if still:
                r["lammps_pending"] = still
                left.setdefault(", ".join(still), []).append(f"{el}/{arm}")
            else:
                del r["lammps_pending"]
                cleared += 1
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"{cleared} candidate arms had nothing outstanding; the key is gone")
    for what, els in sorted(left.items()):
        print(f"  still pending [{what}]: {len(els)} - {' '.join(els)}")


if __name__ == "__main__":
    main()
