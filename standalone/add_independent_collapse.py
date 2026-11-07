#!/usr/bin/env python3
"""
Collapses the library's own screens do not show, found by independent tests.

The library judges a record on its fitted structure: the phonons, the
compression scan, the MD screen, the box-free structure test.  An independent
ground-state test relaxes seven prototypes with the cell free, and the vacancy
test relaxes the cell around a vacancy.  For a few switched records those runs
fall into energies no metal has - lithium's A15 arrangement 339 eV/atom below
its bcc, chromium's switched UG crystal 85 eV down once the cell may relax
round a vacancy - or do not complete.  The record should say so where a reader
choosing it will look.

collapse_independent.json carries the numbers (with the criteria in its
_source); this writes them as `independent_collapse` into the record
of the set they belong to, and nothing else.

    python add_independent_collapse.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = {"ug", "tap", "tap_ug", "rc", "rc_ug", "tap_force", "hard_disp", "tap_nudge"}


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    src = json.load(open(os.path.join(HERE, "collapse_independent.json")))
    #  cleared first, so a record taken off the list loses its note
    for v in lib.values():
        for r in ([v] + list(v.values())) if isinstance(v, dict) else []:
            if isinstance(r, dict):
                r.pop("independent_collapse", None)
    n, absent = 0, []
    for key, v in src.items():
        if key.startswith("_"):
            continue
        el, arm = key.split("|", 1)
        if el not in lib or arm not in ARMS:
            raise SystemExit(f"{key}: no such record")
        rec = lib[el].get(arm)
        #  a clean clone carries no switched or force-matched arm, so the
        #  records these notes belong to are not there to write into
        if not isinstance(rec, dict):
            absent.append(key)
            continue
        rec["independent_collapse"] = v
        n += 1
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"independent_collapse written for {n} records"
          + (f"; {len(absent)} not in this checkout" if absent else ""))


if __name__ == "__main__":
    main()
