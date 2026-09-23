#!/usr/bin/env python3
"""
Fold the self-interstitial formation energies into the library.

Writes lib[el]["interstitial"] from lammps/interstitial_tap.json: the formation
energy of the lowest dumbbell, which orientation that was, and the energy of
each orientation tried, so the page can show the ordering and not only the
winner.  The ordering is the interesting part - it is what a potential has to
get right to describe radiation damage, and a pair potential has no reason to.

**Twelve elements are absent and they are exactly the twelve hexagonal ones.**
The three orientations measured, <100>, <110> and <111>, are cubic directions;
the hexagonal close-packed lattice has its own interstitial sites (octahedral,
tetrahedral, crowdion, and two distinct basal dumbbells) and attempting it with
a cubic site list would report the lowest of the wrong set.  The source file
says so itself and the count is checked here rather than trusted: if a
hexagonal element ever appears in the source, or a cubic one goes missing, this
script says which.

    python add_interstitial.py
    python add_interstitial.py --dry
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "lammps", "interstitial_tap.json")


def main():
    dry = "--dry" in sys.argv
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    src = json.load(open(SRC, encoding="utf-8"))
    note = src.get("_note", "")

    have, lack, odd = [], [], []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or "a0" not in v:
            continue
        rec = src.get(el)
        if not isinstance(rec, dict):
            v.pop("interstitial", None)
            lack.append(el)
            continue
        if v.get("struct") == "hcp":
            odd.append(el)
        v["interstitial"] = dict(
            E_I=rec["E_I"], site=rec["site"], per_site=rec["per_site"],
            set=rec.get("set", "tap"), natoms=rec.get("natoms"), note=note)
        have.append(el)

    #  the absence has to be explainable, so it is checked, not assumed
    wrong = [e for e in lack if lib[e].get("struct") != "hcp"]
    if wrong or odd:
        raise SystemExit(
            "the cubic-only rule no longer describes the source: missing and "
            "not hexagonal %s; present and hexagonal %s" % (wrong, odd))

    if not dry:
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(lib, fh, indent=1, sort_keys=True, default=str)
        os.replace(tmp, path)
    print("interstitial written for %d elements; the %d without it are the "
          "hexagonal ones (%s)%s"
          % (len(have), len(lack), " ".join(lack), "  [dry]" if dry else ""))


if __name__ == "__main__":
    main()
