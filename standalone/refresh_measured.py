#!/usr/bin/env python3
"""
Re-attach every MEASURED field to the library after a refit.

A record in library.json holds two different kinds of thing.  There are the
fitted parameters, which come out of the search, and there are measurements
*of* those parameters - whether the crystal survives 600 K, where its lowest
phonon sits, how far it can be compressed before the lattice stops being the
ground state.  The second kind is only true of the exact numbers it was
measured from.

So a refit invalidates all of it, and the merge that brings new parameters in
deliberately drops `md`, `md_screen`, `min_cm1`, `stable`, `R` and the rest.
That is the safe direction to fail - a missing field renders as "not screened",
a stale one renders as a verdict about a potential that no longer exists.  But
it leaves the page with empty sections until something puts them back, and
until now nothing did: these fields were being written by hand from a scratch
snippet, which is exactly how `tgt["md"] = {...}` once wiped the taper
narrative out of 38 elements without raising an error.

This is that missing step, written down.  It reads the measurement files that
already exist and folds them in, key by key, never replacing a record.

    python refresh_measured.py            # every element, every set
    python refresh_measured.py --dry      # say what would change, write nothing
"""
import json
import os
import sys

import numpy as np

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

#  (library key, pair style, filename suffix)
SETS = {
    "hard": (None, "ugur", ".ugur"),
    "tap": ("tap", "ugur", "_taper.ugur"),
    "ug": ("ug", "ugur/ang", ".ugur.ang"),
    "tap_ug": ("tap_ug", "ugur/ang", "_taper.ugur.ang"),
}
#  the angular sets need the angular tree's latdyn, which cannot be imported
#  here - the two trees define incompatible versions.  angular/dynstab_ug.py
#  writes them out instead; absent that file the fields are simply left alone
#  rather than being filled in with a number from the wrong engine.
UGDYN = os.path.join(ROOT, "angular", "dynstab_ug.json")


def verdict(r):
    if r is None:
        return None
    if r.get("lost"):
        return "DAGILDI"
    if r.get("collapsed"):
        return "COKTU"
    if r.get("T", 300) > 400:
        return "SUPHELI"
    return "saglam"


def main():
    dry = "--dry" in sys.argv
    lib = json.load(open(os.path.join(HERE, "library.json")))
    scr = {}
    p = os.path.join(ROOT, "lammps", "md_screen_all.json")
    if os.path.exists(p):
        scr = json.load(open(p))
    comp = {}
    p = os.path.join(HERE, "compression.json")
    if os.path.exists(p):
        comp = json.load(open(p))
    ugdyn = json.load(open(UGDYN)) if os.path.exists(UGDYN) else {}
    jig = {}
    p = os.path.join(ROOT, "lammps", "jiggle_test.json")
    if os.path.exists(p):
        jig = json.load(open(p))

    n = {"md_screen": 0, "R": 0, "stable": 0, "compression": 0, "jiggle": 0}
    for el in sorted(lib):
        e = refdata.ELEMENTS[el]
        for name, (key, style, suffix) in SETS.items():
            rec = lib[el] if key is None else lib[el].get(key)
            if not rec or "m" not in rec:
                continue

            #  --- the MD verdict -------------------------------------------
            s = scr.get(f"{el}|{name}")
            if s is not None:
                rec["md_screen"] = {
                    "T": s.get("T"), "collapsed": bool(s.get("collapsed")),
                    "lost": s.get("lost"), "state": verdict(s),
                    "style": style, "file": el + suffix,
                }
                n["md_screen"] += 1

            #  --- anisotropy, cubic only ----------------------------------
            cij = rec.get("Cij") or {}
            ce = e["Cij"]
            if e["struct"] in ("fcc", "bcc") and "C11" in cij:
                cp = 0.5 * (cij["C11"] - cij["C12"])
                rec["R"] = cij["C44"] / cp if cp > 0 else None
                rec["R_exp"] = ce["C44"] / (0.5 * (ce["C11"] - ce["C12"]))
                n["R"] += 1

            #  --- lowest phonon -------------------------------------------
            if name in ("hard", "tap"):
                try:
                    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                                    mass=refdata.MASSES[el])
                    fr = L.spectrum(cry, L.Potential.from_record(rec), 8)
                    #  THz -> cm^-1, signed: an imaginary mode is the point
                    mn = float(np.nanmin(np.asarray(fr))) * 33.35641
                    rec["stable"] = bool(mn > -1e-3)
                    rec["min_cm1"] = round(mn, 3)
                    n["stable"] += 1
                except Exception:
                    pass
            else:
                d = ugdyn.get(f"{el}|{name}")
                if d:
                    rec["stable"] = bool(d["stable"])
                    rec["min_cm1"] = d["min_cm1"]
                    n["stable"] += 1

            #  --- finite-amplitude stability -------------------------------
            #  The phonon screen above answers a question at zero amplitude.
            #  This one asks whether the lattice survives a 1e-5 A nudge, which
            #  five bcc records pass the first and fail the second.
            j = jig.get(f"{el}|{name}")
            if j is not None:
                rec["jiggle"] = {"keep": j["keep"], "dE": j["dE"],
                                 "ok": bool(j["ok"])}
                n["jiggle"] += 1

            #  --- compression escape --------------------------------------
            if f"{el}|{name}" in comp:
                c = comp[f"{el}|{name}"]
                rec["compression"] = ({"basin": False} if c is None else
                                      {"basin": True, "x": c["x"],
                                       "barrier": c["barrier"],
                                       "depth": c["depth"],
                                       "reachable": c["reachable"]})
                n["compression"] += 1

    print(f"yeniden yazilan alanlar: {n}")
    if dry:
        print("--dry: nothing written to the file")
        return
    json.dump(lib, open(os.path.join(HERE, "library.json"), "w"),
              indent=1, sort_keys=True)
    print("-> library.json")


if __name__ == "__main__":
    main()
