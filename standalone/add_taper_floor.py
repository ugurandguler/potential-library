#!/usr/bin/env python3
"""
Merge the tapered C44/C' floors into library.json, beside the hard-cutoff ones.

The shipped library is truncated hard and its floor stands as the verdict the
page reports.  These are the same measurement with the cutoff switched off
smoothly, and they are stored alongside rather than replacing it, because the
comparison is the finding.

One distinction has to survive into the record, or the numbers are misread.
For the bcc metals the floor is a measurement: chromium reaches 0.362 with
C44/C11 = 0.078, well clear of any constraint.  For the fcc metals it is not.
Silver, copper and platinum all report 0.023, identical to three figures, and
the reason is that all three sit exactly on `refine_R.SHEAR_MIN_FRAC`: C44 is
1.00 per cent of C11 and C'/C11 is 0.437, so R = 0.01/0.437 = 0.023 whatever
the element.  That is where the search was told to stop, not where the form
stops.  The guard exists for a good reason - a crystal with less shear
resistance than that is not a potential anyone can use - but a number produced
by it is a bound on the search, and saying otherwise would be a claim we have
not measured.

    python add_taper_floor.py <dir with R_floor_*.json>
"""
import glob
import json
import os
import sys

import fit as F
import refdata
import refine_R as R

HERE = os.path.dirname(os.path.abspath(__file__))
GUARD_TOL = 1.05          # within 5 % of the guard counts as sitting on it


def guard_limited(el, rec):
    """is this floor the shear guard rather than the form's own limit?"""
    e = refdata.ELEMENTS[el]
    try:
        m, g, s, C = rec["argmin"]
        F.TAPER = rec.get("taper", 0.85)
        o = F.evaluate(el, e, m, g, C, s)
    except Exception:
        return None
    if not o:
        return None
    c = o["Cij"]
    cp = 0.5 * (c["C11"] - c["C12"])
    frac = min(c["C44"], cp) / c["C11"]
    return bool(frac < GUARD_TOL * R.SHEAR_MIN_FRAC)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    src = sys.argv[1]
    new = {}
    for f in glob.glob(os.path.join(src, "*.json")):
        new.update(json.load(open(f)))
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))

    print(f"{'el':4s}{'R meas':>8s}{'hard':>9s}{'taper':>9s}"
          f"{'verdict':>12s}   baseline?")
    print("-" * 58)
    n = flips = 0
    for el, rec in sorted(new.items()):
        v = lib.get(el)
        if not v or not v.get("reach"):
            continue
        ok = bool(rec["R_exp"] >= rec["R_floor"])
        gl = guard_limited(el, rec)
        v["reach"]["R_floor_taper"] = round(rec["R_floor"], 3)
        v["reach"]["ok_taper"] = ok
        v["reach"]["taper_floor_is_guard"] = gl
        n += 1
        if v["reach"].get("ok") is False and ok:
            flips += 1
        what = ("olculdu" if gl is False
                else "KISIT SINIRI" if gl else "belirsiz")
        print(f"{el:4s}{rec['R_exp']:8.2f}{v['reach']['R_floor']:9.3f}"
              f"{rec['R_floor']:9.3f}"
              f"{('erisilir' if ok else 'FORM DISI'):>12s}   {what}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\ntapered floor folded into {n} elements, {flips} verdicts changed")
    print("birlestirildi:", path)


if __name__ == "__main__":
    main()
