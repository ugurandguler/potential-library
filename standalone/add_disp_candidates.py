#!/usr/bin/env python3
"""
Carry the dispersion-selected candidates BESIDE the fitted records.

Same pattern as the nudge-constrained refits `add_structure.py` writes into
`tap_nudge`, and for the same reason: a solution the objective cannot
distinguish from the shipped one, chosen on a criterion the objective cannot
see, is worth having in front of a reader - and picking one silently would hide
the only interesting thing about it.  Both are shipped, both are labelled.

Written under `hard_disp`: the hard-cut arm, selected by dispersion.

WHAT THEY ARE.  Fits made with today's `fit.py` to the UNCHANGED targets in
refdata.py.  They reach the same 0.000 % elastic residual as the shipped
records - the objective genuinely cannot tell them apart - and a lower error
against the measured phonon dispersion, which never enters the fit:

    Pb  18.30 % -> 8.00 %      Cu  8.50 % -> 6.30 %
    Ag  13.10 % -> 11.90 %     Pd  3.50 % -> 3.20 %

WHY THEY EXIST AT ALL.  Because `fit.py` enforces two constraints today that
the shipped records predate, and both were tightened for the TAPERED sets and
never applied back to the hard-cut library: `E3_OVER_E2_MAX = 0.30` (22 of the
36 shipped records are above it - see the note beside that constant) and
`REQUIRE_COMPRESSION_STABLE`.  All four shipped records here fail both.  The
candidates satisfy both:

    element   E3/E2            compression guard
    Pb        0.676 -> 0.298   FAIL -> PASS
    Cu        0.531 -> 0.248   FAIL -> PASS
    Ag        0.487 -> 0.281   FAIL -> PASS
    Pd        0.496 -> 0.299   FAIL -> PASS

So this is not a new objective.  It is what today's objective returns for the
same targets, and it happens to describe the phonons better.

WHAT WAS CHECKED, AND WHAT WAS NOT.  All four are dynamically stable on the
8^3 U 9^3 union plus the Gamma neighbourhood (fit.evaluate enforces exactly
that mesh), put fcc below bcc and hcp, have no compression basin on
compression.py's own barrier test, and survive the MD screen - where they sit
FURTHER above their static lattice energy than the shipped records do in all
four cases, and closer to 300 K.  Lead is the sharpest: the shipped record ends
0.005 eV/atom BELOW its own static lattice (inside the 0.05 tolerance, so not a
collapse) and runs to 317 K, while the candidate sits 0.035 above it at 300 K.

NOT checked, and it is not an omission: vacancy, surface and stacking-fault
energies.  These are hard-cut records, phi2 does not vanish at the cutoff, and
a relaxation walks into the discontinuity - which is why this library carries
vacancy_tap.json and vacancy_tap_ug.json and no hard-cut equivalent.

    python add_disp_candidates.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

WHY = ("same targets and the same 0.000 % elastic residual as the fitted "
       "record, a lower error against the measured dispersion, and inside "
       "the E3/E2 and compression constraints that fit.py enforces today "
       "and the fitted record fails")
FROM = ("dense_fit.py at 32 restarts under today's constraints; selected "
        "among tied solutions by curve_mae against the measured dispersion")


def main():
    lib_path = os.path.join(HERE, "library.json")
    with open(lib_path) as f:
        lib = json.load(f)
    with open(os.path.join(HERE, "disp_candidates.json")) as f:
        cand = json.load(f)

    n = 0
    for el, rec in sorted(cand.items()):
        if el not in lib or not isinstance(lib[el], dict):
            print("  %s: not in the library" % el)
            continue
        out = dict(rec)
        out["from"] = FROM
        out["why"] = WHY
        lib[el]["hard_disp"] = out
        n += 1
        ev = rec.get("evidence", {})
        print("  %-3s dispersion %.2f -> %.2f %%   E3/E2 %.3f -> %.3f   "
              "rms %.4f %%"
              % (el, ev.get("dispersion_published", float("nan")),
                 ev.get("dispersion_here", float("nan")),
                 ev.get("E3_over_E2_published", float("nan")),
                 ev.get("E3_over_E2_here", float("nan")),
                 rec.get("rms", float("nan"))))

    with open(lib_path, "w") as f:
        json.dump(lib, f, indent=1, sort_keys=True)
    print("\n%d candidates carried beside the fitted records, under "
          "`hard_disp`" % n)
    print("Nothing is replaced.  The fitted record is untouched and is still "
          "what every other")
    print("field in these entries refers to.")


if __name__ == "__main__":
    main()
