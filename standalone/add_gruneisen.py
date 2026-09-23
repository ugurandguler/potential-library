#!/usr/bin/env python3
"""
Fold the Grueneisen parameter and the quasi-harmonic expansion into the library.

Writes lib[el]["gruneisen"] from standalone/gruneisen_tap.json: the
heat-capacity-weighted mode Grueneisen parameter at 298.15 K, the value implied
by the measured expansion, and the three expansion coefficients that follow -
measured, quasi-harmonic prediction from gamma, and the model's own molecular
dynamics.

**The window scan travels with it and is the reason to trust the number.**
gamma is a derivative taken by finite difference, so it depends on the strain
window it was taken over; the source scans four windows from 0.5 % to 3 % and
records the spread.  A gamma whose spread across those windows is a few per
cent is a derivative; one that moves by tens of per cent, or whose modes change
sign inside the window, is an artefact of the step size.  Both are stored and
the page shows which it is, because a quasi-harmonic expansion computed from an
unconverged derivative looks exactly like a converged one.

    python add_gruneisen.py
    python add_gruneisen.py --dry
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "gruneisen_tap.json")
KEEP = ("gamma", "gamma_exp", "gamma_windows", "window_spread", "sign_changes",
        "alpha_V_meas_1e6", "alpha_V_pred_1e6", "alpha_V_model_1e6",
        "B_GPa", "Cv_298", "Vm_m3", "modes", "modes_used")


def main():
    dry = "--dry" in sys.argv
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    src = json.load(open(SRC, encoding="utf-8"))
    note = src.get("_note", "")

    have, lack, wide = [], [], []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or "a0" not in v:
            continue
        rec = src.get(el)
        if not isinstance(rec, dict) or rec.get("gamma") is None:
            v.pop("gruneisen", None)
            lack.append(el)
            continue
        out = {k: rec[k] for k in KEEP if k in rec}
        #  What disqualifies a value is gamma CHANGING SIGN between the 1 % and
        #  the 3 % window - then the quantity is not defined at this
        #  resolution.  The relative spread is stored because it is
        #  informative, but it is not the test and must not become one: gamma
        #  passes through zero for several elements, and a spread divided by a
        #  near-zero gamma reads as hundreds of per cent while the derivative
        #  itself is perfectly well behaved.  The manuscript qualifies the same
        #  five records on the same criterion, and the two must not disagree.
        g, s = out.get("gamma"), out.get("window_spread")
        if g and s is not None:
            out["window_spread_pct"] = 100.0 * s / abs(g)
        if out.get("sign_changes"):
            wide.append(el)
        out["note"] = note
        v["gruneisen"] = out
        have.append(el)

    if not dry:
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(lib, fh, indent=1, sort_keys=True, default=str)
        os.replace(tmp, path)
    print("gruneisen written for %d elements; for %d gamma changes sign "
          "between windows and the page says it is not defined there (%s); "
          "without a value: %s%s"
          % (len(have), len(wide), " ".join(wide) or "none",
             " ".join(lack) or "none", "  [dry]" if dry else ""))


if __name__ == "__main__":
    main()
