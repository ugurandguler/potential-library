#!/usr/bin/env python3
"""
Merge the tapered UG fits into library.json, beside the tapered MAU ones.

These eight are the metals the hard cutoff put out of reach, and they are the
reason the angular factor exists at all: λ₂ and λ₄ were introduced to reach a
C44/C′ the plain φ₂+φ₃ form cannot produce.  The taper reaches them without any
angular term, so the obvious question was whether the angular factor had been
made redundant.  It has not, and the answer is only visible with both runs in
the same record:

    taper alone     Fe, Nb, Ta, V exact; Al 8.89, W 8.79, Cr 14.62, Mo 16.85
    + angular       Al, W, Cr also exact; Mo 16.46 and still the odd one out

so the two mechanisms are complementary rather than alternatives.  What makes
that a statement about the functional form rather than about the search is the
fitted weights, which are stored here for exactly that reason: every element
that gained carries a large λ, and every element that did not has λ driven to
zero by the search itself.

Stored under `tap_ug`.  It does not replace anything - `tap` is the tapered MAU
fit and remains what the MD parameters are built from, since `pair_style ugur`
implements the published φ₃ and not the Legendre factor.

    python add_ug_taper.py <ug_tap.json>
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
#  below this the weights are off, not small: the search turned them off rather
#  than fitting noise into them, which is the check worth recording per element
LAM_OFF = 0.01


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    new = json.load(open(sys.argv[1]))
    #  The control is the ANGULAR=0 arm of the same job, not library.json's
    #  `tap`.  That one comes from the standalone tree, fitted with the
    #  analytic force constants and the stability screen; the angular tree
    #  differentiates numerically and cannot use that screen.  Comparing across
    #  the two would mix the angular term's contribution with the difference
    #  between the trees, so the control travels with the record.
    CTRL = {"Al": 0.0889, "Cr": 0.1462, "Fe": 0.0, "Mo": 0.1685,
            "Nb": 0.0, "Ta": 0.0, "V": 0.0, "W": 0.0879}
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))

    print(f"{'el':4s}{'taper only':>11s}{'+ angular':>10s}{'gain':>9s}"
          f"{'lam2':>9s}{'lam4':>9s}   angular term")
    print("-" * 62)
    n = 0
    for el, rec in sorted(new.items()):
        v = lib.get(el)
        if not v:
            print(f"{el:4s}  not in the library, skipped")
            continue
        before = CTRL.get(el)
        after = rec["score"]
        rec = dict(rec)
        rec["lam_off"] = (abs(rec.get("lam2", 0.0)) < LAM_OFF
                          and abs(rec.get("lam4", 0.0)) < LAM_OFF)
        rec["score_ctrl"] = before
        v["tap_ug"] = rec
        n += 1
        gain = (100 * (before - after)) if before is not None else float("nan")
        what = ("kapali - gerekmedi" if rec["lam_off"]
                else "ACIK - farki bu yapiyor" if gain > 0.5
                else "open but no gain")
        b = f"{100*before:11.2f}" if before is not None else f"{'?':>11s}"
        print(f"{el:4s}{b}{100*after:10.2f}{gain:+9.2f}"
              f"{rec.get('lam2', 0.0):9.3f}{rec.get('lam4', 0.0):9.3f}   {what}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\ntapered UG folded into {n} elements -> {path}")


if __name__ == "__main__":
    main()
