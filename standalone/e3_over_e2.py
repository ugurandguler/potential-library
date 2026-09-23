#!/usr/bin/env python3
"""Re-evaluate every shipped hard-cut record against fit.py's E3/E2 bound.

The bound - the three-body term may correct the pair term but may not cancel
it - was tightened for the tapered sets and never applied back to the hard-cut
library, which had already been fitted.  Nothing about those records changed;
the acceptance test in front of them did.  So re-running the fitter today does
not return the shipped hard-cut record for the ones that are above the bound,
and a reader reproducing the library needs to know which those are before they
read a version difference as a defect.

That count was measured once, on 2026-08-29, and typed into fit.py's comment,
into the manuscript and into the page.  It was typed with a denominator of 36,
because the comment's "within" list had dropped two elements - lutetium and
zirconium, which both comply - while its violator list was complete.  The
count of violators was right, the count of records was not, and three
documents carried it.

Hence this file: the ratio is computed from the shipped parameters rather than
remembered, and the three documents read the number from its output.

    python e3_over_e2.py        # writes e3_over_e2.json

The ratio is |E_coh + E_2| / |E_2| with E_2 the energy the same record gives
with C set to zero, which is how fit.py itself tests it.
"""
import io
import json
import os

import fit
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "e3_over_e2.json")


def ratios():
    """{element: E3/E2 of the shipped hard-cut record}"""
    lib = json.load(io.open(os.path.join(HERE, "library.json"), encoding="utf-8"))
    out = {}
    for el, d in lib.items():
        if not isinstance(d, dict) or not d.get("struct"):
            continue
        e = refdata.ELEMENTS[el]
        #  the pair-only energy of the same record: C = 0, everything else as
        #  shipped, scaled by the record's own D
        p = (d["m"], 1.0, d["alpha"], d["r0"], d["gamma"], 0.0, d["s3"])
        e2 = fit.constraints(el, e, p)["E"] * d["D"]
        out[el] = abs(-e["Ecoh"] - e2) / abs(e2)
    return out


def main():
    r = ratios()
    viol = sorted(el for el, x in r.items() if x > fit.E3_OVER_E2_MAX)
    out = {
        "bound": fit.E3_OVER_E2_MAX,
        "n": len(r),
        "n_violating": len(viol),
        "violating": viol,
        "within": sorted(set(r) - set(viol)),
        "ratio": {el: round(x, 4) for el, x in sorted(r.items())},
        "_source": "standalone/library.json through fit.constraints with C = 0",
    }
    json.dump(out, io.open(OUT, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"{out['n_violating']} of {out['n']} shipped hard-cut records are above "
          f"E3/E2 = {out['bound']}")
    print("  violating:", " ".join(viol))
    print("  within:   ", " ".join(out["within"]))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
