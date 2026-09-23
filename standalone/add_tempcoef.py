#!/usr/bin/env python3
"""
Fold the measured temperature coefficients of the elastic constants into the
library, beside the model's own.

Writes lib[el]["tempcoef"] from standalone/tempcoef_compare.json: for C11, C12
and C44, the measured d(lnC)/dT against two model derivatives, all in 10^-4/K.

**There are two model derivatives and confusing them is the whole difficulty.**
The finite-temperature sweep holds the volume, so what it measures is the
constant-volume derivative, model_V.  An experiment does not hold the volume -
the sample expands - so what it reports also contains the elastic constants'
dependence on volume driven by thermal expansion.  model_P adds that term,
(dlnC/dlnV) * 3 alpha_L, using the MEASURED expansion coefficient rather than
the model's, so that the comparison tests the elastic constants and not the
expansion, which the library already reports separately.  model_P is the number
to compare with experiment; model_V is the one the sweep actually measured, and
the page shows both because the gap between them is the expansion term and is
worth seeing.

    python add_tempcoef.py
    python add_tempcoef.py --dry
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "tempcoef_compare.json")
LB = os.path.join(HERE, "lb29a_tempcoef.json")
#  hexagonal elements have five independent coefficients, not three; a list of
#  the cubic three silently drops Tc13 and Tc33, and thallium is the element
#  that shows it, because those two are the only ones it has
COEF = ("Tc11", "Tc12", "Tc13", "Tc33", "Tc44")


def main():
    dry = "--dry" in sys.argv
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    src = json.load(open(SRC, encoding="utf-8"))
    lb = json.load(open(LB, encoding="utf-8"))
    note = src.get("_note", "")
    source = src.get("_source", "")

    have, lack, n_cmp = [], [], 0
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or "a0" not in v:
            continue
        rec = src.get(el)
        if not isinstance(rec, dict):
            v.pop("tempcoef", None)
            lack.append(el)
            continue
        out = {k: rec[k] for k in COEF if isinstance(rec.get(k), dict)}
        if not out:
            #  an element with no comparable coefficient is kept, with the
            #  reason, rather than dropped: a row missing without explanation
            #  is the defect that a reader has no way to resolve.  Two reasons
            #  exist and they are different - the sweep had too few stable
            #  points, or the source marks every measured coefficient
            #  uncertain and a comparison against a doubtful number is not one
            m = lb.get(el) or {}
            unc = set(m.get("uncertain") or ())
            meas = [k for k in COEF if k in m]
            why = rec.get("skipped")
            if not why and meas and all(k in unc for k in meas):
                why = ("every coefficient the source lists for this element is "
                       "marked uncertain there")
            v["tempcoef"] = dict(none=True, why=why or "no comparable "
                                 "coefficient", system=rec.get("system"),
                                 note=note, source=source)
            lack.append(el)
            continue
        for c in out.values():
            if c.get("experiment") is not None and c.get("model_P") is not None:
                n_cmp += 1
        out.update(system=rec.get("system"), T_max=rec.get("T_max"),
                   T_points=rec.get("T_points"),
                   alpha_exp_1e6=rec.get("alpha_exp_1e6"),
                   note=note, source=source)
        v["tempcoef"] = out
        have.append(el)

    if not dry:
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(lib, fh, indent=1, sort_keys=True, default=str)
        os.replace(tmp, path)
    print("tempcoef written for %d elements, %d coefficients with both a "
          "measured and a constant-pressure model value; without any: %s%s"
          % (len(have), n_cmp, " ".join(lack) or "none",
             "  [dry]" if dry else ""))


if __name__ == "__main__":
    main()
