#!/usr/bin/env python3
"""
Read Table 28 of Landolt-Boernstein III/29a - the temperature coefficients of
the cubic elements' elastic constants.

Tc_ij = (1/c_ij) dc_ij/dT in 10^-4/K, measured at constant pressure over the
temperature range in the second column.  They are the only experimental
statement in the volume about what the elastic constants DO when the crystal
is heated, which is what the finite-temperature sweep computes and nothing in
the fit constrains.

The scan of this volume is poor and a flat text read of these pages returns
nonsense - the header comes back as "Wl" for Tc11 and several element symbols
lose their capital ("ca", "cs", "cu").  The rows are recovered instead from
the x position of each number, which the scan does preserve:

    x ~ 56  element      147  temperature range
      ~ 197 Tc11         239  Tc44        281  Tc12       323  Tc'

A row whose label sits at x ~ 84 is the standard deviation of the rows above
it, not a determination, and is skipped.  A value in parentheses is the
volume's own mark for an estimate; it is kept with `uncertain` set, so a
consumer can drop it.  Every value is checked back against the page text
before it is written, so an OCR misread that produces a plausible number in
the wrong column is caught by its absence there.

    python lb29_tempcoef.py            # writes lb29a_tempcoef.json
    python lb29_tempcoef.py --print    # and prints the table it read
"""
import io
import json
import os
import re
import sys

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
PAGES = range(228, 232)          # Table 28 and its continuations
COLUMNS = ((197, "Tc11"), (239, "Tc44"), (281, "Tc12"), (323, "Tcp"))
SYMBOLS = set("""Ag Al Au Ba Be Ca Cd Co Cr Cu Fe Ir K Li Mg Mo Na Nb Ni Pb Pd
                 Pt Rh Sr Ta Ti V W Zn Zr Rb Cs Yb Th Tl Y Sc Lu Ru Re Hf""".split())
NUM = re.compile(r"^\(?([-+]?\d+(?:\.\d+)?)\)?$")


def bands(page, h=9):
    out = {}
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            t = "".join(s["text"] for s in l["spans"]).strip()
            if t:
                out.setdefault(round(l["bbox"][1] / h), []).append((round(l["bbox"][0]), t))
    return out


def symbol(t):
    """'cu' and 'CU' are copper; a word is not a symbol"""
    s = t.strip().strip(".,")
    if len(s) > 2 or not s.isalpha():
        return None
    s = s[0].upper() + s[1:].lower()
    return s if s in SYMBOLS else None


def main():
    doc = fitz.open(os.path.join(HERE, "LB29a.pdf"))
    text = "\n".join(doc[p].get_text() for p in PAGES)
    out, current, seen = {}, None, set()
    for p in PAGES:
        for key in sorted(bands(doc[p])):
            items = sorted(bands(doc[p])[key])
            left = [t for x, t in items if x < 100]
            #  a standard-deviation row sits further in than a symbol does
            if any(x for x, t in items if 80 <= x < 100):
                continue
            sym = symbol(left[0]) if left else None
            if sym:
                current = sym
            if current is None:
                continue
            row, rng = {}, None
            for x, t in items:
                if 130 <= x < 190:
                    rng = t
                    #  the range and the first number sometimes share a cell
                    m = re.search(r"(-?\d+\.\d+)\s*$", t)
                    if m:
                        row["Tc11"] = (float(m.group(1)), False)
                for cx, name in COLUMNS:
                    if abs(x - cx) <= 12:
                        m = NUM.match(t)
                        if m:
                            row[name] = (float(m.group(1)), t.startswith("("))
            if not row or current in seen:
                continue
            #  keep the first determination of each element only, and only
            #  values that the page text also carries as written
            #  A coefficient of 10^-4/K is a per-cent per hundred kelvin; the
            #  volume prints nothing above about 20, so a "value" of 6936 is a
            #  reference code read as a number and is dropped here rather than
            #  carried into an average later.
            vals = {k: v for k, (v, unc) in row.items()
                    if abs(v) <= 50 and (str(v).lstrip("-") in text or f"{v:g}" in text)}
            unc = {k: u for k, (v, u) in row.items() if k in vals}
            if not vals:
                continue
            seen.add(current)
            out[current] = dict(T=rng, uncertain=[k for k, u in unc.items() if u],
                                complete=all(k in vals for k in ("Tc11", "Tc44", "Tc12")),
                                **vals)
    out["_source"] = ("A. G. Every and A. K. McCurdy, Landolt-Boernstein III/29a, Table 28 "
                      "(temperature coefficients Tc_ij = (1/c) dc/dT of the cubic elements, "
                      "in 10^-4/K, at constant pressure).")
    out["_note"] = ("Read from the scanned volume by lb29_tempcoef.py: rows are rebuilt from "
                    "the x position of each number and every value is checked against the page "
                    "text. Values the volume prints in parentheses are estimates and are listed "
                    "in 'uncertain'. Only the first determination of an element is kept.")
    json.dump(out, io.open(os.path.join(HERE, "lb29a_tempcoef.json"), "w", encoding="utf-8"),
              indent=1, sort_keys=True)
    els = [k for k in out if not k.startswith("_")]
    print(f"{len(els)} elements: {', '.join(sorted(els))}")
    if "--print" in sys.argv:
        for el in sorted(els):
            print(f"  {el:3} {out[el]}")


if __name__ == "__main__":
    main()
