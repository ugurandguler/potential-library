#!/usr/bin/env python3
"""
Read the standard entropy and heat capacity of the elements from the CRC
Handbook, section 5, "Standard Thermodynamic Properties of Chemical Substances".

refdata.THERMO_298 was hand-entered from this table with a "verify before
publishing" warning attached, and ten elements added on 2026-08-03 had no row at
all - which is why their thermodynamics panel had no experimental rings and the
298 K table had blank reference columns.

The table is printed sideways and its column geometry is awkward, but the text
layer comes out in a clean order for the element rows:

    Sym \\n Name \\n 0.0 \\n S° \\n Cp°          then the gas-phase quartet

The leading 0.0 is the enthalpy of formation, which is zero by definition for an
element in its standard state, and that is what makes the row recognisable
without reading column positions.  Elements whose standard state is a gas do not
match, and neither does technetium, whose crystal row the volume leaves blank.

Values are J/(mol K) at 298.15 K.  The volume prints one decimal.

The twenty eight elements already in refdata are the check: this reads them too
and reports any disagreement, so a parser that drifts onto the wrong column
cannot pass silently.

    python crc_thermo.py

Writes crc_thermo.json.
"""
import json
import os
import re

import fitz

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(os.path.dirname(HERE), "CRChandbook.pdf")
OUT = os.path.join(HERE, "crc_thermo.json")
PAGES = range(849, 875)          # the substance table

NAMES = {
    "Ag": "Silver", "Al": "Aluminum", "Au": "Gold", "Ba": "Barium",
    "Be": "Beryllium", "Ca": "Calcium", "Cd": "Cadmium", "Co": "Cobalt",
    "Cr": "Chromium", "Cs": "Cesium", "Cu": "Copper", "Fe": "Iron",
    "Hf": "Hafnium", "Ir": "Iridium", "K": "Potassium", "Li": "Lithium",
    "Lu": "Lutetium", "Mg": "Magnesium", "Mo": "Molybdenum", "Na": "Sodium",
    "Nb": "Niobium", "Ni": "Nickel", "Pb": "Lead", "Pd": "Palladium",
    "Pt": "Platinum", "Rb": "Rubidium", "Re": "Rhenium", "Rh": "Rhodium",
    "Ru": "Ruthenium", "Sc": "Scandium", "Sr": "Strontium", "Ta": "Tantalum",
    "Ti": "Titanium", "Tl": "Thallium", "V": "Vanadium", "W": "Tungsten",
    "Y": "Yttrium", "Yb": "Ytterbium", "Zn": "Zinc", "Zr": "Zirconium",
}
NUM = r"-?\d+\.?\d*"


def scan(doc):
    """{symbol: (S, Cp)} for every element row the table gives"""
    out = {}
    for p in PAGES:
        try:
            text = doc[p].get_text()
        except Exception:                                       # noqa: BLE001
            continue
        for sym, name in NAMES.items():
            if sym in out:
                continue
            #  the symbol, the name, a formation enthalpy of exactly zero, then
            #  the two numbers wanted.  Anchored on the name so that "K" does
            #  not match inside another formula.
            m = re.search(rf"(?:^|\n){re.escape(sym)}\n{name}\n0\.0\n"
                          rf"({NUM})\n({NUM})\n", text)
            if m:
                out[sym] = (float(m.group(1)), float(m.group(2)), p)
    return out


def main():
    if not os.path.exists(PDF):
        raise SystemExit(f"not found: {PDF}")
    doc = fitz.open(PDF)
    got = scan(doc)

    print(f"{'el':4s}{'S CRC':>9s}{'Cp CRC':>9s}{'S ours':>9s}{'Cp ours':>9s}"
          f"{'sayfa':>7s}  ")
    print("-" * 50)
    bad, new = [], []
    for sym in sorted(got):
        s, cp, page = got[sym]
        ours = refdata.THERMO_298.get(sym)
        if ours:
            #  the volume prints one decimal, the table carries two, so a
            #  disagreement below 0.05 is the printing and not a difference
            off = max(abs(s - ours[0]), abs(cp - ours[1]))
            mark = "" if off <= 0.05 else f"  <-- {off:.2f} diff"
            if off > 0.05:
                bad.append(sym)
            print(f"{sym:4s}{s:9.1f}{cp:9.1f}{ours[0]:9.2f}{ours[1]:9.2f}"
                  f"{page:7d}{mark}")
        else:
            new.append(sym)
            print(f"{sym:4s}{s:9.1f}{cp:9.1f}{'-':>9s}{'-':>9s}{page:7d}  YENI")
    missing = sorted(set(NAMES) - set(got))
    json.dump({k: [v[0], v[1]] for k, v in got.items()},
              open(OUT, "w"), indent=1, sort_keys=True)
    print(f"\n{len(got)}/{len(NAMES)} element okundu")
    if missing:
        print(f"not found in the table: {', '.join(missing)}")
    if new:
        print(f"refdata'da olmayan (eklenecek): {', '.join(new)}")
    print("check: " + ("FAILED, differing: " + ", ".join(bad)
                           if bad else "mevcut satirlarin hepsi tutuyor"))
    print("written:", OUT)


if __name__ == "__main__":
    main()
