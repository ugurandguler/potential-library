#!/usr/bin/env python3
"""
Read Table I of Brewer, "The Cohesive Energies of the Elements",
LBL-3720 Rev. (1977) - the source Kittel's table is drawn from.

Our Ecoh values already agree with it to 0.005 eV, so this is not a correction.
What it adds is the UNCERTAINTY on each one, which matters because the fit
imposes the cohesive energy as an exact constraint at every trial point.  The
spread runs from 0.3 % (Ag, Au) to 4.6 % (Ba): pinning a number that is only
known to a few per cent distorts everything the fit is free to move.

The volume lists a separate row per crystal structure - "(ccp)", "(bcc)",
"(hcp)", frequently OCR'd as "(ecp)", "(bee)", "(hep)" - so the row matching
the phase we fit is the one taken.  Values are kcal/gram-atom at 0 K.

    python brewer_extract.py
"""
import json
import os
import re

import fitz

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL_EV = 0.0433641
PAGES = range(6, 12)          # pdf pages 7-12
#  Rows the OCR mangles past recognition, read off rendered images instead.
#  Aluminium comes back as "A1" with a digit and vanadium in lower case; the
#  other four have their value on a neighbouring text line.  kcal/gram-atom.
BY_HAND = {
    "Al": (78.10, 1.0,  "unlabelled"),
    "Li": (37.71, 0.2,  "bcc"),
    "Mo": (157.2, 0.5,  "bcc"),
    "Pb": (46.78, 0.3,  "unlabelled"),
    "Rh": (132.5, 1.0,  "fcc"),
    "V":  (122.4, 2.0,  "bcc"),
}

STRUCT = {"ccp": "fcc", "ecp": "fcc", "eep": "fcc", "cep": "fcc",
          "bcc": "bcc", "bee": "bcc",
          "hcp": "hcp", "hep": "hcp", "dhep": None, "dhcp": None}


def rows(page):
    L = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            t = "".join(s["text"] for s in l["spans"]).strip()
            if t:
                L.append((round(l["bbox"][1], 1), round(l["bbox"][0]), t))
    bands = {}
    for y, x, t in sorted(L):
        bands.setdefault(round(y / 8), []).append((x, t))
    return [sorted(v) for _, v in sorted(bands.items())]


def parse(cells):
    """(symbol, structure or None, value_kcal, uncertainty_kcal)"""
    if not cells or cells[0][0] > 90:
        return None
    label = cells[0][1]
    m = re.match(r"([A-Z][a-z]?)\s*(?:'|\.)?\s*(?:\(([a-z]+)\))?", label)
    if not m:
        return None
    sym, st = m.group(1), STRUCT.get((m.group(2) or "").lower(), False)
    if sym not in refdata.ELEMENTS:
        return None
    for x, t in cells[1:]:
        if x > 200:
            break
        #  Silver reads "68. Ot o. 2" and means 68.0 +- 0.2.  Order matters:
        #  join the split decimals FIRST, then turn a t/T that sits between
        #  digits into the plus-minus.  Matching "t 0." instead eats the
        #  leading zero of the uncertainty and reports 2.0 where the volume
        #  says 0.2 - a factor of ten, and in the direction that would make
        #  the constraint look far less certain than it is.
        s = t.replace("O", "0").replace("o", "0").replace("!", "±")
        s = re.sub(r"(\d)\.\s+(\d)", r"\1.\2", s)
        s = re.sub(r"(?<=\d)\s*[TtI]\s*(?=[\d.])", "±", s)
        mm = re.match(r"\(?(\d+\.?\d*)\)?\s*±?\s*(\d+\.?\d*)?", s)
        if mm and float(mm.group(1)) > 5:
            return sym, st, float(mm.group(1)), \
                   float(mm.group(2)) if mm.group(2) else None
    return None


def main():
    d = fitz.open(os.path.join(HERE, "brewer.pdf"))
    got = {}
    for p in PAGES:
        for cells in rows(d[p]):
            r = parse(cells)
            if not r:
                continue
            sym, st, v, u = r
            want = refdata.ELEMENTS[sym]["struct"]
            #  prefer the row for the phase we fit; fall back to an unlabelled
            #  row, which is what the volume gives for elements with one phase
            if st == want or (st is False and sym not in got):
                got[sym] = {"kcal": v, "unc_kcal": u,
                            "eV": round(v * KCAL_EV, 4),
                            "unc_eV": round(u * KCAL_EV, 4) if u else None,
                            "struct_row": st or "unlabelled"}
    for sym, (v, u, st) in BY_HAND.items():
        got[sym] = {"kcal": v, "unc_kcal": u, "eV": round(v * KCAL_EV, 4),
                    "unc_eV": round(u * KCAL_EV, 4), "struct_row": st,
                    "read": "by hand from the page image"}
    json.dump(got, open(os.path.join(HERE, "brewer_ecoh.json"), "w"), indent=1)
    print(f"{len(got)}/{len(refdata.ELEMENTS)} elements read\n")
    print(f"{'el':4s}{'Brewer eV':>11s}{'+-':>8s}{'unc %':>7s}{'ours':>8s}"
          f"{'diff':>8s}  row")
    worst = []
    for el in sorted(got):
        g = got[el]
        ours = refdata.ELEMENTS[el]["Ecoh"]
        pct = 100 * g["unc_eV"] / g["eV"] if g["unc_eV"] else float("nan")
        worst.append((pct, el))
        print(f"{el:4s}{g['eV']:11.3f}"
              f"{(g['unc_eV'] if g['unc_eV'] else 0):8.3f}{pct:7.1f}"
              f"{ours:8.2f}{g['eV']-ours:+8.3f}  {g['struct_row']}")
    ok = [w for w in worst if w[0] == w[0]]
    ok.sort(reverse=True)
    print("\nleast certain cohesive energies (these are imposed exactly):")
    for pct, el in ok[:6]:
        print(f"   {el:3s} +-{pct:.1f} %")


if __name__ == "__main__":
    main()
