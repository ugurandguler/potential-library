#!/usr/bin/env python3
"""
Read the elements out of Table 11 of Landolt-Boernstein III/29a - hexagonal.

Table 3 (cubic) is already read by lb29_extract.py.  Table 11 covers the
hexagonal materials and is where the candidate metals live: scandium, yttrium,
lutetium, hafnium, technetium, rhenium, ruthenium, thallium.  It also contains
the five hcp metals already in the library - beryllium, cobalt, magnesium,
titanium, zirconium - and those are the check: an extractor that reproduces
values typed in from other sources can be trusted on the ones it is new to.

Four things differ from the cubic table.

**The column order is 11, 33, 44, 12, 13.**  Not the cubic S11 S44 S12 C11 C44
C12, and not the order the constants are usually quoted in.  Taken from the
column headers, which sit at x = 251, 312, 366, 425, 484 on every page.

**Each material occupies two data rows**, compliance s then stiffness c.  There
is a mark column at x ~ 207 that says which is which, but the scan renders it as
"c", "C", "5", "S", "$" or omits it entirely - lutetium has no mark on either
row - so it cannot be relied on.  The rows are taken in order instead, the
second of a pair being the stiffness, and the pair is checked: s11 x C11 is of
order one for any material, so a mis-paired row is rejected rather than reported.

**Every determination may be followed by an s(n=N) spread row**, which is a data
row by any structural test and would shift the pairing by one: rhenium's four
rows are s, s(n=3), c, c(n=3), and pairing them in order gives (s, spread) and
(c, spread), both of which fail the product check, so rhenium vanished
altogether.  A spread row is recognised by the qualifier printed beside it in
the 100 < x < 200 column, which is text where a real sub-row - lutetium's
hydrogen content - is a number.

**Sub-labels must not reset the element.**  Material names start at the left
margin, x < 100; qualifiers sit at x ~ 128.  Treating any left-hand text as a
new material made "s(n=3)" end the block it belonged to.

Osmium is not in the table.  That is not an extraction failure and not something
to fill in from elsewhere: the volume has no complete set for it.

    python lb29_hex.py

Writes lb29a_hex.json.
"""
import json
import os
import re

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "LB29a.pdf")
OUT = os.path.join(HERE, "lb29a_hex.json")

PAGES = range(107, 131)         # Table 11, up to Table 12 (incomplete sets)
COLUMNS = ("C11", "C33", "C44", "C12", "C13")
XCOL = (251.0, 312.0, 366.0, 425.0, 483.0)
XTOL = 22.0
XNAME = 100.0                   # material names start left of this
XQUAL = (100.0, 200.0)          # qualifiers: s(n=3), at%H, alloy composition
DY = 3.5                        # a qualifier belongs to the row it lines up with
YCLUSTER = 5.0                  # numbers of one row share a y to about this

ELEMENTS = {
    "beryllium": "Be", "cobalt": "Co", "hafnium": "Hf", "lutetium": "Lu",
    "magnesium": "Mg", "rhenium": "Re", "ruthenium": "Ru", "scandium": "Sc",
    "technetium": "Tc", "thallium": "Tl", "titanium": "Ti", "yttrium": "Y",
    "zirconium": "Zr",
}
#  Labels that begin with an element name and are something else: alloys
#  ("Beryllium-copper alloys", "Magnesium-lithium alloys") and compounds
#  ("Beryllium oxide", "Titanium boride").  What separates them from the element
#  itself is a second real word after the name - the element rows carry only the
#  symbol and a footnote marker, and those are one or two letters.  Matching the
#  tail exactly instead lost technetium ("Technetium, Tc t&") and yttrium
#  ("Yttrium, Y =)"), whose markers the scan renders as stray characters.
WORD = re.compile(r"[A-Za-z]{3,}")


def items(page):
    """(y, x, text) for every line on the page, sorted top to bottom"""
    out = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            t = "".join(s["text"] for s in l["spans"]).strip()
            if t:
                out.append((l["bbox"][1], l["bbox"][0], t))
    return sorted(out)


def number(t):
    #  "- 43.0" is typeset with a space after the sign, and footnote markers
    #  ride along in the same cell: "110.0 d)", "107.2 4"
    t = re.sub(r"^([-+])\s+", r"\1", t.strip())
    t = re.split(r"\s", t)[0]
    return float(t) if re.fullmatch(r"[-+]?\d+\.?\d*", t) else None


def element_of(label):
    """the element a left-margin label names, or None"""
    low = label.lower()
    for name, sym in ELEMENTS.items():
        if low.startswith(name) and not WORD.search(label[len(name):]):
            return sym
    return None


def scan_page(page):
    """[(y, {column: value})] for the data rows, and the qualifier y positions"""
    L = items(page)
    rows, quals, names = {}, [], []
    for y, x, t in L:
        if x < XNAME:
            names.append((y, t))
            continue
        if XQUAL[0] <= x < XQUAL[1]:
            #  a number here is a real sub-row (lutetium's hydrogen content),
            #  text is a spread or composition qualifier
            if number(t) is None:
                quals.append(y)
            continue
        v = number(t)
        if v is None:
            continue
        for name, xc in zip(COLUMNS, XCOL):
            if abs(x - xc) <= XTOL:
                key = next((k for k in rows if abs(k - y) <= YCLUSTER), y)
                rows.setdefault(key, {})[name] = v
    data = [(y, d) for y, d in sorted(rows.items()) if len(d) >= 4]
    return names, quals, data


def main():
    doc = fitz.open(PDF)
    found = {}
    for p in PAGES:
        names, quals, data = scan_page(doc[p])
        blocks = {}
        for y, d in data:
            if any(abs(y - q) <= DY for q in quals):
                continue                      # s(n=N) spread row
            owner = [t for yy, t in names if yy <= y + DY]
            sym = element_of(owner[-1]) if owner else None
            if sym is None:
                continue
            #  keyed by the label, so a new material - element or not - starts a
            #  new block instead of appending to the last element seen
            blocks.setdefault((owner[-1], sym), []).append(d)
        for (label, sym), rs in blocks.items():
            for i in range(0, len(rs) - 1, 2):
                s, c = rs[i], rs[i+1]
                #  s11 [TPa^-1] x C11 [GPa] is of order one for any material
                prod = s.get("C11", 0) * c.get("C11", 0) * 1e-3
                if not (0.2 <= prod <= 5.0) or not all(k in c for k in COLUMNS):
                    continue
                rec = {k: c[k] for k in COLUMNS}
                rec["B_from_Cij"] = round(
                    (2*(c["C11"] + c["C12"]) + c["C33"] + 4*c["C13"]) / 9.0, 2)
                rec["page"] = p
                rec["label"] = label
                found.setdefault(sym, []).append(rec)

    json.dump(found, open(OUT, "w"), indent=1, sort_keys=True)
    print(f"{len(found)} element, "
          f"{sum(len(v) for v in found.values())} belirleme\n")
    print(f"{'el':4s}{'C11':>8s}{'C12':>8s}{'C13':>8s}{'C33':>8s}{'C44':>8s}"
          f"{'B(Cij)':>9s}  sayfa  etiket")
    print("-" * 76)
    for sym in sorted(found):
        for r in found[sym]:
            print(f"{sym:4s}{r['C11']:8.1f}{r['C12']:8.1f}{r['C13']:8.1f}"
                  f"{r['C33']:8.1f}{r['C44']:8.1f}{r['B_from_Cij']:9.1f}"
                  f"   {r['page']:4d}  {r['label'][:28]}")
    missing = sorted(set(ELEMENTS.values()) - set(found))
    if missing:
        print(f"\nnot found in Table 11: {', '.join(missing)}")
    print("\nwritten:", OUT)


if __name__ == "__main__":
    main()
