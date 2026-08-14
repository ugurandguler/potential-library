#!/usr/bin/env python3
"""
Read Table 3 of Landolt-Boernstein III/29a - the cubic elements.

III/29a (Every and McCurdy, 1992) is the standard compilation of second-order
elastic constants: for each element it reports a weighted mean over the
published measurements with the spread, which is more defensible than any
single paper.  refdata.py was typed in from Kittel, Simmons & Wang and the CRC
tables, so this is the check the page footer has been asking for.

Two traps.  The column order is S11, S44, S12, C11, **C44**, C12 - C44 comes
before C12, which is the easiest thing to get backwards.  And the element name
sits in its own text line at the left margin while the numbers are separate
lines further right, so a row has to be rebuilt from everything sharing a y
band; reading the page as flat text silently pairs a name with the previous
row's numbers, which is how Cobalt first came out as C11 = 128.

    python lb29_extract.py
"""
import json
import os
import re

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
#  Table 3 runs to the footnotes on pdf page 18; page 19 starts
#  Table 4 (alloys), and letting the scan run on pulled a Co-elinvar
#  alloy row in as if it were elemental cobalt.
PAGES = range(11, 18)
SYMBOLS = set("""Ag Al Au Ba Be Ca Cd Co Cr Cu Fe Ir K Li Mg Mo Na Nb Ni Pb Pd
                 Pt Rh Sr Ta Ti V W Zn Zr""".split())
#  x centres of the six numeric columns S11 S44 S12 C11 C44 C12
COLUMNS = (208, 258, 308, 366, 417, 467)
NAMES = {"copper": "Cu", "gild": "Au", "gold": "Au", "cesimn": "Cs",
         "aluminum": "Al", "aluminium": "Al"}


def lines(page):
    out = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            t = "".join(s["text"] for s in l["spans"]).strip()
            if t:
                out.append((l["bbox"][1], l["bbox"][0], t))
    return sorted(out)


def symbol(label):
    #  the separator is a comma in most rows but a full stop in a few
    #  ("Strontium. Sr"), and some carry a phase prefix ("p -Co")
    m = re.match(r"([A-Za-z][A-Za-z\-]*)\s*[,.]\s*(?:[a-z]?\s*-\s*)?([A-Za-z]{1,2})\b",
                 label)
    if not m:
        return None
    name, sym = m.group(1).lower(), m.group(2)
    sym = sym[0].upper() + sym[1:].lower()
    if sym in SYMBOLS:
        return sym
    return NAMES.get(name)


def main():
    d = fitz.open(os.path.join(HERE, "LB29a.pdf"))
    out = {}
    for p in PAGES:
        L = lines(d[p])
        #  Rows are keyed by the y band.  A band with a label in the left
        #  margin starts a new element; a band with numbers and NOTHING in the
        #  margin is a second determination of the same element, which is how
        #  barium's two conflicting entries are stored.
        bands = {}
        for y, x, t in L:
            bands.setdefault(round(y / 9), []).append((x, t))
        current = None
        for key in sorted(bands):
            items = sorted(bands[key])
            label = [t for x, t in items if x < 80]
            #  the "s(n=...)" tag that marks a standard-deviation row is not
            #  always recovered from the PDF, so those rows are recognised by
            #  size instead: their compliances are a small fraction of the
            #  determination above them, while a genuine second determination
            #  is comparable (barium 157 against 123.7, cobalt 0.62 against
            #  8.81)
            if label:
                current = symbol(label[0])
                if not current:
                    continue
            elif current is None:
                continue
            y = min(yy for yy, xx, tt in L
                    if round(yy / 9) == key) if items else 0
            sym, t, x = current, "", 0
            #  Keep only cells that are a bare number, in x order.  Reference
            #  codes ("66Hl", "83V3,81L14") and qualifiers ("78K b,", "k,")
            #  contain letters and are dropped; without that filter the year in
            #  a citation is read as an elastic constant.
            cells = [(xx, tt) for xx, tt in items if xx >= 80]
            #  a few cells are typeset "- 43.0" with a space after the sign,
            #  and dropping those shifts every later column left by one - which
            #  is how Lead first came out as C11 = 14.8 instead of 48.8
            nums = []
            for _, tt in cells:
                tt = re.sub(r"^([-+])\s+", r"\1", tt.strip())
                if re.fullmatch(r"[-+]?\d+\.?\d*", tt):
                    nums.append(float(tt))
                elif nums:
                    break            # the reference codes start here
            if len(nums) < 6:
                continue
            #  Several elements carry more than one determination, and they can
            #  disagree badly: barium is listed twice, at C12 = -0.38 and
            #  C12 = +8.0.  Keep every row - the choice between them is made
            #  downstream, on whether (C11 + 2 C12)/3 matches the measured bulk
            #  modulus, which for barium rules the first one out at once
            #  (2.45 GPa against a measured 9.4-10.3).
            ref = next((tt for _, tt in cells
                        if re.search(r"\d\d[A-Z]", tt)), "")
            prev = out.get(sym, [])
            if prev and nums[0] < 0.35 * prev[-1]["_s11"]:
                #  This is the s(n=N) row: the spread over the measurements
                #  that went into the determination above it.  Attach it there
                #  rather than dropping it - it is the only uncertainty the
                #  volume gives, and the objective has no idea it exists.
                prev[-1]["unc"] = {"C11": nums[3], "C44": nums[4],
                                   "C12": nums[5]}
                continue
            row = {"_s11": nums[0],
                   "C11": nums[3], "C44": nums[4], "C12": nums[5],
                   "B_from_Cij": round((nums[3] + 2*nums[5]) / 3.0, 2),
                   "ref": ref.strip()}
            out.setdefault(sym, []).append(row)
    for v in out.values():
        for r in v:
            r.pop("_s11", None)
    json.dump(out, open(os.path.join(HERE, "lb29a_cubic.json"), "w"), indent=1)
    n = sum(len(v) for v in out.values())
    print(f"{len(out)} cubic elements, {n} determinations\n")
    print(f"{'el':4s}{'C11':>9s}{'C12':>9s}{'C44':>9s}{'B(Cij)':>9s}  ref")
    for s in sorted(out):
        for i, v in enumerate(out[s]):
            print(f"{s if i == 0 else '':4s}{v['C11']:9.2f}{v['C12']:9.2f}"
                  f"{v['C44']:9.2f}{v['B_from_Cij']:9.2f}  {v['ref'][:22]}")


if __name__ == "__main__":
    main()
