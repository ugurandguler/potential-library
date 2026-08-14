#!/usr/bin/env python3
"""
Pull measured phonon frequencies out of Landolt-Boernstein III/13a.

The volume tabulates, per element, neutron-scattering frequencies as (zeta, nu)
down a set of columns, one column per branch and direction.  A plain text
extraction scrambles those columns together; the layout only survives if the
words are re-grouped by their x coordinate, which is what this does.

Two further quirks the parser has to handle: several branches are stacked
inside one column, so a branch ends wherever zeta stops increasing; and the
uncertainty is sometimes glued to the value as "0.80(3)".

Nothing here decides which branch is which - it reports every branch with its
endpoint, and the assignment to H, N, P or X, L is done by hand against the
figures, because an automatic guess is exactly the kind of error that would be
invisible afterwards.

    python lb_extract.py 62          # page 62 (Fe)
"""
import re
import sys

import fitz

BOUNDS = (140, 240, 335, 405)


def branches(page, ytop=290, ybot=640, bounds=BOUNDS):
    words = [w for w in page.get_text("words") if ytop < w[1] < ybot]
    cols = {}
    for w in words:
        k = sum(w[0] > b for b in bounds)
        cols.setdefault(k, []).append((round(w[1] / 4) * 4, w[0], w[4]))
    out = []
    for k in sorted(cols):
        rows = {}
        for y, x, t in sorted(cols[k]):
            rows.setdefault(y, []).append((x, t))
        pairs = []
        for y in sorted(rows):
            s = " ".join(t for _, t in sorted(rows[y]))
            s = re.sub(r"(\d)\((\d+)\)", r"\1 (\2)", s)      # unglue 0.80(3)
            m = re.match(r"(\d*\.?\d+)\s+(\d*\.?\d+)", s)
            if m:
                pairs.append((float(m.group(1)), float(m.group(2))))
        #  a column holds several branches stacked; a new one starts wherever
        #  zeta stops increasing
        cur = []
        for z, v in pairs:
            if cur and z <= cur[-1][0]:
                out.append((k, cur))
                cur = []
            cur.append((z, v))
        if cur:
            out.append((k, cur))
    return out


def main(page_no):
    d = fitz.open("LB13a.pdf")
    for k, br in branches(d[page_no - 1]):
        if len(br) < 3:
            continue
        z0, v0 = br[0]
        z1, v1 = br[-1]
        print(f"  col{k}  {len(br):3d} pts   zeta {z0:5.3f}->{z1:5.3f}   "
              f"nu {v0:6.2f}->{v1:6.2f} THz")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 62)
