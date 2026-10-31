#!/usr/bin/env python3
"""
Fold the published-potential phonon comparison into the library: for each
element with a measured dispersion, the library's own error beside that of
every published EAM and MEAM potential scored the same way.

curve_mae_published.json holds the published potentials (92 files on 29
metals), curve_mae_aitken.json the three MEAM files of Aitken et al. (2021),
which are not in the NIST repository and are kept apart from the 92 so that
the set every other comparison uses stays as it was.  Both were made by
curve_mae_published.py from force constants that lammps/phonon_published.py
took through LAMMPS by finite differences; the same reading - branch labels,
the polarisation test, the measurement volume, symmetry points counted once
- scores the library's own records, so only the force constants differ.

The library's errors are read from curve_mae.py's own table, hard-cut and
switched MAU, so the two columns cannot drift apart from what the page
prints elsewhere.

Writes library.json[el]["published_phonons"]:
    ours_hard, ours_switched      per cent of the highest measured frequency
    published                     [[label, family, pct], ...] best first
    median, n                     over `published`
    aitken                        [[label, pct], ...] (Aitken 2021; Ni Al Cu)
"""
import json
import os
import re
import statistics as st
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def ours():
    """the library's hard-cut and switched errors, from curve_mae.py's table;
    the switched one is None where the arm is absent, as in a fresh clone,
    which cannot rebuild it (refresh.py's header says why)"""
    out = subprocess.run([sys.executable, "curve_mae.py"], cwd=HERE, capture_output=True,
                         text=True, encoding="utf-8",
                         env=dict(os.environ, PYTHONIOENCODING="utf-8")).stdout
    block = out.split("Against the FITTED MODELS")[0]
    row = re.compile(r"^([A-Z][a-z]?)\*?\s+\d+K\s+(\d+)\s+([\d.]+)\s+([\d.]+)%"
                     r"(?:\s+([\d.]+)\s+([\d.]+)%|\s+-)", re.M)
    return {m[1]: (float(m[4]), float(m[6]) if m[6] else None) for m in row.finditer(block)}


def label(key):
    return key.split("|", 1)[1]


def main():
    mine = ours()
    pub = json.load(open(os.path.join(HERE, "curve_mae_published.json")))["scores"]
    ref_path = os.path.join(HERE, "curve_mae_aitken.json")
    ref = json.load(open(ref_path))["scores"] if os.path.exists(ref_path) else {}
    lib = json.load(open(os.path.join(HERE, "library.json")))
    n_el = 0
    for el, d in lib.items():
        d.pop("published_phonons", None)
        if el not in mine:
            continue
        rows = sorted(([label(k), v["family"], round(v["pct"], 2)]
                       for k, v in pub.items() if v["el"] == el), key=lambda r: r[2])
        rr = [[label(k), round(v["pct"], 2)] for k, v in ref.items() if v["el"] == el]
        if not rows and not rr:
            continue
        d["published_phonons"] = dict(
            ours_hard=mine[el][0], ours_switched=mine[el][1], published=rows,
            n=len(rows), median=round(st.median(r[2] for r in rows), 2) if rows else None,
            aitken=rr)
        n_el += 1
    tmp = os.path.join(HERE, "library.json.tmp")
    json.dump(lib, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(HERE, "library.json"))
    print(f"published-potential phonons for {n_el} elements -> library.json")


if __name__ == "__main__":
    main()
