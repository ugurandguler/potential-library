#!/usr/bin/env python3
"""The three ways of comparing this library to a measured dispersion, side by
side, so the question the finite-temperature runs were made to answer can be
read off one table.

    0 K @ a0        the harmonic dispersion at the library's lattice constant
    0 K @ a_meas    the same, at the lattice constant of the crystal the
                    paper actually measured, where the paper states it
    finite T        the dispersion LAMMPS measures by displacement correlation
                    during equilibrium MD at the paper's temperature, in a box
                    the barostat set itself

The first two differ by a volume; the third differs by anharmonicity as well,
and by whatever thermal expansion this potential predicts rather than the one
the crystal has.  A potential can therefore lose on the third column while
being more right about the physics, and can win on it for the wrong reason -
which is why the fourth column is here.  It is the disagreement between two
q-points that the crystal's own symmetry forces to be equal, measured in the
same run, and a difference between columns smaller than it has not been
measured at all.

All three columns identify a branch the same way - L by its eigenvector lying
along the wavevector, not by sorted position - which an earlier pass did not
do: it scored the finite-T column by nearest branch, because phana's option 4
returns no eigenvectors.  Option 9 does, so the concession was withdrawn.
That correction moved every finite-T number the wrong way, which is what a
correction to a flattering rule looks like.

The finite-T numbers come from mdres.json, written on the cluster by
harvest.sh, because the eigenvector files themselves do not survive the link.
A run still in progress is scored on its latest partial binary and marked;
those numbers move as the run continues and are not to be quoted.

All three columns are the MAU arm - the tapered record, shipped as
`<El>_taper.ugur`.  The other three arms (UG, and the two re-cut candidates)
were not run at finite temperature; they have screens and potential files, so
that was a scope decision and not a result about them.

    python three_columns.py [El ...]
"""
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import curve_mae as CM                      # noqa: E402

FULL = 6500000          # the step count a completed run reaches


def zero_k(els, use_a_meas):
    """{el: percent} for the tapered arm, at one choice of volume"""
    CM.USE_A_MEAS = use_a_meas
    buf, keep = io.StringIO(), sys.stdout
    sys.stdout = buf
    try:
        CM.main(tuple(els))
    finally:
        sys.stdout = keep
    out = {}
    for ln in buf.getvalue().splitlines():
        p = ln.split()
        if len(p) >= 7 and p[1].endswith("K") and p[-1].endswith("%"):
            try:
                out[p[0].rstrip("*")] = float(p[-1].rstrip("%"))
            except ValueError:
                continue
    return out


def main():
    els = [a for a in sys.argv[1:] if not a.startswith("-")]
    lib = json.load(open(os.path.join(HERE, "library.json")))
    try:
        md = json.load(open(os.path.join(HERE, "mdres.json")))
    except OSError:
        md = {}
    if not els:
        els = sorted(el for el, v in lib.items()
                     if isinstance(v, dict) and "exp_curve" in v
                     and el in md)
    a0 = zero_k(els, False)
    am = zero_k(els, True)

    print()
    print("%-6s%5s%6s%11s%12s%14s%11s%16s"
          % ("el", "T", "points", "box", "0 K @ a0", "0 K @ a_meas",
             "finite T", "symmetry floor"))
    print("-" * 85)
    part = False
    for el in els:
        v = lib[el]
        m = md.get(el)
        if m and "error" in m:
            m = None
        flag = "*" if v["struct"] == "hcp" else ""
        if m and m.get("steps", 0) < FULL:
            flag += "!"
            part = True
        row = "%-6s%4dK%6s%11s" % (el + flag, v["exp_curve"]["T_K"],
                                   m["n"] if m else "-",
                                   m.get("mesh", "?") if m else "-")
        for tab in (a0, am):
            row += ("%11.1f%%" % tab[el]) if el in tab else "%12s" % "-"
        row += ("%10.1f%%" % m["pct"]) if m else "%11s" % "-"
        row += ("%15.1f%%" % m["sc_pct"]) if m else "%16s" % "-"
        print(row)
    print()
    print("MAU arm throughout - the tapered record, <El>_taper.ugur.  The "
          "UG and re-cut arms")
    print("were not run at finite temperature.")
    print()
    print("Per cent of the highest measured frequency of that element.  The "
          "last column is")
    print("the spread between q and a symmetry image of q in the same MD "
          "measurement: it")
    print("is what the run cannot resolve, and a difference between columns "
          "smaller than")
    print("it has not been measured.")
    if any(lib[el]["struct"] == "hcp" for el in els):
        print()
        print("* hcp: six branches, split by polarisation along q and "
              "perpendicular or parallel")
        print("  to the BASAL PLANE, then acoustic or optic by frequency "
              "order WITHIN a class.")
        print("  Sigma_3 and Sigma_4 narrow to a transverse pair, because "
              "the record refuses")
        print("  to guess which is polarised along c; nu1..nu6 carry no "
              "polarisation and keep")
        print("  the six-way nearest read.  Both columns are read the same "
              "way, and the")
        print("  Cartesian frame matters: LAMMPS writes the box in "
              "restricted triclinic")
        print("  form, which rotates the axes for fcc and bcc, and a "
              "polarisation built in")
        print("  the analytic code's frame picked a different branch from "
              "the analytic one")
        print("  at 70% of copper's points until that was corrected.")
    if part:
        print()
        print("! run not finished; scored on a partial binary and will move.")


if __name__ == "__main__":
    main()
