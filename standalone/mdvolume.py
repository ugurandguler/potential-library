#!/usr/bin/env python3
"""The volume each potential CHOOSES at the measurement temperature.

The finite-temperature dispersion is measured in a box the barostat set, not
in the box the crystal has.  A potential that is 3 % too dense at 296 K will
give frequencies that are too high everywhere, and the finite-T column will
charge it for anharmonicity it never got wrong.  This separates the two.

It is also the one prediction in this whole comparison that the 0 K columns
cannot make at all: thermal expansion is not in a harmonic calculation.

Reads avol.txt - element, mean temperature, mean primitive-cell volume -
written on the cluster by avol.sh, and converts with the structure and axial
ratio known here, so nothing on the cluster has to guess either.

CAUTION on the reference.  refdata's a0 is a ~293 K measurement for most
elements but a 5 K one for the alkalis, so for those the difference below is
mostly REAL thermal expansion and not an error; they are marked.  Where the
dispersion paper states its own lattice constant, that is used instead and is
a true like-for-like comparison.

    python mdvolume.py
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

#  refdata carries these at 5 K, not at room temperature
COLD_A0 = {"Li", "Na", "K", "Rb", "Cs"}


def vprim(struct, a, coa):
    """volume of the primitive cell, the same cell mkprim.py builds"""
    if struct == "bcc":
        return a ** 3 / 2.0
    if struct == "fcc":
        return a ** 3 / 4.0
    return (math.sqrt(3.0) / 2.0) * a * a * (coa * a)


def main():
    lib = json.load(open(os.path.join(HERE, "library.json")))
    path = os.path.join(HERE, "avol.txt")
    if not os.path.exists(path):
        print("avol.txt yok - once kumede: bash ~/avol.sh > avol.txt")
        return 1
    print()
    print("%-5s%5s%7s%11s%11s%9s  %s"
          % ("el", "yapi", "T_MD", "V_MD", "V_ref", "fark", "referans"))
    print("-" * 62)
    for ln in open(path):
        p = ln.split()
        if len(p) < 3 or p[0] not in lib:
            continue
        el, T, V = p[0], float(p[1]), float(p[2])
        v = lib[el]
        ec = v.get("exp_curve", {})
        am, cm = ec.get("a_meas"), ec.get("c_meas")
        if am:
            a = am
            coa = (cm / am) if cm else v.get("c_over_a")
            src = "makalenin a_olcum"
        else:
            a = v["a0"]
            coa = v.get("c_over_a")
            src = "kutuphane a0 (5 K!)" if el in COLD_A0 else "kutuphane a0"
        vr = vprim(v["struct"], a, coa)
        print("%-5s%5s%6.0fK%11.3f%11.3f%8.2f%%  %s"
              % (el, v["struct"], T, V, vr, 100.0 * (V - vr) / vr, src))
    print()
    print("Volumes are per PRIMITIVE cell - one atom for the cubic metals, "
          "two for hcp.")
    print("A row marked (5 K!) compares a room-temperature MD box against a "
          "5 K crystal:")
    print("most of its difference is thermal expansion the potential got "
          "RIGHT, not error.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
