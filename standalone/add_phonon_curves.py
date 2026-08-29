#!/usr/bin/env python3
"""
Attach measured dispersion CURVES to library.json, mapped onto the drawn path.

The segment mapping is done here rather than in the viewer.  Which line of the
Brillouin zone a branch labelled [0zz] lies on, and which way round the drawn
path traverses it, is a fact about crystallography and about
build_library.SC_PATH - not about canvas drawing - and putting it in the page
would mean the answer lived in JavaScript where it cannot be tested.

Stores per element, under "exp_curve":
    T_K    the measurement temperature, which is part of the datum
    ref    the paper, with its table
    segs   {"G|X": [[t, nu, sigma, branch], ...], ...}, t the fraction along
           that segment in the direction a -> b

    python add_phonon_curves.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

import refdata_phonon_curves as C          # noqa: E402

#  The fcc path as build_library draws it.  A segment appears here in the
#  sense the path walks it; Sigma is measured G -> K and drawn K -> G, so the
#  fraction is flipped rather than the data being re-sorted.
#  build_library.SC_PATH, in the sense each path walks its segments
DRAWN = {
    "fcc": [("G", "X"), ("X", "W"), ("W", "K"), ("K", "G"),
            ("G", "L"), ("L", "U"), ("U", "W"), ("W", "L"), ("L", "K"),
            ("U", "X")],
    "bcc": [("G", "H"), ("H", "N"), ("N", "G"), ("G", "P"), ("P", "H"),
            ("P", "N")],
    "hcp": [("G", "M"), ("M", "K"), ("K", "G"), ("G", "A"), ("A", "L"),
            ("L", "H"), ("H", "A"), ("L", "M"), ("K", "H")],
}


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print(f"{'el':4s}{'T':>6s}{'seg':>6s}{'points':>8s}   segments")
    print("-" * 58)
    n_el = 0
    for el in sorted(C.PHONON_CURVE):
        if el not in lib or not isinstance(lib[el], dict):
            print(f"{el:4s}  not in library.json, skipped")
            continue
        rec = C.PHONON_CURVE[el]
        struct = rec.get("struct", "fcc")
        segs, total = {}, 0
        #  A paper may report frequencies only AT the symmetry points, with
        #  no branch in between.  Those cannot be placed by direction - L and
        #  H lie on none of [00z], [z00] or [zz0] - so they are placed by
        #  which segment ends where.  Each measured frequency is listed once,
        #  under the paper's own irreducible-representation label, so nothing
        #  here has to decide how many modes a label stands for.
        at = rec.get("points", {})
        for a, b in DRAWN[struct]:
            pts = C.points_on(el, a, b, struct)
            flip = False
            if not pts:
                pts = C.points_on(el, b, a, struct)   # measured the other way
                flip = bool(pts)
            if not pts:
                continue
            rows = [[round(1.0 - t if flip else t, 6), nu, sig, lab]
                    for t, nu, sig, lab in pts]
            segs[f"{a}|{b}"] = rows
            total += len(rows)
        for a, b in DRAWN[struct]:
            rows = ([[0.0, nu, sig, lab] for nu, sig, lab in at.get(a, ())]
                    + [[1.0, nu, sig, lab] for nu, sig, lab in at.get(b, ())])
            if rows:
                segs.setdefault(f"{a}|{b}", []).extend(rows)
                total += len(rows)
        if not segs:
            lib[el].pop("exp_curve", None)
            continue
        lib[el]["exp_curve"] = {"T_K": rec["T_K"], "ref": rec["ref"],
                                "segs": segs,
                                #  the lattice constant of the crystal that was
                                #  actually in the beam, where the paper states
                                #  it.  refdata's a0 is a ~293 K quantity (5 K
                                #  for the alkalis) and these measurements run
                                #  from 9 to 296 K, so for a few elements the
                                #  comparison was being made at the wrong
                                #  volume - see curve_mae.
                                **({"a_meas": rec["a_meas"]}
                                   if "a_meas" in rec else {}),
                                **({"a_meas_from": rec["a_meas_from"]}
                                   if "a_meas_from" in rec else {}),
                                **({"c_meas": rec["c_meas"]}
                                   if "c_meas" in rec else {}),
                                #  true for a record published as a list of
                                #  points rather than as branches - whether
                                #  those points are the symmetry points
                                #  (beryllium) or a scattered handful that the
                                #  paper tabulated beside a figure (tungsten)
                                "points_only": bool(rec.get("sparse")
                                                    or (at and not rec.get(
                                                        "branches")))}
        n_el += 1
        print(f"{el:4s}{rec['T_K']:5d}K{len(segs):6d}{total:8d}   "
              + ", ".join(sorted(segs)))

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\n{n_el} elements carry a measured dispersion curve")
    print(f"merged into {path}")


if __name__ == "__main__":
    main()
