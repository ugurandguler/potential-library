#!/usr/bin/env python3
"""
Attach measured phonon frequencies to library.json and score them.

This is the only comparison on the page against a measurement rather than
another calculation, and the only one at finite q.  The fit sees elastic
constants alone, which are the q -> 0 limit, so every number here is out of
sample.

Stores per element, under "exp_phonon":
    points   [{name, q, exp[3], ours[3], err_pct}]
    mae      mean absolute error over all branches, per cent
    ref      where the measurements come from

    python add_exp_phonon.py
"""
import json
import os

import numpy as np

import latdyn as L
import refdata
import refdata_phonon as RP

HERE = os.path.dirname(os.path.abspath(__file__))
THZ = 33.35641                      # cm^-1 per THz

Q = {"fcc": {"X": (0.5, 0.0, 0.5), "L": (0.5, 0.5, 0.5)},
     "bcc": {"H": (0.5, -0.5, 0.5), "N": (0.0, 0.0, 0.5)},
     #  hcp has two atoms per cell, so six branches everywhere.  At A the six
     #  collapse to two levels of multiplicity four and two and the measured set
     #  is complete; at Gamma the three acoustic ones are zero by construction
     #  and the comparison is against the three optic modes.
     "hcp": {"A": (0.0, 0.0, 0.5), "G": (0.0, 0.0, 0.0)}}


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print(f"{'el':4s}{'pt':>3s}{'ours (THz)':>24s}{'measured':>24s}{'err %':>8s}")
    print("-" * 64)
    every = []
    for el in sorted(lib):
        v = lib[el]
        st = v["struct"]
        want = RP.for_element(el, st)
        if not want:
            v.pop("exp_phonon", None)
            continue
        cry = L.Crystal(st, v["a0"], v.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(v)
        Phi = L.force_constants(cry, pot)
        pts, errs = [], []
        for name, exp in want.items():
            q = np.array(Q[st][name])
            ours = [float(x) for x in
                    sorted(L.frequencies(cry, pot, q, Phi) * L.CM1 / THZ)]
            if name in RP.OPTIC_ONLY:
                #  drop the acoustic branches, which are zero at Gamma whatever
                #  the potential does and would otherwise dominate the average
                #  with a meaningless zero-against-zero comparison
                ours = ours[-len(exp):]
            e = 100 * np.mean(np.abs(np.array(ours) - np.array(exp))
                              / np.array(exp))
            errs.append(e)
            every.append(e)
            pts.append({"name": name, "q": list(map(float, q)),
                        "exp": [round(float(x), 2) for x in exp],
                        "ours": [round(float(x), 2) for x in ours],
                        "err_pct": round(float(e), 1)})
            print(f"{el:4s}{name:>3s}"
                  f"{str([round(x, 2) for x in ours]):>24s}"
                  f"{str([round(x, 2) for x in exp]):>24s}{e:8.1f}")
        v["exp_phonon"] = {"points": pts,
                           "mae": round(float(np.mean(errs)), 1),
                           "ref": RP.PHONON_EXP[el]["ref"]}

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)

    n = sum(1 for v in lib.values() if "exp_phonon" in v)
    print(f"\n{n} elements have measured frequencies; mean absolute error "
          f"{np.mean(every):.1f} % over {len(every)} points")
    print("merged into", path)

    #  which branch is off, and in which direction - the answer is the whole
    #  story of this potential, so it is worth printing rather than eyeballing
    soft = hard = 0
    for v in lib.values():
        for p in v.get("exp_phonon", {}).get("points", []):
            if p["ours"][0] < p["exp"][0]:
                soft += 1
            if p["ours"][-1] > p["exp"][-1]:
                hard += 1
    tot = len(every)
    print(f"lowest branch too soft in {soft}/{tot} points, "
          f"highest branch too stiff in {hard}/{tot}")


if __name__ == "__main__":
    main()
