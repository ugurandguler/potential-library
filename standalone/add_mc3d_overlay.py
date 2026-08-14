#!/usr/bin/env python3
"""
Merge the Materials Cloud MC3D dispersions into library.json and score them.

Same idea as add_mp_overlay.py: our potential is evaluated at THEIR q-points,
so the residual is physics rather than interpolation.  What this adds is
coverage where Materials Project has none - Ag, Au, Cr, Mo, Nb, Ni, Pb, Pd, Rh,
Ta - and among those are the four bcc metals whose anisotropy the functional
form cannot reach, which until now had no independent phonon reference at all.

Stores under "mc3d" per element:
    q, f, marks, breaks   as fetched (frequencies in cm^-1)
    ours                  our frequencies at the same q
    stats                 rms difference and the two maxima

    python add_mc3d_overlay.py            # every element fetched
    python add_mc3d_overlay.py Nb Mo
"""
import json
import os
import sys

import numpy as np

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


def main(only):
    lib = json.load(open(os.path.join(HERE, "library.json")))
    src = json.load(open(os.path.join(HERE, "mc3d_phonon.json")))

    print(f"{'el':4s}{'branches':>9s}{'q':>6s}{'ours max':>10s}"
          f"{'MC3D max':>10s}{'rms':>8s}{'rel':>7s}")
    print("-" * 56)
    rel = []
    for el in sorted(src):
        if only and el not in only:
            continue
        v = lib.get(el)
        if not v:
            continue                      # Cd and Zn have no fit
        rec = src[el]
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(v)
        Phi = L.force_constants(cry, pot)
        ours = np.array([np.sort(L.frequencies(cry, pot, np.array(q), Phi))
                         * L.CM1 for q in rec["q"]])
        ref = np.array([sorted(r) for r in rec["f"]])
        nb = min(ours.shape[1], ref.shape[1])
        a, b = ours[:, :nb], ref[:, :nb]
        diff = a - b
        rms = float(np.sqrt((diff ** 2).mean()))
        scale = float(np.abs(b).max()) or 1.0
        out = dict(rec)
        out["ours"] = [[round(float(x), 1) for x in a[:, k]]
                       for k in range(nb)]
        out["f"] = [[round(float(x), 1) for x in b[:, k]] for k in range(nb)]
        out["stats"] = {"rms_cm1": round(rms, 2),
                        "rel_pct": round(100 * rms / scale, 1),
                        "ours_max": round(float(a.max()), 1),
                        "ref_max": round(float(b.max()), 1)}
        v["mc3d"] = out
        rel.append(100 * rms / scale)
        print(f"{el:4s}{nb:9d}{len(rec['q']):6d}{a.max():10.1f}{b.max():10.1f}"
              f"{rms:8.1f}{100 * rms / scale:6.1f}%")

    path = os.path.join(HERE, "library.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    if rel:
        print(f"\nmedian relative difference {np.median(rel):.1f} % over "
              f"{len(rel)} elements")
    print("merged into", path)


if __name__ == "__main__":
    main(set(sys.argv[1:]))
