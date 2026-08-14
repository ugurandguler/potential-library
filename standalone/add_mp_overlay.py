#!/usr/bin/env python3
"""
For every element where Materials Project has a phonon band structure, evaluate
OUR dispersion at exactly MP's q-points and store it alongside.

Resampling one curve onto the other's grid would blur the comparison; sharing
the q-list makes it point for point, so the residual is entirely physics rather
than interpolation.  MP's labelled q are the Setyawan-Curtarolo ones, which is
the convention latdyn already uses, so the fractional coordinates transfer
directly with no basis change.

Adds, under library.json -> <el> -> mp.phonon:
    "ours"  same shape as "f": our frequencies in cm^-1 at MP's q
    "stats" branch-resolved and overall agreement

    python add_mp_overlay.py            # every element that has MP phonons
    python add_mp_overlay.py Cu Li
"""
import json, os, sys
import numpy as np
import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


def main(only):
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    els = [e for e in sorted(lib)
           if lib[e].get("mp", {}).get("phonon") and (not only or e in only)]
    if not els:
        raise SystemExit("no element has Materials Project phonon data")

    print(f"{'el':4s}{'branches':>9s}{'q':>6s}{'ours max':>10s}{'MP max':>8s}"
          f"{'RMS':>8s}{'rel':>7s}")
    print("-" * 52)
    for el in els:
        v = lib[el]
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(v["struct"], v["a0"], v.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(v)
        Phi = L.force_constants(cry, pot)
        q = np.array(v["mp"]["phonon"]["q"], dtype=float)
        ours = np.array([np.sort(L.frequencies(cry, pot, qq, Phi))*L.CM1
                         for qq in q])              # (nq, nbranch)
        mp = np.array(v["mp"]["phonon"]["f"], dtype=float)   # (nbranch, nq)
        mp = np.sort(mp, axis=0).T                            # (nq, nbranch)

        nb = min(ours.shape[1], mp.shape[1])
        a, b = ours[:, :nb], mp[:, :nb]
        diff = a - b
        rms = float(np.sqrt((diff**2).mean()))
        scale = float(np.abs(b).max()) or 1.0
        v["mp"]["phonon"]["ours"] = [[round(float(x), 1) for x in a[:, k]]
                                     for k in range(nb)]
        v["mp"]["phonon"]["stats"] = {
            "rms_cm1": round(rms, 2),
            "rel_pct": round(100*rms/scale, 1),
            "ours_max": round(float(a.max()), 1),
            "mp_max": round(float(b.max()), 1),
            "nbranch": nb, "nq": int(q.shape[0]),
        }
        print(f"{el:4s}{nb:9d}{q.shape[0]:6d}{a.max():10.1f}{b.max():8.1f}"
              f"{rms:8.1f}{100*rms/scale:6.1f}%")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\nmerged into {path}")


if __name__ == "__main__":
    main(set(sys.argv[1:]))
