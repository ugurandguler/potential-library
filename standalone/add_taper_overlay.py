#!/usr/bin/env python3
"""
Merge the tapered refit into library.json, beside the shipped parameters.

The library is fitted with the cutoff truncated hard.  That is self-consistent
at a fixed geometry - every elastic constant, phonon and thermodynamic quantity
on the page is sound - but phi2 does not vanish at rcut2, so a neighbour
crossing the sphere changes the energy in one step and the potential cannot be
used for anything that moves.  Switching the cutoff off smoothly fixes that,
and refitting with the switch in place turned out to do something else as well.

Eight cubic metals were reported out of reach: their measured C44/C' sits below
anything the form was found to produce.  With the switch on, all eight fit and
are dynamically stable, four of them exactly - niobium reaches the measured
R = 0.503 where the hard-cutoff floor was 2.190.  A control at matched budget
separates the two possible explanations: shortening the hard cutoff instead
(2.60 -> 2.405 -> 2.210 a0) moves niobium only 20.4 -> 18.3 -> 16.8 per cent,
while the switch reaches 0.00.  So it is the discontinuity, not the range.

This is not a free improvement and the page should not read like one.  The
median over all 38 goes the wrong way, 5.83 -> 7.67 per cent: the hcp metals
lose heavily (7.86 -> 19.81) and so do the alkalis.  Cadmium and zinc, which
fail on axial anisotropy rather than on C44/C', are not rescued either.  The
switch removes one specific limitation and costs accuracy elsewhere.

    python add_taper_overlay.py <dir with dense_*.json>
"""
import glob
import json
import os
import sys

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    src = sys.argv[1]
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))

    new = {}
    for f in glob.glob(os.path.join(src, "*.json")):
        for el, r in json.load(open(f)).items():
            if el != "element" and isinstance(r, dict) and "score" in r:
                new[el] = r
    if not new:
        raise SystemExit(f"no fit record inside {src}")

    print(f"{'el':4s}{'hard':>9s}{'tapered':>9s}{'diff':>9s}"
          f"{'R meas':>9s}{'R taper':>9s}  stable")
    print("-" * 60)
    n_stable = 0
    for el, r in sorted(new.items()):
        v = lib.get(el)
        if not v:
            continue
        e = refdata.ELEMENTS[el]
        c = r["Cij"]
        rec = {"rms": 100.0 * r["score"], "taper": r.get("taper"),
               "m": r["m"], "gamma": r["gamma"], "C": r["C"], "s3": r["s3"],
               "D": r["D"], "alpha": r["alpha"], "r0": r["r0"],
               "alpha3": r["alpha3"], "rcut2": r["rcut2"], "rcut3": r["rcut3"],
               "Cij": {k: c[k] for k in c if c.get(k) is not None}}
        #  the anisotropy is the whole point for the cubic metals, so it is
        #  stored rather than left to be recomputed in the viewer
        if e["struct"] in ("fcc", "bcc"):
            cp = 0.5 * (c["C11"] - c["C12"])
            ce = e["Cij"]
            rec["R"] = c["C44"] / cp if cp > 0 else None
            rec["R_exp"] = ce["C44"] / (0.5 * (ce["C11"] - ce["C12"]))
        pot = L.Potential.from_record(r)
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        mn = min(float(L.spectrum(cry, pot, nq=8).min()),
                 float(L.spectrum(cry, pot, nq=9).min()))
        rec["stable"] = bool(mn > -1e-3)
        rec["min_cm1"] = round(mn, 3)
        n_stable += rec["stable"]
        v["tap"] = rec
        print(f"{el:4s}{v['rms']:9.2f}{rec['rms']:9.2f}"
              f"{v['rms'] - rec['rms']:+9.2f}"
              f"{(rec.get('R_exp') or float('nan')):9.2f}"
              f"{(rec.get('R') or float('nan')):9.2f}"
              f"  {'evet' if rec['stable'] else 'HAYIR'}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\ntapered fit added to {len(new)} elements, "
          f"{n_stable} of them dynamically stable")
    print("birlestirildi:", path)


if __name__ == "__main__":
    main()
