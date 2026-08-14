#!/usr/bin/env python3
"""
Is the fitted crystal the bottom of the well, or a ledge above a deeper one?

Everything the fit constrains lives at one volume.  The cohesive energy, the
lattice constant, the bulk modulus and the elastic constants are all evaluated
at the experimental lattice constant, and the elastic constants are curvature
*at* that point.  Nothing in the objective asks what the energy does when the
crystal is squeezed, and nothing in the two constraints added after the vacancy
work asks it either - those look at the first shell and at the ratio of the
two- and three-body sums, both at equilibrium.

That gap has a specific consequence for this functional form.  The cutoffs are
fixed lengths in angstroms, so compressing the cell pulls further shells inside
them and the number of triplets climbs fast.  phi2 has a repulsive core and
resists; phi3 depends only on r1 + r2 and, where C is negative, does not.  If
the three-body sum wins, the energy turns over and falls, and the fitted
structure is a local minimum sitting beside a bottomless basin.

Molecular dynamics finds that basin.  It is what a "collapse" verdict in the
screen is made of, and it is invisible to every static test the library runs -
the elastic constants, the phonons and the vacancy energy are all properties of
the fitted geometry and all look healthy while the crystal sits on the ledge.

The measure that matters is not how deep the basin is but how easy it is to
reach.  Barium's basin is 208 eV/atom deep and barium is fine in MD, because
reaching it costs 0.17 eV/atom and 600 K supplies about 0.078.  Lithium's basin
costs 0.001 eV and lithium disintegrates.  So this reports, per element:

    x_onset   the least compression at which the crystal is already downhill
    barrier   the highest point between equilibrium and that basin

and calls a set reachable when the barrier falls below k_B T_melt - the same
threshold the fit applies, so the diagnostic and the constraint cannot drift
apart.  That threshold is empirical; the derivation it replaced is documented
beside it below.  The equilibrium reference is the local minimum found near x = 1,
not the energy exactly at the experimental lattice constant, because the fit
does not drive the pressure to zero and a residual few meV would otherwise read
as a basin at x = 0.99 in half the library.

    python compression.py                 # every element, tapered sets
    python compression.py Na Li --set tap
"""
import json
import os
import sys

import numpy as np

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = {"hard": None, "tap": "tap", "ug": "ug", "tap_ug": "tap_ug"}
#  The threshold the fit itself applies, so the diagnostic and the constraint
#  can never drift apart - two places with two thresholds is how that starts.
#
#  It is EMPIRICAL.  The obvious derivation, (3/2) k_B T per atom at the screen
#  temperature, is wrong and was refuted by the test built to check it: barium's
#  switched angular set has a barrier of 0.102 eV against 0.116 and survived
#  900 K, then 1000 and 1100, and collapsed only at 1200 - uniform compression
#  of the whole cell is one collective coordinate and thermal fluctuations do
#  not drive it coherently.  Calibrating against that measured collapse moved
#  the threshold by 1.46, which cancels the 3/2 and leaves k_B T_melt.  Checked
#  against the MD screen on every set that has a basin at all: 30 of 30.

#  a basin has to be deeper than the noise in the energy to count as one
TOL = 5e-3
KB = 8.617333262e-5


def probe(el, rec, lo=0.55, step=0.005):
    """(x_onset, barrier, depth) or (None, None, None) if the well is single"""
    e = refdata.ELEMENTS[el]
    pot = L.Potential.from_record(rec)
    xs = np.arange(1.10, lo, -step)
    E = []
    for x in xs:
        cry = L.Crystal(e["struct"], e["a0"] * x, e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        try:
            E.append(L.energy(cry, pot))
        except Exception:
            E.append(np.nan)
    E = np.array(E)
    ok = np.isfinite(E)
    if not ok.any():
        return None, None, None
    #  the equilibrium the potential actually has, not the one it was given
    i0 = int(np.nanargmin(np.where(ok, E, np.inf)[:len(xs) // 3]))
    E0 = E[i0]
    #  walk inward from there and look for a second, deeper region
    below = [i for i in range(i0 + 1, len(xs))
             if ok[i] and E[i] < E0 - TOL]
    if not below:
        return None, None, None
    i = below[0]
    return float(xs[i]), float(np.nanmax(E[i0:i + 1]) - E0), \
        float(np.nanmin(E[i:]) - E0)


def main():
    args = sys.argv[1:]
    keys = ["tap", "tap_ug"]
    if "--set" in args:
        i = args.index("--set")
        keys = [args[i + 1]]
        del args[i:i + 2]
    lib = json.load(open(os.path.join(HERE, "library.json")))
    scr = {}
    sp = os.path.join(os.path.dirname(HERE), "lammps", "md_screen_all.json")
    if os.path.exists(sp):
        scr = json.load(open(sp))

    def verdict(k):
        r = scr.get(k)
        if r is None:
            return "?"
        if r.get("lost"):
            return "DAGILDI"
        if r.get("collapsed"):
            return "COKTU"
        if r.get("T", 300) > 400:
            return "SUPHELI"
        return "saglam"

    els = args or sorted(lib)
    out = {}
    print("Sikisma kacagi.  Olcut: engel > k_B * T_erime (ampirik, bkz. kod).")
    print("Engel bunun altindaysa kristal erime noktasina dayanmaz.\n")
    for key in keys:
        print(f"=== {key} ===")
        print(f"{'el':4s}{'x_baslangic':>12s}{'engel eV':>10s}"
              f"{'derinlik':>11s}{'ulasilir':>10s}{'MD':>10s}")
        print("-" * 57)
        hit, bad = [], []
        for el in els:
            rec = lib[el] if SETS[key] is None else lib[el].get(SETS[key])
            if not rec or "m" not in rec:
                continue
            x, b, d = probe(el, rec)
            v = verdict(f"{el}|{key}")
            if v in ("COKTU", "DAGILDI"):
                bad.append(el)
            if x is None:
                out[f"{el}|{key}"] = None
                continue
            reach = b < KB * refdata.MELTING.get(el, 0.0)
            if reach:
                hit.append(el)
            out[f"{el}|{key}"] = {"x": x, "barrier": b, "depth": d,
                                  "reachable": bool(reach)}
            print(f"{el:4s}{x:12.3f}{b:10.3f}{d:11.1f}"
                  f"{('EVET' if reach else 'hayir'):>10s}{v:>10s}")
        print(f"\n  reachable at 600 K: {' '.join(hit) or 'none'}")
        print(f"  collapses/disperses in MD: {' '.join(bad) or 'none'}\n")
    #  MERGE.  Running this on a handful of elements to check something used to
    #  replace the whole file, and the loss is invisible: the JSON stays valid
    #  and every element it no longer mentions simply stops appearing on the
    #  page.  Same failure the screen and the library both had.
    pa = os.path.join(HERE, "compression.json")
    all_ = {}
    if os.path.exists(pa):
        try:
            all_ = json.load(open(pa))
        except Exception:
            all_ = {}
    kept = len(all_)
    all_.update(out)
    json.dump(all_, open(pa, "w"), indent=1, sort_keys=True)
    print(f"-> compression.json ({len(all_)} records; {kept} already present, "
          f"{len(out)} tazelendi)")


if __name__ == "__main__":
    main()
