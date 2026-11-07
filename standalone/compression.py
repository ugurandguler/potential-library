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
    python compression.py --merge         # compression.json -> library.json
"""
import json
import os
import sys

import importlib.util

import numpy as np

import refdata

#  THE ANGULAR TREE'S latdyn, for every set (v1.4.7).  This file imported the
#  pair latdyn, which has no Legendre factor, so the UG sets (ug, tap_ug,
#  rc_ug) were scanned as if lam2 = lam4 = 0 with their own D, alpha, r0 - a
#  different potential: tungsten's switched UG record came out at -6.68
#  eV/atom at 0.9 a0 against -6.93 with the factor.  The angular latdyn with
#  lam2 = lam4 = 0 is the pair one exactly, so one engine serves every set.
#  Loaded by path so that refdata above stays this tree's.
_spec = importlib.util.spec_from_file_location(
    "latdyn_angular", os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "angular", "latdyn.py"))
L = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(L)
assert "lam2" in open(L.__file__, encoding="utf-8").read(), L.__file__

#  An arm record is one that carries the potential's five parameters.  Nothing
#  else in an element's record does - baseline_*, elasticT, dyn and mech carry
#  none of them - so this separates arms from payloads without a name list.
#  A NAME LIST IS THE BUG: every one of these files was written with
#  ("tap", "tap_ug"), the re-cut arms arrived and were silently skipped, and
#  `tap_force` arrived and was skipped again.  The next arm will not announce
#  itself either.
_PARAM = ("D", "r0", "alpha", "m", "gamma")


def arms(d):
    """the arm records inside one element's entry, by shape not by name"""
    return [k for k, v in d.items()
            if isinstance(v, dict) and all(p in v for p in _PARAM)]



HERE = os.path.dirname(os.path.abspath(__file__))
#  Display name -> library sub-key.  Only the ALIASES need to be here: "hard"
#  is the element's top-level record, which carries the same five parameters as
#  any arm and has no sub-key of its own, so shape detection cannot name it.
#  Every other arm is looked up by its own name - see SETS.get(key, key) below,
#  which is what lets an arm that did not exist when this file was written be
#  screened without editing this table.
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


#  0.50, not 0.55 (v1.4.7): an independent E(a) scan to 0.50 a0 found
#  basins between the two for Ag (switched MAU and UG), Nb (both) and Pb
#  (switched UG) that this scan, stopping at 0.55, reported as single wells
EV_A3_GPA = 160.21766208


def probe(el, rec, lo=0.50, step=0.005):
    """(x_onset, barrier, depth, P_enter) or four Nones if the well is single

    P_enter (v1.4.7) is the static pressure at which the crystal, held at 0 K,
    falls in: the highest pressure p = -dE/dV the compressed branch sustains
    between equilibrium and the barrier.  Above it the branch has no minimum
    left.  `reachable` asks about heat (barrier against k_B T_melt); this asks
    about load, and the two disagree - barium's re-cut record has a barrier
    five times k_B T_melt and is entered at 1.9 GPa.  Zero where there is no
    barrier at all.  The scan scales the cell uniformly (hcp at fixed c/a),
    which for one-atom cubic cells is exactly what a relaxation at fixed
    pressure does; an independent pressure scan agrees for Ba, Y and Lu and
    stopped short of converging for Na."""
    e = refdata.ELEMENTS[el]
    pot = L.Potential.from_record(rec)
    xs = np.arange(1.10, lo, -step)
    E, V = [], []
    for x in xs:
        cry = L.Crystal(e["struct"], e["a0"] * x, e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        V.append(cry.vol / len(cry.frac))
        try:
            E.append(L.energy(cry, pot))
        except Exception:
            E.append(np.nan)
    E, V = np.array(E), np.array(V)
    ok = np.isfinite(E)
    if not ok.any():
        return None, None, None, None
    #  the equilibrium the potential actually has, not the one it was given
    #  searched over a FIXED window, x >= 0.925: it was the first third of the
    #  scan, which moved when the scan was extended to 0.50 and shifted
    #  lithium's onset from 0.92 to 0.90 with nothing about lithium changed
    i0 = int(np.nanargmin(np.where(ok & (xs >= 0.925 - 1e-9), E, np.inf)))
    E0 = E[i0]
    #  walk inward from there and look for a second, deeper region
    below = [i for i in range(i0 + 1, len(xs))
             if ok[i] and E[i] < E0 - TOL]
    if not below:
        return None, None, None, None
    i = below[0]
    barrier = float(np.nanmax(E[i0:i + 1]) - E0)
    ib = i0 + int(np.nanargmax(E[i0:i + 1]))
    p = -np.gradient(E, V) * EV_A3_GPA
    p_enter = max(0.0, float(np.nanmax(p[i0:ib + 1]))) if barrier > TOL else 0.0
    return float(xs[i]), barrier, float(np.nanmin(E[i:]) - E0), p_enter


def main():
    args = sys.argv[1:]
    keys = None                      # decided from the library, see below
    if "--set" in args:
        i = args.index("--set")
        keys = [args[i + 1]]
        del args[i:i + 2]
    lib = json.load(open(os.path.join(HERE, "library.json")))
    if keys is None:
        #  EVERY arm in the library, not ("tap", "tap_ug").  This is a SCREEN;
        #  an arm it does not look at is an arm nobody screened, and the two
        #  named here were the two that existed when the file was written.
        #  Union over elements, in a stable order, so the table's columns do
        #  not depend on which element happens to be first.
        seen = []
        for d in lib.values():
            if not isinstance(d, dict):
                continue
            for k in arms(d):
                if k not in seen:
                    seen.append(k)
        keys = sorted(seen)
    scr = {}
    sp = os.path.join(os.path.dirname(HERE), "lammps", "md_screen_all.json")
    if os.path.exists(sp):
        scr = json.load(open(sp))

    def verdict(k):
        r = scr.get(k)
        if r is None:
            return "?"
        if r.get("lost"):
            return "disintegrated"
        if r.get("collapsed"):
            return "collapsed"
        if r.get("T", 300) > 400:
            return "suspect"
        return "ok"

    els = args or sorted(lib)
    out = {}
    print("Compression leak.  Criterion: barrier > k_B * T_melt (empirical, see the code).")
    print("Below that barrier the crystal does not survive to its melting point.\n")
    for key in keys:
        print(f"=== {key} ===")
        print(f"{'el':4s}{'x_onset':>12s}{'barr. eV':>10s}"
              f"{'depth':>11s}{'reachable':>10s}{'MD':>10s}")
        print("-" * 57)
        hit, bad = [], []
        for el in els:
            sub = SETS.get(key, key)
            rec = lib[el] if sub is None else lib[el].get(sub)
            if not rec or "m" not in rec:
                continue
            x, b, d, pe = probe(el, rec)
            v = verdict(f"{el}|{key}")
            if v in ("collapsed", "disintegrated"):
                bad.append(el)
            if x is None:
                out[f"{el}|{key}"] = None
                continue
            reach = b < KB * refdata.MELTING.get(el, 0.0)
            if reach:
                hit.append(el)
            out[f"{el}|{key}"] = {"x": x, "barrier": b, "depth": d,
                                  "reachable": bool(reach), "P_enter_GPa": pe}
            print(f"{el:4s}{x:12.3f}{b:10.3f}{d:11.1f}"
                  f"{('YES' if reach else 'no'):>10s}{v:>10s}")
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
          f"{len(out)} refreshed)")


def merge():
    """compression.json -> each record's `compression`, in the shape
    refresh_measured.py writes, and nothing else.  That file also rewrites
    stable, min_cm1, R and the MD verdict, and it looks for the dispersion-
    selected arm under the wrong key and does not know the nudge arm, so their
    compression never reached the library."""
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    comp = json.load(open(os.path.join(HERE, "compression.json")))
    n = 0
    for key, c in comp.items():
        el, arm = key.split("|", 1)
        sub = SETS.get(arm, arm)
        rec = lib.get(el) if sub is None else (lib.get(el) or {}).get(sub)
        if not isinstance(rec, dict):
            continue
        rec["compression"] = ({"basin": False} if c is None else
                              {"basin": True, "x": c["x"], "barrier": c["barrier"],
                               "depth": c["depth"], "reachable": c["reachable"],
                               "P_enter_GPa": c.get("P_enter_GPa")})
        n += 1
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"compression merged into {n} records")


if __name__ == "__main__":
    if "--merge" in sys.argv:
        merge()
    else:
        main()
