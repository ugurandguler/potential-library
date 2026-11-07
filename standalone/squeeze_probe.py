#!/usr/bin/env python3
"""
The two-leg squeeze: does a record have an energy hole inside its repulsive
core?  An isolated linear trimer A-X-B with both legs shortened together.

phi3 = C D f(r1 + r2) with C < 0 runs to minus infinity as the leg sum goes to
zero while phi2 resists only per bond, so squeezing the vertex is where the
three-body term wins.  Neither the uniform-compression scan nor the first-shell
screen produces this distortion, and a defect calculation can reach it: the
self-interstitial and the grain-boundary tests do.

For each record of every set: the lowest energy inside the core (legs at most
1.2 A; a negative energy at bonding distance is cohesion, not a hole), the
longest leg at which it is negative, and the barrier in front of it - the
highest energy at legs longer than the lowest point of the scan.  The scan grid
is that of the probe's first run, 2026-11-03, which this reproduces; it does not resolve barriers that sit between grid points.

    python squeeze_probe.py          # writes squeeze.json and library.json `squeeze`
    python squeeze_probe.py --merge  # squeeze.json -> library.json only, no scan

The energy is written out for the three atoms rather than taken from latdyn's
crystal sums: pair terms A-X, X-B, A-B, and three triplets - the vertex X with
legs (L, L) at 180 degrees, and each end with legs (L, 2L) at 0 degrees - each
phi3(ra + rb) h(cos) S(ra) S(rb), every term inside its own cut-off.  The
angular latdyn is used for every set: with lam2 = lam4 = 0 it is MAU exactly.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "angular"))
import latdyn as L                                   # noqa: E402

assert "lam2" in open(L.__file__, encoding="utf-8").read(), L.__file__
LS = np.round(np.arange(2.6, 0.399, -0.02), 3)
CORE = 1.2
SETS = {"mau": lambda v: v, "ug": lambda v: v.get("ug"), "tap": lambda v: v.get("tap"),
        "tap_ug": lambda v: v.get("tap_ug"), "rc": lambda v: v.get("rc"),
        "rc_ug": lambda v: v.get("rc_ug"), "tap_force": lambda v: v.get("tap_force"),
        #  the side candidates ship as files too (_disp, _nudge_taper)
        "hard_disp": lambda v: v.get("hard_disp"), "tap_nudge": lambda v: v.get("tap_nudge")}


def trimer(pot, Lg):
    def p2(r):
        return float(pot.phi2(r)) if r < pot.rcut2 else 0.0

    def t3(ra, rb, c):
        if not pot.C or not pot.rcut3 or ra >= pot.rcut3 or rb >= pot.rcut3:
            return 0.0
        g = float(pot.phi3(ra + rb)) * float(pot.ang(c))
        if pot.taper is not None:
            g *= float(pot.switch(ra, pot.rcut3, 0)) * float(pot.switch(rb, pot.rcut3, 0))
        return g
    return (2.0 * p2(Lg) + p2(2.0 * Lg)
            + t3(Lg, Lg, -1.0) + 2.0 * t3(Lg, 2.0 * Lg, 1.0))


def scan(rec):
    pot = L.Potential.from_record(rec)
    E = np.array([trimer(pot, x) for x in LS])
    core = LS <= CORE
    neg = [float(x) for x, e in zip(LS[core], E[core]) if e < 0]
    return {"hole": bool(neg), "L_negative": max(neg) if neg else None,
            "min_core_eV": float(E[core].min()),
            "barrier_eV": float(E[:int(np.argmin(E)) + 1].max()) if neg else None}


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    out = {}
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict):
            continue
        for key, get in SETS.items():
            rec = get(v)
            if not isinstance(rec, dict) or "m" not in rec:
                continue
            r = scan(rec)
            out[f"{el}|{key}"] = r
            rec["squeeze"] = r
    for key in SETS:
        rows = [r for k, r in out.items() if k.endswith("|" + key)]
        holes = [r for r in rows if r["hole"]]
        print(f"{key:10s} {len(holes):3d} of {len(rows):3d} with a hole; barrier <= 10 eV: "
              f"{sum(r['barrier_eV'] <= 10 for r in holes)}, <= 50 eV: "
              f"{sum(r['barrier_eV'] <= 50 for r in holes)}")
    json.dump(out, open(os.path.join(HERE, "squeeze.json"), "w"), indent=1, sort_keys=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)


def merge():
    """squeeze.json -> each record's `squeeze`, without scanning.  refresh.py
    runs this form: a clean clone holds only the root and UG records, and a
    scan there would rewrite squeeze.json with those two sets and drop the
    other seven, which the shipped file carries."""
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    src = json.load(open(os.path.join(HERE, "squeeze.json")))
    n = 0
    for key, r in src.items():
        el, arm = key.split("|", 1)
        rec = SETS[arm](lib[el]) if isinstance(lib.get(el), dict) else None
        if isinstance(rec, dict) and "m" in rec:
            rec["squeeze"] = r
            n += 1
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"squeeze merged into {n} of {len(src)} records")


if __name__ == "__main__":
    if "--merge" in sys.argv:
        merge()
    else:
        main()
