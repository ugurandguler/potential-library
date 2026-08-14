#!/usr/bin/env python3
"""
Compare archived fits and price the dynamical-stability constraint.

Read the per-run `dense_<El>_<n>.json` files, NOT `fit.json`.  `fit.py` merges
into `fit.json` keeping the better score per element, so a run that comes out
worse everywhere leaves that file untouched — archiving it captures the previous
run, not the one that just finished.  The dense files are written per run and
are the only per-run record.

The only defensible comparison is between runs that differ in one thing.  The
first constrained run used 20 restarts while the baseline had 60, so part of any
difference was search quality rather than the constraint; the control run repeats
the baseline at 20 restarts so that confound is removed.

For every run it also re-measures dynamical stability from scratch, because a
fit that was never asked for it can easily have imaginary modes.

    python compare_runs.py
"""
import glob
import sys
import json
import os

import numpy as np

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
NQ = 6


def load_run(d):
    """Per-element results for one run, preferring the dense per-run files."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "dense_*.json"))):
        data = json.load(open(p))
        if data.get("failed"):
            out[data["element"]] = None
            continue
        for el, v in data.items():
            if isinstance(v, dict) and "score" in v:
                out[el] = v
    if out:
        return out, "dense"
    p = os.path.join(d, "fit.json")
    if os.path.exists(p):
        data = json.load(open(p))
        return ({el: v for el, v in data.items()
                 if isinstance(v, dict) and "score" in v}, "fit.json")
    return {}, None


def imaginary_fraction(el, v):
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    pot = L.Potential.from_record(v)
    f = L.spectrum(cry, pot, nq=NQ)
    return float((f < -1e-6).mean())


def main():
    runs = {}
    for name in sorted(os.listdir(RUNS)):
        d = os.path.join(RUNS, name)
        if not os.path.isdir(d):
            continue
        r, src = load_run(d)
        if r:
            runs[name] = r
            print(f"   {name:34s} {len(r):3d} elements   (from {src})")

    #  Substring matching picked "2026-07-31_dynstable" over
    #  "2026-07-31_dynstable-8cubed" and silently re-compared the superseded
    #  run, so take the LAST match (runs are date-prefixed, so that is the most
    #  recent) and allow an explicit override.
    def pick(tag, argv_pos):
        if len(sys.argv) > argv_pos:
            name = sys.argv[argv_pos]
            if name not in runs:
                raise SystemExit(f"no such run: {name}")
            return name
        hits = [k for k in runs if tag in k]
        return hits[-1] if hits else None

    ctrl = pick("control", 1)
    cons = pick("dynstable", 2)
    base = next((k for k in runs if "elastic-only" in k), None)
    if not (ctrl and cons):
        raise SystemExit("need both a control and a constrained run")
    print(f"\ncomparing control={ctrl}  constrained={cons}")

    cols = [k for k in (base, ctrl, cons) if k]
    print(f"\n{'el':4s}", end="")
    for k in cols:
        print(f"{k.split('_')[-1][:14]:>16s}", end="")
    print(f"{'cost':>8s}   imaginary modes (control -> constrained)")
    print("-" * 92)

    costs, fixed, broke, missing = [], [], [], []
    for el in sorted(refdata.ELEMENTS):
        if not any(el in runs[k] for k in cols):
            continue
        row = f"{el:4s}"
        for k in cols:
            v = runs[k].get(el, "absent")
            if v == "absent":
                row += f"{'-':>16s}"
            elif v is None:
                row += f"{'FAILED':>16s}"
            else:
                row += f"{v['score'] * 100:16.1f}"

        a, b = runs[ctrl].get(el), runs[cons].get(el)
        tail = ""
        if a is None or b is None:
            tail = "   (no solution)"
        elif a and b:
            c = (b["score"] - a["score"]) * 100
            costs.append((c, el))
            fa, fb = imaginary_fraction(el, a), imaginary_fraction(el, b)
            tail = f"{c:+8.1f}   {100 * fa:5.1f}% -> {100 * fb:5.1f}%"
            if fa > 0 and fb == 0:
                fixed.append(el)
            if fa == 0 and fb > 0:
                broke.append(el)
        else:
            missing.append(el)
            tail = "   (not in both runs)"
        print(row + tail)

    arr = np.array([c for c, _ in costs])
    print(f"\ncost of requiring dynamical stability, matched at 20 restarts, "
          f"{len(costs)} elements:")
    print(f"   median {np.median(arr):+.2f} points, mean {arr.mean():+.2f}, "
          f"worst {max(costs)[0]:+.1f} ({max(costs)[1]})")
    print(f"   unaffected (|cost| < 0.5 points): "
          f"{sum(1 for c, _ in costs if abs(c) < 0.5)}/{len(costs)}")
    print(f"   made dynamically stable: {fixed or 'none'}")
    if broke:
        print(f"   NEWLY unstable (should not happen): {broke}")
    if missing:
        print(f"   skipped, absent from one run: {missing}")

    if base:
        print("\nsearch effort, not the constraint: baseline used 60 restarts, "
              "both new runs 20.")
        d60 = [(runs[ctrl][el]["score"] - runs[base][el]["score"]) * 100
               for el in runs[ctrl]
               if runs[ctrl].get(el) and runs[base].get(el)]
        d60 = np.array(d60)
        print(f"   60 -> 20 restarts, constraint off both times: "
              f"median {np.median(d60):+.2f} points, worst {d60.max():+.1f}")


if __name__ == "__main__":
    main()
