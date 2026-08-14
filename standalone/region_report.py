#!/usr/bin/env python3
"""
Does reachability in the (C', C44) plane predict the fit quality?

cprime_ceiling.py records, for each cubic element, every (C'/B, C44/B) the form
can produce under the hard constraints.  The ratio test in that script is only
necessary, not sufficient: an element also needs the right *magnitude*, so the
honest test is whether the experimental point lies inside the cloud itself.

Distance is measured relative to the experimental point, so it reads as a
fractional error: 0 means reachable, 0.4 means the closest the form can get is
40 % away in the (C', C44) plane.

    python region_report.py
"""
import glob
import json
import os

import numpy as np

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


def fits():
    """best known RMS per element, from the archived runs and the live fit"""
    best = {}
    files = ([os.path.join(HERE, "fit.json")] +
             glob.glob(os.path.join(HERE, "runs", "*", "fit.json")) +
             glob.glob(os.path.join(HERE, "runs", "*", "dense_*.json")))
    for p in files:
        if not os.path.exists(p):
            continue
        for el, v in json.load(open(p)).items():
            if isinstance(v, dict) and "score" in v:
                if el not in best or v["score"] < best[el]:
                    best[el] = v["score"]
    return best


def main():
    p = os.path.join(HERE, "cprime_region.json")
    region = json.load(open(p))
    rms = fits()

    print(f"{'el':4s}{'struct':>7s}{'R_exp':>8s}{'R_min':>8s}"
          f"{'ratio ok':>10s}{'2D gap':>9s}{'RMS %':>9s}")
    print("-" * 55)

    rows = []
    for el in sorted(region, key=lambda e: -region[e]["R_exp"]):
        d = region[el]
        e = refdata.ELEMENTS[el]
        c = e["Cij"]
        B = d["B"]
        target = np.array([0.5 * (c["C11"] - c["C12"]) / B, c["C44"] / B])
        pts = np.array(d["pts"])
        gap = float(np.min(np.linalg.norm(pts - target, axis=1))
                    / np.linalg.norm(target))
        ok = d["R_min"] <= d["R_exp"] <= d["R_max"]
        r = rms.get(el)
        rows.append((el, gap, r, ok))
        print(f"{el:4s}{e['struct']:>7s}{d['R_exp']:8.2f}{d['R_min']:8.2f}"
              f"{('yes' if ok else 'NO'):>10s}{gap:9.2f}"
              f"{(f'{100 * r:.1f}' if r else '-'):>9s}")

    have = [(g, r, el) for el, g, r, _ in rows if r is not None]
    g = np.array([x[0] for x in have])
    r = np.array([x[1] for x in have]) * 100
    print(f"\ncorrelation, 2D gap vs RMS: "
          f"{np.corrcoef(g, np.log10(r + 1))[0, 1]:+.2f}   (n = {len(g)})")

    #  0.4 is where the coarse grid stops being the limiting factor.  Below it
    #  the gap is dominated by how finely the parameter space was sampled - K
    #  reports 0.21 and still fits perfectly - while above it the ordering is
    #  essentially monotone in the RMS.  Every gap here is an upper bound:
    #  a denser grid can only move points closer to the target.
    near = [el for gg, _, el in have if gg < 0.40]
    far = [el for gg, _, el in have if gg >= 0.40]
    for name, group in (("reachable (gap < 0.40)", near),
                        ("out of reach (gap >= 0.40)", far)):
        vals = [100 * rms[el] for el in group]
        print(f"   {name:28s} n={len(group):2d}  median RMS "
              f"{np.median(vals):6.1f} %   {group}")

    print("\nElements that are reachable but still fit badly are search "
          "failures,\nnot form failures - those are the ones worth more "
          "restarts.")


if __name__ == "__main__":
    main()
