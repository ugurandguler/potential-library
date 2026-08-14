#!/usr/bin/env python3
"""
Merge the JARVIS-DFT dispersions into library.json, but only the trustworthy
ones.

JARVIS is a third independent DFT opinion and mostly a good one - against
Materials Cloud MC3D it lands within about ten per cent, which is ordinary
between two independent calculations. But not every entry is usable, and the
failures are not random:

  * **Iron** comes out at 121.5 cm^-1 against a measured 285 (the H point,
    Landolt-Boernstein III/13a). Non-spin-polarised iron has far softer phonons
    than the real ferromagnet, and that is almost certainly what this is.
  * **Zinc** is a factor of two below MC3D while keeping perfectly correct
    symmetry behaviour, so it is not a broken run - it is two calculations
    disagreeing about the one element whose c/a sits 14 % off ideal.

So each entry is checked against whatever else we already have, and anything
more than 25 % away from an existing reference is stored but flagged rather
than silently drawn. Better a visible gap than a wrong curve.

    python add_jarvis_overlay.py
"""
import json
import os

import numpy as np

import build_library as B
import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 0.25            # relative disagreement with an existing reference
DENSITY_FLOOR = 0.25  # a segment sampled this much thinner than the rest is a jump
TOL_X = 1e-9          # path-length comparisons


def existing_max(v, el, exp):
    """highest frequency from any reference we already trust, cm^-1"""
    out = []
    mp = (v.get("mp") or {}).get("phonon", {}).get("stats", {}).get("mp_max")
    if mp:
        out.append(("MP", mp))
    if v.get("mc3d"):
        out.append(("MC3D", max(max(r) for r in v["mc3d"]["f"])))
    if exp:
        out.append(("measured", max(max(p["exp"]) for p in exp["points"])
                    * 33.35641))
    return out


def path_geometry(rec, struct):
    """
    (q, keep, marks, breaks) rebuilt from the path length, or None.

    JARVIS stores the distance along the path and the distance of each
    high-symmetry label, but not the q-vectors.  The label sequence is the same
    path we use, so every stored distance maps back by interpolating between the
    two labels that bracket it - which is what lets our potential be evaluated at
    THEIR q-points rather than against an interpolation of their curve.

    Discontinuities have to be found rather than assumed, and they show up in the
    sampling density: on fcc copper the real segments carry 165-195 points per
    unit path length and K->U carries 6.  Those segments are jumps, the same ones
    our own path declares, and the handful of stray samples inside them are
    dropped - their q-vector would be an interpolation straight across a gap, so
    neither the curve nor the residual there means anything.
    """
    pts = B.SC_POINTS[struct]
    seq = []
    for lab, x in zip(rec["labels"], rec["label_x"]):
        if lab not in pts:
            return None
        if not seq or (lab, x) != seq[-1]:
            seq.append((lab, x))
    x = np.asarray(rec["x"], float)
    segs = [(k, k + 1) for k in range(len(seq) - 1)
            if seq[k + 1][1] > seq[k][1]]
    if not segs:
        return None

    inside = [(x > seq[a][1] + TOL_X) & (x < seq[b][1] - TOL_X)
              for a, b in segs]
    dens = [int(m.sum()) / (seq[b][1] - seq[a][1])
            for (a, b), m in zip(segs, inside)]
    med = float(np.median(dens))
    jump = [d < DENSITY_FLOOR * med for d in dens]

    keep = np.ones(len(x), bool)
    for m, j in zip(inside, jump):
        if j:
            keep &= ~m
    xk = x[keep]

    q = []
    for xv in xk:
        hit = None
        for (a, b), j in zip(segs, jump):
            if j:
                continue
            xa, xb = seq[a][1], seq[b][1]
            if xa - TOL_X <= xv <= xb + TOL_X:
                hit = (seq[a][0], xa, seq[b][0], xb)
                break
        if hit is None:
            return None
        la, xa, lb, xb = hit
        t = (xv - xa) / (xb - xa)
        p, r = np.array(pts[la], float), np.array(pts[lb], float)
        q.append(p + t * (r - p))

    #  A label that opens a run after a jump takes the first kept sample at or
    #  past it; every other label takes the last one at or before it.  At a jump
    #  that puts the two labels on adjacent indices, which is what the viewer
    #  needs to draw them either side of the gap.
    opens = {b for (a, b), j in zip(segs, jump) if j}
    marks, breaks = [], []
    for k, (lab, xl) in enumerate(seq):
        if k == 0 or k in opens:
            w = np.flatnonzero(xk >= xl - TOL_X)
        else:
            w = np.flatnonzero(xk <= xl + TOL_X)
        if not len(w):
            return None
        marks.append([int(w[0] if (k == 0 or k in opens) else w[-1]), lab])
    for (a, b), j in zip(segs, jump):
        if j:
            breaks.append(marks[a][0])
    return np.array(q), keep, marks, sorted(set(breaks))


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    src = json.load(open(os.path.join(HERE, "jarvis_phonon.json")))

    print(f"{'el':4s}{'JARVIS':>9s}{'compared with':>26s}{'ratio':>8s}   verdict")
    print("-" * 62)
    kept = flagged = 0
    for el in sorted(src):
        v = lib.get(el)
        if not v:
            continue                      # Cd and Zn have no fit
        rec = dict(src[el])
        top = max(max(r) for r in rec["f"])
        #  our own dispersion at THEIR q-points, the same treatment MP and MC3D
        #  get, so any residual is physics and not interpolation
        geo = path_geometry(rec, v["struct"])
        if geo is not None:
            q, keep, marks, breaks = geo
            rec["x"] = [xx for xx, k in zip(rec["x"], keep) if k]
            rec["f"] = [[vv for vv, k in zip(b, keep) if k] for b in rec["f"]]
            rec["marks"], rec["breaks"] = marks, breaks
            e = refdata.ELEMENTS[el]
            cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                            mass=refdata.MASSES[el])
            pot = L.Potential.from_record(v)
            Phi = L.force_constants(cry, pot)
            ours = np.array([np.sort(L.frequencies(cry, pot, qq, Phi)) * L.CM1
                             for qq in q])
            ref = np.array([sorted(col) for col in zip(*rec["f"])])
            nb = min(ours.shape[1], ref.shape[1])
            rec["ours"] = [[round(float(x), 1) for x in ours[:, k]]
                           for k in range(nb)]
            rec["f"] = [[round(float(x), 1) for x in ref[:, k]]
                        for k in range(nb)]
            d = ours[:, :nb] - ref[:, :nb]
            rec["stats"] = {"rms_cm1": round(float(np.sqrt((d**2).mean())), 2),
                            "ours_max": round(float(ours.max()), 1),
                            "ref_max": round(float(ref.max()), 1)}
        refs = existing_max(v, el, v.get("exp_phonon"))
        worst = None
        for name, other in refs:
            r = top / other if other else None
            if r and (worst is None or abs(np.log(r)) > abs(np.log(worst[1]))):
                worst = (name, r)
        rec["top_cm1"] = round(float(top), 1)
        #  "cleared" and "never tested" are not the same thing.  Al Ba Co Cr Na
        #  have no other phonon reference at all, so the 25 % screen never ran
        #  on them - they are kept because nothing contradicted them, which is
        #  weaker than agreeing with something, and the viewer says so.
        rec["checked"] = worst is not None
        if worst is None:
            rec["trusted"] = True
            note = "no other reference, kept unchecked"
        elif abs(worst[1] - 1.0) > TOL:
            rec["trusted"] = False
            rec["conflict"] = {"against": worst[0], "ratio": round(worst[1], 2)}
            note = "FLAGGED, not drawn"
            flagged += 1
        else:
            rec["trusted"] = True
            note = "kept"
        if rec["trusted"]:
            kept += 1
        v["jarvis"] = rec
        #  name and value must come from the SAME reference as the ratio
        shown = next((f"{n} {x:.1f}" for n, x in refs if n == worst[0]), "-")             if worst else "-"
        print(f"{el:4s}{top:9.1f}{shown:>26s}"
              f"{(worst[1] if worst else float('nan')):8.2f}   {note}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\n{kept} kept, {flagged} flagged as conflicting")
    print("merged into", path)


if __name__ == "__main__":
    main()
