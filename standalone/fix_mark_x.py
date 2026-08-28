#!/usr/bin/env python3
"""
Put the high-symmetry lines where the path actually crosses those points.

The reference dispersions are stored as a list of q-points, and the viewer
draws one x per INDEX.  A label is carried as an index too, so a vertical line
lands on whichever sample was nearest - and the samples do not have to include
the point they bracket.  On the Materials Cloud path for silver the interior
Gamma is never sampled at all: index 48 sits at (0.0117, 0.0117, 0.0235) on
the way in and index 49 at (0.0096, 0.0096, 0.0096) on the way out.  The line
was drawn at 48, a little to the left of Gamma, where the acoustic branches
have not reached zero yet - so the plot appeared to show a potential whose
frequencies do not vanish at Gamma, which is not what it computes.

fix_mp_path records the same failure for MP's endpoints at a discontinuity;
this is the interior version of it.  137 of the 232 labels across the three
reference sources sit off their own symmetry point.

What is stored is a FRACTIONAL index, `marks_x`, found by minimising the
distance to the symmetry point along the two segments that meet at the label.
Segments that span a discontinuity are skipped, so nothing is interpolated
across a jump.  The original `marks` stays untouched.

    python fix_mark_x.py
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from build_library import SC_POINTS          # noqa: E402

SRC = ("mp", "mc3d", "jarvis")


def vertex_x(q, brk, i, P):
    """fractional index of the symmetry point P, which the path may not sample

    P is a VERTEX: the path arrives along one line and leaves along another,
    so the straight chord between the two samples that bracket it does NOT
    pass through it and looking for the nearest point on that chord finds an
    endpoint instead - which is how the first version of this file moved
    labels by exactly one index and no fractions.  What locates a vertex is
    how far each neighbour still is from it along its own leg: with a roughly
    even sampling step, the crossing sits between them in that ratio.
    """
    d = np.linalg.norm(q - P, axis=1)
    #  LOCALLY.  A symmetry point appears more than once on these paths -
    #  Gamma twice on the fcc one, M twice on the hcp one - so a global
    #  nearest sample can be the other occurrence entirely, half a path away.
    #  The stored index is already right to within a sample or two; this only
    #  refines it.
    lo, hi = max(0, i - 3), min(len(q), i + 4)
    m = lo + int(np.argmin(d[lo:hi]))
    if d[m] < 1e-6:
        return float(m)
    #  Which NEIGHBOUR the point sits between is not decided by which is
    #  closer.  On the Materials Cloud hcp path yttrium's interior Gamma has
    #  samples at |q| = 0.0419 and 0.0105 coming in along [110] and 0.0476
    #  going out along [001]: the nearer neighbour of the nearest sample is
    #  the one BEHIND it, on the same leg, and taking that pair put the label
    #  at 34.80 when Gamma is at 35.18 - on the wrong side of the sample it
    #  was meant to bracket.  Yttrium, beryllium and hafnium were all half an
    #  index out that way.
    #
    #  What brackets a vertex is a pair whose two samples leave it in
    #  DIFFERENT directions.  Two samples on one leg point the same way, so
    #  cos = +1 and the pair is rejected; a pair straddling the vertex turns a
    #  corner (cos ~ 0 here) or reverses if the point is interior (cos = -1).
    cand = []
    for a, b in ((m - 1, m), (m, m + 1)):
        if a < 0 or b >= len(q) or a in brk:
            continue      # the path really is broken there, nothing between
        if d[a] < 1e-12 or d[b] < 1e-12:
            cand.append((-1.0, a, b))
            continue
        cand.append((float(np.dot((q[a] - P) / d[a], (q[b] - P) / d[b])),
                     a, b))
    if not cand:
        return float(m)
    straddle = [c for c in cand if c[0] < 1.0 - 1e-6]
    _, a, b = min(straddle) if straddle else min(cand)
    tot = d[a] + d[b]
    return float(a) if tot == 0 else float(a) + float(d[a] / tot)


def x_index(x, brk, i, X):
    """fractional index at which the path coordinate reaches X

    For JARVIS there are no q-points stored, but the source declares its own
    label positions on its own x axis, which is better than inferring them:
    the answer is read off the source rather than reconstructed.  All that is
    needed is to turn an x into the index the viewer draws by.
    """
    lo, hi = max(0, i - 3), min(len(x) - 1, i + 3)
    for a in range(lo, hi):
        if a in brk:                 #  a -> a+1 crosses a discontinuity
            continue
        xa, xb = float(x[a]), float(x[a + 1])
        if xa == xb:
            continue
        if min(xa, xb) - 1e-12 <= X <= max(xa, xb) + 1e-12:
            return float(a) + (X - xa) / (xb - xa)
    return float(i)


def declared_x(r, i, name):
    """the source's own x for the label `name` nearest sample i, or None"""
    labs, lx = r.get("labels"), r.get("label_x")
    if not labs or not lx or len(labs) != len(lx):
        return None
    here = float(r["x"][i])
    cand = [float(X) for L, X in zip(labs, lx) if L == name]
    return min(cand, key=lambda X: abs(X - here)) if cand else None


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    moved = kept = 0
    worst = []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict):
            continue
        pts = SC_POINTS.get(v.get("struct")) or {}
        for src in SRC:
            r = v.get(src)
            if not isinstance(r, dict) or not r.get("marks"):
                continue
            has_q = bool(r.get("q"))
            if not has_q and not r.get("x"):
                continue
            q = np.array(r["q"], dtype=float) if has_q else None
            brk = set(r.get("breaks") or [])
            out = []
            for i, name in r["marks"]:
                P = pts.get(name)
                n = len(q) if has_q else len(r["x"])
                if i >= n:
                    out.append([float(i), name])
                    continue
                if has_q:
                    if P is None:
                        out.append([float(i), name])
                        continue
                    xf = vertex_x(q, brk, i, np.array(P, dtype=float))
                else:
                    X = declared_x(r, i, name)
                    if X is None:
                        out.append([float(i), name])
                        continue
                    xf = x_index(np.asarray(r["x"], float), brk, i, X)
                out.append([round(xf, 4), name])
                if abs(xf - i) > 1e-3:
                    moved += 1
                    worst.append((abs(xf - i), el, src, name, i, xf))
                else:
                    kept += 1
            r["marks_x"] = out
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    worst.sort(reverse=True)
    print(f"{moved} labels moved, {kept} already on their point")
    for s, el, src, name, i, x in worst[:8]:
        print(f"   {el:3s} {src:8s} {name:2s}  {i} -> {x:.3f}")


if __name__ == "__main__":
    main()
