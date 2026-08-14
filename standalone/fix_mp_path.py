#!/usr/bin/env python3
"""
Repair the high-symmetry labels on the Materials Project dispersions, and
record where their path is genuinely discontinuous.

fetch_mp.py located labels by matching q-point coordinates against
bs.labels_dict.  That silently loses the labels on any branch that follows a
discontinuity, because MP's q-list drops the endpoints there: on the hcp L-M
branch the samples run (0.5, 0, 0.495) to (0.5, 0, 0.005) and never touch L or
M.  Mg therefore came back with ten labels for a twelve-point path, and the plot
ran straight from A to K with no break drawn - looking as if the segment between
them were missing data when in fact it is two segments with a jump between them.

The reliable signal is geometric: a discontinuity is a jump between consecutive
samples far larger than the sampling step.  Those, plus the coordinate matches,
give exactly as many boundaries as `kpath` has names, so the names can simply be
assigned in order.

Adds to each phonon record:
    breaks   indices i where the path jumps between sample i and i+1
    marks    [index, name], now complete

    python fix_mp_path.py
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JUMP = 6.0          # a break is this many times the median sampling step


def boundaries(q, marks):
    """(sorted boundary indices, break indices)"""
    dq = np.linalg.norm(np.diff(q, axis=0), axis=1)
    med = float(np.median(dq))
    brk = [int(i) for i in np.where(dq > JUMP * med)[0]]
    b = {0, len(q) - 1}
    b.update(i for i, _ in marks)
    for i in brk:
        b.update((i, i + 1))
    return sorted(b), brk


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print(f"{'el':4s}{'struct':>7s}{'bounds':>8s}{'kpath':>7s}"
          f"{'breaks':>9s}  labels")
    print("-" * 78)
    fixed = 0
    for el in sorted(lib):
        mp = (lib[el].get("mp") or {}).get("phonon")
        if not mp:
            continue
        q = np.array(mp["q"])
        bnd, brk = boundaries(q, mp.get("marks", []))
        kp = mp.get("kpath") or []
        note = ""
        if kp and len(bnd) == len(kp):
            mp["marks"] = [[int(i), n] for i, n in zip(bnd, kp)]
            fixed += 1
        else:
            #  Be's record carries no kpath at all, and a mismatch would mean
            #  the geometry and the label list disagree - keep what was fetched
            #  rather than invent names.
            note = "  kept fetched labels (no kpath)" if not kp else \
                   "  MISMATCH, labels left alone"
        mp["breaks"] = brk
        print(f"{el:4s}{lib[el]['struct']:>7s}{len(bnd):8d}{len(kp):7d}"
              f"{str(brk):>9s}  {'-'.join(n for _, n in mp['marks'])}{note}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\nrelabelled {fixed} elements; breaks recorded for all")


if __name__ == "__main__":
    main()
