#!/usr/bin/env python3
"""q-points for the finite-temperature dispersion CURVE, one file per element.

The panel shows a number; this is what the number is a summary of.  The curve
has to lie on the same path as the 0 K one already stored in `ld.std`, or the
two cannot be drawn on one axis, so the segments and their order are taken
from there rather than rebuilt - a path assembled twice is a path that will
differ once.

Subsampled to NPT points per segment against ld.std's 60.  At panel width a
segment is about sixty pixels, so sixty points would be one per pixel; twenty
is indistinguishable and keeps the page a hundred kilobytes lighter.

    python emit_ftcurve.py            # -> ftq_<El>.txt and ftq_index.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from build_library import sc_segments          # noqa: E402

NPT = 20


def main():
    lib = json.load(open(os.path.join(HERE, "library.json")))
    md = json.load(open(os.path.join(HERE, "mdres.json")))
    index = {}
    for el in sorted(md):
        v = lib[el]
        std = (v.get("ld") or {}).get("std")
        if not std:
            continue
        ends = {(a, b): (np.asarray(ka, float), np.asarray(kb, float))
                for a, ka, b, kb in sc_segments(v["struct"])}
        qs, segs = [], []
        for sg in std:
            key = (sg["a"], sg["b"])
            if key not in ends:
                #  drawn one way, tabulated the other; the path is the same
                key = (sg["b"], sg["a"])
                if key not in ends:
                    continue
                ka, kb = ends[key][1], ends[key][0]
            else:
                ka, kb = ends[key]
            for i in range(NPT):
                qs.append(ka + (i / (NPT - 1.0)) * (kb - ka))
            segs.append({"a": sg["a"], "b": sg["b"], "len": sg["len"]})
        with open(os.path.join(HERE, "ftq_%s.txt" % el), "w") as fh:
            for q in qs:
                fh.write("%.10f %.10f %.10f\n" % tuple(q))
        index[el] = {"segs": segs, "npt": NPT, "nq": len(qs),
                     "struct": v["struct"]}
        print("%-3s %2d segment x %d = %4d q" % (el, len(segs), NPT, len(qs)))
    json.dump(index, open(os.path.join(HERE, "ftq_index.json"), "w"), indent=0)
    print("\n%d element, ftq_index.json yazildi" % len(index))


if __name__ == "__main__":
    main()
