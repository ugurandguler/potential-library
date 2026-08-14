#!/usr/bin/env python3
"""
Is the tapered drift on the light hcp metals the potential or the timestep?

Titanium switched drifts 17.28 meV/atom/ps at 2 fs and yttrium 5.73, against
0.03 to 0.28 for the cubic set.  Read as a property of the potential that would
say the switch works badly on hcp.  But the four hcp elements order by mass -
Ti 47.9, Y 88.9, Ru 101.1, Re 186.2 - and so does the drift, which is what
integrator error does: a lighter atom vibrates faster and 2 fs stops resolving
it.  Verlet error scales as dt^2, so quartering the step must divide the drift
by sixteen if that is the cause and leave it alone if it is not.
"""
import sys
import numpy as np
import nve_check as N
import json, os
sys.path.insert(0, os.path.join(N.ROOT, "standalone"))
import refdata

lib = json.load(open(os.path.join(N.ROOT, "standalone", "library.json")))
print("tapered run, same physical time (40 ps), different step")
print(f"{'el':4s}{'mass':>7s}{'2 fs':>10s}{'0.5 fs':>10s}{'ratio':>8s}"
      f"{'dt^2 beklenen':>15s}   verdict")
print("-" * 64)
for el in sys.argv[1:] or ["Ti", "Y"]:
    e = refdata.ELEMENTS[el]
    tap = dict(lib[el]["tap"])
    w = tap.get("taper") or 0.85
    d2, T2, err2, _ = N.run(el, tap, e, w, dt=0.002, neq=2000, nrun=20000)
    #  same 40 ps and the same 4 ps of equilibration, four times the steps
    d5, T5, err5, _ = N.run(el, tap, e, w, dt=0.0005, neq=8000, nrun=80000,
                            every=2000)
    if d2 is None or d5 is None:
        print(f"{el:4s}  error: {err2 or err5}")
        continue
    ratio = abs(d2) / abs(d5) if d5 else float("inf")
    verdict = ("integrator hatasi" if 8 <= ratio <= 32
               else "potansiyel" if ratio < 3 else "belirsiz")
    print(f"{el:4s}{refdata.MASSES[el]:7.1f}{d2:10.3f}{d5:10.3f}"
          f"{ratio:8.1f}{16.0:15.1f}   {verdict}")
