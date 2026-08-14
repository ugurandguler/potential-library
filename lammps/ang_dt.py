#!/usr/bin/env python3
"""Cr and W blow up under ugur/ang at 2 fs.  Is that the timestep?

The force is already known to be the gradient of the energy (fd_force.py), so
this is not a code error.  Chromium is the stiffest element in the library and
2 fs is documented as too long for it even without the angular term.  If the
blow-up is the integrator, quartering the step removes it; if it is the
potential, a shorter step only postpones it.
"""
import io, json, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_ang as V
import numpy as np, latdyn as L, refdata, cellfile

lib = json.load(open(os.path.join(V.ROOT, "standalone", "library.json")))
TPL = V.IN_STATIC.replace("""run             0
variable        epa equal pe/atoms
print           "UGUR_E ${{epa}}"
""", """velocity        all create 600.0 12345 mom yes rot yes
fix             1 all nve
timestep        {dt}
thermo          {ev}
thermo_style    custom step temp pe etotal
thermo_modify   norm yes
run             {neq}
reset_timestep  0
run             {nrun}
""")

print("ugur/ang, same 40 ps, different step.  T should be 300 K.\n")
print(f"{'el':4s}{'step':>8s}{'drift':>12s}{'T avg':>8s}   status")
print("-" * 40)
for el in sys.argv[1:] or ["Cr", "W"]:
    rec = lib[el]["tap_ug"]
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    for dt, neq, nrun, ev in ((0.002, 2000, 20000, 500),
                              (0.0005, 8000, 80000, 2000)):
        tpl = TPL.replace("{dt}", str(dt)).replace("{neq}", str(neq)) \
                 .replace("{nrun}", str(nrun)).replace("{ev}", str(ev))
        out = V.run(el, rec, cry, f"dt{int(dt*10000)}", tpl)
        d, T = V.drift(out, dt=dt, nrun=nrun)
        if d is None:
            print(f"{el:4s}{dt*1000:7.1f}f   unreadable")
            continue
        print(f"{el:4s}{dt*1000:7.1f}f{d:12.3f}{T:8.0f}   "
              f"{'ok' if abs(d) < 1.0 and T < 400 else 'BLOWS UP'}")
