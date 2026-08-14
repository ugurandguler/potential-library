#!/usr/bin/env python3
"""
Is the analytic force the derivative of the energy?  Finite differences say so.

Chromium and tungsten blow up under `pair_style ugur/ang` - 5204 K and 10157 K
out of a 600 K start - while aluminium and vanadium conserve to hundredths of a
meV/atom/ps with weights of the same size.  Molecular dynamics cannot tell the
two possible causes apart: a wrong force and a genuinely unstable potential
both look like an explosion.

This can.  Displace one atom by +h and -h along one axis, read the energy from
LAMMPS each time, and compare -(E(+h) - E(-h)) / 2h against the force LAMMPS
reports at h = 0.  Both numbers come from LAMMPS, so what is being tested is
the pair style's own internal consistency: if the force is the gradient of the
energy the potential is implemented correctly, whatever the dynamics then do
with it.  The energy itself is already known to be right - it matches
angular/latdyn to 1e-12 on all four elements.

If they agree, the code is right and the instability is the potential's.  If
they disagree, the angular force term is wrong and the two elements that
survive were lucky.

    python fd_force.py            # the four with nonzero weights
    python fd_force.py Cr
"""
import io
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import cellfile         # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"

POTFILE = """# {el}
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

#  One atom is nudged off its site so the forces are not all zero by symmetry -
#  on a perfect lattice every force vanishes and the test measures nothing.
IN = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       {data}
pair_style      {style}
pair_coeff      * * {pot} {el}
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check no
displace_atoms  all random 0.08 0.08 0.08 4321 units box
group           probe id 1
displace_atoms  probe move {dx} {dy} {dz} units box
run             0
variable        ee equal pe
variable        fx equal fx[1]
variable        fy equal fy[1]
variable        fz equal fz[1]
print           "OUT ${{ee}} ${{fx}} ${{fy}} ${{fz}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def one(d, el, style, dx, dy, dz):
    open(os.path.join(d, "in.fd"), "w").write(IN.format(
        data=f"{el}.data", style=style, pot=f"{el}.ugur", el=el,
        dx=dx, dy=dy, dz=dz))
    out = subprocess.run(["wsl", "-e", "bash", "-lc",
                          f"cd {wsl(d)} && {LMP} -in in.fd 2>&1"],
                         capture_output=True, text=True).stdout
    m = re.search(r"OUT\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)"
                  r"\s+([-\d.eE+]+)", out)
    if not m:
        err = [l for l in out.splitlines() if "ERROR" in l]
        raise SystemExit(f"{el}: {(err or ['no output'])[0][:80]}")
    return [float(x) for x in m.groups()]


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or ["Al", "Cr", "V", "W"]
    h = 1e-4
    print("Analytic force against the numerical energy derivative (h = 1e-4 A)")
    print("Same atom, same configuration; both from LAMMPS.\n")
    print(f"{'el':4s}{'style':>10s}{'axis':>7s}{'F analytic':>14s}"
          f"{'F numeric':>14s}{'rel diff':>13s}   status")
    print("-" * 68)
    bad = 0
    for el in els:
        rec = lib[el].get("tap_ug")
        if not rec:
            print(f"{el:4s}  no tapered UG record")
            continue
        for style, params in (("ugur/ang", rec),
                              ("ugur", {**rec, "lam2": 0.0, "lam4": 0.0})):
            d = os.path.join(HERE, "fdruns", f"{el}_{style.replace('/', '_')}")
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, f"{el}.ugur"), "w").write(
                POTFILE.format(el=el, **{k: params[k] for k in
                                         ("m", "D", "alpha", "r0", "gamma",
                                          "C", "alpha3", "rcut2", "rcut3",
                                          "taper", "lam2", "lam4")}))
            e = refdata.ELEMENTS[el]
            cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                            mass=refdata.MASSES[el])
            box, _ = cellfile.orthogonal_cell(cry)
            rc = max(rec["rcut2"], rec["rcut3"])
            rep = tuple(max(3, int(np.ceil(1.6 * (rc + 2.0) / b)))
                        for b in box)
            cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                                refdata.MASSES[el], rep)
            base = one(d, el, style, 0.0, 0.0, 0.0)
            for ax, i in (("x", 1), ("y", 2), ("z", 3)):
                dp = [0.0, 0.0, 0.0]
                dp[i - 1] = h
                ep = one(d, el, style, *dp)[0]
                dp[i - 1] = -h
                em = one(d, el, style, *dp)[0]
                fnum = -(ep - em) / (2 * h)
                fana = base[i]
                scale = max(abs(fana), abs(fnum), 1e-6)
                rel = abs(fana - fnum) / scale
                ok = rel < 2e-4
                bad += not ok
                print(f"{el if ax == 'x' and style == 'ugur/ang' else '':4s}"
                      f"{style if ax == 'x' else '':>10s}{ax:>7s}"
                      f"{fana:14.6f}{fnum:14.6f}{rel:13.2e}"
                      f"   {'ok' if ok else 'DIFFERS'}")
        print()
    if bad:
        raise SystemExit(f"{bad} bilesende kuvvet enerjinin turevi degil")
    print("kuvvet her yerde enerjinin turevi - uygulama tutarli")


if __name__ == "__main__":
    main()
