#!/usr/bin/env python3
"""
The energy along the volume-conserving tetragonal strain, which is what C'
measures - and what vanadium's finite-temperature collapse is a collapse of.

C' = (C11-C12)/2 is the curvature of E(delta) at delta=0 for the strain
(1+d, 1+d, 1/(1+d)^2).  The finite-temperature runs are NVT at fixed volume, so
thermal expansion cannot be the mechanism: whatever softens C' has to be
already present in this curve, away from the origin.

The distinction the sweep cannot draw, and this can:

  a real elastic constant       E(delta) is a well with one minimum, and the
                                quartic correction is small, so C'(T) falls
                                slowly and smoothly
  a fitted curvature on a       E(delta) curves the right way at the origin and
  ledge                         turns over a few per cent out, so the thermal
                                average samples the far side and C'(T) falls off
                                a cliff at the first finite temperature

This is the same failure as the compression escape and the nudge test, in the
third of the three directions a cubic crystal can be strained.  The objective
measured C11 and C12 at delta=0 and nothing anywhere else.

    python bain.py
    python bain.py V Ta Nb --max 0.2
"""
import io
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
os.environ.setdefault("LMP", f"{HOME}/lammps/src/lmp_serial")
os.environ.setdefault("BASEFILE",
                      os.path.join(ROOT, "standalone", "baselines.json"))

import numpy as np         # noqa: E402
import elastic_T as E      # noqa: E402

#  the shipped potentials are staged in the repo here, not in a LAMMPS install
if os.path.isdir(os.path.join(HERE, "potentials")):
    E.POTDIR = os.path.join(HERE, "potentials")

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       cell.data
pair_style      {style}
{coeff}
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

variable        d equal {d}
variable        sx equal 1.0+v_d
variable        sz equal 1.0/((1.0+v_d)*(1.0+v_d))
change_box      all x scale ${{sx}} y scale ${{sx}} z scale ${{sz}} remap units box
run             0
variable        e equal pe/atoms
print           "BAIN {d} ${{e}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def one(el, tag, d):
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    wd = os.path.join(HERE, "bain", f"{el}_{safe}_{d:+.4f}")
    try:
        style, pot = E.setup(wd, el, tag)
    except Exception as ex:
        return None
    open(os.path.join(wd, "in.bain"), "w").write(
        IN.format(style=style, coeff=E.coeff_line(style, pot, el), d=d))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(wd)} && {E.LMP} -in in.bain > out.txt 2>&1"],
                   capture_output=True, text=True)
    p = os.path.join(wd, "log.lammps")
    if not os.path.exists(p):
        return None
    m = re.search(r"BAIN\s+\S+\s+([-\d.eE+]+)",
                  io.open(p, errors="ignore").read())
    return float(m.group(1)) if m else None


def curve(el, tag, ds):
    nw = max(1, (os.cpu_count() or 4) - 2)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        return list(ex.map(lambda d: one(el, tag, d), ds))


def main():
    argv = sys.argv[1:]
    for opt in ("--max", "--sets"):
        if opt in argv:
            i = argv.index(opt)
            del argv[i + 1]
    args = [a for a in argv if not a.startswith("--")]
    #  the whole library by default, like every other sweep here; a four
    #  element default silently produced an eight-record file where a
    #  seventy-six record one was expected
    els = args or sorted(E.PACK)
    dmax = 0.15
    dump = "--dump" in sys.argv
    if "--max" in sys.argv:
        dmax = float(sys.argv[sys.argv.index("--max") + 1])
    ds = [round(x, 6) for x in np.linspace(-dmax, dmax, 31)]

    out = {}
    for el in els:
        sets = ["tap"]
        if "--sets" in sys.argv:
            sets = sys.argv[sys.argv.index("--sets") + 1].split(",")
        tags = [s for s in sets if s in E.PACK.get(el, {})]
        tags += [] if "--nobase" in sys.argv else ["base|" + f
                 for f, _ in E.BASE.get(el, [])
                 if all(os.path.exists(os.path.join(E.POTDIR, g))
                        for g in f.split("+"))]
        print(f"=== {el} ===")
        for tag in tags:
            e = curve(el, tag, ds)
            if any(v is None for v in e):
                print(f"  {tag}: the run did not finish")
                continue
            e = np.array(e)
            e0 = e[len(ds) // 2]
            #  C' from the curvature at the origin, in GPa.  For this strain
            #  E(d) = E0 + 6 V0 C' d^2 + O(d^3) per unit volume; the prefactor
            #  is not needed to answer the question, so the fit is reported as
            #  a curvature and the sign and the shape are what is read.
            c = np.polyfit(ds[13:18], (e - e0)[13:18], 2)[0]
            #  Where the well stops being a well, on each side separately.
            #  This is the number that decides the finite-temperature
            #  behaviour, and it is not the barrier height: a homogeneous
            #  tetragonal strain is a collective coordinate, so the barrier to
            #  crossing it scales with the number of atoms and the crystal
            #  never crosses.  What the fluctuation formula measures is the
            #  curvature averaged over the thermal strain distribution, whose
            #  width is sigma = sqrt(kT / (V_cell C')) - a few parts in a
            #  thousand.  If the well is narrower than a couple of sigma, the
            #  average curvature is negative and C11-C12 comes out negative
            #  while the crystal sits happily in its basin.
            mid = len(ds) // 2
            iu = next((i - 1 for i in range(mid + 1, len(ds))
                       if e[i] < e[i - 1]), None)
            idn = next((i + 1 for i in range(mid - 1, -1, -1)
                        if e[i] < e[i + 1]), None)
            up = ds[iu] if iu is not None else None
            dn = ds[idn] if idn is not None else None
            #  the barrier is the height of the LOCAL hump the well has to be
            #  climbed over, not the largest energy anywhere on that side.
            #  Reading the global maximum reported chromium's 0.9 meV ledge as
            #  a 364 meV barrier - off by two and a half orders of magnitude,
            #  and in the direction that makes a broken well look sound.
            bar_up = float(e[iu] - e0) if iu is not None else None
            bar_dn = float(e[idn] - e0) if idn is not None else None
            deep = float(e.min() - e0)
            lab = "bizimki (MAU)" if tag == "tap" else tag.replace("base|", "")
            fmt = (lambda x, b: f"{x:+.3f} ({b * 1000:.1f} meV)"
                   if x is not None else "none")
            print(f"  {lab:30s} egrilik {c:8.2f}"
                  f"   tepe(+) {fmt(up, bar_up):>20s}"
                  f"   tepe(-) {fmt(dn, bar_dn):>20s}"
                  f"   en dip {deep * 1000:+8.1f} meV")
            out[f"{el}|{tag}"] = {
                "curvature": float(c), "turn_up": up, "turn_dn": dn,
                "barrier_up": bar_up, "barrier_dn": bar_dn,
                "deepest": deep, "d": list(ds),
                "E": [float(x - e0) for x in e]}
            if dump:
                print("      d      E-E0 (meV/atom)")
                for d, v in zip(ds, e):
                    print(f"    {d:+6.3f}  {1000 * (v - e0):+10.2f}")
        print()

    #  merge, for the same reason struct_rank.py does
    fp = os.path.join(HERE, "bain.json")
    old = {}
    if os.path.exists(fp):
        try:
            old = json.load(open(fp))
        except Exception:
            old = {}
    old.update(out)
    tmp = fp + ".tmp"
    json.dump(old, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, fp)
    print(f"{len(out)} records written; file now holds {len(old)}  -> {fp}")


if __name__ == "__main__":
    main()
