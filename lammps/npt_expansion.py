#!/usr/bin/env python3
"""
Thermal expansion by NPT dynamics, with the published potentials beside it.

This replaces a quasi-harmonic calculation that could not be made to converge.
That attempt (`standalone/expansion.py`) minimises E_static + F_vib over volume
at each temperature, and the shift being looked for - two parts in a thousand
of the lattice parameter over four hundred kelvin - is small enough that the
answer moved with the fitting protocol: over a wide volume window sodium came
back contracting on heating, over a narrow one iron and tungsten did too, and
copper moved by thirteen per cent between them.  A first-order formulation
(fit E and F_vib separately, dV = -(dF_vib/dV)/(d2E/dV2)) is much better
conditioned and removed most of it, but by then the numbers had been through
four protocols and none of them could be checked, because latdyn only knows
this potential's own functional form.  There was no baseline.

That is the deciding argument rather than the numerics.  What made the surface
result mean anything was fifty-one published potentials going through the
identical code and landing where they should.  NPT dynamics can do the same
here: a barostat at zero pressure, the average cell edge against temperature,
no fitting, no mode tracking, and every EAM and MEAM potential in the library
runs through it unchanged.

    alpha = (1/a) da/dT, from a straight line through the measured a(T)

Two things are checked rather than assumed.  The average pressure must come
back at zero, or the barostat has not equilibrated and the volume is not the
equilibrium one.  And the average temperature must match the target, which is
the same guard that caught a collapsed crystal being reported as a drift rate
in the NVE test.

    python npt_expansion.py Cu Al --sets tap
    python npt_expansion.py            # everything with an experimental value
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

LOCAL = not os.path.exists("/arf")
HOME = (subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                       capture_output=True, text=True).stdout.strip()
        if LOCAL else os.path.expanduser("~"))
os.environ.setdefault("LMP", f"{HOME}/lammps/src/lmp_serial")
os.environ.setdefault("BASEFILE",
                      os.path.join(ROOT, "standalone", "baselines.json"))

import numpy as np         # noqa: E402
import latdyn as L         # noqa: E402
import cellfile           # noqa: E402
import refdata            # noqa: E402
import elastic_T as E     # noqa: E402

for _cand in (os.path.join(HERE, "potentials"),
              os.path.expanduser("~/lammps/potentials")):
    if os.path.isdir(_cand):
        E.POTDIR = _cand
        break

#  Below the Debye temperature the expansion is not linear and a straight line
#  through it means nothing, so the grid starts at 200 K.  It stops at 700 K or
#  0.6 of the melting point, whichever is lower, to stay clear of the region
#  where a small perfect crystal starts to matter.
def temps(el):
    tm = refdata.MELTING.get(el, 2000.0)
    hi = min(700.0, 0.6 * tm)
    if hi <= 250.0:                     # the alkalis melt too low for that
        hi = max(0.55 * tm, 120.0)
    lo = min(200.0, 0.25 * hi)
    return [round(lo + (hi - lo) * k / 4.0, 1) for k in range(5)]


IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       cell.data
pair_style      {style}
{coeff}
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

velocity        all create {T2} 12345 mom yes rot yes
#  isotropic barostat: the cell is cubic or hexagonal and stays so, and an
#  anisotropic one would let it drift into a shape the reference is not
fix             1 all npt temp {T} {T} 0.1 iso 0.0 0.0 1.0
timestep        0.001
thermo          500
thermo_style    custom step temp press lx ly lz

run             {neq}
reset_timestep  0

variable        vlx equal lx/{nx}
variable        vly equal ly/{ny}
variable        vlz equal lz/{nz}
variable        vt  equal temp
variable        vp  equal press
fix             av all ave/time 10 {nrep} {nfreq} v_vlx v_vly v_vlz v_vt v_vp &
                ave running
run             {nrun}

variable        ax equal f_av[1]
variable        ay equal f_av[2]
variable        az equal f_av[3]
variable        tt equal f_av[4]
variable        pp equal f_av[5]
print           "RES ax ${{ax}}"
print           "RES ay ${{ay}}"
print           "RES az ${{az}}"
print           "RES T ${{tt}}"
print           "RES P ${{pp}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def run(d, text, name):
    open(os.path.join(d, name), "w").write(text)
    p = os.path.join(d, "log.lammps")
    if os.path.exists(p):
        os.remove(p)
    #  DEVNULL, not capture_output - see the note in elastic_T.run.  This job
    #  survived the deadlock only because a half-hour run leaves the window
    #  open too rarely to be caught in it.
    quiet = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if LOCAL:
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {E.LMP} -in {name} > out.txt 2>&1"],
                       **quiet)
    else:
        subprocess.run(f"cd {d} && {E.LMP} -in {name} > out.txt 2>&1",
                       shell=True, **quiet)
    return io.open(p, errors="ignore").read() if os.path.exists(p) else ""


def one(job):
    el, tag, T = job
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "nptruns", f"{el}_{safe}_{int(T)}K")
    os.makedirs(d, exist_ok=True)
    try:
        style, pot = E.setup(d, el, tag)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    rep = tuple(E.PACK[el]["rep"])
    txt = IN.format(style=style, coeff=E.coeff_line(style, pot, el),
                    T=float(T), T2=2.0 * float(T),
                    nx=rep[0], ny=rep[1], nz=rep[2],
                    #  10 ps to settle the barostat and 30 ps to average
                    #  the cell edge.  Ninety was generous for a quantity this
                    #  well behaved and put a single run at two hours for the
                    #  three-body form, which is what made the local machine
                    #  the wrong place for this.
                    #
                    #  nrep is DERIVED, not chosen.  fix ave/time requires
                    #  nevery*(nrepeat-1) <= nfreq, and halving nfreq from 5000
                    #  to 2500 while leaving nrepeat at 500 violated it -
                    #  10*499 = 4990 against 2500.  Every one of six hundred
                    #  and thirty-five runs on the cluster died on that line
                    #  after the local test had passed with the longer window.
                    #  Tying the two together makes the constraint impossible
                    #  to break by editing one of them.
                    neq=10000, nrun=30000,
                    nrep=(2500 // 10), nfreq=2500)
    lg = run(d, txt, "in.npt")
    r = {m.group(1): float(m.group(2))
         for m in re.finditer(r"RES\s+(\w+)\s+([-\d.eE+]+)", lg)}
    if len(r) < 5:
        err = [l for l in lg.splitlines() if "ERROR" in l]
        return {"error": err[0][:70] if err else "the run did not finish"}
    return {"ax": r["ax"], "ay": r["ay"], "az": r["az"],
            "Tavg": r["T"], "Pavg": r["P"],
            #  the two guards: a barostat that has not settled leaves a
            #  pressure behind, and a thermostat that has not is not at T
            "P_ok": bool(abs(r["P"]) < 500.0),
            "T_ok": bool(abs(r["T"] - T) < 0.15 * T)}


def main():
    argv = sys.argv[1:]
    sets = ["tap"]
    if "--sets" in argv:
        i = argv.index("--sets")
        sets = argv[i + 1].split(",")
        del argv[i:i + 2]
    with_base = "--nobase" not in argv
    els = [a for a in argv if not a.startswith("--")]

    sys.path.insert(0, os.path.join(ROOT, "standalone"))
    import expansion as X          # the experimental table lives there
    els = els or sorted(e for e in E.PACK if e in X.ALPHA_EXP)

    jobs = []
    for el in els:
        for tag in sets:
            if tag in E.PACK.get(el, {}):
                jobs += [(el, tag, T) for T in temps(el)]
        if with_base:
            for fn, _ in E.BASE.get(el, []):
                if all(os.path.exists(os.path.join(E.POTDIR, g))
                       for g in fn.split("+")):
                    jobs += [(el, "base|" + fn, T) for T in temps(el)]
    if not jobs:
        print("0 runs - element or set name not recognised")
        return
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"{len(jobs)} runs, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(one, jobs))

    got = {}
    for (el, tag, T), r in zip(jobs, res):
        got.setdefault(f"{el}|{tag}", {})[str(T)] = r

    out = {}
    print()
    print(f"{'el':4s}{'source':>28s}{'alpha':>9s}{'expt':>8s}{'ratio':>7s}"
          f"{'P ort':>9s}   not")
    print("-" * 72)
    for key, series in sorted(got.items()):
        el, tag = key.split("|", 1)
        pts = [(float(T), r) for T, r in series.items()
               if "error" not in r]
        bad = [T for T, r in pts if not (r["P_ok"] and r["T_ok"])]
        pts = [(T, r) for T, r in pts if r["P_ok"] and r["T_ok"]]
        if len(pts) < 4:
            #  say WHY, not just how many.  Reporting "0 of 5 points" hid a
            #  LAMMPS error message that was sitting in every log file and
            #  named the problem exactly.
            errs = [r["error"] for r in series.values() if "error" in r]
            why = errs[0] if errs else (
                "korumalar eledi" if bad else "sebep bilinmiyor")
            print(f"{el:4s}{tag[:28]:>28s}   yetersiz nokta "
                  f"({len(pts)}/{len(series)}): {why[:46]}")
            #  Record the failure instead of writing nothing.  An absent key
            #  is indistinguishable from a key that was never asked for, and
            #  downstream that let lithium keep a stale quasi-harmonic entry
            #  and show it on the page as a NaN.  A run that died should leave
            #  a record saying it died.
            out[key] = {"failed": why[:120], "points": len(pts),
                        "of": len(series),
                        "alpha_exp_1e6": X.ALPHA_EXP.get(el)}
            continue
        pts.sort()
        Ts = np.array([p[0] for p in pts])
        #  the cube root of the cell volume per unit cell: works for cubic and
        #  hexagonal alike, and is what "the lattice expands" means
        aa = np.array([(p[1]["ax"] * p[1]["ay"] * p[1]["az"]) ** (1 / 3.0)
                       for p in pts])
        slope = float(np.polyfit(Ts, aa, 1)[0])
        a300 = float(np.interp(300.0, Ts, aa))
        alpha = slope / a300 * 1e6
        exp = X.ALPHA_EXP.get(el)
        lab = {"tap": "MAU", "tap_ug": "UG"}.get(tag, tag.replace("base|", ""))
        out[key] = {"T": list(Ts), "a": list(aa), "alpha_1e6": alpha,
                    "alpha_exp_1e6": exp,
                    "ratio": (alpha / exp) if exp else None,
                    "Pavg": [p[1]["Pavg"] for p in pts],
                    "dropped": bad}
        print(f"{el:4s}{lab[:28]:>28s}{alpha:9.1f}{(exp or 0):8.1f}"
              f"{(alpha / exp if exp else 0):7.2f}"
              f"{np.mean([p[1]['Pavg'] for p in pts]):9.1f}"
              f"   {('%d nokta atildi' % len(bad)) if bad else ''}")

    fp = os.path.join(HERE, "npt_expansion.json")
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
    print(f"\n{len(out)} records written; file now holds {len(old)}  -> {fp}")


if __name__ == "__main__":
    main()
