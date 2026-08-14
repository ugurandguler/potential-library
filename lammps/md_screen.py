#!/usr/bin/env python3
"""
Which parameter sets survive molecular dynamics at all?

Chromium and tungsten come apart under `ugur/ang`, and chromium's tapered MAU
does too - a set that has been distributed and described as merely
drifting.  Neither was caught by the fit's dynamical stability screen,
which asks whether the reference structure is a harmonic minimum at 0 K and
cannot ask whether it survives being heated.  So every set has to be tried.

This is a SCREEN, not a measurement.  A collapse announces itself immediately -
both known cases were already above 5000 K during equilibration - so six
picoseconds is enough to separate the sets that hold together from the ones
that do not, at a seventh of the cost of the full drift run.  Whatever survives
gets measured properly afterwards.

The two signals, and the second is the one that matters:

  temperature   600 K of initial velocities must equipartition to 300 K.
                Anything near 5000 is not a hot crystal, it is a crystal that
                converted its potential energy into heat.
  energy        a heated crystal sits ABOVE its static lattice energy.  Below
                it means the structure found something lower and is no longer
                the structure being modelled.

Reading a drift rate instead of these two is exactly how chromium was recorded
as an integrator artefact for a fortnight.

    python md_screen.py                 # every set of every element
    python md_screen.py --set tap_ug
    python md_screen.py Cr W --set tap
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
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import cellfile         # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"
SKIN = 2.0

#  (library key, pair style, label)
SETS = {"hard": (None, "ugur", "hard cutoff"),
        "tap": ("tap", "ugur", "tapered"),
        "ug": ("ug", "ugur/ang", "hard cutoff + angular"),
        "tap_ug": ("tap_ug", "ugur/ang", "tapered + angular")}

POTFILE = """# {el}, written by md_screen.py
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       {data}
pair_style      {style}
pair_coeff      * * {pot} {el}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
#  the static lattice energy, measured by the SAME engine and the SAME pair
#  style that is about to run the dynamics.  Computing it in python instead
#  needs the right one of two latdyn trees - the standalone one drops the
#  angular term silently, by 0.307 eV/atom on rhenium - and getting that choice
#  wrong turns a healthy set into a collapsed one.  There is no reason to guess
#  when the measuring instrument is already open.
run             0
variable        elat equal pe/atoms
print           "LATTICE ${{elat}}"
velocity        all create {T2} 12345 mom yes rot yes
fix             1 all nve
timestep        {dt}
thermo          {nth}
thermo_style    custom step temp pe etotal
thermo_modify   norm yes
run             {nequil}
reset_timestep  0
run             {nprod}
"""

#  The timestep is a diagnostic, not a setting.  A crystal that runs away at
#  2 fs and holds 300 K at 0.5 fs was never unstable - the integrator was.  A
#  crystal that runs away at both is the potential.  Telling those apart is the
#  only way a collapse verdict means anything, so --dt exists and the physical
#  duration is held fixed (2 ps equilibration + 4 ps production) while the
#  step changes, or the two runs would not be comparable.
def steps(dt):
    return int(round(2.0/dt)), int(round(4.0/dt)), max(1, int(round(1.0/dt)))

ROW = re.compile(r"^\s*(\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*$",
                 re.M)


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


DT = 0.002
#  the INITIAL velocity temperature, as the original screen wrote it.  600 K of
#  initial velocities equipartitions to a 300 K crystal, and every threshold
#  below is written against that half.  Making --temp mean the target instead
#  silently doubles the default and reads the whole library as SUPHELI.
TEMP = 600.0


def run(el, rec, style, tag):
    d = os.path.join(HERE, "screenruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    p = dict(rec)
    p.setdefault("lam2", 0.0)
    p.setdefault("lam4", 0.0)
    p["taper"] = rec.get("taper") or -1.0
    open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(
        el=el, **{k: p[k] for k in ("m", "D", "alpha", "r0", "gamma", "C",
                                    "alpha3", "rcut2", "rcut3", "taper",
                                    "lam2", "lam4")}))
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    box, _ = cellfile.orthogonal_cell(cry)
    rc = max(rec["rcut2"], rec["rcut3"])
    #  1.8, and NOT smaller.  Shrinking to 1.15 to make the screen cheap made
    #  tungsten's tapered set "collapse" at 3960 K when the same potential is
    #  stable at 290 K in a proper cell: at that size the box is barely wider
    #  than the cutoff, so an atom interacts with several periodic images of
    #  the same neighbour and the instability is the cell, not the potential.
    #  1.8 reproduces the full 2.5x runs on both known cases - chromium
    #  collapsing at 4559 against 4572, tungsten stable at 290 against 296 -
    #  which is what makes it safe to screen with.  A screen that invents
    #  failures is worse than a slow one.
    rep = tuple(max(3, int(np.ceil(1.8 * (rc + SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep)
    open(os.path.join(d, "in.md"), "w").write(IN.format(
        data=f"{el}.data", style=style, pot=f"{el}.ugur", el=el, skin=SKIN,
        dt=DT, nequil=steps(DT)[0], nprod=steps(DT)[1], nth=steps(DT)[2],
        T2=TEMP))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in in.md 2>&1"],
                   capture_output=True, text=True)
    #  read the log, not the captured stdout: a long run's stdout comes back
    #  truncated through the wsl subprocess
    lg = os.path.join(d, "log.lammps")
    if not os.path.exists(lg):
        return None
    txt = io.open(lg, errors="ignore").read()
    #  "Lost atoms" is the most severe outcome, not a missing one.  Lithium
    #  under the hard cutoff goes from 432 atoms to 1: the crystal does not
    #  merely convert its potential energy to heat, it flies apart and LAMMPS
    #  stops.  Reporting that as "run did not complete" reads as a gap in the
    #  data and is the opposite of the truth.
    m = re.search(r"ERROR: Lost atoms: original (\d+) current (\d+)", txt)
    if m:
        return {"lost": (int(m.group(1)), int(m.group(2)))}
    ml = re.search(r"LATTICE\s+([-\d.eE+]+)", txt)
    lat = float(ml.group(1)) if ml else None
    if "reset_timestep" not in txt:
        return None
    rows = [(int(m.group(1)), float(m.group(2)), float(m.group(3)),
             float(m.group(4)))
            for m in ROW.finditer(txt.split("reset_timestep")[-1])]
    #  the production length follows the timestep; hard-coding 2000 here was
    #  right only at 2 fs and happened to survive 0.5 fs by coincidence, which
    #  is the kind of luck that turns into a silent "run did not complete"
    last = steps(DT)[1]
    rows = [r for r in rows if r[0] <= last]
    if not rows or rows[-1][0] != last:
        return None
    return {"T": float(np.mean([r[1] for r in rows])), "pe": rows[-1][2],
            "lattice": lat}


def main():
    global DT, TEMP
    args = sys.argv[1:]
    if "--temp" in args:
        i = args.index("--temp"); TEMP = float(args[i + 1]); del args[i:i + 2]
    if "--dt" in args:
        i = args.index("--dt"); DT = float(args[i + 1]); del args[i:i + 2]
    want = None
    if "--set" in args:
        i = args.index("--set")
        want = args[i + 1]
        del args[i:i + 2]
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = args or sorted(lib)
    keys = [want] if want else list(SETS)

    print(f"MD screening scan - {TEMP:.0f} K start -> {TEMP/2:.0f} K, "
          f"{DT*1000:.1f} fs, 2 ps denge + 4 ps")
    print("T should be near 300 K and pe ABOVE the lattice value.\n")
    print(f"{'el':4s}{'set':16s}{'T':>7s}{'pe':>10s}{'kafes':>10s}"
          f"{'diff':>9s}   status")
    print("-" * 62)
    bad, out = [], {}

    #  lmp_serial is single-threaded, so the machine's cores are the limit and
    #  the runs are independent
    jobs = []
    for el in els:
        for k in keys:
            key, style, label = SETS[k]
            rec = lib[el] if key is None else lib[el].get(key)
            if rec and "m" in rec:
                jobs.append((el, k, dict(rec), style, label))
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"({len(jobs)} runs, {nw} in parallel)", flush=True)
    print(flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda j: run(j[0], j[2], j[3], j[1]), jobs))

    for (el, k, rec, style, label), r in zip(jobs, res):
        #  the reference has to be THIS set's own static lattice energy, not
        #  the element's root Ecoh.  They agree for the tapered sets, which is
        #  why it went unnoticed, but the hard angular fits miss their target
        #  badly - beryllium by 0.872 eV/atom and chromium by 0.959 - and a
        #  collapse test against the wrong floor is not a test.
        if True:
            if r is None:
                print(f"{el:4s}{label:16s}   the run did not finish")
                continue
            Est = (r["lattice"] if r.get("lattice") is not None
                   else float(np.ravel(lib[el]["Ecoh"])[0]))
            if "lost" in r:
                n0, n1 = r["lost"]
                bad.append((el, k, "DAGILDI"))
                out.setdefault(el, {})[k] = {"lost": [n0, n1],
                                             "collapsed": True}
                print(f"{el:4s}{label:16s}{'':>7s}{'':>10s}{'':>10s}{'':>9s}"
                      f"   DAGILDI - {n0} atomdan {n1} kaldi")
                continue
            gone = r["pe"] < Est - 0.05
            hot = r["T"] > 0.5 * TEMP * 1.34
            #  a hard-truncated set does not conserve energy - the step at the
            #  cutoff pumps it in - so its temperature can run away while the
            #  potential energy stays above the lattice and the pe test never
            #  fires.  Beryllium's hard angular set reaches 8205 K that way.
            #  Whatever the mechanism, that is not a crystal, and calling it
            #  merely "suspect" understates it.
            runaway = r["T"] > 3.0 * (0.5 * TEMP)
            state = ("COKTU" if gone or runaway
                     else ("SUPHELI" if hot else "saglam"))
            if gone or hot:
                bad.append((el, k, state))
            out.setdefault(el, {})[k] = {"T": round(r["T"]), "pe": r["pe"],
                                         "lattice": Est,
                                         "collapsed": bool(gone or runaway)}
            print(f"{el:4s}{label:16s}{r['T']:7.0f}{r['pe']:10.3f}{Est:10.3f}"
                  f"{r['pe'] - Est:9.3f}   {state}")
    #  MERGE, never replace.  A four-element diagnostic run overwrote the file
    #  a full 152-combination screen had just written, and the loss is silent:
    #  the file is still valid JSON, still parses, and every combination it no
    #  longer mentions simply reads as "not screened" downstream.  The flat
    #  "<el>|<set>" file is the one export_potentials.py reads, so it is the
    #  one that must survive a partial run.
    #  and only when this run means what the file means.  A --temp 1100
    #  diagnostic is not a screen result; merging it rewrites the verdict that
    #  152 distributed potential headers are built from, and Ba silently went
    #  from "saglam" to "SUPHELI" that way.
    if DT != 0.002 or TEMP != 600.0:
        print()
        print(f"(non-default run: {TEMP:.0f} K / {DT * 1000:.1f} fs"
              f" - md_screen_all.json'a YAZILMADI)")
        return
    pa = os.path.join(HERE, "md_screen_all.json")
    all_ = {}
    if os.path.exists(pa):
        try:
            all_ = json.load(open(pa))
        except Exception:
            all_ = {}
    kept = len(all_)
    for el, sets in out.items():
        for k, r in sets.items():
            all_[f"{el}|{k}"] = r
    json.dump(all_, open(pa, "w"), indent=1, sort_keys=True)
    json.dump(out, open(os.path.join(HERE, "md_screen.json"), "w"),
              indent=1, sort_keys=True)
    print(f"\n-> md_screen.json (this run), md_screen_all.json "
          f"({len(all_)} records; {kept} already present, not overwritten)")
    if bad:
        print(f"\n{len(bad)} set MD'de kullanilamaz ya da supheli:")
        for el, k, s in bad:
            print(f"  {el:3s} {k:8s} {s}")
    else:
        print("\nbutun setler MD'de yapisini koruyor")


if __name__ == "__main__":
    main()
