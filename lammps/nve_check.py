#!/usr/bin/env python3
"""
Does energy stay conserved in NVE?  Hard cutoff against the switch.

The whole case for the taper is that a truncated potential cannot be used for
dynamics: phi2 does not vanish at rcut2, so a neighbour crossing the sphere
changes the energy in one step.  So far that has been argued from the static
curve - beryllium gains a shell under 0.8 per cent compression and drops
0.16 eV/atom.  This measures it where it actually matters.

The test is deliberately plain.  A perfect crystal, heated to a temperature
where atoms move enough to cross the cutoff, then NVE with no thermostat.  In
NVE the total energy is a constant of the motion; anything else is the
integrator meeting a discontinuity.  Two runs per element, same seed, same
timestep, same everything, and the only difference is the taper flag in the
potential file.

What is reported is the drift of the total energy in meV/atom/ps, which is the
number an MD practitioner judges a potential by.  Below about 1 is fine, tens
is unusable, and the scale is set by the size of the step the potential has at
its cutoff, not by the timestep.

Two things to hold on to when reading the numbers.

The drift tracks the step, quantitatively: over the first eight elements, all
cubic, the correlation between log |phi2(rcut2)| and log(drift with the hard
cutoff) is 0.944.  Palladium's step is 0.0002 eV and it drifts by 5; chromium's
is 0.129 and it drifts by 7256.  That is what makes this a measurement of the
cutoff rather than an observation that one potential happened to behave better.

Adding the four hexagonal metals weakens that number and strengthens the
finding.  It falls to 0.653 over all twelve and 0.787 over the nine whose
crystal survived, with rhenium the outlier - its step is comparable to
niobium's and it drifts forty-five times less.  But three of the four hcp runs
with the hard cutoff do not drift at all in any useful sense: the crystal comes
apart, which the KRISTAL DAGILDI column reports separately for exactly this
reason.  A categorical failure is worse than a large rate, and averaging the
two into one regression understates it.

Below about ten the number is integrator error and not the potential, so do not
read it as one.  Chromium is the stiffest element in the set - Debye 630 K -
and 2 fs is too long for it: with the taper on it drifts -24.6 at 2 fs, +10.6
at 1 fs and -1.37 at 0.5 fs, against -1.54 predicted by dt^2 from the first.
The sign flip in the middle is the giveaway that the measurement is in the
noise there.  Use a shorter timestep for a stiff metal, as one would anyway.

    python nve_check.py            # a few representative elements
    python nve_check.py Cu W Be
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
SKIN = 2.0

#  Ruthenium has the deepest phi2 at its cutoff (-0.346 eV) and is hcp, which is
#  why it heads the list now that the hcp cells are built from latdyn's own
#  Crystal rather than LAMMPS's ideal-c/a `lattice hcp`.  The cubic ones span
#  the range, from tungsten and iron with a large step at the cutoff to
#  palladium and gold with almost none.
DEFAULT = ["Ru", "Re", "Y", "Ti", "W", "Fe", "Cr", "Cu", "Al", "Nb", "Pd", "Au"]

POTFILE = """# {el}, written by nve_check.py from library.json
# Tersoff-style: one line per (centre, leg, leg) triple.  For one element that
# is a single line, and it is written this way so the single-species files and
# the alloy files are the same format - a regression here is a regression in
# what the alloy path will use.
# el1 el2 el3  m D alpha r0 gamma C alpha3 rcut2 rcut3 taper lam2 lam4
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} 0 0
"""

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       {data}
pair_style      ugur
pair_coeff      * * {pot} {el}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes

velocity        all create {T} 12345 mom yes rot yes
fix             1 all nve
timestep        {dt}
thermo          {every}
thermo_style    custom step temp pe etotal
thermo_modify   norm yes

#  equilibrate off the perfect lattice first, so the initial relaxation of the
#  potential energy does not masquerade as non-conservation
run             {neq}
reset_timestep  0
run             {nrun}
"""


#  the thermo table, whose last column is the total energy per atom
ROW = re.compile(r"^\s*(\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*$",
                 re.M)


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def replication(cry, rcut):
    """box copies enough to clear the interaction range on every axis

    hcp needs this more than cubic does: the orthogonal cell is a by a*sqrt3
    by c, so a flat 4x4x4 would be half as thick along x as along y.
    """
    box, _ = cellfile.orthogonal_cell(cry)
    return tuple(max(3, int(np.ceil(2.5 * (rcut + SKIN) / b))) for b in box)


def run(el, rec, e, taper, T=600.0, dt=0.002, neq=2000, nrun=20000,
        every=500):
    tag = "taper" if taper else "hard"
    d = os.path.join(HERE, "nveruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(
        el=el, taper=(taper or -1.0),
        **{k: rec[k] for k in ("m", "D", "alpha", "r0", "gamma",
                               "C", "alpha3", "rcut2", "rcut3")}))
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    rep = replication(cry, max(rec["rcut2"], rec["rcut3"]))
    nat, _ = cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                                 refdata.MASSES[el], rep)
    open(os.path.join(d, "in.nve"), "w").write(
        IN.format(data=f"{el}.data", skin=SKIN, pot=f"{el}.ugur", el=el,
                  T=T, dt=dt, neq=neq, nrun=nrun, every=every))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in in.nve > out.txt 2>&1"],
                   capture_output=True, text=True)

    #  Read the LOG, not the screen.  LAMMPS echoes the input commands to
    #  log.lammps and not to the screen, so the "reset_timestep" guard below
    #  never matched when it was fed stdout: four completed twenty-thousand
    #  step runs were all reported as "did not reach reset_timestep", which is
    #  the failure message for a killed run.  The runs were fine; the reader
    #  was looking at the wrong stream.
    class _R:
        pass
    r = _R()
    lg = os.path.join(d, "log.lammps")
    r.stdout = (io.open(lg, errors="ignore").read()
                if os.path.exists(lg) else "")
    #  Only the second leg: everything after the reset.  The guard matters
    #  more than it looks.  If the run has not reached `reset_timestep` -
    #  killed, still going, out of time - the string is absent and split()
    #  hands back the WHOLE log, so the equilibration rows get read as
    #  production.  Those start at 600 K and fall towards 300, which fits a
    #  large positive slope and an average temperature near 355 K.  That is
    #  exactly how titanium and yttrium were first recorded as drifting by
    #  17.28 and 5.73 meV/atom/ps when they in fact conserve to 0.002; the
    #  giveaway was the temperature, since 600 K of initial velocities must
    #  equipartition to 300.
    if "reset_timestep" not in r.stdout:
        return None, None, "the run never reached reset_timestep (half a day)", None
    tail = r.stdout.split("reset_timestep")[-1]
    #  columns are step, temp, poteng, toteng - in that order, and swapping
    #  two of them here once cost an afternoon: the potential energy fluctuates
    #  by tens of meV while the total is flat, so fitting the wrong column
    #  reports a large drift for a run that is conserving perfectly
    rows = [(int(m.group(1)), float(m.group(2)),
             float(m.group(3)), float(m.group(4)))
            for m in ROW.finditer(tail)]
    rows = [x for x in rows if x[0] <= nrun]
    #  and it has to have finished, not merely started
    if rows and rows[-1][0] != nrun:
        return None, None, f"the run was cut off at step {rows[-1][0]}/{nrun}", None
    if len(rows) < 3:
        err = [l for l in r.stdout.splitlines() if "ERROR" in l]
        return None, None, (err[0][:90] if err else "no thermo line"), None
    #  least-squares slope of total energy against time, meV/atom/ps
    import numpy as np
    ps = np.array([x[0] * dt / 1000.0 for x in rows])
    et = np.array([x[3] * 1000.0 for x in rows])      # toteng, meV/atom
    slope = float(np.polyfit(ps, et, 1)[0])
    Tavg = float(np.mean([x[1] for x in rows]))
    #  A heated crystal sits ABOVE its static lattice energy.  When the
    #  potential energy ends up below it, the run has not drifted - the
    #  structure has collapsed, because the step at the cutoff acts as an
    #  energy source and the system finds something that is not the crystal.
    #  Reporting that as a drift rate understates it, so it is returned
    #  separately and the caller says so.
    return slope, Tavg, None, rows[-1][2]


def main():
    import numpy as np
    import latdyn as L
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or DEFAULT
    print("600 K start, 2 fs, 40 ps NVE - total energy drift")
    print(f"{'el':4s}{'phi2(rc2) eV':>14s}{'SERT':>12s}{'GECISLI':>12s}"
          f"{'ratio':>9s}{'T avg K':>9s}")
    print("-" * 60)
    #  The runs are independent and there are two per element, so doing them
    #  one after another wasted the machine: a single 40 ps run takes half an
    #  hour and twenty-four of them in series is thirteen hours.  They are
    #  dispatched together and the table is assembled afterwards, so the
    #  output is unchanged.
    from concurrent.futures import ThreadPoolExecutor
    jobs, meta = [], []
    for el in els:
        e = refdata.ELEMENTS[el]
        hard = lib[el]
        if not hard.get("tap"):
            print(f"{el:4s}  (no tapered fit)")
            continue
        tap = dict(hard["tap"])
        pot = L.Potential.from_record(hard)
        step = float(pot.phi2(np.array([hard["rcut2"]]), 0)[0])
        meta.append((el, e, hard, tap, step))
        jobs.append((el, hard, e, None))
        jobs.append((el, tap, e, tap.get("taper") or 0.85))
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"({len(jobs)} runs, {nw} in parallel)", flush=True)
    print(flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda j: run(*j), jobs))
    for i, (el, e, hard, tap, step) in enumerate(meta):
        dh, Th, eh, peh = res[2 * i]
        dt_, Tt, et, pet = res[2 * i + 1]
        if dh is None or dt_ is None:
            print(f"{el:4s}  error: {eh or et}")
            continue
        fac = (abs(dh) / abs(dt_)) if dt_ else float("inf")
        Estat = float(np.ravel(hard["Ecoh"])[0])
        gone = [t for t, pe in (("hard", peh), ("tapered", pet))
                if pe is not None and pe < Estat - 0.05]
        print(f"{el:4s}{step:14.4f}{dh:12.3f}{dt_:12.3f}{fac:8.0f}x"
              f"{Th:9.0f}   {('KRISTAL DAGILDI: ' + ', '.join(gone)) if gone else ''}")
    print("\nbirim: meV/atom/ps.  1'in alti iyi, onlarca kullanilamaz.")


if __name__ == "__main__":
    main()
