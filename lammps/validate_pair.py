#!/usr/bin/env python3
"""
Prove `pair_style ugur` inside LAMMPS against standalone/latdyn.py.

validate_kernel.py checks the physics per interaction.  This checks the thing
that stands between that kernel and a simulation: neighbour lists, the triplet
loop, the half/full-list bookkeeping, the virial.  A kernel that is right to
1e-16 and a pair style that double-counts a bond give the same-looking
potential and completely different numbers, so both have to be measured.

Three quantities, each catching a different mistake:

  cohesive energy     a factor of two in the pair sum, or triplets counted
                      twice, shows here and nowhere else
  pressure at a0      the fit pins this to zero; if LAMMPS disagrees the
                      derivative is wrong even when the energy is right.
                      Mind the units: `metal` reports pressure in BAR, and
                      latdyn works in GPa with the opposite sign convention,
                      so the raw numbers differ by -10000 and that is not a
                      bug.  Reading it as one cost an hour.
  elastic constants   the whole point, and the only one that exercises the
                      three-body derivatives properly

Writes the potential file LAMMPS needs from library.json, so the two sides
cannot drift: same parameters, same cutoffs, same taper flag.  The geometry
comes from latdyn's own Crystal through cellfile.py rather than from LAMMPS's
`lattice` command, for the same reason: `lattice hcp` fixes c/a at the ideal
sqrt(8/3) and every hcp element in the library sits below it, beryllium by four
per cent, so the two codes would have been given different crystals.  The
twenty-six cubic elements were validated through the `lattice` path first and
give identical numbers through this one, which is what makes it safe to trust
for the twelve where there is nothing to compare against.

    python validate_pair.py            # every element
    python validate_pair.py Cu Fe Nb
"""
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

#  Written by this script, read by pair_ugur::read_file.  Ten numbers in the
#  order the fit stores them; `taper` is not optional and a negative value
#  means hard truncation, because a potential fitted with the switch on and run
#  with it off is a different potential.
POTFILE = """# {el}, written by validate_pair.py from library.json
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
neigh_modify    delay 0 every 1 check no
run             0
variable        epa equal pe/atoms
variable        prs equal press
print           "UGUR_E ${{epa}}"
print           "UGUR_P ${{prs}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def replication(cry, rcut):
    """enough copies that the box comfortably clears the interaction range

    LAMMPS builds ghosts out to the cutoff whatever the box, so this is not
    about correctness so much as staying well away from the regime where it
    has to.  The old `lattice` path used a flat 4x4x4, which is generous for
    copper and thin for caesium; sizing it from the cutoff treats the two the
    same way, and hcp needs it anyway because the three box lengths differ.
    """
    box, _ = cellfile.orthogonal_cell(cry)
    return tuple(max(3, int(np.ceil(2.5 * (rcut + SKIN) / b))) for b in box)


def lammps(el, rec, e, cry=None):
    d = os.path.join(HERE, "pairruns", el)
    os.makedirs(d, exist_ok=True)
    pf = os.path.join(d, f"{el}.ugur")
    open(pf, "w").write(POTFILE.format(
        el=el, taper=(rec.get("taper") or -1.0),
        **{k: rec[k] for k in ("m", "D", "alpha", "r0", "gamma",
                               "C", "alpha3", "rcut2", "rcut3")}))
    if cry is None:
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
    rep = replication(cry, max(rec["rcut2"], rec["rcut3"]))
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep)
    open(os.path.join(d, "in.check"), "w").write(IN.format(
        data=f"{el}.data", skin=SKIN, pot=f"{el}.ugur", el=el))
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {LMP} -in in.check 2>&1"],
                       capture_output=True, text=True)
    out = r.stdout
    g = {}
    for k, pat in (("E", r"UGUR_E\s+([-\d.eE+]+)"),
                   ("P", r"UGUR_P\s+([-\d.eE+]+)")):
        m = re.search(pat, out)
        if m:
            g[k] = float(m.group(1))
    if not g:
        g["err"] = "\n".join(l for l in out.splitlines()
                             if "ERROR" in l or "error" in l)[:200]
    return g


#  metal units report pressure in bar; latdyn returns dE/d(vol strain)/(3V) in
#  GPa, which is minus the thermodynamic pressure.  Both conversions here.
BAR_PER_GPA = 1.0e4


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or sorted(lib)
    print(f"{'el':4s}{'E latdyn':>12s}{'E lammps':>12s}{'dE':>11s}"
          f"{'P latdyn':>11s}{'P lammps':>11s}  status")
    print("-" * 70)
    bad = []
    for el in els:
        e = refdata.ELEMENTS[el]
        rec = lib[el]
        pot = L.Potential.from_record(rec)
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        E_ref = L.energy(cry, pot)
        g = lammps(el, rec, e, cry)
        if "E" not in g:
            print(f"{el:4s}  LAMMPS error: {g.get('err','?')[:60]}")
            bad.append(el)
            continue
        dE = g["E"] - E_ref
        #  latdyn's own pressure at the same geometry, for the virial check.
        #  h = 1e-5, NOT the 2e-4 the fit uses: at the fitting step the
        #  constraint has pinned the pressure to zero, so comparing against it
        #  measures the constraint rather than the virial.  Converged, latdyn
        #  and LAMMPS agree to three figures - platinum 1.63e-4 against
        #  1.64e-4.  Larger steps are worse, not better: at h = 1e-3 tungsten
        #  reads 876 GPa, because the strain carries neighbours across the hard
        #  cutoff and the discontinuity lands in the pressure.
        h = 1e-5
        nat = len(cry.frac)
        ep = L.energy(cry.strained(np.eye(3) * h), pot)
        em = L.energy(cry.strained(np.eye(3) * (-h)), pot)
        P_ref = ((ep - em) / (2 * h) / (3 * cry.vol / nat)) * L.EV_A3_TO_GPA
        P_lmp = -g["P"] / BAR_PER_GPA
        dP = P_lmp - P_ref
        ok = (abs(dE) < 1e-6 * max(abs(E_ref), 1.0)
              and abs(dP) < 1e-6 + 1e-3 * abs(P_ref))
        print(f"{el:4s}{E_ref:12.6f}{g['E']:12.6f}{dE:11.2e}"
              f"{P_ref:11.2e}{P_lmp:11.2e}  {'ok' if ok else 'DIFFERS'}")
        if not ok:
            bad.append(el)
    print()
    if bad:
        print("uyusmayan:", " ".join(bad))
        raise SystemExit(1)
    print("pair_style ugur gives latdyn.py's energy AND pressure")


if __name__ == "__main__":
    main()
