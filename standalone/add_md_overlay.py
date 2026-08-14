#!/usr/bin/env python3
"""
Attach what the potential does in molecular dynamics to library.json.

The page already warns that the shipped parameters cannot be used for dynamics.
That warning was argued from the static curve; this puts the measurement behind
it.  Two things per element:

  step        phi2 at rcut2.  Not zero, which is the whole problem: a neighbour
              crossing the sphere changes the energy by this much in one step.
              Computed here for all 38.

  drift       energy drift in an NVE run, meV/atom/ps, hard cutoff against the
              switch.  Measured in lammps/nve_check.py for eight elements
              chosen to span the range of `step`, since running all 38 would
              add nothing - the point is the correlation, not the census.

Over those eight, log(step) against log(drift) correlates at 0.944.  That is
what makes the cutoff the cause rather than a coincidence: palladium's step is
0.0002 eV and it drifts by 5, chromium's is 0.129 and it drifts by 7256.

    python add_md_overlay.py
"""
import json
import os

import numpy as np

import latdyn as L

HERE = os.path.dirname(os.path.abspath(__file__))

#  lammps/nve_check.py, 600 K start, 2 fs, 40 ps NVE, meV/atom/ps.
#  Chromium's tapered figure is integrator error and not the potential: it
#  falls to -1.37 at 0.5 fs, against -1.54 predicted by dt^2, and chromium is
#  the stiffest element here.  Anything under about ten is in that regime.
#  The four hexagonal metals were added last and they are the strongest four
#  results in the table.  With the hard cutoff three of them do not merely
#  drift, they destroy the crystal - the mean potential energy ends below the
#  static lattice, which is the signature of a structure that has come apart -
#  and ruthenium, which has the deepest phi2 at its cutoff of any element in
#  the library at -0.346 eV, drifts by 3376 meV/atom/ps against 0.065 with the
#  switch.  A ratio of fifty-two thousand.
NVE = {
    "W":  (581.647, 0.254),
    "Fe": (6343.944, 0.278),
    "Cr": (7256.342, -24.612),
    "Cu": (161.796, 0.131),
    "Al": (724.390, 0.232),
    "Nb": (1643.222, 0.038),
    "Pd": (4.967, 0.065),
    "Au": (70.697, 0.035),
    "Ru": (3376.026, 0.065),
    "Re": (36.552, -0.002),
    "Y":  (41.150, 0.005),
    "Ti": (-142.283, -0.026),
}
#  which of them lost the crystal outright under the hard cutoff, rather than
#  merely drifting; reported separately because a drift rate understates it
NVE_BROKE = ("Ru", "Y", "Ti")
#  log-log correlation of the step at the cutoff against the hard-cutoff
#  drift.  0.944 was measured over the first eight, which happened to be all
#  cubic.  With the four hexagonal metals it is 0.653 over all twelve and 0.787
#  over the nine whose crystal survived - and only the second of those is a
#  correlation between two drift rates, because for a run that came apart the
#  slope is not a drift.  The honest figure to quote is 0.79 over nine.
#  Rhenium is the outlier: its step is comparable to niobium's and it drifts
#  forty-five times less.
NVE_R = 0.787


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    n = 0
    print(f"{'el':4s}{'phi2(rcut2)':>13s}{'NVE hard':>12s}{'NVE taper':>12s}")
    print("-" * 41)
    for el, v in sorted(lib.items()):
        pot = L.Potential.from_record(v)
        step = float(pot.phi2(np.array([v["rcut2"]]), 0)[0])
        #  Merge into whatever is already there.  `v["md"] = md` replaced it,
        #  and md_screen.py's T, collapsed and drift live in the same dict, so
        #  running this after the screen silently emptied three fields on every
        #  element.  That is the same accident that once wiped the page's
        #  cutoff narrative, and it leaves no error behind - the page simply
        #  loses a section.
        md = v.get("md") if isinstance(v.get("md"), dict) else {}
        md["step_eV"] = round(step, 6)
        if el in NVE:
            md["drift_hard"] = NVE[el][0]
            md["drift_taper"] = NVE[el][1]
            md["r"] = NVE_R
            md["hard_broke"] = el in NVE_BROKE
        v["md"] = md
        n += 1
        if el in NVE:
            print(f"{el:4s}{step:13.4f}{NVE[el][0]:12.1f}{NVE[el][1]:12.3f}")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\ncutoff step folded into {n} elements, NVE drift into {len(NVE)}")
    print("birlestirildi:", path)


if __name__ == "__main__":
    main()
