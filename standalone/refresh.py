#!/usr/bin/env python3
"""
Rebuild everything downstream of the fits, in the one order that works.

    python refresh.py                       # rebuild from the current fit.json
    python refresh.py runs/2026-08-01_truba # merge that run in first

Everything the viewer shows is derived from fit.json, so a new search only has
to be merged and then this run: elastic constants, the mechanical analysis,
phonons, thermodynamics and the Materials Project overlay are all recomputed
from the merged parameters.  Nothing is carried over except the fetched MP
reference data, which is a property of the element and not of the fit.

The order is not arbitrary:

  merge_fits      stability-aware best-of; must run before anything reads
                  fit.json, and re-measures stability itself rather than
                  trusting whatever screen a run happened to use
  build_library   elastic tensor, mechanics.py analysis, dispersion, thermo
  add_dynstab     stability flags on the 8^3 U 9^3 union mesh
  add_mp_overlay  our dispersion re-evaluated at MP's own q-points
  add_mc3d_overlay  the same against Materials Cloud MC3D, which covers
                  Ag Au Cr Mo Nb Ni Pb Pd Rh Ta where MP has nothing
  add_exp_phonon  measured frequencies at X/L (fcc) and H/N (bcc), scored
  add_jarvis_overlay  a third DFT opinion, NIST JARVIS-DFT.  Entries more than
                  25 % away from a reference we already hold are stored but
                  flagged rather than drawn.  **Must come after add_exp_phonon**,
                  which supplies the only reference several elements have: with
                  the two the other way round, iron had nothing to be checked
                  against and its non-spin-polarised curve - 121 cm^-1 against a
                  measured 309 - was drawn as trusted
  add_reachability  C44/C' floor verdict; reads R_floor.json and
                  cprime_region.json, which change only if the form does
  fix_mp_path     MP labels and discontinuities - AFTER add_mp_overlay, which
                  rewrites the phonon record
  make_gui        the page

A step that fails stops the run, because a later step reading half-written data
is worse than no rebuild.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

#  Third field: True when the step attaches THIRD-PARTY comparison data that a
#  fresh checkout does not have.  Those datasets are fetched by the fetch_*.py
#  scripts and are not redistributed here, so on a clean clone they are absent -
#  which is not a failure, it is the expected state.  Treating it as one made
#  this script stop at step 3 of 10 for anyone who cloned the repository and
#  followed the README.
#
#  The core steps are NOT optional: they compute from the potential itself, and
#  if one of them fails the page must not be built.
CHAIN = [
    ("build_library.py", "elastic tensor, mechanics, phonons, thermodynamics", False),
    ("add_dynstab.py", "dynamical stability on the union mesh", False),
    ("add_mp_overlay.py", "Materials Project comparison at their q-points", True),
    ("add_mc3d_overlay.py", "Materials Cloud MC3D phonon comparison", True),
    ("add_exp_phonon.py", "measured phonon frequencies, out-of-sample", True),
    ("add_jarvis_overlay.py", "JARVIS-DFT comparison, conflicting entries flagged", True),
    ("add_ug_overlay.py", "UG (angular) results beside MAU, if any exist", True),
    ("add_reachability.py", "is this metal inside the form's reach", False),
    ("fix_mp_path.py", "MP high-symmetry labels and path discontinuities", True),
    ("make_gui.py", "potential.html", False),
]


def run(script, *args, optional=False):
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *args],
                       capture_output=True, text=True, cwd=HERE)
    tail = [l for l in p.stdout.strip().splitlines()
            if l.strip() and "Warning" not in l][-1:]
    print(f"    {time.time() - t0:6.1f}s  {tail[0].strip() if tail else 'done'}")
    if not p.returncode:
        return
    #  What may be skipped is "this step's input is not in this checkout" -
    #  either because the script guarded it and stopped itself (SystemExit) or
    #  because it opened a file that is not there (FileNotFoundError).  Some of
    #  these scripts do the first and some the second, and both mean the same
    #  thing.  Anything ELSE that raises is a bug and must still stop the
    #  chain, or a broken overlay would be forgiven by the mechanism that
    #  forgives a missing dataset.
    tb = "Traceback (most recent call last)" in p.stderr
    missing = "FileNotFoundError" in p.stderr
    if optional and (not tb or missing):
        why = (p.stderr.strip() or p.stdout.strip()).splitlines()
        last = why[-1].strip() if why else "no data"
        if missing:
            last = "missing " + last.split("'")[-2].split(os.sep)[-1]                 if "'" in last else last
        print(f"            skipped - {last}")
        return
    print(p.stdout[-2000:])
    print(p.stderr[-2000:])
    raise SystemExit(f"{script} failed ({p.returncode}); stopping so the "
                     f"page is not built from half-written data")


def main(merge_dirs):
    if merge_dirs:
        print(f"[0/{len(CHAIN)}] merge_fits.py  <- {', '.join(merge_dirs)}")
        run("merge_fits.py", *merge_dirs)
    for i, (script, what, optional) in enumerate(CHAIN, 1):
        print(f"[{i}/{len(CHAIN)}] {script:20s} {what}")
        run(script, optional=optional)
    page = os.path.join(HERE, "potential.html")
    print(f"\n{page}  ({os.path.getsize(page) / 1024:.0f} KB)")


if __name__ == "__main__":
    main(sys.argv[1:])
