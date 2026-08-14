#!/usr/bin/env python3
"""
Export the UG (angular) results in a form the standalone viewer can merge.

The two trees cannot import each other - both define latdyn, and the angular one
carries the Legendre factor - so nothing here is shared by import.  This writes a
plain JSON file and `standalone/add_ug_overlay.py` reads it.  That keeps the
boundary at a file rather than at a sys.path trick.

The record is built by `build_library.record`, the same function that builds the
MAU entries on this side, so a UG element carries exactly what a MAU element
carries: the full elastic tensor, the frozen-ion tensor, cohesive energy, bulk
modulus, residual pressure, the mechanical analysis and its polar sections, the
dispersion along the standard path, and the phonon thermodynamics to 800 K.
Written by hand it carried five elastic constants and nothing else, and the page
showed UG as a two-column table beside a fully worked MAU page.

Added on top of that record, because they are questions only the angular form
raises:
    lam2, lam4        the two Legendre coefficients
    dyn               stability on the 8^3 U 9^3 union mesh, screened with angfc
    exp_phonon        measured frequencies at X/L or H/N - out of sample, since
                      the fit sees only the q -> 0 limit

    python export_ug.py                       # default run directory
    python export_ug.py runs/truba_ang5_on runs/2026-08-03_ug

Writes ug_results.json.
"""
import glob
import json
import os
import sys

import numpy as np

import add_exp_phonon as AEP
import angfc
import build_library as B
import fit as F
import latdyn as L
import refdata
import refdata_phonon as RP

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "ug_results.json")
DEFAULT_RUN = "runs/truba_ang5_on"
TOL = -1e-3          # cm^-1; odd meshes contain Gamma, where zero reads -1e-5
THZ = 33.35641       # cm^-1 per THz

#  The measured q-points, taken from add_exp_phonon rather than copied.  A copy
#  here already went stale once: hcp was added there and not here, so the export
#  died on the first hexagonal element with a KeyError.  One table, one place.
QEXP = AEP.Q

#  fields the merged record does not need: they duplicate the MAU entry exactly
#  (the crystal and the reference data are properties of the element, not of the
#  fit) and would grow library.json for nothing.
DROP = ("exp", "struct", "a0", "c_over_a", "mass", "S298", "Cp298")


#  factors a run could have been made with; the first that reproduces the stored
#  elastic constants is the one it was made with
KNOWN_FACTORS = (1.50, 1.12)
CIJ_TOL = 2e-3          # relative; elastic() is a finite difference, not exact


def repair_rcut3(el, e, p):
    """
    Put back the three-body cutoff the fit was actually made with.

    Until this was fixed, fit.py wrote `rcut3 = dnn * 1.12` into every record
    regardless of RCUT3_OVER_DNN, so the ang5 runs - made at 1.50 - claim a
    first-shell cutoff.  Rebuilding from the record then gives a different
    potential: aluminium's stored C11 = 108.2 GPa comes back as 48.4.  The
    record's own Cij settles it, so the repair is a measurement rather than a
    guess, and it refuses instead of guessing when the evidence is thin.

    Taking the first factor inside a loose tolerance is not enough.  Where the
    three-body coefficient is small the two cutoffs give elastic constants that
    differ by well under a per cent - beryllium 192.9 against 191.9 GPa - so a
    2 % test accepts whichever was tried first and silently relabels a genuine
    first-shell fit.  Compare all of them and take the closest.

    When two radii agree to within the tolerance the measurement cannot decide,
    and then the record's own label stands rather than the export refusing.
    Indistinguishable radii describe the same potential, so which label is
    attached has no physical consequence; the measurement exists to catch a
    label that is wrong, and a label can only be wrong when the two differ -
    which is exactly when the measurement can see it.  Refusing instead cost
    scandium its three best fits on the MAU side before the same rule was put
    into merge_fits.py.
    """
    want = p.get("Cij")
    if not want:
        return p
    idx = {"C11": (0, 0), "C12": (0, 1), "C13": (0, 2),
           "C33": (2, 2), "C44": (3, 3)}
    dev = []
    for f in KNOWN_FACTORS:
        q = dict(p, rcut3=p["dnn"] * f)
        try:
            cry, pot = B.build(el, q, e)
            C, _ = L.elastic(cry, pot)
        except (np.linalg.LinAlgError, ValueError, OverflowError):
            continue
        d = max(abs(float(C[idx[k]]) - v) / max(abs(v), 1.0)
                for k, v in want.items() if k in idx)
        dev.append((d, f, q))
    dev.sort(key=lambda t: t[0])
    if not dev or dev[0][0] > CIJ_TOL:
        raise SystemExit(
            f"{el}: no cutoff in {KNOWN_FACTORS} reproduces the stored elastic "
            f"constants {want}; the record cannot be rebuilt")
    if len(dev) > 1 and dev[1][0] <= CIJ_TOL:
        return p                          # indistinguishable; the label stands
    return dev[0][2]


def candidates(runs):
    """
    {element: [fit, ...] sorted by score} over every seed in every run given.

    Several directories rather than one, because the angular fits arrived in
    batches: fourteen elements in one run and the remaining twenty four later,
    on three cluster nodes.  Reading only the newest would have dropped the
    first fourteen from the page; copying them together into one directory
    would have lost which run each came from, and the run is recorded per
    element in the output.
    """
    out = {}
    for run in runs:
        for f in glob.glob(os.path.join(HERE, run, "dense_*.json")):
            try:
                rec = json.load(open(f))
            except ValueError:
                continue
            for el, v in rec.items():
                if isinstance(v, dict) and "score" in v:
                    out.setdefault(el, []).append(dict(v, _run=run))
    for el in out:
        out[el].sort(key=lambda v: v["score"])
    return out


def pick(el, e, cands):
    """
    Lowest-score fit that is dynamically stable, else the lowest-score fit.

    Taking the lowest score outright put aluminium on the page at 4.75 % with a
    branch at -70 cm^-1, beside a MAU entry at 7.00 % that is stable - an
    improvement only if you do not look at the phonons.  Elastic exactness says
    nothing about stability: this run has four fits at 0.00 % RMS with imaginary
    modes.  Candidates are screened in score order and the first stable one wins,
    so the usual cost is one screen.
    """
    fallback = None
    for p in cands:
        p = repair_rcut3(el, e, p)
        cry, pot = B.build(el, p, e)
        Phi = angfc.force_constants(cry, pot)
        mesh, gam = stability(cry, pot, Phi)
        if mesh > TOL and gam > TOL:
            return p, (mesh, gam), True
        if fallback is None:
            fallback = (p, (mesh, gam))
    return fallback[0], fallback[1], False


def stability(cry, pot, Phi):
    """(min frequency on 8^3 U 9^3, min near Gamma), cm^-1"""
    qs = np.vstack([L.mesh(8), L.mesh(9)])
    mesh = float(L.frequencies_many(cry, pot, qs, Phi).min()) * L.CM1
    gam = float(L.frequencies_many(cry, pot, np.asarray(F.NEAR_GAMMA),
                                   Phi).min()) * L.CM1
    return mesh, gam


def exp_phonon(el, struct, cry, pot, Phi):
    """measured frequencies against ours, or None where none are published"""
    want = RP.for_element(el, struct)
    if not want:
        return None
    pts, errs = [], []
    for name, exp in want.items():
        q = np.array(QEXP[struct][name])
        ours = sorted(float(x) for x in
                      L.frequencies(cry, pot, q, Phi) * L.CM1 / THZ)
        if name in RP.OPTIC_ONLY:
            #  at Gamma the acoustic branches are zero whatever the potential
            #  does, so only the optic ones carry information
            ours = ours[-len(exp):]
        e = 100 * np.mean(np.abs(np.array(ours) - np.array(exp))
                          / np.array(exp))
        errs.append(e)
        pts.append({"name": name, "q": list(map(float, q)),
                    "exp": [round(float(x), 2) for x in exp],
                    "ours": [round(x, 2) for x in ours],
                    "err_pct": round(float(e), 1)})
    return {"points": pts, "mae": round(float(np.mean(errs)), 1),
            "ref": RP.PHONON_EXP[el]["ref"]}


def main():
    runs = sys.argv[1:] or [DEFAULT_RUN]
    cands = candidates(runs)
    if not cands:
        raise SystemExit(f"no fits found under {', '.join(runs)}")

    print(f"{'el':4s}{'rms%':>8s}{'lam2':>8s}{'lam4':>8s}"
          f"{'min freq':>11s}{'near G':>9s}{'phonon':>8s}  status")
    print("-" * 64)
    out = {}
    for el, cs in sorted(cands.items()):
        e = refdata.ELEMENTS[el]
        p, (mesh, gam), stable = pick(el, e, cs)
        d = B.record(el, p, e)
        cry, pot = B.build(el, p, e)
        #  screened with angfc explicitly: latdyn.force_constants does delegate
        #  there when the Legendre coefficients are non-zero, but a silent
        #  fallback to the angle-free constants would call an unstable fit
        #  stable, and that is the one error this screen exists to prevent
        Phi = angfc.force_constants(cry, pot)
        d["dyn"] = {"min_mesh_cm1": round(mesh, 3),
                    "min_near_gamma_cm1": round(gam, 3),
                    "stable": bool(stable),
                    "n_screened": len(cs),
                    "best_score_rms": round(cs[0]["score"] * 100.0, 2)}
        ph = exp_phonon(el, e["struct"], cry, pot, Phi)
        if ph:
            d["exp_phonon"] = ph
        d["run"] = p.get("_run", runs[0])
        for k in DROP:
            d.pop(k, None)
        out[el] = d
        mae = f"{ph['mae']:6.1f}%" if ph else "     -"
        print(f"{el:4s}{d['rms']:8.2f}{d['lam2']:8.3f}{d['lam4']:8.3f}"
              f"{mesh:11.2f}{gam:9.2f}{mae:>8s}  "
              f"{'stable' if stable else 'UNSTABLE'}", flush=True)

    tmp = OUT + ".tmp"
    json.dump(out, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, OUT)
    n = sum(1 for v in out.values() if v["dyn"]["stable"])
    print(f"\n{len(out)} elements, {n} of them dynamically stable")
    print("written:", OUT)


if __name__ == "__main__":
    main()
