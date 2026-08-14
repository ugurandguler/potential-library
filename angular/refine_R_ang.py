#!/usr/bin/env python3
"""
Does the angular term actually lower the C44/C' floor?

This is the question the whole angular extension stands or falls on.  Eight
cubic metals sit below the floor the published form can reach, and no amount of
searching moves them - that was established by refine_R.py.  If phi3 gains an
angular factor

    phi3(r1, r2) -> phi3(r1 + r2) [ 1 + lam2 P2(cos t) + lam4 P4(cos t) ]

the reachable region may open up.  Or it may not, and a fit that merely scores
better could be exploiting slack elsewhere in the objective.  So measure the
floor itself rather than a fitted RMS.

**Seeding decides the answer, and one seed is not enough.**  A 6-dimensional
search started cold is weaker than the 4-dimensional one it is meant to improve
on, and the first attempt at this comparison reported *higher* floors with more
freedom - arithmetically impossible for an infimum, and purely a search
artefact.  So the search starts from the argmin of refine_R.py with
lam2 = lam4 = 0, which bounds the result by the isotropic floor.

That alone is still not enough.  Seeded only from the origin, platinum came back
at 1.705 - above the 1.594 an already-converged angular fit had reached, so the
"floor" sat above a point known to be attainable.  Every converged fit in runs/
is an existence proof and is now fed in as a start as well.  Every number here
is an upper bound on the infimum: it can be lowered by a better search, never
raised, and a value above any achieved fit means the search failed, not that the
region is closed.

    python refine_R_ang.py Nb V Mo Cr Al Ir Pt Rh

Writes R_floor_ang.json.
"""
import glob
import json
import os
import sys

import numpy as np

import dense_fit as DF
import fit as F
import refdata

F.REQUIRE_DYNAMICAL_STABILITY = False

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(os.path.dirname(HERE), "standalone", "R_floor.json")
NDIM = 6
ITERS = 600
#  extra starts around the seed, so a single bad simplex cannot decide it
JITTER = [0.0, 0.15, 0.35, 0.7]
LAM_START = [(0.0, 0.0), (0.6, 0.0), (-0.6, 0.0), (0.0, 0.6), (0.9, -0.5)]


def R_of(el, e, v):
    m, g, s, C, l2, l4 = v
    if not (F.M_MIN <= m <= F.M_MAX) or g < 0 or s <= 0.05 or abs(C) > F.C_MAX:
        return 1e6
    if not (DF.LAM2_MIN <= l2 <= DF.LAM2_MAX
            and DF.LAM4_MIN <= l4 <= DF.LAM4_MAX):
        return 1e6
    #  the angular factor has to stay positive over the whole sphere, otherwise
    #  some triplet geometries get a sign-flipped three-body energy
    if not DF._angular_ok(l2, l4):
        return 1e6
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        try:
            r = F.evaluate(el, e, m, g, C, s, lam2=l2, lam4=l4)
        except (OverflowError, ValueError, FloatingPointError,
                np.linalg.LinAlgError):
            return 1e6
    if not r:
        return 1e6
    c = r["Cij"]
    cp = 0.5 * (c["C11"] - c["C12"])
    if cp <= 1e-6 or c["C44"] <= 1e-6:      # C' <= 0 is instability, not a low R
        return 1e6
    return c["C44"] / cp


def nelder_mead(f, x0, step, iters=ITERS):
    x0 = np.asarray(x0, float)
    sim = [x0] + [x0 + np.eye(NDIM)[i] * step[i] for i in range(NDIM)]
    val = [f(s) for s in sim]
    for _ in range(iters):
        o = np.argsort(val)
        sim = [sim[i] for i in o]
        val = [val[i] for i in o]
        if val[0] > 1e5:
            break
        cen = np.mean(sim[:-1], axis=0)
        xr = cen + (cen - sim[-1])
        fr = f(xr)
        if fr < val[0]:
            xe = cen + 2.0 * (cen - sim[-1])
            fe = f(xe)
            sim[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            sim[-1], val[-1] = xr, fr
        else:
            xc = cen + 0.5 * (sim[-1] - cen)
            fc = f(xc)
            if fc < val[-1]:
                sim[-1], val[-1] = xc, fc
            else:
                for i in range(1, NDIM + 1):
                    sim[i] = sim[0] + 0.5 * (sim[i] - sim[0])
                    val[i] = f(sim[i])
    i = int(np.argmin(val))
    return val[i], sim[i]


def fitted_starts(el):
    """
    Starting points taken from angular fits that already converged.

    Seeding from the lam = 0 argmin bounds the answer by the isotropic floor but
    says nothing about where the *angular* optimum lives, and for the fcc metals
    it is nowhere near: the converged platinum fit sits at lam2 = 0.60,
    lam4 = -1.08, which a simplex stepping 0.35 from the origin does not reach.
    Seeded only from the origin this search returned 1.705 for platinum while a
    completed fit had already achieved 1.594 - a floor above a point known to be
    attainable, which is a search failure and nothing else.  Any converged fit is
    an existence proof, so feed them all in.
    """
    #  every run directory, not a name pattern.  Matching only "angular_*" quietly
    #  skipped the cluster results in runs/truba_ang_on/, and the search then
    #  reported rhodium out of reach at 1.720 while a fit sitting in that very
    #  directory had already hit 1.680 exactly.
    out = []
    for f in glob.glob(os.path.join(HERE, "runs", "*", f"dense_{el}_*.json")):
        try:
            rec = json.load(open(f)).get(el)
        except ValueError:
            continue
        if isinstance(rec, dict) and "m" in rec:
            out.append([rec["m"], rec["gamma"], rec["s3"], rec["C"],
                        rec.get("lam2", 0.0), rec.get("lam4", 0.0)])
    return out


def floor_ang(el, seed):
    e = refdata.ELEMENTS[el]
    m, g, s, C = seed["argmin"]
    best, arg = seed["R_floor"], list(seed["argmin"]) + [0.0, 0.0]
    rng = np.random.default_rng(abs(hash(el)) % (2**32))

    for x0 in fitted_starts(el):
        x0 = np.array(x0, float)
        step = [max(0.3, 0.15 * abs(x0[0])), 0.25, 0.25,
                max(1.0, 0.25 * abs(x0[3])), 0.35, 0.35]
        v, x = nelder_mead(lambda z: R_of(el, e, z), x0, step)
        #  the start itself counts: a converged fit proves its own R attainable
        for cand, val in ((x0, R_of(el, e, x0)), (x, v)):
            if val < best:
                best, arg = float(val), [float(t) for t in cand]

    for j in JITTER:
        for l2, l4 in LAM_START:
            x0 = np.array([m, g, s, C, l2, l4], float)
            if j:
                x0[:4] *= 1.0 + j * rng.standard_normal(4)
            step = [max(0.3, 0.15 * abs(x0[0])), 0.25, 0.25,
                    max(1.0, 0.25 * abs(x0[3])), 0.35, 0.35]
            v, x = nelder_mead(lambda z: R_of(el, e, z), x0, step)
            if v < best:
                best, arg = float(v), [float(t) for t in x]
    return best, arg


def main():
    base = json.load(open(BASE))
    els = sys.argv[1:] or sorted(base)
    #  Each element costs about ten minutes, so this gets run in batches.  Load
    #  what is already there and add to it - writing a fresh dict silently threw
    #  away the previous batch.
    p = os.path.join(HERE, "R_floor_ang.json")
    out = json.load(open(p)) if os.path.exists(p) else {}
    print(f"{'el':4s}{'R_exp':>8s}{'iso floor':>11s}{'ang floor':>11s}"
          f"{'gain':>8s}   {'lam2':>6s}{'lam4':>7s}  verdict")
    print("-" * 74)
    for el in els:
        if el not in base:
            continue
        s = base[el]
        R, arg = floor_ang(el, s)
        iso, exp = s["R_floor"], s["R_exp"]
        was, now = exp >= iso, exp >= R
        note = ("still out of reach" if not now else
                "NOW REACHABLE" if not was else "was already reachable")
        out[el] = {"R_exp": exp, "R_floor_iso": iso, "R_floor_ang": R,
                   "argmin": arg, "reachable_iso": was, "reachable_ang": now,
                   "struct": s["struct"]}
        print(f"{el:4s}{exp:8.3f}{iso:11.3f}{R:11.3f}"
              f"{100*(iso-R)/iso:7.1f}%   {arg[4]:6.2f}{arg[5]:7.2f}  {note}",
              flush=True)
    json.dump(out, open(p, "w"), indent=1, sort_keys=True)
    n = sum(1 for v in out.values()
            if v["reachable_ang"] and not v["reachable_iso"])
    print(f"\n{n} element(s) brought into reach by the angular term")
    print("wrote", p)


if __name__ == "__main__":
    main()
