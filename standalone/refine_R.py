#!/usr/bin/env python3
"""
How low can R = C44/C' actually go?

Both cheap answers are only upper bounds on the infimum, and they disagree in
both directions: the coarse grid in cprime_ceiling.py reached R = 1.57 for W
where a simplex from a sparse start reached 0.806, while for Al the grid reached
1.68 and the simplex stalled at 2.171.  Trusting either alone flips elements
between "the form cannot do this" and "the search was not good enough", which is
the one distinction this whole analysis exists to make.

So: run the full grid, keep the lowest-R points, and start a simplex from each.
Grid-seeded refinement is bounded above by the grid and can only improve on it.

Analytic check.  For a pair potential on a Bravais lattice
    C_ijkl = (1/2V) sum_R [phi'' - phi'/R] R_i R_j R_k R_l / R^2
so C12 = C44 identically (Cauchy).  Summing the fcc nearest-neighbour shell
alone gives C11 = 2 C12, i.e. R = 2 exactly; the bcc nearest-neighbour shell
gives C11 = C12, i.e. C' = 0 and R unbounded, so for bcc the value is set by how
much second-shell weight the potential carries and there is no shell-level floor.
That is why the fcc numbers cluster near 2 and the bcc ones sit near 1.

    python refine_R.py Nb V Cr Mo W Ta Al Cu Ni

Writes R_floor.json.
"""
import json
import os
import sys

import numpy as np

import fit as F
import refdata

#  Measure the floor over the potentials that are actually usable.
#
#  With the three-body sum truncated at the first shell this made little
#  difference, and the flag was off so that the answer described the form's
#  reach rather than the search.  At 1.50 d_nn it changes everything: an fcc
#  metal carries 153 triplets per atom instead of 66, and with that much freedom
#  the minimiser reaches a degenerate corner where C44 goes to zero - gold comes
#  back as C11 = 177.6, C12 = 168.7, C44 = 0.0000.  The ratio C44/C' is then
#  zero, so thirteen of twenty three elements reported a floor of zero and the
#  test stopped discriminating.  A crystal with no shear resistance is not a
#  potential anyone can use; the honest question is the lowest anisotropy the
#  form reaches among potentials that stand up.
#
#  It costs nothing.  Measured over a spread of trial points for niobium:
#  14 ms each with the screen on against 20 ms with it off, because an unstable
#  candidate dies at the near-Gamma check before the rest of the work is done.
F.REQUIRE_DYNAMICAL_STABILITY = (
    os.environ.get("FLOOR_REQUIRE_STABLE", "1") != "0")

HERE = os.path.dirname(os.path.abspath(__file__))
N_SEED = 24                   # lowest-R grid points to refine from
SHEAR_MIN_FRAC = 0.01         # C44 and C' as a fraction of C11; see R_of

M_GRID = [1.2, 1.6, 2.2, 3.0, 4.0, 5.5, 7.5, 10.0, 13.0, 16.0, 20.0]
G_GRID = [0.0, 0.25, 0.5, 1.0, 1.6, 2.4, 3.5, 5.0]
S_GRID = [0.4, 0.7, 1.0, 1.4, 2.0, 3.0]
C_GRID = [-60, -40, -25, -15, -8, -4, -1.5, 0.0, 1.5, 4, 8, 15, 25, 40, 60]


def R_of(el, e, v):
    m, g, s, C = v
    if not (F.M_MIN <= m <= F.M_MAX) or g < 0 or s <= 0.05 or abs(C) > F.C_MAX:
        return 1e6
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        try:
            r = F.evaluate(el, e, m, g, C, s)
        except (OverflowError, ValueError, FloatingPointError,
                np.linalg.LinAlgError):
            return 1e6
    if not r:
        return 1e6
    c = r["Cij"]
    cp = 0.5 * (c["C11"] - c["C12"])
    #  C' <= 0 is not a smaller R, it is an unstable crystal.
    #
    #  An absolute threshold does not say that.  With the screen on and the
    #  three-body sum out to 1.50 d_nn the minimiser still walks into the
    #  degenerate corner and clears 1e-6 by a hair: gold came back as
    #  C11 = 177.616, C12 = 168.692, C44 = 1.018e-6 GPa, i.e. C44/C11 = 6e-9,
    #  and six fcc metals reported a floor of 2e-7 instead of a number.
    #
    #  The dynamical screen does not catch it, and should not be expected to.
    #  C44 -> 0 leaves the crystal marginal, not unstable: the transverse
    #  acoustic branch along [100] goes as q^2 instead of q, so every frequency
    #  on the mesh is real and non-negative.  A search for imaginary modes has
    #  nothing to find.
    #
    #  So the shear constants are held to a fraction of the tensor's own scale.
    #  One per cent of C11 is an order of magnitude below the softest metal
    #  measured - niobium, C44/C11 = 0.116 - so it excludes no physical case,
    #  which matters because low R is exactly what this function looks for and
    #  the bcc metals genuinely have small C44.
    if cp <= SHEAR_MIN_FRAC * c["C11"] or c["C44"] <= SHEAR_MIN_FRAC * c["C11"]:
        return 1e6
    return c["C44"] / cp


def nelder_mead(f, x0, step, iters=400):
    x0 = np.asarray(x0, float)
    sim = [x0] + [x0 + np.eye(4)[i] * step[i] for i in range(4)]
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
                for i in range(1, 5):
                    sim[i] = sim[0] + 0.5 * (sim[i] - sim[0])
                    val[i] = f(sim[i])
    i = int(np.argmin(val))
    return val[i], sim[i]


def floor_R(el):
    e = refdata.ELEMENTS[el]

    def f(v):
        return R_of(el, e, v)

    scored = []
    for m in M_GRID:
        for g in G_GRID:
            for s in S_GRID:
                for C in C_GRID:
                    v = f((m, g, s, C))
                    if v < 1e5:
                        scored.append((v, (m, g, s, C)))
    if not scored:
        return None
    scored.sort(key=lambda t: t[0])
    grid = scored[0][0]
    best, arg = grid, scored[0][1]
    #  Two seed families, because neither dominates.  Grid-seeded starts found
    #  0.902 for Nb where a sparse spread found 1.081, but for W the sparse
    #  spread found 0.806 where grid-seeded stalled at 1.032 - the grid's best-R
    #  points are not always in the basin that goes lowest.  Every floor here is
    #  an upper bound, so run both and keep the smaller.
    sparse = [(m, g, s, C)
              for m in (1.5, 3.0, 6.0, 12.0, 20.0)
              for g in (0.0, 1.0, 3.0)
              for s in (0.5, 1.0, 2.0)
              for C in (-40.0, -8.0, 0.0, 8.0, 40.0)]
    for x0 in [x for _, x in scored[:N_SEED]] + sparse:
        r, x = nelder_mead(f, x0, [1.0, 0.4, 0.3, 8.0])
        if r < best:
            best, arg = r, x
    return grid, best, arg


def main():
    els = sys.argv[1:] or ["Nb", "V", "Cr", "Mo", "W", "Ta", "Al", "Cu", "Ni"]
    #  one element per process is how this is run on all of them at once, so
    #  each writes its own file and merge_R_floor folds them together; a shared
    #  file would be corrupted by concurrent writers
    p = os.path.join(HERE, f"R_floor_{els[0]}.json" if len(els) == 1
                     else "R_floor.json")
    out = json.load(open(p)) if os.path.exists(p) else {}

    print(f"{'el':4s}{'struct':>7s}{'R_exp':>8s}{'R_grid':>9s}"
          f"{'R_floor':>10s}{'margin':>9s}{'reachable?':>12s}")
    print("-" * 60)
    for el in els:
        e = refdata.ELEMENTS[el]
        if e["struct"] == "hcp":
            continue
        c = e["Cij"]
        R_exp = c["C44"] / (0.5 * (c["C11"] - c["C12"]))
        got = floor_R(el)
        if not got:
            print(f"{el:4s}  no valid solution")
            continue
        grid, best, arg = got
        ok = best <= R_exp
        out[el] = dict(R_exp=R_exp, R_grid=grid, R_floor=best,
                       argmin=list(arg), reachable=bool(ok),
                       struct=e["struct"])
        print(f"{el:4s}{e['struct']:>7s}{R_exp:8.2f}{grid:9.2f}{best:10.3f}"
              f"{R_exp / best:9.2f}{('yes' if ok else 'NO'):>12s}")

    tmp = p + ".tmp"
    json.dump(out, open(tmp, "w"), indent=1)
    os.replace(tmp, p)

    bad = [el for el, d in out.items() if not d["reachable"]]
    print(f"\nout of reach for any parameter set: {bad or 'none'}")
    print("margin < 1 means the experimental anisotropy is below everything "
          "the form can produce.")


if __name__ == "__main__":
    main()
