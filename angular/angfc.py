#!/usr/bin/env python3
"""
Force constants for the angular potential, so the angular fits can finally be
screened for dynamical stability.

`latdyn.force_constants` refuses to run when lam2 or lam4 is non-zero, and it is
right to: its three-body blocks assume phi3 depends on r1 + r2 alone, so with an
angular factor present every d/dtheta term would be silently missing.  That left
fourteen fits that reproduce every elastic constant exactly and whose phonons
nobody had ever looked at - which is precisely the trap that produced a lithium
fit scoring 7.6 % with a quarter of its modes imaginary.

The way round it does not need analytic derivatives and does not need a
supercell.  Two facts make it cheap:

  * the two-body force constants are untouched by the angular factor, so
    `latdyn.force_constants` already has them right;
  * a single triplet's energy depends only on its two leg vectors,

        E = phi3(|d1| + |d2|) * h(d1.d2 / |d1||d2|),

    six variables, so its 6x6 Hessian can be differenced numerically and then
    placed into the (i, ja, jb) blocks by the chain rule.  With
    d1 = X_ja - X_i and d2 = X_jb - X_i those blocks are exact algebra, not an
    approximation.

The three-body cutoff is short - 0.97 a on niobium against 2.60 a for the pair
term - so this is 28 triplets per atom and a few thousand scalar evaluations,
seconds per element rather than the hours a displaced 512-atom supercell would
have cost.

Correctness is not assumed.  With lam2 = lam4 = 0 the result must reproduce
`latdyn.force_constants` to numerical precision, and `selftest` checks that.

    python angfc.py            # verify against the analytic version
    python angfc.py --screen   # screen every angular fit for stability
"""
import json
import glob
import os
import sys

import numpy as np

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
#  Displacement for the triplet Hessian, Angstrom.  A second derivative by
#  differences balances round-off (~eps/h^2) against truncation (~h^2), and the
#  scan over Al Nb Rh W bottoms out here: relative error against the analytic
#  force constants is 1e-5 at h = 1e-5, 2e-8 at 3e-4, and back to 6e-6 at 1e-2.
H = 3e-4


def _tri_energy(pot, X):
    """
    Triplet energies for a whole batch: X is (N, 6), columns (d1 | d2).

    Batched deliberately.  One triplet at a time, the 73-point stencil below is
    73 Python calls into phi3 and ang for each of ~28 triplets, and a screened
    evaluate() took 1.4-2.6 s against 10 ms unscreened - a factor of 200, which
    makes a 400-restart search impossible.  Vectorised, each stencil point is a
    single numpy call over every triplet at once.
    """
    d1, d2 = X[:, :3], X[:, 3:]
    r1 = np.sqrt((d1 * d1).sum(1))
    r2 = np.sqrt((d2 * d2).sum(1))
    c = (d1 * d2).sum(1) / (r1 * r2)
    e = pot.phi3(r1 + r2) * pot.ang(c)
    #  The cutoff switch, per leg.  Everything downstream differences THIS
    #  function, so the switch derivatives come out of the stencil on their own
    #  - the standalone tree needed the three-body Hessian blocks generalised
    #  by hand for the same change.
    if getattr(pot, "taper", None) is not None and pot.rcut3:
        e = e * (pot.switch(r1, pot.rcut3, 0) * pot.switch(r2, pot.rcut3, 0))
    return e


def _tri_hessians(pot, X0, h=None):
    """
    (N, 6, 6) Hessians of the triplet energy with respect to (d1, d2).

    Central differences in both variables: the mixed term needs four points and
    the diagonal three, which is the standard stencil and keeps the error at
    O(h^2) rather than O(h).
    """
    #  read the module global at call time, not as a default argument - bound
    #  at definition time, a default made H untunable and a step-size scan came
    #  back with the identical error for every h from 1e-5 to 3e-3
    h = H if h is None else h
    n, N = 6, len(X0)
    Hs = np.zeros((N, n, n))
    f0 = _tri_energy(pot, X0)
    E = np.eye(n) * h
    for a in range(n):
        Hs[:, a, a] = (_tri_energy(pot, X0 + E[a]) - 2.0 * f0
                       + _tri_energy(pot, X0 - E[a])) / (h * h)
        for b in range(a + 1, n):
            v = (_tri_energy(pot, X0 + E[a] + E[b])
                 - _tri_energy(pot, X0 + E[a] - E[b])
                 - _tri_energy(pot, X0 - E[a] + E[b])
                 + _tri_energy(pot, X0 - E[a] - E[b])) / (4.0 * h * h)
            Hs[:, a, b] = Hs[:, b, a] = v
    return Hs


def force_constants(cry, pot):
    """
    Phi[(i, j, R)] = d2E / du_i du_j, angular term included.

    Same key convention and same sign convention as latdyn.force_constants, so
    it is a drop-in replacement for it.
    """
    Phi = {}

    def add(i, j, R, blk):
        k = (i, j, tuple(int(x) for x in R))
        Phi[k] = Phi.get(k, np.zeros((3, 3))) + blk

    #  ---- two-body: unchanged by the angular factor ----
    for (i, j, R, d, r) in L.neighbours(cry, pot.rcut2):
        n = d / r
        B = L._pair_block(r, n, pot.phi2(r, 1), pot.phi2(r, 2))
        add(i, j, R, -B)
        add(i, i, (0, 0, 0), B)

    if not pot.C:
        return Phi

    #  ---- three-body: numerical Hessian per triplet, exact placement ----
    #
    #  E = f(d1, d2),  d1 = X_ja - X_i,  d2 = X_jb - X_i, so with
    #  g1 = df/dd1 and g2 = df/dd2 the forces are dE/dX_ja = g1,
    #  dE/dX_jb = g2, dE/dX_i = -(g1 + g2), and differentiating once more gives
    #  the nine blocks below.  Nothing here is approximated - only H itself is
    #  numerical.
    tri = L.triplets(cry, pot)
    if not tri:
        return Phi
    X0 = np.array([np.concatenate([na * ra, nb_ * rb])
                   for (_, _, _, ra, rb, na, nb_) in tri])
    HS = _tri_hessians(pot, X0)
    for t, (i, (ja, Ra), (jb, Rb), ra, rb, na, nb_) in enumerate(tri):
        Hs = HS[t]
        H11, H12 = Hs[:3, :3], Hs[:3, 3:]
        H21, H22 = Hs[3:, :3], Hs[3:, 3:]

        add(ja, ja, (0, 0, 0), H11)
        add(jb, jb, (0, 0, 0), H22)
        add(ja, jb, Rb - Ra, H12)
        add(jb, ja, Ra - Rb, H21)
        add(i, i, (0, 0, 0), H11 + H12 + H21 + H22)
        add(i, ja, Ra, -(H11 + H21))
        add(ja, i, -Ra, -(H11 + H12))
        add(i, jb, Rb, -(H12 + H22))
        add(jb, i, -Rb, -(H21 + H22))
    return Phi


def _tri_gradients(pot, X0, h=None):
    """(N, 6) gradients of the triplet energy with respect to (d1, d2)"""
    h = H if h is None else h
    n = 6
    G = np.zeros((len(X0), n))
    E = np.eye(n) * h
    for a in range(n):
        G[:, a] = (_tri_energy(pot, X0 + E[a])
                   - _tri_energy(pot, X0 - E[a])) / (2.0 * h)
    return G


def gradient(cry, pot):
    """
    dE/dr_i for every atom, angular term included.

    latdyn._gradient computes the three-body force as phi3'(r1 + r2) times the
    two leg directions, which is right only while phi3 has no angular factor -
    and unlike force_constants it does not refuse, it just returns the wrong
    number.  That matters for hcp and only for hcp: the non-affine correction in
    elastic() needs the forces, and a one-atom cell returns before reaching it,
    which is why the cubic angular results were never affected.
    """
    G = np.zeros((len(cry.frac), 3))
    for (i, j, R, d, r) in L.neighbours(cry, pot.rcut2):
        f = 0.5 * pot.phi2(r, 1) * (d / r)
        G[i] -= f
        G[j] += f
    if not pot.C:
        return G
    tri = L.triplets(cry, pot)
    if not tri:
        return G
    X0 = np.array([np.concatenate([na * ra, nb_ * rb])
                   for (_, _, _, ra, rb, na, nb_) in tri])
    GR = _tri_gradients(pot, X0)
    for t, (i, (ja, _), (jb, _), _, _, _, _) in enumerate(tri):
        g1, g2 = GR[t, :3], GR[t, 3:]
        G[ja] += g1
        G[jb] += g2
        G[i] -= g1 + g2
    return G


#  ---------------------------------------------------------------------------


def _build(el, p):
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    pot = L.Potential(p["m"], p["D"], p["alpha"], p["r0"], p["gamma"],
                      C=p["C"], alpha3=p["alpha3"],
                      rcut2=p["rcut2"], rcut3=p["rcut3"],
                      lam2=p.get("lam2", 0.0), lam4=p.get("lam4", 0.0))
    return cry, pot


def verify(run="runs/truba_ang2_off"):
    """
    With lam = 0 this must reproduce the analytic force constants exactly.

    Uses real converged fits rather than invented probe parameters - an
    arbitrary (m, gamma, C, s3) usually has no solution at all, and the first
    version of this check sailed through reporting PASS having tested nothing.
    """
    best = {}
    for f in glob.glob(os.path.join(HERE, run, "dense_*.json")):
        for el, v in json.load(open(f)).items():
            if isinstance(v, dict) and "score" in v and \
                    (el not in best or v["score"] < best[el]["score"]):
                best[el] = v
    print(f"lam = 0 check against the analytic force constants ({run})")
    print(f"{'el':4s}{'blocks':>8s}{'max |diff|':>13s}{'max |Phi|':>12s}"
          f"{'relative':>11s}")
    ok, n = True, 0
    for el, p in sorted(best.items()):
        if p.get("lam2") or p.get("lam4"):
            continue                      # analytic side cannot do those
        cry, pot = _build(el, p)
        A = L.force_constants(cry, pot)
        N = force_constants(cry, pot)
        keys = set(A) | set(N)
        z = np.zeros((3, 3))
        dmax = max(np.abs(A.get(k, z) - N.get(k, z)).max() for k in keys)
        pmax = max(np.abs(A.get(k, z)).max() for k in keys)
        rel = dmax / pmax
        ok &= rel < 1e-6
        n += 1
        print(f"{el:4s}{len(keys):8d}{dmax:13.2e}{pmax:12.4f}{rel:11.2e}")
    #  a check that tested nothing is a failure, not a pass
    if not n:
        print("FAIL - no lam = 0 fits found, nothing was compared")
        return False
    #  written as a plain conditional rather than inside the f-string: a
    #  multi-line expression in an f-string needs Python 3.12, and the cluster
    #  runs 3.9
    verdict = ("relative agreement better than 1e-6" if ok else
               "the placement algebra or the stencil is wrong")
    print(f"{'PASS' if ok else 'FAIL'} - {n} elements, {verdict}")
    return ok


def screen(run="runs/truba_ang2_on"):
    """dynamical stability of every fit in a run directory"""
    best = {}
    for f in glob.glob(os.path.join(HERE, run, "dense_*.json")):
        for el, v in json.load(open(f)).items():
            if isinstance(v, dict) and "score" in v and \
                    (el not in best or v["score"] < best[el]["score"]):
                best[el] = v
    print(f"stability of {run}")
    print(f"{'el':4s}{'rms%':>7s}{'lam2':>7s}{'lam4':>7s}"
          f"{'min freq 8^3 U 9^3':>20s}{'near Gamma':>12s}  verdict")
    print("-" * 72)
    out = {}
    for el, p in sorted(best.items()):
        cry, pot = _build(el, p)
        Phi = force_constants(cry, pot)
        mins = []
        for n in (8, 9):
            for q in L.mesh(n):
                mins.append(L.frequencies(cry, pot, q, Phi).min())
        m_mesh = float(min(mins)) * L.CM1
        import fit as F
        g = min(L.frequencies(cry, pot, q, Phi).min()
                for q in F.NEAR_GAMMA) * L.CM1
        stable = m_mesh > -1e-3 and g > -1e-3
        out[el] = {"min_mesh_cm1": round(m_mesh, 3),
                   "min_near_gamma_cm1": round(float(g), 3),
                   "stable": bool(stable)}
        print(f"{el:4s}{p['score']*100:7.2f}{p.get('lam2',0):7.3f}"
              f"{p.get('lam4',0):7.3f}{m_mesh:20.2f}{g:12.2f}  "
              f"{'stable' if stable else 'UNSTABLE'}", flush=True)
    n = sum(1 for v in out.values() if v["stable"])
    print(f"\n{n}/{len(out)} dynamically stable")
    p = os.path.join(HERE, "angular_stability.json")
    json.dump(out, open(p, "w"), indent=1, sort_keys=True)
    print("wrote", p)
    return out


if __name__ == "__main__":
    if "--screen" in sys.argv:
        screen()
    else:
        verify()
