#!/usr/bin/env python3
"""
Fit the Akgun-Ugur potential from its own derivatives.  gamma is continuous.

  * gamma is a real number.  Restricting it to integers was never a physical
    requirement - it came from an external input format that truncates the
    r-power - so here it is fitted freely.

  * the elastic constants in the objective are the RELAXED ones, carrying the
    non-affine internal-strain correction of Phys. Rev. Materials 7, 073603
    (2023).  Comparing frozen-ion values with experiment is the wrong quantity
    whenever the primitive cell holds more than one atom, i.e. for every hcp
    element here.

Hard constraints, exact at every trial point:
    cohesive energy, zero hydrostatic pressure at the experimental lattice
    constant, bulk modulus, and - for hcp - zero deviatoric stress so the
    experimental c/a is a genuine stationary point.

Searched: m, gamma, s3 = alpha3/alpha, C.
Predicted and scored: every independent elastic constant.

    python fit.py            # all elements
    python fit.py Pd Cu Mg

Set REQUIRE_DYNAMICAL_STABILITY = False to reproduce the elastic-only fit.
"""
import json, math, os, sys, time
import numpy as np
import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
#  Pair range, overridable so the cutoff can be varied against the taper.
RCUT2_OVER_A = float(os.environ.get("RCUT2_OVER_A", "2.6"))

#  Cutoff taper; see latdyn.Potential.  Default off, so fits made before it
#  still evaluate to their own numbers.  Travels with the record, for the
#  same reason rcut3 does.
TAPER = (float(os.environ["TAPER"]) if os.environ.get("TAPER") else None)

#  see standalone/fit.py for the shell scan behind 1.50
RCUT3_OVER_DNN = float(os.environ.get("RCUT3_OVER_DNN", "1.50"))

R0_MIN_FRAC, R0_MAX_FRAC = 0.70, 1.50
D_MIN_FRAC = 0.05
C_MAX = 60.0
#  Same two constraints the standalone tree now carries, for the same reason.
#  3.0 was effectively no constraint: chromium, molybdenum and tungsten came
#  through the tapered refit at 0.55-0.62, phi2 supplying -10 to -20 eV/atom
#  and phi3 handing back +6 to +11, and those are the three whose vacancy
#  formation energy is NEGATIVE and whose crystals come apart at 600 K.  The
#  angular set fails on exactly the same three, so it needs exactly the same
#  guard - the Legendre factor does not rescue an unphysical radial solution,
#  it just weights it by angle.
E3_OVER_E2_MAX = 0.30

#  In a metal the nearest neighbours carry the binding.  Magnitude and not
#  sign: iron's first shell repels by 0.06 eV/atom against a cohesive energy of
#  4.28 and iron is sound, chromium's by 3.06 against 4.10 and it is not.
REQUIRE_ATTRACTIVE_FIRST_SHELL = True
REQUIRE_COMPRESSION_STABLE = True
#  a basin has to beat the noise in the energy to count as one
COMPRESSION_TOL = 5e-3
KB = 8.617333262e-5
FIRST_SHELL_REPULSION_MAX = 0.15

#  m is the ratio of the repulsive to the attractive decay (alpha/beta in the
#  Malinowska-Adamska transcription of the 1998 paper).  Their Ni table spans
#  1.97 to 10.87, so anything far outside that is an artefact of the optimiser
#  exploiting a flat direction rather than physics: unconstrained, Pd ran to
#  m = 40, which reproduces the elastic constants but is a very hard core and
#  unlikely to transfer.  This is a physical prior, not a numerical bound - the
#  report flags any element that ends up sitting on it.
M_MIN, M_MAX = 1.2, 20.0

KEYS = {"fcc": ["C11", "C12", "C44"], "bcc": ["C11", "C12", "C44"],
        "hcp": ["C11", "C12", "C13", "C33", "C44"]}


def make(el, e, m, D, alpha, r0, gamma, C, s3, lam2=0.0, lam4=0.0):
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    dnn = min(r for (_, _, _, _, r) in L.neighbours(cry, e["a0"] * 1.3))
    #  three-body legs: the first neighbour shell.  For hcp that shell splits
    #  into in-plane and out-of-plane sub-shells; take both.
    rc3 = dnn * RCUT3_OVER_DNN
    pot = L.Potential(m, D, alpha, r0, gamma, C=C, alpha3=s3 * alpha,
                      rcut2=RCUT2_OVER_A * e["a0"], rcut3=rc3, lam2=lam2, lam4=lam4,
                      taper=TAPER)
    return cry, pot, dnn


def constraints(el, e, p, h=2e-4):
    """
    Only what the hard constraints need: energy per atom, hydrostatic pressure,
    bulk modulus, and the deviatoric stress.  Deliberately avoids the full 6x6
    elastic tensor - that costs ~40 energy evaluations and is only needed once,
    at the converged point, not inside every Newton step.
    """
    cry, pot, dnn = make(el, e, *p)
    nat = len(cry.frac)
    E0 = L.energy(cry, pot)
    #  isotropic scaling -> P and B
    ep = L.energy(cry.strained(np.eye(3)*h), pot)
    em = L.energy(cry.strained(np.eye(3)*(-h)), pot)
    dE = (ep - em) / (2*h)                       # dE/d(volumetric strain), per atom
    d2E = (ep - 2*E0 + em) / (h*h)
    v = cry.vol / nat
    P = dE / (3.0*v) * L.EV_A3_TO_GPA            # hydrostatic pressure
    B = d2E / (9.0*v) * L.EV_A3_TO_GPA
    dev = 0.0
    if e["struct"] == "hcp":
        #  uniaxial c strain at fixed a: sigma_zz - sigma_xx
        ez = np.zeros((3, 3)); ez[2, 2] = h
        ex = np.zeros((3, 3)); ex[0, 0] = h
        szz = (L.energy(cry.strained(ez), pot)
               - L.energy(cry.strained(-ez), pot)) / (2*h) / v
        sxx = (L.energy(cry.strained(ex), pot)
               - L.energy(cry.strained(-ex), pot)) / (2*h) / v
        dev = (szz - sxx) * L.EV_A3_TO_GPA
    return dict(cry=cry, pot=pot, dnn=dnn, E=E0, P=P, B=B, dev=dev)


def observables(el, e, p):
    """the full relaxed elastic tensor; only called at a converged point"""
    cry, pot, dnn = make(el, e, *p)
    C, Cb = L.elastic(cry, pot)
    return dict(cry=cry, pot=pot, dnn=dnn,
                Cij={k: C[i, j] for k, (i, j) in
                     (("C11", (0, 0)), ("C12", (0, 1)), ("C13", (0, 2)),
                      ("C33", (2, 2)), ("C44", (3, 3)))},
                Cij_frozen={"C11": Cb[0, 0], "C12": Cb[0, 1], "C44": Cb[3, 3]})


def solve(el, e, m, gamma, C, s3, seed=None, lam2=0.0, lam4=0.0):
    """
    Solve for (alpha, r0) - and D analytically - so that
        E/atom = -Ecoh,  P = 0,  B = B_exp   [ and sigma_zz - sigma_xx = 0 ]
    D is a pure prefactor of the energy, so it drops out of the pressure and of
    the B/Ecoh ratio; that leaves two (three for hcp) equations for two (three)
    unknowns.  hcp promotes m to an unknown to pay for the extra condition.
    """
    hexa = e["struct"] == "hcp"
    Bt, Ec = e["B"], e["Ecoh"]

    def resid(x):
        al, r0 = x[0], x[1]
        mm = x[2] if hexa else m
        if al <= 0.02 or r0 <= 0.3 or not (M_MIN <= mm <= M_MAX):
            return None
        try:
            #  The angular weights belong HERE, not only in evaluate().  Solving
            #  P = 0 and B = B_exp for the angle-free potential and then
            #  measuring the elastic constants with the angular factor switched
            #  on leaves the crystal off its own equilibrium, and the residual
            #  pressure says exactly which structures that ruins: on the first
            #  cluster run every fcc metal came back at P ~ 1e-7 while every bcc
            #  one carried 1.5 to 16 - P2 and P4 average to zero over the fcc
            #  triplet angles and do not over the bcc ones.  Nb V Mo Ta W Cr were
            #  all scored at a lattice constant that was not their own.
            o = constraints(el, e, (mm, 1.0, al, r0, gamma, C, s3, lam2, lam4))
        except (OverflowError, ValueError, FloatingPointError):
            return None
        if o["E"] >= 0:
            return None
        D = Ec / (-o["E"])                      # scale to the cohesive energy
        out = [o["P"]*D / Bt, (o["B"]*D - Bt) / Bt]
        if hexa:
            out.append(o["dev"]*D / Bt)
        return out

    n = 3 if hexa else 2
    #  r0 must be seeded RELATIVE to the nearest-neighbour distance.  Absolute
    #  values in Angstrom only work for the mid-row metals: K sits at
    #  d_nn = 4.61 A and Na at 3.72 A, so a fixed 2.4-3.4 A bracket put the
    #  starting point nowhere near the solution and Newton simply failed - which
    #  is why several elements came out worse than their seed, or not at all.
    cry0 = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"), mass=1.0)
    dnn0 = min(r for (_, _, _, _, r) in L.neighbours(cry0, e["a0"]*1.3))
    seeds = ([seed] if seed else []) + \
            [[a/dnn0*2.75, f*dnn0, m]
             for a in (0.4, 0.8, 1.4, 2.2) for f in (0.88, 1.02, 1.18)]
    for s0 in seeds:
        x = list(s0[:n])
        ok = False
        for _ in range(40):
            F = resid(x)
            if F is None:
                break
            if max(abs(v) for v in F) < 1e-7:
                ok = True
                break
            J = np.zeros((n, n))
            h = 1e-6
            for c in range(n):
                xp = list(x); xp[c] += h
                Fp = resid(xp)
                if Fp is None:
                    ok = False; break
                J[:, c] = (np.array(Fp) - np.array(F)) / h
            else:
                try:
                    x = list(np.array(x) - np.linalg.solve(J, np.array(F)))
                except np.linalg.LinAlgError:
                    break
                continue
            break
        if ok:
            return x
    return None


def evaluate(el, e, m, gamma, C, s3, seed=None, lam2=0.0, lam4=0.0):
    x = solve(el, e, m, gamma, C, s3, seed, lam2, lam4)
    if x is None:
        return None
    hexa = e["struct"] == "hcp"
    al, r0 = x[0], x[1]
    mm = x[2] if hexa else m
    cn = constraints(el, e, (mm, 1.0, al, r0, gamma, C, s3, lam2, lam4))
    D = e["Ecoh"] / (-cn["E"])
    o = observables(el, e, (mm, 1.0, al, r0, gamma, C, s3, lam2, lam4))
    dnn = o["dnn"]
    if not (R0_MIN_FRAC*dnn <= r0 <= R0_MAX_FRAC*dnn):
        return None
    if abs(C) > C_MAX or D <= D_MIN_FRAC*e["Ecoh"]:
        return None
    cij = {k: o["Cij"][k]*D for k in KEYS[e["struct"]]}
    if not born_stable(cij, e["struct"]):
        return None
    if C:
        e2 = constraints(el, e, (mm, 1.0, al, r0, gamma, 0.0, s3))["E"]*D
        if abs(-e["Ecoh"] - e2) > E3_OVER_E2_MAX*abs(e2):
            return None
    if REQUIRE_ATTRACTIVE_FIRST_SHELL:
        cry, pot, _ = make(el, e, mm, D, al, r0, gamma, C, s3)
        nb = [r for (_, _, _, _, r) in L.neighbours(cry, pot.rcut2)]
        if nb:
            r1 = min(nb)
            n1 = sum(1 for r in nb if abs(r - r1) < 1e-9)
            c1 = 0.5 * n1 * float(pot.phi2(np.array([r1]))[0]) / len(cry.frac)
            if c1 > FIRST_SHELL_REPULSION_MAX * abs(e["Ecoh"]):
                return None
    if REQUIRE_COMPRESSION_STABLE:
        cry, pot, _ = make(el, e, mm, D, al, r0, gamma, C, s3)
        if not compression_stable(el, cry, pot):
            return None
    #  cheapest tests first; this one builds the force constants, so it runs last
    #
    #  The angular case used to be exempt, because nothing could compute its
    #  force constants.  It is not exempt any more, and it turned out to matter:
    #  screened afterwards, the unconstrained angular run had four of its
    #  fourteen perfect fits dynamically unstable - Mo at -140, Na -97, Ta -98
    #  and V -239 cm^-1 - every one of them scoring 0.00 % on the elastic
    #  constants.  Exactly the trap that produced a lithium fit at 7.6 % with a
    #  quarter of its modes imaginary.
    if REQUIRE_DYNAMICAL_STABILITY:
        cry, pot, _ = make(el, e, mm, D, al, r0, gamma, C, s3, lam2, lam4)
        if not dynamically_stable(cry, pot):
            return None
    sc, n = 0.0, 0
    for k in KEYS[e["struct"]]:
        t = e["Cij"].get(k)
        if t:
            sc += ((cij[k]-t)/t)**2; n += 1
    return dict(m=mm, gamma=gamma, D=D, alpha=al, r0=r0, C=C, s3=s3,
                alpha3=s3*al, dnn=dnn, rcut2=RCUT2_OVER_A*e["a0"],
                #  must be the cutoff the potential was ACTUALLY built with.
                #  Hard-coding 1.12 here while the potential used the value in
                #  RCUT3_OVER_DNN made every ang5 record claim a first-shell
                #  cutoff it had not used, so anything rebuilding from the record
                #  - the stability screen, the viewer - rebuilt a different
                #  potential and reported its numbers under the fit's name.
                rcut3=dnn*RCUT3_OVER_DNN, taper=TAPER, Cij=cij, B=cn["B"]*D,
                Ecoh=-e["Ecoh"],
                P=cn["P"]*D, lam2=lam2, lam4=lam4,
                score=math.sqrt(sc/n) if n else 0.0, seed=[al, r0, mm])




def compression_stable(el, cry, pot):
    """Is the fitted lattice the crystal, or a ledge above a deeper basin?

    Everything else the fit constrains lives at one volume, and the elastic
    constants are curvature *at* that volume.  Nothing asks what happens when
    the cell is squeezed - and for this form the answer is not benign.  The
    cutoffs are fixed lengths, so compressing the cell pulls further shells
    inside them and the triplet count climbs; phi2's repulsive core resists,
    phi3 depends only on r1 + r2 and where C is negative does not.  Where phi3
    wins, the energy turns over and the fitted structure is a local minimum
    beside a bottomless basin.  Molecular dynamics finds it - sodium reached
    125725 K with a potential energy 14 eV/atom below its own lattice - while
    every static test reports health, because the elastic constants, the
    phonons and the vacancy energy are all properties of the fitted geometry.

    What matters is not how deep the basin is but how hard it is to reach.
    Barium's is 208 eV/atom deep and barium is fine; lithium's costs 0.001 eV
    and lithium disintegrates.  So the test is the barrier, against k_B T_melt.

    That threshold is EMPIRICAL, and the obvious derivation of it is wrong.
    Comparing the barrier with (3/2) k_B T per atom predicted a collapse at
    900 K for a set that survived 1100 K and collapsed at 1200 K - uniform
    compression of the whole cell is one collective coordinate and thermal
    fluctuations do not drive it coherently.  Calibrating against that measured
    collapse moved the threshold by 1.46, which cancels the 3/2 and leaves
    k_B T_melt.  Checked against the MD screen on every parameter set that has
    a basin at all: 30 of 30, with sodium at 0.72 failing and barium at 1.18
    passing on either side of the line.
    """
    tm = refdata.MELTING.get(el)
    if tm is None:
        return True
    lim = KB * tm
    e = refdata.ELEMENTS[el]

    def E(x):
        return L.energy(L.Crystal(e["struct"], e["a0"] * x, e.get("c_over_a"),
                                  mass=refdata.MASSES[el]), pot)

    #  the reference is the potential's OWN equilibrium.  The fit does not
    #  drive the pressure to zero, and a few meV of residual reads as a basin
    #  at x = 0.99 in half the library if the experimental lattice constant is
    #  used instead.
    win = np.arange(1.08, 0.919, -0.01)
    Ew = [E(x) for x in win]
    i0 = int(np.argmin(Ew))
    E0, x0 = Ew[i0], win[i0]
    #  walk inward and stop at whichever comes first - the barrier clearing the
    #  threshold, or the crystal turning downhill before it does
    for x in np.arange(x0 - 0.005, 0.549, -0.005):
        v = E(x)
        if v - E0 >= lim:
            return True
        if v < E0 - COMPRESSION_TOL:
            return False
    return True


#  Dynamical stability as a hard rejection.
#
#  Born stability only constrains the q -> 0 limit, so a potential can reproduce
#  every elastic constant and still have imaginary modes across most of the
#  Brillouin zone - which is what happened to Li (24.7 % imaginary, down to
#  -120 cm^-1), Na, Nb and Al.  Nothing in an elastic-only objective can see
#  that.  A coarse mesh inside the search is enough to reject those candidates;
#  the final answer is re-checked on a finer one.
#  Off for the angular runs: the screen needs force constants, and the
#  analytic derivatives of the angular term are not written yet.  The
#  elastic constants themselves are unaffected - for a one-atom cell they
#  come from strain derivatives of the energy alone.
REQUIRE_DYNAMICAL_STABILITY = True
#  Two-stage screen, and the confirm stage needs meshes of different
#  divisibility.  L.mesh is a half-shifted Monkhorst-Pack grid, so a mesh of n
#  samples (k + 1/2)/n and the sets for n = 4, 8, 16 are nearly nested - they
#  all miss the same places.  Measured minimum frequency in cm-1 for the Nb fit
#  that passed an 8^3 screen:
#
#        3^3     4^3     6^3     8^3     9^3    12^3
#      -44.5   +87.8   +45.8    +0.1   -51.9   -39.1
#
#  The soft mode sits at q = (-1/24, 11/24, -11/24), near the bcc N point, which
#  is where the well-known N-point and 2/3[111] anomalies of Nb, V, W and the
#  alkali metals live.  Every power-of-two mesh steps straight over it.  Al is
#  the mirror image: 9^3 sees -2.5 while 12^3 reports +8.6.  So no single mesh
#  is safe, and the confirm stage takes the union of an even and an odd one.
DYN_NQ = 4                    # 4^3 = 64 q-points, cheap reject
DYN_NQ_CONFIRM = (8, 9)       # 512 + 729, only for survivors
DYN_TOL = -1e-6               # frequencies in THz


#  A uniform mesh has no points closer to Gamma than 1/nq, and that is where
#  this potential fails when C' collapses.  The Al fit that passed 4^3 U 8^3 U
#  9^3 has C' = 0.14 GPa against 23.2 experimentally, which leaves the
#  transverse [110] branch flat enough at Gamma that quartic terms tip it under:
#  it reads +0.89 cm-1 at |q| = 0.02, -3.45 at 0.067, and +7.68 by 0.12 - a
#  narrow negative window that every mesh from 4^3 to 12^3 steps straight over.
#  So sample a Gamma neighbourhood explicitly, along the low-index directions
#  and a few generic ones, and do it in the cheap stage: it is 30 q-points.
_DIRS = np.array([[1., 0, 0], [1., 1, 0], [1., 1, 1], [2., 1, 0], [1., 0, 1],
                  [3., 1, 1], [1., 2, 3], [-1., 0, 1], [1., -1, 2], [2., 3, 1]])
_DIRS = _DIRS / np.linalg.norm(_DIRS, axis=1, keepdims=True)
NEAR_GAMMA = np.vstack([s * _DIRS for s in (0.03, 0.07, 0.13)])


#  Every q at once; one at a time this rebuilt the dynamical matrix from the
#  force-constant dict per point and was 91 % of a screened evaluation.
def _min_frequency(cry, pot, nq, Phi):
    return float(L.frequencies_many(cry, pot, L.mesh(nq), Phi).min())


def _min_near_gamma(cry, pot, Phi):
    return float(L.frequencies_many(cry, pot, NEAR_GAMMA, Phi).min())


def _phi(cry, pot):
    """
    Force constants, whichever implementation can do this potential.

    latdyn's are analytic but assume phi3 depends on r1 + r2 alone, so they
    refuse to run once the angular factor is on.  angfc differences each
    triplet's own 6x6 Hessian and places it by the chain rule; with
    lam2 = lam4 = 0 the two agree to about 1e-8 relative on all fourteen cubic
    metals, so this switch changes nothing for the published form.
    """
    if getattr(pot, "lam2", 0.0) or getattr(pot, "lam4", 0.0):
        import angfc
        return angfc.force_constants(cry, pot)
    return L.force_constants(cry, pot)


def dynamically_stable(cry, pot, nq=DYN_NQ, confirm=DYN_NQ_CONFIRM):
    try:
        Phi = _phi(cry, pot)
        if _min_near_gamma(cry, pot, Phi) < DYN_TOL:
            return False
        if _min_frequency(cry, pot, nq, Phi) < DYN_TOL:
            return False
        return all(_min_frequency(cry, pot, n, Phi) >= DYN_TOL
                   for n in confirm)
    except (OverflowError, ValueError, FloatingPointError,
            np.linalg.LinAlgError):
        return False


def born_stable(c, struct):
    try:
        if struct == "hcp":
            return (c["C11"] > abs(c["C12"]) and c["C44"] > 0 and
                    (c["C11"]+c["C12"])*c["C33"] > 2*c["C13"]**2)
        return (c["C11"] > abs(c["C12"]) and c["C44"] > 0 and
                c["C11"] + 2*c["C12"] > 0)
    except KeyError:
        return False


def warm_start(el):
    """
    Optional starting point, read from seed.json (plain numbers, no code).

    A continuous search must never end up worse than a restricted one, since the
    restricted parameter set is a subset.  Without a seed it did: the coarse grid
    missed Pd's optimum near m = 13, C = -40 and the simplex settled far away.
    Starting from a known-good point makes the improvement monotone by
    construction.  Delete seed.json to fit from scratch.
    """
    path = os.path.join(HERE, "seed.json")
    if not os.path.exists(path):
        return None
    v = json.load(open(path)).get(el)
    if not v or "alpha3" not in v:
        return None
    return (v["m"], v["gamma"], v["alpha3"]/v["alpha"], v["C"])


def optimise(el, e, verbose=False):
    best, seed = None, None
    #  gamma is CONTINUOUS here - that is the whole point
    grids = [(m, g, s3, C)
             for m in (1.4, 2.6, 4.5, 7.5, 11.0, 15.0, 19.0)
             for g in (0.0, 0.5, 1.0, 1.5, 2.0)
             for s3 in (0.5, 1.0, 2.0, 4.0, 6.0)
             for C in (0.0, -1.0, 1.0, -4.0, 4.0, -16.0, 16.0, -40.0, 40.0)]
    ws = warm_start(el)
    if ws:
        grids = [ws] + list(grids)
    for (m, g, s3, C) in grids:
        with np.errstate(over="ignore", invalid="ignore"):
            r = evaluate(el, e, m, g, C, s3, seed)
        if r is None:
            continue
        seed = r["seed"]
        if best is None or r["score"] < best["score"]:
            best = r
    if best is None:
        return None

    #  Nelder-Mead over the continuous parameters (m, gamma, s3, C).  Coordinate
    #  descent got stuck well short of the integer-gamma result, which cannot be
    #  right: integers are a subset of the reals, so a continuous search must do
    #  at least as well.  A simplex handles the strong m-gamma correlation that
    #  axis-by-axis stepping cannot.
    store = {}

    def f(v):
        m_, g_, s_, C_ = v
        if (g_ < 0 or not (M_MIN <= m_ <= M_MAX) or s_ <= 0.05
                or abs(C_) > C_MAX):
            return 1e6
        key = tuple(round(t, 7) for t in v)
        if key in store:
            return store[key][0]
        with np.errstate(over='ignore', invalid='ignore'):
            r = evaluate(el, e, m_, g_, C_, s_, best['seed'])
        val = r["score"] if r else 1e6
        store[key] = (val, r)
        return val

    x0 = np.array([best["m"], best["gamma"], best["s3"], best["C"]])
    step = np.array([max(0.6, 0.15*x0[0]), 0.35, max(0.3, 0.3*x0[2]),
                     max(1.0, 0.4*abs(x0[3]) + 0.5)])
    sim = [x0] + [x0 + np.eye(4)[i]*step[i] for i in range(4)]
    val = [f(s) for s in sim]
    for _ in range(260):
        order = np.argsort(val)
        sim = [sim[i] for i in order]; val = [val[i] for i in order]
        if abs(val[-1] - val[0]) < 1e-9 * max(1.0, abs(val[0])):
            break
        cen = np.mean(sim[:-1], axis=0)
        xr = cen + (cen - sim[-1]); fr = f(xr)
        if fr < val[0]:
            xe = cen + 2.0*(cen - sim[-1]); fe = f(xe)
            sim[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            sim[-1], val[-1] = xr, fr
        else:
            xc = cen + 0.5*(sim[-1] - cen); fc = f(xc)
            if fc < val[-1]:
                sim[-1], val[-1] = xc, fc
            else:
                for i in range(1, 5):
                    sim[i] = sim[0] + 0.5*(sim[i] - sim[0])
                    val[i] = f(sim[i])
    for (v, r) in store.values():
        if r and v < best["score"]:
            best = r
    return best


def main(els):
    #  merge, never replace: running on a subset used to overwrite fit.json with
    #  only those elements and silently drop every other result.
    path = os.path.join(HERE, "fit.json")
    out = {}
    if os.path.exists(path):
        try:
            out = json.load(open(path))
        except ValueError:
            out = {}
    t0 = time.time()
    for i, el in enumerate(els, 1):
        e = refdata.ELEMENTS[el]
        r = optimise(el, e)
        if r is None:
            print(f"  [{i:2d}/{len(els)}] {el:3s} FAILED", flush=True)
            continue
        r.pop("seed", None)
        r["struct"] = e["struct"]
        tol = 1e-3
        r["at_bound"] = [nm for nm, v_, lo, hi in
                         (("m", r["m"], M_MIN, M_MAX),
                          ("C", abs(r["C"]), -1.0, C_MAX),
                          ("gamma", r["gamma"], 0.0, 1e9))
                         if abs(v_-lo) < tol or abs(v_-hi) < tol]
        out[el] = r
        print(f"  [{i:2d}/{len(els)}] {el:3s} {e['struct']}  rms={r['score']*100:5.1f}%"
              f"  gamma={r['gamma']:6.3f}  m={r['m']:6.2f}  C={r['C']:+8.3f}"
              f"{' BOUND:'+','.join(r['at_bound']) if r['at_bound'] else '':>14s}"
              f"  ({time.time()-t0:6.1f}s)", flush=True)
    json.dump(out, open(path, "w"), indent=1, sort_keys=True)
    print(f"wrote fit.json ({len(out)} elements)")


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(refdata.ELEMENTS))
