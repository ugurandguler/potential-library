#!/usr/bin/env python3
"""
Wakabayashi, Scherm and Smith's dipolar fluctuation model for hcp metals.

N. Wakabayashi, R. H. Scherm and H. G. Smith, Phys. Rev. B 25, 5122 (1982),
doi:10.1103/PhysRevB.25.5122.  They measured Ti, Co and Tc, re-analysed Sc,
Y, Zr and Hf, and published the FITTED PARAMETERS rather than the
frequencies - which is why cobalt has no transcribable table anywhere and has
to be reconstructed if it is to be on the page at all.

Their Eq. (3):

    D(q) = P(q) - T*(q) [K + S(q)]^-1 T(q)

P is the ion-ion interaction, T the ion-to-dipole coupling, S the
dipole-dipole coupling and K the self-energy of the dipole.  The paper says
"the parameters in T and S are determined only within the scaling factor of
1/sqrt(K) and 1/K", so K divides out and Table V already carries the scaled
quantities:

    D(q) = P(q) - T*(q) [1 + S(q)]^-1 T(q)

with T in dyn^1/2 cm^-1/2, so T^2 is in dyn/cm like P, and S dimensionless.

Each of P, T and S is the SAME lattice sum with its own constants - Eq. (4),
written for L[001], has one Fourier form and three sets of coefficients.  P
runs to the third neighbour, T and S only to the second.

Table II fixes the tensors, and its "typical atom" coordinates fix the frame:
neighbour 1 at (a/sqrt3, 0, c/2), 2 at (0, a, 0), 3 at (-2a/sqrt3, 0, c/2),
4 at (0, 0, c).  Choosing a1 = a(sqrt3/2, -1/2, 0), a2 = a(sqrt3/2, 1/2, 0),
a3 = c z and the second atom at (a1 + a2)/3 + a3/2 puts all four exactly
there, which the module checks rather than assumes.
"""
import math

import numpy as np

DYNCM_TO_NM = 1.0e-3          # 1 N/m = 10^3 dyn/cm
AMU = 1.66053906660e-27
ANG = 1.0e-10


def lattice(a, c):
    """(a1, a2, a3, basis) in Angstrom, in Wakabayashi's frame"""
    a1 = a * np.array([math.sqrt(3) / 2.0, -0.5, 0.0])
    a2 = a * np.array([math.sqrt(3) / 2.0, 0.5, 0.0])
    a3 = np.array([0.0, 0.0, c])
    return a1, a2, a3, (a1 + a2) / 3.0 + a3 / 2.0


def site_ops():
    """the twelve operations of D3h, as Cartesian 3x3 matrices

    Not D6h: an atom in hcp sits on a three-fold axis, not a six-fold one -
    the six-fold is a screw.  Twelve is what generates each shell exactly
    once, and the check below fails loudly if a shell comes out the wrong
    size.
    """
    t = 2.0 * math.pi / 3.0
    c3 = np.array([[math.cos(t), -math.sin(t), 0.0],
                   [math.sin(t), math.cos(t), 0.0],
                   [0.0, 0.0, 1.0]])
    my = np.diag([1.0, -1.0, 1.0])
    mz = np.diag([1.0, 1.0, -1.0])
    out = []
    R = np.eye(3)
    for _ in range(3):
        for A in (np.eye(3), my):
            for B in (np.eye(3), mz):
                out.append(B @ A @ R)
        R = c3 @ R
    return out


#  Table II.  `inter` says whether the neighbour is on the other sublattice.
TYPICAL = {
    1: dict(inter=True, n=6, r=lambda a, c: (a / math.sqrt(3), 0.0, c / 2.0),
            phi=lambda p: np.array([[p["a"], 0.0, p.get("d", 0.0)],
                                    [0.0, p["b"], 0.0],
                                    [p.get("d", 0.0), 0.0, p["g"]]])),
    2: dict(inter=False, n=6, r=lambda a, c: (0.0, a, 0.0),
            phi=lambda p: np.array([[p["a"], p.get("e", 0.0), 0.0],
                                    [-p.get("e", 0.0), p["b"], 0.0],
                                    [0.0, 0.0, p["g"]]])),
    3: dict(inter=True, n=6,
            r=lambda a, c: (-2.0 * a / math.sqrt(3), 0.0, c / 2.0),
            phi=lambda p: np.array([[p["a"], 0.0, p.get("d", 0.0)],
                                    [0.0, p["b"], 0.0],
                                    [p.get("d", 0.0), 0.0, p["g"]]])),
    4: dict(inter=False, n=2, r=lambda a, c: (0.0, 0.0, c),
            phi=lambda p: np.diag([p["a"], p["a"], p["g"]])),
}


def shell(n, comps, a, c):
    """[(R, Phi, inter)] for one shell, every member of it

    Each member is reached from the typical atom by a site operation S, and
    its tensor is S Phi0 S^T.  Where two operations give the same R they must
    give the same tensor - that is the statement that Phi0 has the site
    symmetry it is supposed to have, and it is asserted.
    """
    spec = TYPICAL[n]
    r0 = np.array(spec["r"](a, c), float)
    phi0 = spec["phi"](comps)
    out = {}
    for S in site_ops():
        v = S @ r0
        key = tuple(np.round(v, 8) + 0.0)
        p = S @ phi0 @ S.T
        if key in out:
            assert np.abs(out[key] - p).max() < 1e-8 * max(
                1.0, np.abs(p).max()), f"shell {n}: inconsistent tensor"
        out[key] = p
    assert len(out) == spec["n"], (
        f"shell {n}: {len(out)} members, expected {spec['n']}")
    return [(np.array(k), v, spec["inter"]) for k, v in out.items()]


def pairs(shells, a, c):
    """[(s, s2, R, Phi)] for both sublattices, R in Angstrom

    The second sublattice is not assumed to repeat the first.  Its
    inter-sublattice neighbours are the reversed vectors with the TRANSPOSED
    tensor, which is the general Phi(ls, l's') = Phi(l's', ls)^T, and its
    self-term is built from its own list.
    """
    out = []
    for n, comps in shells.items():
        for R, Phi, inter in shell(n, comps, a, c):
            if inter:
                out.append((0, 1, R, Phi))
                out.append((1, 0, -R, Phi.T))
            else:
                out.append((0, 0, R, Phi))
                out.append((1, 1, R, Phi))
    return out


def block(q, shells, a, c):
    """the 6x6 lattice sum at q, q in CARTESIAN 2 pi / Angstrom

    Includes the self-term, so the sum vanishes at q = 0 - which is what
    makes the acoustic branches go to zero and what makes Eq. (4) read
    P1 (1 - cos) rather than P1 cos.
    """
    #  SIGN.  Wakabayashi's Eq. (4) reads nu^2 = P1 (1 - cos) + ..., so his
    #  constants are spring constants and the sum is Phi (1 - exp), not
    #  Phi (exp - 1).  Getting it the other way round still gives three
    #  acoustic branches that vanish at Gamma - the sum rule does not care -
    #  and turns every optic branch imaginary, which is how it was caught.
    P = pairs(shells, a, c)
    M = np.zeros((6, 6), complex)
    for s, s2, R, Phi in P:
        M[3*s:3*s+3, 3*s2:3*s2+3] -= Phi * np.exp(1j * float(np.dot(q, R)))
    for s in (0, 1):
        M[3*s:3*s+3, 3*s:3*s+3] += sum(Phi for t, _, _, Phi in P if t == s)
    return M


def dyn(q, P, T, S, a, c):
    """Eq. (3) with K divided out, in dyn/cm"""
    M = block(q, P, a, c)
    if T:
        Tq = block(q, T, a, c)
        Sq = block(q, S, a, c) if S else np.zeros((6, 6), complex)
        M = M - Tq.conj().T @ np.linalg.solve(np.eye(6) + Sq, Tq)
    return 0.5 * (M + M.conj().T)


def freqs_THz(q, P, T, S, a, c, mass_amu):
    w2 = np.linalg.eigvalsh(dyn(q, P, T, S, a, c)
                            * DYNCM_TO_NM / (mass_amu * AMU))
    w = np.sign(w2) * np.sqrt(np.abs(w2))
    return np.sort(w / (2.0 * math.pi) / 1.0e12)


#  Table V.  P in dyn/cm, T in dyn^1/2 cm^-1/2, S dimensionless.  Scandium
#  and yttrium carry footnote a, "Born-von Karman parameters only", so they
#  have no dipolar term at all - which is what makes them a clean test of the
#  ion-ion half on its own.
#
#  Two readings of the printed table are worth stating.  The sixth row of the
#  T block is printed "gamma_1(T)" a second time; it is gamma_2(T), by its
#  position in the alpha_2 beta_2 gamma_2 epsilon_2 sequence.  And delta_1
#  and delta_3, the xz couplings that Table II allows for the out-of-plane
#  shells, are not listed for any element, so they are zero here.
WAK_P = {
    "Sc": {1: {"a": 9281.1, "b": 777.7, "g": 12338.4},
           2: {"a": 3313.0, "b": 15414.2, "g": 1980.4},
           3: {"a": 1302.0, "b": -3377.2, "g": -1297.1}},
    "Y":  {1: {"a": 8355.8, "b": 0.3, "g": 11188.7},
           2: {"a": 2824.3, "b": 12691.5, "g": 900.2},
           3: {"a": -9.9, "b": -1236.0, "g": -921.0}},
    "Ti": {1: {"a": 16733.6, "b": -1541.1, "g": 25375.9},
           2: {"a": 1715.5, "b": 37230.5, "g": 2913.6, "e": 10126.2},
           3: {"a": -3192.9, "b": -2246.2, "g": -3608.5}},
    "Zr": {1: {"a": 12832.2, "b": 248.7, "g": 24098.1},
           2: {"a": 2177.5, "b": 37816.8, "g": 1747.0, "e": 12091.9},
           3: {"a": -4105.6, "b": -2015.6, "g": -3563.4}},
    "Hf": {1: {"a": 12996.0, "b": 5179.9, "g": 28736.3},
           2: {"a": 390.7, "b": 37860.6, "g": 6949.5, "e": -1206.8},
           3: {"a": -3426.5, "b": 239.4, "g": -5114.5}},
    "Co": {1: {"a": 11925.5, "b": -1267.4, "g": 36831.4},
           2: {"a": -1819.8, "b": 40089.8, "g": 2171.9, "e": 1502.5},
           3: {"a": 779.8, "b": 455.2, "g": -3507.4}},
}
WAK_T = {
    "Ti": {1: {"a": 13.47, "b": -6.57, "g": 23.22},
           2: {"a": 6.71, "b": 45.50, "g": 9.83, "e": 15.50}},
    "Zr": {1: {"a": 5.64, "b": 3.51, "g": 22.62},
           2: {"a": 7.57, "b": 43.70, "g": 6.44, "e": 23.00}},
    "Hf": {1: {"a": -1.44, "b": 17.85, "g": 20.74},
           2: {"a": 4.70, "b": 36.46, "g": 13.89, "e": -10.85}},
    #  cobalt's column has only the zz entries filled: the paper set the
    #  in-plane dipolar parameters to zero because "no obvious anomalies were
    #  detected in the branches with the polarization vectors parallel to the
    #  hexagonal plane"
    "Co": {1: {"a": 0.0, "b": 0.0, "g": 21.70},
           2: {"a": 0.0, "b": 0.0, "g": 9.79}},
}
WAK_S = {
    "Ti": {1: {"a": 0.0, "b": 0.0, "g": -0.0422},
           2: {"a": 0.0, "b": 0.0, "g": 0.0029}},
    "Zr": {1: {"a": 0.0, "b": 0.0, "g": -0.0325}},
    "Hf": {1: {"a": 0.0, "b": 0.0, "g": -0.0528}},
    "Co": {1: {"a": 0.0, "b": 0.0, "g": -0.0534}},
}
CO_REF = ("N. Wakabayashi, R. H. Scherm and H. G. Smith, Phys. Rev. B "
          "<b>25</b>, 5122 (1982), Table V; the dipolar fluctuation model "
          "fitted to their own neutron measurements, made at room "
          "temperature under the conditions they quote as 295 K for titanium")

#  hcp symmetry points in the fractional reciprocal basis of lattice()
HCP_PTS = {"G": (0, 0, 0), "M": (0.5, 0, 0), "K": (1/3.0, -1/3.0, 0),
           "A": (0, 0, 0.5), "L": (0.5, 0, 0.5), "H": (1/3.0, -1/3.0, 0.5)}


def reciprocal(a, c):
    """rows b1, b2, b3 in 2 pi / Angstrom"""
    a1, a2, a3, _ = lattice(a, c)
    V = float(np.dot(a1, np.cross(a2, a3)))
    return np.array([np.cross(a2, a3), np.cross(a3, a1),
                     np.cross(a1, a2)]) * 2.0 * math.pi / V


def check_ladder(verbose=True):
    """the model against 464 measured frequencies it never saw

    Five of the seven elements in Table V already carry a measured curve on
    this page, from other people's papers, so the reconstruction can be
    tested before it is trusted for cobalt, which carries none.  Scandium and
    yttrium test the ion-ion half alone; titanium, zirconium and hafnium test
    what the dipolar term adds.
    """
    import refdata
    import refdata_phonon_curves as C
    rows = []
    for el in ("Sc", "Y", "Ti", "Zr", "Hf"):
        rec = C.PHONON_CURVE.get(el)
        if not rec or "branches" not in rec:
            continue
        a = refdata.ELEMENTS[el]["a0"]
        c = a * refdata.ELEMENTS[el]["c_over_a"]
        m = refdata.MASSES[el]
        B = reciprocal(a, c)
        ov = rec.get("seg_override") or {}
        ns = rec.get("nu_scale", 1.0)
        axis = {"[00z]": B[2], "[z00]": B[0], "[zz0]": B[0] - B[1]}
        err_p, err_d, top = [], [], 0.0
        for lab, pts in rec["branches"].items():
            br = lab[lab.find("["):]
            if br not in axis:
                continue
            base = C.SEGMENT["hcp"][br][0][3]
            zs = base / ov.get(br, C.SEGMENT["hcp"][br])[0][3]
            for row in pts:
                q = row[0] * zs * axis[br]
                nu = row[1] * ns
                top = max(top, nu)
                err_p.append(min(abs(freqs_THz(q, WAK_P[el], None, None,
                                               a, c, m) - nu)))
                err_d.append(min(abs(freqs_THz(q, WAK_P[el], WAK_T.get(el),
                                               WAK_S.get(el), a, c, m) - nu)))
        rows.append((el, len(err_d), float(np.sqrt(np.mean(np.square(err_p)))),
                     float(np.sqrt(np.mean(np.square(err_d)))),
                     100.0 * float(np.mean(err_d)) / top, el in WAK_T))
    if verbose:
        print("  Wakabayashi's model against measured curves it never saw:")
        print(f"    {'el':4s}{'points':>8s}{'P only':>12s}{'full DFM':>12s}"
              f"{'  of top branch':>16s}")
        for el, n, rp, rd, pct, has in rows:
            tag = "" if has else "   (no dipolar term published)"
            print(f"    {el:4s}{n:8d}{rp:10.3f} THz{rd:9.3f} THz"
                  f"{pct:13.1f}%{tag}")
    ok = all(pct < 3.0 for _, _, _, _, pct, _ in rows)
    if verbose and not ok:
        print("    FAILED")
    return ok
