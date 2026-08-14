#!/usr/bin/env python3
"""
Lattice dynamics for the Akgun-Ugur potential.

    phi2(r)     = D/(m-1) (r0/r)^g [ e^{m a (r0-r)} - m e^{a (r0-r)} ]
    phi3(r1,r2) = C D/(m-1) (r0/x)^g [ e^{m a3 (r0-x)} - m e^{a3 (r0-x)} ],  x = r1+r2

Everything here is a derivative of that expression, so gamma stays a real
number - the integer restriction came from an input format, not physics.

WHAT IS COMPUTED
  * analytic first and second Cartesian derivatives of the energy
  * force-constant matrices  Phi(0k, lk')  ->  dynamical matrix  D(q)
  * phonon frequencies on a path or a Monkhorst-Pack mesh
  * elastic constants with the non-affine (internal relaxation) correction
  * harmonic thermodynamics: zero-point energy, F, S, Cv

ELASTIC CONSTANTS
Following Born & Huang and, for many-body potentials with internal relaxation,
Grieser/Pastewka, Phys. Rev. Materials 7, 073603 (2023) [arXiv:2302.08754]:

    C_abgd = C^Born_abgd - (1/V) Xi_{i,ab} (H^+)_{ij} Xi_{j,gd}

    C^Born  : second derivative of the energy under affine strain, ions frozen
    H       : Gamma-point force-constant (Hessian) matrix, 3N x 3N
    Xi      : d(force on atom i) / d(strain), the strain-displacement coupling
    H^+     : Moore-Penrose pseudo-inverse; the three translations are null
              modes and must be projected out, not inverted

For one atom per primitive cell (fcc, bcc) symmetry forces Xi = 0, so the Born
term is already exact.  For hcp it is not: that correction is the difference
between frozen-ion and relaxed elastic constants.

THREE-BODY DERIVATIVES
phi3 is a function of x = r_ij + r_ik alone (i is the vertex), so with
n_j, n_k the unit vectors from i to j and k:

    dx/dr_j = n_j,   dx/dr_k = n_k,   dx/dr_i = -(n_j + n_k)

    d2/dr_j dr_j = g'' n_j n_j + (g'/r_ij)(I - n_j n_j)
    d2/dr_j dr_k = g'' n_j n_k                 <- the bond-bond cross term,
                                                  this is what breaks Cauchy
Units: eV, Angstrom, amu.  Frequencies in THz and cm^-1.
"""
import math
import numpy as np

EV_A3_TO_GPA = 160.21766208
#  sqrt(eV / (amu A^2)) -> rad/s ;  then /2pi -> Hz
_W = math.sqrt(1.602176634e-19 / (1.66053906660e-27 * 1e-20))
THZ = _W / (2.0 * math.pi) / 1e12          # eigenvalue^(1/2) -> THz
CM1 = 33.35641                             # THz -> cm^-1
KB = 8.617333262e-5                        # eV/K
HBAR_THZ = 4.135667696e-3 / (2*math.pi) * 2*math.pi   # h in eV/THz = 4.1357e-3


# --------------------------------------------------------------------------
#  potential and its radial derivatives
# --------------------------------------------------------------------------
class Potential:
    """phi2 and phi3 share (m, r0, gamma); alpha3 sets the three-body range."""

    def __init__(self, m, D, alpha, r0, gamma, C=0.0, alpha3=None,
                 rcut2=None, rcut3=None, lam2=0.0, lam4=0.0, taper=None):
        self.m, self.D, self.alpha, self.r0, self.gamma = m, D, alpha, r0, gamma
        self.C = C
        self.alpha3 = alpha if alpha3 is None else alpha3
        self.rcut2, self.rcut3 = rcut2, rcut3
        #  Cutoff switch; see standalone/latdyn.Potential for why it exists and
        #  what it costs.  None reproduces hard truncation exactly, and is the
        #  default so that the fits made before it still evaluate to their own
        #  numbers.
        self.taper = taper
        #  Angular extension.  phi3 in the published form depends only on
        #  x = r1 + r2, so it cannot see a distortion that changes bond angles
        #  at fixed bond lengths - which is exactly what C' measures, and why
        #  C44/C' has a floor the form cannot go under.  Multiplying by a
        #  Legendre factor fixes that with one parameter:
        #
        #      phi3(r1, r2, theta) = phi3_radial(r1 + r2) * [1 + lam2 P2(cos)]
        #
        #  P2 averages to zero over the sphere, so the isotropic part - and
        #  with it the cohesive energy, the pressure and the bulk modulus that
        #  the fit holds exactly - is left where it was; only the anisotropy
        #  moves.  lam2 = 0 reproduces the published potential exactly, so this
        #  is a strict generalisation and the old library is a special case.
        self.lam2 = lam2
        #  A second Legendre order costs one more parameter and, like P2,
        #  averages to zero over the sphere, so it too leaves the cohesive
        #  energy and the bulk modulus constraints undisturbed.
        self.lam4 = lam4

    def ang(self, c):
        """angular factor at cos(theta); 1 when lam2 = lam4 = 0"""
        if not (self.lam2 or self.lam4):
            return np.ones_like(np.asarray(c, dtype=float))
        c = np.asarray(c, dtype=float)
        c2 = c * c
        h = 1.0 + self.lam2 * 0.5 * (3.0 * c2 - 1.0)
        if self.lam4:
            h = h + self.lam4 * 0.125 * (35.0 * c2 * c2 - 30.0 * c2 + 3.0)
        return h

    def _f(self, r, al, order):
        """value / d/dr / d2dr2 of (r0/r)^g [e^{m al (r0-r)} - m e^{al (r0-r)}]/(m-1)"""
        m, r0, g = self.m, self.r0, self.gamma
        u = al * (r0 - r)
        e1, em = np.exp(u), np.exp(m * u)
        A = (em - m * e1) / (m - 1.0)                     # bracket
        dA = (-m*al*em + m*al*e1) / (m - 1.0)             # d bracket / dr
        d2A = (m*m*al*al*em - m*al*al*e1) / (m - 1.0)
        P = (r0 / r) ** g                                 # prefactor
        dP = -g * P / r
        d2P = g * (g + 1.0) * P / (r * r)
        if order == 0:
            return P * A
        if order == 1:
            return P * dA + dP * A
        return P * d2A + 2.0 * dP * dA + d2P * A

    @classmethod
    def from_record(cls, rec):
        """rebuild the potential a stored fit was actually made with"""
        return cls(rec["m"], rec["D"], rec["alpha"], rec["r0"], rec["gamma"],
                   C=rec.get("C", 0.0), alpha3=rec.get("alpha3"),
                   rcut2=rec.get("rcut2"), rcut3=rec.get("rcut3"),
                   lam2=rec.get("lam2", 0.0), lam4=rec.get("lam4", 0.0),
                   taper=rec.get("taper"))

    def switch(self, r, rc, order=0):
        """quintic S: 1 below taper*rc, 0 at rc, S' = S'' = 0 at both ends"""
        r = np.asarray(r, dtype=float)
        if self.taper is None or not rc:
            z = np.zeros_like(r)
            return np.ones_like(r) if order == 0 else z
        r_on = self.taper * rc
        w = rc - r_on
        t = np.clip((r - r_on) / w, 0.0, 1.0)
        if order == 0:
            return 1.0 - t**3 * (10.0 - 15.0 * t + 6.0 * t * t)
        inside = (r > r_on) & (r < rc)
        if order == 1:
            return np.where(inside, -30.0 * t * t * (1.0 - t) ** 2 / w, 0.0)
        return np.where(inside,
                        -60.0 * t * (1.0 - 3.0 * t + 2.0 * t * t) / (w * w), 0.0)

    def phi2(self, r, order=0):
        v = self.D * self._f(r, self.alpha, order)
        if self.taper is None or not self.rcut2:
            return v
        S = self.switch(r, self.rcut2, 0)
        if order == 0:
            return v * S
        f0 = self.D * self._f(r, self.alpha, 0)
        S1 = self.switch(r, self.rcut2, 1)
        if order == 1:
            return v * S + f0 * S1
        f1 = self.D * self._f(r, self.alpha, 1)
        S2 = self.switch(r, self.rcut2, 2)
        return v * S + 2.0 * f1 * S1 + f0 * S2

    def phi3(self, x, order=0):
        if self.C == 0.0:
            return np.zeros_like(np.asarray(x, dtype=float))
        return self.C * self.D * self._f(x, self.alpha3, order)


# --------------------------------------------------------------------------
#  crystal
# --------------------------------------------------------------------------
class Crystal:
    """primitive cell: lattice (3x3 rows), fractional basis, masses"""

    STRUCT = {
        "fcc": (np.array([[0., .5, .5], [.5, 0., .5], [.5, .5, 0.]]),
                np.array([[0., 0., 0.]])),
        "bcc": (np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]),
                np.array([[0., 0., 0.]])),
    }

    def __init__(self, struct, a0, c_over_a=None, mass=1.0):
        self.struct, self.a0, self.c_over_a = struct, a0, c_over_a
        if struct == "hcp":
            s3 = math.sqrt(3.0) / 2.0
            self.lat = a0 * np.array([[1., 0., 0.], [-.5, s3, 0.],
                                      [0., 0., c_over_a]])
            self.frac = np.array([[1/3., 2/3., .25], [2/3., 1/3., .75]])
        else:
            L, F = self.STRUCT[struct]
            self.lat, self.frac = a0 * L.copy(), F.copy()
        self.mass = np.full(len(self.frac), mass)

    @property
    def pos(self):
        return self.frac @ self.lat

    @property
    def vol(self):
        return abs(np.linalg.det(self.lat))

    def strained(self, eps):
        """new Crystal with lattice (I+eps) applied; fractional coords kept"""
        c = object.__new__(Crystal)
        c.struct, c.a0, c.c_over_a = self.struct, self.a0, self.c_over_a
        c.lat = self.lat @ (np.eye(3) + eps).T
        c.frac, c.mass = self.frac.copy(), self.mass.copy()
        return c


#  The neighbour list depends only on the geometry, never on the potential
#  parameters.  During a fit the geometry is fixed apart from a handful of small
#  strains, so without this cache the list was being rebuilt hundreds of
#  thousands of times and dominated the entire run time.
_NB_CACHE = {}
_NB_MAX = 512


def _nb_key(cry, rcut):
    return (cry.lat.tobytes(), cry.frac.tobytes(), round(float(rcut), 10))


def neighbours(cry, rcut):
    """
    All (i, j, R, dvec, r) with |dvec| <= rcut, where dvec = r_j + R - r_i.
    R is the lattice translation of the image of j.  Cached.
    """
    key = _nb_key(cry, rcut)
    hit = _NB_CACHE.get(key)
    if hit is not None:
        return hit
    out = _neighbours_build(cry, rcut)
    if len(_NB_CACHE) > _NB_MAX:
        _NB_CACHE.clear()
    _NB_CACHE[key] = out
    return out


def _neighbours_build(cry, rcut):
    lat, pos = cry.lat, cry.pos
    n = len(pos)
    #  how many cells to search in each direction
    inv = np.linalg.inv(lat)
    reach = [int(math.ceil(rcut * np.linalg.norm(inv[:, k]))) + 1
             for k in range(3)]
    out = []
    for i1 in range(-reach[0], reach[0] + 1):
        for i2 in range(-reach[1], reach[1] + 1):
            for i3 in range(-reach[2], reach[2] + 1):
                R = np.array([i1, i2, i3]) @ lat
                for i in range(n):
                    for j in range(n):
                        d = pos[j] + R - pos[i]
                        r = math.sqrt(d @ d)
                        if 1e-8 < r <= rcut:
                            out.append((i, j, np.array([i1, i2, i3]), d, r))
    return out


# --------------------------------------------------------------------------
#  energy, forces, force constants
# --------------------------------------------------------------------------
_DIST_CACHE = {}


def _distances(cry, pot):
    """(pair distances, triplet sums) as arrays - cached alongside the list"""
    key = (_nb_key(cry, pot.rcut2), round(float(pot.rcut3 or 0.0), 10))
    hit = _DIST_CACHE.get(key)
    if hit is not None:
        return hit
    rs = np.array([r for (_, _, _, _, r) in neighbours(cry, pot.rcut2)])
    if pot.rcut3:
        tr = list(triplets(cry, pot))
        xs = np.array([ra + rb for (_, _, _, ra, rb, _, _) in tr])
        cs = np.array([float(na @ nb) for (_, _, _, _, _, na, nb) in tr])
        #  the legs, not only their sum: the three-body cutoff applies to each
        #  leg, so the switch multiplies per leg and needs them separately
        legs = np.array([[ra, rb] for (_, _, _, ra, rb, _, _) in tr],
                        dtype=float)
        if not len(legs):
            legs = np.empty((0, 2))
    else:
        xs = cs = np.empty(0)
        legs = np.empty((0, 2))
    if len(_DIST_CACHE) > _NB_MAX:
        _DIST_CACHE.clear()
    _DIST_CACHE[key] = (rs, xs, cs, legs)
    return rs, xs, cs, legs


def energy(cry, pot):
    """energy per atom (eV)"""
    rs, xs, cs, legs = _distances(cry, pot)
    e = 0.5 * pot.phi2(rs).sum()
    if pot.C and len(xs):
        g = pot.phi3(xs) * pot.ang(cs)
        if pot.taper is not None:
            #  E3 = phi3(ra + rb) h(cos) S(ra) S(rb); the switch is per leg
            #  because the cutoff is, so it cannot be folded into phi3
            g = g * (pot.switch(legs[:, 0], pot.rcut3, 0)
                     * pot.switch(legs[:, 1], pot.rcut3, 0))
        e += g.sum()
    return e / len(cry.frac)


def triplets(cry, pot):
    """
    (i, ja, jb, ra, rb, na, nb) - vertex i with two neighbours inside rcut3,
    unordered pairs (each pair once).  ja/jb carry (atom index, image vector).
    """
    nb = [x for x in neighbours(cry, pot.rcut3)]
    by_i = {}
    for (i, j, R, d, r) in nb:
        by_i.setdefault(i, []).append((j, R, d, r))
    out = []
    for i, lst in by_i.items():
        for a in range(len(lst)):
            ja, Ra, da, ra = lst[a]
            for b in range(a + 1, len(lst)):
                jb, Rb, db, rb = lst[b]
                out.append((i, (ja, Ra), (jb, Rb), ra, rb, da / ra, db / rb))
    return out


def _pair_block(r, n, f1, f2):
    """d2 phi / d r_j d r_j for a pair with unit vector n, given phi' and phi''"""
    nn = np.outer(n, n)
    return f2 * nn + (f1 / r) * (np.eye(3) - nn)


def force_constants(cry, pot):
    """
    Phi[(i, j, R)] = d2 E / d u_i d u_j   (3x3 blocks), including on-site terms
    at R = 0, i == j.  Returned as a dict keyed by (i, j, R-tuple).
    """
    #  The blocks below assume phi3 depends on r1 + r2 alone, so with an angular
    #  factor every d/dtheta term would be missing.  angfc differences each
    #  triplet's own 6x6 Hessian instead; with lam = 0 the two agree to ~1e-8
    #  relative, so this hands off only where it has to.  Imported here rather
    #  than at the top because angfc imports this module.
    if getattr(pot, "lam2", 0.0) or getattr(pot, "lam4", 0.0):
        import angfc
        return angfc.force_constants(cry, pot)

    Phi = {}

    def add(i, j, R, blk):
        k = (i, j, tuple(int(x) for x in R))
        Phi[k] = Phi.get(k, np.zeros((3, 3))) + blk

    #  ---- two-body ----
    for (i, j, R, d, r) in neighbours(cry, pot.rcut2):
        n = d / r
        f1, f2 = pot.phi2(r, 1), pot.phi2(r, 2)
        B = _pair_block(r, n, f1, f2)
        #  The neighbour list holds both (i,j,R) and (j,i,-R), so this loop sees
        #  every bond twice - but the two visits fill DIFFERENT keys, Phi(i,j,R)
        #  and Phi(j,i,-R).  Each key must therefore receive the full -B, not
        #  half of it.  (Halving here made every force constant a factor of two
        #  small; the elastic constants still came out right because they are
        #  computed from strain derivatives of the energy, not from Phi.)
        add(i, j, R, -B)
        add(i, i, (0, 0, 0), B)

    #  ---- three-body: phi3 = g(x), x = r_ij + r_ik ----
    if pot.C:
        for (i, (ja, Ra), (jb, Rb), ra, rb, na, nb_) in triplets(cry, pot):
            x = ra + rb
            g1, g2 = pot.phi3(x, 1), pot.phi3(x, 2)
            naa, nbb = np.outer(na, na), np.outer(nb_, nb_)
            Pa = (np.eye(3) - naa) / ra
            Pb = (np.eye(3) - nbb) / rb
            #  blocks in the (a, b, i) basis
            Baa = g2 * naa + g1 * Pa
            Bbb = g2 * nbb + g1 * Pb
            Bab = g2 * np.outer(na, nb_)
            #  a-a, b-b, a-b
            add(ja, ja, (0, 0, 0), Baa); add(i, i, (0, 0, 0), Baa)
            add(i, ja, Ra, -Baa);        add(ja, i, -Ra, -Baa)
            add(jb, jb, (0, 0, 0), Bbb); add(i, i, (0, 0, 0), Bbb)
            add(i, jb, Rb, -Bbb);        add(jb, i, -Rb, -Bbb)
            add(ja, jb, Rb - Ra, Bab);   add(jb, ja, Ra - Rb, Bab.T)
            add(i, ja, Ra, -Bab.T);      add(ja, i, -Ra, -Bab)
            add(i, jb, Rb, -Bab);        add(jb, i, -Rb, -Bab.T)
            add(i, i, (0, 0, 0), Bab + Bab.T)
    return Phi


def check_force_constants(cry, pot, tol=1e-7):
    """
    Necessary conditions on Phi.  Both are exactly the checks that would have
    caught the factor-of-two pair bug immediately:

      acoustic sum rule   sum_j Phi(i,j) = 0   (rigid translation costs nothing)
      hermiticity         Phi(i,j,R) = Phi(j,i,-R)^T
    """
    Phi = force_constants(cry, pot)
    n = len(cry.frac)
    msgs = []
    for i in range(n):
        S = np.zeros((3, 3))
        for (a, b, R), blk in Phi.items():
            if a == i:
                S += blk
        if np.abs(S).max() > tol:
            msgs.append(f"acoustic sum rule violated for atom {i}: "
                        f"max |sum_j Phi| = {np.abs(S).max():.3e}")
    for (i, j, R), blk in Phi.items():
        other = Phi.get((j, i, tuple(-np.array(R))))
        if other is None:
            msgs.append(f"missing transpose partner for ({i},{j},{R})")
        elif np.abs(blk - other.T).max() > tol:
            msgs.append(f"Phi({i},{j},{R}) != Phi({j},{i},-R)^T, "
                        f"max diff {np.abs(blk-other.T).max():.3e}")
    return msgs


def dynamical(cry, pot, q, Phi=None):
    """D(q) for q in fractional reciprocal coordinates; returns 3N x 3N complex"""
    Phi = force_constants(cry, pot) if Phi is None else Phi
    n = len(cry.frac)
    D = np.zeros((3*n, 3*n), dtype=complex)
    for (i, j, R), blk in Phi.items():
        ph = np.exp(2j * math.pi * np.dot(q, R))
        D[3*i:3*i+3, 3*j:3*j+3] += blk * ph
    M = np.repeat(cry.mass, 3)
    D /= np.sqrt(np.outer(M, M))
    return 0.5 * (D + D.conj().T)


def frequencies(cry, pot, q, Phi=None):
    """phonon frequencies in THz (negative = imaginary mode)"""
    w2 = np.linalg.eigvalsh(dynamical(cry, pot, q, Phi))
    return np.sign(w2) * np.sqrt(np.abs(w2)) * THZ


def _flat(Phi):
    """
    Force constants as arrays instead of a dict: (blocks, i, j, R).

    Cached, because the caller almost always builds Phi once and then asks for
    hundreds of q-points from it.

    **The cache holds a reference to Phi itself, and checks it.**  Keying on
    id() alone is wrong: CPython reuses an id as soon as the object behind it is
    collected, so a freed force-constant dict hands its id to the next one of
    the same length and the cache returns another element's blocks.  That is not
    hypothetical - it silently gave zirconium and potassium frequencies wrong by
    2 and 12 THz while the first four elements of the same loop agreed to 1e-15.
    Holding the reference keeps the id alive and the identity test catches any
    remaining collision.
    """
    hit = _FLAT_CACHE.get(id(Phi))
    if hit is not None and hit[0] is Phi:
        return hit[1]
    keys = list(Phi.keys())
    out = (np.array([Phi[k] for k in keys]),                 # (B, 3, 3)
           np.array([k[0] for k in keys]),                   # (B,)
           np.array([k[1] for k in keys]),                   # (B,)
           np.array([k[2] for k in keys], float))            # (B, 3)
    if len(_FLAT_CACHE) > 64:
        _FLAT_CACHE.clear()
    _FLAT_CACHE[id(Phi)] = (Phi, out)
    return out


_FLAT_CACHE = {}


def frequencies_many(cry, pot, qs, Phi=None):
    """
    Frequencies at many q at once, shape (len(qs), 3N).  Same numbers as
    calling frequencies() in a loop, and about sixty times faster.

    The loop version rebuilds D(q) from the dict for every single q, which for
    a 2.6a pair cutoff means a few hundred Python iterations and a scalar
    np.exp per q-point.  On platinum that is 398 000 iterations for the
    8^3 U 9^3 stability mesh: 2.17 s, of which almost none is arithmetic.
    Here the blocks are flattened once, all the phases are formed in one
    matrix product, and eigvalsh is handed the whole stack - 35 ms.

    The remaining loop runs over BLOCKS, not q-points, so its length is set by
    the cutoff and not by how finely the zone is sampled.
    """
    Phi = force_constants(cry, pot) if Phi is None else Phi
    B, I, J, R = _flat(Phi)
    qs = np.atleast_2d(np.asarray(qs, float))
    n = len(cry.frac)
    ph = np.exp(2j * np.pi * (qs @ R.T))                     # (Q, B)
    D = np.zeros((len(qs), 3 * n, 3 * n), complex)
    for b in range(len(B)):
        D[:, 3*I[b]:3*I[b]+3, 3*J[b]:3*J[b]+3] += ph[:, b, None, None] * B[b]
    M = np.repeat(cry.mass, 3)
    D /= np.sqrt(np.outer(M, M))
    D = 0.5 * (D + np.conj(np.transpose(D, (0, 2, 1))))
    w2 = np.linalg.eigvalsh(D)
    return np.sign(w2) * np.sqrt(np.abs(w2)) * THZ


# --------------------------------------------------------------------------
#  elastic constants
# --------------------------------------------------------------------------
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def _eps(v, h):
    """Voigt index -> strain tensor with engineering shear"""
    e = np.zeros((3, 3))
    a, b = VOIGT[v]
    if a == b:
        e[a, a] = h
    else:
        e[a, b] = e[b, a] = h / 2.0
    return e


def stress(cry, pot, h=1e-5):
    """Cauchy stress in GPa via central differences of the energy under strain"""
    nat = len(cry.frac)
    s = np.zeros(6)
    for v in range(6):
        ep = energy(cry.strained(_eps(v, h)), pot) * nat
        em = energy(cry.strained(_eps(v, -h)), pot) * nat
        s[v] = (ep - em) / (2 * h) / cry.vol
    return s * EV_A3_TO_GPA


def elastic(cry, pot, h=2e-4, relax=True):
    """
    Elastic constants in GPa.

    Returns (C_relaxed, C_born).  For one atom per cell the two are identical
    by symmetry; the difference is the non-affine correction.
    """
    nat = len(cry.frac)
    V = cry.vol

    #  ---- Born term: energy curvature under affine strain, ions frozen ----
    E0 = energy(cry, pot) * nat
    Cb = np.zeros((6, 6))
    ep_cache = {}
    for v in range(6):
        for s in (+1, -1):
            ep_cache[(v, s)] = energy(cry.strained(_eps(v, s*h)), pot) * nat
    for v in range(6):
        Cb[v, v] = (ep_cache[(v, 1)] - 2*E0 + ep_cache[(v, -1)]) / (h*h)
    for v in range(6):
        for w in range(v+1, 6):
            epp = energy(cry.strained(_eps(v, h) + _eps(w, h)), pot) * nat
            emm = energy(cry.strained(_eps(v, -h) + _eps(w, -h)), pot) * nat
            Cb[v, w] = Cb[w, v] = (epp + emm
                                   - ep_cache[(v, 1)] - ep_cache[(w, 1)]
                                   - ep_cache[(v, -1)] - ep_cache[(w, -1)]
                                   + 2*E0) / (2*h*h)
    Cb *= EV_A3_TO_GPA / V

    if not relax or nat == 1:
        return Cb.copy(), Cb

    #  ---- non-affine correction ----
    #  Xi[i alpha, v] = d(force on atom i)/d(strain v), by central difference
    Xi = np.zeros((3*nat, 6))
    for v in range(6):
        fp = -_gradient(cry.strained(_eps(v, h)), pot)
        fm = -_gradient(cry.strained(_eps(v, -h)), pot)
        Xi[:, v] = (fp - fm).ravel() / (2*h)
    #  Gamma-point Hessian, mass-independent.  Built here rather than at the
    #  top of the function: a one-atom cell returns above and never needs it,
    #  which is what lets the angular term be tested on the cubic metals with
    #  the energy alone.
    Phi = force_constants(cry, pot)
    H = np.zeros((3*nat, 3*nat))
    for (i, j, R), blk in Phi.items():
        H[3*i:3*i+3, 3*j:3*j+3] += blk
    H = 0.5 * (H + H.T)
    #  pseudo-inverse: project out the three rigid translations
    Hp = np.linalg.pinv(H, rcond=1e-8)
    Cr = Cb - (Xi.T @ Hp @ Xi) * EV_A3_TO_GPA / V
    return 0.5 * (Cr + Cr.T), Cb


def _gradient(cry, pot):
    """dE/dr_i for every atom, shape (nat, 3), eV/A"""
    #  The three-body force below is phi3'(r1 + r2) along the two legs, which is
    #  the whole story only while phi3 is angle-free.  Unlike force_constants
    #  this used to return a wrong number rather than refuse - harmless for the
    #  cubic metals, because elastic() returns before it for a one-atom cell,
    #  but silently wrong for hcp.
    if getattr(pot, "lam2", 0.0) or getattr(pot, "lam4", 0.0):
        import angfc
        return angfc.gradient(cry, pot)
    nat = len(cry.frac)
    G = np.zeros((nat, 3))
    for (i, j, R, d, r) in neighbours(cry, pot.rcut2):
        #  ordered pairs, factor 1/2 for double counting
        G[i] -= 0.5 * pot.phi2(r, 1) * (d / r)
        G[j] += 0.5 * pot.phi2(r, 1) * (d / r)
    if pot.C:
        for (i, (ja, _), (jb, _), ra, rb, na, nb_) in triplets(cry, pot):
            g1 = pot.phi3(ra + rb, 1)
            G[ja] += g1 * na
            G[jb] += g1 * nb_
            G[i] -= g1 * (na + nb_)
    return G


# --------------------------------------------------------------------------
#  Brillouin-zone sampling and harmonic thermodynamics
# --------------------------------------------------------------------------
def mesh(n):
    """Monkhorst-Pack fractional q points, shape (n^3, 3)"""
    g = (np.arange(n) + 0.5) / n - 0.5
    return np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T


def spectrum(cry, pot, nq=8):
    """all phonon frequencies (THz) on an nq^3 mesh, shape (nq^3, 3N)"""
    Phi = force_constants(cry, pot)
    return frequencies_many(cry, pot, mesh(nq), Phi)


def thermo(freqs, T, nat):
    """
    Harmonic thermodynamics per atom from frequencies in THz.
    Returns dict with zpe (eV), F (eV), S and Cv in J/(mol K).
    """
    h_eV_THz = 4.135667696e-3                 # Planck constant in eV/THz
    w = np.asarray(freqs).ravel()
    w = w[w > 1e-6]                           # drop acoustic zeros / imaginary
    nq = len(freqs)
    zpe = 0.5 * h_eV_THz * w.sum() / (nq * nat)
    out = {"zpe": zpe}
    if T <= 0:
        return {**out, "F": zpe, "S": 0.0, "Cv": 0.0}
    x = h_eV_THz * w / (KB * T)
    x = np.clip(x, 1e-12, 700.0)
    F = (0.5*h_eV_THz*w + KB*T*np.log1p(-np.exp(-x))).sum() / (nq*nat)
    S = KB * ((x/(np.exp(x)-1.0)) - np.log1p(-np.exp(-x))).sum() / (nq*nat)
    Cv = KB * ((x*x*np.exp(x))/(np.exp(x)-1.0)**2).sum() / (nq*nat)
    J = 96485.33212                           # eV/(mol) -> J/mol
    return {**out, "F": F, "S": S*J, "Cv": Cv*J}
