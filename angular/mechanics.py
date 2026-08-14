#!/usr/bin/env python3
"""
Everything the elastic tensor alone can tell you, worked out from the tensor
definitions rather than from any tabulated special case.

The single object is the compliance tensor s_ijkl = inv(C) in Voigt form,
un-contracted with the usual factors: a Voigt index above 3 carries a 1/2 per
occurrence, so

    s_ijkl = S_mn / (f(m) f(n)),      f(m) = 1 for m < 3, 2 otherwise.

From it (Nye, *Physical Properties of Crystals*, ch. VIII):

    1/E(u)      = s_ijkl u_i u_j u_k u_l           Young's modulus
    beta(u)     = s_ijkk u_i u_j                   linear compressibility
    1/G(u,v)    = 4 s_ijkl u_i v_j u_k v_l         shear, u perp v
    nu(u,v)     = -(s_ijkl u_i u_j v_k v_l) E(u)   Poisson's ratio

The averages are Voigt (uniform strain), Reuss (uniform stress) and their Hill
mean; the anisotropy index is A^U = 5 G_V/G_R + B_V/B_R - 6, which is zero only
for an isotropic solid (Ranganathan and Ostoja-Starzewski, PRL 101, 055504).

Nothing here knows about the potential - it takes a 6x6 matrix in GPa.
"""
import numpy as np

HBAR = 1.054571817e-34
KB = 1.380649e-23
AMU = 1.66053907e-27


def _voigt_pairs():
    return [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def compliance_tensor(C):
    """s_ijkl (3,3,3,3) in 1/GPa from the 6x6 stiffness in GPa"""
    S = np.linalg.inv(np.asarray(C, dtype=float))
    pair = _voigt_pairs()
    f = lambda m: 1.0 if m < 3 else 2.0
    s = np.zeros((3, 3, 3, 3))
    for m, (i, j) in enumerate(pair):
        for n, (k, l) in enumerate(pair):
            v = S[m, n] / (f(m) * f(n))
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    s[a, b, c, d] = v
    return s, S


def _sphere(n):
    """n roughly equidistributed unit vectors (Fibonacci lattice)"""
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.sin(phi) * np.cos(theta),
                     np.sin(phi) * np.sin(theta),
                     np.cos(phi)], axis=1)


def _perp(u):
    """two unit vectors completing an orthonormal frame with u"""
    a = np.zeros_like(u)
    small = np.argmin(np.abs(u), axis=-1)
    a[np.arange(len(u)), small] = 1.0
    v = np.cross(u, a)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    w = np.cross(u, v)
    return v, w


def young(s, u):
    return 1.0 / np.einsum("ijkl,ni,nj,nk,nl->n", s, u, u, u, u)


def compressibility(s, u):
    return np.einsum("ijkk,ni,nj->n", s, u, u)


def shear_and_poisson(s, u, nphi=72):
    """min/max over the transverse direction v, for each u"""
    v0, w0 = _perp(u)
    ang = np.linspace(0, np.pi, nphi, endpoint=False)
    Eu = young(s, u)
    gmin = np.full(len(u), np.inf); gmax = np.full(len(u), -np.inf)
    nmin = np.full(len(u), np.inf); nmax = np.full(len(u), -np.inf)
    for a in ang:
        v = np.cos(a) * v0 + np.sin(a) * w0
        g = 1.0 / (4.0 * np.einsum("ijkl,ni,nj,nk,nl->n", s, u, v, u, v))
        nu = -np.einsum("ijkl,ni,nj,nk,nl->n", s, u, u, v, v) * Eu
        gmin = np.minimum(gmin, g); gmax = np.maximum(gmax, g)
        nmin = np.minimum(nmin, nu); nmax = np.maximum(nmax, nu)
    return gmin, gmax, nmin, nmax


def averages(C):
    C = np.asarray(C, dtype=float)
    S = np.linalg.inv(C)
    Bv = (C[0, 0] + C[1, 1] + C[2, 2]
          + 2 * (C[0, 1] + C[0, 2] + C[1, 2])) / 9.0
    Gv = (C[0, 0] + C[1, 1] + C[2, 2] - (C[0, 1] + C[0, 2] + C[1, 2])
          + 3 * (C[3, 3] + C[4, 4] + C[5, 5])) / 15.0
    Br = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2]
                + 2 * (S[0, 1] + S[0, 2] + S[1, 2]))
    Gr = 15.0 / (4 * (S[0, 0] + S[1, 1] + S[2, 2])
                 - 4 * (S[0, 1] + S[0, 2] + S[1, 2])
                 + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
    B, G = 0.5 * (Bv + Br), 0.5 * (Gv + Gr)
    return dict(B_V=Bv, B_R=Br, B_H=B, G_V=Gv, G_R=Gr, G_H=G,
                E_H=9 * B * G / (3 * B + G),
                nu_H=(3 * B - 2 * G) / (2 * (3 * B + G)),
                A_U=5 * Gv / Gr + Bv / Br - 6.0,
                pugh=B / G, cauchy=float(C[0, 1] - C[3, 3]))


def sound(B, G, rho):
    """longitudinal, transverse and Debye-averaged speeds, m/s; rho in kg/m^3"""
    vt = np.sqrt(G * 1e9 / rho)
    vl = np.sqrt((B + 4 * G / 3) * 1e9 / rho)
    vm = (1.0 / 3.0 * (2.0 / vt ** 3 + 1.0 / vl ** 3)) ** (-1.0 / 3.0)
    return float(vl), float(vt), float(vm)


def analyse(C, mass_amu=None, volume_A3=None, natoms=1, nu=4000):
    """
    Full report.  mass/volume are per unit cell and only feed the sound speeds
    and the Debye temperature; leave them out to skip those.
    """
    C = np.asarray(C, dtype=float)
    out = averages(C)
    s, S = compliance_tensor(C)

    u = _sphere(nu)
    E = young(s, u)
    beta = compressibility(s, u)
    gmin, gmax, nmin, nmax = shear_and_poisson(s, u)

    def rng(name, arr, vecs=u):
        i, j = int(np.argmin(arr)), int(np.argmax(arr))
        out[name + "_min"] = float(arr[i])
        out[name + "_max"] = float(arr[j])
        out[name + "_min_dir"] = [round(float(x), 4) for x in vecs[i]]
        out[name + "_max_dir"] = [round(float(x), 4) for x in vecs[j]]
        out[name + "_aniso"] = (float(arr[j] / arr[i])
                                if arr[i] > 0 else None)

    rng("E", E)
    rng("beta", beta)
    rng("G", gmin)
    out["G_max"] = float(gmax.max())
    out["G_max_dir"] = [round(float(x), 4) for x in u[int(np.argmax(gmax))]]
    out["G_aniso"] = float(gmax.max() / gmin.min()) if gmin.min() > 0 else None
    out["nu_min"] = float(nmin.min())
    out["nu_max"] = float(nmax.max())

    #  eigenvalues of C: all positive is the general Born criterion, and the
    #  smallest one says how close to instability the crystal sits
    ev = np.linalg.eigvalsh(C)
    out["eigenvalues"] = [float(x) for x in ev]
    out["born_stable"] = bool(ev.min() > 0)

    if mass_amu and volume_A3:
        rho = mass_amu * AMU / (volume_A3 * 1e-30)          # kg/m^3
        vl, vt, vm = sound(out["B_H"], out["G_H"], rho)
        n = natoms / (volume_A3 * 1e-30)                     # atoms / m^3
        out.update(rho=rho, v_l=vl, v_t=vt, v_m=vm,
                   debye=float(HBAR / KB * (6 * np.pi ** 2 * n) ** (1 / 3) * vm))
    return out


def plane_curves(C, plane, n=120):
    """
    E, beta and the shear / Poisson envelopes around one coordinate plane.
    plane is 'xy', 'xz' or 'yz'.  The angles are a plain linspace over the full
    turn, so they are not stored - the viewer regenerates them and only the
    curves travel.
    """
    s, _ = compliance_tensor(C)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    z = np.zeros_like(t)
    u = {"xy": np.stack([np.cos(t), np.sin(t), z], 1),
         "xz": np.stack([np.cos(t), z, np.sin(t)], 1),
         "yz": np.stack([z, np.cos(t), np.sin(t)], 1)}[plane]
    gmin, gmax, nmin, nmax = shear_and_poisson(s, u)
    return dict(E=[round(float(x), 2) for x in young(s, u)],
                beta=[round(float(x), 6) for x in compressibility(s, u)],
                G_min=[round(float(x), 2) for x in gmin],
                G_max=[round(float(x), 2) for x in gmax],
                nu_min=[round(float(x), 3) for x in nmin],
                nu_max=[round(float(x), 3) for x in nmax])


if __name__ == "__main__":
    #  checks that do not depend on any tabulated result
    #  1. an isotropic tensor must give A^U = 0 and E, G, nu independent of u
    E0, nu0 = 200.0, 0.30
    lam = E0 * nu0 / ((1 + nu0) * (1 - 2 * nu0))
    mu = E0 / (2 * (1 + nu0))
    Ciso = np.zeros((6, 6))
    Ciso[:3, :3] = lam
    Ciso[0, 0] = Ciso[1, 1] = Ciso[2, 2] = lam + 2 * mu
    Ciso[3, 3] = Ciso[4, 4] = Ciso[5, 5] = mu
    r = analyse(Ciso, nu=800)
    print(f"isotropic: A_U={r['A_U']:.2e}  E {r['E_min']:.3f}-{r['E_max']:.3f} "
          f"(input {E0})  nu {r['nu_min']:.4f}-{r['nu_max']:.4f} (input {nu0})")
    print(f"           G_H={r['G_H']:.3f} (input {mu:.3f})  "
          f"B_H={r['B_H']:.3f} (input {lam + 2 * mu / 3:.3f})")

    #  2. cubic: E along [100] must be 1/S11 and the Zener ratio must match
    C = np.zeros((6, 6))
    c11, c12, c44 = 169.0, 122.0, 75.3          # a strongly anisotropic cubic
    C[:3, :3] = c12
    C[0, 0] = C[1, 1] = C[2, 2] = c11
    C[3, 3] = C[4, 4] = C[5, 5] = c44
    s, S = compliance_tensor(C)
    E100 = young(s, np.array([[1.0, 0, 0]]))[0]
    print(f"cubic: E[100]={E100:.4f}  1/S11={1 / S[0, 0]:.4f}   "
          f"Zener={2 * c44 / (c11 - c12):.4f}")
    #  3. the linear compressibilities along any orthonormal triad must sum to
    #     1/B_R - a theorem, so a genuine check on the contraction factors
    r = analyse(C, nu=2000)
    tot = compressibility(s, np.eye(3)).sum()
    print(f"       sum beta(x,y,z) = {tot:.6f}   1/B_R = {1 / r['B_R']:.6f}"
          f"   (differ by {abs(tot - 1 / r['B_R']):.1e})")
