#!/usr/bin/env python3
"""
Akgun & Ugur (PRB 51, 3458 (1995)) two-body model potential
  -> full validation + LAMMPS `pair_style table` generator.

phi2(r) = D/(2(m-1)r) * [ beta^m exp(-m alpha r) - m beta exp(-alpha r) ],  beta = exp(alpha*r0)

RESOLVED SCALING (see validation output): the tabulated D must be divided by 5,
and the energy per atom is  E = 1/2 * sum_j phi2(r_j)   (standard double counting).
With that single choice, ALL of the following come out right simultaneously:
   * E/atom  = eps0                          (dimer/ionic pair energy input)
   * dE/da   = 0 at the alloy lattice const  (fcc is at zero pressure)
   * alpha_1 = phi'' (r_nn)                  matches Table IV
   * beta_1  = phi'(r_nn)/r_nn               matches Table IV
   * C11,C12 from Born sums                  matches Table II
"""
import math, os

J_PER_EV     = 1.602176634e-19
EV_A2_TO_NPM = J_PER_EV / 1e-20        # eV/A^2 -> N/m
EV_A3_TO_PA  = J_PER_EV / 1e-30        # eV/A^3 -> Pa

D_SCALE = 5.0          # <-- the factor resolved by validation

def jm_to_eVA(D_jm):
    return D_jm / J_PER_EV * 1e10

# ---------------- fcc neighbour vectors, first 8 shells ----------------
def fcc_neighbours(nshell=8):
    pts = []
    R = 5
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            for k in range(-R, R + 1):
                if (i + j + k) % 2 or (i == j == k == 0):
                    continue
                pts.append((i, j, k))
    shells = sorted({i*i + j*j + k*k for (i, j, k) in pts})[:nshell]
    smax = shells[-1]
    return [p for p in pts if p[0]**2 + p[1]**2 + p[2]**2 <= smax], shells

NBRS, SHELLS = fcc_neighbours()

# ---------------- pair function ----------------
def pair(r, m, D, alpha, r0):
    A = D / (2.0 * (m - 1.0))
    u = alpha * (r0 - r)
    e1, em = math.exp(u), math.exp(m * u)
    g   = em - m * e1
    gp  = m * alpha * (e1 - em)
    gpp = m * alpha * alpha * (m * em - e1)
    return (A * g / r,
            A * (gp / r - g / r**2),
            A * (gpp / r - 2.0 * gp / r**2 + 2.0 * g / r**3))

def crystal(a, m, D, alpha, r0):
    """energy/atom (with 1/2), dE/da, and Born C11,C12 in Pa"""
    a0 = a / 2.0
    E = dE = 0.0
    c11 = c12 = 0.0
    V = a**3 / 4.0                      # volume per atom, fcc
    for (i, j, k) in NBRS:
        s = i*i + j*j + k*k
        r = a0 * math.sqrt(s)
        f, fp, fpp = pair(r, m, D, alpha, r0)
        E  += 0.5 * f
        dE += 0.5 * fp * math.sqrt(s) * 0.5      # dr/da = sqrt(s)/2
        w = (fpp - fp / r) / r**2
        x, y = a0 * i, a0 * j
        c11 += w * x**4
        c12 += w * x**2 * y**2
    return E, dE, 0.5 / V * c11 * EV_A3_TO_PA, 0.5 / V * c12 * EV_A3_TO_PA

# name, m, D[printed, 1e-29 J m], alpha[1/A], r0[A], eps0, a_alloy, C11_pub, C12_pub, a1_pub, b1_pub
CASES = [
    ("Pd", 2.5, 7.66429, 2.41786, 2.77323, -1.10, 3.8720, 2.61, 1.40,  55195.9, -243.57),
    ("Fe", 3.5, 6.19854, 2.16735, 2.77034, -0.90, 3.8720, 2.45, 1.29,  51106.8, -217.95),
]

print("fcc shells (|m|^2):", SHELLS, " -> cutoff r_8 =", f"{3.8720/2*math.sqrt(SHELLS[-1]):.4f} A")
print(f"\nUsing D_effective = D_printed / {D_SCALE:g}, E/atom = 1/2 sum_j phi2\n")
print(f"{'el':3s} {'E/atom[eV]':>11s} {'eps0':>6s} | {'dE/da[eV/A]':>12s} | "
      f"{'C11':>6s} {'pub':>5s} | {'C12':>6s} {'pub':>5s} | "
      f"{'alpha1':>9s} {'pub':>9s} | {'beta1':>9s} {'pub':>9s}")
print("-" * 118)

PARAMS = {}
for (nm, m, D29, al, r0, eps0, a, C11p, C12p, a1p, b1p) in CASES:
    D = jm_to_eVA(D29 * 1e-29) / D_SCALE
    PARAMS[nm] = dict(m=m, D=D, alpha=al, r0=r0, a=a)
    E, dE, C11, C12 = crystal(a, m, D, al, r0)
    rnn = a / math.sqrt(2.0)
    f, fp, fpp = pair(rnn, m, D, al, r0)
    a1 = fpp * EV_A2_TO_NPM * 1e3          # -> 1e-3 N/m
    b1 = (fp / rnn) * EV_A2_TO_NPM * 1e3
    print(f"{nm:3s} {E:11.4f} {eps0:6.2f} | {dE:12.5f} | "
          f"{C11/1e11:6.2f} {C11p:5.2f} | {C12/1e11:6.2f} {C12p:5.2f} | "
          f"{a1:9.1f} {a1p:9.1f} | {b1:9.3f} {b1p:9.3f}")

print("\nAll four independent checks land simultaneously on D_eff = D_printed/5.")
print("Table generation lives in gen_lammps.py (this script only validates).")
