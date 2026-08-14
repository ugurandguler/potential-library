#!/usr/bin/env python3
"""
How much does the 3-body term matter?  i.e. what does a pair-only LAMMPS run lose?

Pair-only dynamical matrix (central pair potential, monatomic fcc, 8 shells):
    Phi_ab(r) = (alpha - beta) x_a x_b / r^2 + beta delta_ab ,  alpha=phi''(r), beta=phi'(r)/r
    D_ab(q)   = sum_j Phi_ab(r_j) [1 - cos(q.r_j)]
    M w^2     = eig(D)

3-body part, exactly as published (Akgun&Ugur eq.12, after Mishra et al.):
    D^m_aa = 4 g [4 - 2C_2i - C_i(C_j + C_k)]
    D^m_ab = 4 g [C_i(C_j + C_k) - 2]          C_i = cos(pi q_i), C_2i = cos(2 pi q_i)
"""
import math

J_PER_EV, AMU = 1.602176634e-19, 1.66053906660e-27
EV_A2_TO_NPM = J_PER_EV / 1e-20

# ---- validated potential (D_printed/5), see akgun_ugur_lammps.py ----
def jm(D29): return D29 * 1e-29 / J_PER_EV * 1e10 / 5.0
POT = {"Pd": dict(m=2.5, D=jm(7.66429), al=2.41786, r0=2.77323),
       "Fe": dict(m=3.5, D=jm(6.19854), al=2.16735, r0=2.77034)}
A_ALLOY, X_FE = 3.8720, 0.10
GAMMA = -0.623772                      # N/m, alloy, from the paper
MASS  = ((1 - X_FE) * 106.42 + X_FE * 55.845) * AMU

def pair(r, m, D, al, r0):
    A = D / (2.0 * (m - 1.0)); u = al * (r0 - r)
    e1, em = math.exp(u), math.exp(m * u)
    g, gp, gpp = em - m*e1, m*al*(e1-em), m*al*al*(m*em-e1)
    return A*g/r, A*(gp/r - g/r**2), A*(gpp/r - 2*gp/r**2 + 2*g/r**3)

def fcc_nbrs(nshell=8):
    pts, R = [], 5
    for i in range(-R, R+1):
        for j in range(-R, R+1):
            for k in range(-R, R+1):
                if (i+j+k) % 2 or (i==j==k==0): continue
                pts.append((i, j, k))
    smax = sorted({i*i+j*j+k*k for i,j,k in pts})[:nshell][-1]
    return [p for p in pts if p[0]**2+p[1]**2+p[2]**2 <= smax]
NBRS = fcc_nbrs()

# concentration-averaged radial/tangential force constants, in N/m
def alloy_fc(r):
    a_ = b_ = 0.0
    for nm, w in (("Pd", 1-X_FE), ("Fe", X_FE)):
        p = POT[nm]
        _, fp, fpp = pair(r, p['m'], p['D'], p['al'], p['r0'])
        a_ += w * fpp * EV_A2_TO_NPM
        b_ += w * (fp/r) * EV_A2_TO_NPM
    return a_, b_

def jacobi3(A):
    a = [row[:] for row in A]
    for _ in range(60):
        off = max(((i, j) for i in range(3) for j in range(i+1, 3)),
                  key=lambda t: abs(a[t[0]][t[1]]))
        i, j = off
        if abs(a[i][j]) < 1e-18: break
        th = 0.5 * math.atan2(2*a[i][j], a[j][j]-a[i][i])
        c, s = math.cos(th), math.sin(th)
        for k in range(3):
            aik, ajk = a[i][k], a[j][k]
            a[i][k], a[j][k] = c*aik - s*ajk, s*aik + c*ajk
        for k in range(3):
            aki, akj = a[k][i], a[k][j]
            a[k][i], a[k][j] = c*aki - s*akj, s*aki + c*akj
    return sorted(a[k][k] for k in range(3))

def omega(q, with_3body):
    """q in units of 2pi/a ; returns 3 frequencies in 1e13 rad/s"""
    a0 = A_ALLOY / 2.0
    D = [[0.0]*3 for _ in range(3)]
    for (mi, ni, li) in NBRS:
        v = (mi, ni, li)
        s = mi*mi + ni*ni + li*li
        r = a0 * math.sqrt(s)
        al_, be_ = alloy_fc(r)
        phase = 1.0 - math.cos(math.pi * (q[0]*mi + q[1]*ni + q[2]*li))
        for A in range(3):
            for B in range(3):
                x = a0 * v[A]; y = a0 * v[B]
                Pab = (al_ - be_) * x * y / r**2 + (be_ if A == B else 0.0)
                D[A][B] += Pab * phase
    if with_3body:
        C  = [math.cos(math.pi * q[i]) for i in range(3)]
        C2 = [math.cos(2*math.pi * q[i]) for i in range(3)]
        for A in range(3):
            j, k = (A+1) % 3, (A+2) % 3
            D[A][A] += 4*GAMMA * (4 - 2*C2[A] - C[A]*(C[j]+C[k]))
        for A in range(3):
            for B in range(3):
                if A == B: continue
                j, k = (A+1) % 3, (A+2) % 3
                D[A][B] += 4*GAMMA * (C[A]*(C[j]+C[k]) - 2)
    ev = jacobi3(D)
    out = []
    for e in ev:
        w2 = e / MASS
        out.append(math.copysign(math.sqrt(abs(w2)), w2) / 1e13)
    return out

PTS = [("X  [100]",  (1.0, 0.0, 0.0)),
       ("[0.5 0 0]", (0.5, 0.0, 0.0)),
       ("K  [110]",  (0.75, 0.75, 0.0)),
       ("[0.5.5 0]", (0.5, 0.5, 0.0)),
       ("L  [111]",  (0.5, 0.5, 0.5)),
       ("[.25x3]",   (0.25, 0.25, 0.25))]

print("Pd-10%Fe phonon frequencies (1e13 rad/s)   [Fig.1 y-axis tops out at 5]")
print(f"{'q point':11s} | {'pair only (dashed)':>26s} | {'pair+3body (solid)':>26s} | softening")
print("-" * 96)
for nm, q in PTS:
    p1, p2 = omega(q, False), omega(q, True)
    s = [f"{100*(1-b/a):.0f}%" if a > 1e-6 else "-" for a, b in zip(p1, p2)]
    print(f"{nm:11s} | " + " ".join(f"{v:8.3f}" for v in p1) +
          " | " + " ".join(f"{v:8.3f}" for v in p2) + " | " + " ".join(f"{t:>4s}" for t in s))

print("\nNote: pair-only == the DASHED curves of Fig.1; pair+3body == the SOLID curves")
print("(the solid ones are what the papers show agreeing with experiment).")
