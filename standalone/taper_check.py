#!/usr/bin/env python3
"""
Verify the cutoff taper: is it continuous, and are its derivatives right?

The taper is the change that makes the potential usable in molecular dynamics,
and it touches the energy, the forces and the force constants at once.  Three
independent checks, because each catches something the others cannot:

  continuity     energy against lattice constant across a shell crossing.
                 Untapered, beryllium drops 0.16 eV/atom in one step when a
                 shell enters rcut2 - that single number is the reason for all
                 of this, so it is the first thing measured.
  gradient       analytic forces against a finite difference of the energy on a
                 displaced cell.  A taper that is continuous but whose
                 derivative is wrong gives a smooth energy and wrong dynamics,
                 which is worse than an obvious jump.
  force constants  analytic Phi against a finite difference of the gradient,
                 plus the acoustic sum rule.  The three-body blocks are where
                 the leg switches make the algebra genuinely different - the
                 term stops being a function of x = ra + rb alone - so this is
                 the check that earns its keep.

    python taper_check.py [element ...]
"""
import sys

import numpy as np

import latdyn as L
import refdata

TAPER = 0.85
DEFAULT = ["Be", "Ru", "Rh", "Fe", "Mg", "Cu", "Nb", "Ti"]


def build(el, taper):
    import json
    import os
    lib = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "library.json")))
    d, e = lib[el], refdata.ELEMENTS[el]
    pot = L.Potential(d["m"], d["D"], d["alpha"], d["r0"], d["gamma"],
                      C=d["C"], alpha3=d["alpha3"], rcut2=d["rcut2"],
                      rcut3=d["rcut3"], taper=taper)
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    return cry, pot, e


def biggest_step(el, taper):
    """largest kink in E(a) as the lattice is compressed.

    A plain step size will not do: E genuinely varies with a, so a smooth curve
    already changes by about 0.01 eV per 0.1 per cent of a0 and would be
    reported as a fault.  What a truncation leaves behind is a *discontinuity*,
    so the second difference is the thing to look at - small and slowly varying
    for a smooth curve, a spike where a shell crosses the sphere.
    """
    cry, pot, e = build(el, taper)
    vals = []
    f = 1.02
    while f > 0.90:
        c = L.Crystal(e["struct"], e["a0"] * f, e.get("c_over_a"),
                      mass=refdata.MASSES[el])
        vals.append(L.energy(c, pot))
        f -= 0.001
    v = np.array(vals)
    return float(np.abs(v[2:] - 2 * v[1:-1] + v[:-2]).max())


def displaced(cry, amp=0.06, seed=0):
    """same crystal with the basis pushed off site, so forces are non-trivial"""
    rng = np.random.default_rng(seed)
    c = L.Crystal(cry.struct, cry.a0, cry.c_over_a, mass=cry.mass[0])
    if len(c.frac) == 1:
        #  a one-atom cell can only translate rigidly, which every force law
        #  gets right; use a cell with a basis to test anything at all
        return None
    c.frac = c.frac + rng.normal(0.0, amp, c.frac.shape)
    return c


def _shift(c, i, k, h):
    """copy of c with atom i moved h angstrom along CARTESIAN direction k.

    pos = frac @ lat, so a displacement of h e_k needs frac += h (lat^-1)[k].
    Adding h to frac[i, k] instead moves the atom along lattice row k, which
    for hcp is neither of unit length nor parallel to e_k - the first version
    of this check did exactly that and reported the taper as broken when what
    was broken was the check.
    """
    cc = L.Crystal(c.struct, c.a0, c.c_over_a, mass=c.mass[0])
    cc.frac = c.frac.copy()
    cc.frac[i] += h * np.linalg.inv(cc.lat)[k]
    return cc


def grad_error(el, taper, h=1e-5):
    cry, pot, _ = build(el, taper)
    c = displaced(cry)
    if c is None:
        return None
    nat = len(c.frac)
    G = L._gradient(c, pot)
    worst = 0.0
    for i in range(nat):
        for k in range(3):
            ep = L.energy(_shift(c, i, k, +h), pot) * nat
            em = L.energy(_shift(c, i, k, -h), pot) * nat
            worst = max(worst, abs((ep - em) / (2 * h) - G[i, k]))
    return worst


def fc_error(el, taper, h=2e-4):
    """analytic Phi against a finite difference of the gradient"""
    cry, pot, _ = build(el, taper)
    c = displaced(cry)
    if c is None:
        return None, None
    Phi = L.force_constants(c, pot)
    nat = len(c.frac)
    #  on-site block only: it needs no image bookkeeping and is the sum of
    #  everything else by the acoustic rule, so an error anywhere shows here
    worst = 0.0
    for i in range(nat):
        for k in range(3):
            cols = [L._gradient(_shift(c, i, k, s * h), pot) for s in (+1, -1)]
            num = (cols[0] - cols[1]) / (2 * h)
            ana = np.zeros((nat, 3))
            for (a, b, R), blk in Phi.items():
                if b == i:
                    ana[a] += blk[:, k]
            worst = max(worst, np.abs(num - ana).max())
    #  acoustic sum rule
    asr = 0.0
    for i in range(nat):
        s = np.zeros((3, 3))
        for (a, b, R), blk in Phi.items():
            if a == i:
                s += blk
        asr = max(asr, np.abs(s).max())
    return worst, asr


def main():
    els = sys.argv[1:] or DEFAULT
    print(f"taper = {TAPER} of the cutoff\n")
    print(f"{'el':4s}{'step OFF':>11s}{'step ON':>10s}{'|dE/dx| err':>13s}"
          f"{'|Phi| err':>11s}{'ASR':>11s}")
    print("-" * 60)
    bad = []
    for el in els:
        off = biggest_step(el, None)
        on = biggest_step(el, TAPER)
        g = grad_error(el, TAPER)
        f, asr = fc_error(el, TAPER)
        print(f"{el:4s}{off:11.4f}{on:10.4f}"
              f"{(g if g is not None else float('nan')):13.2e}"
              f"{(f if f is not None else float('nan')):11.2e}"
              f"{(asr if asr is not None else float('nan')):11.2e}")
        if on > 0.02 * max(off, 1e-9) and on > 1e-3:
            bad.append(f"{el}: taper did not remove the step ({on:.4f})")
        for name, v, tol in (("gradient", g, 1e-6), ("force constants", f, 1e-4),
                             ("acoustic sum rule", asr, 1e-7)):
            if v is not None and v > tol:
                bad.append(f"{el}: {name} off by {v:.2e}")
    print()
    if bad:
        for b in bad:
            print("FAIL ", b)
        raise SystemExit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
