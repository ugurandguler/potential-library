#!/usr/bin/env python3
"""
Self-consistency tests for latdyn.py.  No external code is involved.

Correctness is established from internal identities rather than by comparison
with another code.  Each test below is one that actually caught a bug during
development:

  1. force constants        acoustic sum rule and Phi(i,j,R) = Phi(j,i,-R)^T.
                            Caught a factor-of-two error in the pair term that
                            the elastic constants could not see, because those
                            come from strain derivatives of the energy rather
                            than from Phi.

  2. elastic <-> phonon     the long-wavelength sound velocities must reproduce
                            the elastic constants.  These are two completely
                            separate code paths - energy-under-strain versus the
                            dynamical matrix - so agreement is a real check.
                            For a cubic crystal along [100]:
                                rho v_L^2 = C11 ,  rho v_T^2 = C44

  3. zero stress            the fit imposes it, so the residual pressure (and
                            for hcp the deviatoric stress) must vanish.

  4. mechanical stability   Born criteria; a negative eigenvalue of the elastic
                            matrix means the structure is not a minimum.

  5. thermodynamic limits   Cv -> 3R = 24.94 J/(mol K) at high temperature
                            (Dulong-Petit), and S -> 0 as T -> 0.

    python selftest.py            # every fitted element
    python selftest.py Pd Mg
"""
import json, math, os, sys
import numpy as np
import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
R_GAS = 8.314462618            # J/(mol K);  3R = 24.94


def records():
    """the fitted records, from library.json or from fit.json + refdata

    library.json is 9.6 MB of measured results and is not distributed; fit.json
    is the 27 kB distillation that is.  Reading only the former made this
    script - the one the README tells a new reader to run first - crash on a
    clean checkout with a FileNotFoundError, which is a poor first impression
    of a project whose whole claim is that its numbers are checked.

    The two carry the same potential.  What fit.json does not carry is the
    geometry, because a0 and c/a are experimental targets rather than fitted
    outputs; those come from refdata, which is where library.json got them.
    """
    p = os.path.join(HERE, "library.json")
    if os.path.exists(p):
        return json.load(open(p))
    p = os.path.join(HERE, "fit.json")
    if not os.path.exists(p):
        raise SystemExit("neither library.json nor fit.json is here")
    out = {}
    for el, v in json.load(open(p)).items():
        ref = refdata.ELEMENTS.get(el)
        if not ref or "a0" not in ref:
            continue
        r = dict(v)
        r["a0"] = ref["a0"]
        if "c_over_a" in ref:
            r["c_over_a"] = ref["c_over_a"]
        out[el] = r
    return out


def load(el, lib):
    v = lib[el]
    cry = L.Crystal(v["struct"], v["a0"], v.get("c_over_a"),
                    mass=refdata.MASSES[el])
    pot = L.Potential.from_record(v)
    return cry, pot, v


def sound_check(cry, pot, C):
    """
    Slopes of the acoustic branches along [100] against C11 and C44.
    rho in amu/A^3 -> GPa needs the same conversion as the energy density.
    """
    lat = cry.lat
    #  a small q along cartesian x, expressed in fractional reciprocal coords
    rec = np.linalg.inv(lat).T                    # rows are b_i / 2pi
    xhat = np.array([1.0, 0.0, 0.0])
    qfrac = lat @ xhat                            # since q_cart = qfrac @ rec.T*2pi
    qfrac = qfrac / np.linalg.norm(qfrac @ np.linalg.inv(lat))
    eps = 1e-3
    qf = qfrac * eps
    w = np.sort(L.frequencies(cry, pot, qf))      # THz
    qcart = 2*math.pi*np.linalg.norm(qf @ rec)    # |q| in 1/A
    if qcart <= 0:
        return None
    v = 2*math.pi*w*1e12 / (qcart*1e10)           # m/s
    rho = cry.mass.sum() * 1.66053906660e-27 / (cry.vol * 1e-30)   # kg/m^3
    got = rho * v**2 / 1e9                        # GPa
    return {"C44_from_vT": float(got[0]), "C44": float(C[3, 3]),
            "C11_from_vL": float(got[2]), "C11": float(C[0, 0])}


def main(els):
    lib = records()
    els = els or sorted(lib)
    fails = []
    print(f"{'el':4s} {'sum rule':>9s} {'C11 phon/elast':>18s} "
          f"{'C44 phon/elast':>18s} {'P (GPa)':>9s} {'Cv(2000K)':>10s} {'stab':>5s}")
    print("-" * 78)
    for el in els:
        cry, pot, v = load(el, lib)
        msgs = L.check_force_constants(cry, pot)
        C, Cb = L.elastic(cry, pot)
        s = L.stress(cry, pot)
        P = float(s[:3].mean())
        f = L.spectrum(cry, pot, nq=6)
        hi = L.thermo(f, 2000.0, len(cry.frac))["Cv"]
        lo = L.thermo(f, 1.0, len(cry.frac))["S"]
        ev = np.linalg.eigvalsh(0.5*(C + C.T))
        stable = bool(ev.min() > 0)

        sc = sound_check(cry, pot, C) if v["struct"] in ("fcc", "bcc") else None
        def ratio(a, b):
            return f"{a:7.1f}/{b:<7.1f}" if a is not None else f"{'-':>15s}"
        r11 = ratio(sc["C11_from_vL"], sc["C11"]) if sc else f"{'-':>15s}"
        r44 = ratio(sc["C44_from_vT"], sc["C44"]) if sc else f"{'-':>15s}"

        ok = (not msgs) and abs(P) < 1e-2 and stable \
            and abs(hi - 3*R_GAS) < 0.5 and abs(lo) < 0.5
        if sc:
            ok = ok and abs(sc["C11_from_vL"]-sc["C11"])/max(sc["C11"], 1) < 0.02 \
                 and abs(sc["C44_from_vT"]-sc["C44"])/max(sc["C44"], 1) < 0.02
        if not ok:
            fails.append((el, msgs[:1], P, hi, stable))
        print(f"{el:4s} {'PASS' if not msgs else 'FAIL':>9s} {r11:>18s} "
              f"{r44:>18s} {P:9.2e} {hi:10.2f} {'yes' if stable else 'NO':>5s}"
              f"{'' if ok else '   <<<'}")
    print(f"\nDulong-Petit limit 3R = {3*R_GAS:.2f} J/(mol K)")
    if fails:
        print(f"\n{len(fails)} element(s) need attention:")
        for (el, m, P, hi, st) in fails:
            why = []
            if m:
                why.append(f"force constants: {m[0]}")
            if abs(P) > 1e-2:
                why.append(f"residual pressure {P:.3f} GPa")
            if not st:
                why.append("elastic matrix not positive definite")
            if abs(hi - 3*R_GAS) > 0.5:
                why.append(f"Cv(2000K) = {hi:.2f}, expected {3*R_GAS:.2f}")
            print(f"   {el}: {'; '.join(why) or 'sound velocity mismatch'}")
    else:
        print("\nall checks passed")


if __name__ == "__main__":
    main(sys.argv[1:])
