#!/usr/bin/env python3
"""
Thermal expansion: does the form work where coordination does NOT change?

The vacancy test settled what this potential cannot do.  Removing an atom
changes the coordination, the bond energy should respond, nothing in phi2 +
phi3 makes it respond, and the formation energy comes out two to three times
too large for a reason no refitting can touch.

Thermal expansion asks the opposite question.  The atoms keep every neighbour
they had and simply sit further apart, so the quantity depends on the
anharmonicity of the same bonds the elastic constants already probe - the
third derivative of a curve whose second derivative was fitted.  It is
genuinely out of sample and it is inside the form's reach in a way defect
energies are not.  If it works, there is a real domain here: lattice dynamics
at finite temperature, thermal transport, anything driven by force constants
rather than by defects.  If it does not, the usable scope ends at the elastic
tensor and that is worth knowing precisely.

The method is the quasi-harmonic approximation rather than NPT dynamics, which
is both cheaper and the standard way this is compared with experiment:

    F(V, T) = E_static(V) + F_vib(V, T)

minimised over V at each temperature, with F_vib from the phonon spectrum
computed at that volume.  The expansion is then alpha = (1/a) da/dT.

The whole thing rests on the frequencies shifting with volume - if they did
not, F_vib would not depend on V and there would be no expansion at all.  That
shift is the Grueneisen parameter, and it is reported alongside because a
sensible alpha with an absurd gamma would mean the two errors had cancelled.

    python expansion.py                # tapered set, a few elements
    python expansion.py Cu Al Ag --nq 6
"""
import json
import os
import sys

import numpy as np

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
KB = 8.617333262e-5

#  Linear thermal expansion coefficient at 25 C, 1e-6 / K, from the CRC
#  Handbook of Chemistry and Physics tables.  Quoted because without a
#  reference the number below is uncalibrated, and this project has already
#  produced one uncalibrated negative result that was worth nothing until a
#  baseline was attached to it.  For the alkalis the tabulated values carry
#  several per cent of spread between compilations; for the refractories they
#  are good to about one per cent.
#
#  The edition, in the publisher's own recommended form:
#
#    David R. Lide, ed., CRC Handbook of Chemistry and Physics, Internet
#    Version 2005, <http://www.hbcpnetbase.com>, CRC Press, Boca Raton, FL,
#    2005.
#
#  Note that this is NOT the same edition as the melting points in refdata.py,
#  which name the 97th edition (2016).  Both are recorded rather than
#  reconciled: they are different tables and there is no reason a value should
#  have been copied from one to the other, but anyone comparing the two should
#  know they are not from the same printing.  CRC's figures do move between
#  editions - refdata.py's cohesive-energy note says so - which is why the
#  spread quoted above is the right precision to read these at.
ALPHA_EXP = {
    "Ag": 18.9, "Al": 23.1, "Au": 14.2, "Ba": 20.6, "Be": 11.3, "Ca": 22.3,
    "Cd": 30.8, "Co": 13.0, "Cr": 4.9, "Cs": 97.0, "Cu": 16.5, "Fe": 11.8,
    "Hf": 5.9, "Ir": 6.4, "K": 83.0, "Li": 46.0, "Lu": 9.9, "Mg": 24.8,
    "Mo": 4.8, "Na": 71.0, "Nb": 7.3, "Ni": 13.4, "Pb": 28.9, "Pd": 11.8,
    "Pt": 8.8, "Rb": 90.0, "Re": 6.2, "Rh": 8.2, "Ru": 6.4, "Sc": 10.2,
    "Sr": 22.5, "Ta": 6.3, "Ti": 8.6, "Tl": 29.9, "V": 8.4, "W": 4.5,
    "Y": 10.6, "Yb": 26.3, "Zn": 30.2, "Zr": 5.7,
}


def spectra(el, rec, xs, nq):
    """phonon spectrum and static energy at each scaled lattice parameter

    Computed once per volume rather than once per volume AND temperature.  The
    frequencies do not depend on T, only on the geometry, so recomputing them
    inside the temperature loop was five times the work for the same answer.
    """
    e = refdata.ELEMENTS[el]
    pot = L.Potential.from_record(rec)
    out = []
    for x in xs:
        cry = L.Crystal(e["struct"], e["a0"] * x, e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        fr = L.spectrum(cry, pot, nq)
        #  imaginary modes mean the structure is not a minimum at this volume;
        #  carried as None rather than dropped, so the caller can say so
        bad = (np.isnan(fr).any() or (np.asarray(fr) < -1e-6).any())
        out.append((None if bad else np.asarray(fr),
                    L.energy(cry, pot), len(cry.frac)))
    return out


def a_of_T(el, sp, xs, T):
    """equilibrium lattice parameter at temperature T

    First-order quasi-harmonic, and deliberately not a minimisation of the
    total free energy.  The shift being measured is two parts in a thousand of
    the lattice parameter, while E_static and F_vib are each of order an
    electronvolt, so fitting their sum and looking for its minimum subtracts
    two large numbers to get a small one and the answer depends on the fitting
    window.  It did: over a wide window sodium contracted on heating, over a
    narrow one iron and tungsten did too, and copper moved by thirteen per cent
    between the two.  None of those were physical - the heat-capacity weighted
    Grueneisen parameter was positive throughout, which forbids contraction.

    Instead the two pieces are fitted separately, which is well conditioned
    because each is smooth on its own:

        dE/dV + dF_vib/dV = 0   =>   dV = -(dF_vib/dV) / (d2E/dV2)

    E_static(V) carries no phonon noise at all and gives the curvature; F_vib
    is nearly linear in V over this range and gives the slope.  Nothing large
    is ever subtracted from anything large.
    """
    e = refdata.ELEMENTS[el]
    ok = [k for k, (fr, _, _) in enumerate(sp) if fr is not None]
    if len(ok) < 5:
        return None
    x = np.array([xs[k] for k in ok])
    v = x ** 3                                   # volume, in units of a0^3
    est = np.array([sp[k][1] for k in ok])       # static energy per atom
    fvib = np.array([L.thermo(sp[k][0], T, sp[k][2])["F"] for k in ok])

    ce = np.polyfit(v, est, 3)
    d1e, d2e = np.polyder(ce), np.polyder(ce, 2)
    #  the static minimum, from the static curve alone
    roots = [r.real for r in np.roots(d1e)
             if abs(r.imag) < 1e-9 and v.min() <= r.real <= v.max()
             and np.polyval(d2e, r.real) > 0]
    if not roots:
        return None
    v0 = min(roots, key=lambda z: np.polyval(ce, z))
    k = float(np.polyval(d2e, v0))               # d2E/dV2 > 0
    if k <= 0:
        return None
    slope = float(np.polyval(np.polyder(np.polyfit(v, fvib, 2)), v0))
    vT = v0 - slope / k
    if not (v.min() <= vT <= v.max()):
        return None
    return float(e["a0"] * vT ** (1.0 / 3.0))


def gruneisen(sp, xs, T=300.0):
    """Grueneisen parameter, both the plain mean and the heat-capacity weighted

    The weighted one is the thermodynamic quantity, and the distinction is not
    cosmetic.  alpha = gamma_w C_v / (3 B V), so the SIGN of the expansion is
    the sign of the weighted gamma - not of the mean.  Mode Grueneisen
    parameters can have mixed signs, and when the low-frequency ones are
    negative they dominate the weighting and the lattice contracts on heating.
    That is the mechanism of real negative-expansion materials, and reading the
    mean instead of the weighted average led to calling a negative alpha here
    "impossible" when it is merely wrong.
    """
    good = [(x, fr) for x, (fr, _, _) in zip(xs, sp) if fr is not None]
    if len(good) < 3:
        return None, None
    lv = np.array([3.0 * np.log(x) for x, _ in good])
    lw = np.array([np.log(np.asarray(fr)[np.asarray(fr) > 1e-6].mean())
                   for _, fr in good])
    g_mean = float(-np.polyfit(lv, lw, 1)[0])

    #  mode by mode, on the reference volume's branch ordering
    ref = min(range(len(good)), key=lambda k: abs(good[k][0] - 1.0))
    n = min(len(np.asarray(fr).ravel()) for _, fr in good)
    F = np.array([np.sort(np.asarray(fr).ravel())[:n] for _, fr in good])
    g_mode, w = [], []
    THZ_TO_K = 47.9924  # 1 THz -> K in hbar w / kB
    for m in range(n):
        f = F[:, m]
        if (f <= 1e-6).any():
            continue
        g_mode.append(-np.polyfit(lv, np.log(f), 1)[0])
        #  Einstein heat capacity of this mode as the weight
        u = THZ_TO_K * float(F[ref, m]) / T
        w.append(u * u * np.exp(u) / (np.exp(u) - 1.0) ** 2 if u < 500 else 0.0)
    if not g_mode or sum(w) <= 0:
        return g_mean, None
    g_w = float(np.dot(g_mode, w) / sum(w))
    return g_mean, g_w


def main():
    args = sys.argv[1:]
    nq = 6
    if "--nq" in args:
        i = args.index("--nq"); nq = int(args[i + 1]); del args[i:i + 2]
    which = "tap"
    if "--set" in args:
        i = args.index("--set"); which = args[i + 1]; del args[i:i + 2]
    lib = json.load(open(os.path.join(HERE, "library.json")))
    els = args or ["Cu", "Al", "Ag", "Au", "Ni", "Pb", "Mg", "W"]

    TS = [100.0, 200.0, 300.0, 400.0, 500.0]
    print(f"Thermal expansion, {which} set, quasi-harmonic, {nq}^3 mesh")
    print("Coordination does not change - the form should represent this.\n")
    print(f"{'el':4s}{'a(100K)':>10s}{'a(300K)':>10s}{'a(500K)':>10s}"
          f"{'alpha':>10s}{'expt':>9s}{'ratio':>7s}{'gamma':>7s}{'gamma_w':>9s}"
          f"   not")
    print("-" * 80)
    out = {}
    #  Narrow, and deliberately so.  The expansion being measured is two parts
    #  in a thousand of the lattice parameter over four hundred kelvin, while
    #  a scan of plus or minus six per cent reaches, for sodium, into the
    #  compression basin documented in section 10b - and a quartic fitted
    #  across that is pulled by a region the answer has nothing to do with.
    #  Sodium and tantalum both came back CONTRACTING on heating with a
    #  positive heat-capacity-weighted Grueneisen parameter, which cannot
    #  happen, and this is why.
    xs = np.linspace(0.97, 1.05, 9)
    for el in els:
        rec = lib[el] if which == "hard" else lib[el].get(which)
        if not rec:
            print(f"{el:4s}  no {which} record"); continue
        try:
            sp = spectra(el, rec, xs, nq)
        except Exception as ex:
            print(f"{el:4s}  error: {type(ex).__name__} {str(ex)[:40]}"); continue
        aa = [a_of_T(el, sp, xs, T) for T in TS]
        if any(a is None for a in aa):
            print(f"{el:4s}  hesaplanamadi (hayali kip?)"); continue
        aa = np.array(aa)
        #  alpha near 300 K from a straight line through the whole set; the
        #  curvature over 100-500 K is small enough that this is the number
        #  tables quote
        slope = np.polyfit(TS, aa, 1)[0]
        i300 = TS.index(300.0)
        alpha = slope / aa[i300] * 1e6
        g, g_w = gruneisen(sp, xs)
        exp = ALPHA_EXP.get(el)
        #  how many of the scanned volumes had an imaginary mode: a record
        #  that is only a minimum over part of the range is a different
        #  animal from one that is a minimum throughout, and alpha alone
        #  does not say which
        nbad = sum(1 for fr, _, _ in sp if fr is None)
        out[el] = {"a": list(aa), "T": TS, "alpha_1e6": float(alpha),
                   "gruneisen": g, "gruneisen_weighted": g_w,
                   "alpha_exp_1e6": exp,
                   "n_imaginary": nbad, "n_volumes": len(sp),
                   "ratio": (alpha / exp) if exp else None}
        note = "" if 0 < alpha < 200 else "SUPHELI"
        if nbad:
            note = (note + f"  {nbad}/{len(sp)} hacimde hayali kip").strip()
        print(f"{el:4s}{aa[0]:10.4f}{aa[i300]:10.4f}{aa[-1]:10.4f}"
              f"{alpha:10.1f}{(exp or 0):9.1f}"
              f"{(alpha / exp if exp else 0):7.2f}"
              f"{(g if g is not None else 0):7.2f}"
              f"{(g_w if g_w is not None else 0):9.2f}"
              f"   {note}")

    #  merge, never overwrite: a run over a few elements must not erase a run
    #  over all of them, which has happened to four other files in this project
    fp = os.path.join(HERE, f"expansion_{which}.json")
    old = {}
    if os.path.exists(fp):
        try:
            old = json.load(open(fp))
        except Exception:
            old = {}
    old.update(out)
    out = old
    json.dump(out, open(os.path.join(HERE, f"expansion_{which}.json"), "w"),
              indent=1, sort_keys=True)
    print(f"\n-> expansion_{which}.json")


if __name__ == "__main__":
    main()
