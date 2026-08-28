#!/usr/bin/env python3
"""
Is the alkali contraction physics, or is it the thermostat?

npt_expansion returns a strongly NEGATIVE thermal expansion for the re-cut
alkalis - Na -39, K -42, Rb -37, Cs -37 in 1e-6/K, against experiment near
+80.  Three explanations were checked and all three are ruled out: those
potentials are dynamically STABLE along the whole symmetry path (the
published arm is the one with imaginary modes at H-N), their Grueneisen
parameters are positive everywhere with not one negative mode out of 213,
and their zero-kelvin energy curve has a single minimum with no denser basin
to fall into.  Quasi-harmonic theory built from these very potentials
predicts POSITIVE expansion.  So the molecular dynamics and the lattice
dynamics disagree, and one of them is not describing the potential.

The suspect is the coupling constant, not the potential.  npt_expansion
fixes the thermostat damping at 0.1 ps for every element.  That is the usual
hundred timesteps, and for copper it is about one period of the highest
phonon - 7.5 THz is 0.13 ps.  Caesium's highest phonon is 1.07 THz, a period
of 0.93 ps, so the same 0.1 ps thermostat is seven times FASTER than the
motion it is meant to be sampling.  A Nose-Hoover chain driven faster than
the vibration it thermostats does not simply add noise; it can shift the
average volume the barostat then equilibrates to.

This runs caesium both ways, with copper alongside as a control - copper's
damping is already correct on this argument, so if copper moves too, the
explanation is something else.

    python npt_damp_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import npt_expansion as M          # noqa: E402

#  Same template, but the two damping constants scaled to the material and a
#  longer settle.  Pdamp is kept at ten times Tdamp, as npt_expansion has it.
#  KEEP THE ORIGINAL.  The first version of this file set `M.IN = M.IN` for
#  the baseline, which does nothing: once the slow template had been
#  installed for caesium it stayed installed, so copper's "current damping"
#  run silently used the slow one too and the two copper columns came out
#  identical to every digit.  That is what an identical control means here -
#  not agreement, but the same run reported twice.
ORIG = M.IN
SLOW = M.IN.replace("temp {T} {T} 0.1 iso 0.0 0.0 1.0",
                    "temp {T} {T} {td} iso 0.0 0.0 {pd}")


def alpha(el, tag, td=None, pd=None, label=""):
    ts = M.temps(el)
    a, out = [], []
    for T in ts:
        M.IN = ORIG if td is None else SLOW.replace(
            "{td}", str(td)).replace("{pd}", str(pd))
        r = M.one((el, tag, T))
        if "error" in r:
            print(f"  {el} {label} {T:6.1f} K  HATA: {r['error']}")
            return None
        a.append(r["ax"])
        out.append((T, r["ax"], r["Pavg"], r["Tavg"]))
    #  least squares slope of a(T) / a(T0)
    n = len(ts)
    mt = sum(ts) / n
    ma = sum(a) / n
    sl = (sum((ts[i] - mt) * (a[i] - ma) for i in range(n))
          / sum((t - mt) ** 2 for t in ts))
    al = sl / a[0] * 1e6
    for T, ax, P, Tv in out:
        print(f"  {el:3s} {label:14s} {T:6.1f} K   a={ax:.5f}   "
              f"P={P:7.2f}   T={Tv:6.1f}")
    print(f"  {el:3s} {label:14s} alpha = {al:+7.1f} e-6/K")
    return al


def main():
    base = dict(td=None, pd=None)
    runs = [("Cs", "rc", base, "0.1 ps (mevcut)"),
            ("Cs", "rc", dict(td=1.0, pd=10.0), "1.0 ps (esler)"),
            ("Cu", "rc", base, "0.1 ps (mevcut)"),
            ("Cu", "rc", dict(td=1.0, pd=10.0), "1.0 ps (kontrol)")]
    got = {}
    for el, tag, kw, lab in runs:
        print(f"\n=== {el} {tag}, sonum {lab} ===", flush=True)
        got[(el, lab)] = alpha(el, tag, label=lab, **kw)
    print("\n" + "=" * 58)
    for k, v in got.items():
        print(f"  {k[0]:3s} {k[1]:18s} {v if v is None else round(v,1)}")
    print("\nCs isaret degistirir ve Cu yerinde kalirsa sebep sonum sabitidir.")
    print("Ikisi de kaymazsa sebep baska yerde.")


if __name__ == "__main__":
    main()
