#!/usr/bin/env python3
"""
Which switched set to recommend for molecular dynamics, element by element.

The shell-gap re-cut (`rc`, shipped as <El>_recut.ugur) moves the pair cutoff
into a gap between neighbour shells.  It exists for 23 elements, and the
evidence has been one-sided for the ones it suits: the intrinsic stacking
fault comes out NEGATIVE in all 15 fcc/hcp tapered records that also have a
re-cut, and POSITIVE in all 15 re-cut ones, closer to experiment in all 8
that have a compiled value (Rosengaard and Skriver 1993).  It is also a per-element question, not a global
switch: the same re-cut makes sodium, potassium, rubidium, caesium and lithium
expand NEGATIVELY on heating and stiffen instead of softening.

So the recommendation is computed here, from what the library already holds,
by a rule written down in one place - not typed into the page as a list, which
is how the re-cut dispersion table went stale for a release.

An element's re-cut is RECOMMENDED for MD when all of these hold:

  cold    stable on mesh AND path; the fitted structure is its ground state;
          the 300 K molecular-dynamics screen holds the crystal
  warm    thermal expansion positive
          C11 does not STIFFEN between the lowest molecular-dynamics
          temperature and 300 K - the sign of dC11/dT is the screen that
          exposes the alkali failure
          300 K elastic constants, against III/29a's room-temperature values,
          no more than WORSE_PTS points worse than the published switched set
  cold    measured dispersion, where scored, no more than WORSE_PTS points
          worse than the published switched set

Otherwise the published switched set, <El>_taper.ugur, stays the
recommendation, and the reason is recorded.  Nothing is deleted: both files
keep shipping, and the field only says which one to use.

THE FALLBACK IS CHECKED TOO.  Falling back to the switched set used to be
unconditional, and for lithium that recommended a record with no barrier
against collapse: squeezed uniformly, its energy falls without limit
(-434 eV/atom by 0.92 a0).  fit.py rejects exactly that - a compression basin
whose barrier is below k_B T_melt, the threshold calibrated against the MD
screen on 30 of 30 parameter sets - but the shipped records were fitted before
the constraint existed, and a constrained refit of lithium's switched set has
no solution in either form (MAU or UG).  So when the re-cut is rejected and the
switched set fails that same test, the recommendation is "none": no record of
this element is fit for molecular dynamics, and the files say so.

THE STIFFENING CONDITION USES MOLECULAR-DYNAMICS POINTS ONLY.  Its first
version divided C11 at 300 K by the elasticT table's T = 0 row, and that row is
not a measurement: it is the fit's static target (barium reads exactly 12.600
GPa in all four arms).  The jump from that frozen-ion constant to the first
finite-temperature point is the non-affine internal relaxation the
fluctuation term supplies - physics, not a temperature dependence.  It lowers
barium's switched C11 by 12 % at 50 K and raises its re-cut C11 by 3 %, and
that jump was the whole of the "1 % stiffening" that first kept barium out.
Between MD points barium's re-cut softens by 2 % from 50 to 300 K.  The
corrected condition changed no other element's verdict: the alkalis still
stiffen x2.4-3.1 and ytterbium x1.055 between MD points.

Writes lib[el]["md_recommended"] = {"set": "rc" | "tap" | "none", "why": [...],
"rule": RULE, "facts": {...}}.  Elements with no re-cut get nothing.

    python recommend_recut.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import refdata                                       # noqa: E402

KB = 8.617333262e-5                                  # eV/K, as fit.py

WORSE_PTS = 2.0
C11_RATIO_MAX = 1.00
RULE = ("stable, ground state, MD screen; expansion > 0; C11 not stiffening "
        "between the lowest MD temperature and 300 K; 300 K elastic error and "
        f"measured dispersion no more than {WORSE_PTS:g} points worse than the "
        "switched set")

#  The same numbers as RC_DISP in make_gui.py and export_potentials.py:
#  curve_mae.py, arms rc and tap, rerun 2026-09-13.
DISP = {"Ag": (13.0, 12.7), "Au": (4.2, 4.4), "Ba": (8.3, 9.6),
        "Ca": (8.6, 8.3), "Cs": (6.4, 12.9), "Cu": (9.3, 10.0),
        "Ir": (9.0, 9.9), "K": (11.8, 13.9), "Li": (7.0, 25.9),
        "Mg": (8.3, 6.6), "Na": (5.8, 16.2), "Ni": (12.1, 13.4),
        "Pb": (12.9, 18.3), "Pd": (3.5, 3.6), "Pt": (5.0, 4.8),
        "Rb": (7.6, 11.8), "W": (17.2, 17.8), "Yb": (8.5, 8.6)}


def near(pts, T):
    pts = [p for p in pts if not p.get("above_melt")]
    return min(pts, key=lambda p: abs(p["T"] - T))


def err300(v, el, arm):
    et = (v.get("elasticT") or {}).get(arm)
    if not et:
        return None
    p = near(et["pts"], 300)
    ref = refdata.ELEMENTS[el]["Cij"]
    ks = [k for k in ("C11", "C12", "C13", "C33", "C44")
          if k in ref and p.get(k) is not None]
    return 100.0 * sum(abs(p[k] - ref[k]) / abs(ref[k]) for k in ks) / len(ks)


def c11_md(v, arm):
    """(T_low, T_high, C11(T_high)/C11(T_low)) over molecular-dynamics points.

    T = 0 is excluded: it is the static target, not an MD point."""
    et = (v.get("elasticT") or {}).get(arm)
    if not et:
        return None
    pts = sorted((p for p in et["pts"]
                  if not p.get("above_melt") and p["T"] > 0),
                 key=lambda p: p["T"])
    if len(pts) < 2:
        return None
    lo = pts[0]
    hi = min(pts[1:], key=lambda p: abs(p["T"] - 300))
    return lo["T"], hi["T"], hi["C11"] / lo["C11"]


def decide(el, v):
    r, t = v["rc"], v["tap"]
    why = []
    if r.get("stable") is not True:
        why.append("imaginary modes on the mesh or path")
    if (r.get("ground") or {}).get("ok") is not True:
        why.append("the fitted structure is not its ground state")
    md = r.get("md_screen") or {}
    if md.get("collapsed") or md.get("lost"):
        why.append("does not hold the crystal in the MD screen")
    a = (r.get("expansion") or {}).get("alpha_1e6")
    if a is None or a <= 0:
        why.append(f"thermal expansion {a:.0f}e-6/K, negative" if a is not None
                   else "no thermal expansion measured")
    c = c11_md(v, "rc")
    if c is None:
        why.append("no elastic constants against temperature")
    elif c[2] > C11_RATIO_MAX:
        why.append(f"C11 stiffens from {c[0]:.0f} K to {c[1]:.0f} K "
                   f"(x{c[2]:.3f})")
    er, et = err300(v, el, "rc"), err300(v, el, "tap")
    if er is not None and et is not None and er > et + WORSE_PTS:
        why.append(f"300 K elastic error {er:.1f} % against {et:.1f} % for "
                   "the switched set")
    d = DISP.get(el)
    if d and d[0] > d[1] + WORSE_PTS:
        why.append(f"measured dispersion {d[0]:.1f} % against {d[1]:.1f} %")
    facts = {"alpha_1e6": a,
             "c11_md": ({"T_low": c[0], "T_high": c[1], "ratio": c[2]}
                        if c else None),
             "e300_rc": er, "e300_tap": et, "disp": d,
             "isf_rc": (r.get("stacking") or {}).get("isf"),
             "isf_tap": (t.get("stacking") or {}).get("isf")}
    if not why:
        return "rc", why, facts
    #  the fallback has to pass fit.py's compression test itself
    cp = t.get("compression") or {}
    tm = refdata.MELTING.get(el)
    if cp.get("basin") and tm and (cp.get("barrier") or 0.0) < KB * tm:
        facts["tap_compression"] = {"barrier_eV": cp.get("barrier"),
                                    "kT_melt_eV": KB * tm, "depth_eV": cp.get("depth")}
        why = why + ["and the switched set itself has no barrier against collapse "
                     f"under compression ({cp.get('barrier', 0):.3f} eV against "
                     f"k_B T_melt {KB * tm:.3f} eV)"]
        return "none", why, facts
    return "tap", why, facts


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    rec, keep, none = [], [], []
    for el in sorted(lib):
        v = lib[el]
        if not (isinstance(v, dict) and isinstance(v.get("rc"), dict)
                and "m" in v["rc"] and isinstance(v.get("tap"), dict)):
            continue
        s, why, facts = decide(el, v)
        v["md_recommended"] = {"set": s, "why": why, "rule": RULE,
                               "facts": facts}
        (rec if s == "rc" else none if s == "none" else keep).append((el, why))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"re-cut RECOMMENDED for MD ({len(rec)}): "
          + " ".join(e for e, _ in rec))
    print(f"switched set stays ({len(keep)}):")
    for el, why in keep:
        print(f"   {el:3s} {'; '.join(why)}")
    print(f"NO record fit for MD ({len(none)}):")
    for el, why in none:
        print(f"   {el:3s} {'; '.join(why)}")


if __name__ == "__main__":
    main()
