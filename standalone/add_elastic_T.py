#!/usr/bin/env python3
"""
Fold the finite-temperature elastic sweeps into the library, with the derived
mechanical properties already computed.

The page plots; it does not do physics.  Voigt-Reuss-Hill averaging has a
different closed form for cubic and hexagonal symmetry, the bulk modulus has no
Voigt-Reuss spread for cubic at all, and Chen's hardness has an exponent that
is easy to mistype - none of that belongs in a template literal inside an HTML
file where it cannot be tested.  So it is done here and the page receives
numbers.

Two flags travel with every point, because both change what a curve means:

  above_melt   the run sits above the element's melting point.  A perfect
               108-atom crystal has nowhere to nucleate from and will happily
               superheat, so those points describe a metastable solid.  The
               hexagonal identity C66 = (C11-C12)/2 confirms it independently:
               below melting it holds to 1.4 %, above it degrades to 67 %.
  nudge_bad    the parameter set does not hold its lattice against a 1e-5 A
               displacement (see lammps/jiggle_test.py).  Thermal motion at any
               temperature here is four orders of magnitude larger, so whatever
               was measured is a property of some other structure.  Plotted,
               but never plotted as if it were sound.

    python add_elastic_T.py            # every source it can find
    python add_elastic_T.py --dry
"""
import json
import math
import os
import sys

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
#  (file, whether its grid follows the melting point)
#  one file per cluster job, because three of them run at once and a shared
#  read-then-write merge loses whichever finishes second
SOURCES = [
    (os.path.join(ROOT, "lammps", "elastic_T.json"), True),
    (os.path.join(ROOT, "lammps", "elastic_T_re.json"), True),
    (os.path.join(ROOT, "lammps", "elastic_T_newbase.json"), True),
    (os.path.join(ROOT, "lammps", "elastic_T_meam.json"), True),
    (os.path.join(ROOT, "lammps", "elastic_T_sheng.json"), True),
    (os.path.join(ROOT, "lammps", "hcp_born.json"), False),
]
KS = ("C11", "C12", "C13", "C33", "C44", "C66")
_lp = os.path.join(HERE, "baseline_labels.json")
LABELS = json.load(open(_lp)) if os.path.exists(_lp) else {}
_bp = os.path.join(HERE, "baselines.json")
BASELINES = json.load(open(_bp)) if os.path.exists(_bp) else {}


def cubic(c):
    """closed-form VRH for cubic symmetry; B has no Voigt-Reuss spread"""
    c11, c12, c44 = c["C11"], c["C12"], c["C44"]
    B = (c11 + 2 * c12) / 3.0
    gv = (c11 - c12 + 3 * c44) / 5.0
    den = 4 * c44 + 3 * (c11 - c12)
    gr = 5 * c44 * (c11 - c12) / den if abs(den) > 1e-9 else float("nan")
    return B, B, gv, gr


def hexag(c):
    """closed-form VRH for hexagonal symmetry"""
    c11, c12, c13, c33, c44 = (c["C11"], c["C12"], c["C13"],
                               c["C33"], c["C44"])
    c66 = 0.5 * (c11 - c12)
    bv = (2 * (c11 + c12) + 4 * c13 + c33) / 9.0
    gv = (c11 + c12 + 2 * c33 - 4 * c13 + 12 * c44 + 12 * c66) / 30.0
    C2 = (c11 + c12) * c33 - 2 * c13 * c13
    d1 = c11 + c12 + 2 * c33 - 4 * c13
    br = C2 / d1 if abs(d1) > 1e-9 else float("nan")
    d2 = 3 * bv * c44 * c66 + C2 * (c44 + c66)
    gr = 2.5 * C2 * c44 * c66 / d2 if abs(d2) > 1e-9 else float("nan")
    return bv, br, gv, gr


def derive(struct, c):
    bv, br, gv, gr = hexag(c) if struct == "hcp" else cubic(c)
    B, G = 0.5 * (bv + br), 0.5 * (gv + gr)
    out = {"B": B, "G": G}
    if not (B > 0 and G > 0):
        #  a negative average is not a modulus; it is the Born criteria being
        #  violated, and squaring it into a hardness would hide that
        out.update({k: None for k in ("E", "nu", "BG", "AU", "Hv")})
        return out
    out["E"] = 9 * B * G / (3 * B + G)
    out["nu"] = (3 * B - 2 * G) / (2 * (3 * B + G))
    out["BG"] = B / G
    out["AU"] = (5 * gv / gr + bv / br - 6.0) if gr else None
    k = G / B
    out["Hv"] = 2 * (k * k * G) ** 0.585 - 3
    return out


def stable(struct, c):
    """Born criteria - a curve that crosses zero here has stopped being elastic"""
    if struct == "hcp":
        return (c["C11"] > abs(c["C12"])
                and (c["C11"] + c["C12"]) * c["C33"] > 2 * c["C13"] ** 2
                and c["C44"] > 0)
    return (c["C11"] - c["C12"] > 0 and c["C44"] > 0
            and c["C11"] + 2 * c["C12"] > 0)


def main():
    dry = "--dry" in sys.argv
    lib = json.load(open(os.path.join(HERE, "library.json")))
    jig = {}
    p = os.path.join(ROOT, "lammps", "jiggle_test.json")
    if os.path.exists(p):
        jig = json.load(open(p))

    got = {}
    for path, melt_grid in SOURCES:
        if not os.path.exists(path):
            continue
        raw = json.load(open(path))
        for key, series in raw.items():
            el, tag = key.split("|", 1)
            if el not in lib:
                continue
            struct = refdata.ELEMENTS[el]["struct"]
            tm = refdata.MELTING.get(el, 1e9)
            pts = []
            for Ts, r in series.items():
                if not r or "error" in r:
                    continue
                if any(r.get(k) is None for k in ("C11", "C12", "C44")):
                    continue
                c = {k: r.get(k, 0.0) for k in KS}
                if struct != "hcp":
                    c["C33"], c["C13"] = c["C11"], c["C12"]
                d = derive(struct, c)
                pts.append({"T": float(Ts), **{k: c[k] for k in KS}, **d,
                            "above_melt": bool(float(Ts) > tm),
                            "born_ok": bool(stable(struct, c)),
                            "Tavg": r.get("Tavg")})
            if not pts:
                continue
            pts.sort(key=lambda x: x["T"])
            nb = jig.get(f"{el}|{tag}")
            rec = {"pts": pts, "grid": "Tmelt" if melt_grid else "fixed",
                   "Tmelt": tm if tm < 1e8 else None,
                   "nudge_bad": (nb is not None and nb.get("ok") is False)}
            #  hcp_born.py wrote the baseline tag as a bare "base" without
            #  the file name, so the curve was labelled, literally, "base".
            #  The element has exactly one shipped potential in that sweep, so
            #  it can be resolved; elastic_T.py carries the name in the tag.
            if tag == "base" and el in BASELINES and len(BASELINES[el]) == 1:
                tag = "base|" + BASELINES[el][0][0]
            if tag.startswith("base|"):
                #  the file name is not a label.  "Mg_mm" tells a reader
                #  nothing; the attribution read out of the file's own CITATION
                #  header does, and it is the reference the curve is being
                #  measured against, so it should be nameable.
                fn = tag.split("|", 1)[1]
                meta = LABELS.get(fn, {})
                rec["label"] = meta.get("label", fn)
                rec["file"] = fn
                rec["citation"] = meta.get("citation", "")
                rec["kind"] = "base"
            else:
                rec["label"] = {"tap": "MAU", "tap_ug": "UG"}.get(tag, tag)
                rec["kind"] = "ours"
            #  a Tmelt-scaled run supersedes a fixed-grid one for the same
            #  series; the fixed grid only exists for the ruthenium comparison
            slot = got.setdefault(el, {})
            if tag not in slot or melt_grid:
                slot[tag] = rec
            slot.pop("base", None) if tag != "base" and "base" in slot and                 slot.get("base", {}).get("label") == "base" else None

    #  AFLOW's AEL is a single density-functional point, not a curve: Voigt-
    #  Reuss-Hill B and G at 0 K and nothing else.  It is attached separately
    #  and drawn as a marker, because presenting one DFT number as if it were
    #  a temperature series would be a lie about what it is.  Entries in the
    #  wrong phase are carried but flagged - an elastic modulus belongs to a
    #  structure, and AFLOW's barium is hcp where this library's is bcc.
    #  WITHDRAWN FROM THE PUBLISHED PAGE - do not re-enable without reading
    #  this.  AFLOW's terms are "free for scientific, academic and
    #  non-commercial purposes.  Any other use is prohibited"
    #  (aflowlib.duke.edu).  The library page is distributed under CC-BY 4.0,
    #  which permits commercial use, so embedding AFLOW values in it would
    #  grant a right we do not hold.  Materials Project (CC-BY 4.0), MC3D
    #  (CC-BY 4.0) and JARVIS (US Government work) are all compatible and stay.
    #
    #  The cost of dropping it was measured first: 22 elements carried an AFLOW
    #  record, only 14 were usable, and Materials Project already supplies
    #  elastic data for 33.  AFLOW added K and Na and nothing else.
    #
    #  Fetching it for our own checking is still fine - that is the use AFLOW
    #  permits.  What may not happen is the result going into the page.

    n = 0
    for el, series in got.items():
        lib[el]["elasticT"] = series
        n += len(series)
    print(f"{len(got)} element, {n} egri")
    bad = [f"{el}/{t}" for el, s in got.items() for t, r in s.items()
           if r["nudge_bad"]]
    print(f"failed the nudge test and flagged: {' '.join(bad) or 'none'}")
    am = sum(1 for s in got.values() for r in s.values()
             for q in r["pts"] if q["above_melt"])
    nb = sum(1 for s in got.values() for r in s.values()
             for q in r["pts"] if not q["born_ok"])
    print(f"erime ustu nokta: {am}   Born olcutunu ihlal eden nokta: {nb}")
    if dry:
        print("--dry: nothing written")
        return
    json.dump(lib, open(os.path.join(HERE, "library.json"), "w"),
              indent=1, sort_keys=True)
    print("-> library.json")


if __name__ == "__main__":
    main()
