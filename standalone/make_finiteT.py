#!/usr/bin/env python3
"""finiteT.json - the table the interface's finite-temperature panel draws.

Written to a file rather than embedded in make_gui.py so that a new run can
be folded in by re-running this and rebuilding the page, without touching the
renderer.  Caesium landed on 2026-08-25 and did exactly that: its 8^3 box
could not resolve its own answer - a floor of 11.2 % under a difference of
8.3 - and the 21^3 box brings the floor to 0.8 %, which turns the row from
unresolved into a measured gain.  It is the clearest case in the set for the
floor being computed rather than eyeballed.

WHICH ARM.  Every number here is the MAU arm - the tapered record,
`lib[el]["tap"]`, shipped as `<El>_taper.ugur`.  That is not the arm the
parameter-set selector calls MAU: its "MAU" is the HARD-TRUNCATED root
record.  The two are different potentials and the panel says so, because a
reader comparing the curve on screen with a number in this table would
otherwise be comparing two different things.  The switched arm is used here
because it is the only one that can be run at temperature at all - a hard cut
leaves the pair energy discontinuous and copper drifts 350 meV/atom/ns.

THE FLOOR IS THE POINT.  Every q is scored twice, at q and at a symmetry
image of q, and the spread between them is what the run cannot resolve.  A
difference between columns smaller than that has not been measured, so the
verdict is computed FROM the floor rather than left to the reader's eye.  The
page counts them from this file rather than quoting a number.

NVT, NOT NPT, SINCE 2026-10-01.  fix phonon reads unwrapped CARTESIAN
positions, so under a barostat the affine box strain enters <u u*> on the
q-lines along the primitive axes and softens the longitudinal branch there:
every hcp element and every bcc one was hit (hcp LA and LO split by 15-25 % at
A, where symmetry makes them degenerate), the fcc rows were not (no scored or
drawn q lies on an axis).  The twelve affected elements were re-run in NVT at
the MEAN cell of their own zero-pressure NPT run, same seed, steps and mesh;
mdres.json keeps the NPT numbers under `npt_2026_08`.  The leak had flattered
the hcp scores and inflated their floors two- to thirteen-fold.

FIVE ELEMENTS LEFT THE STUDY at the same time: the switched records of Cs, Fe,
Nb, Rb and W do not stay body-centred cubic at the measurement's temperature
(the structure gate, add_gate.py), so what was measured was not a bcc
spectrum.  They are kind "notbcc" below.

    python make_finiteT.py        # -> finiteT.json
"""
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import curve_mae as CM                      # noqa: E402

#  Elements outside the study, and why.  A blank cell reads as a failure, so
#  every one of them carries its own sentence.  Four groups, and they are not
#  the same kind of absence: a missing measurement is not a failed potential.
OUT = {
    "cold": ("Al", "Be", "K", "Pt"),
    "screen": ("Li",),
    "model": ("Co", "Cr", "Mo", "Rh", "Sr"),
    "none": (),
    #  Ir, Ru and Re were "pending" and "gamma" from 2026-09-13 until their
    #  runs landed on 2026-09-14; they are rows now.  The two kinds stay so a
    #  later entry can use them.
    "pending": (),
    "gamma": (),
    #  vanadium, 2026-09-15: a MEASURED dispersion (Bosak et al. 2008, x-ray)
    #  replaced the Colella model as its reference, and no fix phonon run has
    #  been made against it - not waiting on a job, simply not run
    "unrun": ("V",),
    #  2026-10-01: the switched record leaves bcc at the run temperature - the
    #  spectrum that was measured belongs to a different structure
    "notbcc": ("Cs", "Fe", "Nb", "Rb", "W"),
}
OUTSIDE = {el for g in OUT.values() for el in g}


#  rows whose run measured one Gamma frequency, scored on that alone at 0 K
#  as well (rhenium; its [001] dispersion entered the library 2026-10-01 and
#  is not what this run measured)
GAMMA_ONLY = ("Re",)


def zero_k(els, use_a_meas):
    """{el: percent} for the tapered arm at one choice of volume"""
    CM.USE_A_MEAS, CM.NEAREST = use_a_meas, False
    CM.ONLY_POINTS = set(GAMMA_ONLY)
    buf, keep = io.StringIO(), sys.stdout
    sys.stdout = buf
    try:
        CM.main(tuple(els))
    finally:
        sys.stdout = keep
    out = {}
    for ln in buf.getvalue().splitlines():
        p = ln.split()
        if len(p) >= 7 and p[1].endswith("K") and p[-1].endswith("%"):
            try:
                out[p[0].rstrip("*")] = float(p[-1].rstrip("%"))
            except ValueError:
                continue
    return out


def main():
    lib = json.load(open(os.path.join(HERE, "library.json")))
    md = json.load(open(os.path.join(HERE, "mdres.json")))
    els = sorted(md)
    a0, am = zero_k(els, False), zero_k(els, True)

    rows = {}
    for el in els:
        m = md[el]
        if "error" in m or el in OUTSIDE:
            continue
        if el not in a0 or el not in am:
            #  curve_mae prints "-" in the MAU column when the tapered arm is
            #  absent, and the parse above then returns nothing.  A clean clone
            #  is always in that position: `tap` comes from add_taper_overlay
            #  reading the dense_*.json search output, which is not published.
            #  Saying so beats KeyError, which was what this did.
            raise SystemExit(
                "no 0 K column for %s: library.json has no tapered arm.\n"
                "  `tap` is written by add_taper_overlay.py from the "
                "dense_*.json\n  search output, which is not part of the "
                "public repository.\n"
                "  The finished table ships as finiteT.json and the page "
                "already carries it." % el)
        best = min(a0[el], am[el])
        d = best - m["pct"]
        fl = m["sc_pct"]
        rows[el] = {
            "n": m["n"], "mesh": m.get("mesh", "?"), "steps": m.get("steps", 0),
            "T": lib[el]["exp_curve"]["T_K"],
            "a0": round(a0[el], 1), "am": round(am[el], 1),
            "ft": round(m["pct"], 1), "floor": round(fl, 1),
            "d": round(d, 1),
            "v": "flat" if abs(d) < fl else ("gain" if d > 0 else "loss"),
            #  which ensemble the spectrum was sampled in: "nvt" at the mean
            #  NPT cell (re-run 2026-10-01), "npt" otherwise (fcc and rhenium,
            #  which the barostat leak does not reach)
            "box": m.get("box", "npt"),
        }
        #  a row scored differently from the rest says how, in its own words
        #  (rhenium: one Gamma frequency, floored by its degenerate split)
        if m.get("note"):
            rows[el]["note"] = m["note"]
        #  what the barostat run had said, so the page can show both
        old = m.get("npt_2026_08")
        if old:
            rows[el]["npt"] = {"ft": round(old["pct"], 1),
                               "floor": round(old["sc_pct"], 1)}

    why = {}
    for kind, group in OUT.items():
        for el in group:
            v = lib.get(el)
            if not isinstance(v, dict) or "m" not in v:
                continue
            why[el] = {"kind": kind,
                       "T": (v.get("exp_curve") or {}).get("T_K")}

    out = {"rows": rows, "out": why}
    p = os.path.join(HERE, "finiteT.json")
    json.dump(out, open(p, "w"), indent=0, sort_keys=True)
    ver = {}
    for r in rows.values():
        ver[r["v"]] = ver.get(r["v"], 0) + 1
    print("finiteT.json: %d elements measured, %d out of scope" % (len(rows),
                                                                len(why)))
    print("   verdicts:", ver)
    print("   %d bytes" % os.path.getsize(p))


if __name__ == "__main__":
    main()
