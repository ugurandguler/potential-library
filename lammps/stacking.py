#!/usr/bin/env python3
"""
The generalised stacking fault energy on (111), and the prediction it tests.

The ground-state result said seventy-two of seventy-six records prefer a
structure other than the fitted one, and for the face-centred metals the
preferred one is hexagonal close-packed.  An intrinsic stacking fault IS two
layers of hexagonal stacking inside a face-centred crystal, so that finding has
an immediate and falsifiable consequence:

    gamma_isf ~ 2 (E_hcp - E_fcc) / A

which, from the energies already measured, comes out NEGATIVE for twelve of
this library's thirteen face-centred records - copper at -63 mJ/m^2 against an
experimental +45, aluminium -124 against +166, platinum -32 against +322.  A
negative intrinsic fault energy means the perfect crystal is unstable against
spontaneous faulting: every close-packed plane would rather shear than stay.

That estimate is the first-neighbour one and it is not a measurement, which is
what this script is for.  The geometry, checked rather than assumed:

  the Shockley partial is a/6<112> and has to lie IN the fault plane, so it is
  a/6[11-2] and not a/6[112] - the latter is not perpendicular to [111].  In
  the oriented cell (x along [1-10], y along [11-2], z along [111]) that
  partial points along y, with length a/sqrt(6) against a y period of
  a*sqrt(6)/2, so the fault sits at exactly one third of the cell edge.

The cell is periodic, so the box tilt is shifted along with the atoms: without
that the wrap becomes a second interface carrying the opposite displacement,
which in a close-packed stacking is the forbidden arrangement and costs far
more than the fault itself.  Relaxation is allowed along the normal only;
letting the in-plane coordinates move lets the fault slide back out.

The published potentials go through the same code, and NIST publishes its own
900-point gamma-surface for most of them, so this has both an internal control
and an external one.

    python stacking.py Cu Al --sets tap
    python stacking.py                  # every fcc element, both sets
"""
import io
import json
import os
import re
import subprocess
import sys
#  PROCESSES, not threads.  Every job here spawns LAMMPS, and spawning a
#  subprocess means fork(), and fork() from one thread of a many-threaded
#  process copies the other threads' locks in whatever state they were in.
#  If any of them held the allocator lock at that instant the child inherits
#  a lock nobody will ever release.
#
#  This is what stopped the cluster job twice.  A stack dump of the second
#  attempt showed it exactly: one worker parked in subprocess._execute_child
#  at the fork, thirty-odd others parked in oriented_cell, which is numpy and
#  allocates.  Not a slow job - a dead one, at fourteen runs of eleven
#  hundred and fifty, with the interpreter reported as sleeping.
#
#  Removing the captured pipes (see surface.run) was a real fix for a real
#  second bug and did not touch this one.  The runs here are sub-second, so
#  the fork window opens hundreds of times a minute; the thermal expansion
#  job survives on the same machinery only because half an hour passes
#  between its forks.  A pool of processes has one thread each and cannot
#  reach the condition at all.
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))

LOCAL = not os.path.exists("/arf")
HOME = (subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                       capture_output=True, text=True).stdout.strip()
        if LOCAL else os.path.expanduser("~"))
os.environ.setdefault("LMP", f"{HOME}/lammps/src/lmp_serial")
os.environ.setdefault("BASEFILE",
                      os.path.join(ROOT, "standalone", "baselines.json"))

import numpy as np         # noqa: E402
import refdata            # noqa: E402
import elastic_T as E     # noqa: E402
import surface as SF      # noqa: E402  (oriented cell and data writer)

EV_A2_TO_MJ_M2 = 16021.766208
SKIN = 2.0
#  from 0 to the partial, in fractions of the y edge; the intrinsic fault is
#  the last point and the unstable fault is the maximum on the way
#  A FULL period of the y edge, not just up to the partial.  The stacking in
#  this cell advances by -y/3 going up, so shifting by +y/3 runs against the
#  sequence and produces the forbidden arrangement rather than the intrinsic
#  fault - a monotonic climb to 502 mJ/m^2 with no minimum anywhere.  Scanning
#  the whole period finds the fault wherever it actually is instead of assuming
#  which direction is which.
FRACS = [round(k / 24.0, 6) for k in range(25)]      # 0 .. 1

#  Scanning the whole period is right for finding the fault, but the two
#  numbers people quote live in the FIRST THIRD of it and taking them from the
#  whole scan gets both wrong:
#
#    gamma_isf is the fault, at one partial = 1/3 of the edge.  Reading the
#    LAST point instead makes it a whole lattice translation, which is
#    identically zero - copper came back at -0.0 mJ/m^2 for every potential in
#    the library, ours and the published ones alike, and a column of zeros
#    looks like a converged answer rather than a wrong question.
#
#    gamma_usf is the barrier ON THE WAY to that fault.  The maximum over the
#    whole period is the 2/3 point, where the shift stacks a layer directly on
#    top of its neighbour - the arrangement close packing forbids.  That is a
#    real energy but it is not the unstable stacking fault: it read 796 for
#    Mishin copper against a published 162.
#
#  So the window is [0, 1/3], and the rest of the scan is kept for two checks
#  it can make that the window cannot: the curve must close (gamma at a full
#  period is zero to within a fraction of a mJ/m^2, or the geometry is wrong)
#  and it must be symmetric about 2/3.
ISF_FRAC = 1.0 / 3.0


def reduce_curve(fr, gm):
    """the physical numbers, from the raw scan"""
    g = dict(zip([round(f, 6) for f in fr], gm))
    fr = sorted(g)
    i3 = min(fr, key=lambda f: abs(f - ISF_FRAC))
    if abs(i3 - ISF_FRAC) > 1e-3:
        return {"error": "1/3 noktasi taranmamis"}
    win = [f for f in fr if f <= i3 + 1e-9]
    fu = max(win, key=lambda f: g[f])
    out = {"isf": g[i3], "usf": g[fu], "usf_frac": fu,
           "back_barrier": g[fu] - g[i3]}
    #  A stacking fault is a MINIMUM of this curve.  If the point at one third
    #  is higher than its neighbours the shift went the wrong way and what has
    #  been measured is the forbidden arrangement, which is a perfectly good
    #  number for a completely different question.  Both times that happened
    #  the answer looked ordinary: no error, no warning, just a positive fault
    #  energy where a negative one belonged.
    nb = [f for f in fr if abs(f - i3) < 1.5 / 24.0 + 1e-9 and f != i3]
    out["is_minimum"] = bool(nb and all(g[f] >= g[i3] for f in nb))
    #  period closure: a shift of one whole edge is a lattice translation
    if abs(fr[-1] - 1.0) < 1e-6:
        out["closure"] = g[fr[-1]]
        out["closure_ok"] = bool(abs(g[fr[-1]]) < 1.0)
    return out

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       cell.data
pair_style      {style}
{coeff}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes

#  Everything above the mid-plane slides by the partial, and the box tilt
#  follows it.  Without the tilt the periodic wrap becomes a second interface
#  with the OPPOSITE displacement, which for a close-packed stacking is the
#  forbidden arrangement - two identical layers face to face - and costs far
#  more than the fault being measured.  Halving the sum then gave 419 mJ/m^2
#  where the reference potential's published value is 44.  Matching the tilt
#  to the shift makes the wrap continuous and leaves exactly one fault.
change_box      all triclinic
region          top block INF INF INF INF {zmid} INF units box
group           top region top
displace_atoms  top move 0.0 {dy} 0.0 units box
change_box      all yz delta {dy} units box

#  Relax along the normal only.  With the in-plane coordinates free the fault
#  simply slides back out and every answer is zero.
fix             1 all setforce 0.0 0.0 NULL
min_style       cg
minimize        1e-14 1e-12 20000 200000

variable        e equal pe
variable        aa equal lx*ly
print           "EFAULT ${{e}}"
print           "AREA ${{aa}}"
"""

PERFECT = IN.replace('''region          top block INF INF INF INF {zmid} INF units box
group           top region top
displace_atoms  top move 0.0 {dy} 0.0 units box

''', "").replace("EFAULT", "EPERFECT")


def one(job):
    el, tag, frac = job
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "sfruns", f"{el}_{safe}_{frac:.4f}")
    os.makedirs(d, exist_ok=True)
    try:
        style, pot = E.setup(d, el, tag)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    a0 = float(E.PACK[el]["a0"])
    struct = refdata.ELEMENTS[el]["struct"]
    #  Both close-packed structures are done by the same code, because the
    #  geometry works out to the same thing.  The oriented fcc cell has y
    #  along [11-2] with a period a*sqrt(6)/2 and the Shockley partial is
    #  a/sqrt(6); the orthogonal hcp cell has y along [01-10] with a period
    #  a*sqrt(3) and the basal partial is a/sqrt(3).  In both the partial is
    #  exactly ONE THIRD of the y edge, so the same scan, the same tilt and
    #  the same 1/3 read-off apply unchanged.
    #
    #  What differs is only what the fault means.  In fcc it is a slab of
    #  hexagonal stacking inside a cubic crystal; in hcp it is a slab of cubic
    #  stacking inside a hexagonal one.  The prediction therefore flips sign
    #  with the structure - 2(E_hcp-E_fcc)/A for fcc, 2(E_fcc-E_hcp)/A for hcp
    #  - which is what makes running both a stronger test than running either.
    def build(reps):
        if struct == "hcp":
            coa = float(E.PACK[el].get("c_over_a")
                        or refdata.ELEMENTS[el]["c_over_a"])
            return SF.hcp_cell(a0, coa, tuple(reps), "z")
        return SF.oriented_cell("fcc", a0, "111", tuple(reps))

    _, base = build((1, 1, 1))
    rc = 8.6
    reps = [max(1, int(np.ceil(2.4 * (rc + SKIN) / b))) for b in base]
    #  thick enough that the two faults do not see each other
    reps[2] = max(reps[2], int(np.ceil(2.0 * 2.0 * (rc + SKIN) / base[2])), 6)
    coords, box = build(reps)
    SF.write_cell(os.path.join(d, "cell.data"), coords, box,
                  refdata.MASSES[el])
    coeff = E.coeff_line(style, pot, el)
    #  The mid-plane goes between two layers, not through one - and WHICH two
    #  matters for hcp in a way it does not for fcc.
    #
    #  Face-centred stacking is ABCABC, and the three interfaces A|B, B|C and
    #  C|A are related by a cyclic relabelling, so a partial shift produces the
    #  same intrinsic fault whichever one is cut.  Hexagonal stacking is ABAB,
    #  where A|B and B|A are related by REVERSING the in-plane direction, so
    #  cutting the wrong one flips which way the partial has to point.
    #
    #  The number of layers depends on how many cells were stacked to make the
    #  slab thick enough, which varies with c - so the parity of the middle
    #  layer, and with it the answer, changed from element to element.  Seven
    #  of twelve hexagonal elements came back with the fault energy replaced by
    #  the forbidden arrangement: rhenium at +2531 mJ/m^2, titanium at +805.
    #  They were caught by the local-minimum test rather than by inspection.
    zs = np.unique(np.round(coords[:, 2], 6))
    mi = len(zs) // 2
    if struct == "hcp" and mi % 2 == 0:
        mi -= 1
    mid = float(zs[mi] - 0.25 * (zs[1] - zs[0]))
    #  The shift is a fraction of ONE oriented cell's y edge, not of the
    #  repeated box.  Taking it from box[1] made "one third" equal two whole
    #  lattice periods in a six-cell box, so the displacement was a lattice
    #  translation and every fault energy came back as exactly zero - the
    #  curve had zeros at 0, 1/6 and 1/3, which is 0, 1 and 2 whole periods.
    #  base[1]/3 is the Shockley partial a/sqrt(6); that identity is asserted
    #  rather than trusted.
    part = base[1] / 3.0
    want = a0 / np.sqrt(3.0 if struct == "hcp" else 6.0)
    assert abs(part - want) < 1e-6, (
        f"kismi vektor boyu tutmuyor: {part:.6f} != {want:.6f} ({struct})")
    #  The partial has a direction as well as a length, and the two structures
    #  need opposite ones in these cells.  Getting it wrong does not fail, it
    #  produces a curve that climbs to a MAXIMUM where the fault should be:
    #  magnesium came back at +316 mJ/m^2 at one third with its real fault,
    #  -31, sitting at two thirds.  The same mistake was made once already on
    #  the cubic side.  The read-off is checked afterwards rather than trusted
    #  - see reduce_curve, which now requires the fault to be a local minimum.
    sgn = -1.0 if struct == "hcp" else 1.0
    txt = (PERFECT if frac == 0.0 else IN).format(
        style=style, coeff=coeff, skin=SKIN, zmid=mid,
        dy=sgn * frac * 3.0 * part)
    lg = SF.run(d, txt, "in.sf")
    key = "EPERFECT" if frac == 0.0 else "EFAULT"
    e = SF.grab(lg, key)
    area = SF.grab(lg, "AREA")
    if e is None or area is None:
        err = [l for l in lg.splitlines() if "ERROR" in l]
        return {"error": err[0][:70] if err else "the run did not finish"}
    return {"E": e, "area": area, "atoms": len(coords),
            "struct": struct, "plane": "0001" if struct == "hcp" else "111"}


def redo(fp):
    """Recompute isf/usf from the raw scan already in a results file.

    The cluster job was launched before the window was corrected, so it writes
    the old fields.  The scan itself is unaffected - the same 25 points, the
    same energies - so the fix is a re-read rather than a re-run.
    """
    d = json.load(open(fp))
    n = 0
    print(f"{'record':40s}{'isf':>9s}{'usf':>9s}{'closure':>9s}")
    print("-" * 67)
    for k, v in sorted(d.items()):
        if "frac" not in v or "gamma" not in v:
            continue
        r = reduce_curve(v["frac"], v["gamma"])
        if "error" in r:
            print(f"{k[:40]:40s}   {r['error']}")
            continue
        v.update(r)
        n += 1
        print(f"{k[:40]:40s}{r['isf']:9.1f}{r['usf']:9.1f}"
              f"{r.get('closure', float('nan')):9.2f}"
              f"{'  KAPANMA BOZUK' if r.get('closure_ok') is False else ''}"
              f"{'  1/3 MINIMUM DEGIL' if r.get('is_minimum') is False else ''}")
    tmp = fp + ".tmp"
    json.dump(d, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, fp)
    print(f"\n{n} records reduced again -> {fp}")


def main():
    argv = sys.argv[1:]
    if "--reduce" in argv:
        i = argv.index("--reduce")
        fp = (argv[i + 1] if len(argv) > i + 1 and
              not argv[i + 1].startswith("--")
              else os.path.join(HERE, "stacking.json"))
        return redo(fp)
    sets = ["tap"]
    if "--sets" in argv:
        i = argv.index("--sets")
        sets = argv[i + 1].split(",")
        del argv[i:i + 2]
    with_base = "--nobase" not in argv
    els = [a for a in argv if not a.startswith("--")]
    els = els or sorted(E.PACK)
    #  Say what is being left out and why, rather than filtering in silence.
    #  A stacking fault in the sense measured here is a slab of the OTHER
    #  close-packed stacking, and body-centred cubic has no such thing: its
    #  {110} and {112} gamma-lines rise to a maximum and come back down with
    #  no minimum in between, so there is an unstable fault and no intrinsic
    #  one.  Running the same scan on those thirteen elements would return a
    #  number that looks like the others and is not the same quantity.
    skipped = sorted(e for e in els
                     if refdata.ELEMENTS[e]["struct"] not in ("fcc", "hcp"))
    els = [e for e in els
           if refdata.ELEMENTS[e]["struct"] in ("fcc", "hcp")]
    if skipped:
        print(f"skipped {len(skipped)} body-centred elements "
              f"(intrinsic stacking fault undefined): {' '.join(skipped)}")

    jobs = []
    for el in els:
        for tag in sets:
            if tag in E.PACK.get(el, {}):
                jobs += [(el, tag, f) for f in FRACS]
        if with_base:
            for fn, _ in E.BASE.get(el, []):
                if all(os.path.exists(os.path.join(E.POTDIR, g))
                       for g in fn.split("+")):
                    jobs += [(el, "base|" + fn, f) for f in FRACS]
    if not jobs:
        print("0 runs - fcc element or set name not recognised")
        return
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"{len(jobs)} runs, {nw} parallel processes", flush=True)
    with ProcessPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(one, jobs, chunksize=1))

    got = {}
    for (el, tag, f), r in zip(jobs, res):
        got.setdefault(f"{el}|{tag}", {})[f"{f:.6f}"] = r

    out = {}
    print()
    print(f"{'el':4s}{'source':>28s}{'gamma_isf':>11s}{'gamma_usf':>11s}"
          f"{'expt':>8s}   note")
    print("-" * 68)
    for key, series in sorted(got.items()):
        el, tag = key.split("|", 1)
        if any("error" in r for r in series.values()):
            e = next(r["error"] for r in series.values() if "error" in r)
            print(f"{el:4s}{tag[:28]:>28s}   {e[:40]}")
            continue
        f0 = series["0.000000"]
        area = f0["area"]
        #  ONE fault now, not two
        g = {float(k): (v["E"] - f0["E"]) / area * EV_A2_TO_MJ_M2
             for k, v in series.items()}
        fr = sorted(g)
        rec = {"frac": fr, "gamma": [g[f] for f in fr], "atoms": f0["atoms"]}
        rec.update(reduce_curve(fr, [g[f] for f in fr]))
        out[key] = rec
        isf, usf = rec.get("isf"), rec.get("usf")
        lab = {"tap": "MAU", "tap_ug": "UG"}.get(tag, tag.replace("base|", ""))
        note = "ISF NEGATIF" if isf is not None and isf < 0 else ""
        if rec.get("closure_ok") is False:
            note = (note + " KAPANMA BOZUK %.1f" % rec["closure"]).strip()
        if rec.get("is_minimum") is False:
            note = (note + " 1/3 MINIMUM DEGIL - YON YANLIS").strip()
        print(f"{el:4s}{lab[:28]:>28s}{isf:11.1f}{usf:11.1f}"
              f"{'':>8s}   {note}")

    fp = os.path.join(HERE, "stacking.json")
    old = {}
    if os.path.exists(fp):
        try:
            old = json.load(open(fp))
        except Exception:
            old = {}
    old.update(out)
    tmp = fp + ".tmp"
    json.dump(old, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, fp)
    print(f"\n{len(out)} records written; file now holds {len(old)}  -> {fp}")


if __name__ == "__main__":
    main()
