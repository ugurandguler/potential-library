#!/usr/bin/env python3
"""
Merge every independent search into one best-of result.

fit.py and dense_fit.py optimise exactly the same objective, but each starts
from different points and the constraint solver converges to different branches
depending on where it starts.  Neither run dominates the other: the dense local
search beat the plain fit on Ta, V, W, Au and Mg, and lost to it on Li, Mo and
Ti.  They are independent samples of one landscape, so the right answer per
element is the best anyone found.

"Best" is not the lowest score, though.  Some runs imposed dynamical stability
and some did not, and an unconstrained fit can win on score while having
imaginary modes across a quarter of the zone (Li: 7.6 % with 24 % imaginary,
against 37.0 % for the stable one).  So stable candidates outrank unstable ones
and scores are only compared within a group.  Stability is re-measured here on
the 8^3 U 9^3 union screen rather than trusted from whatever the run used.

One further rule: **every winner must have been fitted at the same three-body
cutoff.**  Ranking on score alone let ten elements into the library at the old
first-shell truncation, because at 1.12 they scored a little better - chromium
49.6 % against 55.3 %, rhodium 7.24 against 7.39 - and the library then mixed two
different potentials.  A cutoff picked per element is a free discrete parameter,
not a prescription, and it also breaks the MAU-against-UG comparison, which only
means anything when both sides truncate at the same radius.  So candidates are
filtered to RCUT3_OVER_DNN and the rejected ones are reported.

The cutoff cannot be read off the record, either.  Until fit.py was fixed it
wrote 1.12 into every record whatever the run used, so `true_rcut3` recovers it
by rebuilding at each candidate radius and keeping the one that reproduces the
record's own elastic constants.

    python merge_fits.py            # fit.json + dense_*.json in this directory
    python merge_fits.py extra/     # also pull in results from another directory

Writes fit.json in place (backing up the previous one) and reports where each
winner came from.
"""
import glob, json, os, shutil, sys

import numpy as np

import build_library as BL
import fit as F
import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = float(os.environ.get("RCUT3_OVER_DNN", "1.50"))
KNOWN_FACTORS = (1.50, 1.12)
CIJ_TOL = 2e-3
IDX = {"C11": (0, 0), "C12": (0, 1), "C13": (0, 2),
       "C33": (2, 2), "C44": (3, 3)}


def true_rcut3(el, rec, cache={}):
    """
    The cutoff factor the record was really fitted at.

    Its own Cij is the evidence: rebuilding at the wrong radius moves aluminium
    from C11 = 108.2 GPa to 48.4.  Where the three-body coefficient is small the
    two radii give the same constants - beryllium 192.9 against 191.9 - and
    there the measurement cannot decide.

    In that case the record's own label stands.  Refusing instead threw away
    scandium's three best fits, all of them made at 1.50 by a job script that
    set the radius explicitly, because its C = -0.062 makes the three-body term
    too weak for the two radii to be told apart.  That is the point: when the
    radii cannot be distinguished from the elastic constants they describe the
    same potential, so which label is attached has no physical consequence.
    The measurement exists to catch a label that is *wrong*, and a wrong label
    is only wrong when the two differ - which is exactly when the measurement
    can see it.
    """
    want = rec.get("Cij")
    if not want or "dnn" not in rec:
        return None
    key = (el, rec["m"], rec["gamma"], rec["C"], rec["alpha3"], rec["r0"],
           rec["D"], rec["alpha"])
    if key in cache:
        return cache[key]
    e = refdata.ELEMENTS[el]
    dev = []
    for f in KNOWN_FACTORS:
        q = dict(rec, rcut3=rec["dnn"] * f)
        try:
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                cry, pot = BL.build(el, q, e)
                C, _ = L.elastic(cry, pot)
        except (OverflowError, ValueError, FloatingPointError,
                np.linalg.LinAlgError):
            continue
        dev.append((max(abs(float(C[IDX[k]]) - v) / max(abs(v), 1.0)
                        for k, v in want.items() if k in IDX), f))
    dev.sort()
    if not dev or dev[0][0] > CIJ_TOL:
        out = None                       # no radius reproduces the record
    elif len(dev) > 1 and dev[1][0] <= CIJ_TOL:
        out = rec["rcut3"] / rec["dnn"]   # indistinguishable; the label stands
    else:
        out = dev[0][1]
    cache[key] = out
    return out


def sources(extra_dirs):
    """(label, {element: record}) for every result file we can find"""
    out = []
    p = os.path.join(HERE, "fit.json")
    if os.path.exists(p):
        out.append(("plain", json.load(open(p))))
    dirs = [HERE] + [d if os.path.isabs(d) else os.path.join(HERE, d)
                     for d in extra_dirs]
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "dense_*.json"))):
            try:
                data = json.load(open(f))
            except ValueError:
                continue
            if "failed" in data:
                continue
            #  same basenames appear in every run directory, so keep the
            #  directory in the label or the origin report is ambiguous
            tag = os.path.join(os.path.basename(d) or ".",
                               os.path.basename(f)[:-5])
            out.append((tag, data))
    return out


def stable(el, rec, cache={}):
    """dynamically stable on the union screen, 4^3 reject then 8^3 and 9^3"""
    key = (el, rec["m"], rec["gamma"], rec["C"], rec["alpha3"], rec["r0"])
    if key in cache:
        return cache[key]
    e = refdata.ELEMENTS[el]
    try:
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(rec)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            ok = F.dynamically_stable(cry, pot)
    except (KeyError, OverflowError, ValueError, FloatingPointError,
            np.linalg.LinAlgError):
        ok = False
    cache[key] = ok
    return ok


def main(extra_dirs):
    srcs = sources(extra_dirs)
    if not srcs:
        raise SystemExit("no result files found")

    #  Lowest score alone is the wrong rule once runs differ in whether they
    #  imposed dynamical stability.  Li's unconstrained fit scores 7.6 % and has
    #  24 % imaginary modes; the constrained one scores 37.0 % and is real.  A
    #  plain best-of hands the library the unphysical one.  So: rank stable
    #  candidates ahead of unstable ones, and only compare scores within a
    #  group.  Elements where nothing stable exists are kept but flagged.
    best, origin = {}, {}
    rejected, ambiguous = {}, {}
    for label, data in srcs:
        for el, rec in data.items():
            if not isinstance(rec, dict) or "score" not in rec:
                continue
            f = true_rcut3(el, rec)
            if f is None:
                ambiguous[el] = ambiguous.get(el, 0) + 1
                continue
            if abs(f - TARGET) > 1e-6:
                rejected[el] = min(rejected.get(el, 9e9), rec["score"])
                continue
            #  the record's own label may be the pre-fix 1.12; the measured one
            #  is what everything downstream rebuilds from
            rec = dict(rec, rcut3=rec["dnn"] * f)
            rec = dict(rec, dyn_stable=stable(el, rec))
            cur = best.get(el)
            better = (cur is None or
                      (rec["dyn_stable"], -rec["score"]) >
                      (cur["dyn_stable"], -cur["score"]))
            if better:
                best[el], origin[el] = rec, label

    if rejected:
        print(f"three-body cutoff fixed at {TARGET:.2f} d_nn; dropped fits made "
              f"at another radius:")
        for el in sorted(rejected):
            here = best.get(el, {}).get("score")
            got = f"{here*100:.2f} %" if here is not None else "nothing kept"
            print(f"   {el:3s} best rejected {rejected[el]*100:6.2f} %  ->  {got}")
        print()
    if ambiguous:
        print("cutoff not recoverable from the record (the two radii give the "
              "same elastic constants), so these candidates were skipped:")
        print("   " + ", ".join(f"{el}x{n}" for el, n in sorted(ambiguous.items())))
        print()

    unstable = sorted(el for el, r in best.items() if not r["dyn_stable"])
    if unstable:
        print(f"no dynamically stable fit found for: {', '.join(unstable)}")
        print("   (kept the best-scoring one, flagged dyn_stable = false)\n")

    dst = os.path.join(HERE, "fit.json")
    if os.path.exists(dst):
        shutil.copy(dst, dst + ".bak")
    json.dump(best, open(dst, "w"), indent=1, sort_keys=True)

    from collections import Counter
    tally = Counter(origin.values())
    print(f"merged {len(srcs)} result sets -> {len(best)} elements")
    print("winner came from:")
    for k, n in tally.most_common():
        print(f"   {k:22s} {n:3d} element(s)")
    print()
    print(f"{'el':4s}{'best':>8s}   from")
    for el in sorted(best, key=lambda e: best[e]["score"]):
        print(f"{el:4s}{best[el]['score']*100:7.1f}%   {origin[el]}")
    missing = sorted({e for _, d in srcs for e in d} - set(best))
    if missing:
        print(f"\nno usable fit: {', '.join(missing)}")


if __name__ == "__main__":
    main(sys.argv[1:])
