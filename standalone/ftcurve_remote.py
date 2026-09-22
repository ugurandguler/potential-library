#!/usr/bin/env python3
"""The measured dispersion along the path, for every element, on the cluster.

Runs next to the binaries because the eigenvector files never travelled well;
only the curve comes home, and it comes as one file fetched in pieces.

INSIDE THE FIRST MESH CELL nothing was measured.  The mesh is fixed by the
box, and between Gamma and the first mesh point phana is interpolating across
a gap - badly enough, for the softest branches, to return negative
frequencies.  There omega = v |q| holds exactly, so the value is taken from a
larger q along the same direction and divided by the ratio.  That is the
physics rather than a fit to it and it cannot return a negative number.

    For a CUBIC crystal every branch is acoustic and the substitution is
    valid for all of them.  For hcp it is valid only for the three acoustic
    branches: an optic branch does not go to zero at Gamma and scaling it by
    |q| is simply wrong - lutetium's measured TOperp is 3.670 THz there, the
    interpolant returns 3.551 and is right, the acoustic scaling would return
    1.060.  So hcp keeps the interpolated optic branches and substitutes only
    the lowest three.  merge_ftcurve.py reshapes them afterwards.

ONE MESH INTERVAL, NOT TWO.  The earlier run used two, on the grounds that
thallium still went negative a little past the first mesh point.  Two throws
away the mesh point AT one interval, which is a measurement - and where the
mesh is coarse that is the whole curve.  Niobium's mesh is 8 and P sits at
(1/4,1/4,1/4), exactly two intervals from Gamma in every direction, so a
radius of two replaced the entire Gamma-P branch with a straight line and
discarded the one real point in between.  The radius is therefore one
interval, and the points that actually break are repaired individually:

    pass 1   everything strictly inside the first mesh cell
    pass 2   any point still carrying a negative frequency, scaled in from
             the next shell out along its own direction

which repairs thallium without flattening niobium.  MEASURED, not assumed:
with ten ASR passes the second pass fires for NO element - all 24 come back
with zero negative frequencies from the first pass alone.  So the thallium
negative that motivated the two-interval radius was an ASR artefact, and the
radius was paying for it in every other element.  Both passes are recorded
in `sub` so the page can draw them dotted; the difference between "nothing was
measured here" and "the interpolation failed here" is not one the reader needs,
but the difference between measured and not is.
"""
import glob
import json
import os
import subprocess
import sys

#  This runs next to the fix-phonon binaries, wherever those are, so nothing
#  here is an absolute path.  Set the three variables to match the layout on
#  the machine that holds the runs:
#
#    FTC_DIR    holds ftq_index.json, the ftq_<El>.txt q-lists, and receives
#               ftcurve.json                        (default: this directory)
#    FTC_RUNS   the parent of the fix-phonon run directories, each named
#               fp_<El>_n21, fp_<El>_n10 or fp_<El> (default: one level up)
#    PHANA      LAMMPS' tools/phonon executable     (default: found on PATH)
M = os.environ.get("FTC_DIR", ".")
RUNS = os.environ.get("FTC_RUNS", "..")
PHANA = os.environ.get("PHANA", "phana")


def binary(el):
    for d in ("fp_%s_n21" % el, "fp_%s_n10" % el, "fp_%s" % el,
              "fixphonon_Cs" if el == "Cs" else None,
              "fixphonon" if el == "Ba" else None):
        if not d:
            continue
        p = os.path.join(RUNS, d)
        b = sorted(glob.glob(p + "/%sPhonon.bin.*" % el),
                   key=lambda x: int(x.rsplit(".", 1)[1]))
        if b:
            return b[-1], int(open(p + "/map.in").readline().split()[0])
    return None, None


def run(b, qs):
    #  Ten ASR passes, not one.  One is enough for a ONE-ATOM cell and not for
    #  two: the acoustic sum rule couples the sublattices, and hafnium's
    #  acoustic branches at Gamma come out at -1.223 THz after one pass and
    #  -0.072 after five.  It does not move the scored numbers - Hf 17.7
    #  either way, no scored point sits at Gamma - but it is what a drawn
    #  curve runs into.
    if not qs:
        return []
    inp = "10\n1\n4\n" + "".join("%.10f %.10f %.10f\n" % tuple(q) for q in qs) \
        + "q\n0\n"
    r = subprocess.run([PHANA, b], input=inp, capture_output=True, text=True)
    out, take = [], False
    for ln in r.stdout.splitlines():
        if "vibrational frequencies at this q-point" in ln:
            take = True
            continue
        if take:
            try:
                out.append(sorted(float(x) for x in ln.split()))
            except ValueError:
                pass
            take = False
    return out


def lift(b, qs, idx, radius, hcp, f):
    """Replace f[i] for i in idx by the value carried in from `radius`.

    Acoustic branches follow omega = v|q| and are scaled by the ratio; optic
    branches do not go to zero at Gamma and must not be scaled, so they are
    taken from the shell as they stand.
    """
    if not idx:
        return
    scaled, ratio = [], []
    for i in idx:
        s = radius / max(abs(v) for v in qs[i])
        ratio.append(s)
        scaled.append([v * s for v in qs[i]])
    g = run(b, scaled)
    if len(g) != len(idx):
        print("   shell %d/%d blocks, this pass skipped" % (len(g), len(idx)))
        return
    for j, i in enumerate(idx):
        up = [x / ratio[j] for x in g[j]]
        f[i] = (up[:3] + g[j][3:]) if hcp else up


def main():
    index = json.load(open(os.path.join(M, "ftq_index.json")))
    out = {}
    for el in sorted(index):
        b, N = binary(el)
        if not b:
            print("%-3s binary MISSING" % el)
            continue
        qs = [list(map(float, l.split())) for l in open(os.path.join(M, "ftq_%s.txt" % el))]
        f = run(b, qs)
        if len(f) != len(qs):
            print("%-3s %d/%d blocks, skipped" % (el, len(f), len(qs)))
            continue
        hcp = index[el]["struct"] == "hcp"
        step = 1.0 / N
        rad = [max(abs(v) for v in q) for q in qs]

        first = [i for i in range(len(qs)) if 1e-9 < rad[i] < step]
        lift(b, qs, first, step, hcp, f)

        #  whatever still comes back negative, one shell further out
        broke = [i for i in range(len(qs))
                 if rad[i] > 1e-9 and i not in first and min(f[i]) < 0]
        lift(b, qs, broke, 2.0 * step, hcp, f)

        #  At Gamma the three acoustic frequencies are ZERO, exactly, by the
        #  same sum rule phana is iterating towards.  Ten passes get hafnium
        #  to -0.072 THz and not to zero, so the value is set rather than
        #  approached: this is a constraint, not a measurement, and a curve
        #  that dips below the axis at Gamma would be reporting an
        #  instability the crystal does not have.
        for i, q in enumerate(qs):
            if rad[i] < 1e-9:
                f[i] = [0.0, 0.0, 0.0] + f[i][3:]

        neg = sum(1 for r in f for v in r if v < 0)
        sub = sorted(set(first) | set(broke))
        out[el] = {"f": [[round(v, 3) for v in r] for r in f],
                   "sub": sub, "mesh": N, "neg": neg}
        print("%-3s N=%-3d %4d q, first cell %2d, repaired %2d, negative %d"
              % (el, N, len(f), len(first), len(broke), neg))
    p = os.path.join(M, "ftcurve.json")
    json.dump(out, open(p, "w"), separators=(",", ":"))
    print("\n%s  %d bytes" % (p, os.path.getsize(p)))


if __name__ == "__main__":
    sys.exit(main())
