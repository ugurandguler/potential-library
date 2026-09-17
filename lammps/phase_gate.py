#!/usr/bin/env python3
"""
Is there actually a solid in this cell?  The gate every melting run must pass.

WHY THIS EXISTS.  On 2026-09-05 a two-phase melting run reported a tantalum
melting point of 3438 K.  The cell contained **no solid at any point** - it had
melted in the first equilibration stage, before the two-phase construction, and
nothing downstream noticed.  The number looked entirely reasonable and was
reported as a held-out win before it was checked.

Nothing detected it because the run's own diagnostic, `compute SolidTemp`, is
defined on a group fixed by **initial position**, not by structure.  Two halves
of a uniform liquid report two temperatures that differ only by noise - which
is exactly what the log showed, and exactly what a real coexistence also shows.
**A temperature cannot tell you a phase.**  A structure factor can.

    S(G) = |sum_j exp(i G.r_j)| / N

with G a reciprocal-lattice vector of the crystal.  S = 1 for a perfect
lattice, ~1/sqrt(N) for a liquid.  Cheap, and it answers the actual question.

TWO THINGS THAT MAKE IT WRONG IF SKIPPED.

**G must be commensurate with the periodic box.**  Take G = 2*pi/a with a
guessed lattice constant and a box whose length is not an exact multiple of it,
and the sum picks up a spurious phase that drives S toward zero for a PERFECT
crystal.  Here G is built from `round(L/a)` whole wavelengths per box edge, so
it is exact by construction.

**Always print the control.**  A perfect crystal must come back at S = 1.000.
If it does not, the metric is broken and every other number on the page is
meaningless.  That control is what turned "both halves look liquid" from a
suspicion into a result.

    python phase_gate.py C2_both.dump --halves --data perfect.data --a 3.3052
    python phase_gate.py *.dump --halves --data perfect.data --a 3.3052

`a` is the lattice constant of the crystal the cell was built from; the run
writes `perfect.data` before it heats anything, which is the control.  See
in.coexistence_bracket for how the verdicts below are read as a melting point.
"""
import glob
import io
import os
import sys

import numpy as np

SOLID_MIN = 0.30      # over this, crystalline; a liquid of 1000 atoms is ~0.03
LIQUID_MAX = 0.15


def read_dump(path):
    """positions and box lengths from a LAMMPS text dump"""
    L = io.open(path, encoding="utf-8", errors="replace").read().splitlines()
    n = int(L[3])
    box = [list(map(float, L[i].split()[:2])) for i in (5, 6, 7)]
    head = L[8].split()[2:]
    ix, iy, iz = head.index("x"), head.index("y"), head.index("z")
    pos = np.array([[float(r.split()[k]) for k in (ix, iy, iz)]
                    for r in L[9:9 + n]])
    lo = np.array([b[0] for b in box])
    return pos - lo, np.array([b[1] - b[0] for b in box])


def read_data(path):
    """the same, from a LAMMPS data file - used for the control"""
    L = io.open(path, encoding="utf-8", errors="replace").read().splitlines()
    box = [list(map(float, l.split()[:2])) for l in L
           if any(t in l for t in ("xlo", "ylo", "zlo"))]
    i = [k for k, l in enumerate(L) if l.strip().startswith("Atoms")][0]
    pos = []
    for l in L[i + 1:]:
        f = l.split()
        if len(f) >= 5:
            try:
                pos.append([float(f[2]), float(f[3]), float(f[4])])
            except ValueError:
                pass
    lo = np.array([b[0] for b in box])
    return np.array(pos) - lo, np.array([b[1] - b[0] for b in box])


def sofg(pos, boxlen, a):
    """S(G) with G forced commensurate with the box - see the docstring"""
    n = np.maximum(np.round(boxlen / a), 1).astype(int)
    best = 0.0
    for axis in ((1, 1, 0), (1, 0, 1), (0, 1, 1), (2, 0, 0)):
        G = 2 * np.pi * np.array([axis[k] * n[k] / boxlen[k] for k in range(3)])
        if not G.any():
            continue
        best = max(best, abs(np.exp(1j * (pos @ G)).sum()) / len(pos))
    return best


def verdict(s):
    return ("KATI" if s > SOLID_MIN else
            "SIVI" if s < LIQUID_MAX else "ARADA")


def main():
    args = sys.argv[1:]
    a = 3.4
    if "--a" in args:
        i = args.index("--a"); a = float(args[i + 1]); del args[i:i + 2]
    halves = "--halves" in args
    if halves:
        args.remove("--halves")
    control = None
    if "--data" in args:
        i = args.index("--data"); control = args[i + 1]; del args[i:i + 2]
    files = [f for p in args for f in sorted(glob.glob(p))]
    if not files:
        print(__doc__)
        return 1

    #  THE CONTROL FIRST.  Without it the numbers below mean nothing.
    if control and os.path.exists(control):
        pos, L = read_data(control)
        s = max(sofg(pos, L, x) for x in np.arange(a * 0.9, a * 1.15, 0.01))
        print("KONTROL  %-22s S = %.3f   %s"
              % (os.path.basename(control), s,
                 "metrik saglam" if s > 0.9 else
                 "METRIK BOZUK - asagidaki hicbir sayiya guvenme"))
        print()
    else:
        print("KONTROL YOK.  --data ile kusursuz bir baslangic hucresi ver;")
        print("aksi halde dusuk S in sivi mi yoksa hatali G mi oldugunu"
              " bilemezsin.\n")

    hdr = ("%-22s %9s %9s   %s" % ("dosya", "alt", "ust", "hukum")) if halves \
        else ("%-22s %9s   %s" % ("dosya", "S", "hukum"))
    print(hdr)
    print("-" * len(hdr))
    for f in files:
        pos, L = read_dump(f)
        grid = np.arange(a * 0.9, a * 1.15, 0.01)
        if halves:
            zc = L[2] / 2
            lo_, hi_ = pos[pos[:, 2] < zc], pos[pos[:, 2] >= zc]
            s1 = max(sofg(lo_, L, x) for x in grid)
            s2 = max(sofg(hi_, L, x) for x in grid)
            v = ("IKI FAZ" if abs(s1 - s2) > 0.2 else
                 "HEPSI " + verdict(max(s1, s2)))
            print("%-22s %9.3f %9.3f   %s" % (os.path.basename(f), s1, s2, v))
        else:
            s = max(sofg(pos, L, x) for x in grid)
            print("%-22s %9.3f   %s" % (os.path.basename(f), s, verdict(s)))
    print("\nS ~ 1 kristal, S ~ 1/sqrt(N) sivi.  Birlikte-varolma kosusunda iki")
    print("yari FARKLI olmali; ikisi de dusukse ortada erime noktasi yoktur.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
