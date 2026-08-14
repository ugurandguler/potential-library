#!/usr/bin/env python3
"""
Prove ugurpot.h against standalone/latdyn.py, which is the reference.

The C file is meant to be the one place the physics lives, with the LAMMPS pair
style and any Python binding as thin wrappers.  That is only worth anything if
it demonstrably agrees with the implementation every published number here came
from, so this compares them directly - value and derivative, both terms, with
the cutoff switch on and off - over every element in the library.

Neighbour bookkeeping is deliberately not tested: it is LAMMPS's, it is already
trusted, and mixing it in would hide a kernel error behind a list error.  What
is tested is each interaction.

    python validate_kernel.py            # every element, taper on and off
    python validate_kernel.py Cu Fe
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import latdyn as L          # noqa: E402

LIB = os.path.join(ROOT, "standalone", "library.json")

#  Scaled by the magnitude of the term over the whole sample, not by the value
#  at each point.  Both sides evaluate the same closed form and differ only in
#  the order of the arithmetic, so where the switch has taken a quantity five
#  orders below its own scale the pointwise ratio measures cancellation rather
#  than disagreement: magnesium's worst leg pair sits 0.1 % inside rcut3 and
#  differs by 9.4e-19 eV/A on a term of size 3.5e-3, which read as 3.8e-11
#  pointwise and is nothing.
TOL = 1e-13                 # of the largest value the term reaches


def wsl(path):
    """windows path -> the path WSL sees"""
    p = os.path.abspath(path).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def build():
    cmd = (f"cd {wsl(HERE)} && gcc -O2 -std=c99 -o kernel_probe "
           f"kernel_probe.c -lm")
    r = subprocess.run(["wsl", "-e", "bash", "-lc", cmd],
                       capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("derleme basarisiz:\n" + r.stdout + r.stderr)


def probe(par, lines):
    head = " ".join(f"{x:.17g}" for x in par)
    inp = head + "\n" + "\n".join(lines) + "\n"
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(HERE)} && ./kernel_probe"],
                       input=inp, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("probe hatasi:\n" + r.stderr)
    return [[float(x) for x in ln.split()]
            for ln in r.stdout.strip().splitlines()]


def worst_scaled(pairs):
    """max |C - python| over the sample, divided by the sample's own scale"""
    if not pairs:
        return 0.0
    scale = max(max(abs(x), abs(y)) for x, y in pairs)
    if scale <= 0.0:
        return 0.0
    return max(abs(x - y) for x, y in pairs) / scale


def check(el, d, taper):
    pot = L.Potential.from_record(dict(d, taper=taper))
    rc2, rc3 = d["rcut2"], d["rcut3"]
    par = [d["m"], d["D"], d["alpha"], d["r0"], d["gamma"], d["C"],
           d["alpha3"], rc2, rc3, taper if taper else -1.0]

    #  sample right through the switch window as well as the ordinary range,
    #  because that is where an implementation can differ and still look right
    r2s = np.concatenate([np.linspace(0.6 * d["dnn"], rc2 * 0.999, 240),
                          np.linspace(0.80 * rc2, 0.999 * rc2, 160)])
    legs = [(a, b)
            for a in np.linspace(0.6 * d["dnn"], 0.999 * rc3, 26)
            for b in np.linspace(0.6 * d["dnn"], 0.999 * rc3, 26)]

    lines = [f"2 {r:.17g}" for r in r2s]
    lines += [f"3 {a:.17g} {b:.17g}" for a, b in legs]
    out = probe(par, lines)

    #  each quantity judged against its own scale, so a small derivative is not
    #  held to the tolerance of a large energy or the other way round
    v2, g2, v3, g3 = [], [], [], []
    for r, row in zip(r2s, out[:len(r2s)]):
        v2.append((row[0], float(pot.phi2(np.array([r]), 0)[0])))
        g2.append((row[1], float(pot.phi2(np.array([r]), 1)[0])))
    for (a, b), row in zip(legs, out[len(r2s):]):
        E, d1, d2, _, _, _ = pot.leg3(a, b)
        v3.append((row[0], float(E)))
        g3.append((row[1], float(d1)))
        g3.append((row[2], float(d2)))
    return max(worst_scaled(v2), worst_scaled(g2),
               worst_scaled(v3), worst_scaled(g3))


def main():
    build()
    lib = json.load(open(LIB))
    els = sys.argv[1:] or sorted(lib)
    print(f"{'el':4s}{'taper off':>14s}{'taper 0.85':>14s}")
    print("-" * 32)
    bad = []
    for el in els:
        d = lib[el]
        a = check(el, d, None)
        b = check(el, d, 0.85)
        print(f"{el:4s}{a:14.2e}{b:14.2e}")
        if a > TOL or b > TOL:
            bad.append(f"{el}: {max(a, b):.2e}")
    print()
    if bad:
        for x in bad:
            print("FAIL ", x)
        raise SystemExit(1)
    print(f"ugurpot.h agrees with latdyn.py to {TOL:.0e} relative")


if __name__ == "__main__":
    main()
