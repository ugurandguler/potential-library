#!/usr/bin/env python3
"""
Write a LAMMPS data file straight from latdyn's own Crystal.

`validate_pair.py` builds its cells with LAMMPS's `lattice` command, which is
fine for fcc and bcc and wrong for every hcp metal in the library: `lattice hcp`
fixes c/a at the ideal sqrt(8/3), and all twelve of ours sit below it, beryllium
by four per cent.  Feeding LAMMPS the ideal ratio and comparing against latdyn
at the real one would measure the difference between two crystals, not two
codes.

Rather than hand-typing a basis for the orthogonal hcp cell - which is the sort
of thing that looks right and puts one atom in the wrong place - the cell is
derived from the primitive vectors latdyn already uses.  For each structure
there is an integer matrix T with

    T @ lat_primitive = a diagonal box

so the supercell is orthogonal by construction and holds |det T| times the
primitive basis.  The atoms are then enumerated from the primitive lattice
itself and wrapped, so the only geometry anywhere is latdyn's.

    fcc   T = [[-1,1,1],[1,-1,1],[1,1,-1]]   det 4   conventional cube
    bcc   T = [[0,1,1],[1,0,1],[1,1,0]]      det 2   conventional cube
    hcp   T = [[1,0,0],[1,2,0],[0,0,1]]      det 2   a x a*sqrt3 x c

All three come out diagonal, so no triclinic box and no tilt factors are
needed - which matters, because the natural hexagonal cell has xy = -a/2 sitting
exactly on the skew limit LAMMPS enforces.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "standalone"))
import latdyn as L      # noqa: E402

T_SUPER = {
    "fcc": np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
    "bcc": np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]),
    "hcp": np.array([[1, 0, 0], [1, 2, 0], [0, 0, 1]]),
}
TOL = 1e-9


def orthogonal_cell(cry):
    """(box 3-vector, cartesian positions) for the diagonal supercell"""
    T = T_SUPER[cry.struct]
    sup = T @ cry.lat                       # rows are the supercell vectors
    off = sup - np.diag(np.diag(sup))
    if np.abs(off).max() > TOL * max(1.0, np.abs(sup).max()):
        raise ValueError(f"{cry.struct}: supercell is not diagonal\n{sup}")
    box = np.diag(sup).copy()
    if (box <= 0).any():
        raise ValueError(f"{cry.struct}: non-positive box {box}")

    #  Enumerate primitive translations far enough out to cover the box, then
    #  keep whatever lands inside it.  |det T| tells us how many to expect, so a
    #  miscount is caught here rather than showing up as a wrong energy.
    n = int(abs(round(np.linalg.det(T))))
    inv = np.linalg.inv(sup)
    rng = range(-3, 4)
    pos = []
    for i in rng:
        for j in rng:
            for k in rng:
                shift = np.array([i, j, k]) @ cry.lat
                for b in cry.pos:
                    f = (b + shift) @ inv
                    f = f - np.floor(f + TOL)       # wrap, seam atoms to 0
                    f[np.abs(f - 1.0) < TOL] = 0.0
                    pos.append(f)
    pos = np.array(pos)

    keep = []
    for f in pos:
        if not any(np.abs(f - g).max() < 1e-6 for g in keep):
            keep.append(f)
    want = n * len(cry.frac)
    if len(keep) != want:
        raise ValueError(f"{cry.struct}: {len(keep)} atoms, expected {want}")
    return box, np.array(keep) @ sup


def write_data(path, cry, mass, rep=(1, 1, 1), triclinic=False):
    """LAMMPS data file for cry, replicated rep times along each axis

    `triclinic` writes a zero `xy xz yz` line.  The box is still orthogonal -
    the tilts are zero - but LAMMPS decides once, at read time, whether a box
    can ever be sheared, and refuses `change_box ... xy` on one that was
    declared orthogonal.  Anything that applies a shear strain later needs
    this set from the start.
    """
    box, pos = orthogonal_cell(cry)
    nx, ny, nz = rep
    big = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                big.append(pos + np.array([i, j, k]) * box)
    pos = np.vstack(big)
    box = box * np.array(rep, float)

    with open(path, "w") as f:
        f.write(f"# {cry.struct} a0={cry.a0} c/a={cry.c_over_a} "
                f"written by cellfile.py from latdyn.Crystal\n\n")
        f.write(f"{len(pos)} atoms\n1 atom types\n\n")
        for lo, hi, ax in ((0.0, box[0], "x"), (0.0, box[1], "y"),
                           (0.0, box[2], "z")):
            f.write(f"{lo:.17g} {hi:.17g} {ax}lo {ax}hi\n")
        if triclinic:
            f.write("0.0 0.0 0.0 xy xz yz\n")
        f.write(f"\nMasses\n\n1 {mass:.17g}\n\nAtoms # atomic\n\n")
        for i, p in enumerate(pos, 1):
            f.write(f"{i} 1 {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}\n")
    return len(pos), box


def selftest():
    """the supercell must hold the same energy per atom as the primitive cell"""
    import refdata
    print(f"{'el':4s}{'struct':>6s}{'atoms':>6s}{'E primitive':>13s}"
          f"{'E supercell':>13s}{'diff':>11s}")
    print("-" * 53)
    import json
    lib = json.load(open(os.path.join(os.path.dirname(HERE),
                                      "standalone", "library.json")))
    bad = 0
    for el in sorted(lib):
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        pot = L.Potential.from_record(lib[el])
        E1 = L.energy(cry, pot)
        box, pos = orthogonal_cell(cry)
        #  same crystal, described by a different cell: build it as a Crystal
        #  and let latdyn compute the energy the ordinary way
        sup = object.__new__(L.Crystal)
        sup.struct, sup.a0, sup.c_over_a = e["struct"], e["a0"], e.get("c_over_a")
        sup.lat = np.diag(box)
        sup.frac = pos @ np.linalg.inv(np.diag(box))
        sup.mass = np.full(len(pos), refdata.MASSES[el])
        E2 = L.energy(sup, pot)
        d = E2 - E1
        if abs(d) > 1e-9 * max(abs(E1), 1.0):
            bad += 1
        print(f"{el:4s}{e['struct']:>6s}{len(pos):6d}{E1:13.6f}{E2:13.6f}"
              f"{d:11.2e}")
    print()
    print("mismatched:" if bad else "every cell gives the primitive cell's energy",
          bad if bad else "")
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if selftest() else 0)
