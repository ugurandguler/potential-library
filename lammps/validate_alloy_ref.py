#!/usr/bin/env python3
"""
The alloy path against an independent implementation, on a genuinely mixed cell.

`validate_alloy.py` proved the multi-species code is self-consistent: declare
one element twice and nothing moves.  That catches structural errors and cannot
catch a wrong convention, because every entry in a twin file holds the same
numbers.  This is the other half - two implementations, one ordered binary,
different numbers on the unlike bonds.

`standalone/alloy.py` is the reference.  It reads the same file LAMMPS reads,
so the only thing under test is what the two codes do with identical
parameters, and it is written as a plain double loop so that when they disagree
the fast one is the suspect.

Three structures, chosen because they put different bonds in different places:

  L1_0   alternating (001) planes in fcc - like and unlike bonds in the same
         coordination shell, so the (i,j,j) lookup has to be right per bond
  B2     caesium-chloride, bcc - every first neighbour is unlike and every
         second neighbour is like, which separates the two cutoffs cleanly
  L1_2   A3B in fcc - the only one of the three where a centre atom sees a
         genuinely mixed triplet, B on one leg and A on the other, which is
         the case the centre rule was written for

The last one is the point.  L1_0 and B2 can both be got right by code that
never looks at the second leg.

    python validate_alloy_ref.py                # Cu-Ni
    python validate_alloy_ref.py Ni Al
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import alloy            # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"
SKIN = 2.0

IN = """units           metal
boundary        p p p
atom_style      atomic
read_data       {data}
pair_style      ugur
pair_coeff      * * {pot} {els}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check no
run             0
variable        epa equal pe/atoms
print           "UGUR_E ${{epa}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def make_cell(kind, a, c_over_a=1.0):
    """(lattice rows, fractional positions, species index per atom)"""
    if kind == "L1_0":
        lat = np.diag([a, a, a * c_over_a])
        frac = np.array([[0., 0., 0.], [.5, .5, 0.],
                         [.5, 0., .5], [0., .5, .5]])
        sp = [0, 0, 1, 1]                    # z = 0 is A, z = 1/2 is B
    elif kind == "B2":
        lat = np.diag([a, a, a])
        frac = np.array([[0., 0., 0.], [.5, .5, .5]])
        sp = [0, 1]
    elif kind == "L1_2":
        lat = np.diag([a, a, a])
        frac = np.array([[0., 0., 0.], [.5, .5, 0.],
                         [.5, 0., .5], [0., .5, .5]])
        sp = [1, 0, 0, 0]                    # B at the corner, A on the faces
    else:
        raise ValueError(kind)
    cry = object.__new__(L.Crystal)
    cry.struct, cry.a0, cry.c_over_a = kind, a, c_over_a
    cry.lat, cry.frac = lat, frac
    cry.mass = np.ones(len(frac))
    return cry, sp


def write_data(path, cry, sp, masses, rep):
    box = np.diag(cry.lat).copy()
    pos = cry.frac @ cry.lat
    out, ids = [], []
    for i in range(rep[0]):
        for j in range(rep[1]):
            for k in range(rep[2]):
                out.append(pos + np.array([i, j, k]) * box)
                ids.extend(sp)
    pos = np.vstack(out)
    box = box * np.array(rep, float)
    with open(path, "w") as f:
        f.write(f"# {cry.struct}\n\n{len(pos)} atoms\n"
                f"{len(masses)} atom types\n\n")
        for lo, hi, ax in ((0., box[0], "x"), (0., box[1], "y"),
                           (0., box[2], "z")):
            f.write(f"{lo:.17g} {hi:.17g} {ax}lo {ax}hi\n")
        f.write("\nMasses\n\n")
        for t, m in enumerate(masses, 1):
            f.write(f"{t} {m:.17g}\n")
        f.write("\nAtoms # atomic\n\n")
        for n, (t, p) in enumerate(zip(ids, pos), 1):
            f.write(f"{n} {t + 1} {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}\n")
    return len(pos)


def main():
    els = sys.argv[1:] or ["Cu", "Ni"]
    A, B = els[0], els[1]
    d = os.path.join(HERE, "refruns", f"{A}{B}")
    os.makedirs(d, exist_ok=True)
    pot_path = os.path.join(d, f"{A}{B}.ugur.alloy")
    with open(pot_path, "w") as f:
        f.write(subprocess.run(
            [sys.executable, os.path.join(HERE, "make_alloy_file.py"),
             "--set", "tap", A, B],
            capture_output=True, text=True).stdout)
    pot = alloy.AlloyPotential.from_file(pot_path, [A, B])

    a0 = 0.5 * (refdata.ELEMENTS[A]["a0"] + refdata.ELEMENTS[B]["a0"])
    masses = [refdata.MASSES[A], refdata.MASSES[B]]

    print(f"{A}-{B}, tapered parameters, the same file given to both codes\n")
    print(f"{'struct':7s}{'atoms':>6s}{'E alloy.py':>14s}{'E lammps':>14s}"
          f"{'diff':>12s}   status")
    print("-" * 58)
    bad = 0
    for kind in ("L1_0", "B2", "L1_2"):
        cry, sp = make_cell(kind, a0, 0.95 if kind == "L1_0" else 1.0)
        ref = alloy.energy(cry, pot, sp)
        rcut = max(pot.rcut2, pot.rcut3)
        box = np.diag(cry.lat)
        rep = tuple(max(3, int(np.ceil(2.0 * (rcut + SKIN) / b))) for b in box)
        dd = os.path.join(d, kind)
        os.makedirs(dd, exist_ok=True)
        nat = write_data(os.path.join(dd, "cell.data"), cry, sp, masses, rep)
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cp {wsl(pot_path)} {wsl(dd)}/"], capture_output=True)
        open(os.path.join(dd, "in.check"), "w").write(IN.format(
            data="cell.data", skin=SKIN, pot=os.path.basename(pot_path),
            els=f"{A} {B}"))
        out = subprocess.run(["wsl", "-e", "bash", "-lc",
                              f"cd {wsl(dd)} && {LMP} -in in.check 2>&1"],
                             capture_output=True, text=True).stdout
        m = re.search(r"UGUR_E\s+([-\d.eE+]+)", out)
        if not m:
            e = [l for l in out.splitlines() if "ERROR" in l]
            print(f"{kind:7s}  LAMMPS: {(e or ['no output'])[0][:60]}")
            bad += 1
            continue
        got = float(m.group(1))
        diff = got - ref
        ok = abs(diff) < 1e-9 * max(abs(ref), 1.0)
        bad += not ok
        print(f"{kind:7s}{nat:6d}{ref:14.8f}{got:14.8f}{diff:12.2e}"
              f"   {'ok' if ok else 'DIFFERS'}")
    print()
    if bad:
        raise SystemExit(f"{bad} structures did not match")
    print("the multi-species kernel matches the independent implementation,")
    print("karisik ucluler dahil (L1_2)")


if __name__ == "__main__":
    main()
