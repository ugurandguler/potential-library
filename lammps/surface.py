#!/usr/bin/env python3
"""
Surface energies, per facet, against density functional theory and experiment.

This is the test the vacancy result asks for.  Removing one atom costs 0.7-0.8
of the cohesive energy in this form where a real metal pays 0.2-0.35, because
the energy of a bond here does not depend on how many other bonds an atom has:
the neighbours of a hole gain nothing back.  A surface is the same physics
spread over a plane, so the prediction written down before running this was
that the surface energies would come out high by a similar factor, two to
three.  Whether they do decides whether six pathologies and a vacancy reduce to
one missing ingredient or to several.

    gamma = (E_slab - N E_bulk) / (2 A)

with the bulk relaxed first so that E_bulk belongs to the potential's own
lattice constant and not to the fitted one - otherwise a residual pressure in
the slab is counted as surface energy.  The slab is periodic in the plane and
carries vacuum along the normal, the atoms relax and the cell does not, and the
factor of two is the two faces.

Geometry.  Cubic facets are cut by orienting LAMMPS's own lattice, which is
what that command exists for.  The hexagonal ones need no reorientation at all:
the orthogonal hcp cell is a by a*sqrt(3) by c with x along [2-1-10] and y along
[01-10], so cutting normal to z, y and x gives (0001), (10-10) and (11-20) -
exactly the three facets the reference database resolves.

    python surface.py                    # every element with a reference
    python surface.py Pd Cu --sets tap,tap_ug
"""
import io
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))

#  windows+wsl here, a plain shell on the cluster
LOCAL = not os.path.exists("/arf")
HOME = (subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                       capture_output=True, text=True).stdout.strip()
        if LOCAL else os.path.expanduser("~"))
os.environ.setdefault("LMP", f"{HOME}/lammps/src/lmp_serial")
os.environ.setdefault("BASEFILE",
                      os.path.join(ROOT, "standalone", "baselines.json"))

import numpy as np         # noqa: E402
import latdyn as L         # noqa: E402
import cellfile           # noqa: E402
import refdata            # noqa: E402
import elastic_T as E     # noqa: E402

for _cand in (os.path.join(HERE, "potentials"),
              os.path.expanduser("~/lammps/potentials")):
    if os.path.isdir(_cand):
        E.POTDIR = _cand
        break

EV_A2_TO_J_M2 = 16.021766208
VACUUM = 20.0        # A, comfortably past the longest cutoff in the library
SKIN = 2.0

#  orient vectors for the cubic facets: z is the surface normal, x and y span
#  the plane, all three mutually orthogonal as LAMMPS requires
CUBIC = {
    "100": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    "110": ((-1, 1, 0), (0, 0, 1), (1, 1, 0)),
    "111": ((1, -1, 0), (1, 1, -2), (1, 1, 1)),
}
#  which facets to report, in the order the reference resolves them
WANT = {"fcc": ("111", "100", "110"),
        "bcc": ("110", "100", "111"),
        "hcp": ("0001", "10-10", "11-20")}
#  the hcp facets are the three faces of the orthogonal cell
HCP_AXIS = {"0001": "z", "10-10": "y", "11-20": "x"}

BULK = """units           metal
boundary        p p p
atom_style      atomic
{build}
pair_style      {style}
{coeff}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
fix             1 all box/relax iso 0.0 vmax 0.001
minimize        1e-14 1e-12 20000 200000
variable        e equal pe/atoms
variable        n equal atoms
variable        lxx equal lx
print           "EBULK ${{e}}"
print           "NAT ${{n}}"
print           "LX ${{lxx}}"
"""

SLAB = """units           metal
boundary        p p p
atom_style      atomic
{build}
pair_style      {style}
{coeff}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes

#  vacuum along the normal, added after the atoms exist so the slab keeps the
#  bulk spacing it was built with; boundary stays periodic because the gap is
#  wider than the cutoff and the images no longer see each other
change_box      all z delta 0.0 {vac} units box

variable        a0 equal lx*ly
variable        n0 equal atoms
run             0
variable        e0 equal pe
print           "AREA ${{a0}}"
print           "NAT ${{n0}}"
print           "EUNREL ${{e0}}"

min_style       cg
min_modify      dmax 0.05 line quadratic
minimize        1e-14 1e-12 20000 200000
variable        e1 equal pe
print           "ESLAB ${{e1}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def run(d, text, name):
    open(os.path.join(d, name), "w").write(text)
    p = os.path.join(d, "log.lammps")
    #  never read a log the current run did not write: the nudge filter was
    #  caught doing exactly that, and a stale log is indistinguishable from a
    #  fresh one once it is on disk
    if os.path.exists(p):
        os.remove(p)
    #  No pipes.  `capture_output=True` deadlocks this function when it is
    #  called from many threads at once: each call creates a pipe for the
    #  child's output, and a pipe's write end is inherited by whatever other
    #  children happen to be forked while it is open, so a thread waits for an
    #  end-of-file that a sibling's child is holding shut.  With sub-second
    #  runs and thirty-eight workers it fires almost immediately - the
    #  stacking-fault job on the cluster finished sixty-two runs in fifteen
    #  seconds and then sat for twenty-three minutes with the interpreter
    #  spinning at 121 % and a queue of unreaped children.  The thermal
    #  expansion job on the same code did not hit it, because its runs take
    #  half an hour and the window is never open.
    #
    #  Nothing was being read from those pipes in any case: LAMMPS writes to
    #  out.txt inside the shell and the answers are read back from log.lammps.
    quiet = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if LOCAL:
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {E.LMP} -in {name} > out.txt 2>&1"],
                       **quiet)
    else:
        subprocess.run(f"cd {d} && {E.LMP} -in {name} > out.txt 2>&1",
                       shell=True, **quiet)
    return io.open(p, errors="ignore").read() if os.path.exists(p) else ""


def grab(log, key):
    m = re.search(rf"^{key}\s+([-\d.eE+]+)", log, re.M)
    return float(m.group(1)) if m else None


#  Cubic lattice vectors in units of a: what a "lattice point" is.
PRIM = {"fcc": [(0.5, 0.5, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5)],
        "bcc": [(0.5, 0.5, 0.5), (-0.5, 0.5, 0.5), (0.5, -0.5, 0.5)]}


def period(struct, h):
    """shortest lattice translation along the integer direction h, in units of a

    Not the length of h itself: along [111] of an fcc lattice the shortest
    translation is a(1,1,1), because a(1,1,1)/2 is not a lattice vector, while
    along [1-10] it is a(1,-1,0)/2.  Getting this wrong makes a cell that is not
    periodic, and the slab then has a fault running through it that is counted
    as surface energy.
    """
    h = np.array(h, float)
    u = h / np.linalg.norm(h)
    P = np.array(PRIM[struct])
    best = None
    R = 4
    for n1 in range(-R, R + 1):
        for n2 in range(-R, R + 1):
            for n3 in range(-R, R + 1):
                v = n1 * P[0] + n2 * P[1] + n3 * P[2]
                if not v.any():
                    continue
                #  parallel to u, and pointing along it
                if np.linalg.norm(np.cross(v, u)) > 1e-9:
                    continue
                d = float(np.dot(v, u))
                if d > 1e-9 and (best is None or d < best):
                    best = d
    return best


def oriented_cell(struct, a0, facet, reps):
    """atoms of an orthorhombic supercell with facet's normal along z

    Built here rather than with LAMMPS's `lattice ... orient`, which hung on
    the fcc (111) orientation without producing an error - the log stopped
    after the banner.  Doing it explicitly also makes the cell auditable: the
    atom count has to equal the volume divided by the atomic volume, and that
    is checked below rather than assumed.
    """
    x, y, z = CUBIC[facet]
    ax = np.array([period(struct, x) * np.array(x) / np.linalg.norm(x),
                   period(struct, y) * np.array(y) / np.linalg.norm(y),
                   period(struct, z) * np.array(z) / np.linalg.norm(z)])
    L = np.linalg.norm(ax, axis=1) * a0                     # cell lengths
    U = (ax.T / np.linalg.norm(ax, axis=1)).T               # unit frame
    P = np.array(PRIM[struct]) * a0
    #  Generate lattice points over a range certain to cover the box, then keep
    #  the ones that land inside it.  The range has to follow the REPEATED box,
    #  not one cell: taking it from L alone produced 489 atoms where 1800 were
    #  needed, which the audit in cubic_build caught.  A quarter-filled slab
    #  would have looked like an enormous surface energy.
    n = int(np.ceil(max(L * np.array(reps, float)) / (a0 * 0.5))) + 3
    g = np.arange(-n, n + 1)
    N1, N2, N3 = np.meshgrid(g, g, g, indexing="ij")
    pts = (N1.ravel()[:, None] * P[0] + N2.ravel()[:, None] * P[1]
           + N3.ravel()[:, None] * P[2])
    frac = pts @ U.T                                        # into the frame
    tol = 1e-7
    keep = np.ones(len(frac), bool)
    for k in range(3):
        keep &= (frac[:, k] > -tol) & (frac[:, k] < L[k] * reps[k] - tol)
    return frac[keep], L * np.array(reps, float)


def write_cell(path, coords, box, mass):
    out = ["# surface.py", "", f"{len(coords)} atoms", "1 atom types", ""]
    for hi, t in zip(box, ("x", "y", "z")):
        out.append(f"0.0 {hi:.12f} {t}lo {t}hi")
    out += ["", "Masses", "", f"1 {mass:.6f}", "", "Atoms # atomic", ""]
    for i, c in enumerate(coords, 1):
        out.append(f"{i} 1 {c[0]:.12f} {c[1]:.12f} {c[2]:.12f}")
    open(path, "w").write("\n".join(out) + "\n")


def cubic_build(el, d, a0, facet, reps):
    struct = refdata.ELEMENTS[el]["struct"]
    coords, box = oriented_cell(struct, a0, facet, reps)
    #  the audit: atoms must equal volume over atomic volume, or the cell has
    #  either lost a layer or double-counted one
    vat = a0 ** 3 / (4 if struct == "fcc" else 2)
    want = box.prod() / vat
    if abs(len(coords) - want) > 0.5:
        raise RuntimeError(f"{facet}: {len(coords)} atom, beklenen {want:.1f}")
    write_cell(os.path.join(d, "cell.data"), coords, box, refdata.MASSES[el])
    return "read_data       cell.data\n"


#  the four-atom orthogonal hcp cell, in fractions of (a, a*sqrt3, c)
HCP_BASIS = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0),
             (0.5, 1.0 / 6.0, 0.5), (0.0, 2.0 / 3.0, 0.5)]


def hcp_cell(a0, coa, reps, axis):
    """orthogonal hcp supercell with the chosen face turned to lie normal to z

    The cell is a by a*sqrt(3) by c with x along [2-1-10] and y along [01-10],
    so its three faces are already (11-20), (10-10) and (0001) - the same three
    the reference database resolves.  Since the cell is orthogonal, pointing a
    different one along z is a permutation of the coordinates and of the box,
    not a rotation, which is why no reorientation machinery is needed here.
    """
    c = a0 * float(coa)
    L0 = np.array([a0, a0 * np.sqrt(3.0), c])
    cells = []
    for i in range(reps[0]):
        for j in range(reps[1]):
            for k in range(reps[2]):
                for f in HCP_BASIS:
                    cells.append((np.array(f) + (i, j, k)) * L0)
    pts = np.array(cells)
    box = L0 * np.array(reps, float)
    perm = {"z": (0, 1, 2), "y": (2, 0, 1), "x": (1, 2, 0)}[axis]
    return pts[:, perm], box[list(perm)]


def hcp_build(el, d, a0, coa, facet, reps):
    axis = HCP_AXIS[facet]
    coords, box = hcp_cell(a0, coa, reps, axis)
    #  same audit as the cubic cells: atoms against volume over atomic volume
    vat = (np.sqrt(3.0) / 2.0) * a0 * a0 * (a0 * float(coa)) / 2.0
    want = box.prod() / vat
    if abs(len(coords) - want) > 0.5:
        raise RuntimeError(f"{facet}: {len(coords)} atom, beklenen {want:.1f}")
    write_cell(os.path.join(d, "cell.data"), coords, box, refdata.MASSES[el])
    return "read_data       cell.data\n"


def one(job):
    el, tag, facet = job
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "surfruns", f"{el}_{safe}_{facet}")
    os.makedirs(d, exist_ok=True)
    struct = refdata.ELEMENTS[el]["struct"]
    try:
        style, pot = E.setup(d, el, tag)      # writes the potential file
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    coeff = E.coeff_line(style, pot, el)
    rc = 8.6      # the longest cutoff in the library, with room to spare
    a0 = float(E.PACK[el]["a0"])

    if struct == "hcp":
        coa = float(E.PACK[el].get("c_over_a") or 1.633)
        #  repeats chosen in the cell's own frame, then permuted with it, so
        #  the slab is thick along whichever axis becomes the normal
        base = np.array([a0, a0 * np.sqrt(3.0), a0 * coa])
        need = 2.4 * (rc + SKIN)
        reps = [max(1, int(np.ceil(need / b))) for b in base]
        ax = HCP_AXIS[facet]
        k = {"x": 0, "y": 1, "z": 2}[ax]
        reps[k] = max(reps[k], int(np.ceil(2.6 * (rc + SKIN) / base[k])), 4)
        build_bulk = hcp_build(el, d, a0, coa, facet, tuple(reps))
        build_slab = build_bulk
        normal = "z"
    else:
        #  repeats of the oriented cell, not of the cubic one: the oriented
        #  cell's own edges differ per facet, so the counts are chosen from
        #  its measured lengths rather than from a0
        struct = refdata.ELEMENTS[el]["struct"]
        _, base = oriented_cell(struct, a0, facet, (1, 1, 1))
        need = 2.4 * (rc + SKIN)
        reps = [max(1, int(np.ceil(need / b))) for b in base]
        #  and thicker along the normal, so the middle of the slab is bulk
        reps[2] = max(reps[2], int(np.ceil(2.6 * (rc + SKIN) / base[2])), 4)
        build_bulk = cubic_build(el, d, a0, facet, tuple(reps))
        build_slab = build_bulk
        normal = "z"

    if normal != "z":
        return {"error": f"normal {normal} henuz desteklenmiyor"}

    lb = run(d, BULK.format(build=build_bulk, style=style, coeff=coeff,
                            skin=SKIN), "in.bulk")
    ebulk = grab(lb, "EBULK")
    if ebulk is None:
        err = [l for l in lb.splitlines() if "ERROR" in l]
        return {"error": err[0][:80] if err else "bulk kosmadi"}

    ls = run(d, SLAB.format(build=build_slab, style=style, coeff=coeff,
                            skin=SKIN, vac=VACUUM), "in.slab")
    area, nat = grab(ls, "AREA"), grab(ls, "NAT")
    eslab, eunrel = grab(ls, "ESLAB"), grab(ls, "EUNREL")
    if None in (area, nat, eslab):
        err = [l for l in ls.splitlines() if "ERROR" in l]
        return {"error": err[0][:80] if err else "dilim kosmadi"}
    g = (eslab - nat * ebulk) / (2.0 * area)
    g_un = (eunrel - nat * ebulk) / (2.0 * area)
    return {"gamma": g * EV_A2_TO_J_M2, "gamma_unrelaxed": g_un * EV_A2_TO_J_M2,
            "E_bulk": ebulk, "atoms": int(nat), "area": area}


def main():
    argv = sys.argv[1:]
    sets = ["tap"]
    if "--sets" in argv:
        i = argv.index("--sets")
        sets = argv[i + 1].split(",")
        del argv[i:i + 2]
    els = [a for a in argv if not a.startswith("--")]

    #  next to the script as well as in standalone/, because on the cluster the
    #  tree is flat: looking only in ../standalone found nothing, the element
    #  list came back empty and the job reported zero runs in two seconds.  The
    #  guard below turned that into a visible failure instead of a silent one.
    ref = {}
    for p in (os.path.join(ROOT, "standalone", "surface_ref.json"),
              os.path.join(HERE, "surface_ref.json")):
        if os.path.exists(p):
            ref = json.load(open(p))
            break
    #  every element with parameters, not only those with a reference: the
    #  ordering test needs no reference, and filtering on one quietly dropped
    #  cobalt, rhenium and ytterbium from the whole sweep
    els = els or sorted(E.PACK)

    jobs = [(el, tag, f) for el in els for tag in sets
            if tag in E.PACK.get(el, {})
            for f in WANT[refdata.ELEMENTS[el]["struct"]]]
    #  the published potentials through the identical code.  Without them a
    #  ratio of three and a half says nothing about whose fault it is: the
    #  slab thickness, the vacuum, the bulk reference and the relaxation are
    #  all choices, and a baseline that comes out at its own published value
    #  is what turns them from assumptions into a measurement.
    if "--nobase" not in sys.argv:
        for el in els:
            for fn, _ in E.BASE.get(el, []):
                if all(os.path.exists(os.path.join(E.POTDIR, g))
                       for g in fn.split("+")):
                    jobs += [(el, "base|" + fn, f)
                             for f in WANT[refdata.ELEMENTS[el]["struct"]]]
    if not jobs:
        print("0 runs - element or set name not recognised")
        return
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"{len(jobs)} runs, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(one, jobs))

    out = {}
    for (el, tag, f), r in zip(jobs, res):
        out.setdefault(f"{el}|{tag}", {})[f] = r

    print()
    print(f"{'el':4s}{'set':8s}{'facet':8s}{'ours':>9s}{'DFT':>8s}{'ratio':>7s}"
          f"{'expt':>8s}{'ratio':>7s}   atoms")
    print("-" * 66)
    ratios = []
    for key, fs in out.items():
        el, tag = key.split("|", 1)
        R = ref.get(el, {})
        for f, r in fs.items():
            if "error" in r:
                print(f"{el:4s}{tag:8s}{f:8s}   {r['error'][:44]}")
                continue
            dft = (R.get("facets") or {}).get(f)
            exp = R.get("tyson")
            rd = (r["gamma"] / dft) if dft else None
            re_ = (r["gamma"] / exp) if exp else None
            if rd:
                ratios.append(rd)
            print(f"{el:4s}{tag:8s}{f:8s}{r['gamma']:9.3f}"
                  f"{(dft or 0):8.3f}{(rd or 0):7.2f}"
                  f"{(exp or 0):8.2f}{(re_ or 0):7.2f}   {r['atoms']}")
    if ratios:
        a = np.array(ratios)
        print(f"\nratio to DFT: median {np.median(a):.2f}, "
              f"range {a.min():.2f}-{a.max():.2f}, {len(a)} facets")
        print("Prediction recorded beforehand: 2-3 times too high.")
    fp = os.path.join(HERE, "surface.json")
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
    print(f"{len(out)} records written; file now holds {len(old)}  -> {fp}")


if __name__ == "__main__":
    main()
