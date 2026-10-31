#!/usr/bin/env python3
"""
Force constants of PUBLISHED potentials, for the dispersion comparison.

The library's own phonons come from analytic second derivatives (latdyn).  A
published EAM or MEAM file has no such code here, so its force constants are
taken from LAMMPS by finite differences: one atom of a supercell is moved by
+h and -h along x, y and z, and the force on every atom gives one row of
Phi(0, j) = -dF_j / du_0.  The blocks are keyed exactly as latdyn keys its
own, (i, j, R) in the primitive cell of latdyn.Crystal, so curve_mae scores
them through the same branch-label reading as the library's records - the
comparison differs in the force constants and in nothing else.

Two things decide whether that is true, and both are checked rather than
assumed:

  * the supercell must hold every pair the force constants couple.  In an
    embedded-atom model atom 0 moves the density at every neighbour k, and
    k's embedding couples 0 to every neighbour of k, so Phi(0, j) reaches out
    to TWICE the cut-off.  The box edge is therefore at least 4 r_c + 1 A, so
    that minimum image never folds two of those j onto one;

  * --check runs the library's own switched copper record through LAMMPS
    (pair_style ugur) the same way, and its frequencies must reproduce
    latdyn's analytic ones (curve_mae_published.py --check compares them).

The cell is built at the volume curve_mae scores the library at - the
measurement's lattice constant where the record gives one - so a published
potential and ours are evaluated at the same geometry.

Runs in WSL, where the LAMMPS build lives:
    python3 phonon_published.py            # every job -> phonon_published.json
    python3 phonon_published.py Cu W       # only these elements
    python3 phonon_published.py --check    # our Cu|tap record through LAMMPS
"""
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STAND = os.path.join(HERE, "..", "standalone")
#  the NIST downloads of the hold-out comparison (one folder per NIST id, with
#  a MANIFEST.json naming each id's files and pair style); not redistributed,
#  so set HOLDOUT_BASELINES to a local copy, or those potentials are skipped
UG2 = os.environ.get("HOLDOUT_BASELINES", "")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/src/lmp_serial"))
POTS = os.path.expanduser("~/lammps/potentials")
#  the NIST downloads the baselines were taken from live in a second folder
SEARCH = [POTS, os.path.expanduser("~/nistpot")]
WORK = os.path.expanduser("~/phonon_published")
OUT = os.path.join(HERE, "phonon_published.json")
H = 0.005              # A, central-difference step
SKIN = 0.3

sys.path.insert(0, STAND)
import latdyn as L     # noqa: E402

LIB = json.load(open(os.path.join(STAND, "library.json")))


def family(style):
    return "MEAM" if "meam" in style else ("EAM" if style.startswith("eam") else "other")


def coeff(style, files, el, syms=None):
    """pair_coeff for each flavour (see elastic_T.coeff_line: a funcfl .eam
    takes no element name, and writing one makes LAMMPS stop silently).  A
    multi-element MEAM file lists every element of its library, in the
    library's order, and maps the one atom type onto the element wanted"""
    if style == "eam":
        return f"1 1 {files[0]}"
    if style == "meam":
        return f"* * {files[0]} {' '.join(syms or [el])} {files[1]} {el}"
    return f"* * {files[0]} {el}"


def find(name):
    for d in SEARCH:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)


#  withdrawn by NIST (18 Nov 2025: phi(r) wrong at small r); ipr4 replaces it
SUPERSEDED = {"mo.fs.eam.alloy"}
EXTRA = [("Mo", "1987--Ackland-G-J--Mo--LAMMPS--ipr4", "eam/fs",
          [os.path.join(UG2, "1987--Ackland-G-J--Mo--LAMMPS--ipr4",
                        "Ackland_Thetford_1987_Mo.eam.fs")], None)] if UG2 else []


#  the MEAM potentials of Aitken et al., Phys. Rev. B 103, 094116 (2021): the
#  LAMMPS files of that paper's supplementary material, which are not in the
#  NIST repository and are not redistributed here (put them in AITKEN_DIR).
#  Kept out of OUT and scored on their own (--aitken), so the published set
#  the other comparisons use stays as it was
AITKEN_DIR = os.path.join(HERE, "aitken2021_meam")
AITKEN = [(el, f"Aitken2021-MEAM-{el}", "meam",
            [os.path.join(AITKEN_DIR, f"{el}_library.meam"), os.path.join(AITKEN_DIR, f"{el}.meam")],
            None)
           for el in ("Ni", "Al", "Cu")]
AITKEN_OUT = os.path.join(HERE, "phonon_aitken.json")


def jobs():
    """(el, label, style, [absolute files]) - the published potentials the
    library is already compared with elsewhere: the baselines of the
    elastic, expansion and equation-of-state screens, and those of the
    hold-out study, without repeats"""
    out, seen = [], set()
    for f in sorted(glob.glob(os.path.join(STAND, "baselines*.json"))):
        for el, rows in json.load(open(f)).items():
            for fn, style in rows:
                if fn in SUPERSEDED:
                    continue
                files = [find(x) for x in fn.split("+")]
                key = (el, os.path.basename(files[-1]))
                if key in seen:
                    continue
                seen.add(key)
                out.append((el, fn, style, files, None))
    if not UG2 or not os.path.exists(os.path.join(UG2, "MANIFEST.json")):
        print("HOLDOUT_BASELINES not set: the hold-out comparison's NIST potentials "
              "are skipped", flush=True)
        man = {}
    else:
        man = json.load(open(os.path.join(UG2, "MANIFEST.json")))
    for pid, m in man.items():
        #  a multi-element file is also a potential for each of its elements;
        #  the binaries and the high-entropy MEAMs fit the pure metals too
        syms = eval(m["symbols"])
        style = m["pair_style"]
        names = [x[0] for x in m["files"]
                 if not x[0].lower().endswith((".txt", ".pdf", ".json", ".md", ".in"))]
        for el in syms:
            if style == "meam":
                lib = next(n for n in names if "library" in n.lower())
                par = next(n for n in names if n != lib)
                files = [os.path.join(UG2, pid, lib), os.path.join(UG2, pid, par)]
            else:
                files = [os.path.join(UG2, pid, names[0])]
            key = (el, os.path.basename(files[-1]))
            if key in seen:
                continue
            seen.add(key)
            out.append((el, pid, style, files, syms if len(syms) > 1 else None))
    out = [j for j in out if os.path.basename(j[3][-1]) not in SUPERSEDED]
    return out + EXTRA


def geometry(el):
    """(Crystal, a, c/a) exactly as curve_mae builds it for scoring"""
    v = LIB[el]
    ec = v["exp_curve"]
    a = ec.get("a_meas", v["a0"])
    coa = v.get("c_over_a")
    if "c_meas" in ec:
        coa = ec["c_meas"] / a
    return L.Crystal(v["struct"], a, coa), a, coa


def supercell(cry, edge):
    """orthogonal supercell of at least `edge` A on every side; returns
    (box lengths, positions, [(basis index, integer R)] per atom)"""
    st, a = cry.struct, cry.a0
    if st == "hcp":
        cell = np.array([a, math.sqrt(3) * a, a * cry.c_over_a])
    else:
        cell = np.array([a, a, a])
    n = np.ceil(edge / cell).astype(int)
    box = cell * n
    inv = np.linalg.inv(cry.lat)
    #  every lattice point R + basis k whose position falls inside the box
    span = int(np.ceil(np.abs(box @ np.abs(inv)).max())) + 2
    rng = np.arange(-span, span + 1)
    R = np.array(np.meshgrid(rng, rng, rng, indexing="ij")).reshape(3, -1).T
    pos, lab = [], []
    for k, fk in enumerate(cry.frac):
        p = (R + fk) @ cry.lat
        ok = np.all((p > -1e-6) & (p < box - 1e-6), axis=1)
        pos.append(p[ok])
        lab += [(k, tuple(r)) for r in R[ok]]
    pos = np.vstack(pos)
    want = {"fcc": 4, "bcc": 2, "hcp": 4}[st] * int(np.prod(n))
    assert len(pos) == want, (len(pos), want)
    return box, pos, lab


def write_data(path, box, pos):
    with open(path, "w") as fh:
        fh.write(f"phonon supercell\n\n{len(pos)} atoms\n1 atom types\n\n")
        for lo, hi, ax in zip((0, 0, 0), box, "xyz"):
            fh.write(f"0.0 {hi:.10f} {ax}lo {ax}hi\n")
        fh.write("\nMasses\n\n1 1.0\n\nAtoms # atomic\n\n")
        for i, p in enumerate(pos, 1):
            fh.write(f"{i} 1 {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")


def lammps(d, text):
    open(os.path.join(d, "in.ph"), "w").write(text)
    subprocess.run(f"cd {d} && {LMP} -in in.ph -log log.ph > out.txt 2>&1", shell=True)
    return open(os.path.join(d, "log.ph")).read()


def head(style, cf):
    return (f"units metal\natom_style atomic\nboundary p p p\nread_data cell.data\n"
            f"pair_style {style}\npair_coeff {cf}\nneighbor {SKIN} bin\n"
            f"thermo_style custom step pe\n")


def cutoff(d, style, cf, cry):
    box, pos, _ = supercell(cry, 14.0)
    write_data(os.path.join(d, "cell.data"), box, pos)
    log = lammps(d, head(style, cf) + "run 0 post no\n")
    m = re.search(r"master list distance cutoff = ([\d.]+)", log)
    if not m:
        raise RuntimeError("no cutoff in log: " + log[-600:])
    return float(m[1]) - SKIN


def read_forces(path, n):
    F = np.zeros((n, 3))
    lines = open(path).read().splitlines()
    i = lines.index("ITEM: ATOMS id fx fy fz") + 1
    for ln in lines[i:i + n]:
        a = ln.split()
        F[int(a[0]) - 1] = [float(x) for x in a[1:4]]
    return F


def meam_library_ints(src, dst):
    """copy a MEAM library file with z, ielement and ibar written as integers

    Current LAMMPS reads those three as integers and stops on "12." or
    "72.000000000000" (pair_meam.cpp: "Not a valid integer number"), which
    published files written by older versions or by other codes carry - the
    Hf entry of Huang et al. (2021) and all three of Aitken et al. (2021).
    Every entry is 19 values whatever its line breaks (Huang's Hf is one
    line): elt lat z ielement atwt / alpha b0 b1 b2 b3 alat esub asub /
    t0 t1 t2 t3 rozero ibar.  Only a value that is a whole number is
    rewritten, and the rewritten file is checked to hold the same numbers."""
    head, toks = [], []
    for line in open(src, encoding="utf-8", errors="replace"):
        body = line.split("#", 1)[0]
        if not body.strip():
            if not toks:
                head.append(line.rstrip("\n"))
            continue
        toks += body.split()
    if len(toks) % 19:
        raise ValueError(f"{src}: {len(toks)} values, not a multiple of 19")

    def whole(x):
        v = float(x)
        if v != int(v):
            raise ValueError(f"{src}: integer field {x} is not a whole number")
        return str(int(v))
    out = []
    for k in range(0, len(toks), 19):
        e = list(toks[k:k + 19])
        e[2], e[3], e[18] = whole(e[2]), whole(e[3]), whole(e[18])
        assert all(float(a) == float(b) for a, b in zip(e[2:], toks[k + 2:k + 19]))
        out += ["  ".join(e[0:5]), "  ".join(e[5:13]), "  ".join(e[13:19])]
    open(dst, "w").write("\n".join(head + out) + "\n")


def force_constants(el, label, style, files, d, syms=None):
    os.makedirs(d, exist_ok=True)
    local = []
    for k, f in enumerate(files):
        if style == "meam" and k == 0:
            meam_library_ints(f, os.path.join(d, os.path.basename(f)))
        else:
            shutil.copy(f, d)
        local.append(os.path.basename(f))
    cf = coeff(style, local, el, syms) if style != "ugur" else f"* * {local[0]} {el}"
    cry, a, coa = geometry(el)
    rc = cutoff(d, style, cf, cry)
    box, pos, lab = supercell(cry, 4 * rc + 1.0)
    write_data(os.path.join(d, "cell.data"), box, pos)
    n = len(pos)
    #  the atom of each basis site that sits at R = 0
    centre = {k: lab.index((k, (0, 0, 0))) for k in range(len(cry.frac))}
    text = head(style, cf)
    for k, i0 in centre.items():
        text += f"group d{k} id {i0 + 1}\n"
    text += "run 0 post no\nwrite_dump all custom f_ref.dump id fx fy fz modify sort id\n"
    for k in centre:
        for ax, vec in enumerate(("1 0 0", "0 1 0", "0 0 1")):
            for sgn, tag in ((1, "p"), (-1, "m")):
                mv = " ".join(str(sgn * H * float(x)) for x in vec.split())
                back = " ".join(str(-sgn * H * float(x)) for x in vec.split())
                text += (f"displace_atoms d{k} move {mv} units box\nrun 0 post no\n"
                         f"write_dump all custom f_{k}{ax}{tag}.dump id fx fy fz modify sort id\n"
                         f"displace_atoms d{k} move {back} units box\n")
    log = lammps(d, text)
    if "ERROR" in log or not os.path.exists(os.path.join(d, f"f_{len(centre) - 1}2m.dump")):
        raise RuntimeError("LAMMPS: " + log[-800:])
    inv = np.linalg.inv(cry.lat)
    Phi = {}
    for k0, i0 in centre.items():
        rows = []
        for ax in range(3):
            Fp = read_forces(os.path.join(d, f"f_{k0}{ax}p.dump"), n)
            Fm = read_forces(os.path.join(d, f"f_{k0}{ax}m.dump"), n)
            rows.append(-(Fp - Fm) / (2 * H))           # (n, 3): d2E/du0_ax du_j
        rows = np.stack(rows, 1)                           # (n, 3 ax, 3 b)
        for j in range(n):
            blk = rows[j]
            if np.abs(blk).max() < 1e-7:
                continue
            dvec = pos[j] - pos[i0]
            dvec -= box * np.round(dvec / box)             # minimum image
            r_abs = cry.frac[k0] @ cry.lat + dvec
            f = r_abs @ inv
            hit = None
            for k, fk in enumerate(cry.frac):
                Rf = f - fk
                if np.abs(Rf - np.round(Rf)).max() < 1e-5:
                    hit = (k, tuple(int(x) for x in np.round(Rf)))
            assert hit is not None, (el, label, j)
            key = (k0, hit[0], hit[1])
            Phi[key] = Phi.get(key, np.zeros((3, 3))) + blk
    #  the acoustic sum rule, enforced on the on-site block (the finite
    #  difference leaves it at the 1e-6 level; it is reported, not hidden)
    asr = 0.0
    for k0 in centre:
        s = sum((b for (i, j, R), b in Phi.items() if i == k0), np.zeros((3, 3)))
        asr = max(asr, float(np.abs(s).max()))
        Phi[(k0, k0, (0, 0, 0))] -= s
    return dict(el=el, label=label, style=style, family=family(style), a=a, c_over_a=coa,
                rc=rc, n_atoms=n, asr_before=asr,
                phi=[[i, j, list(R), b.tolist()] for (i, j, R), b in Phi.items()])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    os.makedirs(WORK, exist_ok=True)
    if "--check" in sys.argv:
        #  our own switched copper record, through the same finite differences
        sys.path.insert(0, HERE)
        rec = dict(LIB["Cu"]["tap"])
        rec.setdefault("lam2", 0.0)
        rec.setdefault("lam4", 0.0)
        d = os.path.join(WORK, "check_Cu_tap")
        os.makedirs(d, exist_ok=True)
        #  the potential file in the form elastic_T.POTFILE writes it
        src = os.path.join(WORK, "Cu.ugur")
        open(src, "w").write("# Cu|tap, written by phonon_published.py --check\nCu Cu Cu "
                             + " ".join(f"{rec[k]:.17g}" for k in
                                        ("m", "D", "alpha", "r0", "gamma", "C", "alpha3",
                                         "rcut2", "rcut3", "taper", "lam2", "lam4")) + "\n")
        r = force_constants("Cu", "Cu|tap via LAMMPS", "ugur", [src], d)
        json.dump(r, open(os.path.join(HERE, "phonon_published_check.json"), "w"))
        print("check written:", r["n_atoms"], "atoms, rc", r["rc"], "ASR", r["asr_before"])
        return
    if "--aitken" in sys.argv:
        out = {}
        for el, label, style, files, syms in AITKEN:
            key = f"{el}|{label}"
            out[key] = force_constants(el, label, style, files,
                                       os.path.join(WORK, re.sub(r"[^\w.+-]", "_", key)), syms)
            r = out[key]
            print(f"{key:40s} rc {r['rc']:.2f}  {r['n_atoms']:6d} atoms  "
                  f"ASR {r['asr_before']:.1e}", flush=True)
        json.dump(out, open(AITKEN_OUT, "w"))
        return
    done = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for k in [k for k in done if k.split("|", 1)[1].split("+")[-1] in SUPERSEDED
              or any(x in k for x in ("ipr3",))]:
        del done[k]
    for el, label, style, files, syms in jobs():
        if el not in LIB or "exp_curve" not in LIB[el]:
            continue                     # no measured dispersion to score against
        if args and el not in args:
            continue
        key = f"{el}|{label}"
        if key in done and "phi" in done[key] and not args:
            continue
        d = os.path.join(WORK, re.sub(r"[^\w.+-]", "_", key))
        try:
            done[key] = force_constants(el, label, style, files, d, syms)
            r = done[key]
            print(f"{key:60s} {r['family']:5s} rc {r['rc']:.2f}  {r['n_atoms']:6d} atoms  "
                  f"ASR {r['asr_before']:.1e}", flush=True)
        except Exception as ex:
            done[key] = dict(el=el, label=label, style=style, family=family(style),
                             error=str(ex)[-500:])
            print(f"{key:60s} FAILED {str(ex)[-160:]}", flush=True)
        json.dump(done, open(OUT, "w"))


if __name__ == "__main__":
    main()
