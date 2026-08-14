#!/usr/bin/env python3
"""
The EAM comparison on hcp, which the cubic recipe could not do.

`eam_compare.py` stops at fcc and bcc for a specific reason: LAMMPS's ELASTIC
example relaxes the box with `fix box/relax aniso`, which lets the three box
lengths move independently.  On the orthogonal hcp cell those lengths are a,
a*sqrt(3) and c, so `aniso` is free to break the ratio between the first two -
that is, to leave the basal plane no longer hexagonal - and every constant
measured afterwards belongs to a crystal that is not the one being studied.
`couple xy` is the fix: x and y scale together, so a relaxes as one number
while c relaxes separately, which is exactly the two degrees of freedom hcp
has.

The cells come from `cellfile.py`, for the same reason as everywhere else here:
`lattice hcp` fixes c/a at the ideal sqrt(8/3) and all twelve of ours sit below
it.

What is compared is the five independent constants.  The recipe returns the
whole 6x6, so C66 is left over and hexagonal symmetry requires it to equal
(C11-C12)/2 - a free check on each side separately, and one the EAM sets have
to pass too.

Only magnesium and zirconium have an EAM set shipped with LAMMPS; titanium's
files are MEAM and RANN, which are different pair styles and a different claim.
Two elements is thin and is reported as two elements.

    python eam_hcp.py            # Mg and Zr
    python eam_hcp.py Mg
"""
import json
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
import cellfile         # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"
POT = f"{HOME}/lammps/potentials"

#  what LAMMPS actually ships for hcp elements in this library
SETS = {
    "Mg": [("Mg_mm.eam.fs", "eam/fs", ["Mg"])],
    "Zr": [("Zr_mm.eam.fs", "eam/fs", ["Zr"])],
}

INIT = """variable up equal 1.0e-5
variable atomjiggle equal 1.0e-5
units           metal
variable cfac equal 1.0e-4
variable cunits string GPa
variable etol equal 0.0
variable ftol equal 1.0e-12
variable maxiter equal 1000
variable maxeval equal 10000
variable dmax equal 1.0e-2
boundary        p p p
read_data       {data}
neigh_modify every 1 delay 0 check yes
"""

POTMOD = """pair_style      {style}
pair_coeff      {coeff}
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes
min_style       cg
min_modify      dmax ${{dmax}} line quadratic
thermo          1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
"""

POTFILE = """# {el}, written by eam_hcp.py from library.json
# Tersoff-style: one line per (centre, leg, leg) triple.  For one element that
# is a single line, and it is written this way so the single-species files and
# the alloy files are the same format - a regression here is a regression in
# what the alloy path will use.
# el1 el2 el3  m D alpha r0 gamma C alpha3 rcut2 rcut3 taper lam2 lam4
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} 0 0
"""

WANT = ("C11", "C33", "C12", "C13", "C44", "C66")
LXYZ = re.compile(r"^\s*\d+\s+\S+\s+\S+\s+\S+(?:\s+\S+){6}\s+"
                  r"([\d.]+)\s+([\d.]+)\s+([\d.]+)", re.M)


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def couple_xy(path):
    """rewrite the ELASTIC driver so relaxation keeps the basal plane hexagonal"""
    txt = open(path).read()
    out = txt.replace("fix 3 all box/relax  aniso 0.0",
                      "fix 3 all box/relax  x 0.0 y 0.0 z 0.0 couple xy")
    if out == txt:
        raise SystemExit("no box/relax line in in.elastic")
    open(path, "w").write(out)


def run(el, tag, style, coeff, relax, rcut):
    d = os.path.join(HERE, "hcpruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    for f in ("in.elastic", "displace.mod"):
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cp {POT}/../examples/ELASTIC/{f} {wsl(d)}/"],
                       capture_output=True)
    if relax:
        couple_xy(os.path.join(d, "in.elastic"))
    else:
        #  ours is pinned to the experimental cell by construction, so relaxing
        #  it would measure a different crystal from the one it was fitted to.
        #  Read before opening for write: `open(p, "w")` truncates immediately,
        #  so doing both in one expression hands back an empty file.
        p = os.path.join(d, "in.elastic")
        txt = open(p).read()
        #  and `unfix 3` has to go with it, or LAMMPS stops on a fix it was
        #  never given
        open(p, "w").write(txt
                           .replace("fix 3 all box/relax  aniso 0.0",
                                    "# box fixed: fitted geometry")
                           .replace("unfix 3", "# unfix 3"))

    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    box, _ = cellfile.orthogonal_cell(cry)
    rep = tuple(max(2, int(np.ceil(1.2 * (rcut + 1.0) / b))) for b in box)
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep, triclinic=True)
    open(os.path.join(d, "init.mod"), "w").write(
        INIT.format(data=f"{el}.data"))
    open(os.path.join(d, "potential.mod"), "w").write(
        POTMOD.format(style=style, coeff=coeff))
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {LMP} -in in.elastic 2>&1"],
                       capture_output=True, text=True)
    got = {}
    for k in WANT:
        m = re.search(rf"{k}all\s*=\s*([-\d.eE+]+)", r.stdout)
        if m:
            got[k] = float(m.group(1))
    if not got:
        err = [l for l in r.stdout.splitlines() if "ERROR" in l]
        return {"err": err[0][:80] if err else "output unreadable"}
    m = LXYZ.findall(r.stdout)
    if m:
        lx, ly, lz = (float(x) for x in m[-1])
        got["a"] = lx / rep[0]
        got["c_over_a"] = (lz / rep[2]) / (lx / rep[0])
    return got


def err(got, exp):
    ks = [k for k in ("C11", "C33", "C12", "C13", "C44") if k in got]
    return 100.0 * (sum(((got[k] - exp[k]) / exp[k]) ** 2
                        for k in ks) / len(ks)) ** 0.5


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or sorted(SETS)
    print("hcp esneklik sabitleri (GPa) - EAM kutusunu gevsetiyor (couple xy),")
    print("ours is fixed at the experimental cell, because that is where it was fitted.\n")
    for el in els:
        e = refdata.ELEMENTS[el]
        exp = e["Cij"]
        rows = []
        for fn, style, names in SETS.get(el, []):
            coeff = f"* * {POT}/{fn} " + " ".join(names)
            rows.append((fn, run(el, fn.split(".")[0], style, coeff,
                                 True, 6.0)))
        for tag, rec in (("ours, hard", lib[el]),
                         ("ours, tapered", lib[el].get("tap"))):
            if not rec:
                continue
            d = os.path.join(HERE, "hcpruns", f"{el}_{tag.split()[-1]}")
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(
                el=el, taper=(rec.get("taper") or -1.0),
                **{k: rec[k] for k in ("m", "D", "alpha", "r0", "gamma",
                                       "C", "alpha3", "rcut2", "rcut3")}))
            rows.append((tag, run(el, tag.split()[-1], "ugur",
                                  f"* * {el}.ugur {el}",
                                  False, max(rec["rcut2"], rec["rcut3"]))))

        print(f"=== {el} ===")
        print(f"{'':16s}" + "".join(f"{k:>9s}" for k in
                                    ("C11", "C33", "C12", "C13", "C44"))
              + f"{'RMS %':>9s}{'a0 err':>9s}{'C66 check':>13s}")
        print(f"{'experiment':16s}" + "".join(f"{exp[k]:9.1f}" for k in
                                         ("C11", "C33", "C12", "C13", "C44")))
        for name, got in rows:
            if "err" in got:
                print(f"{name:16s}  error: {got['err']}")
                continue
            da = (f"{100*(got['a']-e['a0'])/e['a0']:+8.2f}%"
                  if "a" in got else f"{'fixed':>9s}")
            #  hexagonal symmetry, imposed nowhere
            chk = (f"{got['C66'] - 0.5*(got['C11']-got['C12']):+12.3f}"
                   if "C66" in got else f"{'-':>13s}")
            print(f"{name:16s}"
                  + "".join(f"{got.get(k, float('nan')):9.1f}" for k in
                            ("C11", "C33", "C12", "C13", "C44"))
                  + f"{err(got, exp):9.2f}{da}{chk}")
        print()


if __name__ == "__main__":
    main()
