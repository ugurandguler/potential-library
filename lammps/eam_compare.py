#!/usr/bin/env python3
r"""
Elastic constants from EAM, for comparison with this potential.

The question a referee asks first is why anyone should use this form when
LAMMPS already ships EAM for most metals.  Answering it needs EAM's own numbers
computed the same way ours are, not quoted from papers that used different
targets, so they are computed here with LAMMPS's own ELASTIC recipe.

Three things make the comparison fair, and each costs something:

  Each potential at its own equilibrium.  `in.elastic` relaxes the box before
  straining it, so EAM is measured where EAM wants to be.  Ours is pinned to
  the experimental lattice constant by construction - a hard constraint, not a
  fitted target - so EAM's lattice-constant error is reported alongside and
  belongs in any honest table.

  Every set that ships, not one.  Copper alone has five.  Picking the weakest
  would flatter us, so all of them are run and the best is reported with the
  file that produced it; the rest are kept in the JSON.

  Both of ours.  The library ships hard-truncated parameters; the tapered
  refit is a second set.  Quoting only whichever is better per element would
  be the same cheat in the other direction, so both columns are printed.

Morse is deliberately absent.  A central pair potential gives C12 = C44
exactly on a cubic Bravais lattice, so the best it can do is a closed-form
bound rather than a measurement.  Running it would be theatre.

hcp is not covered yet: the recipe relaxes a cubic box and c/a needs its own
treatment.

    python eam_compare.py            # every element with a potential
    python eam_compare.py Cu Fe
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import refdata          # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
#  LAMMPS does not expand ~ in a potential path, so it is resolved here.
LMP = f"{HOME}/lammps/src/lmp_serial"
POT = f"{HOME}/lammps/potentials"

#  element -> [(file, pair_style, element names)]
#
#  `eam` is the single-element funcfl format and its pair_coeff is `i j file`
#  with no element names; `eam/alloy` and `eam/fs` are setfl and take
#  `* * file El...`.  Mixing the two is a hard error rather than a silent one,
#  which is the only good thing about it.
CHOICE = {
    "Ag": [("Ag_u3.eam", "eam", [])],
    "Al": [("Al_zhou.eam.alloy", "eam/alloy", ["Al"]),
           ("Al_jnp.eam", "eam", []),
           ("Al_mm.eam.fs", "eam/fs", ["Al"])],
    "Au": [("Au_u3.eam", "eam", [])],
    "Cu": [("Cu_mishin1.eam.alloy", "eam/alloy", ["Cu"]),
           ("Cu_zhou.eam.alloy", "eam/alloy", ["Cu"]),
           ("Cu_u3.eam", "eam", []),
           ("Cu_u6.eam", "eam", []),
           ("Cu_smf7.eam", "eam", [])],
    "Fe": [("Fe_mm.eam.fs", "eam/fs", ["Fe"])],
    "Ni": [("Ni_u3.eam", "eam", []), ("Ni_smf7.eam", "eam", [])],
    "Pd": [("Pd_u3.eam", "eam", [])],
    "Pt": [("Pt_u3.eam", "eam", [])],
    "Ta": [("CuTa.eam.alloy", "eam/alloy", ["Ta"])],
    "V":  [("VFe_mm.eam.fs", "eam/fs", ["V"])],
    "W":  [("W_zhou.eam.alloy", "eam/alloy", ["W"])],
}

LAT = {"fcc": "fcc", "bcc": "bcc"}

INIT = """variable up equal 1.0e-6
variable atomjiggle equal 1.0e-5
units           metal
variable cfac equal 1.0e-4
variable cunits string GPa
variable etol equal 0.0
variable ftol equal 1.0e-10
variable maxiter equal 200
variable maxeval equal 2000
variable dmax equal 1.0e-2
boundary        p p p
lattice         {lat} {a0}
region          box prism 0 1 0 1 0 1 0 0 0
create_box      1 box
create_atoms    1 box
mass 1 {mass}
neigh_modify every 1 delay 0 check yes
"""

POTMOD = """pair_style      {style}
pair_coeff      {coeff}
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes
min_style       cg
min_modify      dmax ${{dmax}} line quadratic
thermo          1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
"""

CIJ = (("C11", r"C11all\s*=\s*([-\d.eE+]+)"),
       ("C12", r"C12all\s*=\s*([-\d.eE+]+)"),
       ("C44", r"C44all\s*=\s*([-\d.eE+]+)"))
LX = r"^\s*\d+\s+\S+\s+\S+\s+\S+(?:\s+\S+){6}\s+([\d.]+)"


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def one(el, e, fn, style, names, d):
    coeff = (f"1 1 {POT}/{fn}" if style == "eam"
             else f"* * {POT}/{fn} " + " ".join(names))
    open(os.path.join(d, "init.mod"), "w").write(
        INIT.format(lat=LAT[e["struct"]], a0=e["a0"], mass=refdata.MASSES[el]))
    open(os.path.join(d, "potential.mod"), "w").write(
        POTMOD.format(style=style, coeff=coeff))
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {LMP} -in in.elastic 2>&1"],
                       capture_output=True, text=True)
    out = r.stdout
    got = {}
    for k, pat in CIJ:
        m = re.search(pat, out)
        if m:
            got[k] = float(m.group(1))
    if len(got) < 3:
        return None
    m = re.findall(LX, out, re.M)
    got["a0"] = float(m[-1]) if m else None
    got["file"] = fn
    c = e["Cij"]
    errs = [(got[k] - c[k]) / c[k] for k in ("C11", "C12", "C44")]
    got["rms"] = 100.0 * (sum(x * x for x in errs) / 3) ** 0.5
    return got


def run(el):
    """best of every set that ships for this element"""
    e = refdata.ELEMENTS[el]
    if e["struct"] not in LAT or el not in CHOICE:
        return None, []
    d = os.path.join(HERE, "eamruns", el)
    os.makedirs(d, exist_ok=True)
    for f in ("in.elastic", "displace.mod"):
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cp ~/lammps/examples/ELASTIC/{f} {wsl(d)}/"],
                       capture_output=True)
    tried = [g for g in (one(el, e, fn, st, nm, d)
                         for fn, st, nm in CHOICE[el]) if g]
    if not tried:
        return None, []
    return min(tried, key=lambda g: g["rms"]), tried


def main():
    els = sys.argv[1:] or sorted(CHOICE)
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    out = {}
    head = ("el", "best EAM", "n", "a0 %", "EAM", "OURS hard", "OURS taper")
    print("{:4s}{:>24s}{:>4s}{:>8s}{:>8s}{:>10s}{:>11s}".format(*head))
    print("-" * 69)
    for el in els:
        g, tried = run(el)
        if not g:
            print(f"{el:4s}  (no output)")
            continue
        e = refdata.ELEMENTS[el]
        da = (100.0 * (g["a0"] - e["a0"]) / e["a0"]
              if g.get("a0") else float("nan"))
        g["tried"] = {t["file"]: round(t["rms"], 2) for t in tried}
        out[el] = g
        v = lib[el]
        tap = (v.get("tap") or {}).get("rms", float("nan"))
        print(f"{el:4s}{g['file']:>24s}{len(tried):4d}{da:8.2f}"
              f"{g['rms']:8.2f}{v['rms']:10.2f}{tap:11.2f}")
    json.dump(out, open(os.path.join(HERE, "eam_results.json"), "w"),
              indent=1, sort_keys=True)
    print("\nyazildi:", os.path.join(HERE, "eam_results.json"))


if __name__ == "__main__":
    main()
