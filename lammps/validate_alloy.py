#!/usr/bin/env python3
"""
Does the multi-species path give the same answer as the single-species one?

There is no independent implementation of a mixed-species ugur potential to
check against, so the check is built out of one that exists: **declare the same
element twice**.  Write a two-element file in which A and B carry identical
parameters on all eight triples, take the crystal that was already validated,
and label its atoms alternately type 1 and type 2.  Physically nothing has
changed, so the energy, the pressure and the forces have to be bit-for-bit what
the one-type run gives.

That sounds weak and is not.  Everything the alloy path added runs in this
test and nothing else does: `map` from type to element, `elem3param` indexing,
`pair_par(i,j)` picking the two-body entry off the (i,j,j) line, the per-leg
cutoffs rcut3_AB and rcut3_AC, and the rule that a triple takes its radial
shape from the centre.  A wrong index in any of them reaches a different entry
and, because every entry here holds the same numbers, still returns the right
answer only if the indexing is right for a reason rather than by luck - so any
*structural* error shows up, while an error that merely picks the wrong element
cannot hide.

The second test is the one that would catch that last class: entry (A,B,B) is
given a deliberately wrong value and the run must then disagree.  If it still
agrees, the code is not reading that entry at all.

    python validate_alloy.py            # Cu, then Fe
    python validate_alloy.py Nb
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
SKIN = 2.0

LINE = ("{a} {b} {c} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} "
        "{C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} 0 0")

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
variable        prs equal press
print           "UGUR_E ${{epa}}"
print           "UGUR_P ${{prs}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def write_twin(path, el, rec, spoil=None):
    """A and B identical; `spoil` multiplies C on the (A,B,B) entry only"""
    p = {k: rec[k] for k in ("m", "D", "alpha", "r0", "gamma",
                             "C", "alpha3", "rcut2", "rcut3")}
    p["taper"] = rec.get("taper") or -1.0
    out = [f"# twin test for {el}: two labels, one element"]
    for a in ("A", "B"):
        for b in ("A", "B"):
            for c in ("A", "B"):
                q = dict(p)
                if spoil is not None and (a, b, c) == ("A", "B", "B"):
                    q["C"] = p["C"] * spoil
                out.append(LINE.format(a=a, b=b, c=c, **q))
    open(path, "w").write("\n".join(out) + "\n")


def write_data_two_types(path, cry, mass, rep):
    """same cell as cellfile.write_data, atoms alternating between two types"""
    cellfile.write_data(path, cry, mass, rep)
    lines = open(path).read().splitlines()
    out, n = [], 0
    for ln in lines:
        if ln.strip() == "1 atom types":
            out.append("2 atom types")
            continue
        if ln.startswith("1 ") and "Masses" in "\n".join(out[-4:]):
            out.append(ln)
            out.append(ln.replace("1 ", "2 ", 1))
            continue
        parts = ln.split()
        if len(parts) == 5 and parts[1] == "1" and parts[0].isdigit():
            n += 1
            out.append(f"{parts[0]} {1 + (n % 2)} "
                       f"{parts[2]} {parts[3]} {parts[4]}")
            continue
        out.append(ln)
    open(path, "w").write("\n".join(out) + "\n")
    return n


def run(el, rec, cry, tag, spoil=None):
    d = os.path.join(HERE, "alloyruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    write_twin(os.path.join(d, "twin.ugur"), el, rec, spoil)
    rcut = max(rec["rcut2"], rec["rcut3"])
    box, _ = cellfile.orthogonal_cell(cry)
    rep = tuple(max(3, int(np.ceil(2.0 * (rcut + SKIN) / b))) for b in box)
    n = write_data_two_types(os.path.join(d, f"{el}.data"), cry,
                             refdata.MASSES[el], rep)
    open(os.path.join(d, "in.check"), "w").write(IN.format(
        data=f"{el}.data", skin=SKIN, pot="twin.ugur", els="A B"))
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {wsl(d)} && {LMP} -in in.check 2>&1"],
                       capture_output=True, text=True)
    g = {"nat": n}
    for k, pat in (("E", r"UGUR_E\s+([-\d.eE+]+)"),
                   ("P", r"UGUR_P\s+([-\d.eE+]+)")):
        m = re.search(pat, r.stdout)
        if m:
            g[k] = float(m.group(1))
    if "E" not in g:
        err = [l for l in r.stdout.splitlines() if "ERROR" in l]
        g["err"] = err[0][:90] if err else "output unreadable"
    return g


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or ["Cu", "Fe", "Ti"]
    print("Twin test: the same element labelled as two types.\n")
    print(f"{'el':4s}{'atom':>7s}{'E one type':>13s}{'E two types':>13s}"
          f"{'diff':>11s}{'C spoiled':>13s}   status")
    print("-" * 68)
    bad = 0
    for el in els:
        e = refdata.ELEMENTS[el]
        rec = lib[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        ref = L.energy(cry, pot=L.Potential.from_record(rec))
        g = run(el, rec, cry, "twin")
        if "err" in g:
            print(f"{el:4s}  error: {g['err']}")
            bad += 1
            continue
        d = g["E"] - ref
        #  and now break one entry: if the answer does not move, that entry is
        #  never read and the agreement above meant nothing
        s = run(el, rec, cry, "spoil", spoil=1.5)
        moved = abs(s.get("E", g["E"]) - g["E"])
        ok = abs(d) < 1e-9 * max(abs(ref), 1.0) and moved > 1e-6
        bad += not ok
        print(f"{el:4s}{g['nat']:7d}{ref:13.8f}{g['E']:13.8f}{d:11.2e}"
              f"{moved:13.2e}   {'ok' if ok else 'FAILED'}")
    print()
    if bad:
        print(f"{bad} element uyusmadi")
        raise SystemExit(1)
    print("the multi-species path gives the single-species energy, and")
    print("bozulan girdi sonucu degistiriyor - yani gercekten okunuyor")


if __name__ == "__main__":
    main()
