#!/usr/bin/env python3
"""
Fold the structure gate into the library: does each record KEEP its crystal
structure at temperature once the box is free to change shape?

The molecular-dynamics screen (lammps/md_screen.py) cannot ask that.  It runs
constant-energy dynamics in a box whose shape is fixed, and a Bain shear or a
collapse to close packing IS a change of box shape - so the screen is
structurally blind to the failures that matter most for the body-centred
metals.  Found 2026-09-30: niobium's switched record, which passes the screen,
leaves bcc within four picoseconds at 296 K and settles 200 meV/atom below its
own static lattice.  The finite-temperature dispersion of five elements had
been measured on crystals that were no longer bcc.

The gate (lammps/gate/, run on TRUBA, jobs 6428949-51):

    0 K       box/relax with all six cell degrees of freedom -> E0
    heat      NPT with a fully triclinic barostat at T = min(300 K, 0.6 Tm),
              40 ps, the last 20 averaged
    quench    minimise again, box free
    CNA       the fraction of atoms in the reference structure, fixed cutoff
              1.207 x (bcc lattice parameter, or nearest-neighbour distance
              for fcc and hcp), during the run and after the quench

    FAIL      the quenched crystal lies more than 2 meV/atom BELOW E0, or
              fewer than 90 % of its atoms are still in the reference
              structure after the quench

An energy below the static lattice cannot be a heated version of that
lattice; it is something else, and the CNA says what fraction is left.

HARD-CUT records are measured and stored but carry no verdict.  A hard
truncation leaves the pair energy discontinuous, so the quench can land a
neighbour shell on the other side of the cut and report an energy jump of
electron-volts with the structure intact (hafnium's hard + angular record:
100 % hcp, -3.2 eV).  Those records are not for dynamics anyway, and the page
says so where they are listed.

    python add_gate.py
    python add_gate.py --dry
"""
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "lammps", "gate", "gate.json")

#  gate.json names arms by potential-file suffix; the library by record key.
#  "hard" is the element's top-level record itself.
KEY = {"hard": None, "hard_ug": "ug", "tap": "tap", "tap_ug": "tap_ug",
       "rc": "rc", "rc_ug": "rc_ug", "tap_force": "tap_force",
       "nudge": "tap_nudge", "disp": "hard_disp"}
HARD = {"hard", "hard_ug", "disp"}
#  NOT "gate": that name is taken twice already - tap_force["gate"] is the
#  three-body ceiling verdict ("pass"/"reject") and model_curve["gate"] the
#  quantity a reconstruction was checked against.  A first version of this
#  file cleared "gate" everywhere and deleted both (caught by a key-path diff).
FIELD = "structure_gate"


def main():
    dry = "--dry" in sys.argv
    if not os.path.exists(SRC):
        print(f"{SRC} is missing - nothing folded in")
        return
    lib = json.load(open(os.path.join(HERE, "library.json")))
    gate = json.load(open(SRC))

    #  clear first: a record should carry a gate result only if this file has one
    for d in lib.values():
        if isinstance(d, dict):
            d.pop(FIELD, None)
            for sub in d.values():
                if isinstance(sub, dict):
                    sub.pop(FIELD, None)

    tally, fails, missing = Counter(), {}, []
    for name, g in sorted(gate.items()):
        el, arm = g["el"], g["arm"]
        d = lib.get(el)
        key = KEY.get(arm, "?")
        rec = d if key is None else (d or {}).get(key)
        if not isinstance(rec, dict):
            missing.append(name)
            continue
        out = {"T": g["T"], "struct": g["struct"]}
        if g["verdict"] == "NORUN":
            out["verdict"] = "crashed"
            out["error"] = g.get("err", "")
        else:
            out.update({
                "frac_md": round(g["fmd"], 3),      # reference fraction, last 20 ps
                "frac_quench": round(g["fq"], 3),
                "dE_md_meV": round(g["dE_md"], 1),  # <PE> - E0, per atom
                "dE_quench_meV": round(g["dE_q"], 1),
                "kT15_meV": round(g["kT15"], 1),
            })
            out["verdict"] = ("na_hard" if arm in HARD else
                              ("fail" if g["verdict"] == "FAIL" else "pass"))
        rec[FIELD] = out
        tally[(arm, out["verdict"])] += 1
        if out["verdict"] in ("fail", "crashed"):
            fails.setdefault(arm, []).append(el)

    arms = sorted({a for a, _ in tally})
    for a in arms:
        row = {v: n for (b, v), n in tally.items() if b == a}
        print(f"  {a:10s} {row}")
    for a, els in sorted(fails.items()):
        print(f"  FAIL {a:10s} {' '.join(sorted(els))}")
    if missing:
        print(f"  no library record for {len(missing)}: {' '.join(missing)}")
    if dry:
        print("--dry: nothing written")
        return
    tmp = os.path.join(HERE, "library.json.tmp")
    json.dump(lib, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(HERE, "library.json"))
    print("-> library.json")


if __name__ == "__main__":
    main()
