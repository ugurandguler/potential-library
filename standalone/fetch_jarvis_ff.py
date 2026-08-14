#!/usr/bin/env python3
"""
Pull classical-potential elastic constants from JARVIS-FF (NIST).

Why this rather than running LAMMPS ourselves.  `lammps/eam_compare.py` does
run it, and covers eleven elements, but two things about that are unsatisfying:
the numbers are ours, so "you ran EAM badly" is an available objection; and
where several parameter sets ship, one has to be chosen.  That choice is not
small - aluminium reads 5.58 % with Mendelev's set and 26.54 % with Zhou's, so
picking the wrong one flatters us by a factor of five.

JARVIS-FF removes both.  Every publicly available potential is run through one
standardised workflow by a third party, so the comparison is against the whole
published distribution rather than against a set we picked, and the elastic
constants are not ours to get wrong.  Copper has 73 entries, nickel 62,
aluminium 50.

    python fetch_jarvis_ff.py            # fetch, check, write jarvis_ff.json
    python fetch_jarvis_ff.py --check    # convention check only

What is NOT covered, and it is a real gap: sixteen of our thirty-eight have no
entry at all - Ba, Be, Ca, Cs, Ir, K, Li, Lu, Rb, Re, Rh, Sc, Sr, Tl, Y, Yb.
Alkalis, alkaline earths and rare earths, which is where the classical-potential
literature is thin to begin with.  That absence is itself worth reporting.
"""
import json
import os
import re
import sys

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "jarvis_ff.json")

#  Which family a potential belongs to, from its name.  The comparison has to
#  be like with like: MAU has no angular term and belongs against EAM, while UG
#  carries one and belongs against MEAM and ADP.  Anything machine-learned is
#  kept separate - it is a different kind of object and beating it or losing to
#  it says nothing about a seven-parameter analytic form.
FAMILY = (
    (r"\bmeam\b", "MEAM"),
    (r"\badp\b", "ADP"),
    (r"\beam\b|\bfs\b|finnis", "EAM"),
    (r"tersoff|\bsw\b|stillinger|vashishta|comb|edip|bop", "bond-order"),
    (r"snap|gap|nnp|mlip|ace|mtp|alignn|\bml\b", "machine-learned"),
    (r"lj|lennard|morse|buck", "pair"),
)


def family(func):
    f = (func or "").lower()
    for pat, name in FAMILY:
        if re.search(pat, f):
            return name
    return "other"


def cubic_from_tensor(et):
    """C11, C12, C44 from a JARVIS-FF elastic record, GPa, or None"""
    raw = et.get("raw_et_tensor") if isinstance(et, dict) else None
    if not raw or len(raw) < 6:
        return None
    try:
        c11 = float(raw[0][0])
        c12 = float(raw[0][1])
        c44 = float(raw[3][3])
    except (TypeError, ValueError, IndexError):
        return None
    if not all(abs(x) < 5000 for x in (c11, c12, c44)):
        return None
    return c11, c12, c44


def main():
    from jarvis.db.figshare import data
    lib = json.load(open(os.path.join(HERE, "library.json")))
    d = data("jff")

    out = {}
    for r in d:
        el = r.get("formula")
        if el not in lib:
            continue
        c = cubic_from_tensor(r.get("elastic_tensor_data") or {})
        if not c:
            continue
        out.setdefault(el, []).append(
            {"func": r.get("func", ""), "family": family(r.get("func")),
             "C11": c[0], "C12": c[1], "C44": c[2], "jid": r.get("jid")})

    #  Convention check.  Their tensor is computed by their own workflow and
    #  there is no point comparing against it until it is known to mean the
    #  same thing.  The test that matters: does the BEST potential for a
    #  well-studied element land near experiment?  If the whole distribution is
    #  offset, the convention differs and the comparison is void.
    print("RULE CHECK - the potential closest to experiment for each element")
    print(f"{'el':4s}{'n':>4s}{'en iyi rms %':>14s}{'aile':>16s}  potansiyel")
    print("-" * 74)
    for el in ("Cu", "Ni", "Al", "Fe", "Ta"):
        rows = out.get(el, [])
        if not rows:
            continue
        e = refdata.ELEMENTS[el]["Cij"]
        for row in rows:
            errs = [(row[k] - e[k]) / e[k] for k in ("C11", "C12", "C44")]
            row["rms"] = 100.0 * (sum(x * x for x in errs) / 3) ** 0.5
        best = min(rows, key=lambda x: x["rms"])
        print(f"{el:4s}{len(rows):4d}{best['rms']:14.2f}{best['family']:>16s}"
              f"  {best['func'][:34]}")

    if "--check" in sys.argv:
        return
    for el, rows in out.items():
        e = refdata.ELEMENTS[el]["Cij"]
        for row in rows:
            if "rms" not in row:
                errs = [(row[k] - e[k]) / e[k] for k in ("C11", "C12", "C44")]
                row["rms"] = 100.0 * (sum(x * x for x in errs) / 3) ** 0.5
    json.dump(out, open(OUT, "w"), indent=1, sort_keys=True)
    n = sum(len(v) for v in out.values())
    print(f"\n{len(out)} elements, {n} potential records -> {OUT}")
    print("kapsanmayan:",
          " ".join(sorted(e for e in lib if e not in out)))


if __name__ == "__main__":
    main()
