#!/usr/bin/env python3
"""
Audit the hand-entered experimental table in refdata.py.

The fit is exact; its output is only as good as these targets, and every value
was typed in by hand from memory of standard references.  Three independent
checks, none of which needs an external source:

  1. INTERNAL IDENTITY.  For a cubic crystal the bulk modulus is not an extra
     measurement, it is B = (C11 + 2 C12)/3 exactly.  For hexagonal the Voigt
     average is B = [2(C11 + C12) + 4 C13 + C33]/9.  If the entered B and the
     entered Cij disagree, one of them is wrong - this cannot be argued away as
     "different experiments".

  2. MECHANICAL STABILITY.  Born criteria on the experimental Cij themselves.
     A real crystal satisfies them; a violation means a typo.

  3. DFT CROSS-CHECK.  Materials Project at 0 K versus these room-temperature
     numbers.  Not an identity - DFT and experiment legitimately differ by
     several per cent, and MP has its own failures - but a 50 % gap is a flag.

Prints a ranked list of everything worth re-reading in the source.

    python check_refdata.py
"""
import json, os
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


def b_from_cij(c, struct):
    """the bulk modulus implied by the elastic constants"""
    if struct == "hcp":
        need = ("C11", "C12", "C13", "C33")
        if not all(k in c for k in need):
            return None
        return (2*(c["C11"] + c["C12"]) + 4*c["C13"] + c["C33"]) / 9.0
    if not all(k in c for k in ("C11", "C12")):
        return None
    return (c["C11"] + 2*c["C12"]) / 3.0


def stability(c, struct):
    bad = []
    if struct == "hcp":
        if not c["C11"] > abs(c["C12"]):
            bad.append("C11 <= |C12|")
        if not (c["C11"] + c["C12"])*c["C33"] > 2*c["C13"]**2:
            bad.append("(C11+C12)C33 <= 2 C13^2")
        if c["C44"] <= 0:
            bad.append("C44 <= 0")
    else:
        if not c["C11"] > abs(c["C12"]):
            bad.append("C11 <= |C12|")
        if c["C11"] + 2*c["C12"] <= 0:
            bad.append("C11 + 2 C12 <= 0")
        if c["C44"] <= 0:
            bad.append("C44 <= 0")
    return bad


def main():
    lib = {}
    p = os.path.join(HERE, "library.json")
    if os.path.exists(p):
        lib = json.load(open(p))

    print("1. bulk modulus: entered value vs the one implied by the entered Cij")
    print(f"{'el':4s}{'str':5s}{'B entered':>11s}{'B from Cij':>12s}{'diff':>9s}")
    print("-"*42)
    flags = []
    for el, e in sorted(refdata.ELEMENTS.items()):
        bc = b_from_cij(e["Cij"], e["struct"])
        if bc is None:
            continue
        d = 100*(e["B"] - bc)/bc
        mark = ""
        if abs(d) > 8:
            mark = "  <<<"
            flags.append((abs(d), el, f"B {e['B']:.1f} vs {bc:.1f} from Cij ({d:+.0f}%)"))
        print(f"{el:4s}{e['struct']:5s}{e['B']:11.1f}{bc:12.1f}{d:8.1f}%{mark}")

    print("\n2. Born stability of the experimental constants themselves")
    any_bad = False
    for el, e in sorted(refdata.ELEMENTS.items()):
        bad = stability(e["Cij"], e["struct"])
        if bad:
            any_bad = True
            flags.append((999, el, "unstable: " + "; ".join(bad)))
            print(f"   {el}: {'; '.join(bad)}   {e['Cij']}")
    if not any_bad:
        print("   all pass")

    print("\n3. Materials Project (0 K DFT) vs the entered room-temperature values")
    print(f"{'el':4s}{'C11 exp':>9s}{'MP':>8s}{'diff':>8s}   "
          f"{'C44 exp':>9s}{'MP':>8s}{'diff':>8s}")
    print("-"*54)
    for el in sorted(refdata.ELEMENTS):
        me = lib.get(el, {}).get("mp", {}).get("elastic")
        if not me:
            continue
        e = refdata.ELEMENTS[el]["Cij"]
        row = f"{el:4s}"
        worst = 0.0
        for k in ("C11", "C44"):
            if k in e and k in me:
                d = 100*(me[k]-e[k])/e[k]
                worst = max(worst, abs(d))
                row += f"{e[k]:9.0f}{me[k]:8.0f}{d:7.0f}%"
            else:
                row += f"{'-':>9s}{'-':>8s}{'-':>8s}"
            row += "   " if k == "C11" else ""
        print(row + ("  <<<" if worst > 40 else ""))
        if worst > 40:
            flags.append((worst, el, f"DFT differs by {worst:.0f}% - check both"))

    print("\n" + "="*60)
    if flags:
        print("worth re-reading in the source, most suspicious first:")
        for _, el, why in sorted(flags, reverse=True):
            print(f"   {el:3s} {why}")
    else:
        print("nothing flagged")


if __name__ == "__main__":
    main()
