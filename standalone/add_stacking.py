#!/usr/bin/env python3
"""
Fold the generalised stacking fault into the library, with its two controls.

This is the test that turns a ranking into a consequence.  The ground-state
result said seventy-two of seventy-six records prefer a structure other than
the fitted one, and for the face-centred metals the preferred one is
hexagonal.  An intrinsic stacking fault IS a slab of hexagonal stacking inside
a face-centred crystal, so that finding has a falsifiable prediction attached:

    gamma_isf ~ 2 (E_hcp - E_fcc) / A     must come out NEGATIVE

and a negative intrinsic fault energy is not a small error, it is a change of
sign.  In a real face-centred metal a fault costs energy and heals; at a
negative value the faulted crystal is LOWER than the perfect one, so any fault
that forms stays, and the partial dislocations bounding it repel without
limit.  Copper is the case to quote: experiment +45 mJ/m^2, ours -63.

Two controls travel with it, because a wrong sign is exactly the kind of
result that deserves suspicion of the machinery first.

  the published potentials go through the identical code.  Five of five
  copper baselines come back positive, from 13.6 to 86.3 mJ/m^2 - so the code
  does not manufacture negative fault energies.

  NIST computes the same quantity independently, by a different geometry (a
  relaxed free-surface slab, sliced; ours is a fully periodic cell with the
  box tilt matched to the shift).  Where the potential is confirmed to be the
  same file, the two agree to four significant figures.

That last qualification is load-bearing and is enforced rather than assumed -
see nist_match.py.  Cu_zhou.eam.alloy would otherwise have contributed a 4x
disagreement that says nothing about either calculation, because the file
LAMMPS ships under that name is not the one NIST hosts.

    python add_stacking.py
    python add_stacking.py --dry
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nist_match                                     # noqa: E402
import refdata                                        # noqa: E402

EV_A2_TO_MJ_M2 = 16021.766208

SRC = os.path.join(HERE, "..", "lammps", "stacking.json")
EXP = {                       # mJ/m^2, intrinsic fault, room temperature
    "Ag": 16.0, "Al": 166.0, "Au": 32.0, "Cu": 45.0,
    "Ni": 125.0, "Pd": 180.0, "Pt": 322.0, "Rh": 750.0, "Ir": 480.0,
    "Pb": 30.0,
}


#  the prototype and plane our fault lives on, per crystal structure
PLANE = {"fcc": ("A1", "111"), "hcp": ("A3", "0001")}


def nist_isf(props, el, fn, struct="fcc"):
    """NIST's intrinsic and unstable fault for one baseline file

    Keyed on the prototype that IS this element's crystal.  NIST computes
    every structure it can relax for every potential, so an aluminium file
    carries an A3 basal fault as well as its A1 one; comparing our hexagonal
    result against their cubic entry, or the reverse, would be comparing two
    different quantities that happen to have the same name.
    """
    want = PLANE.get(struct)
    if not want:
        return None, None
    pre, plane = want
    r = (props.get(f"{el}|{fn}") or {}).get("stacking") or {}
    for proto, planes in r.items():
        if proto.split("--")[0] != pre:
            continue
        d = planes.get(plane)
        if not d:
            continue
        isf = d.get("E_isf")
        usf = next((v for k, v in d.items() if k.startswith("E_usf")), None)
        return isf, usf
    return None, None


def predicted(rank, a0, struct="fcc"):
    """gamma_isf expected from the bulk hcp-fcc energy difference alone

    The fault is two layers of the OTHER close-packed stacking, so if a bond's
    energy really does not know how many neighbours an atom has, the fault
    should cost whatever those two layers cost in bulk:

        fcc:  gamma ~ 2 (E_hcp - E_fcc) / A
        hcp:  gamma ~ 2 (E_fcc - E_hcp) / A

    The area is per atom in the fault plane, and it is NOT the same expression
    for the two structures even though both planes are triangular lattices.
    The lattice constant means different things: in fcc the in-plane spacing
    is a/sqrt(2), giving sqrt(3)/4 a^2, while in hcp it is a itself, giving
    sqrt(3)/2 a^2 - twice as much. Carrying the cubic formula over to the
    hexagonal elements made every hcp prediction exactly twice too large,
    which looked like a physical shortfall (magnesium predicted -65 against a
    measured -31) rather than an arithmetic one. It was caught because the
    cubic side agreed to 1 % and there was no reason for the hexagonal side to
    be systematically half.

    This shares no code with the shift calculation - one relaxes bulk
    crystals, the other slides half a slab - so agreement is a real check.
    """
    if not rank or "hcp" not in rank or "fcc" not in rank:
        return None
    dE = ((rank["fcc"]["E"] - rank["hcp"]["E"]) if struct == "hcp"
          else (rank["hcp"]["E"] - rank["fcc"]["E"]))
    area = (math.sqrt(3.0) / 2.0 if struct == "hcp"
            else math.sqrt(3.0) / 4.0) * a0 * a0
    return 2.0 * dE / area * EV_A2_TO_MJ_M2


def main():
    dry = "--dry" in sys.argv
    lib = json.load(open(os.path.join(HERE, "library.json")))
    if not os.path.exists(SRC):
        print(f"{SRC} is missing")
        return
    st = json.load(open(SRC))
    props = json.load(open(os.path.join(HERE, "nist_props.json")))
    ok = nist_match.confirmed()
    rp = os.path.join(HERE, "..", "lammps", "struct_rank.json")
    rank = json.load(open(rp)) if os.path.exists(rp) else {}

    n = 0
    ours, base_rows, checks, preds = [], [], [], []
    for key, v in sorted(st.items()):
        if "isf" not in v:
            continue
        el, tag = key.split("|", 1)
        rec = {"isf": v["isf"], "usf": v["usf"],
               "back_barrier": v.get("back_barrier"),
               "frac": v["frac"], "gamma": v["gamma"],
               "closure": v.get("closure"), "atoms": v.get("atoms")}
        if tag.startswith("base|"):
            fn = tag[5:]
            ni, nu = nist_isf(props, el, fn,
                              refdata.ELEMENTS.get(el, {}).get("struct", "fcc"))
            #  only where the file has been confirmed to be the same one
            if ni is not None and f"{el}|{fn}" in ok:
                rec["nist_isf"], rec["nist_usf"] = ni, nu
                checks.append((abs(v["isf"] - ni), el, fn, v["isf"], ni))
            base_rows.append((el, fn, v["isf"], v["usf"]))
            #  same shape the surface section uses, so the page reads them
            #  the same way
            lib.setdefault(el, {}).setdefault("baseline_stacking", {})[fn] = rec
        else:
            if el not in lib or not isinstance(lib[el].get(tag), dict):
                continue
            rec["exp"] = EXP.get(el)
            a0 = (lib[el].get("a0") or {})
            a0 = a0 if isinstance(a0, float) else lib[el].get("a0")
            #  NOT `st` - that name already holds the whole results file two
            #  scopes up.  Shadowing one like it silently dropped every
            #  surface reference once already in this project.
            struct = refdata.ELEMENTS.get(el, {}).get("struct", "fcc")
            p = predicted(rank.get(key), float(a0), struct) if a0 else None
            if p is not None:
                rec["predicted"] = p
                preds.append((p, v["isf"], key))
            lib[el][tag]["stacking"] = rec
            ours.append((v["isf"], el, tag))
        n += 1

    print(f"{n} records ({len(ours)} ours, {len(base_rows)} baseline)")

    neg = [x for x in ours if x[0] < 0]
    print(f"\nOURS: intrinsic stacking fault in {len(neg)}/{len(ours)} records "
          f"NEGATIF")
    if ours:
        ours.sort()
        print(f"  range {ours[0][0]:.1f} ({ours[0][1]}) .. "
              f"{ours[-1][0]:.1f} ({ours[-1][1]}) mJ/m2")
    #  Say what the baselines actually did, not what would be convenient.
    #  Nineteen of twenty are positive and one - Al_jnp, at -4.1 mJ/m^2 - is
    #  not, so "the code cannot produce a negative fault energy" would be
    #  false.  What the baselines establish is weaker and still sufficient:
    #  the code does not produce negatives BY ITSELF, and the one it does
    #  produce is a tenth the size of ours and a known weakness of that
    #  particular aluminium potential.
    pos = [x for x in base_rows if x[2] > 0]
    negb = [x for x in base_rows if x[2] <= 0]
    print(f"TABAN : {len(pos)}/{len(base_rows)} kayitta POZITIF")
    for el, fn, isf, _ in negb:
        print(f"   negative baseline: {el} {fn} {isf:.1f} mJ/m2")

    if preds:
        sign = sum(1 for p, m, _ in preds if (p < 0) == (m < 0))
        rat = sorted(m / p for p, m, _ in preds if p)
        print(f"\nPREDICTION from the ground-state difference ({len(preds)} records):")
        print(f"  isaret uyusmasi {sign}/{len(preds)}")
        print(f"  measured/predicted: median {rat[len(rat)//2]:.2f}, "
              f"range {rat[0]:.2f}-{rat[-1]:.2f}")
        #  Carry the library-wide figures INTO each record so the page states
        #  them from data.  The surface section once carried a hand-typed
        #  "0.8-1.1" that had drifted to a true 0.48-1.97, and it flattered
        #  our own comparison for weeks before anybody re-derived it.
        summ = {"n": len(preds), "sign_ok": sign,
                "ratio_median": rat[len(rat) // 2],
                "ratio_lo": rat[0], "ratio_hi": rat[-1]}
        for _, _, key in preds:
            el, tag = key.split("|", 1)
            lib[el][tag]["stacking"]["library"] = summ

    if checks:
        checks.sort()
        d = [c[0] for c in checks]
        print(f"\nNIST comparison ({len(checks)} confirmed files):")
        print(f"  difference: median {d[len(d)//2]:.2f}, largest {d[-1]:.2f} mJ/m2")
        for dd, el, fn, o, ni in checks[-3:]:
            print(f"   {el:3s} {fn[:30]:30s} ours {o:8.2f}  NIST {ni:8.2f}")
    else:
        print("\nNIST comparison: no confirmed match")

    if dry:
        print("\n--dry: nothing written")
        return
    tmp = os.path.join(HERE, "library.json.tmp")
    json.dump(lib, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(HERE, "library.json"))
    print("\n-> library.json")


if __name__ == "__main__":
    main()
