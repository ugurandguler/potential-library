#!/usr/bin/env python3
"""
Surface energies to test against: density functional per facet, and experiment.

Written before the calculation that will be compared with it, and with a
prediction recorded, because this project has already produced one uncalibrated
negative result and it was worth nothing until a baseline was attached.

**The prediction.** The vacancy formation energy came out at 0.7-0.8 of the
cohesive energy where a real metal gives 0.2-0.35, and the reason is structural
rather than a fitting failure: in a pair-like form the energy of a bond does not
depend on how many other bonds an atom has, so removing an atom costs the full
sum of its bonds and the neighbours gain nothing back.  A surface is the same
physics - atoms with fewer neighbours - so the surface energies should come out
high by a similar factor, roughly two to three.  If they do, the six pathologies
and the vacancy all reduce to one missing ingredient.  If they do not, there is
something else going on and it is worth knowing.

Two references, because they fail differently:

  DFT      Materials Project's surface database (Crystalium; Tran, Xu, Radhakrishnan,
           Winston, Sun, Persson, Ong, Sci. Data 3, 160080 (2016)).  Per facet,
           which is what a potential should be judged on - the anisotropy is a
           stronger test than the average, since any form can be scaled to hit
           one number.  Systematically low by 10-30 % against experiment, in the
           way GGA surface energies always are.

  experiment  Tyson and Miller, Surf. Sci. 62, 267 (1977): liquid-metal surface
           tension measured against temperature and extrapolated to the melting
           point, then to the solid.  One number per element, no facet
           resolution, and the extrapolation is the weak part.  Quoted because
           it is the standard experimental compilation and because it brackets
           the DFT from the other side.

    python fetch_surface.py
    python fetch_surface.py Pd Cu Al
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import fetch_mp      # noqa: E402  (key handling; the key is never printed)
import refdata       # noqa: E402

#  Tyson & Miller, Surf. Sci. 62, 267 (1977), Table 1: gamma_s at the melting
#  point, J/m^2, from liquid surface tension extrapolated to the solid.  These
#  are the numbers everyone compares potentials against; they carry an
#  uncertainty the paper puts at a few per cent for the well-measured metals and
#  considerably more for the alkalis and the refractories.
TYSON = {
    "Ag": 1.25, "Al": 1.14, "Au": 1.50, "Ba": 0.38, "Be": 1.63, "Ca": 0.49,
    "Co": 2.52, "Cr": 2.30, "Cs": 0.10, "Cu": 1.79, "Fe": 2.42, "Hf": 2.20,
    "Ir": 3.05, "K": 0.13, "Li": 0.52, "Lu": 0.90, "Mg": 0.79, "Mo": 2.91,
    "Na": 0.26, "Nb": 2.66, "Ni": 2.45, "Pb": 0.59, "Pd": 2.00, "Pt": 2.49,
    "Rb": 0.11, "Re": 3.63, "Rh": 2.70, "Ru": 3.05, "Sc": 1.20, "Sr": 0.41,
    "Ta": 2.90, "Ti": 2.10, "Tl": 0.58, "V": 2.55, "W": 3.27, "Y": 1.13,
    "Yb": 0.50, "Zr": 1.91,
}

#  the facets worth having per structure: the close-packed one, which is the
#  lowest, and the two that bracket it.  A potential that gets the ordering
#  wrong is more interesting than one that gets the magnitude wrong.
#  space group of the structure this library fits, so the reference is the same
#  crystal as the calculation: fcc Fm-3m, bcc Im-3m, hcp P6_3/mmc
SG = {"fcc": 225, "bcc": 229, "hcp": 194}

FACETS = {
    "fcc": [(1, 1, 1), (1, 0, 0), (1, 1, 0)],
    "bcc": [(1, 1, 0), (1, 0, 0), (1, 1, 1)],
    "hcp": [(0, 0, 1), (1, 0, 0), (1, 1, 0)],
}


def main():
    els = sys.argv[1:] or sorted(refdata.ELEMENTS)
    ids = json.load(open(os.path.join(HERE, "jarvis_ids.json"))) \
        if os.path.exists(os.path.join(HERE, "jarvis_ids.json")) else {}
    out = {}
    p = os.path.join(HERE, "surface_ref.json")
    if os.path.exists(p):
        try:
            out = json.load(open(p))
        except Exception:
            out = {}

    from mp_api.client import MPRester
    with MPRester(fetch_mp.api_key()) as mpr:
        for el in els:
            if el not in refdata.ELEMENTS:
                continue
            struct = refdata.ELEMENTS[el]["struct"]
            try:
                #  Chosen by SPACE GROUP first and hull distance second.  Hull
                #  distance alone returns whichever polymorph is lowest in MP's
                #  own calculation, and that is not always this library's
                #  structure: silver came back as P6_3/mmc, sodium and
                #  strontium as hexagonal, titanium as P6/mmm.  A surface
                #  energy belongs to a crystal, so comparing those with an fcc
                #  or bcc calculation would be comparing two different metals -
                #  the same mistake the AFLOW entries were caught making.
                docs = mpr.materials.summary.search(
                    formula=el, fields=["material_id", "symmetry",
                                        "energy_above_hull", "structure"])
                docs = [d for d in docs if d.energy_above_hull is not None]
                if not docs:
                    #  no density functional entry, but the experimental value
                    #  is a reference in its own right and ytterbium has one
                    if TYSON.get(el):
                        out[el] = {"struct": struct, "tyson": TYSON[el],
                                   "mp_id": None, "spacegroup": None}
                        print(f"{el:3s} no MP record; experiment "
                              f"kept against {TYSON[el]:.2f}")
                    else:
                        print(f"{el:3s} no MP record")
                    continue
                want_sg = SG[struct]
                right = [d for d in docs
                         if d.symmetry and d.symmetry.number == want_sg]
                if not right:
                    docs.sort(key=lambda d: d.energy_above_hull)
                    sgn = docs[0].symmetry.number if docs[0].symmetry else "?"
                    print(f"{el:3s} no record for {struct} (space group {want_sg}) "
                          f"missing; the lowest is {sgn} - skipped")
                    out.pop(el, None)
                    continue
                right.sort(key=lambda d: d.energy_above_hull)
                mid = str(right[0].material_id)
                sg = right[0].symmetry.symbol
                hull = right[0].energy_above_hull
                sp = mpr.materials.surface_properties.search(
                    material_ids=[mid],
                    fields=["material_id", "weighted_surface_energy",
                            "surface_anisotropy", "surfaces"])
            except Exception as ex:
                print(f"{el:3s} error: {type(ex).__name__} {str(ex)[:70]}")
                continue
            rec = {"mp_id": mid, "spacegroup": sg, "struct": struct,
                   "hull": hull, "tyson": TYSON.get(el)}
            if sp:
                d = sp[0]
                rec["weighted"] = d.weighted_surface_energy
                rec["anisotropy"] = d.surface_anisotropy
                rec["facets"] = {
                    "".join(str(i) for i in s.miller_index): s.surface_energy
                    for s in (d.surfaces or [])}
            out[el] = rec
            if not rec.get("facets"):
                #  No facet-resolved DFT for this entry.  The record is KEPT if
                #  it still carries the experimental value, because the test
                #  that matters - whether the close-packed face is the cheapest
                #  - is a universal ordering and needs no per-element
                #  reference at all.  Dropping these silently removed cobalt,
                #  rhenium and ytterbium from the sweep entirely, which is not
                #  the same as having nothing to say about them.
                rec.pop("weighted", None)
                rec.pop("anisotropy", None)
                if not rec.get("tyson"):
                    print(f"{el:3s} {mid:12s} {sg:8s} neither DFT nor experiment - skipped")
                    out.pop(el, None)
                    continue
                out[el] = rec
                print(f"{el:3s} {mid:12s} {sg:8s} no facet DFT; "
                      f"kept against experiment {rec['tyson']:.2f}")
                continue
            f = rec.get("facets") or {}
            want = ["".join(str(i) for i in m) for m in FACETS[struct]]
            got = " ".join(f"{k}={f[k]:.3f}" for k in want if k in f)
            print(f"{el:3s} {mid:12s} {sg:8s} agirlikli "
                  f"{(rec.get('weighted') or 0):5.3f}  expt "
                  f"{(rec.get('tyson') or 0):5.2f}   {got}")

    tmp = p + ".tmp"
    json.dump(out, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, p)
    nf = sum(1 for r in out.values() if r.get("facets"))
    print(f"\n{len(out)} elements, {nf} of them facet-resolved, "
          f"{sum(1 for r in out.values() if r.get('tyson'))} with experiment")
    print("->", p)


if __name__ == "__main__":
    main()
