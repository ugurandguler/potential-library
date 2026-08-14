#!/usr/bin/env python3
"""
Merge the UG (angular) results into library.json so the viewer can put MAU and
UG side by side.

The two trees define incompatible versions of latdyn, so nothing is imported
across the boundary: angular/export_ug.py writes ug_results.json and this reads
it.  If that file is absent the step is a no-op and the page simply shows MAU
alone, which is the correct behaviour before any angular run exists.

**The comparison is only meaningful when both sides used the same three-body
cutoff.**  Extending it moves niobium from 216.7 % to 24.8 % with no angular
term at all, so a comparison across truncations measures the truncation.  Each
merged record therefore carries the cutoff it was fitted with, and the viewer
refuses to draw the comparison as a verdict when the two differ.

What "the same" means had to be settled by measurement rather than by comparing
floats.  The two radii disagreed in the fifth decimal for seven hcp elements -
hafnium at 4.687988 against 4.688042 - because their lattice constants were
rounded differently between the run and the reference table, and the difference
scattered in both directions, so it was rounding and not a data revision.  An
absolute tolerance of 1e-6 called those incomparable, and the viewer duly
dropped UG from the potential plot for Hf, Lu, Re, Ru, Sc, Tl and Y.

The quantity that actually matters is which neighbours the three-body sum
covers, so that is what is tested: the two radii are the same cutoff when the
same neighbours fall inside them.  For all seven that is 36 neighbours and 153
triplets at either radius - the same sum, to the atom.  A shell gap is a few per
cent and lattice-constant rounding is 1e-5, so nothing sits close to the line.

    python add_ug_overlay.py
"""
import json
import os

import latdyn as L
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), "angular", "ug_results.json")


def same_cutoff(el, r_mau, r_ug):
    """Do the two radii enclose the same neighbours?"""
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    lo, hi = sorted((r_mau, r_ug))
    #  anything strictly between the two radii is a neighbour one sum has and
    #  the other does not; if there is none, the two sums are identical
    return not any(lo < r <= hi for (_, _, _, _, r) in L.neighbours(cry, hi))


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))

    if not os.path.exists(SRC):
        print("ug_results.json is missing; UG comparison skipped "
              "(angular/export_ug.py has not been run)")
        return

    ug = json.load(open(SRC))
    print(f"{'el':4s}{'MAU rms':>10s}{'UG rms':>9s}{'gain':>9s}"
          f"{'MAU kar.':>10s}{'UG kar.':>9s}{'rcut3':>18s}")
    print("-" * 66)

    merged = comparable = 0
    for el, u in sorted(ug.items()):
        v = lib.get(el)
        if not v:
            continue
        same_cut = same_cutoff(el, v["rcut3"], u["rcut3"])
        u = dict(u)
        #  a comparison across different truncations measures the truncation,
        #  so the flag travels with the data rather than living in the viewer
        u["comparable"] = bool(same_cut)
        u["rcut3_mau"] = v["rcut3"]
        v["ug"] = u
        merged += 1
        comparable += same_cut
        m_stab = (v.get("dyn") or {}).get("stable")
        cut = ("same" if same_cut
               else f"{v['rcut3']:.2f} vs {u['rcut3']:.2f}")
        print(f"{el:4s}{v['rms']:9.2f}%{u['rms']:8.2f}%"
              f"{v['rms'] - u['rms']:+9.2f}"
              f"{str(m_stab):>10s}{str(u['dyn']['stable']):>9s}{cut:>18s}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\nUG added to {merged} elements, {comparable} of them at the same "
          f"cutoff radius")
    if comparable < merged:
        print("WARNING: where the cutoff radii differ, the difference carries "
              "the cutoff's effect as well as the angular term's")
    print("merged:", path)


if __name__ == "__main__":
    main()
