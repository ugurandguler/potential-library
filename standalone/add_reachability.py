#!/usr/bin/env python3
"""
Attach the reachability verdict to library.json.

Holding the bulk modulus fixed pins C11 + 2 C12, so a cubic tensor has exactly
two free constants left, C' = (C11 - C12)/2 and C44.  Scanning the parameter
space and recording every (C'/B, C44/B) the form can produce turns "why does
this element fit badly" into geometry: either the measured point is inside the
reachable set or it is not.

Two numbers go in, both computed elsewhere:

    R_floor.json        the lowest C44/C' any parameter set reaches, from
                        refine_R.py (grid, then simplex from two seed families,
                        keeping the smaller - each is only an upper bound)
    cprime_region.json  the reachable point cloud, from cprime_ceiling.py

The ratio test is necessary but not sufficient - Ta's ratio is reachable while
its point is not - so the distance to the cloud is carried too and is the one
that correlates with the fit quality.

    python add_reachability.py
"""
import json
import os

import numpy as np

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


def R_of_fit(rec):
    """C44/C' of a fitted tensor, or None if it is not a usable cubic one.

    build_library stores the constants as [ours, measured] pairs for the
    records the viewer draws and as bare numbers elsewhere, so both are read.
    """
    if not rec:
        return None
    try:
        c11, c12, c44 = rec["C11"], rec["C12"], rec["C44"]
    except (KeyError, TypeError):
        return None
    if isinstance(c11, list):
        c11, c12, c44 = c11[0], c12[0], c44[0]
    cp = 0.5 * (c11 - c12)
    return c44 / cp if cp > 0 else None


def stale_floors(floors):
    """True while R_floor.json is still the one measured at rcut3 = 1.12 d_nn.

    The library moved to 1.50 and the floors did not follow, so the verdict on
    the page would be the old form's.  Rather than let that pass silently the
    page carries a note, and the note clears itself: the check is a content
    comparison against the kept 1.12 file, so writing genuinely new floors
    makes it false without anyone remembering to.
    """
    old = os.path.join(HERE, "R_floor.rcut112.json")
    if not os.path.exists(old):
        return False
    ref = json.load(open(old))
    common = set(floors) & set(ref)
    if not common:
        return False
    return all(abs(floors[el]["R_floor"] - ref[el]["R_floor"]) < 1e-9
               for el in common)


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    floors = json.load(open(os.path.join(HERE, "R_floor.json")))
    region = json.load(open(os.path.join(HERE, "cprime_region.json")))
    stale = stale_floors(floors)
    if stale:
        print("R_floor.json is still the 1.12 d_nn measurement; the page will "
              "say so until it is replaced")

    print(f"{'el':4s}{'R_exp':>8s}{'floor':>8s}{'margin':>8s}{'2D gap':>8s}"
          f"{'R_MAU':>8s}{'R_UG':>8s}{'RMS %':>8s}  verdict")
    print("-" * 76)
    rows = []
    for el in sorted(lib):
        v = lib[el]
        f = floors.get(el)
        if not f:
            v.pop("reach", None)
            continue
        c = refdata.ELEMENTS[el]["Cij"]
        #  What the fits in the library actually reach.  The floor from
        #  refine_R.py is a search result and therefore an upper bound on the
        #  infimum, never a proof; a fit that lands lower is a better upper
        #  bound and replaces it.  Chromium needed this - its MAU fit sits at
        #  R = 2.338 where the search stopped at 2.384.
        r_mau = R_of_fit(v)
        r_ug = R_of_fit(v.get("ug"))
        floor = f["R_floor"]
        from_fit = r_mau is not None and r_mau < floor
        if from_fit:
            floor = r_mau
        reachable = bool(f["R_exp"] >= floor)
        gap = None
        if el in region:
            B = region[el]["B"]
            t = np.array([0.5 * (c["C11"] - c["C12"]) / B, c["C44"] / B])
            pts = np.array(region[el]["pts"])
            gap = float(np.min(np.linalg.norm(pts - t, axis=1))
                        / np.linalg.norm(t))
        v["reach"] = {"R_exp": round(f["R_exp"], 3),
                      "R_floor": round(floor, 3),
                      "margin": round(f["R_exp"] / floor, 3),
                      "gap": round(gap, 3) if gap is not None else None,
                      "ok": reachable,
                      "from_fit": from_fit,
                      "R_mau": round(r_mau, 3) if r_mau is not None else None,
                      "R_ug": round(r_ug, 3) if r_ug is not None else None,
                      #  the angular factor is a wider form, so its own fits
                      #  are entitled to sit below a floor measured without it
                      "ug_below": bool(r_ug is not None and r_ug < floor),
                      "stale": stale}
        rows.append((el, v["reach"], v["rms"]))
        print(f"{el:4s}{f['R_exp']:8.2f}{floor:8.3f}"
              f"{f['R_exp'] / floor:8.2f}"
              f"{(gap if gap is not None else float('nan')):8.2f}"
              f"{(r_mau if r_mau is not None else float('nan')):8.3f}"
              f"{(r_ug if r_ug is not None else float('nan')):8.3f}"
              f"{v['rms']:8.1f}  "
              f"{'reachable' if reachable else 'OUT OF REACH'}"
              f"{'  (floor from the fit)' if from_fit else ''}"
              f"{'  UG below' if v['reach']['ug_below'] else ''}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)

    ok = [r for _, r, _ in rows if r["ok"]]
    no = [r for _, r, _ in rows if not r["ok"]]
    mo = np.median([s for _, r, s in rows if r["ok"]])
    mn = np.median([s for _, r, s in rows if not r["ok"]])
    print(f"\nreachable   {len(ok):2d} elements, median RMS {mo:5.1f} %")
    print(f"out of reach{len(no):3d} elements, median RMS {mn:5.1f} %")
    print("merged into", path)


if __name__ == "__main__":
    main()
