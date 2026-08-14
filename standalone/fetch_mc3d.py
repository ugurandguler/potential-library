#!/usr/bin/env python3
"""
Pull phonon dispersions from the Materials Cloud MC3D database.

MC3D itself is a structural database - its paper does not mention phonons at
all - but the same web application serves a separate contribution,
`supercon_phonon-vis`, holding full dispersions from an electron-phonon
superconductivity study (PBEsol, Quantum ESPRESSO).  That is where the useful
data is, and it covers the metals Materials Project does not:

    MP gives us      Be Ca Cu Ir K Li Mg Pt Sr Ti V W Zr
    MC3D adds        Ag Au Cr Mo Nb Ni Pb Pd Rh Ta

The additions matter because Nb, Cr, Mo and Ta are exactly the bcc metals whose
anisotropy the potential cannot reach, and until now nothing independent
described their phonons.

The API is not documented anywhere; these routes were read out of the site's
own JavaScript bundle:

    /mc3d/dataset-index/<id>                       which contributions exist
    /mc3d/pbesol-v1/supercon_phonon-vis/<id>       the dispersion itself

Their q-points are in the same primitive fractional convention we use - H is
(1/2, -1/2, 1/2), P is (1/4, 1/4, 1/4), N is (0, 0, 1/2), all identical to
build_library.SC_POINTS - so our potential can be evaluated at their q-points
directly, with no conversion.

One caveat worth carrying: this data was produced for electron-phonon coupling,
so the q-meshes and smearing were chosen for that, not for a phonon benchmark.
Treat it as a second DFT opinion, not as a measurement.

    python fetch_mc3d.py            # everything we have an id for
    python fetch_mc3d.py Nb Mo

Writes mc3d_phonon.json; add_mc3d_overlay.py merges it into library.json.
"""
import json
import os
import subprocess
import sys

API = "https://mcxd-api.materialscloud.org/mc3d"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "mc3d_phonon.json")

#  Found by scanning /mc3d/pbesol-v2/overview for single-element formulas and
#  keeping the entry whose space group matches the phase we fit: 225 for fcc,
#  229 for bcc, 194 for hcp.  Elements absent here either have no MC3D entry in
#  our phase or carry no phonon contribution (Al, Ba, Ca, Co, Fe, Li, Mg, Na,
#  Sr) - Al is the notable gap, and the one place the measured
#  Landolt-Boernstein values matter most.
IDS = {
    "Ag": "mc3d-75785", "Au": "mc3d-66164", "Be": "mc3d-40378",
    "Cd": "mc3d-73497", "Cr": "mc3d-49144", "Cu": "mc3d-24185",
    "Ir": "mc3d-74998", "K": "mc3d-55121", "Mo": "mc3d-10871",
    "Nb": "mc3d-10833", "Ni": "mc3d-4988", "Pb": "mc3d-48868",
    "Pd": "mc3d-36406", "Pt": "mc3d-31294", "Rh": "mc3d-4421",
    "Ta": "mc3d-47329", "Ti": "mc3d-21616", "V": "mc3d-45322",
    "W": "mc3d-38286", "Zn": "mc3d-50813", "Zr": "mc3d-53433",
    #  added 2026-08-03 with the ten new elements, found the same way: the
    #  pbesol-v2 overview scanned for a single-element formula whose space group
    #  is the phase we fit.  Ytterbium and lutetium have no MC3D entry at all,
    #  and rhenium has entries only at space groups 225, 213 and 139 - none of
    #  them the hcp phase - so those three get no MC3D column.
    "Rb": "mc3d-71910", "Cs": "mc3d-71072", "Sc": "mc3d-12651",
    "Y": "mc3d-43809", "Hf": "mc3d-54508", "Ru": "mc3d-57338",
    "Tl": "mc3d-37604",
}
THIN = 3            # keep every third q-point; 350 -> 117, plenty for a curve


def get(url):
    p = subprocess.run(["curl", "-s", "-L", "--max-time", "60", url],
                       capture_output=True, text=True, timeout=90)
    return json.loads(p.stdout)


def main():
    els = sys.argv[1:] or sorted(IDS)
    store = json.load(open(OUT)) if os.path.exists(OUT) else {}
    print(f"{'el':4s}{'mc3d id':>14s}{'q':>6s}{'branches':>10s}"
          f"{'max cm-1':>10s}   path")
    print("-" * 74)
    for el in els:
        mid = IDS.get(el)
        if not mid:
            print(f"{el:4s}  no MC3D id")
            continue
        try:
            d = get(f"{API}/pbesol-v1/supercon_phonon-vis/{mid}")
        except Exception as exc:                          # noqa: BLE001
            print(f"{el:4s}  fetch failed: {type(exc).__name__}")
            continue
        if "eigenvalues" not in d:
            print(f"{el:4s}  no dispersion in the record")
            continue
        q = d["qpoints"][::THIN]
        f = d["eigenvalues"][::THIN]
        #  labels come as [index, name]; a name like "H|P" marks the
        #  discontinuity, exactly as in the Materials Project records
        marks, breaks = [], []
        for idx, name in d["highsym_qpts"]:
            j = idx // THIN
            if "|" in name:
                a, b = name.split("|")
                marks.append([j, a])
                breaks.append(j)
                marks.append([j + 1, b])
            else:
                marks.append([j, name])
        store[el] = {
            "id": mid, "method": "pbesol, Quantum ESPRESSO (supercon_phonon)",
            "q": [[round(float(x), 6) for x in v] for v in q],
            "f": [[round(float(x), 1) for x in row] for row in f],
            "marks": marks, "breaks": breaks,
            "source": ("Materials Cloud MC3D, supercon_phonon-vis "
                       "contribution (PBEsol, Quantum ESPRESSO)"),
        }
        top = max(max(r) for r in f)
        print(f"{el:4s}{mid:>14s}{len(q):6d}{len(f[0]):10d}{top:10.1f}   "
              f"{'-'.join(n for _, n in marks)}")

    tmp = OUT + ".tmp"
    json.dump(store, open(tmp, "w"))
    os.replace(tmp, OUT)
    print(f"\n{len(store)} elements in {OUT}")


if __name__ == "__main__":
    main()
