#!/usr/bin/env python3
"""
The measured-dispersion error of PUBLISHED potentials, by the library's rules.

lammps/phonon_published.py takes each published EAM and MEAM file through
LAMMPS and writes its force constants in latdyn's own keying; this hands them
to curve_mae, so every point is read with the same branch labels, the same
eigenvector test for L and T, the same symmetry-point de-duplication and the
same normalisation by the highest measured frequency as the library's own
records.  Only the force constants differ.

    python curve_mae_published.py --check    our Cu|tap through LAMMPS, against
                                             latdyn's analytic frequencies
    python curve_mae_published.py            every published potential ->
                                             curve_mae_published.json
    python curve_mae_published.py --aitken   the MEAM files of Aitken et al.
                                             (2021) -> curve_mae_aitken.json
"""
import json
import os
import statistics as st
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LMPDIR = os.path.join(HERE, "..", "lammps")
sys.path.insert(0, HERE)
import curve_mae as CM                      # noqa: E402
import latdyn as L                          # noqa: E402

OUT = os.path.join(HERE, "curve_mae_published.json")


def phi_of(r):
    return {(i, j, tuple(R)): np.array(b) for i, j, R, b in r["phi"]}


def check():
    """the library's switched copper record two ways: analytic (latdyn) and
    finite differences through LAMMPS; the frequencies and the score must agree"""
    r = json.load(open(os.path.join(LMPDIR, "phonon_published_check.json")))
    lib = json.load(open(os.path.join(HERE, "library.json")))
    v = lib["Cu"]
    cry = L.Crystal("fcc", r["a"], None, mass=63.546)
    pot = L.Potential.from_record(v["tap"])
    qs = L.mesh(6)
    fa = L.frequencies_many(cry, pot, qs)
    fl = L.frequencies_many(cry, None, qs, phi_of(r))
    print(f"Cu|tap on a 6^3 mesh: max |analytic - LAMMPS| = {np.abs(fa - fl).max():.4f} THz "
          f"(top {fa.max():.3f} THz)")
    CM.EXTERNAL = lambda el: [("LAMMPS", phi_of(r))] if el == "Cu" else []
    got = CM.main(("Cu",))
    CM.EXTERNAL = None
    print("curve_mae through LAMMPS:", got[0].get("LAMMPS"))
    print("curve_mae analytic (MAU):  rerun `python curve_mae.py Cu` and compare the MAU column")


def main(src="phonon_published.json", dst=OUT):
    raw = json.load(open(os.path.join(LMPDIR, src)))
    ok = {k: r for k, r in raw.items() if "phi" in r}
    by_el = {}
    for k, r in ok.items():
        by_el.setdefault(r["el"], []).append((k, r))
    CM.EXTERNAL = lambda el: [(k, phi_of(r)) for k, r in by_el.get(el, [])]
    rows = CM.main(tuple(sorted(by_el)))
    out = {}
    for row in rows:
        for k, r in by_el.get(row["el"], []):
            if k in row:
                out[k] = dict(el=row["el"], family=r["family"], style=r["style"],
                              mae_THz=row[k][0], pct=row[k][1], n=row["n"],
                              hcp=row.get("hcp", False), near=row.get("near", False))
    failed = {k: r.get("error", "")[-200:] for k, r in raw.items() if "phi" not in r}
    json.dump(dict(scores=out, failed=failed), open(dst, "w"), indent=1)
    print(f"{len(out)} published potentials scored over {len(by_el)} elements; "
          f"{len(failed)} failed in LAMMPS")
    for fam in ("EAM", "MEAM"):
        v = [s["pct"] for s in out.values() if s["family"] == fam]
        if v:
            print(f"  {fam}: {len(v)} potentials, mean {st.mean(v):.1f}%, median {st.median(v):.1f}%")


if __name__ == "__main__":
    if "--check" in sys.argv:
        check()
    elif "--aitken" in sys.argv:
        #  Aitken et al.'s MEAM files, scored the same way, kept apart
        main("phonon_aitken.json", os.path.join(HERE, "curve_mae_aitken.json"))
    else:
        main()
