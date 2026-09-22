#!/usr/bin/env python3
"""
Merge the re-ranked candidates into a COPY of library.json as two more sets,
so the interface can be looked at before anything is decided.

`add_taper_overlay.py` is the template - same fields, same shape - with three
differences that matter:

  * latdyn comes from angular/, because 23 of the 46 candidate records carry
    lam2 and lam4 and standalone/latdyn.py has no angular terms at all.  Built
    with the wrong module they would be evaluated as if the Legendre factor
    were absent, which is quietly wrong rather than loudly wrong.
  * stability is measured on the mesh AND the Setyawan-Curtarolo path, at a
    physical -1 cm^-1 threshold.  `add_taper_overlay` uses `mn > -1e-3` on the
    mesh alone; that is how iridium came to be published as stable with a mode
    at -8.1 cm^-1 half way along K-Gamma.
  * the ground-state verdict is carried over from refit/dyn_candidates.json
    rather than recomputed, so the page cannot disagree with REFIT section 8b.

Only the LAMMPS-free quantities are here.  The four panels that need LAMMPS -
elasticT, bain, the gamma line and thermal expansion - are not, and the record
says so rather than showing an empty chart.

    python add_candidates.py           # writes library.json in THIS tree
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "angular"))

import latdyn as L                                  # noqa: E402
assert "lam2" in open(L.__file__, encoding="utf-8").read(), \
    f"wrong latdyn: {L.__file__} carries no angular terms"

sys.path.append(HERE)
import refdata                                      # noqa: E402
from build_library import sc_segments, NQ_PATH      # noqa: E402

CAND = os.path.join(ROOT, "refit")
SETS = {"mau": "rc", "ug": "rc_ug"}
IMAG_TOL_CM1 = -1.0


def worst_frequency(el, p):
    """lowest frequency over the mesh and the standard path, in cm^-1"""
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    pot = L.Potential.from_record(p)
    Phi = L.force_constants(cry, pot)
    qs = []
    for (_, ka, __, kb) in sc_segments(e["struct"]):
        qs.extend(ka + np.linspace(0, 1, NQ_PATH)[:, None] * (kb - ka))
    f = np.concatenate([L.spectrum(cry, pot, nq=n).ravel() for n in (8, 9)]
                       + [L.frequencies_many(cry, pot, np.array(qs),
                                             Phi).ravel()]) * L.CM1
    return float(f.min())


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    ground = {}
    #  THE GROUND-STATE VERDICT IS OPTIONAL AND NOW SAYS SO WHEN IT IS ABSENT.
    #  `dyn_candidates.json` is in NEITHER tree - it was never carried over
    #  from the screening work - so this used to skip in silence and every
    #  candidate came out with an empty `ground` column.  The verdicts already
    #  in library.json are untouched; only a REBUILD loses them, and a rebuild
    #  that drops a column without saying so is the failure this project keeps
    #  finding.  Checked 2026-09-08.
    gp = os.path.join(CAND, "dyn_candidates.json")
    if os.path.exists(gp):
        for k, v in json.load(open(gp)).items():
            el, arm = k.split("|")
            ground[(el, arm.lower())] = v["verdict"]
    else:
        print("WARNING: %s missing - the 'ground' column will stay EMPTY. The "
              "verdicts already in the records are unaffected; only this rebuild "
              "cannot produce them." % os.path.relpath(gp, ROOT), flush=True)

    print(f"{'el':4s}{'set':7s}{'rms':>8s}{'shipped':>9s}{'min cm-1':>10s}"
          f"{'ground':>9s}  stable")
    print("-" * 56)
    n = 0
    for f in sorted(glob.glob(os.path.join(CAND, "fit_*_force_*.json"))):
        el, _, arm = os.path.basename(f)[4:-5].split("_")
        if arm not in SETS or el not in lib or el not in refdata.ELEMENTS:
            continue
        r = json.load(open(f))["params"]
        e = refdata.ELEMENTS[el]
        c = r["Cij"]
        rec = {"rms": 100.0 * r["score"], "taper": r.get("taper"),
               "m": r["m"], "gamma": r["gamma"], "C": r["C"], "s3": r["s3"],
               "D": r["D"], "alpha": r["alpha"], "r0": r["r0"],
               "alpha3": r["alpha3"], "rcut2": r["rcut2"],
               "rcut3": r["rcut3"], "dnn": r["dnn"], "struct": e["struct"],
               "Ecoh": r.get("Ecoh"), "B": r.get("B"), "P": r.get("P"),
               "lam2": r.get("lam2", 0.0), "lam4": r.get("lam4", 0.0),
               "at_bound": r.get("at_bound", []),
               "Cij": {k: c[k] for k in c if c.get(k) is not None},
               #  said explicitly rather than left as a missing key, so the
               #  viewer can state it instead of drawing an empty panel
               "lammps_pending": ["elasticT", "bain", "gamma", "expansion",
                                  "stacking", "surface", "md_screen",
                                  "jiggle"]}
        if e["struct"] in ("fcc", "bcc"):
            cp = 0.5 * (c["C11"] - c["C12"])
            ce = e["Cij"]
            rec["R"] = c["C44"] / cp if cp > 0 else None
            rec["R_exp"] = ce["C44"] / (0.5 * (ce["C11"] - ce["C12"]))
        mn = worst_frequency(el, r)
        rec["min_cm1"] = round(mn, 3)
        rec["stable"] = bool(mn > IMAG_TOL_CM1)
        g = ground.get((el, arm))
        if g:
            rec["ground"] = {"want": g["want"], "lowest": g["winner_raw"],
                             "rel": g["gap_meV"],
                             "ok": g["winner_raw"] == g["want"]}
        lib[el][SETS[arm]] = rec
        n += 1
        gtxt = ("-" if not g else
                ("ok" if rec["ground"]["ok"] else rec["ground"]["lowest"]))
        print(f"{el:4s}{SETS[arm]:7s}{rec['rms']:8.2f}"
              f"{(lib[el].get('tap') or {}).get('rms', float('nan')):9.2f}"
              f"{mn:10.2f}{gtxt:>9s}  "
              f"{'yes' if rec['stable'] else 'NO'}")

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"\n{n} candidate records merged into {path}")
    print("this tree only - the published library.json is untouched")


if __name__ == "__main__":
    main()
