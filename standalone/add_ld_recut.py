#!/usr/bin/env python3
"""
Give the re-cut candidate arms a dispersion, so the phonon panel can draw them.

build_library computes `ld` for the shipped arms only - the top-level record
and, through angular/export_ug.py, the `ug` one.  The candidates were merged
in later from cluster runs and never had one, which is why they are missing
from the phonon plot while they appear everywhere else.

This does for `rc` exactly what build_library does for the top-level arm: the
same Setyawan-Curtarolo path, the same NQ_PATH points per panel, the same
`len` from the reciprocal metric, so the two are drawn on one x mapping and
stay in register.  Frequencies are stored in cm^-1 and rounded to 0.1, as
there.

`rc_ug` is NOT done here.  It needs the angular latdyn, and the two versions
of that module must not meet in one process - see add_ug_overlay.  It goes
through the same door the shipped angular arm uses.

    python add_ld_recut.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import latdyn as L                                   # noqa: E402
import refdata                                       # noqa: E402
from build_library import sc_segments, recip, NQ_PATH, NQ_MESH, T0, DT, NT  # noqa: E402

assert "lam2" not in open(L.__file__, encoding="utf-8").read(), \
    "the ANGULAR latdyn is on the path; this script is for the pair arm only"


def build(el, e, rec):
    cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    pot = L.Potential.from_record(rec)
    Phi = L.force_constants(cry, pot)
    B_ = recip(e["struct"], e.get("c_over_a") or 1.0)

    def sample(ka, kb, npts):
        ks = ka[None, :] + np.linspace(0, 1, npts)[:, None] * (kb - ka)[None, :]
        br = np.array([L.frequencies(cry, pot, q, Phi) * L.CM1 for q in ks])
        return [[round(float(v), 1) for v in br[:, b]]
                for b in range(br.shape[1])]

    std = [{"a": la, "b": lb, "n": NQ_PATH,
            "len": float(np.linalg.norm((kb - ka) @ B_)),
            "branches": sample(ka, kb, NQ_PATH)}
           for (la, ka, lb, kb) in sc_segments(e["struct"])]
    nat = len(cry.frac)
    f = L.spectrum(cry, pot, nq=NQ_MESH)
    th = [dict({"T": T0 + n * DT}, **{k: L.thermo(f, T0 + n * DT, nat)[k]
                                      for k in ("zpe", "F", "S", "Cv")})
          for n in range(NT)]
    return {"std": std, "thermo": th, "maxfreq": float(f.max() * L.CM1)}


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    n = 0
    print(f"{'el':4s}{'max cm-1':>10s}{'most neg.':>12s}   path")
    print("-" * 46)
    for el in sorted(lib):
        v = lib[el]
        rec = v.get("rc")
        if not isinstance(rec, dict) or "m" not in rec:
            continue
        e = refdata.ELEMENTS[el]
        ld = build(el, e, rec)
        rec["ld"] = ld
        lo = min(min(min(b) for b in p["branches"]) for p in ld["std"])
        n += 1
        print(f"{el:4s}{ld['maxfreq']:10.1f}{lo:12.1f}   "
              + "-".join([ld["std"][0]["a"]] + [p["b"] for p in ld["std"]]))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print()
    print(f"{n} candidate arms now carry a dispersion -> {path}")


if __name__ == "__main__":
    main()
