#!/usr/bin/env python3
"""
What the elastic constants DO when the crystal is heated, against measurement.

The fit is at one temperature: the constants are targets at 0 K in a static
lattice, matched to room-temperature measurements.  Their temperature
derivative is not fitted anywhere, and Landolt-Boernstein III/29a tabulates it
for most of the library (lb29a_tempcoef.json, read from the page images).

Two derivatives are not the same quantity, and the difference is the whole
care of this script:

  Tc(V)   what the finite-temperature sweep measures.  elastic_T.py holds the
          cell at its 0 K volume, so its C_ij(T) is the constant-volume
          derivative: pure anharmonicity of the vibrations.

  Tc(P)   what the experiment measures.  A heated crystal also expands, and
          the elastic constants fall with volume, which is most of the effect
          in a simple metal.

So the expansion term is added back here,

    Tc(P) = Tc(V) + (dlnC/dlnV) * 3 alpha_L,

with dlnC/dlnV computed analytically at +/-1 % of the lattice constant and
alpha_L the MEASURED linear expansion coefficient - not the model's own, which
is 32 % low and would fold that error into this comparison.

    python tempcoef_compare.py            # writes tempcoef_compare.json
    python tempcoef_compare.py --print
"""
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import latdyn as L      # noqa: E402
import refdata          # noqa: E402

CUBIC = (("Tc11", "C11", (0, 0)), ("Tc44", "C44", (3, 3)), ("Tc12", "C12", (0, 1)))
HEX = CUBIC + (("Tc33", "C33", (2, 2)), ("Tc13", "C13", (0, 2)))
ARM = "tap"


def dln_dlnV(el, rec, keys):
    """d ln C_ij / d ln V at the fitted volume, analytic, +/-1 % in a"""
    e = refdata.ELEMENTS[el]
    C = {}
    for s in (0.99, 1.0, 1.01):
        cry = L.Crystal(e["struct"], s * e["a0"], e.get("c_over_a"), mass=1.0)
        C[s] = L.elastic(cry, L.Potential.from_record(rec))[0]
    dlnV = 3 * np.log(1.01 / 0.99)
    return {k: float(((C[1.01][i, j] - C[0.99][i, j]) / C[1.0][i, j]) / dlnV)
            for k, _, (i, j) in keys}


def main():
    lib = json.load(io.open(os.path.join(HERE, "library.json"), encoding="utf-8"))
    ref = json.load(io.open(os.path.join(HERE, "lb29a_tempcoef.json"), encoding="utf-8"))
    out = {}
    for el, r in ref.items():
        if el.startswith("_") or el not in lib:
            continue
        rec = lib[el].get(ARM)
        pts = ((lib[el].get("elasticT") or {}).get(ARM) or {})
        if not isinstance(rec, dict) or not pts:
            continue
        keys = HEX if lib[el]["struct"] == "hcp" else CUBIC
        #  the measured coefficients are quoted over ranges that end at or
        #  below room temperature, so the model slope is taken over the same
        #  span: every grid point up to 300 K that is still Born stable
        use = [p for p in pts["pts"] if p["T"] <= 300.0 and p["born_ok"]]
        if len(use) < 3:
            out[el] = {"skipped": f"only {len(use)} stable points at or below 300 K"}
            continue
        T = np.array([p["T"] for p in use])
        alpha = ((rec.get("expansion") or {}).get("alpha_exp_1e6"))
        dln = dln_dlnV(el, rec, keys)
        row = {"T_points": len(use), "T_max": float(T.max()),
               "alpha_exp_1e6": alpha, "system": r.get("system")}
        for name, key, _ in keys:
            if name not in r or name in r.get("uncertain", []):
                continue
            C = np.array([p[key] for p in use], float)
            slope = float(np.polyfit(T, C, 1)[0])
            tcV = 1e4 * slope / float(C[0])
            row[name] = dict(
                model_V=tcV,
                model_P=(tcV + 1e4 * dln[name] * 3 * alpha * 1e-6) if alpha else None,
                dlnC_dlnV=dln[name], experiment=r[name])
        out[el] = row
    out["_note"] = ("Tc in 10^-4/K. model_V is the constant-volume derivative the "
                    "finite-temperature sweep measures; model_P adds (dlnC/dlnV) * 3 alpha_L "
                    "with the MEASURED expansion coefficient, which is what the experiment "
                    "reports. Coefficients the volume prints as estimates are skipped.")
    out["_source"] = ("library.json elasticT." + ARM + " (lammps/elastic_T.py) + latdyn.elastic "
                      "at +/-1 % of a + lb29a_tempcoef.json")
    json.dump(out, io.open(os.path.join(HERE, "tempcoef_compare.json"), "w", encoding="utf-8"),
              indent=1, sort_keys=True)
    els = [k for k in out if not k.startswith("_") and "skipped" not in out[k]]
    print(f"{len(els)} elements compared")
    if "--print" in sys.argv:
        print(f"{'el':4}{'quantity':9}{'model V':>9}{'model P':>9}{'experiment':>12}")
        for el in sorted(els):
            for name in ("Tc11", "Tc33", "Tc44", "Tc12", "Tc13"):
                v = out[el].get(name)
                if v:
                    print(f"{el:4}{name:9}{v['model_V']:9.2f}"
                          f"{(v['model_P'] if v['model_P'] is not None else float('nan')):9.2f}"
                          f"{v['experiment']:12.2f}")


if __name__ == "__main__":
    main()
