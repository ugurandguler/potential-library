#!/usr/bin/env python3
"""
The Grueneisen parameter: the phonon side of the thermal expansion.

The library's thermal expansion is 32 per cent low and nine records contract
on heating.  Expansion is not a separate property from the spectrum: in the
quasi-harmonic picture it is the spectrum's volume dependence,

    gamma_i = -d ln omega_i / d ln V ,     gamma = sum_i c_i gamma_i / sum_i c_i

with c_i the mode heat capacity, and then

    alpha_V = gamma C_V / (B_T V) .

So a record can be asked where its expansion comes from: whether the modes
stiffen too little under compression (gamma too small) or the rest of the
chain is at fault.  Both sides of that identity are computable here and the
measured value is in the library, which makes it a closed test rather than a
statement of belief.

    gamma      from frequencies on a Monkhorst-Pack mesh at +/-1 % of the
               lattice constant, weighted by each mode's own heat capacity at
               298.15 K - analytic, no molecular dynamics
    alpha_V    predicted from that gamma with the record's own C_V and bulk
               modulus, against the measured alpha_V
    gamma_exp  from the MEASURED alpha_V, B and C_p at 298 K

    python gruneisen.py                # every element, the switched arm
    python gruneisen.py Cu W --set tap
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

KB = 8.617333262e-5            # eV/K
THZ_TO_EV = 4.135667696e-3     # h * 1 THz in eV
NA = 6.02214076e23
T0 = 298.15
MESH = 8


WINDOWS = (0.005, 0.01, 0.02, 0.03)


def mode_gamma(el, rec, mass, n=MESH, T=T0, d=0.01):
    """heat-capacity-weighted mode Grueneisen parameter at T

    d is the half-width of the volume window the derivative is taken over.
    It is an argument because the answer depends on it: for the body-centred
    refractories the frequencies are not smooth in volume on the per-cent
    scale - the same wiggle the equation-of-state scan shows - and gamma
    changes sign between a one and a three per cent window.  main() therefore
    reports the whole set rather than one number."""
    e = refdata.ELEMENTS[el]
    pot = L.Potential.from_record(rec)
    q = L.mesh(n)
    f = {}
    for s in (1 - d, 1.0, 1 + d):
        cry = L.Crystal(e["struct"], s * e["a0"], e.get("c_over_a"), mass=mass)
        f[s] = np.asarray(L.frequencies_many(cry, pot, q))      # THz, (nq, nbranch)
    w0 = f[1.0]
    dlnV = 3 * np.log((1 + d) / (1 - d))
    with np.errstate(divide="ignore", invalid="ignore"):
        g = -(np.log(np.abs(f[1 + d])) - np.log(np.abs(f[1 - d]))) / dlnV
    #  the acoustic modes at Gamma are zero and their gamma is undefined; so
    #  is any imaginary branch, which is not a mode of a stable crystal
    ok = (w0 > 0.05) & np.isfinite(g) & (f[1 + d] > 0) & (f[1 - d] > 0)
    x = THZ_TO_EV * w0[ok] / (KB * T)
    c = KB * x ** 2 * np.exp(x) / (np.exp(x) - 1.0) ** 2          # per mode, eV/K
    return float(np.sum(c * g[ok]) / np.sum(c)), int(ok.sum()), int(ok.size)


def molar_volume(el):
    e = refdata.ELEMENTS[el]
    a = e["a0"]
    if e["struct"] == "fcc":
        v = a ** 3 / 4
    elif e["struct"] == "bcc":
        v = a ** 3 / 2
    else:
        v = np.sqrt(3) / 4 * a ** 3 * e.get("c_over_a", 1.633)
    return v * 1e-30 * NA          # m^3/mol


def main():
    argv = sys.argv[1:]
    arm = "tap"
    if "--set" in argv:
        i = argv.index("--set")
        arm = argv[i + 1]
        del argv[i:i + 2]
    lib = json.load(io.open(os.path.join(HERE, "library.json"), encoding="utf-8"))
    els = [a for a in argv if not a.startswith("--")] or sorted(lib)
    out = {}
    print(f"{'el':4}{'gamma':>8}{'gamma exp':>11}{'alpha model':>13}{'alpha meas':>12}"
          f"{'alpha from g':>14}")
    for el in els:
        rec = lib[el].get(arm)
        if not isinstance(rec, dict):
            continue
        exp = (rec.get("expansion") or {})
        aL_meas = exp.get("alpha_exp_1e6")
        aL_model = exp.get("alpha_1e6")
        B = lib[el]["B"][1] if isinstance(lib[el].get("B"), list) else lib[el].get("B")
        cv = None
        for r in sorted((lib[el].get("ld") or {}).get("thermo", []), key=lambda r: r["T"]):
            if r["T"] >= T0:
                cv = r["Cv"]                      # J/(mol K)
                break
        if not (aL_meas and B and cv):
            continue
        gw = {d: mode_gamma(el, rec, lib[el]["mass"], d=d)[0] for d in WINDOWS}
        g, nok, ntot = mode_gamma(el, rec, lib[el]["mass"])
        Vm = molar_volume(el)
        #  alpha_V = gamma C_V / (B V), with B in Pa and V in m^3/mol
        a_pred = g * cv / (B * 1e9 * Vm)
        g_exp = (3 * aL_meas * 1e-6) * (B * 1e9) * Vm / cv
        out[el] = dict(gamma=g, gamma_exp=g_exp, modes_used=nok, modes=ntot,
                       gamma_windows={str(k): v for k, v in gw.items()},
                       window_spread=max(gw.values()) - min(gw.values()),
                       sign_changes=min(gw.values()) < 0 < max(gw.values()),
                       alpha_V_pred_1e6=1e6 * a_pred, alpha_V_meas_1e6=3 * aL_meas,
                       alpha_V_model_1e6=3 * aL_model if aL_model else None,
                       Cv_298=cv, B_GPa=B, Vm_m3=Vm)
        print(f"{el:4}{g:8.2f}{g_exp:11.2f}"
              f"{(3 * aL_model if aL_model else float('nan')):13.1f}{3 * aL_meas:12.1f}"
              f"{1e6 * a_pred:14.1f}")
    out["_note"] = ("gamma is the heat-capacity-weighted mode Grueneisen parameter at 298.15 K "
                    f"on an {MESH}^3 mesh, from frequencies at +/-1 % of the lattice constant; "
                    "gamma_exp is alpha_V B V / C_V with the measured alpha_V and B and the "
                    "record's own C_V; alpha_V_pred is the quasi-harmonic prediction from "
                    "gamma, alpha_V_model the molecular-dynamics value, both in 10^-6/K.")
    out["_source"] = f"library.json {arm} records + latdyn frequencies; arm = {arm}"
    json.dump(out, io.open(os.path.join(HERE, f"gruneisen_{arm}.json"), "w", encoding="utf-8"),
              indent=1, sort_keys=True)
    print(f"\n{len(out) - 2} elements -> gruneisen_{arm}.json")


if __name__ == "__main__":
    main()
