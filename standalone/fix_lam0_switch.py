#!/usr/bin/env python3
"""
v1.4.7: the elastic constants of the switched angle-free records, recomputed
with the per-leg three-body switch.

angular/latdyn.py, which add_candidates.py and the UG searches use, built the
three-body force constants and gradient from phi3'(r1 + r2) alone, without the
S(r1) S(r2) of the switch that its own energy() applies (fixed 2026-11-04 by
porting standalone/latdyn.Potential.leg3).  The gradient enters the internal
relaxation of the hcp cell, so every switched record with lam2 = lam4 = 0
carried hcp elastic constants from an unswitched three-body term: Co, Mg and
Re of the re-cut set, and Be and Zr of tap_ug (V too, but it is cubic, and a
one-atom cell has no internal relaxation).  Found by an independent elastic test,
2026-11-01; the true values are what that test, standalone/latdyn and a direct
finite-strain probe all give.

The parameters are NOT refitted.  The re-cut fits were made with the old code,
so the parameters were chosen to fit its numbers, and the honest record of the
shipped parameters is their true constants and the error those carry against
experiment - which goes up.  Rewritten per record, only where a value moves:

  Cij      from standalone/latdyn.elastic (= the fixed fork, bit for bit)
  rms      the relative RMS against refdata over the measured constants,
           which is how every stored rms is defined (checked on all of them)
  score    = rms / 100 where the record carries the search's score
  mech, mech_planes   as add_mech_recut.py + add_plane_d.py + refresh_mech.py
           build them, for the arms that carry them (rc, rc_ug)

`dyn`/`stable` come from add_dynstab_arms.py and tap_ug's `min_cm1` from
angular/dynstab_ug.py: rerun those, with the fixed fork, rather than here.

    python fix_lam0_switch.py [--dry]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import latdyn as L                                   # noqa: E402
import mechanics                                     # noqa: E402
import refdata                                       # noqa: E402
from add_mech_recut import tensor                    # noqa: E402

assert "leg3" in open(L.__file__, encoding="utf-8").read(), L.__file__
ARMS = ("rc", "rc_ug", "tap_ug")
KEYS = {"hcp": ("C11", "C12", "C13", "C33", "C44"), "fcc": ("C11", "C12", "C44"),
        "bcc": ("C11", "C12", "C44")}
IDX = {"C11": (0, 0), "C12": (0, 1), "C13": (0, 2), "C33": (2, 2), "C44": (3, 3)}
TOL = 1e-3                    # GPa; the two codes differentiate to ~1e-6, the bug moves tens


def rel_rms(cij, el):
    ce = refdata.ELEMENTS[el]["Cij"]
    ks = [k for k in cij if ce.get(k)]
    return 100.0 * float(np.sqrt(np.mean([((cij[k] - ce[k]) / ce[k]) ** 2 for k in ks])))


def mech_of(el, cij):
    """the mech block and its four plane sections, built the way the chain does"""
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"), mass=refdata.MASSES[el])
    nat = len(cry.frac)
    C = tensor(e["struct"], cij)
    m = mechanics.analyse(C, mass_amu=refdata.MASSES[el] * nat, volume_A3=cry.vol, natoms=nat)
    pl = {p: mechanics.plane_curves(C, p) for p in ("xy", "xz", "yz", "d")}
    md = mechanics.analyse(C, mass_amu=refdata.MASSES[el] * nat, volume_A3=cry.vol,
                           natoms=nat, nu=400)
    for k in ("debye", "debye_iso", "v_m", "v_m_iso"):
        if k in md:
            m[k] = md[k]
    return m, pl


def same(a, b):
    return json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)


def main(dry):
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    changed, checked = [], 0
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or el not in refdata.ELEMENTS:
            continue
        e = refdata.ELEMENTS[el]
        for arm in ARMS:
            rec = v.get(arm)
            if not isinstance(rec, dict) or "m" not in rec or rec.get("taper") is None:
                continue
            if rec.get("lam2", 0.0) or rec.get("lam4", 0.0):
                continue                         # angfc path, switched already
            cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                            mass=refdata.MASSES[el])
            C, _ = L.elastic(cry, L.Potential.from_record(rec))
            new = {k: float(C[IDX[k]]) for k in KEYS[e["struct"]] if k in rec["Cij"]}
            d = max(abs(new[k] - rec["Cij"][k]) for k in new)
            checked += 1
            if d < TOL:
                #  the control: an untouched record must come back as stored,
                #  mech included, or this script is not the chain it imitates
                if isinstance(rec.get("mech"), dict):
                    m, pl = mech_of(el, rec["Cij"])
                    if not same(m, rec["mech"]) or not same(pl, rec.get("mech_planes")):
                        raise SystemExit(f"{el}/{arm}: mech not reproduced from its own Cij")
                continue
            if e["struct"] != "hcp":
                raise SystemExit(f"{el}/{arm}: a cubic tensor moved by {d:.3g} GPa - "
                                 "a one-atom cell has no internal relaxation, look first")
            old = dict(rec["Cij"])
            rms_old = rec.get("rms")
            if rms_old is not None and abs(rel_rms(old, el) - rms_old) > 1e-9:
                raise SystemExit(f"{el}/{arm}: stored rms is not the relative RMS of its Cij")
            rec["Cij"] = {**rec["Cij"], **new}
            rec["rms"] = rel_rms(rec["Cij"], el)
            if "score" in rec and rms_old is not None and abs(rec["score"] - rms_old / 100) < 1e-12:
                rec["score"] = rec["rms"] / 100.0
            if isinstance(rec.get("mech"), dict):
                rec["mech"], rec["mech_planes"] = mech_of(el, rec["Cij"])
            changed.append((el, arm))
            print(f"{el:3} {arm:7} " + "  ".join(f"{k} {old[k]:.1f}->{new[k]:.1f}" for k in new)
                  + f"   rms {rms_old:.2f}->{rec['rms']:.2f} %"
                  + (f"   B_H {rec['mech']['B_H']:.1f}" if isinstance(rec.get("mech"), dict) else ""))
    print(f"\n{checked} switched angle-free records checked, {len(changed)} changed: "
          + ", ".join(f"{e}/{a}" for e, a in changed))
    if dry or not changed:
        return
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print("library.json written")


if __name__ == "__main__":
    main("--dry" in sys.argv)
