"""ug against tap, on every element with a measured dispersion.

The question is whether copper is a special case or whether the angular
hard-cut arm is generally the one that reproduces dispersion.  Both scored
the same way - polarisation, extended-zone frame, the labelled reading.

Run twice: --ang for the angular records (ug, tap_ug) with the angular
latdyn on the path, plain for the isotropic ones (root, tap).
"""
import json, os, sys
import numpy as np
STAND = r"C:\Users\Admin\Desktop\ugurpotential_screen\standalone"
ANG = r"C:\Users\Admin\Desktop\ugurpotential_screen\angular"
A = "--ang" in sys.argv
sys.argv = [sys.argv[0]] + (["--ug"] if A else [])
sys.path.insert(0, STAND)
if A: sys.path.insert(0, ANG)
import latdyn as L, refdata, curve_mae as CM
from build_library import sc_segments
angfc = __import__("angfc") if A else None
KEYS = (("ug", "ug"), ("tap_ug", "tap_ug")) if A else ((None, "kok"), ("tap", "tap"))

def one(el, v, rec):
    ec = v["exp_curve"]
    cry = L.Crystal(v["struct"], v["a0"], v.get("c_over_a"),
                    mass=refdata.MASSES[el])
    pot = L.Potential.from_record(rec)
    Phi = angfc.force_constants(cry, pot) if A else None
    segs = {(a, b): (np.asarray(ka, float), np.asarray(kb, float))
            for a, ka, b, kb in sc_segments(v["struct"])}
    R = L.reciprocal(cry); hcp = v["struct"] == "hcp"
    err, e = [], []
    for sk, pts in ec["segs"].items():
        a, b = sk.split("|")
        if (a, b) not in segs: continue
        ka, kb = segs[(a, b)]
        for row in pts:
            lab, nu = row[3], row[1]
            q = ka + row[0]*(kb-ka); nrm = None
            if not hcp:
                line = CM.ext_line(v["struct"], el, lab[lab.find("["):], (a,b))
                if line is not None:
                    u,zlo,zhi,fl = line; t = 1.0-row[0] if fl else row[0]
                    q = (zlo+t*(zhi-zlo))*np.asarray(u,float); nc = q@R
                    if np.linalg.norm(nc) > 1e-9: nrm = nc/np.linalg.norm(nc)
            f, vec = L.modes_many(cry, pot, q.reshape(1,3), Phi)
            if hcp:
                val, how = CM.pick_hcp(f[0], vec[0], (q@R), lab)
                got = (float(min(val, key=lambda x: abs(x-nu))) if how=="cift"
                       else float(val) if val is not None
                       else float(min(f[0], key=lambda x: abs(x-nu))))
            elif nrm is not None:
                got = CM.pick_polar(f[0], vec[0], nrm, lab)
            else:
                got = CM.pick(f[0], lab, (a,b))
            err.append(got-nu); e.append(nu)
    e = np.array(e); d = np.array(err)
    return 100*np.abs(d).mean()/e.max(), 100*d.mean()/e.max(), len(e)

lib = json.load(open(os.path.join(STAND, "library.json")))
els = sorted(el for el, v in lib.items()
             if isinstance(v, dict) and "exp_curve" in v)
print("%-4s%6s" % ("el", "n") + "".join("%18s" % n for _, n in KEYS))
for el in els:
    v = lib[el]; out = "%-4s" % el; n0 = 0
    for k, name in KEYS:
        rec = v if k is None else v.get(k)
        if not rec or "m" not in rec: out += "%18s" % "-"; continue
        try:
            mae, sg, n0 = one(el, v, rec)
            out += "%12.1f%%%+5.1f" % (mae, sg)
        except Exception as ex:
            out += "%18s" % ("hata")
    print("%-4s%6d%s" % (el, n0, out[4:]))
