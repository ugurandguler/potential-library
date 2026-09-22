#!/usr/bin/env python3
"""Score one element's finite-temperature dispersion where the data is.

mdcurve_mae.py does this locally, and would be the only copy if the cluster
link could carry an eigenvector file.  It cannot: anything over about 8 kB
arrives empty, and phana's output for niobium alone is 120 kB.  So the
picking runs next to the file and only the answer travels - one line per
element, plus the picked frequency at each point, which is small enough to
come back and is what a plot would need later.

Everything it needs is in mdq_<El>.json, written by `mdcurve_mae.py emit`:
the measured frequency, the branch label, and the wavevector direction each
label is a polarisation with respect to.  The branch a label names is decided
by its eigenvector, exactly as the 0 K column decides it.

    python3 mdscore_remote.py El ev_El.dat mdq_El.json
"""
import json
import sys

import numpy as np


def pick_polar(f3, vec, n, label):
    """the branch a polarisation label names, decided by the eigenvector

    Copied from curve_mae.pick_polar rather than imported, because the tree
    it lives in is not on the cluster.  If one is ever changed the other must
    be: the whole point of this column is that it is picked the same way.
    """
    w = np.abs(np.asarray(vec).conj().T @ np.asarray(n, float))
    lo = int(np.argmax(w))
    head = label.split("[")[0]
    if head == "L":
        return float(f3[lo])
    tr = sorted(float(f3[k]) for k in range(len(f3)) if k != lo)
    if head == "T1":
        return tr[0]
    if head == "T2":
        return tr[-1]
    return 0.5 * (tr[0] + tr[-1])


#  ---- lifted verbatim from curve_mae.py ----------------------------------
#  Copied, not imported, because the tree curve_mae lives in is not on the
#  cluster.  The whole point of this column is that it is read the same way
#  as the 0 K one, so if one copy changes the other must; the lift is done
#  by script from the original text for exactly that reason.
#  the c axis of the hcp setting latdyn and mkprim.py both use
CHAT = np.array([0.0, 0.0, 1.0])


def hcp_classes(vec, qc):
    """(class index per mode, whether perp and par are distinguishable)

    Splits the six branches by polarisation three ways at once: along q,
    perpendicular to the BASAL PLANE, and parallel to it.  The records say
    in as many words that this is what the papers' perp and par mean, so it
    is not an interpretation of the labels.

    Only the MAGNITUDE of each atom's polarisation along each axis is used,
    never the relative phase of the two atoms, so the answer does not depend
    on whether the phase convention puts the basis offset inside the
    dynamical matrix.  latdyn and phana need not agree about that, and this
    way it never has to be checked that they do.

    Returns None when q is at the zone centre, where the directions
    degenerate and there is nothing to classify.
    """
    n = np.linalg.norm(qc)
    if n < 1e-9:
        return None
    qh = np.asarray(qc, float) / n
    p = CHAT - np.dot(CHAT, qh) * qh
    np_ = np.linalg.norm(p)
    #  q along c: the two transverse directions are degenerate by symmetry
    #  and the paper does not label them separately there either
    split = np_ > 1e-6
    if split:
        p = p / np_
        t = np.cross(qh, p)
    else:
        #  any two directions orthogonal to q will do; they are equivalent
        p = np.array([qh[1], -qh[0], 0.0])
        if np.linalg.norm(p) < 1e-9:
            p = np.array([1.0, 0.0, 0.0])
        p = p / np.linalg.norm(p)
        t = np.cross(qh, p)
    axes = (qh, p, t)
    cls = []
    for m in range(vec.shape[1]):
        e = np.asarray(vec[:, m]).reshape(-1, 3)
        w = [float(np.sum(np.abs(e @ d) ** 2)) for d in axes]
        cls.append(int(np.argmax(w)))
    return cls, split


def pick_hcp(f6, vec, qc, label):
    """the branch an hcp label names, decided by the eigenvectors

    Six branches, and the labels split them three ways at once: L or T with
    respect to q, perp or par with respect to the basal plane, and acoustic
    or optic.  The first two come out of the eigenvector.  The third is then
    the frequency ORDER within a polarisation class, which is what acoustic
    and optic mean once the class is fixed - and is why this does not repeat
    the mistake the cubic rows made, where sorted position was used across
    classes rather than within one.

    Sigma_3 and Sigma_4 - zirconium - are the two transverse REPRESENTATIONS,
    and the record deliberately does not say which of them is polarised along
    c, because the paper labels them by representation and guessing would put
    a real datum under a made-up name.  Those points are taken to the nearer
    of the two transverse branches OF THE RIGHT acoustic or optic character:
    a two-way choice, not the six-way one they had before.

    nu1..nu6 are zone-boundary frequencies the paper lists by index with no
    polarisation at all, and they keep the six-way nearest read.

    Returns (value, how) with how in {"label", "pair", "nearest"} so the
    caller can report how much of the comparison is actually labelled.
    """
    f6 = np.asarray(f6, float)
    head = label.split("[")[0]
    got = hcp_classes(vec, qc)
    if got is None or head.startswith("nu"):
        return None, "nearest"
    cls, split = got
    #  members of each class, in frequency order: first acoustic, then optic
    band = {}
    for k in (0, 1, 2):
        idx = sorted((j for j in range(len(f6)) if cls[j] == k),
                     key=lambda j: f6[j])
        band[k] = idx
    L_, P_, T_ = band[0], band[1], band[2]

    def one(idx, which):
        if len(idx) < 2:
            return None
        return float(f6[idx[0] if which == "A" else idx[-1]])

    if head in ("LA", "LO"):
        v = one(L_, head[1])
        return (v, "label") if v is not None else (None, "nearest")
    which = head[1] if len(head) > 1 else "A"
    if head in ("TAperp", "TOperp"):
        v = one(P_, which)
        return (v, "label") if v is not None else (None, "nearest")
    if head in ("TApar", "TOpar"):
        v = one(T_, which)
        return (v, "label") if v is not None else (None, "nearest")
    if head in ("TA", "TO"):
        #  q along c, where the two transverse classes are degenerate; the
        #  paper does not separate them there either
        vs = [x for x in (one(P_, which), one(T_, which)) if x is not None]
        if not vs:
            return None, "nearest"
        return float(np.mean(vs)), "label" if not split else "pair"
    if head.startswith("TA_") or head.startswith("TO_"):
        vs = [x for x in (one(P_, which), one(T_, which)) if x is not None]
        if not vs:
            return None, "nearest"
        return vs, "pair"
    return None, "nearest"
#  ---- end of lift --------------------------------------------------------


def parse_ev(path):
    """[(frequencies, eigenvectors with modes as columns)] in file order"""
    out, f, vecs, cur = [], [], [], None
    for ln in open(path, encoding="utf-8", errors="replace"):
        s = ln.strip()
        if s.startswith("# q-point:"):
            if f:
                out.append((np.array(f), np.array(vecs).T))
            f, vecs, cur = [], [], None
        elif s.startswith("# frequency"):
            f.append(float(s.rsplit(":", 1)[1]))
            cur = []
            vecs.append(cur)
        elif s and not s.startswith("#") and cur is not None:
            g = [float(x) for x in s.split(":")[0].split()[1:]]
            cur.extend(complex(g[k], g[k + 1]) for k in range(0, len(g), 2))
    if f:
        out.append((np.array(f), np.array(vecs).T))
    return out


def main():
    el, evfile, qfile = sys.argv[1], sys.argv[2], sys.argv[3]
    meta = json.load(open(qfile))
    rows, n = meta["idx"], meta["n"]
    ev = parse_ev(evfile)
    if len(ev) != 2 * n:
        print(json.dumps({"el": el,
                          "error": "%d q blocks, %d expected" % (len(ev),
                                                                 2 * n)}))
        return
    got, err, scatter, near = [], [], [], []
    #  how each hcp point was resolved: fully labelled, narrowed to a
    #  transverse pair, or left to the six-way nearest read
    how_n = {}
    nu_max = 0.0
    for k, r in enumerate(rows):
        f, vec = ev[k]
        fr = ev[k + n][0]
        #  q and a symmetry image of q must give the same frequencies for any
        #  correct force constants.  These are measured, with no symmetry
        #  imposed, so they do not - and the difference is the resolution of
        #  this run, reported next to the error rather than hidden in it.
        scatter.append(float(np.abs(np.sort(f) - np.sort(fr)).max()))
        if r["hcp"]:
            val, how = pick_hcp(f, vec, r.get("qc"), r["label"])
            how_n[how] = how_n.get(how, 0) + 1
            if how == "pair":
                g = float(min(val, key=lambda x: abs(x - r["nu"])))
            elif val is None:
                g = float(min(f, key=lambda x: abs(x - r["nu"])))
            else:
                g = float(val)
        elif r["nrm"] is None:
            g = float(min(f, key=lambda x: abs(x - r["nu"])))
        else:
            g = pick_polar(f, vec[:3], r["nrm"], r["label"])
        got.append(round(g, 4))
        err.append(abs(g - r["nu"]))
        #  the same point scored by nearest branch as well, on purpose.  It
        #  is the reading most of the literature uses, it flatters by 0.1 to
        #  14.8 points, and carrying both is what made a frame bug visible:
        #  a published EAM cannot be 1.8 % out by one rule and 9.3 % out by
        #  the other, and it was not - the polarised column was wrong.
        near.append(abs(float(min(f, key=lambda x: abs(x - r["nu"])))
                        - r["nu"]))
        nu_max = max(nu_max, r["nu"])
    print(json.dumps({
        "el": el, "n": len(err), "nu_max": round(nu_max, 4),
        "mae": round(float(np.mean(err)), 4),
        "pct": round(100.0 * float(np.mean(err)) / nu_max, 2),
        "sc": round(float(np.mean(scatter)), 4),
        "sc_max": round(float(max(scatter)), 4),
        "sc_pct": round(100.0 * float(np.mean(scatter)) / nu_max, 2),
        "near_pct": round(100.0 * float(np.mean(near)) / nu_max, 2),
        "how": how_n,
        "got": got,
    }))


if __name__ == "__main__":
    main()
