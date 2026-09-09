#!/usr/bin/env python3
"""Mean absolute error of the FINITE-TEMPERATURE dispersion against measurement.

curve_mae.py scores the 0 K harmonic dispersion the library computes from the
force constants.  This scores the dispersion LAMMPS measures by displacement
correlation during equilibrium MD at the temperature the neutron experiment
was run at, which is the spectrum the experiment actually sees: renormalised
by anharmonicity and by the crystal's own thermal expansion.

WHICH ARM.  Every finite-temperature run in this study is the MAU arm - the
tapered record, `lib[el]["tap"]`, shipped as `<El>_taper.ugur` and run under
`pair_style ugur`.  Not UG (that would need `ugur/ang` and the `.ugur.ang`
files), and not the `rc` / `rc_ug` re-cut candidates.  All four arms have MD
stability screens and LAMMPS files, so this was a choice rather than a
constraint: MAU is the primary record and four arms would have cost four
times the cluster time.  Whatever this column concludes is therefore a
statement about MAU, not about the library.  The 0 K columns it is compared
against are read from the same arm.

It is the same comparison in every other respect, deliberately.  The q-points
are the ones curve_mae builds, including the extended-zone remapping that puts
the eigenvectors in the frame the paper's L/T labels refer to; the branch a
label names is decided by its EIGENVECTOR, not by sorted position; hcp keeps
the nearest-branch reading and is named as a lower bound.  Anything that
differs between the two columns is then the physics, not the bookkeeping.

An earlier pass scored this column by nearest branch because phana's menu
option 4 returns frequencies only.  Option 9 returns eigenvectors as well and
costs the same, so that concession was unnecessary and is withdrawn.

Two stages, because the frequencies come from the cluster:

    python mdcurve_mae.py emit  El [El ...]   # -> mdq_El.txt, mdq_El.json
    (run phana option 9 on the cluster with mdq_El.txt -> ev_El.dat)
    python mdcurve_mae.py score El [El ...]   # -> the table

SYMMETRY GATE.  The measured force constants have no symmetry imposed on
them, so two q-points that a cubic crystal must make identical need not come
out identical.  Every q is therefore emitted a second time, rotated by a
proper symmetry operation of the lattice; the spread between the two is the
noise floor of this measurement and is reported next to the error.  It is not
small: tantalum's two N-point entries differ by 12 % on the top branch.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import latdyn as L                                  # noqa: E402
from build_library import sc_segments               # noqa: E402
from curve_mae import ext_line, pick_polar          # noqa: E402


def rot_cart(struct):
    """a proper rotation of the lattice, in Cartesian coordinates

    Cubic: the 3-fold about [111], which permutes the axes.  Hexagonal: the
    3-fold about c.  Both are symmetries of the crystal, so the frequencies
    at q and at R q are equal for any correct set of force constants.
    """
    if struct == "hcp":
        c, s = math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3)
        return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    return np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])


def lammps_frame(cry):
    """the lattice as LAMMPS holds it: rotated into restricted triclinic form

    mkprim.py writes the box with a1 along x and a2 in the xy plane, because
    that is the only triclinic form LAMMPS accepts.  That is a ROTATION of the
    Cartesian frame, and phana returns its eigenvectors in the rotated frame
    while latdyn works in the unrotated one.  Projecting one onto a direction
    computed in the other is meaningless: measured on copper it put the
    longitudinal branch on a different mode from the analytic calculation at
    70 % of the points, with an average of 2.2 THz between the two branches,
    and it inflated the whole finite-temperature column by about seven points.

    fcc and bcc are rotated; hcp is not, because a1 is already along x and a2
    already lies in the xy plane, so the hcp rows were never affected.
    """
    a1, a2, a3 = np.asarray(cry.lat, float)
    ax = np.linalg.norm(a1)
    u1 = a1 / ax
    bx = float(np.dot(a2, u1))
    by = math.sqrt(float(np.dot(a2, a2)) - bx * bx)
    cx = float(np.dot(a3, u1))
    cy = (float(np.dot(a2, a3)) - bx * cx) / by
    cz = math.sqrt(float(np.dot(a3, a3)) - cx * cx - cy * cy)
    return np.array([[ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]])


def reference(v, el):
    """(kind, record) - the measured curve if there is one, else the model

    `exp_curve` is preferred and `model_curve` is accepted, and the caller is
    told which so the difference reaches the printed table.  A model reference
    is a weaker test - the reference is itself a fit - and `curve_mae` already
    says so in its own output.

    Vanadium is why this exists rather than a refusal: it scatters neutrons
    almost entirely incoherently, which is what makes it the standard neutron
    calibrant and what makes a coherent measurement of it impossible, so its
    only reference is Colella and Batterman's Born-von Karman model.  That
    model is gated on the elastic constants it was NOT fitted to.
    """
    if "exp_curve" in v:
        return "exp", v["exp_curve"]
    if "model_curve" in v:
        return "model", v["model_curve"]
    raise AssertionError("%s: ne olculen ne model egrisi var" % el)


def build(el, lib):
    """(q rows, index rows, n) for one element, in phana's fractional units"""
    v = lib[el]
    kind, ec = reference(v, el)
    segs = {(a, b): (ka, kb) for a, ka, b, kb in sc_segments(v["struct"])}
    #  The MD box was built at the library a0 and then let go under NPT, so
    #  the volume is the potential's own at this temperature; the a_meas
    #  correction curve_mae applies has nothing to correct here.
    cry = L.Crystal(v["struct"], v["a0"], v.get("c_over_a"))
    hcp = v["struct"] == "hcp"
    R = L.reciprocal(cry)
    Rinv = np.linalg.inv(R)
    #  Cartesian quantities that will be compared with phana's eigenvectors
    #  must be built in the frame LAMMPS rotated the box into, not in
    #  latdyn's.  The fractional wavevectors themselves are frame-free and
    #  stay as they are.
    RL = 2.0 * np.pi * np.linalg.inv(lammps_frame(cry)).T
    rot = rot_cart(v["struct"])

    qs, idx = [], []
    for sk, pts in ec["segs"].items():
        a, b = sk.split("|")
        if (a, b) not in segs:
            continue
        ka, kb = np.asarray(segs[(a, b)][0], float), \
            np.asarray(segs[(a, b)][1], float)
        for i, row in enumerate(pts):
            #  A model row is [t, nu1, nu2, nu3] with no branch label; a
            #  measured row is [t, nu, sigma, label].  The model side keeps
            #  the whole sorted triple and is compared sorted, because
            #  assigning a label it does not carry is how molybdenum's
            #  transverse branch gets read as longitudinal between P and H.
            if kind == "model":
                nus = sorted(float(x) for x in row[1:4])
                #  DROP THE GAMMA ENDPOINTS.  Three of vanadium's 240 model
                #  rows are identically zero - the Gamma end of each segment
                #  that touches it - and the acoustic modes are zero there by
                #  construction on BOTH sides.  Scoring them adds three
                #  perfect agreements nothing earned, and phana has nothing
                #  useful to say at q = 0 either.
                if max(nus) < 1e-6:
                    continue
                q = (np.asarray(segs[(a, b)][0], float)
                     + row[0] * (np.asarray(segs[(a, b)][1], float)
                                 - np.asarray(segs[(a, b)][0], float)))
                qs.append(np.asarray(q, float))
                idx.append({"seg": sk, "i": i, "label": None, "nu": nus[0],
                            "nus": nus, "nrm": None, "hcp": hcp,
                            "qc": (np.asarray(q, float) @ RL).tolist()})
                continue
            nu, lab = row[1], row[3]
            q = ka + row[0] * (kb - ka)
            nrm = None
            if not hcp:
                line = ext_line(v["struct"], el, lab[lab.find("["):], (a, b))
                if line is not None:
                    u, zlo, zhi, fl = line
                    t = 1.0 - row[0] if fl else row[0]
                    q = (zlo + t * (zhi - zlo)) * np.asarray(u, float)
                    nc = q @ RL
                    if np.linalg.norm(nc) > 1e-9:
                        nrm = (nc / np.linalg.norm(nc)).tolist()
            qs.append(np.asarray(q, float))
            #  the Cartesian wavevector, which the hcp reading needs to say
            #  what "along q" and "perpendicular to the basal plane" mean
            idx.append({"seg": sk, "i": i, "label": lab, "nu": nu,
                        "nrm": nrm, "hcp": hcp,
                        "qc": (np.asarray(q, float) @ RL).tolist()})
    #  the rotated copies, appended as one block so the file stays in order
    n = len(qs)
    for q in list(qs):
        qs.append((q @ R @ rot.T) @ Rinv)
    return qs, idx, n


def parse_ev(path, want=None):
    """[(frequencies, eigenvectors with modes as columns)] in file order

    `want` is the number of q blocks the caller asked phana for.  Give it and
    a TRUNCATED file raises instead of being scored.

    That is not hypothetical.  On 2026-09-08 an abandoned line-range download
    left `ev_Cu_hard.dat` holding 308 of its 2,576 lines and `ev_Ni_hard.dat`
    768 of 2,184 - files that parse perfectly, contain real numbers, and would
    have produced a confident wrong answer from a fifth of the data.  Nothing
    in the format says how long it should be, so the count has to come from
    the caller, which knows: it wrote the q list.
    """
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
            #  "id  re im re im re im  : |e|" - one row per atom in the cell
            g = [float(x) for x in s.split(":")[0].split()[1:]]
            cur.extend(complex(g[k], g[k + 1]) for k in range(0, len(g), 2))
    if f:
        out.append((np.array(f), np.array(vecs).T))
    if want is not None and len(out) != want:
        raise AssertionError(
            "%s KESILMIS ya da eksik: %d q blogu var, %d bekleniyor. "
            "Puanlanmadi." % (os.path.basename(path), len(out), want))
    return out


def score(el, lib, tag=None):
    """tag lets several FITS of one element be scored without collision

    The study had one run per element, so `ev_<El>.dat` was unambiguous.  This
    batch has one per FIT - `Mo_lib` and `Mo_hard` are two different
    potentials for molybdenum - and naming both `ev_Mo.dat` would mean the
    second silently overwrote the first and the table reported one number
    twice.  So the eigenvector file is `ev_<tag>.dat` when a tag is given.
    """
    meta = json.load(open(os.path.join(HERE, "mdq_%s.json" % el)))
    rows, n = meta["idx"], meta["n"]
    #  the q list this element was emitted with says how many blocks phana
    #  was asked for: n real points plus n rotated copies
    with open(os.path.join(HERE, "mdq_%s.json" % el)) as fh:
        want = 2 * int(json.load(fh)["n"])
    ev = parse_ev(os.path.join(HERE, "ev_%s.dat" % (tag or el)), want)
    assert len(ev) == 2 * n, "%s: %d q blocks, %d expected" % (el, len(ev),
                                                               2 * n)
    err, nu_max, scatter = [], 0.0, []
    for k, r in enumerate(rows):
        f, vec = ev[k]
        fr = ev[k + n][0]
        scatter.append(float(np.abs(np.sort(f) - np.sort(fr)).max()))
        if r.get("nus") is not None:
            #  model reference: no labels on either side, so both are sorted
            ours = sorted(float(x) for x in f)[:len(r["nus"])]
            err.append(float(np.mean(np.abs(np.asarray(ours)
                                            - np.asarray(r["nus"])))))
            nu_max = max(nu_max, max(r["nus"]))
            continue
        if r["hcp"] or r["nrm"] is None:
            got = float(min(f, key=lambda x: abs(x - r["nu"])))
        else:
            got = pick_polar(f, vec[:3], r["nrm"], r["label"])
        err.append(abs(got - r["nu"]))
        nu_max = max(nu_max, r["nu"])
    mae = float(np.mean(err))
    kind, ec = reference(lib[el], el)
    return {"el": el, "tag": tag or el, "n": len(err), "T": ec["T_K"],
            "kind": kind,
            "mae": mae, "pct": 100.0 * mae / nu_max, "hcp": rows[0]["hcp"],
            "sc": float(np.mean(scatter)), "sc_max": float(max(scatter)),
            "sc_pct": 100.0 * float(np.mean(scatter)) / nu_max}


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "score"
    els = sys.argv[2:]
    lib = json.load(open(os.path.join(HERE, "library.json")))
    if mode == "emit":
        for el in els:
            qs, idx, n = build(el, lib)
            with open(os.path.join(HERE, "mdq_%s.txt" % el), "w") as fh:
                for q in qs:
                    fh.write("%.10f %.10f %.10f\n" % tuple(q))
            json.dump({"idx": idx, "n": n},
                      open(os.path.join(HERE, "mdq_%s.json" % el), "w"))
            print("%s: %d nokta + %d donmus kopya = %d q" % (el, n, n,
                                                             len(qs)))
        return

    out = []
    for spec in els:
        el, _, tag = spec.partition(":")
        try:
            out.append(score(el, lib, tag or None))
        except (OSError, AssertionError) as e:
            print("%s: %s" % (spec, e))
    head = "%-10s%5s%5s%10s%7s%16s%6s" % ("kayit", "T", "pts", "MAE", "",
                                         "simetri sacilmasi", "")
    print()
    print(head)
    print("-" * 62)
    for r in out:
        tag = (r["tag"] + ("*" if r["hcp"] else "")
               + ("~" if r.get("kind") == "model" else ""))
        print("%-10s%4dK%5d%10.3f%6.1f%%%12.3f%5.1f%%  (en buyuk %.2f)"
              % (tag, r["T"], r["n"], r["mae"], r["pct"], r["sc"],
                 r["sc_pct"], r["sc_max"]))
    print()
    print("THz, and as a percentage of the highest measured frequency of that")
    print("element.  The scatter column is the same quantity evaluated between")
    print("q and a symmetry image of q: it is what this measurement cannot")
    print("resolve, and an error below it is not a result.")
    if any(r.get("kind") == "model" for r in out):
        print()
        print("~ scored against a FITTED MODEL, not a measurement, and both")
        print("  sides compared sorted because the model carries no branch")
        print("  labels.  A weaker test: the reference is itself a fit.")
    if any(r["hcp"] for r in out):
        print()
        print("* hcp: distance to the NEAREST branch, a lower bound, as in "
              "curve_mae.")


if __name__ == "__main__":
    main()
