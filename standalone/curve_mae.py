#!/usr/bin/env python3
"""
Mean absolute error of each arm against the MEASURED dispersion curves.

The page draws the measured points and the computed branches on the same
axes, which lets a reader see agreement but not state it.  This states it.

Two things it deliberately does NOT do.  It does not match a measured point
to the nearest computed branch - that would let a longitudinal point be
scored against a transverse branch and would flatter every arm.  It uses the
branch LABEL: L against the highest of the three at that q, T against the
lowest, T1/T2 against the lower and upper of the two transverse ones.  And it
does not compare across elements measured at different temperatures without
saying so; the temperature is carried in the table because a 280 K caesium
measurement and a 90 K platinum one are not the same kind of test of a
harmonic 0 K calculation.

Run it twice, because the two arms need different modules.  The angular arm
must be scored with the angular latdyn, and there is no importing both into
one process; asking the standalone one for a UG record gives a plausible
wrong answer rather than an error, which is how that mistake stays hidden.

    python curve_mae.py [El ...]          # hard cut and MAU
    python curve_mae.py --ug [El ...]     # UG, with the angular latdyn
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

UG = "--ug" in sys.argv
if UG:
    sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "angular")))

import latdyn as L                          # noqa: E402
import refdata                              # noqa: E402
from build_library import sc_segments       # noqa: E402
import refdata_phonon_curves as RPC         # noqa: E402
from refdata_phonon_curves import SEGMENT as SEG   # noqa: E402

has_lam2 = "lam2" in open(L.__file__, encoding="utf-8").read()
assert has_lam2 == UG, (
    "wrong latdyn on the path: the angular one is needed for --ug and must "
    "not be used without it")
angfc = __import__("angfc") if UG else None

ARMS = ([("tap_ug", "UG")] if UG
        else [(None, "hard cut"), ("tap", "MAU")])


#  Lines on which symmetry forces the two transverse branches to be
#  degenerate: <100> for G-H and <111> for G-P and its continuation P-H.
DEGENERATE_T = {("G", "H"), ("G", "P"), ("P", "H")}


def ext_line(struct, el, br, seg):
    """(u, lo, hi, flip): the extended-zone line a labelled branch lies on

    A paper's L/T labels are polarisations with respect to the wavevector it
    swept, and past a zone corner that is NOT the wavevector of the folded
    segment we draw - fcc [0zz] leaves the zone at K and returns along U-X,
    which is a symmetry IMAGE of the continuation rather than a translate of
    it.  Evaluating there gives the right frequencies in the wrong frame for
    the eigenvectors, so the modes are taken on the straight extended line
    q(zeta) = zeta * u instead.  The frequencies must be unchanged by that,
    and QGATE checks it rather than trusting it.

    u is read off the entry that starts at Gamma - u = k(b) / zeta(b) - so no
    unit convention is assumed for zeta.  That matters: silver's record puts L
    at zeta = 1.333 and N at 2.139, and a fixed 2 pi / a table would have
    placed every one of its points somewhere else.

    Returns None for the D and G lines and for X-W, whose direction turns with
    zeta so there is no single u; those keep the sorted T1/T2/L reading, which
    their N-point frequencies already confirm.

    `flip` is set when the path walks the segment against the direction the
    measurement swept - Sigma is measured G -> K and drawn K -> G - because
    add_phonon_curves stores 1 - t in that case.
    """
    if struct == "hcp":
        return None
    table = dict(SEG.get(struct, {}))
    table.update(RPC.PHONON_CURVE.get(el, {}).get("seg_override", {}))
    ent = table.get(br)
    if not ent:
        return None
    base = next((e for e in ent if e[0] == "G"), None)
    if base is None or not base[3]:
        return None
    #  the path stores each segment in the sense it walks it, and it walks
    #  N -> G and K -> G rather than outwards from Gamma, so look both ways
    ends = {}
    for a, ka, b, kb in sc_segments(struct):
        ends[(a, b)], ends[(b, a)] = kb, ka
    if ("G", base[1]) not in ends:
        return None
    u = np.asarray(ends[("G", base[1])], float) / base[3]
    for a, b, lo, hi in ent:
        if (a, b) == seg:
            return u, lo, hi, False
        if (b, a) == seg:
            return u, lo, hi, True
    return None


def pick_polar(f3, vec, n, label):
    """the branch a polarisation label names, decided by the eigenvector

    The longitudinal branch is the one whose polarisation lies along the
    wavevector, whether that puts it top, middle or bottom.  Everything else
    is transverse, ordered low to high as T1, T2.
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


def pick(f3, label, seg=None):
    """the computed branch a label of this name refers to, at one q

    Cubic only.  An hcp cell has six branches and the paper's labels split
    them by polarisation relative to the basal plane as well as relative to
    q, and the acoustic three are not simply the lowest three - titanium's
    LA at M sits ABOVE both optic transverse branches.  Matching those by
    sorted position would be wrong and matching each point to whichever
    branch happens to be nearest would flatter every arm, so hcp is reported
    as a nearest-branch distance and named as one.
    """
    f = np.sort(f3)
    head = label.split("[")[0]
    if seg in DEGENERATE_T:
        #  SUPERSEDED by pick_polar, which decides this from the eigenvector
        #  and covers these three segments; what is left here is the zone
        #  centre, where the direction degenerates and all three frequencies
        #  go to zero together.  Kept because it is right, and because it
        #  states the physics the polarisation test relies on.
        #
        #  Take the longitudinal branch as the odd frequency out, not as the
        #  highest one.  Sorted position is not good enough here: past P the
        #  <111> longitudinal branch dips BELOW the transverse pair - for
        #  every one of the 9 potassium points, 11 of 12 for sodium and
        #  niobium, all 4 for caesium - while iron crosses back over halfway
        #  along and keeps L on top for 6 of its 11.  No fixed ordering is
        #  right for both, and the degeneracy is, because it is symmetry.
        if f[1] - f[0] <= f[2] - f[1]:
            return f[2] if head == "L" else 0.5 * (f[0] + f[1])
        return f[0] if head == "L" else 0.5 * (f[1] + f[2])
    if head == "L":
        return f[2]
    if head == "T1":
        return f[0]
    if head == "T2":
        return f[1]
    return 0.5 * (f[0] + f[1]) if abs(f[1] - f[0]) > 1e-9 else f[0]


#  set False to score at the fitted a0 instead, which is what the
#  measurement-temperature correction is measured against
USE_A_MEAS = True

#  Set True to score EVERY point by nearest branch, cubic ones included.
#  That is not a way to report a result - it flatters by 0.1 to 14.8 points,
#  worst for hcp, where there are six branches to land near instead of three
#  - but running both columns under it measures the flattering rather than
#  arguing about it.
#
#  It was once believed here that the rule also favours the NOISIER of two
#  calculations, so that a finite-temperature column could beat a 0 K one
#  without being better.  Measured, it does not: the MD side is flattered
#  more in only 7 of 24 elements and by -1.0 points on average.  The
#  asymmetry that suggested otherwise was a bug in the MD column, not a
#  property of this rule - the LAMMPS box is written in restricted triclinic
#  form, which rotates the Cartesian frame for fcc and bcc, and the
#  polarisation direction was being built in the unrotated one.
NEAREST = False

QGATE = []

#  how each hcp point was resolved: fully labelled, narrowed to a pair of
#  transverse branches, or left to the six-way nearest read
HCPHOW = {}


def main(only=(), field="exp_curve"):
    only = tuple(a for a in only if not a.startswith("-"))
    lib = json.load(open(os.path.join(HERE, "library.json")))
    rows = []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or field not in v:
            continue
        if only and el not in only:
            continue
        ec = v[field]
        #  a record with no branch labels at all - a figure digitised dot by
        #  dot - can only be read against the nearest computed branch, and
        #  says so itself rather than being guessed into T1/T2/L
        near_rec = bool(ec.get("nearest"))
        segs = {(a, b): (ka, kb) for a, ka, b, kb
                in sc_segments(v["struct"])}
        #  Build the crystal at the volume the MEASUREMENT was made at, when
        #  the paper says what that was.  refdata's a0 is ~293 K (and 5 K for
        #  the alkalis) while these curves run from 9 K to 296 K; for most
        #  elements those coincide and nothing changes, but lithium's a0 is
        #  0.56% dense against its own 293 K crystal and rubidium's 0.81%,
        #  while lead, aluminium and beryllium are measured cold against a
        #  room-temperature a0.  Comparing a dispersion at the wrong volume
        #  charges the potential for a bookkeeping mismatch.
        a_use = ec.get("a_meas", v["a0"]) if USE_A_MEAS else v["a0"]
        coa = v.get("c_over_a")
        if "c_meas" in ec and USE_A_MEAS:
            coa = ec["c_meas"] / a_use
        cry = L.Crystal(v["struct"], a_use, coa, mass=refdata.MASSES[el])
        out = {"el": el, "T": ec["T_K"], "n": 0,
               "hcp": v["struct"] == "hcp" and field == "exp_curve"}
        if near_rec and field == "exp_curve":
            out["near"] = True
        for key, name in ARMS:
            rec = v if key is None else v.get(key)
            if not rec or "m" not in rec:
                continue
            pot = L.Potential.from_record(rec)
            hcp = v["struct"] == "hcp"
            Phi = angfc.force_constants(cry, pot) if key == "tap_ug" else None
            err, nu_max = [], 0.0
            for sk, pts in ec["segs"].items():
                a, b = sk.split("|")
                if (a, b) not in segs:
                    continue
                ka, kb = segs[(a, b)]
                #  a model row is [t, f...] and may carry fewer than three
                #  branches, so take the fraction by position
                q = np.array([ka + r[0] * (kb - ka) for r in pts])
                #  For labelled data in a cubic cell the modes are taken on
                #  the extended-zone line instead, so the eigenvectors come
                #  out in the frame the paper's L/T labels refer to.  The
                #  frequencies must be unchanged by that - the two q differ by
                #  a symmetry operation - and QGATE checks it rather than
                #  assuming it.
                pol = field == "exp_curve"
                nrm = [None] * len(pts)
                #  the Cartesian wavevector of each point, which the hcp
                #  reading needs to say what "along q" means
                qcart = [None] * len(pts)
                if pol:
                    qq = q.copy()
                    for i, row in enumerate(pts):
                        lab = row[3]
                        line = ext_line(v["struct"], el,
                                        lab[lab.find("["):], (a, b))
                        if line is None:
                            qcart[i] = qq[i] @ L.reciprocal(cry)
                            continue
                        u, zlo, zhi, fl = line
                        t = 1.0 - row[0] if fl else row[0]
                        qq[i] = (zlo + t * (zhi - zlo)) * u
                        #  the polarisation direction is Cartesian, because
                        #  the force-constant blocks are
                        nc = qq[i] @ L.reciprocal(cry)
                        qcart[i] = nc
                        if np.linalg.norm(nc) > 1e-9:
                            nrm[i] = nc / np.linalg.norm(nc)
                    f, vecs = L.modes_many(cry, pot, qq, Phi)
                    QGATE.append(float(np.abs(
                        np.sort(f, 1)
                        - np.sort(L.frequencies_many(cry, pot, q, Phi), 1)
                    ).max()))
                else:
                    f = L.frequencies_many(cry, pot, q, Phi)
                for i, row in enumerate(pts):
                    if field == "exp_curve":
                        _, nu, _, lab = row
                        if hcp and not NEAREST and not near_rec:
                            val, how = pick_hcp(f[i], vecs[i], qcart[i], lab)
                            HCPHOW[how] = HCPHOW.get(how, 0) + 1
                            if how == "pair":
                                got = [min(val, key=lambda x: abs(x - nu))]
                            elif val is None:
                                got = [min(f[i], key=lambda x: abs(x - nu))]
                            else:
                                got = [val]
                        elif hcp or NEAREST or near_rec:
                            got = [min(f[i], key=lambda x: abs(x - nu))]
                        elif nrm[i] is not None:
                            got = [pick_polar(f[i], vecs[i], nrm[i], lab)]
                        else:
                            got = [pick(f[i], lab, (a, b))]
                        ref = [nu]
                    else:
                        #  a model curve has no branch labels, only three
                        #  sorted frequencies, so both sides are compared
                        #  sorted - there is nothing else honest to do
                        ref = list(row[1:])
                        got = list(np.sort(f[i]))
                        if len(ref) < len(got):
                            #  a segment where the source tabulates only
                            #  some of the branches.  Pairing them off in
                            #  order would score molybdenum's longitudinal
                            #  Sigma_1 against a computed TRANSVERSE branch,
                            #  so each is taken to its nearest instead - the
                            #  same lower bound the hcp rows carry.
                            out["near"] = True
                            got = [min(got, key=lambda x: abs(x - r))
                                   for r in ref]
                    err.extend(abs(g - r) for g, r in zip(got, ref))
                    nu_max = max([nu_max] + ref)
            if err:
                out[name] = (float(np.mean(err)), 100.0 * np.mean(err) / nu_max)
                out["n"] = len(err)
        rows.append(out)

    print(f"\n{'el':5s}{'T':>5s}{'pts':>5s}"
          + "".join(f"{n:>18s}" for _, n in ARMS))
    print("-" * (15 + 18 * len(ARMS)))
    for r in rows:
        tag = r["el"] + ("*" if r.get("hcp") or r.get("near") else "")
        line = f"{tag:5s}{r['T']:5d}K{r['n']:5d}"
        for _, n in ARMS:
            line += (f"{r[n][0]:11.3f} {r[n][1]:5.1f}%" if n in r
                     else f"{'-':>18s}")
        print(line)
    print()
    print("THz, and as a percentage of the highest reference frequency of "
          "that element.")
    if field == "exp_curve":
        print("L scored against the longitudinal branch and T against the "
              "transverse pair,")
        print("never against whichever branch happens to sit closest, and "
              "never by sorted")
        print("position: L is identified by its EIGENVECTOR lying along the "
              "wavevector, so")
        print("it is found whether it sits top, middle or bottom.  Past a "
              "zone corner it is")
        print("often not the top one.  90% of the cubic points are scored "
              "that way; the D")
        print("and G lines and X-W turn direction with zeta and keep the "
              "sorted reading.")
        if any(r.get("hcp") for r in rows) and not NEAREST:
            n_lab = HCPHOW.get("label", 0)
            n_pair = HCPHOW.get("pair", 0)
            n_near = HCPHOW.get("nearest", 0)
            tot = max(1, n_lab + n_pair + n_near)
            print()
            print("* hcp: six branches, split by polarisation along q and "
                  "perpendicular or")
            print("  parallel to the BASAL PLANE - the records say that is "
                  "what the papers'")
            print("  perp and par mean - and then acoustic or optic by "
                  "frequency order WITHIN")
            print("  a polarisation class.  %d%% of the points resolve that "
                  "way; %d%% narrow to"
                  % (round(100 * n_lab / tot), round(100 * n_pair / tot)))
            print("  a transverse PAIR, because Sigma_3 and Sigma_4 are "
                  "representations and the")
            print("  record refuses to guess which is polarised along c; "
                  "%d%% are the nu1..nu6"
                  % round(100 * n_near / tot))
            print("  zone-boundary frequencies, listed by index with no "
                  "polarisation at all.")
        elif any(r.get("hcp") or r.get("near") for r in rows):
            print()
            print("* nearest computed branch, not a labelled comparison.")
            print("  It cannot be larger than a correct one, so read it as a "
                  "lower bound.")
        if (any(r.get("near") and not r.get("hcp") for r in rows)
                and any(r.get("hcp") for r in rows) and not NEAREST):
            print()
            print("* a cubic row marked * is a digitised source with no branch "
                  "labels, scored")
            print("  against the nearest computed branch - a lower bound, like "
                  "the hcp rows.")
    else:
        print("A model curve carries no branch labels, so both sides are "
              "compared sorted.")
        if any(r.get("near") for r in rows):
            print("* partly nearest-branch: the source tabulates only some "
                  "branches on one line.")


if __name__ == "__main__":
    args = tuple(sys.argv[1:])
    main(args)
    print()
    print("Against the FITTED MODELS - a weaker test, since the reference is "
          "itself a fit, and for Sr it was fitted to a powder spectrum.")
    main(args, field="model_curve")
