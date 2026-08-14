#!/usr/bin/env python3
"""
Assemble library.json entirely from latdyn.py.

Reads fit.json (continuous gamma) and computes, for every element:
    relaxed and frozen-ion elastic constants
    phonon dispersion on the Setyawan-Curtarolo path
    phonon thermodynamics on a Monkhorst-Pack mesh

Every number is a derivative of the potential itself, so the r-power gamma stays
a real number instead of being truncated by an external input format.

    python fit.py && python build_library.py && python make_gui.py
"""
import json, math, os, sys
import numpy as np
import latdyn as L
import mechanics
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
T0, DT, NT = 25.0, 25.0, 32               # 25 K -> 800 K
NQ_MESH = 8                               # thermodynamics mesh
NQ_PATH = 60                              # points per dispersion panel

#  ---------------------------------------------------------------------------
#  High-symmetry paths.
#
#  STANDARD: the Setyawan-Curtarolo / HPKOT convention (Comput. Mater. Sci. 49,
#  299 (2010); 128, 140 (2017)) that seekpath implements and that Materials
#  Project, phononwebsite and most published band structures use.  Using it means
#  our dispersion can be overlaid directly on theirs.  The fractional
#  coordinates below are for exactly the primitive cells latdyn.Crystal builds:
#      fcc  a1=(0,1/2,1/2) a2=(1/2,0,1/2) a3=(1/2,1/2,0)
#      bcc  a1=(-1/2,1/2,1/2) a2=(1/2,-1/2,1/2) a3=(1/2,1/2,-1/2)
#      hcp  hexagonal, 2 atoms
#  Verified by mapping each point back to Cartesian: fcc X -> (0,1,0) 2pi/a,
#  K -> (3/4,3/4,0), L -> (1/2,1/2,1/2), W -> (1,1/2,0), U -> (1/4,1,1/4).
#
#  PAPERS: three separate branches out of Gamma, the layout of figure 1 in the
#  Akgun-Ugur papers.  Kept because it is what those figures show.
SC_POINTS = {
    "fcc": {"G": (0, 0, 0), "X": (0.5, 0.0, 0.5), "L": (0.5, 0.5, 0.5),
            "W": (0.5, 0.25, 0.75), "U": (0.625, 0.25, 0.625),
            "K": (0.375, 0.375, 0.75)},
    "bcc": {"G": (0, 0, 0), "H": (0.5, -0.5, 0.5), "P": (0.25, 0.25, 0.25),
            "N": (0.0, 0.0, 0.5)},
    "hcp": {"G": (0, 0, 0), "M": (0.5, 0.0, 0.0), "K": (1/3., 1/3., 0.0),
            "A": (0.0, 0.0, 0.5), "L": (0.5, 0.0, 0.5),
            "H": (1/3., 1/3., 0.5)},
}
#  Segments; a break (discontinuity) is a "|" entry.  These are the exact
#  sequences Materials Project's pheasy runs sample, read back off the fetched
#  band structures, so the two dispersion views share one x axis and an element
#  with MP data and one without can be compared by eye.  Note fcc: the
#  Setyawan-Curtarolo path breaks at
#  W|L, but MP's q-list runs a real straight segment from W to L (51 sampled
#  points), so keeping the break here put a gap in our plot where theirs has
#  data.  The remaining "|" are genuine discontinuities in both - no data
#  exists there and both views draw a gap.
SC_PATH = {
    "fcc": ["G", "X", "W", "K", "G", "L", "U", "W", "L", "K", "|", "U", "X"],
    "bcc": ["G", "H", "N", "G", "P", "H", "|", "P", "N"],
    "hcp": ["G", "M", "K", "G", "A", "L", "H", "A", "|", "L", "M",
            "|", "K", "H"],
}


def sc_segments(struct):
    """[(label_a, k_a, label_b, k_b), ...] for the standard path"""
    pts, seq = SC_POINTS[struct], SC_PATH[struct]
    segs, i = [], 0
    while i < len(seq) - 1:
        a, b = seq[i], seq[i+1]
        if a == "|":
            i += 1
            continue
        if b == "|":
            i += 2
            continue
        segs.append((a, np.array(pts[a]), b, np.array(pts[b])))
        i += 1
    return segs


def recip(struct, c_over_a):
    if struct == "fcc":
        return np.array([[-1., 1, 1], [1, -1, 1], [1, 1, -1]])
    if struct == "bcc":
        return np.array([[0., 1, 1], [1, 0, 1], [1, 1, 0]])
    s3 = math.sqrt(3.0)
    return np.array([[1., 1/s3, 0], [0, 2/s3, 0], [0, 0, 1/c_over_a]])


def build(el, p, e):
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    #  lam2/lam4 travel with the parameter record.  A MAU fit simply has no
    #  such keys and gets zero, which is the same potential it was fitted with;
    #  omitting them here would have built every UG entry without its angular
    #  factor and reported MAU numbers under a UG label.
    pot = L.Potential(p["m"], p["D"], p["alpha"], p["r0"], p["gamma"],
                      C=p["C"], alpha3=p["alpha3"],
                      rcut2=p["rcut2"], rcut3=p["rcut3"],
                      lam2=p.get("lam2", 0.0), lam4=p.get("lam4", 0.0))
    return cry, pot


def record(el, p, e=None):
    """
    Everything the viewer shows for one element, from the parameters alone.

    Split out of main() so that export_ug.py builds UG entries through exactly
    this code path.  The alternative - a second, shorter record written by hand
    on the angular side - is how the page ended up showing UG with elastic
    constants and nothing else, while MAU carried mechanics, thermodynamics and
    a dispersion curve.  One producer, one shape.
    """
    e = e or refdata.ELEMENTS[el]
    cry, pot = build(el, p, e)
    nat = len(cry.frac)

    msgs = L.check_force_constants(cry, pot)
    C, Cb = L.elastic(cry, pot)
    s = L.stress(cry, pot)
    exp = e["Cij"]
    keys = ["C11", "C12", "C44"] + (["C13", "C33"]
                                    if e["struct"] == "hcp" else [])
    idx = {"C11": (0, 0), "C12": (0, 1), "C13": (0, 2),
           "C33": (2, 2), "C44": (3, 3)}

    d = dict(struct=e["struct"], a0=e["a0"], c_over_a=e.get("c_over_a"),
             mass=refdata.MASSES[el], exp=exp,
             m=p["m"], gamma=p["gamma"], D=p["D"], alpha=p["alpha"],
             r0=p["r0"], alpha3=p["alpha3"], s3=p["s3"], C=p["C"],
             lam2=p.get("lam2", 0.0), lam4=p.get("lam4", 0.0),
             dnn=p["dnn"], rcut2=p["rcut2"], rcut3=p["rcut3"],
             ntrip=len(L.triplets(cry, pot)) // nat,
             rms=p["score"]*100.0,
             at_bound=p.get("at_bound", []),
             Ecoh=[L.energy(cry, pot), -e["Ecoh"]],
             P_resid=float(s[:3].mean()),
             B=[float((C[0, 0]+C[1, 1]+C[2, 2]
                       + 2*(C[0, 1]+C[0, 2]+C[1, 2]))/9.0), e["B"]],
             fc_check=msgs)
    for k in keys:
        if k in exp:
            d[k] = [float(C[idx[k]]), exp[k]]
    d["frozen"] = {k: float(Cb[idx[k]]) for k in keys}

    #  the full tensor, and everything it implies on its own
    d["Cfull"] = [[round(float(x), 4) for x in row] for row in C]
    try:
        d["mech"] = mechanics.analyse(
            C, mass_amu=refdata.MASSES[el]*nat, volume_A3=cry.vol,
            natoms=nat)
        d["mech_planes"] = {pl: mechanics.plane_curves(C, pl)
                            for pl in ("xy", "xz", "yz")}
    except np.linalg.LinAlgError:
        d["mech"] = None              # singular C: the Born check says why

    #  ---- dispersion ----
    Phi = L.force_constants(cry, pot)
    B_ = recip(e["struct"], e.get("c_over_a") or 1.0)

    def sample(ka, kb, npts):
        ks = ka[None, :] + np.linspace(0, 1, npts)[:, None]*(kb - ka)[None, :]
        br = L.frequencies_many(cry, pot, ks, Phi) * L.CM1
        return [[round(float(v), 1) for v in br[:, b]]
                for b in range(br.shape[1])]

    std = []
    for (la, ka, lb, kb) in sc_segments(e["struct"]):
        std.append({"a": la, "b": lb, "n": NQ_PATH,
                    "len": float(np.linalg.norm((kb - ka) @ B_)),
                    "branches": sample(ka, kb, NQ_PATH)})

    #  ---- thermodynamics ----
    f = L.spectrum(cry, pot, nq=NQ_MESH)
    th = []
    for n in range(NT):
        T = T0 + n*DT
        t = L.thermo(f, T, nat)
        th.append({"T": T, "zpe": t["zpe"], "F": t["F"],
                   "S": t["S"], "Cv": t["Cv"]})
    d["ld"] = {"std": std, "thermo": th, "maxfreq": float(f.max()*L.CM1)}
    t = refdata.THERMO_298.get(el)
    if t:
        d["S298"], d["Cp298"] = t
    return d


def main(only=None):
    fit = json.load(open(os.path.join(HERE, "fit.json")))
    prev = {}
    path = os.path.join(HERE, "library.json")
    if os.path.exists(path):
        try:
            prev = json.load(open(path))
        except ValueError:
            pass

    out = {}
    for i, (el, p) in enumerate(sorted(fit.items()), 1):
        if only and el not in only:
            if el in prev:
                out[el] = prev[el]
            continue
        d = record(el, p)
        #  keep externally fetched reference data across rebuilds; it is a
        #  property of the element, not of the fit
        if el in prev and "mp" in prev[el]:
            d["mp"] = prev[el]["mp"]
        out[el] = d
        msgs = d["fc_check"]
        bad = " FC-CHECK:" + str(msgs[:1]) if msgs else ""
        print(f"  [{i:2d}/{len(fit)}] {el:3s} rms={d['rms']:5.1f}%  "
              f"gamma={p['gamma']:6.3f}  max={d['ld']['maxfreq']:6.1f} cm-1  "
              f"S298={d['ld']['thermo'][11]['S']:6.2f}{bad}", flush=True)


    json.dump(out, open(path, "w"), indent=1, sort_keys=True)
    print(f"wrote {path} ({len(out)} elements)")


if __name__ == "__main__":
    main(set(sys.argv[1:]) or None)
