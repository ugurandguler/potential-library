#!/usr/bin/env python3
"""
Attach Born-von Karman dispersion curves, built from published force
constants, to library.json.

Stored under "model_curve", NOT "exp_curve", and the viewer draws them as
lines rather than as the open circles it uses for measured phonons.  The
distinction is the whole point: most elements on this page are compared
against individually measured frequencies, and these three are compared
against a model, because for each of them a table of frequencies is not what
the paper published.  Putting the two in the same field would make the page
say something it cannot support.

    Sr  no single crystal could be grown, so the phonons were never measured
        one at a time at all - the model is fitted to a powder spectrum
    Cr  measured with single crystals, but the paper tabulates the fitted
        force constants and plots the frequencies
    Rh  the same, a 24-parameter model the paper says reproduces its own
        measurements to 0.031 THz

Every model is checked before anything is written, and the check is not the
same for all three - see bvk_model.  A model that fails its check is not
written and stops the run.

    python add_bvk_curve.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bvk_model as B                      # noqa: E402
import fourier_model as F                  # noqa: E402
import dfm_model as DFM                    # noqa: E402
import refdata                             # noqa: E402

#  High-symmetry points in CARTESIAN units of 2 pi / a, which is the
#  convention bvk_model works in (R in units of the conventional cube edge),
#  and the segments in the sense build_library.SC_PATH walks them.
CART = {
    "fcc": {"G": (0.0, 0.0, 0.0), "X": (1.0, 0.0, 0.0), "W": (1.0, 0.5, 0.0),
            "K": (0.75, 0.75, 0.0), "L": (0.5, 0.5, 0.5),
            "U": (1.0, 0.25, 0.25)},
    "bcc": {"G": (0.0, 0.0, 0.0), "H": (1.0, 0.0, 0.0),
            "N": (0.5, 0.5, 0.0), "P": (0.5, 0.5, 0.5)},
}
DRAWN = {
    "fcc": [("G", "X"), ("X", "W"), ("W", "K"), ("K", "G"), ("G", "L"),
            ("L", "U"), ("U", "W"), ("W", "L"), ("L", "K"), ("U", "X")],
    "bcc": [("G", "H"), ("H", "N"), ("N", "G"), ("G", "P"), ("P", "H"),
            ("P", "N")],
}
NQ = 40

#  `why` says why THIS element is a model rather than a measurement, because
#  the reason is not the same for all three and the page used to state
#  strontium's reason for whatever element was on screen.  `gate` says what
#  the model was actually checked against.
POWDER = ("No single crystal of this element large enough for triple-axis "
          "spectroscopy has been grown, so its phonons were reached instead "
          "by fitting to time-of-flight spectra from a <em>powder</em>")
TABULATED = ("It was measured with single crystals, but what the paper "
             "tabulates is the fitted force constants rather than a list of "
             "frequencies")

MODELS = {
    "Sr": {"struct": "fcc", "shells": B.SR_293, "nat": 4, "T_K": 293,
           "ref": B.SR_REF, "check": B.check, "why": POWDER,
           "want": B.SR_ELASTIC, "src": "the same paper",
           "gate": "the elastic constants"},
    "Cr": {"struct": "bcc", "shells": B.CR_300, "nat": 2, "T_K": 300,
           "ref": B.CR_REF, "check": B.check_cr, "why": TABULATED,
           "want": None, "src": "measured",
           "gate": (f"the {B.CR_H_THZ} THz the paper quotes for the zone "
                    "boundary along [001], which this reconstruction "
                    "reproduces to 0.02 THz.  Its c&prime; lands within "
                    "0.2% of the measured value; its c<sub>44</sub> sits "
                    "about 30% above, which is a property of the published "
                    "fit and not of the reconstruction")},
    "Rh": {"struct": "fcc", "shells": B.RH_297, "nat": 4, "T_K": 297,
           "ref": B.RH_REF, "check": B.check_rh, "why": TABULATED,
           "want": B.RH_ELASTIC, "src": "the same paper",
           "gate": "the elastic constants"},
    #  V carries its own builder because its constants are a general TENSOR
    #  per shell, the form written for niobium, rather than the component
    #  dictionary the other three use.
    "V": {"struct": "bcc", "shells": B.V_T7, "nat": 2, "T_K": 296,
          "fc": lambda: B.tensor_force_constants(B.V_T7),
          "ref": B.V_REF, "check": B.check_V,
          "why": ("Vanadium scatters neutrons almost entirely "
                  "<em>incoherently</em> &mdash; which is what makes it the "
                  "standard neutron calibrant, and what makes coherent "
                  "measurements on it hard &mdash; so its dispersion was "
                  "taken by x-ray thermal diffuse scattering, and the paper "
                  "publishes the fitted force constants rather than a list "
                  "of frequencies"),
          "want": None, "src": "measured",
          "gate": ("the elastic constants, which this column was not fitted "
                   "to &mdash; the paper reports its elastic-constant "
                   "constrained variant as a separate and poorer fit.  The "
                   "reconstruction returns 228.3, 119.1 and 42.6 GPa against "
                   "measured 230.0 &plusmn; 5.0, 120.0 &plusmn; 4.0 and 43.1 "
                   "&plusmn; 0.4")},
}


#  Molybdenum is published as a Fourier series per branch per direction
#  rather than as a force-constant tensor, so it is generated separately -
#  see fourier_model.  Sigma_4 is not in that table, so [zz0] carries two
#  branches where the others carry three, and the row is simply shorter.
FOURIER = {
    "Mo": {"struct": "bcc", "T_K": 296, "coeffs": F.MO_296, "ref": F.MO_REF,
           "check": F.check,
           "on": {("G", "H"): ("[z00]", 0.0, 1.0),
                  ("N", "G"): ("[zz0]", 0.5, 0.0),
                  ("G", "P"): ("[zzz]", 0.0, 0.5),
                  ("P", "H"): ("[zzz]", 0.5, 1.0)},
           "why": ("The paper plots its measured frequencies and tabulates "
                   "the Fourier coefficients of each branch instead of the "
                   "frequencies themselves"),
           "gate": ("the Landolt-B&ouml;rnstein values this page already "
                    "carries for H and N, a different compilation that was "
                    "used to set nothing in the reconstruction: it returns "
                    "5.49&ndash;5.55 THz at H against their 5.52, and 8.15 "
                    "and 4.58 THz at N against their 8.14 and 4.56")},
}


#  Cobalt: the dipolar fluctuation model.  Its own dispersion was measured
#  and PLOTTED, never tabulated, so this is the only route onto the page, and
#  it is only taken because five of the other six elements in the same table
#  already carry measured curves here and can test the reconstruction first -
#  see dfm_model.check_ladder.
DFM_DRAWN = [("G", "M"), ("M", "K"), ("K", "G"), ("G", "A"), ("A", "L"),
             ("L", "H"), ("H", "A"), ("L", "M"), ("K", "H")]


def dfm_curves(lib):
    import refdata
    if not DFM.check_ladder():
        raise SystemExit("Co: the DFM ladder failed; nothing written")
    el = "Co"
    if el not in lib or not isinstance(lib[el], dict):
        print(f"{el:4s}  not in library.json, skipped")
        return
    a = refdata.ELEMENTS[el]["a0"]
    c = a * refdata.ELEMENTS[el]["c_over_a"]
    mass = refdata.MASSES[el]
    Brec = DFM.reciprocal(a, c)
    segs = {}
    for x, z in DFM_DRAWN:
        qa = np.array(DFM.HCP_PTS[x], float) @ Brec
        qz = np.array(DFM.HCP_PTS[z], float) @ Brec
        rows = []
        for k in range(NQ):
            t = k / (NQ - 1.0)
            f = DFM.freqs_THz(qa + t * (qz - qa), DFM.WAK_P[el],
                              DFM.WAK_T[el], DFM.WAK_S[el], a, c, mass)
            rows.append([round(t, 5)] + [round(float(v), 4) for v in f])
        segs[f"{x}|{z}"] = rows
    lib[el]["model_curve"] = {
        "T_K": 295, "kind": "dipolar fluctuation",
        "ref": DFM.CO_REF, "segs": segs, "shells": 3,
        "why": ("Cobalt absorbs neutrons strongly and scatters them weakly, "
                "so the only full measurement of hcp cobalt plots its "
                "frequencies and tabulates the fitted model instead.  "
                "Landolt-B&ouml;rnstein's cobalt table is no help either: it "
                "is fcc Co<sub>0.92</sub>Fe<sub>0.08</sub>, an alloy in the "
                "wrong phase"),
        "gate": ("the same model, with the same code, against the measured "
                 "curves this page already carries for the five other "
                 "elements of that table &mdash; 464 frequencies that "
                 "entered nothing.  Scandium and yttrium, which have no "
                 "dipolar term published, come out at 0.9% and 1.3% of "
                 "their top branch; titanium, zirconium and hafnium at 1.4%, "
                 "1.1% and 1.4%, each about four times better than the "
                 "ion-ion half alone, which is the dipolar term doing what "
                 "it is supposed to do"),
    }
    n = sum(len(v) for v in segs.values())
    print(f"{el:4s}{295:5d}K{3:8d}{len(segs):6d}{n:8d}")


def fourier_curves(lib):
    import refdata
    for el, m in sorted(FOURIER.items()):
        if not m["check"]():
            raise SystemExit(f"{el}: check failed; nothing written")
        if el not in lib or not isinstance(lib[el], dict):
            print(f"{el:4s}  not in library.json, skipped")
            continue
        mass = refdata.MASSES[el]
        segs = {}
        for (a, z), (d, lo, hi) in m["on"].items():
            rows = []
            for k in range(NQ):
                t = k / (NQ - 1.0)
                f = F.branches_at(m["coeffs"], d, lo + t * (hi - lo), mass)
                rows.append([round(t, 5)] + [round(x, 4) for x in f])
            segs[f"{a}|{z}"] = rows
        lib[el]["model_curve"] = {
            "T_K": m["T_K"], "kind": "Fourier interplanar force-constant",
            "ref": m["ref"], "segs": segs, "why": m["why"], "gate": m["gate"],
        }
        n = sum(len(v) for v in segs.values())
        print(f"{el:4s}{m['T_K']:5d}K{'-':>8s}{len(segs):6d}{n:8d}")


#  Axially symmetric models: two constants per shell rather than a tensor.
#  Tungsten's third-shell bending constant is illegible in the scan, so the
#  generator does not take a value - it takes whichever sign bvk_model.check_W
#  finds against the eight frequencies of the same paper's Table 1, which
#  entered the model nowhere.  If neither sign passes, nothing is written.
AS_CART = {"bcc": {"G": (0.0, 0.0, 0.0), "H": (0.0, 0.0, 1.0),
                   "N": (0.5, 0.5, 0.0), "P": (0.5, 0.5, 0.5)}}
AS_DRAWN = {"bcc": [("G", "H"), ("H", "N"), ("N", "G"), ("G", "P"),
                    ("P", "H"), ("P", "N")]}


def as_curves(lib):
    import refdata
    sh, rms, which = B.check_W()
    if sh is None:
        print("W: check failed; nothing written")
        return
    el = "W"
    mass = refdata.MASSES[el]
    fc = B.as_force_constants(sh)
    segs = {}
    for a, z in AS_DRAWN["bcc"]:
        qa = np.array(AS_CART["bcc"][a]); qz = np.array(AS_CART["bcc"][z])
        rows = []
        for k in range(NQ):
            t = k / (NQ - 1.0)
            f = B.freqs_THz(fc, qa + t * (qz - qa), mass)
            rows.append([round(t, 5)] + [round(float(x), 4) for x in f])
        segs[f"{a}|{z}"] = rows
    lib[el]["model_curve"] = {
        "T_K": 296, "kind": "third-neighbour axially symmetric",
        "ref": B.W_REF, "segs": segs, "shells": 3, "neighbours": 26,
        "why": ("Chen and Brockhouse plot tungsten's dispersion and tabulate "
                "only eight points of it, so the curve here is their force "
                "constants rather than their measurement"),
        "gate": (f"the eight frequencies of the same paper's Table 1, which "
                 f"the model never saw - it reproduces them to "
                 f"{rms:.2f} THz rms, and it was that comparison, not the "
                 f"scan, that settled the illegible sign of the third-shell "
                 f"bending constant ({which})"),
    }
    print(f"{el:4s}{296:5d}K{'-':>8s}{len(segs):6d}{sum(len(v) for v in segs.values()):8d}")


def main():
    for el, m in sorted(MODELS.items()):
        if not m["check"]():
            raise SystemExit(f"{el}: check failed; nothing written")

    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    print()
    print(f"{'el':4s}{'T':>6s}{'shells':>8s}{'segs':>6s}{'points':>8s}")
    print("-" * 32)
    for el, m in sorted(MODELS.items()):
        if el not in lib or not isinstance(lib[el], dict):
            print(f"{el:4s}  not in library.json, skipped")
            continue
        mass = refdata.MASSES[el]
        fc = m["fc"]() if "fc" in m else B.force_constants(m["struct"],
                                                           m["shells"])
        pts = CART[m["struct"]]
        segs = {}
        for a, z in DRAWN[m["struct"]]:
            qa, qz = np.array(pts[a]), np.array(pts[z])
            rows = []
            for k in range(NQ):
                t = k / (NQ - 1.0)
                f = B.freqs_THz(fc, qa + t * (qz - qa), mass)
                rows.append([round(t, 5)] + [round(float(x), 4) for x in f])
            segs[f"{a}|{z}"] = rows
        c11, c12, c44, cp = B.elastic_GPa(fc, float(
            refdata.ELEMENTS[el]["a0"]), mass, m["nat"])
        want = m["want"]
        if want is None:
            w = refdata.ELEMENTS[el]["Cij"]
            want = {"c11": w["C11"], "c44": w["C44"],
                    "cp": 0.5 * (w["C11"] - w["C12"])}
        lib[el]["model_curve"] = {
            "T_K": m["T_K"], "kind": "Born-von K&aacute;rm&aacute;n",
            "ref": m["ref"], "segs": segs, "shells": len(m["shells"]),
            "neighbours": len(fc), "why": m["why"], "gate": m["gate"],
            #  two decimals, because strontium's c-prime is 2.45 GPa and
            #  one decimal turns that into 2.4 against a published 2.5 -
            #  a disagreement created entirely by the rounding
            "check": {"c11": round(c11, 2), "c44": round(c44, 2),
                      "cp": round(cp, 2), "published_c11": want["c11"],
                      "published_c44": want["c44"],
                      "published_cp": round(want["cp"], 2),
                      "src": m["src"]},
        }
        n = sum(len(v) for v in segs.values())
        print(f"{el:4s}{m['T_K']:5d}K{len(m['shells']):8d}"
              f"{len(segs):6d}{n:8d}")

    fourier_curves(lib)
    as_curves(lib)
    dfm_curves(lib)

    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print()
    print(f"merged into {path}")


if __name__ == "__main__":
    main()
