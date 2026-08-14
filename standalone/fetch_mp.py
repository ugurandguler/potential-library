#!/usr/bin/env python3
"""
Pull Materials Project reference data for the elements in library.json and merge
it in under a "mp" key.

What is fetched, per element (the ground-state entry matching our structure):
    mp-id, DFT lattice constant, energy above hull
    elastic tensor  -> C11, C12, C44 (+ C13, C33 for hcp), K and G
    phonon band structure, if MP has one, resampled onto OUR standard path

Why the elastic tensor matters: our fit targets experiment, and experiment for
elastic constants is room-temperature while DFT is 0 K.  Having both columns
side by side separates "the potential is wrong" from "the reference is a
different quantity".

The API key is read from MP_API_KEY or from mp_key.txt next to this file; it is
never written into the output.

    python fetch_mp.py            # everything in library.json
    python fetch_mp.py Pd Cu Mg
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

#  our structure -> the spacegroup MP should report for the elemental ground state
SG = {"fcc": 225, "bcc": 229, "hcp": 194}


def api_key():
    k = os.environ.get("MP_API_KEY")
    if k:
        return k.strip()
    p = os.path.join(HERE, "mp_key.txt")
    if os.path.exists(p):
        return open(p).read().strip()
    raise SystemExit("no API key: set MP_API_KEY or create mp_key.txt")


def pick(docs, struct):
    """the lowest-energy entry whose symmetry matches the structure we fitted"""
    want = SG[struct]
    same = [d for d in docs
            if getattr(getattr(d, "symmetry", None), "number", None) == want]
    pool = same or list(docs)
    if not pool:
        return None
    return min(pool, key=lambda d: (getattr(d, "energy_above_hull", None)
                                    if getattr(d, "energy_above_hull", None)
                                    is not None else 9e9))


#  MP gives the elastic tensor twice: `ieee_format`, rotated into the standard
#  orientation but ROUNDED TO INTEGERS, and `raw`, unrounded but in the
#  calculation's own frame.  For a stiff metal the rounding is nothing - copper
#  is 186 against 185.949 - but for the alkali metals it is the whole signal.
#  Rubidium's ieee tensor is 3, 3, 3 / 2, 2, 2, so C11 = C12 exactly, C' = 0 and
#  the crystal reads as mechanically unstable; unrounded it is 3.058, 2.663,
#  1.889, giving C' = 0.198 GPa against a measured 0.26.  The data was fine and
#  the rounding destroyed it.
#
#  So: average `raw` into the cubic or hexagonal form ourselves, and only fall
#  back to `ieee_format` when raw is missing.  The averaging is checked - the
#  components symmetry says are equal have to agree - because raw is not
#  guaranteed to be in the standard orientation, and averaging a rotated tensor
#  would quietly produce numbers that are not elastic constants at all.
RAW_SPREAD_MAX = 0.15       # relative disagreement allowed among equivalents
C_MAX_GPA = 5000.0          # diamond, the stiffest solid there is, has C11 1080


def _sym_average(C, struct):
    """(constants, worst relative spread among components symmetry equates)"""
    def avg(*v):
        v = [float(x) for x in v]
        m = sum(v) / len(v)
        sp = (max(v) - min(v)) / max(abs(m), 1e-9)
        return m, sp
    if struct == "hcp":
        c11, s1 = avg(C[0, 0], C[1, 1])
        c13, s2 = avg(C[0, 2], C[1, 2])
        c44, s3 = avg(C[3, 3], C[4, 4])
        out = {"C11": c11, "C12": float(C[0, 1]), "C13": c13,
               "C33": float(C[2, 2]), "C44": c44}
        return out, max(s1, s2, s3)
    c11, s1 = avg(C[0, 0], C[1, 1], C[2, 2])
    c12, s2 = avg(C[0, 1], C[0, 2], C[1, 2])
    c44, s3 = avg(C[3, 3], C[4, 4], C[5, 5])
    return {"C11": c11, "C12": c12, "C44": c44}, max(s1, s2, s3)


def tensor_problem(c, struct):
    """
    Why this tensor is unusable, or None if it describes a solid.

    Born stability rather than a list of bad material ids, because the failures
    are of different kinds and a list would only catch the ones already seen.
    Potassium's raw tensor has one diverged strain, C11 = 599445 GPa against
    C22 = C33 = 4.12, which the IEEE symmetrisation then spreads over all three
    diagonal entries.  Sodium's is not broken at all - MP simply finds bcc
    sodium mechanically unstable at 0 K, with C' = -0.08 GPa, which is a real
    result about the phase rather than a bad calculation and is worth saying so
    on the page rather than hiding.
    """
    big = max(abs(v) for v in c.values())
    if not all(v == v for v in c.values()):
        return "the tensor contains non-finite entries"
    if big > C_MAX_GPA:
        return (f"largest component {big:.0f} GPa, above {C_MAX_GPA:.0f}; "
                f"no elemental solid is that stiff")
    if struct == "hcp":
        cp = 0.5 * (c["C11"] - c["C12"])
        ok = (c["C44"] > 0 and cp > 0 and c["C33"] > 0 and
              (c["C11"] + c["C12"]) * c["C33"] > 2 * c["C13"] ** 2)
    else:
        cp = 0.5 * (c["C11"] - c["C12"])
        ok = c["C44"] > 0 and cp > 0 and (c["C11"] + 2 * c["C12"]) > 0
    if not ok:
        return (f"not Born stable: C' = {cp:.2f}, C44 = {c['C44']:.2f} GPa "
                f"- MP's own tensor says this phase does not stand up")
    return None


def elastic_from(doc, struct):
    """MP elasticity doc -> the constants we compare against"""
    t = getattr(doc, "elastic_tensor", None)
    if t is None:
        return None
    src, C = "raw", getattr(t, "raw", None)
    if C is None:
        src, C = "ieee_format", getattr(t, "ieee_format", None)
    if C is None:
        return None
    out, spread = _sym_average(np.array(C, dtype=float), struct)
    if src == "raw" and spread > RAW_SPREAD_MAX:
        #  not in the standard orientation, so averaging it is meaningless;
        #  the rounded tensor is still better than a wrong one
        ie = getattr(t, "ieee_format", None)
        if ie is None:
            return {"rejected": f"raw tensor components that symmetry equates "
                                f"disagree by {spread*100:.0f} %"}
        src = "ieee_format"
        out, spread = _sym_average(np.array(ie, dtype=float), struct)
    bad = tensor_problem(out, struct)
    if bad:
        return {"rejected": bad}
    out = {k: round(v, 2) for k, v in out.items()}
    out["tensor_source"] = src
    for name, attr in (("K", "k_vrh"), ("G", "g_vrh")):
        v = getattr(doc, attr, None)
        if v is None and isinstance(getattr(doc, "bulk_modulus", None), dict):
            v = doc.bulk_modulus.get("vrh")
        if v is not None:
            out[name] = float(v)
    return out


def main(only):
    lib = json.load(open(os.path.join(HERE, "library.json")))
    els = [e for e in sorted(lib) if not only or e in only]
    from mp_api.client import MPRester

    got = {"entry": 0, "elastic": 0, "phonon": 0}
    with MPRester(api_key()) as mpr:
        for i, el in enumerate(els, 1):
            struct = lib[el]["struct"]
            rec = {}
            try:
                docs = mpr.materials.summary.search(
                    chemsys=el, num_elements=1,
                    fields=["material_id", "symmetry", "energy_above_hull",
                            "structure", "formula_pretty"])
            except Exception as exc:                      # noqa: BLE001
                print(f"  {el:3s} summary query failed: {exc}")
                continue
            doc = pick(docs, struct)
            if doc is None:
                print(f"  {el:3s} no entry")
                continue
            mid = str(doc.material_id)
            rec["mp_id"] = mid
            rec["e_above_hull"] = (float(doc.energy_above_hull)
                                   if doc.energy_above_hull is not None else None)
            sg = getattr(getattr(doc, "symmetry", None), "number", None)
            rec["spacegroup"] = sg
            rec["matches_structure"] = (sg == SG[struct])
            st = getattr(doc, "structure", None)
            if st is not None:
                lat = st.lattice
                rec["a_dft"] = float(lat.abc[0])
                rec["volume_per_atom"] = float(lat.volume/len(st))
            got["entry"] += 1

            try:
                ed = mpr.materials.elasticity.search(
                    material_ids=[mid],
                    fields=["elastic_tensor", "bulk_modulus", "shear_modulus"])
                if ed:
                    e = elastic_from(ed[0], struct)
                    if e and "rejected" in e:
                        rec["elastic_rejected"] = e["rejected"]
                        e = None
                    if e:
                        rec["elastic"] = e
                        got["elastic"] += 1
            except Exception as exc:                      # noqa: BLE001
                rec["elastic_error"] = str(exc)[:120]

            #  ---- phonon band structure ----
            #  Two traps here.  The convenience accessor
            #  MPRester.get_phonon_bandstructure_by_material_id() hard-codes
            #  phonon_method="dfpt", but MP's current elemental phonons were all
            #  produced with "pheasy", so it reports "no data" for entries that
            #  plainly have some.  Go through the route directly and ask for the
            #  method the metadata says was used.
            #
            #  We store MP's own q-points, not just the frequencies: our
            #  dispersion is then evaluated at exactly those q, so the two curves
            #  share an x-axis point for point instead of being resampled onto
            #  each other.  MP's labelled points are the Setyawan-Curtarolo ones
            #  (verified: fcc X=(1/2,0,1/2), W=(1/2,1/4,3/4), K=(3/8,3/8,3/4)),
            #  which is the convention this library already uses.
            try:
                pdocs = mpr.materials.phonon.search(
                    material_ids=[mid], fields=["phonon_method", "symmetry"])
                if pdocs:
                    meth = getattr(pdocs[0], "phonon_method", None) or "pheasy"
                    meth = getattr(meth, "value", None) or str(meth)
                    bs = mpr.materials.phonon.get_bandstructure_from_material_id(
                        mid, phonon_method=meth)
                    fr = np.array(bs.frequencies, dtype=float)          # THz
                    qq = np.array([getattr(p, "frac_coords", p)
                                   for p in bs.qpoints], dtype=float)
                    step = max(1, len(qq)//360)                          # thin out
                    idx = list(range(0, len(qq), step))
                    if idx[-1] != len(qq)-1:
                        idx.append(len(qq)-1)
                    labels = {}
                    for k, v in bs.labels_dict.items():
                        name = "G" if "Gamma" in k else k
                        labels[name] = [float(x) for x in np.array(v, dtype=float)]
                    marks = []
                    for j, i0 in enumerate(idx):
                        for name, c in labels.items():
                            if np.allclose(qq[i0], c, atol=2e-3):
                                if not marks or marks[-1][1] != name:
                                    marks.append([j, name])
                                break
                    rec["phonon"] = {
                        "method": meth,
                        "q": [[round(float(x), 6) for x in qq[i]] for i in idx],
                        "f": [[round(float(fr[b, i])*33.35641, 1) for i in idx]
                              for b in range(fr.shape[0])],   # THz -> cm^-1
                        "marks": marks,
                        "kpath": [("G" if "Gamma" in s else s) for s in
                                  (bs.kpath or [])],
                    }
                    got["phonon"] += 1
            except Exception as exc:                      # noqa: BLE001
                rec["phonon_error"] = f"{type(exc).__name__}: {str(exc)[:110]}"

            lib[el]["mp"] = rec
            tag = (f"{mid:12s} sg={sg}"
                   f"{'' if rec['matches_structure'] else ' (DIFFERENT STRUCTURE)'}"
                   f"  a={rec.get('a_dft', float('nan')):.4f}"
                   f"  elastic={'yes' if 'elastic' in rec else 'no'}"
                   f"  phonon={'yes' if 'phonon' in rec else 'no'}")
            print(f"  [{i:2d}/{len(els)}] {el:3s} {tag}", flush=True)

    dst = os.path.join(HERE, "library.json")
    tmp = dst + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, dst)
    print(f"\nentries {got['entry']}/{len(els)}, "
          f"elastic tensors {got['elastic']}, phonon band structures {got['phonon']}")


if __name__ == "__main__":
    main(set(sys.argv[1:]))
