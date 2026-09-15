#!/usr/bin/env python3
"""
Fold the surface energies into the library, with the two things worth reading
off them computed here rather than in a template literal.

**Magnitude** is the obvious one and the weaker one.  Any form can be rescaled
to hit a number, and the prediction written down before the calculation - that
these would come out two to three times high, for the same reason the vacancy
does - was directionally right and quantitatively too optimistic: the measured
factor is nearer three and a half, and aluminium is worse than that.

**Ordering** is the one that cannot be rescaled away.  The usual rule is that
the close-packed face is the cheapest - gamma(111) < gamma(100) < gamma(110)
for fcc - but the reference contradicts that rule for 17 of the 35 elements it
covers, so each record is scored against its OWN element's reference ordering,
and orderings the reference resolves by less than ORDER_MARGIN are not scored
at all.  On the records that can be judged, published potentials get 22 of 32
and ours get 3 of 38.  A form that inverts a resolvable ordering is not
describing a surface, whatever its magnitude.

The comparison is only worth anything because the published potentials went
through the identical slab thickness, vacuum, bulk reference and relaxation.
Five copper potentials come back at 0.82-1.08 of the DFT value with the correct
ordering, and Mishin 2001 reproduces its own published 1.24 J/m^2 for Cu(111).
Whatever the numbers below say, they are not saying it about the method.

    python add_surface.py
    python add_surface.py --dry
"""
import json
import statistics
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

import refdata      # noqa: E402

#  the ordering a real metal has, per structure, cheapest first
RIGHT = {"fcc": ["111", "100", "110"],
         "bcc": ["110", "100", "111"],
         "hcp": ["0001", "10-10", "11-20"]}

#  reference facets closer together than this are treated as a tie
ORDER_MARGIN = 0.05


def summarise(struct, gam, ref):
    """what to say about one record's three facets"""
    fs = {f: v for f, v in gam.items() if v is not None}
    if len(fs) < 2:
        return None
    order = sorted(fs, key=fs.get)
    rf = (ref or {}).get("facets") or {}
    #  RIGHT is a rule of thumb and the reference calculation contradicts it
    #  for 17 of the 35 elements it covers: close packing does not always
    #  win.  Where the reference has every facet the record does, ITS order
    #  is the target and the rule is only the fallback.
    #
    #  But the target has to carry its own margin.  In 19 of those 35 the two
    #  closest reference facets are within 5 % of each other - Ca 0.6 %,
    #  Zr 0.7 %, Rh 0.8 % - and demanding that a potential resolve a gap that
    #  narrow asks for an accuracy neither it nor the reference has.  Those
    #  orderings are recorded as unresolved rather than scored, so a record is
    #  neither credited nor blamed for a coin toss.
    ref_ord = [f for f in sorted(rf, key=rf.get) if f in fs]
    margin = None
    if len(ref_ord) == len(fs) and len(ref_ord) > 1:
        want, basis = ref_ord, "reference"
        v = [rf[f] for f in ref_ord]
        margin = min((v[i + 1] - v[i]) / v[i] for i in range(len(v) - 1))
    else:
        want, basis = [f for f in RIGHT[struct] if f in fs], "rule"
    resolved = margin is None or margin >= ORDER_MARGIN
    out = {"gamma": fs, "order": order,
           "order_ok": bool(order == want) if resolved else None,
           "order_want": want, "order_basis": basis,
           "order_margin": margin, "order_resolved": bool(resolved),
           #  spread as a fraction of the mean: a form with no facet
           #  anisotropy at all is a different failure from one with the
           #  wrong anisotropy, and the number separates them
           "spread": (max(fs.values()) - min(fs.values()))
           / (sum(fs.values()) / len(fs))}
    rat = {f: fs[f] / rf[f] for f in fs if rf.get(f)}
    if rat:
        out["ratio_dft"] = rat
        out["ratio_dft_median"] = statistics.median(rat.values())
    if (ref or {}).get("tyson"):
        out["ratio_exp_median"] = (statistics.median(fs.values())
                                   / ref["tyson"])
    if rf:
        ro = sorted(rf, key=rf.get)
        out["ref_order"] = [f for f in ro if f in fs][:3]
    return out


def main():
    dry = "--dry" in sys.argv
    lib = json.load(open(os.path.join(HERE, "library.json")))

    def load(p):
        return json.load(open(p)) if os.path.exists(p) else {}

    surf = load(os.path.join(ROOT, "lammps", "surface.json"))
    ref = load(os.path.join(HERE, "surface_ref.json"))
    labels = load(os.path.join(HERE, "baseline_labels.json"))

    nours = nbase = nbad = ntie = 0
    for key, fs in surf.items():
        el, tag = key.split("|", 1)
        if el not in lib:
            continue
        struct = refdata.ELEMENTS[el]["struct"]
        gam = {f: (r.get("gamma") if isinstance(r, dict) and "error" not in r
                   else None) for f, r in fs.items()}
        s = summarise(struct, gam, ref.get(el))
        if not s:
            continue
        s["unrelaxed"] = {f: r.get("gamma_unrelaxed")
                          for f, r in fs.items() if "error" not in r}
        if tag.startswith("base|"):
            fn = tag.split("|", 1)[1]
            s["label"] = (labels.get(fn) or {}).get("label", fn)
            s["file"] = fn
            lib[el].setdefault("baseline_surface", {})[fn] = s
            nbase += 1
        elif tag in lib[el] and isinstance(lib[el][tag], dict):
            lib[el][tag]["surface"] = s
            nours += 1
            if s["order_ok"] is None:
                ntie += 1
            elif not s["order_ok"]:
                nbad += 1
        else:
            continue

    #  NIST's own calculation of the same quantity for the same potential,
    #  attached to the baseline it belongs to.  This is not a property of the
    #  element - it is evidence that the slab machinery is right - so it goes
    #  beside the baseline rather than into a section of its own, and the page
    #  states it once per element rather than arguing it again each time.
    nist = load(os.path.join(HERE, "nist_props.json"))
    PROTO = {"fcc": "A1--Cu--fcc", "bcc": "A2--W--bcc", "hcp": "A3--Mg--hcp"}
    nn = 0
    for key, r in nist.items():
        el, fn = key.split("|", 1)
        b = ((lib.get(el) or {}).get("baseline_surface") or {}).get(fn)
        if not b or not r.get("surface"):
            continue
        #  NOT named `ref`: that is the surface-reference dictionary in the
        #  enclosing scope, and shadowing it here emptied every element's
        #  reference without any error being raised.  The script's own count
        #  line said "referans 0 element" and that is the only reason it was
        #  noticed.
        nref = r["surface"].get(PROTO[refdata.ELEMENTS[el]["struct"]])
        if not nref:
            continue
        a0 = refdata.ELEMENTS[el]["a0"]
        #  same crystal, or the comparison is between two different metals
        if abs(nref["a"] - a0) / a0 > 0.05:
            b["nist"] = {"skipped": "different lattice constant",
                         "a_nist": nref["a"], "a_ours": a0}
            continue
        #  Is the NIST record consistent with ITSELF?  The close-packed face
        #  of fcc, hcp and double-hcp differ only in how the layers stack, so
        #  their surface energies must be close.  Ruthenium's record gives
        #  0.1036 J/m2 for both the fcc (111) and the hcp (0001) - the same
        #  number to four figures, which is not physics - while its own
        #  double-hcp entry gives 2.8376, next to our 3.0345 and an
        #  experimental 3.05.  Comparing against the broken half of that
        #  record produced a "2828 % disagreement" on the ruthenium page,
        #  which read as our calculation failing when the opposite is true.
        CP = {"A1--Cu--fcc": "111", "A3--Mg--hcp": "0001",
              "A3'--alpha-La--double-hcp": "0001"}
        cp = [pr["gamma"][f] for k, f in CP.items()
              if (pr := r["surface"].get(k)) and f in pr["gamma"]
              and pr["gamma"][f] > 0]
        if len(cp) > 1 and max(cp) / min(cp) > 1.5:
            b["nist"] = {"unreliable": ("the NIST record is internally inconsistent: "
                                        "close-packed surface values "
                                        + " / ".join(f"{v:.4f}" for v in cp)),
                         "close_packed": cp, "impl": r.get("impl")}
            continue
        d = {f: nref["gamma"][f] for f in b["gamma"] if f in nref["gamma"]}
        if not d:
            continue
        err = {f: 100.0 * (b["gamma"][f] - v) / v for f, v in d.items()}
        b["nist"] = {"gamma": d, "err_pct": err,
                     "worst_pct": max(abs(x) for x in err.values()),
                     "impl": r.get("impl")}
        nn += 1
    print(f"baselines linked to the NIST comparison: {nn}")

    #  attach a reference that has EITHER facet-resolved DFT or the
    #  experimental value.  Requiring facets dropped cobalt's and rhenium's
    #  measured 2.52 and 3.63 J/m2 from the page, which is the second time
    #  "no DFT" has been allowed to mean "nothing to compare against".
    nref = 0
    for el, r in ref.items():
        if el in lib and (r.get("facets") or r.get("tyson")):
            lib[el]["surface_ref"] = r
            nref += 1

    print(f"ours {nours} records ({nbad} with the facet ordering wrong, "
          f"{ntie} where the reference facets are too close to call), "
          f"baseline {nbase} records, reference {nref} elements")
    #  EVERY arm that produced a surface record, not a list of two names.
    #  Written as ("tap", "tap_ug") this summary silently omitted rc, rc_ug and
    #  tap_force - so the "ratio to DFT" line quoted below has been describing
    #  half the library while reading as though it described all of it.
    #  Keyed on the payload, so an arm with no surface run simply does not
    #  appear, which is what a statistic over results should do.
    ours = [sub["surface"] for el in lib if isinstance(lib[el], dict)
            for k, sub in lib[el].items()
            if isinstance(sub, dict) and not k.startswith("baseline")
            and isinstance(sub.get("surface"), dict)]
    rr = [s["ratio_dft_median"] for s in ours if s.get("ratio_dft_median")]
    if rr:
        rr.sort()
        print(f"ratio to DFT, our records: median {statistics.median(rr):.2f}, "
              f"range {rr[0]:.2f}-{rr[-1]:.2f}, {len(rr)} records")
    bs = [s for el in lib for s in (lib[el].get("baseline_surface") or {}).values()]
    rb = [s["ratio_dft_median"] for s in bs if s.get("ratio_dft_median")]
    if rb:
        rb.sort()
        nb_ok = sum(1 for s in bs if s.get("order_ok"))
        nb_tie = sum(1 for s in bs if s.get("order_ok") is None)
        print(f"ratio to DFT, published baselines: median {statistics.median(rb):.2f}, "
              f"range {rb[0]:.2f}-{rb[-1]:.2f}, {len(rb)} records; "
              f"siralamasi dogru olan {nb_ok}/{len(bs)}"
              f" ({nb_tie} ayirt edilemez)")
    if dry:
        print("--dry: nothing written")
        return
    tmp = os.path.join(HERE, "library.json.tmp")
    json.dump(lib, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(HERE, "library.json"))
    print("-> library.json")


if __name__ == "__main__":
    main()
