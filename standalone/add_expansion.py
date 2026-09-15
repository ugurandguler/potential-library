#!/usr/bin/env python3
"""
Fold the thermal expansion into the library.

This is the one test on the list that lives INSIDE the regime the form
describes.  The vacancy, the surface, the ground-state ranking and the
stacking fault all ask what happens when an atom's coordination changes.
Thermal expansion does not: the atoms keep every neighbour and simply sit
further apart, so what is measured is the anharmonicity of the same bonds
whose curvature was fitted - the third derivative of a curve whose second
derivative is a target.

The number comes from NPT dynamics (`lammps/npt_expansion.py`), not from the
quasi-harmonic calculation that preceded it.  That one could not be made to
converge: the shift being looked for is two parts in a thousand of the lattice
parameter over four hundred kelvin, and the answer moved with the fitting
protocol.  More decisively, it had no baseline - `latdyn` only knows this
potential's own functional form, so there was nothing to check it against.  A
barostat at zero pressure has one: every published EAM and MEAM potential runs
through the identical code.

The Grueneisen parameter is kept from the harmonic spectra, because it is what
keeps the expansion honest.  Thermal expansion exists only because the
frequencies shift with volume; an alpha that happens to look reasonable with
the wrong gamma is two errors cancelling.  A metal sits between about 1 and 3.

    python add_expansion.py
    python add_expansion.py --dry
"""
import json
import statistics
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refdata                                        # noqa: E402
import expansion as X                                # noqa: E402

#  An arm record is one that carries the potential's five parameters.  Nothing
#  else in an element's record does - baseline_*, elasticT, dyn and mech carry
#  none of them - so this separates arms from payloads without a name list.
#  A NAME LIST IS THE BUG: every one of these files was written with
#  ("tap", "tap_ug"), the re-cut arms arrived and were silently skipped, and
#  `tap_force` arrived and was skipped again.  The next arm will not announce
#  itself either.
_PARAM = ("D", "r0", "alpha", "m", "gamma")


def arms(d):
    """the arm records inside one element's entry, by shape not by name"""
    return [k for k, v in d.items()
            if isinstance(v, dict) and all(p in v for p in _PARAM)]



SRC = os.path.join(HERE, "..", "lammps", "npt_expansion.json")
#  the harmonic Grueneisen still comes from the quasi-harmonic run; only the
#  expansion coefficient itself was replaced
QHA = {"tap": "expansion_tap.json", "tap_ug": "expansion_tap_ug.json"}


def main():
    dry = "--dry" in sys.argv
    lib = json.load(open(os.path.join(HERE, "library.json")))
    if not os.path.exists(SRC):
        print(f"{SRC} is missing")
        return
    npt = json.load(open(SRC))

    gam = {}
    for tag, fn in QHA.items():
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            for el, r in json.load(open(p)).items():
                if r.get("gruneisen") is not None:
                    gam[f"{el}|{tag}"] = {
                        k: r[k] for k in ("gruneisen", "gruneisen_weighted")
                        if k in r}

    #  Clear first, then write.  Merging into whatever is already there keeps
    #  records from runs that no longer exist: lithium's quasi-harmonic entry
    #  survived the switch to NPT and sat in the library as a NaN, rendered on
    #  the page as though it were a measurement.  A record should exist only
    #  if this run produced it.
    cleared = 0
    for el, d in lib.items():
        if not isinstance(d, dict):
            continue
        #  Every sub-record that carries one, not a fixed list of arm names.
        #  The list said ("tap", "tap_ug") and the re-ranked candidates arrived
        #  under "rc" and "rc_ug", so their results would have survived a
        #  rerun that no longer produced them - which is the exact failure the
        #  note above describes, one arm later.
        for t, sub in list(d.items()):
            if isinstance(sub, dict) and "expansion" in sub:
                del sub["expansion"]
                cleared += 1
        d.pop("baseline_expansion", None)

    n = 0
    ours, base = [], []
    for key, v in sorted(npt.items()):
        el, tag = key.split("|", 1)
        rec = {k: v[k] for k in ("T", "a", "alpha_1e6", "alpha_exp_1e6",
                                 "ratio", "Pavg", "dropped", "failed",
                                 "points", "of") if k in v}
        #  The measured slope spans the whole grid; the experimental value is
        #  at 25 C.  alpha follows C_V, so the two are not the same quantity
        #  and the ratio carried a factor <C_V>_grid / C_V(298) belonging to
        #  neither the potential nor the experiment.  See expansion.py's
        #  window_factor for why, and for the test that says it is a
        #  conversion rather than a fudge - it improves the published EAM and
        #  MEAM baselines MORE than it improves ours.
        #
        #  The raw ratio is kept beside it: this changes a published number,
        #  and a reader should be able to see what it was.
        w = X.window_factor(el, v.get("T") or [])
        if w and v.get("ratio") and not v.get("failed"):
            rec["window"] = round(w, 4)
            rec["ratio_raw"] = v["ratio"]
            rec["ratio"] = v["ratio"] / w
            v = dict(v, ratio=rec["ratio"])
        #  a run that died carries its reason and no numbers
        if v.get("failed"):
            rec.pop("alpha_1e6", None)
            rec.pop("ratio", None)
        if tag.startswith("base|"):
            lib.setdefault(el, {}).setdefault("baseline_expansion",
                                              {})[tag[5:]] = rec
            if v.get("ratio"):
                base.append((v["ratio"], el, tag[5:]))
        else:
            if el not in lib or not isinstance(lib[el].get(tag), dict):
                continue
            rec.update(gam.get(key, {}))
            lib[el][tag]["expansion"] = rec
            if v.get("ratio"):
                ours.append((v["ratio"], el, tag))
        n += 1

    print(f"{n} records ({len(ours)} ours, {len(base)} baseline), "
          f"{cleared} stale records cleared")

    #  A non-finite number is not a measurement.  Say so loudly rather than
    #  writing it and letting toFixed render it as "NaN" on the page.
    bad = [k for k, v in npt.items() if not v.get("failed")
           and (v.get("alpha_1e6") is None
                or v["alpha_1e6"] != v["alpha_1e6"])]
    if bad:
        print(f"  SAYI OLMAYAN {len(bad)}: {bad[:6]}")
    died = [(k, v["failed"]) for k, v in sorted(npt.items()) if v.get("failed")]
    if died:
        print(f"  KOSU COKEN {len(died)}:")
        for k, why in died:
            print(f"     {k:14s} {why[:58]}")

    #  and the gap that matters: parameters exist, no result came back
    #  the same question, over whichever arms this element actually has
    #  TWO QUESTIONS, AND THE OLD LIST CONFLATED THEM.  "Parameters exist, no
    #  result came back" is one thing for an arm this sweep covers - that is a
    #  run that failed - and quite another for an arm the sweep has never
    #  touched, which is simply not measured yet.  Reported together over every
    #  arm in the library the second buries the first: the hard-cut `ug` record
    #  alone contributes one line per element, 38 of them, and none is news.
    #
    #  So the covered set is read from the sweep's OWN OUTPUT - the arms that
    #  appear in npt for at least one element - and everything else is counted
    #  on one line instead of listed.  No name list either way.
    covered = {k.split("|", 1)[1] for k in npt}
    gone = sorted(f"{el}|{t}" for el, d in lib.items()
                  if isinstance(d, dict)
                  for t in arms(d)
                  if t in covered
                  and "expansion" not in d[t] and f"{el}|{t}" not in npt)
    never = sorted({t for d in lib.values() if isinstance(d, dict)
                    for t in arms(d) if t not in covered})
    if never:
        print(f"  bu taramanin hic olcmedigi kol: {' '.join(never)}")
    if gone:
        print(f"  SONUC DONMEYEN {len(gone)}: {' '.join(gone)}")
        print("   (the run collapsed or the guards rejected it - see the npt output)")

    def report(name, rows):
        if not rows:
            return None
        r = sorted(x[0] for x in rows)
        med = statistics.median(r)
        neg = [x for x in rows if x[0] < 0]
        lo = [x for x in rows if 0 <= x[0] < 0.75]
        hi = [x for x in rows if x[0] > 1.25]
        print(f"\n{name}: ratio to experiment, median {med:.2f}, "
              f"range {r[0]:.2f}..{r[-1]:.2f}")
        print(f"   NEGATIVE expansion {len(neg)}/{len(rows)}"
              + (f"  {[f'{e}|{t}' for _, e, t in neg]}" if neg else ""))
        print(f"   %25'ten dusuk {len(lo)}, %25'ten yuksek {len(hi)}")
        return med

    m_ours = report("BIZIM", ours)
    m_base = report("TABAN", base)
    #  the baselines are the scale.  Without them "median 0.68" is a number
    #  with nothing to be 0.68 of - the same argument the surface and stacking
    #  fault sections rest on.
    if m_ours and m_base:
        print(f"\nscale: published potentials through the same code {m_base:.2f}, "
              f"biz {m_ours:.2f}")

    #  which structures the negatives belong to: a failure that sorts by
    #  crystal is a different statement from one scattered at random
    neg = [(e, t) for r, e, t in ours if r < 0]
    if neg:
        by = {}
        for e, t in neg:
            by.setdefault(refdata.ELEMENTS.get(e, {}).get("struct", "?"),
                          []).append(e)
        print("   negatiflerin yapisi: "
              + ", ".join(f"{k} {sorted(set(v))}" for k, v in sorted(by.items())))

    #  every arm with an expansion record, keyed on the payload - the two-name
    #  list here omitted the re-cut arms whose expansion runs this same file
    #  had just written, so the Gruneisen spread was quoted over a subset
    gs = sorted(v["gruneisen"] for el in lib if isinstance(lib[el], dict)
                for k, sub in lib[el].items()
                if isinstance(sub, dict) and not k.startswith("baseline")
                for v in [sub.get("expansion") or {}]
                if isinstance(v, dict) and v.get("gruneisen") is not None)
    if gs:
        odd = [g for g in gs if not 0.5 < g < 4.0]
        print(f"\nGruneisen: median {statistics.median(gs):.2f}, "
              f"range {gs[0]:.2f}..{gs[-1]:.2f}; "
              f"fiziksel araligin disinda {len(odd)}/{len(gs)}")

    if dry:
        print("\n--dry: nothing written")
        return
    tmp = os.path.join(HERE, "library.json.tmp")
    json.dump(lib, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(HERE, "library.json"))
    print("\n-> library.json")


if __name__ == "__main__":
    main()
