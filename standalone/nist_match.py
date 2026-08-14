#!/usr/bin/env python3
"""
Which of our baseline potentials NIST's records actually describe.

The index that maps a potential file to a NIST entry works by filename, and
for the files that came from NIST that is exact.  For the rest - the ones
LAMMPS ships under its own names, and the MEAM pairs keyed on a shared library
file - the mapping was done by hand, by author and year.  Each alias was
checked by fetching its page, and that check was too weak: it confirmed the
page EXISTS, not that it describes the file we hold.

The strong check was available all along.  We compute surface energies for
every baseline with our own code, and NIST publishes surface energies for
every potential it hosts.  If the two describe the same potential they agree
to five or six figures - not approximately, exactly, because it is the same
arithmetic on the same tabulated functions.  If they do not, they are two
different potentials wearing the same author's name.

Run over the thirty-seven matched records the split is not a gradient, it is
two populations:

    twenty-eight agree to better than 1 per cent - median 0.089, worst 0.721
    nine disagree by 1 to 2828 per cent

Copper is the clean demonstration.  Mishin, u3 and smf7 agree to six figures;
Cu_zhou.eam.alloy is out by 6 per cent on every facet, and its NIST intrinsic
stacking fault is 22.4 mJ/m^2 where ours is 86.3.  Neither number is wrong -
they are answers about different potentials.  LAMMPS's Zhou files are built
from the analytic parameters in the 2004 paper, and none of the three
implementations NIST hosts reproduces them.

So the gate here is agreement, and anything that fails it is dropped from the
comparison rather than quietly averaged into it.  A tolerance of 1 per cent is
far looser than the agreeing population needs and far tighter than any
disagreeing one survives, so nothing sits on the boundary.

One of the nine is excluded for a different reason and is labelled as such:
ruthenium's A3 entry reports 0.10 J/m^2 for the basal plane where the
reference is 2.9, a record that does not describe a crystal rather than a file
we failed to identify.  Saying "no match" there would have blamed our own
potential file for someone else's broken relaxation.

    python nist_match.py            # the table
    import nist_match; nist_match.confirmed()
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refdata                                       # noqa: E402

TOL = 1.0                                   # per cent, on any facet
REF = {}                                    # element -> {facet: DFT gamma}
PROTO = {"fcc": "A1", "bcc": "A2", "hcp": "A3"}


#  the densest plane of each prototype, and the same plane in the reference
CLOSE = {"A1": "111", "A2": "110", "A3": "0001", "A3'": "0001"}
WILD = 5.0                                  # ratio to the reference


def _broken(gamma, pre, el):
    """is NIST's record for this prototype describing a crystal at all?

    Ruthenium is the case that forced this.  Its A3 entry reports 0.104 J/m^2
    for the basal plane and its A3' entry 2.838 at almost the same lattice
    constant - a factor of twenty-seven between two stackings that differ only
    in period.  Ours agrees with the second, so it is the first that is
    broken, and calling that a failed match would have said something false
    about our own file.  Titanium's A3' entry is negative outright and
    beryllium's A7 entry is all zeros; same fault, easier to see.

    The obvious test - hold the prototypes of one element against each other -
    is not safe, because it assumes the broken ones are in the minority AND
    that they fail downwards.  Sodium breaks both assumptions at once: three
    of its four prototypes sit at 40 J/m^2, a hundred times too HIGH, and the
    one good record was the one voted out.  So the anchor is external instead,
    the density-functional reference this project already compares surfaces
    against, and the tolerance is deliberately enormous: a factor of five is
    far beyond any disagreement a real potential produces - ours, the worst in
    the library, is 2.9 - and far inside every failure seen here, the nearest
    of which is off by twenty-seven.
    """
    f = CLOSE.get(pre)
    ref = REF.get(el) or {}
    if not f or f not in gamma:
        return False
    #  the reference writes hexagonal facets with signs
    r = ref.get(f) or ref.get(f.replace("10-10", "1010"))
    if not r or r <= 0:
        return False
    g = gamma[f]
    return bool(g <= 0 or g / r > WILD or r / g > WILD)


def _nist_gamma(rec, struct, el):
    """the surface energies of the prototype that IS this element's crystal"""
    pre = PROTO.get(struct)
    if not pre or not rec.get("surface"):
        return None, False
    for proto, d in rec["surface"].items():
        #  A3' is double-hcp and is a different crystal from A3
        if proto.split("--")[0] == pre:
            g = d.get("gamma") or {}
            return g, _broken(g, pre, el)
    return None, False


def compare(surface_path=None, props_path=None):
    """[(worst deviation %, element, filename, n facets, ok)]"""
    sp = surface_path or os.path.join(HERE, "..", "lammps", "surface.json")
    pp = props_path or os.path.join(HERE, "nist_props.json")
    ours = json.load(open(sp))
    nist = json.load(open(pp))
    if not REF:
        rp = os.path.join(HERE, "surface_ref.json")
        if os.path.exists(rp):
            for e, d in json.load(open(rp)).items():
                #  the facets are nested; reading the top level instead picks
                #  up the anisotropy and leaves every element with no facets
                #  at all, which makes the broken-record test silently never
                #  fire - the failure mode this whole module exists to catch
                f = d.get("facets") if isinstance(d, dict) else None
                if isinstance(f, dict):
                    REF[e] = {k: v for k, v in f.items()
                              if isinstance(v, (int, float))}
        if not REF:
            raise RuntimeError("surface_ref.json is unreadable - "
                               "the corrupted-record test will not run")

    rows = []
    for key, v in sorted(ours.items()):
        if "|base|" not in key:
            continue
        el, fn = key.split("|base|")
        rec = nist.get(f"{el}|{fn}")
        if not rec:
            continue
        nd, broken = _nist_gamma(
            rec, refdata.ELEMENTS.get(el, {}).get("struct"), el)
        if not nd:
            continue
        devs = []
        for facet, o in v.items():
            if not isinstance(o, dict) or "gamma" not in o:
                continue
            n = nd.get(facet)
            if n is None or n <= 0:
                continue
            devs.append(abs(o["gamma"] - n) / n * 100.0)
        if not devs:
            continue
        worst = max(devs)
        rows.append({"el": el, "file": fn, "worst_pct": worst,
                     "facets": len(devs), "nist_broken": broken,
                     "ok": bool(worst <= TOL and not broken)})
    rows.sort(key=lambda r: -r["worst_pct"])
    return rows


def confirmed(**kw):
    """{"El|file"} for the records whose NIST entry is the same potential"""
    return {f"{r['el']}|{r['file']}" for r in compare(**kw) if r["ok"]}


def main():
    rows = compare()
    print(f"{'dev %':>9s}  {'el':3s} {'file':34s}{'facet':>6s}  status")
    print("-" * 68)
    for r in rows:
        note = ("NIST KAYDI BOZUK - disarida" if r["nist_broken"] else
                "eslesme DOGRULANDI" if r["ok"] else
                "ESLESME YOK - disarida")
        print(f"{r['worst_pct']:9.3f}  {r['el']:3s} {r['file'][:34]:34s}"
              f"{r['facets']:6d}  {note}")
    ok = [r for r in rows if r["ok"]]
    nb = [r for r in rows if r["nist_broken"]]
    print(f"\n{len(ok)}/{len(rows)} eslesme dogrulandi "
          f"(esik %{TOL:g}, herhangi bir fasette)")
    print(f"  eslesmeyen {len(rows) - len(ok) - len(nb)}, "
          f"NIST record corrupt {len(nb)}"
          + (f" ({', '.join(r['el'] for r in nb)})" if nb else ""))
    if ok:
        w = sorted(r["worst_pct"] for r in ok)
        print(f"worst deviation among the confirmed: median {w[len(w)//2]:.3f} %, "
              f"en fazla %{w[-1]:.3f}")
    out = os.path.join(HERE, "nist_match.json")
    json.dump({r["el"] + "|" + r["file"]: r for r in rows},
              open(out, "w"), indent=1, sort_keys=True)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
