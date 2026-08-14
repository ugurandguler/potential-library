#!/usr/bin/env python3
"""
NIST's own computed properties for the published potentials we compare against.

The NIST Interatomic Potentials Repository does not only host potential files -
it runs a fixed battery of calculations on every potential it hosts and
publishes the results.  That makes it something this project has not had
before: an INDEPENDENT calculation of the same quantity, for the same
potential, by different code.

The reason to want it is a specific claim.  Our surface energies say that
seventy-two of seventy-six records get the facet ordering wrong while
forty-eight of fifty-one published potentials get it right, and the whole
weight of that rests on our slab machinery being correct.  So far the argument
for that has been internal: the published potentials land where a published
potential should.  NIST's numbers close it from outside.  For Mishin's 2001
copper the two agree to four significant figures -

    facet     ours      NIST
    (111)   1.2395    1.2395
    (100)   1.3453    1.3453
    (110)   1.4755    1.4754

- and this fetches the rest so that one agreement becomes fifty.

Two traps, both hit while working this out.

The per-potential CSVs carry SEVERAL crystal prototypes.  Keying the surface
energies by facet alone lets the A15 and double-hcp rows overwrite the fcc
ones, which turned an exact agreement into apparent errors of 1, 5 and 12 per
cent.  Filter on the prototype.

The thermal expansion column is VOLUMETRIC.  Compared directly with the
tabulated linear coefficient it makes a good EAM potential look 2.9 times too
large; divided by three, Mishin copper gives 15.7 against an experimental 16.5.
Verified here against the volume column rather than assumed.

    python fetch_nist_props.py
    python fetch_nist_props.py --only Cu,Al
"""
import ast
import csv
import io
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
MIRROR = ("https://raw.githubusercontent.com/lmhale99/potentials-library/"
          "master/potential_LAMMPS.csv")
ENTRY = "https://www.ctcms.nist.gov/potentials/entry"
UA = {"User-Agent": "ugurpotential-baseline-check/1.0"}


def get(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return fh.read().decode("utf-8", "replace")


#  Matching by filename only works for the files that came FROM NIST.  Those
#  that ship with LAMMPS carry LAMMPS's names, and the MEAM pairs are keyed on
#  a library file shared by many potentials, so neither is in the index under
#  the name we hold.  These are matched by hand, by author and year, and every
#  one was confirmed by fetching its page before being written down - an alias
#  that quietly points at the wrong potential would be worse than no alias.
ALIAS = {
    "Cu_mishin1.eam.alloy": "2001--Mishin-Y--Cu-1--LAMMPS--ipr1",
    "Cu_zhou.eam.alloy": "2004--Zhou-X-W--Cu--LAMMPS--ipr2",
    "Al_zhou.eam.alloy": "2004--Zhou-X-W--Al--LAMMPS--ipr2",
    "W_zhou.eam.alloy": "2004--Zhou-X-W--W--LAMMPS--ipr2",
    "Al_mm.eam.fs": "2008--Mendelev-M-I--Al--LAMMPS--ipr1",
}


def index():
    """filename -> (implementation id, potential id, elements)"""
    p = os.path.join(HERE, "nist_lammps_index.csv")
    if not os.path.exists(p):
        io.open(p, "w", encoding="utf-8").write(get(MIRROR))
    out = {}
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        if r.get("status") != "active":
            continue
        try:
            arts = ast.literal_eval(r["artifacts"]) or []
            els = ast.literal_eval(r["elements"]) or []
        except Exception:
            continue
        for a in arts:
            fn = a.get("filename")
            if fn:
                out.setdefault(fn, []).append((r["id"], r["potid"], els))
    return out


def props(impl, potid, el):
    """the computed properties NIST publishes for one implementation"""
    base = f"{ENTRY}/{potid}/{impl}"
    out = {"impl": impl, "potid": potid}

    #  free surfaces.  The prototype column matters: the same facet label
    #  appears once per crystal structure and keying on the label alone lets
    #  the wrong structure overwrite the right one.
    try:
        rows = list(csv.DictReader(io.StringIO(
            get(f"{base}/freesurface.{el}.csv"))))
        #  the lattice parameter travels with the facets.  NIST computes
        #  surfaces for every prototype it can relax, including ones that
        #  relaxed somewhere absurd - ruthenium's A3 entry sits at a lattice
        #  constant that gives 0.10 J/m2 for the basal plane against an
        #  experimental 3.05 - so a comparison has to check that both sides
        #  are describing the same crystal before it means anything.
        proto = {}
        for r in rows:
            d = proto.setdefault(r["prototype"], {"a": float(r["a"]),
                                                  "gamma": {}})
            d["gamma"][r["surface"]] = float(r["gamma_fs"]) / 1000.0
        out["surface"] = proto
    except Exception as ex:
        out["surface_error"] = f"{type(ex).__name__}"

    #  stacking faults, which this project has predicted but never measured
    try:
        rows = list(csv.DictReader(io.StringIO(
            get(f"{base}/stackingfault.{el}.csv"))))
        out["stacking"] = [{k: r.get(k) for k in
                            ("prototype", "stackingfault_id",
                             "E_isf", "E_usf")} for r in rows]
    except Exception:
        pass

    #  thermal expansion.  The published column is VOLUMETRIC; the linear
    #  coefficient is a third of it, and that was checked against the volume
    #  column rather than taken on trust.
    for name, key in (("alpha", "alpha_V"), ("V", "volume")):
        try:
            rows = list(csv.DictReader(io.StringIO(
                get(f"{base}/phonon.{el}.{name}.csv"))))
            col = next(c for c in rows[0]
                       if c and c not in ("", "temperature"))
            d = {}
            for r in rows:
                try:
                    d[float(r["temperature"])] = float(r[col])
                except (ValueError, TypeError):
                    pass
            if d:
                out[key] = {"column": col,
                            "T": sorted(d), "y": [d[t] for t in sorted(d)]}
        except Exception:
            pass
    return out


def main():
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))

    base = json.load(open(os.path.join(HERE, "baselines.json")))
    idx = index()
    print(f"NIST index: {len(idx)} file names")

    jobs, missing = [], []
    for el, lst in sorted(base.items()):
        if only and el not in only:
            continue
        for fn, _ in lst:
            first = fn.split("+")[0]
            hits = idx.get(first)
            #  the MEAM pairs are keyed on the element file, not the shared
            #  library one, so try the second half too
            if not hits and "+" in fn:
                hits = idx.get(fn.split("+")[1])
            if not hits and first in ALIAS:
                want = ALIAS[first]
                hits = [h for hs in idx.values() for h in hs if h[0] == want]
            if not hits:
                missing.append((el, fn))
                continue
            #  prefer an implementation whose element list contains ours
            hit = next((h for h in hits if el in h[2]), hits[0])
            jobs.append((el, fn, hit[0], hit[1]))

    print(f"eslesen: {len(jobs)}   eslesmeyen: {len(missing)}")
    if missing:
        print("  eslesmeyenler:", [f"{e}/{f}" for e, f in missing][:8])

    out = {}
    p = os.path.join(HERE, "nist_props.json")
    if os.path.exists(p):
        try:
            out = json.load(open(p))
        except Exception:
            out = {}

    def one(j):
        el, fn, impl, potid = j
        try:
            return (el, fn, props(impl, potid, el))
        except Exception as ex:
            return (el, fn, {"error": f"{type(ex).__name__}: {ex}"})

    with ThreadPoolExecutor(max_workers=6) as ex:
        for el, fn, r in ex.map(one, jobs):
            out[f"{el}|{fn}"] = r
            _s = r.get("surface") or {}
            n = len((_s.get(next(iter(_s), "")) or {}).get("gamma", {}))
            print(f"  {el:3s} {fn[:34]:34s} surface {n:2d} facets"
                  f"{'  stacking' if r.get('stacking') else ''}"
                  f"{'  expansion' if r.get('alpha_V') else ''}"
                  f"{'  ' + r['error'][:30] if 'error' in r else ''}")

    tmp = p + ".tmp"
    json.dump(out, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, p)
    ns = sum(1 for r in out.values() if r.get("surface"))
    nk = sum(1 for r in out.values() if r.get("stacking"))
    na = sum(1 for r in out.values() if r.get("alpha_V"))
    print(f"\n{len(out)} records: surface {ns}, stacking {nk}, expansion {na}"
          f"  -> {p}")


if __name__ == "__main__":
    main()
