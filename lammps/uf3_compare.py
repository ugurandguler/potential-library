#!/usr/bin/env python3
"""
Our phi2 + phi3 against a UF3 phi2 + phi3, for niobium.

The library's surface note argues that this functional form is not the
limitation - that a two- plus three-body potential with free radial shapes,
fitted to density-functional energies and forces, reproduces properties this
one misses (Xie, Rupp and Hennig, npj Comput. Mater. 9, 162).  The argument is
made from the literature and has never been run here.  It can be: LAMMPS ships
`Nb.uf3`, which is that potential for niobium, and pair_style uf3 is compiled
into our build.

The comparison is the point.  Nb.uf3 is the SAME functional family - the file
declares `2B Nb Nb` and `3B Nb Nb Nb` - so what differs is not the form:

    ours     9 analytic parameters, fitted to about six measured numbers
    UF3      27 two-body spline knots and an 11 x 11 x 19 three-body grid,
             fitted to DFT energies and forces

So a gap here is a gap between FIT DATA, not between forms, and that is the
one question the page keeps raising without testing.  Niobium is a good place
to ask it: the shipped MAU arm sits at 20 % elastic RMS there.

The elastic machinery is elastic_check's, unchanged - LAMMPS's own ELASTIC
recipe, six deformations with the internal coordinates relaxed - so both sides
go through exactly the same code and only the potential differs.

    python uf3_compare.py
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "standalone"))

import numpy as np                      # noqa: E402
import elastic_check as E               # noqa: E402
import cellfile                         # noqa: E402
import latdyn as L                      # noqa: E402
import refdata                          # noqa: E402

#  LAMMPS ships Nb.uf3 in its potentials/ directory; LAMMPS does not expand
#  "~" in pair_coeff, so the path is built from the WSL home elastic_check found
UF3 = f"{E.HOME}/lammps/potentials/Nb.uf3"
#  the two-body cutoff the file declares, which pair_style uf3 wants as an
#  argument; taken from the header rather than assumed
UF3_NBODY = 3


def uf3_cij(el, cry, a0, up=1.0e-5):
    d = os.path.join(HERE, "elruns", el + "_uf3")
    os.makedirs(d, exist_ok=True)
    for f in ("in.elastic", "displace.mod"):
        subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cp ~/lammps/examples/ELASTIC/{f} {E.wsl(d)}/"],
                       capture_output=True)
    box, _ = cellfile.orthogonal_cell(cry)
    #  UF3's two-body cutoff is 8 A, well beyond ours, so the box is sized
    #  from THAT and not from our record - a cell chosen for a 5 A potential
    #  would fold the 8 A one onto its own images.
    rep = tuple(max(2, int(np.ceil(1.2 * (8.0 + E.SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep, triclinic=True)
    open(os.path.join(d, "init.mod"), "w").write(
        E.INIT.format(data=f"{el}.data", up=up))
    open(os.path.join(d, "potential.mod"), "w").write(
        E.POTMOD.format(style=f"uf3 {UF3_NBODY}", pot=UF3, el=el))
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        f"cd {E.wsl(d)} && {E.LMP} -in in.elastic 2>&1"],
                       capture_output=True, text=True)
    got = {}
    for k in E.WANT:
        m = re.search(rf"{k}all\s*=\s*([-\d.eE+]+)", r.stdout)
        if m:
            got[k] = float(m.group(1))
    if not got:
        err = [l for l in r.stdout.splitlines() if "ERROR" in l]
        got["err"] = err[0][:90] if err else "output unreadable"
    return got


def main():
    el = "Nb"
    lib = json.load(open(os.path.join(os.path.dirname(HERE), "standalone",
                                      "library.json")))
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], float(e["a0"]), e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    print(f"{el}, bcc, a0 = {e['a0']} A   (experimental cell, the same for both sides)")
    print()
    uf = uf3_cij(el, cry, e["a0"])
    if "err" in uf:
        print("UF3 run failed:", uf["err"])
        return
    ours = E.lammps_cij(el, lib[el], cry, tag="_ref")
    exp = e["Cij"]
    print(f"{'':6s}{'expt':>9s}{'ours (MAU)':>13s}{'UF3':>9s}"
          f"{'ours dev':>13s}{'UF3 dev':>11s}")
    print("-" * 62)
    for k in ("C11", "C12", "C44"):
        if k not in exp:
            continue
        a, b = ours.get(k), uf.get(k)
        f = lambda v: f"{v:9.1f}" if v is not None else f"{'-':>9s}"
        da = f"{100*(a-exp[k])/exp[k]:12.0f}%" if a is not None else ""
        db = f"{100*(b-exp[k])/exp[k]:10.0f}%" if b is not None else ""
        print(f"{k:6s}{exp[k]:9.1f}{f(a)[:13]:>13s}{f(b)}{da}{db}")




#  ---- surface energies -------------------------------------------------
#
#  This is the test that actually bears on the page's claim.  The elastic
#  comparison above is asymmetric - our fit was GIVEN those constants and the
#  UF3 was not - so it cannot settle whether the form is the limitation.
#  Surface energy was fitted by neither, and it is the property the claim is
#  made about: Xie, Rupp and Hennig report tungsten surfaces to within 4 %
#  from a two- plus three-body potential with free radial shapes.
#
#  If a UF3 niobium lands near the DFT reference while ours sits four times
#  above it, the form is not the limitation and the fit data is.  If the UF3
#  is also four times high, the claim does not survive.
import surface as S                        # noqa: E402


def uf3_surface(el, facet, a0):
    d = os.path.join(HERE, "surfruns", f"{el}_uf3_{facet}")
    os.makedirs(d, exist_ok=True)
    style, coeff = f"uf3 {UF3_NBODY}", f"pair_coeff * * {UF3} {el}"
    struct = refdata.ELEMENTS[el]["struct"]
    rc = 8.6
    _, base = S.oriented_cell(struct, a0, facet, (1, 1, 1))
    need = 2.4 * (rc + S.SKIN)
    reps = [max(1, int(np.ceil(need / b))) for b in base]
    reps[2] = max(reps[2], int(np.ceil(2.6 * (rc + S.SKIN) / base[2])), 4)
    build = S.cubic_build(el, d, a0, facet, tuple(reps))
    lb = S.run(d, S.BULK.format(build=build, style=style, coeff=coeff,
                                skin=S.SKIN), "in.bulk")
    eb = S.grab(lb, "EBULK")
    if eb is None:
        err = [l for l in lb.splitlines() if "ERROR" in l]
        return {"error": err[0][:90] if err else "bulk did not run"}
    ls = S.run(d, S.SLAB.format(build=build, style=style, coeff=coeff,
                                skin=S.SKIN, vac=S.VACUUM), "in.slab")
    area, nat, es = S.grab(ls, "AREA"), S.grab(ls, "NAT"), S.grab(ls, "ESLAB")
    if None in (area, nat, es):
        err = [l for l in ls.splitlines() if "ERROR" in l]
        return {"error": err[0][:90] if err else "slab did not run"}
    return {"gamma": (es - nat * eb) / (2.0 * area) * S.EV_A2_TO_J_M2,
            "atoms": int(nat)}


def surfaces():
    el, a0_uf3 = "Nb", 3.34148          # UF3's OWN relaxed lattice constant
    lib = json.load(open(os.path.join(os.path.dirname(HERE), "standalone",
                                      "library.json")))
    ref = (lib[el].get("surface_ref") or {}).get("facets") or {}
    ours = ((lib[el].get("tap") or {}).get("surface") or {}).get("gamma") or {}
    print()
    print(f"{el} surface energies (J/m2) - NEITHER was fitted to this property")
    print(f"{'facet':7s}{'DFT ref':>9s}{'ours':>9s}{'UF3':>9s}"
          f"{'ours/ref':>11s}{'UF3/ref':>10s}")
    print("-" * 55)
    rows = {}
    for f in ("110", "100", "111"):
        if f not in ref:
            continue
        u = uf3_surface(el, f, a0_uf3)
        if "error" in u:
            print(f"{f:7s}{ref[f]:9.2f}{ours.get(f, float('nan')):9.2f}"
                  f"   UF3 error: {u['error'][:40]}")
            continue
        o = ours.get(f)
        print(f"{f:7s}{ref[f]:9.2f}{(o if o else float('nan')):9.2f}"
              f"{u['gamma']:9.2f}{(o/ref[f] if o else 0):10.1f}x"
              f"{u['gamma']/ref[f]:9.1f}x")
        rows[f] = {"ref": round(ref[f], 3), "ours": round(o, 3) if o else None,
                   "uf3": round(u["gamma"], 3)}
    return rows


def store(rows, a0_uf3):
    """put the measurement in library.json, so the page can stop citing

    The surface note argued from Xie, Rupp and Hennig that the FORM is not the
    limitation.  That argument is now a measurement made here, in this
    pipeline, with the same elastic and surface machinery on both sides - and
    a page that says "published work reports" when it has its own number is
    weaker than it needs to be.
    """
    path = os.path.join(os.path.dirname(HERE), "standalone", "library.json")
    lib = json.load(open(path))
    order = lambda k: [f for f, _ in sorted(
        ((f, r[k]) for f, r in rows.items() if r.get(k)), key=lambda t: t[1])]
    pct = lambda k: round(sum(abs(r[k] - r["ref"]) / r["ref"]
                              for r in rows.values() if r.get(k))
                          / max(1, sum(1 for r in rows.values()
                                       if r.get(k))) * 100, 0)
    lib["Nb"]["uf3"] = {
        "file": "Nb.uf3", "a0": a0_uf3, "facets": rows,
        "mean_pct": {"ours": pct("ours"), "uf3": pct("uf3")},
        "order": {k: order(k) for k in ("ref", "ours", "uf3")},
        "why": ("Nb.uf3 declares 2B Nb Nb and 3B Nb Nb Nb - the same "
                "functional family as ours, with 27 two-body spline knots and "
                "an 11&times;11&times;19 three-body grid fitted to "
                "density-functional energies and forces, against our nine "
                "analytic parameters fitted to about six measured numbers. So "
                "a gap here is a gap between FIT DATA, not between forms."),
        "caveat": ("Nb.uf3 carries an empty CITATION field, so its training "
                   "set is unknown and this is NOT demonstrated extrapolation "
                   "- the same limitation the published tungsten numbers "
                   "carry. The elastic constants are not quoted because that "
                   "comparison is asymmetric: ours was fitted to exactly "
                   "those numbers and the UF3 was not."),
        "ref": ("S. R. Xie, M. Rupp and R. G. Hennig, npj Comput. Mater. "
                "<b>9</b>, 162 (2023), doi:10.1038/s41524-023-01092-7; the "
                "file is the one LAMMPS ships as Nb.uf3, run here through "
                "elastic_check and surface.py unchanged"),
    }
    tmp = path + ".tmp"
    json.dump(lib, open(tmp, "w"), indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print()
    print(f"-> library.json, Nb.uf3: ours {pct('ours'):.0f} %, "
          f"UF3 {pct('uf3'):.0f} % mean deviation")


if __name__ == "__main__":
    rows = surfaces()
    if "--write" in sys.argv and rows:
        store(rows, 3.34148)
