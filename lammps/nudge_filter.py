#!/usr/bin/env python3
"""
Pick, from each element's pool of fitted solutions, the best one that holds its
lattice.

The nudge test needs LAMMPS, and LAMMPS cannot go inside the search: Nelder-
Mead evaluates thousands of candidates per restart and a minimisation each time
would cost more than the whole fit.  So the search runs unchanged and keeps
every distinct solution it finds (dense_fit.py, N_POOL), and the constraint is
applied here, afterwards, to a few dozen candidates instead of thousands.

That is only worth doing because the pool is close together.  Chromium's best
eight solutions span 63.1 to 66.6 per cent, so if the winner fails and the
fourth passes, the constraint costs three points of RMS rather than the
element.  If none of them passes, that is the answer - the same answer lithium
gave when 900 restarts found no solution satisfying the compression
constraint, and it is a statement about the functional form rather than about
the search.

What "holds its lattice" means is measured, not judged: displace every atom by
1e-5 A, minimise with the cell fixed, and ask how much of the displacement
survives.  Sound records keep 0.01-0.3 of it; unsound ones keep 8000-16000 and
settle 0.2-1.3 meV/atom lower.  Four orders of magnitude separate the two.

It is measured along FIVE directions, not one, and a record counts as sound
only if it returns along all of them.  One direction is one sample: iron's best
candidate comes back under seed 87287 and keeps eight to ten thousand times the
displacement under the other four, so the single-direction test had been
certifying a record that holds its lattice along a line.  As a detector of
broken records one direction is plenty - they fail every direction by four
orders of magnitude - but certifying is a different job.

    python nudge_filter.py <dir with dense_*.json>
    python nudge_filter.py <dir> --set tap_ug
"""
import io
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import cellfile         # noqa: E402

LOCAL = not os.path.exists("/arf")      # windows+wsl here, plain shell on TRUBA
#  Ask WSL where home is only where WSL exists.  This ran unconditionally at
#  import, so on the cluster the module died on its first line - after the fit
#  had already spent thirty-seven minutes.  The launcher resubmits whenever
#  nudge_picked.json is missing, so it repeated that six times before anyone
#  read the traceback.  The pools survive in dense_*.json, which is the only
#  reason it cost time rather than work.
HOME = (subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                       capture_output=True, text=True).stdout.strip()
        if LOCAL else os.path.expanduser("~"))
LMP = os.environ.get("LMP", f"{HOME}/lammps/src/lmp_serial")
SKIN = 2.0
AMP = 1.0e-5
KEEP_MAX = 10.0
DE_MAX = -5.0e-5

POTFILE = """# {el}, written by nudge_filter.py
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

IN = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       cell.data
pair_style      {style}
pair_coeff      * * pot.ugur {el}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
min_modify      dmax 1.0e-3 line quadratic
compute         disp all displace/atom
displace_atoms  all random {amp} {amp} {amp} {seed} units box
variable        e1 equal pe/atoms
run             0
print           "E_NUDGED ${{e1}}"
minimize        0 1e-12 100000 1000000
variable        e2 equal pe/atoms
print           "E_RELAXED ${{e2}}"
variable        d2 atom c_disp[4]*c_disp[4]
compute         sd all reduce sum v_d2
variable        rms equal sqrt(c_sd/atoms)
run             0
print           "RMS ${{rms}}"
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def shell(d, cmd):
    if LOCAL:
        return subprocess.run(["wsl", "-e", "bash", "-lc",
                               f"cd {wsl(d)} && {cmd}"],
                              capture_output=True, text=True)
    return subprocess.run(f"cd {d} && {cmd}", shell=True,
                          capture_output=True, text=True)


#  One displacement direction is one sample, and for iron that turned out to be
#  the whole answer: its best candidate returns to its lattice under seed 87287
#  and keeps eight to ten thousand times the displacement under four others.
#  The single-direction test is a good detector - broken records fail every
#  direction by four orders of magnitude - but it cannot certify a record as
#  sound, which is exactly what it was being used for here.
SEEDS = (87287, 11311, 40213, 90007, 55501)


def test_many(el, rec, style, tag, seeds=SEEDS):
    """the nudge test along several directions; sound means all of them"""
    keeps, des = [], []
    for s in seeds:
        v = test(el, rec, style, f"{tag}_s{s}", seed=s)
        keeps.append(v["keep"] if v else None)
        des.append(v["dE"] if v else None)
    good = [k for k in keeps if k is not None]
    npass = sum(1 for i, k in enumerate(keeps)
                if k is not None and (k < KEEP_MAX or des[i] > DE_MAX))
    return {"keep": (max(good) if good else None),
            "keep_all": keeps, "dE_all": des, "seeds": list(seeds),
            "n_pass": npass, "n": len(seeds),
            "ok": bool(good and len(good) == len(seeds)
                       and npass == len(seeds))}


def test(el, rec, style, tag, seed=87287):
    d = os.path.join(HERE, "nudgefilter", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    p = dict(rec)
    p.setdefault("lam2", 0.0)
    p.setdefault("lam4", 0.0)
    p["taper"] = rec.get("taper") or -1.0
    p.setdefault("alpha3", rec.get("s3", 1.0) * rec["alpha"])
    open(os.path.join(d, "pot.ugur"), "w").write(POTFILE.format(
        el=el, **{k: p[k] for k in ("m", "D", "alpha", "r0", "gamma", "C",
                                    "alpha3", "rcut2", "rcut3", "taper",
                                    "lam2", "lam4")}))
    e = refdata.ELEMENTS[el]
    cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                    mass=refdata.MASSES[el])
    box, _ = cellfile.orthogonal_cell(cry)
    rc = max(rec["rcut2"], rec["rcut3"])
    rep = tuple(max(3, int(np.ceil(1.8 * (rc + SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, "cell.data"), cry,
                        refdata.MASSES[el], rep)
    open(os.path.join(d, "in.jig"), "w").write(
        IN.format(style=style, el=el, skin=SKIN, amp=AMP, seed=seed))
    lg = os.path.join(d, "log.lammps")
    #  Remove the previous log BEFORE running.  Without this a run that fails
    #  to start leaves the previous invocation's log in place and it is read as
    #  a fresh result, and because the candidate pools are regenerated on every
    #  submission, candidate 008 of one run is a different parameter set from
    #  candidate 008 of the next.  This was found while chasing iron's bad
    #  verdict and is NOT what caused it - fixing it left the verdict
    #  unchanged, which is how the real cause (one displacement direction is
    #  one sample; see SEEDS) came to light.  Fixed anyway.
    if os.path.exists(lg):
        os.remove(lg)
    shell(d, f"{LMP} -in in.jig > o.txt 2>&1")
    if not os.path.exists(lg):
        return None
    t = io.open(lg, errors="ignore").read()
    g = {}
    for k in ("E_NUDGED", "E_RELAXED", "RMS"):
        m = re.search(rf"{k}\s+([-\d.eE+]+)", t)
        if m:
            g[k] = float(m.group(1))
    if "RMS" not in g:
        return None
    keep = g["RMS"] / AMP
    dE = g["E_RELAXED"] - g["E_NUDGED"]
    return {"keep": keep, "dE": dE,
            "ok": bool(keep < KEEP_MAX or dE > DE_MAX)}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    style = "ugur"
    if "--set" in sys.argv:
        if sys.argv[sys.argv.index("--set") + 1] == "tap_ug":
            style = "ugur/ang"
    src = args[0] if args else "."
    out = {}
    for f in sorted(os.listdir(src)):
        if not (f.startswith("dense_") and f.endswith(".json")):
            continue
        try:
            raw = json.load(open(os.path.join(src, f)))
        except Exception:
            continue
        for el, rec in raw.items():
            if not isinstance(rec, dict) or "pool" not in rec:
                continue
            out.setdefault(el, []).extend(rec["pool"])
    if not out:
        print("no dense_*.json carrying a pool was found")
        return
    print(f"{'el':4s}{'cand':>6s}{'best rms':>12s}{'first pass':>11s}"
          f"{'bedeli':>9s}   secilen")
    print("-" * 62)
    picked = {}
    for el in sorted(out):
        #  distinct again: the same solution can appear in several seeds
        seen, cand = set(), []
        for r in sorted(out[el], key=lambda x: x["score"]):
            sig = tuple(round(float(r[k]), 6)
                        for k in ("m", "gamma", "D", "alpha", "r0", "C", "s3"))
            if sig in seen:
                continue
            seen.add(sig)
            cand.append(r)
        nw = max(1, (os.cpu_count() or 4) - 2)
        with ThreadPoolExecutor(max_workers=nw) as ex:
            res = list(ex.map(
                lambda ir: test_many(el, ir[1], style, f"{ir[0]:03d}"),
                list(enumerate(cand))))
        best0 = cand[0]["score"] * 100
        hit = next(((r, v) for r, v in zip(cand, res) if v and v["ok"]), None)
        if hit is None:
            print(f"{el:4s}{len(cand):6d}{best0:12.2f}{'-':>11s}{'-':>9s}"
                  f"   GECEN YOK")
            continue
        r, v = hit
        rank = cand.index(r) + 1
        print(f"{el:4s}{len(cand):6d}{best0:12.2f}{r['score']*100:11.2f}"
              f"{r['score']*100 - best0:+9.2f}   {rank}. sirada, "
              f"kalan/sarsma {v['keep']:.2f} ({v['n_pass']}/{v['n']} yon)")
        picked[el] = dict(r, nudge=v, rank=rank, best_rms=best0)
    json.dump(picked, open(os.path.join(src, "nudge_picked.json"), "w"),
              indent=1, sort_keys=True)
    print(f"\na passing solution was found for {len(picked)}/{len(out)} elements"
          f"  -> {os.path.join(src, 'nudge_picked.json')}")


if __name__ == "__main__":
    main()
