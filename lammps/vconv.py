#!/usr/bin/env python3
"""
Does vanadium's finite-temperature collapse survive longer sampling?

The Born stress-fluctuation method divides the stress covariance by kT, so the
statistical error in that term grows as the temperature falls.  Vanadium's Born
violations sit between 0.05 and 0.60 T_melt and vanish above, which is the
shape a sampling artefact would have; potassium and rubidium, whose first
finite-T point lands at 17 and 16 K, return C11 of -880 and -1432, which is
what that artefact looks like when it is unmistakable.  Tantalum, same
structure, same reduced temperature, is clean - so it is not universal.

The only way to separate the two is to sample longer and watch whether the
answer moves.  A real elastic constant converges; a covariance that has not
been sampled wanders.  Each element is run at the standard length and at three
and eight times it, with everything else - cell, seed, thermostat, delta -
identical to the sweep.

Vanadium's two published baselines are included for the same reason: they went
through this pipeline unchanged, so if they hold their lattice where ours
collapses, the instability belongs to the potential and not to the method.

    python vconv.py
    python vconv.py --els V,Ta --mults 1,3
"""
import io
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
#  on the cluster the two trees are copied side by side and this is a no-op;
#  here latdyn and cellfile still live in standalone/
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
os.environ.setdefault("LMP", f"{HOME}/lammps/src/lmp_serial")
#  baselines.json lives in standalone/ here; os.path.join keeps an absolute
#  second argument, which is what elastic_T does with BASEFILE
os.environ.setdefault("BASEFILE",
                      os.path.join(ROOT, "standalone", "baselines.json"))

import elastic_T as E      # noqa: E402


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


#  elastic_T.py is written for the cluster, where the shell is the shell.  Here
#  LAMMPS lives inside WSL, so the one call that reaches it is redirected and
#  nothing else about the module changes.
def run_wsl(d, text, name):
    open(os.path.join(d, name), "w").write(text)
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {E.LMP} -in {name} > out.txt 2>&1"],
                   capture_output=True, text=True)
    p = os.path.join(d, "log.lammps")
    return io.open(p, errors="ignore").read() if os.path.exists(p) else ""


#  C_ij is a sum of three terms and only one of them can be noisy.  The Born
#  term is a numerical second derivative of the potential, the kinetic term is
#  an exact function of N, T and V, and the fluctuation term is the sampled
#  stress covariance divided by kT.  Printing them apart says immediately
#  whether a low C11 is the potential going soft or the covariance not having
#  converged.
DECOMP = """
variable        B11s equal v_B[1]
variable        B12s equal v_B[7]
variable        B44s equal v_B[4]
print           "DEC B11 ${B11s}"
print           "DEC B12 ${B12s}"
print           "DEC B44 ${B44s}"
print           "DEC F11 ${F11}"
print           "DEC F12 ${F12}"
print           "DEC F44 ${F44}"
print           "DEC K   ${kfac}"
"""


def one(el, tag, T, mult):
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "vconv", f"{el}_{safe}_{T}K_x{mult}")
    try:
        style, pot = E.setup(d, el, tag)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    text = E.MD.format(T=float(T), style=style, el=el, delta=E.DELTA,
                       coeff=E.coeff_line(style, pot, el)) + DECOMP
    #  the single change that matters: production length.  Equilibration is
    #  left alone, so the starting state is the sweep's starting state.
    old = "variable        nrun      equal 100*${nthermo}"
    assert old in text
    text = text.replace(old, f"variable        nrun      equal "
                             f"{100 * mult}*${{nthermo}}")
    lg = E.cached(d) or run_wsl(d, text, "in.md")
    r = {m.group(1): float(m.group(2))
         for m in re.finditer(r"RES\s+(\w+)\s+([-\d.eE+]+)", lg)}
    if len(r) < 11:
        return {"error": "the run did not finish"}
    dec = {m.group(1): float(m.group(2))
           for m in re.finditer(r"DEC\s+(\w+)\s+([-\d.eE+]+)", lg)}
    return {"C11": (r["C11"] + r["C22"]) / 2, "C12": r["C12"],
            "C44": (r["C44"] + r["C55"]) / 2, "Tavg": r["Tavg"], **dec}


def main():
    els = ["V", "Ta"]
    mults = [1, 3, 8]
    T = 300
    for flag, conv in (("--els", None), ("--mults", int), ("--T", int)):
        if flag in sys.argv:
            v = sys.argv[sys.argv.index(flag) + 1]
            if flag == "--els":
                els = v.split(",")
            elif flag == "--mults":
                mults = [int(x) for x in v.split(",")]
            else:
                T = int(v)

    jobs = []
    for el in els:
        jobs += [(el, "tap", T, m) for m in mults]
        #  every shipped potential for the same element, at the longest length
        #  only: the question they answer is whether the pipeline is sound, not
        #  whether they converge.  Skipped when the file is not staged locally,
        #  since the main sweep already ran them through this pipeline.
        for fn, _ in E.BASE.get(el, []):
            if all(os.path.exists(os.path.join(E.POTDIR, f))
                   for f in fn.split("+")):
                jobs.append((el, "base|" + fn, T, max(mults)))

    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"{T} K, {len(jobs)} runs, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(lambda j: one(*j), jobs))

    print()
    print(f"{'el':4s}{'source':>26s}{'sample':>7s}{'C11':>9s}{'C12':>9s}"
          f"{'C44':>8s}{'C11-C12':>10s}{'Tolc':>7s}  Born")
    print("-" * 82)
    last = None
    for (el, tag, T_, m), r in zip(jobs, res):
        if el != last:
            print("-" * 82) if last else None
            last = el
        lab = {"tap": "bizimki (MAU)"}.get(tag, tag.replace("base|", ""))[:26]
        if "error" in r:
            print(f"{el:4s}{lab:>26s}{m:6d}x   {r['error']}")
            continue
        d1 = r["C11"] - r["C12"]
        ok = "OK" if (d1 > 0 and r["C44"] > 0
                      and r["C11"] + 2 * r["C12"] > 0) else "IHLAL"
        print(f"{el:4s}{lab:>26s}{m:6d}x{r['C11']:9.1f}{r['C12']:9.1f}"
              f"{r['C44']:8.1f}{d1:10.1f}{r['Tavg']:7.0f}  {ok}")

    print("\n--- terimlere ayrilmis (GPa) ---")
    print(f"{'el':4s}{'source':>26s}{'sample':>7s}"
          f"{'B11':>9s}{'F11':>9s}{'B12':>9s}{'F12':>9s}{'kin':>8s}")
    print("-" * 72)
    for (el, tag, T_, m), r in zip(jobs, res):
        if "error" in r or "B11" not in r:
            continue
        lab = {"tap": "bizimki (MAU)"}.get(tag, tag.replace("base|", ""))[:26]
        print(f"{el:4s}{lab:>26s}{m:6d}x{r['B11']:9.1f}{r['F11']:9.1f}"
              f"{r['B12']:9.1f}{r['F12']:9.1f}{r['K']:8.2f}")
    print("\nsample: multiple of the baseline run (150 000 steps).")
    print("B = Born terimi (potansiyelin ikinci turevi), F = dalgalanma "
          "terimi (gerilme kovaryansi / kT), kin = kinetik terim.")


if __name__ == "__main__":
    main()
