#!/usr/bin/env python3
"""
Prove `pair_style ugur/ang` - the angular form - against the angular tree.

Two things have to be right and they fail differently.

The **energy** is checked against `angular/latdyn.py`, which is the
implementation the UG fits were made with.  That catches a wrong Legendre
polynomial, a wrong sign on lam4, a factor in h(cos).

The **force** is checked by conservation.  h(cos theta) is the first thing in
this project that puts a force component across the legs rather than along
them, and an error there is invisible in the energy: the static number stays
right while the dynamics quietly stop conserving.  So each element is also run
in NVE, and the drift has to land where the angle-free form lands - hundredths
of a meV/atom/ps, not tens.

The angular tree defines its own `latdyn` and must never share a sys.path with
the standalone one, so it is called in a subprocess rather than imported.

    python validate_ang.py          # the tap_ug elements with nonzero weights
    python validate_ang.py V W
"""
import io
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import refdata          # noqa: E402
import cellfile         # noqa: E402

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"
SKIN = 2.0

POTFILE = """# {el}, UG form, written by validate_ang.py
# el1 el2 el3  m D alpha r0 gamma C alpha3 rcut2 rcut3 taper lam2 lam4
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

IN_STATIC = """units           metal
boundary        p p p
atom_style      atomic
read_data       {data}
pair_style      ugur/ang
pair_coeff      * * {pot} {el}
neighbor        {skin} bin
neigh_modify    delay 0 every 1 check no
run             0
variable        epa equal pe/atoms
print           "UGUR_E ${{epa}}"
"""

IN_NVE = IN_STATIC.replace("""run             0
variable        epa equal pe/atoms
print           "UGUR_E ${{epa}}"
""", """velocity        all create 600.0 12345 mom yes rot yes
fix             1 all nve
timestep        0.002
thermo          500
thermo_style    custom step temp pe etotal
thermo_modify   norm yes
run             2000
reset_timestep  0
run             20000
""")

#  the angular tree computes the reference; it is run out of process because
#  its latdyn is not the same module as the one imported above
REF = r"""
import json, sys
sys.path.insert(0, r"{ang}")
import latdyn as A, refdata as R
rec = json.loads(sys.argv[1])
e = R.ELEMENTS[sys.argv[2]]
cry = A.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                mass=R.MASSES[sys.argv[2]])
pot = A.Potential(**{{k: rec[k] for k in
                     ("m","D","alpha","r0","gamma","C","alpha3",
                      "rcut2","rcut3","taper","lam2","lam4")}})
print("REF %.12f" % A.energy(cry, pot))
"""


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def reference(el, rec):
    src = REF.format(ang=os.path.join(ROOT, "angular"))
    r = subprocess.run([sys.executable, "-c", src, json.dumps(rec), el],
                       capture_output=True, text=True,
                       cwd=os.path.join(ROOT, "angular"))
    m = re.search(r"REF\s+([-\d.eE+]+)", r.stdout)
    if not m:
        return None, (r.stderr.strip().splitlines() or ["?"])[-1][:90]
    return float(m.group(1)), None


def run(el, rec, cry, tag, template):
    d = os.path.join(HERE, "angruns", f"{el}_{tag}")
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, f"{el}.ugur"), "w").write(POTFILE.format(
        el=el, **{k: rec[k] for k in
                  ("m", "D", "alpha", "r0", "gamma", "C", "alpha3",
                   "rcut2", "rcut3", "taper", "lam2", "lam4")}))
    rcut = max(rec["rcut2"], rec["rcut3"])
    box, _ = cellfile.orthogonal_cell(cry)
    rep = tuple(max(3, int(np.ceil(2.0 * (rcut + SKIN) / b))) for b in box)
    cellfile.write_data(os.path.join(d, f"{el}.data"), cry,
                        refdata.MASSES[el], rep)
    open(os.path.join(d, "in.run"), "w").write(template.format(
        data=f"{el}.data", skin=SKIN, pot=f"{el}.ugur", el=el))
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in in.run 2>&1"],
                   capture_output=True, text=True)
    #  Read the log rather than the captured stdout.  A 30-minute run's stdout
    #  came back truncated through the wsl subprocess - the thermo table was
    #  complete on disk and missing from what Python saw, so the drift read as
    #  "could not parse" on a run that had finished perfectly.  The file is the
    #  authority.
    lg = os.path.join(d, "log.lammps")
    return io.open(lg, errors="ignore").read() if os.path.exists(lg) else ""


ROW = re.compile(r"^\s*(\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*$",
                 re.M)


def drift(out, dt=0.002, nrun=20000):
    if "reset_timestep" not in out:
        return None, None
    rows = [(int(m.group(1)), float(m.group(2)), float(m.group(3)),
             float(m.group(4)))
            for m in ROW.finditer(out.split("reset_timestep")[-1])]
    rows = [x for x in rows if x[0] <= nrun]
    if not rows or rows[-1][0] != nrun or len(rows) < 3:
        return None, None
    ps = np.array([x[0] * dt / 1000.0 for x in rows])
    et = np.array([x[3] * 1000.0 for x in rows])
    return float(np.polyfit(ps, et, 1)[0]), float(np.mean([x[1] for x in rows]))


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    els = sys.argv[1:] or [e for e in sorted(lib)
                           if lib[e].get("tap_ug")
                           and not lib[e]["tap_ug"].get("lam_off")]
    print("pair_style ugur/ang, UG formu: enerji angular/latdyn'e karsi,")
    print("kuvvet NVE korunumuyla (600 K, 2 fs, 40 ps).\n")
    print(f"{'el':4s}{'lam2':>8s}{'lam4':>8s}{'E ref':>12s}{'E lammps':>12s}"
          f"{'dE':>11s}{'drift':>10s}{'T':>6s}   status")
    print("-" * 74)
    bad = 0
    for el in els:
        rec = lib[el].get("tap_ug")
        if not rec:
            print(f"{el:4s}  no tapered UG record")
            continue
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        ref, err = reference(el, rec)
        if ref is None:
            print(f"{el:4s}  referans hatasi: {err}")
            bad += 1
            continue
        out = run(el, rec, cry, "static", IN_STATIC)
        m = re.search(r"UGUR_E\s+([-\d.eE+]+)", out)
        if not m:
            e0 = [l for l in out.splitlines() if "ERROR" in l]
            print(f"{el:4s}  LAMMPS: {(e0 or ['no output'])[0][:60]}")
            bad += 1
            continue
        got = float(m.group(1))
        d, T = drift(run(el, rec, cry, "nve", IN_NVE))
        ok = (abs(got - ref) < 1e-8 * max(abs(ref), 1.0)
              and d is not None and abs(d) < 1.0)
        bad += not ok
        print(f"{el:4s}{rec['lam2']:8.3f}{rec['lam4']:8.3f}{ref:12.6f}"
              f"{got:12.6f}{got-ref:11.2e}"
              f"{(f'{d:10.4f}' if d is not None else f'{chr(63):>10s}')}"
              f"{(f'{T:6.0f}' if T is not None else f'{chr(63):>6s}')}"
              f"   {'ok' if ok else 'FAILED'}")
    print()
    if bad:
        raise SystemExit(f"{bad} element uyusmadi")
    print("the angular form matches angular/latdyn on energy, and the force")
    print("enerjinin turevi - NVE korunumu bunu gosteriyor")


if __name__ == "__main__":
    main()
