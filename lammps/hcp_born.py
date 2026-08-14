#!/usr/bin/env python3
"""
Finite-temperature elastic constants of every hexagonal element in the library.

This is the ruthenium calculation from the group's own paper (Guler, Ugur,
Guler, Ugur, Eur. Phys. J. Plus 139:372, 2024) applied to all twelve, using
LAMMPS's own examples/ELASTIC_T/BORN_MATRIX recipe unchanged - the Born matrix
by numerical differentiation, plus the stress-fluctuation and kinetic terms:

    C_ij = <B_ij> - (V/kT) cov(s_i, s_j) + (N kT/V) * Kronecker

108 atoms (3x3x3 of the orthogonal four-atom hcp cell), 0 to 1200 K in 200 K
steps, nthermo = 1500, nrun = 100*nthermo, which is what the paper's "150000
steps with 1500 intervals" means.

Three things are deliberate.

**Only the tapered sets.**  A hard-truncated potential does not conserve energy
- the step at the cutoff pumps it in - and has no business in a 150 ps run.

**Published potentials go through the identical pipeline.**  LAMMPS ships EAM/FS
for magnesium and zirconium and a spline MEAM for titanium.  Running those here,
with the same cell, the same recipe and the same post-processing, is what turns
a number into a comparison; without it the result is uncalibrated, which is the
mistake the nickel-aluminium work already made once.  Three of twelve is not
coverage, but it is a scale.

**The experimental cell, not a relaxed one.**  Every set was fitted at the
experimental lattice constant, so that is where its elastic constants mean what
they were fitted to mean.  Relaxing first would compare each potential at its
own volume and confound the comparison with the equation of state.

Temperatures above an element's melting point are run and reported, but flagged:
a perfect 108-atom crystal has nowhere to nucleate from and will happily
superheat, so those rows describe a metastable solid, not the material.

    python3 hcp_born.py            # everything
    python3 hcp_born.py Mg Zr
"""
import io
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import cellfile         # noqa: E402

LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/src/lmp_serial"))
POTDIR = os.path.expanduser("~/lammps/potentials")
PACK = json.load(open(os.path.join(HERE, "hcp_pack.json")))
TEMPS = [0, 200, 400, 600, 800, 1000, 1200]
DELTA = 1.0e-6
REP = (3, 3, 3)

#  our two switched sets
OURS = {"tap": ("ugur", "ugur"), "tap_ug": ("ugur/ang", "ugur.ang")}
#  and the published potentials LAMMPS ships for three of these elements
BASE = {
    "Mg": ("eam/fs", "Mg_mm.eam.fs", "EAM/FS, Sun et al."),
    "Zr": ("eam/fs", "Zr_mm.eam.fs", "EAM/FS, Mendelev-Ackland"),
    "Ti": ("meam/spline", "Ti.meam.spline", "spline MEAM, Hennig et al."),
}

POTFILE = """# {el}, written by hcp_born.py
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

HEAD = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       cell.data
pair_style      {style}
pair_coeff      * * {pot} {el}
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
"""

STATIC = HEAD + """min_style       cg
minimize        0 1e-12 10000 100000
compute         virial all pressure NULL virial
compute         born all born/matrix numdiff {delta} virial
variable        bfac equal 1.0e-4*1.6021765e6/vol
run             0
variable        i loop 21
label           lp
variable        v equal c_born[${{i}}]*${{bfac}}
print           "BORN ${{i}} ${{v}}"
variable        v delete
next            i
jump            SELF lp
"""

MD = HEAD + """
variable        temp   equal {T}
variable        nthermo   equal 1500
variable        nevery    equal 10
variable        neveryborn equal 100
variable        delta     equal {delta}
variable        nfreq     equal ${{nthermo}}
variable        nrepeat   equal floor(${{nfreq}}/${{nevery}})
variable        nrepeatborn equal floor(${{nfreq}}/${{neveryborn}})
variable        nequil    equal 10*${{nthermo}}
variable        nrun      equal 100*${{nthermo}}

run             0
variable        elat equal pe/atoms
print           "LATTICE ${{elat}}"

velocity        all create ${{temp}} 87287 mom yes rot yes
timestep        0.001
fix             4 all nve
fix             5 all langevin ${{temp}} ${{temp}} 0.1 123457
thermo          ${{nthermo}}
thermo_style    custom step temp pe press
run             ${{nequil}}
reset_timestep  0

compute         stress all pressure thermo_temp
variable        s1 equal c_stress[1]
variable        s2 equal c_stress[2]
variable        s3 equal c_stress[3]
variable        s4 equal c_stress[6]
variable        s5 equal c_stress[5]
variable        s6 equal c_stress[4]
variable        s11 equal v_s1*v_s1
variable        s22 equal v_s2*v_s2
variable        s33 equal v_s3*v_s3
variable        s44 equal v_s4*v_s4
variable        s55 equal v_s5*v_s5
variable        s66 equal v_s6*v_s6
variable        s12 equal v_s1*v_s2
variable        s13 equal v_s1*v_s3
variable        s14 equal v_s1*v_s4
variable        s15 equal v_s1*v_s5
variable        s16 equal v_s1*v_s6
variable        s23 equal v_s2*v_s3
variable        s24 equal v_s2*v_s4
variable        s25 equal v_s2*v_s5
variable        s26 equal v_s2*v_s6
variable        s34 equal v_s3*v_s4
variable        s35 equal v_s3*v_s5
variable        s36 equal v_s3*v_s6
variable        s45 equal v_s4*v_s5
variable        s46 equal v_s4*v_s6
variable        s56 equal v_s5*v_s6

variable        mytemp equal temp
variable        mype equal pe/atoms
fix             avt all ave/time ${{nevery}} ${{nrepeat}} ${{nfreq}} v_mytemp ave running
fix             avpe all ave/time ${{nevery}} ${{nrepeat}} ${{nfreq}} v_mype ave running
fix             avs all ave/time ${{nevery}} ${{nrepeat}} ${{nfreq}} v_s1 v_s2 v_s3 v_s4 v_s5 v_s6 ave running
fix             avssq all ave/time ${{nevery}} ${{nrepeat}} ${{nfreq}} &
                v_s11 v_s22 v_s33 v_s44 v_s55 v_s66 &
                v_s12 v_s13 v_s14 v_s15 v_s16 &
                v_s23 v_s24 v_s25 v_s26 &
                v_s34 v_s35 v_s36 &
                v_s45 v_s46 &
                v_s56 ave running

variable        boltz equal 8.617343e-5
variable        nktv2p equal 1.6021765e6
variable        vkt equal vol/(${{boltz}}*${{temp}})/${{nktv2p}}
variable        ffac equal 1.0e-4*${{vkt}}
variable        F11 equal -(f_avssq[1]-f_avs[1]*f_avs[1])*${{ffac}}
variable        F22 equal -(f_avssq[2]-f_avs[2]*f_avs[2])*${{ffac}}
variable        F33 equal -(f_avssq[3]-f_avs[3]*f_avs[3])*${{ffac}}
variable        F44 equal -(f_avssq[4]-f_avs[4]*f_avs[4])*${{ffac}}
variable        F55 equal -(f_avssq[5]-f_avs[5]*f_avs[5])*${{ffac}}
variable        F66 equal -(f_avssq[6]-f_avs[6]*f_avs[6])*${{ffac}}
variable        F12 equal -(f_avssq[7]-f_avs[1]*f_avs[2])*${{ffac}}
variable        F13 equal -(f_avssq[8]-f_avs[1]*f_avs[3])*${{ffac}}
variable        F23 equal -(f_avssq[12]-f_avs[2]*f_avs[3])*${{ffac}}

compute         virial all pressure NULL virial
compute         born all born/matrix numdiff ${{delta}} virial
fix             avborn all ave/time ${{neveryborn}} ${{nrepeatborn}} ${{nfreq}} c_born[*] ave running
variable        bfac equal 1.0e-4*${{nktv2p}}/vol
variable        B vector f_avborn*${{bfac}}
variable        kfac equal 1.0e-4*${{nktv2p}}*atoms*${{boltz}}*${{temp}}/vol

thermo_style    custom step temp pe press f_avt f_avpe
thermo_modify   norm no
run             ${{nrun}}

variable        C11 equal v_F11+v_B[1]+4.0*${{kfac}}
variable        C22 equal v_F22+v_B[2]+4.0*${{kfac}}
variable        C33 equal v_F33+v_B[3]+4.0*${{kfac}}
variable        C44 equal v_F44+v_B[4]+2.0*${{kfac}}
variable        C55 equal v_F55+v_B[5]+2.0*${{kfac}}
variable        C66 equal v_F66+v_B[6]+2.0*${{kfac}}
variable        C12 equal v_F12+v_B[7]
variable        C13 equal v_F13+v_B[8]
variable        C23 equal v_F23+v_B[12]
print           "RES C11 ${{C11}}"
print           "RES C22 ${{C22}}"
print           "RES C33 ${{C33}}"
print           "RES C44 ${{C44}}"
print           "RES C55 ${{C55}}"
print           "RES C66 ${{C66}}"
print           "RES C12 ${{C12}}"
print           "RES C13 ${{C13}}"
print           "RES C23 ${{C23}}"
variable        tavg equal f_avt
variable        peavg equal f_avpe
print           "RES Tavg ${{tavg}}"
print           "RES PEavg ${{peavg}}"
"""


def setup(d, el, tag):
    """write cell.data and the potential file; return (pair_style, pot name)"""
    os.makedirs(d, exist_ok=True)
    p = PACK[el]
    cry = L.Crystal("hcp", p["a0"], p["c_over_a"], mass=p["mass"])
    cellfile.write_data(os.path.join(d, "cell.data"), cry, p["mass"], REP)
    if tag in OURS:
        style, ext = OURS[tag]
        name = f"{el}.{ext}"
        open(os.path.join(d, name), "w").write(
            POTFILE.format(el=el, **{k: p[tag][k] for k in
                                     ("m", "D", "alpha", "r0", "gamma", "C",
                                      "alpha3", "rcut2", "rcut3", "taper",
                                      "lam2", "lam4")}))
        return style, name
    style, fn, _ = BASE[el]
    #  copied in rather than referenced by path, so the run directory is a
    #  complete record of what was run
    open(os.path.join(d, fn), "wb").write(
        open(os.path.join(POTDIR, fn), "rb").read())
    return style, fn


def run(d, text, name):
    open(os.path.join(d, name), "w").write(text)
    subprocess.run(f"cd {d} && {LMP} -in {name} > out.txt 2>&1",
                   shell=True, capture_output=True, text=True)
    p = os.path.join(d, "log.lammps")
    return io.open(p, errors="ignore").read() if os.path.exists(p) else ""


def one(job):
    el, tag, T = job
    d = os.path.join(HERE, "runs", f"{el}_{tag}_{T}K")
    try:
        style, pot = setup(d, el, tag)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    if T == 0:
        lg = run(d, STATIC.format(style=style, pot=pot, el=el, delta=DELTA),
                 "in.static")
        v = {int(m.group(1)): float(m.group(2))
             for m in re.finditer(r"BORN\s+(\d+)\s+([-\d.eE+]+)", lg)}
        if len(v) < 21:
            return {"error": "born/matrix cikmadi"}
        return {"C11": 0.5 * (v[1] + v[2]), "C33": v[3],
                "C44": 0.5 * (v[4] + v[5]), "C66": v[6],
                "C12": v[7], "C13": 0.5 * (v[8] + v[12]),
                "Tavg": 0.0, "static": True}
    lg = run(d, MD.format(T=float(T), style=style, pot=pot, el=el,
                          delta=DELTA), "in.md")
    r = {m.group(1): float(m.group(2))
         for m in re.finditer(r"RES\s+(\w+)\s+([-\d.eE+]+)", lg)}
    if len(r) < 11:
        return {"error": "the run did not finish"}
    ml = re.search(r"LATTICE\s+([-\d.eE+]+)", lg)
    lat = float(ml.group(1)) if ml else None
    out = {"C11": 0.5 * (r["C11"] + r["C22"]), "C33": r["C33"],
           "C44": 0.5 * (r["C44"] + r["C55"]), "C66": r["C66"],
           "C12": r["C12"], "C13": 0.5 * (r["C13"] + r["C23"]),
           "Tavg": r["Tavg"], "PEavg": r["PEavg"], "lattice": lat,
           "static": False}
    #  the same two guards the MD screen uses: a run that lost its structure
    #  still prints numbers, and they are numbers about a different solid
    out["hot"] = bool(r["Tavg"] > 1.5 * T)
    out["collapsed"] = bool(lat is not None and r["PEavg"] < lat - 0.05)
    out["above_melt"] = bool(T > PACK[el]["Tmelt"])
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    els = args or sorted(PACK)
    jobs = [(el, tag, T) for el in els for tag in OURS for T in TEMPS]
    jobs += [(el, "base", T) for el in els if el in BASE for T in TEMPS]
    nw = max(1, int(os.environ.get("SLURM_CPUS_ON_NODE",
                                   os.cpu_count() or 4)) - 2)
    print(f"hcp finite-T elasticity: {len(jobs)} runs, {nw} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(one, jobs))
    out = {}
    for (el, tag, T), r in zip(jobs, res):
        out.setdefault(f"{el}|{tag}", {})[str(T)] = r
    json.dump(out, open(os.path.join(HERE, "hcp_born.json"), "w"),
              indent=1, sort_keys=True)
    bad = sum(1 for r in res if r is None or "error" in r)
    print(f"-> hcp_born.json  ({len(res) - bad} ok, {bad} failed)",
          flush=True)


if __name__ == "__main__":
    main()
