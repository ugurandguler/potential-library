#!/usr/bin/env python3
"""
hcp Ru at finite temperature, under the conditions of the group's own paper.

E. Guler, S. Ugur, M. Guler, G. Ugur, "Molecular dynamics exploration of the
temperature-dependent elastic, mechanical and anisotropic properties of hcp
ruthenium", Eur. Phys. J. Plus 139:372 (2024) computed C_ij(T) from 0 to 1200 K
with a Finnis-Sinclair potential, by the Born stress-fluctuation method with
numerical derivatives, in a 108-atom cell.  The question here is what the same
calculation gives when the potential is replaced by ours.

The protocol is not reconstructed by guesswork.  The paper's "150000 steps with
1500 intervals" is exactly nrun = 100*nthermo, nthermo = 1500 of LAMMPS's own
examples/ELASTIC_T/BORN_MATRIX template, so that template is what is used here,
term for term:

    C_ij = <B_ij>  -  (V/kT) cov(sigma_i, sigma_j)  +  (N kT/V) * (Kronecker)

with the Born term from `compute born/matrix numdiff`, which is the only route
available for this pair style - pair_ugur has no analytic born_matrix method,
and neither did the Finnis-Sinclair style the paper used.

Two cells are run, deliberately, because the paper's own numbers disagree with
each other.  Its Table 1 quotes experiment as a = 2.705, c = 4.281 (c/a =
1.582), and its method section states the simulated cell as a = 2.720,
c = 4.441 (c/a = 1.633, the ideal value).  Those are not the same crystal:
the second has a 3.7 % longer c axis, and C33 and C13 in hcp are strongly c/a
dependent.  Reporting one without the other would hide which of the two any
agreement or disagreement belongs to.

Only the tapered sets are run.  Ru's hard-truncated MAU and UG both fail the
MD screen outright - Ru has the deepest step at the cutoff in the library,
-0.346 eV - and a hard-truncated potential has no business in a 150 ps NVT run.

    python ru_born.py --check          # T=0 Born term vs our own static C_ij
    python ru_born.py --delta-scan     # numdiff sensitivity, as the docs require
    python ru_born.py                  # the sweep, 0-1200 K in 200 K steps
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

HOME = subprocess.run(["wsl", "-e", "bash", "-lc", "echo $HOME"],
                      capture_output=True, text=True).stdout.strip()
LMP = f"{HOME}/lammps/src/lmp_serial"
RUNS = os.path.join(HERE, "ruruns")

MASS = refdata.MASSES["Ru"]
TEMPS = [0, 200, 400, 600, 800, 1000, 1200]

#  the two cells, kept apart on purpose - see the docstring
CELLS = {
    "paper":  dict(a=2.720, coa=4.441 / 2.720,
                   note="method section of the paper, c/a = 1.633 (ideal)"),
    "expt":   dict(a=2.71, coa=1.5793,
                   note="experiment as our refdata holds it, c/a = 1.579"),
}

#  (label, library key, pair style, file extension)
SETS = {
    "MAU_tap": ("tap", "ugur", "ugur"),
    "UG_tap": ("tap_ug", "ugur/ang", "ugur.ang"),
}

POTFILE = """# Ru, written by ru_born.py
Ru Ru Ru {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

#  ---------------------------------------------------------------- static ---
#  At T = 0 there is no fluctuation term and no kinetic term, so C = B alone.
#  That makes this the check on the whole chain: if the Born term at the
#  minimised structure does not reproduce the C_ij our own static code gives
#  for the same cell, something in the units, the ordering or the pair style's
#  virial is wrong, and no finite-temperature number afterwards means anything.
STATIC = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       Ru.data
pair_style      {style}
pair_coeff      * * Ru.{ext} Ru
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
min_style       cg
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
variable        pea equal pe/atoms
print           "BORN_PE ${{pea}}"
"""

#  ------------------------------------------------------------------- MD ---
#  This is examples/ELASTIC_T/BORN_MATRIX rewritten for a read_data cell and
#  our pair style.  The three terms, the running averages, the sampling
#  intervals and the unit constants are the template's, unchanged.
MD = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       Ru.data

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

pair_style      {style}
pair_coeff      * * Ru.{ext} Ru
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

velocity        all create ${{temp}} 87287 mom yes rot yes
timestep        0.001
fix             4 all nve
fix             5 all langevin ${{temp}} ${{temp}} 0.1 123457

thermo          ${{nthermo}}
thermo_style    custom step temp pe press
run             ${{nequil}}

reset_timestep  0

#  --- stress fluctuation term F ---
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
variable        mypress equal press
variable        mype equal pe/atoms
fix             avt all ave/time ${{nevery}} ${{nrepeat}} ${{nfreq}} v_mytemp ave running
fix             avp all ave/time ${{nevery}} ${{nrepeat}} ${{nfreq}} v_mypress ave running
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

#  --- Born term ---
compute         virial all pressure NULL virial
compute         born all born/matrix numdiff ${{delta}} virial
fix             avborn all ave/time ${{neveryborn}} ${{nrepeatborn}} ${{nfreq}} c_born[*] ave running
variable        bfac equal 1.0e-4*${{nktv2p}}/vol
variable        B vector f_avborn*${{bfac}}

#  --- kinetic term ---
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


def wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:]


def setup(d, setname, cell):
    """write the data file and the potential file into directory d"""
    os.makedirs(d, exist_ok=True)
    key, style, ext = SETS[setname]
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    rec = dict(lib["Ru"][key])
    rec.setdefault("lam2", 0.0)
    rec.setdefault("lam4", 0.0)
    rec["taper"] = rec.get("taper") or -1.0
    open(os.path.join(d, f"Ru.{ext}"), "w").write(POTFILE.format(
        **{k: rec[k] for k in ("m", "D", "alpha", "r0", "gamma", "C",
                               "alpha3", "rcut2", "rcut3", "taper",
                               "lam2", "lam4")}))
    cry = L.Crystal("hcp", cell["a"], cell["coa"], mass=MASS)
    #  3x3x3 of the orthogonal 4-atom cell = 108 atoms, the paper's size
    cellfile.write_data(os.path.join(d, "Ru.data"), cry, MASS, (3, 3, 3))
    return style, ext, cry


def run(d, text, name="in.ru"):
    open(os.path.join(d, name), "w").write(text)
    subprocess.run(["wsl", "-e", "bash", "-lc",
                    f"cd {wsl(d)} && {LMP} -in {name} > out.txt 2>&1"],
                   capture_output=True, text=True)
    p = os.path.join(d, "log.lammps")
    return io.open(p, errors="ignore").read() if os.path.exists(p) else ""


#  ------------------------------------------------------------------------
def born_static(d, style, ext, delta):
    """C_ij from the Born term alone at the minimised structure"""
    lg = run(d, STATIC.format(style=style, ext=ext, delta=delta), "in.static")
    vals = {}
    for m in re.finditer(r"BORN\s+(\d+)\s+([-\d.eE+]+)", lg):
        vals[int(m.group(1))] = float(m.group(2))
    if len(vals) < 21:
        return None
    return vals


def check(delta=1.0e-6):
    """T=0 Born term against our own static elastic constants"""
    print("T = 0 check: the Born term alone should give C_ij.")
    print("Compared against latdyn's static C_ij in the same cell.\n")
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    for cname, cell in CELLS.items():
        for sname in SETS:
            key, style, ext = SETS[sname]
            d = os.path.join(RUNS, f"chk_{cname}_{sname}")
            style, ext, cry = setup(d, sname, cell)
            v = born_static(d, style, ext, delta)
            if v is None:
                print(f"{cname:6s} {sname:8s}  born/matrix cikmadi")
                continue
            pot = L.Potential.from_record(lib["Ru"][key])
            #  the comparison must be against the FROZEN-ION curvature: that is
            #  what compute born/matrix measures.  The relaxed value is printed
            #  beside it because the gap between them is the non-affine term,
            #  which at T > 0 is what the stress-fluctuation term supplies -
            #  so its size says how much of the answer the T = 0 row is missing.
            Crel, Cb = L.elastic(cry, pot)
            print(f"--- {cname} (a={cell['a']:.3f}, c/a={cell['coa']:.4f})"
                  f"  {sname} ---")
            print(f"{'':6s}{'LAMMPS Born':>13s}{'latdyn Born':>13s}"
                  f"{'diff %':>9s}{'latdyn relaxed':>15s}")
            #  ortho cell: x and y are both basal, z = c.  Voigt order is the
            #  same on both sides (4=yz, 5=xz, 6=xy), checked in the sources.
            pairs = [("C11", v[1], Cb[0, 0], Crel[0, 0]),
                     ("C22", v[2], Cb[1, 1], Crel[1, 1]),
                     ("C33", v[3], Cb[2, 2], Crel[2, 2]),
                     ("C44", v[4], Cb[3, 3], Crel[3, 3]),
                     ("C55", v[5], Cb[4, 4], Crel[4, 4]),
                     ("C66", v[6], Cb[5, 5], Crel[5, 5]),
                     ("C12", v[7], Cb[0, 1], Crel[0, 1]),
                     ("C13", v[8], Cb[0, 2], Crel[0, 2]),
                     ("C23", v[12], Cb[1, 2], Crel[1, 2])]
            for nm, got, want, rel in pairs:
                e = 100 * (got - want) / abs(want) if abs(want) > 1e-9 else 0.0
                print(f"{nm:6s}{got:13.1f}{want:13.1f}{e:9.2f}{rel:15.1f}")
            print(f"{'':6s}C66 vs (C11-C12)/2: {v[6]:.1f} vs "
                  f"{0.5*(v[1]-v[7]):.1f}   (altigen ozdesligi)")
            print()


def delta_scan():
    """the sensitivity study the docs say is the only way to pick delta"""
    print("numdiff delta taramasi (T=0, paper hucresi, MAU_tap).")
    print("Dokuman: cok kucukse gurultu, cok buyukse yuksek mertebe.\n")
    print(f"{'delta':>10s}{'C11':>10s}{'C33':>10s}{'C44':>10s}{'C12':>10s}")
    print("-" * 50)
    for dl in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
        d = os.path.join(RUNS, f"scan_{dl:.0e}")
        style, ext, cry = setup(d, "MAU_tap", CELLS["paper"])
        v = born_static(d, style, ext, dl)
        if v is None:
            print(f"{dl:10.0e}   cikmadi"); continue
        print(f"{dl:10.0e}{v[1]:10.1f}{v[3]:10.1f}{v[4]:10.1f}{v[7]:10.1f}")


def one_md(job):
    cname, sname, T, delta = job
    d = os.path.join(RUNS, f"{cname}_{sname}_{T}K")
    style, ext, cry = setup(d, sname, CELLS[cname])
    if T == 0:
        v = born_static(d, style, ext, delta)
        if v is None:
            return None
        return {"C11": 0.5 * (v[1] + v[2]), "C33": v[3],
                "C44": 0.5 * (v[4] + v[5]), "C66": v[6],
                "C12": v[7], "C13": 0.5 * (v[8] + v[12]),
                "Tavg": 0.0, "static": True}
    lg = run(d, MD.format(T=float(T), style=style, ext=ext, delta=delta))
    if "RES C11" not in lg:
        return None
    r = {m.group(1): float(m.group(2))
         for m in re.finditer(r"RES\s+(\w+)\s+([-\d.eE+]+)", lg)}
    if len(r) < 11:
        return None
    return {"C11": 0.5 * (r["C11"] + r["C22"]), "C33": r["C33"],
            "C44": 0.5 * (r["C44"] + r["C55"]), "C66": r["C66"],
            "C12": r["C12"], "C13": 0.5 * (r["C13"] + r["C23"]),
            "Tavg": r["Tavg"], "PEavg": r["PEavg"], "static": False}


def sweep(delta=1.0e-6):
    jobs = [(c, s, T, delta) for c in CELLS for s in SETS for T in TEMPS]
    nw = max(1, (os.cpu_count() or 4) - 2)
    print(f"hcp Ru, 108 atom, NVT, Born dalgalanma yontemi "
          f"({len(jobs)} runs, {nw} in parallel)\n")
    with ThreadPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(one_md, jobs))
    out = {}
    for (c, s, T, _), r in zip(jobs, res):
        out.setdefault(f"{c}|{s}", {})[str(T)] = r
    json.dump(out, open(os.path.join(HERE, "ru_born.json"), "w"),
              indent=1, sort_keys=True)
    for k in sorted(out):
        c, s = k.split("|")
        print(f"=== {c} hucresi, {s} ===  ({CELLS[c]['note']})")
        print(f"{'T':>6s}{'C11':>9s}{'C33':>9s}{'C44':>9s}{'C12':>9s}"
              f"{'C13':>9s}{'C66':>9s}{'(C11-C12)/2':>13s}{'Tolc':>8s}")
        for T in TEMPS:
            r = out[k][str(T)]
            if r is None:
                print(f"{T:6d}   hesaplanamadi"); continue
            print(f"{T:6d}{r['C11']:9.1f}{r['C33']:9.1f}{r['C44']:9.1f}"
                  f"{r['C12']:9.1f}{r['C13']:9.1f}{r['C66']:9.1f}"
                  f"{0.5*(r['C11']-r['C12']):13.1f}{r['Tavg']:8.0f}")
        print()
    print("-> ru_born.json")


if __name__ == "__main__":
    os.makedirs(RUNS, exist_ok=True)
    if "--check" in sys.argv:
        check()
    elif "--delta-scan" in sys.argv:
        delta_scan()
    else:
        sweep()


#  --------------------------------------------------------------- report ---
#  Table 1 of the paper, "This work" column: hcp Ru at T = 0 K, P = 0 GPa,
#  Finnis-Sinclair, in the a = 2.720 / c = 4.441 cell.  Experiment is the
#  column the paper itself quotes, ref. [6] there.
PAPER_T0 = {"C11": 597.1, "C12": 162.2, "C13": 174.0, "C33": 720.2,
            "C44": 239.7, "C66": 217.4}
EXPT = {"C11": 576.0, "C12": 187.0, "C13": 167.0, "C33": 641.0,
        "C44": 189.0, "C66": 194.5}
#  and the cells the three of us settle on
CELL_REF = {"FS (paper)": (2.720, 4.441), "experiment [24]": (2.705, 4.281)}


def report():
    """our numbers beside the paper's and experiment's"""
    out = json.load(open(os.path.join(HERE, "ru_born.json")))
    ks = ("C11", "C12", "C13", "C33", "C44", "C66")
    print("hcp Ru, T = 0 K.  Makalenin Tablo 1'i ile yan yana.\n")
    hdr = f"{'':6s}{'expt':>8s}{'FS(pap.)':>10s}"
    for cn in CELLS:
        for sn in SETS:
            hdr += f"{sn + '/' + cn:>14s}"
    print(hdr)
    print("-" * len(hdr))
    for k in ks:
        row = f"{k:6s}{EXPT[k]:8.1f}{PAPER_T0[k]:10.1f}"
        for cn in CELLS:
            for sn in SETS:
                r = out.get(f"{cn}|{sn}", {}).get("0")
                row += f"{r[k]:14.1f}" if r else f"{'-':>14s}"
        print(row)
    #  Table 1 does not stop at the stiffness constants, so neither does this
    print()
    mk = ("B", "G", "E", "B/G", "nu", "AU")
    hdr = f"{'':6s}{'expt':>8s}{'FS(pap.)':>10s}"
    for cn in CELLS:
        for sn in SETS:
            hdr += f"{sn + '/' + cn:>14s}"
    print(hdr)
    print("-" * len(hdr))
    he, hp = hill(EXPT), hill(PAPER_T0)
    for k in mk:
        row = f"{k:6s}{he[k]:8.2f}{hp[k]:10.2f}"
        for cn in CELLS:
            for sn in SETS:
                r = out.get(f"{cn}|{sn}", {}).get("0")
                row += f"{hill(r)[k]:14.2f}" if r else f"{'-':>14s}"
        print(row)

    #  a single number for "how close", against both references
    print()
    for cn in CELLS:
        for sn in SETS:
            r = out.get(f"{cn}|{sn}", {}).get("0")
            if not r:
                continue
            for nm, ref in (("experiment", EXPT), ("paper", PAPER_T0)):
                rms = np.sqrt(np.mean([((r[k] - ref[k]) / ref[k]) ** 2
                                       for k in ks])) * 100
                print(f"{sn}/{cn:6s} vs {nm:7s} rms {rms:6.2f} %")

    #  and the temperature dependence, which is what the paper is actually
    #  about - its Table 1 is one row, its figures are the result
    print()
    for k in sorted(out):
        cn, sn = k.split("|")
        print(f"=== {sn}, {cn} hucresi ===")
        print(f"{'T':>6s}{'C11':>9s}{'C12':>9s}{'C13':>9s}{'C33':>9s}"
              f"{'C44':>9s}{'C66':>9s}{'B':>9s}{'G':>9s}{'Tolc':>7s}")
        for T in TEMPS:
            r = out[k].get(str(T))
            if r is None:
                print(f"{T:6d}   hesaplanamadi")
                continue
            h = hill(r)
            print(f"{T:6d}{r['C11']:9.1f}{r['C12']:9.1f}{r['C13']:9.1f}"
                  f"{r['C33']:9.1f}{r['C44']:9.1f}{r['C66']:9.1f}"
                  f"{h['B']:9.1f}{h['G']:9.1f}{r['Tavg']:7.0f}")
        print()


#  ------------------------------------------------- mechanical averages ---
#  The paper reports Hill averages, not only C_ij, so a comparison that stops
#  at the stiffness constants leaves half its Table 1 unanswered.  These are
#  the standard hexagonal Voigt-Reuss-Hill expressions; they are checked
#  against the paper's own row in `hill_selftest`, which is the only way to be
#  sure the convention matches the one it used.
def hill(C):
    c11, c12, c13, c33, c44 = (C["C11"], C["C12"], C["C13"],
                               C["C33"], C["C44"])
    c66 = 0.5 * (c11 - c12)
    Bv = (2 * (c11 + c12) + 4 * c13 + c33) / 9.0
    Gv = (c11 + c12 + 2 * c33 - 4 * c13 + 12 * c44 + 12 * c66) / 30.0
    C2 = (c11 + c12) * c33 - 2 * c13 * c13
    Br = C2 / (c11 + c12 + 2 * c33 - 4 * c13)
    Gr = 2.5 * (C2 * c44 * c66) / (3 * Bv * c44 * c66 + C2 * (c44 + c66))
    B, G = 0.5 * (Bv + Br), 0.5 * (Gv + Gr)
    return {"B": B, "G": G, "E": 9 * B * G / (3 * B + G),
            "B/G": B / G, "nu": (3 * B - 2 * G) / (2 * (3 * B + G)),
            "AU": 5 * Gv / Gr + Bv / Br - 6.0, "C66": c66}


def hill_selftest():
    """do these formulas reproduce the paper's own averages from its own C_ij?"""
    got = hill(PAPER_T0)
    want = {"B": 324.7, "G": 232.4, "E": 563.0, "B/G": 1.4,
            "nu": 0.21, "AU": 0.03}
    print("Hill formullerinin makalenin kendi satirini yeniden vermesi:")
    print(f"{'':6s}{'ours':>9s}{'paper':>9s}")
    for k in ("B", "G", "E", "B/G", "nu", "AU"):
        print(f"{k:6s}{got[k]:9.2f}{want[k]:9.2f}")
