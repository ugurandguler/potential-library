#!/usr/bin/env python3
"""
Elastic constants against temperature, for the whole library and for every
published potential LAMMPS ships alongside it.

This is the ruthenium paper's calculation (Guler, Ugur, Guler, Ugur, Eur. Phys.
J. Plus 139:372, 2024) generalised: LAMMPS's own examples/ELASTIC_T/BORN_MATRIX
recipe unchanged, the Born matrix by numerical differentiation plus the
stress-fluctuation and kinetic terms,

    C_ij = <B_ij> - (V/kT) cov(s_i, s_j) + (N kT/V) * Kronecker

with the derived quantities the paper reports - Voigt-Reuss-Hill B, G, E, B/G,
Poisson's ratio, the universal anisotropy A_U and Chen's hardness - computed
afterwards from the tensor.

Three choices that are not the paper's, each because copying it would have
produced something meaningless here.

**The temperature grid follows the melting point, not a fixed 0-1200 K.**
Melting points in this library run from caesium at 302 K to tungsten at 3695.
A fixed grid puts every alkali point above melting and never lifts tungsten
past a third of its range, so the curves would not be comparable between
elements and half of them would describe a superheated solid.  Sampling at
fractions of T_melt makes every element's curve span the same physical range.
300 K is added to all of them as a common reference, and it is skipped where
it exceeds 0.9 T_melt.

**Every available baseline is run, not one per element.**  LAMMPS ships five
EAM sets for copper and three for aluminium, and they do not agree: on the
static elastic constants aluminium spans 5.6 to 26.5 per cent depending which
is chosen.  Picking one is picking an answer.

**Records that fail the nudge test are run and reported, but marked.**  Five
bcc tapered sets do not hold their lattice against a 1e-5 A displacement, and
thermal motion at any temperature here is four orders of magnitude larger, so
whatever is measured is a property of some other structure.  Hiding them would
be worse - the curves are informative about what the potential actually does -
but presenting them beside the sound ones without a mark would be dishonest.

    python3 elastic_T.py                # everything
    python3 elastic_T.py Cu Al --sets tap
"""
import io
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np      # noqa: E402
import latdyn as L      # noqa: E402
import cellfile         # noqa: E402

LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/src/lmp_serial"))
POTDIR = os.path.expanduser("~/lammps/potentials")
PACK = json.load(open(os.path.join(HERE, "elastic_T_pack.json")))
#  BASEFILE lets a follow-up job run only the potentials added since the
#  main sweep started, instead of repeating four hundred runs to reach the
#  twenty new ones
BASE = json.load(open(os.path.join(HERE,
                     os.environ.get("BASEFILE", "baselines.json"))))
DELTA = 1.0e-6
#  fractions of the melting point; 0.05 stands in for "cold" without asking a
#  thermostat to hold a temperature the Langevin damping cannot resolve
FRACS = (0.05, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90)
OURS = {"tap": ("ugur", "ugur"), "tap_ug": ("ugur/ang", "ugur.ang")}

#  A candidate parameter set can be tested exactly like a shipped one by
#  pointing PICKS at a nudge_picked.json.  It joins PACK under the name "pick",
#  so bain.py, struct_rank.py and this sweep all reach it through the same
#  path as tap and tap_ug - which is the point.  A candidate judged by a
#  private code path is a candidate judged by a different measurement.
_picks = os.environ.get("PICKS")
if _picks and os.path.exists(_picks):
    OURS["pick"] = ("ugur", "ugur")
    for _el, _r in json.load(open(_picks)).items():
        if _el not in PACK:
            continue
        _q = dict(_r)
        _q.setdefault("alpha3", _r.get("s3", 1.0) * _r["alpha"])
        _q.setdefault("lam2", 0.0)
        _q.setdefault("lam4", 0.0)
        if not _q.get("taper"):
            _q["taper"] = -1.0
        PACK[_el]["pick"] = _q

POTFILE = """# {el}, written by elastic_T.py
{el} {el} {el} {m:.17g} {D:.17g} {alpha:.17g} {r0:.17g} {gamma:.17g} {C:.17g} {alpha3:.17g} {rcut2:.17g} {rcut3:.17g} {taper:.17g} {lam2:.17g} {lam4:.17g}
"""

HEAD = """units           metal
boundary        p p p
atom_style      atomic
atom_modify     map array
read_data       cell.data
pair_style      {style}
{coeff}
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


def temps(el):
    tm = PACK[el]["Tmelt"]
    ts = sorted({0} | {round(f * tm) for f in FRACS}
                | ({300} if 300 < 0.9 * tm else set()))
    return ts


def setup(d, el, tag):
    os.makedirs(d, exist_ok=True)
    p = PACK[el]
    cry = L.Crystal(p["struct"], p["a0"], p.get("c_over_a"), mass=p["mass"])
    #  the same cell for every potential of a given element, so a difference
    #  between curves is the potential and never the geometry
    cellfile.write_data(os.path.join(d, "cell.data"), cry, p["mass"],
                        tuple(p["rep"]))
    if tag in OURS:
        style, ext = OURS[tag]
        name = f"{el}.{ext}"
        open(os.path.join(d, name), "w").write(
            POTFILE.format(el=el, **{k: p[tag][k] for k in
                                     ("m", "D", "alpha", "r0", "gamma", "C",
                                      "alpha3", "rcut2", "rcut3", "taper",
                                      "lam2", "lam4")}))
        return style, name
    fn, style = tag.split("|", 1)[1], None
    for f, s in BASE[el]:
        if f == fn:
            style = s
    #  MEAM needs two files and a four-argument pair_coeff, and NIST ships
    #  every entry under the same two names, so they are stored with an entry
    #  prefix and carried here as "library+parameter" in one field.
    for one_f in fn.split("+"):
        open(os.path.join(d, one_f), "wb").write(
            open(os.path.join(POTDIR, one_f), "rb").read())
    return style, fn


def coeff_line(style, pot, el):
    """pair_coeff, which is NOT the same for every EAM flavour.

    A plain .eam file is funcfl and takes `pair_coeff I J file` with no element
    name; eam/alloy and eam/fs take setfl and want `pair_coeff * * file El`.
    Writing the setfl form for a funcfl file does not warn - LAMMPS errors out
    and the run simply produces nothing, which reads downstream as "the
    potential failed" rather than "the input was wrong".  Ten of the nineteen
    baselines are funcfl.
    """
    if style == "eam":
        return f"pair_coeff      1 1 {pot}"
    if style == "meam":
        lib, par = pot.split("+")
        return f"pair_coeff      * * {lib} {el} {par} {el}"
    return f"pair_coeff      * * {pot} {el}"


def run(d, text, name):
    open(os.path.join(d, name), "w").write(text)
    #  DEVNULL rather than capture_output: a captured pipe's write end is
    #  inherited by children forked from other threads while it is open, so a
    #  thread can wait forever for an end-of-file a sibling is holding shut.
    #  It is not hypothetical - it stopped the stacking-fault job dead after
    #  sixty-two runs.  Nothing here reads that pipe anyway; LAMMPS writes to
    #  out.txt and the answers come back from log.lammps.
    subprocess.run(f"cd {d} && {LMP} -in {name} > out.txt 2>&1", shell=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p = os.path.join(d, "log.lammps")
    return io.open(p, errors="ignore").read() if os.path.exists(p) else ""


def cached(d):
    """a finished run's log, if this directory already holds one

    Sixteen hours was not enough for this sweep and the first attempt would
    have been killed with 141 runs complete and nothing saved, because results
    were only written at the end.  Every one of those runs left a complete log
    behind; reading it back costs milliseconds and turns a resubmission into a
    resume rather than a restart.
    """
    p = os.path.join(d, "log.lammps")
    if not os.path.exists(p):
        return None
    t = io.open(p, errors="ignore").read()
    if "RES PEavg" in t or t.count("BORN ") >= 21:
        return t
    return None


def one(job):
    el, tag, T = job
    safe = tag.replace("|", "_").replace("/", "-").replace(".", "")
    d = os.path.join(HERE, "eTruns", f"{el}_{safe}_{T}K")
    try:
        style, pot = setup(d, el, tag)
    except Exception as ex:
        return {"error": f"setup: {ex}"}
    if T == 0:
        lg = cached(d) or run(d, STATIC.format(
            style=style, el=el, delta=DELTA,
            coeff=coeff_line(style, pot, el)), "in.static")
        v = {int(m.group(1)): float(m.group(2))
             for m in re.finditer(r"BORN\s+(\d+)\s+([-\d.eE+]+)", lg)}
        if len(v) < 21:
            return {"error": "born/matrix cikmadi"}
        return {"C11": (v[1] + v[2]) / 2, "C33": v[3],
                "C44": (v[4] + v[5]) / 2, "C66": v[6],
                "C12": v[7], "C13": (v[8] + v[12]) / 2,
                "Tavg": 0.0, "static": True}
    lg = cached(d) or run(d, MD.format(
        T=float(T), style=style, el=el, delta=DELTA,
        coeff=coeff_line(style, pot, el)), "in.md")
    r = {m.group(1): float(m.group(2))
         for m in re.finditer(r"RES\s+(\w+)\s+([-\d.eE+]+)", lg)}
    if len(r) < 11:
        return {"error": "the run did not finish"}
    ml = re.search(r"LATTICE\s+([-\d.eE+]+)", lg)
    lat = float(ml.group(1)) if ml else None
    return {"C11": (r["C11"] + r["C22"]) / 2, "C33": r["C33"],
            "C44": (r["C44"] + r["C55"]) / 2, "C66": r["C66"],
            "C12": r["C12"], "C13": (r["C13"] + r["C23"]) / 2,
            "Tavg": r["Tavg"], "PEavg": r["PEavg"], "lattice": lat,
            "static": False,
            "hot": bool(r["Tavg"] > 1.3 * T),
            "collapsed": bool(lat is not None and r["PEavg"] < lat - 0.05)}


def main():
    #  Parse the option and its VALUE together.  Filtering on "--" alone left
    #  the value of --sets in the element list, so `elastic_T.py --sets base`
    #  asked for an element called "base", found none, planned zero runs and
    #  reported success.  The two jobs that happened to also name elements were
    #  unaffected, which is why it took a job that did nothing to notice.
    argv = sys.argv[1:]
    want = None
    if "--sets" in argv:
        i = argv.index("--sets")
        want = argv[i + 1].split(",")
        del argv[i:i + 2]
    args = [a for a in argv if not a.startswith("--")]
    els = args or sorted(PACK)
    if any(a not in PACK for a in args):
        bad = [a for a in args if a not in PACK]
        print(f"bilinmeyen element: {' '.join(bad)}", flush=True)
        return
    jobs = []
    for el in els:
        for tag in OURS:
            if want and tag not in want:
                continue
            if tag in PACK[el]:
                jobs += [(el, tag, T) for T in temps(el)]
        for f, _s in BASE.get(el, []):
            if want and "base" not in want:
                continue
            jobs += [(el, f"base|{f}", T) for T in temps(el)]
    nw = max(1, int(os.environ.get("SLURM_CPUS_ON_NODE",
                                   os.cpu_count() or 4)) - 2)
    print(f"Finite-T elasticity: {len(jobs)} runs, {nw} in parallel", flush=True)
    #  Each job writes its own file.  Three of these run at once and a shared
    #  read-then-write merge loses whichever finishes second.
    p = os.path.join(HERE, os.environ.get("OUTFILE", "elastic_T.json"))
    out = {}
    if os.path.exists(p):
        try:
            out = json.load(open(p))
        except Exception:
            out = {}

    def save():
        #  through a temporary file and renamed, so a job killed mid-write
        #  leaves the previous complete file rather than a truncated one
        tmp = p + ".tmp"
        json.dump(out, open(tmp, "w"), indent=1, sort_keys=True)
        os.replace(tmp, p)

    res = []
    with ThreadPoolExecutor(max_workers=nw) as ex:
        futs = {ex.submit(one, j): j for j in jobs}
        for k, fut in enumerate(as_completed(futs), 1):
            el, tag, T = futs[fut]
            try:
                r = fut.result()
            except Exception as ex_:
                r = {"error": str(ex_)[:80]}
            res.append(r)
            out.setdefault(f"{el}|{tag}", {})[str(T)] = r
            if k % 20 == 0:
                save()
                print(f"  {k}/{len(jobs)} saved", flush=True)
    save()
    bad = sum(1 for r in res if r is None or "error" in r)
    print(f"-> {os.path.basename(p)}  ({len(res) - bad} ok, {bad} failed, "
          f"{len(out)} records in total)", flush=True)


if __name__ == "__main__":
    main()
