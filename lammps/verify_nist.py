#!/usr/bin/env python3
"""
Check every downloaded NIST potential before it is allowed to be a baseline.

A file that does not load, or loads and returns nonsense, is worse than no
baseline at all: it becomes a grey curve on a plot that a reader will take as
the published state of the art.  So each one is asked for the static elastic
constants at the experimental lattice constant and compared with experiment.

The comparison is deliberately loose.  These potentials were fitted at their
own equilibrium volumes, not at the experimental one, so a real disagreement of
tens of per cent is expected and fine - the plot measures every potential in the
same cell on purpose, and that is a stated choice, not an error.  What this
rejects is the other thing: a file that fails to parse, a pair_style mismatch,
or a tensor that violates the Born criteria and therefore is not describing a
solid at all.
"""
import io, json, os, re, subprocess, sys
sys.path.insert(0, os.path.abspath("../standalone")); sys.path.insert(0, ".")
import numpy as np, latdyn as L, refdata, cellfile
from concurrent.futures import ThreadPoolExecutor

HOME=subprocess.run(["wsl","-e","bash","-lc","echo $HOME"],capture_output=True,text=True).stdout.strip()
LMP=f"{HOME}/lammps/src/lmp_serial"
def wsl(p):
    p=os.path.abspath(p).replace("\\","/"); return "/mnt/"+p[0].lower()+p[2:]

IN="""units metal
boundary p p p
atom_style atomic
atom_modify map array
read_data cell.data
pair_style {ps}
{coeff}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
min_style cg
minimize 0 1e-12 10000 100000
compute virial all pressure NULL virial
compute born all born/matrix numdiff 1.0e-6 virial
variable bfac equal 1.0e-4*1.6021765e6/vol
run 0
variable i loop 21
label lp
variable v equal c_born[${{i}}]*${{bfac}}
print "BORN ${{i}} ${{v}}"
variable v delete
next i
jump SELF lp
"""

def one(job):
    el, pid, ps, fn = job
    d=os.path.join("nistchk", f"{el}_{pid[:22]}"); os.makedirs(d, exist_ok=True)
    e=refdata.ELEMENTS[el]
    cry=L.Crystal(e["struct"], e["a0"], e.get("c_over_a"), mass=refdata.MASSES[el])
    box,_=cellfile.orthogonal_cell(cry)
    rep=tuple(max(3,int(np.ceil(1.6*8.0/b))) for b in box)
    cellfile.write_data(os.path.join(d,"cell.data"), cry, refdata.MASSES[el], rep)
    coeff = f"pair_coeff 1 1 {fn}" if ps=="eam" else f"pair_coeff * * {fn} {el}"
    open(os.path.join(d,"in.x"),"w").write(IN.format(ps=ps, coeff=coeff))
    subprocess.run(["wsl","-e","bash","-lc",
        f"cp ~/nistpot/{fn} {wsl(d)}/ 2>/dev/null; cd {wsl(d)} && {LMP} -in in.x > o.txt 2>&1"],
        capture_output=True)
    lg=os.path.join(d,"log.lammps")
    if not os.path.exists(lg): return (el,pid,ps,fn,None,"no log")
    t=io.open(lg,errors="ignore").read()
    v={int(m.group(1)):float(m.group(2)) for m in re.finditer(r"BORN\s+(\d+)\s+([-\d.eE+]+)",t)}
    if len(v)<21:
        err=[l for l in t.splitlines() if "ERROR" in l]
        return (el,pid,ps,fn,None, err[0][:60] if err else "born cikmadi")
    c={"C11":0.5*(v[1]+v[2]),"C12":v[7],"C13":0.5*(v[8]+v[12]),
       "C33":v[3],"C44":0.5*(v[4]+v[5])}
    return (el,pid,ps,fn,c,None)

cand=json.load(open("../standalone/nist_candidates.json"))
have=json.load(open("../standalone/baselines.json"))
jobs=[]
for el in sorted(cand):
    if el in have: continue
    for x in cand[el]:
        for f in x["files"]:
            jobs.append((el, x["id"], x["pair_style"], f["name"]))
nw=max(1,(os.cpu_count() or 4)-2)
print(f"{len(jobs)} aday sinaniyor\n")
with ThreadPoolExecutor(max_workers=nw) as ex: res=list(ex.map(one, jobs))

good={}
print(f"{'el':4s}{'pair_style':12s}{'file':26s}{'C11':>8s}{'C12':>8s}{'C44':>8s}   verdict")
print("-"*84)
for el,pid,ps,fn,c,err in res:
    ex_=refdata.ELEMENTS[el]["Cij"]
    if c is None:
        print(f"{el:4s}{ps:12s}{fn[:25]:26s}{'':>24s}   DUSTU: {err}"); continue
    st=refdata.ELEMENTS[el]["struct"]
    born = ((c["C11"]>abs(c["C12"]) and (c["C11"]+c["C12"])*c["C33"]>2*c["C13"]**2
             and c["C44"]>0) if st=="hcp" else
            (c["C11"]-c["C12"]>0 and c["C44"]>0 and c["C11"]+2*c["C12"]>0))
    dev=100*abs(c["C11"]-ex_["C11"])/ex_["C11"]
    ok = born and dev<200
    if ok: good.setdefault(el,[]).append([fn,ps])
    print(f"{el:4s}{ps:12s}{fn[:25]:26s}{c['C11']:8.1f}{c['C12']:8.1f}{c['C44']:8.1f}"
          f"   {'accept' if ok else 'REJECT'}  (expt C11={ex_['C11']:.0f}, dev {dev:.0f} %"
          f"{', Born IHLAL' if not born else ''})")
json.dump(good, open("../standalone/nist_verified.json","w"), indent=1, sort_keys=True)
print(f"\nkabul edilen: {sum(len(v) for v in good.values())} potansiyel, {len(good)} element")
