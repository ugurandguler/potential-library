#!/usr/bin/env python3
"""
Check the MEAM and ADP candidates before they are allowed to be baselines.

Two things make these harder than the EAM family and both are quiet failures.

MEAM needs TWO files and a four-argument pair_coeff, `* * library.meam El
param.meam El`, and the element name appears twice; get either wrong and LAMMPS
either errors or - worse - silently uses a different element from the library
file.  And NIST ships every MEAM entry under the same two names, library.meam
and <El>.meam, so three chromium entries would overwrite each other on
download.  They are stored here with a per-entry prefix for that reason.

ADP is a single file and behaves like EAM, but the pair_style is different and
the package is separate.
"""
import io, json, os, re, subprocess, sys
sys.path.insert(0, os.path.abspath("../standalone")); sys.path.insert(0, ".")
import numpy as np, latdyn as L, refdata, cellfile
from concurrent.futures import ThreadPoolExecutor

HOME=subprocess.run(["wsl","-e","bash","-lc","echo $HOME"],capture_output=True,text=True).stdout.strip()
LMP=f"{HOME}/lammps/src/lmp_serial"
def wsl(p):
    p=os.path.abspath(p).replace("\\","/"); return "/mnt/"+p[0].lower()+p[2:]

BODY = """units metal
boundary p p p
atom_style atomic
atom_modify map array
read_data cell.data
pair_style @PS@
@COEFF@
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
variable v equal c_born[${i}]*${bfac}
print "BORN ${i} ${v}"
variable v delete
next i
jump SELF lp
"""

def coeff(ps, el, files):
    if ps == "meam":
        lib = [f for f in files if "library" in f][0]
        par = [f for f in files if "library" not in f][0]
        return f"pair_coeff * * {lib} {el} {par} {el}"
    return f"pair_coeff * * {files[0]} {el}"

def one(c):
    el, ps, files = c["el"], c["ps"], c["files"]
    d=os.path.join("meamchk", f"{el}_{c['tag']}"); os.makedirs(d, exist_ok=True)
    e=refdata.ELEMENTS[el]
    cry=L.Crystal(e["struct"], e["a0"], e.get("c_over_a"), mass=refdata.MASSES[el])
    box,_=cellfile.orthogonal_cell(cry)
    rep=tuple(max(3,int(np.ceil(1.6*8.0/b))) for b in box)
    cellfile.write_data(os.path.join(d,"cell.data"), cry, refdata.MASSES[el], rep)
    open(os.path.join(d,"in.x"),"w").write(
        BODY.replace("@PS@",ps).replace("@COEFF@",coeff(ps,el,files)))
    cp=" ".join(f"cp ~/lammps/potentials/{f} {wsl(d)}/;" for f in files)
    subprocess.run(["wsl","-e","bash","-lc",
        f"{cp} cd {wsl(d)} && {LMP} -in in.x > o.txt 2>&1"], capture_output=True)
    p=os.path.join(d,"log.lammps")
    t=io.open(p,errors="ignore").read() if os.path.exists(p) else ""
    v={int(m.group(1)):float(m.group(2)) for m in re.finditer(r"BORN\s+(\d+)\s+([-\d.eE+]+)",t)}
    if len(v)<21:
        err=[l for l in t.splitlines() if "ERROR" in l]
        return c, None, (err[0][:58] if err else "born cikmadi")
    return c, {"C11":0.5*(v[1]+v[2]),"C12":v[7],"C13":0.5*(v[8]+v[12]),
               "C33":v[3],"C44":0.5*(v[4]+v[5])}, None

plan=json.load(open("../standalone/meam_plan.json"))
with ThreadPoolExecutor(max_workers=max(1,(os.cpu_count() or 4)-2)) as ex:
    res=list(ex.map(one, plan))
print(f"{'el':4s}{'style':6s}{'record':26s}{'C11':>8s}{'C12':>8s}{'C44':>8s}  verdict")
print("-"*80)
ok=[]
for c,v,err in res:
    el=c["el"]; ex_=refdata.ELEMENTS[el]["Cij"]; st=refdata.ELEMENTS[el]["struct"]
    if v is None:
        print(f"{el:4s}{c['ps']:6s}{c['id'][:25]:26s}{'':>24s}  DUSTU: {err}"); continue
    born=((v["C11"]>abs(v["C12"]) and (v["C11"]+v["C12"])*v["C33"]>2*v["C13"]**2 and v["C44"]>0)
          if st=="hcp" else (v["C11"]-v["C12"]>0 and v["C44"]>0 and v["C11"]+2*v["C12"]>0))
    dev=100*abs(v["C11"]-ex_["C11"])/ex_["C11"]
    good=born and dev<60
    if good: ok.append(c)
    print(f"{el:4s}{c['ps']:6s}{c['id'][:25]:26s}{v['C11']:8.1f}{v['C12']:8.1f}{v['C44']:8.1f}"
          f"  {'accept' if good else 'REJECT'} (expt {ex_['C11']:.0f}, {dev:.0f} %"
          f"{', Born IHLAL' if not born else ''})")
json.dump(ok, open("../standalone/meam_verified.json","w"), indent=1)
print(f"\nkabul: {len(ok)}")
