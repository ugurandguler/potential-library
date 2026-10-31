#!/usr/bin/env python3
"""Structure gate over every shipped record.  build | collect"""
import os, re, sys, json, glob, shutil
#  where the gate runs (default: this directory) and where the expansion
#  sweep's run directories are, whose cell.data files the gate starts from
G = os.environ.get("GATE_DIR", os.path.dirname(os.path.abspath(__file__))) + "/"
NPT = os.environ.get("GATE_CELLS", os.path.join(G, "..", "nptruns")) + "/"
M = json.load(open(G + "meta.json"))
ARMS = [("_force_taper.ugur", "tap_force"), ("_nudge_taper.ugur", "nudge"), ("_disp.ugur", "disp"),
        ("_taper.ugur.ang", "tap_ug"), ("_taper.ugur", "tap"), ("_recut.ugur.ang", "rc_ug"),
        ("_recut.ugur", "rc"), (".ugur.ang", "hard_ug"), (".ugur", "hard")]
def arm_of(f):
    el = re.match(r"([A-Z][a-z]?)", f).group(1)
    for suf, a in ARMS:
        if f == el + suf: return el, a
    return el, None
def build():
    jobs = []
    for f in sorted(os.listdir(G + "pots")):
        el, a = arm_of(f)
        if a is None: print("skip", f); continue
        cells = sorted(glob.glob(NPT + f"{el}_tap_[0-9]*K/cell.data"))
        if not cells: print("NO CELL", el, a); continue
        m = M[el]; d = G + f"runs/{el}_{a}/"
        os.makedirs(d, exist_ok=True)
        shutil.copy(cells[0], d + "cell.data"); shutil.copy(G + "pots/" + f, d + f)
        st = "ugur/ang" if f.endswith(".ang") else "ugur"
        open(d + "pair.in", "w").write(f"pair_style {st}\npair_coeff * * {f} {el}\n")
        jobs.append(f"{d} {m['T']} {m['rc']} {m['ct']}")
    open(G + "jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
def grab(t, k):
    m = re.findall(r"^%s .*" % k, t, re.M)
    return dict(re.findall(r"(\w+)=([-\d.eE+]+)", m[-1])) if m else None
def collect():
    out = {}
    for d in sorted(glob.glob(G + "runs/*/")):
        name = d.rstrip("/").split("/")[-1]; el, a = name.split("_", 1)
        t = open(d + "log.gate").read() if os.path.exists(d + "log.gate") else ""
        e0, md, q = grab(t, "E0"), grab(t, "MD"), grab(t, "QUENCH")
        r = dict(el=el, arm=a, T=M[el]["T"], struct=M[el]["struct"])
        if not (e0 and md and q):
            err = re.findall(r"ERROR.*", t); r["verdict"] = "NORUN"; r["err"] = err[0][:80] if err else "unfinished"
        else:
            N = float(e0["N"]); kT = 8.617333e-5 * M[el]["T"]
            r.update(E0=float(e0["pe_atom"]), f0=float(e0["nref"]) / N, fmd=float(md["nref"]) / N,
                     fq=float(q["nref"]) / N, dE_md=1e3 * (float(md["pe_atom"]) - float(e0["pe_atom"])),
                     kT15=1.5e3 * kT, dE_q=1e3 * (float(q["pe_atom"]) - float(e0["pe_atom"])))
            bad = r["dE_q"] < -2.0 or r["fq"] < 0.9
            r["verdict"] = "FAIL" if bad else ("PASS" if r["f0"] > 0.99 else "CHECK")
        out[name] = r
    json.dump(out, open(G + "gate.json", "w"), indent=1, sort_keys=True)
    from collections import Counter
    print(Counter(v["verdict"] for v in out.values()))
    for k, v in sorted(out.items()):
        if v["verdict"] != "PASS":
            print("%-16s %-5s %s" % (k, v["verdict"], {x: round(v[x], 3) for x in ("f0", "fmd", "fq", "dE_md", "kT15", "dE_q") if x in v} if "f0" in v else v.get("err")))
if __name__ == "__main__":
    {"build": build, "collect": collect}[sys.argv[1]]()
