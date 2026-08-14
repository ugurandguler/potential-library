#!/usr/bin/env python3
"""
Intensive multi-start search for one element.  Written for the cluster: one
element per core, results merged afterwards.

The plain fit does a coarse grid then a single simplex.  That is enough for most
elements but leaves V, W, Nb and Pt worse than their seed, because the
constraint solver converges to a different branch depending on where it starts
and a single simplex cannot escape it.  Here the same objective is attacked from
many independent starting points, so the reported optimum is the best of a large
sample rather than the first one found.

    python dense_fit.py Pd [n_restarts] [seed]

Writes dense_<el>.json.  Never returns something worse than the warm start.
"""
import json, math, os, random, sys, time
import numpy as np
import fit as F
import refdata

HERE = os.path.dirname(os.path.abspath(__file__))


#  Two extra search directions, lam2 and lam4, the Legendre weights of the
#  angular factor.  Both zero reproduces the published phi3 exactly, so the
#  6-parameter search contains the 4-parameter one and cannot do worse - as
#  long as the starting points include lam = 0, which run() ensures.
#  Set ANGULAR = 0 to search the published four parameters only.  Needed for a
#  fair comparison: the archived baseline was fitted WITH the dynamical
#  stability screen and the angular runs cannot use it (no analytic derivatives
#  for the angular term), so part of any gain would just be the missing
#  constraint.  Rerunning the published form here, screen off and everything
#  else identical, isolates what the angular term actually buys.
ANGULAR = os.environ.get("ANGULAR", "1") != "0"

LAM2_MIN, LAM2_MAX = -0.95, 1.90
LAM4_MIN, LAM4_MAX = -1.20, 1.20
NDIM = 6 if ANGULAR else 4


def _angular_ok(lam2, lam4, n=41):
    """the angular factor must stay positive, or phi3 changes sign"""
    c = np.linspace(-1.0, 1.0, n)
    h = (1.0 + lam2*0.5*(3*c**2 - 1)
         + lam4*0.125*(35*c**4 - 30*c**2 + 3))
    return h.min() > 0.05


def simplex(el, e, x0, step, budget=520):
    """Nelder-Mead on (m, gamma, s3, C, lam2, lam4)"""
    store = {}

    def f(v):
        if ANGULAR:
            m_, g_, s_, C_, l2_, l4_ = v
        else:
            (m_, g_, s_, C_), l2_, l4_ = v, 0.0, 0.0
        if (g_ < 0 or not (F.M_MIN <= m_ <= F.M_MAX) or s_ <= 0.05
                or abs(C_) > F.C_MAX
                or not (LAM2_MIN <= l2_ <= LAM2_MAX)
                or not (LAM4_MIN <= l4_ <= LAM4_MAX)
                or not _angular_ok(l2_, l4_)):
            return 1e6
        key = tuple(round(t, 6) for t in v)
        if key in store:
            return store[key][0]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            r = F.evaluate(el, e, m_, g_, C_, s_, lam2=l2_, lam4=l4_)
        val = r["score"] if r else 1e6
        store[key] = (val, r)
        return val

    x0 = np.asarray(x0, dtype=float)
    sim = [x0] + [x0 + np.eye(NDIM)[i]*step[i] for i in range(NDIM)]
    val = [f(s) for s in sim]
    for _ in range(budget):
        order = np.argsort(val)
        sim = [sim[i] for i in order]; val = [val[i] for i in order]
        if abs(val[-1] - val[0]) < 1e-10 * max(1.0, abs(val[0])):
            break
        cen = np.mean(sim[:-1], axis=0)
        xr = cen + (cen - sim[-1]); fr = f(xr)
        if fr < val[0]:
            xe = cen + 2.0*(cen - sim[-1]); fe = f(xe)
            sim[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            sim[-1], val[-1] = xr, fr
        else:
            xc = cen + 0.5*(sim[-1] - cen); fc = f(xc)
            if fc < val[-1]:
                sim[-1], val[-1] = xc, fc
            else:
                for i in range(1, NDIM+1):
                    sim[i] = sim[0] + 0.5*(sim[i] - sim[0])
                    val[i] = f(sim[i])
    best = (1e6, None)
    for (v, r) in store.values():
        if r and v < best[0]:
            best = (v, r)
    return best


def run(el, n_restarts=48, rseed=0):
    e = refdata.ELEMENTS[el]
    rng = random.Random(rseed)
    best = (1e6, None)
    t0 = time.time()

    def keep(res):
        nonlocal best
        if res[1] and res[0] < best[0]:
            best = res

    #  1. the warm start, so the answer can never be worse than the seed
    ws = F.warm_start(el)
    if ws:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            r = F.evaluate(el, e, ws[0], ws[1], ws[3], ws[2])
        if r:
            keep((r["score"], r))
        #  from the published-form optimum with the angular weights at zero,
        #  which guarantees the 6-parameter answer is at least as good
        keep(simplex(el, e, [ws[0], ws[1], ws[2], ws[3]] + ([0.0, 0.0] if ANGULAR else []),
                     [1.5, 0.4, 0.5, max(1.0, 0.4*abs(ws[3]))] + ([0.3, 0.3] if ANGULAR else [])))

    #  2. a systematic grid, coarse but wide
    for m in (1.4, 2.2, 3.2, 4.6, 6.5, 9.0, 12.5, 16.0, 19.5):
        for g in (0.0, 0.4, 0.8, 1.2, 1.8, 2.5):
            for s3 in (0.4, 0.8, 1.5, 3.0, 5.0):
                for C in (0.0, -0.8, 0.8, -3.0, 3.0, -12.0, 12.0, -45.0, 45.0):
                    with np.errstate(over="ignore", invalid="ignore",
                                     divide="ignore"):
                        r = F.evaluate(el, e, m, g, C, s3)
                    if r:
                        keep((r["score"], r))

    #  3. many independent simplex restarts from random points
    for k in range(n_restarts):
        x0 = [rng.uniform(F.M_MIN, F.M_MAX),
              rng.uniform(0.0, 2.6),
              math.exp(rng.uniform(math.log(0.15), math.log(6.0))),
              rng.choice([1, -1]) * math.exp(rng.uniform(math.log(0.05),
                                                         math.log(F.C_MAX)))]
        #  half the restarts begin with no angular term and half with a random
        #  one, so neither region of the space is starved
        if ANGULAR:
            x0 += ([0.0, 0.0] if k % 2 == 0 else
                   [rng.uniform(-0.9, 1.5), rng.uniform(-1.0, 1.0)])
        step = [2.0, 0.5, 0.6, max(1.0, 0.35*abs(x0[3]))]
        if ANGULAR:
            step += [0.35, 0.35]
        keep(simplex(el, e, x0, step))

    if best[1] is None:
        return None
    r = dict(best[1])
    r.pop("seed", None)
    r["struct"] = e["struct"]
    tol = 1e-3
    r["at_bound"] = [nm for nm, v_, lo, hi in
                     (("m", r["m"], F.M_MIN, F.M_MAX),
                      ("C", abs(r["C"]), -1.0, F.C_MAX),
                      ("gamma", r["gamma"], 0.0, 1e9),
                      ("lam2", r.get("lam2", 0.0), LAM2_MIN, LAM2_MAX),
                      ("lam4", r.get("lam4", 0.0), LAM4_MIN, LAM4_MAX))
                     if abs(v_-lo) < tol or abs(v_-hi) < tol]
    r["search"] = {"restarts": n_restarts, "seconds": round(time.time()-t0, 1)}
    return r


if __name__ == "__main__":
    el = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 48
    rs = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    res = run(el, n, rs)
    out = os.path.join(HERE, f"dense_{el}_{rs}.json")
    if res is None:
        json.dump({"element": el, "failed": True}, open(out, "w"))
        print(f"{el}: FAILED")
    else:
        json.dump({el: res}, open(out, "w"), indent=1, sort_keys=True)
        print(f"{el}: rms={res['score']*100:.2f}%  gamma={res['gamma']:.4f}  "
              f"m={res['m']:.3f}  C={res['C']:+.3f}  "
              f"{res['search']['seconds']:.0f}s")
