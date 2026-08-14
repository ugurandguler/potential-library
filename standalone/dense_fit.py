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

#  how many distinct solutions to carry out of the search, for the
#  post-filter to choose from
N_POOL = int(os.environ.get("N_POOL", 40))

HERE = os.path.dirname(os.path.abspath(__file__))


def simplex(el, e, x0, step, budget=320):
    """Nelder-Mead on (m, gamma, s3, C); returns (best_score, best_record)"""
    store = {}

    def f(v):
        m_, g_, s_, C_ = v
        if (g_ < 0 or not (F.M_MIN <= m_ <= F.M_MAX) or s_ <= 0.05
                or abs(C_) > F.C_MAX):
            return 1e6
        key = tuple(round(t, 6) for t in v)
        if key in store:
            return store[key][0]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            r = F.evaluate(el, e, m_, g_, C_, s_)
        val = r["score"] if r else 1e6
        store[key] = (val, r)
        return val

    x0 = np.asarray(x0, dtype=float)
    sim = [x0] + [x0 + np.eye(4)[i]*step[i] for i in range(4)]
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
                for i in range(1, 5):
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
    #  Every distinct solution found, not only the winner.  The nudge test
    #  (lammps/jiggle_test.py) needs LAMMPS and cannot live inside the search -
    #  Nelder-Mead evaluates thousands of candidates and a minimisation each
    #  time would cost more than the whole fit.  So the search runs unchanged
    #  and the constraint is applied afterwards, which needs something to
    #  choose from: the best solution is often the one that fails, and the
    #  second or fifth may be sound at a small cost in RMS.
    #
    #  "Distinct" means the parameters differ, not the score: chromium and
    #  copper both have several solutions at identical RMS with entirely
    #  different parameters, and keeping one of each would defeat the purpose.
    pool = []
    KEY = ("m", "gamma", "D", "alpha", "r0", "C", "s3")

    def keep(res):
        nonlocal best, pool
        if not res[1]:
            return
        if res[0] < best[0]:
            best = res
        r = res[1]
        sig = tuple(round(float(r[k]), 6) for k in KEY)
        if not any(sig == p_[0] for p_ in pool):
            pool.append((sig, res[0], r))
            pool.sort(key=lambda x: x[1])
            del pool[N_POOL:]

    #  1. the warm start, so the answer can never be worse than the seed
    ws = F.warm_start(el)
    if ws:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            r = F.evaluate(el, e, ws[0], ws[1], ws[3], ws[2])
        if r:
            keep((r["score"], r))
        keep(simplex(el, e, [ws[0], ws[1], ws[2], ws[3]],
                     [1.5, 0.4, 0.5, max(1.0, 0.4*abs(ws[3]))]))

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
        step = [2.0, 0.5, 0.6, max(1.0, 0.35*abs(x0[3]))]
        keep(simplex(el, e, x0, step))

    #  the pool travels with the answer so the filter has something to work on
    if best[1] is not None:
        best[1]["pool"] = [dict(r, score=sc) for _sig, sc, r in pool]
    if best[1] is None:
        return None
    r = dict(best[1])
    r.pop("seed", None)
    r["struct"] = e["struct"]
    tol = 1e-3
    r["at_bound"] = [nm for nm, v_, lo, hi in
                     (("m", r["m"], F.M_MIN, F.M_MAX),
                      ("C", abs(r["C"]), -1.0, F.C_MAX),
                      ("gamma", r["gamma"], 0.0, 1e9))
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
