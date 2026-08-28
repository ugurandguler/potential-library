#!/usr/bin/env python3
"""Fold the finite-temperature CURVES into finiteT.json.

The curve was generated on ld.std's own segments, in ld.std's order, so the
two share one x axis without either being re-derived.  That is the whole
reason emit_ftcurve.py read the path out of the library instead of rebuilding
it - but "was generated from" is a claim about a script that ran on another
machine three steps ago, so it is CHECKED here rather than assumed: same
number of segments, same endpoints, same order.  An element that fails the
check is dropped and falls back to the number-only panel it had before.

THE OPTIC BRANCHES NEAR GAMMA ARE RESHAPED HERE.  ftcurve_remote.py holds them
flat at their value on the substitution shell, which was the right call for a
scored number - no scored point sits inside that radius - and the wrong one
for a drawn curve.  Gamma is a mesh point, so its frequencies are measured
rather than interpolated: titanium's highest branch is 12.126 THz there and
10.312 at the first point outside the substituted cell along Gamma-A, and
both of those are measurements.  Holding flat throws the Gamma value away and
opens a step of nearly 2 THz in the middle of the picture.

An optic branch is flat at Gamma to first order - it has no term linear in q -
so the leading behaviour is

    omega(q) = omega(0) + [omega(shell) - omega(0)] (q/q_shell)^2

which is exact at both ends, has zero slope at Gamma as it must, and closes
the step.  It is still not a measurement and the page still draws it dotted.

The acoustic branches are left exactly as the cluster produced them: there
omega = v|q| is already exact at both ends and there is nothing to close.

Reads ftcurve.json fresh every time, so re-running after a new harvest is
safe - the reshaping is never applied twice to its own output.

    python merge_ftcurve.py           # after fetching ftcurve.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def reshape_optic(el, f, subs, npt, mesh):
    """q^2 through Gamma and the shell, for the optic branches only.

    Returns the number of points touched.  A cubic cell has three branches and
    all of them are acoustic, so it returns zero and nothing moves.
    """
    if len(f[0]) <= 3 or not subs:
        return 0
    qs = [[float(x) for x in ln.split()]
          for ln in open(os.path.join(HERE, "ftq_%s.txt" % el))]
    assert len(qs) == len(f), "%s: %d q, %d satir" % (el, len(qs), len(f))
    step = 2.0 / mesh
    #  the substituted points of one segment share a Gamma end and a shell
    #  value; group them so each segment is closed with its own pair
    by_seg = {}
    for i in subs:
        by_seg.setdefault(i // npt, []).append(i)
    n = 0
    for si, idx in by_seg.items():
        gam = [j for j in range(si * npt, (si + 1) * npt)
               if max(abs(v) for v in qs[j]) < 1e-9]
        if not gam:
            continue                      # not a Gamma-adjacent segment
        #  COPIES, both of them.  f is rewritten in place below, and idx[0] is
        #  one of the rows being rewritten - binding the row itself made the
        #  shell value follow the first assignment and every later point was
        #  then closed against a target that had already moved.
        f0 = list(f[gam[0]])
        shell = list(f[idx[0]])           # held flat, so any of them will do
        for i in idx:
            r = max(abs(v) for v in qs[i]) / step
            for b in range(3, len(f[i])):
                f[i][b] = round(f0[b] + (shell[b] - f0[b]) * r * r, 3)
            n += 1
    return n


def main():
    lib = json.load(open(os.path.join(HERE, "library.json")))
    ft = json.load(open(os.path.join(HERE, "finiteT.json")))
    cur = json.load(open(os.path.join(HERE, "ftcurve.json")))
    ix = json.load(open(os.path.join(HERE, "ftq_index.json")))

    out, bad = {}, []
    for el in sorted(cur):
        if el not in ft["rows"]:
            bad.append((el, "finiteT satiri yok"))
            continue
        std = ((lib.get(el) or {}).get("ld") or {}).get("std")
        segs, npt = ix[el]["segs"], ix[el]["npt"]
        if not std or len(std) != len(segs):
            bad.append((el, "segment sayisi %s vs %s"
                        % (len(std or []), len(segs))))
            continue
        if any((s["a"], s["b"]) != (t["a"], t["b"])
               for s, t in zip(std, segs)):
            bad.append((el, "yol sirasi farkli"))
            continue
        f = [list(r) for r in cur[el]["f"]]
        if len(f) != len(segs) * npt:
            bad.append((el, "%d q, %d beklenen" % (len(f), len(segs) * npt)))
            continue
        #  ONE RUN, or none.  binary() on the cluster takes the largest box it
        #  finds, and fix phonon writes its binary as it goes - so caesium's
        #  curve came off the 21^3 job while it was still running, next to a
        #  table of numbers from the finished 8^3 one.  A panel whose picture
        #  and whose numbers are two different measurements is worse than a
        #  panel with no picture, and the mesh is what tells them apart.
        want = int(str(ft["rows"][el]["mesh"]).split("x")[0])
        if cur[el]["mesh"] != want:
            bad.append((el, "egri %d^3, tablo %s - ayri kosular"
                        % (cur[el]["mesh"], ft["rows"][el]["mesh"])))
            continue
        subs = cur[el]["sub"]
        nopt = reshape_optic(el, f, subs, npt, cur[el]["mesh"])
        out[el] = {"f": f, "sub": subs, "npt": npt, "mesh": cur[el]["mesh"]}
        print("%-3s %2d x %-2d = %4d q, %2d degistirilmis, %d dal%s"
              % (el, len(segs), npt, len(f), len(subs), len(f[0]),
                 ", %d optik yeniden sekillendirildi" % nopt if nopt else ""))

    for el, why in bad:
        print("%-3s ATLANDI: %s" % (el, why))

    ft["curve"] = out
    p = os.path.join(HERE, "finiteT.json")
    json.dump(ft, open(p, "w"), separators=(",", ":"), sort_keys=True)
    print("\n%d egri, %s  %d bayt" % (len(out), p, os.path.getsize(p)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
