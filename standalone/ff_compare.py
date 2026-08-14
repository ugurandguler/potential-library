#!/usr/bin/env python3
"""
Where this potential sits in the published distribution, per element.

The earlier comparison picked one EAM set per element and ran it here.  This
compares against every classical potential JARVIS-FF has for that element -
73 for copper, 62 for nickel - computed by NIST's own workflow rather than by
us.  Two objections disappear with it: that we ran EAM badly, and that we chose
a weak set to beat.

What is reported per element is the position rather than a winner.  The best
published potential, the median of them, and where our two parameter sets fall
in the ranking.  "Better than 60 of 73" is a statement about the distribution
and survives a referee; "we beat EAM" does not.

The comparison is by family where the name allows it.  MAU has no angular term
and belongs against EAM; UG carries one and belongs against MEAM and ADP.
Machine-learned potentials are excluded from the ranking - they are a different
kind of object with orders of magnitude more parameters, and beating or losing
to one says nothing about a seven-parameter analytic form.

    python ff_compare.py            # every element JARVIS-FF covers
    python ff_compare.py Cu Fe
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP = {"machine-learned"}


def main():
    ff = json.load(open(os.path.join(HERE, "jarvis_ff.json")))
    lib = json.load(open(os.path.join(HERE, "library.json")))
    els = sys.argv[1:] or sorted(ff)

    print("JARVIS-FF dagilimina karsi konum  (elastik RMS %, dusuk iyi)")
    print(f"{'el':4s}{'struct':5s}{'n':>4s}{'best':>8s}{'median':>9s}"
          f"{'OURS hard':>10s}{'rank':>8s}{'OURS taper':>11s}{'rank':>8s}")
    print("-" * 68)
    tot = {"hard": [], "tap": [], "best": [], "med": []}
    for el in els:
        rows = [r for r in ff.get(el, []) if r["family"] not in SKIP]
        if not rows or el not in lib:
            continue
        rms = sorted(r["rms"] for r in rows)
        n = len(rms)
        best, med = rms[0], rms[n // 2]
        v = lib[el]
        hard = v["rms"]
        tap = (v.get("tap") or {}).get("rms")
        rank_h = sum(1 for x in rms if x < hard) + 1
        rank_t = (sum(1 for x in rms if x < tap) + 1) if tap is not None else None
        print(f"{el:4s}{v['struct']:5s}{n:4d}{best:8.2f}{med:9.2f}"
              f"{hard:10.2f}{f'{rank_h}/{n+1}':>8s}"
              f"{(tap if tap is not None else float('nan')):11.2f}"
              f"{(f'{rank_t}/{n+1}' if rank_t else '-'):>8s}")
        tot["hard"].append(rank_h / (n + 1))
        if rank_t:
            tot["tap"].append(rank_t / (n + 1))
        tot["best"].append(best)
        tot["med"].append(med)

    import statistics as st
    print()
    print(f"median published potential  : {st.median(tot['med']):6.2f} %")
    print(f"median BEST published       : {st.median(tot['best']):6.2f} %")
    print(f"ours, hard cutoff,  median percentile : "
          f"{100*st.median(tot['hard']):5.0f}. (kucuk = ustte)")
    print(f"ours, tapered,      median percentile : "
          f"{100*st.median(tot['tap']):5.0f}.")


if __name__ == "__main__":
    main()
