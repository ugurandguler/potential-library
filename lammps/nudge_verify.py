#!/usr/bin/env python3
"""
Re-measure the nudge test on the candidate records, here, against several
displacement directions.

One displacement direction is one sample, and iron proved it.  The filter
originally certified an iron solution on the strength of a single direction;
re-measured along five, the best candidate it could offer returned to its
lattice under seed 87287 and kept eight to ten thousand times the displacement
under the other four.  One in five.  A record that holds its lattice does so
whichever way it is pushed - a record that holds it along a line does not hold
it, and the only way to tell the two apart is to ask more than once.

The filter now tests five directions itself, and this stays as the independent
check on it: same measurement, different code path, run wherever the library
is being assembled rather than wherever the fit happened.  Two of tonight's
corrections came from exactly that, so it earns its place.

With five directions the separation is unambiguous: the sound candidates keep
0.01-0.10 on every one, the shipped records keep 8000-16000 on every one, and
nothing sits between.

    python nudge_verify.py                 # every element in nudge_picked.json
    python nudge_verify.py Fe --seeds 9
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nudge_filter as NF      # noqa: E402

#  arbitrary and fixed, so the measurement is reproducible
SEEDS = (87287, 11311, 40213, 90007, 55501, 33199, 71723, 20261, 60649)


def main():
    argv = sys.argv[1:]
    nseed = 5
    if "--seeds" in argv:
        i = argv.index("--seeds")
        nseed = int(argv[i + 1])
        del argv[i:i + 2]
    els = [a for a in argv if not a.startswith("--")]

    picks = json.load(open(os.path.join(HERE, "nudge_picked.json")))
    lib = json.load(open(os.path.join(os.path.dirname(HERE),
                                      "standalone", "library.json")))
    seeds = SEEDS[:nseed]
    out = {}
    #  The seed is passed as an argument, not patched into the input template.
    #  It used to be `NF.IN.replace("87287", str(s))`, which worked only while
    #  the template carried that number literally; once nudge_filter grew its
    #  own seed parameter the replace matched nothing, every "direction" was
    #  the same direction, and five identical numbers were reported as five
    #  independent agreements.  A verification that cannot fail is not one.
    if "{seed}" not in NF.IN:
        raise SystemExit("nudge_filter.IN artik {seed} tasimiyor - "
                         "the check may not be changing direction")
    print(f"{'el':4s}{'record':8s}" + "".join(f"{s:>11d}" for s in seeds)
          + "   passing")
    print("-" * (12 + 11 * len(seeds) + 9))
    for el in sorted(els or picks):
        for lab, rec in (("shipped", lib[el]["tap"]), ("candidate", picks[el])):
            keeps, des = [], []
            for s in seeds:
                v = NF.test(el, rec, "ugur", f"verify{s}_{lab}", seed=s)
                keeps.append(v["keep"] if v else None)
                des.append(v["dE"] if v else None)
            good = [k for k in keeps if k is not None]
            npass = sum(1 for k in good if k < NF.KEEP_MAX)
            print(f"{el:4s}{lab:8s}"
                  + "".join(f"{k:11.2f}" if k is not None else f"{'-':>11s}"
                            for k in keeps)
                  + f"{npass:>6d}/{len(seeds)}")
            out[f"{el}|{lab}"] = {
                "seeds": list(seeds), "keep": keeps, "dE": des,
                "n_pass": npass, "n": len(seeds),
                #  every seed, not a majority: a lattice that holds only
                #  sometimes does not hold
                "ok": bool(good and npass == len(seeds))}
        print()
    p = os.path.join(HERE, "nudge_verify.json")
    old = {}
    if os.path.exists(p):
        try:
            old = json.load(open(p))
        except Exception:
            old = {}
    old.update(out)
    tmp = p + ".tmp"
    json.dump(old, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, p)
    print(f"{len(out)} measurements written; file now holds {len(old)}  -> {p}")


if __name__ == "__main__":
    main()
