#!/usr/bin/env python3
"""
Give the switched arm, `tap`, its own 0 K dispersion.

The finite-temperature panel measures `<El>_taper.ugur` - the switched record,
the only one that can be run at temperature - and draws it against a dashed
"0 K, harmonic" curve.  That dashed curve was `d.ld.std`, the TOP-LEVEL
hard-cut record, because `tap` had no dispersion of its own.  So the picture
set a switched run against a hard-cut calculation, and the gap a reader took
for temperature was partly truncation.  It is not small: hard and switched
differ at 0 K by more than the finite-temperature difference itself in 15 of
the 24 elements (yttrium 8.4 against 19.8 %, titanium 12.1 against 19.4 %).
The table beside the picture was always right - make_finiteT scores `tap` -
only the drawn reference was the wrong record.

This is add_ld_recut's `build`, unchanged, applied to `tap`: the same
Setyawan-Curtarolo path, the same NQ_PATH points and `len`, so the curve stays
in register with every other panel.  Checked before it was trusted: the same
function applied to the top-level record reproduces its stored `ld.std`
exactly (Cu, Nb, Mg, largest difference 0.0 cm^-1).

Every element with a `tap` record gets one, not only those in the finite-T
study, so that a run added later does not have to remember this step.

`tap.ld` is read by the finite-temperature panel only.  The parameter-set
selector's "MAU" resolves to the top-level record, not to `tap`, and the
thermodynamics panel drops `tap` from its extra curves explicitly.

    python add_ld_tap.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import refdata                                       # noqa: E402
from add_ld_recut import build                       # noqa: E402


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    n = 0
    print(f"{'el':4s}{'max cm-1':>10s}{'most neg.':>12s}{'root max':>10s}")
    print("-" * 36)
    for el in sorted(lib):
        v = lib[el]
        rec = v.get("tap") if isinstance(v, dict) else None
        if not isinstance(rec, dict) or "m" not in rec:
            continue
        ld = build(el, refdata.ELEMENTS[el], rec)
        rec["ld"] = ld
        lo = min(min(min(b) for b in p["branches"]) for p in ld["std"])
        n += 1
        print(f"{el:4s}{ld['maxfreq']:10.1f}{lo:12.1f}"
              f"{(v.get('ld') or {}).get('maxfreq', float('nan')):10.1f}")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print()
    print(f"{n} switched arms now carry a dispersion -> {path}")


if __name__ == "__main__":
    main()
