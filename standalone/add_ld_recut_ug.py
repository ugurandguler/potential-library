#!/usr/bin/env python3
"""
Merge the angular re-cut candidate's dispersion into library.json.

angular/export_recut_ug_ld.py computes it on the other side of the boundary -
the angular latdyn cannot be imported here - and writes recut_ug_ld.json.
This only reads that file, exactly as add_ug_overlay reads ug_results.json.
If the file is absent the step is a no-op and the page draws the pair
candidate alone, which is the correct behaviour before any angular run.

    python add_ld_recut_ug.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "angular", "recut_ug_ld.json")


def main():
    path = os.path.join(HERE, "library.json")
    if not os.path.exists(SRC):
        print(f"{SRC} missing - skipped")
        return
    src = json.load(open(SRC))
    lib = json.load(open(path))
    n, skipped = 0, []
    for el, ld in sorted(src.items()):
        rec = (lib.get(el) or {}).get("rc_ug")
        if not isinstance(rec, dict) or "m" not in rec:
            skipped.append(el)
            continue
        rec["ld"] = ld
        n += 1
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print(f"{n} angular candidate arms carry a dispersion"
          + (f"; skipped {' '.join(skipped)}" if skipped else ""))
    print(f"merged into {path}")


if __name__ == "__main__":
    main()
