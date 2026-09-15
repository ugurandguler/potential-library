#!/usr/bin/env python3
"""
Attach the electronic heat-capacity coefficient to each element's record.

Writes lib[el]["gamma_e"] (mJ mol^-1 K^-2) from refdata_electronic.GAMMA, for
the elements that table covers, and removes the key from any element it does
not - so a value dropped from the source cannot linger in the library.  The
page reads it to show gamma*T beside the lattice term; see
refdata_electronic.py for the source and for why gamma*T is an upper bound at
room temperature.

    python add_electronic.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import refdata_electronic as RE                     # noqa: E402


def main():
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    have, missing = [], []
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict) or "Cp298" not in v:
            continue
        if el in RE.GAMMA:
            v["gamma_e"] = RE.GAMMA[el]
            have.append(el)
        else:
            v.pop("gamma_e", None)
            missing.append(el)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(lib, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    print("gamma_e written for %d elements; not in the source: %s"
          % (len(have), " ".join(missing) or "none"))


if __name__ == "__main__":
    main()
