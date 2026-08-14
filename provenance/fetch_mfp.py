#!/usr/bin/env python3
"""
Pull phonon dispersions from the Chiral Phonon Materials Database
(materialsfingerprint.com, Zhejiang University; Nature Physics, 2026).

The catalogue is built to classify phonon angular momentum, but it computed a
full dispersion for all 11,614 compounds to do so, and 1,075 of those are
elemental.  It reaches 25 of our 30 - including **Zn and Cd**, the two the
potential cannot fit at all, for which an independent dispersion is worth
having even though there is no fit to compare it against.

Three undocumented endpoints, read out of the site's own JavaScript:

    POST /api/search                       the whole catalogue, ignores filters
    GET  /api/material?id=<id>             formula, space group, attachments
    GET  /api/attachment?kind=dat&id=<id>  the dispersion itself

The .dat is columns: path length, frequency in **meV**, then angular-momentum
and group-velocity components.  meV, not THz and not cm^-1: 1 meV = 0.2417989
THz = 8.065544 cm^-1.

Several ICSD entries exist per element - 27 for zinc alone - differing in the
source structure.  The one taken is whichever has the most atoms in the
primitive cell matching what we fit (1 for cubic, 2 for hcp) and, among those,
the first listed; a lattice check against refdata would be better and is worth
adding if any of these are ever used quantitatively.

    python fetch_mfp.py            # everything matched
    python fetch_mfp.py Zn Cd

Writes mfp_phonon.json.
"""
import json
import os
import subprocess
import sys

import refdata

API = "https://materialsfingerprint.com/api"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "mfp_phonon.json")
SG = {"fcc": 225, "bcc": 229, "hcp": 194}
NAT = {"fcc": 1, "bcc": 1, "hcp": 2}
MEV_CM1 = 8.065544
THIN = 4                      # the .dat is dense; keep every fourth point


def get(url, timeout=90):
    p = subprocess.run(["curl", "-s", "-L", "--max-time", str(timeout), url],
                       capture_output=True, text=True, timeout=timeout + 30)
    return p.stdout


def catalogue():
    path = os.path.join(HERE, "mfp_all.json")
    if not os.path.exists(path):
        raise SystemExit("mfp_all.json missing; fetch it with\n"
                         "  curl -X POST -H 'Content-Type: application/json' "
                         "-d '{}' " + API + "/search > mfp_all.json")
    return json.load(open(path))


def main():
    els = sys.argv[1:] or sorted(refdata.ELEMENTS)
    cat = catalogue()
    store = json.load(open(OUT)) if os.path.exists(OUT) else {}

    print(f"{'el':4s}{'id':>16s}{'atoms':>7s}{'points':>8s}"
          f"{'branches':>10s}{'max cm-1':>10s}")
    print("-" * 56)
    for el in els:
        e = refdata.ELEMENTS[el]
        want_sg, want_nat = SG[e["struct"]], NAT[e["struct"]]
        cand = [r for r in cat
                if r.get("cnt_elements") == 1
                and r.get("chemistry_formula") == el
                and r.get("space_group") == want_sg
                and r.get("cnt_atoms") == want_nat]
        if not cand:
            print(f"{el:4s}  no entry in the right structure")
            continue
        rec = cand[0]
        raw = get(f"{API}/attachment?kind=dat&id={rec['id']}")
        rows = []
        for line in raw.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 2:
                continue
            try:
                rows.append((float(f[0]), float(f[1])))
            except ValueError:
                continue
        if not rows:
            print(f"{el:4s}  {rec['id']}  no numeric rows")
            continue
        #  the file lists one branch after another down the same path, so a
        #  branch ends wherever the path length stops increasing
        branches, cur = [], []
        for x, v in rows:
            if cur and x < cur[-1][0] - 1e-9:
                branches.append(cur)
                cur = []
            cur.append((x, v))
        if cur:
            branches.append(cur)
        store[el] = {
            "id": rec["id"], "icsd": rec.get("ICSD"), "mpid": rec.get("mpid"),
            "space_group": rec["space_group"], "natoms": rec["cnt_atoms"],
            "source": ("Chiral Phonon Materials Database, "
                       "materialsfingerprint.com (Zhejiang University)"),
            "x": [round(x, 5) for x, _ in branches[0][::THIN]],
            "f": [[round(v * MEV_CM1, 1) for _, v in b[::THIN]]
                  for b in branches],
        }
        top = max(v for b in branches for _, v in b) * MEV_CM1
        print(f"{el:4s}{rec['id']:>16s}{rec['cnt_atoms']:>7d}"
              f"{len(branches[0]):>8d}{len(branches):>10d}{top:>10.1f}")

    tmp = OUT + ".tmp"
    json.dump(store, open(tmp, "w"))
    os.replace(tmp, OUT)
    print(f"\n{len(store)} elements in {OUT}")


if __name__ == "__main__":
    main()
