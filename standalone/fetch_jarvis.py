#!/usr/bin/env python3
"""
Pull phonon dispersions from JARVIS-DFT (NIST).

The `dft_3d` JSON dataset carries elastic tensors but no phonon bands, and the
documentation says "Gamma-point only", which is what the JSON exposes.  The
per-material XML is a different story: it holds the full band structure under
`phonon_bandstructure_frequencies`, sampled along the same high-symmetry path we
use.  So JARVIS is a usable third opinion, but only through the XML.

Checked against Materials Cloud MC3D before being trusted:

    Cu  JARVIS 238.7   MC3D 266.5 cm-1   ratio 0.90
    Nb  JARVIS 212.0   MC3D 224.5 cm-1   ratio 0.94

Six to ten per cent between two independent DFT calculations is ordinary.  Zinc
is the exception at ratio 0.53 - its symmetry behaviour is right (three acoustic
zeros at Gamma, E2g doubly degenerate, fourfold degeneracy at K), so it is not a
broken run, but two calculations disagreeing by a factor of two on the one
element whose c/a is 14 % off ideal is itself worth knowing.

Two traps.  `spg_number` in `dft_3d` is a **string**, so comparing it with an
int silently matches nothing.  And the frequency array is branch-major:
reshape to (nbranches, npoints), not the other way round.

Frequencies are cm^-1, the same unit the rest of the library stores.

    python fetch_jarvis.py            # every element we can match
    python fetch_jarvis.py Fe Ni Be

Writes jarvis_phonon.json.
"""
import json
import os
import re
import subprocess
import sys

import numpy as np

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "jarvis_phonon.json")
IDS = os.path.join(HERE, "jarvis_ids.json")
XML = "https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/{}.xml"
SG = {"fcc": "225", "bcc": "229", "hcp": "194"}
#  elements with no JARVIS entry in the phase we fit, confirmed by a full scan
#  of dft_3d rather than assumed; listed so a stale-cache check does not
#  re-download the dataset on every run because of them
NO_JARVIS = set()
THIN = 40                      # the path is sampled at ~10,000 points


def ids():
    #  The cache is keyed by nothing, so it silently goes stale when refdata
    #  grows: ten elements were added to the library and every one of them came
    #  back "no JVASP entry" because the file had been built when there were
    #  thirty.  Rebuild whenever an element is missing from it - the elements
    #  genuinely absent from JARVIS in our phase are recorded so that their
    #  absence is not mistaken for staleness a second time.
    if os.path.exists(IDS):
        cached = json.load(open(IDS))
        if not (set(refdata.ELEMENTS) - set(cached) - set(NO_JARVIS)):
            return cached
    from jarvis.db.figshare import data
    d = data("dft_3d")
    out = {}
    for r in d:
        el = r.get("formula", "")
        if el in refdata.ELEMENTS and \
                str(r.get("spg_number")) == SG[refdata.ELEMENTS[el]["struct"]]:
            out.setdefault(el, []).append(r["jid"])
    json.dump(out, open(IDS, "w"), indent=1)
    return out


def field(text, tag):
    m = re.search(rf"<{tag}>(.*?)</{tag}>", text, re.S)
    return m.group(1) if m else ""


def numbers(text, tag):
    return np.array([float(x) for x in
                     re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?",
                                field(text, tag))])


def main():
    els = sys.argv[1:] or sorted(refdata.ELEMENTS)
    jid = ids()
    store = json.load(open(OUT)) if os.path.exists(OUT) else {}

    print(f"{'el':4s}{'jid':>14s}{'branches':>10s}{'points':>8s}"
          f"{'max cm-1':>10s}{'min':>9s}")
    print("-" * 56)
    for el in els:
        cand = jid.get(el)
        if not cand:
            print(f"{el:4s}  no JVASP entry in our structure")
            continue
        got = None
        for j in cand:
            try:
                t = subprocess.run(["curl", "-s", "-L", "--max-time", "90",
                                    XML.format(j)], capture_output=True,
                                   text=True, timeout=120).stdout
            except Exception:                                  # noqa: BLE001
                continue
            f = numbers(t, "phonon_bandstructure_frequencies")
            d = numbers(t, "phonon_bandstructure_distances")
            if len(f) and len(d) and len(f) % len(d) == 0:
                got = (j, f, d, field(t, "phonon_bandstructure_labels"),
                       numbers(t, "phonon_bandstructure_label_points"))
                break
        if not got:
            print(f"{el:4s}  no band structure in the XML")
            continue
        j, f, d, labels, lpts = got
        nb = len(f) // len(d)
        F = f.reshape(nb, -1)
        store[el] = {
            "jid": j, "branches": nb,
            "source": "JARVIS-DFT (NIST), per-material XML",
            "x": [round(float(x), 5) for x in d[::THIN]],
            "f": [[round(float(v), 1) for v in F[b][::THIN]]
                  for b in range(nb)],
            "labels": [s.replace("\\Gamma", "G").strip()
                       for s in labels.strip().strip("'").split(",")],
            #  JARVIS gives the path length but not the q-vectors.  Keeping the
            #  distance of each label lets them be rebuilt by interpolating
            #  between the high-symmetry points, which is what allows our
            #  potential to be evaluated at THEIR q-points - without that the
            #  comparison would measure interpolation rather than physics.
            "label_x": [round(float(x), 6) for x in lpts],
        }
        print(f"{el:4s}{j:>14s}{nb:>10d}{len(d):>8d}{F.max():>10.1f}"
              f"{F.min():>9.1f}")

    tmp = OUT + ".tmp"
    json.dump(store, open(tmp, "w"))
    os.replace(tmp, OUT)
    print(f"\n{len(store)} elements in {OUT}")


if __name__ == "__main__":
    main()
