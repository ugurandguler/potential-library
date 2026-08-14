#!/usr/bin/env python3
"""
Multi-species energy, as an independent reference for `pair_style ugur`.

Everything else in this project is validated by two implementations agreeing.
The alloy path had no second implementation, so what had been checked was that
the code is self-consistent - the twin test, the regression - and not that it
computes the intended potential on a mixed structure.  This is the second
implementation.

It reads **the same file LAMMPS reads**, on purpose.  A reference that took its
parameters from library.json and applied its own mixing would be checking the
mixing rule twice and the kernel not at all; reading the file means the only
thing being compared is what the two codes do with identical numbers.

The conventions it has to match, all of them decisions rather than mathematics
and all recorded in lammps/ALLOYS.md:

  * two-body parameters come from the (i,j,j) entry, Tersoff's rule
  * the switch acts per leg, so S_AB(r1) uses rcut3 of the A-B pair
  * a triple takes C and alpha3 from its own (i,j,k) entry and the radial
    shape - D, m, r0, gamma - from the CENTRE's own (i,i,i) entry, because
    (A,B,C) and (A,C,B) carry different pair columns and a triplet must not
    depend on which leg the loop reached first

    python alloy.py CuNi.ugur.alloy Cu Ni        # self-test on a written cell
"""
import math
import os
import sys

import numpy as np

import latdyn as L

HERE = os.path.dirname(os.path.abspath(__file__))
FIELDS = ("m", "D", "alpha", "r0", "gamma", "C", "alpha3",
          "rcut2", "rcut3", "taper", "lam2", "lam4")


class AlloyPotential:
    """the Tersoff-format file, parsed, with the lookup rules applied"""

    def __init__(self, entries, elements):
        self.elements = list(elements)
        self.idx = {e: i for i, e in enumerate(self.elements)}
        n = len(self.elements)
        self.tri = {}
        for (a, b, c), v in entries.items():
            if a in self.idx and b in self.idx and c in self.idx:
                self.tri[(self.idx[a], self.idx[b], self.idx[c])] = v
        missing = [(a, b, c) for a in range(n) for b in range(n)
                   for c in range(n) if (a, b, c) not in self.tri]
        if missing:
            raise ValueError(f"dosyada eksik ucluler: {missing[:4]}")
        tapers = {v["taper"] for v in self.tri.values()}
        if len(tapers) > 1:
            raise ValueError(f"the file carries different taper values: {tapers}")
        self.rcut2 = max(v["rcut2"] for v in self.tri.values())
        self.rcut3 = max(v["rcut3"] for v in self.tri.values())

    @classmethod
    def from_file(cls, path, elements):
        entries = {}
        for line in open(path):
            line = line.split("#")[0].strip()
            if not line:
                continue
            t = line.split()
            if len(t) < 3 + 10:
                continue
            vals = [float(x) for x in t[3:3 + 12]] if len(t) >= 15 else \
                   [float(x) for x in t[3:3 + 10]] + [0.0, 0.0]
            entries[(t[0], t[1], t[2])] = dict(zip(FIELDS, vals))
        if not entries:
            raise ValueError(f"{path}: no readable line")
        return cls(entries, elements)

    def pair(self, i, j):
        """(i,j,j) - where the two-body parameters live"""
        return self.tri[(i, j, j)]

    def triple(self, i, j, k):
        """C and alpha3 from the triple, radial shape from the centre"""
        p = dict(self.pair(i, i))
        t = self.tri[(i, j, k)]
        p["C"] = t["C"]
        p["alpha3"] = t["alpha3"]
        p["lam2"] = t["lam2"]
        p["lam4"] = t["lam4"]
        return p

    #  ---- the two terms, written to mirror latdyn.Potential exactly --------
    @staticmethod
    def _as_pot(p, rcut3=None):
        q = L.Potential(m=p["m"], D=p["D"], alpha=p["alpha"], r0=p["r0"],
                        gamma=p["gamma"], C=p["C"], alpha3=p["alpha3"],
                        rcut2=p["rcut2"],
                        rcut3=(p["rcut3"] if rcut3 is None else rcut3),
                        taper=(p["taper"] if p["taper"] > 0 else None))
        return q

    def phi2(self, i, j, r):
        return self._as_pot(self.pair(i, j)).phi2(np.asarray([r]))[0]

    def e3(self, i, j, k, r1, r2, c):
        """one triplet, including both switches and the angular factor"""
        p = self.triple(i, j, k)
        q = self._as_pot(p)
        g = float(q.phi3(np.asarray([r1 + r2]))[0])
        if p["taper"] > 0:
            g *= float(q.switch(np.asarray([r1]), self.pair(i, j)["rcut3"], 0)[0])
            g *= float(q.switch(np.asarray([r2]), self.pair(i, k)["rcut3"], 0)[0])
        if p["lam2"] or p["lam4"]:
            c2 = c * c
            h = 1.0 + p["lam2"] * 0.5 * (3.0 * c2 - 1.0)
            if p["lam4"]:
                h += p["lam4"] * 0.125 * (35.0 * c2 * c2 - 30.0 * c2 + 3.0)
            g *= h
        return g


def energy(cry, pot, species):
    """energy per atom (eV); `species` gives an element index per basis atom

    Written as a plain double loop rather than vectorised.  It is the slow
    way and it is the right way here: this exists to be obviously correct so
    that when it disagrees with the fast implementation the fast one is the
    suspect.
    """
    sp = list(species)
    if len(sp) != len(cry.frac):
        raise ValueError(f"{len(sp)} tur, {len(cry.frac)} atom")

    e = 0.0
    for (i, j, R, d, r) in L.neighbours(cry, pot.rcut2):
        if r < pot.pair(sp[i], sp[j])["rcut2"]:
            e += 0.5 * pot.phi2(sp[i], sp[j], r)

    #  triplets: same enumeration as latdyn.triplets, but the cutoff that
    #  decides membership is the per-bond one, so the neighbour list is built
    #  at the widest and filtered per leg
    nb = {}
    for (i, j, R, d, r) in L.neighbours(cry, pot.rcut3):
        nb.setdefault(i, []).append((j, d, r))
    for i, lst in nb.items():
        for a in range(len(lst)):
            ja, da, ra = lst[a]
            if ra >= pot.pair(sp[i], sp[ja])["rcut3"]:
                continue
            for b in range(a + 1, len(lst)):
                jb, db, rb = lst[b]
                if rb >= pot.pair(sp[i], sp[jb])["rcut3"]:
                    continue
                c = float(np.dot(da, db) / (ra * rb))
                e += pot.e3(sp[i], sp[ja], sp[jb], ra, rb, c)
    return e / len(sp)


def selftest():
    """a one-element alloy file must reproduce the single-species energy"""
    import json
    import refdata
    lib = json.load(open(os.path.join(HERE, "library.json")))
    tmp = os.path.join(HERE, "_selftest.ugur")
    bad = 0
    print(f"{'el':4s}{'single-species':>13s}{'alloy.py':>13s}{'diff':>11s}")
    print("-" * 41)
    for el in ("Cu", "Fe", "Ti", "Mg"):
        rec = lib[el]
        vals = [rec[k] for k in FIELDS[:9]] + [rec.get("taper") or -1.0, 0.0, 0.0]
        open(tmp, "w").write(f"{el} {el} {el} " +
                             " ".join(f"{v:.17g}" for v in vals) + "\n")
        pot = AlloyPotential.from_file(tmp, [el])
        e = refdata.ELEMENTS[el]
        cry = L.Crystal(e["struct"], e["a0"], e.get("c_over_a"),
                        mass=refdata.MASSES[el])
        ref = L.energy(cry, L.Potential.from_record(rec))
        got = energy(cry, pot, [0] * len(cry.frac))
        d = got - ref
        bad += abs(d) > 1e-10 * max(abs(ref), 1.0)
        print(f"{el:4s}{ref:13.8f}{got:13.8f}{d:11.2e}")
    os.remove(tmp)
    print()
    print("alloy.py matches latdyn on a single element" if not bad
          else f"{bad} element uyusmadi")
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if selftest() else 0)
