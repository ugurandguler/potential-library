#!/usr/bin/env python3
"""
Measured phonon frequencies at high-symmetry points, in THz.

Nothing thermal or dynamical enters the fit - the objective sees only elastic
constants, which are the q -> 0 limit - so these are a genuine out-of-sample
test, and a harder one than the Materials Project comparison because they are
measurements rather than another calculation.

Values are stored as the SORTED triple of the three branch frequencies at each
point, not as named branches.  The named labels in the literature are not
consistent about which transverse branch is "T1": for Nb the N-point set is
N_L 5.66, N_T1 3.93, N_T2 5.07, so T1 < T2, while for Mo it is 8.14, 5.73, 4.56,
so T2 < T1.  Comparing sorted sets sidesteps a labelling convention that carries
no physics.

Degeneracies are written out in full: at fcc X and L the two transverse branches
are degenerate, and at bcc H all three branches are.

The q-points are the ones already in build_library.SC_POINTS:
    fcc  X = (1/2, 0, 1/2)     L = (1/2, 1/2, 1/2)
    bcc  H = (1/2, -1/2, 1/2)  N = (0, 0, 1/2)
"""

#  ---------------------------------------------------------------------------
#  Sources.  There IS a single authoritative compilation for the phonon
#  dispersion of the elements - Landolt-Boernstein New Series III/13a - and
#  every value below comes from it, so the set is internally consistent: one
#  compiler, one set of conventions, no mixing of measurement temperatures
#  between elements.
#
#  It is a printed volume, so the values here are taken from the portion of it
#  tabulated by Savrasov, which covers X_L, X_T, L_L, L_T for the fcc metals
#  and H, N_L, N_T1, N_T2 for the bcc ones.  V appears there as a calculation
#  with no measurement quoted, so V is absent below.  That tabulation is why
#  the coverage is seven elements rather than all twenty-eight; extending it
#  means going to III/13a itself or to the original scattering papers, and each
#  addition should carry its own reference in this dict rather than inheriting
#  this one.
LB_DIRECT = "Landolt-B&ouml;rnstein III/13a (Springer 1981)"

LB = ("Landolt-B&ouml;rnstein III/13a (Springer 1981), inelastic neutron "
      "scattering; values as tabulated in S. Y. Savrasov and D. Y. Savrasov, "
      "Phys. Rev. B <b>54</b>, 16487 (1996), Table I")

PHONON_EXP = {
    #  fcc: X and L, transverse branches doubly degenerate
    "Al": {"X": [5.78, 5.78, 9.69], "L": [4.19, 4.19, 9.69], "ref": LB},
    #  The 100 K measurement, and the temperature belongs in the label: lead
    #  has a Debye temperature near 105 K and the library's lattice constant
    #  is the 293 K one, so part of the disagreement on this element is that
    #  mismatch rather than the potential.  Both values are the zeta = 1 and
    #  zeta = 0.867 rows of the source's own tables.
    "Pb": {"X": [0.89, 0.89, 1.86], "L": [0.89, 0.89, 2.18],
           "ref": LB + ", at 100 K (B. N. Brockhouse, T. Arase, G. Caglioti, "
                       "K. R. Rao and A. D. B. Woods, Phys. Rev. <b>128</b>, "
                       "1099 (1962), Tables I and II)"},
    "Cu": {"X": [5.13, 5.13, 7.25], "L": [3.42, 3.42, 7.30], "ref": LB},
    #  This is the 120 K column of Miller and Brockhouse, not a room-
    #  temperature figure: their table gives 120, 296, 673 and 853 K, and
    #  4.64/6.72 at X and 3.34/7.02 at L match the 120 K column exactly.  The
    #  generic reference above does not say so, and a phonon frequency without
    #  its temperature cannot be compared with anything - palladium's [00z]L
    #  at X moves 6.72 -> 6.70 -> 6.47 THz from 120 to 673 K.
    "Pd": {"X": [4.64, 4.64, 6.72], "L": [3.34, 3.34, 7.02],
           "ref": LB + ", at 120 K (A. P. Miller and B. N. Brockhouse, "
                       "Can. J. Phys. <b>49</b>, 704 (1971), Table 3)"},
    #  Read off the papers themselves rather than through a compilation, so
    #  each carries its own page and temperature.  Both are the zeta = 1 point
    #  of [00z] and the zeta = 0.5 point of [zzz] in refdata_phonon_curves,
    #  which is where the whole branch lives.
    "Au": {"X": [2.75, 2.75, 4.61], "L": [1.86, 1.86, 4.70],
           "ref": "J. W. Lynn, H. G. Smith and R. M. Nicklow, Phys. Rev. B "
                  "<b>8</b>, 3493 (1973), Table I, 296 K"},
    "Pt": {"X": [3.84, 3.84, 5.80], "L": [2.90, 2.90, 5.85],
           "ref": "D. H. Dutton, B. N. Brockhouse and A. P. Miller, Can. J. "
                  "Phys. <b>50</b>, 2915 (1972), Table 1, 90 K"},
    #  Zone-boundary values read off the same tables as the full branches in
    #  refdata_phonon_curves: X is zeta = 1 on [00z], L is zeta = 0.5 on
    #  [zzz].  Each carries its own paper and temperature rather than
    #  inheriting a compilation's.
    "Ca": {"X": [3.63, 3.63, 4.52], "L": [2.36, 2.36, 4.61],
           "ref": "C. Stassis et al., Phys. Rev. B <b>27</b>, 3303 (1983), "
                  "Table I, 295 K"},
    "Yb": {"X": [1.85, 1.85, 2.40], "L": [1.16, 1.16, 2.30],
           "ref": "C. Stassis et al., Phys. Rev. B <b>26</b>, 4106 (1982), "
                  "Table I, 295 K"},
    "Ni": {"X": [6.27, 6.27, 8.55], "L": [4.24, 4.24, 8.88],
           "ref": "R. J. Birgeneau, J. Cordes, G. Dolling and A. D. B. Woods, "
                  "Phys. Rev. <b>136</b>, A1359 (1964), Table I, 296 K"},
    #  bcc: H triply degenerate, N has three distinct branches
    "Nb": {"H": [6.49, 6.49, 6.49], "N": [3.93, 5.07, 5.66], "ref": LB},
    "Ta": {"H": [5.03, 5.03, 5.03], "N": [2.63, 4.35, 4.35], "ref": LB},
    "Mo": {"H": [5.52, 5.52, 5.52], "N": [4.56, 5.73, 8.14], "ref": LB},
    #  H is the zeta = 1 end of [z00], where the three branches meet; the
    #  table stops at 0.95 on L and reaches 1.0 on T, and they agree there
    #  (2.13 against 2.15), which is the degeneracy showing up in the data
    #  rather than being assumed.  N is zeta = 0.5 on [zz0].
    "Ba": {"H": [2.15, 2.15, 2.15], "N": [0.75, 1.53, 2.30],
           "ref": "J. Mizuki, Y. Chen, K.-M. Ho and C. Stassis, Phys. Rev. B "
                  "<b>32</b>, 666 (1985), Table I, 295 K"},
    #  Read directly off III/13a rather than through the tabulation above, so
    #  these carry their own page reference.  The H values are the safest
    #  numbers in the volume: the [00zeta] L and T columns are listed
    #  separately and must meet at zeta = 1, and they do, which checks the
    #  column assignment without relying on the OCR of the headers.
    "Fe": {"H": [8.56, 8.56, 8.56], "N": [4.53, 6.45, 9.26],
           "ref": LB_DIRECT + " p. 62 (Minkiewicz et al. 1967, 295 K)"},
    "Na": {"H": [3.58, 3.58, 3.58],
           "ref": LB_DIRECT + " p. 107 (Woods et al., [00&zeta;] L at "
                              "&zeta; = 1)"},
    #  Read off Table 2 and cross-checked against Fig. 3 on the facing page,
    #  which plots all four measured temperatures: at 12 K the three branches
    #  meet at H, and at N they separate into 0.34, 0.96 and 1.50 THz, which is
    #  what the figure shows.  The check matters because the table interleaves
    #  four temperature columns and the OCR does not keep them apart.
    "Rb": {"H": [1.385, 1.385, 1.385], "N": [0.34, 0.96, 1.50],
           "ref": LB_DIRECT + " p. 123-125 (Copley and Brockhouse 1973, 12 K)"},

    #  ---- hcp ----------------------------------------------------------
    #  Two points, both complete without any symmetry assignment.
    #
    #  A = (0, 0, 1/2).  The little group there forces TA and TO to meet and LA
    #  and LO to meet, so the six branches collapse to two levels with
    #  multiplicity four and two.  Zirconium's table prints the degeneracy
    #  outright - TA and TO are both 1.81 THz at zeta = 0.5 - which is the
    #  volume checking itself.
    #
    #  Gamma.  The three optic modes of hcp are E2g, doubly degenerate, and B1g.
    #  Where the volume gives them under two labels they agree across
    #  directions: zirconium's Delta5 and Sigma4 both read 2.56 at zeta = 0, and
    #  Sigma1 reads 2.55.  Out of sample in the strongest sense - the fit sees
    #  the acoustic slopes at Gamma and nothing else.
    "Ti": {"A": [3.05, 3.05, 3.05, 3.05, 5.73, 5.73],
           "ref": LB_DIRECT + " Table 2 Ti, [79St2], 295 K "
                              "([00&zeta;] TA and LA at &zeta; = 0.5)"},
    "Zr": {"A": [1.81, 1.81, 1.81, 1.81, 4.15, 4.15],
           "G": [2.56, 2.56, 4.16],
           "ref": LB_DIRECT + " Table 2 Zr, [78St1], 295 K"},
    "Y":  {"G": [2.68, 2.68, 4.64],
           "ref": LB_DIRECT + " Table 2 Y, [70Si2], 295 K"},
    "Mg": {"G": [3.70, 3.70, 7.25],
           "ref": LB_DIRECT + " Table 2 Mg, [66Sq1, 68Py1] "
                              "(&Gamma;<sub>5</sub><sup>+</sup> and "
                              "&Gamma;<sub>3</sub><sup>+</sup>)"},
    #  80 K rather than room temperature; the volume gives the shift to 300 K
    #  for the E2g mode as 0.08 THz, so the temperature costs under a per cent
    "Be": {"G": [13.73, 13.73, 20.28],
           "ref": LB_DIRECT + " Table 2 Be, [76St1], 80 K"},
}

#  Elements looked for in III/13a and deliberately left out:
#
#  Ba  the volume states outright, p. 24: "The phonon dispersion curves have
#      not been measured in barium."  That is why no database carries it - the
#      measurement does not exist, so its absence here is not a gap to fill.
#  Co  measured, but the volume covers both alpha-Co (hcp) and beta-Co (fcc)
#      and the three data sets sit at 77-300 K and 472-721 K across the two
#      phases.  We fit hcp, so picking numbers here needs the phase of each
#      data set settled first; guessing would put a fcc measurement against an
#      hcp calculation.
NOT_MEASURED = {"Ba": "no measurement exists (Landolt-Boernstein III/13a p. 24)"}

#  Which points are usable per structure.  For hcp this is A and the zone
#  centre, not M: at M only four of the six branches are resolved in the volume
#  (magnesium gives M1+, M2-, M3-, M4+ and nothing for the other two), and
#  comparing an incomplete set against six computed branches would need an
#  irreducible-representation assignment we do not do.  A and Gamma are complete
#  and need none - see the comments on the hcp rows.
POINTS = {"fcc": ("X", "L"), "bcc": ("H", "N"), "hcp": ("A", "G")}

#  At Gamma the three acoustic branches are zero by construction, so a Gamma
#  entry holds the three OPTIC frequencies and the comparison takes the three
#  highest computed branches.  Nothing else needs special handling.
OPTIC_ONLY = ("G",)


def for_element(el, struct):
    """{point: sorted [THz]} for the points that belong to this structure"""
    rec = PHONON_EXP.get(el)
    if not rec:
        return {}
    return {p: sorted(rec[p]) for p in POINTS.get(struct, ()) if p in rec}


if __name__ == "__main__":
    for el, rec in sorted(PHONON_EXP.items()):
        pts = {k: v for k, v in rec.items() if k != "ref"}
        print(f"{el:3s} " + "  ".join(
            f"{k}={[round(x, 2) for x in sorted(v)]}" for k, v in pts.items()))
    print(f"\n{len(PHONON_EXP)} elements; source: {LB}")
