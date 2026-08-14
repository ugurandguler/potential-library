#!/usr/bin/env python3
"""
Experimental reference data used as FIT TARGETS for the Akgun-Ugur pair potential.

>>> PROVENANCE WARNING <<<
These are standard textbook / handbook room-temperature values, entered by hand:
  a0, Ecoh, B      : C. Kittel, "Introduction to Solid State Physics" (tables 1, 3, 4)
  Cij              : G. Simmons & H. Wang, "Single Crystal Elastic Constants and
                     Calculated Aggregate Properties", 2nd ed. (MIT Press, 1971)
                     and the CRC Handbook.
VERIFY EVERY ROW against the cited sources before using the resulting library in
published work.  The fitter is exact; its output is only as good as this table.
Editing a row here and re-running fit.py is all that is needed.

Units:  a0, (c/a) in Angstrom / dimensionless ; Ecoh in eV/atom (POSITIVE magnitude)
        B and Cij in GPa
Cij keys: cubic -> C11, C12, C44 ; hcp -> C11, C12, C13, C33, C44
"""

# structure, a0, c_over_a, Ecoh, B, Cij
#  Elastic constants for the CUBIC elements come from Landolt-Boernstein New
#  Series III/29a (Every and McCurdy, Springer 1992), Table 3, read off the
#  volume itself with lb29_extract.py and checked row by row against rendered
#  images of the page.  That table gives a weighted mean over the published
#  measurements, so it is more defensible than any single paper and, unlike a
#  mixture of compilations, it is internally consistent.
#
#  Several elements are listed twice with determinations that disagree.  The
#  one kept is whichever satisfies (C11 + 2 C12)/3 = measured bulk modulus:
#      Ba  84B8 (B = 9.53) over 85M4, whose C12 = -0.38 gives B = 2.45 against
#          a measured 9.4-10.3
#      Ir  66M1 (B = 355) over 65P1 (B = 373); measured is 355
#      Ca  86H6 (B = 18.3) over 83S9 (B = 21.4)
#      Sr  84B8 (B = 12.0) over 85M2 (B = 8.8)
#      Cu  the room-temperature row over the 4 K one
#  Adopting the volume also removed a long-standing inconsistency in Cr, whose
#  hand-entered bulk modulus of 190.1 GPa disagreed with its own Cij by 18 %;
#  the III/29a constants give 160.7, which is the accepted value.
#
#  The hexagonal elements come from Table 11 of the same volume, whose column
#  order is C11 C33 C44 C12 C13 - not the cubic order.  Agreement with the
#  previous hand-entered values was already good (most within a few per cent,
#  and inside the quoted spread), so this is for consistency of source rather
#  than a correction.  One number deserves care: **Be C13 = 6 with an
#  uncertainty of 9**, i.e. it is not really determined.  It still enters the
#  RMS with full weight, which is worth remembering when reading beryllium's
#  fit quality.
#
#  Beryllium is therefore NOT taken from III/29a.  Its C13 = 6 there is an
#  outlier: resonant ultrasound spectroscopy twelve years later (A. Migliori,
#  H. Ledbetter, D. J. Thoma and T. W. Darling, J. Appl. Phys. 95, 2436 (2004))
#  gives 14, as does the older hand-entered value, and Luo et al., J. Appl.
#  Phys. 111, 053503 (2012) compute 19.1.  The same bulk-modulus test that
#  settled barium settles this too - against a measured B of 116.8 GPa the
#  Voigt combination [2(C11+C12) + 4 C13 + C33]/9 gives 117.1 for Migliori,
#  114.5 for the old values, and 111.7 for III/29a.  So Be uses Migliori 2004.

#  Cohesive energies are Brewer, "The Cohesive Energies of the Elements",
#  LBL-3720 Rev. (1977), Table I, converted from kcal/gram-atom at
#  0.0433641 eV per kcal.  This is the source Kittel's table is drawn from, so
#  the hand-entered values already agreed with it to 0.005 eV everywhere except
#  tantalum, which was 8.10 against Brewer's 8.0657 and has been corrected.
#  The volume lists a row per crystal structure, and the row for the phase we
#  fit is the one taken.  Six rows defeated the OCR - aluminium comes back as
#  "A1" with a digit, vanadium in lower case, and four have their value on a
#  neighbouring text line - so those were read off page images; see
#  brewer_extract.py, which marks them.

#  ---------------------------------------------------------------------------
#  Ten elements added 2026-08-03, screened before being let in: every one fits
#  and every one is dynamically stable on the 8^3 U 9^3 union mesh.
#      a0, c/a   Kittel, "Introduction to Solid State Physics" 8th ed., Table 3
#                of chapter 1 (ICSD).  Three figures for most, which is coarser
#                than the rest of this table; Kittel's Table 4 gives the
#                nearest-neighbour distance independently and agrees with every
#                row, so they are right to the precision printed and no better.
#                Rubidium and caesium are its 5 K values, as Li, Na and K
#                already are.
#      B         Kittel Table 3 of chapter 3, in 10^12 dyn/cm^2 = 100 GPa - the
#                same table the rest of the B column comes from.  Substituted
#                below by the value the elastic constants imply, as for every
#                other element; for Sc, Y and Lu the two differ by 13-28 %,
#                because Gschneidner's 1964 rare-earth moduli are older than
#                the III/29a constants.
#      Ecoh      Brewer LBL-3720 Rev. (1977), via brewer_extract.py
#      Cij       III/29a Table 3 (cubic) and Table 11 (hexagonal).  Rubidium is
#                listed twice there and the row kept is 67G4, whose B = 2.90
#                matches the measured 3.1 to 6 %; the other gives 2.61, off by
#                16 %.  Caesium is also listed twice and the kept row is the
#                one that matches.
#  Technetium and osmium were looked at and left out; see candidates.py.

ELEMENTS = {
    # ---------------- fcc ----------------
    "Al": dict(struct="fcc", a0=4.0495, Ecoh=3.3867, B=72.2,
               Cij=dict(C11=108.0, C12=62.0, C44=28.3)),
    "Ni": dict(struct="fcc", a0=3.5240, Ecoh=4.4405, B=186.0,
               Cij=dict(C11=247.0, C12=153.0, C44=122.0)),
    "Cu": dict(struct="fcc", a0=3.6150, Ecoh=3.4865, B=137.0,
               Cij=dict(C11=169.0, C12=122.0, C44=75.3)),
    "Pd": dict(struct="fcc", a0=3.8907, Ecoh=3.8941, B=180.8,
               Cij=dict(C11=221.0, C12=171.0, C44=70.8)),
    "Ag": dict(struct="fcc", a0=4.0857, Ecoh=2.9488, B=100.7,
               Cij=dict(C11=122.0, C12=92.0, C44=45.5)),
    "Pt": dict(struct="fcc", a0=3.9239, Ecoh=5.8411, B=278.3,
               Cij=dict(C11=347.0, C12=251.0, C44=76.5)),
    "Au": dict(struct="fcc", a0=4.0782, Ecoh=3.8143, B=173.2,
               Cij=dict(C11=191.0, C12=162.0, C44=42.2)),
    "Pb": dict(struct="fcc", a0=4.9502, Ecoh=2.0286, B=43.0,
               Cij=dict(C11=48.8, C12=41.4, C44=14.8)),
    "Rh": dict(struct="fcc", a0=3.8034, Ecoh=5.7457, B=270.4,
               Cij=dict(C11=413.0, C12=194.0, C44=184.0)),
    "Ir": dict(struct="fcc", a0=3.8390, Ecoh=6.9426, B=355.0,
               Cij=dict(C11=580.0, C12=242.0, C44=256.0)),
    "Ca": dict(struct="fcc", a0=5.5884, Ecoh=1.843, B=15.2,
               Cij=dict(C11=22.8, C12=16.0, C44=14.0)),
    "Sr": dict(struct="fcc", a0=6.0849, Ecoh=1.7216, B=11.6,
               Cij=dict(C11=15.3, C12=10.3, C44=9.9)),
    "Yb": dict(struct="fcc", a0=5.48, Ecoh=1.605, B=13.3,
               Cij=dict(C11=18.6, C12=10.4, C44=17.7)),

    # ---------------- bcc ----------------
    "Fe": dict(struct="bcc", a0=2.8665, Ecoh=4.28, B=168.0,
               Cij=dict(C11=230.0, C12=135.0, C44=117.0)),
    "Cr": dict(struct="bcc", a0=2.8839, Ecoh=4.0979, B=190.1,
               Cij=dict(C11=348.0, C12=67.0, C44=100.0)),
    "Mo": dict(struct="bcc", a0=3.1470, Ecoh=6.8168, B=272.5,
               Cij=dict(C11=465.0, C12=163.0, C44=109.0)),
    "W":  dict(struct="bcc", a0=3.1652, Ecoh=8.8983, B=323.2,
               Cij=dict(C11=523.0, C12=203.0, C44=160.0)),
    "V":  dict(struct="bcc", a0=3.0240, Ecoh=5.3078, B=161.9,
               Cij=dict(C11=230.0, C12=120.0, C44=43.1)),
    "Nb": dict(struct="bcc", a0=3.3008, Ecoh=7.567, B=170.2,
               Cij=dict(C11=245.0, C12=132.0, C44=28.4)),
    "Ta": dict(struct="bcc", a0=3.3058, Ecoh=8.0657, B=200.0,
               Cij=dict(C11=264.0, C12=158.0, C44=82.6)),
    "Li": dict(struct="bcc", a0=3.4910, Ecoh=1.6353, B=11.6,
               Cij=dict(C11=13.4, C12=11.3, C44=9.6)),
    "Na": dict(struct="bcc", a0=4.2906, Ecoh=1.1127, B=6.8,
               Cij=dict(C11=7.59, C12=6.33, C44=4.3)),
    "K":  dict(struct="bcc", a0=5.3280, Ecoh=0.9341, B=3.7,
               Cij=dict(C11=3.71, C12=3.15, C44=1.88)),
    "Ba": dict(struct="bcc", a0=5.0280, Ecoh=1.895, B=9.4,
               Cij=dict(C11=12.6, C12=8.0, C44=9.5)),
    "Rb": dict(struct="bcc", a0=5.585, Ecoh=0.852, B=3.1,
               Cij=dict(C11=3.25, C12=2.73, C44=1.98)),
    "Cs": dict(struct="bcc", a0=6.045, Ecoh=0.804, B=2.0,
               Cij=dict(C11=2.47, C12=2.06, C44=1.48)),

    # ---------------- hcp ----------------
    "Mg": dict(struct="hcp", a0=3.2094, c_over_a=1.6236, Ecoh=1.5047, B=35.4,
               Cij=dict(C11=59.3, C12=25.7, C13=21.4, C33=61.5, C44=16.4)),
    "Ti": dict(struct="hcp", a0=2.9506, c_over_a=1.5873, Ecoh=4.8481, B=105.1,
               Cij=dict(C11=160.0, C12=90.0, C13=66.0, C33=181.0, C44=46.5)),
    "Zr": dict(struct="hcp", a0=3.2320, c_over_a=1.5931, Ecoh=6.2531, B=83.3,
               Cij=dict(C11=144.0, C12=74.0, C13=67.0, C33=166.0, C44=33.4)),
    "Co": dict(struct="hcp", a0=2.5071, c_over_a=1.6228, Ecoh=4.3928, B=191.4,
               Cij=dict(C11=295.0, C12=159.0, C13=111.0, C33=335.0, C44=71.0)),
    "Be": dict(struct="hcp", a0=2.2858, c_over_a=1.5677, Ecoh=3.3174, B=116.8,
               Cij=dict(C11=293.6, C12=26.8, C13=14.0, C33=356.7, C44=162.2)),
    "Zn": dict(struct="hcp", a0=2.6649, c_over_a=1.8563, Ecoh=1.346, B=59.8,
               Cij=dict(C11=165.0, C12=31.1, C13=50.0, C33=61.8, C44=39.6)),
    "Cd": dict(struct="hcp", a0=2.9793, c_over_a=1.8859, Ecoh=1.1591, B=46.7,
               Cij=dict(C11=114.1, C12=41.0, C13=40.3, C33=49.9, C44=19.0)),
    "Sc": dict(struct="hcp", a0=3.31, c_over_a=1.5921, Ecoh=3.898, B=43.5,
               Cij=dict(C11=99.3, C12=39.7, C13=29.4, C33=107.0, C44=27.7)),
    "Y":  dict(struct="hcp", a0=3.65, c_over_a=1.5699, Ecoh=4.371, B=36.6,
               Cij=dict(C11=77.9, C12=29.2, C13=20.0, C33=76.9, C44=24.3)),
    "Lu": dict(struct="hcp", a0=3.50, c_over_a=1.5857, Ecoh=4.432, B=41.1,
               Cij=dict(C11=86.2, C12=32.0, C13=28.0, C33=80.9, C44=26.8)),
    "Hf": dict(struct="hcp", a0=3.19, c_over_a=1.5831, Ecoh=6.435, B=109.0,
               Cij=dict(C11=181.0, C12=77.0, C13=66.0, C33=197.0, C44=55.7)),
    "Re": dict(struct="hcp", a0=2.76, c_over_a=1.6159, Ecoh=8.031, B=372.0,
               Cij=dict(C11=616.0, C12=273.0, C13=206.0, C33=683.0, C44=161.0)),
    "Ru": dict(struct="hcp", a0=2.71, c_over_a=1.5793, Ecoh=6.739, B=320.8,
               Cij=dict(C11=563.0, C12=188.0, C13=168.0, C33=624.0, C44=181.0)),
    "Tl": dict(struct="hcp", a0=3.46, c_over_a=1.5954, Ecoh=1.882, B=35.9,
               Cij=dict(C11=40.8, C12=35.4, C13=29.0, C33=52.8, C44=7.3)),
}

#  Standard molar entropy S° and heat capacity Cp° at 298.15 K, J/(mol.K).
#  Source: CRC Handbook of Chemistry and Physics, standard thermodynamic tables.
#  >>> Same provenance warning as above: hand-entered, verify before publishing. <<<
#
#  These are NOT fit targets.  They are here purely so the phonon thermodynamics
#  computed from the potential can be scored against experiment - a genuine test,
#  since nothing thermal enters the fit.
#
#  NB: the calculation gives Cv, the table gives Cp.  For metals near 300 K
#  Cp - Cv = T V alpha^2 B is roughly 1-2 J/(mol.K), so a calculated Cv a little
#  below the tabulated Cp is the expected behaviour, not an error.
THERMO_298 = {
    "Al": (28.30, 24.20), "Ni": (29.87, 26.07), "Cu": (33.15, 24.44),
    "Pd": (37.57, 25.98), "Ag": (42.55, 25.35), "Pt": (41.63, 25.86),
    "Au": (47.49, 25.42), "Pb": (64.81, 26.44), "Rh": (31.51, 24.98),
    "Ir": (35.48, 25.10), "Ca": (41.59, 25.93), "Sr": (55.00, 26.40),
    "Fe": (27.28, 25.10), "Cr": (23.77, 23.35), "Mo": (28.66, 24.06),
    "W":  (32.64, 24.27), "V":  (28.94, 24.89), "Nb": (36.40, 24.60),
    "Ta": (41.51, 25.36), "Li": (29.12, 24.86), "Na": (51.30, 28.23),
    "K":  (64.68, 29.60), "Ba": (62.80, 28.07), "Mg": (32.67, 24.87),
    "Ti": (30.72, 25.06), "Zr": (39.00, 25.36), "Co": (30.04, 24.81),
    "Be": (9.50, 16.44),  "Zn": (41.63, 25.39), "Cd": (51.80, 26.02),
}

#  ---------------------------------------------------------------------------
#  Make the bulk modulus consistent with the elastic constants.
#
#  For a cubic crystal B = (C11 + 2 C12)/3 is an identity, not a separate
#  measurement; for hexagonal the Voigt average is
#  B = [2(C11 + C12) + 4 C13 + C33]/9.  The B values above were taken from one
#  compilation and the Cij from another, and for five elements they contradict
#  each other by more than 10 %:
#
#      Cr  190.1 entered vs 152.3 implied   (+25 %)   - literature commonly 160
#      Zn   59.8          vs  75.1          (-20 %)
#      Cd   46.7          vs  57.6          (-19 %)
#      Ca   15.2          vs  18.2          (-17 %)
#      Zr   83.3          vs  95.4          (-13 %)
#
#  That contradiction was being fed straight into the fit, which holds B as a
#  hard constraint while scoring the Cij - so for those elements it was being
#  pulled toward two incompatible targets at once.  Deriving B from the Cij
#  removes the conflict and keeps the target set self-consistent.  The originally
#  entered numbers are preserved as B_literature for comparison.
#  Uncertainty on the cohesive energy, eV, from the same table.  It is not used
#  by the fit yet and that is the point: Ecoh is imposed EXACTLY at every trial
#  point, while Brewer knows it to between 0.2 % (Pt) and 4.6 % (Ba).  Barium
#  and strontium reach RMS 0.0 on their elastic constants, but they do it while
#  pinned to a cohesive energy uncertain at the per-cent level, so some of that
#  zero is false precision.  An uncertainty-weighted objective would use this.
ECOH_UNC = {
    "Ag": 0.0087,
    "Al": 0.0434,
    "Au": 0.013,
    "Ba": 0.0867,
    "Be": 0.065,
    "Ca": 0.0173,
    "Cd": 0.0065,
    "Co": 0.026,
    "Cr": 0.0434,
    "Cu": 0.013,
    "Fe": 0.013,
    "Ir": 0.065,
    "Li": 0.0087,
    "Mg": 0.013,
    "Mo": 0.0217,
    "Na": 0.0043,
    "Nb": 0.1735,
    "Ni": 0.0217,
    "Pb": 0.013,
    "Pd": 0.0217,
    "Pt": 0.013,
    "Rh": 0.0434,
    "Sr": 0.0434,
    "Ti": 0.0217,
    "V": 0.0867,
    "W": 0.0434,
    "Zr": 0.0434,
    #  added 2026-08-03, Brewer's own spread; see brewer_ecoh.json
    "Rb": 0.002, "Cs": 0.002, "Yb": 0.013, "Sc": 0.043, "Y": 0.030,
    "Lu": 0.009, "Hf": 0.043, "Re": 0.065, "Ru": 0.043, "Tl": 0.013,
}


#  Spread on the elastic constants, GPa, from the s(n=N) rows of III/29a
#  Table 3 - the scatter over the measurements that went into each weighted
#  mean.  Only the cubic elements, and only the twelve whose adopted row has
#  more than one determination behind it; a single measurement gets no spread
#  quoted (Ir, Pt, Rh, K, Li, Ba, Ca, Ni, Sr).  Beryllium is absent because it
#  does not use III/29a at all.  Hexagonal spreads are in the volume too but
#  Table 11 has a different column order and has not been extracted yet.
#
#  Like ECOH_UNC this is not used by the fit.  It is here so that an
#  uncertainty-weighted objective becomes a change of a few lines rather than a
#  data-collection exercise: at the moment a 5 GPa error in Cr's C12, which is
#  one standard deviation, is scored the same as a 5 GPa error in W's C11,
#  which is five.
CIJ_UNC = {
    "Ag": {"C11": 2.0, "C12": 3.0, "C44": 1.0},
    "Al": {"C11": 2.0, "C12": 2.0, "C44": 0.2},
    "Au": {"C11": 2.0, "C12": 3.0, "C44": 0.8},
    "Cr": {"C11": 4.0, "C12": 5.0, "C44": 0.5},
    "Cu": {"C11": 1.5, "C12": 1.8, "C44": 0.6},
    "Fe": {"C11": 5.0, "C12": 4.0, "C44": 1.0},
    "Na": {"C11": 0.15, "C12": 0.13, "C44": 0.09},
    "Nb": {"C11": 5.0, "C12": 5.0, "C44": 0.3},
    "Pb": {"C11": 1.0, "C12": 1.2, "C44": 0.3},
    "Ta": {"C11": 5.0, "C12": 5.0, "C44": 0.6},
    "V": {"C11": 5.0, "C12": 4.0, "C44": 0.4},
    "W": {"C11": 1.0, "C12": 1.0, "C44": 1.0},
}


def _bulk_from_cij(c, struct):
    if struct == "hcp":
        return (2*(c["C11"] + c["C12"]) + 4*c["C13"] + c["C33"]) / 9.0
    return (c["C11"] + 2*c["C12"]) / 3.0

for _el, _e in ELEMENTS.items():
    _e["B_literature"] = _e["B"]
    _e["B"] = _bulk_from_cij(_e["Cij"], _e["struct"])


# atomic masses (amu) - needed for the dynamical matrix
#  Melting points, K.  CRC Handbook of Chemistry and Physics, 97th ed. (2016),
#  section 4, "Physical Constants of Inorganic Compounds" / elemental table.
#
#  These are NOT fitting targets and never enter the objective.  They set the
#  temperature a parameter set has to remain a crystal up to, which is what the
#  compression-escape constraint compares its barrier against.  A potential
#  whose lattice is only metastable at half its own melting point is not usable
#  for dynamics however well it reproduces the elastic tensor.
MELTING = {
    "Li": 453.7, "Na": 371.0, "K": 336.7, "Rb": 312.5, "Cs": 301.6,
    "Be": 1560.0, "Mg": 923.0, "Ca": 1115.0, "Sr": 1050.0, "Ba": 1000.0,
    "Al": 933.5, "Ti": 1941.0, "V": 2183.0, "Cr": 2180.0, "Fe": 1811.0,
    "Co": 1768.0, "Ni": 1728.0, "Cu": 1357.8, "Zn": 692.7, "Y": 1799.0,
    "Zr": 2128.0, "Nb": 2750.0, "Mo": 2896.0, "Ru": 2607.0, "Rh": 2237.0,
    "Pd": 1828.1, "Ag": 1234.9, "Cd": 594.2, "Sc": 1814.0, "Hf": 2506.0,
    "Ta": 3290.0, "W": 3695.0, "Re": 3459.0, "Ir": 2719.0, "Pt": 2041.4,
    "Au": 1337.3, "Tl": 577.0, "Pb": 600.6, "Lu": 1925.0, "Yb": 1097.0,
}


MASSES = {
    "Al": 26.9815, "Ni": 58.6934, "Cu": 63.546, "Pd": 106.42, "Ag": 107.868,
    "Pt": 195.084, "Au": 196.967, "Pb": 207.2, "Rh": 102.906, "Ir": 192.217,
    "Ca": 40.078, "Sr": 87.62, "Fe": 55.845, "Cr": 51.996, "Mo": 95.95,
    "W": 183.84, "V": 50.9415, "Nb": 92.906, "Ta": 180.948, "Li": 6.94,
    "Na": 22.9898, "K": 39.098, "Ba": 137.327, "Mg": 24.305, "Ti": 47.867,
    "Zr": 91.224, "Co": 58.9332, "Be": 9.0122, "Zn": 65.38, "Cd": 112.414,
    #  added 2026-08-03, IUPAC standard atomic weights (CIAAW 2021)
    "Rb": 85.4678, "Cs": 132.90545, "Yb": 173.045, "Sc": 44.955908,
    "Y": 88.90584, "Lu": 174.9668, "Hf": 178.486, "Re": 186.207,
    "Ru": 101.07, "Tl": 204.38,
}
