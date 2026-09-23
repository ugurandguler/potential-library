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
#                Rubidium and caesium are its 5 K values.  So is nothing
#                else: see the correction below - lithium is at 78 K and
#                sodium and potassium are at room temperature.
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
#
#  ---------------------------------------------------------------------------
#  TEMPERATURE.  Read this before comparing anything here with a calculation.
#
#  The three anchors come from three different temperatures, and the fit that
#  uses them is STATIC - a zero-kelvin lattice sum.  Nothing here is wrong in
#  itself; what was wrong is that none of it was written down.
#
#      Ecoh      0 K          Brewer LBL-3720 Rev. (1977)
#      a0        ~293 K       Kittel Table 3, ICSD - EXCEPT Li at 78 K and
#                             Rb, Cs at 5 K.  Na and K are ~293 K like the
#                             rest; see the correction below.
#      Cij       ~300 K       III/29a section 1.1.9, in the volume's own words:
#                             "Unless otherwise stated, all elastic constants
#                             are given at room temperature, RT, (= 300 K)."
#
#  So a static calculation is being matched to room-temperature stiffnesses and
#  a room-temperature lattice, while its cohesive energy is a 0 K quantity.
#  The library is therefore fitted slightly too soft and slightly too expanded
#  against a true 0 K reference.  III/29a Table 28 gives the temperature
#  coefficients, and the size of that "slightly" is not uniform:
#
#      W    C44 x1.02 at 0 K        Cu   C11 x1.06  C44 x1.11
#      Ta   C11 x1.03  C44 x1.08    Al   C11 x1.10  C44 x1.15
#      Ba   C44 x1.46               Na   C11 x1.21  C44 x1.67
#
#  The refractory metals barely notice; sodium's C44 target is SIXTY PER CENT
#  of its zero-kelvin value.  C44 moves about twice as far as C11 everywhere,
#  so the anisotropy the fit is asked to reach is distorted as well as the
#  magnitude.
#
#  THREE ELEMENTS CARRY AN EXPLICIT TEMPERATURE in Table 3 and the row kept
#  here is not always the room-temperature one:
#
#      Li   two rows: RT (13.4/11.3/9.6) and 195 K (13.9/11.7/9.85).
#           The RT row is the one used.  Correct.
#      Rb   two rows: RT (2.96/2.44/1.60) and ~80 K (3.25/2.73/1.98).
#           The 80 K row is the one used.
#      Cs   two rows: 78 K (2.47/2.06/1.48) and 280 K (1.60/0.99/1.44).
#           The 78 K row is the one used - C11 is 54 % above the 280 K value
#           and C12 is 108 % above it.
#
#  Rb and Cs were NOT chosen for their temperature; they were chosen because
#  their bulk modulus matched Kittel's, and Kittel's alkali entries are
#  themselves low-temperature.  The reasoning was accidental and the outcome
#  is defensible: their a0 IS Kittel's 5 K value, so a low-temperature Cij is
#  more consistent with a static fit than the room-temperature row would be.
#  Rubidium and caesium are, by accident, the two most internally consistent
#  records in this table.
#
#  CORRECTED 2026-08-29, and the correction is not cosmetic.  The note above
#  used to say all five alkalis carried Kittel's 5 K lattice constant.  Three
#  of the five do not.  Kittel's Table 3 is a periodic-table graphic that
#  nothing extracts, so it was read off a rendering; it states its own
#  convention - "The data given are at room temperature for the most common
#  form, or at the stated temperature in deg K", crediting the Inorganic
#  Crystal Structure Database - and stamps Li 78 K, Na 5 K,
#  K 5 K, Rb 5 K, Cs 5 K.  Against the values actually in this file:
#
#      Li  3.491  = Kittel's 78 K value        -> 78 K, not 5 K
#      Rb  5.585  = Kittel's 5 K value         -> 5 K, as claimed
#      Cs  6.045  = Kittel's 5 K value         -> 5 K, as claimed
#      Na  4.2906 against Kittel's 4.225       -> NOT Kittel's; ~293 K
#      K   5.328  against Kittel's 5.225       -> NOT Kittel's; ~293 K
#
#  So SODIUM AND POTASSIUM ARE NOT THE INCONSISTENT RECORDS this header used
#  to name.  They are the opposite: lattice constant and elastic constants
#  both at room temperature, consistent with each other and both simply in
#  the wrong place, which is the ordinary case rather than the pathological
#  one.  Lithium is genuinely inconsistent, by 78 K rather than by 5 K.
#
#  The earlier claim came from Kittel's CHAPTER 6 Table 1, free electron
#  Fermi surface parameters, whose header does say "except for Na, K, Rb, Cs
#  at 5 K and Li at 78 K".  That is a different table about a different
#  quantity.  The note below already flagged that it proved nothing about a0;
#  it understated the problem.
#
#  ZERO-POINT ENERGY, which is a separate systematic and larger than the
#  temperature one for the light elements.
#
#  A measured cohesive energy is the work to take the crystal apart FROM ITS
#  ZERO-POINT STATE.  A classical potential has no zero-point motion, so its
#  static well should be deeper than the measured value by exactly that
#  energy: the target for D is E_coh + E_ZPE, not E_coh.  Every element in
#  this table is therefore fitted a little too shallow.
#
#  Using the Debye model, E_ZPE = (9/8) k_B theta_D:
#
#      Be   4.2 %     Mg 2.6 %    Li 2.0 %    Cr 1.5 %    Na 1.4 %
#      Al   1.2 %     Fe 1.1 %    Ni 1.0 %    W  0.44 %   Au 0.42 %
#      Pt   0.40 %    Nb 0.35 %   Ta 0.29 %
#
#  CAVEAT, NOW CLOSED (2026-08-29).  These per cents were quoted rather than
#  subtracted because the Debye temperatures behind them were entered by hand
#  and were in no file.  They are now in THETA_D below, from Stewart, Rev.
#  Sci. Instrum. 54, 1 (1983), Table I, and they reproduce every figure here
#  to within 0.13 points.  The hand values were right; they are now sourced.
#  What has NOT changed is that nothing is subtracted: see the closing note.
#
#  Looked for, 2026-08-12, and not found in a usable form:
#
#    CRC (the 2000 printing in the tree) has no single elemental table.  Its
#    superconductive-elements table gives theta_D for Al Ir Mo Nb Pb Ru Ta Ti
#    Tl V W Zn Zr, and its rare-earth table adds Sc Y Yb Lu.  Where those
#    overlap the hand values they agree to about 5 % (worst: Pb 9 %, Ta 7 %),
#    which is far inside what a per-cent-level estimate needs - so the sizing
#    above stands.  But the table MISSES Be, Li, Mg, Na and K, which are
#    exactly the five elements where the correction is largest.
#
#    One column trap in it: rhenium's theta_D cell reads 4.5, which is the
#    neighbouring gamma column.  Rhenium's is about 430 K.  A blind extraction
#    would have written 4.5 into the library.
#
#    Kittel has the complete table (ch. 5, Table 1) but in this PDF it is a
#    periodic-table graphic; nothing extracts.  Same for the lattice-constant
#    table of ch. 1, which is why the note below could not be checked either.
#
#  OPEN, and flagged rather than fixed: the a0 note above says Li, Na, K, Rb
#  and Cs are all Kittel's 5 K values.  Kittel's ch. 6 Table 1 states its own
#  convention as "room temperature, except for Na, K, Rb, Cs at 5 K and Li at
#  78 K".  That is a different table from the one a0 comes from, so it proves
#  nothing about a0 - but it is a reason to check rather than to assume that
#  lithium sits with the other four.
#
#  The ordering is what matters and it is robust: the correction is largest
#  exactly where the library already struggles - beryllium, magnesium, lithium
#  - and negligible for the refractory metals, which are its cleanest records.
#
#  Nothing is changed here.  Changing a target means refitting the element,
#  and refitting on this account means refitting the library; that is a
#  decision about scope, not a correction, and it is recorded here in this
#  header rather than made silently in the data below.
#
#  ---------------------------------------------------------------------------
#  AND THE DECISION WAS TAKEN, 2026-08-29: THE TARGETS STAY WHERE THEY ARE.
#  Not for want of doing the work.  All three anchors were carried to 0 K and
#  the fit was tried against them, and it does not help.
#
#  WHAT WAS MEASURED.  The elastic anchor moved with the temperature
#  coefficients of Landolt-Boernstein III/29a - Table 28 for the cubic
#  elements, Table 35 for the hexagonal ones - which cover 37 of the 40.
#  Cobalt, scandium and ytterbium have no coefficient anywhere in the volume.
#  The lattice constant moved by the Debye-shaped expansion integral under the
#  CRC's alpha(25 C) with THETA_D below, from each element's OWN reference
#  temperature (78 K for Li, 5 K for Rb and Cs, ~293 K for the rest).  The
#  cohesive energy moved by (9/8) k_B theta_D.  How far:
#
#      volume            Cs +0.0 %    ...   Na +5.1 %,  K +6.5 %
#      bulk modulus      Ba -2.0 %    ...   Na +19.3 %, Cr +27.5 %
#      cohesive energy   Ta +0.3 %    ...   Mg +2.6 %,  Be +4.3 %
#
#  Every corrected tensor is still positive definite; caesium is closest, with
#  a smallest eigenvalue of 0.29 GPa against a C' target that falls from 0.2 to
#  0.1 GPa.  So the corrected set is reachable in principle.
#
#  WHAT IT COST TO CHECK, AND WHAT CAME BACK.  Four elements - copper as the
#  control, lead, silver, and sodium as a labelled probe - were fitted twice,
#  to the old targets and the new, at the same search budget and seed, with and
#  without the cutoff taper, and scored against the measured dispersion, which
#  never enters the fit.  Two improved, six got worse, and THE CONTROL GOT
#  WORSE IN BOTH CUTOFFS.  Sodium is the instructive one: its targets are not
#  reachable at all (the published tapered sodium sits at 31 % residual after
#  400 restarts) and the corrected targets halve that to 16 % - better
#  reachability, worse dispersion.  Fitting those targets more closely is
#  fitting the wrong thing more closely.
#
#  SO THE CORRECTION IS REAL AND IT IS NOT A LEVER.  The shift clears the
#  published uncertainty on the target it corrects for 37 of 40 elements, so
#  "it is inside the noise" is not available as a reason to ignore it.  It is
#  simply not what stands between this library and a better one.  That is a
#  result about the objective, not a failure of the measurement, and it is the
#  question the whole exercise was opened to answer.
#
#  READ THE FIGURES ABOVE AS ERROR BARS, NOT AS CORRECTIONS.  They say how far
#  each target sits from the 0 K quantity the static fit actually computes,
#  which is worth knowing when a residual is being interpreted - sodium's C44
#  target is sixty per cent of its zero-kelvin value and no fit can be judged
#  against it as if it were exact.  They are not a better place to aim.
#
#  Two cautions that came out of the same work, for anyone who does apply them:
#    * C' is a DIFFERENCE, and six elements have a coefficient for one of c11
#      and c12 and not the other (Nb, Pt, Re, Ti, Tl, Zr).  Correcting one and
#      not the other moves the anisotropy for a reason that is not physics -
#      thallium's C' falls 41 % on that account alone.  Both or neither.
#    * chromium's three coefficients are all bracketed in the volume, its
#      Tc12 = (-15) is the largest magnitude in the cubic table, it alone moves
#      chromium's bulk modulus by 28 %, and chromium's Neel transition at 311 K
#      sits just above the range they were measured over.
#  ---------------------------------------------------------------------------
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
#
#  Source, in the publisher's recommended form and the SAME printing as the
#  expansion coefficients of expansion.py:
#
#    David R. Lide, ed., CRC Handbook of Chemistry and Physics, Internet
#    Version 2005, <http://www.hbcpnetbase.com>, CRC Press, Boca Raton, FL,
#    section 5, "Standard Thermodynamic Properties of Chemical Substances".
#
#  VERIFIED 2026-09-23, which closes the hand-entered warning this block used
#  to carry.  crc_thermo.py reads all 38 element rows out of the PDF in the
#  repository; 36 agree with the values below to the one decimal the volume
#  prints, and two did not, by more than rounding:
#
#      Al  Cp  was 24.20, volume prints 24.4
#      Sr  Cp  was 26.40, volume prints 26.8
#
#  Both are corrected below.  The extra digit on the other entries is kept -
#  it agrees with the volume wherever the volume can resolve it - but nothing
#  in the library needs it, and a reader checking against the printed table
#  will find one decimal.
#
#  These are NOT fit targets.  They are here purely so the phonon thermodynamics
#  computed from the potential can be scored against experiment - a genuine test,
#  since nothing thermal enters the fit.
#
#  NB: the calculation gives Cv, the table gives Cp.  For metals near 300 K
#  Cp - Cv = T V alpha^2 B is roughly 1-2 J/(mol.K), so a calculated Cv a little
#  below the tabulated Cp is the expected behaviour, not an error.
THERMO_298 = {
    "Al": (28.30, 24.40), "Ni": (29.87, 26.07), "Cu": (33.15, 24.44),
    "Pd": (37.57, 25.98), "Ag": (42.55, 25.35), "Pt": (41.63, 25.86),
    "Au": (47.49, 25.42), "Pb": (64.81, 26.44), "Rh": (31.51, 24.98),
    "Ir": (35.48, 25.10), "Ca": (41.59, 25.93), "Sr": (55.00, 26.80),
    "Fe": (27.28, 25.10), "Cr": (23.77, 23.35), "Mo": (28.66, 24.06),
    "W":  (32.64, 24.27), "V":  (28.94, 24.89), "Nb": (36.40, 24.60),
    "Ta": (41.51, 25.36), "Li": (29.12, 24.86), "Na": (51.30, 28.23),
    "K":  (64.68, 29.60), "Ba": (62.80, 28.07), "Mg": (32.67, 24.87),
    "Ti": (30.72, 25.06), "Zr": (39.00, 25.36), "Co": (30.04, 24.81),
    "Be": (9.50, 16.44),  "Zn": (41.63, 25.39), "Cd": (51.80, 26.02),
    #  Added 2026-08-03 with the ten new elements, read off the volume itself
    #  by crc_thermo.py (section 5, "Standard Thermodynamic Properties of
    #  Chemical Substances") rather than typed in - the pages are recorded in
    #  crc_thermo.json.  One decimal, which is what that table prints.
    #
    #  That read also covers the twenty eight rows above, and they agree, but
    #  not all to the last digit: this edition prints Al Cp 24.4 against the
    #  24.20 here, Ba S 62.5 against 62.80, Sr Cp 26.8 against 26.40, and the
    #  rest within 0.1.  The older rows are left as they are.  They match the
    #  commonly quoted values and CRC's own editions differ between printings;
    #  overwriting a settled column on the strength of one copy would be a
    #  change, not a correction.  It is recorded here so the difference is
    #  known rather than discovered again.
    "Rb": (76.8, 31.1), "Cs": (85.2, 32.2), "Yb": (59.9, 26.7),
    "Sc": (34.6, 25.5), "Y":  (44.4, 26.5), "Lu": (51.0, 26.9),
    "Hf": (43.6, 25.7), "Re": (36.9, 25.5), "Ru": (28.5, 24.1),
    "Tl": (64.2, 26.3),
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
    #  Cadmium and zinc, from the s(n=6) and s(n=5) spread rows of III/29a
    #  Table 11 - the volume's own mean-of-determinations and its scatter, read
    #  positionally because those two pages scan badly.  Zinc's C12 is the
    #  outlier at 20 %, which Ledbetter's independent NBS compilation
    #  (J. Phys. Chem. Ref. Data 6, 1181 (1977), 100 references) puts at 17 %;
    #  he also notes zinc's S12 is positive where every other hexagonal metal
    #  measured has it negative, so the scatter is not simply carelessness.
    #  Relevant because these are the two elements with no solution: over half
    #  the hexagonal target for zinc is known to worse than 10 %.
    "Cd": {"C11": 4.0, "C33": 1.9, "C44": 1.7, "C12": 3.7, "C13": 2.0},
    "Zn": {"C11": 7.0, "C33": 6.8, "C44": 1.0, "C12": 6.2, "C13": 6.0},
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



#  Debye temperatures at 0 K, in K.
#
#  Source: G. R. Stewart, "Measurement of low-temperature specific heat",
#  Rev. Sci. Instrum. 54, 1 (1983), TABLE I - specific heat parameters gamma
#  and theta_D for the elements, laid out as a periodic table.  Its own
#  footnote settles the one thing that matters for a zero-point energy: "The
#  Debye temperatures quoted are the low-temperature values, T << theta_D."
#  That is the T -> 0 quantity, stated by the source rather than assumed.
#  All forty were read off the paper's page at 4x.
#
#  THIS CLOSES THE CAVEAT IN THE HEADER ABOVE, which recorded that the Debye
#  temperatures behind its per-cent figures were entered by hand and were in
#  no file.  They are here now, and they reproduce those figures: across a
#  fifteen-fold range the largest disagreement is 0.13 points - Mg 2.60
#  against 2.60, Nb 0.35 against 0.35, Ta 0.29 against 0.29, Be 4.20 against
#  4.33.  The hand values were right.
#
#  Checked two ways.  The CRC's superconductive-elements table agrees on the
#  fifteen elements it covers, and its rare-earth table agrees on Y, Lu and
#  Yb - but it MISSES Be, Li, Mg, Na and K, which are the five where the
#  correction is largest, and it prints theta_D = 4.5 K for rhenium, which a
#  4x rendering confirms is the printed page and not the OCR (rhenium is
#  416 K, as Stewart has it).  A source can be wrong in the shape of data.
#  Independently, an elastic Debye temperature computed from the Cij above,
#  carried to 0 K and averaged over the Christoffel slowness surface, runs at
#  1.00 of these values with a 5 % spread.
#
#  STEWART STARS THE VALUES HE DERIVED FROM ELASTIC CONSTANTS rather than
#  from a specific heat.  Fe (477*) and Yb (118*) are the two here, so for
#  those two an elastic calculation is not an independent check of this
#  column - it is the same method twice.
#
#  These are NOT fit targets.  Nothing in the fit uses them.
THETA_D = {
    "Ag": 227.3, "Al": 433.0, "Au": 162.3, "Ba": 111.0, "Be": 1481.0,
    "Ca": 229.0, "Cd": 210.0, "Co": 460.0, "Cr": 606.0, "Cs": 40.5,
    "Cu": 347.0, "Fe": 477.0, "Hf": 252.0, "Ir": 420.0, "K": 91.1,
    "Li": 344.0, "Lu": 183.0, "Mg": 403.0, "Mo": 423.0, "Na": 156.5,
    "Nb": 276.0, "Ni": 477.0, "Pb": 105.0, "Pd": 271.0, "Pt": 237.0,
    "Rb": 56.5, "Re": 416.0, "Rh": 512.0, "Ru": 555.0, "Sc": 346.0,
    "Sr": 147.0, "Ta": 245.0, "Ti": 420.0, "Tl": 78.5, "V": 399.0,
    "W": 383.0, "Y": 248.0, "Yb": 118.0, "Zn": 329.0, "Zr": 290.0,
}
#  derived from elastic constants in Stewart's own table, not calorimetry
THETA_D_FROM_ELASTIC = {"Fe", "Yb"}

#  Melting points, K.  CRC Handbook of Chemistry and Physics, 97th ed.
#  (2016), section 4, "Physical Constants of Inorganic Compounds" /
#  elemental table.
#
#  THE EDITION IS PART OF THE CITATION.  Checked 2026-09-07 against the CRC
#  copy in this repository, which is the Internet Version 2005 - a DIFFERENT
#  printing - and the two disagree:
#      Cu  ours 1357.8 K   2005 ed. 1084.62 C = 1357.77 K    +0.03
#      Mo  ours 2896.0 K   2005 ed. 2622    C = 2895.15 K    +0.85
#      Ta  ours 3290.0 K   2005 ed. 3007    C = 3280.15 K    +9.85
#      W   ours 3695.0 K   2005 ed. 3414    C = 3687.15 K    +7.85
#  Copper agrees to 0.03 K and tantalum is ten kelvin out, so a reader who
#  checks against the PDF on disk will find a mismatch and conclude the table
#  is wrong.  It is not; it is the 97th.  The same edition drift is already
#  recorded for CRC's expansion coefficients in expansion.py, where both
#  printings are kept rather than reconciled.
#
#  These are NOT fitting targets and never enter the objective.  They set the
#  temperature a parameter set has to remain a crystal up to, which is what
#  the compression-escape constraint compares its barrier against.  A
#  potential whose lattice is only metastable at half its own melting point is
#  not usable for dynamics however well it reproduces the elastic tensor.
#
#  The finite-temperature melting work quotes tantalum's 3290 K as the
#  benchmark its force-matched record misses by about thirty per cent, so this
#  is the citation that number needs.  Ten kelvin of edition drift is nothing
#  against a ~1000 K discrepancy.
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


#  Atomic masses, amu - needed for the dynamical matrix, where they set the
#  frequency scale as sqrt(1/M).
#
#  IUPAC standard atomic weights (CIAAW 2021).  The attribution used to sit
#  INSIDE the table, on the row of the ten added 2026-08-03, which is why an
#  audit reading the lines above the table found nothing; it covers all forty.
#  This is a later revision than the 2005 CRC printing on disk, which prints
#  Mo 95.94(2) against our 95.95 - 1 part in 10^4, so 5 parts in 10^5 on a
#  frequency, far under any error bar here.
MASSES = {
    "Al": 26.9815, "Ni": 58.6934, "Cu": 63.546, "Pd": 106.42, "Ag": 107.868,
    "Pt": 195.084, "Au": 196.967, "Pb": 207.2, "Rh": 102.906, "Ir": 192.217,
    "Ca": 40.078, "Sr": 87.62, "Fe": 55.845, "Cr": 51.996, "Mo": 95.95,
    "W": 183.84, "V": 50.9415, "Nb": 92.906, "Ta": 180.948, "Li": 6.94,
    "Na": 22.9898, "K": 39.098, "Ba": 137.327, "Mg": 24.305, "Ti": 47.867,
    "Zr": 91.224, "Co": 58.9332, "Be": 9.0122, "Zn": 65.38, "Cd": 112.414,
    #  added 2026-08-03
    "Rb": 85.4678, "Cs": 132.90545, "Yb": 173.045, "Sc": 44.955908,
    "Y": 88.90584, "Lu": 174.9668, "Hf": 178.486, "Re": 186.207,
    "Ru": 101.07, "Tl": 204.38,
}
