#!/usr/bin/env python3
"""
Measured phonon dispersion CURVES - whole branches, not just the zone-boundary
point.

refdata_phonon.py carries three frequencies at X and at L per element, which is
what fits in a table.  This carries the branches themselves, so the measured
points can be drawn along the dispersion panel wherever the path passes them.
It is the only comparison on the page against a measurement rather than another
calculation, and the fit sees none of it: the objective is elastic constants,
the q -> 0 limit alone.

Each point is [zeta, nu, sigma], nu and sigma in THz, exactly as printed.
Nothing is converted, averaged or symmetrised here.

TEMPERATURE IS PART OF THE DATUM, not a footnote.  Frequencies soften
measurably with it - palladium's zone-boundary [00z]L moves 6.72 -> 6.70 ->
6.47 THz between 120 and 673 K in the source below - so a curve without its
temperature cannot be compared with anything.  The three sets here sit at three
different temperatures because that is what was measured, and each says so.
Gold is the only one whose full dispersion exists at room temperature at all;
platinum's does not, which is why every density-functional paper compares
against 90 K.

Branch labels stay in each source's own notation.  Gold writes the Sigma
direction [zz0]; platinum and palladium write it [0zz].  They are the same
line, and the mapping to the drawn path is made where the drawing happens
rather than by rewriting what the papers printed.
"""

#  fcc high-symmetry lines.  zeta is the reduced wave vector q a / 2 pi.
#      [00z]  Delta,  Gamma -> X   (X at zeta = 1)
#      [zzz]  Lambda, Gamma -> L   (L at zeta = 0.5)
#      [zz0] / [0zz]  Sigma, Gamma -> K -> X   (K at zeta = 0.75)
#      [1z0]  Z,      X -> W       (W at zeta = 0.5)

AU_REF = ("J. W. Lynn, H. G. Smith and R. M. Nicklow, "
          "Phys. Rev. B <b>8</b>, 3493 (1973), Table I; "
          "coherent inelastic neutron scattering at 296 K")
PT_REF = ("D. H. Dutton, B. N. Brockhouse and A. P. Miller, "
          "Can. J. Phys. <b>50</b>, 2915 (1972), Tables 1 and 2; "
          "inelastic neutron scattering at 90 K")
PD_REF = ("A. P. Miller and B. N. Brockhouse, "
          "Can. J. Phys. <b>49</b>, 704 (1971), Table 3, 296 K column; "
          "inelastic neutron scattering")
CA_REF = ("C. Stassis, J. Zarestky, D. K. Misemer, H. L. Skriver, "
          "B. N. Harmon and R. M. Nicklow, Phys. Rev. B <b>27</b>, 3303 "
          "(1983), Table I; inelastic neutron scattering at 295 K")
YB_REF = ("C. Stassis, C.-K. Loong, B. N. Harmon, S. H. Liu and R. M. "
          "Nicklow, Phys. Rev. B <b>26</b>, 4106 (1982), Table I; "
          "inelastic neutron scattering at 295 K")
NI_REF = ("R. J. Birgeneau, J. Cordes, G. Dolling and A. D. B. Woods, "
          "Phys. Rev. <b>136</b>, A1359 (1964), Table I; "
          "inelastic neutron scattering at 296 K")
CU_REF = ("E. C. Svensson, B. N. Brockhouse and J. M. Rowe, "
          "Phys. Rev. <b>155</b>, 619 (1967), Table I; "
          "inelastic neutron scattering at 296 K")
BA_REF = ("J. Mizuki, Y. Chen, K.-M. Ho and C. Stassis, "
          "Phys. Rev. B <b>32</b>, 666 (1985), Table I; "
          "inelastic neutron scattering on single crystals at 295 K")
CS_REF = ("J. Mizuki and C. Stassis, "
          "Phys. Rev. B <b>34</b>, 5890 (1986), Table I, 280 K.  Caesium "
          "melts at 302 K, so this was measured at 0.93 of the melting "
          "temperature and the paper attributes the softening it sees to "
          "cubic anharmonicity - a harmonic 0 K calculation is not expected "
          "to land on it")

Y_REF = ("S. K. Sinha, T. O. Brun, L. D. Muhlestein and J. Sakurai, "
         "Phys. Rev. B <b>1</b>, 2430 (1970), Table I, 295 K.  The paper is "
         "a scan with no text layer, so this was transcribed from the page "
         "images")

ZR_REF = ("R. E. Schmunk, H. F. Bezdek and L. Finegold, "
          "phys. stat. sol. <b>42</b>, 275 (1970), Table 1, room "
          "temperature.  Energies are printed in meV and zeta as "
          "q/q<sub>max</sub>, so the zone boundary is at 1 rather than at "
          "1/2; both are handled on the way in and the table below holds "
          "the printed numbers")

SC_REF = ("N. Wakabayashi, S. K. Sinha and F. H. Spedding, "
          "Phys. Rev. B <b>4</b>, 2398 (1971), Table II, 295 K.  The two "
          "basal directions are printed with different q units - 4 pi / "
          "(sqrt3 a) along Gamma-M and 4 pi / a along Gamma-K-M - and both "
          "put the zone boundary where the other tables here put it")

BE_REF = ("R. Stedman, Z. Amilius, R. Pauli and O. Sundin, "
          "J. Phys. F <b>6</b>, 157 (1976), Table 1, 80 K.  This is the "
          "measurement Landolt-B&ouml;rnstein quotes for beryllium as "
          "[76St1]; the value already carried for Gamma comes from that "
          "compilation and agrees with this table exactly")

HF_REF = ("C. Stassis, D. Arch, O. D. McMasters and B. N. Harmon, "
          "Phys. Rev. B <b>24</b>, 730 (1981), Table I, 295 K column.  The "
          "paper also measured 800 K and 1300 K")

TI_REF = ("C. Stassis, D. Arch and B. N. Harmon, "
          "Phys. Rev. B <b>19</b>, 181 (1979), Table I, 295 K column.  The "
          "paper also measured 773 K and 1054 K; only the room-temperature "
          "column is carried, to keep this set comparable")

MO_REF = ("A. D. B. Woods and S. H. Chen, Solid State Commun. <b>2</b>, 233 "
          "(1964), Table 1, 296 K.  The companion paper to the tungsten "
          "record above - same spectrometer at Chalk River, same method, "
          "same temperature - and like it, the dispersion is PLOTTED and "
          "only eleven modes are tabulated, so this is a handful of "
          "measurements rather than a curve.  Frequencies are in "
          "10<sup>12</sup> c/s, which is THz.  Two caveats travel with it. "
          "The paper reports a pronounced anomaly near H, which it argues is "
          "a Kohn anomaly tied to the Fermi surface; a short-ranged "
          "potential cannot produce one, so the error at H is not a fitting "
          "failure and should not be read as one.  And between P and H the "
          "TRANSVERSE branch runs ABOVE the longitudinal - 7.08 against 6.28 "
          "- so a comparison that assigns branches by sorted position rather "
          "than by the paper's own labels gets molybdenum backwards there.")


W_REF = ("S. H. Chen and B. N. Brockhouse, Solid State Commun. <b>2</b>, 73 "
         "(1964), Table 1, room temperature.  The paper plots its dispersion "
         "and tabulates only eight specific points, so this is a handful of "
         "measurements and not a curve.  Frequencies are in "
         "10<sup>12</sup> c/s, which is THz, and the coordinates are in "
         "units of 2 pi / a with a = 3.165 &Aring;")

MG_REF = ("P. K. Iyengar, G. Venkataraman, P. R. Vijayaraghavan and A. P. Roy, "
          "<i>Lattice dynamics of magnesium</i>, in <i>Inelastic Scattering "
          "of Neutrons</i> (IAEA, Vienna, 1965), Vol. I, p. 153, Table IV, "
          "room temperature; frequencies in 10<sup>12</sup> c/s against "
          "q/q<sub>max</sub>.  The three K values and the 5.32 THz at M come "
          "from G. L. Squires, Proc. Phys. Soc. <b>88</b>, 919 (1966), "
          "Table 3, which this table does not reach")

NB_REF = ("Y. Nakagawa and A. D. B. Woods, in <i>Lattice Dynamics</i> "
          "(Pergamon 1965), p. 39, Table 1, 296 K.  Frequencies are printed "
          "in 10<sup>12</sup> c/s, which is THz.  The paper also measured "
          "[zz1] and [1/2 1/2 z]; neither lies on the standard path and "
          "neither is carried here")

AG_REF = ("W. Drexel, Z. Physik <b>255</b>, 281 (1972), Table 1, room "
          "temperature.  q is the MAGNITUDE in &Aring;<sup>-1</sup> and "
          "omega is in 10<sup>13</sup> rad/s, so each direction has its own "
          "zone-boundary value - 1.538 along [001], 1.631 along [011] and "
          "1.332 along [111] - and all three are handled by the ranges "
          "below.  The values are interpolated between neighbouring "
          "detectors and allowed up to 3 degrees off the symmetry direction; "
          "the paper compares itself with a triple-axis measurement and puts "
          "the worst disagreement at about 4 %, in the [011] T1 branch")

FE_REF = ("V. J. Minkiewicz, G. Shirane and R. Nathans, "
          "Phys. Rev. <b>162</b>, 528 (1967), Tables I and II, 295 K.  "
          "Energies are in meV and q is the CARTESIAN COMPONENT in "
          "&Aring;<sup>-1</sup> rather than a reduced coordinate, so the "
          "zone boundary sits at 2 pi / a = 2.192 for [00z] and [zzz] and "
          "at half that for [zz0]; both are handled by the ranges below and "
          "the printed numbers are kept.  The paper's pi and Lambda blocks "
          "lie on a line the standard path does not follow and are not "
          "carried")

AL_REF = ("R. Stedman and G. Nilsson, Phys. Rev. <b>145</b>, 492 (1966), "
          "<b>80 K</b>, measured by neutron spectrometry; their numbers are "
          "tabulated in Landolt-B&ouml;rnstein III/13a (Schober and "
          "Dederichs 1981), Table 2 Al, where they are marked unpublished")

LI_REF = ("M. M. Beg and M. Nielsen, Phys. Rev. B <b>14</b>, 4266 (1976), "
          "Table I, <b>293 K</b>.  Coherent inelastic neutron scattering on "
          "<sup>7</sup>Li; energies in meV")

RB_REF = ("J. R. D. Copley and B. N. Brockhouse, Can. J. Phys. <b>51</b>, 657 "
          "(1973), Table 1, <b>120 K</b> column.  Frequencies in THz")

NA_REF = ("A. D. B. Woods, B. N. Brockhouse, R. H. March, A. T. Stewart "
          "and R. Bowers, Phys. Rev. <b>128</b>, 1112 (1962), Table III, "
          "<b>90 K</b>.  Frequencies in 10<sup>12</sup> cps.  Sodium "
          "undergoes a martensitic transformation below about 40 K, which is "
          "why the measurement is not colder")

K_REF = ("R. A. Cowley, A. D. B. Woods and G. Dolling, "
         "Phys. Rev. <b>150</b>, 487 (1966), Table I, <b>9 K</b>.  "
         "Frequencies in 10<sup>12</sup> cps, zeta the reduced coordinate.  "
         "Nine kelvin is 0.027 of potassium's melting point, so this is very "
         "nearly the harmonic limit and is a CLEANER comparison for a 0 K "
         "lattice-dynamics calculation than a room-temperature measurement "
         "would be.  The same group's Phys. Rev. <b>180</b>, 755 (1969) "
         "measures the widths and shifts up to 299 K but tabulates no "
         "frequencies")

LU_REF = ("J. Pleschiutschnig, O. Blaschko and W. Reichardt, "
          "Phys. Rev. B <b>41</b>, 975 (1990), Table I, 295 K.  Each of its "
          "three blocks carries its own q unit - 2 pi / a along Gamma-K-M, "
          "2 pi / c along Gamma-A and 2 pi / (a sqrt3) along Gamma-M - and "
          "all three are mapped by seg_override rather than by rescaling the "
          "printed numbers")

TL_REF = ("T. G. Worlton and R. E. Schmunk, Phys. Rev. B <b>3</b>, 4115 "
          "(1971), <b>77 K</b>, tabulated in Landolt-B&ouml;rnstein III/13a "
          "(Schober and Dederichs 1981), Table 2 Tl.  Abscissa "
          "q/q<sub>max</sub>; only &Delta; and &Sigma; were measured")

TA_REF = ("A. D. B. Woods, Phys. Rev. <b>136</b>, A781 (1964), Table I, "
          "296 K.  Frequencies are printed in 10<sup>12</sup> cps, which is "
          "THz")

PB_REF = ("B. N. Brockhouse, T. Arase, G. Caglioti, K. R. Rao and "
          "A. D. B. Woods, Phys. Rev. <b>128</b>, 1099 (1962), Tables I "
          "and II; inelastic neutron scattering at 100 K")

SQ2 = 1.41421356
RADS13_THZ = 10.0 / 6.283185307   # 10^13 rad/s -> THz
MEV_THZ = 1.0 / 4.135667696      # h in meV / THz
SQ3 = 1.73205081

PHONON_CURVE = {

    #  ---- gold, 296 K, all nine branches -------------------------------
    "Au": {"T_K": 296, "ref": AU_REF, "branches": {
        "L[00z]": [[0.052, 0.50, 0.07], [0.092, 0.75, 0.06], [0.116, 1.00, 0.05],
                   [0.179, 1.50, 0.05], [0.206, 1.65, 0.04], [0.293, 2.25, 0.04],
                   [0.30, 2.32, 0.04], [0.40, 2.85, 0.06], [0.47, 3.30, 0.08],
                   [0.552, 3.70, 0.08], [0.60, 3.88, 0.06], [0.70, 4.17, 0.06],
                   [0.80, 4.40, 0.06], [0.90, 4.58, 0.06], [1.00, 4.61, 0.05]],
        "T[00z]": [[0.10, 0.41, 0.04], [0.15, 0.58, 0.03], [0.20, 0.75, 0.03],
                   [0.25, 0.91, 0.03], [0.30, 1.09, 0.03], [0.40, 1.46, 0.03],
                   [0.50, 1.84, 0.03], [0.60, 2.16, 0.03], [0.70, 2.41, 0.04],
                   [0.80, 2.61, 0.05], [0.90, 2.73, 0.05], [1.00, 2.75, 0.04]],
        "L[zz0]": [[0.058, 0.75, 0.08], [0.074, 1.00, 0.07], [0.109, 1.25, 0.08],
                   [0.12, 1.50, 0.10], [0.17, 2.00, 0.08], [0.201, 2.40, 0.06],
                   [0.28, 3.00, 0.05], [0.35, 3.50, 0.05], [0.45, 3.85, 0.10],
                   [0.50, 3.96, 0.06], [0.60, 3.84, 0.05], [0.66, 3.70, 0.05],
                   [0.70, 3.55, 0.06], [0.75, 3.34, 0.05], [0.80, 3.15, 0.05],
                   [0.90, 2.83, 0.05], [1.00, 2.75, 0.04]],
        "T1[zz0]": [[0.10, 0.31, 0.02], [0.15, 0.48, 0.02], [0.20, 0.63, 0.04],
                    [0.25, 0.83, 0.03], [0.30, 1.02, 0.03], [0.35, 1.23, 0.03],
                    [0.40, 1.44, 0.03], [0.45, 1.64, 0.05], [0.50, 1.79, 0.04],
                    [0.55, 1.96, 0.05], [0.60, 2.08, 0.05], [0.65, 2.18, 0.05],
                    [0.70, 2.34, 0.05], [0.75, 2.43, 0.05], [0.80, 2.55, 0.08],
                    [0.90, 2.70, 0.06], [0.95, 2.74, 0.06], [1.00, 2.75, 0.04]],
        "T2[zz0]": [[0.10, 0.56, 0.03], [0.20, 1.10, 0.03], [0.30, 1.63, 0.04],
                    [0.40, 2.15, 0.05], [0.50, 2.71, 0.05], [0.60, 3.25, 0.05],
                    [0.70, 3.77, 0.07], [0.80, 4.26, 0.07], [0.90, 4.53, 0.07],
                    [1.00, 4.61, 0.05]],
        "L[zzz]": [[0.035, 0.60, 0.08], [0.06, 1.00, 0.08], [0.12, 1.80, 0.09],
                   [0.176, 2.40, 0.10], [0.23, 3.00, 0.15], [0.26, 3.40, 0.15],
                   [0.30, 3.75, 0.06], [0.35, 4.14, 0.06], [0.40, 4.45, 0.05],
                   [0.45, 4.64, 0.05], [0.50, 4.70, 0.04]],
        "T[zzz]": [[0.10, 0.51, 0.03], [0.15, 0.76, 0.03], [0.20, 1.00, 0.02],
                   [0.25, 1.23, 0.02], [0.30, 1.44, 0.02], [0.35, 1.63, 0.03],
                   [0.40, 1.74, 0.03], [0.45, 1.85, 0.04], [0.50, 1.86, 0.04]],
        "pi[1z0]": [[0.00, 4.61, 0.05], [0.10, 4.57, 0.05], [0.20, 4.46, 0.05],
                    [0.30, 4.23, 0.05], [0.40, 3.92, 0.05], [0.50, 3.63, 0.05],
                    [0.60, 3.30, 0.05], [0.70, 3.06, 0.06], [0.80, 2.90, 0.07],
                    [0.90, 2.80, 0.07], [1.00, 2.75, 0.04]],
        "Lambda[1z0]": [[0.00, 2.75, 0.04], [0.10, 2.74, 0.05], [0.20, 2.73, 0.05],
                        [0.30, 2.70, 0.04], [0.40, 2.63, 0.03], [0.50, 2.63, 0.03],
                        [0.60, 2.64, 0.03], [0.70, 2.69, 0.04], [0.80, 2.73, 0.05],
                        [0.90, 2.74, 0.05], [1.00, 2.75, 0.04]],
    }},

    #  ---- platinum, 90 K.  Table 1 for eight branches, Table 2 for T1 ---
    #  Measured at 90 K and, until 2026-08-29, scored at refdata's ~293 K
    #  a0.  a(90 K) from the Debye-shaped expansion integral under the CRC's
    #  alpha(25 C) with Stewart's theta_D; the integral runs ~12 % high for
    #  the softest metals and platinum is not one of them.
    "Pt": {"T_K": 90, "a_meas": 3.9177, "a_meas_from": "expansion integral", "ref": PT_REF, "branches": {
        "T[00z]": [[0.10, 0.50, 0.03], [0.15, 0.75, 0.02], [0.20, 1.00, 0.02],
                   [0.25, 1.23, 0.02], [0.30, 1.48, 0.03], [0.35, 1.71, 0.02],
                   [0.40, 1.98, 0.02], [0.45, 2.22, 0.02], [0.50, 2.45, 0.05],
                   [0.55, 2.70, 0.03], [0.60, 2.93, 0.03], [0.65, 3.09, 0.03],
                   [0.70, 3.30, 0.03], [0.80, 3.57, 0.04], [0.90, 3.76, 0.04],
                   [1.00, 3.84, 0.05]],
        "L[00z]": [[0.20, 1.99, 0.05], [0.30, 2.79, 0.04], [0.40, 3.52, 0.05],
                   [0.50, 4.19, 0.05], [0.60, 4.77, 0.03], [0.70, 5.18, 0.03],
                   [0.80, 5.56, 0.04], [0.90, 5.73, 0.05], [1.00, 5.80, 0.08]],
        "T[zzz]": [[0.10, 0.76, 0.02], [0.15, 1.06, 0.03], [0.20, 1.40, 0.03],
                   [0.25, 1.72, 0.03], [0.30, 2.08, 0.03], [0.35, 2.38, 0.04],
                   [0.40, 2.63, 0.05], [0.50, 2.90, 0.03]],
        "L[zzz]": [[0.10, 1.76, 0.05], [0.15, 2.64, 0.05], [0.20, 3.44, 0.06],
                   [0.30, 4.77, 0.04], [0.40, 5.60, 0.04], [0.50, 5.85, 0.05]],
        "T2[0zz]": [[0.10, 0.70, 0.02], [0.15, 1.05, 0.02], [0.20, 1.38, 0.02],
                    [0.25, 1.73, 0.03], [0.30, 2.06, 0.03], [0.35, 2.36, 0.04],
                    [0.40, 2.71, 0.03], [0.45, 3.10, 0.03], [0.50, 3.47, 0.03],
                    [0.55, 3.85, 0.03], [0.60, 4.20, 0.03], [0.70, 4.87, 0.03],
                    [0.80, 5.36, 0.05], [0.90, 5.66, 0.07], [1.00, 5.80, 0.08]],
        "L[0zz]": [[0.10, 1.49, 0.04], [0.15, 2.14, 0.03], [0.20, 2.82, 0.04],
                   [0.25, 3.40, 0.05], [0.30, 3.94, 0.05], [0.325, 4.20, 0.08],
                   [0.40, 4.77, 0.07], [0.50, 4.98, 0.07], [0.60, 4.95, 0.07],
                   [0.75, 4.30, 0.07], [0.90, 3.89, 0.06], [1.00, 3.84, 0.05]],
        "Lambda[1z0]": [[0.10, 3.72, 0.10], [0.30, 3.50, 0.10], [0.40, 3.36, 0.05],
                        [0.50, 3.25, 0.08]],
        "pi[1z0]": [[0.10, 5.70, 0.09], [0.20, 5.58, 0.07], [0.30, 5.28, 0.07],
                    [0.40, 4.95, 0.07], [0.50, 4.65, 0.07], [0.60, 4.44, 0.08],
                    [0.75, 4.03, 0.05]],
        #  Table 2, 90 K column - the branch with the kink at zeta ~ 0.33
        #  that the paper is largely about.
        "T1[0zz]": [[0.100, 0.580, 0.02], [0.150, 0.835, 0.02], [0.200, 1.050, 0.02],
                    [0.250, 1.255, 0.02], [0.275, 1.345, 0.02], [0.300, 1.430, 0.02],
                    [0.325, 1.510, 0.03], [0.350, 1.540, 0.03], [0.375, 1.650, 0.03],
                    [0.400, 1.770, 0.03], [0.425, 1.895, 0.02], [0.450, 2.020, 0.02],
                    [0.500, 2.255, 0.02], [0.550, 2.510, 0.02], [0.600, 2.810, 0.02],
                    [0.650, 3.080, 0.02], [0.700, 3.27, 0.03], [0.800, 3.57, 0.03],
                    [0.900, 3.76, 0.04], [1.000, 3.84, 0.05]],
    }},

    #  ---- palladium, the 296 K column of a four-temperature table ------
    "Pd": {"T_K": 296, "ref": PD_REF, "branches": {
        "T[00z]": [[0.1, 0.66, 0.04], [0.2, 1.29, 0.03], [0.25, 1.63, 0.02],
                   [0.3, 1.93, 0.03], [0.4, 2.555, 0.03], [0.5, 3.16, 0.03],
                   [0.6, 3.65, 0.03], [0.7, 4.10, 0.05], [0.75, 4.28, 0.06],
                   [0.8, 4.36, 0.04], [0.9, 4.54, 0.06], [1.0, 4.56, 0.06]],
        "L[00z]": [[0.1, 1.18, 0.06], [0.2, 2.21, 0.03], [0.25, 2.72, 0.06],
                   [0.3, 3.15, 0.05], [0.4, 4.01, 0.06], [0.5, 4.78, 0.06],
                   [0.6, 5.48, 0.07], [0.7, 6.06, 0.07], [0.8, 6.39, 0.09],
                   [0.9, 6.66, 0.09], [1.0, 6.70, 0.09]],
        "T[zzz]": [[0.1, 0.84, 0.03], [0.15, 1.22, 0.03], [0.2, 1.59, 0.04],
                   [0.25, 2.01, 0.04], [0.3, 2.42, 0.04], [0.35, 2.76, 0.04],
                   [0.4, 3.02, 0.05], [0.45, 3.18, 0.07], [0.5, 3.21, 0.08]],
        "L[zzz]": [[0.1, 2.13, 0.07], [0.15, 3.12, 0.09], [0.2, 4.00, 0.11],
                   [0.25, 4.84, 0.10], [0.3, 5.55, 0.13], [0.35, 6.16, 0.12],
                   [0.4, 6.65, 0.13], [0.45, 6.85, 0.10], [0.5, 6.86, 0.13]],
        "T2[0zz]": [[0.1, 0.88, 0.02], [0.2, 1.84, 0.02], [0.3, 2.70, 0.03],
                    [0.4, 3.56, 0.03], [0.5, 4.43, 0.04], [0.6, 5.15, 0.05],
                    [0.7, 5.80, 0.06], [0.75, 6.05, 0.07], [0.8, 6.32, 0.07],
                    [0.9, 6.64, 0.10]],
        "T1[0zz]": [[0.1, 0.575, 0.02], [0.2, 1.055, 0.02], [0.3, 1.475, 0.03],
                    [0.4, 2.135, 0.03], [0.5, 2.835, 0.03], [0.6, 3.420, 0.03],
                    [0.7, 3.90, 0.04], [0.75, 4.09, 0.04], [0.8, 4.29, 0.04],
                    [0.875, 4.44, 0.07], [0.9, 4.47, 0.04]],
    }},

    #  ---- calcium, 295 K, nine branches --------------------------------
    #  Written [0z1] here for the X -> W line that gold and platinum write
    #  [1z0]; the papers disagree about which index carries the 1 and the
    #  line is the same one.
    "Ca": {"T_K": 295, "ref": CA_REF, "branches": {
        "L[00z]": [[0.2, 1.50, 0.04], [0.3, 2.28, 0.06], [0.4, 3.17, 0.08],
                   [0.5, 3.66, 0.10], [0.6, 4.13, 0.08], [0.7, 4.35, 0.08],
                   [0.8, 4.40, 0.06], [0.9, 4.59, 0.08], [1.0, 4.52, 0.08]],
        "T[00z]": [[0.15, 0.83, 0.04], [0.2, 1.07, 0.10], [0.3, 1.69, 0.10],
                   [0.4, 2.28, 0.04], [0.5, 2.54, 0.05], [0.6, 2.92, 0.06],
                   [0.7, 3.29, 0.04], [0.8, 3.42, 0.04], [0.9, 3.61, 0.03],
                   [1.0, 3.63, 0.06]],
        "L[zzz]": [[0.1, 1.60, 0.08], [0.15, 2.32, 0.06], [0.2, 2.90, 0.10],
                   [0.25, 3.34, 0.06], [0.3, 3.90, 0.14], [0.4, 4.51, 0.12],
                   [0.5, 4.61, 0.08]],
        "T[zzz]": [[0.2, 1.26, 0.04], [0.3, 1.81, 0.05], [0.4, 2.16, 0.04],
                   [0.5, 2.36, 0.06]],
        "L[zz0]": [[0.1, 1.29, 0.04], [0.15, 1.78, 0.04], [0.2, 2.44, 0.03],
                   [0.3, 3.22, 0.04], [0.4, 3.94, 0.04], [0.5, 4.14, 0.10],
                   [0.6, 4.19, 0.05], [0.7, 4.04, 0.04], [0.8, 3.84, 0.06],
                   [0.9, 3.68, 0.06]],
        "T1[zz0]": [[0.2, 0.80, 0.04], [0.3, 1.15, 0.04], [0.4, 1.64, 0.04],
                    [0.5, 2.01, 0.08], [0.6, 2.63, 0.03], [0.7, 3.02, 0.03],
                    [0.8, 3.14, 0.10], [0.9, 3.30, 0.08], [1.0, 3.63, 0.06]],
        "T2[zz0]": [[0.3, 2.42, 0.09], [0.4, 2.91, 0.03], [0.5, 3.50, 0.07],
                    [0.6, 4.14, 0.04], [0.7, 4.39, 0.04], [0.8, 4.61, 0.03],
                    [0.9, 4.70, 0.03], [1.0, 4.62, 0.04]],
        "Lambda[0z1]": [[0.1, 3.64, 0.06], [0.2, 3.50, 0.06], [0.3, 3.29, 0.06],
                        [0.4, 3.09, 0.06], [0.5, 2.89, 0.09]],
        "pi[0z1]": [[0.2, 4.62, 0.04], [0.4, 4.42, 0.06], [0.6, 4.16, 0.08],
                    [0.8, 3.79, 0.03]],
    }},

    #  ---- ytterbium, 295 K, nine branches ------------------------------
    #  The same group and the same instruments as calcium a year earlier,
    #  and the paper is explicit that the two look alike because both metals
    #  are divalent with similar bands.  It is also the source of the one
    #  fact that closes the fcc set: "the dispersion curves of fcc Sr have
    #  not yet been measured".
    "Yb": {"T_K": 295, "ref": YB_REF, "branches": {
        "L[00z]": [[0.15, 0.44, 0.04], [0.2, 0.60, 0.03], [0.3, 0.95, 0.03],
                   [0.4, 1.20, 0.04], [0.5, 1.55, 0.04], [0.6, 1.80, 0.04],
                   [0.7, 2.05, 0.04], [0.8, 2.24, 0.12], [0.9, 2.33, 0.05],
                   [1.0, 2.40, 0.06]],
        "T[00z]": [[0.15, 0.385, 0.02], [0.2, 0.535, 0.02], [0.3, 0.835, 0.02],
                   [0.4, 1.10, 0.04], [0.5, 1.35, 0.04], [0.6, 1.54, 0.04],
                   [0.7, 1.70, 0.06], [0.8, 1.75, 0.04], [0.9, 1.80, 0.05],
                   [1.0, 1.85, 0.05]],
        "L[zzz]": [[0.1, 0.65, 0.04], [0.15, 1.04, 0.05], [0.2, 1.40, 0.05],
                   [0.3, 1.95, 0.04], [0.4, 2.23, 0.06], [0.5, 2.30, 0.04]],
        "T[zzz]": [[0.1, 0.28, 0.03], [0.15, 0.45, 0.03], [0.2, 0.63, 0.03],
                   [0.25, 0.80, 0.04], [0.3, 0.95, 0.05], [0.35, 1.05, 0.04],
                   [0.4, 1.10, 0.04], [0.45, 1.14, 0.04], [0.5, 1.16, 0.04]],
        "L[0zz]": [[0.1, 0.50, 0.03], [0.15, 0.75, 0.03], [0.2, 1.05, 0.03],
                   [0.3, 1.57, 0.04], [0.4, 1.80, 0.04], [0.5, 2.00, 0.04],
                   [0.6, 2.05, 0.04], [0.7, 2.00, 0.05], [0.8, 1.95, 0.05],
                   [0.9, 1.88, 0.06], [1.0, 1.85, 0.04]],
        "T1[0zz]": [[0.15, 0.30, 0.03], [0.2, 0.425, 0.02], [0.3, 0.64, 0.03],
                    [0.4, 0.86, 0.03], [0.5, 1.10, 0.04], [0.6, 1.35, 0.03],
                    [0.7, 1.55, 0.04], [0.8, 1.70, 0.05], [0.9, 1.83, 0.05],
                    [1.0, 1.85, 0.07]],
        "T2[0zz]": [[0.1, 0.38, 0.02], [0.15, 0.58, 0.02], [0.2, 0.80, 0.03],
                    [0.3, 1.18, 0.04], [0.4, 1.56, 0.04], [0.5, 1.85, 0.04],
                    [0.6, 2.05, 0.04], [0.7, 2.20, 0.05], [0.8, 2.30, 0.04],
                    [0.9, 2.40, 0.06], [1.0, 2.40, 0.05]],
        "Lambda[0z1]": [[0.1, 1.86, 0.04], [0.2, 1.80, 0.05], [0.3, 1.70, 0.07],
                        [0.4, 1.63, 0.06], [0.5, 1.65, 0.05]],
        "pi[0z1]": [[0.1, 2.40, 0.06], [0.2, 2.38, 0.05], [0.4, 2.30, 0.05],
                    [0.6, 2.10, 0.03], [0.8, 1.90, 0.07]],
    }},

    #  ---- nickel, 296 K, nine branches ---------------------------------
    #  Ferromagnetic, and the paper says so: no one-magnon groups were seen,
    #  and the magnon and phonon branches cross only below zeta = 0.1.  So
    #  these are phonons, not a mixture, which matters for a potential with
    #  no spin in it at all.
    "Ni": {"T_K": 296, "ref": NI_REF, "branches": {
        "L[00z]": [[0.1, 1.71, 0.10], [0.2, 3.12, 0.09], [0.3, 4.42, 0.11],
                   [0.4, 5.58, 0.12], [0.5, 6.54, 0.13], [0.6, 7.34, 0.12],
                   [0.7, 7.94, 0.14], [0.8, 8.34, 0.13], [0.85, 8.50, 0.24],
                   [0.9, 8.56, 0.18], [0.95, 8.65, 0.20], [1.0, 8.55, 0.13]],
        "T[00z]": [[0.2, 2.03, 0.04], [0.3, 2.99, 0.06], [0.4, 3.83, 0.06],
                   [0.5, 4.49, 0.08], [0.6, 5.12, 0.09], [0.7, 5.67, 0.11],
                   [0.75, 5.83, 0.12], [0.8, 6.01, 0.12], [0.85, 6.07, 0.12],
                   [0.9, 6.23, 0.13], [0.95, 6.24, 0.14], [1.0, 6.27, 0.10]],
        "L[zzz]": [[0.1, 3.05, 0.05], [0.15, 4.39, 0.07], [0.2, 5.60, 0.09],
                   [0.25, 6.61, 0.10], [0.3, 7.44, 0.14], [0.35, 8.14, 0.16],
                   [0.4, 8.53, 0.17], [0.425, 8.61, 0.24], [0.45, 8.79, 0.18],
                   [0.475, 8.58, 0.25], [0.5, 8.88, 0.17]],
        "T[zzz]": [[0.1, 1.33, 0.04], [0.15, 1.89, 0.05], [0.2, 2.47, 0.05],
                   [0.25, 2.99, 0.05], [0.3, 3.37, 0.05], [0.35, 3.76, 0.05],
                   [0.375, 3.90, 0.07], [0.4, 4.02, 0.06], [0.425, 4.10, 0.08],
                   [0.45, 4.26, 0.06], [0.475, 4.24, 0.08], [0.5, 4.24, 0.06]],
        "L[0zz]": [[0.1, 2.34, 0.09], [0.2, 4.44, 0.11], [0.3, 6.08, 0.10],
                   [0.4, 7.25, 0.17], [0.5, 7.63, 0.20], [0.55, 7.69, 0.27],
                   [0.6, 7.68, 0.18], [0.65, 7.47, 0.19], [0.7, 7.39, 0.17],
                   [0.75, 7.30, 0.24], [0.8, 6.85, 0.17], [0.85, 6.74, 0.25],
                   [0.9, 6.54, 0.17], [0.95, 6.36, 0.15]],
        "T1[0zz]": [[0.2, 1.96, 0.05], [0.3, 2.81, 0.08], [0.4, 3.62, 0.09],
                    [0.5, 4.36, 0.08], [0.6, 4.98, 0.10], [0.7, 5.59, 0.12],
                    [0.8, 5.97, 0.13], [0.9, 6.26, 0.15]],
        "T2[0zz]": [[0.1, 1.28, 0.05], [0.2, 2.76, 0.10], [0.3, 4.14, 0.14],
                    [0.4, 5.50, 0.18], [0.5, 6.15, 0.12], [0.6, 6.85, 0.20],
                    [0.65, 7.22, 0.18], [0.7, 7.67, 0.15], [0.75, 7.93, 0.23],
                    [0.8, 8.13, 0.14], [0.9, 8.52, 0.17]],
        "Lambda[0z1]": [[0.5, 6.21, 0.15], [0.6, 6.20, 0.12], [0.7, 6.40, 0.14],
                        [0.8, 6.36, 0.15], [0.9, 6.32, 0.16]],
        "pi[0z1]": [[0.1, 8.52, 0.20], [0.2, 8.39, 0.15], [0.3, 8.16, 0.16],
                    [0.4, 7.83, 0.14], [0.5, 7.49, 0.14], [0.6, 7.11, 0.13],
                    [0.7, 6.80, 0.12], [0.8, 6.47, 0.12], [0.9, 6.40, 0.11]],
    }},

    #  ---- copper, 296 K, nine branches ---------------------------------
    #  Table I, the Chalk River set.  The paper's Table II is a finer scan of
    #  the low-zeta part of T1[0zz] taken on a different spectrometer at
    #  McMaster; the two agree within their combined errors where they
    #  overlap (1.31 against 1.35 at zeta = 0.2, 2.00 against 2.03 at 0.3),
    #  and only Table I is carried here so that one branch is not drawn twice
    #  from two instruments.  The paper also warns that its quoted errors are
    #  probably overestimates by up to a factor of two, so the bars here are
    #  conservative rather than tight.
    "Cu": {"T_K": 296, "ref": CU_REF, "branches": {
        "T[00z]": [[0.15, 1.17, 0.04], [0.2, 1.56, 0.04], [0.25, 1.92, 0.04],
                   [0.275, 2.12, 0.04], [0.3, 2.30, 0.04], [0.35, 2.64, 0.04],
                   [0.4, 3.01, 0.04], [0.45, 3.30, 0.05], [0.5, 3.62, 0.04],
                   [0.55, 3.88, 0.05], [0.6, 4.15, 0.05], [0.65, 4.34, 0.05],
                   [0.7, 4.54, 0.05], [0.75, 4.73, 0.06], [0.8, 4.86, 0.07],
                   [0.9, 5.02, 0.07], [1.0, 5.08, 0.08]],
        "L[00z]": [[0.15, 1.90, 0.09], [0.2, 2.42, 0.07], [0.25, 3.02, 0.08],
                   [0.3, 3.56, 0.06], [0.4, 4.47, 0.07], [0.5, 5.32, 0.07],
                   [0.6, 6.05, 0.08], [0.65, 6.37, 0.08], [0.7, 6.60, 0.08],
                   [0.75, 6.77, 0.12], [0.8, 6.99, 0.13], [0.85, 7.14, 0.14],
                   [0.9, 7.17, 0.12], [1.0, 7.19, 0.12]],
        "T2[0zz]": [[0.1, 1.11, 0.03], [0.15, 1.69, 0.04], [0.2, 2.27, 0.04],
                    [0.25, 2.82, 0.04], [0.3, 3.37, 0.04], [0.4, 4.30, 0.05],
                    [0.5, 5.07, 0.06], [0.6, 5.71, 0.06], [0.65, 6.04, 0.06],
                    [0.7, 6.31, 0.07], [0.75, 6.54, 0.10], [0.8, 6.80, 0.11],
                    [0.9, 7.13, 0.15], [1.0, 7.19, 0.12]],
        "T1[0zz]": [[0.2, 1.35, 0.04], [0.3, 2.03, 0.04], [0.4, 2.70, 0.04],
                    [0.5, 3.34, 0.04], [0.6, 3.89, 0.05], [0.7, 4.34, 0.05],
                    [0.75, 4.55, 0.05], [0.8, 4.75, 0.07], [0.9, 5.03, 0.08],
                    [1.0, 5.08, 0.08]],
        "L[0zz]": [[0.1, 2.03, 0.10], [0.2, 3.70, 0.08], [0.3, 5.11, 0.07],
                   [0.4, 5.97, 0.08], [0.5, 6.36, 0.10], [0.6, 6.38, 0.12],
                   [0.7, 5.91, 0.10], [0.75, 5.73, 0.08], [0.8, 5.51, 0.07],
                   [0.9, 5.19, 0.07], [1.0, 5.08, 0.08]],
        "L[zzz]": [[0.05, 1.24, 0.06], [0.075, 1.86, 0.07], [0.1, 2.46, 0.07],
                   [0.125, 2.99, 0.06], [0.15, 3.59, 0.06], [0.2, 4.54, 0.06],
                   [0.25, 5.43, 0.07], [0.3, 6.14, 0.07], [0.35, 6.67, 0.08],
                   [0.4, 7.06, 0.10], [0.45, 7.25, 0.13], [0.5, 7.40, 0.13]],
        "T[zzz]": [[0.075, 0.79, 0.04], [0.1, 1.01, 0.05], [0.125, 1.23, 0.06],
                   [0.15, 1.47, 0.06], [0.2, 1.87, 0.06], [0.25, 2.29, 0.05],
                   [0.3, 2.66, 0.06], [0.35, 2.97, 0.06], [0.4, 3.17, 0.07],
                   [0.45, 3.34, 0.07], [0.5, 3.37, 0.07]],
        "pi[0z1]": [[0.0, 7.19, 0.12], [0.1, 7.17, 0.15], [0.2, 7.07, 0.14],
                    [0.3, 6.80, 0.09], [0.4, 6.44, 0.09], [0.5, 6.10, 0.08],
                    [0.6, 5.77, 0.08], [0.7, 5.47, 0.08], [0.8, 5.27, 0.08],
                    [0.9, 5.13, 0.08], [1.0, 5.08, 0.08]],
        "Lambda[0z1]": [[0.0, 5.08, 0.08], [0.1, 5.03, 0.08], [0.2, 4.99, 0.07],
                        [0.3, 4.97, 0.08], [0.4, 4.89, 0.09], [0.5, 4.89, 0.08]],
    }},

    #  ---- barium, 295 K, BCC, single crystals --------------------------
    #  The one entry here that is not face-centred, and the reason it exists:
    #  Buchenau et al. reached barium only through a Born-von Karman fit to
    #  polycrystalline time-of-flight spectra, because crystals could not be
    #  grown.  Mizuki et al. grew them.  So barium is a measurement and not a
    #  model, and the model is not carried here at all.
    #  The anomaly the paper is about: between zeta = 0.2 and 0.8 along
    #  [z00] the LONGITUDINAL branch lies BELOW the transverse, which is not
    #  how a simple metal behaves and which they attribute to d-state
    #  occupation.
    "Ba": {"T_K": 295, "struct": "bcc", "ref": BA_REF, "branches": {
        "L[z00]": [[0.1, 0.28, 0.03], [0.15, 0.39, 0.05], [0.2, 0.56, 0.07],
                   [0.3, 0.82, 0.07], [0.4, 1.05, 0.07], [0.45, 1.16, 0.05],
                   [0.5, 1.22, 0.12], [0.55, 1.38, 0.03], [0.6, 1.44, 0.12],
                   [0.7, 1.74, 0.05], [0.8, 2.01, 0.05], [0.85, 2.09, 0.05],
                   [0.9, 2.10, 0.04], [0.95, 2.13, 0.05]],
        "T[z00]": [[0.1, 0.31, 0.02], [0.15, 0.50, 0.03], [0.2, 0.65, 0.02],
                   [0.3, 0.98, 0.04], [0.4, 1.24, 0.03], [0.45, 1.40, 0.03],
                   [0.5, 1.52, 0.04], [0.55, 1.64, 0.05], [0.6, 1.73, 0.03],
                   [0.65, 1.85, 0.04], [0.7, 1.93, 0.05], [0.75, 1.98, 0.05],
                   [0.8, 2.05, 0.05], [0.9, 2.15, 0.05], [1.0, 2.15, 0.07]],
        "T1[zz0]": [[0.1, 0.27, 0.02], [0.15, 0.38, 0.02], [0.2, 0.50, 0.03],
                    [0.25, 0.59, 0.03], [0.3, 0.65, 0.03], [0.35, 0.71, 0.03],
                    [0.4, 0.76, 0.04], [0.45, 0.76, 0.04], [0.5, 0.75, 0.04]],
        "T2[zz0]": [[0.1, 0.47, 0.02], [0.15, 0.69, 0.02], [0.2, 0.88, 0.03],
                    [0.3, 1.23, 0.03], [0.4, 1.44, 0.03], [0.5, 1.53, 0.04]],
        "L[zz0]": [[0.1, 0.60, 0.04], [0.2, 1.25, 0.03], [0.3, 1.81, 0.04],
                   [0.4, 2.10, 0.05], [0.5, 2.30, 0.07]],
        "T[zzz]": [[0.15, 0.615, 0.02], [0.2, 0.80, 0.04], [0.3, 1.16, 0.03],
                   [0.4, 1.47, 0.03], [0.45, 1.54, 0.03], [0.5, 1.72, 0.05],
                   [0.6, 1.90, 0.04], [0.7, 1.99, 0.03], [0.8, 2.05, 0.10],
                   [0.9, 2.08, 0.08]],
        "L[zzz]": [[0.1, 0.78, 0.04], [0.2, 1.56, 0.05], [0.3, 2.12, 0.08],
                   [0.35, 2.10, 0.09], [0.4, 1.97, 0.06], [0.5, 1.72, 0.03],
                   [0.52, 1.60, 0.04], [0.575, 1.40, 0.04], [0.667, 1.11, 0.04],
                   [0.79, 1.40, 0.05], [0.8, 1.36, 0.05], [0.9, 1.90, 0.06]],
    }},

    #  ---- lead, 100 K --------------------------------------------------
    #  The one set here whose zeta is |q| a / 2pi rather than the component
    #  along the direction.  The paper says so itself - Table II is headed
    #  "the zone boundary is at aq/2pi = 1.061", and 1.061 is 3 sqrt(2)/4,
    #  which is K.  So L sits at sqrt(3)/2 = 0.866 and the [zzz] rows run to
    #  0.867 rather than to 0.5.  The numbers are left exactly as printed and
    #  the segment ranges are moved instead; converting them would have made
    #  the table impossible to check against the paper by eye.
    #
    #  100 K, and that matters more here than anywhere else in this file:
    #  lead is soft, its Debye temperature is about 105 K, and the library's
    #  lattice constant is the 293 K one.  Some of the disagreement on this
    #  element is that mismatch rather than the potential.
    #  Same group and the same spectrometer as barium, and the same seven
    #  branches, so the bcc segment table above needs nothing added for it.
    #  The internal check that the labels are right: L and T meet at 0.70 THz
    #  at zeta = 1/2, which is P, where bcc branches must be degenerate, and
    #  both run up to the 0.895 THz that T[z00] reaches at H.
    #  Beryllium is the one element here whose paper reports frequencies only
    #  AT the symmetry points - Table 1 is a list of points, not of branches -
    #  so there is no dispersion curve to draw and the page says so.
    #
    #  Each measured frequency appears ONCE, under the paper's own
    #  irreducible-representation label.  Listing it once per degenerate mode
    #  instead would mean deciding that K5 stands for two modes and A3 for
    #  four, and that decision is group theory rather than measurement.  It
    #  is not needed to draw the point, so it is not made here.
    #
    #  These are NOT merged into refdata_phonon.  That set is deliberately
    #  one compiler and the same four points for every element, which is what
    #  makes its per-element error comparable across the table; giving
    #  beryllium six points and everything else four would quietly end that.
    #  "the lattice constants determined by us at 80 K were a = 2.2825 +/-
    #  0.0005 A and c = 3.579 +/- 0.001 A, which may be compared with the
    #  values for 300 K quoted by a modern reference work, 2.2858 and 3.5843".
    #  That 300 K pair is what refdata carries, to the digit - a check on its
    #  own note saying a0 is a ~293 K quantity.
    "Be": {"T_K": 80, "struct": "hcp", "ref": BE_REF,
           "a_meas": 2.2825, "c_meas": 3.5790, "points": {
        "G": [(20.28, 0.04, "G3+"), (13.73, 0.03, "G5+")],
        "M": [(17.79, 0.04, "M1-"), (17.45, 0.05, "M2-"), (16.73, 0.03, "M3-"),
              (12.41, 0.03, "M3+"), (16.90, 0.03, "M4-"), (11.86, 0.03, "M4+")],
        "K": [(18.80, 0.05, "K1"), (15.97, 0.05, "K3"), (14.65, 0.05, "K5"),
              (14.59, 0.03, "K6")],
        "A": [(14.90, 0.03, "A1"), (10.72, 0.04, "A3")],
        "L": [(19.04, 0.04, "L1+"), (15.07, 0.03, "L1-"), (14.23, 0.03, "L2")],
        "H": [(19.10, 0.04, "H1"), (17.30, 0.04, "H3+"), (12.30, 0.03, "H3-")],
    }},
    #  The worst of the five, and it runs the other way: measured at 280 K
    #  and scored at refdata's a0, which for caesium is Kittel's **5 K**
    #  value - a crystal 7.5 % of volume too small.  a(280 K) from the same
    #  expansion integral, carried UP rather than down.
    "Cs": {"T_K": 280, "a_meas": 6.1962, "a_meas_from": "expansion integral", "struct": "bcc", "ref": CS_REF,
           "branches": {
        "L[z00]": [[0.2, 0.29, 0.03], [0.3, 0.49, 0.04], [0.4, 0.65, 0.04],
                   [0.5, 0.73, 0.04], [0.6, 0.88, 0.04], [0.8, 0.93, 0.04]],
        "T[z00]": [[0.1, 0.16, 0.01], [0.15, 0.18, 0.01], [0.2, 0.26, 0.03],
                   [0.3, 0.34, 0.02], [0.4, 0.46, 0.02], [0.5, 0.58, 0.03],
                   [0.55, 0.62, 0.03], [0.6, 0.65, 0.05], [0.7, 0.75, 0.03],
                   [0.8, 0.82, 0.05], [0.85, 0.85, 0.05], [1.0, 0.895, 0.05]],
        "T1[zz0]": [[0.1, 0.097, 0.01], [0.15, 0.121, 0.01],
                    [0.2, 0.155, 0.01], [0.25, 0.174, 0.01],
                    [0.3, 0.208, 0.01], [0.35, 0.206, 0.02],
                    [0.4, 0.225, 0.02], [0.45, 0.230, 0.02],
                    [0.5, 0.220, 0.03]],
        "L[zz0]": [[0.1, 0.31, 0.01], [0.15, 0.40, 0.01], [0.2, 0.53, 0.06],
                   [0.25, 0.70, 0.04], [0.3, 0.82, 0.02], [0.35, 0.91, 0.03],
                   [0.4, 1.00, 0.03], [0.5, 1.07, 0.04]],
        "T2[zz0]": [[0.1, 0.21, 0.01], [0.15, 0.28, 0.01], [0.2, 0.34, 0.02],
                    [0.25, 0.39, 0.02], [0.3, 0.46, 0.02], [0.4, 0.56, 0.02],
                    [0.5, 0.58, 0.03]],
        "T[zzz]": [[0.1, 0.16, 0.01], [0.15, 0.22, 0.01], [0.2, 0.28, 0.02],
                   [0.3, 0.41, 0.03], [0.4, 0.54, 0.06], [0.5, 0.70, 0.06],
                   [0.55, 0.76, 0.07], [0.60, 0.80, 0.08], [0.7, 0.90, 0.06],
                   [0.8, 0.89, 0.06], [0.9, 0.91, 0.04]],
        "L[zzz]": [[0.1, 0.32, 0.02], [0.15, 0.53, 0.03], [0.2, 0.65, 0.03],
                   [0.3, 0.90, 0.04], [0.35, 0.91, 0.05], [0.4, 0.83, 0.04],
                   [0.45, 0.79, 0.02], [0.5, 0.70, 0.06], [0.55, 0.54, 0.04],
                   [0.6, 0.33, 0.04], [0.667, 0.30, 0.03], [0.75, 0.39, 0.02],
                   [0.8, 0.54, 0.02], [0.9, 0.79, 0.02]],
    }},
    #  The same two structural checks the caesium table passes: L and T meet
    #  at 3.78 THz at zeta = 1/2 along [zzz], which is P, and all four
    #  branches of the upper block meet at 5.03 THz at zeta = 1, which is H.
    #  Entries printed as "..." in the table are absent here, not zero.
    #  The first hcp element here, so six branches rather than three, and the
    #  paper's own labels are kept: perp and par are the polarisations
    #  perpendicular and parallel to the basal plane, which is a distinction
    #  a cubic lattice does not have.
    #  The structural check: LA and LO both reach 5.73 THz at zeta = 1/2
    #  along [00z], which is A, where hcp branches must meet in pairs.
    #  Two things differ from titanium and both are handled by fields rather
    #  than by editing the numbers: energies in meV (nu_scale) and zeta
    #  normalised to the zone boundary (seg_override).  Delta is the
    #  Gamma-A line and Sigma the Gamma-M line, which is the paper's own
    #  labelling and is why [0110] appears here as [z00].
    #
    #  The check that both blocks describe the same crystal: at q = 0 the
    #  optic modes come out at 18.96 meV and 10.77 meV from the Delta
    #  columns and at 18.96 and 10.77 from the Sigma columns.
    #
    #  Sigma_3 and Sigma_4 are the two transverse representations.  Which of
    #  them is polarised along c is NOT recorded here, because the paper
    #  labels them by representation and guessing would put a real datum
    #  under a made-up name.
    #  This table validates the [zz0] folding above, and it does it with
    #  measurements rather than with an argument.  M is reached twice - once
    #  along [z00] at zeta = 1/2 and once along [zz0] at zeta = 1/2, which is
    #  the far end of the folded part - and the six frequencies agree:
    #      along [z00]   2.08  2.35  3.50  3.62  3.85  4.27
    #      along [zz0]   2.08  2.32  3.50  3.62  3.85  4.27
    #  If [zz0] beyond K had been drawn as more of Gamma-K instead, those two
    #  sets would have had no reason to match.
    #
    #  The perpendicular/parallel labelling is checked the same way rather
    #  than read off a blurred subscript: at Gamma the optic mode polarised
    #  along c is 3.1 THz (it is LO along [00z], where c IS the longitudinal
    #  direction) and the basal E2g pair is 2.7 THz.  So the [z00] optic
    #  branch starting at 3.1 is the perpendicular one and the one starting
    #  at 2.7 is parallel - which the branches then confirm by joining up
    #  across the page break.
    #
    #  LA[zz0] also has a point at zeta = 0.054 printed with NO error.  It is
    #  left out rather than given a zero error bar, which would draw it as
    #  the most precise point on the plot.
    "Hf": {"T_K": 295, "struct": "hcp", "ref": HF_REF, "branches": {
        "TA[00z]": [[0.1, 0.42, 0.01], [0.2, 0.90, 0.02], [0.3, 1.21, 0.03],
                    [0.4, 1.58, 0.03], [0.5, 1.90, 0.10]],
        "TO[00z]": [[0.0, 2.70, 0.16], [0.1, 2.74, 0.12], [0.2, 2.54, 0.12],
                    [0.3, 2.46, 0.12], [0.4, 2.16, 0.10], [0.5, 1.90, 0.10]],
        "LA[00z]": [[0.1, 0.76, 0.02], [0.2, 1.60, 0.04], [0.3, 2.10, 0.10],
                    [0.4, 2.60, 0.10], [0.5, 3.00, 0.15]],
        "LO[00z]": [[0.0, 3.10, 0.15], [0.05, 3.15, 0.15], [0.1, 3.10, 0.2],
                    [0.2, 3.40, 0.15], [0.25, 3.58, 0.15], [0.3, 3.70, 0.15],
                    [0.35, 3.56, 0.15], [0.4, 3.50, 0.10], [0.45, 3.40, 0.10],
                    [0.5, 3.00, 0.15]],
        "LA[z00]": [[0.1, 1.31, 0.02], [0.2, 2.36, 0.10], [0.3, 3.29, 0.10],
                    [0.4, 3.80, 0.15], [0.5, 3.85, 0.15]],
        "TApar[z00]": [[0.15, 1.20, 0.05], [0.2, 1.41, 0.04], [0.3, 1.73, 0.03],
                       [0.4, 2.03, 0.15], [0.5, 2.08, 0.06]],
        "TAperp[z00]": [[0.1, 0.79, 0.03], [0.15, 1.18, 0.02],
                        [0.2, 1.48, 0.02], [0.25, 1.78, 0.03],
                        [0.3, 1.91, 0.04], [0.4, 2.16, 0.12],
                        [0.5, 2.35, 0.15]],
        "TOperp[z00]": [[0.0, 3.10, 0.15], [0.05, 3.12, 0.12],
                        [0.15, 3.21, 0.2], [0.2, 3.23, 0.15], [0.25, 3.37, 0.2],
                        [0.3, 3.52, 0.1], [0.4, 3.45, 0.1], [0.5, 3.50, 0.12]],
        "TOpar[z00]": [[0.0, 2.70, 0.16], [0.1, 2.70, 0.2], [0.15, 2.78, 0.12],
                       [0.2, 2.90, 0.15], [0.3, 3.40, 0.15], [0.35, 3.50, 0.15],
                       [0.4, 3.55, 0.15], [0.5, 3.62, 0.12]],
        "LO[z00]": [[0.0, 2.70, 0.16], [0.05, 2.75, 0.10], [0.1, 2.95, 0.10],
                    [0.15, 3.25, 0.10], [0.2, 3.45, 0.10], [0.25, 3.70, 0.10],
                    [0.3, 3.80, 0.10], [0.5, 4.27, 0.05]],
        "TAperp[zz0]": [[0.1, 1.37, 0.07], [0.2, 2.28, 0.15],
                        [0.333, 3.13, 0.15], [0.5, 3.50, 0.15]],
        "TApar[zz0]": [[0.05, 0.61, 0.03], [0.1, 1.39, 0.02],
                       [0.15, 1.87, 0.05], [0.2, 2.56, 0.06], [0.25, 2.97, 0.1],
                       [0.33, 3.32, 0.16], [0.4, 2.90, 0.06], [0.5, 2.08, 0.06]],
        "LA[zz0]": [[0.1, 2.27, 0.2], [0.15, 2.63, 0.10], [0.2, 2.90, 0.10],
                    [0.25, 3.15, 0.10], [0.3, 3.40, 0.15], [0.35, 3.50, 0.15],
                    [0.4, 3.52, 0.15], [0.45, 3.58, 0.12], [0.5, 3.62, 0.12]],
        "LO[zz0]": [[0.0, 2.70, 0.16], [0.15, 3.60, 0.15], [0.25, 3.72, 0.12],
                    [0.5, 3.85, 0.15]],
        "TOperp[zz0]": [[0.0, 3.10, 0.15], [0.1, 3.16, 0.12], [0.2, 3.57, 0.12],
                        [0.4, 2.82, 0.10], [0.5, 2.32, 0.1]],
        "TOpar[zz0]": [[0.0, 2.70, 0.16], [0.15, 3.30, 0.1], [0.4, 4.03, 0.15],
                       [0.5, 4.27, 0.05]],
    }},
    #  Same group and the same table layout as scandium, and the same M-point
    #  check: 4.14 THz and 2.67 THz each appear once in the [z00] column at
    #  zeta = 1/2 and once in the [zz0] column at zeta = 1/2.  T2 has a point
    #  printed at 0.333, which is K.
    "Y": {"T_K": 295, "struct": "hcp", "ref": Y_REF, "branches": {
        "LA[00z]": [[0.2, 1.37, 0.02], [0.225, 1.54, 0.02], [0.25, 1.71, 0.02],
                    [0.275, 1.90, 0.02], [0.3, 2.05, 0.02], [0.325, 2.21, 0.02],
                    [0.35, 2.35, 0.02], [0.375, 2.52, 0.02], [0.4, 2.66, 0.02],
                    [0.425, 2.82, 0.03], [0.45, 2.96, 0.03], [0.475, 3.05, 0.03],
                    [0.5, 3.20, 0.04]],
        "LO[00z]": [[0.0, 4.64, 0.05], [0.05, 4.65, 0.06], [0.1, 4.60, 0.05],
                    [0.2, 4.40, 0.05], [0.225, 4.27, 0.05], [0.25, 4.28, 0.05],
                    [0.275, 4.24, 0.04], [0.3, 4.10, 0.05], [0.325, 4.07, 0.05],
                    [0.35, 3.88, 0.04], [0.375, 3.79, 0.03], [0.4, 3.84, 0.04],
                    [0.45, 3.66, 0.02], [0.475, 3.42, 0.03]],
        "TA[00z]": [[0.25, 0.95, 0.02], [0.275, 1.05, 0.02], [0.4, 1.67, 0.04],
                    [0.425, 1.80, 0.05]],
        "TO[00z]": [[0.0, 2.68, 0.03], [0.1, 2.69, 0.02], [0.2, 2.48, 0.03],
                    [0.25, 2.47, 0.03], [0.3, 2.27, 0.05], [0.35, 2.23, 0.05],
                    [0.4, 2.20, 0.02], [0.5, 1.96, 0.02]],
        "LA[z00]": [[0.1, 1.19, 0.02], [0.15, 1.93, 0.05], [0.2, 2.48, 0.06],
                    [0.25, 3.19, 0.06], [0.3, 3.45, 0.06], [0.35, 3.80, 0.05],
                    [0.4, 3.83, 0.04], [0.45, 3.91, 0.04], [0.5, 4.02, 0.05]],
        "LO[z00]": [[0.0, 2.68, 0.03], [0.05, 2.71, 0.05], [0.075, 2.86, 0.03],
                    [0.1, 2.94, 0.02], [0.125, 3.12, 0.03], [0.15, 3.23, 0.05],
                    [0.175, 3.33, 0.03], [0.2, 3.53, 0.05], [0.225, 3.66, 0.04],
                    [0.25, 3.72, 0.08], [0.275, 3.83, 0.04], [0.3, 3.78, 0.05],
                    [0.35, 3.99, 0.05], [0.4, 4.08, 0.05], [0.45, 4.10, 0.05]],
        "TOperp[z00]": [[0.0, 4.64, 0.05], [0.1, 4.65, 0.06], [0.2, 4.58, 0.06],
                        [0.3, 4.42, 0.04], [0.4, 4.28, 0.04], [0.5, 4.14, 0.04]],
        "TAperp[z00]": [[0.3, 2.25, 0.05], [0.4, 2.45, 0.03], [0.5, 2.67, 0.04]],
        "TApar[z00]": [[0.3, 2.03, 0.02], [0.4, 2.29, 0.03], [0.5, 2.30, 0.05]],
        "TOpar[z00]": [[0.3, 3.53, 0.04], [0.4, 3.90, 0.05], [0.5, 4.04, 0.05]],
        "TAperp[zz0]": [[0.2, 2.36, 0.02], [0.302, 3.19, 0.02],
                        [0.406, 3.73, 0.03], [0.5, 4.14, 0.04]],
        "TOperp[zz0]": [[0.0, 4.64, 0.05], [0.1, 4.53, 0.05], [0.2, 4.31, 0.03],
                        [0.25, 4.06, 0.02], [0.3, 3.69, 0.02],
                        [0.333, 3.43, 0.02], [0.36, 3.25, 0.02],
                        [0.4, 2.91, 0.02], [0.5, 2.67, 0.04]],
    }},
    "Zr": {"T_K": 295, "struct": "hcp", "ref": ZR_REF, "nu_scale": MEV_THZ,
           "seg_override": {"[00z]": [("G", "A", 0.0, 1.0)],
                            "[z00]": [("G", "M", 0.0, 1.0)]},
           "branches": {
        "LA[00z]": [[0.2, 3.75, 0.06], [0.4, 7.90, 0.06], [0.6, 11.46, 0.08],
                    [0.8, 14.75, 0.08], [1.0, 16.88, 0.09]],
        "LO[00z]": [[0.0, 18.96, 0.04], [0.2, 18.84, 0.05], [0.4, 19.02, 0.05],
                    [0.6, 18.97, 0.05], [0.8, 18.66, 0.05], [1.0, 16.88, 0.09]],
        "TA[00z]": [[0.2, 1.95, 0.07], [0.3, 2.92, 0.07], [0.4, 3.84, 0.03],
                    [0.5, 4.52, 0.04], [0.6, 5.23, 0.07], [0.8, 6.72, 0.03],
                    [1.0, 7.77, 0.03]],
        "TO[00z]": [[0.0, 10.77, 0.05], [0.2, 10.67, 0.06], [0.4, 10.02, 0.05],
                    [0.6, 9.47, 0.06], [0.8, 8.69, 0.05], [1.0, 7.77, 0.03]],
        "LA[z00]": [[0.2, 6.52, 0.05], [0.4, 13.03, 0.08], [0.6, 17.22, 0.05],
                    [0.8, 20.87, 0.11], [1.0, 23.28, 0.14]],
        "LO[z00]": [[0.0, 10.77, 0.05], [0.2, 13.18, 0.12], [0.4, 17.28, 0.15],
                    [0.6, 20.84, 0.09], [0.8, 22.18, 0.09], [1.0, 22.87, 0.07]],
        "TA_S3[z00]": [[0.2, 3.60, 0.03], [0.3, 5.14, 0.05], [0.4, 6.66, 0.04],
                       [0.5, 7.90, 0.07], [0.6, 9.16, 0.16], [0.8, 10.74, 0.05],
                       [1.0, 11.13, 0.06]],
        "TO_S3[z00]": [[0.0, 18.96, 0.04], [0.2, 18.92, 0.06],
                       [0.4, 19.02, 0.08], [0.6, 19.65, 0.10],
                       [0.8, 19.78, 0.06], [1.0, 19.85, 0.08]],
        "TA_S4[z00]": [[0.2, 3.92, 0.03], [0.3, 5.23, 0.05], [0.4, 6.65, 0.03],
                       [0.5, 7.78, 0.04], [0.6, 8.35, 0.06], [0.8, 9.30, 0.04],
                       [1.0, 9.77, 0.04]],
        "TO_S4[z00]": [[0.0, 10.77, 0.05], [0.2, 12.25, 0.03],
                       [0.4, 15.13, 0.03], [0.6, 17.70, 0.04],
                       [0.8, 19.17, 0.06], [1.0, 19.78, 0.08]],
    }},
    #  Thallium is measured on two lines only, and III/13a says so outright:
    #  "only the dispersion in the symmetry directions Delta and Sigma has
    #  been measured".  So this record carries Gamma-A and Gamma-M and
    #  nothing else, which is honest rather than thin - there is no more.
    #
    #  The abscissa is q/q_max, 1 at the zone boundary, not the reduced
    #  component this file uses, where M and A both sit at 0.5.  The
    #  seg_override below carries that factor of two; the frequencies are
    #  untouched.
    #
    #  Check the transcription did not use: at A the table gives Delta_1(LA)
    #  and Delta_2(LO) both as 1.84, and Delta_5(TO) and Delta_6(TA) both as
    #  0.76.  That is the hcp A-point sticking, and it shows up independently
    #  in four branches, so the two columns are lined up correctly.
    #
    #  The error bars are large - 0.36 THz on 1.28, twenty-eight per cent -
    #  and they are the source's own.  Read this element's agreement with
    #  that in mind rather than as a tight test.
    #  Measured at 77 K, scored at a ~293 K a0 until 2026-08-29.  Only a is
    #  corrected: the expansion integral has one polycrystalline alpha and
    #  cannot say how c/a moves, so the axial ratio is left alone.
    "Tl": {"T_K": 77, "a_meas": 3.4385, "a_meas_from": "expansion integral", "struct": "hcp", "ref": TL_REF,
           "seg_override": {"[z00]": [("G", "M", 0.0, 1.0)],
                            "[00z]": [("G", "A", 0.0, 1.0)]},
           "branches": {
        "LA[z00]": [[0.4, 1.21, 0.19], [0.6, 1.68, 0.10], [0.8, 2.07, 0.12],
                    [1.0, 2.16, 0.12]],
        "LO[z00]": [[0.0, 1.28, 0.36], [0.2, 1.38, 0.36], [0.4, 1.63, 0.10],
                    [0.6, 1.99, 0.07], [0.8, 2.25, 0.07], [1.0, 2.32, 0.07]],
        "TAperp[z00]": [[0.4, 0.57, 0.12], [0.6, 0.85, 0.07], [0.8, 1.02, 0.07],
                        [1.0, 1.03, 0.10]],
        "TOperp[z00]": [[0.0, 2.84, 0.12], [0.2, 2.82, 0.07], [0.4, 2.76, 0.07],
                        [0.6, 2.67, 0.07], [0.8, 2.55, 0.12], [1.0, 2.40, 0.12]],
        "TApar[z00]": [[0.4, 0.41, 0.19], [0.6, 0.65, 0.15], [0.8, 0.89, 0.12],
                       [1.0, 0.97, 0.15]],
        "TOpar[z00]": [[0.2, 1.28, 0.36], [0.4, 1.31, 0.36], [0.6, 1.39, 0.12],
                       [0.8, 1.50, 0.10], [1.0, 1.58, 0.07]],
        "LA[00z]": [[0.4, 0.80, 0.12], [0.6, 1.20, 0.15], [0.8, 1.50, 0.12],
                    [1.0, 1.84, 0.10]],
        "LO[00z]": [[0.2, 2.78, 0.07], [0.4, 2.61, 0.07], [0.6, 2.27, 0.17],
                    [0.8, 2.04, 0.10], [1.0, 1.84, 0.10]],
        "TA[00z]": [[0.4, 0.30, 0.10], [0.6, 0.45, 0.12], [0.8, 0.58, 0.22],
                    [1.0, 0.76, 0.15]],
        "TO[00z]": [[0.2, 1.12, 0.22], [0.4, 1.11, 0.22], [0.6, 0.97, 0.12],
                    [0.8, 0.85, 0.10], [1.0, 0.76, 0.15]],
    }},
    "Ti": {"T_K": 295, "struct": "hcp", "ref": TI_REF, "branches": {
        "TA[00z]": [[0.1, 0.74, 0.03], [0.2, 1.38, 0.03], [0.3, 2.03, 0.03],
                    [0.4, 2.56, 0.03], [0.5, 3.05, 0.03]],
        "TO[00z]": [[0.0, 4.10, 0.15], [0.1, 4.07, 0.15], [0.2, 3.95, 0.10],
                    [0.3, 3.79, 0.05], [0.4, 3.43, 0.05]],
        "LA[00z]": [[0.1, 1.45, 0.05], [0.2, 2.74, 0.06], [0.3, 3.91, 0.05],
                    [0.4, 4.97, 0.06], [0.5, 5.73, 0.05]],
        "LO[00z]": [[0.0, 5.54, 0.15], [0.1, 5.67, 0.15], [0.2, 5.75, 0.20],
                    [0.3, 6.06, 0.20], [0.4, 6.28, 0.12], [0.5, 5.73, 0.05]],
        "TAperp[z00]": [[0.1, 1.29, 0.03], [0.2, 2.46, 0.05],
                        [0.3, 3.32, 0.10], [0.4, 3.78, 0.10], [0.5, 3.82, 0.1]],
        "TApar[z00]": [[0.1, 1.13, 0.03], [0.2, 2.12, 0.02], [0.3, 2.83, 0.04],
                       [0.4, 3.15, 0.03], [0.5, 3.40, 0.04]],
        "LA[z00]": [[0.05, 1.16, 0.03], [0.1, 2.45, 0.1], [0.2, 4.45, 0.1],
                    [0.3, 5.80, 0.10], [0.4, 6.58, 0.2], [0.5, 7.10, 0.06]],
        "TOperp[z00]": [[0.1, 5.50, 0.2], [0.2, 5.62, 0.15], [0.3, 5.86, 0.15],
                        [0.4, 6.12, 0.12], [0.5, 6.06, 0.15]],
        "TOpar[z00]": [[0.1, 4.64, 0.08], [0.2, 5.45, 0.04], [0.3, 6.24, 0.06],
                       [0.4, 6.62, 0.08], [0.5, 6.95, 0.12]],
        "LO[z00]": [[0.1, 4.70, 0.06], [0.2, 6.26, 0.04], [0.4, 7.59, 0.04],
                    [0.5, 7.69, 0.04]],
        "TAperp[zz0]": [[0.05, 1.07, 0.04], [0.1, 2.10, 0.03],
                        [0.2, 3.70, 0.04], [0.3, 4.82, 0.06],
                        [0.4, 5.93, 0.04]],
        "TApar[zz0]": [[0.05, 1.09, 0.04], [0.1, 2.26, 0.05], [0.2, 4.75, 0.04],
                       [0.33, 6.00, 0.02], [0.4, 5.09, 0.02]],
        "LA[zz0]": [[0.1, 3.70, 0.08], [0.2, 4.96, 0.08], [0.33, 6.23, 0.03],
                    [0.4, 6.65, 0.03]],
        "TOperp[zz0]": [[0.1, 5.88, 0.04], [0.2, 5.90, 0.02], [0.3, 5.44, 0.02],
                        [0.4, 4.50, 0.04]],
        "TOpar[zz0]": [[0.1, 5.41, 0.04], [0.15, 6.03, 0.06], [0.2, 6.60, 0.04],
                       [0.25, 6.90, 0.03], [0.3, 7.10, 0.2], [0.33, 7.01, 0.04]],
        "LO[zz0]": [[0.1, 5.73, 0.1], [0.2, 7.30, 0.1], [0.3, 6.59, 0.1],
                    [0.4, 6.60, 0.1]],
    }},
    #  The same M-point check as hafnium, and it passes twice: the [zz0]
    #  branch reaching zeta = 1/2 gives 3.97 THz and 6.23 THz, and both
    #  numbers appear again in the [z00] column at its own zeta = 1/2.  The
    #  paper says as much in its text - it calls these the "Gamma K M"
    #  branches - so the line really does run through K and on to M.
    #  T2 even has a point printed at zeta = 0.333, which is K itself.
    "Sc": {"T_K": 295, "struct": "hcp", "ref": SC_REF, "branches": {
        "LA[00z]": [[0.105, 1.05, 0.02], [0.128, 1.31, 0.02],
                    [0.155, 1.59, 0.02], [0.175, 1.66, 0.03],
                    [0.1875, 1.94, 0.03], [0.2, 1.98, 0.03],
                    [0.2125, 2.10, 0.02], [0.225, 2.25, 0.01],
                    [0.2375, 2.29, 0.02], [0.25, 2.38, 0.02],
                    [0.2625, 2.65, 0.02], [0.275, 2.72, 0.02],
                    [0.282, 2.80, 0.02], [0.2875, 2.84, 0.02],
                    [0.293, 2.95, 0.02], [0.3, 2.97, 0.02],
                    [0.3125, 3.15, 0.02], [0.325, 3.24, 0.02],
                    [0.3375, 3.34, 0.03], [0.35, 3.43, 0.02],
                    [0.375, 3.69, 0.03], [0.3875, 3.76, 0.04],
                    [0.4, 3.89, 0.04], [0.425, 4.16, 0.02],
                    [0.45, 4.23, 0.02], [0.5, 4.74, 0.03]],
        "LO[00z]": [[0.0, 6.91, 0.03], [0.05, 6.89, 0.03], [0.1, 6.83, 0.04],
                    [0.12, 6.69, 0.04], [0.15, 6.66, 0.04], [0.172, 6.51, 0.03],
                    [0.2, 6.58, 0.03], [0.257, 6.29, 0.06], [0.272, 6.08, 0.03],
                    [0.29, 6.09, 0.02], [0.3, 6.01, 0.02], [0.34, 5.78, 0.02],
                    [0.39, 5.50, 0.02], [0.44, 5.15, 0.05], [0.491, 4.75, 0.02]],
        "TA[00z]": [[0.1, 0.62, 0.01], [0.15, 0.89, 0.01], [0.2, 1.17, 0.01],
                    [0.25, 1.47, 0.01], [0.3, 1.83, 0.04], [0.35, 2.05, 0.01],
                    [0.41, 2.38, 0.04], [0.5, 2.87, 0.02]],
        "TO[00z]": [[0.0, 4.04, 0.04], [0.05, 4.04, 0.07], [0.1, 4.14, 0.04],
                    [0.15, 4.09, 0.04], [0.2, 3.92, 0.04], [0.25, 3.83, 0.05],
                    [0.3, 3.68, 0.04], [0.35, 3.51, 0.03], [0.4, 3.35, 0.02]],
        "LA[z00]": [[0.09, 1.75, 0.02], [0.1, 2.00, 0.02], [0.125, 2.43, 0.01],
                    [0.15, 2.86, 0.01], [0.175, 3.37, 0.02], [0.2, 3.77, 0.02],
                    [0.225, 4.20, 0.02], [0.25, 4.52, 0.02], [0.275, 4.87, 0.01],
                    [0.3, 5.22, 0.02], [0.35, 5.70, 0.02], [0.5, 6.21, 0.1]],
        "LO[z00]": [[0.0, 4.04, 0.04], [0.05, 4.08, 0.04], [0.1, 4.41, 0.03],
                    [0.15, 4.71, 0.03], [0.2, 5.22, 0.02], [0.25, 5.66, 0.02],
                    [0.3, 5.93, 0.04], [0.35, 5.92, 0.04], [0.4, 6.05, 0.02],
                    [0.5, 6.23, 0.05]],
        "TAperp[z00]": [[0.1, 1.07, 0.01], [0.15, 1.59, 0.01], [0.2, 2.10, 0.01],
                        [0.3, 3.06, 0.02], [0.4, 3.75, 0.02], [0.45, 4.08, 0.03],
                        [0.5, 3.97, 0.02]],
        "TOperp[z00]": [[0.0, 6.91, 0.03], [0.05, 6.70, 0.03], [0.1, 6.58, 0.06],
                        [0.2, 6.59, 0.06], [0.25, 6.51, 0.03], [0.3, 6.38, 0.05],
                        [0.4, 6.39, 0.04], [0.45, 6.11, 0.06], [0.5, 6.23, 0.04]],
        "TApar[z00]": [[0.15, 1.60, 0.01], [0.2, 2.12, 0.01], [0.25, 2.58, 0.01],
                       [0.3, 2.96, 0.01], [0.35, 3.30, 0.03], [0.4, 3.49, 0.03],
                       [0.45, 3.69, 0.04], [0.5, 3.57, 0.04]],
        "TOpar[z00]": [[0.0, 4.04, 0.04], [0.05, 4.11, 0.03], [0.1, 4.30, 0.03],
                       [0.15, 4.55, 0.03], [0.2, 4.88, 0.03], [0.25, 5.15, 0.04],
                       [0.3, 5.49, 0.02], [0.35, 5.66, 0.03], [0.4, 5.90, 0.02],
                       [0.5, 6.11, 0.04]],
        "TAperp[zz0]": [[0.1, 1.83, 0.01], [0.15, 2.79, 0.01], [0.2, 3.61, 0.01],
                        [0.25, 4.32, 0.02], [0.3, 5.00, 0.02], [0.4, 5.86, 0.02],
                        [0.45, 6.16, 0.07], [0.5, 6.23, 0.04]],
        "TOperp[zz0]": [[0.0, 6.91, 0.03], [0.05, 6.87, 0.05], [0.1, 6.85, 0.05],
                        [0.2, 6.41, 0.04], [0.25, 6.08, 0.05], [0.3, 5.67, 0.03],
                        [0.333, 5.35, 0.03], [0.4, 4.68, 0.06],
                        [0.45, 4.19, 0.04], [0.5, 3.97, 0.02]],
    }},
    #  M is reached TWICE - at zeta = 1 of the Gamma-K-M block and at zeta = 1
    #  of the Gamma-M block - and the six frequencies are the same numbers:
    #  1.84 1.99 3.06 3.23 3.26 3.39.  That confirms the transcription and
    #  the two different q units at once.  A is reached once and pairs up as
    #  hcp requires: LA = LO = 2.41 and TA = TO = 1.47.
    #
    #  ONE POINT IS LEFT OUT.  The paper prints Delta_5(TA) = 1.90 THz at
    #  zeta = 0.3, which is above its own value at A (1.47) and above
    #  Delta_6(TO) at the same q (1.84).  Along Gamma-A the transverse
    #  acoustic branch rises to meet the optic one at A and cannot pass it,
    #  so that entry contradicts the paper's own labelling and hexagonal
    #  symmetry.  It is dropped rather than plotted as a measurement our
    #  potential would then be scored against; the other two points of that
    #  branch are kept as printed.
    #  Two bcc checks, both from the table itself: all four branches reach
    #  6.49 THz at zeta = 1, which is H and triply degenerate, and along
    #  [zzz] the L and T branches are both 5.04 at zeta = 1/2, which is P.
    #
    #  The table mixes two techniques.  Rows printed as "0.031 +- 0.004  0.50"
    #  carry the error on ZETA - they were taken at constant energy transfer -
    #  and rows printed as "0.10  1.63 +- 0.06" carry it on the frequency.
    #  Only the second kind is kept: this structure has one error field and it
    #  belongs to nu, so storing a zeta error there would put a horizontal
    #  uncertainty on a vertical bar.  That drops a handful of points nearest
    #  Gamma on the acoustic branches and nothing else.
    #  Eight points, not a dispersion - and the paper's own two labels are the
    #  check on the mapping: it calls (0,0,1) "all modes degenerate", which is
    #  H, and (1/2,1/2,1/2) "L and T degenerate", which is P.  Both come out
    #  at 5.50 THz.
    #
    #  The (1/2,1/2,0)L entry carries the authors' own warning - "this value
    #  may possibly be in serious error because of poorly formed neutron
    #  groups" - and is kept with it recorded here rather than dropped: it is
    #  their measurement and their caveat, and hiding it would leave the
    #  longitudinal branch at N with nothing at all.
    #
    #  T1 and T2 cross between zeta = 0.4 and N: T1 rises 4.12 -> 4.40 while
    #  T2 falls 4.30 -> 4.15.  Transcribed as printed.
    #  Woods and Chen's molybdenum, the companion to the tungsten record
    #  below: same journal volume, same group, same 296 K.  Eleven tabulated
    #  modes against tungsten's eight.
    #
    #  ZETA IS THE PAPER'S OWN and is not folded here.  Two of the [zzz]
    #  points sit past P - 0.7 and 0.93 - and SEGMENT already sends that half
    #  of Lambda to P -> H.  Folding them onto Gamma -> P by hand is exactly
    #  what would have drawn barium's second half on top of its first.
    #
    #  The two degenerate modes are listed under BOTH labels, as tungsten's
    #  are: the paper says "all polarizations degenerate" at (0,0,1) and at
    #  (1/2,1/2,1/2), and filing each under one label only would make the
    #  other branch look unmeasured there.
    "Mo": {"T_K": 296, "struct": "bcc", "ref": MO_REF, "sparse": True,
           "branches": {
        "L[z00]": [[0.6, 7.61, 0.10], [1.0, 5.51, 0.08]],
        "T[z00]": [[0.8, 5.97, 0.10], [1.0, 5.51, 0.08]],
        "L[zz0]": [[0.5, 8.22, 0.10]],
        "T1[zz0]": [[0.5, 5.75, 0.10]],
        "T2[zz0]": [[0.4, 4.85, 0.08], [0.5, 4.59, 0.08]],
        "L[zzz]": [[0.4, 7.26, 0.10], [0.5, 6.53, 0.10], [0.93, 6.28, 0.10]],
        "T[zzz]": [[0.5, 6.53, 0.10], [0.7, 7.08, 0.10]],
    }},
    "W": {"T_K": 296, "struct": "bcc", "ref": W_REF, "sparse": True,
          "branches": {
        "L[z00]": [[0.7, 6.30, 0.07], [1.0, 5.50, 0.1]],
        "T[z00]": [[1.0, 5.50, 0.1]],
        "L[zz0]": [[0.5, 6.75, 0.15]],
        "T1[zz0]": [[0.4, 4.12, 0.03], [0.5, 4.40, 0.05]],
        "T2[zz0]": [[0.4, 4.30, 0.05], [0.5, 4.15, 0.05]],
        "L[zzz]": [[0.5, 5.50, 0.05]],
        "T[zzz]": [[0.5, 5.50, 0.05]],
    }},
    #  Like tungsten and beryllium: a handful of points, not a dispersion.
    #  The paper's own consistency is visible in it - nu_2 and nu_3 are both
    #  5.65 at K, which is the degenerate pair hexagonal symmetry requires
    #  there - and nu_0 at Gamma, 3.70 THz, is the E2g mode Raman scattering
    #  also measures.
    #
    #  T_K is 295 by INFERENCE.  See MG_REF: the paper never states one.  It
    #  is carried as a number because the page needs one, and the reason it is
    #  a guess is carried with it rather than left out.
    #  This REPLACES a symmetry-point record.  What was here before was
    #  Squires' eight interpolated frequencies at Gamma, K and M - useful,
    #  but not a dispersion.  Iyengar and co-workers measured along the two
    #  symmetry lines and tabulated 64 frequencies, so magnesium now carries
    #  a curve like the rest of the set.  Squires' K values are kept below as
    #  `points`, because this table never goes to K, and so is his M value at
    #  5.32, the one branch there that he reached and this table did not.
    #
    #  The abscissa is q/q_max, 1 at the zone boundary, so the override
    #  carries a factor of two against the reduced component used here.
    #  Frequencies are in 10^12 c/s, which is THz, and are untouched.
    #
    #  Temperature: the paper says "at room temperature" and gives no number,
    #  so 295 K is this file's room-temperature convention and not a quoted
    #  value - the same caveat Squires' record carried.
    #
    #  FOUR CHECKS THE TRANSCRIPTION DID NOT USE.  Within the table, both
    #  [0001] pairs stick at A - LO and LA both reach 5.20, TO and TA both
    #  2.94 - which is the hcp A-point degeneracy, and the two lines agree at
    #  Gamma, where [0001]TO and [01-10]LO both read 3.75 and [0001]LO and
    #  [01-10]TO(perp) both read 7.30.  Against SQUIRES, who is a different
    #  experiment entirely, the three branches both of them reached at M agree
    #  to 0.02 THz or better: 6.88 against 6.88, 6.58 against 6.60, 3.70
    #  against 3.71.  Put his 5.32 and their 4.15 beside those and M has all
    #  six of its branches.
    "Mg": {"T_K": 295, "struct": "hcp", "ref": MG_REF,
           "seg_override": {"[00z]": [("G", "A", 0.0, 1.0)],
                            "[z00]": [("G", "M", 0.0, 1.0)]},
           "points": {"K": [(6.05, 0.15, "nu1"), (5.65, 0.15, "nu2"),
                            (5.65, 0.15, "nu3")],
                      "M": [(5.32, 0.10, "nu6")]},
           "branches": {
        "LO[00z]": [[0.00, 7.30, 0.15], [0.10, 7.34, 0.16], [0.20, 7.13, 0.14],
                    [0.30, 7.07, 0.14], [0.40, 7.03, 0.14], [0.50, 6.78, 0.13],
                    [0.60, 6.56, 0.13], [0.70, 6.20, 0.12], [0.80, 5.92, 0.12],
                    [0.90, 5.615, 0.11], [1.00, 5.20, 0.10]],
        "LA[00z]": [[0.50, 2.70, 0.08], [0.60, 3.20, 0.09], [0.70, 3.70, 0.11],
                    [1.00, 5.20, 0.10]],
        "TO[00z]": [[0.00, 3.75, 0.08], [0.30, 3.83, 0.08], [0.40, 3.60, 0.07],
                    [0.50, 3.71, 0.07], [0.60, 3.58, 0.06], [0.70, 3.46, 0.05],
                    [0.80, 3.34, 0.03], [0.90, 3.15, 0.03], [1.00, 2.94, 0.03]],
        "TA[00z]": [[0.50, 1.55, 0.02], [0.70, 2.03, 0.02], [0.80, 2.29, 0.02],
                    [0.85, 2.66, 0.03], [0.95, 2.82, 0.03], [1.00, 2.94, 0.03]],
        "LO[z00]": [[0.00, 3.75, 0.08], [0.10, 3.92, 0.08], [0.21, 4.20, 0.08],
                    [0.30, 4.70, 0.09], [0.40, 5.26, 0.11], [0.52, 5.69, 0.12],
                    [0.806, 6.69, 0.13], [0.90, 6.88, 0.14], [1.00, 6.88, 0.14]],
        "LA[z00]": [[0.30, 3.20, 0.07], [0.60, 5.58, 0.11], [0.69, 5.96, 0.06],
                    [0.78, 6.25, 0.06], [0.914, 6.46, 0.13], [1.00, 6.60, 0.13]],
        "TOperp[z00]": [[0.00, 7.30, 0.15], [0.10, 7.28, 0.15], [0.20, 7.18, 0.14],
                        [0.30, 7.08, 0.10], [0.60, 6.51, 0.13], [0.90, 6.25, 0.08],
                        [1.00, 6.12, 0.08]],
        "TApar[z00]": [[0.50, 2.56, 0.05], [0.60, 2.96, 0.06], [0.65, 3.28, 0.07],
                       [0.70, 3.40, 0.07], [0.90, 3.82, 0.08], [1.00, 3.71, 0.07]],
        "TAperp[z00]": [[0.50, 2.80, 0.03], [0.60, 3.30, 0.03], [0.70, 3.58, 0.04],
                        [0.80, 3.90, 0.04], [0.90, 4.10, 0.04], [1.00, 4.15, 0.04]],
    }},
    "Nb": {"T_K": 296, "struct": "bcc", "ref": NB_REF, "branches": {
        "L[z00]": [[0.10, 1.63, 0.06], [0.15, 2.39, 0.08], [0.20, 3.10, 0.08],
                   [0.25, 3.73, 0.08], [0.30, 4.35, 0.08], [0.35, 4.82, 0.06],
                   [0.375, 5.05, 0.08], [0.40, 5.30, 0.08], [0.425, 5.46, 0.08],
                   [0.45, 5.61, 0.08], [0.475, 5.70, 0.08], [0.50, 5.73, 0.10],
                   [0.525, 5.90, 0.12], [0.55, 5.88, 0.10], [0.60, 5.77, 0.12],
                   [0.65, 5.55, 0.08], [0.70, 5.49, 0.12], [0.725, 5.56, 0.10],
                   [0.75, 5.63, 0.08], [0.775, 5.67, 0.10], [0.80, 5.72, 0.12],
                   [0.825, 5.98, 0.14], [0.85, 6.03, 0.12], [0.875, 6.11, 0.12],
                   [0.90, 6.35, 0.12], [0.95, 6.40, 0.12], [1.00, 6.49, 0.10]],
        "T[z00]": [[0.143, 0.73, 0.04], [0.193, 0.93, 0.06], [0.213, 1.05, 0.04],
                   [0.233, 1.14, 0.04], [0.243, 1.19, 0.06], [0.253, 1.25, 0.06],
                   [0.273, 1.37, 0.06], [0.293, 1.51, 0.06], [0.30, 1.54, 0.06],
                   [0.313, 1.63, 0.04], [0.333, 1.78, 0.04], [0.343, 1.87, 0.06],
                   [0.353, 1.97, 0.06], [0.393, 2.28, 0.06], [0.40, 2.27, 0.08],
                   [0.443, 2.82, 0.08], [0.493, 3.38, 0.08], [0.50, 3.55, 0.08],
                   [0.60, 4.57, 0.10], [0.65, 5.04, 0.10], [0.75, 5.90, 0.10],
                   [0.80, 6.23, 0.12], [0.90, 6.42, 0.10], [0.95, 6.47, 0.10],
                   [1.00, 6.49, 0.10]],
        "L[zz0]": [[0.05, 1.13, 0.06], [0.075, 1.62, 0.06], [0.10, 2.12, 0.06],
                   [0.125, 2.55, 0.08], [0.15, 3.00, 0.08], [0.20, 3.83, 0.08],
                   [0.25, 4.32, 0.08], [0.30, 4.86, 0.08], [0.338, 5.05, 0.08],
                   [0.35, 5.12, 0.08], [0.40, 5.26, 0.12], [0.45, 5.47, 0.16],
                   [0.50, 5.66, 0.12]],
        "T1[zz0]": [[0.075, 0.78, 0.06], [0.10, 1.04, 0.06], [0.125, 1.31, 0.06],
                    [0.15, 1.57, 0.04], [0.20, 2.15, 0.04], [0.25, 2.71, 0.06],
                    [0.30, 3.24, 0.06], [0.35, 3.57, 0.06], [0.40, 3.80, 0.06],
                    [0.45, 3.88, 0.10], [0.50, 3.93, 0.06]],
        "T2[zz0]": [[0.10, 0.73, 0.08], [0.125, 0.98, 0.08], [0.15, 1.19, 0.10],
                    [0.20, 1.76, 0.10], [0.25, 2.39, 0.08], [0.30, 3.17, 0.10],
                    [0.35, 3.82, 0.10], [0.40, 4.53, 0.10], [0.45, 4.92, 0.10],
                    [0.50, 5.07, 0.10]],
        "L[zzz]": [[0.10, 2.46, 0.08], [0.15, 3.51, 0.08], [0.20, 4.34, 0.08],
                   [0.25, 5.04, 0.08], [0.275, 5.42, 0.08], [0.30, 5.71, 0.10],
                   [0.325, 5.90, 0.08], [0.35, 5.91, 0.08], [0.375, 5.77, 0.08],
                   [0.40, 5.51, 0.10], [0.425, 5.29, 0.08], [0.45, 5.21, 0.08],
                   [0.475, 5.12, 0.08], [0.50, 5.04, 0.08], [0.525, 4.89, 0.08],
                   [0.55, 4.70, 0.08], [0.575, 4.43, 0.08], [0.60, 4.23, 0.10],
                   [0.625, 3.91, 0.06], [0.65, 3.68, 0.06], [0.675, 3.51, 0.06],
                   [0.70, 3.47, 0.10], [0.725, 3.66, 0.08], [0.75, 3.80, 0.08],
                   [0.775, 3.95, 0.08], [0.80, 4.17, 0.08], [0.85, 4.99, 0.12],
                   [0.90, 5.80, 0.10], [0.95, 6.29, 0.10], [1.00, 6.49, 0.10]],
        "T[zzz]": [[0.15, 1.71, 0.06], [0.20, 2.11, 0.06], [0.25, 2.62, 0.06],
                   [0.30, 3.23, 0.08], [0.35, 3.88, 0.06], [0.40, 4.39, 0.06],
                   [0.45, 4.82, 0.08], [0.50, 5.04, 0.08], [0.55, 5.07, 0.08],
                   [0.60, 5.10, 0.10], [0.625, 5.05, 0.20], [0.65, 5.11, 0.10],
                   [0.70, 5.20, 0.10], [0.75, 5.45, 0.12], [0.775, 5.60, 0.16],
                   [0.80, 5.80, 0.16], [0.85, 6.12, 0.10], [0.90, 6.30, 0.12],
                   [0.95, 6.41, 0.10], [1.00, 6.49, 0.10]],
    }},
    #  q is in inverse angstroms, so the mapping is fixed by the lattice
    #  constant: zeta = q a / 2 pi with a = 2.8665 A, which puts 2 pi / a at
    #  2.192.  The table confirms it without being asked - the [zz0] columns
    #  stop at 1.096, exactly half of that, which is N; the [zzz] columns run
    #  to 2.192, which along that line is H.
    #
    #  And the check that comes free: at q = 2.192 all four branches that
    #  reach it - L and T of [00z], L and T of [zzz] - give 35.4 meV, which
    #  is H and triply degenerate.
    #  TWO ENTRIES ARE DROPPED, both plainly corrupt in print rather than
    #  merely surprising: T2[011] at q = 0.950 reads "0.123" where its
    #  neighbours are 1.16 and 1.73, and T[111] at q = 0.860 reads "10.8"
    #  where its neighbours are 1.08 and 1.36.  Each is its own value with a
    #  decimal point in the wrong place, and each would draw an order of
    #  magnitude away from its own branch.  They are left out rather than
    #  repaired, because repairing a printed number is guessing however
    #  obvious it looks.
    #
    #  The T1[011] branch runs past K, to 1.800 against K at 1.631.  Those
    #  points are kept in the table and dropped by the range, the same way
    #  the other fcc records handle Sigma beyond K.
    "Ag": {"T_K": 296, "struct": "fcc", "ref": AG_REF,
           "nu_scale": RADS13_THZ,
           "seg_override": {"[00z]": [("G", "X", 0.0, 1.5378)],
                            "[zz0]": [("G", "K", 0.0, 1.6311)],
                            "[zzz]": [("G", "L", 0.0, 1.3318)]},
           "branches": {
        "L[00z]": [[0.400, 1.39, 0.10], [0.800, 2.30, 0.15],
                   [1.119, 2.96, 0.15], [1.189, 2.85, 0.15],
                   [1.522, 3.19, 0.20]],
        "T[00z]": [[0.350, 0.73, 0.05], [0.500, 1.05, 0.05],
                   [0.770, 1.40, 0.10], [1.350, 2.01, 0.15],
                   [1.420, 2.10, 0.15]],
        "L[zz0]": [[0.300, 1.00, 0.05], [0.500, 1.70, 0.10],
                   [0.950, 2.49, 0.15], [1.700, 2.30, 0.15],
                   [1.910, 2.22, 0.15], [2.007, 2.14, 0.10],
                   [2.010, 2.23, 0.07], [2.102, 2.08, 0.10],
                   [2.139, 2.08, 0.10]],
        "T2[zz0]": [[0.700, 0.87, 0.05], [0.810, 1.06, 0.10],
                    [0.880, 1.16, 0.10], [1.450, 1.73, 0.15]],
        "T1[zz0]": [[0.280, 0.58, 0.05], [0.460, 0.96, 0.10],
                    [0.520, 1.00, 0.10], [0.753, 1.56, 0.08],
                    [0.840, 1.62, 0.15], [0.894, 1.71, 0.09],
                    [0.981, 1.85, 0.08], [1.040, 1.83, 0.15],
                    [1.081, 2.05, 0.07], [1.225, 2.21, 0.20],
                    [1.400, 2.47, 0.20], [1.480, 2.58, 0.20],
                    [1.520, 2.75, 0.15], [1.580, 2.77, 0.15],
                    [1.640, 2.82, 0.15], [1.750, 2.88, 0.15],
                    [1.800, 2.88, 0.15]],
        "L[zzz]": [[0.400, 1.44, 0.10], [0.530, 2.00, 0.15],
                   [0.600, 2.02, 0.15], [0.880, 2.78, 0.25],
                   [1.067, 2.78, 0.15], [1.086, 3.09, 0.15]],
        "T[zzz]": [[0.330, 0.49, 0.05], [0.400, 0.56, 0.05],
                   [0.530, 0.80, 0.05], [0.650, 0.92, 0.10],
                   [0.740, 1.08, 0.10], [1.117, 1.40, 0.05],
                   [1.180, 1.36, 0.15]],
    }},
    "Fe": {"T_K": 295, "struct": "bcc", "ref": FE_REF, "nu_scale": MEV_THZ,
           "seg_override": {"[z00]": [("G", "H", 0.0, 2.192)],
                            "[zz0]": [("G", "N", 0.0, 1.096)],
                            "[zzz]": [("G", "P", 0.0, 1.096),
                                      ("P", "H", 1.096, 2.192)]},
           "branches": {
        "L[z00]": [[0.4, 13.8, 0.3], [0.6, 19.8, 0.4], [0.8, 25.0, 0.4],
                   [1.0, 29.3, 0.4], [1.2, 32.2, 0.3], [1.4, 34.2, 0.4],
                   [1.6, 35.5, 0.7], [1.8, 35.8, 0.5], [2.0, 36.0, 0.6],
                   [2.192, 35.4, 0.2]],
        "T[z00]": [[0.2, 5.0, 0.2], [0.4, 10.2, 0.05], [0.6, 15.0, 0.1],
                   [0.8, 19.3, 0.2], [1.0, 23.3, 0.2], [1.096, 25.1, 0.2],
                   [1.315, 28.7, 0.2], [1.534, 31.4, 0.3], [1.753, 33.4, 0.5],
                   [1.973, 35.4, 0.5], [2.192, 35.4, 0.2]],
        "L[zz0]": [[0.141, 7.5, 0.3], [0.212, 11.5, 0.3], [0.230, 12.9, 0.2],
                   [0.248, 13.8, 0.2], [0.265, 14.8, 0.2], [0.283, 15.6, 0.3],
                   [0.301, 16.6, 0.2], [0.318, 17.4, 0.2], [0.336, 18.3, 0.2],
                   [0.354, 18.9, 0.3], [0.424, 22.3, 0.3], [0.495, 25.5, 0.3],
                   [0.566, 28.6, 0.2], [0.636, 31.3, 0.2], [0.707, 33.6, 0.2],
                   [0.743, 34.4, 0.3], [0.778, 35.5, 0.3], [0.813, 36.1, 0.3],
                   [0.884, 37.1, 0.3], [0.955, 38.2, 0.3], [1.025, 38.3, 0.4],
                   [1.096, 38.3, 0.5]],
        "T1[zz0]": [[1.096, 18.5, 0.3]],
        "T2[zz0]": [[0.219, 7.9, 0.2], [0.329, 11.9, 0.2], [0.438, 15.4, 0.2],
                    [0.548, 18.5, 0.2], [0.566, 19.0, 0.2], [0.636, 20.8, 0.2],
                    [0.707, 22.4, 0.2], [0.778, 23.9, 0.2], [0.849, 25.1, 0.2],
                    [0.919, 26.1, 0.2], [0.990, 26.5, 0.2], [1.061, 26.7, 0.2],
                    [1.096, 26.7, 0.2]],
        "L[zzz]": [[0.231, 15.5, 0.4], [0.346, 22.4, 0.5], [0.462, 28.0, 0.5],
                   [0.577, 32.7, 0.5], [0.693, 33.9, 0.4], [0.751, 34.7, 0.4],
                   [0.808, 34.6, 0.3], [0.924, 33.8, 0.2], [1.039, 31.4, 0.3],
                   [1.097, 29.8, 0.2], [1.155, 28.5, 0.1], [1.212, 27.3, 0.1],
                   [1.270, 25.9, 0.1], [1.328, 24.7, 0.2], [1.386, 23.8, 0.1],
                   [1.501, 23.7, 0.4], [1.617, 24.9, 0.5], [1.730, 27.7, 0.2],
                   [1.846, 30.4, 0.2], [1.961, 33.2, 0.2], [2.077, 34.9, 0.3],
                   [2.192, 35.4, 0.2]],
        "T[zzz]": [[0.231, 7.9, 0.4], [0.346, 12.0, 0.2], [0.462, 16.1, 0.2],
                   [0.577, 19.9, 0.3], [0.693, 23.4, 0.4], [0.808, 25.8, 0.2],
                   [0.924, 27.8, 0.2], [1.039, 29.4, 0.2], [1.155, 30.4, 0.2],
                   [1.270, 31.8, 0.2], [1.386, 32.9, 0.3], [1.501, 33.4, 0.3],
                   [1.617, 34.5, 0.3], [1.732, 34.4, 0.4], [1.848, 34.9, 0.6],
                   [1.963, 35.1, 0.4], [2.079, 35.6, 0.4], [2.192, 35.4, 0.2]],
    }},
    #  Two segments here that no other bcc element in this set reaches.
    #  [1/2 1/2 z] is the D line, N -> P, and [zz1] is the G line, H -> N:
    #  at zeta = 1/2 it arrives at (1/2, 1/2, 1), and subtracting the
    #  reciprocal-lattice vector (1, 0, 1) turns that into (-1/2, 1/2, 0),
    #  which is an N point.  Both go in as seg_override, because SEGMENT
    #  carries only the three lines through Gamma.
    #
    #  The paper's own footnote settles the longitudinal [zzz] branch: F is
    #  "a continuation of the Lambda branch", so Lambda_1 below zeta = 1/2 and
    #  F_1 above it are one branch through P - which is how this file already
    #  treats [zzz] for tantalum and caesium.
    #  This one is a MEASUREMENT, not a model: 5.225 is Kittel's own stamped
    #  5 K value for potassium, and 5 K and 9 K are the same volume to well
    #  inside anything here can resolve.  It was being scored at refdata's
    #  5.328, a room-temperature value, 5.8 % of volume too large.
    "K": {"T_K": 9, "a_meas": 5.225, "a_meas_from": "Kittel 5 K row", "struct": "bcc", "ref": K_REF,
          "seg_override": {"[hhz]": [("N", "P", 0.0, 0.5)],
                           "[zz1]": [("H", "N", 0.0, 0.5)]},
          "branches": {
        "L[z00]": [[0.15, 0.66, 0.04], [0.20, 0.89, 0.025], [0.25, 1.08, 0.05],
                   [0.30, 1.21, 0.03], [0.35, 1.35, 0.04], [0.40, 1.53, 0.03],
                   [0.45, 1.63, 0.03], [0.50, 1.74, 0.035], [0.55, 1.81, 0.04],
                   [0.60, 1.89, 0.03], [0.65, 1.96, 0.06], [0.70, 1.965, 0.04],
                   [0.75, 2.085, 0.03], [0.80, 2.11, 0.025], [0.85, 2.15, 0.03],
                   [0.90, 2.19, 0.025], [0.95, 2.21, 0.025], [1.00, 2.21, 0.02]],
        "T[z00]": [[0.15, 0.52, 0.03], [0.20, 0.68, 0.025], [0.25, 0.82, 0.03],
                   [0.30, 0.99, 0.02], [0.35, 1.15, 0.03], [0.40, 1.28, 0.025],
                   [0.45, 1.44, 0.03], [0.50, 1.57, 0.02], [0.55, 1.69, 0.02],
                   [0.60, 1.79, 0.025], [0.65, 1.89, 0.02], [0.70, 1.99, 0.02],
                   [0.75, 2.07, 0.03], [0.80, 2.11, 0.025], [0.85, 2.19, 0.04],
                   [0.90, 2.20, 0.03], [0.95, 2.21, 0.04]],
        "L[zz0]": [[0.10, 0.74, 0.02], [0.15, 1.08, 0.03], [0.20, 1.41, 0.02],
                   [0.25, 1.69, 0.05], [0.30, 1.94, 0.03], [0.35, 2.08, 0.05],
                   [0.40, 2.25, 0.035], [0.45, 2.33, 0.06], [0.50, 2.40, 0.04]],
        "T[zz0]": [[0.15, 0.67, 0.03], [0.20, 0.93, 0.03], [0.25, 1.11, 0.03],
                   [0.30, 1.23, 0.03], [0.35, 1.36, 0.03], [0.40, 1.44, 0.025],
                   [0.45, 1.49, 0.05], [0.50, 1.50, 0.025]],
        "L[zzz]": [[0.05, 0.53, 0.03], [0.10, 0.95, 0.02], [0.15, 1.34, 0.04],
                   [0.20, 1.695, 0.02], [0.25, 1.95, 0.04], [0.30, 2.10, 0.025],
                   [0.35, 2.19, 0.05], [0.40, 2.15, 0.025], [0.45, 2.03, 0.06],
                   [0.55, 1.55, 0.04], [0.58, 1.38, 0.05], [0.60, 1.24, 0.03],
                   [0.62, 1.17, 0.03], [0.65, 1.04, 0.025], [0.68, 1.005, 0.025],
                   [0.70, 1.02, 0.03], [0.72, 1.06, 0.03], [0.75, 1.20, 0.03],
                   [0.78, 1.36, 0.03], [0.80, 1.47, 0.02], [0.85, 1.77, 0.03],
                   [0.90, 1.99, 0.03], [0.95, 2.19, 0.04]],
        "T[zzz]": [[0.20, 0.77, 0.05], [0.25, 1.00, 0.06], [0.30, 1.20, 0.06],
                   [0.40, 1.53, 0.05], [0.45, 1.68, 0.04], [0.50, 1.785, 0.02],
                   [0.55, 1.91, 0.04], [0.60, 2.00, 0.04], [0.65, 2.05, 0.06],
                   [0.70, 2.09, 0.03], [0.75, 2.12, 0.04], [0.80, 2.16, 0.03],
                   [0.85, 2.19, 0.03], [0.90, 2.18, 0.025], [0.95, 2.22, 0.03]],
        "D1[hhz]": [[0.10, 1.52, 0.04], [0.20, 1.62, 0.04], [0.30, 1.69, 0.04],
                    [0.40, 1.72, 0.04]],
        "D3[hhz]": [[0.10, 2.35, 0.06], [0.20, 2.30, 0.035], [0.30, 2.20, 0.05],
                    [0.40, 2.01, 0.04]],
        "D4[hhz]": [[0.10, 0.70, 0.03], [0.20, 0.97, 0.03], [0.30, 1.27, 0.03],
                    [0.40, 1.55, 0.04]],
        "G1[zz1]": [[0.05, 2.22, 0.05], [0.10, 2.14, 0.04], [0.15, 2.01, 0.06],
                    [0.20, 1.78, 0.03], [0.25, 1.55, 0.04], [0.30, 1.31, 0.04],
                    [0.35, 1.07, 0.04], [0.40, 0.81, 0.03], [0.45, 0.61, 0.03],
                    [0.50, 0.53, 0.02]],
        "G3[zz1]": [[0.05, 2.18, 0.05], [0.10, 2.13, 0.04], [0.15, 2.10, 0.035],
                    [0.20, 2.02, 0.035], [0.25, 1.895, 0.03], [0.30, 1.78, 0.025],
                    [0.35, 1.70, 0.04], [0.40, 1.57, 0.03], [0.45, 1.53, 0.04]],
    }},
    #  Sodium reaches every one of the six bcc segments, as potassium does.
    #  Two conventions are worth stating because the table states them itself:
    #  its zeta columns are headed aq/2pi for [00z], (1/sqrt2)(aq/2pi) for
    #  [zz0] and (1/sqrt3)(aq/2pi) for [zzz], and all three of those ARE the
    #  reduced component this file uses - the headings only say how to get it
    #  from |q|, so nothing is rescaled.
    #
    #  The paper's own text supplies the checks: "at the point (0,0,1) the
    #  longitudinal and transverse branches are degenerate by symmetry" and
    #  the same "at the point (1/2,1/2,1/2)".  Both hold in the numbers - all
    #  four branches reach 3.58 at H, and L and T meet at 2.88 at P - and the
    #  three N frequencies 0.93, 2.56 and 3.82 appear independently in the
    #  [zz0], [1/2 1/2 z] and [zz1] blocks.
    #
    #  THE D LINE IS FOLDED.  [1/2 1/2 z] runs N (zeta = 0) through P (1/2) to
    #  N' (1), and N' is an N point, so the second half is the same segment
    #  travelled back.  Points above 1/2 are stored at 1 - zeta.  That is a
    #  symmetry operation and not a rescaling; drawing them unfolded would put
    #  two copies of one segment on top of each other, and dropping them would
    #  throw away thirteen of the twenty-one measurements on that line.
    #  The fold checks itself twice.  The Lambda branch is measured on both
    #  halves and lands on 2.92 +/- 0.07 and 2.80 +/- 0.08 at the same folded
    #  abscissa - the same number to within the quoted errors.  And all three
    #  folded branches meet at 2.88 at P, which is the L-T degeneracy the
    #  paper states for (1/2,1/2,1/2).
    #  Measured at 90 K, scored at refdata's ~293 K a0 until 2026-08-29.
    #  a(90 K) from the expansion integral, whose prediction for sodium was
    #  checked against Kittel's 5 K value: measured 1.529 % against predicted
    #  1.714 %, so this carries about 12 % slack in the direction of too
    #  large a correction.
    "Na": {"T_K": 90, "a_meas": 4.2337, "a_meas_from": "expansion integral", "struct": "bcc", "ref": NA_REF,
           "seg_override": {"[hhz]": [("N", "P", 0.0, 0.5)],
                            "[zz1]": [("H", "N", 0.0, 0.5)]},
           "branches": {
        "L[z00]": [[0.20, 1.43, 0.07], [0.30, 1.94, 0.06], [0.40, 2.44, 0.05],
                   [0.45, 2.68, 0.1], [0.50, 2.78, 0.06], [0.55, 2.91, 0.07],
                   [0.60, 3.01, 0.07], [0.65, 3.19, 0.07], [0.68, 3.14, 0.1],
                   [0.70, 3.24, 0.06], [0.72, 3.25, 0.08], [0.74, 3.31, 0.08],
                   [0.75, 3.36, 0.1], [0.76, 3.36, 0.07], [0.80, 3.44, 0.05],
                   [0.85, 3.53, 0.06], [0.90, 3.55, 0.05], [1.00, 3.58, 0.04]],
        "T[z00]": [[0.18, 0.97, 0.03], [0.20, 1.09, 0.04], [0.22, 1.18, 0.04],
                   [0.24, 1.29, 0.04], [0.26, 1.42, 0.04], [0.28, 1.52, 0.04],
                   [0.30, 1.64, 0.03], [0.32, 1.74, 0.05], [0.34, 1.83, 0.04],
                   [0.36, 1.94, 0.05], [0.38, 2.07, 0.05], [0.40, 2.17, 0.04],
                   [0.42, 2.25, 0.05], [0.44, 2.32, 0.04], [0.50, 2.59, 0.05],
                   [0.60, 2.96, 0.03], [0.70, 3.23, 0.04], [0.75, 3.35, 0.04],
                   [0.80, 3.45, 0.05], [0.90, 3.57, 0.06], [1.00, 3.58, 0.04]],
        "L[zz0]": [[0.10, 1.25, 0.04], [0.20, 2.32, 0.03], [0.25, 2.77, 0.05],
                   [0.30, 3.17, 0.05], [0.35, 3.46, 0.06], [0.40, 3.67, 0.05],
                   [0.45, 3.75, 0.09], [0.50, 3.82, 0.07]],
        "T1[zz0]": [[0.14, 0.43, 0.03], [0.21, 0.61, 0.03], [0.28, 0.76, 0.03],
                    [0.35, 0.87, 0.03], [0.42, 0.92, 0.04], [0.50, 0.93, 0.02]],
        "T2[zz0]": [[0.15, 1.16, 0.04], [0.20, 1.52, 0.04], [0.25, 1.81, 0.03],
                    [0.28, 1.97, 0.03], [0.30, 2.09, 0.03], [0.35, 2.27, 0.04],
                    [0.40, 2.47, 0.04], [0.45, 2.52, 0.06], [0.50, 2.56, 0.05]],
        "L[zzz]": [[0.10, 1.53, 0.05], [0.20, 2.72, 0.06], [0.25, 3.16, 0.06],
                   [0.30, 3.38, 0.06], [0.35, 3.44, 0.05], [0.40, 3.42, 0.06],
                   [0.45, 3.22, 0.06], [0.50, 2.88, 0.04], [0.55, 2.48, 0.05],
                   [0.60, 2.06, 0.04], [0.62, 1.90, 0.04], [0.64, 1.78, 0.04],
                   [0.65, 1.74, 0.03], [0.66, 1.71, 0.04], [0.68, 1.67, 0.04],
                   [0.70, 1.68, 0.03], [0.72, 1.78, 0.05], [0.725, 1.77, 0.05],
                   [0.74, 1.89, 0.05], [0.75, 1.94, 0.04], [0.76, 2.04, 0.04],
                   [0.78, 2.22, 0.05], [0.80, 2.43, 0.04], [0.82, 2.64, 0.07],
                   [0.84, 2.82, 0.09], [0.85, 2.87, 0.05], [0.90, 3.28, 0.06],
                   [1.00, 3.58, 0.04]],
        "T[zzz]": [[0.20, 1.28, 0.06], [0.30, 1.92, 0.06], [0.40, 2.47, 0.05],
                   [0.50, 2.88, 0.04], [0.60, 3.21, 0.06], [0.70, 3.42, 0.06],
                   [0.72, 3.44, 0.09], [0.74, 3.46, 0.09], [0.75, 3.46, 0.05],
                   [0.76, 3.44, 0.09], [0.78, 3.44, 0.09], [0.80, 3.48, 0.05],
                   [0.82, 3.48, 0.09], [0.84, 3.54, 0.09], [0.86, 3.52, 0.09],
                   [0.90, 3.56, 0.05], [1.00, 3.58, 0.04]],
        "T2[hhz]": [[0.0, 2.56, 0.05], [0.125, 2.62, 0.06], [0.25, 2.74, 0.06],
                    [0.375, 2.92, 0.07], [0.50, 2.88, 0.04],
                    [0.375, 2.80, 0.08], [0.25, 2.73, 0.07], [0.0, 2.56, 0.05]],
        #  The Pi branch as tabulated runs N -> P -> N', and its two ends are
        #  NOT the same mode: 3.82 is L at N and 0.93 is T1 there.  Folded
        #  back it is therefore two branches of the N-P segment, not one, and
        #  they are split here so each is scored against its own polarisation.
        "L[hhz]": [[0.0, 3.82, 0.07], [0.25, 3.57, 0.08], [0.50, 2.88, 0.04]],
        "T1[hhz]": [[0.0, 0.93, 0.02], [0.02, 0.93, 0.03], [0.04, 0.96, 0.04],
                    [0.06, 1.01, 0.03], [0.08, 1.07, 0.03], [0.10, 1.14, 0.03],
                    [0.12, 1.23, 0.04], [0.14, 1.28, 0.04], [0.25, 1.81, 0.05],
                    [0.375, 2.38, 0.06], [0.50, 2.88, 0.04]],
        "T1[zz1]": [[0.0, 3.58, 0.04], [0.125, 3.33, 0.07], [0.25, 2.50, 0.05],
                    [0.375, 1.49, 0.04], [0.50, 0.93, 0.02]],
        "T2[zz1]": [[0.0, 3.58, 0.04], [0.125, 3.48, 0.07], [0.25, 3.14, 0.06],
                    [0.375, 2.75, 0.06], [0.50, 2.56, 0.05]],
    }},
    #  Copley and Brockhouse tabulate four temperatures, 12, 85, 120 and
    #  205 K, and the 120 K column is taken here although 12 K is the colder
    #  and so the fairer comparison with a 0 K lattice sum.  The reason is
    #  coverage: along Gamma-N the 12 K column holds two points on T2 and a
    #  single one on T1, against nine and five at 120 K, and the [zz1] and
    #  [1/2 1/2 z] lines are measured almost entirely at 120 K.  It is also
    #  the column the authors fit their own Born-von Karman model to.
    #
    #  The cost is stated rather than hidden.  Where both columns have the
    #  same mode the softening from 12 K to 120 K is 4.7% at H, 2.7% at P,
    #  2.3% on L at N, 5.9% on T1 at N and 7.8% on T2 at N - so this record
    #  sits a few per cent low against a harmonic calculation, in one
    #  direction, by that much.  Rubidium melts at 312 K, so 120 K is 0.38 of
    #  the melting point.
    #
    #  The paper supplies its own degeneracy check and it is used: "the
    #  phonon frequency at symmetry point H is given for [00z]T and L,
    #  z = 1.0, [zzz]T and L, z = 1.0, and [zz1]A, p2, and p1, z = 0.0".
    #  All five of those read 1.32.  At P the L and T branches meet at 1.10,
    #  and the three frequencies at N - 0.32, 0.885 and 1.465 - appear
    #  independently in the [zz0], [1/2 1/2 z] and [zz1] blocks.
    #
    #  The D line is folded as for sodium, and here the paper makes the
    #  reason explicit: its pi branch starts at 1.465, which is L at N, and
    #  ends at 0.32, which is T1 there.  Those are different modes, so the
    #  two halves are two branches of the N-P segment and are split as such.
    "Rb": {"T_K": 120, "struct": "bcc", "ref": RB_REF,
           #  "a is the lattice constant, 5.63 A at 120 K", their own
           #  words about their own crystal.  III/13a says 5.69 A at
           #  120 K for the same element, one per cent away and the
           #  same size as this whole correction - the crystal that was
           #  in the beam wins.
           "a_meas": 5.63,
           "seg_override": {"[hhz]": [("N", "P", 0.0, 0.5)],
                            "[zz1]": [("H", "N", 0.0, 0.5)]},
           "branches": {
        "T[z00]": [[0.10, 0.215, 0.01], [0.20, 0.40, 0.01], [0.25, 0.47, 0.02],
                   [0.30, 0.57, 0.01], [0.40, 0.755, 0.02], [0.50, 0.895, 0.025],
                   [0.60, 1.065, 0.015], [0.80, 1.27, 0.015], [0.90, 1.30, 0.015],
                   [1.00, 1.32, 0.02]],
        "L[z00]": [[0.20, 0.50, 0.025], [0.25, 0.62, 0.02], [0.30, 0.72, 0.02],
                   [0.35, 0.82, 0.02], [0.40, 0.88, 0.025], [0.45, 0.98, 0.02],
                   [0.50, 1.075, 0.03], [0.60, 1.225, 0.05], [0.65, 1.23, 0.04],
                   [0.70, 1.23, 0.04], [0.75, 1.27, 0.04], [0.80, 1.275, 0.03],
                   [0.90, 1.305, 0.02], [1.00, 1.32, 0.02]],
        "T[zzz]": [[0.10, 0.235, 0.03], [0.15, 0.35, 0.03], [0.20, 0.48, 0.02],
                   [0.30, 0.71, 0.03], [0.40, 0.92, 0.03], [0.45, 1.02, 0.03],
                   [0.50, 1.10, 0.02], [0.55, 1.20, 0.04], [0.60, 1.26, 0.04],
                   [0.70, 1.325, 0.02], [0.80, 1.32, 0.03], [0.90, 1.33, 0.05],
                   [1.00, 1.32, 0.02]],
        "L[zzz]": [[0.10, 0.525, 0.03], [0.20, 1.01, 0.03], [0.30, 1.28, 0.03],
                   [0.35, 1.31, 0.025], [0.40, 1.305, 0.025], [0.45, 1.24, 0.04],
                   [0.50, 1.10, 0.02], [0.55, 0.935, 0.03], [0.60, 0.72, 0.03],
                   [0.65, 0.64, 0.03], [0.70, 0.60, 0.03], [0.75, 0.685, 0.025],
                   [0.80, 0.87, 0.03], [0.90, 1.20, 0.04], [0.95, 1.285, 0.03],
                   [1.00, 1.32, 0.02]],
        "T1[zz0]": [[0.10, 0.11, 0.025], [0.20, 0.20, 0.02], [0.30, 0.265, 0.02],
                    [0.40, 0.315, 0.02], [0.50, 0.32, 0.025]],
        "T2[zz0]": [[0.10, 0.285, 0.02], [0.15, 0.395, 0.02], [0.20, 0.525, 0.015],
                    [0.25, 0.66, 0.03], [0.30, 0.735, 0.02], [0.35, 0.79, 0.015],
                    [0.40, 0.85, 0.015], [0.45, 0.89, 0.015], [0.50, 0.885, 0.02]],
        "L[zz0]": [[0.10, 0.41, 0.02], [0.20, 0.82, 0.03], [0.30, 1.185, 0.025],
                   [0.40, 1.415, 0.02], [0.50, 1.465, 0.02]],
        "T2[hhz]": [[0.00, 0.885, 0.02], [0.10, 0.895, 0.025], [0.15, 0.92, 0.035],
                    [0.20, 0.945, 0.04], [0.25, 0.975, 0.03], [0.35, 1.035, 0.03],
                    [0.40, 1.08, 0.06], [0.50, 1.10, 0.02]],
        "L[hhz]": [[0.00, 1.465, 0.02], [0.20, 1.42, 0.02], [0.30, 1.34, 0.02],
                   [0.40, 1.25, 0.025], [0.50, 1.10, 0.02]],
        "T1[hhz]": [[0.00, 0.32, 0.025], [0.10, 0.405, 0.015], [0.20, 0.59, 0.015],
                    [0.30, 0.77, 0.015], [0.40, 0.965, 0.015], [0.50, 1.10, 0.02]],
        "T1[zz1]": [[0.00, 1.32, 0.02], [0.10, 1.27, 0.02], [0.20, 1.04, 0.03],
                    [0.30, 0.78, 0.03], [0.40, 0.49, 0.025], [0.45, 0.37, 0.03],
                    [0.50, 0.32, 0.025]],
        "T2[zz1]": [[0.00, 1.32, 0.02], [0.10, 1.29, 0.025], [0.20, 1.205, 0.02],
                    [0.30, 1.07, 0.02], [0.40, 0.94, 0.03], [0.50, 0.885, 0.02]],
        "L[zz1]": [[0.00, 1.32, 0.02], [0.10, 1.335, 0.02], [0.20, 1.36, 0.05],
                   [0.30, 1.42, 0.03], [0.40, 1.47, 0.03], [0.50, 1.465, 0.02]],
    }},
    #  Lithium is measured at 293 K and that is not a choice - below about
    #  80 K it leaves the bcc phase for a close-packed 9R structure, so there
    #  is no cold bcc dispersion to have.  293 K is 0.65 of the melting point,
    #  the warmest reference in this set after caesium, and the paper's own
    #  Table IV says what that costs: the same phonons measured at 110 K and
    #  at 293 K soften by 1.9% (36.3 -> 35.6), 3.6% (33.4 -> 32.2), 4.4%
    #  (9.0 -> 8.6), 4.7% (19.1 -> 18.2) and 5.3% (23.4 -> 22.15), rising to
    #  11% (5.7 -> 5.05) for the lowest-frequency mode they tracked.
    #
    #  ISOTOPE.  This is 7Li, mass 7.016, while the library computes with the
    #  natural mass 6.94.  Frequencies go as 1/sqrt(m), so these numbers are
    #  0.54% below what natural lithium would give.  That is left uncorrected
    #  because it is far inside the error bars and inside the thermal shift
    #  above, but it is a real offset and it is in one direction.
    #
    #  The three zeta columns are headed aq/2pi, (1/sqrt2)(aq/2pi) and
    #  (1/sqrt3)(aq/2pi), which all reduce to the reduced component used here,
    #  so nothing is rescaled but the meV -> THz conversion.
    #
    #  Checks the reconstruction did not use: all four branches reach 36.0 at
    #  H, and L crosses BELOW T exactly at P - 30.3 against 27.7 at zeta
    #  0.451, then 25.2 against 29.3 at 0.548.  Interpolated to P the two give
    #  27.7 and 28.5, equal to within one error bar, which is the L-T
    #  degeneracy there.  That crossing is also a fifth independent sighting
    #  of the branch-order inversion past P that curve_mae.pick handles.
    "Li": {"T_K": 293, "struct": "bcc", "ref": LI_REF, "nu_scale": MEV_THZ,
           #  Table III of the same paper: 3.4900 at 110 K, 3.5105 at 293 K.
           #  Our a0 is the 5 K value, so the comparison was being made
           #  half a per cent too dense.
           "a_meas": 3.5105,
           "branches": {
        "L[z00]": [[0.112, 7.6, 0.3], [0.168, 10.9, 0.5], [0.196, 12.3, 0.5],
                   [0.223, 14.2, 0.5], [0.279, 17.0, 0.5], [0.335, 19.9, 0.5],
                   [0.363, 20.5, 0.5], [0.447, 25.1, 1.0], [0.503, 27.5, 1.0],
                   [0.559, 27.1, 1.0], [0.614, 28.0, 0.7], [0.670, 30.0, 0.8],
                   [0.782, 32.2, 0.8], [0.950, 35.6, 0.6], [1.00, 36.0, 1.0]],
        "T[z00]": [[0.112, 5.1, 0.3], [0.168, 7.9, 0.4], [0.223, 10.6, 0.4],
                   [0.279, 13.2, 0.5], [0.335, 16.5, 0.6], [0.391, 20.0, 0.8],
                   [0.447, 23.1, 1.0], [0.559, 27.5, 0.6], [0.614, 29.7, 0.6],
                   [0.782, 33.3, 1.0], [1.00, 36.0, 1.0]],
        "L[zz0]": [[0.04, 4.5, 0.3], [0.08, 8.6, 0.3], [0.119, 12.8, 0.6],
                   [0.158, 16.8, 0.6], [0.237, 24.8, 0.6], [0.316, 30.8, 0.6],
                   [0.356, 33.0, 0.5], [0.375, 34.0, 0.5], [0.395, 34.5, 0.5],
                   [0.435, 35.8, 0.6], [0.50, 36.5, 0.5]],
        "T1[zz0]": [[0.08, 2.2, 0.3], [0.138, 3.7, 0.3], [0.198, 4.9, 0.4],
                    [0.257, 6.0, 0.5], [0.316, 7.2, 0.5], [0.356, 7.3, 0.6],
                    [0.407, 8.0, 0.6], [0.50, 8.6, 0.6]],
        "T2[zz0]": [[0.08, 5.0, 0.3], [0.119, 7.5, 0.4], [0.158, 10.3, 0.5],
                    [0.198, 13.0, 0.5], [0.237, 14.8, 0.5], [0.316, 18.1, 0.5],
                    [0.356, 19.6, 0.5], [0.435, 21.5, 0.5], [0.50, 22.2, 0.5]],
        "L[zzz]": [[0.097, 13.0, 1.0], [0.161, 20.8, 1.0], [0.226, 27.5, 0.6],
                   [0.290, 31.8, 1.0], [0.355, 33.3, 1.2], [0.451, 30.3, 0.8],
                   [0.548, 25.2, 0.8], [0.645, 16.5, 0.6], [0.677, 15.0, 0.5],
                   [0.710, 14.2, 1.0], [0.742, 17.0, 0.8], [0.839, 27.3, 1.0],
                   [1.00, 36.0, 1.0]],
        "T[zzz]": [[0.097, 6.4, 0.5], [0.129, 8.4, 0.5], [0.161, 9.0, 0.5],
                   [0.194, 11.7, 0.6], [0.226, 12.7, 0.6], [0.258, 15.2, 0.6],
                   [0.290, 18.3, 0.8], [0.323, 20.0, 1.0], [0.355, 23.5, 1.2],
                   [0.419, 26.0, 1.0], [0.451, 27.7, 0.8], [0.548, 29.3, 0.8],
                   [0.710, 33.5, 0.6], [0.774, 34.0, 0.8], [1.00, 36.0, 1.0]],
    }},
    #  The paper this comes from tabulates nothing.  Stedman and Nilsson
    #  present aluminium as dispersion curves and contour maps and say so:
    #  "errors have not been indicated for points on dispersion curves, but
    #  the size of the points corresponds to an error of 0.03".  Their NUMBERS
    #  reached print in Landolt-Boernstein III/13a, Table 2 Al, headed
    #  "unpublished data of [66St1]" - so the measurement is theirs and the
    #  tabulation is Schober and Dederichs'.  Both are credited in AL_REF.
    #
    #  80 K is taken over the 300 K column that sits beside it.  III/13a puts
    #  the softening between them at 7.5% on average, and a further 15% up to
    #  930 K, so the choice is worth about that much and it is the colder one.
    #
    #  T1 AND T2 ARE SWAPPED against the source, deliberately.  III/13a calls
    #  the HIGHER transverse branch T1 - at zeta 0.5 its T1 is 6.37 and its T2
    #  4.65 - while every other record here, and pick_polar in curve_mae, take
    #  T1 as the lower of the pair.  The data below is relabelled to this
    #  file's convention rather than the convention silently disagreeing.
    #
    #  Sigma runs to zeta = 1, which is past K at 0.75, so the tail folds onto
    #  U-X.  The check that it is folded right is at the end point: [0zz]L
    #  reads 5.78 there while [00z]L reads 9.69 at the same X.  Those are not
    #  in conflict - they are longitudinal with respect to DIFFERENT arrival
    #  directions, [011] against [001], and at X the mode polarised along
    #  [001] has no component along [011].  curve_mae identifies the branch
    #  from its eigenvector and lands on 5.78, which is the source's own L.
    "Al": {"T_K": 80, "struct": "fcc", "ref": AL_REF,
           #  The paper gives no lattice constant but does give the unit it
           #  measures q in: "2 pi / (the side of the unit cube), i.e.
           #  1.561 A^-1 at 80 K, 1.555 A^-1 at 300 K".  The RATIO of those
           #  is what is used - 0.99616 - applied to our own 293 K a0, since
           #  their absolute values are rounded to four figures and land
           #  0.2% off the accepted room-temperature value.
           "a_meas": 4.0495 * (1.555 / 1.561),
           "seg_override": {"[0zz]": [("G", "K", 0.0, 0.75),
                                      ("U", "X", 0.75, 1.0)]},
           "branches": {
        "L[00z]": [[0.175, 2.85, 0.03], [0.200, 3.20, 0.03], [0.225, 3.60, 0.03],
                   [0.250, 4.04, 0.03], [0.275, 4.38, 0.03], [0.300, 4.74, 0.03],
                   [0.401, 5.98, 0.03], [0.448, 6.51, 0.03], [0.504, 7.08, 0.03],
                   [0.600, 7.93, 0.03], [0.700, 8.53, 0.03], [0.750, 8.85, 0.08],
                   [0.800, 9.04, 0.05], [0.850, 9.37, 0.10], [0.900, 9.53, 0.05],
                   [0.950, 9.68, 0.06], [1.000, 9.69, 0.05]],
        "T[00z]": [[0.202, 1.75, 0.02], [0.301, 2.58, 0.02], [0.398, 3.29, 0.02],
                   [0.497, 4.06, 0.02], [0.599, 4.62, 0.02], [0.699, 5.05, 0.02],
                   [0.797, 5.43, 0.02], [0.898, 5.76, 0.03], [0.988, 5.79, 0.03]],
        "L[zzz]": [[0.101, 2.96, 0.03], [0.151, 4.36, 0.03], [0.200, 5.55, 0.03],
                   [0.252, 6.72, 0.03], [0.275, 7.26, 0.03], [0.301, 7.80, 0.03],
                   [0.325, 8.28, 0.03], [0.350, 8.66, 0.03], [0.400, 9.20, 0.03],
                   [0.449, 9.57, 0.03], [0.498, 9.69, 0.10]],
        "T[zzz]": [[0.151, 2.15, 0.02], [0.204, 2.83, 0.03], [0.250, 3.37, 0.02],
                   [0.300, 3.76, 0.02], [0.353, 4.04, 0.02], [0.400, 4.20, 0.03],
                   [0.425, 4.19, 0.03], [0.450, 4.19, 0.03], [0.500, 4.19, 0.03]],
        "L[0zz]": [[0.071, 1.62, 0.03], [0.108, 2.47, 0.03], [0.141, 3.25, 0.03],
                   [0.175, 4.03, 0.03], [0.216, 4.84, 0.03], [0.251, 5.47, 0.03],
                   [0.283, 6.02, 0.03], [0.322, 6.56, 0.03], [0.354, 6.97, 0.03],
                   [0.386, 7.35, 0.03], [0.424, 7.72, 0.03], [0.457, 8.09, 0.03],
                   [0.494, 8.39, 0.03], [0.530, 8.56, 0.03], [0.562, 8.67, 0.08],
                   [0.601, 8.63, 0.05], [0.634, 8.53, 0.05], [0.672, 8.20, 0.05],
                   [0.706, 7.94, 0.03], [0.745, 7.64, 0.03], [0.780, 7.27, 0.03],
                   [0.815, 6.96, 0.03], [0.861, 6.54, 0.03], [0.892, 6.29, 0.03],
                   [0.928, 6.02, 0.03], [0.955, 5.92, 0.03], [0.998, 5.78, 0.03]],
        #  the source's T2 - the LOWER branch
        "T1[0zz]": [[0.212, 2.42, 0.03], [0.283, 3.18, 0.02], [0.354, 3.74, 0.02],
                    [0.423, 4.20, 0.02], [0.496, 4.65, 0.02], [0.566, 4.98, 0.02],
                    [0.633, 5.28, 0.02], [0.707, 5.54, 0.02], [0.779, 5.71, 0.03],
                    [0.847, 5.83, 0.03], [0.915, 5.84, 0.03], [0.991, 5.76, 0.03]],
        #  the source's T1 - the HIGHER branch
        "T2[0zz]": [[0.144, 1.75, 0.02], [0.180, 2.23, 0.02], [0.215, 2.70, 0.02],
                    [0.251, 3.23, 0.02], [0.284, 3.64, 0.02], [0.316, 4.07, 0.02],
                    [0.354, 4.55, 0.02], [0.386, 4.95, 0.02], [0.424, 5.46, 0.02],
                    [0.457, 5.87, 0.03], [0.496, 6.37, 0.03], [0.528, 6.72, 0.03],
                    [0.564, 7.10, 0.03], [0.598, 7.42, 0.05], [0.632, 7.77, 0.03],
                    [0.672, 8.05, 0.03], [0.710, 8.32, 0.03], [0.745, 8.61, 0.03],
                    [0.780, 8.88, 0.03], [0.823, 9.17, 0.05], [0.853, 9.28, 0.08],
                    [0.887, 9.47, 0.03], [0.921, 9.57, 0.03], [0.958, 9.66, 0.03],
                    [1.000, 9.69, 0.05]],
    }},
    "Lu": {"T_K": 295, "struct": "hcp", "ref": LU_REF,
           "seg_override": {"[zz0]": [("G", "K", 0.0, 2.0 / 3.0),
                                      ("K", "M", 2.0 / 3.0, 1.0)],
                            "[00z]": [("G", "A", 0.0, 0.5)],
                            "[z00]": [("G", "M", 0.0, 1.0)]},
           "branches": {
        "LA[zz0]": [[0.08, 0.70, 0.02], [0.1, 0.86, 0.02], [0.2, 1.60, 0.05],
                    [0.3, 2.18, 0.01], [0.4, 2.45, 0.01], [0.5, 2.69, 0.01],
                    [0.7, 2.77, 0.01], [0.8, 2.83, 0.01], [0.9, 2.98, 0.01],
                    [1.0, 3.06, 0.01]],
        "TOpar[zz0]": [[0.0, 2.04, 0.01], [0.2, 2.34, 0.01], [0.3, 2.63, 0.01],
                       [0.35, 2.78, 0.02], [0.4, 2.88, 0.05], [0.5, 3.14, 0.01],
                       [0.6, 3.20, 0.02], [0.8, 3.27, 0.01], [1.0, 3.23, 0.01]],
        "TOperp[zz0]": [[0.0, 3.63, 0.02], [0.1, 3.54, 0.01], [0.2, 3.52, 0.01],
                        [0.4, 3.26, 0.02], [0.6, 2.92, 0.01], [0.7, 2.55, 0.01],
                        [0.9, 2.09, 0.01], [1.0, 1.99, 0.01]],
        "TAperp[zz0]": [[0.2, 0.96, 0.01], [0.3, 1.40, 0.01], [0.4, 1.82, 0.02],
                        [0.5, 2.21, 0.02], [0.6, 2.52, 0.02], [0.7, 2.84, 0.03],
                        [0.8, 3.07, 0.03], [0.9, 3.26, 0.05], [1.0, 3.26, 0.05]],
        "LO[zz0]": [[0.0, 2.04, 0.01], [0.1, 2.16, 0.02], [0.2, 2.46, 0.02],
                    [0.3, 2.82, 0.02], [0.4, 3.10, 0.01], [0.6, 3.08, 0.01],
                    [0.7, 3.08, 0.01], [0.8, 3.18, 0.01], [1.0, 3.39, 0.01]],
        "TApar[zz0]": [[0.15, 0.76, 0.01], [0.2, 1.03, 0.01], [0.3, 1.50, 0.01],
                       [0.4, 1.94, 0.01], [0.45, 2.15, 0.01], [0.5, 2.34, 0.02],
                       [0.7, 2.71, 0.01], [0.8, 2.36, 0.01], [0.9, 1.95, 0.01],
                       [1.0, 1.84, 0.01]],
        "LA[00z]": [[0.1, 0.55, 0.01], [0.2, 1.05, 0.01], [0.3, 1.54, 0.01],
                    [0.4, 2.00, 0.01], [0.5, 2.41, 0.01]],
        "LO[00z]": [[0.0, 3.63, 0.02], [0.1, 3.57, 0.01], [0.2, 3.45, 0.01],
                    [0.3, 3.17, 0.01], [0.4, 2.83, 0.01], [0.5, 2.41, 0.01]],
        "TA[00z]": [[0.4, 1.20, 0.01], [0.5, 1.47, 0.01]],
        "TO[00z]": [[0.0, 2.04, 0.01], [0.1, 2.02, 0.02], [0.2, 1.95, 0.01],
                    [0.3, 1.84, 0.01], [0.4, 1.68, 0.01], [0.5, 1.47, 0.01]],
        "LA[z00]": [[0.2, 1.14, 0.01], [0.3, 1.60, 0.01], [0.4, 1.96, 0.01],
                    [0.5, 2.31, 0.01], [0.6, 2.64, 0.01], [0.7, 2.95, 0.01],
                    [0.8, 3.16, 0.01], [1.0, 3.23, 0.01]],
        "LO[z00]": [[0.0, 2.04, 0.01], [0.2, 2.28, 0.01], [0.4, 2.72, 0.01],
                    [0.6, 3.03, 0.01], [0.9, 3.33, 0.01], [1.0, 3.39, 0.01]],
        "TAperp[z00]": [[0.2, 0.53, 0.01], [0.3, 0.81, 0.01], [0.4, 1.12, 0.01],
                        [0.5, 1.32, 0.01], [0.6, 1.56, 0.01], [0.7, 1.76, 0.01],
                        [0.8, 1.95, 0.01], [0.9, 2.10, 0.01], [1.0, 1.99, 0.01]],
        "TOperp[z00]": [[0.0, 3.63, 0.02], [0.1, 3.67, 0.01], [0.2, 3.63, 0.02],
                        [0.3, 3.59, 0.01], [0.4, 3.57, 0.02], [0.5, 3.55, 0.01],
                        [0.6, 3.45, 0.01], [0.7, 3.37, 0.01], [0.8, 3.31, 0.01],
                        [0.9, 3.18, 0.01], [1.0, 3.26, 0.05]],
        "TApar[z00]": [[0.2, 0.55, 0.01], [0.3, 0.85, 0.01], [0.4, 1.09, 0.01],
                       [0.5, 1.31, 0.01], [0.6, 1.51, 0.01], [0.7, 1.67, 0.01],
                       [0.8, 1.79, 0.02], [0.9, 1.83, 0.01], [1.0, 1.84, 0.01]],
        "TOpar[z00]": [[0.0, 2.04, 0.01], [0.2, 2.18, 0.01], [0.4, 2.41, 0.01],
                       [0.6, 2.72, 0.01], [0.8, 2.94, 0.01], [1.0, 3.06, 0.01]],
    }},
    "Ta": {"T_K": 296, "struct": "bcc", "ref": TA_REF, "branches": {
        "L[z00]": [[0.10, 1.23, 0.06], [0.20, 2.28, 0.05], [0.30, 3.10, 0.05],
                   [0.40, 3.82, 0.06], [0.50, 4.20, 0.08], [0.60, 4.25, 0.08],
                   [0.70, 4.45, 0.08], [0.80, 4.68, 0.08], [0.90, 4.98, 0.10],
                   [1.00, 5.03, 0.07]],
        "T[z00]": [[0.20, 1.28, 0.04], [0.30, 1.88, 0.04], [0.40, 2.60, 0.04],
                   [0.50, 3.37, 0.05], [0.60, 4.03, 0.07], [0.70, 4.63, 0.08],
                   [0.80, 4.85, 0.08], [0.90, 5.03, 0.08], [1.00, 5.03, 0.07]],
        "L[zzz]": [[0.10, 2.24, 0.05], [0.20, 3.75, 0.07], [0.30, 4.38, 0.10],
                   [0.40, 4.30, 0.08], [0.50, 3.78, 0.06], [0.60, 3.03, 0.06],
                   [0.70, 2.70, 0.06], [0.80, 3.70, 0.10], [0.90, 4.78, 0.10],
                   [1.00, 5.03, 0.07]],
        "T[zzz]": [[0.10, 0.98, 0.04], [0.20, 1.73, 0.04], [0.30, 2.48, 0.04],
                   [0.40, 3.28, 0.05], [0.50, 3.78, 0.06], [0.60, 3.80, 0.10],
                   [0.70, 4.20, 0.15], [0.80, 4.48, 0.15], [0.90, 4.90, 0.10],
                   [1.00, 5.03, 0.07]],
        "L[zz0]": [[0.10, 1.65, 0.10], [0.20, 3.10, 0.10], [0.25, 3.51, 0.07],
                   [0.30, 3.95, 0.08], [0.35, 4.15, 0.08], [0.40, 4.24, 0.08],
                   [0.41, 4.28, 0.08], [0.45, 4.32, 0.08], [0.50, 4.35, 0.08]],
        "T1[zz0]": [[0.50, 2.63, 0.08]],
        "T2[zz0]": [[0.10, 0.96, 0.04], [0.20, 1.95, 0.04], [0.30, 2.97, 0.05],
                    [0.35, 3.51, 0.05], [0.40, 3.97, 0.06], [0.45, 4.24, 0.06],
                    [0.50, 4.35, 0.06]],
    }},
    "Pb": {"T_K": 100, "ref": PB_REF,
           #  "a is the cubic lattice constant (4.924 A for Pb at 100 K)"
           "a_meas": 4.9240,
           "seg_override": {"[zzz]": [("G", "L", 0.0, SQ3 / 2)],
                            "[zz0]": [("G", "K", 0.0, 0.75 * SQ2)]},
           "branches": {
        "T[z00]": [[0.20, 0.47, 0.04], [0.30, 0.73, 0.03], [0.40, 0.90, 0.02],
                   [0.45, 0.96, 0.03], [0.50, 1.04, 0.02], [0.60, 1.115, 0.02],
                   [0.70, 1.10, 0.03], [0.80, 1.03, 0.03], [0.90, 0.95, 0.02],
                   [1.00, 0.89, 0.02]],
        "L[z00]": [[0.20, 0.87, 0.03], [0.30, 1.27, 0.03], [0.40, 1.61, 0.03],
                   [0.45, 1.71, 0.04], [0.50, 1.83, 0.04], [0.55, 1.91, 0.03],
                   [0.60, 2.00, 0.03], [0.65, 2.07, 0.03], [0.70, 2.14, 0.04],
                   [0.75, 2.16, 0.02], [0.80, 2.15, 0.02], [0.85, 2.14, 0.02],
                   [0.90, 2.05, 0.03], [0.95, 1.94, 0.04], [1.00, 1.86, 0.03]],
        "T[zzz]": [[0.19, 0.35, 0.04], [0.26, 0.44, 0.03], [0.35, 0.55, 0.03],
                   [0.40, 0.61, 0.02], [0.45, 0.66, 0.02], [0.50, 0.73, 0.02],
                   [0.55, 0.77, 0.02], [0.60, 0.79, 0.02], [0.65, 0.805, 0.02],
                   [0.70, 0.835, 0.02], [0.75, 0.86, 0.02], [0.80, 0.88, 0.02],
                   [0.867, 0.89, 0.02]],
        "L[zzz]": [[0.143, 0.82, 0.05], [0.25, 1.24, 0.04], [0.332, 1.51, 0.04],
                   [0.433, 1.76, 0.03], [0.519, 1.91, 0.03], [0.563, 1.97, 0.03],
                   [0.606, 2.00, 0.03], [0.649, 2.05, 0.03], [0.693, 2.085, 0.02],
                   [0.736, 2.09, 0.02], [0.779, 2.08, 0.03], [0.823, 2.16, 0.03],
                   [0.867, 2.185, 0.02]],
        "L[zz0]": [[0.212, 0.99, 0.04], [0.284, 1.26, 0.04], [0.325, 1.40, 0.04],
                   [0.35, 1.48, 0.03], [0.375, 1.545, 0.02], [0.40, 1.605, 0.02],
                   [0.425, 1.64, 0.02], [0.45, 1.665, 0.02], [0.475, 1.675, 0.03],
                   [0.50, 1.715, 0.02], [0.525, 1.76, 0.02], [0.55, 1.795, 0.02],
                   [0.60, 1.85, 0.02], [0.65, 1.92, 0.02], [0.70, 2.01, 0.02],
                   [0.778, 2.10, 0.02], [0.85, 2.09, 0.02], [0.90, 2.04, 0.03],
                   [0.99, 1.925, 0.04], [1.061, 1.75, 0.04]],
        "T2[zz0]": [[0.20, 0.53, 0.03], [0.30, 0.75, 0.03], [0.40, 0.95, 0.03],
                    [0.50, 1.17, 0.03], [0.60, 1.37, 0.03], [0.70, 1.57, 0.03],
                    [0.814, 1.78, 0.03], [0.914, 1.92, 0.04], [1.00, 2.02, 0.03]],
        "T1[zz0]": [[0.141, 0.21, 0.03], [0.283, 0.41, 0.03], [0.424, 0.56, 0.03],
                    [0.566, 0.76, 0.03], [0.707, 0.91, 0.03], [0.848, 1.12, 0.03],
                    [0.959, 1.20, 0.03], [1.056, 1.25, 0.03]],
    }},
}


#  --- Ru and Ir, DIGITISED from figures (2026-09-13) -----------------------
#  Neither source tabulates a single frequency; both publish the measurement
#  only as dots in a figure.  The values below are READ from those figures and
#  are kept in meV exactly as read, converted by nu_scale on the way out like
#  the other energy-printed records.  The zeta of each custom bracket key runs
#  from 0 at the first named point to 1 at the second, as the figure draws it.
#  They carry no branch labels, so both are scored by nearest branch.
RU_REF = ("R. Heid, L. Pintschovius, W. Reichardt and K.-P. Bohnen, "
          "Phys. Rev. B <b>61</b>, 12059 (2000), Fig. 1: inelastic neutron "
          "scattering on a single crystal at room temperature.  <b>Digitised from "
          "the figure</b> - the paper tabulates no frequencies - 283 points on the "
          "drawn path; the three low-symmetry panels lie off it.  Reading "
          "uncertainty 0.06 meV (median; 90%% under 0.14), measured by reading the "
          "same &Gamma;&ndash;M data again from its Fig. 2.  The paper puts its error "
          "bars at about the symbol size, carried here as 0.37 meV")
IR_REF = ("R. Heid, K.-P. Bohnen, K. Felix, K. M. Ho and W. Reichardt, "
          "J. Phys.: Condens. Matter <b>10</b>, 7967 (1998), Fig. 1, filled circles: "
          "first results of an inelastic neutron scattering study at Saclay, room "
          "temperature (cited there as Reichardt, to be published).  The open "
          "symbols in that figure are theory and are not used.  <b>Digitised from "
          "the figure</b>, which is vector in the PDF, so the coordinates are "
          "exact; no error bars are given.  85 distinct points; 4 along L&ndash;X are off "
          "the drawn path and not used.  Checked against S. Wang et al., Phys. "
          "Rev. B <b>112</b>, 184305 (2025), Fig. 9, 300 K: &Gamma;&ndash;X agrees "
          "to 0.35 meV (median); on &Gamma;&ndash;L the transverse branch of the 2025 "
          "data lies up to 1.8 meV higher - a difference between the two "
          "experiments, not a reading error")
PHONON_CURVE["Ru"] = {
    "T_K": 295, "struct": "hcp", "ref": RU_REF, "nu_scale": MEV_THZ,
    "nearest": True,
    "digitised": ("These points were read from Fig. 1 of the source, which "
                  "tabulates no frequencies.  Reading uncertainty 0.06&nbsp;meV, "
                  "measured by reading the same &Gamma;&ndash;M data from its Fig. 2."),
    "seg_override": {
        "[AH]": [("A", "H", 0.0, 1.0)],
        "[GA]": [("G", "A", 0.0, 1.0)],
        "[GK]": [("G", "K", 0.0, 1.0)],
        "[HL]": [("H", "L", 0.0, 1.0)],
        "[KH]": [("K", "H", 0.0, 1.0)],
        "[KM]": [("K", "M", 0.0, 1.0)],
        "[LA]": [("L", "A", 0.0, 1.0)],
        "[MG]": [("M", "G", 0.0, 1.0)],
        "[ML]": [("M", "L", 0.0, 1.0)],
    },
    "branches": {
        "D[AH]": [[0.1455, 17.63, 0.37], [0.1467, 26.85, 0.37], [0.2973, 20.34, 0.37], [0.2988, 19.64, 0.37], [0.3004, 28.72, 0.37], [0.4475, 20.16, 0.37], [0.4475, 31.1, 0.37], [0.4481, 23.66, 0.37], [0.5959, 26.18, 0.37], [0.5982, 31.92, 0.37], [0.5994, 19.93, 0.37], [0.7489, 27.85, 0.37], [0.7498, 20.23, 0.37], [0.7511, 31.18, 0.37], [0.8989, 21.18, 0.37], [0.9002, 30.19, 0.37], [0.9012, 29.19, 0.37], [1.0, 21.65, 0.37], [1.0, 29.85, 0.37]],
        "D[GA]": [[0.0978, 1.83, 0.37], [0.1034, 3.44, 0.37], [0.1034, 32.65, 0.37], [0.1952, 3.69, 0.37], [0.196, 32.2, 0.37], [0.1969, 6.84, 0.37], [0.2, 23.47, 0.37], [0.2976, 31.64, 0.37], [0.3938, 30.61, 0.37], [0.4011, 7.17, 0.37], [0.4013, 13.28, 0.37], [0.4023, 22.68, 0.37], [0.4901, 29.7, 0.37], [0.4942, 16.26, 0.37], [0.5871, 28.71, 0.37], [0.5932, 10.34, 0.37], [0.5944, 19.0, 0.37], [0.5966, 21.23, 0.37], [0.6971, 21.43, 0.37], [0.6987, 27.77, 0.37], [0.7885, 27.25, 0.37], [0.7917, 13.59, 0.37], [0.7922, 23.48, 0.37], [0.7938, 19.19, 0.37], [0.8952, 25.27, 0.37], [0.8952, 26.84, 0.37]],
        "D[GK]": [[0.0013, 32.77, 0.37], [0.0021, 23.69, 0.37], [0.075, 5.13, 0.37], [0.076, 24.11, 0.37], [0.0772, 2.97, 0.37], [0.1493, 24.18, 0.37], [0.1503, 32.33, 0.37], [0.151, 10.33, 0.37], [0.151, 25.01, 0.37], [0.1513, 5.92, 0.37], [0.2255, 26.12, 0.37], [0.2257, 24.98, 0.37], [0.2262, 15.17, 0.37], [0.2984, 25.75, 0.37], [0.2991, 31.19, 0.37], [0.2992, 11.1, 0.37], [0.2998, 27.25, 0.37], [0.3, 12.01, 0.37], [0.3001, 19.15, 0.37], [0.3729, 30.19, 0.37], [0.3736, 28.19, 0.37], [0.3745, 22.16, 0.37], [0.3745, 26.26, 0.37], [0.4465, 15.67, 0.37], [0.4465, 24.13, 0.37], [0.4471, 17.14, 0.37], [0.4477, 27.34, 0.37], [0.4478, 28.93, 0.37], [0.5277, 25.38, 0.37], [0.5278, 17.47, 0.37], [0.5278, 28.7, 0.37], [0.5285, 19.22, 0.37], [0.6004, 29.68, 0.37], [0.601, 25.63, 0.37], [0.6012, 26.32, 0.37], [0.6014, 18.88, 0.37], [0.6014, 27.98, 0.37], [0.602, 21.16, 0.37], [0.6758, 27.42, 0.37], [0.6758, 30.41, 0.37], [0.6763, 22.54, 0.37], [0.6765, 19.92, 0.37], [0.6765, 25.66, 0.37], [0.7487, 31.09, 0.37], [0.7491, 26.77, 0.37], [0.7495, 24.31, 0.37], [0.7498, 23.56, 0.37], [0.7499, 20.6, 0.37], [0.7501, 25.37, 0.37], [0.824, 21.18, 0.37], [0.824, 26.5, 0.37], [0.8246, 24.72, 0.37], [0.8246, 31.51, 0.37], [0.8981, 22.91, 0.37], [0.8981, 25.21, 0.37], [0.8996, 31.95, 0.37], [0.8997, 21.86, 0.37], [0.8998, 26.3, 0.37], [0.9985, 26.19, 0.37], [0.9994, 22.43, 0.37], [0.9994, 32.24, 0.37]],
        "D[HL]": [[0.0986, 29.79, 0.37], [0.0992, 21.41, 0.37], [0.2434, 29.8, 0.37], [0.2466, 20.93, 0.37], [0.3926, 20.47, 0.37], [0.3937, 30.43, 0.37], [0.3939, 27.19, 0.37], [0.7009, 20.15, 0.37], [0.7032, 31.39, 0.37], [0.7045, 24.03, 0.37], [0.9989, 20.01, 0.37], [1.0, 22.53, 0.37], [1.0, 31.34, 0.37], [1.0, 31.85, 0.37]],
        "D[KH]": [[0.0, 22.38, 0.37], [0.0, 26.17, 0.37], [0.0005, 32.23, 0.37], [0.2072, 22.83, 0.37], [0.2076, 26.3, 0.37], [0.3319, 31.91, 0.37], [0.3351, 23.0, 0.37], [0.3353, 27.32, 0.37], [0.3369, 26.73, 0.37], [0.4992, 22.56, 0.37], [0.5005, 27.81, 0.37], [0.5039, 31.42, 0.37], [0.505, 28.87, 0.37], [0.6646, 21.99, 0.37], [0.6655, 28.71, 0.37], [0.6667, 30.81, 0.37], [0.6683, 29.62, 0.37], [0.801, 29.87, 0.37], [0.9784, 29.84, 0.37], [0.9815, 21.72, 0.37]],
        "D[KM]": [[0.0996, 22.77, 0.37], [0.1015, 25.42, 0.37], [0.1021, 26.15, 0.37], [0.1027, 32.25, 0.37], [0.1042, 22.15, 0.37], [0.25, 23.7, 0.37], [0.25, 31.64, 0.37], [0.2511, 25.8, 0.37], [0.3926, 21.98, 0.37], [0.3955, 25.56, 0.37], [0.3966, 30.64, 0.37], [0.3967, 22.91, 0.37], [0.3992, 21.51, 0.37], [0.5477, 25.81, 0.37], [0.5487, 19.71, 0.37], [0.55, 29.56, 0.37], [0.6953, 21.95, 0.37], [0.6955, 18.22, 0.37], [0.6955, 25.92, 0.37], [0.6955, 26.41, 0.37], [0.6964, 22.64, 0.37], [0.6966, 28.75, 0.37], [0.85, 25.97, 0.37], [0.85, 28.17, 0.37], [0.8513, 17.25, 0.37], [0.9982, 26.17, 0.37], [0.9982, 26.72, 0.37], [0.9982, 27.28, 0.37], [0.9983, 27.84, 0.37], [0.9996, 21.83, 0.37], [1.0, 16.88, 0.37], [1.0, 22.39, 0.37]],
        "D[LA]": [[0.1003, 31.61, 0.37], [0.2021, 20.39, 0.37], [0.2021, 22.81, 0.37], [0.2021, 31.92, 0.37], [0.3008, 32.56, 0.37], [0.3999, 32.49, 0.37], [0.4011, 20.61, 0.37], [0.4038, 22.66, 0.37], [0.5013, 22.14, 0.37], [0.5016, 21.07, 0.37], [0.5025, 20.3, 0.37], [0.5035, 31.0, 0.37], [0.6028, 20.46, 0.37], [0.605, 29.57, 0.37], [0.7005, 19.18, 0.37], [0.7009, 19.71, 0.37], [0.7018, 27.86, 0.37], [0.8016, 17.7, 0.37], [0.8019, 26.92, 0.37], [0.8071, 18.26, 0.37], [0.9043, 26.41, 0.37], [0.9995, 16.51, 0.37], [0.9996, 26.31, 0.37]],
        "D[MG]": [[0.0962, 26.48, 0.37], [0.0963, 22.86, 0.37], [0.0966, 16.74, 0.37], [0.0966, 27.03, 0.37], [0.0971, 27.57, 0.37], [0.0983, 21.23, 0.37], [0.1972, 16.65, 0.37], [0.1978, 26.69, 0.37], [0.1979, 20.19, 0.37], [0.1979, 29.23, 0.37], [0.1984, 27.3, 0.37], [0.1998, 23.85, 0.37], [0.2844, 27.62, 0.37], [0.2982, 25.56, 0.37], [0.2992, 30.21, 0.37], [0.2999, 16.87, 0.37], [0.3, 19.11, 0.37], [0.316, 27.51, 0.37], [0.3476, 27.41, 0.37], [0.3496, 30.4, 0.37], [0.3792, 27.3, 0.37], [0.3991, 17.57, 0.37], [0.3994, 30.23, 0.37], [0.3998, 16.48, 0.37], [0.4108, 27.19, 0.37], [0.4951, 29.4, 0.37], [0.4976, 26.59, 0.37], [0.498, 24.96, 0.37], [0.4983, 14.9, 0.37], [0.4987, 15.45, 0.37], [0.5039, 28.93, 0.37], [0.5979, 12.93, 0.37], [0.5989, 21.58, 0.37], [0.6, 26.03, 0.37], [0.6003, 30.37, 0.37], [0.6008, 27.85, 0.37], [0.6996, 16.99, 0.37], [0.6997, 10.15, 0.37], [0.7005, 31.39, 0.37], [0.7013, 25.38, 0.37], [0.7018, 26.55, 0.37], [0.8, 24.68, 0.37], [0.8005, 25.16, 0.37], [0.8013, 6.9, 0.37], [0.8016, 11.8, 0.37], [0.8021, 32.02, 0.37], [0.8984, 5.95, 0.37], [0.899, 24.1, 0.37], [0.899, 32.4, 0.37], [0.8991, 3.43, 0.37], [1.0, 23.66, 0.37], [1.0, 32.29, 0.37], [1.0, 32.8, 0.37]],
        "D[ML]": [[0.0, 22.38, 0.37], [0.0015, 16.79, 0.37], [0.0038, 26.49, 0.37], [0.2026, 26.62, 0.37], [0.2061, 23.88, 0.37], [0.2071, 27.15, 0.37], [0.2076, 21.11, 0.37], [0.2083, 16.94, 0.37], [0.209, 29.07, 0.37], [0.3236, 26.33, 0.37], [0.3313, 20.2, 0.37], [0.3313, 27.98, 0.37], [0.3345, 24.25, 0.37], [0.3363, 29.87, 0.37], [0.4038, 17.85, 0.37], [0.4774, 19.28, 0.37], [0.5025, 25.76, 0.37], [0.5038, 23.53, 0.37], [0.5061, 31.25, 0.37], [0.5075, 29.58, 0.37], [0.5376, 19.32, 0.37], [0.5977, 19.35, 0.37], [0.6579, 19.39, 0.37], [0.6615, 22.45, 0.37], [0.663, 25.16, 0.37], [0.6667, 32.01, 0.37], [0.668, 30.62, 0.37], [0.801, 24.2, 0.37], [0.8025, 20.84, 0.37], [0.8031, 32.26, 0.37], [0.8049, 21.46, 0.37], [0.806, 30.82, 0.37], [0.9995, 19.98, 0.37], [0.9995, 31.46, 0.37], [1.0, 22.54, 0.37]],
    },
}
PHONON_CURVE["Ir"] = {
    "T_K": 295, "struct": "fcc", "ref": IR_REF, "nu_scale": MEV_THZ,
    "nearest": True,
    "digitised": ("These points were read from Fig. 1 of the source, drawn as "
                  "vectors in the PDF, so their coordinates are exact.  The source "
                  "gives no branch labels, so they are scored against the nearest "
                  "computed branch, a lower bound.  A second experiment (Wang et "
                  "al. 2025, 300&nbsp;K) agrees on &Gamma;&ndash;X to 0.35&nbsp;meV and "
                  "lies up to 1.8&nbsp;meV higher on the &Gamma;&ndash;L transverse branch."),
    #  Sigma past K returns along U-X, as for aluminium; the half of [0z1]
    #  beyond W is the image of X-W, since (0,z,1) is equivalent to (0,1-z,1).
    #  The figure draws each point of a DEGENERATE transverse branch twice, two
    #  filled circles at identical coordinates (31 such positions along Gamma-X
    #  and Gamma-L); that is one measurement, so it is kept once.
    "seg_override": {
        "[0zz]": [("G", "K", 0.0, 0.75), ("U", "X", 0.75, 1.0)],
        "[0z1w]": [("W", "X", 0.5, 1.0)],
    },
    "branches": {
        "D[00z]": [[0.1019, 3.585, 0.0], [0.1542, 5.393, 0.0], [0.2066, 6.985, 0.0], [0.2066, 10.673, 0.0], [0.2502, 8.431, 0.0], [0.3025, 10.022, 0.0], [0.3025, 15.411, 0.0], [0.3548, 11.18, 0.0], [0.4072, 12.192, 0.0], [0.4072, 19.1, 0.0], [0.4595, 13.097, 0.0], [0.5075, 13.928, 0.0], [0.5075, 21.885, 0.0], [0.5555, 14.652, 0.0], [0.6078, 15.556, 0.0], [0.6078, 23.332, 0.0], [0.6602, 16.243, 0.0], [0.7081, 16.677, 0.0], [0.7081, 24.417, 0.0], [0.7561, 17.328, 0.0], [0.8084, 17.762, 0.0], [0.8608, 17.907, 0.0], [0.9131, 18.124, 0.0], [0.9655, 18.196, 0.0]],
        "D[zzz]": [[0.0498, 5.249, 0.0], [0.1019, 5.466, 0.0], [0.1019, 9.95, 0.0], [0.1497, 7.925, 0.0], [0.1497, 13.856, 0.0], [0.2017, 9.516, 0.0], [0.2017, 17.871, 0.0], [0.2495, 10.384, 0.0], [0.2495, 20.294, 0.0], [0.3016, 11.397, 0.0], [0.3016, 22.536, 0.0], [0.3536, 12.626, 0.0], [0.3536, 23.476, 0.0], [0.4014, 13.639, 0.0], [0.4014, 24.489, 0.0], [0.4535, 14.145, 0.0], [0.4535, 25.357, 0.0]],
        "D[0zz]": [[0.0473, 4.055, 0.0], [0.0988, 4.308, 0.0], [0.0988, 5.032, 0.0], [0.0988, 7.852, 0.0], [0.1502, 6.334, 0.0], [0.1502, 7.563, 0.0], [0.1502, 11.903, 0.0], [0.2016, 8.286, 0.0], [0.2016, 10.022, 0.0], [0.2016, 14.579, 0.0], [0.2531, 10.239, 0.0], [0.2531, 16.894, 0.0], [0.296, 11.758, 0.0], [0.296, 18.919, 0.0], [0.3474, 13.205, 0.0], [0.3474, 20.872, 0.0], [0.3988, 14.218, 0.0], [0.3988, 22.174, 0.0], [0.4503, 14.941, 0.0], [0.4503, 23.223, 0.0], [0.4931, 15.303, 0.0], [0.4931, 23.332, 0.0], [0.5446, 15.52, 0.0], [0.5446, 22.753, 0.0], [0.596, 15.375, 0.0], [0.596, 21.523, 0.0], [0.6475, 15.664, 0.0], [0.6475, 20.728, 0.0], [0.6989, 16.098, 0.0], [0.6989, 20.149, 0.0], [0.7461, 16.641, 0.0], [0.7461, 19.679, 0.0], [0.7932, 20.366, 0.0], [0.8447, 17.762, 0.0], [0.8447, 20.221, 0.0], [0.8961, 18.268, 0.0], [0.8961, 19.136, 0.0], [0.9475, 18.413, 0.0], [0.9475, 18.63, 0.0], [0.9947, 18.124, 0.0], [0.9947, 18.341, 0.0], [0.9947, 24.344, 0.0]],
        "D[0z1]": [[0.0089, 18.124, 0.0], [0.0089, 24.344, 0.0]],
    },
}

#  Re: ONE frequency, the Raman-active E2g transverse optic mode at Gamma.
#  Stored in THz directly because a points record is not passed through
#  nu_scale; 120.9 cm-1 x 0.0299792458 = 3.6245 THz.
RE_REF = ("H. Olijnyk, A. P. Jephcoat and K. Refson, Europhys. Lett. <b>53</b>, "
          "504 (2001), Table I: the transverse-optical &Gamma;-point (E<sub>2g</sub>) "
          "mode by high-pressure Raman spectroscopy, &nu;<sub>0</sub> = 120.9 "
          "&plusmn; 0.3 cm<sup>&minus;1</sup>, the zero-pressure coefficient of a "
          "fit to data up to 60 GPa (Table II gives 121 cm<sup>&minus;1</sup> at "
          "0 GPa), here 3.6245 THz.  <b>One frequency, not a dispersion.</b>  The "
          "temperature is not stated; pressure was calibrated by ruby "
          "fluorescence and the measurement is taken as room temperature.  "
          "Checked on ruthenium in the same paper: its 190.6 cm<sup>&minus;1</sup> "
          "(23.63 meV) agrees with the &Gamma;-point transverse optic branch of the "
          "neutron data of Heid et al. 2000, read from their figure at 23.68 meV")
PHONON_CURVE["Re"] = {"T_K": 295, "struct": "hcp", "ref": RE_REF, "points": {
    "G": [(3.6245, 0.009, "G5+")],
}}


V_REF = ("A. Bosak, M. Hoesch, D. Antonangeli, D. L. Farber, I. Fischer and "
         "M. Krisch, Phys. Rev. B <b>78</b>, 020301(R) (2008), Fig. 2(a): inelastic "
         "<b>x-ray</b> scattering on beamline ID28 of the ESRF, 3.0 meV energy "
         "resolution.  Vanadium scatters neutrons almost entirely incoherently, "
         "which is why its measured dispersion is an x-ray one.  <b>Digitised from "
         "the figure</b>, which is vector in the PDF, so the coordinates are exact: "
         "86 points, every marker in the panel, on the branches the paper labels L, T, "
         "T<sub>1</sub> and T<sub>2</sub> (polarisations in its Table I).  No error "
         "bars are given.  The temperature is not stated and the measurement is "
         "taken as room temperature.  Checked against a number the reading did not "
         "use: the initial slopes of L and T along [100], L along [110] and L and T "
         "along [111] agree with the sound velocities of the elastic constants to "
         "within 3.5&nbsp;%")
PHONON_CURVE["V"] = {
    "T_K": 295, "struct": "bcc", "ref": V_REF, "nu_scale": MEV_THZ,
    "digitised": ("These points were read from Fig. 2(a) of the source, drawn as "
                  "vectors in the PDF, so their coordinates are exact.  The source "
                  "labels every branch, so they are scored by polarisation like a "
                  "tabulated measurement.  The two transverse branches along "
                  "&Gamma;&ndash;N are nearly degenerate up to &xi; &asymp; 0.15, and "
                  "the paper reports an anomaly in each."),
    #  Two markers sit on the H line: the [z00] scan ends there and the [zzz]
    #  scan starts there.  H is a point where the branches are degenerate, so
    #  one goes to each line and neither assignment changes a score.
    "branches": {
        "L[z00]": [[0.050, 4.34, 0.0], [0.100, 8.45, 0.0],
                  [0.200, 15.77, 0.0], [0.300, 22.07, 0.0],
                  [0.400, 26.76, 0.0], [0.501, 29.29, 0.0],
                  [0.550, 29.74, 0.0], [0.601, 29.46, 0.0],
                  [0.634, 29.45, 0.0], [0.668, 29.58, 0.0],
                  [0.701, 30.01, 0.0], [0.734, 30.07, 0.0],
                  [0.767, 30.33, 0.0], [0.800, 30.96, 0.0],
                  [0.851, 31.33, 0.0], [0.900, 31.84, 0.0],
                  [1.000, 31.98, 0.0]],
        "T[z00]": [[0.050, 1.82, 0.0], [0.100, 3.63, 0.0],
                  [0.200, 6.90, 0.0], [0.300, 10.67, 0.0],
                  [0.400, 15.87, 0.0], [0.501, 21.28, 0.0],
                  [0.601, 26.17, 0.0], [0.701, 29.60, 0.0],
                  [0.800, 30.98, 0.0], [0.900, 32.01, 0.0],
                  [1.000, 32.66, 0.0]],
        "L[zzz]": [[0.050, 6.89, 0.0], [0.100, 13.77, 0.0],
                  [0.200, 24.32, 0.0], [0.300, 29.35, 0.0],
                  [0.349, 29.31, 0.0], [0.400, 27.62, 0.0],
                  [0.425, 26.60, 0.0], [0.450, 25.99, 0.0],
                  [0.475, 25.31, 0.0], [0.499, 25.34, 0.0],
                  [0.524, 25.03, 0.0], [0.549, 23.89, 0.0],
                  [0.599, 21.77, 0.0], [0.699, 20.19, 0.0],
                  [0.799, 24.25, 0.0], [0.899, 30.23, 0.0],
                  [1.000, 32.45, 0.0]],
        "T[zzz]": [[0.100, 6.71, 0.0], [0.200, 11.97, 0.0],
                  [0.300, 18.01, 0.0], [0.400, 22.99, 0.0],
                  [0.499, 25.44, 0.0], [0.549, 25.36, 0.0],
                  [0.599, 25.13, 0.0], [0.632, 25.50, 0.0],
                  [0.666, 26.15, 0.0], [0.699, 26.74, 0.0],
                  [0.749, 28.09, 0.0], [0.799, 29.37, 0.0],
                  [0.899, 31.60, 0.0]],
        "L[zz0]": [[0.050, 5.95, 0.0], [0.100, 11.29, 0.0],
                  [0.150, 16.47, 0.0], [0.200, 20.23, 0.0],
                  [0.250, 23.35, 0.0], [0.300, 25.86, 0.0],
                  [0.350, 27.20, 0.0], [0.400, 28.20, 0.0],
                  [0.425, 28.36, 0.0], [0.450, 28.04, 0.0],
                  [0.475, 28.06, 0.0], [0.500, 28.22, 0.0]],
        "T1[zz0]": [[0.051, 2.86, 0.0], [0.100, 5.71, 0.0],
                   [0.200, 10.89, 0.0], [0.300, 16.05, 0.0],
                   [0.400, 19.17, 0.0], [0.500, 20.03, 0.0]],
        "T2[zz0]": [[0.051, 2.81, 0.0], [0.100, 5.26, 0.0],
                   [0.150, 8.14, 0.0], [0.200, 11.32, 0.0],
                   [0.251, 14.69, 0.0], [0.300, 18.73, 0.0],
                   [0.350, 22.50, 0.0], [0.400, 25.97, 0.0],
                   [0.451, 27.54, 0.0], [0.500, 27.82, 0.0]],
    }}


#  Which drawn segment each branch lies on, and over what range of the paper's
#  zeta.  The path is build_library.SC_PATH["fcc"], G-X-W-K-G-L-U-W-L-K.
#  Sigma is drawn as K->G, so a [zz0] branch runs backwards along it and only
#  the part with zeta <= 0.75 is on the drawn path at all; the rest of the
#  branch, from K out to X, is measured but never passed by this path and is
#  dropped rather than folded onto a segment it does not belong to.
#  Keyed by structure, because the same bracket means a different line in a
#  different lattice: [zzz] runs Gamma -> L in fcc and Gamma -> P -> H in bcc.
#  A direction may cover more than one drawn segment, which is why the value
#  is a list.
SEGMENT = {
    "fcc": {
        "[00z]": [("G", "X", 0.0, 1.0)],
        "[zzz]": [("G", "L", 0.0, 0.5)],
        "[zz0]": [("G", "K", 0.0, 0.75)],
        "[0zz]": [("G", "K", 0.0, 0.75)],
        "[1z0]": [("X", "W", 0.0, 0.5)],
        "[0z1]": [("X", "W", 0.0, 0.5)],   # calcium's name for the same line
    },
    "bcc": {
        "[z00]": [("G", "H", 0.0, 1.0)],
        "[zz0]": [("G", "N", 0.0, 0.5)],
        #  Lambda reaches P at zeta = 1/2 and (1,1,1) at zeta = 1, and
        #  (1,1,1) - (1,1,0) = (0,0,1) is H, so the far half of the branch is
        #  the P -> H segment rather than more of Gamma -> P.  Barium was
        #  measured out to 0.9 and folding that onto Gamma -> P would have
        #  drawn the second half of the branch on top of the first.
        "[zzz]": [("G", "P", 0.0, 0.5), ("P", "H", 0.5, 1.0)],
    },
    #  hcp.  M sits at zeta = 1/2 along [z00] and A at 1/2 along [00z], but
    #  [zz0] reaches K at 1/3 and KEEPS GOING - Stassis measured out to 0.4.
    #  Those points are outside the first zone.  Folding by subtracting b1
    #  sends zeta = 1/3 to a K point and zeta = 1/2 to an M point, and the
    #  segment between them has the same length as the standard K-M edge, so
    #  the far part of the branch belongs on M-K and not on more of Gamma-K.
    #  Dropping it on Gamma-K would have drawn the branch back over itself.
    "hcp": {
        "[00z]": [("G", "A", 0.0, 0.5)],
        "[z00]": [("G", "M", 0.0, 0.5)],
        "[zz0]": [("G", "K", 0.0, 1.0 / 3.0), ("K", "M", 1.0 / 3.0, 0.5)],
    },
}


def direction(label):
    """the bracketed part of a branch label: T1[0zz] -> [0zz]"""
    i = label.find("[")
    return label[i:] if i >= 0 else None


def points_on(el, seg_a, seg_b, struct="fcc"):
    """[(t, nu, sigma, label)] for every measured point on segment a -> b

    t is the fraction along that segment, so a caller needs only the two
    endpoints' x on the drawn path.  Points outside the segment are dropped,
    never extrapolated.  A segment given in the reverse sense to the table
    (K->G against Sigma's G->K) is handled by the caller, which knows which
    way it is drawing.
    """
    rec = PHONON_CURVE.get(el)
    if not rec:
        return []
    out = []
    #  A paper may measure zeta as |q| a / 2pi rather than as the component
    #  along the direction, which puts L at sqrt(3)/2 and K at 3 sqrt(2)/4
    #  instead of at 1/2 and 3/4.  Lead does exactly that and says so - its
    #  Table II is headed "the zone boundary is at aq/2pi = 1.061".  The
    #  values stay as printed and the RANGE moves, so nothing is converted
    #  here and the table can still be checked against the paper by eye.
    table = dict(SEGMENT.get(struct, {}))
    table.update(rec.get("seg_override", {}))
    #  A paper may tabulate ENERGIES rather than frequencies.  The numbers
    #  are kept exactly as printed so the table can still be checked against
    #  the paper by eye, and the conversion happens here instead.
    k = rec.get("nu_scale", 1.0)
    for label, pts in rec.get("branches", {}).items():
        for a_, b_, lo, hi in table.get(direction(label), ()):
            if (a_, b_) != (seg_a, seg_b):
                continue
            for z, nu, sig in pts:
                if lo <= z <= hi:
                    out.append(((z - lo) / (hi - lo), round(nu * k, 4),
                                round(sig * k, 4), label))
    return out


def summary():
    for el in sorted(PHONON_CURVE):
        r = PHONON_CURVE[el]
        n = sum(len(v) for v in r["branches"].values())
        print(f"{el:3s} {r['T_K']:4d} K  {len(r['branches'])} branches  "
              f"{n:3d} points")


if __name__ == "__main__":
    summary()
