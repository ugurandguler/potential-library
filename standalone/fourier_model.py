#!/usr/bin/env python3
"""
Dispersion from published INTERPLANAR force constants (a Fourier series).

A third way a paper can publish its dispersion, after a list of frequencies
and a Born-von Karman tensor model: as the Fourier coefficients of each
branch along each symmetry direction,

    4 pi^2 M nu^2 = sum_n Phi_n [ 1 - cos(n pi zeta / zeta_max) ]

which is Eq. (1) of Powell, Martel and Woods.  The Phi_n are combinations of
Born-von Karman constants, but the series is fitted branch by branch, so it
reproduces one line of the zone at a time and says nothing off it.  That is
why this is stored as a model curve and drawn as a line: it is an
interpolation of measured points, not the points.

TWO THINGS HAD TO BE PINNED DOWN AND NEITHER WAS GUESSED.

The unit.  The table is headed only "Fourier coefficients"; taking the
printed numbers as dyn/cm puts the zone boundary at 0.17 THz, which is
wrong by a factor of a thousand in energy.  They are 10^3 dyn/cm.

zeta_max.  The paper says only that it "is a function of the particular
branch".  For Sigma the series must end at N, because its value there is a
maximum rather than zero and the [zz0] line returns to a reciprocal-lattice
point at zeta = 1.  For Lambda the branch is labelled Lambda_1(F_1): Lambda
is Gamma-P and F is P-H, so the series spans BOTH and ends at H, not at P.
The evidence is that its endpoint agrees with the Delta endpoint, and Delta
ends at H where all three branches are degenerate - if Lambda stopped at P
there would be no reason for those two numbers to agree, and for tantalum
they differ by 25%.

check() tests the whole reconstruction against Landolt-Boernstein, which is
a different compilation and was not used to set anything above.

    python fourier_model.py
"""
import math

AMU_G = 1.66054e-24
PHI_UNIT = 1.0e3          # the table is in 10^3 dyn/cm

#  Powell, Martel and Woods, Phys. Rev. 171, 727 (1968), Table I, the
#  Molybdenum block, 296 K.  Sigma_4 is not tabulated, so [zz0] carries two
#  of its three branches and the third is simply absent.
MO_296 = {
    "L[z00]":  [119.19, 117.42, -12.99, 16.94, -8.04, 3.15, -3.48],
    "T[z00]":  [114.38, 14.49, -9.46, 5.71, -8.04, 6.04],
    "L[zzz]":  [45.40, 78.38, 63.56, 29.15, -1.66, 2.28, 4.38, 6.94, -6.32,
                -0.27, -4.09, 5.28, -3.00, 2.62, -0.22, 2.52, -1.18, 0.52,
                -1.80, 1.53],
    "T[zzz]":  [112.32, 85.01, -0.42, 22.47, -11.91, -0.88, -2.56, -0.08,
                -0.05, 1.01, -1.77, 2.09],
    "L[zz0]":  [198.31, 11.13, 6.35, -2.30, 4.03, -2.30],
    "T1[zz0]": [76.41, 4.89, -8.07, 5.42, -2.44],
}
#  transverse branches along [z00] and [zzz] are doubly degenerate in bcc
MO_WEIGHT = {"T[z00]": 2, "T[zzz]": 2}
ZMAX = {"[z00]": 1.0, "[zzz]": 1.0, "[zz0]": 0.5}

MO_REF = ("R. M. Powell, P. Martel and A. D. B. Woods, "
          "Phys. Rev. <b>171</b>, 727 (1968), Table I, molybdenum block, "
          "296 K.  The paper plots its measured frequencies and tabulates "
          "the Fourier coefficients of each branch, so what is drawn here "
          "is their series rather than their points")
#  Landolt-Boernstein, through the library's own reference set, for the check
MO_LB = {"H": 5.52, "N_L": 8.14, "N_T": 4.56}


def nu_THz(phi, zeta, zeta_max, mass_amu):
    s = sum(p * (1.0 - math.cos((n + 1) * math.pi * zeta / zeta_max))
            for n, p in enumerate(phi))
    if s <= 0.0:
        return 0.0
    return math.sqrt(s * PHI_UNIT
                     / (4.0 * math.pi ** 2 * mass_amu * AMU_G)) / 1.0e12


def direction(label):
    return label[label.find("["):]


def branches_at(coeffs, label_dir, zeta, mass_amu):
    """[(frequency, branch label)] on one direction, degeneracies expanded"""
    out = []
    for lab, phi in coeffs.items():
        if direction(lab) != label_dir:
            continue
        f = nu_THz(phi, zeta, ZMAX[label_dir], mass_amu)
        out.extend([f] * MO_WEIGHT.get(lab, 1))
    return sorted(out)


def check():
    import refdata
    m = refdata.MASSES["Mo"]
    h = [nu_THz(MO_296[b], 1.0, 1.0, m)
         for b in ("L[z00]", "T[z00]", "L[zzz]", "T[zzz]")]
    n_l = nu_THz(MO_296["L[zz0]"], 0.5, 0.5, m)
    n_t = nu_THz(MO_296["T1[zz0]"], 0.5, 0.5, m)
    print()
    print("Mo, Fourier series against Landolt-Boernstein, which set nothing "
          "above")
    print(f"{'':12s}{'from the series':>18s}{'Landolt-B.':>14s}")
    print(f"{'H':12s}{min(h):8.2f} -{max(h):6.2f}{MO_LB['H']:14.2f}")
    print(f"{'N, L':12s}{n_l:18.2f}{MO_LB['N_L']:14.2f}")
    print(f"{'N, T':12s}{n_t:18.2f}{MO_LB['N_T']:14.2f}")
    #  H is where the three branches must coincide, so their spread is
    #  itself a test of the reconstruction
    ok = (max(h) - min(h) < 0.12
          and abs(sum(h) / 4 - MO_LB["H"]) < 0.10
          and abs(n_l - MO_LB["N_L"]) < 0.15
          and abs(n_t - MO_LB["N_T"]) < 0.15)
    print(f"the four branches agree at H to {max(h) - min(h):.3f} THz, "
          "where bcc symmetry makes them degenerate")
    print("reconstruction " + ("CONFIRMED" if ok else "DOES NOT MATCH - do "
          "not use this model until it does"))
    return ok


if __name__ == "__main__":
    check()
