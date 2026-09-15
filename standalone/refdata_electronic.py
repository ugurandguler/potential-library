"""
Electronic heat-capacity coefficients, gamma, for the thermodynamics panel.

The panel compares the model's C_v at 298 K with the CRC's tabulated C_p.
Two terms separate them that no interatomic potential contains: the lattice
term C_p - C_v = 9 alpha^2 B V_m T, which the page already computes, and the
electronic heat capacity gamma*T of the conduction electrons.  This module is
the source for the second.

SOURCE, read from the rendered page because the PDF page carries no text
layer: C. Kittel, Introduction to Solid State Physics, 8th ed. (Wiley),
Chapter 6, Table 2, p. 146, "Experimental and free electron values of
electronic heat capacity constant gamma of metals", from compilations by
N. Phillips and N. Pearlman.  The row used is "Observed gamma in
mJ mol^-1 K^-2".  36 of the library's 38 elements are in it; lutetium and
ytterbium are not, and they are left out rather than estimated.

WHAT THE NUMBER IS AND IS NOT.  gamma is measured at low temperature, as the
intercept of C/T against T^2 (Kittel eq. 37).  At those temperatures it
carries the electron-phonon mass enhancement, which fades once T is well above
the Debye temperature.  So gamma*T at 298 K is an UPPER bound on the
electronic term, and it overshoots most where gamma is largest - Nb, V, Sc,
Pd, Y, Ta, Pt come out 0.8 to 2.3 J/(mol K) past the tabulated C_p.  No
enhancement factors are applied here: none are in this source, and a number
without a source is not quoted.

A separate module rather than a table in refdata.py: refdata is imported by
the fitting code, and nothing that fits a potential should depend on a value
that only the thermodynamics panel uses.
"""

GAMMA_SOURCE = ("C. Kittel, <i>Introduction to Solid State Physics</i>, 8th ed. "
                "(Wiley), Table 2, p. 146 (from compilations by N. Phillips and "
                "N. Pearlman): observed &gamma;, mJ mol<sup>&minus;1</sup> "
                "K<sup>&minus;2</sup>")

#  mJ mol^-1 K^-2, Kittel 8th ed. Table 2, "Observed gamma" row
GAMMA = {
    "Li": 1.63, "Be": 0.17, "Na": 1.38, "Mg": 1.3, "Al": 1.35,
    "K": 2.08, "Ca": 2.9, "Sc": 10.7, "Ti": 3.35, "V": 9.26, "Cr": 1.40,
    "Fe": 4.98, "Co": 4.73, "Ni": 7.02, "Cu": 0.695,
    "Rb": 2.41, "Sr": 3.6, "Y": 10.2, "Zr": 2.80, "Nb": 7.79, "Mo": 2.0,
    "Ru": 3.3, "Rh": 4.9, "Pd": 9.42, "Ag": 0.646,
    "Cs": 3.20, "Ba": 2.7, "Hf": 2.16, "Ta": 5.9, "W": 1.3, "Re": 2.3,
    "Ir": 3.1, "Pt": 6.8, "Au": 0.729, "Tl": 1.47, "Pb": 2.98,
}
