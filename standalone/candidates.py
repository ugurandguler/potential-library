#!/usr/bin/env python3
"""
The twelve elements considered for the library on 2026-08-03, and what happened
to them.

TEN WERE ACCEPTED and now live in refdata.py: Rb, Cs, Yb, Sc, Y, Lu, Hf, Re, Ru
and Tl.  Every one fits and every one is dynamically stable on the 8^3 U 9^3
union mesh - Yb, Rb and Cs exactly, Tl at 0.36 %, Re at 3.17 %, the worst being
Lu at 13.75 %.  Their rows are in refdata.py rather than here, so that this file
does not become a second table to keep in step with the first.

TWO WERE NOT.

**Technetium.**  Three soft inputs in a row.  Its bulk modulus is one of the two
Kittel prints in parentheses, meaning an estimate rather than a measurement; its
cohesive energy is the least certain of the twelve at +-2.5 %, and the fit
imposes Ecoh exactly at every trial point; and its elastic constants have
C12 = C13 = 199 exactly, which is what a single determination with an assumed
relation looks like rather than an average over measurements.

**Osmium.**  III/29a has no complete elastic set for it, so there is nothing to
fit against.  Filling the gap from another compilation would break the one
property that makes the elastic column defensible - every constant in it comes
from a single compiler with a single set of conventions.

Their reference data is kept below.  If a complete osmium determination turns up,
or a better technetium one, the row is ready and only needs moving.

Sources are the same as for the accepted ten, and are documented in refdata.py:
Kittel Table 3 of chapter 1 for a0 and c/a, Kittel Table 3 of chapter 3 for B,
Brewer LBL-3720 Rev. (1977) for Ecoh, Landolt-Boernstein III/29a for Cij.
"""

#  B in GPa, a0 in Angstrom, Ecoh in eV/atom (positive magnitude).  Values
#  Kittel prints in parentheses are estimates, not measurements - both of these
#  bulk moduli are among them.
NOT_ADOPTED = {
    "Tc": dict(struct="hcp", a0=2.74, c_over_a=4.40/2.74, Ecoh=6.851, B=297.0,
               B_is_estimate=True, ecoh_unc=0.173,
               Cij=dict(C11=433.0, C12=199.0, C13=199.0, C33=470.0, C44=177.0),
               why="B is an estimate; Ecoh +-2.5 %, the least certain of the "
                   "twelve; C12 = C13 exactly suggests one determination"),
    "Os": dict(struct="hcp", a0=2.74, c_over_a=4.32/2.74, Ecoh=8.170, B=418.0,
               B_is_estimate=True, ecoh_unc=0.039, Cij=None,
               why="no complete elastic set in III/29a"),
}

MASSES = {"Tc": 97.90721, "Os": 190.23}


if __name__ == "__main__":
    import refdata
    took = [e for e in ("Rb","Cs","Yb","Sc","Y","Lu","Hf","Re","Ru","Tl")
            if e in refdata.ELEMENTS]
    print(f"accepted {len(took)}/10: {', '.join(took)}")
    print(f"alinmayan {len(NOT_ADOPTED)}: " +
          ", ".join(f"{k} ({v['why']})" for k, v in NOT_ADOPTED.items()))
