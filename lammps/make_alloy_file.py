#!/usr/bin/env python3
"""
Write a multi-element ugur potential file from the pure-element library.

The format is LAMMPS's Tersoff convention: one line per ordered triple
(centre, leg, leg), N**3 lines for N elements, and two-body parameters read
only from the lines where the two legs agree.  See ALLOYS.md for why the
problem has that shape - phi3 depends on the legs through r1 + r2 and does not
factorise, so C and alpha3 belong to a triple and cannot be built out of pair
quantities.

The mixing rules here are conventional and **none of them is validated**.  They
are written into the header of every file produced, so that whoever is holding
one does not have to come back here to find out what it assumes.  The intended
use is as a starting point for a fit against alloy data, with the refitted
unlike entries written straight back into the same file.

    python make_alloy_file.py Cu Ni > CuNi.ugur.alloy
    python make_alloy_file.py --set tap Cu Ni Al

`--set tap` uses the tapered parameters, which is what anything going into MD
should use; the default is the hard-truncated reference set.

WHICH SETS ARE AVAILABLE DEPENDS ON THE LIBRARY YOU HAVE.  This reads
`standalone/library.json`, and the published tree does not ship one - it ships
`fit.json` and `refresh.py` builds the rest.  That chain rebuilds the 0 K
analytic layer only, so a freshly rebuilt library carries the hard-truncated
sets and NOT `tap`, `rc`, `tap_force` or the rest; asking for one of those
gets a refusal naming it rather than a file.  The arms live in the library
that ships inside the page, and `UG_LIB` points this at another tree's copy.

An ANGULAR set - `ug`, `tap_ug`, `rc_ug` - is refused outright: there is no
mixing rule for lam2/lam4 and writing them as zero would produce a file that
runs and is not that arm.

    python make_alloy_file.py --correct --set tap Cu Ni
    python make_alloy_file.py --chi-d 1.03 --chi-r0 0.99 --set tap Cu Ni

`--correct` applies the measured cross correction for the system where one
exists; `--chi-d` and `--chi-r0` give it explicitly.  Both are off by default.
Off is the honest default and also the dangerous one: with no correction the
mixing rule puts the formation energy on the wrong side of zero for every
compound measured so far, so an uncorrected file will not reproduce that the
alloy forms at all.  CORRECTIONS below carries what is known.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

NL = chr(10)

PAIR = ("m", "D", "alpha", "r0", "gamma")
CUT = ("rcut2", "rcut3")


def geometric(vals):
    """|product|**(1/n), which is the geometric mean where signs allow"""
    p = 1.0
    for v in vals:
        p *= abs(v)
    return p ** (1.0 / len(vals))


def arithmetic(vals):
    return sum(vals) / len(vals)


#  Measured cross corrections, refit/alloy_chi3.py.  chi_D multiplies the
#  cross well depth and chi_r0 the cross bond length; 1.0 for both is the
#  uncorrected rule below.  Each was solved compound by compound against the
#  DFT formation energy AND the excess volume, then taken at the median, and
#  the spread over the compounds is quoted because it is what says whether one
#  pair of numbers is a description or a curve fit:
#
#    system   chi_D   spread   chi_r0   spread   compounds
#    Cu-Ni    1.029    0.4 %   0.993     2.4 %   3
#    Au-Cu    1.180    8.2 %   0.925     9.2 %   5
#    Al-Ni    1.324   23.9 %   0.926     5.0 %   9
#
#  Cu-Ni is the one that is really described by two numbers - 0.4 per cent over
#  its compounds - and its chi_r0 is 1 to within 0.7 per cent, so it is a
#  one-parameter correction.  Al-Ni is NOT: a 24 per cent spread in chi_D means
#  no single depth serves its nine compounds and the cross terms have to be
#  fitted properly.  Applying the Al-Ni entry buys mechanical stability - AlNi
#  B2 goes from C11 = 84.2 below C12 = 118.5, which is Born-unstable, to
#  299/180/121 - while overshooting DFT's 204/134/113 badly.  Stability, not
#  accuracy.  Read ALLOYS.md before using any of them.
CORRECTIONS = {
    frozenset(("Cu", "Ni")): (1.029, 0.993),
    frozenset(("Au", "Cu")): (1.180, 0.925),
    frozenset(("Al", "Ni")): (1.324, 0.926),
}


def mix_pair(a, b, chi_d=1.0, chi_r0=1.0):
    """parameters for the A-B bond from the two pure elements"""
    out = {}
    #  a well depth multiplies, a length adds: the Lorentz-Berthelot split,
    #  which is a convention and not a derivation
    out["D"] = chi_d * geometric([a["D"], b["D"]])
    out["r0"] = chi_r0 * arithmetic([a["r0"], b["r0"]])
    out["alpha"] = arithmetic([a["alpha"], b["alpha"]])
    #  m and gamma are exponents with no dimensional argument either way
    out["m"] = arithmetic([a["m"], b["m"]])
    out["gamma"] = arithmetic([a["gamma"], b["gamma"]])
    #  cutoffs average so that the switch window stays the same fraction of
    #  the range for every bond
    for k in CUT:
        out[k] = arithmetic([a[k], b[k]])
    return out


def mix_triple(recs):
    """C and alpha3 for a triple, given (centre, leg, leg) records

    The sign of C is taken from the centre rather than from the product.  C
    changes sign across this library, so a geometric mean is undefined where
    the product is negative, and the centre is the atom whose bonds are being
    bent.  It is a convention; ALLOYS.md says so and says it can be argued
    with.
    """
    cs = [r["C"] for r in recs]
    sign = -1.0 if recs[0]["C"] < 0 else 1.0
    return {"C": sign * geometric(cs),
            "alpha3": arithmetic([r["alpha3"] for r in recs])}


HEADER = """\
# Ugur potential, {n} elements: {els}
# Generated by make_alloy_file.py from the {which} parameter set.
#
# FORMAT - LAMMPS Tersoff convention.  One line per ordered triple
# (centre, leg, leg).  Two-body parameters are read ONLY from the lines where
# the two legs agree; elsewhere they are repeated for readability and ignored.
# Three-body parameters (C, alpha3) are read from every line.
#
#   el1 el2 el3   m  D  alpha  r0  gamma  C  alpha3  rcut2  rcut3  taper \
lam2 lam4
#
# UNLIKE ENTRIES ARE NOT FITTED.  They come from these rules and nothing has
# tested them:
#   D_AB      = sqrt(D_A D_B)          r0_AB     = (r0_A + r0_B)/2
#   alpha_AB  = (alpha_A + alpha_B)/2  m, gamma  = arithmetic mean
#   rcut_AB   = arithmetic mean        taper     = must be identical
#   C_ABC     = |C_A C_B C_C|^(1/3), sign of the centre's C
#   alpha3_ABC= arithmetic mean of the three
#
# Only C and alpha3 are read from a triple line.  phi3's radial shape - D, m,
# r0, gamma - is taken from the CENTRE element's own (A,A,A) entry, because a
# triple must be symmetric under swapping its legs and the (A,B,C) and (A,C,B)
# lines carry different pair columns.  The pair columns therefore serve phi2
# alone and are not doing double duty.
# Like entries (el1 = el2 = el3) are the fitted pure-element parameters and are
# exact.  Treat an alloy built from this file as a starting point for a fit,
# not as a prediction.  See ALLOYS.md.
"""


def main():
    args = sys.argv[1:]
    which = "hard-truncated"
    if "--set" in args:
        i = args.index("--set")
        which = "tapered" if args[i + 1] == "tap" else args[i + 1]
        del args[i:i + 2]
    #  The correction is OFF unless asked for, and it is off by default on
    #  purpose: the uncorrected file is what the mixing rules say, and a
    #  generator that silently applied a fitted number would make the two
    #  impossible to tell apart afterwards.
    chi_d = chi_r0 = 1.0
    chi_src = "none"
    for flag, name in (("--chi-d", "d"), ("--chi-r0", "r")):
        if flag in args:
            i = args.index(flag)
            val = float(args[i + 1])
            if name == "d":
                chi_d = val
            else:
                chi_r0 = val
            chi_src = "given on the command line"
            del args[i:i + 2]
    auto = "--correct" in args
    if auto:
        args.remove("--correct")
    if len(args) < 2:
        raise SystemExit("two or more elements are required, for example: "
                         "python make_alloy_file.py Cu Ni")

    #  UG_LIB lets this read another tree's library, which is how an arm that
    #  lives outside this one is generated without copying the file in.
    lib = json.load(open(os.environ.get(
        "UG_LIB", os.path.join(ROOT, "standalone", "library.json"))))
    recs = {}
    for el in args:
        if el not in lib:
            raise SystemExit(f"{el} is not in the library")
        #  ANY ARM, BY NAME.  This line used to read
        #      r = lib[el]["tap"] if which == "tapered" else lib[el]
        #  so every name except `tap` fell through to the top-level hard-cut
        #  record while the header went on saying which set had been asked
        #  for.  Measured: `--set tap_ug Re Na` wrote "from the tap_ug
        #  parameter set" above m = 10.897, D = 0.5781, taper off - the
        #  hard-cut numbers - where tap_ug is m = 1.2048, D = 1.0175,
        #  taper 0.85.  A header that misnames its own contents is worse than
        #  no header, and the hard-cut set is the one whose own file says
        #  "DO NOT use for molecular dynamics".
        if which == "tapered":
            r = lib[el].get("tap")
        elif which == "hard-truncated":
            r = lib[el]
        else:
            r = lib[el].get(which)
        if not r:
            raise SystemExit(f"no {which} record for {el}")
        recs[el] = r

    #  AN ANGULAR ARM CANNOT BE MIXED, and asking for one used to write a file
    #  with lam2 = lam4 = 0 that ran perfectly well and was not that arm - the
    #  UG radial parameters carrying MAU's angle-free physics, with no error
    #  and no warning.
    #
    #  The kernel is not the limitation: pair_ugur.cpp reads lam2/lam4 per
    #  TRIPLE and refuses them under a non-angular style.  What is missing is a
    #  rule for combining them, and there is no obvious one - lambda weights
    #  the angular factor of a (centre, leg, leg) triple, so it belongs to the
    #  triple rather than to a pair, and whether it should follow the centre
    #  alone or all three is a physical question, not an arithmetic one.
    #
    #  It matters across the library: only 7 of the 38 records carry
    #  |lam2| + |lam4| below 0.01, and the largest are near 3.1 (Re, Na, Li).
    ang = {e: r for e, r in recs.items()
           if abs(r.get("lam2") or 0.0) + abs(r.get("lam4") or 0.0) > 1e-12}
    if ang:
        carried = ", ".join("%s(lam2=%.4g, lam4=%.4g)"
                            % (e, r.get("lam2") or 0.0, r.get("lam4") or 0.0)
                            for e, r in sorted(ang.items()))
        raise SystemExit(
            ("%s is an angular set and there is no mixing rule for lam2/lam4."
             % which) + NL
            + "  carried by: " + carried + NL
            + "Generating it would write lam2 = lam4 = 0, which runs and is"
              " NOT this arm." + NL
            + "Use a non-angular set for a multi-element file: tap, rc,"
              " tap_force," + NL
            + "hard_disp, tap_nudge, or the hard-truncated default.")

    tapers = {r.get("taper") or -1.0 for r in recs.values()}
    if len(tapers) > 1:
        raise SystemExit(f"the elements carry different tapers: {tapers}; "
                         "pencere modelin parcasi, bag basina parametre degil")
    taper = tapers.pop()

    if auto:
        key = frozenset(args)
        if key not in CORRECTIONS:
            raise SystemExit(
                "--correct: no measured correction for %s; the ones that "
                "exist are %s.  Give --chi-d/--chi-r0 explicitly, or leave "
                "them off and treat the file as the starting point it is."
                % (" ".join(args),
                   ", ".join("-".join(sorted(k)) for k in CORRECTIONS)))
        chi_d, chi_r0 = CORRECTIONS[key]
        chi_src = "--correct, the measured median for this system"

    out = [HEADER.format(n=len(args), els=" ".join(args), which=which)]
    if chi_d != 1.0 or chi_r0 != 1.0:
        out.append(
            "# CROSS CORRECTION APPLIED: chi_D = %g on the cross well depth,\n"
            "# chi_r0 = %g on the cross bond length (%s).  The unlike-pair D\n"
            "# and r0 columns below are therefore NOT what the mixing rules\n"
            "# in this header give; multiply them out to recover those.  See\n"
            "# ALLOYS.md for what these numbers were measured against and,\n"
            "# more importantly, for what they were not.\n#"
            % (chi_d, chi_r0, chi_src))
    else:
        out.append(
            "# No cross correction: this is the uncorrected mixing rule, and\n"
            "# for every system measured so far that rule puts the formation\n"
            "# energy on the WRONG SIDE OF ZERO.  --correct applies the\n"
            "# measured one where it exists.  ALLOYS.md has the numbers.\n#")
    for a in args:
        for b in args:
            pb = (mix_pair(recs[a], recs[b], chi_d, chi_r0)
                  if a != b else dict(recs[a]))
            for c in args:
                tri = (mix_triple([recs[a], recs[b], recs[c]])
                       if not (a == b == c) else
                       {"C": recs[a]["C"], "alpha3": recs[a]["alpha3"]})
                #  lam2/lam4 are written as ZERO, and that is not a
                #  placeholder waiting for a rule - see the refusal in main().
                #  Reaching here means the arm carries no angular term.
                v = [pb["m"], pb["D"], pb["alpha"], pb["r0"], pb["gamma"],
                     tri["C"], tri["alpha3"], pb["rcut2"], pb["rcut3"],
                     taper, 0.0, 0.0]
                out.append(f"{a:<3s} {b:<3s} {c:<3s} "
                           + " ".join(f"{x:.12g}" for x in v))
    print("\n".join(out))


if __name__ == "__main__":
    main()
