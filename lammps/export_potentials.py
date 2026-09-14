#!/usr/bin/env python3
"""
Write the distributable potential files - the thing a LAMMPS user picks up.

Everything else here generates a file into a scratch directory, uses it and
forgets it.  That is right for a validator and wrong for a library: what has
been demonstrated so far is that the parameters and the pair style work, not
that anyone else can run them.  This produces the set.

Four of them, because there are four fitted potentials and they are not
interchangeable:

    potentials/<El>.ugur              hard truncation, the reference set
    potentials/<El>_taper.ugur        switched - use this one for MD
    potentials/<El>.ugur.ang          hard truncation, with the angular factor
    potentials/<El>_taper.ugur.ang    switched, with the angular factor

The naming follows LAMMPS's own convention on both axes, and both halves of it
matter.  The **extension tracks the pair style**, as `.eam` / `.eam.alloy` /
`.eam.fs` do, so a user reading the filename knows which style the file needs
without opening it.  The **variant goes in the stem**, as `Cu_u3` and
`Cu_mishin1` do, which is what keeps the hard and switched sets from colliding:
they were in separate directories under identical names, and a user copying
both into one place would have silently overwritten one with the other.  Two
potentials with the same name is a worse trap than two with a confusing one.

Every file carries its own provenance in the header, including which pair style
it needs and whether it is switched, because the one class of mistake this
project has already made twice is running a potential under a truncation it was
not fitted with.  A file that says `taper -1` and one that says `taper 0.85`
look identical at a glance and are different potentials.

The `ug/` files require `pair_style ugur/ang`; `pair_style ugur` refuses them
rather than dropping the angular term silently.

    python export_potentials.py
"""
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "standalone"))
import refdata          # noqa: E402

OUT = os.path.join(HERE, "potentials")
def phi2_at_cut(rec):
    """phi2 at the cutoff, in eV - the size of the step a hard cut leaves

    Written here rather than taken from the library because it is a property
    of the record's own parameters and has to follow them if they change.
    Evaluated just inside rcut2, since at rcut2 exactly the pair term is
    defined to be zero by the cutoff rather than by the function.

    The STANDALONE latdyn is used even for the angular records, which is
    normally the mistake curve_mae warns about.  It is safe here and only
    here: the angular factor multiplies the three-body term, phi2 is the
    same function in both modules, and phi2 is all this asks for.
    """
    import latdyn as _L
    return float(_L.Potential.from_record(rec).phi2(rec["rcut2"] - 1e-6))

KEYS = ("m", "D", "alpha", "r0", "gamma", "C", "alpha3", "rcut2", "rcut3")

#  (label, library key, pair style, filename suffix, description)
SETS = (
    ("mau", None, "ugur", ".ugur",
     "hard truncation, the reference parameters"),
    ("mau_taper", "tap", "ugur", "_taper.ugur",
     "cutoff switched off over its outer 15 %"),
    ("ug", "ug", "ugur/ang", ".ugur.ang",
     "hard truncation, with the angular factor"),
    ("ug_taper", "tap_ug", "ugur/ang", "_taper.ugur.ang",
     "switched, with the angular factor"),
    #  The shell-gap re-cut.  Same pair styles and the same 0.85 taper as the
    #  two switched sets above, fitted at a different cutoff.  Where
    #  standalone/recommend_recut.py finds it passes every cold and warm
    #  screen, <El>_recut.ugur is the RECOMMENDED set for molecular dynamics
    #  and the description below is overridden for that element.  It breaks
    #  lithium, sodium, potassium, rubidium and caesium once warm (negative
    #  expansion, C11 stiffening), and the molybdenum and tungsten records
    #  failed selection outright - an earlier version of this comment said
    #  the re-cut repaired those two.
    ("rc", "rc", "ugur", "_recut.ugur",
     "CANDIDATE, re-cut, switched"),
    ("rc_ug", "rc_ug", "ugur/ang", "_recut.ugur.ang",
     "CANDIDATE, re-cut, switched, with the angular factor"),
    #  The dispersion-selected records.  Same targets and the same 0.000 %
    #  elastic residual as the reference set above, chosen among the solutions
    #  that residual cannot distinguish by which one describes the measured
    #  phonons better.  They ship for the reason the re-cut candidates do: a
    #  reader who wants to check the comparison on the page needs the file it
    #  was made with.  Hard truncation, so the same MD warning applies.
    ("disp", "hard_disp", "ugur", "_disp.ugur",
     "CANDIDATE, dispersion-selected, hard truncation"),
    #  The nudge-constrained refits.  Carried in library.json since the
    #  nudge filter ran and quoted on the page ever since, but never written
    #  out - so one candidate set was downloadable and the other was not,
    #  which is the inconsistency rather than the export.  None is withdrawn
    #  and all five pass the displacement test.
    ("nudge", "tap_nudge", "ugur", "_nudge_taper.ugur",
     "CANDIDATE, nudge-constrained, switched"),
    #  The force-matched arm.  Fitted to DFT forces rather than to the
    #  elastic constants alone, which is what repairs the surfaces - copper
    #  goes from 3.4x the DFT surface energy to 1.2x with the elastic
    #  constants intact.  Switched, taper 0.85, same pair style as the
    #  shipped tapered set.
    #
    #  A CANDIDATE, and by the library's own gate rather than by caution:
    #  four of the seven exceed the three-body ceiling that standalone/fit.py
    #  enforces with no tolerance, so `fit.py` would not have kept them, and
    #  vanadium is not a harmonic minimum at all.  Each file says which.
    ("force", "tap_force", "ugur", "_force_taper.ugur",
     "CANDIDATE, force-matched, switched"),
)

#  Which published force database each force-matched element was fitted to.
#  Taken from refit/DATA_REFERENCES.md, whose own header records that two of
#  its citations were WRONG when checked against Crossref: the ColabFit dataset
#  names say "PRM2019" and the paper is Phys. Rev. Materials 4, 093802 (2020).
#  A dataset's NAME is not its citation, so these are the resolved ones.
FORCE_SRC = {
    "Mo": "Byggmastar, Nordlund, Djurabekova, Phys. Rev. Materials 4, 093802 (2020)",
    "Nb": "Byggmastar, Nordlund, Djurabekova, Phys. Rev. Materials 4, 093802 (2020)",
    "Ta": "Byggmastar, Nordlund, Djurabekova, Phys. Rev. Materials 4, 093802 (2020)",
    "V": "Byggmastar, Nordlund, Djurabekova, Phys. Rev. Materials 4, 093802 (2020)",
    "W": "Byggmastar, Hamedani, Nordlund, Djurabekova, Phys. Rev. B 100, 144105 (2019)",
    "Cu": "Fellman, Byggmastar, Granberg, Nordlund, Djurabekova, "
          "Phys. Rev. Materials 9, 053807 (2025)",
    "Ni": "Fellman, Byggmastar, Granberg, Nordlund, Djurabekova, "
          "Phys. Rev. Materials 9, 053807 (2025)",
}

HEADER = """\
# {el} - Ugur interatomic potential, {what}
#
# pair_style {style}
# pair_coeff * * {fn} {el}
#
# Truncation: {trunc}
{warn}
#
{fitted}
#
# Format is LAMMPS's Tersoff convention - one line per ordered triple
# (centre, leg, leg) - so that this file and a multi-element one are read the
# same way.  For a single element that is one line.
#
# el1 el2 el3  m D alpha r0 gamma C alpha3 rcut2 rcut3 taper lam2 lam4
{el} {el} {el} {vals}
"""


def main():
    lib = json.load(open(os.path.join(ROOT, "standalone", "library.json")))
    #  Which sets survive molecular dynamics at all.  Written into the file
    #  itself, not only the README: a user who copies one file out of a
    #  directory takes the header with it and leaves the README behind, and a
    #  potential that turns a crystal into a 4500 K liquid should not be
    #  distributed with nothing but a filename to warn them.
    jig = {}
    jp = os.path.join(HERE, "jiggle_test.json")
    if os.path.exists(jp):
        for k, v in json.load(open(jp)).items():
            el_, st_ = k.split("|")
            jig[(el_, st_)] = v
    scr = {}
    sp = os.path.join(HERE, "md_screen_all.json")
    if os.path.exists(sp):
        for k, v in json.load(open(sp)).items():
            el, st = k.split("|")
            scr[(el, st)] = v
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    print(f"{'set':12s}{'element':>9s}{'skipped':>9s}   what")
    print("-" * 58)
    total = 0
    #  counted in the loop rather than inferred from filenames afterwards:
    #  ".ugur" is a suffix of nothing but "_taper.ugur" is, and getting that
    #  arithmetic wrong put a "0 elements" in the README for a set of 38
    written = {}
    os.makedirs(OUT)
    for name, key, style, suffix, what in SETS:
        n, skipped = 0, []
        for el in sorted(lib):
            rec = lib[el] if key is None else lib[el].get(key)
            if not rec or any(k not in rec for k in KEYS):
                skipped.append(el)
                continue
            taper = rec.get("taper") or -1.0
            vals = [rec[k] for k in KEYS] + [taper,
                                             rec.get("lam2", 0.0),
                                             rec.get("lam4", 0.0)]
            #  The first header line names what the file is.  For an element
            #  whose re-cut is recommended, "CANDIDATE" on that line would
            #  contradict the paragraph below it, and the switched file must
            #  say where MD users should go instead.
            what_el = what
            if (lib[el].get("md_recommended") or {}).get("set") == "rc":
                if name == "rc":
                    what_el = ("RECOMMENDED for molecular dynamics, re-cut, "
                               "switched")
                elif name == "mau_taper":
                    what_el = (what + " - for MD on " + el + " the re-cut "
                               + el + "_recut.ugur is recommended instead")
            #  What each truncation is for, with the measurement behind it.
            #
            #  The hard sets reproduce the MEASURED dispersion better - 9.6 %
            #  mean over the 32 elements that carry a neutron curve against
            #  12.5 % for the switched ones, and better in 25 of the 32 - and
            #  they cannot be run at temperature.  Both halves are measured
            #  and neither is a preference.  The discontinuity is what does
            #  it, so each file carries its own rather than an average: it
            #  runs from 0.25 meV for palladium to 123 meV for yttrium, a
            #  factor of five hundred, and "not for dynamics" is not equally
            #  true across that range.
            jump = abs(phi2_at_cut(rec)) * 1000.0        # meV
            kT296 = 8.617333e-5 * 296 * 1000.0           # meV
            if taper <= 0:
                trunc = (
                    "hard - phi2 does not vanish at rcut2, it stops at "
                    f"{-abs(jump):.3g} meV ({jump / kT296:.3g} of k_B T at "
                    "296 K).  USE for static properties and for lattice "
                    "dynamics: measured against neutron dispersion the hard "
                    "sets average 9.6 % over the 32 elements that have one, "
                    "against 12.5 % for the switched sets, and are closer in "
                    "25 of them.  DO NOT use for molecular dynamics: "
                    "copper, whose step is 8.4 meV, drifts 350 meV/atom/ns in "
                    "NVE and climbs from 296 K to over 1100 K in 200 ps")
                if name == "disp":
                    trunc += (
                        ".  CANDIDATE, not part of the published library: the "
                        "same targets and the same elastic residual as the "
                        "reference set, chosen among the solutions that "
                        "residual cannot tell apart by which describes the "
                        "measured dispersion better.  It also satisfies the "
                        "E3/E2 and compression constraints the fitter "
                        "enforces today and the reference record does not.  "
                        "Read the section on the page before using it")
            elif name == "nudge":
                trunc = (
                    f"switched from {taper:g} of each cutoff to the cutoff, "
                    "quintic, C2 - energy is conserved in MD.  CANDIDATE, not "
                    "part of the published library.  This is the best "
                    "solution in the same pool as the shipped switched record "
                    "that survives a 1e-5 A displacement, which the shipped "
                    "one does not; it costs elastic accuracy to get there and "
                    "the page carries both errors side by side.  Read them "
                    "before using this")
            elif name in ("rc", "rc_ug"):
                #  A different cutoff, so the 12.5 % measured on the published
                #  switched sets is not theirs.  curve_mae.py, rerun
                #  2026-09-13: over the 17 elements that have a re-cut record
                #  and a measured dispersion and were not rejected, the re-cut
                #  reaches 8.3 % against 11.2 % for the published switched
                #  set.  Most of that gap is lithium, sodium, potassium,
                #  rubidium and caesium, the elements it BREAKS warm; without
                #  them it is 8.6 % against 9.2 % over 12, better by more
                #  than half a point in 5.  So the dispersion is not the
                #  reason to use it.  What recommend_recut.py weighs is the
                #  warm screens, per element, and the file says the result.
                mr = lib[el].get("md_recommended") or {}
                base = (f"switched from {taper:g} of each cutoff to the "
                        "cutoff, quintic, C2 - energy is conserved in MD.  ")
                if name == "rc" and mr.get("set") == "rc":
                    trunc = base + (
                        f"RECOMMENDED for molecular dynamics on {el}, in place"
                        f" of {el}_taper.ugur, which still ships.  It passes "
                        "every screen of the rule - " + mr.get("rule", "")
                        + " - which standalone/recommend_recut.py applies to "
                        "every element that has a re-cut, not by hand")
                elif name == "rc":
                    why = "; ".join(mr.get("why") or [])
                    trunc = base + (
                        "CANDIDATE, not the recommended set"
                        + (f": {why}" if why else "")
                        + ".  For MD on this element use "
                        f"{el}_taper.ugur")
                else:
                    trunc = base + (
                        "CANDIDATE: the angular re-cut is never the "
                        "recommended set.  Its dispersion advantage is mostly "
                        "the alkalis, which it breaks warm; read the notes "
                        "below before using it")
            else:
                trunc = (
                    f"switched from {taper:g} of each cutoff to the cutoff, "
                    "quintic, C2 - energy is conserved in MD: measured drift "
                    "0.4 meV/atom/ns for copper, 876x less than its hard "
                    "twin.  USE for molecular dynamics, which no hard set "
                    "can do.  It pays for that on the dispersion: the "
                    "switched sets average 12.5 % against the hard sets' "
                    "9.6 %, because switching the pair term off over the "
                    "outer 15 % moves the force constants at the largest "
                    "separations, which is the short-wavelength end, and "
                    "nothing in the fit sees it")
            #  the measured verdict, which is not the same question as whether
            #  the fit reproduced its targets
            #  the export labels the sets mau/mau_taper/ug/ug_taper and the
            #  screen labels them hard/tap/ug/tap_ug; without this the lookup
            #  silently misses and every file claims "not screened"
            #  The candidate arms keep their screen inside their own record
            #  rather than in md_screen_all.json, so they are looked up there.
            #  The four shipped sets keep the old path untouched - reading
            #  rec first for those would change what they report.
            if name in ("rc", "rc_ug", "disp", "nudge", "force"):
                md = rec.get("md_screen")
            else:
                md = scr.get((el, {"mau": "hard", "mau_taper": "tap",
                                   "ug": "ug", "ug_taper": "tap_ug"}[name]))
            if md is None:
                warn = "# MD: not screened."
            elif md.get("lost"):
                n0, n1 = md["lost"]
                warn = (f"# MD: DO NOT USE.  At 600 K this structure"
                        f" disintegrates - {n0} atoms become {n1} and LAMMPS"
                        f" stops.")
            elif md.get("collapsed"):
                warn = (f"# MD: DO NOT USE.  At 600 K the crystal collapses:"
                        f" it reaches {md['T']} K and its potential energy"
                        f" falls below the static lattice.  Static properties"
                        f" from this set are still valid.")
            elif md.get("T", 300) > 400:
                warn = (f"# MD: SUSPECT.  At 600 K this runs at {md['T']} K"
                        f" where equipartition gives 300; something is"
                        f" releasing energy.  Check before trusting it.")
            else:
                warn = (f"# MD: structure screened - holds its shape at"
                        f" 600 K ({md['T']} K, as equipartition requires)."
                        f"  This tested the SHAPE, not energy conservation;"
                        f" for that see the truncation line above.")
            #  What three warm screens say about the candidates, which no
            #  cold screen can see.  The alkalis pass every 0 K test - they
            #  hit every fitted target and are dynamically stable along the
            #  whole symmetry path, MORE stable than the published arm, which
            #  has imaginary modes at H-N for K, Rb and Cs - and then fail
            #  every warm one.
            #  The nudge test, which is a different question from the MD
            #  screen and from the phonon screen: both of those can pass while
            #  the lattice fails to survive a 1e-5 A displacement.  Five bcc
            #  records do exactly that, and a user running a defect or an
            #  interface calculation would meet it immediately.
            #  same split as the screen above: the candidates carry their
            #  own jiggle result, the shipped four keep the old lookup
            if name in ("rc", "rc_ug", "disp", "nudge", "force"):
                jg = rec.get("jiggle")
            else:
                jg = jig.get((el, {"mau": "hard", "mau_taper": "tap",
                                   "ug": "ug", "ug_taper": "tap_ug"}[name]))
            nl = chr(10)
            if jg is None:
                jwarn = "#" + nl + "# Nudge test: not run."
            elif jg.get("ok"):
                jwarn = ("#" + nl
                         + "# Nudge test: passed - displaced by 1e-5 A and"
                         " relaxed, the lattice returns.")
            else:
                jwarn = ("#" + nl
                         + "# Nudge test: FAILED.  Displace every atom by"
                         f" 1e-5 A and relax, and the crystal keeps"
                         f" {jg['keep']:.0f}x that displacement and settles"
                         f" {abs(jg['dE']) * 1000:.2f} meV/atom LOWER.  The"
                         " reference lattice is not this potential's minimum,"
                         " so the elastic constants above are the constants of"
                         " a structure it does not hold.  Static reference"
                         " only; do not use for defects, surfaces or dynamics.")
            #  The symmetry-path stability screen.  It is a DIFFERENT question
            #  from the 600 K run above: that one asks whether the crystal
            #  survives being heated, this one whether the reference lattice
            #  is a harmonic minimum at all.  They do not cover each other -
            #  caesium's tapered set holds its structure at 600 K, passes the
            #  nudge test, and still carries a mode at -1.5 cm^-1 between two
            #  symmetry points.  A half-shifted Monkhorst-Pack mesh cannot see
            #  a mode confined to a line however fine it is made, which is why
            #  the fit's own check calls these records stable.
            dyn = rec.get("dyn") or {}
            if dyn.get("stable") is False:
                mn = dyn.get("min_path_cm1")
                if mn is None:
                    mn = dyn.get("most_neg_cm1")
                warn = (warn + nl + "#" + nl
                        + "# LATTICE STABILITY: this set is NOT a harmonic"
                        " minimum." + nl
                        + f"# Along the Setyawan-Curtarolo path its frequencies"
                        f" reach {mn} cm^-1," + nl
                        + f"# with {100 * dyn.get('imag_frac', 0.0):.2f} % of"
                        " the sampled modes imaginary.  That is a separate"
                        + nl + "# failure from the 600 K screen above and is"
                        " not covered by it.")
            elif dyn:
                warn = (warn + nl + "#" + nl
                        + "# Lattice stability: screened on the 8^3 and 9^3"
                        " meshes AND along the" + nl
                        + "# symmetry path; no imaginary modes.")
            #  THE THREE-BODY CEILING.  `standalone/fit.py` refuses any
            #  solution whose three-body sum exceeds 0.30 of the two-body one,
            #  with no tolerance, because that ratio is what the compression
            #  escape and the vacancy failures were traced to.  The
            #  force-matched arm was fitted without it and five of its seven
            #  records land above the line - by 0.0002 to 0.0012, small
            #  numbers against a bound that is enforced exactly.  A reader
            #  holding one of these files is holding a set the library's own
            #  acceptance test would have rejected, and has to be told so
            #  here rather than on a web page they may never have opened.
            if rec.get("gate") == "reject":
                warn = (warn + nl + "#" + nl
                        + "# THREE-BODY CEILING: this set is OVER it."
                        + nl
                        + f"# E3/E2 = {rec.get('ratio')} against the bound of"
                        f" 0.30, exceeding it by {rec.get('violation')}."
                        + nl
                        + "# standalone/fit.py applies that bound with no"
                        " tolerance, so this solution" + nl
                        + "# would not have been accepted by the fit that"
                        " produced the shipped sets.")
            warn = warn + nl + jwarn
            #  What the measured dispersion says about this candidate,
            #  which is NOT a second copy of the finite-temperature warning
            #  further down - it is the other half of the same fact.  The
            #  re-cut reaches 8.8 % over the 18 elements that have both a
            #  candidate record and a neutron curve, against 11.6 % for the
            #  published switched arm on the same 18, closer in 14 of them.
            #  Its largest gains are the alkalis, and the alkalis are exactly
            #  what it breaks warm: a 0 K dispersion and an elastic constant
            #  are both properties of the curvature AT the minimum, so an arm
            #  can have that neighbourhood right and the shape of the well
            #  away from it wrong.
            #  The same table as RC_DISP in standalone/make_gui.py, and it
            #  must stay the same: this copy had drifted older than the page's
            #  and told K_recut.ugur it was worse than the published arm, when
            #  the current scoring says better.  Recomputed 2026-09-13 with
            #  curve_mae.py, arms rc and tap.  Ir is nearest-branch against a
            #  digitised figure, so both of its numbers are lower bounds.
            RC_DISP = {"Ag": (13.0, 12.7), "Au": (4.2, 4.4), "Ba": (8.3, 9.6),
                       "Ca": (8.6, 8.3), "Cs": (6.4, 12.9), "Cu": (9.3, 10.0),
                       "Ir": (9.0, 9.9), "K": (11.8, 13.9), "Li": (7.0, 25.9),
                       "Mg": (8.3, 6.6), "Na": (5.8, 16.2), "Ni": (12.1, 13.4),
                       "Pb": (12.9, 18.3), "Pd": (3.5, 3.6), "Pt": (5.0, 4.8),
                       "Rb": (7.6, 11.8), "W": (17.2, 17.8), "Yb": (8.5, 8.6)}
            warm = ""
            if name in ("rc", "rc_ug") and el in RC_DISP:
                mine, pub = RC_DISP[el]
                verdict = ("better" if mine < pub - 0.5 else
                           "worse" if mine > pub + 0.5 else "the same")
                warm = ("#" + nl
                        + f"# Measured dispersion: {mine:g} % against"
                        f" {pub:g} % for the published switched" + nl
                        + f"# arm - {verdict}." + nl)
                if el == "Ir":
                    warm += ("# Both are nearest-branch scores against points"
                             " read from a figure," + nl
                             + "# so both are lower bounds." + nl)
                if el in ("Na", "K", "Rb", "Cs", "Li"):
                    warm += ("# Read that beside the finite-temperature note"
                             " below rather than instead of" + nl
                             + "# it.  Both are the same fact: a dispersion at"
                             " 0 K is a property of the" + nl
                             + "# curvature AT the minimum, and this arm has"
                             " that neighbourhood right" + nl
                             + "# while the shape of the well away from it is"
                             " wrong." + nl)
                #  Recomputed 2026-09-13 from library.json, per arm.  The
                #  sentence this replaces put both arms under "41 records"
                #  and quoted 300 K elastic constants 3.4 -> 2.7 %, 7 better
                #  and 3 worse, which no slice of the current data
                #  reproduces; the expansion figures did, for `rc` alone.
                elif name == "rc":
                    warm += ("# Over the 21 accepted re-cut records: thermal"
                             " expansion error 30.7 % -> 25.6 %" + nl
                             + "# median, 10 better and 10 worse; elastic"
                             " constants near 300 K against the" + nl
                             + "# room-temperature values 5.2 % -> 3.4 %"
                             " median, 11 better and 10 worse." + nl)
                else:
                    warm += ("# Over the 20 accepted angular re-cut records:"
                             " thermal expansion error" + nl
                             + "# 35.0 % -> 34.7 % median, 11 better and 8"
                             " worse; elastic constants near 300 K" + nl
                             + "# against the room-temperature values 3.7 %"
                             " -> 3.4 % median, 9 better and 11 worse." + nl)
                warm = warm.rstrip(nl)
            warn = warn + (nl + warm if warm else "")
            #  The 600 K screen asks whether the crystal survives, which is
            #  not the same question as whether it behaves.  Sodium's
            #  candidate holds its structure and still expands the wrong way
            #  and stiffens as it is heated, so the screen alone reads as
            #  reassurance it has not earned.  This adds the two finite-
            #  temperature failures to the header, where somebody who copies
            #  the file will meet them.
            if name in ("rc", "rc_ug"):
                #  Four of these records did not pass selection at all -
                #  molybdenum and tungsten in both arms, on a mode at about
                #  -13 cm^-1 on the symmetry path.  They are kept as the
                #  control the others are read against.  Without this line the
                #  file says CANDIDATE, reports the 600 K screen and the nudge
                #  test as passed, and never mentions that it was rejected.
                gr = rec.get("ground") or {}
                if rec.get("stable") is False or gr.get("ok") is False:
                    bits = []
                    if rec.get("stable") is False:
                        bits.append("it carries a mode at"
                                    f" {rec.get('min_cm1', 0.0):.1f} cm^-1 on"
                                    " the symmetry path")
                    if gr.get("ok") is False:
                        bits.append(f"it puts {gr.get('lowest')} lowest,"
                                    f" {abs(gr.get('rel', 0)):.0f} meV/atom"
                                    f" below {gr.get('want')}")
                    warn = (warn + nl + "#" + nl
                            + "# REJECTED: this record did NOT pass selection -"
                            + nl + "# " + "; ".join(bits) + "." + nl
                            + "# It is kept as the control the accepted"
                            " records are read against," + nl
                            + "# not as a potential.  Do not use it.")
                ft = []
                xp = rec.get("expansion") or {}
                a, ae = xp.get("alpha_1e6"), xp.get("alpha_exp_1e6")
                if a is not None and ae and a < 0 < ae:
                    ft.append(f"thermal expansion comes out NEGATIVE"
                              f" ({a:.0f}e-6/K against a measured {ae:.0f})")
                et = (lib[el].get("elasticT") or {}).get(name) or {}
                pts = [q for q in et.get("pts", [])
                       if q.get("T") is not None
                       and q["T"] <= 0.7 * (et.get("Tmelt") or 1e9)]
                if len(pts) >= 4 and pts[0].get("C11"):
                    rise = (pts[-1]["C11"] - pts[0]["C11"]) / pts[0]["C11"]
                    if rise > 0.05:
                        ft.append(f"C11 RISES {100 * rise:.0f} % between 0 K"
                                  f" and {pts[-1]['T']:.0f} K instead of"
                                  " softening")
                if ft:
                    warn = (warn + nl + "#" + nl
                            + "# FINITE TEMPERATURE: this candidate passes the"
                            " screens above and fails" + nl
                            + "# below room temperature - "
                            + ("; ".join(ft)) + "." + nl
                            + "# Use it for static properties only, if at"
                            " all.")
            #  Per set, and named.  The force-matched files say which
            #  published database they were fitted to, because a reader
            #  holding one has to be able to check it without the web page.
            if name == "force":
                fitted = ("# Fitted to: DFT FORCES on displaced supercells of "
                          + el + ", from the" + nl
                          + "# database cited here, with the cohesive energy,"
                          " lattice constant and" + nl
                          + "# elastic constants kept as anchors." + nl
                          + "# Source: "
                          + FORCE_SRC.get(el, "see refit/DATA_REFERENCES.md")
                          + ".")
            else:
                fitted = ("# Fitted to: cohesive energy, lattice constant,"
                          " bulk modulus and the elastic" + nl
                          + "# constants of " + el + " at its experimental"
                          " lattice constant.  Nothing else.")
            open(os.path.join(OUT, el + suffix), "w").write(HEADER.format(
                fitted=fitted,
                el=el, fn=el + suffix, style=style, what=what_el, trunc=trunc,
                warn=warn,
                vals=" ".join(f"{v:.17g}" for v in vals)))
            n += 1
        total += n
        written[name] = n
        print(f"{name:12s}{n:9d}{len(skipped):9d}   {what}")
        if skipped:
            print(f"{'':12s}skipped: {' '.join(skipped)}")
    rec_els = sorted(e for e in lib if isinstance(lib[e], dict)
                     and (lib[e].get("md_recommended") or {}).get("set")
                     == "rc")
    #  a README beside them, because a directory of bare numbers is not a
    #  distribution
    open(os.path.join(OUT, "README"), "w").write(f"""\
Ugur interatomic potential - LAMMPS parameter files
===================================================

    <El>.ugur              {written['mau']} elements, pair_style ugur
    <El>_taper.ugur        {written['mau_taper']} elements, pair_style ugur
    <El>.ugur.ang          {written['ug']} elements, pair_style ugur/ang
    <El>_taper.ugur.ang    {written['ug_taper']} elements, pair_style ugur/ang

and the re-cut sets.  <El>_recut.ugur is the RECOMMENDED set for molecular
dynamics on {len(rec_els)} elements and a candidate for the rest;
<El>_recut.ugur.ang is always a candidate.  Recommended for:

    {' '.join(rec_els)}

    <El>_recut.ugur        {written['rc']} elements, pair_style ugur
    <El>_recut.ugur.ang    {written['rc_ug']} elements, pair_style ugur/ang

and the CANDIDATE sets, which are not part of the published library:

    <El>_disp.ugur         {written['disp']} elements, pair_style ugur
    <El>_nudge_taper.ugur  {written['nudge']} elements, pair_style ugur
    <El>_force_taper.ugur  {written['force']} elements, pair_style ugur

They are here because the page compares them with the published sets and a
reader who wants to check that comparison needs the files it was made with.
Every one of them carries, in its own header, what it is and what is known
against it.  Read that before using one.  In one line each:

  _recut        a different pair cutoff.  Where it is recommended it passes
                every cold and warm screen and gets the sign of the
                intrinsic stacking fault right, which the _taper sets get
                wrong.  It BREAKS lithium, sodium, potassium, rubidium and
                caesium once the crystal is warm, and its molybdenum and
                tungsten records failed selection.
  _disp         the same targets and the same elastic residual as the plain
                set, chosen among the solutions that residual cannot tell
                apart by which describes the measured phonon dispersion
                better.  Hard truncation, so the MD prohibition below applies.
  _nudge_taper  the best solution in the same pool as the plain _taper set
                that survives a 1e-5 A displacement, which the shipped one
                does not.  It pays for that in elastic accuracy.
  _force_taper  fitted to DFT FORCES on displaced supercells, not to the
                elastic constants alone.  That is what repairs the surfaces:
                copper falls from 3.4x the DFT surface energy to 1.2x with
                the elastic constants intact.  The price is stated per file -
                four of the seven are over the three-body ceiling that the
                fit enforces without tolerance, and vanadium is not a
                harmonic minimum.  Seven elements only: Cu, Mo, Nb, Ni, Ta,
                V, W.

The extension says which pair style the file needs, the way .eam / .eam.alloy /
.eam.fs do.  The stem says which parameter set it is: a plain name is
hard-truncated, _taper is switched.  Those two are different potentials, fitted
under different truncations, and must not be swapped - each file states which
it is in its own header.

Which one to use
----------------
For molecular dynamics, use a *_taper set - or, for the {len(rec_els)} elements
named above, <El>_recut.ugur, which is recommended there instead.

The hard-truncated potentials do not conserve energy: phi2 does not vanish at
the cutoff, so a neighbour crossing it changes the energy in one step.
Measured, that is a drift of 5 to
7256 meV/atom/ps against 0.03 to 0.28 for the switched sets, and on three of
the hexagonal metals the crystal does not merely drift but comes apart.

For static elastic constants and phonons the hard-truncated sets are the
reference and are better for the hexagonal metals and the alkalis.  Both are
provided for that reason.

The angular sets (ug/) add h(cos theta) = 1 + lam2 P2 + lam4 P4 to the
three-body term.  They need pair_style ugur/ang; pair_style ugur refuses a file
with nonzero weights rather than ignoring the term.

What these were fitted to
-------------------------
Cohesive energy, lattice constant, bulk modulus and elastic constants, at the
experimental lattice constant.  Nothing else - no defect energies, no surface
energies, no melting behaviour, no liquid structure.  Do not assume
transferability to any of those.

Alloys
------
See ALLOYS.md.  A multi-element file can be generated with make_alloy_file.py.
The unlike-pair and unlike-triple entries it writes come from conventional
mixing rules, and those rules have now been measured rather than assumed.

What the measurement says.  On MoNbTaVW the mixed potential reproduces alloy
DFT forces about as well as it reproduces the pure elements it was fitted to -
0.91x the single-element control - so the mixing rule costs nothing on forces.
That does not make the alloy usable.  Under a free relaxation the bcc crystal
leaves bcc entirely, axis ratio 1.412 +/- 0.004 against the Bain value of
sqrt(2), on every seed, while all five constituents stay cubic under the
identical relaxation.  Of the ten binaries of those five elements, five reach
the full Bain value, two sit between, and three stay cubic or nearly so; only
Nb-Ta, at 1.0005, is cleanly cubic.  Use a generated file for forces at a
fixed cell shape, not for anything that lets the cell relax.

Which sets can be mixed.  An angular set - ug, tap_ug, rc_ug - is refused,
because there is no mixing rule for lam2 and lam4.  A set an element does not
have is refused by name; the force-matched set covers only Cu, Mo, Nb, Ni, Ta,
V and W.  That leaves the hard-truncated, tapered and re-cut sets.

One thing to know first.  A library rebuilt from the shipped fit.json by
standalone/refresh.py carries only the hard-truncated sets, so --set tap -
which the generator recommends for anything going into molecular dynamics -
will report that there is no tapered record.  The tapered parameters are in
the library embedded in the published page, docs/index.html; the public tree
ships fit.json and no library.json at all.
""")
    print(f"\n{total} files -> {OUT}")


if __name__ == "__main__":
    main()
