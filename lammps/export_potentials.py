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
    #  The shell-gap re-cut candidates.  NOT part of the shipped library and
    #  named so they cannot be mistaken for it: same pair styles and the same
    #  0.85 taper as the two switched sets above, but fitted at a different
    #  cutoff.  Before using either one, read the verdict column on the page -
    #  they repair molybdenum and tungsten and they break the alkalis, whose
    #  elastic constants come out three times too stiff at 300 K and whose
    #  thermal expansion comes out negative.
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
)

HEADER = """\
# {el} - Ugur interatomic potential, {what}
#
# pair_style {style}
# pair_coeff * * {fn} {el}
#
# Truncation: {trunc}
{warn}
#
# Fitted to: cohesive energy, lattice constant, bulk modulus and the elastic
# constants of {el} at its experimental lattice constant.  Nothing else.
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
            #  What each truncation is for, with the measurement behind it.
            #
            #  The hard sets reproduce the MEASURED dispersion better - 9.5 %
            #  mean over the 29 elements that carry a neutron curve against
            #  12.7 % for the switched ones, and better in 23 of the 29 - and
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
                    "sets average 9.5 % over the 29 elements that have one, "
                    "against 12.7 % for the switched sets, and are closer in "
                    "23 of them.  DO NOT use for molecular dynamics: "
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
                #  A different cutoff, so the 12.7 % measured on the shipped
                #  switched sets is not theirs.  Measured on the 17 elements
                #  that have both a re-cut record and a neutron curve, the
                #  re-cut arm reaches 8.5 % against 11.6 % for the shipped
                #  switched arm on the same 17, and is closer in 13 of them.
                #
                #  That is not a recommendation, and the reason is the whole
                #  point of these files.  The biggest gains are the alkalis -
                #  lithium 25.2 % to 5.8, caesium 16.7 to 3.9, rubidium 13.5
                #  to 5.4 - and those are exactly the elements the re-cut
                #  BREAKS once the crystal is warm.  Both facts have the same
                #  cause: a dispersion at 0 K and an elastic constant are
                #  properties of the curvature AT the minimum, and the re-cut
                #  gets that neighbourhood right while getting the shape of
                #  the well away from it wrong.
                trunc = (
                    f"switched from {taper:g} of each cutoff to the cutoff, "
                    "quintic, C2 - energy is conserved in MD.  CANDIDATE, not "
                    "part of the published library.  Against measured neutron "
                    "dispersion this arm reaches 7.9 % over the 16 elements "
                    "that have both a candidate record and a neutron curve "
                    "and were not rejected outright, against 11.2 % for the "
                    "published switched arm on the same 16.  That gap is "
                    "almost all alkali: sodium, potassium, rubidium and "
                    "caesium are where the arm gains most and are also what "
                    "it breaks warm.  Leave them out and it is 8.3 % against "
                    "8.9 % over 11 elements, better in 7 - a wash.  There is "
                    "no measured reason to prefer this arm for its dispersion")
            else:
                trunc = (
                    f"switched from {taper:g} of each cutoff to the cutoff, "
                    "quintic, C2 - energy is conserved in MD: measured drift "
                    "0.4 meV/atom/ns for copper, 876x less than its hard "
                    "twin.  USE for molecular dynamics, which no hard set "
                    "can do.  It pays for that on the dispersion: the "
                    "switched sets average 12.7 % against the hard sets' "
                    "9.5 %, because switching the pair term off over the "
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
            if name in ("rc", "rc_ug", "disp", "nudge"):
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
            if name in ("rc", "rc_ug", "disp", "nudge"):
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
            warn = warn + nl + jwarn
            #  What the measured dispersion says about this candidate,
            #  which is NOT a second copy of the finite-temperature warning
            #  further down - it is the other half of the same fact.  The
            #  re-cut reaches 8.5 % over the 17 elements that have both a
            #  candidate record and a neutron curve, against 11.6 % for the
            #  published switched arm on the same 17, closer in 13 of them.
            #  Its largest gains are the alkalis, and the alkalis are exactly
            #  what it breaks warm: a 0 K dispersion and an elastic constant
            #  are both properties of the curvature AT the minimum, so an arm
            #  can have that neighbourhood right and the shape of the well
            #  away from it wrong.
            RC_DISP = {"Ag": (13.0, 12.7), "Au": (4.2, 4.4), "Ba": (8.3, 9.6),
                       "Ca": (8.6, 8.3), "Cs": (3.9, 16.7), "Cu": (9.3, 10.0),
                       "K": (16.6, 11.4), "Li": (5.8, 25.2), "Mg": (8.3, 6.6),
                       "Na": (9.4, 15.4), "Ni": (12.1, 13.4), "Pb": (10.6, 15.9),
                       "Pd": (3.5, 3.6), "Pt": (4.8, 4.7), "Rb": (5.4, 13.5),
                       "W": (17.2, 17.8), "Yb": (8.5, 8.6)}
            warm = ""
            if name in ("rc", "rc_ug") and el in RC_DISP:
                mine, pub = RC_DISP[el]
                verdict = ("better" if mine < pub - 0.5 else
                           "worse" if mine > pub + 0.5 else "the same")
                warm = ("#" + nl
                        + f"# Measured dispersion: {mine:g} % against"
                        f" {pub:g} % for the published switched" + nl
                        + f"# arm - {verdict}." + nl)
                if el in ("Na", "K", "Rb", "Cs", "Li"):
                    warm += ("# Read that beside the finite-temperature note"
                             " below rather than instead of" + nl
                             + "# it.  Both are the same fact: a dispersion at"
                             " 0 K is a property of the" + nl
                             + "# curvature AT the minimum, and this arm has"
                             " that neighbourhood right" + nl
                             + "# while the shape of the well away from it is"
                             " wrong." + nl)
                else:
                    warm += ("# On the 41 accepted candidate records the"
                             " re-cut is close to break-even" + nl
                             + "# overall: thermal expansion 30.7 % -> 25.6 %"
                             " median, 10 better and 10" + nl
                             + "# worse; 300 K elastic constants 3.4 % ->"
                             " 2.7 %, 7 better and 3 worse." + nl)
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
            open(os.path.join(OUT, el + suffix), "w").write(HEADER.format(
                el=el, fn=el + suffix, style=style, what=what, trunc=trunc,
                warn=warn,
                vals=" ".join(f"{v:.17g}" for v in vals)))
            n += 1
        total += n
        written[name] = n
        print(f"{name:12s}{n:9d}{len(skipped):9d}   {what}")
        if skipped:
            print(f"{'':12s}skipped: {' '.join(skipped)}")
    #  a README beside them, because a directory of bare numbers is not a
    #  distribution
    open(os.path.join(OUT, "README"), "w").write(f"""\
Ugur interatomic potential - LAMMPS parameter files
===================================================

    <El>.ugur              {written['mau']} elements, pair_style ugur
    <El>_taper.ugur        {written['mau_taper']} elements, pair_style ugur
    <El>.ugur.ang          {written['ug']} elements, pair_style ugur/ang
    <El>_taper.ugur.ang    {written['ug_taper']} elements, pair_style ugur/ang

and the CANDIDATE sets, which are not part of the published library:

    <El>_recut.ugur        {written['rc']} elements, pair_style ugur
    <El>_recut.ugur.ang    {written['rc_ug']} elements, pair_style ugur/ang
    <El>_disp.ugur         {written['disp']} elements, pair_style ugur
    <El>_nudge_taper.ugur  {written['nudge']} elements, pair_style ugur

They are here because the page compares them with the published sets and a
reader who wants to check that comparison needs the files it was made with.
Every one of them carries, in its own header, what it is and what is known
against it.  Read that before using one.  In one line each:

  _recut        a different pair cutoff.  Repairs molybdenum and tungsten and
                BREAKS the alkalis once the crystal is warm - their elastic
                constants come out three times too stiff at 300 K and their
                thermal expansion negative.
  _disp         the same targets and the same elastic residual as the plain
                set, chosen among the solutions that residual cannot tell
                apart by which describes the measured phonon dispersion
                better.  Hard truncation, so the MD prohibition below applies.
  _nudge_taper  the best solution in the same pool as the plain _taper set
                that survives a 1e-5 A displacement, which the shipped one
                does not.  It pays for that in elastic accuracy.

The extension says which pair style the file needs, the way .eam / .eam.alloy /
.eam.fs do.  The stem says which parameter set it is: a plain name is
hard-truncated, _taper is switched.  Those two are different potentials, fitted
under different truncations, and must not be swapped - each file states which
it is in its own header.

Which one to use
----------------
For molecular dynamics, use a *_taper set.  The hard-truncated potentials do
not conserve energy: phi2 does not vanish at the cutoff, so a neighbour
crossing it changes the energy in one step.  Measured, that is a drift of 5 to
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
See ALLOYS.md.  A multi-element file can be generated with make_alloy_file.py,
but the unlike-pair and unlike-triple entries it writes come from conventional
mixing rules that have not been tested against alloy data.
""")
    print(f"\n{total} files -> {OUT}")


if __name__ == "__main__":
    main()
