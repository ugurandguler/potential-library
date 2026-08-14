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
            trunc = ("hard - phi2 is cut at rcut2 and does NOT vanish there, "
                     "so this set is for static properties, not dynamics"
                     if taper <= 0 else
                     f"switched from {taper:g} of each cutoff to the cutoff, "
                     "quintic, C2 - energy is conserved in MD")
            #  the measured verdict, which is not the same question as whether
            #  the fit reproduced its targets
            #  the export labels the sets mau/mau_taper/ug/ug_taper and the
            #  screen labels them hard/tap/ug/tap_ug; without this the lookup
            #  silently misses and every file claims "not screened"
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
                warn = (f"# MD: screened - holds its structure at 600 K"
                        f" ({md['T']} K, as equipartition requires).")
            #  The nudge test, which is a different question from the MD
            #  screen and from the phonon screen: both of those can pass while
            #  the lattice fails to survive a 1e-5 A displacement.  Five bcc
            #  records do exactly that, and a user running a defect or an
            #  interface calculation would meet it immediately.
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
            warn = warn + nl + jwarn
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
