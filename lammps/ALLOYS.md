# Alloys: the decision, the file format, and what is not yet decided

## The problem the published method does not solve

The original work treats an alloy by the **virtual-crystal** route: take the
alloy's own lattice constant, cohesive energy, ionic energy and bulk modulus as
inputs, and solve for a single effective potential. There is no mixing rule for
the parameters, and none is needed — including for the three-body term, which
is otherwise the hardest part of extending a many-body potential past one
species, since a triplet can have a centre of one element and legs of two
others.

That is a complete method for what it was built for and it is the wrong tool
here. A virtual crystal has one kind of atom, so in molecular dynamics it
cannot show local chemistry, short-range order, segregation, antisite defects
or anything else that makes an alloy an alloy. It would reproduce the elastic
constants it was fitted to and nothing else.

So a mixing rule has to be chosen. This file records what was chosen and why,
because the choice is a modelling decision and not an implementation detail.

## What has to be indexed by what

Write the two terms with species labels. For a pair of atoms of elements
A and B at separation r:

    E2 = phi2_AB(r)

and for a triplet with **centre** of element A and legs of elements B and C at
distances r1 and r2 from it:

    E3 = phi3_ABC(r1 + r2) * h_ABC(cos theta) * S_AB(r1) * S_AC(r2)

Three things follow, and the third is the one that dictates the file format.

1. `phi2` needs a parameter set per **unordered pair** (A,B): m, D, alpha, r0,
   gamma, plus its cutoff rcut2_AB.
2. The switch acts **per leg**, so its cutoff belongs to a pair, not to the
   triplet: S_AB uses rcut3_AB.
3. `phi3` is a function of x = r1 + r2. It does **not** factorise into
   something belonging to leg AB times something belonging to leg AC. So C and
   alpha3 belong to the **ordered triple** (A; B, C), symmetric under B↔C, and
   there is no way to build them from pair quantities without inventing a rule.

Point 3 is why this is a Tersoff-shaped problem rather than a Lennard-Jones
shaped one, and why a `pair_coeff` line with one file per element cannot work.

## The format

Exactly LAMMPS's Tersoff convention, because it is the format that already
exists for this shape of problem and the one a LAMMPS user will expect:

    # el1 el2 el3   m  D  alpha  r0  gamma  C  alpha3  rcut2  rcut3  taper  lam2  lam4
    Cu  Cu  Cu      ...
    Cu  Cu  Ni      ...
    ...

with `pair_coeff * * ugur.alloy Cu Ni`, and one line for every ordered triple —
N³ lines for N elements, 8 for a binary, 27 for a ternary.

The reading rule is Tersoff's, unchanged:

- **Two-body parameters are taken only from lines where el2 == el3.** The
  `Cu Ni Ni` line carries the Cu–Ni pair parameters. On every other line the
  two-body columns are ignored and are conventionally written as the same
  values, so a file stays readable.
- **Three-body parameters are taken from every line**, indexed by
  (centre, leg, leg).

Following the existing convention exactly means a user who has written a
Tersoff or Stillinger-Weber file already knows how to write this one, and it
means the entries can be *fitted* to alloy data rather than generated, which
is the whole point of not hard-coding a mixing rule in the C++.

`lam2` and `lam4` are the UG Legendre weights. They are in the format because
leaving them out would mean a second incompatible format later; `pair_style
ugur` reads them and requires them to be zero, since it implements the
published angle-free phi3. A separate style will read them properly.

## The default mixing rule, and its status

A file has to come from somewhere before anyone has fitted an alloy, so
`make_alloy_file.py` generates one from the pure-element library. The rules are
the conventional ones and each is written into the file's own header so that no
one has to read this document to know what they are holding:

| quantity | rule | why |
| --- | --- | --- |
| D_AB | sqrt(D_A · D_B) | geometric, as for a well depth |
| r0_AB | (r0_A + r0_B)/2 | arithmetic, as for a length |
| alpha_AB | (alpha_A + alpha_B)/2 | arithmetic; it is an inverse length |
| m_AB, gamma_AB | arithmetic mean | no better argument exists |
| rcut2_AB, rcut3_AB | arithmetic mean | keeps the switch window fractional |
| taper | must be identical | a mixed window is a different model |
| C_ABC | cube root of C_A·C_B·C_C, sign of the centre's | geometric, extended |
| alpha3_ABC | arithmetic mean of the three | as for alpha |

**None of this is validated.** These rules reproduce the pure-element
parameters exactly when all three labels agree — which is a consistency
requirement, not evidence — and beyond that they are guesses of the kind that
happen to be standard. An alloy potential built this way should be treated as a
starting point for a fit, not as a prediction. The intended path is: generate,
then refit the unlike-pair and unlike-triple entries against alloy data, and
the format exists so that the refitted numbers can simply be written back in.

The sign convention for C_ABC deserves a note. C changes sign between elements
in this library, so a geometric mean is undefined where the product is
negative. The rule taken is |C_A·C_B·C_C|^(1/3) carrying the sign of the
centre's own C, because the centre is the atom whose bonds are being bent. It
is a convention and it is arbitrary; it is written here so that it can be
argued with.

## What is deliberately not done

- **No fitting of alloys.** Nothing here produces an alloy potential that has
  been tested against alloy data.
- **No claim about transferability.** The pure-element fits reproduce elastic
  constants they were fitted to; nothing follows about a mixed environment.
- **`pair_style ugur` accepts more than one element.** An earlier version of
  this section said it refuses them "until the multi-species kernel is
  validated"; that validation was done — the twin test agrees to 1e-11 over six
  elements and three structures, and the deliberate-corruption half proves the
  unlike entries are really being read — and the refusal was removed. The
  sentence stayed behind and was wrong. **The kernel being correct is not the
  same as an alloy being right**, and the next point is the one that matters.

## The measured failure, stated plainly

`make_alloy_file.py` will generate an alloy file for any pair of elements in
the library, and `pair_style ugur` will run it. **The one case that was
measured against alloy data came out mechanically unstable**: Ni₃Al returns
C₁₁ = −19.9 GPa. The Ni–Al and Au–Cu hold-out tests both failed.

That is a property of the default mixing rules, which are the conventional
guesses tabulated above and have never been fitted to anything. So the path
from this repository to a working alloy potential does not exist yet: what
exists is a validated kernel, a file format, and a generator whose output is
a starting point for a fit that has not been done.

Treat an alloy file from `make_alloy_file.py` as an initial guess, never as a
prediction, and run the elastic constants before anything else — the instability
above was one command away from being missed.

### A second system, and it separates two things the sentence above runs together

Ni₃Al's C₁₁ = −19.9 GPa says the mixing rule fails. It does not say *where*.
A second measurement does, on a system with nothing in common with Ni–Al: the
equimolar refractory alloy **MoNbTaVW**, whose five elements are exactly the
five bcc elements this library force-matched, against 4,417 DFT configurations
from Byggmästar, Nordlund and Djurabekova (*Phys. Rev. B* **104**, 104101,
2021; data doi:10.23729/1e6d0215-d26b-4f5f-8f5b-df575efa6594). χ = 1 — the
rules exactly as tabulated above, not one free parameter.

**On forces the rule costs nothing.** Matched on each frame's own DFT force
magnitude, so that quiet frames are not compared against violent ones, the
relative force error of the five-element alloy equals that of the pure elements
the campaign actually fitted — 0.91× the unary control. Bin by bin, its
correlation with the reference forces runs 0.68–0.81 against the control's
0.66, then 0.82–0.85 against 0.78, then 0.87–0.89 against 0.88: above the
control in the two lower bins and level with it in the highest. The ~50 % relative force error is
the *form's* own baseline, present with one element and unchanged with five.

**On structure it fails completely.** Relax an equimolar random bcc solid
solution and let the cell change shape, and it leaves bcc: the cubic cell goes
to an axis ratio of **1.41 = √2, the Bain path**, common-neighbour analysis
reports fcc and no bcc, and the energy drops 0.107 eV/atom. Every seed tested
does the same. Real MoNbTaVW is single-phase bcc — O. N. Senkov, G. B. Wilks, J. M. Scott and D. B. Miracle, *Mechanical properties of Nb25Mo25Ta25W25 and V20Nb20Mo20Ta20W20 refractory high entropy alloys*, **Intermetallics 19**, 698 (2011) — so this is a
verdict on the potential. **All five constituent elements, put through the
byte-identical relaxation, stay perfectly bcc** — cell ratio 1.0000, CNA 100 %
bcc — so it is the mixing that produces it, not the individual fits.

**Everything cheap passed.** The lattice parameter comes out within 0.3 % of the
3.216 Å reported for this alloy — a figure taken from the secondary literature
and **not yet checked against Senkov's own page**, so it is quoted as reported
rather than as measured; the MD screen holds at 300 K in every seed. Both are
computed with the cell shape held fixed, and a Bain transformation *is* a
change of cell shape — so neither could have seen it. Passing them is a
constraint, not evidence.

**What this changes about the plan.** The failure is not in the cross
parameters, so fitting them attacks a term that is not costing anything: the
gradients are already as good as the form manages anywhere. What the form gets
wrong is which structure sits at the minimum, and forces do not determine that.
Measure the relaxed cell shape — not just the elastic constants — before
trusting any alloy file from this generator, and read a good force score as
evidence about gradients only.

### What is in this release and what is not

Shipped, because it is validated and a reader has to be able to check it:
`pair_ugur.cpp` and `pair_ugur.h` (multi-species), `alloy.py` (the independent
implementation the kernel is checked against), `make_alloy_file.py`, and
`validate_alloy.py` / `validate_alloy_ref.py`, which are the twin and
corruption tests themselves.

Held back: the layer that fits an alloy to alloy data — `alloy_fit.py`,
`alloy_geom.py`, `alloy_holdout.py`, `fetch_alloy_ref.py` and the Ni–Al EAM
comparisons. That work is unfinished and its only measured outcome is the
instability above. It will appear in a later release when there is a result
to go with it. Nothing shipped here depends on it.

## What the mixing rule actually gets wrong (measured 2026-08-29)

The section above says the rules are unvalidated. They have now been measured,
against the 17 ordered compounds in `standalone/alloy_ref.json` — nine Al–Ni,
five Au–Cu, three Cu–Ni — and the result is worse and more specific than
"unvalidated".

**With no correction, all 17 formation energies come out positive.** Every
Al–Ni intermetallic that DFT puts at −0.27 to −0.66 eV/atom is predicted not to
form. That is a sign error, not a calibration error, and it is the Ni₃Al
C₁₁ = −19.9 GPa failure seen from a second direction.

The test that says why: put one multiplicative correction on the cross pair
well depth, `chi_D`, and solve it separately for each compound. If the rule had
the right shape the compounds would agree on one number.

**One knob is not enough, and the way it fails is the finding.** The `chi_D`
that fixes a compound's formation energy makes its *volume* worse — Al–Ni ends
19 to 37 % too large. Deepening a bond should contract a structure. It expands
because `r0_AB` is the arithmetic mean, that mean is too long, and a stronger
bond obeys the wrong length more strictly.

**A confound had to be removed before that could be believed.** The pure
elements here are fitted to *experimental* lattice constants; the compound
reference cells are PBE. The two rulers disagree by −1.4 % (Ni) to +2.3 % (Au)
in nearest-neighbour distance, so demanding that the potential reproduce the
DFT cell exactly makes the free length knob absorb the difference. The second
condition is therefore the **excess** volume, measured against the elements on
each ruler — the same construction the formation energy already uses. Correcting
this moved Cu–Ni's `chi_r0` from 0.960 to 0.993 and tightened its `chi_D`
spread from 7.4 % to 0.4 %, so it was not a detail.

With that done, solving both knobs per compound against formation energy and
excess volume gives:

| system | `chi_D` | spread | `chi_r0` | spread | compounds |
| --- | --- | --- | --- | --- | --- |
| Cu–Ni | 1.029 | **0.4 %** | 0.993 | 2.4 % | 3 |
| Au–Cu | 1.180 | 8.2 % | 0.925 | 9.2 % | 5 |
| Al–Ni | 1.324 | **23.9 %** | 0.926 | 5.0 % | 9 |

**The spread is the result, not the value.** Cu–Ni is genuinely described by
these numbers — 0.4 % over its compounds — and its `chi_r0` is 1 to within
0.7 %, so it is a one-parameter correction: the cross bond is 2.9 % deeper than
the geometric mean and its length rule was right all along. Al–Ni is not
described by them: a 24 % spread in `chi_D` means no single depth serves its
nine compounds, and its cross terms have to be fitted rather than corrected.

**The stability gate.** AlNi (B2) uncorrected has C₁₁ = 84.2 below
C₁₂ = 118.5 — Born-unstable. Corrected it returns 299/180/121, which is stable
and badly overshoots DFT's 204/134/113. **Stability is bought; accuracy is
not.** For Cu–Ni both L₁₂ prototypes are Born-stable with and without the
correction, and the correction barely moves the elastic constants
(183/125/84 → 180/129/84).

Cu–Ni gets no comparison against DFT elastic constants because the reference
carries none for it, and its three stored cells are not the L₁₂ prototypes they
look like — a rebuild check caught a 21.5 meV mismatch. All three sit 41 to
59 meV/atom above the hull, which is the correct physics: **Cu–Ni has no stable
ordered compound**, so what is stored are hypothetical orderings.

`--correct` applies the table above where it exists. It is off by default,
because an uncorrected file is what the mixing rules say and a generator that
silently applied a fitted number would make the two impossible to tell apart
afterwards. The measurement is `refit/alloy_chi3.py`, with `alloy_chi.py`,
`alloy_chi2.py`, `alloy_chi_stab.py` and `alloy_cuni.py` beside it.

**What this does not do.** It corrects two numbers per system against formation
energies and volumes. It does not touch the three-body cross term, it has not
been tested on any property outside the ones it was solved against, and for
Al–Ni it does not work at all. The path to a real alloy potential is still a
DFT force database — but the target is now narrower: the length rule is
adequate and the **energy scale** is what has to be learned.
