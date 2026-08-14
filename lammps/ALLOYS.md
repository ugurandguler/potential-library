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
