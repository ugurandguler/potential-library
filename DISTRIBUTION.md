# Getting the potential into molecular dynamics

This began as a planning note written before any of it existed. Most of what it
proposed has since been built and is in this repository, so it has been rewritten
to describe what is here and what is still open, rather than what might be done.

## The forces

Molecular dynamics needs forces, and `latdyn._gradient` computes them
analytically for both terms. Verified against finite differences on a displaced
two-atom Mg cell: the largest disagreement was **2.6 × 10⁻¹⁰ eV/Å**.

The three-body force is unusually simple, and for a reason worth stating:

    g1 = phi3'(r1 + r2)
    F_ja += g1 * n_a
    F_jb += g1 * n_b
    F_i  -= g1 * (n_a + n_b)

`phi3` depends only on `r1 + r2`, so there are no angular derivatives at all.
Stillinger-Weber and Tersoff need pages of algebra here; this needs three lines.
The same property that limits the potential's accuracy — it cannot reach a
C44/C′ ratio below its floor — is what makes it trivial to implement.

## What was built: one kernel, thin wrappers

The rule was not to write a molecular dynamics code. Ninety per cent of that
work is neighbour lists, integrators, thermostats, barostats, boundaries,
restart files and parallelism — none of it specific to this potential, all of it
debugged elsewhere over thirty years. So the physics went into one
dependency-free C kernel with thin wrappers around it:

    lammps/ugurpot.h         the physics, header-only C99, libm and nothing else
      |- lammps/pair_ugur.cpp    pair_style ugur and ugur/ang, thin
      |- standalone/latdyn.py    the independent reference implementation

`latdyn.py` is analytic and lattice-based; `pair_ugur.cpp` is a
molecular-dynamics kernel. They are written from the same equations but share no
code, and they **agree on energy and pressure to 10⁻¹¹** across the library
(`lammps/validate_pair.py`, `lammps/validate_kernel.py`). That is the asset:
most people writing a pair style have nothing to check against, and here a
disagreement anywhere else is a result rather than a bug in one of them.

Potential files for both arms are in `lammps/potentials/` — `.ugur` for MAU,
`.ugur.ang` for UG. A Python or ASE binding over `ugurpot.h` has not been
written; nothing prevents it.

## Licensing

Settled: **GPL-2.0** across the tree, forced by the pair style, because LAMMPS
is GPL-2.0. The fitted parameters and the library page are CC-BY 4.0. See
`LICENSE.md` and `COPYING`; the reasoning for taking one licence across
everything rather than splitting it is written there.

## Worth considering: OpenKIM

[openkim.org](https://openkim.org) is a code-agnostic standard for interatomic
models, supported directly by LAMMPS. Publishing as a KIM model would give:

- a **DOI**, so the model is citable independently of the papers
- use from LAMMPS, ASE and others through the same interface
- **automated verification** — lattice constants, elastic constants, phonon
  stability, run by KIM itself

That last point matters here. Part of the validation chain built by hand in this
repository is something KIM does as standard, and "which elements are
dynamically stable and which are not" would appear on the model page. We already
know that answer and have no interest in hiding it.

This one is genuinely still open.

## Alloys: the format is decided, the fit is not

The multi-species kernel exists and is validated — the twin test agrees to
10⁻¹¹ over six elements and three structures, and its deliberate-corruption half
proves the unlike entries are really read. The file format is LAMMPS's Tersoff
convention, one line per ordered triple.

The question this note once listed as undecided — which parameter set applies to
a triplet with centre A and legs B and C — **has been decided**: φ₃'s radial
shape comes from the centre element's own (A,A,A) entry, because φ₃ depends on
its legs through r₁ + r₂ and does not factorise into per-leg pieces. The
reasoning, and the wrong first version that read it off the triple line, are in
`lammps/ALLOYS.md`.

What is still open is the part that matters: **no alloy has been fitted**. The
default mixing rules are the conventional guesses, and the one case measured
against alloy data came out mechanically unstable (Ni₃Al, C₁₁ = −19.9 GPa). The
fitting layer is held back from this release for that reason. Treat a generated
alloy file as a starting point for a fit, never as a prediction.

## Which elements are fit to use

Reproducing elastic constants does not make a potential safe for molecular
dynamics — a fit can match every C_ij and still have imaginary modes, which in
MD melts or collapses the crystal. Six tests sit between a good RMS and a usable
potential, and two of them apply to nearly everything:

  - **The tetragonal well.** For 30 of 76 records the energy along the
    volume-conserving tetragonal strain turns over within a few per cent of
    zero, over a hump of one or two meV/atom. Where that hump is small the
    elastic constants go Born-unstable at a few per cent of the melting point
    even though the crystal itself sits perfectly still. `lammps/bain.py`, and
    it costs seconds.
  - **The ground state.** 72 of 76 records put a structure other than the
    fitted one lowest — usually hcp, by tens to hundreds of meV/atom. Eighteen
    of nineteen published potentials tested the same way get this right. It
    does not invalidate the fitted tensor, but it does mean the crystal has
    somewhere to go. `lammps/struct_rank.py`.

**Practical rule before any MD run:** check the sign of the vacancy formation
energy, then run `jiggle_test.py`, `bain.py` and `struct_rank.py`. All four are
minutes of work and between them they catch every failure this project found the
expensive way.

A sensible starting set is the fcc noble metals and the alkalis — Cu, Ag, Au,
Pd, Pt, Pb, Ca, Sr, Ba, K, Na — which are both low-RMS and stable, and which
come through both tests in far better shape than the body-centred metals, though
not untouched:

  - Nine of their twenty-two records have a tetragonal well that never turns
    over. The rest turn over only far out, at strains of 0.04 to 0.11 rather
    than vanadium's 0.008, and the weakest hump in the whole set is lead's
    2.7 meV/atom against vanadium's 0.94 and iron's 0.66. None of them shows a
    finite-temperature Born violation.
  - Their ground state is still nominally wrong, by 0.4 to 31.5 meV/atom rather
    than the hundreds the body-centred metals show. For fcc against hcp — two
    close-packed structures that differ only in stacking — a few tens of meV is
    inside what any pair-like form can claim to resolve, and it is the
    body-centred cases, where the competing structure is genuinely different,
    that should be read as a failure.

And the standing limit, which no amount of MD hygiene fixes: these parameters
reproduce the elastic tensor at the fitted volume and do not transfer to
coordination changes. The README's table says by how much.
