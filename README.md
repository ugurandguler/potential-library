# Uğur and Güler interatomic potential library

[Prof. Dr. Gökay Uğur](https://avesis.gazi.edu.tr/gokay) ·
[Prof. Dr. Şule Uğur](https://avesis.gazi.edu.tr/suleugur) — Gazi University
· [Prof. Dr. Melek Güler](https://avesis.hacibayram.edu.tr/melek.guler) ·
[Prof. Dr. Emre Güler](https://avesis.hacibayram.edu.tr/guler.emre) — Ankara
Hacı Bayram Veli University

A two-body plus three-body potential for 38 metals, refitted against
experimental elastic constants, lattice constants and cohesive energies, with
a LAMMPS pair style, a self-contained interactive library page, and a written
record of what the fits do and do not reproduce.

**The interactive library is at `docs/index.html`.** It is a single
self-contained file: no server, no external requests, no scripts loaded from
anywhere. Open it in a browser, or serve it through GitHub Pages.

## Read this before using the parameters

The library reproduces **the elastic tensor at the experimental lattice
constant** as well as tabulated EAM does, and that is what it was fitted to. It
is a statement about the second derivative of the energy at fixed coordination.

It does **not** transfer to coordination changes. Measured against
density-functional theory and against 51 published potentials run through the
identical code:

| quantity | how it does |
| --- | --- |
| elastic constants C_ij | median RMS 1.45 % (UG), 5.83 % (MAU) |
| vacancy formation energy | **2 to 3× too large**; negative for Cr, Mo, W |
| surface energies | **2.9× DFT**, with the facet ordering wrong in 72 of 76 records |
| intrinsic stacking fault | negative in all 25 records where it is defined |
| thermal expansion | 31 % low; 9 records contract on heating |

**So: do not use these parameters for defect energies, surface energies,
diffusion barriers or melting.** Use them for elastic and vibrational
properties near the fitted volume.

Every row above is a measurement, not an estimate, and the evidence is in this
repository: open `docs/index.html` and pick an element. Each page carries the
vacancy, the surfaces facet by facet, the stacking fault and the thermal
expansion, with density-functional values and 51 published EAM and MEAM
potentials — run through this same code — beside them. The drivers that
produced those numbers are in `lammps/`.

Two of these limits are limits of the *functional form*, established by
scanning the parameter space rather than inferred from one fit: Cd and Zn's
axial anisotropy, and the C₄₄/C′ floor for eight cubic metals. The rest are
limits of the **fitting data** rather than of the form. These parameters were
fitted to seven derived numbers of one bulk structure at one volume and never
saw a vacancy or a surface. Published work on the same two- plus three-body
form, with free radial shapes and fitted instead to density-functional energies
and forces on configurations that contain those defects, reproduces them to
within a few per cent (Xie, Rupp and Hennig, *npj Comput. Mater.* **9**, 162
(2023)). Whether that carries over to this parameterisation is under test and
will be reported separately.

## Three potentials, one hierarchy

| | name | parameters | note |
| --- | --- | --- | --- |
| **AU** | Akgün–Uğur | D, C, r₀, α | the published form |
| **MAU** | modified Akgün–Uğur | D, C, r₀, α, α₃ | this work; φ₃ carries its own decay constant |
| **UG** | Uğur–Güler | D, C, r₀, α, α₃, λ₂, λ₄ | angular generalisation; φ₃ sees the bond angle |

AU ⊂ MAU ⊂ UG — each contains the previous one exactly as a special case
(s₃ = 1 recovers AU, λ₂ = λ₄ = 0 recovers MAU), so none of them replaces
another.

## Layout

```
standalone/   the production tree.  latdyn.py computes elastic constants,
              phonons and thermodynamics from the potential's own analytic
              derivatives; fit.py finds the parameters.
angular/      the UG branch: the Legendre factor on phi3, plus angfc.py, which
              supplies force constants for it so the fits can be screened for
              dynamical stability.
lammps/       pair_ugur.cpp, the pair style, and the drivers that measured
              everything above - vacancy.py, surface.py, stacking.py,
              npt_expansion.py, elastic_T.py.  potentials/ holds our own
              .ugur files; the published baselines are not redistributed.
docs/         index.html, the interactive library.  GitHub Pages serves from
              here.
provenance/   how the reference numbers were read out of the source volumes.
              Not runnable from a clean checkout - the volumes are publisher
              copyright and are not in this repository.
```

Two independent implementations of the same potential exist on purpose.
`standalone/latdyn.py` is analytic and lattice-based; `lammps/pair_ugur.cpp`
is a molecular-dynamics kernel. They agree on energy and pressure to 10⁻¹¹,
which is what makes a disagreement anywhere else a real result rather than a
bug in one of them.

## Running it

```
cd standalone
python refresh.py            # rebuild everything downstream of fit.json
python fit.py Cu Pd Mg       # refit selected elements
python selftest.py           # verify latdyn against known limits
```

`refresh.py` runs the whole chain in the one order that works and stops on the
first failure, so the page is never built from half-written data.

`fit.py` has the same property and it is worth knowing before you use it: it
**overwrites `standalone/fit.json`**, and its own search is short. The shipped
parameters came from `dense_fit.py` at 400 restarts, merged in with
`merge_fits.py`; a quick `fit.py` run can land worse. On a test run Cu and Pd
came back exact either way, but magnesium went from 13.1 to 15.0 % RMS — it is
hexagonal, its γ sits on a bound, and it is the kind of element a short search
loses. Use `dense_fit.py <el> <restarts>` when the number matters.

One thing to know before running `refresh.py`: it **overwrites `standalone/potential.html`
with what your checkout can compute**. The third-party comparison overlays —
Materials Project, Materials Cloud, JARVIS — are fetched data that is not
redistributed here, so those steps skip and say so, and the rebuilt page comes
out around half the size of the shipped one with the density-functional
comparisons missing. The shipped `docs/index.html` is the complete version;
keep a copy if you want it back.

For the LAMMPS side, build `pair_ugur.cpp` into LAMMPS as a normal pair style,
then:

```
cd lammps
python validate_pair.py      # pair_style ugur against latdyn.py
python validate_kernel.py    # the shared kernel against latdyn.py
```

## What is checked

Nothing in the chain is trusted without a test that could fail: the acoustic
sum rule to machine precision, the Dulong–Petit limit of the heat capacity, the
frozen and relaxed elastic constants agreeing for one-atom cells, analytic
forces against finite differences (2.6 × 10⁻¹⁰ eV/Å), the numerical angular
force constants against the analytic ones (~10⁻⁸), and the pair style against
the analytic code on energy and pressure (10⁻¹¹). `python selftest.py` and
`lammps/validate_*.py` run them.

The comparisons against published potentials are checked from outside as well:
the NIST Interatomic Potentials Repository publishes its own numbers for the
potentials it hosts, and over the 74 surface facets belonging to potentials
confirmed to be the same file, the median difference from the numbers here is
0.033 %. Identity is established by numerical agreement, not by filename — a
distinction that cost a whole comparison to learn.

## Reference data

Experimental targets are traceable to primary compilations — Landolt-Börnstein
III/29a for elastic constants, III/13a for measured phonons, Brewer LBL-3720
Rev. (1977) for cohesive energies — and are cited in `refdata.py` at the point
of use. The compilations themselves are publisher copyright and are not
redistributed here.

The three anchors are not at one temperature: cohesive energies are 0 K,
lattice constants about 293 K, elastic constants about 300 K, and zero-point
energy is not subtracted anywhere. The size of that inconsistency is measured
element by element in the header of `refdata.py`. Nothing
has been changed on its account, because changing a target means refitting the
library, and that is a decision about scope rather than a correction.

Calculated phonon dispersions from the Materials Project, Materials Cloud MC3D
and JARVIS-DFT are drawn for comparison only; none enters the fit, and each
carries its own citation and licence terms.

## Citation

The functional form is not ours to claim — it is İ. Akgün and G. Uğur,
*Phys. Rev. B* **51**, 3458 (1995); *Nuovo Cimento D* **19**, 779 (1997);
*Nuovo Cimento D* **20**, 1549 (1998). Cite those for the potential.

The parameters here are a fresh fit and are not the published ones. Cite this
repository for them — <https://github.com/ugurandguler/potential-library> — and
please state which release you used. That request is not a formality: the
numbers move as the fits improve, and a reader who cannot tell which release a
result came from cannot reproduce it. `CITATION.cff` carries the machine-
readable form.

## Licence

**GPL-2.0** for the code, forced by the LAMMPS pair style and taken across the
tree for consistency; **CC-BY 4.0** for the fitted parameters and the library
page. The reference data are other people's measurements and are cited as
theirs. See `LICENSE.md` and `COPYING`.
