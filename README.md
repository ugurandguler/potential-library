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

`docs/index.html` is the published copy. The file it is made from is
`standalone/potential.html`, which `make_gui.py` writes — and **two trees write
a file of that name**, so the page names itself in its own footer, "screen
tree" or "published tree". As of 1.1.0 the two carry the same panels — the
screen tree's renderer and library were promoted into this one for the release
— so the stamp now says only which tree built the file you have open. It is
still the first line to read when a change seems missing: twice a change has
been reported as absent when the file on disk already had it and the browser
was showing an older copy.

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
| surface energies | **2.9× DFT**; of the 38 records whose facet ordering can be decided at all, 3 come out right |
| intrinsic stacking fault | negative in 45 of the 50 records where it is defined |
| thermal expansion | 31 % low; 9 records contract on heating |

**So: do not use these parameters for defect energies, surface energies,
diffusion barriers or melting.** Use them for elastic and vibrational
properties near the fitted volume.

The surface row carries a criterion, not just a count. Most disagreements about
facet ordering are near-ties: where the closest two faces differ by less than
five per cent neither this potential nor the reference resolves them, so the
verdict is three-way — right, wrong, **undecidable** — with the threshold in
`standalone/add_surface.py`. Of our 76 records, 38 are decidable and 3 of those
are right. The same change cost the published potentials more than it cost us,
because the close-packed-first rule had been flattering them too: they go from
48 of 51 to 22 of 32. The conclusion survives the threshold — at two per cent
it is 27/39 against 5/54, at ten per cent 6/7 against 3/14.

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
(2023)).

That is no longer only a citation. LAMMPS ships `Nb.uf3`, which is that
potential for niobium, and `lammps/uf3_compare.py` runs it and ours through the
same elastic and surface machinery with only the potential differing. On
niobium surface energies — which neither was fitted to — the mean deviation is
**277 % for ours against 10 % for the UF3**, and the ordering goes the same
way: the reference gives 110 < 100 < 111, the UF3 reproduces it exactly, ours
gives 100 < 111 < 110. The form reaches both the magnitude and the ordering and
this parameterisation reaches neither. `Nb.uf3` carries an empty CITATION
field, so its training set is unknown and this is not demonstrated
extrapolation.

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
doc/          make_doc.py and make_doc_tr.py, which generate the English and
              Turkish notes from facts.json.  The .docx files are output, not
              source - edit the generators.
refit/        the same nine parameters refitted to density-functional forces.
              Written up separately in REFIT.md.
```

A second working tree, `ugurpotential_screen/`, sits beside this one. It is
where the next release is built, and everything it held for 1.1.0 — the
shell-gap re-cut candidate arms, the (11̄0) polar panel, the reference-path
label fix and the finite-temperature dispersion study — was promoted into this
tree for the release.

**Promotion is not a copy.** `make_gui.py` and `library.json` normally do not
travel: the renderers differ on purpose, and copying the published one over the
screen one destroyed a day's work once. Carrying them the other way, for a
release, is safe only after proving the source is a superset of the
destination, and the proof is not the file size. For 1.1.0 a key-PATH diff of
the two `library.json` files — not a key-name diff — found exactly one field
the screen tree lacked, `dyn/path`, present on all 38 elements and nested
inside `dyn` where a shallow comparison could not see it. It was merged back,
and the result was checked to carry every path either file had and to leave
every fitted parameter and every target byte-identical. Four producer modules
were likewise kept from this tree rather than overwritten, because the release
work had moved them ahead of the screen copies.

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

One thing to know before running `refresh.py`: it **overwrites
`standalone/potential.html` with what your checkout can compute**, which is
less than the shipped page. This was measured, not estimated — the repository
was cloned into an empty directory and run. All 129 Python files compile,
`selftest.py` passes every check, and `refresh.py` completes its fifteen steps
in about a hundred seconds. The page it writes is about 3.3 MB against the
shipped 7 MB, and what is missing is data that is not ours to redistribute or
is too large to carry:

| missing from a rebuild | why |
| --- | --- |
| Materials Project, MC3D, JARVIS overlays | fetched data. Set `MP_API_KEY` for the first; the other two need their own downloads |
| the switched arms `tap` and `tap_ug` | written by `add_taper_overlay.py` from the `dense_*.json` search output, which is gigabytes and is not carried |
| the re-cut candidate arms | `add_candidates.py` reads `refit/`, which is not published |
| *nothing else* | the finite-temperature panel **is** rebuilt: `finiteT.json` ships and the page reads it, so all 24 rows and 24 curves appear. Only regenerating that file needs the switched arm, and `make_finiteT.py` says so rather than failing obscurely |

Every one of those steps says what it could not do and continues, rather than
failing. **The shipped `docs/index.html` is the complete page** and carries all
of it, including the finite-temperature table as `finiteT.json`; keep a copy
before rebuilding if you want it back.

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
III/29a for elastic constants, Brewer LBL-3720 Rev. (1977) for cohesive
energies — and are cited in `refdata.py` at the point of use.

The **measured phonon** comparison outgrew a single compilation. It began as
III/13a's zone-boundary points and is now a set of dispersion curves along the
symmetry path: **35 of the 38 elements carry one** — 29 from tabulated
measured frequencies, 6 reconstructed from published force constants — with
only iridium, rhenium and ruthenium still without. Seven reconstructions are
drawn in all, because tungsten has both: Chen and Brockhouse tabulate a dozen
of its points and plot the rest, so it is counted among the 29 measured and
carries their force constants for the remainder of the path. Each is
transcribed from its own source's table and
cited there; a reconstruction is drawn as a dashed line, never with the open
circles used for individually measured points, and carries the reason its
paper published a fit instead of frequencies. Three of them needed a model
form of their own — a Fourier interplanar series for molybdenum, an axially
symmetric model for tungsten, a dipolar fluctuation model for cobalt — and
each was tested against measured curves already on the page before it was
trusted. The compilations themselves are publisher copyright and are not
redistributed here.

The three anchors are not at one temperature: cohesive energies are 0 K,
lattice constants about 293 K, elastic constants about 300 K, and zero-point
energy is not subtracted anywhere. The size of that inconsistency is measured
element by element in the header of `refdata.py`. Nothing
has been changed on its account, because changing a target means refitting the
library, and that is a decision about scope rather than a correction.

### The dispersion at the temperature it was measured

Every phonon comparison above scores a 0 K harmonic calculation against a
measurement made at 9 to 296 K, which charges the potential for a temperature
difference it was never asked to reproduce. LAMMPS' `fix phonon` can measure
the spectrum a molecular-dynamics trajectory actually has, so it was run for
each element at its own paper's temperature. **24 of the 38 carry that
measurement**; each element page shows the curve beside the 0 K one and the
neutron points, and the 14 without it carry a sentence saying which kind of
absence it is.

Two things make the table readable. The numbers are the **switched** arm
(`<El>_taper.ugur`), not the hard-truncated root the parameter selector calls
MAU — a hard cut leaves the pair energy discontinuous and copper drifts
350 meV per atom per nanosecond under it, so the switched arm is the only one
that can be run at temperature at all. And every run is scored twice, at q and
at a symmetry image of q, so it reports **its own noise floor**; the verdict is
computed from that rather than left to the reader, because 12 of the 24
differences are smaller than it.

The result is **12 gains, no losses, and 12 differences too small to read**.
The gains are mostly one to three points, which is the size the effect ought to
be: copper's real 49 K to 298 K shift is 0.5 per cent on average, and two
independent 1967 neutron measurements of copper disagree by the same 0.5 per
cent. Getting the temperature right is worth about as much as the disagreement
between two laboratories, and a result claiming more than that has another
explanation.

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
