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

**Before anything else, if you want to run this in LAMMPS:** `pair_style ugur`
is not a LAMMPS package. Copy `lammps/pair_ugur.cpp`, `lammps/pair_ugur.h` and
`lammps/ugurpot.h` into your LAMMPS `src/` and rebuild — `lammps/README.md`
section 0 has the two commands and the one-line check that it took. For the
Python side you need Python 3 and **numpy**, and that is the whole requirement
for the analysis and for rebuilding the page. Three optional packages appear
and none is needed to use the library: `esprima`, which `make_gui.py` uses to
parse the page's own JavaScript before writing it and skips with a message if
absent; `mp_api` and `jarvis`, used only by the `fetch_*.py` scripts that
download reference data already carried here; and `fitz`, only in
`provenance/`, for publisher PDFs that are deliberately not shipped.

The library reproduces **the elastic tensor at the experimental lattice
constant** as well as tabulated EAM does, and that is what it was fitted to. It
is a statement about the second derivative of the energy at fixed coordination.

It does **not** transfer to coordination changes. Measured against
density-functional theory and against 51 published potentials run through the
identical code:

| quantity | how it does |
| --- | --- |
| elastic constants C_ij | median RMS 1.45 % (UG), 5.83 % (MAU) |
| ground state of the model | the fitted structure is the lowest of bcc, fcc and hcp in **0 of 38** switched MAU records and 4 of 38 switched UG; the structure the model prefers sits a median 29 meV/atom lower (chromium 412). The re-cut arm is the exception, 21 of 23, and its selection rule screens for it. Nineteen published potentials run through the same code: 18 right. The hard-cut arms cannot be asked — their pair term does not vanish at the cutoff and the three structures hold different neighbour counts inside it, so the differences reach 19 eV/atom (iridium) and are truncation, not energy. `lammps/struct_rank.py`, in `library.json` as `<arm>.ground` |
| vacancy formation energy | **2.25× the median of the published classical potentials** in JARVIS-FF (not DFT) for the switched arms, over the 21 elements that have one, and above it in 18 of the 21. The hard-cut reference sets are closer on the median (1.05 MAU, 1.55 UG) and **turn negative**: 13 elements in MAU, 17 in UG. The switched arms give no negative anywhere. `lammps/vacancy.py`, all 38 elements and all four sets, in `lammps/vacancy_*.json`; **not carried in `library.json`, so this is the one row not on the page** |
| surface energies | **2.9× DFT**; of the 38 records whose facet ordering can be decided at all, 3 come out right |
| intrinsic stacking fault | negative in 45 of the 50 records where it is defined |
| thermal expansion | 32 % low; 9 records contract on heating |

**So: do not use these parameters for defect energies, surface energies,
diffusion barriers or melting.** Use them for elastic and vibrational
properties near the fitted volume.

**And check the finite-temperature panel before running one warm.** Passing the
0 K screen and holding the crystal in the molecular-dynamics screen is not the
same as staying stable on heating: the switched records of Cs, Fe, K, Mo, Rb, V
and W lose C′ = (C₁₁ − C₁₂)/2 between 5 and 10 per cent of their melting
point, measured with the Born stress-fluctuation method. The page marks those
points hollow and now says so in the panel note.

The surface row carries a criterion, not just a count. Most disagreements about
facet ordering are near-ties: where the closest two faces differ by less than
five per cent neither this potential nor the reference resolves them, so the
verdict is three-way — right, wrong, **undecidable** — with the threshold in
`standalone/add_surface.py`. Of the 76 records in the two shipped arms, 38 are
decidable and 3 of those are right. **That is the count for `tap` and `tap_ug`;
the producer now summarises every arm** and prints 121 records, 63 decidable
and 9 right, because the re-cut arms carry surface runs too. The row above is
the shipped library's figure and the producer's is the wider one — the same
measurement over a larger set, not a disagreement. The same change cost the published potentials more than it cost us,
because the close-packed-first rule had been flattering them too: they go from
48 of 51 to 22 of 32. The conclusion survives the threshold — at two per cent
it is 27/39 against 5/54, at ten per cent 6/7 against 3/14.

Every row above is a measurement, not an estimate, and the evidence is in this
repository: open `docs/index.html` and pick an element. Each page carries the
vacancy, the surfaces facet by facet, the stacking fault and the thermal
expansion, with density-functional values and 51 published EAM and MEAM
potentials — run through this same code — beside them. The drivers that
produced those numbers are in `lammps/`.

Two limits were established by scanning the parameter space rather than
inferred from one fit. Cd and Zn's axial anisotropy has no solution in any of
the three forms as parameterised here. The C₄₄/C′ floor that keeps eight cubic
metals out of reach of the MAU form is a property of the **hard cut-off**, not of
the form (the angular factor already lowers it for Al, Cr, Nb, V and W): with
the switched cut-off all eight are within reach, and aluminium, iron, niobium,
tantalum and vanadium fit exactly — but only aluminium and tantalum keep that
fit through the later checks (iron and niobium leave the bcc lattice when the
atoms are displaced and relaxed, and vanadium's C′ turns negative at 5 % of its
melting point). Chromium, molybdenum and tungsten would need a repulsive first
neighbour shell, which the fit forbids. The rest are
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

## Force matching: what it repairs, and what it does not

The rows above are properties of the **shipped** records, which are fitted to
experimental anchors — the cohesive energy, the lattice constant, the bulk
modulus and the elastic constants — and to nothing else. A separate line of
work asks what happens if the same functional form is fitted to **DFT forces**
while those anchors are held. Eleven such fits exist, all converged —
the last two, tantalum and vanadium, landed on 2026-09-08. **They are
not in the shipped library**; they are an experiment about the form, reported
here because two of the caveats above turn out to be reachable.

| fit | element | E_vac | E_int | stacking fault | dispersion |
| --- | --- | --- | --- | --- | --- |
| unconstrained | Nb | **−1.973** | **−2.808** | — | 21.2 % |
| constrained | Nb | **+3.674** | **+2.012** | — | 12.6 % |
| constrained | Cu | 1.847 | 3.137 | −83.4 | 5.4 % |
| + stacking hinge | Cu | 1.859 | 3.153 | −15.0 | 5.3 % |
| constrained | Ni | 2.299 | 3.798 | −239.5 | 7.7 % |
| + stacking hinge | Ni | 2.368 | 3.939 | **+39.9** | 7.4 % |

Three things are worth reading off it.

**A negative vacancy energy travels with a negative interstitial energy, and
one constraint repairs both.** Niobium's unconstrained fit prefers a hole *and*
prefers an extra atom crammed in. Bounding the three-body term against the pair
term — the same ceiling the shipped fitting script already rejects on — turns
both positive at once. The interstitial had never been computed anywhere in
this project, while the training data carries thirty-two interstitial
configurations; the fit was being shown them and nothing was scoring the
result.

**The stacking fault can be given the right sign, and it is nearly free.** A
fifth penalty asking for the hexagonal stacking to sit above the cubic one
takes nickel's fault from −239.5 to +39.9 mJ/m² against a reference of +125 —
the first record in this work with the right sign — while the cohesive energy,
the elastic constants, the vacancy and the interstitial stay where they were
and the dispersion gets slightly better. The sign survives at 300 K, which had
to be checked separately: the molecular-dynamics screen runs a *perfect*
crystal and so contains no fault to test.

**The dispersion improves everywhere it can be checked**, from 21.2 to 12.6 per
cent for niobium, 13.4 to 7.4 for nickel and 10.0 to 5.3 for copper, on curves
that are in no objective. What does not improve is the facet ordering of the
body-centred fits, which stays wrong.

Seven of them ship, as the `tap_force` records (Cu, Mo, Nb, Ni, Ta, V, W).
Three of those - niobium, tantalum and vanadium - sit strictly inside every
constraint; the other four sit on the three-body ceiling E3/E2 = 0.30 and
exceed it by at most 0.0002 in that ratio, which the shipped fitting script
rejects on with no tolerance. That is a decision still to be taken, and
the ceiling should not be moved to accommodate results produced under it.

## The melting point, measured for the first time

The caveat above says not to use these parameters for melting. That is now a
number rather than a caution. Tantalum's force-matched record melts at
**2625 ± 125 K against an experimental 3290** — about twenty per cent low, at
zero pressure. (Release 1.4.3 said 2250 ± 250 K; see below for why that moved.)

It was measured wrongly first, and the failure is worth recording because it
produced a plausible number. A two-phase coexistence run returned 3438 K, which
is void: there was never a solid phase in the cell. Velocities created at twice
the target temperature melt everything before the two-phase construction
begins, and the run's own solid-temperature diagnostic cannot see it because
the group it measures is fixed by *initial position* — two halves of a uniform
liquid report two indistinguishable temperatures. Every melting run here is now
gated on the static structure factor per half-cell, printed beside a
perfect-crystal control that must return 1.000.

**The protocol is in the repository.** `lammps/in.coexistence_bracket` builds
the solid at 1200 K — never at the melting point, where it nucleates at random
— melts half of it at 6000 K, brings **both** halves to a trial temperature and
holds them there. `lammps/phase_gate.py` then reads which phase won:

```
lmp -in in.coexistence_bracket -var TTRY 2500
python phase_gate.py *.dump --halves --data perfect.data --a 3.3052
```

| trial T | S per half after 120 ps | what happened | conclusion |
| --- | --- | --- | --- |
| 2000 K | 0.760 / 0.802 | the liquid half froze | Tm above 2000 |
| 2250 K | 0.847 / 0.857 | the liquid half froze | Tm above 2250 |
| 2500 K | 0.751 / 0.708 | the liquid half froze | Tm above 2500 |
| 2750 K | 0.068 / 0.098 | the solid half melted | Tm below 2750 |

**What changed the number, and what is not known.** The copy of
`in.coexistence_bracket` in 1.4.3 and 1.4.4 held the *volume* instead of the
pressure in its last three stages, which leaves the cell at about 2 GPa near
2500 K and raises the apparent melting point, while the numbers in its own
comment came from zero-pressure runs; it now runs every stage after the first
at zero pressure. A 40 ps hold does not decide: at every trial temperature from
2000 to 2750 K the two halves still differ when it ends, so the default is now
120 ps, which decides all four. The earlier 2250 ± 250 K came from
zero-pressure runs started from a different cell, one of which ended at 2250 K
with both halves disordered; the released file reproduces neither that run nor
its number at either hold time, and why is not known.

**That is the whole measurement**: at a fixed temperature and zero pressure,
whether the solid melts or the liquid freezes *is* the answer, so it needs no
interface-settling stage and no temperature reading — the fragile part, and the part that produced
the void 3438 K. `lammps/in.coexistence` carries that stage for anyone who
wants the settling temperature as well. Both files default to tantalum's
force-matched record and take the element, potential file, lattice constant and
mass as command-line variables.

## The same improvement, twice — and it does not add up

Two separate things reduce the error against a **measured** dispersion. Force
matching moves the 0 K harmonic spectrum. Measuring the spectrum instead by
displacement correlation during equilibrium molecular dynamics, at the
temperature the neutron experiment was actually run at, gives the renormalised
spectrum rather than the harmonic one.

The natural expectation is that they add: one repairs the form, the other
supplies the anharmonicity no harmonic form can carry. **They do not add. They
substitute.**

| | library 0 K | library at T | fit 0 K | fit at T | gain at T |
| --- | --- | --- | --- | --- | --- |
| V\* | 14.5 | 14.7 | 23.6 | **12.1** | **+11.5** |
| Nb | 12.3 | 12.8 | 12.6 | **10.6** | +2.0 |
| Ta | 14.8 | 12.9 | 12.2 | **9.7** | +2.5 |
| Cu | 10.0 | 8.4 | 5.4 | **4.1** | +1.3 |
| Mo | 14.8 | 9.2 | 9.6 | 12.4 | −2.8 |
| Ni | 13.4 | 11.9 | 7.7 | 10.8 | −3.1 |
| W | 17.8 | 10.3 | 8.9 | 14.6 | −5.7 |

Per cent of the highest measured frequency, all against the same curves and
the same q-points. \* vanadium is scored against a fitted model, not a
measurement — see below.

Ordered by how much force matching gained at 0 K, **the gain at temperature
falls steadily and changes sign**, with a correlation of −0.96 over the seven.
The crossover sits near five points. Copper is face-centred and falls between
two body-centred elements in that ordering, so this is not a structure effect.
Whichever lever is pulled first takes the improvement; the other moves back.

Vanadium is the case that makes it plain, because it is the only one at the
far end. Force matching made its 0 K dispersion **9.1 points worse**, and the
thermal measurement recovered 11.5 — ending below the library's own figure. A
fit that looks like a failure at zero kelvin is the best vanadium record here
once it is measured at the temperature its reference was taken at. **A
0 K-only reading would have discarded it.**

What this does **not** say is what the shared content is. Anharmonic
renormalisation and the effect of fitting to strained and defected
configurations are both plausible, and this measurement does not separate
them. It also gives no floor: four body-centred elements suggested a common
limit near 9–10 per cent and copper's 4.1 removed it.

Vanadium's row in this comparison is scored against a Born–von Kármán model
(Colella and Batterman 1970), gated on the elastic constants it was *not*
fitted to. Vanadium scatters neutrons almost entirely incoherently, which is
what makes it the standard neutron calibrant, so no neutron dispersion exists;
a measured one does, by inelastic x-ray scattering (Bosak *et al.* 2008), and
since September 2026 the library's 0 K comparison uses it. Tungsten has twelve measured points and molybdenum fifteen against
niobium's hundred and thirty-eight, so those two rows are thin.

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
              npt_expansion.py, elastic_T.py, and the melting measurement
              (in.coexistence_bracket, in.coexistence, phase_gate.py).
              potentials/ holds our own .ugur files; the published baselines
              are not redistributed.
docs/         index.html, the interactive library.  GitHub Pages serves from
              here.  forms.html beside it is a separate study page: the
              published form against two embedding variants on five bcc
              metals and their alloys, not offered for use (make_forms.py).
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
is a molecular-dynamics kernel. They agree on the energy to 10⁻¹¹ eV/atom and on
the pressure to 10⁻⁶ GPa (`lammps/validate_pair.py`),
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

There is a second thing to know before re-fitting. The admissibility test in
`fit.py` — the ceiling on the three-body share of the cohesive energy — was
tightened *after* the hard-cut library had been fitted, and it was derived from
the switched sets. Re-evaluating a shipped record through the current file
therefore returns nothing for **22 of the 36 hard-cut fits** (the list is in
`fit.py`, measured 2026-08-29). No parameter changed; the acceptance test in
front of them did. Without this note a reader re-fitting has no way to tell a
version difference from a defect.

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
| *nothing else* | the finite-temperature panel **is** rebuilt: `finiteT.json` ships and the page reads it, so all 27 rows and 26 curves appear (rhenium has one measured frequency, not a curve). Only regenerating that file needs the switched arm, and `make_finiteT.py` says so rather than failing obscurely |

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
the analytic code on energy (10⁻¹¹ eV/atom) and pressure (10⁻⁶ GPa). `python selftest.py` and
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
energies — and are cited in `refdata.py` at the point of use. The electronic
heat-capacity coefficients the thermodynamics panel subtracts come from
Kittel's *Introduction to Solid State Physics*, 8th ed., Table 2, and are
cited in `refdata_electronic.py`; lutetium and ytterbium are not in that table
and carry no electronic term.

The **measured phonon** comparison outgrew a single compilation. It began as
III/13a's zone-boundary points and is now a set of dispersion curves along the
symmetry path, and **all 38 elements carry one**: 27 dispersions transcribed
from tabulated measured frequencies; 3 digitised from published figures —
iridium and ruthenium by neutrons, vanadium by x-rays — and flagged as such;
4 records of individual measured frequencies rather than a dispersion
(beryllium, molybdenum, rhenium, tungsten); and 4 reconstructed only from
published force constants (cobalt, chromium, rhodium, strontium). Seven
reconstructions are drawn in all, because molybdenum, tungsten and vanadium
carry one beside their measured points. Each is transcribed from its own
source and cited there; a reconstruction is drawn as a dashed line, a
digitised point as a diamond, and only a tabulated measurement gets the open
circle, and a reconstruction carries the reason its paper published a fit
instead of frequencies. Three of them needed a model
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
each element at its own paper's temperature. **27 of the 38 carry that
measurement**; each element page shows the curve beside the 0 K one and the
neutron points, and the 11 without it carry a sentence saying which kind of
absence it is. Iridium, ruthenium and rhenium joined the original 24 in
September 2026: the first two scored against points digitised from published
figures, rhenium against its single Raman frequency at Γ.

Two things make the table readable. The numbers are the **switched** arm
(`<El>_taper.ugur`), not the hard-truncated root the parameter selector calls
MAU — a hard cut leaves the pair energy discontinuous and copper drifts
176 meV per atom per nanosecond under it (40 ps of NVE from 600 K,
`lammps/nve_check.py`), so the switched arm is the only one
that can be run at temperature at all. And every run is scored twice, at q and
at a symmetry image of q, so it reports **its own noise floor**; the verdict is
computed from that rather than left to the reader, because 15 of the 27
differences are smaller than it.

The result is **12 gains, no losses, and 15 differences too small to read**.
The three newest are all in the last group: iridium 9.9 to 9.4 per cent
against a floor of 1.3, ruthenium 9.0 to 8.4 against 2.9, rhenium 12.0 to 12.7
against 2.4.
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
