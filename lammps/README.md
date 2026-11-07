# The library in LAMMPS: `pair_style ugur` and `ugur/ang`

`pair_ugur.cpp`, `pair_ugur.h` and `ugurpot.h` implement the two- plus three-body form of
İ. Akgün and G. Uğur (Phys. Rev. B 51, 3458 (1995)) with the refitted parameters of this library:
`ugur` for the MAU form and `ugur/ang` for the UG form with the angular factor. Every file in
`potentials/` is read by one of the two.

---

## 1. Building the pair style

Nothing in `potentials/` can be read by a stock LAMMPS: `pair_style ugur` is
not a LAMMPS package, it is the three files here. Copy them into the source
tree and rebuild.

```
cp pair_ugur.cpp pair_ugur.h ugurpot.h  <lammps>/src/
cd <lammps>/src && make serial          # or: make mpi, or your usual target
```

CMake builds pick the new files up on the next configure:

```
cd <lammps>/build && cmake ../cmake && cmake --build .
```

No package needs to be enabled — `make yes-...` is not involved, because these
are plain pair styles and not a package. Check it took:

```
<lammps>/src/lmp_serial -h | grep -i ugur
```

which must list `ugur` and `ugur/ang`. `ugur` is the two- plus three-body form;
`ugur/ang` adds the angular factor and reads the `.ugur.ang` files. A file's
extension says which style it needs, and the header of every file in
`potentials/` repeats the two lines to paste into an input script.

**Python side**, for the analysis and the page rather than for MD: Python 3
and **numpy**. There is no scipy, matplotlib, pandas or requests anywhere in
this distribution. Four packages are optional and none is needed to use the
library: `esprima` (`make_gui.py` parses the page's JavaScript with it before
writing, and says so rather than skipping silently if it is missing), `mp_api`
and `jarvis` (only the `fetch_*.py` downloaders, whose output is already
carried here), and `fitz` (only `provenance/`, for publisher PDFs that are
deliberately not shipped).

---

## 2. A first run

The two lines every potential file's header gives, on a small copper crystal:

```
units metal
boundary p p p
lattice fcc 3.615
region b block 0 4 0 4 0 4
create_box 1 b
create_atoms 1 box
mass 1 63.546
pair_style ugur
pair_coeff * * potentials/Cu_taper.ugur Cu
velocity all create 600 4928 mom yes rot yes
fix 1 all nve
thermo 200
run 1000
```

With a stock LAMMPS (the `stable` branch) built as above, the 256 atoms run at constant energy and
the total energy stays at −872.76 eV to within a few meV over the 1 ps (1000 steps at the
default 1 fs) (measured 2026-11-03 on a
fresh clone of this repository). For the angular form, write `pair_style ugur/ang` and read
`Cu_taper.ugur.ang` instead.

---

## 3. Which file to use

`potentials/README` explains every file stem, and the `README.md` at the top of the repository has
the table of parameter sets (hard MAU, switched MAU, re-cut, ...). In short: the hard-cut files
(`El.ugur`, `El.ugur.ang`) are for static properties and lattice dynamics, not for molecular
dynamics; for MD use the file the element's page recommends - the re-cut `El_recut.ugur` for 15
elements, the switched `El_taper.ugur` for most others - and read its header, which says what it
was screened for and where it fails.

---

## 4. History

`README_1995_table.md` documents the first step of this project: the two-body part of the
published 1995 Pd–Fe potential as a `pair_style table` (`AkgunUgur1995.table` and its inputs),
with the correction of the printed D by a factor of five. It predates `pair_style ugur`, and its
section on the three-body term is superseded.
