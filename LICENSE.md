# Licensing

Three different things are bundled here and they are not all ours to license.
They are kept apart deliberately.

---

## 1. The code — **GPL-2.0**

`standalone/`, `angular/` and `lammps/` are licensed under the **GNU General
Public License, version 2**. The full text is in [`COPYING`](COPYING).

The choice is forced at one point and taken everywhere for consistency:
`lammps/pair_ugur.cpp` is a LAMMPS pair style, LAMMPS is GPL-2.0, and a pair
style distributed for it has to be GPL-2.0 as well. Splitting the repository —
a permissive licence for the fitting pipeline and GPL for the kernel — would
have left every file needing a note saying which one it was under. One licence
across the tree costs a little freedom and removes that question entirely.

    Copyright (C) 2026  Gökay Uğur, Şule Uğur, Melek Güler and Emre Güler
    Gazi University and Ankara Hacı Bayram Veli University

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

## 2. The fitted parameters and the library page

`standalone/fit.json`, `lammps/potentials/*.ugur`, `docs/index.html` and the
derived tables are the scientific contribution of this work, and they are data
rather than program text. **CC-BY 4.0** applies to them: reuse and
redistribution, including modification, on condition of attribution — which is
what the citation request in the README asks for anyway.

## 3. Reference data — **not ours**

The measured quantities used as fitting targets belong to their original
sources. Individual physical constants are facts and are cited normally at the
point of use in `refdata.py`; the compilations they come from are publisher
copyright and are **not** redistributed in this repository.

| data | source | terms |
| --- | --- | --- |
| elastic constants | Landolt-Börnstein III/29a (Springer, 1992) | publisher copyright; values cited, volume not redistributed |
| measured phonons | Landolt-Börnstein III/13a (Springer, 1981) | as above |
| cohesive energies | Brewer, LBL-3720 Rev. (1977) | US national laboratory report |
| beryllium elastic constants | Migliori *et al.*, J. Appl. Phys. **95**, 2436 (2004) | cited |
| surface energies (experiment) | Tyson and Miller, Surf. Sci. **62**, 267 (1977) | cited |
| tungsten DFT energies and forces | Szlachta, Bartók and Csányi, Phys. Rev. B **90**, 104108 (2014) | used in a refit reported separately; dataset not redistributed here |

Calculated values drawn for comparison come from third parties and keep their
own terms:

| source | used for | terms |
| --- | --- | --- |
| Materials Project | phonons, surface energies | CC-BY 4.0; attribution required |
| Materials Cloud MC3D | phonons | **CC-BY 4.0**, from the record's own licence (Materials Cloud Archive, doi:10.24435/materialscloud:rw-t0); attribution required |
| JARVIS-DFT (NIST) | phonons, force-field properties | US Government work; citation requested |
| NIST Interatomic Potentials Repository | independent check of our baseline numbers | US Government work; citation requested |

None of these enters the fit. All of them are drawn for comparison only, and
where one conflicts with a measurement by more than 25 % it is withheld from
the plot rather than shown.

### LAMMPS

`lammps/pair_ugur.cpp` and `pair_ugur.h` are LAMMPS pair styles. Their
file-reading and parameter-indexing structure follows LAMMPS's own Tersoff
style (LAMMPS's own `src/MANYBODY/pair_tersoff.cpp`, by Aidan Thompson,
SNL - a file in the LAMMPS distribution, not in this one), and both carry the
LAMMPS copyright header for that reason. LAMMPS is GPL-2.0 and
Copyright (2003) Sandia Corporation; this is why the code here is GPL-2.0 and
not something more permissive.

### Published potentials used as baselines are not redistributed either

Fifty-one published EAM, MEAM and ADP potentials were run through this
project's own code to produce the baseline columns. Those potential files stay
out of this repository, and `lammps/potentials/` carries only the `.ugur` and
`.ugur.ang` sets that are ours. Fetch the baselines from the NIST Interatomic
Potentials Repository or the LAMMPS distribution, where each one is identified
by more than its filename — which matters, because `Cu_zhou.eam.alloy` as
LAMMPS ships it is a different potential from every implementation NIST hosts
under that author and year, and it took a surface-energy comparison to notice.

---

## A note on the reference data

Systematically extracting a substantial part of a commercial compilation and
republishing it is a different act from citing individual measured values in a
paper, even though the values themselves are facts. This repository does the
latter: the numbers appear as fitting targets alongside the fitted results, in
the way any paper's table would present them. If the intention ever changes —
for instance publishing a downloadable replica of the tables — that is worth a
question to the university library or research office first.
