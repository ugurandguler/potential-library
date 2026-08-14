# Provenance — how the reference numbers were obtained

These scripts are not part of the pipeline and **will not run from a clean
checkout**. They read source PDFs that are publisher copyright and are therefore
not in this repository. They are kept because they record exactly where each
experimental target came from, which the numbers themselves cannot.

| script | source | what it produced |
| --- | --- | --- |
| `lb29_extract.py` | Landolt-Börnstein III/29a, Tables 3 and 11 | elastic constants, cubic and hexagonal |
| `lb_extract.py` | Landolt-Börnstein III/13a | measured phonon frequencies at zone-boundary points |
| `brewer_extract.py` | Brewer, LBL-3720 Rev. (1977), Table I | cohesive energies |
| `fetch_mfp.py` | materialsfingerprint.com | **rejected** — see below |

Everything they produced now lives in `standalone/refdata.py` and
`standalone/refdata_phonon.py`, with the source cited at the point of use.

## Why these are worth keeping

Reading tables out of a scanned volume is not a clean operation and it went
wrong repeatedly. Every one of the following was caught by rendering the page as
an image and reading it back, not by the extraction code:

- cobalt came out as 128 instead of 242, from flat text pairing the wrong
  columns;
- lead as 14.8 instead of 48.8, because the source prints `- 43.0` with a space
  after the minus sign;
- molybdenum was dropped entirely, because its S44 sat on the following line;
- rows giving the standard deviation of several determinations were read as if
  they were determinations;
- Brewer's silver uncertainty was parsed as ±2.0 instead of ±0.2, because the
  OCR renders `±0.2` as `Ot o. 2` and the pattern ate the leading zero;
- six Brewer rows were missed entirely because aluminium OCRs as `A1` and
  vanadium in lower case.

Anyone re-deriving these numbers should expect the same class of problem. The
scripts are the record of which values were machine-read, which were read by eye
and marked as such, and which were corrected by hand.

## The rejected source

`fetch_mfp.py` targets materialsfingerprint.com, which publishes phonon
dispersions. It is kept as a record of a negative result rather than as a tool:
of 25 entries checked, 15 carried imaginary modes as deep as −966 cm⁻¹ for
elements where Materials Cloud MC3D returns exactly zero. The data were not
used.

## The one exception to the compilations

Beryllium does not use III/29a. That volume gives C13 = 6 ± 9 GPa, an
uncertainty of 150 %; resonant ultrasound spectroscopy gives 14, and the
bulk-modulus identity settles it — against a measured 116.8 GPa the hexagonal
Voigt combination returns 117.1 for the Migliori values and 111.7 for III/29a.
See `standalone/refdata.py`, which explains this at the point of use.
