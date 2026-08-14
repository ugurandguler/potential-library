# The Akgün–Uğur model potential → LAMMPS

Sources: **İ. Akgün and G. Uğur, Phys. Rev. B 51, 3458 (1995)** (Pd–10%Fe)
and **Il Nuovo Cimento D 19, 779 (1997)** (Fe–28%Pd).

```
phi2(r) = D / (2(m-1) r) * [ beta^m exp(-m*alpha*r) - m*beta*exp(-alpha*r) ],   beta = exp(alpha*r0)
```

---

## 1. IMPORTANT: the D of Table III must be divided by 5

The `D` values as printed, used with the stated unit (10⁻²⁹ J·m) and the
stated formula, do not reproduce the paper's **own** numbers — they inflate
every one of them by exactly a factor of five. The effective value is:

```
D_eff = D_printed / 5
```

with the energy per atom carrying the usual double-counting factor,
`E = ½ Σ_j phi2(r_j)`.

That single correction puts **four independent** results of the paper in place
at once (below). Without it the energy, the pressure, the force constants and
the elastic constants all come out 5× wrong — while the simulation still
"runs", which is to say it fails silently.

| element | D_printed (10⁻²⁹ J·m) | D_eff (eV·Å) | α (Å⁻¹) | r₀ (Å) | m |
|---|---|---|---|---|---|
| Pd | 7.66429 | 0.956870 | 2.41786 | 2.77323 | 2.5 |
| Fe | 6.19854 | 0.773764 | 2.16735 | 2.77034 | 3.5 |

## 2. Verification (all with D_eff, with nothing left free to choose)

| test | computed | paper | source |
|---|---|---|---|
| E/atom, Pd | −1.0999 eV | −1.10 (ε₀) | Table I input |
| E/atom, Fe | −0.8999 eV | −0.90 (ε₀) | Table I input |
| dE/da at a = 3.8720 Å | ~6×10⁻⁵ eV/Å | 0 (equilibrium condition) | Eq. (2) |
| C₁₁ / C₁₂, Pd | 2.62 / 1.40 | 2.61 / 1.40 | Table II |
| C₁₁ / C₁₂, Fe | 2.45 / 1.30 | 2.45 / 1.29 | Table II |
| α₁, Pd | 55192.9 | 55195.9 | Table IV |
| β₁, Pd | −243.458 | −243.570 | Table IV |
| α₁, Fe | 51103.5 | 51106.8 | Table IV |
| β₁, Fe | −217.837 | −217.950 | Table IV |

End to end through LAMMPS (`in.validate_Pd`):

```
PE/atom (a=3.8720 Å) = -1.09988 eV      pressure = -7.2 bar (~0)
a after box/relax    =  3.871995 Å      (target 3.8720 Å)
```

**The 1997 paper (Fe–28%Pd) does NOT pass the same tests.** The α ratios are
5.18/4.90 and the β ratios 6.32/6.59 — there is no clean common factor; and
the lattice sum has its minimum at 3.755 Å rather than the reported 3.750 Å.
Either the printed (D, α, r₀) or the Table IV values are corrupted (the scan
also shows 7↔1 digit confusions). Only the **1995 parameter set** is used here.

## 3. Files

| file | contents |
|---|---|
| `AkgunUgur1995.table` | `pair_style table` data: `Pd_Pd`, `Fe_Fe`, `PdFe_mean`, `Pd_Fe` |
| `in.validate_Pd` | static check: pure fcc Pd, E/atom and pressure |
| `in.alloy_mean` | average-crystal Pd–10%Fe, NVT at 300 K |
| `validate.py` | produces the verification table above |
| `gen_lammps.py` | regenerates the table and the inputs (to change range or resolution) |
| `phonon_check.py` | two-body vs two- plus three-body phonon frequencies |

Usage:

```
lmp -in in.validate_Pd
lmp -in in.alloy_mean
```

Table entries:

- `Pd_Pd`, `Fe_Fe` — the pure elements, fitted at the **alloy's** lattice
  constant.
- `PdFe_mean` — the average crystal, `0.9*phi_Pd + 0.1*phi_Fe`. Averaging the
  pair function reproduces the paper's Eq. (10) force-constant average
  **exactly**, because differentiation is linear. Mass = 101.3625 amu. This is
  the recommended route for the alloy.
- `Pd_Fe` — **NOT IN THE PAPERS.** They give no cross Pd–Fe potential (the
  average-crystal model does not need one). This entry was made up with an
  arithmetic mean of `m, α, r₀` and a geometric mean of `D`. It is there as a
  convenience for anyone wanting to build a genuine two-type random alloy; it
  has no physical basis, and should not be used in a publication.

The cutoff of 7.80 Å sits just beyond the 8th fcc shell (7.7440 Å) — the paper
carries the two-body interaction out to the 8th neighbour. The energy at the
cutoff is ~10⁻⁶ eV and no shift is applied.

## 4. The three-body term is NOT included, and why

```
phi3(r1,r2) = C*D / (2(m-1)(r1+r2)) * { beta^m exp[-m*alpha*(r1+r2)] - m*beta*exp[-alpha*(r1+r2)] }
```

It cannot be carried into LAMMPS as it stands, for three separate reasons:

1. **It has no angular dependence.** It depends only on `(r1+r2)`. Neither
   Stillinger–Weber, nor Tersoff, nor MEAM — no standard LAMMPS three-body
   style — has that functional form.

2. **The papers never use it as a force field.** They extract a single scalar
   from `phi3`: `γ = phi3''`. That γ then goes into the phenomenological
   force-constant matrix of Mishra *et al.*
   (`D^m_αα = 4γ[4−2C₂ᵢ−Cᵢ(Cⱼ+C_k)]` and so on). Those matrix elements are
   **not** the analytic dynamical matrix of `phi3`; they are an ansatz specific
   to the fcc lattice. So there is no consistent (energy, force) pair — there
   is nothing to tabulate.

3. **The triplet counting rule is lattice-specific.** "In fcc, a first
   neighbour counts as the common nearest neighbour of the second and third
   neighbours" is not a cutoff radius; it is a statement about ideal fcc
   geometry. It cannot be applied to a thermally displaced or disordered MD box
   without generalising it, and how you generalise it changes the answer.

**What is lost in practice?** `phonon_check.py` estimates the three-body
contribution from the paper's Eq. (12): a 0–10 % stiffening of the frequencies,
largest in the [111] and [110] transverse branches. So a pair-only LAMMPS gives
the **dashed** curves of Fig. 1 — not the **solid** ones the paper shows
agreeing with experiment. (The index convention of Eq. (12) is ambiguous, so
those percentages are indicative rather than exact.)

If you want it anyway, there are two routes:

- `exp[-α(r1+r2)] = exp(-αr1)·exp(-αr2)` factorises; only the `1/(r1+r2)`
  prefactor spoils the separation. Replace it with `1/(2r̄)` and the term
  becomes exactly separable, at which point the `Σ_triplets f(r1)f(r2)`
  structure is equivalent to an EAM-type embedding function with `F(ρ)=ρ²` →
  it can be built as a `hybrid/overlay` of two `eam/alloy` tables. That is an
  approximation, not the paper's result.
- Or write a `pair_style` in C++ for the true triplet form. The formula is
  simple; the hard part is deciding the counting rule of point (3).

## 5. Other limitations

- The parameters were fitted at **a single volume** (a = 3.8720 Å).
  Transferability to other volumes is untested; this is not suitable for
  pressure work.
- The input data (ε₀, k) come from **dimer** spectroscopy (Morse 1986;
  Lin/Strauss/Kant 1969). So the cohesive energy is not a solid-state value:
  −1.10 eV for Pd, against a true cohesive energy of −3.89 eV.
- A purely central pair interaction forces the Cauchy relation `C₁₂ = C₄₄`.
  The paper knowingly steps around this and uses `C₄₄ = ⅓(2C₁₁−C₁₂)`
  (Milstein–Rasky), which is **not valid** in LAMMPS. The C₄₄ you measure from
  LAMMPS will equal C₁₂ (~1.40×10¹¹ N/m²), not the paper's 1.27.
- The thermal pressure is high: ~3 GPa at fixed volume at 300 K, i.e. a thermal
  expansion coefficient about 1.5× experiment. Typical of simple pair
  potentials.
- LAMMPS warns "force inconsistent with -dE/dr" at 1 of the 15000 points of
  `PdFe_mean`. E and F come from the same analytic formula; the warning comes
  from a finite-difference check at an inflection point (LAMMPS itself notes it
  "should only be flagged at inflection points"). Harmless.
