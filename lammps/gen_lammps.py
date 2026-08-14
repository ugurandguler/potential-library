#!/usr/bin/env python3
"""Generate LAMMPS pair_style table files + input scripts for the
Akgun & Ugur (PRB 51, 3458 (1995)) two-body model potential."""
import math, os

J_PER_EV = 1.602176634e-19
OUT = r"C:\Users\Admin\Desktop\ugurpotential\lammps"
os.makedirs(OUT, exist_ok=True)

def jm(D29):
    """printed D [1e-29 J m] -> effective D [eV A], including the /5 correction"""
    return D29 * 1e-29 / J_PER_EV * 1e10 / 5.0

POT = {
    "Pd": dict(m=2.5, D=jm(7.66429), al=2.41786, r0=2.77323, mass=106.42),
    "Fe": dict(m=3.5, D=jm(6.19854), al=2.16735, r0=2.77034, mass=55.845),
}
A_ALLOY, X_FE = 3.8720, 0.10
RLO, RHI, NPT = 1.00, 7.80, 15000       # r_8 = 7.7440 A

def pair(r, m, D, al, r0):
    A = D / (2.0 * (m - 1.0)); u = al * (r0 - r)
    e1, em = math.exp(u), math.exp(m * u)
    g  = em - m * e1
    gp = m * al * (e1 - em)
    return A * g / r, A * (gp / r - g / r**2)

def phi_of(kind):
    """returns callable r -> (E, dE/dr)"""
    if kind in POT:
        p = POT[kind]
        return lambda r: pair(r, p['m'], p['D'], p['al'], p['r0'])
    if kind == "mean":                   # mean-crystal (VCA): linear mix of the pair FUNCTIONS
        def f(r):                        # -> reproduces eq.(10) force-constant averaging exactly
            e = d = 0.0
            for nm, w in (("Pd", 1-X_FE), ("Fe", X_FE)):
                p = POT[nm]; ee, dd = pair(r, p['m'], p['D'], p['al'], p['r0'])
                e += w*ee; d += w*dd
            return e, d
        return f
    if kind == "cross":                  # NOT from the papers - Lorentz-Berthelot style mixing
        a, b = POT["Pd"], POT["Fe"]
        m  = 0.5*(a['m'] + b['m']); al = 0.5*(a['al'] + b['al'])
        r0 = 0.5*(a['r0'] + b['r0']); D = math.sqrt(a['D'] * b['D'])
        return lambda r: pair(r, m, D, al, r0)
    raise KeyError(kind)

ENTRIES = [("Pd_Pd", "Pd"), ("Fe_Fe", "Fe"), ("PdFe_mean", "mean"), ("Pd_Fe", "cross")]

tbl = os.path.join(OUT, "AkgunUgur1995.table")
with open(tbl, "w") as fh:
    fh.write("# Akgun & Ugur two-body model potential -- Phys. Rev. B 51, 3458 (1995)\n"
             "#   phi(r) = D/(2(m-1)r) * [ beta^m exp(-m alpha r) - m beta exp(-alpha r) ],"
             " beta = exp(alpha*r0)\n"
             "# Units: metal (eV, Angstrom).  Truncate after the 8th fcc shell: r_8 = 7.7440 A.\n"
             "#\n"
             "# IMPORTANT: the effective D is the PRINTED table-III value divided by 5.\n"
             "#   With D_eff (and E/atom = 1/2 sum_j phi) the paper's own numbers are reproduced:\n"
             "#     E/atom  = eps0 (-1.10 eV Pd / -0.90 eV Fe),  dE/da = 0 at a = 3.8720 A,\n"
             "#     alpha_1, beta_1 = table IV to 5 digits,  C11/C12 = table II.\n"
             "#   Using D exactly as printed makes every one of those wrong by 5x.\n"
             "#\n"
             "# Parameters were fitted AT a = 3.8720 A (the Pd-10%Fe lattice constant).\n"
             "# They are not transferable to very different volumes.\n#\n")
    for key, kind in ENTRIES:
        if kind in POT:
            p = POT[kind]
            fh.write(f"# {key:10s} m={p['m']:<4} D_eff={p['D']:.6f} eV*A  "
                     f"alpha={p['al']} 1/A  r0={p['r0']} A\n")
        elif kind == "mean":
            fh.write(f"# {key:10s} mean-crystal  {1-X_FE:.2f}*phi_Pd + {X_FE:.2f}*phi_Fe "
                     f"(paper's eq.10 for the 2-body part)\n")
        else:
            fh.write(f"# {key:10s} *** NOT IN THE PAPERS *** ad-hoc mixing, see README\n")
    for key, kind in ENTRIES:
        f = phi_of(kind)
        fh.write(f"\n{key}\nN {NPT} R {RLO} {RHI}\n\n")
        for n in range(NPT):
            r = RLO + (RHI - RLO) * n / (NPT - 1)
            e, de = f(r)
            fh.write(f"{n+1} {r:.10f} {e:.12e} {-de:.12e}\n")

m_mean = (1-X_FE)*POT['Pd']['mass'] + X_FE*POT['Fe']['mass']
print(f"wrote {tbl}  ({os.path.getsize(tbl)/1e6:.1f} MB)")
print(f"mean-crystal mass = {m_mean:.4f} amu")
for key, kind in ENTRIES:
    e, _ = phi_of(kind)(RHI)
    print(f"  {key:10s} phi(7.80 A) = {e:+.3e} eV")

# ---------------- input script: static validation of pure fcc Pd ----------------
val = """# Validate the Akgun-Ugur (PRB 51, 3458) pair potential: pure fcc Pd at a = 3.8720 A
# Expected:  PE/atom = -1.0999 eV   and   Pxx=Pyy=Pzz ~ 0  (the fit imposes both)
units           metal
atom_style      atomic
boundary        p p p

lattice         fcc 3.8720
region          box block 0 5 0 5 0 5
create_box      1 box
create_atoms    1 box
mass            1 106.42

pair_style      table linear 15000
pair_coeff      1 1 AkgunUgur1995.table Pd_Pd 7.80

thermo_style    custom step pe etotal press pxx pyy pzz vol
run             0

variable        peat equal pe/atoms
print           "PE per atom  = ${peat} eV      (target -1.0999)"
print           "Pressure     = $(press) bar    (target ~0)"

# equilibrium lattice constant from a quick volume relaxation
fix             1 all box/relax iso 0.0 vmax 0.001
min_style       cg
minimize        1e-12 1e-14 2000 20000
variable        arel equal (vol/count(all)*4)^(1/3)
print           "relaxed a    = ${arel} Ang      (target 3.8720)"
print           "relaxed PE/atom = $(pe/atoms) eV"
"""
p = os.path.join(OUT, "in.validate_Pd"); open(p, "w").write(val); print("wrote", p)

# ---------------- input script: mean-crystal alloy MD ----------------
alloy = f"""# Mean-crystal (VCA) Pd-10%Fe with the Akgun-Ugur two-body potential.
# This is the 2-body half of the papers' model: the pair function itself is
# concentration-averaged, which reproduces their eq.(10) force-constant average exactly.
# The 3-body term of the papers is NOT included -- see README.
units           metal
atom_style      atomic
boundary        p p p

lattice         fcc {A_ALLOY}
region          box block 0 6 0 6 0 6
create_box      1 box
create_atoms    1 box
mass            1 {m_mean:.4f}                 # (1-x)M_Pd + x M_Fe, x = {X_FE}

pair_style      table linear {NPT}
pair_coeff      1 1 AkgunUgur1995.table PdFe_mean 7.80

velocity        all create 300.0 12345 mom yes rot yes
fix             1 all nvt temp 300.0 300.0 0.1
timestep        0.002
thermo          200
thermo_style    custom step temp pe etotal press
run             5000
"""
p = os.path.join(OUT, "in.alloy_mean"); open(p, "w").write(alloy); print("wrote", p)
