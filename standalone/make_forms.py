"""Build the form-comparison page: a study, not an arm of the library.

The interface's five arms are all the SAME functional form - phi2 + phi3 with an
optional angular factor - differing in cutoff, taper and fitting target. The two
forms compared here replace phi3 with an embedding, so none of the producers can
evaluate them, `pair_ugur` cannot run them, and `export_potentials.py` cannot ship
them. Wiring them in as arms would put an entry in the menu that nobody can
download or reproduce, which is the "silently dead arm" failure in a new guise.

They are also not offered for use. So this is a separate page that states what
was measured and what it means, and it touches nothing the library ships.

The page has to STAND ALONE: the project notes are .docx files and are not
distributed, so pointing a reader at a section number is pointing at nothing.
Whatever a reader needs in order to judge a number - the objective, the anchor
weight and why it differs per element, how the hold-out was chosen, how each
screen is defined, how the intervals were computed - belongs on the page.

Every number below is typed in from the study's result files, and each table
says what was run to produce it. Nothing here is recomputed by the page, so a
changed result has to be carried into this file by hand.

    python make_forms.py            ->  forms.html beside this script
"""
import datetime, io, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ the forms
#  The parameters live beside this script rather than inside it: forms_params
#  .json is the same file the calculations read, so the page cannot drift from
#  what was computed.

PARAMS = json.load(open(os.path.join(HERE, "forms_params.json")))
ELEMENTS = ["W", "Mo", "Nb", "Ta", "V"]

#  which parameters each form actually has, in the order they are displayed
PARKEYS = {
    "mau": ["D", "alpha", "r0", "m", "gamma", "C", "alpha3", "rcut2", "rcut3"],
    "ug":  ["D", "alpha", "r0", "m", "gamma", "C", "alpha3", "rcut2", "rcut3",
            "lam2", "lam4"],
    "fs":  ["D", "alpha", "r0", "m", "gamma", "A", "beta", "rcut2", "rcutd"],
    "ug2": ["D", "alpha", "r0", "m", "gamma", "A", "beta", "rcut2", "rcutd",
            "t2", "t4"],
}
PARLABEL = {"D": "D", "alpha": "&alpha;", "r0": "r<sub>0</sub>", "m": "m",
            "gamma": "&gamma;", "C": "C", "alpha3": "&alpha;<sub>3</sub>",
            "rcut2": "r<sub>c2</sub>", "rcut3": "r<sub>c3</sub>",
            "rcutd": "r<sub>cd</sub>", "lam2": "&lambda;<sub>2</sub>",
            "lam4": "&lambda;<sub>4</sub>", "A": "A", "beta": "&beta;",
            "t2": "t<sub>2</sub>", "t4": "t<sub>4</sub>"}


def shape(x, m, r0, gamma, al):
    """the shape function all four forms share"""
    u = al * (r0 - x)
    return (r0 / x) ** gamma * (math.exp(m * u) - m * math.exp(u)) / (m - 1.0)


def switch(r, rc, taper=0.85):
    if r >= rc:
        return 0.0
    ron = taper * rc
    if r <= ron:
        return 1.0
    t = (r - ron) / (rc - ron)
    return 1.0 - t ** 3 * (10.0 - 15.0 * t + 6.0 * t * t)


def phi2(r, p):
    return p["D"] * shape(r, p["m"], p["r0"], p["gamma"], p["alpha"]) \
        * switch(r, p["rcut2"])


def dens(r, p):
    return math.exp(-p["beta"] * (r - p["r0"])) * switch(r, p["rcutd"])


def legendre(c, l):
    return 0.5 * (3 * c * c - 1) if l == 2 else \
        0.125 * (35 * c ** 4 - 30 * c * c + 3)


# ------------------------------------------------------------- tiny SVG plots

def svg(curves, xlim, ylim, xlabel, ylabel, w=440, h=230, zero=True):
    """polyline plot, no script and no dependency"""
    L, R, T, B = 46, 12, 12, 34
    x0, x1 = xlim
    y0, y1 = ylim

    def px(x):
        return L + (x - x0) / (x1 - x0) * (w - L - R)

    def py(y):
        return T + (1 - (y - y0) / (y1 - y0)) * (h - T - B)

    out = [f'<svg viewBox="0 0 {w} {h}" class="plot" role="img">']
    out.append(f'<rect x="{L}" y="{T}" width="{w-L-R}" height="{h-T-B}" '
               f'fill="var(--surface)" stroke="var(--line-2)"/>')
    for k in range(5):
        x = x0 + k * (x1 - x0) / 4.0
        out.append(f'<line x1="{px(x):.1f}" y1="{h-B}" x2="{px(x):.1f}" '
                   f'y2="{h-B+4}" stroke="var(--line)"/>')
        out.append(f'<text x="{px(x):.1f}" y="{h-B+15}" class="tk" '
                   f'text-anchor="middle">{x:g}</text>')
    for k in range(3):
        y = y0 + k * (y1 - y0) / 2.0
        out.append(f'<line x1="{L-4}" y1="{py(y):.1f}" x2="{L}" '
                   f'y2="{py(y):.1f}" stroke="var(--line)"/>')
        out.append(f'<text x="{L-7}" y="{py(y)+3:.1f}" class="tk" '
                   f'text-anchor="end">{y:g}</text>')
    if zero and y0 < 0 < y1:
        out.append(f'<line x1="{L}" y1="{py(0):.1f}" x2="{w-R}" y2="{py(0):.1f}" '
                   f'stroke="var(--line)" stroke-dasharray="3 3"/>')
    for label, colour, dash, pts in curves:
        d = " ".join(f"{px(x):.1f},{py(max(min(y, y1), y0)):.1f}" for x, y in pts)
        da = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<polyline points="{d}" fill="none" stroke="{colour}" '
                   f'stroke-width="1.8"{da}/>')
    out.append(f'<text x="{(L+w-R)/2:.0f}" y="{h-4}" class="ax" '
               f'text-anchor="middle">{xlabel}</text>')
    out.append(f'<text x="12" y="{(T+h-B)/2:.0f}" class="ax" '
               f'text-anchor="middle" transform="rotate(-90 12 {(T+h-B)/2:.0f})">'
               f'{ylabel}</text>')
    out.append("</svg>")
    key = " &middot; ".join(
        f'<span style="color:{c}">&#9473;&#9473;</span> {lab}'
        for lab, c, _, _ in curves)
    return "".join(out) + f'<p class="key">{key}</p>'


C_MAU, C_UG, C_FS, C_UG2 = ("var(--phi2)", "var(--phi3)", "var(--good)",
                            "var(--bad)")


def plots_html(el="W"):
    p = {f: PARAMS[f][el] for f in ("mau", "ug", "fs", "ug2")}
    #  1) the pair term
    rr = [1.8 + i * (6.6 - 1.8) / 160 for i in range(161)]
    c1 = [("MAU", C_MAU, "", [(r, phi2(r, p["mau"])) for r in rr]),
          ("UG", C_UG, "6 3", [(r, phi2(r, p["ug"])) for r in rr]),
          ("FS", C_FS, "2 2", [(r, phi2(r, p["fs"])) for r in rr]),
          ("UG-2", C_UG2, "5 2 1 2", [(r, phi2(r, p["ug2"])) for r in rr])]
    lo = min(min(y for _, y in c[3]) for c in c1)
    plot1 = svg(c1, (1.8, 6.6), (max(-1.5, lo * 1.1), 1.5),
                "r (&#8491;)", "&phi;&#8322; (eV)")
    #  2) density and embedding
    c2 = [("FS density f(r)", C_FS, "2 2", [(r, dens(r, p["fs"])) for r in rr]),
          ("UG-2 density f(r)", C_UG2, "5 2 1 2",
           [(r, dens(r, p["ug2"])) for r in rr])]
    #  the axes follow the data: the density is not bounded by one, and rho-bar
    #  runs to roughly a first shell of it.  Hard-coding 1.6 and 12 clipped both
    #  curves flat against the frame.
    fmax = max(y for _, _, _, pts in c2 for _, y in pts)
    plot2 = svg(c2, (1.8, 6.6), (0, math.ceil(fmax * 1.15)),
                "r (&#8491;)", "f (r)", zero=False)
    rnn = min(((r, phi2(r, p["mau"])) for r in rr), key=lambda t: t[1])[0]
    rho_max = math.ceil(8 * max(dens(rnn, p["fs"]), dens(rnn, p["ug2"])) / 10) * 10
    rho = [0.02 + i * (rho_max - 0.02) / 160 for i in range(161)]
    c3 = [("FS  F = &minus;A&#8730;&rho;&#772;", C_FS, "2 2",
           [(x, -p["fs"]["A"] * math.sqrt(x)) for x in rho]),
          ("UG-2", C_UG2, "5 2 1 2",
           [(x, -p["ug2"]["A"] * math.sqrt(x)) for x in rho])]
    fmin = min(y for _, _, _, pts in c3 for _, y in pts)
    plot3 = svg(c3, (0, rho_max), (math.floor(fmin), 0), "&rho;&#772;", "F (eV)")
    #  3) the angular objects
    cc = [-1 + i * 2.0 / 160 for i in range(161)]
    ugp, u2p = p["ug"], p["ug2"]
    c4 = [("UG  h(cos&theta;)", C_UG, "6 3",
           [(c, 1 + ugp["lam2"] * legendre(c, 2) + ugp["lam4"] * legendre(c, 4))
            for c in cc]),
          ("UG-2  1 + t&#8322;P&#8322; + t&#8324;P&#8324;", C_UG2, "5 2 1 2",
           [(c, 1 + u2p["t2"] * legendre(c, 2) + u2p["t4"] * legendre(c, 4))
            for c in cc])]
    ys = [y for _, _, _, pts in c4 for _, y in pts]
    plot4 = svg(c4, (-1, 1), (min(-1, min(ys) * 1.1), max(3, max(ys) * 1.1)),
                "cos &theta;", "angular weight")
    return (f'<div class="grid">'
            f'<figure>{plot1}<figcaption>The pair term of all four forms for '
            f'tungsten. They share one shape function; what differs is what sits '
            f'beside it.</figcaption></figure>'
            f'<figure>{plot2}<figcaption>The density each embedding sums over a '
            f'neighbourhood. UG-2&rsquo;s is shorter-ranged than the plain '
            f'embedding&rsquo;s.</figcaption></figure>'
            f'<figure>{plot3}<figcaption>The embedding itself, &minus;A&radic;&rho;&#772;, '
            f'over the range a first shell of eight neighbours reaches. Its '
            f'curvature is what makes the energy respond to coordination rather '
            f'than to each bond separately.</figcaption></figure>'
            f'<figure>{plot4}<figcaption>The angular objects, and they are not the '
            f'same thing: UG multiplies its three-body term by h(cos&theta;), '
            f'while UG-2&rsquo;s t&#8322; and t&#8324; weight partial densities '
            f'inside the square root.</figcaption></figure>'
            f'</div>')


def params_html():
    out = []
    for form, name in (("mau", "MAU"), ("ug", "UG"), ("fs", "FS"),
                       ("ug2", "UG-2")):
        keys = PARKEYS[form]
        head = "".join(f"<th>{PARLABEL[k]}</th>" for k in keys)
        rows = []
        for el in ELEMENTS:
            p = PARAMS[form][el]
            cells = "".join(f"<td>{p[k]:.4g}</td>" for k in keys)
            rows.append(f"<tr><td class='k'>{el}</td>{cells}</tr>")
        out.append(f"<h3>{name} &middot; {len(keys)} parameters</h3>"
                   f"<table><thead><tr><th>Element</th>{head}</tr></thead>"
                   f"<tbody>{''.join(rows)}</tbody></table>")
    return "".join(out)


# ---------------------------------------------------------------- references

REFS = [
 ("The density-functional data", [
  ("Tungsten &mdash; the W-14 set, 280 training and 280 hold-out configurations.",
   "W. J. Szlachta, A. P. Bart&oacute;k and G. Cs&aacute;nyi, "
   "<i>Accuracy and transferability of Gaussian approximation potential models "
   "for tungsten</i>, Phys. Rev. B <b>90</b>, 104108 (2014).",
   "10.1103/PhysRevB.90.104108",
   "CASTEP 6.01, PBE, ultrasoft, 600 eV, k-spacing 0.015 &Aring;<sup>&minus;1</sup>, "
   "Gaussian smearing 0.1 eV."),
  ("Molybdenum, niobium, tantalum and vanadium &mdash; 464 training "
   "configurations each, subsampled by configuration type.",
   "J. Byggm&auml;star, K. Nordlund and F. Djurabekova, "
   "<i>Gaussian approximation potentials for body-centered-cubic transition "
   "metals</i>, Phys. Rev. Materials <b>4</b>, 093802 (2020).",
   "10.1103/PhysRevMaterials.4.093802",
   "Data via ColabFit Exchange, CC BY 4.0. VASP, PBE, ENCUT 500 eV, "
   "KSPACING 0.15 &Aring;<sup>&minus;1</sup>, ISMEAR 1, SIGMA 0.1 eV. "
   "The ColabFit dataset names say PRM2019; the paper is 2020."),
  ("The alloys &mdash; every binary, ternary, quaternary and quinary "
   "configuration used here, and the independent test file.",
   "J. Byggm&auml;star, K. Nordlund and F. Djurabekova, "
   "<i>Modeling refractory high-entropy alloys with efficient machine-learned "
   "interatomic potentials: Defects and segregation</i>, "
   "Phys. Rev. B <b>104</b>, 104101 (2021); arXiv:2106.03369.",
   "10.1103/PhysRevB.104.104101",
   "Files <code>db_HEA_v2.xyz</code> and <code>testset_HEA_v2_all.xyz</code>, "
   "from the Fairdata deposit <i>Data and tabGAP v.2 potential files for "
   "Mo-Nb-Ta-V-W alloys</i> (J. Byggm&auml;star, 2022), "
   "<a href=\'https://doi.org/10.23729/1e6d0215-d26b-4f5f-8f5b-df575efa6594\'>"
   "doi:10.23729/1e6d0215-d26b-4f5f-8f5b-df575efa6594</a>."),
 ]),
 ("The functional forms", [
  ("MAU and UG &mdash; the published form and its Legendre angular factor.",
   "&#304;. Akg&uuml;n and G. U&#287;ur, Phys. Rev. B <b>51</b>, 3458 (1995); "
   "Nuovo Cimento D <b>19</b>, 779 (1997); Nuovo Cimento D <b>20</b>, 1549 (1998).",
   "10.1103/PhysRevB.51.3458",
   "The 1998 paper carries the five-parameter version with the "
   "(r<sub>0</sub>/r)<sup>&gamma;</sup> prefactor used here."),
  ("FS &mdash; the square-root embedding.",
   "M. W. Finnis and J. E. Sinclair, <i>A simple empirical N-body potential for "
   "transition metals</i>, Philos. Mag. A <b>50</b>, 45 (1984).",
   "10.1080/01418618408244210",
   "The right citation for the square root, but their own pair and density "
   "functions are cubic splines, which these are not."),
  ("FS &mdash; what the form in these files actually is.",
   "V. Rosato, M. Guillop&eacute; and B. Legrand, <i>Thermodynamical and "
   "structural properties of f.c.c. transition metals using a simple "
   "tight-binding model</i>, Philos. Mag. A <b>59</b>, 321 (1989).",
   "10.1080/01418618908205062",
   "An exponential density under a square root is the second-moment "
   "tight-binding form, not Finnis and Sinclair&rsquo;s."),
  ("FS &mdash; the same form for transition metals and their alloys.",
   "F. Cleri and V. Rosato, <i>Tight-binding potentials for transition metals "
   "and alloys</i>, Phys. Rev. B <b>48</b>, 22 (1993).",
   "10.1103/PhysRevB.48.22",
   "LAMMPS&rsquo;s <code>pair_style eam/fs</code> means something else again: "
   "there the letters stand for element-pair-specific density functions."),
  ("The embedding idea.",
   "M. S. Daw and M. I. Baskes, <i>Embedded-atom method: Derivation and "
   "application to impurities, surfaces, and other defects in metals</i>, "
   "Phys. Rev. B <b>29</b>, 6443 (1984).",
   "10.1103/PhysRevB.29.6443", ""),
  ("UG-2 &mdash; the precedent for angular partial densities.",
   "M. I. Baskes, <i>Modified embedded-atom potentials for cubic materials and "
   "impurities</i>, Phys. Rev. B <b>46</b>, 2727 (1992).",
   "10.1103/PhysRevB.46.2727",
   "UG-2&rsquo;s &rho;&#772;<sup>2</sup> = &rho;<sub>0</sub><sup>2</sup> + "
   "t<sub>2</sub>Q<sub>2</sub> + t<sub>4</sub>Q<sub>4</sub> follows the same "
   "construction, with Q<sub>l</sub> a sum of squares by the addition theorem."),
 ]),
 ("Comparison targets and method", [
  ("The ceiling of the two- and three-body class.",
   "S. Pozdnyakov, A. R. Oganov, E. Mazhnik, A. Mazitov and I. Kruglov, "
   "<i>Fast general two- and three-body interatomic potential</i>, "
   "Phys. Rev. B <b>107</b>, 125160 (2023).",
   "10.1103/PhysRevB.107.125160",
   "Free splines for &phi;<sub>2</sub> and &phi;<sub>3</sub>, fitted by linear "
   "regression at a cost independent of the parameter count."),
  ("Why a higher body order is not automatically progress.",
   "S. Chong, T. Jiang, M. Domina, F. Bigi, F. Grasselli, J. Lee and "
   "M. Ceriotti, <i>Resolving the body-order paradox of machine learning "
   "interatomic potentials</i>, J. Chem. Phys. <b>164</b>, 064121 (2026).",
   "10.1063/5.0303302",
   "On hydrogen clusters with machine-learned models, so it frames the "
   "reasoning rather than proving anything about these forms."),
  ("The density-functional elastic constants of the quaternary.",
   "X.-G. Li, C. Chen, H. Zheng, Y. Zuo and S. P. Ong, <i>Complex "
   "strengthening mechanisms in the NbMoTaW multi-principal element alloy</i>, "
   "npj Comput. Mater. <b>6</b>, 70 (2020).",
   "10.1038/s41524-020-0339-0",
   "Their Table 1 &mdash; a DFT reference rather than another potential."),
  ("The mechanical stability criteria.",
   "F. Mouhat and F.-X. Coudert, <i>Necessary and sufficient elastic stability "
   "conditions in various crystal systems</i>, Phys. Rev. B <b>90</b>, 224104 "
   "(2014).",
   "10.1103/PhysRevB.90.224104", ""),
 ]),
]


# --------------------------------------------------------------------- data
#  Each table records where it came from. The page prints that line under it.

CEILING = dict(
    caption="Hold-out force error, eV/Å, at each element's knee anchor weight. "
            "The spline column is a free two- and three-body potential in the "
            "manner of Pozdnyakov et al., between about 270 and 950 parameters "
            "against nine or eleven, fitted by linear regression on the same "
            "training set.",
    source="Free-spline fits on the same training sets at four to six grid "
           "resolutions per element. For each resolution the ridge is chosen on "
           "a validation slice of the training set, the resolution is then "
           "chosen on the same slice, and the hold-out is scored once. For the four "
           "forms, each element was fitted in up to three searches at its knee weight "
           "(cold, warm-started, and a converged 42-start run); the fit reported is the one "
           "with the lowest TRAINING cost, and the hold-out plays no part in choosing it. "
           "UG-2 lands within 0.01 eV/Å across searches (tantalum 0.15); the other forms "
           "move by up to 0.23, so differences below about 0.1 eV/Å are not resolved. Forces verified against finite "
           "differences to 10<sup>−9</sup> eV/Å &mdash; and the fitting "
           "machinery recovers a known potential to 0.0002 eV/Å on a synthetic "
           "hold-out.",
    cols=["Element", "MAU", "UG", "FS", "UG-2", "free splines", "UG-2 / splines"],
    rows=[["W (240)", "0.409", "0.396", "0.299", "0.294", "0.091", "3.2×"],
          ["Mo (380)", "0.431", "0.400", "0.308", "0.308", "0.096", "3.2×"],
          ["Nb (272)", "0.314", "0.328", "0.407", "0.236", "0.076", "3.1×"],
          ["Ta (385)", "0.397", "0.374", "0.432", "0.248", "0.082", "3.0×"],
          ["V (248)", "0.306", "0.308", "0.238", "0.248", "0.069", "3.6×"]],
    best=[None, None, None, None, None, 5, None])

ANCHOR = dict(
    caption="The same four forms refitted with the experimental anchors switched "
            "off (W_EXP = 0), same data and quota, eight cold starts. Hold-out "
            "force error in eV/Å.",
    source="The same fitting procedure with the anchor weight set to zero, eight "
           "cold starts per form, on the same data and the same subsample.",
    cols=["Fit", "MAU", "UG", "FS", "UG-2"],
    rows=[["Nb at its knee (272)", "0.314", "0.328", "0.407", "0.236"],
          ["Nb with no anchors", "0.192", "0.188", "0.207", "0.195"],
          ["Mo at its knee (380)", "0.431", "0.400", "0.308", "0.308"],
          ["Mo with no anchors", "0.291", "0.287", "0.307", "0.298"]],
    best=[None, None, None, None, None])

PRICE = dict(
    caption="Read as a price: what each form gives up in force accuracy in order "
            "to satisfy the experimental anchors.",
    source="Derived from the table above.",
    cols=["Element", "MAU", "UG", "FS", "UG-2"],
    rows=[["Nb", "+64 %", "+74 %", "+97 %", "+21 %"],
          ["Mo", "+48 %", "+40 %", "+0 %", "+3 %"]],
    best=[None, None, None, None, None])

BINARIES = dict(
    caption="Five-fold cross-validation of the ten bcc binaries: every one of a "
            "pair's 39 configurations held out exactly once. Hold-out force "
            "error in eV/Å. The last column's confidence interval excludes zero "
            "everywhere except the two marked n.s.",
    source="Five folds over each pair’s 39 configurations, two hundred fits in "
           "all, giving 390 independent hold-out configurations. Only the five "
           "cross-pair parameters are free; the single-element fits stay frozen "
           "at their knee weights.",
    cols=["Pair", "Group", "MAU", "UG", "FS", "UG-2", "UG-2 vs FS"],
    rows=[["Mo-Nb", "6-5", "0.568", "0.561", "0.713", "0.572", "−19.8 %"],
          ["Mo-Ta", "6-5", "0.549", "0.602", "0.648", "0.574", "−11.4 %"],
          ["Mo-V", "6-5", "0.595", "0.516", "0.551", "0.532", "−3.4 %"],
          ["Mo-W", "6-6", "1.127", "1.063", "0.697", "0.664", "−4.8 %"],
          ["Nb-Ta", "5-5", "0.627", "0.663", "0.522", "0.463", "−11.1 %"],
          ["Nb-V", "5-5", "0.462", "0.408", "0.526", "0.590", "+12.2 %"],
          ["Nb-W", "5-6", "1.292", "1.200", "1.265", "0.701", "−44.6 %"],
          ["Ta-V", "5-5", "0.626", "0.585", "0.451", "0.422", "−6.3 %"],
          ["Ta-W", "5-6", "1.269", "1.223", "0.710", "0.689", "−3.0 % n.s."],
          ["V-W", "5-6", "1.144", "1.037", "0.681", "0.675", "−0.9 % n.s."],
          ["pooled", "—", "0.885", "0.839", "0.710", "0.595", "−16.1 %"]],
    best=[None, None, None, None, None, None, None])

TRANSFER = dict(
    caption="Zero free parameters: the single-element fits frozen and the binary "
            "cross terms used exactly as fitted, applied to configurations with "
            "three, four and five species from the authors' own test file. "
            "Force error in eV/Å; the DFT column is the error of predicting no "
            "force at all.",
    source="240 configurations with three, four and five species, taken from the "
           "dataset authors’ own independent test file, which no fit in this "
           "work has touched. Nothing is refitted for this table.",
    cols=["Set", "n", "DFT rms", "MAU", "UG", "FS", "UG-2"],
    rows=[["crystals, 3/4/5 species", "120", "0.956", "0.655", "0.618", "0.685", "0.482"],
          ["five-species HEA", "40", "0.868", "0.686", "0.637", "0.678", "0.494"],
          ["all non-liquid", "200", "1.126", "0.695", "0.654", "0.674", "0.510"],
          ["W-Mo-Ta, liquid included", "80", "2.697", "1.224", "1.193", "1.099", "1.194"]],
    best=[None, None, None, None, None, None, None])

EXTERNAL = dict(
    caption="Equiatomic MoNbTaVW at 0 K against the table published by the group "
            "whose data this work uses. Elastic constants in GPa, mixing energy "
            "in meV per atom. Still zero free parameters.",
    source="128-atom random cells at their own relaxed volume, three chemical "
           "arrangements, positions relaxed at every strain. Reference: "
           "Byggmästar, Nordlund and Djurabekova, Phys. Rev. B 104, 104101 "
           "(2021), Table II.",
    cols=["", "a (Å)", "E<sub>mix</sub>", "B", "C<sub>11</sub>",
          "C<sub>12</sub>", "C<sub>44</sub>"],
    rows=[["reference (tabGAP)", "3.195", "−41.85", "210.4", "382.2", "124.5", "47.5"],
          ["MAU", "3.16", "+1556", "234", "321", "191", "112"],
          ["UG", "3.17", "+1380", "214", "279", "182", "99"],
          ["FS", "3.20", "+973", "180", "217", "161", "65"],
          ["UG-2", "3.14", "−160", "239", "330", "194", "116"]],
    best=[None, None, None, None, None, None, None])

STABILITY = dict(
    caption="Mechanical stability of the ten binaries with the fitted cross term, "
            "counted over pair-and-seed runs out of twenty. Bain: the cubic cell "
            "is a local minimum and the face-centred end is uphill. Born: the "
            "relaxed C11−C12, C44 and B are all positive.",
    source="54-atom random 50:50 solutions, two chemical arrangements per pair, "
           "reference volume taken from the alloy’s own energy minimum rather "
           "than an experimental lattice table.",
    cols=["Screen", "MAU", "UG", "FS", "UG-2"],
    rows=[["Bain, unfitted mixing rule", "18/20", "20/20", "8/20", "13/20"],
          ["Bain, fitted cross term", "18/20", "17/20", "14/20", "19/20"],
          ["Bain, violations deeper than 5 meV", "2", "2", "5", "0"],
          ["Born criteria, fitted cross term", "19/20", "20/20", "18/20", "20/20"]],
    best=[None, None, None, None, None])

DEFECTS = dict(
    caption="Point defects in the five-species alloy, zero free parameters, "
            "against the published values. Energies in eV.",
    source="128-atom random equiatomic cells, three arrangements, one removed "
           "site per species. Migration energy is the maximum of a constrained "
           "path — the hopping atom pinned at nine fractions — because in a "
           "random alloy the midpoint of a jump is not a saddle. A collapsed "
           "vacancy (a neighbour falls into the site) is excluded from both "
           "means. A path lying entirely below the relaxed vacancy gives a "
           "negative 'barrier', which means the vacancy was not a minimum; "
           "those have their own row and are excluded from the barrier mean. "
           "Reference: Byggmästar et al. (2021): the vacancy value is their DFT "
           "average, the migration range is their tabGAP's own mean over some "
           "1100 nudged-elastic-band barriers, from Nb (1.06) to W (1.56).",
    cols=["Quantity", "published", "MAU", "UG", "FS", "UG-2"],
    rows=[["mean vacancy formation", "3.3 (DFT)", "1.11", "1.66", "1.16 (10)", "1.97"],
          ["mean migration barrier", "1.06–1.56 (tabGAP)", "0.21", "0.13", "0.39", "0.54"],
          ["cases with no barrier at all", "—", "7/14", "9/15", "1/8", "3/15"],
          ["negative barrier: vacancy not a minimum", "—", "1/15", "0", "2/10", "0"],
          ["vacancy collapsed", "—", "0", "0", "5/15", "0"]],
    best=[None, None, None, None, None, None])

TABLES = [("The ceiling of the form class", CEILING),
          ("What the anchors cost", ANCHOR),
          ("The same numbers as a price", PRICE),
          ("Ten bcc binaries, cross-validated", BINARIES),
          ("Transfer to three, four and five species", TRANSFER),
          ("Mechanical stability", STABILITY),
          ("Against a published table", EXTERNAL),
          ("Point defects against published values", DEFECTS)]

FORMS_DEF = [
    ("MAU", "phi2 + phi3(r1 + r2)", "9",
     "the published Akgün–Uğur form."),
    ("UG", "MAU × h(cosθ) = 1 + λ₂P₂ + λ₄P₄", "11",
     "the published form with its Legendre angular factor."),
    ("FS", "phi2 + F(ρ) = −A√ρ,  ρ = Σ exp(−β(r − r₀))", "9",
     "the three-body term replaced by a second-moment embedding. Called FS in "
     "the working files for Finnis–Sinclair, but the exponential density with a "
     "square-root band term is the second-moment tight-binding form of Rosato, "
     "Guillopé and Legrand and of Cleri and Rosato; LAMMPS's pair_style eam/fs "
     "means something else again."),
    ("UG-2", "FS with ρ̄² = ρ₀² + t₂Q₂ + t₄Q₄", "11",
     "the embedding with Legendre partial densities, Q_l a sum of squares by "
     "the addition theorem."),
]


def refs_html():
    out = []
    for group, items in REFS:
        out.append(f"<h3>{group}</h3><ol class=\'refs\'>")
        for what, cite, doi, note in items:
            link = f" <a href=\'https://doi.org/{doi}\'>doi:{doi}</a>" if doi else ""
            tail = f"<br><span class=\'rn\'>{note}</span>" if note else ""
            out.append(f"<li><b>{what}</b><br>{cite}{link}{tail}</li>")
        out.append("</ol>")
    return "".join(out)


def table_html(t):
    head = "".join(f"<th>{c}</th>" for c in t["cols"])
    body = []
    for r in t["rows"]:
        cells = []
        for i, v in enumerate(r):
            cls = ' class="k"' if i == 0 else ""
            cells.append(f"<td{cls}>{v}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table>'
            f'<p class="cap">{t["caption"]}</p>'
            f'<p class="src">{t["source"]}</p>')


CSS = """
:root{
  --paper:#E9EDF2; --surface:#FBFCFD; --sunk:#E1E7EE;
  --ink:#131A24; --ink-2:#43546A; --ink-3:#6E8098;
  --line:#C6D0DB; --line-2:#DCE3EB;
  --phi2:#2E5C8A; --phi3:#B4622F;
  --good:#2F7250; --mid:#8A6A18; --bad:#9E3B2E;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --paper:#0F151C; --surface:#172029; --sunk:#111921;
  --ink:#E4EAF0; --ink-2:#A8B7C7; --ink-3:#77899C;
  --line:#2C3A49; --line-2:#222E3A;
  --phi2:#79ADDD; --phi3:#E29A63;
  --good:#5FB98A; --mid:#C9A44A; --bad:#D9705F;
}}
:root[data-theme="dark"]{
  --paper:#0F151C; --surface:#172029; --sunk:#111921;
  --ink:#E4EAF0; --ink-2:#A8B7C7; --ink-3:#77899C;
  --line:#2C3A49; --line-2:#222E3A;
  --phi2:#79ADDD; --phi3:#E29A63;
  --good:#5FB98A; --mid:#C9A44A; --bad:#D9705F;
}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",sans-serif;
  line-height:1.55;font-size:15px}
.wrap{max-width:1000px;margin:0 auto;padding:32px 16px 80px}
header{border-bottom:1px solid var(--line);padding-bottom:18px;margin-bottom:26px}
h1{font-size:1.5rem;margin:0 0 6px;letter-spacing:-.01em}
h2{font-size:1.15rem;margin:38px 0 10px;color:var(--ink);
   border-left:3px solid var(--phi3);padding-left:10px}
h3{font-size:1rem;margin:26px 0 8px;color:var(--ink-2)}
p{margin:10px 0;color:var(--ink-2)}
.lede{color:var(--ink);font-size:1.02rem}
.warn{background:var(--sunk);border:1px solid var(--line);
  border-left:3px solid var(--bad);padding:12px 14px;border-radius:5px;margin:18px 0}
.warn b{color:var(--ink)}
table{border-collapse:collapse;width:100%;margin:16px 0 6px;
  background:var(--surface);font-variant-numeric:tabular-nums;font-size:.92rem}
th,td{border:1px solid var(--line-2);padding:6px 9px;text-align:right}
th{background:var(--sunk);color:var(--ink);font-weight:600;text-align:right}
th:first-child,td.k{text-align:left;color:var(--ink)}
tbody tr:last-child td{font-weight:600}
.cap{font-size:.86rem;color:var(--ink-2);margin:6px 0 2px}
.src{font-size:.78rem;color:var(--ink-3);margin:0 0 4px}
dl{margin:14px 0}
dt{font-weight:600;color:var(--ink);margin-top:12px}
dd{margin:2px 0 0 0;color:var(--ink-2);font-size:.94rem}
code{background:var(--sunk);padding:1px 5px;border-radius:3px;font-size:.88em}
footer{margin-top:50px;padding-top:16px;border-top:1px solid var(--line);
  font-size:.82rem;color:var(--ink-3)}
a{color:var(--phi2)}
.eq{background:var(--sunk);border:1px solid var(--line-2);border-radius:5px;
  padding:9px 12px;margin:10px 0;color:var(--ink);font-size:.98rem;
  overflow-x:auto}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));
  gap:18px;margin:16px 0}
figure{margin:0}
figcaption{font-size:.84rem;color:var(--ink-2);margin-top:2px}
svg.plot{width:100%;height:auto;display:block}
svg .tk{font-size:9px;fill:var(--ink-3)}
svg .ax{font-size:10px;fill:var(--ink-2)}
.key{font-size:.8rem;color:var(--ink-2);margin:2px 0 0}

ol.refs{margin:10px 0 18px;padding-left:22px}
ol.refs li{margin:12px 0;color:var(--ink-2);font-size:.92rem}
ol.refs li b{color:var(--ink);font-weight:600}
.rn{color:var(--ink-3);font-size:.86rem}

.toggle{position:fixed;top:12px;right:14px;background:var(--surface);
  border:1px solid var(--line);color:var(--ink-2);border-radius:6px;
  padding:5px 10px;font-size:.8rem;cursor:pointer}
@media (max-width:640px){.wrap{padding:22px 16px 60px}table{font-size:.82rem}
  th,td{padding:5px 6px}}
"""

JS = """
(function(){
  var b=document.querySelector('.toggle');
  function cur(){return document.documentElement.getAttribute('data-theme')||'';}
  b.addEventListener('click',function(){
    var n = cur()==='dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme',n);
    try{localStorage.setItem('forms-theme',n);}catch(e){}
  });
  try{var s=localStorage.getItem('forms-theme');
      if(s){document.documentElement.setAttribute('data-theme',s);}}catch(e){}
})();
"""


def build():
    defs = "".join(
        f"<dt>{n} &middot; {p} parameters</dt>"
        f"<dd><code>{f}</code><br>{d}</dd>"
        for n, f, p, d in FORMS_DEF)
    sections = "".join(f"<h2>{title}</h2>{table_html(t)}" for title, t in TABLES)
    today = datetime.date.today().isoformat()
    refs = refs_html()
    plots = plots_html()
    params = params_html()
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Four Forms Compared</title>
<style>{CSS}</style></head>
<body><button class="toggle">theme</button><div class="wrap">
<header>
  <h1>Four forms compared</h1>
  <p class="lede">What happens to a &phi;<sub>2</sub>&nbsp;+&nbsp;&phi;<sub>3</sub>
  potential for bcc transition metals when the three-body term is replaced by a
  second-moment embedding, and what happens when Legendre partial densities are
  added to it. Five elements, ten binaries, and a five-species alloy.</p>
</header>

<div class="warn">
  <b>These are not potentials offered for use.</b> They are a controlled
  comparison of functional forms at matched parameter counts: the same data, the
  same objective, the same hold-out. They cover five of the library's
  thirty-eight elements, they are not implemented in <code>pair_ugur</code>, and
  no parameter file is shipped for them. Nothing on this page changes the
  library. The surface energies are still wrong by tens of per cent, the facet
  ordering is still wrong for tungsten and molybdenum, and the mixing energies
  below show what they cannot do at all.
</div>

<h2>The question, and the forms</h2>

<p>The published form for these metals is a pair term plus a three-body term.
Two of its answers are wrong in opposite directions: <b>the vacancy comes out too
expensive and the migration barrier far too cheap</b>, and both are properties of
what the energy does when an atom&rsquo;s <i>coordination</i> changes rather than
when a single bond stretches. That is the job the three-body term is doing. So the
experiment is to take it out, put a second-moment embedding in its place, change
nothing else, and see which half of the form was buying what.</p>

<p>All four forms are built on <b>one shape function</b>, so the pair term is
common to all of them and the comparison is not confounded by it:</p>

<p class="eq">f(x) = (r<sub>0</sub>/x)<sup>&gamma;</sup> &middot;
[ e<sup>m&alpha;(r<sub>0</sub>&minus;x)</sup> &minus;
m&thinsp;e<sup>&alpha;(r<sub>0</sub>&minus;x)</sup> ] / (m &minus; 1)</p>

<p>Every term is cut off by the same quintic switch, which is smooth to second
order at both ends and starts at 85&nbsp;per&nbsp;cent of the cutoff:</p>

<p class="eq">S(r) = 1 &minus; t<sup>3</sup>(10 &minus; 15t + 6t<sup>2</sup>),
&nbsp; t = (r &minus; 0.85&thinsp;r<sub>c</sub>) /
(0.15&thinsp;r<sub>c</sub>), clipped to [0,&nbsp;1]</p>

<h3>MAU &mdash; the published form, 9 parameters</h3>
<p class="eq">E = &frac12;&sum;<sub>ij</sub> D f(r<sub>ij</sub>) S(r<sub>ij</sub>)
&nbsp;+&nbsp; &sum;<sub>i</sub> &sum;<sub>j&lt;k</sub>
C&thinsp;D f(r<sub>ij</sub>+r<sub>ik</sub>; &alpha;<sub>3</sub>)
S(r<sub>ij</sub>) S(r<sub>ik</sub>)</p>
<p>The three-body term takes the <i>sum of the two legs</i> as its argument, and
uses the same shape function with its own range &alpha;<sub>3</sub>. Free: D,
&alpha;, r<sub>0</sub>, m, &gamma;, C, &alpha;<sub>3</sub> and the two cutoffs.</p>

<h3>UG &mdash; the published angular factor, 11 parameters</h3>
<p class="eq">&hellip; &times; h(cos&theta;<sub>jik</sub>) = 1 +
&lambda;<sub>2</sub>P<sub>2</sub>(cos&theta;) +
&lambda;<sub>4</sub>P<sub>4</sub>(cos&theta;)</p>
<p>The same three-body term multiplied by a Legendre factor in the angle at the
central atom. Adds &lambda;<sub>2</sub> and &lambda;<sub>4</sub>.</p>

<h3>FS &mdash; the second-moment embedding, 9 parameters</h3>
<p class="eq">E = &frac12;&sum;<sub>ij</sub> D f(r<sub>ij</sub>) S(r<sub>ij</sub>)
&nbsp;+&nbsp; &sum;<sub>i</sub> F(&rho;&#772;<sub>i</sub>), &nbsp;
F(&rho;&#772;) = &minus;A&radic;&rho;&#772;</p>
<p class="eq">&rho;&#772;<sub>i</sub> = &rho;<sub>0,i</sub> =
&sum;<sub>j</sub> f<sub>ij</sub>, &nbsp;
f<sub>ij</sub> = e<sup>&minus;&beta;(r<sub>ij</sub>&minus;r<sub>0</sub>)</sup>
S(r<sub>ij</sub>; r<sub>cd</sub>)</p>
<p>The three-body term is gone. In its place each atom sums a density over its
neighbours and the energy is the square root of it &mdash; the second moment of a
tight-binding band. The pair term is unchanged, and the parameter count is the
same nine: A and &beta; and a density cutoff replace C, &alpha;<sub>3</sub> and
r<sub>c3</sub>.</p>

<h3>UG-2 &mdash; the embedding with Legendre partial densities, 11 parameters</h3>
<p class="eq">&rho;&#772;<sub>i</sub><sup>2</sup> =
&rho;<sub>0,i</sub><sup>2</sup> + t<sub>2</sub>Q<sub>2,i</sub> +
t<sub>4</sub>Q<sub>4,i</sub>, &nbsp;
Q<sub>l,i</sub> = &sum;<sub>jk</sub> f<sub>ij</sub> f<sub>ik</sub>
P<sub>l</sub>(cos&theta;<sub>jik</sub>)</p>
<p>The density gains angular content. <b>Q<sub>l</sub> is never negative</b>: by
the spherical-harmonic addition theorem the double sum is a sum of squares, which
is what keeps the square root real without a floor being imposed by hand. Setting
t<sub>2</sub> = t<sub>4</sub> = 0 gives back FS exactly, so the two differ by two
parameters and nothing else &mdash; the same relation UG has to MAU.</p>

<p>Note what is <i>not</i> the same between the two angular constructions. UG
multiplies an energy by an angular factor, so a single triplet can lower the
energy or raise it. UG-2 puts the angle inside the density, under a square root
that is already concave in coordination. That is why the two behave differently
even though both are &ldquo;Legendre corrections&rdquo;.</p>

<h2>What the forms look like</h2>
{plots}

<h2>The fitted parameters</h2>
<p>These are the single-element parameter sets the alloy, transfer and
point-defect tables were built on, as read from <code>forms_params.json</code>
beside the page&rsquo;s own build script. The ceiling table reports, for each
element and form, the search with the lowest training cost. In thirteen of
the twenty cases that is the same fit; in seven (tungsten's UG-2, molybdenum's
MAU, UG and UG-2, vanadium's MAU, FS and UG-2) it is a different search of the
same objective, and every set behind that table is in
<code>forms_params_ceiling.json</code> beside it, with the search it came
from. Lengths are in &aring;ngstr&ouml;m and energies in
electronvolts; the switch turns on at 0.85 of each cutoff. <b>They are published so the numbers can
be checked, not because they are recommended:</b> see the warning at the top of
this page.</p>
{params}


<h2>What the comparison established</h2>
<p><b>The limit was the parameter count, not the form.</b> A free spline
&phi;<sub>2</sub>&nbsp;+&nbsp;&phi;<sub>3</sub> potential fitted by linear
regression sits three to three and a half times below the best of these forms on
held-out forces, and further below the rest, consistently across five elements and two groups of the periodic
table. So no claim of the kind &ldquo;a two- and three-body form cannot do
this&rdquo; is supported. What survives is a statement about compactness.</p>

<p><b>The advantage is about the anchors, not about forces.</b> Refit the same
four forms with the experimental anchors switched off and they come within ten
per cent of one another, with UG-2 third rather than first. The embedding with
Legendre densities satisfies the anchors nearly for free on both elements
tested (a fifth on niobium, 3 % on molybdenum); the angle-free and angular
three-body forms pay two fifths to three quarters, and the plain embedding
pays all of it on niobium and almost nothing on molybdenum. Carrying the anchors is
what an empirical potential is for, so this is the narrower and more useful
claim.</p>

<p><b>In alloys the added densities earn their place against the plain
embedding.</b> With the cross term fitted to binary forces only, UG-2 beats it in
seven of ten binaries, ties in two and loses in one. Against all four forms it is
the best in six; in the other four (Mo-Nb, Mo-Ta, Mo-V, Nb-V) one of the two
three-body forms is ahead of it. The cross terms then
transfer, unchanged and with no free parameter, to three, four and five species.
Mechanical stability agrees: with a fitted cross term UG-2 has no violation
deeper than relaxation noise anywhere.</p>

<p><b>Chemical energetics do not transfer.</b> Where the density-functional
energy differences between arrangements of the same composition are small, every
form errs by four to six times the signal, and the two embedding forms are
anti-correlated with the reference. The mixing energy of the five-species alloy
shows the same thing plainly: three of the four forms predict that the alloy
cannot form at all. No ordering, segregation or phase-stability claim follows
from any of this work.</p>

<h2>How these numbers were produced</h2>

<h3>The fits</h3>
<p>All four forms were fitted to the same data with the same objective, the same
hold-out and <b>matched parameter counts</b> &mdash; nine against nine, eleven
against eleven &mdash; so that the only thing that changes between columns is the
form. The objective is a least-squares sum of three blocks: <b>forces</b> at
weight 1, <b>energies per atom</b> at weight 20 with one free constant per atom
solved exactly rather than fitted, and the <b>experimental anchors</b> &mdash;
cohesive energy, lattice constant and the elastic constants &mdash; at a weight
written W<sub>exp</sub>.</p>

<p><b>The number in parentheses beside each element is its anchor weight.</b> It
is not a free choice. The fitting script&rsquo;s default is 40, at which every
form sacrifices the anchors and misses C<sub>11</sub> and C<sub>44</sub> by tens
of GPa &mdash; and at which the ranking of the forms is not the same as at the
knee. The knee measured on the tungsten set is 240, and carrying it to each
training set by the block-scale rule
W<sub>exp</sub> = 240&nbsp;&radic;(n<sub>F</sub>&sigma;<sub>F</sub><sup>2</sup> /
2.245&times;10<sup>4</sup>) gives W 240, Nb 272, Mo 380, Ta 385, V 248. Every
single-element number on this page is at its own knee, except where the table
says the anchors were switched off.</p>

<h3>The hold-out</h3>
<p>Tungsten was fitted to 280 configurations of the W-14 set and scored on 280
more. Molybdenum, niobium, tantalum and vanadium were fitted to 464
configurations each, subsampled by configuration type with a fixed quota and a
fixed seed; for those four <b>the hold-out is not a random split</b> but three
categories the fit never sees &mdash; the (112) surface, the C15 structure and
the di-vacancy. A random split would measure interpolation; withholding whole
categories measures what the fit does with a structure it has not met.</p>

<h3>The alloys</h3>
<p>For every binary the five single-element fits stay <b>frozen</b> at their knee
weights and only the five cross-pair parameters are free &mdash; the pair
depth, its range, its minimum, and the two exponents. The binaries use all 39
configurations of each pair that are neither liquid nor dimer, five-fold
cross-validated so that every configuration is held out exactly once. The
ternary fit frees fifteen to seventeen parameters across three unlike pairs and
is scored on the dataset authors&rsquo; own test file. <b>The transfer table
fits nothing at all</b>: the cross terms are used exactly as the binaries left
them, on configurations with three, four and five species from a file no fit
here has touched.</p>

<h3>The screens</h3>
<p><b>Bain</b> is the volume-conserving tetragonal path through a 54-atom random
50:50 solution with the positions relaxed at every c/a and two chemical
arrangements per pair; a pair counts as cubic when c/a = 1 is a local minimum
<i>and</i> the face-centred end of the path is uphill, and the reference volume
is the alloy&rsquo;s own energy minimum rather than an experimental lattice
table. That path probes only the tetragonal shear, so a form can pass it with a
negative C<sub>44</sub>; the <b>Born</b> row therefore computes relaxed
C<sub>11</sub>&minus;C<sub>12</sub>, C<sub>44</sub> and B from three strain modes
and requires all three to be positive.</p>

<p><b>Point defects</b> are computed in a 128-atom random equiatomic cell at its
own relaxed volume, one removed site per species, three arrangements. The
migration energy is the <b>maximum of a constrained path</b> &mdash; the hopping
atom pinned at nine fractions of the way to the vacancy while everything else
relaxes &mdash; because in a random alloy the midpoint of a jump is not a saddle:
there is no symmetry to make it one, and pinning an atom there lets the rest
rearrange into a state degenerate with the vacancy. The number is therefore a
lower bound on the true barrier.</p>

<h3>Intervals, and what &ldquo;n.s.&rdquo; means</h3>
<p>Every percentage carries a 95&nbsp;per&nbsp;cent interval from a paired
bootstrap over hold-out <i>configurations</i> &mdash; 4000 resamples, the two
forms scored on the same resampled set each time, so the comparison is paired
rather than two independent error bars. <b>n.s.</b> marks an interval that
contains zero. Pooled root-mean-square weights the hardest configurations most,
which is why the per-pair rows matter: a group average can be carried almost
entirely by one difficult pair.</p>

<h3>Verification</h3>
<p>Analytic forces were checked against finite differences for every form and
agree to better than 10<sup>&minus;8</sup> eV/&Aring; throughout, and to
10<sup>&minus;9</sup> for the free-spline one; the multi-species evaluator
reproduces the single-species one to 10<sup>&minus;14</sup>; and the spline
fitting machinery recovers a known potential to 0.0002 eV/&Aring; on a synthetic
hold-out. Each table below states what was computed for it and on how much
data.</p>

{sections}

<h2>What is not claimed</h2>
<p>None of these fits is a usable potential. The alloy results speak about forces
and mechanical stability and not about formation energies or phase diagrams. The
liquid results show the gain disappearing, so nothing here carries into melting
or transport. And no claim of the form &ldquo;alloys can be predicted from the
pure elements&rdquo; is supported: the cross term was <i>fitted</i>, to alloy
data, one binary at a time. What the transfer test shows is narrower and still
worth having &mdash; that a cross term fitted on binaries carries to three, four
and five species without refitting.</p>

<h2>References</h2>
<p>Every dataset, every functional form and every number compared against, with a
resolvable identifier. Compiled against <code>refit/DATA_REFERENCES.md</code>,
whose own header records that two citations were wrong when checked through
Crossref &mdash; a dataset&rsquo;s <i>name</i> is not its citation.</p>
{refs}

<div class="warn">
  <b>These datasets are not one ruler.</b> The tungsten set is CASTEP with
  ultrasoft pseudopotentials at 600 eV; the other four elements and all the
  alloys are VASP at 500 eV with different smearing. Numbers from different
  sources are not directly comparable and the source has to travel with the
  number. The library&rsquo;s own fitting targets are experimental &mdash;
  cohesive energies, lattice constants and elastic constants &mdash; so a fit
  that also carries anchors mixes an experimental frame with a
  density-functional one by construction.
</div>

<h2>Acknowledgement</h2>
<p>The numerical calculations reported in this work were partially performed at
TÜBİTAK ULAKBİM, High Performance and Grid Computing Center
(TRUBA resources).</p>

<footer>
  Generated {today} by <code>make_forms.py</code>, which carries the numbers and
  their provenance as data, and reads the parameters from
  <code>forms_params.json</code> beside it.
</footer>
</div><script>{JS}</script></body></html>
"""
    p = os.path.join(HERE, "forms.html")
    io.open(p, "w", encoding="utf-8").write(html)
    print(f"wrote {p}  ({len(html)/1024:.1f} kB)")
    print(f"tables: {len(TABLES)}  forms: {len(FORMS_DEF)}")


if __name__ == "__main__":
    build()
