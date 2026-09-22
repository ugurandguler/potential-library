#!/usr/bin/env python3
"""
Build the standalone viewer.  Everything shown is computed by latdyn.py from the
potential's own derivatives - no external code is involved anywhere in the chain.

    python fit.py && python build_library.py && python make_gui.py

Source is deliberately pure ASCII; Greek letters go in as HTML entities.
"""
import datetime as _dt
import json, math, os

import refdata
import refdata_electronic

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "library.json")))
#  The finite-temperature table, built by make_finiteT.py from the
#  cluster results.  Kept in its own file so a new run folds in by
#  re-running that and rebuilding, without touching the renderer.
#  Absent on a tree that has not run the study - the panel then does
#  not appear, rather than appearing empty.
try:
    FT = json.load(open(os.path.join(HERE, "finiteT.json")))
except OSError:
    FT = {"rows": {}, "out": {}}

POS = {
    "Li": (2, 1), "Be": (2, 2),
    "Na": (3, 1), "Mg": (3, 2), "Al": (3, 13),
    "K": (4, 1), "Ca": (4, 2), "Ti": (4, 4), "V": (4, 5), "Cr": (4, 6),
    "Fe": (4, 8), "Co": (4, 9), "Ni": (4, 10), "Cu": (4, 11), "Zn": (4, 12),
    "Sr": (5, 2), "Zr": (5, 4), "Nb": (5, 5), "Mo": (5, 6), "Rh": (5, 9),
    "Pd": (5, 10), "Ag": (5, 11), "Cd": (5, 12),
    "Ba": (6, 2), "Ta": (6, 5), "W": (6, 6), "Ir": (6, 9), "Pt": (6, 10),
    "Au": (6, 11), "Pb": (6, 14),
    #  added 2026-08-03.  Sc and Y take group 3 of their own periods; Hf follows
    #  the lanthanides in period 6, as in the textbook layout.
    "Rb": (5, 1), "Cs": (6, 1), "Sc": (4, 3), "Y": (5, 3),
    "Hf": (6, 4), "Re": (6, 7), "Ru": (5, 8), "Tl": (6, 13),
    #  Yb and Lu are the last two lanthanides and have no slot in the main body,
    #  so they sit on a row of their own with row 7 left empty between - which
    #  is how a printed table separates the f block, and it keeps the grid at
    #  fourteen columns instead of widening the whole page for two elements.
    "Yb": (8, 13), "Lu": (8, 14),
}
NAMES = {
    "Li": "Lithium", "Be": "Beryllium", "Na": "Sodium", "Mg": "Magnesium",
    "Al": "Aluminium", "K": "Potassium", "Ca": "Calcium", "Ti": "Titanium",
    "V": "Vanadium", "Cr": "Chromium", "Fe": "Iron", "Co": "Cobalt",
    "Ni": "Nickel", "Cu": "Copper", "Zn": "Zinc", "Sr": "Strontium",
    "Zr": "Zirconium", "Nb": "Niobium", "Mo": "Molybdenum", "Rh": "Rhodium",
    "Pd": "Palladium", "Ag": "Silver", "Cd": "Cadmium", "Ba": "Barium",
    "Ta": "Tantalum", "W": "Tungsten", "Ir": "Iridium", "Pt": "Platinum",
    "Au": "Gold", "Pb": "Lead",
    "Rb": "Rubidium", "Cs": "Caesium", "Sc": "Scandium", "Y": "Yttrium",
    "Hf": "Hafnium", "Re": "Rhenium", "Ru": "Ruthenium", "Tl": "Thallium",
    "Yb": "Ytterbium", "Lu": "Lutetium",
}
for el, v in DATA.items():
    v["pos"] = POS.get(el)
    v["name"] = NAMES.get(el, el)
    #  The SYMBOL, which is otherwise only the dictionary key and so is not
    #  reachable from a record once one has been handed to a function.  Three
    #  places had reached for d.name instead and got "Copper": the download
    #  table named files that do not exist, and two per-element lookups
    #  silently missed every time.
    v["sym"] = el

REQUIRED = ("m", "gamma", "D", "alpha", "r0", "alpha3", "C", "dnn",
            "rcut2", "rcut3", "a0", "rms", "struct", "pos", "ld")


def validate(data):
    bad = []
    for el, v in sorted(data.items()):
        for k in REQUIRED:
            if k not in v or v[k] is None:
                bad.append(f"{el}: missing {k}")
            elif isinstance(v[k], float) and not math.isfinite(v[k]):
                bad.append(f"{el}: {k} is {v[k]}")
        if v.get("ld") and not v["ld"].get("std"):
            bad.append(f"{el}: no dispersion data")
        try:
            r, m, al, r0, g = v["dnn"], v["m"], v["alpha"], v["r0"], v["gamma"]
            u = al*(r0 - r)
            val = v["D"]*(r0/r)**g*(math.exp(m*u) - m*math.exp(u))/(m - 1)
            if not math.isfinite(val):
                bad.append(f"{el}: phi2(d_nn) is not finite")
        except (KeyError, TypeError, OverflowError, ZeroDivisionError) as exc:
            bad.append(f"{el}: phi2(d_nn) failed - {exc}")
    return bad


_bad = validate(DATA)
if _bad:
    raise SystemExit("library.json is incomplete, refusing to build a page "
                     "with blank plots:\n  " + "\n  ".join(_bad[:20]))

HTML = r"""<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Ugur and Guler - Interatomic Potential Library</title>
<style>
:root{
  --paper:#E9EDF2; --surface:#FBFCFD; --sunk:#E1E7EE;
  --ink:#131A24; --ink-2:#43546A; --ink-3:#6E8098;
  --line:#C6D0DB; --line-2:#DCE3EB;
--ref:#55606D;
  --phi2:#2E5C8A; --phi3:#B4622F;
  --good:#2F7250; --mid:#8A6A18; --bad:#9E3B2E; --focus:#B4622F;
  --cell:#C2726C;
}
@media (prefers-color-scheme:dark){:root{
  --paper:#0F151C; --surface:#172029; --sunk:#111921;
  --ink:#E4EAF0; --ink-2:#A8B7C7; --ink-3:#77899C;
  --line:#2C3A49; --line-2:#222E3A;
--ref:#9AA7B4;
  --phi2:#79ADDD; --phi3:#E29A63;
  --good:#5FB98A; --mid:#C9A44A; --bad:#D9705F;
  --cell:#E0968E;
}}
:root[data-theme="dark"]{
  --paper:#0F151C; --surface:#172029; --sunk:#111921;
  --ink:#E4EAF0; --ink-2:#A8B7C7; --ink-3:#77899C;
  --line:#2C3A49; --line-2:#222E3A;
--ref:#9AA7B4;
  --phi2:#79ADDD; --phi3:#E29A63;
  --good:#5FB98A; --mid:#C9A44A; --bad:#D9705F;
  --cell:#E0968E;
}
:root[data-theme="light"]{
  --paper:#E9EDF2; --surface:#FBFCFD; --sunk:#E1E7EE;
  --ink:#131A24; --ink-2:#43546A; --ink-3:#6E8098;
  --line:#C6D0DB; --line-2:#DCE3EB;
--ref:#55606D;
  --phi2:#2E5C8A; --phi3:#B4622F;
  --good:#2F7250; --mid:#8A6A18; --bad:#9E3B2E;
  --cell:#C2726C;
}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:15px;
  line-height:1.55}
.disp{font-family:"Bahnschrift","DIN Alternate","Roboto Condensed",system-ui,
  sans-serif;font-weight:600;letter-spacing:.01em}
.mono{font-family:ui-monospace,"Cascadia Mono",Consolas,"SF Mono",monospace;
  font-variant-numeric:tabular-nums}
.wrap{max-width:1180px;margin:0 auto;padding:28px 22px 60px}
header{border-bottom:2px solid var(--ink);padding-bottom:14px;margin-bottom:22px}
/*  ---- mark -------------------------------------------------------------
    Inline SVG, not a file: the page makes no external request of any kind, so
    a linked image would break it offline and from a file:// URL.

    The mark is the object the potential computes - one triplet.  An apex atom
    with two neighbours, the two legs whose lengths enter phi3 as
    x = r_ij + r_ik, and the angle at the apex that the Legendre factor acts on.
    The two neighbours ARE the initials, each leg carrying its own letter's
    colour: U in the pair colour, G in the three-body colour, the same key
    every plot below uses.

    Two earlier drafts are worth not repeating.  Discs behind the letters turned
    the whole thing into a face - two eyes and a mouth - and a wide-radius angle
    arc was the smile; the arc belongs close to the vertex, which is both
    correct and unfaceable.  And with the bonds stopping far short, the letters
    floated free of the geometry instead of being part of it.                 */
.mark{width:clamp(132px,19vw,214px);height:auto;flex:0 0 auto}
.mark .plate{fill:var(--sunk)}
/*  The cell is drawn in its own red rather than in --bad, which means
    unstable everywhere else on this page: a reader who has learned that
    red flags a problem should not meet it on the logo.                  */
.mark .edge{stroke:var(--cell);fill:none;stroke-linecap:round}
.mark .legU{stroke:var(--phi2);stroke-width:3;stroke-linecap:round;fill:none}
.mark .legG{stroke:var(--phi3);stroke-width:3;stroke-linecap:round;fill:none}
/*  Arc and label are deliberately NOT the same colour, and which gets which
    follows from what each one is.  The arc spans from the U leg to the G
    leg, so wearing either leg's colour would tie the angle to one side of
    a thing that belongs to both; it is drawn as geometry, in neutral ink,
    thick enough to read at this size.  h(theta) is not geometry - it is the
    factor multiplying phi3 - so it takes the three-body colour.          */
.mark .arc{stroke:var(--ink-2);stroke-width:2.3;fill:none;stroke-linecap:round;
  stroke-linejoin:round;opacity:0.85}
/*  All three atoms in one colour, because in an elemental crystal they are one
    species: colouring the neighbours to match their letters would say the
    triplet is made of two different elements, which it is not.  The legs
    already carry the letter colours, so U and G stay attached to their own
    sides.  The apex is drawn a little larger only because it is the vertex the
    angle belongs to.                                                         */
.mark .apex{fill:var(--ink)}
.mark .node{fill:var(--ink)}
/*  Handwritten, and the fallback chain matters more than usual: a script face
    that is missing degrades to whatever the system calls cursive, which on
    Windows is Comic Sans.  Segoe Script and Brush Script MT are the two that
    are actually present on the machines this is read on, so they lead.       */
.mark .lt{font-family:"Segoe Script","Brush Script MT","Bradley Hand",
  "Snell Roundhand","Apple Chancery",cursive;font-weight:700;font-size:18px;
  text-anchor:middle;dominant-baseline:central}
/*  The angle carries the one thing UG adds, so it is named.  Set in an
    italic serif rather than the script the initials use: those are a monogram,
    this is a function, and letting them share a face would blur the two.
    h(theta) and not the whole h = 1 + lam2 P2 + lam4 P4 - at 132 px the box is
    64 units across and anything longer stops being readable and starts being
    texture.                                                                 */
.mark .ang{font-family:"Cambria Math","Latin Modern Math",Georgia,
  "Times New Roman",serif;font-style:italic;font-weight:600;font-size:7px;
  fill:var(--phi3);text-anchor:middle;dominant-baseline:central;opacity:0.9}
.mark .u{fill:var(--phi2)}
.mark .g{fill:var(--phi3)}
.brand{display:flex;align-items:center;justify-content:space-between;
  gap:26px;margin-bottom:10px;flex-wrap:wrap}
.brand>div{flex:1 1 440px;min-width:0}
/*  ---- authors ---------------------------------------------------------- */
.authors{display:flex;flex-wrap:wrap;gap:4px 18px;margin:7px 0 0;
  font-size:13px;color:var(--ink-2)}
.authors a{color:var(--ink);text-decoration:none;
  border-bottom:1px solid var(--line)}
.authors a:hover{border-bottom-color:var(--phi3);color:var(--phi3)}
h1{margin:0 0 2px;text-wrap:balance}
/*  Two things on one line that are not the same thing: whose potential it is,
    and what the page is.  Sized and weighted apart so the eye takes the names
    first, joined by a dot in the three-body colour that the mark beside it and
    every plot below already use.  clamp() rather than a fixed size, because at
    30px the line wrapped mid-name on a narrow window.                        */
/*  Both halves at one size: the names and what the page is are set in the same
    type, and only weight, case and colour tell them apart.  The clamp is the
    size the second half already rendered at, so nothing on the line grew or
    shrank except the names coming down to meet it.                           */
.title{display:block;font-size:clamp(17px,2.7vw,27px);line-height:1.15}
.title .who{font-weight:700;letter-spacing:.01em}
.title .dot{color:var(--phi3);font-weight:400;margin:0 .34em;
  font-size:.72em;vertical-align:.12em}
.title .what{font-weight:400;letter-spacing:.055em;
  color:var(--ink-2);text-transform:uppercase;white-space:nowrap}
.sub{color:var(--ink-2);max-width:70ch;margin:0}
/*  A sub-heading, not a paragraph: this labels a block that belongs to the
    section above it rather than starting a new one.  Without a rule of its own
    it inherited .sub's margin:0 and read as body text.                       */
h4.sub{font-size:13.5px;font-weight:600;letter-spacing:.01em;
  margin:18px 0 8px;color:var(--ink-2)}
/*  .sub inside a heading is a footnote, not a heading: without this it
    inherits the h3 size and weight and reads as a second title. */
h3 .sub{display:block;font-size:12.5px;font-weight:400;line-height:1.45;margin-top:5px;max-width:78ch}
.eyebrow{font-size:11px;letter-spacing:.14em;text-transform:uppercase;
  color:var(--ink-3);margin:0 0 6px}
/*  ---- equations -------------------------------------------------------
    Built from spans rather than a maths library on purpose: the page makes
    no external request of any kind, so it works offline, behind a firewall
    and from a file:// URL.  A CDN script would cost all three.            */
.forms{margin:16px 0 4px;border:1px solid var(--line);background:var(--surface)}
.forms>summary{cursor:pointer;padding:11px 14px;font-size:13px;
  color:var(--ink);list-style:none;display:flex;align-items:center;gap:10px}
.forms>summary::-webkit-details-marker{display:none}
.forms>summary::after{content:"\25BE";margin-left:auto;color:var(--ink-3);
  transition:transform .15s}
.forms[open]>summary::after{transform:rotate(180deg)}
.forms .body{padding:2px 14px 14px;border-top:1px solid var(--line-2)}
.fcard{padding:11px 0;border-top:1px dashed var(--line-2)}
.fcard:first-child{border-top:0}
.tag{display:inline-block;font:600 11px/1.6 ui-monospace,Menlo,Consolas,monospace;
  letter-spacing:.06em;padding:1px 7px;border:1px solid currentColor}
.tag.au{color:var(--ink-3)} .tag.mau{color:var(--phi2)}
.tag.ug{color:var(--phi3)}
.inlinesel{font:inherit;font-size:13px;font-weight:400;margin-left:10px;
  padding:1px 6px;background:var(--surface);color:var(--ink);
  border:1px solid var(--line);vertical-align:middle}
.eq{font-family:Cambria,Georgia,'Times New Roman',serif;font-size:15px;
  text-align:center;margin:9px 0;color:var(--ink);overflow-x:auto;
  line-height:2.4}
.eq i{font-style:italic}
.fr{display:inline-block;vertical-align:-0.6em;text-align:center;margin:0 3px}
.fr>span{display:block;padding:0 5px}
.fr>span:first-child{border-bottom:1px solid currentColor}
.sg{display:inline-block;vertical-align:-0.4em;text-align:center;margin:0 3px}
.sg>b{display:block;font:400 1.45em/.85 Cambria,Georgia,serif}
.sg>u{display:block;font-size:.6em;text-decoration:none;color:var(--ink-2)}
.legend{display:flex;flex-wrap:wrap;gap:16px;margin:14px 0 4px;font-size:12.5px;
  color:var(--ink-2)}
/*  The row is a flex line, but each entry's TEXT must stay one inline run:
    a flex container blockifies its children, so <sub> inside a flex item stops
    being a subscript and lands beside the symbol as its own box - which is how
    the legend came to read "phi 2 pair" on three lines.  The marker and the
    label are the two flex items; everything inside the label is inline.       */
.legend span{display:flex;align-items:center;gap:6px}
.legend span > .lbl{display:inline}
.mk{width:11px;height:11px;border:1.5px solid currentColor;flex:none}
.mk.good{background:currentColor;color:var(--good)}
.mk.mid{background:linear-gradient(135deg,currentColor 50%,transparent 50%);
  color:var(--mid)}
.mk.bad{color:var(--bad)}
.ptable{display:grid;grid-template-columns:repeat(14,1fr);gap:5px;margin:18px 0 8px}
.cell{grid-column:var(--c);grid-row:var(--r);background:var(--surface);
  border:1px solid var(--line);padding:6px 4px 5px;text-align:center;
  cursor:pointer;position:relative;transition:border-color .12s,transform .12s;
  min-height:56px}
.cell:hover{border-color:var(--phi3);transform:translateY(-1px)}
/*  A cell showing a set the selector did not ask for, because this element
    has no such set.  Without a mark the table reads as though every colour
    came from the selected fit; the tooltip said otherwise but a tooltip is
    not the picture.  Dashed, and the tier bar dimmed.  */
.cell.fellback{border-style:dashed}
.cell.fellback .tier{opacity:.35}
.cell:focus-visible{outline:2px solid var(--focus);outline-offset:2px}
.cell[aria-pressed="true"]{border-color:var(--phi3);border-width:2px;
  background:var(--sunk)}
.cell .sym{font-size:19px;display:block;line-height:1.1}
.cell .st{font-size:9.5px;letter-spacing:.08em;color:var(--ink-3);
  text-transform:uppercase}
.cell .tier{position:absolute;top:4px;right:4px;width:8px;height:8px;
  border:1.5px solid currentColor}
.tier.good{background:currentColor;color:var(--good)}
.tier.mid{background:linear-gradient(135deg,currentColor 50%,transparent 50%);
  color:var(--mid)}
.tier.bad{color:var(--bad)}
/*  square marker, top left, so it cannot be confused with the round RMS tier
    dot at top right */
.cell .ugdot{position:absolute;top:4px;left:4px;width:7px;height:7px;
  background:var(--phi3);opacity:.85}
.panel{background:var(--surface);border:1px solid var(--line);
  padding:20px 20px 22px;margin-top:22px}
.phead{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap;
  border-bottom:1px solid var(--line-2);padding-bottom:12px;margin-bottom:16px}
.phead h2{font-size:26px;margin:0}
.phead .meta{color:var(--ink-2);font-size:13.5px}
.rms{margin-left:auto;font-size:13px;padding:2px 9px;border:1px solid currentColor}
.rms.good{color:var(--good)} .rms.mid{color:var(--mid)} .rms.bad{color:var(--bad)}
h3{font-size:11px;letter-spacing:.14em;text-transform:uppercase;
  color:var(--ink-3);margin:24px 0 9px;font-weight:600}
h3:first-of-type{margin-top:0}
.pars{display:grid;grid-template-columns:repeat(auto-fit,minmax(112px,1fr));
  gap:1px;background:var(--line-2);border:1px solid var(--line-2)}
.par{background:var(--surface);padding:8px 10px}
.par .k{font-size:11px;color:var(--ink-3);display:block}
.par .v{font-size:15px}
.par.p2{border-top:2px solid var(--phi2)}
.par.p3{border-top:2px solid var(--phi3)}
.par.flag{background:var(--sunk)}
.grid2{display:grid;grid-template-columns:1.15fr .85fr;gap:26px;align-items:start}
@media(max-width:820px){.grid2{grid-template-columns:1fr}}
canvas{width:100%;height:auto;display:block;background:var(--sunk);
  border:1px solid var(--line-2)}
.plotnote{font-size:12px;color:var(--ink-3);margin:7px 0 0}
.dim{color:var(--ink-3);font-size:11px}
.swatch{display:inline-block;width:20px;height:2px;vertical-align:middle;
  margin-right:5px}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th,td{text-align:right;padding:5px 8px;border-bottom:1px solid var(--line-2)}
th:first-child,td:first-child{text-align:left}
th{font-size:11px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600}
td.err{width:1%;white-space:nowrap}
/*  The finite-temperature table.  Narrower than the wide comparison tables
    above it because it is five rows of one number each, and the verdict row
    is set apart by a rule rather than by a colour block: the colour says
    which way it went, the rule says that it is a conclusion and not another
    measurement.  A row marked mid is one whose difference is smaller than
    the run could resolve, which is the commonest case and must not read as
    a failure.  */
table.ft{width:auto;min-width:min(100%,430px);margin:11px 0 3px;
  font-size:13px}
table.ft td{padding:5px 14px 5px 0}
table.ft td:last-child{text-align:right;padding-right:0;
  font-family:ui-monospace,Consolas,monospace;
  font-variant-numeric:tabular-nums;white-space:nowrap}
table.ft tr:last-child td{border-bottom:none;border-top:1px solid var(--ink-3);
  padding-top:8px}
table.ft tr.ok td:last-child{color:var(--good)}
table.ft tr.bad td:last-child{color:var(--bad)}
table.ft tr.warn td:last-child{color:var(--mid)}
.bar{display:inline-block;height:9px;background:var(--phi2);opacity:.5;
  vertical-align:middle;margin-left:6px}
.cols3{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:1px;background:var(--line-2);border:1px solid var(--line-2)}
.note{font-size:12px;color:var(--ink-3);margin:9px 0 0;max-width:76ch}
.warn{border-left:3px solid var(--mid);padding:9px 12px;background:var(--sunk);
  font-size:12.5px;color:var(--ink-2);margin:12px 0 0}
.gen{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:10px}
select,button{font:inherit;color:inherit;background:var(--surface);
  border:1px solid var(--line);padding:5px 10px;cursor:pointer}
select:focus-visible,button:focus-visible{outline:2px solid var(--focus);
  outline-offset:1px}
button.primary{border-color:var(--phi3);color:var(--phi3)}
button.primary:hover{background:var(--phi3);color:var(--surface)}
pre{background:var(--sunk);border:1px solid var(--line-2);padding:12px 14px;
  margin:0;overflow-x:auto;font-size:12.5px;line-height:1.5;max-height:320px}
footer{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);
  color:var(--ink-3);font-size:12.5px;max-width:80ch}
@media(prefers-reduced-motion:reduce){*{transition:none!important}}
</style>

<div class="wrap">
<header>
  <div class="brand">
    <div>
      <p class="eyebrow">&phi;<sub>2</sub> + &phi;<sub>3</sub> interatomic
        potential</p>
      <h1 class="disp title"><span class="who">U&#286;UR and
        G&Uuml;LER</span><span class="dot">&#9679;</span><span
        class="what">Interatomic Potential Library</span></h1>
      <p class="authors">
        <a href="https://avesis.gazi.edu.tr/gokay" target="_blank"
          rel="noopener noreferrer">Prof. Dr. G&ouml;kay U&#286;UR</a>
        <a href="https://avesis.gazi.edu.tr/suleugur" target="_blank"
          rel="noopener noreferrer">Prof. Dr. &#350;ule U&#286;UR</a>
        <a href="https://avesis.hacibayram.edu.tr/melek.guler" target="_blank"
          rel="noopener noreferrer">Prof. Dr. Melek G&Uuml;LER</a>
        <a href="https://avesis.hacibayram.edu.tr/guler.emre" target="_blank"
          rel="noopener noreferrer">Prof. Dr. Emre G&Uuml;LER</a>
      </p>
    </div>
    <svg class="mark" viewBox="0 0 64 64" role="img" id="mk"
         aria-label="U and G as two neighbours of the atom at the centre
                     of a body-centred cubic cell, turning in three
                     dimensions">
      <rect class="plate" x="1" y="1" width="62" height="62" rx="9"/>
      <line class="edge" id="mkE0" x1="20.55" y1="40.23" x2="18.44" y2="53.84" style="stroke-width:0.99;opacity:0.46"/>
      <line class="edge" id="mkE1" x1="20.55" y1="40.23" x2="19.52" y2="18.62" style="stroke-width:0.86;opacity:0.39"/>
      <line class="edge" id="mkE2" x1="20.55" y1="40.23" x2="43.45" y2="40.23" style="stroke-width:0.70;opacity:0.30"/>
      <line class="edge" id="mkE3" x1="18.44" y1="53.84" x2="16.97" y2="29.28" style="stroke-width:1.44;opacity:0.71"/>
      <line class="edge" id="mkE4" x1="18.44" y1="53.84" x2="45.56" y2="53.84" style="stroke-width:1.29;opacity:0.63"/>
      <line class="edge" id="mkE5" x1="19.52" y1="18.62" x2="16.97" y2="29.28" style="stroke-width:1.31;opacity:0.64"/>
      <line class="edge" id="mkE6" x1="19.52" y1="18.62" x2="44.48" y2="18.62" style="stroke-width:1.01;opacity:0.47"/>
      <line class="edge" id="mkE7" x1="16.97" y1="29.28" x2="47.03" y2="29.28" style="stroke-width:1.60;opacity:0.80"/>
      <line class="edge" id="mkE8" x1="43.45" y1="40.23" x2="45.56" y2="53.84" style="stroke-width:0.99;opacity:0.46"/>
      <line class="edge" id="mkE9" x1="43.45" y1="40.23" x2="44.48" y2="18.62" style="stroke-width:0.86;opacity:0.39"/>
      <line class="edge" id="mkE10" x1="45.56" y1="53.84" x2="47.03" y2="29.28" style="stroke-width:1.44;opacity:0.71"/>
      <line class="edge" id="mkE11" x1="44.48" y1="18.62" x2="47.03" y2="29.28" style="stroke-width:1.31;opacity:0.64"/>
      <polyline class="arc" id="mkArc" points="27.45,29.34 28.12,28.94 28.93,28.57 29.86,28.26 30.9,28.06 32.0,27.98 33.1,28.06 34.14,28.26 35.07,28.57 35.88,28.94 36.55,29.34"/>
      <line class="legU" id="mkLU" x1="32.0" y1="35.5" x2="19.52" y2="18.62" style="stroke-width:2.69"/>
      <line class="legG" id="mkLG" x1="32.0" y1="35.5" x2="44.48" y2="18.62" style="stroke-width:2.69"/>
      <circle class="apex" id="mkA" cx="32.0" cy="35.5" r="3.90"/>
      <circle class="node" id="mkJ" cx="19.52" cy="18.62" r="2.88"/>
      <circle class="node" id="mkK" cx="44.48" cy="18.62" r="2.88"/>
      <text class="lt u" id="mkU" x="15.52" y="13.21" style="font-size:13.0px">U</text>
      <text class="lt g" id="mkG" x="48.48" y="13.21" style="font-size:13.0px">G</text>
      <text class="ang" id="mkH" x="32.0" y="20.6">h(&#952;)</text>
    </svg>
  </div>
  <p class="sub">A two-body plus three-body potential for metals, refitted with
  a continuous r-power and evaluated entirely from the potential's own
  derivatives. <strong>__NOK__ of the __NTRIED__ metals attempted</strong> yield
  a fit. Cohesive
  energy, zero pressure and bulk modulus are exact constraints, so every elastic
  constant, phonon and thermodynamic quantity below is a prediction rather than
  a fitted quantity.</p>
  <p class="sub" style="margin-top:9px;font-size:13.5px">
  The functional form is that of
  <strong>&#304;. Akg&uuml;n and G. U&#287;ur</strong>,
  <i>Phys. Rev. B</i> <strong>51</strong>, 3458 (1995);
  <i>Nuovo Cimento D</i> <strong>19</strong>, 779 (1997);
  <i>Nuovo Cimento D</i> <strong>20</strong>, 1549 (1998) &mdash; the last being
  the five-parameter version carrying the (r<sub>0</sub>/r)<sup>&gamma;</sup>
  prefactor used here. The parameters on this page are a fresh
  fit, not the published ones.</p>
  <details class="forms">
    <summary><strong>The three forms</strong>
      <span style="color:var(--ink-3)">AU &sub; MAU &sub; UG &mdash; each
        contains the previous one exactly</span></summary>
    <div class="body">

      <div class="fcard">
        <span class="tag au">AU</span>
        <strong style="margin-left:8px">Akg&#252;n&ndash;U&#287;ur</strong>
        <span style="color:var(--ink-3);font-size:12.5px">&mdash; the published
          form; <i>D</i>, <i>C</i>, <i>r</i><sub>0</sub>, &alpha;</span>
        <div class="eq">
          &phi;<sub>2</sub>(<i>r<sub>ij</sub></i>) =
          <span class="fr"><span><i>D</i></span><span>2(<i>m</i>&minus;1)</span></span>
          <span class="sg"><b>&Sigma;</b><u><i>i</i>&ne;<i>j</i></u></span>
          (<i>r</i><sub>0</sub>/<i>r<sub>ij</sub></i>)<sup>&gamma;</sup>
          [ e<sup><i>m</i>&alpha;(<i>r</i><sub>0</sub>&minus;<i>r<sub>ij</sub></i>)</sup>
          &minus; <i>m</i>&thinsp;e<sup>&alpha;(<i>r</i><sub>0</sub>&minus;<i>r<sub>ij</sub></i>)</sup> ]
        </div>
        <div class="eq">
          &phi;<sub>3</sub> =
          <span class="fr"><span><i>CD</i></span><span>2(<i>m</i>&minus;1)</span></span>
          <span class="sg"><b>&Sigma;</b><u><i>j</i>&ne;<i>k</i></u></span>
          <span class="sg"><b>&Sigma;</b><u><i>i</i></u></span>
          (<i>r</i><sub>0</sub>/<i>x</i>)<sup>&gamma;</sup>
          [ e<sup><i>m</i>&alpha;(<i>r</i><sub>0</sub>&minus;<i>x</i>)</sup>
          &minus; <i>m</i>&thinsp;e<sup>&alpha;(<i>r</i><sub>0</sub>&minus;<i>x</i>)</sup> ],
          &nbsp; <i>x</i> = <i>r<sub>ij</sub></i> + <i>r<sub>ik</sub></i>
        </div>
        <p class="note">&phi;<sub>3</sub> depends on <i>r<sub>ij</sub></i> +
          <i>r<sub>ik</sub></i> only, so it cannot see the angle at the central
          atom. That single fact sets what the form can and cannot reach.</p>
      </div>

      <div class="fcard">
        <span class="tag mau">MAU</span>
        <strong style="margin-left:8px">modified Akg&#252;n&ndash;U&#287;ur</strong>
        <span style="color:var(--ink-3);font-size:12.5px">&mdash; this library;
          adds &alpha;<sub>3</sub></span>
        <div class="eq">
          &phi;<sub>3</sub> =
          <span class="fr"><span><i>CD</i></span><span>2(<i>m</i>&minus;1)</span></span>
          <span class="sg"><b>&Sigma;</b><u><i>j</i>&ne;<i>k</i></u></span>
          <span class="sg"><b>&Sigma;</b><u><i>i</i></u></span>
          (<i>r</i><sub>0</sub>/<i>x</i>)<sup>&gamma;</sup>
          [ e<sup><i>m</i>&alpha;<sub>3</sub>(<i>r</i><sub>0</sub>&minus;<i>x</i>)</sup>
          &minus; <i>m</i>&thinsp;e<sup>&alpha;<sub>3</sub>(<i>r</i><sub>0</sub>&minus;<i>x</i>)</sup> ]
        </div>
        <p class="note">The three-body term gets its own decay constant, searched
          through <i>s</i><sub>3</sub> = &alpha;<sub>3</sub>/&alpha;.
          <i>s</i><sub>3</sub> = 1 recovers AU exactly. &phi;<sub>2</sub> is
          unchanged.</p>
      </div>

      <div class="fcard">
        <span class="tag ug">UG</span>
        <strong style="margin-left:8px">U&#287;ur&ndash;G&uuml;ler</strong>
        <span style="color:var(--ink-3);font-size:12.5px">&mdash; angular
          generalisation; adds &lambda;<sub>2</sub>, &lambda;<sub>4</sub></span>
        <div class="eq">
          &phi;<sub>3</sub><sup>UG</sup> =
          &phi;<sub>3</sub><sup>MAU</sup> &middot;
          <i>h</i>(cos&thinsp;&theta;<sub><i>jik</i></sub>)
        </div>
        <div class="eq">
          <i>h</i>(cos&thinsp;&theta;) = 1 +
          &lambda;<sub>2</sub><i>P</i><sub>2</sub>(cos&thinsp;&theta;) +
          &lambda;<sub>4</sub><i>P</i><sub>4</sub>(cos&thinsp;&theta;)
        </div>
        <p class="note">&theta;<sub><i>jik</i></sub> is the angle between the two
          legs at the central atom. <i>P</i><sub>2</sub> and
          <i>P</i><sub>4</sub> average to zero over the sphere, so the cohesive
          energy, pressure and bulk modulus &mdash; imposed exactly as
          constraints &mdash; are untouched and only the anisotropy moves.
          &lambda;<sub>2</sub> = &lambda;<sub>4</sub> = 0 recovers MAU.</p>
      </div>

    </div>
  </details>

  <div class="legend">
    <span><i class="mk good"></i><span class="lbl">RMS &le; 20%</span></span>
    <span><i class="mk mid"></i><span class="lbl">20-35%</span></span>
    <span><i class="mk bad"></i><span class="lbl">&gt; 35% &mdash; refit before use</span></span>
    <span><span class="lbl" style="color:var(--ink-3)">&mdash; of the parameter
      set selected below; dashed cell = that element has no such set, so its
      colour is MAU</span></span>
    <span><i class="mk" style="background:var(--phi3);border-color:var(--phi3);
      border-radius:0"></i><span class="lbl">UG fit available &mdash; __NUG__
      metals carry the angular term beside MAU</span></span>
    <span style="margin-left:auto"><i class="swatch"
      style="background:var(--phi2)"></i><span class="lbl">&phi;<sub>2</sub>
      pair</span></span>
    <span><i class="swatch" style="background:var(--phi3)"></i><span
      class="lbl">&phi;<sub>3</sub> three-body</span></span>
  </div>
</header>

<div class="ptable" id="pt" role="group" aria-label="Element selection"></div>
<div id="plotfail" style="display:none;margin:10px 0;padding:8px 10px;border:1px solid var(--bad);color:var(--bad);font-size:13px"></div>
<div class="panel" id="panel"></div>

<footer>
  <strong style="color:var(--ink)">A separate study</strong><br>
  <a href="forms.html">forms.html</a> compares the published form with two
  variants that replace or extend its three-body term by an embedding, on five
  bcc metals and on their alloys. Those variants are not offered for use and are
  not in the menu above; that page says why, and what they do and do not show.
  <br><br>
  <strong style="color:var(--ink)">How to cite</strong><br>
  The <em>functional form</em> is not ours to claim &mdash; it is
  &#304;. Akg&uuml;n and G. U&#287;ur, <i>Phys. Rev. B</i>
  <strong>51</strong>, 3458 (1995); <i>Nuovo Cimento D</i> <strong>19</strong>,
  779 (1997); <i>Nuovo Cimento D</i> <strong>20</strong>, 1549 (1998), the last
  being the five-parameter version used here. Cite those for the potential.
  The <em>parameters on this page</em> are a fresh fit and are not the published
  ones; cite this page for them, and please say which build you used, since the
  numbers move as the fits improve.
  <br><br>
  <strong style="color:var(--ink)">Acknowledgement</strong><br>
  The numerical calculations reported in this work were partially performed at
  T&Uuml;B&#304;TAK ULAKB&#304;M, High Performance and Grid Computing Center
  (TRUBA resources).
  <br><br>
  The reference data are other people's measurements and should be cited as
  theirs, not as ours: elastic constants from
  <strong>Landolt-B&ouml;rnstein III/29a</strong> (Every and McCurdy, Springer
  1992), except beryllium from <strong>Migliori <i>et al.</i></strong>,
  <i>J. Appl. Phys.</i> <strong>95</strong>, 2436 (2004); cohesive energies from
  <strong>Brewer</strong>, LBL-3720 Rev. (1977); measured phonon frequencies
  from <strong>Landolt-B&ouml;rnstein III/13a</strong> (Schober and Dederichs,
  Springer 1981); standard entropies and heat capacities at 298.15 K from the
  <strong>CRC Handbook of Chemistry and Physics</strong> thermodynamic tables;
  Debye temperatures from <strong>G. R. Stewart</strong>,
  <i>Rev. Sci. Instrum.</i> <strong>54</strong>, 1 (1983), Table I, whose own
  footnote states they are the low-temperature values,
  T&nbsp;&Lt;&nbsp;&theta;<sub>D</sub>; linear thermal expansion at 25&nbsp;&deg;C
  from the same CRC Handbook, whose expansion column is
  <strong>Touloukian</strong>, <i>Thermophysical Properties of Matter</i>,
  vol. 12. The calculated dispersions drawn for comparison come from the
  <strong>Materials Project</strong> (Jain <i>et al.</i>, <i>APL Materials</i>
  <strong>1</strong>, 011002 (2013); the phonon entries are computed with
  <i>pheasy</i> and are labelled with that method in the plot selector), the
  <strong>Materials Cloud</strong> archive (Talirz <i>et al.</i>,
  <i>Sci. Data</i> <strong>7</strong>, 299 (2020)) &mdash; specifically its
  <code>supercon_phonon-vis</code> contribution under MC3D, PBEsol dispersions
  from an electron&ndash;phonon study rather than from MC3D's own structural
  data &mdash; and <strong>JARVIS-DFT</strong> (Choudhary <i>et al.</i>,
  <i>npj Comput. Mater.</i> <strong>6</strong>, 173 (2020)), and carry their own
  citation requirements. None of the three is treated as an arbiter: they are
  independent calculations, they disagree with each other by up to a factor of
  two on some elements, and where one of them conflicts with a measurement by
  more than 25 % it is withheld from the plot rather than drawn.
  <br><br>
  Everything shown here is computed from the potential's own derivatives; no
  external code enters the chain. What that means in practice is set out below.
  <br><br>
  Elastic constants include the non-affine internal-strain correction, which for
  a two-atom cell such as hcp is not a small effect. Phonons come from the
  dynamical matrix built from the same second derivatives; the force constants
  are checked against the acoustic sum rule and Hermitian symmetry on every
  build. Nothing thermal enters the fit, so the entropy and heat capacity shown
  are predictions. Elastic constants for the cubic elements are taken from
  <strong>Landolt-B&ouml;rnstein New Series III/29a</strong> (Every and
  McCurdy, 1992), Table 3, which reports a weighted mean over the published
  measurements; where that table lists two determinations the one consistent
  with the measured bulk modulus was kept, which is what decides barium
  (its other entry gives B = 2.45 GPa against a measured 9.5). Cohesive
  hexagonal constants come from Table 11 of the same volume. Beryllium is the
  one exception: III/29a gives C<sub>13</sub> = 6 &plusmn; 9, an outlier against
  14 from resonant ultrasound (Migliori <i>et al.</i>, <i>J. Appl. Phys.</i>
  <strong>95</strong>, 2436 (2004)) and 19.1 from DFT, and the same
  bulk-modulus test decides it. Cohesive energies are Brewer,
  <i>The Cohesive Energies of the Elements</i>, LBL-3720 Rev. (1977), Table I,
  the source Kittel's table is drawn from. Measured phonon frequencies, where shown,
  come from <strong>Landolt-B&ouml;rnstein New Series III/13a</strong> (Springer
  1981), the standard compilation of inelastic-neutron-scattering dispersion
  data for the elements, read from the tabulation in Savrasov and Savrasov,
  <i>Phys. Rev. B</i> <strong>54</strong>, 16487 (1996). Using one compiler for
  every element keeps the measurement conventions consistent; it also limits the
  comparison to the seven metals that tabulation covers.

  <br><br><strong>The three anchors are not at one temperature, and the fit
  that uses them is static.</strong> Cohesive energies are 0 K; lattice
  constants are room temperature (except <b>lithium at 78 K and rubidium and
  caesium at 5 K</b> &mdash; sodium and potassium are room temperature like the
  rest, which this page said otherwise until 2026-08-29); and III/29a states
  its own convention plainly &mdash; <i>&ldquo;unless
  otherwise stated, all elastic constants are given at room temperature, RT,
  (= 300 K)&rdquo;</i>. So a zero-kelvin lattice sum is being matched to
  room-temperature stiffnesses at a room-temperature volume. The library is
  therefore built a little too soft and a little too expanded against a true
  0 K reference, and how little depends strongly on the element: III/29a's own
  temperature coefficients put tungsten's C<sub>44</sub> only 2 % higher at 0 K,
  copper's 11 %, aluminium's 15 % &mdash; and sodium's <strong>67 %</strong>.
  C<sub>44</sub> moves about twice as far as C<sub>11</sub> everywhere, so the
  anisotropy the fit is asked to reach is distorted as well as the magnitude.
  Three elements carry an explicit temperature in Table 3: lithium's
  room-temperature row is the one used, while rubidium's &asymp;80 K and
  caesium's 78 K rows are &mdash; which, with their lattice constants genuinely
  at 5 K, makes them by accident the two most internally consistent records
  here.
  <br><br><strong>And the fit was tried against the corrected anchors, in
  2026.</strong> All three were carried to 0 K &mdash; the elastic constants by
  III/29a's Tables 28 and 35, the lattice constant by a Debye-shaped expansion
  integral, the cohesive energy by (9/8)&nbsp;k<sub>B</sub>&theta;<sub>D</sub>
  &mdash; and four elements were refitted to them and scored against the
  measured dispersion, which never enters the fit. Two improved, six got worse,
  and the control got worse in both cutoff treatments. <b>So the correction is
  real and it is not a lever</b>: read the percentages above as error bars on
  how far each target sits from the 0 K quantity the static fit computes, not
  as a better place to aim. Nothing has been changed on this account.

  <br><br><strong>A second systematic sits beside it and is larger for the
  light elements.</strong> A measured cohesive energy is the work to take the
  crystal apart <em>from its zero-point state</em>; a classical potential has
  no zero-point motion, so its static well should be deeper than the measured
  value by exactly that energy. Every element here is therefore fitted a little
  too shallow &mdash; by about 4.2 % for beryllium, 2.6 % for magnesium and
  2.0 % for lithium, against 0.3 % for tantalum. Largest, again, where the
  library already struggles. The Debye temperatures behind those figures are
  not in the reference data; they were entered by hand to size the effect, so
  the ordering is sound and the third digit is not.
  <p class="plotnote" style="margin-top:14px;opacity:.7">
  <b>__TREE__ tree</b> &middot; built __BUILT__ &middot; __NSET__ elements
  with more than one selectable parameter set. Two trees produce a page of
  this name and, since 1.1.0, they carry the same panels &mdash; the screen
  tree's renderer and library were promoted into the published one for the
  release &mdash; so this line says only which tree wrote the file you have
  open. If the timestamp is older than the change you are looking for, the
  browser is showing a cached copy &mdash; reload with the cache
  bypassed.</p>
</footer>
</div>

<script>
/*  A runtime error in this script is otherwise invisible: the browser writes
    it to a console nobody has open, and the page just stops - half-drawn, or
    with the logo frozen because the animation never got its first frame.
    Anything thrown from here on says so on the page itself.  The esprima
    check in make_gui.py catches SYNTAX errors before the file is written;
    this is the other half.  */
addEventListener("error", e => {
  let b = document.getElementById("jserr");
  if(!b){
    b = document.createElement("div");
    b.id = "jserr";
    b.style.cssText = "position:fixed;left:0;right:0;bottom:0;z-index:9999;"
      + "background:#a00;color:#fff;font:12px ui-monospace,monospace;"
      + "padding:8px 12px;white-space:pre-wrap;max-height:40vh;overflow:auto";
    (document.body||document.documentElement).appendChild(b);
  }
  b.textContent += `script error: ${e.message}
    at ${e.filename||"page"}`
                 + `:${e.lineno}:${e.colno}
`;
});

const DATA = __DATA__;
const FT = __FT__;
const $ = (s,r=document)=>r.querySelector(s);
const fmt=(x,n=4)=>(x===null||x===undefined||isNaN(x))?"&mdash;":Number(x).toFixed(n);
const tier=r=>r<=20?"good":r<=35?"mid":"bad";
/*  A measured value printed exactly as it is stored, with no padding and no
    rounding.  The two places that show the 298 K experiments used toFixed(1)
    and toFixed(2), so the same number appeared as 28.3 in one table and 28.30
    in the other, and 34.6 was padded to 34.60 as though a digit had been
    measured that was not.  The table's own precision varies - the older rows
    carry two decimals and the ones read out of the CRC volume carry one - and
    showing that is the point, not hiding it. */
const exact=v=>(v===null||v===undefined||isNaN(v))?"&mdash;":String(v);

function fshape(r,m,al,r0,g){
  const u=al*(r0-r);
  return Math.pow(r0/r,g)*(Math.exp(m*u)-m*Math.exp(u))/(m-1);
}
const phi2=(r,d)=>d.D*fshape(r,d.m,d.alpha,d.r0,d.gamma);
const phi3=(x,d)=>d.C*d.D*fshape(x,d.m,d.alpha3,d.r0,d.gamma);
/*  The angular factor of UG, h = 1 + lam2 P2 + lam4 P4, swept over the apex
    angle.  Returns [min, max] - and that is why UG's phi3 cannot be drawn as a
    curve the way MAU's can: at one leg-sum x it is a whole interval, one value
    per geometry.  Drawing a single line for it would be a picture of a
    potential that does not exist. */
function hRange(u){
  let lo=Infinity, hi=-Infinity;
  for(let i=0;i<=200;i++){
    const c1=-1+2*i/200;
    const h=1+u.lam2*0.5*(3*c1*c1-1)+u.lam4*0.125*(35*Math.pow(c1,4)-30*c1*c1+3);
    if(h<lo)lo=h; if(h>hi)hi=h;
  }
  return [lo,hi];
}

/*  Opens on palladium, and it is worth saying why it is not iron any more.

    Iron was chosen when the only test that mattered was the fit, and the
    reason given was that palladium had no angular fit so the MAU/UG comparison
    was invisible on the first screen.  Both halves of that have expired: the
    angular arm now covers all thirty-eight elements, and iron has since failed
    three of the tests added afterwards - it does not hold its lattice against
    a 1e-5 A displacement, its C11-C12 goes negative at five per cent of the
    melting point, and its ground state is hcp rather than bcc.  Opening on it
    showed the library at its worst by accident.

    Palladium is the record that comes through everything: 0.00 per cent in
    both arms, holds its lattice along all five displacement directions, no
    Born violation at any temperature, a tetragonal well that never turns over,
    and - alone with platinum, rhenium and lithium - a ground state the angular
    term gets right.  It also has phonons and measured S and Cp, so no view on
    the page is empty.

    What it is not best at, so that this is a choice and not a claim: it has
    one published potential to be compared against where copper has five, and
    the smallest step at the cutoff in the library, so the drift measurement
    that argues for the taper is least dramatic here of anywhere.  Both of
    those sections sit far below the fold. */
let cur="Pd";

/*  This block must stay ABOVE the table: buildTable() runs as soon as it is
    defined and reads PAR_SET, and `let` in the temporal dead zone is a
    ReferenceError, not undefined.  Declared 600 lines lower down it killed the
    whole script at load - no table, no plots, nothing - because everything
    after the throw never ran.  */
/*  Which parameter set the panels describe.  This was a boolean while there
    were two sets; it is a key so that adding one is a data change.  PAR_UG
    is kept as a derived value because several panels ask only "is this the
    angular arm", which is true of both UG sets.  */
let PAR_SET = "mau";
/*  Four selectable parameter sets, not two.  The re-cut candidates are full
    records - their own m, D, alpha, gamma, Cij, dispersion, expansion and
    screens - so the whole panel can stand on one of them, which is the only
    way to read the candidate against the published arm rather than beside
    fragments of it.  */
/*  THE SELECTOR IS A DIFFERENT AXIS FROM ARMS, and adding an arm to one
    does not add it to the other: ARMS decides which CURVES a panel draws,
    SETKEY and setsOf decide which PARAMETER SET the tables show.  Note also
    that `mau` here resolves to `d`, the TOP-LEVEL hard-cut record, not to
    d.tap - the selector's "MAU" and the ARMS entry "tap" are different
    records.  */
const SETKEY = {mau:d=>d, ug:d=>d.ug, rc:d=>d.rc, rc_ug:d=>d.rc_ug,
                tap_force:d=>d.tap_force};
const parSet = d => ((SETKEY[PAR_SET]||SETKEY.mau)(d)) || d;

Object.defineProperty(window,"PAR_UG",{get:()=>PAR_SET==="ug"});

const pt=$("#pt");
/*  All four sets of every element, in display order, for whichever of them
    that element actually has.  */
/*  The one place the arms are named.  Every plot and every table used to
    carry its own [["tap",d.tap],["tap_ug",d.tap_ug]], written out by hand in
    more than eighty places, so adding a set meant editing all of them and
    forgetting one was a certainty rather than a risk - the same "changed it
    here, left it there" failure that put a stale caption under a marker whose
    meaning had changed.  Adding a set is now a data question: give an element
    the key and it appears wherever it has something to show.
    Solid line = no angular term, dashed = angular.  */
/*  rc and rc_ug are the SHELL-GAP RE-CUT candidates, not part of the
    published library.  They are here because the comparison is the whole
    point: they pass every cold screen - lattice constant, cohesive energy,
    C11, C12, C44, B, and dynamical stability along the whole symmetry path,
    more of it than the published arm - and then fail the warm ones for a
    specific set of elements.  Keeping them off the page would leave that
    finding invisible.  */
/*  tap_force is the FORCE-MATCHED arm: the same functional form fitted to
    DFT forces and energies while the experimental anchors are held.  It is
    not a better fit of the same thing - it is a different reference frame,
    which is why it sits beside the others rather than replacing them, and
    why `refit/ARM_NAMING.md` refused the name `tap_dft`.  Seven of the
    thirty-eight elements have one.  */
const ARMS = [["tap", "MAU"], ["tap_ug", "UG"],
              ["rc", "re-cut"], ["rc_ug", "re-cut UG"],
              ["tap_force", "force-matched"]];
const ARM_DASH = {tap: [], tap_ug: [6, 3], rc: [2, 2], rc_ug: [5, 2, 1, 2],
                  tap_force: [1, 3]};
const ARM_ANG  = {tap: false, tap_ug: true, rc: false, rc_ug: true,
                  tap_force: false};
/*  Four of the forty-five candidate records failed selection outright -
    Mo and W in both arms, with a mode near -13 cm^-1 on the symmetry path and
    molybdenum also putting fcc 209 meV/atom below bcc.  They are kept as the
    control the accepted ones are read against, and they say so in their own
    label, because a rejected fit drawn beside an accepted one with the same
    styling is a trap.  */
/*  Rejection means DYNAMICALLY unstable, and nothing else.  Testing
    ground.ok as well looked equivalent - the four rejected candidates fail
    both - and is not: 72 of the 76 SHIPPED tapered records also predict the
    wrong ground state, which is a known finding of this library and not a
    reason to stamp them rejected.  The two criteria agree only on the
    candidates, which is exactly the case that hides the difference.  */
const armBad = r => !!(r && r.stable === false);
/*  A force-matched record can fail for a DIFFERENT reason, and the two must
    not wear the same word.  `armBad` means the crystal is dynamically
    unstable.  Over the ceiling means `standalone/fit.py`'s acceptance test -
    a strict three-body bound with no tolerance - would not have kept the fit;
    four of the seven sit 0.0002-0.0013 eV/atom above it while niobium,
    tantalum and vanadium sit strictly inside with a violation of exactly
    zero.  The record carries the verdict in `gate`, precomputed, so this page
    does not reimplement the rule and cannot drift from it.  */
const armOverCeiling = r => !!(r && r.gate === "reject");
/*  An unstable record is not always a REJECTED one.  REJECTED is a selection
    verdict and only the re-cut candidates went through selection; the four
    that failed are kept as controls.  A PUBLISHED arm with an imaginary mode
    was shipped, and its potential file says so in a LATTICE STABILITY
    warning - calling it rejected would conflate "never selected" with
    "shipped with a known defect".  Since 2026-09-13 `stable` means mesh AND
    path on every record (unify_stable.py); before that the labels read a
    mesh-only verdict on 16 records and said nothing while the red box below
    them said dynamically unstable.  ON THE PATH when the mesh alone is clean,
    because that is the case a reader would otherwise not believe.  */
const CANDIDATE_ARMS = ["rc", "rc_ug"];
const armStab = (r, k) => !armBad(r) ? ""
  : CANDIDATE_ARMS.indexOf(k) >= 0 ? " — REJECTED"
  : ((r.dyn||{}).min_mesh_cm1 !== undefined && r.dyn.min_mesh_cm1 >= -1)
    ? " — UNSTABLE ON THE PATH" : " — DYNAMICALLY UNSTABLE";
const armMark = (r, k) => armBad(r) ? armStab(r, k)
                   : armOverCeiling(r) ? " — OVER CEILING" : "";
/*  [key, record, label] for every arm this element has that satisfies `has` */
const armRecs = (d, has) =>
  ARMS.map(a => [a[0], d[a[0]], a[1] + armMark(d[a[0]], a[0])])
      .filter(o => o[1] && has(o[1]));

/*  ARM_DASH in words, for the caption.  The two have to agree, so they sit
    together: a pattern here that no longer matches the one the plot uses is
    worse than no legend at all.  */
const ARM_STYLE = {tap:"solid", tap_ug:"dashed",
                   rc:"dotted", rc_ug:"dash-dot", tap_force:"fine dots"};

/*  The legend for a panel, built from the SAME armRecs call that panel draws
    from.  Written by hand it named MAU and UG in four captions while the
    plots drew up to four curves, which is the failure the ARMS table exists
    to prevent.  */
const armKey = recs => recs.map(o =>
  `<i class="swatch" style="background:var(--${ARM_ANG[o[0]]?"phi3":"phi2"})"></i>`
  + o[2] + ` (${ARM_STYLE[o[0]]||"solid"})`).join("&nbsp;&middot;&nbsp;");

const setsOf = d => [["mau","MAU",d],["ug","UG",d.ug],
                     ["rc","re-cut",d.rc],["rc_ug","re-cut UG",d.rc_ug],
                     ["tap_force","force-matched",d.tap_force]]
                    /*  the rms filter is why a new arm needs a NUMERIC rms:
                        without one the record is dropped here silently, the
                        page still builds, and the arm simply never appears in
                        the menu.  */
                    .filter(o=>o[2]&&typeof o[2].rms==="number");

/*  The table used to colour every cell by d.rms - the hard-cut MAU fit -
    whatever set was selected, and it was built once so it could not follow a
    change.  That put the library's WORST arm on the front page as though it
    were the library: iridium read 10.65 % while its tapered and angular sets
    are at or under 6.26, and niobium read 20.35 against a tapered 0.00.  It
    is the same defect as a stability badge that
    disagrees with the dispersion drawn beside it - a verdict measured from
    something other than what is on screen - so the cell now follows the
    selection, and falls back to MAU for an element that has no such set.  The
    title carries every set the element has, because one number cannot say
    that the spread is a factor of six.  */
function buildTable(){
  pt.innerHTML="";
  Object.keys(DATA).sort((a,b)=>DATA[a].pos[0]-DATA[b].pos[0]
    ||DATA[a].pos[1]-DATA[b].pos[1]).forEach(el=>{
    const d=DATA[el], b=document.createElement("button");
    b.className="cell"; b.style.setProperty("--r",d.pos[0]-1);
    b.style.setProperty("--c",d.pos[1]);
    b.setAttribute("aria-pressed",el===cur); b.dataset.el=el;
    const sets=setsOf(d);
    const shown=(sets.find(o=>o[0]===PAR_SET)||sets[0]);
    /*  The warning triangle has to belong to the set the cell is showing.  It
        read the root record whatever was selected, so a candidate arm that is
        unstable where the published one is not - or the other way round -
        went unmarked.  */
    const unst=((shown[2].dyn)||d.dyn||{}).stable===false;
    const rms=shown[2].rms, fellBack=shown[0]!==PAR_SET;
    if(fellBack) b.classList.add("fellback");
    /*  The dot marks "has a UG fit".  That was informative when fourteen
        elements had one; all thirty-eight do now, so it marks nothing.  Kept
        rather than removed because a legend entry that stops appearing is a
        change with no reader behind it.  */
    const hasUG=!!d.ug;
    b.innerHTML=`<i class="tier ${tier(rms)}"></i>
      <span class="sym disp">${el}${unst?'<sup style="color:var(--bad)">&#9888;</sup>':""}</span>
      ${hasUG?'<i class="ugdot" title="UG comparison available"></i>':""}
      <span class="st">${d.struct}</span>`;
    /*  Newlines via fromCharCode, not an escape: this JS lives inside a
        Python string, so a backslash-n written here is turned into a real
        line break before it ever reaches the browser and the string literal
        it sits in ends early.  export_potentials.py uses chr(10) for the same
        reason.  */
    const NL=String.fromCharCode(10);
    b.title=`${d.name} - ${d.struct}`+NL
      +sets.map(o=>`${o[1]}: ${o[2].rms.toFixed(2)} %`).join(NL)
      +(fellBack?NL+`(no ${PAR_SET} set - showing ${shown[1]})`:"")
      +(unst?NL+`dynamically unstable: ${
          ((((shown[2].dyn)||d.dyn||{}).imag_frac||0)*100).toFixed(1)
        } % imaginary modes`:"");
    b.onclick=()=>{cur=el;
      pt.querySelectorAll(".cell").forEach(c=>
        c.setAttribute("aria-pressed",c.dataset.el===el));
      render();};
    pt.appendChild(b);
  });
}
buildTable();

function setup(id,ratio=0.62){
  const cv=$(id); if(!cv) return null;
  const dpr=devicePixelRatio||1, W=cv.clientWidth, H=Math.round(W*ratio);
  cv.width=W*dpr; cv.height=H*dpr; cv.style.height=H+"px";
  const c=cv.getContext("2d"); c.scale(dpr,dpr); c.clearRect(0,0,W,H);
  const cs=getComputedStyle(document.documentElement);
  return {c,W,H,
    p2:cs.getPropertyValue("--phi2").trim(), p3:cs.getPropertyValue("--phi3").trim(),
    lin:cs.getPropertyValue("--line-2").trim(),
    ink:cs.getPropertyValue("--ink-3").trim()};
}

/*  A solid frame with inward ticks on all four sides.

    The dispersion panel had one and the other two did not, so the potential
    curve and the thermodynamics floated on bare gridlines with no axes at all.
    Shared rather than copied, because three plots drifting apart in style is
    how that happened.                                                        */
function frame(c, L, T, pw, ph, ink, xt, yt){
  c.strokeStyle=ink; c.lineWidth=1.2;
  c.strokeRect(L, T, pw, ph);
  (xt||[]).forEach(x=>{
    c.beginPath(); c.moveTo(x, T+ph); c.lineTo(x, T+ph-5); c.stroke();
    c.beginPath(); c.moveTo(x, T);    c.lineTo(x, T+5);    c.stroke();});
  (yt||[]).forEach(y=>{
    c.beginPath(); c.moveTo(L, y);      c.lineTo(L+5, y);      c.stroke();
    c.beginPath(); c.moveTo(L+pw, y);   c.lineTo(L+pw-5, y);   c.stroke();});
}


/*  Elastic constants and the mechanical properties derived from them, against
    temperature.  The physics is done in add_elastic_T.py and this only draws:
    Voigt-Reuss-Hill has a different closed form for cubic and hexagonal, and a
    template literal inside an HTML file is the worst place to keep a formula
    that cannot be unit-tested.

    Three things are drawn differently on purpose, because they change what a
    curve means:
      - a point above the melting point is hollow.  A perfect small crystal has
        nowhere to nucleate from and superheats happily, so those describe a
        metastable solid.  Two independent signals agree with the flag: the
        hexagonal identity C66 = (C11-C12)/2 degrades from 1.4 % to 67 % there,
        and every Born-criterion violation in the set sits above melting.
      - a series whose lattice fails the nudge test is drawn in the warning
        colour and says so.  Thermal motion is ten thousand times the 1e-5 A
        that already destroys it, so the curve is a property of some other
        structure.  Hiding it would be worse; drawing it as if it were sound
        would be dishonest.
      - published potentials are grey.  They are the calibration, not the
        subject, and they should not compete for attention with the two sets
        this library ships.                                                   */
let ET_Q = "C11";
const ET_QS = [["C11","C11"],["C12","C12"],["C13","C13"],["C33","C33"],
               ["C44","C44"],["C66","C66"],["B","B"],["G","G"],["E","E"],
               ["nu","ν"],["BG","B/G"],["AU","Aᵁ"],["Hv","Hᵥ"]];
const ET_UNIT = q => (q==="nu"||q==="BG"||q==="AU") ? "" : " (GPa)";

function drawElasticT(d){
  const et=d.elasticT; if(!et) return;
  const s=setup("#elasT",0.60); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const cs=getComputedStyle(document.documentElement);
  const bad=cs.getPropertyValue("--bad").trim();
  const ref=cs.getPropertyValue("--ref").trim();
  //  R is wide because the curve labels live in that margin
  const L=58,R=104,T=12,B=30, pw=W-L-R, ph=H-T-B;

  /*  hexagonal-only constants simply do not exist for a cubic element, and an
      empty frame is more confusing than a missing option                     */
  const keys=Object.keys(et);
  const series=keys.map(k=>({k, r:et[k]}))
    .filter(o=>o.r.pts.some(q=>q[ET_Q]!==null&&q[ET_Q]!==undefined&&isFinite(q[ET_Q])));
  if(!series.length){ c.fillStyle=ink; c.font="12px system-ui";
    c.fillText("no data for this quantity",L,T+ph/2); return; }

  let tmax=0, lo=Infinity, hi=-Infinity;
  series.forEach(o=>o.r.pts.forEach(q=>{
    const v=q[ET_Q]; if(v===null||v===undefined||!isFinite(v)) return;
    tmax=Math.max(tmax,q.T); lo=Math.min(lo,v); hi=Math.max(hi,v);}));
  if(!(hi>lo)){ hi=lo+1; }
  const pad=(hi-lo)*0.12; lo-=pad; hi+=pad;
  if(lo>0&&lo<(hi-lo)) lo=0;
  const X=t=>L+(tmax?t/tmax:0)*pw, Y=v=>T+ph-(v-lo)/(hi-lo)*ph;

  c.strokeStyle=lin;c.lineWidth=1;
  c.font="11px ui-monospace,Consolas,monospace";c.fillStyle=ink;
  for(let k=0;k<=4;k++){const v=lo+(hi-lo)*k/4,y=Y(v);
    c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
    c.textAlign="right";
    c.fillText(Math.abs(v)>=100?v.toFixed(0):v.toFixed(2),L-6,y+4);}
  for(let k=0;k<=4;k++){const t=tmax*k/4;
    c.textAlign="center";c.fillText(String(Math.round(t)),X(t),H-17);}
  c.textAlign="center";c.fillText("T (K)",L+pw/2,H-3);
  frame(c,L,T,pw,ph,ink,[0,1,2,3,4].map(k=>X(tmax*k/4)),
        [0,1,2,3,4].map(k=>Y(lo+(hi-lo)*k/4)));

  /*  the melting point, where it is on the plot: past it the solid is
      metastable and the curve is describing something the material is not   */
  const tm=series[0].r.Tmelt;
  if(tm&&tm<tmax){ c.save();c.setLineDash([3,3]);c.strokeStyle=ink;
    c.beginPath();c.moveTo(X(tm),T);c.lineTo(X(tm),T+ph);c.stroke();
    c.restore(); c.textAlign="center";c.fillStyle=ink;
    c.fillText("Tₘ",X(tm),T+11); }

  series.forEach(o=>{
    const ours=o.r.kind==="ours";
    /*  Colour says angular or not and the dash says which arm, both from the
        tables every other panel reads.  Written out by hand here it gave the
        re-cut arms the same colour AND the same dash, so two potentials were
        one line, and gave the non-angular re-cut the angular colour.  */
    const col = o.r.nudge_bad ? bad
              : (ours ? (ARM_ANG[o.k] ? p3 : p2) : ref);
    const pts=o.r.pts.filter(q=>isFinite(q[ET_Q]));
    /*  Line style carries the identity, colour only the warning.  The three
        colours here sit within 1.4 of each other in luminance, so on a
        greyscale print or to a colour-blind reader they are one curve drawn
        three times; the dash pattern is what survives that.  It also matches
        the potential panel, where MAU is already solid and UG dashed.       */
    const dash = o.r.kind==="base" ? [1.5,2.5] : (ARM_DASH[o.k]||[]);
    c.strokeStyle=col; c.lineWidth=ours?2:1.3; c.save();
    c.setLineDash(dash);
    c.beginPath(); pts.forEach((q,i)=>{const x=X(q.T),y=Y(q[ET_Q]);
      i?c.lineTo(x,y):c.moveTo(x,y);}); c.stroke(); c.restore();
    pts.forEach(q=>{const x=X(q.T),y=Y(q[ET_Q]);
      c.beginPath(); c.arc(x,y,ours?3:2.2,0,7);
      if(q.above_melt||!q.born_ok){ c.strokeStyle=col; c.lineWidth=1.2;
        c.stroke(); } else { c.fillStyle=col; c.fill(); }});
  });

  /*  AFLOW's density-functional point, where there is one and it belongs to
      this element's structure.  Only B and G exist there, so it appears on
      those two panels and nowhere else, as a hollow ring at T = 0 - a single
      first-principles value, not a curve, and drawn so it cannot be mistaken
      for one.  It answers a different question from the published potentials:
      not "how do we compare with another model" but "how close is the model
      to the underlying physics".                                             */
  /*  Inert since the AFLOW withdrawal: d.aflow is stripped when DATA is
      serialised, for the licence reason recorded in fetch_aflow.py and
      add_elastic_T.py.  The drawing is kept rather than deleted because the
      question it answers - how close is the model to the underlying physics,
      as opposed to another model - is still the right question, and whatever
      replaces AFLOW will want this marker back.                              */
  const af=d.aflow;
  if(af&&af.usable&&(ET_Q==="B"||ET_Q==="G")&&af[ET_Q]!=null){
    const y=Y(af[ET_Q]);
    if(y>=T&&y<=T+ph){
      c.strokeStyle=ink; c.lineWidth=1.4; c.setLineDash([]);
      c.beginPath(); c.arc(X(0),y,5,0,7); c.stroke();
      c.beginPath(); c.moveTo(X(0)-8,y); c.lineTo(X(0)+8,y); c.stroke();
      c.font="10px system-ui"; c.textAlign="left"; c.fillStyle=ink;
      c.fillText("AFLOW (DFT)", X(0)+10, y-7);
    }
  }

  /*  Labels on the curves themselves.  A legend underneath makes the reader
      hold three colours in their head and look away from the plot to spend
      them, and with a published potential beside two of ours there is a real
      question - which line is the baseline - that a swatch answers badly.
      Placed at each curve's last point, pushed apart where they would
      collide, and drawn in the curve's own colour.                          */
  const lab=[];
  series.forEach(o=>{
    const pts=o.r.pts.filter(q=>isFinite(q[ET_Q]));
    if(!pts.length) return;
    const q=pts[pts.length-1];
    const txt=o.r.label;
    /*  THE SAME TABLE THE CURVE USED, not a second copy written by hand.
        This line said `o.k==="tap"?p2:p3`, which gave every arm but `tap` the
        angular colour - so the re-cut arm and the force-matched arm, neither
        of which carries an angular term, were labelled in one colour and
        drawn in another.  The comment above the curve records that exact bug
        being fixed there; it was left standing here, which is what a table
        written out twice does.                                              */
    lab.push({y:Y(q[ET_Q]), x:X(q.T),
              col:(o.r.nudge_bad?bad
                   :(o.r.kind==="ours"?(ARM_ANG[o.k]?p3:p2):ref)),
              txt:txt+(o.r.nudge_bad?" ⚠":"")});
  });
  lab.sort((a,b)=>a.y-b.y);
  for(let i=1;i<lab.length;i++)
    if(lab[i].y-lab[i-1].y<12) lab[i].y=lab[i-1].y+12;
  const over=lab.length?lab[lab.length-1].y-(T+ph):0;
  if(over>0) lab.forEach(l=>l.y-=over);
  c.font="11px system-ui,-apple-system,sans-serif";
  c.textAlign="left";
  /*  FIT THE LABEL TO THE MARGIN.  Vertical collision and vertical overflow
      were both handled; horizontal width was never measured, and 32 of the
      179 labels in this library do not fit the 98 px that R leaves - a
      published potential's is the long one, "Mendelev, unattributed
      (EAM/FS)" at about 164 px.  The canvas clips at its own edge, so the
      reader saw a truncated NAME with nothing to say it had been cut.
      Shortened in the order that loses least: the full name, then without its
      trailing parenthetical - the model class is already carried by the dash
      pattern and by the caption - then a hard cut with an ellipsis, which at
      least admits that something is missing.                                */
  const room=R-8;
  const fit=t=>{
    if(c.measureText(t).width<=room) return t;
    const noparen=t.replace(/\s*\([^()]*\)\s*$/,"");
    if(noparen!==t&&c.measureText(noparen).width<=room) return noparen;
    let u=noparen;
    while(u.length>1&&c.measureText(u+"…").width>room) u=u.slice(0,-1);
    return u+"…";
  };
  lab.forEach(l=>{ c.fillStyle=l.col; c.fillText(fit(l.txt), W-R+6, l.y+4); });

  c.save();c.translate(14,T+ph/2);c.rotate(-Math.PI/2);
  c.textAlign="center";c.fillStyle=ink;
  c.fillText((ET_QS.find(a=>a[0]===ET_Q)||["",ET_Q])[1]+ET_UNIT(ET_Q),0,0);
  c.restore();
}

/*  The pair term and the three-body term, one panel each.

    They were one panel and it stopped working the moment UG joined.  The two
    terms differ by an order of magnitude - tungsten's pair well is 0.93 eV deep
    while its three-body term spans 0.23 - so on shared axes phi3 lay flat on
    zero and was unreadable, and the UG pair well, twice the depth of MAU's,
    pushed the frame until nothing else had room.  Two panels, one x axis, each
    scaled to what it holds.

    Within a panel both potentials share the y axis, because that is the whole
    point: MAU solid, UG dashed, and the difference is the picture.            */
/*  The tetragonal energy path, whose curvature at the origin IS C' =
    (C11-C12)/2.  It is on the page because it is the cheapest honest answer to
    a question the fitted numbers cannot be asked: the curvature at zero strain
    is a target of the fit and comes out right, and it says nothing whatever
    about whether the well it belongs to is a well.

    Reading it: the curve should rise on both sides.  Where it turns over, the
    crystal has a downhill direction under a strain of that size, and across the
    seven cubic metals here the height of that first hump orders the
    finite-temperature failures exactly - iron 0.7 meV fails at five per cent of
    the melting point, tantalum 3.2 meV never fails.  A published potential run
    through the identical scan rises monotonically to twelve per cent.         */
function drawBain(d){
  const recs=armRecs(d,r=>r.bain&&r.bain.E&&r.bain.E.length);
  if(!recs.length) return;
  const s=setup("#bain",0.52); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=58,R=16,T=12,B=30, pw=W-L-R, ph=H-T-B;

  let lo=Infinity, hi=-Infinity, dmax=0;
  recs.forEach(o=>o[1].bain.E.forEach((v,i)=>{
    const y=v*1000; if(!isFinite(y)) return;
    lo=Math.min(lo,y); hi=Math.max(hi,y);
    dmax=Math.max(dmax,Math.abs(o[1].bain.d[i]));}));
  /*  clipped, because the far ends of the path run to hundreds of meV and
      would flatten the only part anyone needs to see                        */
  hi=Math.min(hi,80); lo=Math.max(lo,-140);
  if(!(hi>lo)){hi=lo+1;}
  const X=x=>L+(x+dmax)/(2*dmax)*pw, Y=v=>T+ph-(v-lo)/(hi-lo)*ph;

  c.strokeStyle=lin;c.lineWidth=1;
  c.font="11px ui-monospace,Consolas,monospace";c.fillStyle=ink;
  for(let k=0;k<=4;k++){const v=lo+(hi-lo)*k/4,y=Y(v);
    c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
    c.textAlign="right";c.fillText(v.toFixed(0),L-6,y+4);}
  for(let k=0;k<=4;k++){const x=-dmax+2*dmax*k/4;
    c.textAlign="center";c.fillText(x.toFixed(2),X(x),H-17);}
  c.textAlign="center";c.fillText("tetragonal strain δ",L+pw/2,H-3);
  c.save();c.translate(13,T+ph/2);c.rotate(-Math.PI/2);c.textAlign="center";
  c.fillText("E − E₀  (meV/atom)",0,0);c.restore();
  frame(c,L,T,pw,ph,ink,[0,1,2,3,4].map(k=>X(-dmax+2*dmax*k/4)),
        [0,1,2,3,4].map(k=>Y(lo+(hi-lo)*k/4)));

  //  the zero line is the whole reading: below it, bcc is not the minimum
  c.save();c.setLineDash([2,3]);c.strokeStyle=ink;
  c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();c.restore();

  recs.forEach(o=>{
    const b=o[1].bain, col=ARM_ANG[o[0]]?p3:p2;
    c.strokeStyle=col;c.lineWidth=2;c.save();
    c.setLineDash(ARM_DASH[o[0]]||[]);
    c.beginPath();
    b.d.forEach((x,i)=>{const y=Y(Math.max(lo,Math.min(hi,b.E[i]*1000)));
      i?c.lineTo(X(x),y):c.moveTo(X(x),y);});
    c.stroke();c.restore();
    //  mark the ledge, if there is one
    if(b.turn_up!==null&&b.turn_up!==undefined){
      const x=X(b.turn_up), y=Y(Math.max(lo,Math.min(hi,b.barrier_up*1000)));
      c.fillStyle=col;c.beginPath();c.arc(x,y,3.5,0,7);c.fill();
      c.beginPath();c.arc(x,y,6.5,0,7);c.strokeStyle=col;c.lineWidth=1;c.stroke();
    }
  });
}

/*  The generalised stacking fault curve.  The whole reading is where the
    curve sits at one third of the period against the zero line: above it a
    fault costs energy and heals, below it the faulted crystal is the cheaper
    one and the fault is permanent.  So the zero line is drawn heavier than
    the grid, and the two points that have names are marked.                */
function drawGamma(d){
  /*  frac as well as gamma: the curve is drawn from g.frac.forEach with
      g.gamma[i] read inside it, so a record carrying one and not the other
      takes the whole panel down with "reading 'forEach' of undefined".  Every
      record in the library has both, which is exactly why the guard was
      written against gamma alone and the mismatch went unnoticed: the first
      record to arrive half-written from a separate run would have found it,
      and a guard that only holds because the data happens to be complete is
      not a guard.  */
  const ok=r=>r&&r.gamma&&r.gamma.length&&r.frac&&r.frac.length;
  const recs=armRecs(d,r=>ok(r.stacking));
  const bases=Object.entries(d.baseline_stacking||{}).filter(o=>ok(o[1]));
  if(!recs.length&&!bases.length) return;
  const s=setup("#gamma",0.52); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=58,R=16,T=12,B=30, pw=W-L-R, ph=H-T-B;

  let lo=0, hi=0;
  const all=recs.map(o=>o[1].stacking).concat(bases.map(o=>o[1]));
  all.forEach(g=>g.gamma.forEach(v=>{if(isFinite(v)){
    lo=Math.min(lo,v); hi=Math.max(hi,v);}}));
  if(!(hi>lo)) hi=lo+1;
  const pad=(hi-lo)*0.06; lo-=pad; hi+=pad;
  const X=x=>L+x*pw, Y=v=>T+ph-(v-lo)/(hi-lo)*ph;

  c.strokeStyle=lin;c.lineWidth=1;
  c.font="11px ui-monospace,Consolas,monospace";c.fillStyle=ink;
  for(let k=0;k<=4;k++){const v=lo+(hi-lo)*k/4,y=Y(v);
    c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
    c.textAlign="right";c.fillText(v.toFixed(0),L-6,y+4);}
  [[0,"0"],[1/6,"1/6"],[1/3,"1/3"],[2/3,"2/3"],[1,"1"]].forEach(t=>{
    c.textAlign="center";c.fillText(t[1],X(t[0]),H-17);});
  c.textAlign="center";
  c.fillText("shift along [11̄2̄], in periods",L+pw/2,H-3);
  c.save();c.translate(13,T+ph/2);c.rotate(-Math.PI/2);c.textAlign="center";
  c.fillText("γ  (mJ/m²)",0,0);c.restore();
  frame(c,L,T,pw,ph,ink,[0,1/6,1/3,2/3,1].map(X),
        [0,1,2,3,4].map(k=>Y(lo+(hi-lo)*k/4)));

  //  zero, and the partial.  A curve that ends below zero at 1/3 is the
  //  finding; a curve that does not close at 1 is a broken cell.
  c.save();c.strokeStyle=ink;c.globalAlpha=0.55;c.lineWidth=1.5;
  c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();
  c.setLineDash([2,3]);c.globalAlpha=0.4;
  c.beginPath();c.moveTo(X(1/3),T);c.lineTo(X(1/3),T+ph);c.stroke();
  c.restore();

  bases.forEach(o=>{
    const g=o[1];
    c.strokeStyle=lin;c.lineWidth=1;c.beginPath();
    g.frac.forEach((x,i)=>{const y=Y(g.gamma[i]);
      i?c.lineTo(X(x),y):c.moveTo(X(x),y);});
    c.stroke();});

  recs.forEach(o=>{
    const g=o[1].stacking, col=ARM_ANG[o[0]]?p3:p2;
    c.strokeStyle=col;c.lineWidth=2;c.save();
    c.setLineDash(ARM_DASH[o[0]]||[]);
    c.beginPath();
    g.frac.forEach((x,i)=>{const y=Y(g.gamma[i]);
      i?c.lineTo(X(x),y):c.moveTo(X(x),y);});
    c.stroke();c.restore();
    [[1/6,g.usf],[1/3,g.isf]].forEach(m=>{
      if(m[1]===null||m[1]===undefined) return;
      c.fillStyle=col;c.beginPath();c.arc(X(m[0]),Y(m[1]),3.5,0,7);c.fill();});
  });
}

/*  Thermal expansion: the cell edge against temperature, from NPT dynamics.
    The published potentials are drawn in thin grey through the identical
    barostat, because "our slope is 0.7 of experiment" means nothing until it
    is known what a working potential's slope does in the same code.  Each
    curve is normalised to its own value at the lowest temperature, since what
    is being compared is a slope and not a lattice constant.               */
function drawExpansion(d){
  const recs=armRecs(d,r=>r.expansion&&r.expansion.T&&r.expansion.T.length>1);
  const bases=Object.entries(d.baseline_expansion||{})
    .filter(o=>o[1]&&o[1].T&&o[1].T.length>1);
  if(!recs.length&&!bases.length) return;
  const s=setup("#expan",0.52); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=62,R=16,T=12,B=30, pw=W-L-R, ph=H-T-B;

  const series=recs.map(o=>({g:o[1].expansion,col:ARM_ANG[o[0]]?p3:p2,
                             dash:ARM_DASH[o[0]]||[],wide:2}))
    .concat(bases.map(o=>({g:o[1],col:lin,dash:[],wide:1})));
  let t0=Infinity,t1=-Infinity,lo=Infinity,hi=-Infinity;
  const norm=g=>g.a.map(x=>100*(x/g.a[0]-1));
  series.forEach(o=>{o.y=norm(o.g);
    o.g.T.forEach(t=>{t0=Math.min(t0,t);t1=Math.max(t1,t);});
    o.y.forEach(v=>{lo=Math.min(lo,v);hi=Math.max(hi,v);});});
  /*  The experimental slope, drawn from the same starting point.  Taken from
      `series`, whose entries carry a .g - NOT from `recs`, whose entries are
      still the ["tap", record] pairs they were built from.  Reading recs[0].g
      is undefined, and asking undefined for a property throws, which killed
      this function before it drew a single line and left an empty canvas with
      nothing else on the page disturbed.  A blank plot is the one failure
      mode a build cannot see.                                              */
  const exp=series.length?series[0].g.alpha_exp_1e6:null;
  if(exp) hi=Math.max(hi,100*exp*1e-6*(t1-t0));
  if(!(hi>lo)){hi=lo+0.1;}
  const pad=(hi-lo)*0.08; lo-=pad; hi+=pad;
  const X=t=>L+(t-t0)/(t1-t0||1)*pw, Y=v=>T+ph-(v-lo)/(hi-lo)*ph;

  c.strokeStyle=lin;c.lineWidth=1;
  c.font="11px ui-monospace,Consolas,monospace";c.fillStyle=ink;
  for(let k=0;k<=4;k++){const v=lo+(hi-lo)*k/4,y=Y(v);
    c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
    c.textAlign="right";c.fillText(v.toFixed(2),L-6,y+4);}
  for(let k=0;k<=4;k++){const t=t0+(t1-t0)*k/4;
    c.textAlign="center";c.fillText(t.toFixed(0),X(t),H-17);}
  c.textAlign="center";c.fillText("T (K)",L+pw/2,H-3);
  c.save();c.translate(13,T+ph/2);c.rotate(-Math.PI/2);c.textAlign="center";
  c.fillText("Δa/a  (%)",0,0);c.restore();
  frame(c,L,T,pw,ph,ink,[0,1,2,3,4].map(k=>X(t0+(t1-t0)*k/4)),
        [0,1,2,3,4].map(k=>Y(lo+(hi-lo)*k/4)));

  //  zero: below it the crystal is CONTRACTING on heating
  c.save();c.strokeStyle=ink;c.globalAlpha=0.5;c.lineWidth=1.5;
  c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();c.restore();

  if(exp){
    c.save();c.setLineDash([4,4]);c.strokeStyle=ink;c.globalAlpha=0.75;
    c.lineWidth=1.5;c.beginPath();
    c.moveTo(X(t0),Y(0));c.lineTo(X(t1),Y(100*exp*1e-6*(t1-t0)));
    c.stroke();c.restore();
  }
  series.forEach(o=>{
    c.strokeStyle=o.col;c.lineWidth=o.wide;c.save();c.setLineDash(o.dash);
    c.beginPath();
    o.g.T.forEach((t,i)=>{const y=Y(o.y[i]);
      i?c.lineTo(X(t),y):c.moveTo(X(t),y);});
    c.stroke();c.restore();});
}

/*  Every plot goes through here.  A draw function that throws leaves its
    canvas blank and disturbs nothing else on the page, which is the one
    failure this project cannot see from a build: the HTML is written, the
    checker passes, and a section is quietly empty.  It happened once, to the
    thermal expansion plot, from reading a property off the wrong array.  Now
    the canvas says so.                                                     */
function plots(d){
  const jobs=[["dispersion",drawDisp],["finite-T dispersion",drawFT],
              ["thermo",drawThermo],
              ["polar",drawPolar],["elasticT",drawElasticT],
              ["bain",drawBain],["gamma",drawGamma],
              ["expansion",drawExpansion],["main",draw]];
  jobs.forEach(([name,fn])=>{
    try{ fn(d); }
    catch(e){
      console.error("plot "+name+" failed:",e);
      const box=$("#plotfail");
      if(box){ box.style.display="block";
        box.textContent="A plot failed to draw ("+name+"): "+e.message
          +" — the numbers above are unaffected."; }
    }});
}

function draw(d){
  const s=setup("#plot",0.92); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const u=ugData(d);          //  hRange now travels with each extra series
  const L=52,R=12,GAP=26,TOP=10,BOT=30;
  const ph=(H-TOP-BOT-GAP)/2, pw=W-L-R;
  /*  Where the window starts is data, not a constant.  It was fixed at
      0.78 d_nn, which suits an element whose phi2 minimum sits near d_nn but
      cuts the repulsive wall off for the few whose well lies further in:
      vanadium's minimum is at 0.818 d_nn and its wall at 0.494, so the panel
      opened past the wall and showed a bare rising curve that looked nothing
      like the others.  The left edge now backs off until phi2 has climbed one
      well-depth above zero, which puts the wall at the same visual height for
      every element and is what makes them comparable at a glance.  It is
      clamped at 0.45 d_nn, and 32 of 38 elements are unaffected.

      A fit with no wall to find keeps the old edge; see coreless() below.    */
  function wallOf(dd){
    const lo0=0.45*d.dnn, hi0=2.05*d.dnn, M=400;
    let vmin=Infinity, imin=0, v=[];
    for(let i=0;i<M;i++){const r=lo0+(hi0-lo0)*i/(M-1);
      const y=phi2(r,dd); v.push([r,y]); if(y<vmin){vmin=y;imin=i;}}
    for(let i=imin;i>=0;i--) if(v[i][1]>=Math.abs(vmin)) return v[i][0];
    return null;                      // monotone inward: no wall exists
  }
  const walls=[wallOf(d)].concat(u?[wallOf(u)]:[])
      .concat(["rc","rc_ug"].filter(k=>d[k]&&typeof d[k].D==="number")
                            .map(k=>wallOf(d[k])))
      .filter(w=>w!==null);
  const r0=walls.length?Math.max(0.45*d.dnn,
              Math.min(0.78*d.dnn, Math.min.apply(null,walls)-0.04*d.dnn))
            :0.78*d.dnn,
        r1=2.05*d.dnn, N=380;
  /*  Every record drawn beside the reference.  A record carrying lam2 or
      lam4 is drawn as a BAND, not a line: with the angular factor phi3 is an
      interval at each leg-sum, one value per apex angle, and a single curve
      would be a picture of a potential that does not exist.  That is why UG
      has always been a band here.  */
  const extras=[];
  if(u) extras.push({rec:u, label:"UG", dash:[5,3]});
  /*  The re-cut candidates belong here more than anywhere else on the page:
      what makes them a different potential is where phi2 is cut - rcut3 is
      unchanged for every element - and the parameters that go with it are not
      small adjustments: copper's m goes 13.78 to 1.87 and its alpha 0.285 to
      0.833.  Everywhere else on the page the candidate is read through a
      consequence; here it is the curve itself.  */
  [["rc","re-cut",[2,2]],["rc_ug","re-cut UG",[6,2,1,2]]].forEach(a=>{
    const r=d[a[0]];
    if(r && typeof r.D==="number")
      extras.push({rec:r, label:a[1]+(armBad(r)?" — REJECTED":""),
                   dash:a[2]});
  });
  extras.forEach(e=>{
    e.ang=!!(e.rec.lam2||e.rec.lam4);
    e.hr=e.ang?hRange(e.rec):null;
    e.p2=[]; e.lo=[]; e.hi=[];
  });
  const xs=[],y2=[],y3=[];
  for(let i=0;i<N;i++){const r=r0+(r1-r0)*i/(N-1);
    xs.push(r); y2.push(phi2(r,d)); y3.push(phi3(2*r,d));
    extras.forEach(e=>{
      const b=phi3(2*r,e.rec);
      e.p2.push(phi2(r,e.rec));
      if(e.ang){e.lo.push(Math.min(b*e.hr[0],b*e.hr[1]));
                e.hi.push(Math.max(b*e.hr[0],b*e.hr[1]));}
      else {e.lo.push(b); e.hi.push(b);}
    });}
  const X=r=>L+(r-r0)/(r1-r0)*pw;
  c.font="11px ui-monospace,Consolas,monospace";

  /*  a panel: its own vertical range, clipped so the repulsive wall - which
      rises without bound - cannot decide the framing                          */
  /*  The vertical range comes from the data, not from the depth of a well.
      Tying the ceiling to |lo| assumes there IS a well: it works for phi2 and
      fails for phi3, which is one-signed for most elements - tungsten's is
      positive everywhere, so a range built around a well spent half the panel
      below zero and clipped the band at a tenth of its height.  Instead the
      ceiling is the largest value from 0.9 d_nn outwards, which leaves the
      repulsive wall to run off the top where it belongs, and the floor is the
      true minimum.  */
  const iRef=xs.findIndex(r=>r>=0.9*d.dnn);
  /*  Each panel carries its own key.  A single one at the top described line
      styles, and in the lower panel UG is not a line at all but a band, so the
      key there explained nothing about the largest thing on the plot. */
  const key=(T,items)=>{
    c.font="11px ui-monospace,Consolas,monospace";c.textAlign="left";
    items.forEach(([t,kind,col],i)=>{
      const yy=T+13+i*14, x=W-R-92;
      if(kind==="band"){
        c.save();c.globalAlpha=0.20;c.fillStyle=col;
        c.fillRect(x,yy-5,20,10);c.restore();
      }else{
        c.strokeStyle=col;c.lineWidth=2;c.setLineDash(kind==="dash"?[5,3]:[]);
        c.beginPath();c.moveTo(x,yy);c.lineTo(x+20,yy);c.stroke();c.setLineDash([]);
      }
      c.fillStyle=ink;c.fillText(t,x+25,yy+4);});
  };
  const panel=(T,sets,label,cuts)=>{
    let lo=0,hi=0;
    sets.forEach(([ys])=>ys.forEach((v,i)=>{
      if(v<lo)lo=v;
      if(i>=iRef&&v>hi)hi=v;}));
    /*  where there is a real well, keep the wall from taking more than half
        the panel; where the term is one-signed - phi3 usually is - there is no
        well to scale against and the cap must not fire */
    if(lo < -0.05*Math.max(hi,1e-9)) hi=Math.min(hi,Math.abs(lo)*1.5);
    const pad=0.08*((hi-lo)||1);
    lo-=pad; hi+=pad;
    if(hi-lo<1e-6){hi=lo+1e-6;}
    const Y=v=>T+ph-(v-lo)/(hi-lo)*ph;
    c.strokeStyle=lin;c.lineWidth=1;c.fillStyle=ink;
    for(let k=0;k<=3;k++){const v=lo+(hi-lo)*k/3,y=Y(v);
      c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
      c.textAlign="right";c.fillText(v.toFixed(2),L-6,y+4);}
    c.strokeStyle=ink;c.save();c.setLineDash([3,3]);
    c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();
    c.beginPath();c.moveTo(X(d.dnn),T);c.lineTo(X(d.dnn),T+ph);c.stroke();
    c.restore();
    sets.forEach(([ys,col,w,dash,band])=>{
      if(band){
        c.save();c.globalAlpha=0.20;c.fillStyle=col;c.beginPath();
        ys.forEach((v,i)=>{const y=Y(Math.min(Math.max(v,lo),hi));
          i?c.lineTo(X(xs[i]),y):c.moveTo(X(xs[i]),y);});
        band.forEach((v,i)=>{const j=N-1-i;
          c.lineTo(X(xs[j]),Y(Math.min(Math.max(band[j],lo),hi)));});
        c.closePath();c.fill();c.restore();
        return;
      }
      c.strokeStyle=col;c.lineWidth=w;c.setLineDash(dash||[]);c.beginPath();
      let on=false;
      ys.forEach((v,i)=>{if(v>hi||v<lo){on=false;return;}
        const x=X(xs[i]),y=Y(v); on?c.lineTo(x,y):(c.moveTo(x,y),on=true);});
      c.stroke();c.setLineDash([]);});
    /*  A cut that falls inside the window is marked on the axis.  Most do
        not - copper's rcut2 is 8.3 A against a window that ends at 5.2 - but
        where the re-cut moved it a long way inwards it is the whole story:
        caesium goes from 15.72 A to 7.30, and 7.30 is on screen.  */
    (cuts||[]).forEach(cu=>{
      if(cu.r<=r0||cu.r>=r1) return;
      const x=X(cu.r);
      c.save();c.strokeStyle=cu.col;c.globalAlpha=0.8;c.lineWidth=1.2;
      c.setLineDash(cu.dash&&cu.dash.length?cu.dash:[3,3]);
      c.beginPath();c.moveTo(x,T);c.lineTo(x,T+ph);c.stroke();c.restore();
      c.save();c.fillStyle=ink;c.globalAlpha=0.75;c.textAlign="center";
      c.font="10px ui-monospace,Consolas,monospace";
      c.fillText(cu.tag,x,T+ph-4);c.restore();
    });
    frame(c,L,T,pw,ph,ink,[0,1,2,3,4].map(k=>X(r0+(r1-r0)*k/4)),
          [0,1,2,3].map(k=>Y(lo+(hi-lo)*k/3)));
    c.fillStyle=ink;c.textAlign="left";c.fillText(label,L+7,T+13);
  };

  const cuts2=[{r:d.rcut2,col:p2,dash:[],tag:"cut"}].concat(
      extras.map(e=>({r:e.rec.rcut2,col:p2,dash:e.dash,tag:"cut"})))
    .filter(cu=>typeof cu.r==="number");
  panel(TOP, [[y2,p2,2.2]].concat(extras.map(e=>[e.p2,p2,1.6,e.dash])),
        "φ₂(r)", cuts2);
  if(extras.length)
    key(TOP, [["MAU","line",p2]].concat(extras.map(e=>[e.label,"dash",p2])));
  //  bands first so the lines stay readable on top of them
  panel(TOP+ph+GAP,
        extras.filter(e=>e.ang).map(e=>[e.hi,p3,0,null,e.lo])
          .concat([[y3,p3,1.8]])
          .concat(extras.filter(e=>!e.ang).map(e=>[e.lo,p3,1.6,e.dash])),
        "φ₃(r,r)");
  if(extras.length)
    key(TOP+ph+GAP, [["MAU","line",p3]].concat(extras.map(e=>
      [e.ang?e.label+", all θ":e.label, e.ang?"band":"dash", p3])));

  /*  one axis label for both panels: the quantity is the same in each and the
      unit was previously carried inside the panel captions, where it read as
      part of the function name */
  c.save();c.translate(14,TOP+ph+GAP/2);c.rotate(-Math.PI/2);
  c.textAlign="center";c.fillStyle=ink;c.fillText("Energy (eV)",0,0);c.restore();
  c.fillStyle=ink;c.textAlign="center";
  /*  d_nn shares the tick row, so a tick that lands under it is dropped
      rather than drawn through it.  The marker is the more useful of the two
      - it says where the nearest neighbour sits - and one missing tick out of
      five costs nothing, while two labels on the same pixels cost the row. */
  const xdnn=X(d.dnn);
  for(let k=0;k<=4;k++){const r=r0+(r1-r0)*k/4;
    if(Math.abs(X(r)-xdnn)<24) continue;
    c.fillText(r.toFixed(2),X(r),H-18);}
  c.textAlign="center";c.fillText("d_nn",xdnn,H-18);
  c.fillText("r (A)",L+(W-L-12)/2,H-4);
}

let DISP_MODE = null;

/*  Draw the UG result on top of the MAU one, everywhere it exists.

    There is no global switch any more.  It lived in the UG section far below
    the plots it governed, and it did nothing at all while the dispersion was in
    one of its DFT-comparison modes - which is the default for most elements, so
    the usual experience was a phonon panel with no UG in it and no visible
    reason why.  The directional response and the thermodynamics now show both
    unconditionally, and the dispersion offers UG as one of its own view modes,
    beside "vs Materials Project" where the choice belongs.

    UG used to appear as a five-row table of elastic constants beside a MAU
    entry that carried mechanics, thermodynamics, a dispersion curve and a
    measured-phonon comparison, which made the two look like different kinds of
    result rather than the same calculation with one extra term.  They are
    produced by the same code now, so the page draws them on one set of axes:
    the difference between two curves is the angular contribution, and reading
    it off a plot is the whole point.

    Off by default - the single-potential view is the common case - and only
    offered where d.ug exists and was fitted at the same three-body cutoff. */
const ugData = d => (d.ug && d.ug.comparable) ? d.ug : null;
/*  The selector's "two potentials together" modes.  UG needs `comparable`
    because it is a different FORM fitted to the same targets and the flag says
    whether the two can be put on one axis; the re-cut candidates are the same
    form at a different cutoff, so having a dispersion on the standard path is
    the whole condition.  */
const SECOND = {
  ug:       d => ugData(d),
  recut:    d => (d.rc    && d.rc.ld    && d.rc.ld.std)    ? d.rc    : null,
  recut_ug: d => (d.rc_ug && d.rc_ug.ld && d.rc_ug.ld.std) ? d.rc_ug : null,
};
const SECOND_LABEL = {ug: "UG", recut: "re-cut", recut_ug: "re-cut UG"};
const secondOf = (d, m) => (SECOND[m] ? SECOND[m](d) : null);

/*  Does phi2 turn over on the way in, or fall for ever?  A pair term with no
    minimum has no repulsive core: two atoms lower their energy by merging.
    One fit of the 76 in this library does this - niobium's MAU - and the
    symptom is visible on the plot, so the plot says what it is rather than
    leaving the reader to wonder why that panel looks wrong.                 */
function coreless(d){
  for(let i=0;i<40;i++){
    const r=(0.28+0.012*i)*d.dnn;
    if(phi2(r,d) > 0) return false;
  }
  return true;
}

/*  Which parameter set the page is showing.

    One state for the chip row at the top and the export block at the bottom,
    because they must never disagree: every parameter differs between the two
    fits - iron's D is 0.345 eV for MAU and 0.284 for UG - so showing lambda
    beside MAU's D would describe a potential that was never fitted.  The row
    switches as a whole or not at all. */

/* frequencies are stored in cm-1 and plotted in THz */
const CM1_PER_THZ=33.35641;
const THZ=v=>v/CM1_PER_THZ;

/* a round tick spacing giving roughly `want` intervals over `span` */
function niceStep(span,want){
  const raw=span/Math.max(want,1);
  const p=Math.pow(10,Math.floor(Math.log10(raw||1)));
  for(const m of [1,2,2.5,5,10]) if(raw<=m*p) return m*p;
  return 10*p;
}

function drawDisp(d){
  /*  The "this potential only" view follows the parameter set on screen, so
      the dispersion belongs to the parameters printed above it.  The DFT
      modes cannot: their `ours` array was resampled at the reference's own
      q-points for the ROOT record only, and no such resampling exists for
      the other arms - so those views stay MAU whatever is selected, and the
      note under the panel says so rather than letting a reader assume.  */
  const Pd=parSet(d);
  const g=(Pd.ld&&Pd.ld.std)?Pd.ld:(d.ld||{});
  /* Two independent DFT references, and they cover different metals: MP has
     13 of ours, MC3D adds Ag Au Cr Mo Nb Ni Pb Pd Rh Ta - including the bcc
     transition metals whose anisotropy the form cannot reach, which had no
     comparison at all before.  Both are drawn the same way, our curve
     evaluated at THEIR q-points, so the residual is physics not interpolation. */
  /*  JARVIS is a third DFT source and mostly agrees to within ten per cent,
      but a few of its runs disagree badly - iron by a factor of 2.6, which is
      what a non-spin-polarised ferromagnet looks like.  add_jarvis_overlay.py
      compares each one against whatever reference already exists and clears
      only those within 25 %; the rest stay in the data, flagged, and are not
      offered here.  */
  const REF={mp:(d.mp||{}).phonon, mc3d:d.mc3d,
             jarvis:(d.jarvis&&d.jarvis.trusted)?d.jarvis:null};
  const has=k=>REF[k]&&REF[k].ours;
  /*  Al Ba Co Cr Na have JARVIS and nothing else, so it has to be reachable as
      a default too - otherwise their only comparison sits in the data unseen. */
  if(SECOND[DISP_MODE] && !secondOf(d,DISP_MODE)) DISP_MODE=null;  // no such arm
  if(DISP_MODE===null||!(DISP_MODE==="std"||SECOND[DISP_MODE]||has(DISP_MODE)))
    /*  Where a MEASURED curve exists it is the better default, and it is only
        drawn on the standard path - the DFT modes replace that path with the
        reference's own.  Opening on a DFT comparison for an element that has
        neutron data would hide the stronger test behind a selector.  The same
        applies to a MODEL curve: cobalt has no measured points (vanadium had
        none until its x-ray dispersion was entered on 2026-09-15) and opened
        on JARVIS, with its reconstructed dispersion invisible until someone
        changed the selector.  Either kind of reference counts.  */
    DISP_MODE = ((d.exp_curve && d.exp_curve.segs)
              || (d.model_curve && d.model_curve.segs)) ? "std"
              : (has("mp") ? "mp"
              : (has("mc3d") ? "mc3d" : (has("jarvis") ? "jarvis" : "std")));
  /*  The DFT references are stored at THEIR q-points and UG at ours, so the
      two cannot share an x axis; "ug" is therefore a view of the standard path
      with the second potential drawn on it, not a reference mode. */
  const withUG = secondOf(d, DISP_MODE);
  const MODE = withUG ? "std" : DISP_MODE;
  const mpp = REF[MODE];
  const s=setup("#disp",0.66); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=58,R=16,T=20,B=34, pw=W-L-R, ph=H-T-B;
  const G=String.fromCharCode(915);
  const lab=t=>t==="G"?G:t;
  const GAP=10;          /* width of a discontinuity in the path, both views */

  /* Solid frame on all four sides.  Ticks on the left and the right, numbers
     only on the left.  Returns the y mapping - imaginary modes are stored as
     negative frequencies, and the range opens downwards to keep them inside
     the box instead of drawing them off the bottom of the canvas. */
  const axes=(loRaw,hiRaw)=>{
    /*  A mode at -0.1 cm^-1 is arithmetic, not an instability, and letting it
        open the axis downwards costs a fifth of the panel and puts a negative
        tick under a dispersion that has nothing negative in it.  Anything
        shallower than 1 cm^-1 is floored to zero; anything deeper is drawn,
        because a real imaginary branch must be visible.  */
    if(loRaw > -1.0/CM1_PER_THZ) loRaw = 0;
    const step=niceStep(hiRaw-Math.min(0,loRaw),5);
    const hi=Math.ceil(hiRaw/step)*step || step;
    const lo=loRaw<0?-Math.ceil(-loRaw/step)*step:0;
    const dec=(step<1||Math.abs(step-Math.round(step))>1e-9)?1:0;
    const Y=v=>T+ph-(v-lo)/(hi-lo)*ph;
    c.font="11px ui-monospace,Consolas,monospace";
    for(let k=Math.round(lo/step);k<=Math.round(hi/step);k++){
      const v=k*step, y=Y(v);
      c.strokeStyle=ink;c.lineWidth=1.2;
      c.beginPath();c.moveTo(L,y);c.lineTo(L+5,y);c.stroke();
      c.beginPath();c.moveTo(W-R,y);c.lineTo(W-R-5,y);c.stroke();
      c.fillStyle=ink;c.textAlign="right";
      c.fillText(v.toFixed(dec),L-8,y+4);
    }
    if(lo<0){            /* omega = 0, solid like the rest of the frame */
      c.strokeStyle=ink;c.lineWidth=1;
      c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();
    }
    c.strokeStyle=ink;c.lineWidth=1.2;c.strokeRect(L,T,pw,ph);
    c.fillStyle=ink;c.textAlign="center";
    c.save();c.translate(15,T+ph/2);c.rotate(-Math.PI/2);
    c.fillText("Frequency (THz)",0,0);c.restore();
    return Y;
  };

  /* measured frequencies, drawn where the plot passes that point.  Collected
     while the verticals are drawn, because only then is the x of each label
     known, and both path layouts share this. */
  const expPts=(d.exp_phonon||{}).points||[];
  const seen={};
  const marker=(x,name,Y)=>{
    if(seen[name])return; seen[name]=1;
    const p=expPts.find(p=>p.name===name); if(!p)return;
    c.strokeStyle=ink;c.lineWidth=1.4;
    p.exp.forEach(f=>{c.beginPath();c.arc(x,Y(f),3.4,0,2*Math.PI);c.stroke();});
  };

  /* dashed vertical through a high-symmetry point, inward ticks top and bottom */
  const vline=(x,name)=>{
    c.save();c.strokeStyle=ink;c.lineWidth=1;c.globalAlpha=0.55;
    c.setLineDash([5,4]);
    c.beginPath();c.moveTo(x,T);c.lineTo(x,T+ph);c.stroke();c.restore();
    c.strokeStyle=ink;c.lineWidth=1.2;
    c.beginPath();c.moveTo(x,T+ph);c.lineTo(x,T+ph-6);c.stroke();
    c.beginPath();c.moveTo(x,T);c.lineTo(x,T+6);c.stroke();
    if(name!==undefined){
      c.textAlign="center";c.fillStyle=ink;c.fillText(lab(name),x,T+ph+17);}
  };

  /* ---- shared q-points with Materials Project ---- */
  if(MODE!=="std"){
    const ours=mpp.ours, mp=mpp.f, n=ours[0].length;
    let hi=0, lo=0;
    [ours,mp].forEach(set=>set.forEach(b=>b.forEach(v=>{
      if(v>hi)hi=v; if(v<lo)lo=v;})));
    const Y=axes(THZ(lo),THZ(hi));
    /* MP's q-list jumps at the genuine discontinuities of the path - hcp A|L
       and M|K, bcc H|P, fcc K|U.  Leave the same gap the single-potential view
       leaves, instead of joining across the jump with a straight line. */
    const brk=new Set(mpp.breaks||[]);
    const unit=(pw-GAP*brk.size)/Math.max(n-1-brk.size,1);
    const xs=new Array(n); xs[0]=L;
    for(let i=1;i<n;i++) xs[i]=xs[i-1]+(brk.has(i-1)?GAP:unit);
    /*  A label's index is the nearest SAMPLE, and the reference path does not
        have to sample the point it brackets - silver's interior Gamma is
        never visited, the list steps from (0.0117,0.0117,0.0235) straight to
        (0.0096,0.0096,0.0096).  Drawing the line on the nearest index put it
        a little to the left of Gamma, where the acoustic branches have not
        reached zero, and the plot then seemed to show frequencies that do
        not vanish there.  marks_x carries the fractional index computed in
        fix_mark_x.py; xat interpolates the x between the two samples.  */
    const MX = mpp.marks_x || mpp.marks || [];
    const xat = xf => {const j=Math.floor(xf), t=xf-j;
      return (j+1<xs.length) ? xs[j]+t*(xs[j+1]-xs[j])
                             : xs[Math.min(j,xs.length-1)];};
    MX.forEach(([xf,name])=>vline(xat(xf),name));
    const curve=(set,col,w,dash)=>{
      c.strokeStyle=col;c.lineWidth=w;c.setLineDash(dash);
      set.forEach(b=>{c.beginPath();
        b.forEach((v,i)=>{const x=xs[i],y=Y(THZ(v));
          if(i===0||brk.has(i-1))c.moveTo(x,y);else c.lineTo(x,y);});
        c.stroke();});
      c.setLineDash([]);
    };
    curve(mp,p3,1.4,[4,3]);
    curve(ours,p2,1.8,[]);
    MX.forEach(([xf,name])=>marker(xat(xf),name,Y));
    return;
  }

  /* ---- our curve on the same path, no MP reference available ---- */
  const segs = g.std||[];
  if(!segs.length) return;
  const u=withUG||null, usegs=(u&&u.ld)?u.ld.std:null;
  let hi=0, lo=0;
  const span=set=>set.forEach(p=>p.branches.forEach(b=>b.forEach(v=>{
    if(v===null)return; if(v>hi)hi=v; if(v<lo)lo=v;})));
  span(segs); if(usegs) span(usegs);
  /*  The measured points and a fitted model are drawn on this axis as well,
      so they set its range too.  Scaled on the computed branches alone,
      niobium's highest measured frequencies (6.49 THz, 6.59 with the error
      bar) sat above an axis that stopped at 6, outside the frame.  */
  const onPath=[];
  if(d.exp_curve&&d.exp_curve.segs) Object.values(d.exp_curve.segs).forEach(
    v=>v.forEach(pt=>onPath.push(pt[1]+(pt[2]||0))));
  if(d.model_curve&&d.model_curve.segs) Object.values(d.model_curve.segs).forEach(
    v=>v.forEach(row=>row.slice(1).forEach(x=>{if(x!=null)onPath.push(x);})));
  expPts.forEach(p=>(p.exp||[]).forEach(f=>onPath.push(f)));
  onPath.forEach(f=>{const v=f*CM1_PER_THZ; if(v>hi)hi=v;});
  const Y=axes(THZ(lo),THZ(hi));
  const brk=segs.map((p,i)=>i<segs.length-1 && p.b!==segs[i+1].a);
  const tot=segs.reduce((a,p)=>a+p.len,0);
  const avail=pw-GAP*brk.filter(Boolean).length;
  /*  UG is sampled on the same path with the same number of points, so the two
      share one x mapping; drawing it on its own segment lengths would put the
      curves out of register wherever a segment length differed.  */
  const panel=(p,X,col,w,dash)=>{
    c.strokeStyle=col;c.lineWidth=w;c.setLineDash(dash);
    p.branches.forEach(b=>{c.beginPath();let on=false;
      b.forEach((v,k)=>{if(v===null){on=false;return;}
        const x=X(k),y=Y(THZ(v));
        if(!on){c.moveTo(x,y);on=true;}else c.lineTo(x,y);});
      c.stroke();});
    c.setLineDash([]);
  };
  let x0=L;
  segs.forEach((p,i)=>{
    const w=avail*p.len/tot, n=p.n;
    const X=k=>x0+(n>1?k/(n-1):0)*w;
    vline(x0,lab(p.a));  marker(x0,p.a,Y);
    if(i===segs.length-1||brk[i]){vline(x0+w,lab(p.b)); marker(x0+w,p.b,Y);}
    if(usegs&&usegs[i]) panel(usegs[i],k=>x0+(usegs[i].n>1?k/(usegs[i].n-1):0)*w,
                              p3,1.5,[5,3]);
    panel(p,X,p2,1.8,[]);
    /*  Measured points along this segment, with their quoted error.  Drawn
        after the curves so they sit on top, and as open circles so a
        measurement is never mistaken for a computed line.  The fraction t
        already runs in the direction this segment is drawn - the flip for
        Sigma, measured Gamma->K and drawn K->Gamma, is done in
        add_phonon_curves.py where the path is defined.  */
    /*  A fitted model, drawn as a LINE and never as the circles used for
        measured phonons.  Strontium, chromium and rhodium have one, for
        different reasons - strontium because a powder spectrum plus a
        Born-von Karman fit is all anyone has for it, the other two because
        their papers tabulate the fitted force constants rather than the
        frequencies.  Giving any of them the same symbol as gold's 114
        individually measured frequencies would say something the data
        cannot support; the reason is carried per element in
        model_curve.why.  */
    const mcs=(d.model_curve&&d.model_curve.segs)
              ?d.model_curve.segs[p.a+"|"+p.b]:null;
    if(mcs&&mcs.length){
      c.save();c.strokeStyle=ink;c.globalAlpha=0.55;c.lineWidth=1.3;
      c.setLineDash([4,3]);
      /*  not always three: molybdenum's source tabulates two of the three
          [zz0] branches, so the row is shorter there and the missing one is
          drawn as missing rather than faked  */
      for(let br=1;br<mcs[0].length;br++){
        c.beginPath();
        mcs.forEach((row,k)=>{const x=x0+row[0]*w, y=Y(row[br]);
          k?c.lineTo(x,y):c.moveTo(x,y);});
        c.stroke();}
      c.restore();
    }
    const ecs=(d.exp_curve&&d.exp_curve.segs)?d.exp_curve.segs[p.a+"|"+p.b]:null;
    if(ecs){
      c.save();c.strokeStyle=ink;c.lineWidth=1.1;c.globalAlpha=0.85;
      ecs.forEach(pt=>{
        const x=x0+pt[0]*w, y=Y(pt[1]), e=pt[2]||0;
        if(e>0){const y1=Y(pt[1]-e), y2=Y(pt[1]+e);
          c.beginPath();c.moveTo(x,y1);c.lineTo(x,y2);c.stroke();}
        /*  a point READ from a figure is drawn as a diamond, never with
            the circle a tabulated measurement gets  */
        if(d.exp_curve.digitised){c.beginPath();c.moveTo(x,y-3.3);
          c.lineTo(x+3.3,y);c.lineTo(x,y+3.3);c.lineTo(x-3.3,y);
          c.closePath();c.stroke();}
        else{c.beginPath();c.arc(x,y,2.6,0,2*Math.PI);c.stroke();}});
      c.restore();
    }
    x0 += w + (brk[i]?GAP:0);
  });
  if(usegs) legend(c,L+8,T+14,ink,p2,p3,
                   SECOND_LABEL[DISP_MODE]||"UG",
                   {mau:"MAU",ug:"UG",rc:"re-cut",rc_ug:"re-cut UG"}[PAR_SET]
                   ||"MAU");
}

/*  two-entry key, drawn inside the frame so it cannot be clipped */
function legend(c,x,y,ink,p2,p3,second,first){
  /*  The second curve is not always UG - it is whichever arm the selector
      names - so the key has to be told, or it will label a re-cut candidate
      as the angular potential.  */
  c.font="11px ui-monospace,Consolas,monospace";c.textAlign="left";
  [[first||"MAU",p2,[]],[second||"UG",p3,[5,3]]]
    .forEach(([t,col,dash],i)=>{
    const yy=y+i*14;
    c.strokeStyle=col;c.lineWidth=2;c.setLineDash(dash);
    c.beginPath();c.moveTo(x,yy);c.lineTo(x+20,yy);c.stroke();c.setLineDash([]);
    c.fillStyle=ink;c.fillText(t,x+26,yy+4);
  });
}

/* Directional response, one polar panel per coordinate plane.

   Radius is the value itself, so the shape is the anisotropy: a circle means
   isotropic in that plane.  Poisson's ratio can be negative, so its radius is
   drawn from a floor at min(0, nu_min) and the zero circle is marked. */
const MPROP={
  E:   {keys:["E"],             unit:"GPa",  name:"Young's modulus E"},
  beta:{keys:["beta"],          unit:"1/TPa",name:"linear compressibility",
        scale:1000},
  G:   {keys:["G_min","G_max"], unit:"GPa",  name:"shear modulus G"},
  nu:  {keys:["nu_min","nu_max"],unit:"",    name:"Poisson's ratio"},
};
let MPROP_SEL="E";
/*  Which arm the panel draws BESIDE the published one.  It used to be UG and
    only UG, hard-coded, which was fine while UG was the only other fit; with
    the re-cut candidates on the page a fixed second curve would silently show
    one comparison and label it as the comparison.  */
let MARM = "ug";
/*  every arm of this element that has directional data, as [key, label].
    MAU is in the list: once the main curve follows the selection, a reader
    looking at a re-cut candidate needs the published arm as the thing to
    compare it against, and leaving it out made that the one comparison the
    panel could not draw.  */
const marmRec = (d, k) => k === "mau" ? d
                : k === "ug" ? ((d.ug && d.ug.comparable) ? d.ug : null)
                : d[k];
const marmOpts = d => [["mau", "MAU"], ["ug", "UG"]].concat(
    ARMS.filter(a => a[0] !== "tap")
        .map(a => [a[0], a[1] + armStab(d[a[0]], a[0])]))
  .filter((o, i, all) => all.findIndex(x => x[0] === o[0]) === i)
  .filter(o => {const r = marmRec(d, o[0]); return r && r.mech_planes;});

function drawPolar(d){
  /*  The MAIN curve is the set on screen, not the root record.  It was
      d.mech_planes unconditionally, so selecting a re-cut candidate changed
      the Hill table and the anisotropy beside this panel while the panel
      itself went on drawing MAU - the two disagreeing about the same
      element, side by side.  */
  const Pp=parSet(d);
  const pl=(Pp.mech_planes||d.mech_planes); if(!pl) return;
  const s=setup("#polar",0.36); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const spec=MPROP[MPROP_SEL], sc=spec.scale||1;
  /*  The fourth panel is not a coordinate plane and that is the point.  The
      x-y, x-z and y-z sections between them contain [100], [110] and [101]
      and NO member of <111> - which is where Young's modulus is extremal for
      a cubic crystal, so three panels understate the stiffest direction
      badly and silently.  Measured on this library's own tensors: Rb 2.1
      against 4.8 GPa, W 195 against 401, Cu 130 against 191.  Nothing about
      the three panels looks wrong.  "d" is the plane spanned by [110] and
      [001]; its normal is [110] x [001] = [1 -1 0], so it is (11-0) with the
      bar on the second digit, and <111> sits inside it at 35.26 degrees.  */
  const planes=["xy","xz","yz"].concat(pl["d"]?["d"]:[]);
  /*  Drawn whenever a comparable UG fit exists, not behind the overlay switch:
      the table beside this panel now shows both potentials unconditionally, and
      a plot that disagreed with the table next to it would be a trap.  The
      switch stays for the dispersion and the thermodynamics, where the second
      curve changes how the first is read. */
  /*  The comparison curve is any arm EXCEPT the one already drawn as the main
      one, or the panel would put a potential beside itself.  */
  const opts=marmOpts(d).filter(o=>o[0]!==PAR_SET);
  if(!opts.some(o=>o[0]===MARM)) MARM = opts.length ? opts[0][0] : null;
  const u = MARM ? marmRec(d, MARM) : null;
  const upl = (u && u.mech_planes) ? u.mech_planes : null;
  /* one radial scale for all three panels, so they are comparable - and the
     same scale for both potentials, or the shapes could not be compared */
  let vmax=-Infinity, vmin=Infinity;
  [pl,upl].forEach(src=>{ if(!src) return;
    planes.forEach(p=>spec.keys.forEach(k=>(src[p]||{})[k]&&src[p][k].forEach(v=>{
      const x=v*sc; if(x>vmax)vmax=x; if(x<vmin)vmin=x;})));});
  const floor=Math.min(0,vmin), span=(vmax-floor)||1;
  const pad=26, w=W/planes.length, R=Math.min(w,H)/2-pad;
  c.font="11px ui-monospace,Consolas,monospace";
  planes.forEach((p,ip)=>{
    const cx=w*(ip+0.5), cy=H/2;
    const rad=v=>R*(v*sc-floor)/span;
    /* rings at the ends of the scale, plus zero when the range crosses it */
    c.strokeStyle=lin;c.lineWidth=1;
    [floor,floor+span/2,vmax].forEach(v=>{
      c.beginPath();c.arc(cx,cy,R*(v-floor)/span,0,2*Math.PI);c.stroke();});
    if(floor<0){c.strokeStyle=ink;c.save();c.setLineDash([3,3]);
      c.beginPath();c.arc(cx,cy,R*(-floor)/span,0,2*Math.PI);c.stroke();
      c.restore();}
    c.strokeStyle=lin;
    c.beginPath();c.moveTo(cx-R,cy);c.lineTo(cx+R,cy);
    c.moveTo(cx,cy-R);c.lineTo(cx,cy+R);c.stroke();
    const loop=(arr,col,w,dash)=>{
      const n=arr.length;
      c.strokeStyle=col;c.lineWidth=w;c.setLineDash(dash);
      c.beginPath();
      for(let i=0;i<=n;i++){
        const t=2*Math.PI*(i%n)/n, r=rad(arr[i%n]);
        const x=cx+r*Math.cos(t), y=cy-r*Math.sin(t);
        i?c.lineTo(x,y):c.moveTo(x,y);
      }
      c.stroke();c.setLineDash([]);
    };
    spec.keys.forEach((k,ik)=>{
      /*  min and max already use solid/dashed, so UG cannot be a third dash
          pattern without becoming unreadable: it is the same pattern, lighter
          and thinner, which reads as "the same quantity, other potential". */
      if(upl&&upl[p]&&upl[p][k]){
        c.save();c.globalAlpha=0.85;
        loop(upl[p][k],p3,1.2,ik?[4,3]:[2,3]);c.restore();
      }
      loop(pl[p][k],ik?p3:p2,1.7,ik?[4,3]:[]);
    });
    c.fillStyle=ink;c.textAlign="center";
    /*  fillText does not decode HTML entities - an &ndash; here once reached
        the canvas as the literal text "&ndas" - so the bar is the real
        combining macron U+0304, and it sits on the digit it belongs to.  */
    c.fillText(p==="d" ? "(11̄0) plane" : p[0]+"-"+p[1]+" plane",cx,H-7);
  });
  c.fillStyle=ink;c.textAlign="left";
  c.fillText(`${(vmax).toFixed(vmax<10?2:0)} ${spec.unit}`.trim()+" outer ring"
             +(upl?"   (thin = "
                +((marmOpts(d).find(o=>o[0]===MARM)||[0,"UG"])[1])+")":""),
             6,12);
}

function drawThermo(d){
  /*  same rule as the dispersion: the solid curve is the set on screen  */
  const Pt=parSet(d);
  const g=(Pt.ld&&Pt.ld.thermo)?Pt.ld:(d.ld||{}); if(!g.thermo) return;
  const s=setup("#thermo"); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=48,R=10,T=10,B=26, pw=W-L-R, ph=H-T-B;
  const th=g.thermo.filter(r=>r.S!==undefined&&r.Cv!==undefined);
  if(!th.length) return;
  /*  Every arm that has a thermodynamic curve, not just UG.  S(T) and Cv(T)
      come from the same dispersion the panel above draws, so an arm whose
      phonons differ shows it here - which is the point for the re-cut
      candidates, whose whole disagreement with the published arm is thermal.  */
  const extra = armRecs(d, r => r.ld && r.ld.thermo)
      .filter(o => o[0] !== "tap" && o[0] !== PAR_SET)
      .map(o => [o[2], o[1].ld.thermo.filter(
                          r => r.S !== undefined && r.Cv !== undefined)])
      .filter(o => o[1].length)
      .concat((()=>{const u=ugData(d);
        return (u&&u.ld&&u.ld.thermo)?[["UG", u.ld.thermo.filter(
          r=>r.S!==undefined&&r.Cv!==undefined)]]:[];})());
  const uth = extra.length ? extra[0][1] : null;
  const Tm=th[th.length-1].T;
  let hi=Math.max(...th.map(r=>Math.max(r.S,r.Cv)), d.S298||0, d.Cp298||0,
                  ...extra.flatMap(o=>o[1].map(r=>Math.max(r.S,r.Cv))), 0);
  hi=Math.ceil(hi/10)*10;
  const X=t=>L+t/Tm*pw, Y=v=>T+ph-v/hi*ph;
  c.strokeStyle=lin;c.lineWidth=1;
  c.font="11px ui-monospace,Consolas,monospace";c.fillStyle=ink;
  for(let k=0;k<=4;k++){const v=hi*k/4,y=Y(v);
    c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
    c.textAlign="right";c.fillText(String(Math.round(v)),L-6,y+4);}
  for(let k=0;k<=4;k++){const t=Tm*k/4;
    c.textAlign="center";c.fillText(String(Math.round(t)),X(t),H-16);}
  c.textAlign="center";c.fillText("T (K)",L+pw/2,H-3);
  frame(c, L, T, pw, ph, ink,
        [0,1,2,3,4].map(k=>X(Tm*k/4)),
        [0,1,2,3,4].map(k=>Y(hi*k/4)));
  c.save();c.translate(13,T+ph/2);c.rotate(-Math.PI/2);
  c.textAlign="center";c.fillStyle=ink;
  c.fillText("J / (mol K)",0,0);c.restore();
  //  Dulong-Petit limit
  c.strokeStyle=ink;c.setLineDash([2,3]);
  c.beginPath();c.moveTo(L,Y(24.94));c.lineTo(W-R,Y(24.94));c.stroke();
  c.setLineDash([]);
  c.textAlign="right";c.fillText("3R",W-R-2,Y(24.94)-3);
  const line=(rows,key,col,w,dash)=>{c.strokeStyle=col;c.lineWidth=w;
    c.setLineDash(dash);c.beginPath();
    rows.forEach((r,i)=>{const x=X(r.T),y=Y(r[key]);i?c.lineTo(x,y):c.moveTo(x,y);});
    c.stroke();c.setLineDash([]);};
  const DASHES=[[5,3],[2,2],[6,2,1,2],[1,3]];
  extra.forEach((o,i)=>{const dl=DASHES[i%DASHES.length];
    line(o[1],"S",p2,1.3,dl); line(o[1],"Cv",p3,1.3,dl);});
  line(th,"S",p2,2,[]); line(th,"Cv",p3,2,[]);
  if(extra.length){c.fillStyle=ink;c.textAlign="left";
    c.font="11px ui-monospace,Consolas,monospace";
    extra.forEach((o,i)=>{
      const yy=T+12+i*13, dl=DASHES[i%DASHES.length];
      c.strokeStyle=ink;c.lineWidth=1.3;c.setLineDash(dl);
      c.beginPath();c.moveTo(L+6,yy-4);c.lineTo(L+26,yy-4);c.stroke();
      c.setLineDash([]);
      c.fillText(o[0],L+32,yy);});}
  const ring=(t,v,col)=>{if(!v)return;c.strokeStyle=col;c.lineWidth=2;
    c.beginPath();c.arc(X(t),Y(v),4.5,0,6.284);c.stroke();};
  ring(298,d.S298,p2); ring(298,d.Cp298,p3);
}

/*  pre defaults to "exp" because most of these tiles carry a measured value
    underneath; the mechanical ones carry the UG number instead, and labelling
    a second calculation "exp" would be a lie in a place nobody would check. */
const cell3=(k,v,ref,pre)=>`<div class="par"><span class="k">${k}</span>
  <span class="v mono">${v}</span>${ref?`<span class="k">${pre||"exp"} ${ref}</span>`:""}</div>`;
function tv(g,key){
  if(!g.thermo) return "&mdash;";
  let best=g.thermo[0];
  for(const r of g.thermo) if(Math.abs(r.T-298)<Math.abs(best.T-298)) best=r;
  return best[key]!==undefined?best[key].toFixed(key==="F"?4:2):"&mdash;";
}
/*  relaxed | frozen-ion | DFT (Materials Project) | experiment
    The error column is taken against EXPERIMENT, because that is what the fit
    targets.  MP is a 0 K DFT number and experiment is room temperature, so the
    two reference columns are different quantities, not competing measurements -
    where they disagree with each other, our deviation from either is not by
    itself evidence that the potential is wrong.                              */
/*  MAU against UG for one element.

    Everything MAU carries, UG carries too - the same producer builds both
    records - so this is not a summary but the whole comparison: parameters,
    elastic constants, the mechanical averages the tensor implies, the phonon
    spectrum and the measured frequencies.  The plots above take the overlay.

    Two things are shown together that are easy to present separately and
    misleadingly: the elastic error AND the dynamical stability.  UG reaches
    0.00 % on several cubic metals, and an elastic-only table would read as a
    clean win regardless of what the phonons do.                              */
/*  The LAMMPS layer.  The page described the physics of the cutoff at length
    and never mentioned that a pair style and a set of files exist, which left
    a reader knowing why the switch matters and not knowing they could run it.

    The verdict column is the reason this is here rather than in a README.  A
    fit reproducing its targets says nothing about whether the crystal survives
    being heated - chromium's switched set reproduces its elastic constants and
    turns into a 4559 K liquid - and the two questions had never been asked
    separately.  A user copies a file out of a directory and takes its header
    with it; they do not take the README.  So the warning lives in both.        */
/*  The argument above used to end at the citation.  It no longer has to:
    LAMMPS ships Nb.uf3, which is that potential for niobium, and it went
    through this page's own elastic and surface machinery unchanged, so only
    the potential differed.  A page that says "published work reports" when it
    has its own number is weaker than it needs to be.  */
function uf3Note(nb, here){
  const u = nb && nb.uf3; if(!u) return "";
  const f = k => ["110","100","111"].map(x=>(u.facets[x]||{})[k])
                   .filter(v=>v!=null).map(v=>v.toFixed(2)).join(" / ");
  const ord = k => u.order[k].join(" &lt; ");
  return `<br><br><strong>That is measured here, not only cited.</strong>
  ${here?"For this element":"For niobium"} the two potentials went through
  the same code and only the potential differed &mdash; 110 / 100 / 111 in
  J&nbsp;m<sup>&minus;2</sup>: reference ${f("ref")}, ours ${f("ours")},
  the UF3 ${f("uf3")}. Mean deviation <strong>${u.mean_pct.ours.toFixed(0)} %
  for ours against ${u.mean_pct.uf3.toFixed(0)} % for the UF3</strong>, on a
  property neither was fitted to. The ordering goes the same way: the
  reference gives ${ord("ref")}, the UF3 reproduces ${ord("uf3")} exactly,
  ours gives ${ord("ours")}. ${u.caveat}</p><p class="note">${u.why}
  <br><span style="opacity:.75">${u.ref}</span>`;
}

/*  What to make of the re-cut candidate for THIS element, which is not the
    same answer twice.  The arm looks strong in aggregate - 8.3 % against the
    published switched arm's 11.2 % over the 17 elements that have both and
    were not rejected - and that number is misleading on its own: almost all
    of the gap is lithium, sodium, potassium, rubidium and caesium, and the
    last four are exactly the elements the arm breaks once the crystal is
    warm.  Take those five out and it is 8.6 against 9.2 over 12 elements,
    better by more than half a point in 5.  A wash.  (Recomputed 2026-09-13,
    after the 2026-08-29 lattice-constant correction and with iridium.)

    So the page says it per element rather than quoting the aggregate.  The
    dispersion numbers are the measured comparison, arm against published
    switched arm, as a percentage of that element's highest measured
    frequency; the warm classes come from three independent screens - thermal
    expansion, 300 K elastic constants, and the sign of dC11/dT.  */
/*  curve_mae.py with the arms rc and tap, 2026-09-13.  Li, Pb and Rb had been
    left at their values from before the 2026-08-29 lattice-constant
    correction; a copy of this table in lammps/export_potentials.py was older
    still.  Iridium's points are read from a figure and carry no branch
    labels, so both of its numbers are nearest-branch lower bounds.  */
const RC_DISP = {Ag:[13.0,12.7], Au:[4.2,4.4], Ba:[8.3,9.6], Ca:[8.6,8.3],
                 Cs:[6.4,12.9], Cu:[9.3,10.0], Ir:[9.0,9.9], K:[11.8,13.9],
                 Li:[7.0,25.9], Mg:[8.3,6.6], Na:[5.8,16.2], Ni:[12.1,13.4],
                 Pb:[12.9,18.3], Pd:[3.5,3.6], Pt:[5.0,4.8], Rb:[7.6,11.8],
                 W:[17.2,17.8], Yb:[8.5,8.6]};
const RC_NEAR = ["Ir"];
/*  Distinct measured frequencies in a points-only record.  One Gamma point is
    drawn at every place the path passes Gamma, so counting drawn marks says
    three where rhenium has one.  */
function nMeasured(ec){
  return new Set(Object.values(ec.segs).flatMap(
    v=>v.map(p=>String(p[3])+"|"+p[1]))).size;
}
/*  C_p - C_v = 9 alpha^2 B V_m T for the crystal as measured: the CRC linear
    expansion coefficient (carried in tap.expansion), the bulk modulus of the
    elastic constants the fit targets (B[1]), the molar volume from a0.  It
    turns the tabulated C_p into the quantity the curve is.  null where any
    input is missing.  */
function cpLat(d){
  const al=((d.tap||{}).expansion||{}).alpha_exp_1e6,
        Bx=Array.isArray(d.B)?d.B[1]:null;
  if(!al||!Bx||!d.a0||!d.Cp298) return null;
  const a=d.a0, v=d.struct==="fcc"?a*a*a/4:d.struct==="bcc"?a*a*a/2
        :Math.sqrt(3)/4*a*a*a*(d.c_over_a||1.633);
  return 9*Math.pow(al*1e-6,2)*Bx*1e9*v*1e-30*6.02214076e23*298.15;
}
/*  Lithium added 2026-09-13.  Its re-cut expands at -37e-6/K against a
    measured +46 and its C11 rises x4.6 across its grid - the alkali failure,
    and it had been left off this list.  Which element the re-cut is
    RECOMMENDED for is not a list typed here: recommend_recut.py computes
    `md_recommended` from the library by a rule stated in the record, which is
    what keeps it from going stale the way RC_DISP did.  */
const RC_BROKEN = ["Na","K","Rb","Cs","Li"];  /* warm screens: catastrophic  */
const RC_SUSPECT = ["Yb"];                   /* the same failure, milder    */
const RCREC = d => !!(d && d.md_recommended && d.md_recommended.set === "rc");
/*  The warm failure in this element's own numbers.  The sentence used to quote
    fixed ranges - "-36 to -42", "84-94 % against 12-18 %" - which no longer
    matched the stored data and could not cover lithium, whose grid stops at
    68 K.  */
function rcWarm(d){
  const et = (d.elasticT||{}).rc, x = (d.rc||{}).expansion||{},
        f = (d.md_recommended||{}).facts||{}, bits = [];
  if(x.alpha_1e6 != null && x.alpha_1e6 <= 0)
    bits.push(`thermal expansion comes out NEGATIVE,
      ${x.alpha_1e6.toFixed(0)}&times;10<sup>&minus;6</sup>/K${
      x.alpha_exp_1e6 != null ? ` against a measured +${x.alpha_exp_1e6.toFixed(0)}` : ""}`);
  if(et && et.pts && et.pts.length > 1){
    /*  Molecular-dynamics points only.  The table's T = 0 row is the fit's
        static target, and the jump from it to the first finite-temperature
        point is internal relaxation, not a temperature dependence - which
        is what had made barium look as if it stiffened.  */
    const pts = et.pts.filter(p => !p.above_melt && p.T > 0);
    if(pts.length > 1){
      const p0 = pts.reduce((a,p) => p.T < a.T ? p : a);
      const p3 = pts.reduce((a,p) => Math.abs(p.T-300) < Math.abs(a.T-300) ? p : a);
      if(p3 !== p0 && p3.C11 > p0.C11)
        bits.push(`C<sub>11</sub> rises &times;${(p3.C11/p0.C11).toFixed(2)} from
          ${p0.T.toFixed(0)} to ${p3.T.toFixed(0)}&nbsp;K instead of softening,
          so dC<sub>11</sub>/dT has the wrong SIGN`);
    }
  }
  if(f.e300_rc != null && f.e300_tap != null)
    bits.push(`the elastic constants near 300&nbsp;K land
      ${f.e300_rc.toFixed(0)}&nbsp;% off the room-temperature values, against
      ${f.e300_tap.toFixed(0)}&nbsp;% for the published switched set`);
  return bits.length ? bits.join("; ") : "see the screens below";
}
/*  Why an element has no candidate at all.  It is not that the fit was tried
    and failed - a cutoff has to sit in a GAP between neighbour shells, and
    in two of the three structures there is nowhere to put one.  Shells in
    neighbour-distance units: fcc 1.000 1.414 1.732, bcc 1.000 1.155 1.633,
    hcp 1.000 1.019 1.428.  A 0.019 gap is not a place for a cutoff; it is
    two shells with a cutoff on top of both.  */
const RC_NOGAP = {
  bcc: "1.000 and 1.155 &mdash; the first two shells nearly touch",
  hcp: "1.000 and 1.019 &mdash; the densest of the three, and c/a-dependent",
  fcc: "1.000, 1.414 and 1.732 &mdash; wide and evenly spaced"};

function forceNote(d){
  if(d.tap_force) return "";
  return `<p class="plotnote"><strong>No force-matched record exists for
    ${d.name}.</strong> That is a gap in what has been computed, not a fit
    that failed: force matching needs a density-functional force dataset for
    the element and a fitting run of its own, and this arm covers
    <strong>seven of the thirty-eight</strong> &mdash; copper, nickel,
    niobium, tantalum, tungsten, molybdenum and vanadium. Where a panel shows
    nothing for it, nothing has been measured.</p>`;
}

/*  A force-matched record that IS present still needs a sentence, because
    four of the seven would not pass the library's own acceptance test.  The
    ceiling is not moved to admit them and the violation travels in the
    record; this says so where the numbers are read rather than in a caption
    somewhere else.  */
/*  What the arm has been measured ON, read off the record.  A fixed list
    would be a sentence that starts lying the day the surface run is folded
    in; this one follows the data.  */
const FORCE_PANELS = [["Cij","the elastic constants"],
                      ["dyn","dynamical stability"],
                      ["compression","the compression screen"],
                      ["surface","surface energy"],
                      ["stacking","the stacking fault"],
                      ["ground","the ground state"],
                      ["bain","the Bain path"],
                      ["expansion","thermal expansion"],
                      ["md_screen","the MD screen"],
                      ["jiggle","the nudge test"]];

function forceCoverageNote(d){
  const r = d.tap_force; if(!r) return "";
  const has = FORCE_PANELS.filter(p=>r[p[0]]!=null).map(p=>p[1]);
  if(d.elasticT&&d.elasticT.tap_force) has.push("elastic constants against temperature");
  const lacks = FORCE_PANELS.filter(p=>r[p[0]]==null).map(p=>p[1]);
  if(!lacks.length||!has.length) return "";
  const list=a=>a.length<2?a[0]:a.slice(0,-1).join(", ")+" and "+a[a.length-1];
  return `<p class="plotnote"><strong>The force-matched arm is measured on
    ${list(has)}, and on nothing else here.</strong> It carries no
    ${list(lacks)}, so those panels draw ${d.sym}'s other arms and not this
    one. Nothing failed: those runs have not been done for this arm, and an
    absent curve means an absent measurement, never a rejected result.</p>`;
}

function forceGateNote(d){
  const r = d.tap_force;
  if(!r) return "";
  if(r.gate === "pass")
    return `<p class="plotnote">${d.sym}'s force-matched record sits
      <strong>inside</strong> the three-body ceiling &mdash;
      e<sub>3</sub>/e<sub>2</sub> = ${r.ratio.toFixed(4)} with a violation of
      exactly zero &mdash; so the shipped fitting script would have kept it.</p>`;
  return `<p class="plotnote"><strong>${d.sym}'s force-matched record is over
    the three-body ceiling.</strong> e<sub>3</sub>/e<sub>2</sub> =
    ${r.ratio.toFixed(4)} against a limit of 0.30, a violation of
    ${r.violation.toFixed(4)} eV/atom, and <code>standalone/fit.py</code>
    rejects on that comparison with no tolerance. It is shown because the arm
    is a separate reference frame rather than a candidate for the library, and
    because niobium, tantalum and vanadium sit strictly inside the same
    ceiling &mdash; so 0.30 is not a wall every element is driven into, and
    moving it to admit these four would delete that.</p>`;
}

function recutNote(d){
  const has = d.rc || d.rc_ug;
  const el = d.sym, nm = d.name;
  if(!has){
    return `<p class="plotnote"><strong>No re-cut candidate exists for
      ${nm}.</strong> That is a property of its crystal rather than a failed
      fit: a re-cut moves the pair cutoff into a GAP between neighbour shells,
      and ${d.struct} shells sit at ${RC_NOGAP[d.struct]||"distances that leave no usable gap"} in neighbour-distance units. A gap of two per cent
      is not a place to put a cutoff; it is two shells with a cutoff on top of
      both. Twenty-five of the thirty-eight elements have one &mdash; every
      fcc metal, eight of thirteen bcc, four of twelve hcp.</p>`;
  }
  const rej = armBad(d.rc) || armBad(d.rc_ug);
  if(rej){
    return `<p class="plotnote"><strong>${el}'s re-cut records failed
      selection.</strong> A mode near &minus;13 cm<sup>&minus;1</sup> on the
      symmetry path${el==="Mo"?", and fcc 209 meV/atom below bcc":""}. They
      are kept as the control the accepted candidates are read against, and
      are not to be used for anything.</p>`;
  }
  const dsp = RC_DISP[el];
  const cmp = dsp ? (dsp[0] < dsp[1] - 0.5 ? "better than"
                   : dsp[0] > dsp[1] + 0.5 ? "worse than" : "the same as") : null;
  const nPts = (d.exp_curve&&d.exp_curve.points_only)?nMeasured(d.exp_curve):0;
  const num = dsp ? `Against the measured dispersion it reaches
      <strong>${dsp[0].toFixed(1)}&nbsp;%</strong> where the published
      switched arm reaches ${dsp[1].toFixed(1)}&nbsp;% &mdash; ${cmp} it.${
      RC_NEAR.indexOf(el)>=0?` Both are nearest-branch scores against points
      read from a figure, so both are lower bounds.`:""}`
    : nPts ? `${el} has ${nPts===1?"one measured frequency":nPts+" measured frequencies"},
      not a dispersion, so that comparison is not made here.`
    : `${el} has no measured dispersion, so that comparison cannot be made
      here.`;
  if(RC_BROKEN.indexOf(el) >= 0){
    return `<p class="plotnote"><strong>${el}: do not use the re-cut at finite
      temperature.</strong> ${num} And it fails the warm screens:
      ${rcWarm(d)}. The dispersion and the failure are one fact, not two:
      a dispersion at 0&nbsp;K and an elastic constant are both properties of
      the curvature AT the minimum, and this arm has that neighbourhood right
      while the shape of the well away from it is wrong. Static properties
      only, and with care.</p>`;
  }
  if(RC_SUSPECT.indexOf(el) >= 0){
    return `<p class="plotnote"><strong>${el}: treat the re-cut as
      suspect.</strong> ${num} It shows a milder form of the failure that
      rules the alkalis out &mdash; check thermal expansion and the 300&nbsp;K
      elastic constants before using it warm.</p>`;
  }
  const mr = d.md_recommended || {}, f = mr.facts || {};
  if(RCREC(d)){
    return `<p class="plotnote"><strong>${el}: the re-cut is the RECOMMENDED
      set for molecular dynamics</strong> &mdash; <code>${el}_recut.ugur</code>
      in place of <code>${el}_taper.ugur</code>, which still ships. It passes
      every cold and warm screen${f.isf_rc != null && f.isf_tap != null
        ? `, and it gets the sign of the intrinsic stacking fault right:
          ${f.isf_rc.toFixed(0)} mJ&nbsp;m<sup>&minus;2</sup> against
          ${f.isf_tap.toFixed(0)} for the switched set` : ""}${
      f.e300_rc != null && f.e300_tap != null
        ? `; near 300&nbsp;K its elastic constants are
          ${f.e300_rc.toFixed(1)}&nbsp;% off the room-temperature values
          against ${f.e300_tap.toFixed(1)}&nbsp;%` : ""}. ${num} Not chosen by
      hand: <code>recommend_recut.py</code> applies one rule to every element
      that has a re-cut &mdash; ${mr.rule}. The angular re-cut,
      <code>${el}_recut.ugur.ang</code>, is not part of the recommendation.</p>`;
  }
  return `<p class="plotnote"><strong>${el}: the re-cut is usable, and it is
    not recommended.</strong> ${mr.why && mr.why.length
      ? `It misses the rule on ${mr.why.join("; ")}.` : ""} ${num} The
    published switched set, <code>${el}_taper.ugur</code>, stays the set to
    use for molecular dynamics.</p>`;
}

/*  The finite-temperature dispersion, measured rather than computed.

    Everything else on this page that carries a phonon is a 0 K harmonic
    calculation from the force constants.  This is the spectrum LAMMPS
    measures by displacement correlation during equilibrium molecular
    dynamics at the temperature the neutron experiment was actually run at,
    which is the thing the experiment sees: renormalised by anharmonicity and
    by whatever thermal expansion the potential predicts.

    WHICH ARM, first, because it is not the one the selector shows.  Every
    number here is the SWITCHED record, <El>_taper.ugur, while the
    parameter-set selector's "MAU" is the hard-truncated root record.  The
    switched arm is used because it is the only one that can be run at
    temperature: a hard cut leaves the pair energy discontinuous and copper
    drifts 350 meV/atom/ns in constant-energy dynamics, reaching 1100 K from
    a 296 K start inside 200 ps.

    THE FLOOR decides what may be read.  fix phonon imposes no symmetry on
    what it inverts, so two wavevectors the crystal forces to be equal do not
    come out equal; every q is scored twice and the spread is the run's
    resolution.  The verdict column is computed from it, because twelve of
    the twenty-four differences are smaller than it and a reader subtracting
    two printed numbers would rank the library on noise.                    */
const FT_WHY = {
  cold: T => `measured at ${T} K, far below the Debye temperature &mdash; `
           + `classical dynamics describes the wrong physics there, and a `
           + `finite-temperature run would be answering a different question`,
  screen: () => `the molecular-dynamics screen fails: this record does not `
           + `hold its structure when heated, so there is nothing to measure `
           + `a spectrum on`,
  /*  The element's OWN reason, the same string the dispersion panel above
      prints, because the six elements in this group are in six different
      positions and one sentence about all of them is unfair to most.  No
      single crystal of strontium large enough for triple-axis spectroscopy
      has ever been grown; molybdenum's frequencies WERE measured and the
      paper simply tabulated Fourier coefficients instead of them.  Saying
      the same thing about both told the reader neither.
      And the closing clause is the general rule rather than a special
      pleading for this panel: a reconstruction is drawn everywhere in this
      library and scored nowhere, at 0 K as much as at temperature.  */
  model: (T, d) => {
    const w = d && d.model_curve && d.model_curve.why;
    const lead = w ? (/^[A-Z][a-z]/.test(w) ? w[0].toLowerCase() + w.slice(1) : w)
                   : `the reference here is a reconstruction`;
    //  a full stop, not another clause: strontium's own reason already ends
    //  in "so its phonons were reached by fitting to a powder", and a second
    //  "so" in the same sentence makes the reader read the join twice
    return lead + `. What exists for it is therefore a RECONSTRUCTED `
         + `dispersion rather than individually measured frequencies, and a `
         + `reconstruction is drawn throughout this library, never scored, `
         + `at 0 K as much as at temperature`;
  },
  /*  entered after the study ran, with a run queued: the absence is of a
      RESULT, not of a measurement, and the sentence must say which  */
  pending: T => `a measured dispersion for this element was entered after the `
           + `study ran, read from a published figure; a run at the `
           + `measurement's own ${T} K has been submitted and its result is `
           + `not in yet`,
  gamma: T => `the only measurement is one frequency, the Raman-active optic `
           + `mode at &Gamma;, so there is no dispersion to set a curve `
           + `against; a small-cell run at ${T} K that measures that one `
           + `frequency has been submitted and its result is not in yet`,
  unrun: T => `a measured dispersion for this element was entered after the `
           + `study ran, read from a published figure, and no run at the `
           + `measurement's ${T} K has been made against it yet`,
  none: () => `no dispersion reference exists for this element, measured or `
           + `reconstructed`,
};
/*  kinds that are not outside the study, only not in it yet  */
const FT_WAIT = k => k === "pending" || k === "gamma" || k === "unrun";
const FT_V = {gain: ["gain", "ok"], loss: ["loss", "bad"],
              flat: ["below the floor", "warn"]};

/*  The finite-temperature dispersion.  finiteTBlock is the panel; this is its
    picture.  The dashed 0 K curve is `tap.ld`, the SAME switched record the
    run measures (add_ld_tap.py).  It used to be `d.ld`, the top-level hard-cut
    record, which set a switched run against a hard-cut calculation - and at
    0 K those two differ by more than the temperature effect in 15 of the 24
    elements.  ld.std is cm^-1 and phana returns THz, so the harmonic curve is
    converted on the way in and the axis is THz throughout - the same axis the
    dispersion panel uses.                                                   */
function drawFT(d){
  const cv=(FT.curve||{})[d.sym], segs=((d.tap||{}).ld||{}).std;
  if(!cv||!segs||!segs.length) return;
  const s=setup("#ftdisp",0.60); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=58,R=16,T=20,B=34, pw=W-L-R, ph=H-T-B;
  const G=String.fromCharCode(915), lab=t=>t==="G"?G:t;
  const GAP=10, npt=cv.npt, sub=new Set(cv.sub||[]);

  /*  One range over BOTH curves.  Two auto-scaled axes would make a softening
      look like an artefact of the scale, which is the one thing this plot
      exists to show.  The MEASURED curve never opens it downwards - the
      substitution near Gamma removed its last negative frequency.  The 0 K
      reference can: it is the switched record's own harmonic dispersion, and
      for iron, niobium, tungsten, rubidium and caesium that record has an
      imaginary mode on the path.  That is a property of the record, not a drawing
      fault, so the axis opens under the same rule as the dispersion panel:
      shallower than 1 cm^-1 is arithmetic and floored to zero, deeper is
      drawn.                                                                 */
  let hi=0, lo=0;
  segs.forEach(p=>p.branches.forEach(b=>b.forEach(v=>{
    if(v===null) return;
    if(THZ(v)>hi) hi=THZ(v); if(THZ(v)<lo) lo=THZ(v);})));
  if(lo > -1.0/CM1_PER_THZ) lo = 0;
  cv.f.forEach(r=>r.forEach(v=>{if(v>hi)hi=v;}));
  /*  and the measured points drawn on top, with their error bars  */
  if(d.exp_curve&&d.exp_curve.segs) Object.values(d.exp_curve.segs).forEach(
    v=>v.forEach(pt=>{const t=pt[1]+(pt[2]||0); if(t>hi)hi=t;}));
  const step=niceStep(hi-lo,5);
  hi=Math.ceil(hi/step)*step||step;
  lo=lo<0?-Math.ceil(-lo/step)*step:0;
  const dec=(step<1||Math.abs(step-Math.round(step))>1e-9)?1:0;
  const Y=v=>T+ph-(v-lo)/(hi-lo)*ph;

  c.font="11px ui-monospace,Consolas,monospace";
  if(lo<0){            /* omega = 0, solid like the rest of the frame */
    c.strokeStyle=ink;c.lineWidth=1;
    c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();
  }
  for(let k=Math.round(lo/step);k<=Math.round(hi/step);k++){
    const v=k*step, y=Y(v);
    c.strokeStyle=ink;c.lineWidth=1.2;
    c.beginPath();c.moveTo(L,y);c.lineTo(L+5,y);c.stroke();
    c.beginPath();c.moveTo(W-R,y);c.lineTo(W-R-5,y);c.stroke();
    c.fillStyle=ink;c.textAlign="right";c.fillText(v.toFixed(dec),L-8,y+4);
  }
  c.strokeStyle=ink;c.lineWidth=1.2;c.strokeRect(L,T,pw,ph);
  c.fillStyle=ink;c.textAlign="center";
  c.save();c.translate(15,T+ph/2);c.rotate(-Math.PI/2);
  c.fillText("Frequency (THz)",0,0);c.restore();

  const vline=(x,name)=>{
    c.save();c.strokeStyle=ink;c.lineWidth=1;c.globalAlpha=0.55;
    c.setLineDash([5,4]);
    c.beginPath();c.moveTo(x,T);c.lineTo(x,T+ph);c.stroke();c.restore();
    c.strokeStyle=ink;c.lineWidth=1.2;
    c.beginPath();c.moveTo(x,T+ph);c.lineTo(x,T+ph-6);c.stroke();
    c.beginPath();c.moveTo(x,T);c.lineTo(x,T+6);c.stroke();
    c.textAlign="center";c.fillStyle=ink;c.fillText(lab(name),x,T+ph+17);
  };

  /*  Same segment layout as the dispersion panel: width in proportion to the
      segment's length in reciprocal space, and a gap rather than a straight
      line across a genuine discontinuity of the path.                       */
  const brk=segs.map((p,i)=>i<segs.length-1&&p.b!==segs[i+1].a);
  const tot=segs.reduce((a,p)=>a+p.len,0);
  const avail=pw-GAP*brk.filter(Boolean).length;
  let x0=L;
  segs.forEach((p,i)=>{
    const w=avail*p.len/tot;
    vline(x0,p.a);
    if(i===segs.length-1||brk[i]) vline(x0+w,p.b);

    /*  0 K first and dashed, so the finite-temperature curve sits on top of
        it.  ld.std carries sixty points a segment against the curve's twenty
        and keeps all sixty here - it costs nothing, it is already in the
        page, and drawing the reference at the coarser sampling would blunt
        exactly the sharp features the comparison is about.                  */
    /*  p3 dashed, not the gridline colour.  --line-2 is #DCE3EB, which is
        what a gridline is meant to be and leaves the reference curve
        invisible against the panel; the page's own grammar for "the thing
        being compared against" is already p3 dashed, in the dispersion
        panel directly above.                                              */
    const n=p.n;
    c.save();c.setLineDash([5,3]);c.strokeStyle=p3;c.lineWidth=1.3;
    c.globalAlpha=0.85;
    p.branches.forEach(b=>{c.beginPath();let on=false;
      b.forEach((v,k)=>{if(v===null){on=false;return;}
        const x=x0+(n>1?k/(n-1):0)*w, y=Y(THZ(v));
        if(on){c.lineTo(x,y);}else{c.moveTo(x,y);on=true;}});
      c.stroke();});
    c.restore();

    /*  The measured curve, drawn one interval at a time so the dash can
        change mid-branch: inside the first mesh cell nothing was resolved
        and the value there is the leading behaviour of the branch rather
        than a measurement.  Solid for measured, dotted for substituted, and an
        interval touching a substituted end is dotted - it is the join that
        is uncertain, not one of its two ends.                               */
    const nb=cv.f[0].length;
    c.strokeStyle=p2;c.lineWidth=1.8;
    for(let b=0;b<nb;b++)
      for(let k=1;k<npt;k++){
        const ia=i*npt+k-1, ib=i*npt+k;
        const va=cv.f[ia][b], vb=cv.f[ib][b];
        if(va===null||vb===null||va===undefined||vb===undefined) continue;
        c.save();c.setLineDash((sub.has(ia)||sub.has(ib))?[2,2.5]:[]);
        c.beginPath();
        c.moveTo(x0+(k-1)/(npt-1)*w,Y(va));
        c.lineTo(x0+k/(npt-1)*w,Y(vb));
        c.stroke();c.restore();
      }

    /*  and the neutron points on top, as open circles with their error bar,
        the same symbol the dispersion panel uses - a measurement must never
        be drawn with the same mark as a computed line.  The fraction t
        already runs in the direction this segment is drawn.                 */
    const ecs=(d.exp_curve&&d.exp_curve.segs)?d.exp_curve.segs[p.a+"|"+p.b]:null;
    if(ecs){
      c.save();c.strokeStyle=ink;c.lineWidth=1.1;c.globalAlpha=0.85;
      ecs.forEach(pt=>{
        const x=x0+pt[0]*w, y=Y(pt[1]), e=pt[2]||0;
        if(e>0){c.beginPath();c.moveTo(x,Y(pt[1]-e));c.lineTo(x,Y(pt[1]+e));
          c.stroke();}
        /*  a point READ from a figure is drawn as a diamond, never with
            the circle a tabulated measurement gets  */
        if(d.exp_curve.digitised){c.beginPath();c.moveTo(x,y-3.3);
          c.lineTo(x+3.3,y);c.lineTo(x,y+3.3);c.lineTo(x-3.3,y);
          c.closePath();c.stroke();}
        else{c.beginPath();c.arc(x,y,2.6,0,2*Math.PI);c.stroke();}});
      c.restore();
    }
    x0+=w+(brk[i]?GAP:0);
  });

  /*  A key inside the frame.  The caption below says the same thing in words,
      but a figure that has to be read together with a paragraph underneath it
      is a figure that will be read wrong once it is on its own.             */
  const rT=(FT.rows||{})[d.sym];
  const keys=[[p2,[],(rT?rT.T:"")+" K, molecular dynamics"],
              [p3,[5,3],"0 K, harmonic, same record"]];
  /*  On an opaque patch of the panel's own background.  Titanium's highest
      branch runs through this corner, and a key drawn straight onto the
      curves is a key that has to be deciphered.                            */
  c.textAlign="left";
  const kw=Math.max(...keys.map(k=>c.measureText(k[2]).width))+44;
  c.fillStyle=getComputedStyle(document.documentElement)
              .getPropertyValue("--surface").trim();
  c.fillRect(L+1,T+1,kw,keys.length*14+8);
  keys.forEach(([col,dash,txt],k)=>{
    const yy=T+14+k*14;
    c.strokeStyle=col;c.lineWidth=2;c.setLineDash(dash);
    c.beginPath();c.moveTo(L+10,yy);c.lineTo(L+30,yy);c.stroke();
    c.setLineDash([]);c.fillStyle=ink;c.fillText(txt,L+36,yy+4);
  });
}


/*  Counted from finiteT.json, never typed: these sentences said
    "twenty-four" and "twelve" for a study that has since grown to 27
    elements, 15 of them below the floor.  */
function ftCount(){
  const rows = Object.values(FT.rows||{});
  return {rows: rows.length, out: Object.keys(FT.out||{}).length,
          all: rows.length + Object.keys(FT.out||{}).length,
          flat: rows.filter(r => r.v === "flat").length};
}

function finiteTBlock(d){
  const el = d.sym, nm = d.name,
        r = (FT.rows||{})[el], o = (FT.out||{})[el];
  if(!r && !o) return "";
  if(!r){
    return `<h3>The dispersion at the temperature it was measured</h3>
      <p class="plotnote"><strong>${nm} ${FT_WAIT(o.kind)
        ?"is not in this study yet":"is outside this study"}</strong> &mdash;
      ${(FT_WHY[o.kind]||FT_WHY.none)(o.T, d)}. ${FT_WAIT(o.kind)
      ?`${ftCount().rows} of the ${ftCount().all} elements carry a
      finite-temperature measurement so far; ${o.kind==="unrun"
      ?"no run has been made for this one yet":"this one is waiting on its run"}.`
      :`${ftCount().rows} of the ${ftCount().all}
      elements carry a finite-temperature measurement; this is one of the
      ${ftCount().out} that do not, and that is a fact about the reference data
      rather than about the potential.`}</p>`;
  }
  const [word, cls] = FT_V[r.v];
  const cv = (FT.curve||{})[el], hasCurve = !!cv && !!(((d.tap||{}).ld||{}).std);
  const amRow = (r.am !== r.a0) ? `<tr><td>0 K, at the paper's own lattice
      constant</td><td>${r.am.toFixed(1)} %</td></tr>` : "";
  return `<h3>The dispersion at the temperature it was measured</h3>
  ${hasCurve?`<canvas id="ftdisp"></canvas>
  <p class="plotnote">The same path and the same axis as the dispersion above.
    <strong>The dashed curve is the 0&nbsp;K harmonic dispersion of this same
    switched record</strong>, not the hard-cut curve drawn above, so the two
    temperatures are read against each other with nothing else changed.${(()=>{
      const m=Math.min(...d.tap.ld.std.flatMap(p=>p.branches.flat()));
      return m<=-1?` At 0&nbsp;K this record has an <strong>imaginary mode on
        the path</strong>, down to ${m.toFixed(1)}&nbsp;cm<sup>&minus;1</sup>,
        which is why the dashed curve dips below zero.`:"";})()} ${
    d.exp_curve&&d.exp_curve.digitised
      ?`The diamonds are the measurement itself, read from the published
        figure because the paper tabulates no frequencies,`
      :d.exp_curve&&d.exp_curve.points_only
      ?`The circles are the measured frequencies &mdash; single points, not a
        dispersion &mdash;`
      :`The open circles with their error bars are the neutron measurement
        itself,`}
    at ${d.exp_curve?d.exp_curve.T_K:r.T}&nbsp;K.
    <br><strong>The dotted stretch beside &Gamma; was not measured.</strong>
    The mesh is fixed by the box &mdash; ${r.mesh} cells &mdash; and inside
    the first mesh cell there is nothing to interpolate between, so each
    branch is carried in from the first mesh shell along its
    own direction by its own leading behaviour: the three acoustic branches
    by &omega;&nbsp;=&nbsp;v|q|, which is exact as
    q&nbsp;&rarr;&nbsp;0${cv.f[0].length>3?`, and the optic branches by
    &omega;&nbsp;=&nbsp;&omega;(&Gamma;)&nbsp;+&nbsp;&Delta;(q/q<sub>s</sub>)<sup>2</sup>,
    which is flat at &Gamma; as an optic branch must be. &Gamma; itself is a
    mesh point and keeps its measured value`:""}. That is the physics rather
    than a fit to it, but it is not a measurement, and drawing it solid would
    pass it off as one.${d.exp_curve&&d.exp_curve.ref?`
    <br>Source: ${d.exp_curve.ref}.`:""}</p>`:""}
  <table class="ft"><tbody>
    <tr><td>0 K harmonic, at the library's lattice constant</td>
        <td>${r.a0.toFixed(1)} %</td></tr>
    ${amRow}
    <tr><td><strong>measured by molecular dynamics at ${r.T} K</strong></td>
        <td><strong>${r.ft.toFixed(1)} %</strong></td></tr>
    <tr><td>the run's own symmetry floor</td>
        <td>${r.floor.toFixed(1)} %</td></tr>
    <tr class="${cls}"><td>difference, against that floor</td>
        <td>${r.d>0?"+":""}${r.d.toFixed(1)} &mdash; <strong>${word}</strong></td></tr>
  </tbody></table>
  <p class="plotnote">Per cent of the highest measured frequency, over
    ${r.n} point${r.n===1?"":"s"}. ${r.mesh} cells, ${(r.steps/1e6).toFixed(1)} million steps
    of 2 fs under a barostat, so the volume is the potential's own at that
    temperature and not one imposed on it.
    <strong>This is the switched record</strong>
    (<code>${el}_taper.ugur</code>), which is NOT what the parameter-set
    selector above calls MAU &mdash; that is the hard-truncated root record.
    The switched arm is the only one that can be run at temperature at all:
    a hard cut leaves the pair energy discontinuous, and copper drifts
    350 meV per atom per nanosecond under it.${r.note?`
    <br>${r.note}`:""}${RCREC(d)?`
    <br><strong>For ${el} the recommended set for molecular dynamics is now
    the re-cut</strong>, <code>${el}_recut.ugur</code>, and this study has not
    run it: what is measured here is the switched file that still ships
    beside it.`:""}
    ${r.v==="flat"?`<br>The difference here is smaller than what the run can
      resolve, so it is <strong>not a result</strong> in either direction.
      ${ftCount().flat} of the ${ftCount().rows} elements are in that position;
      the floor is printed so that the reader does not have to guess which.`:""}</p>`;
}


function lammpsBlock(d){
  /*  The "what" column carries the USE, not just the name.  Which truncation
      a reader wants is not a matter of taste and the two answers point
      opposite ways: measured against neutron dispersion the hard sets average
      9.7 % over the 33 elements that have a curve and the switched ones
      12.6 %, hard being closer in 25 of the 33 - and no hard set can be run
      at temperature, because phi2 does not vanish at the cutoff.  Copper
      drifts 350 meV per atom per nanosecond in NVE with the hard set and 0.4
      with the switched one, a factor of 876.  Each file's own header carries
      its own discontinuity, which runs from 0.25 meV for palladium to 123 for
      yttrium, so "not for dynamics" is not equally true across the library.  */
  const SETS = [
    ["hard", d,        d.sym+".ugur",             "ugur",
     "hard truncation &mdash; static properties and lattice dynamics, "
     + "<strong>not molecular dynamics</strong>"],
    ["tap",  d.tap,    d.sym+"_taper.ugur",       "ugur",
     (RCREC(d)
       ? "switched &mdash; runs MD, but for " + d.sym + " the re-cut below is "
         + "<strong>recommended instead</strong>"
       : "switched &mdash; <strong>the set to use for MD</strong>; "
         + "about 3 points worse on the measured dispersion")],
    ["ug",   d.ug,     d.sym+".ugur.ang",         "ugur/ang",
     "hard + angular &mdash; static properties and lattice dynamics, "
     + "<strong>not molecular dynamics</strong>"],
    ["tap_ug", d.tap_ug, d.sym+"_taper.ugur.ang", "ugur/ang",
     "switched + angular &mdash; <strong>MD</strong>; about 3 points worse "
     + "on the measured dispersion"],
    /*  The re-cut ships under its own stem.  Where recommend_recut.py finds
        it passes every cold and warm screen it is the RECOMMENDED set for
        molecular dynamics; elsewhere it is a candidate, shipped so the
        comparison on this page can be checked with the file it was made
        with.  The angular re-cut is never the recommendation.  */
    ["rc",    d.rc,    d.sym+"_recut.ugur",     "ugur",
     (RCREC(d)
       ? "re-cut, switched &mdash; <strong>recommended for MD</strong> for "
         + d.sym
       : "re-cut candidate, not recommended"
         + (armBad(d.rc)?" — REJECTED":""))],
    ["rc_ug", d.rc_ug, d.sym+"_recut.ugur.ang", "ugur/ang",
     "re-cut + angular, candidate"
     + (armBad(d.rc_ug)?" — REJECTED":"")],
    /*  The two other candidate sets ship for the same reason the re-cut ones
        do: a reader who wants to check a comparison made on this page needs
        the file it was made with.  Neither is part of the published library
        and both rows say so.  */
    ["hard_disp", d.hard_disp, d.sym+"_disp.ugur", "ugur",
     "dispersion-selected candidate, not published &mdash; hard truncation, "
     + "<strong>not molecular dynamics</strong>"],
    ["tap_nudge", d.tap_nudge, d.sym+"_nudge_taper.ugur", "ugur",
     "nudge-constrained candidate, not published &mdash; switched, holds its "
     + "lattice where the shipped one does not"]
  ].filter(r=>r[1]);
  if(!SETS.length) return "";
  const V = r => (r[1] && r[1].md_screen) || null;
  const bad = SETS.filter(r=>{const v=V(r); return v && v.collapsed;});
  const rows = SETS.map(r=>{
    const v = V(r);
    let verdict, cls;
    if(!v){ verdict = "not screened"; cls = ""; }
    else if(v.lost){ verdict = `<strong>do not use</strong> &mdash; the crystal
      disintegrates, ${v.lost[0]} atoms become ${v.lost[1]}`; cls = "bad"; }
    else if(v.collapsed){ verdict = `<strong>do not use</strong> &mdash; collapses
      to ${v.T} K, energy below the static lattice`; cls = "bad"; }
    else if(v.T > 400){ verdict = `suspect &mdash; runs at ${v.T} K where
      equipartition gives 300`; cls = "warn"; }
    else verdict = `holds its structure (${v.T} K)`;
    return `<tr${cls?` class="${cls}"`:""}><td><code>${r[2]}</code></td>
      <td><code>${r[3]}</code></td><td>${r[4]}</td><td>${verdict}</td></tr>`;
  }).join("");
  const rcnote = recutNote(d);
  const fnote = forceNote(d) + forceGateNote(d) + forceCoverageNote(d);
  return `<h3>Running this in LAMMPS</h3>
  ${bad.length?`<div class="warn" style="border-left-color:var(--bad)">
    <strong>${bad.length===SETS.length?"None of the":
      bad.length===1?"One of the":`${bad.length} of the`} parameter set${
      bad.length===1?"":"s"} for ${d.name} can be used for molecular
    dynamics.</strong> They reproduce the elastic constants they were fitted
    to and the crystal does not survive being heated to 600 K: the potential
    energy ends up below the static lattice, which means the structure has
    become something else. Static properties from those sets are unaffected.
    The 0 K stability screen inside the fit cannot see this &mdash; it asks
    whether the reference structure is a harmonic minimum, not whether it
    survives thermal excitation.</div>`:""}
  <table class="tbl"><thead><tr><th>file</th><th>pair_style</th>
    <th>parameters</th><th>molecular dynamics</th></tr></thead>
    <tbody>${rows}</tbody></table>
  ${rcnote}
  ${fnote}
  <p class="note">Files are in <code>lammps/potentials/</code>. The extension
  says which pair style the file needs, the way <code>.eam</code> and
  <code>.eam.alloy</code> do; the stem says which parameter set it is, a plain
  name being hard-truncated and <code>_taper</code> switched. Those are
  different potentials fitted under different truncations and must not be
  swapped. The verdict is a six-picosecond screen at 600 K: temperature against
  the 300 K equipartition requires, and the sign of the potential energy
  against the static lattice. It is not a drift measurement &mdash; reading a
  drift rate is how chromium was recorded as an integrator artefact while its
  crystal was in fact coming apart.</p>`;
}

function ugBlock(d){
  const u=d.ug, ms=(d.dyn||{}).stable, us=(u.dyn||{}).stable;
  const um=u.mech||{}, dm=d.mech||{};
  //  dg === "exact" prints the stored value untouched, for measured numbers
  //  whose precision is the source's and not ours to invent
  const cell=(v,dg)=>v===undefined||v===null||isNaN(v)?"&mdash;"
                     :(dg==="exact"?String(v)
                                   :(+v).toFixed(dg===undefined?1:dg));
  /*  dg formats the two calculated columns, dgExp the measured one, because
      they are not the same kind of number.  Giving the whole row "exact" so
      that the experiment kept its own precision printed the calculated entropy
      as 46.670359764216904 - a float's full width, in a table read by eye. */
  const row=(k,mv,uv,ev,dg,dgExp)=>`<tr><td>${k}</td>
      <td class="mono">${cell(mv,dg)}</td>
      <td class="mono">${cell(uv,dg)}</td>
      <td class="mono" style="color:var(--ink-3)">${cell(ev,
        dgExp===undefined?dg:dgExp)}</td></tr>`;
  const keys=d.struct==="hcp"?["C11","C12","C13","C33","C44"]:["C11","C12","C44"];
  const flag=s=>s===undefined?"&mdash;":
    (s?`<span style="color:var(--good)">stable</span>`
      :`<span style="color:var(--bad)">unstable</span>`);
  const par=(k,mv,uv,dg)=>`<tr><td>${k}</td>
      <td class="mono">${cell(mv,dg)}</td><td class="mono">${cell(uv,dg)}</td>
      <td></td></tr>`;
  const uph=u.exp_phonon, dph=d.exp_phonon;
  const gain=d.rms-u.rms;
  return `
  <h3>MAU against UG</h3>
  <p class="plotnote">${u.comparable===false
    ? `The two fits used different three-body cutoffs, so UG is kept out of the
       plots above: the difference between the curves would be the cutoff as
       much as the angular term.`
    : `Both potentials are already drawn above &mdash; the directional response,
       the thermodynamics and the mechanical table show them side by side, and
       the dispersion offers <em>MAU and UG together</em> in its view selector.`}
  </p>
  <div class="grid2">
    <div>
      <table><thead><tr><th>elastic (GPa)</th>
        <th><span class="tag mau">MAU</span></th>
        <th><span class="tag ug">UG</span></th>
        <th>exp</th></tr></thead><tbody>
        ${keys.map(k=>row(k, d[k]?d[k][0]:null, u[k]?u[k][0]:null,
                          d[k]?d[k][1]:null)).join("")}
        <tr><td><strong>RMS</strong></td>
          <td class="mono"><strong>${d.rms.toFixed(2)}%</strong></td>
          <td class="mono"><strong>${u.rms.toFixed(2)}%</strong></td>
          <td></td></tr>
        ${row("B", d.B?d.B[0]:null, u.B?u.B[0]:null, d.B?d.B[1]:null)}
        ${row("E<sub>coh</sub> (eV)", d.Ecoh?d.Ecoh[0]:null,
              u.Ecoh?u.Ecoh[0]:null, d.Ecoh?d.Ecoh[1]:null, 3)}
        <tr><td>dynamically</td><td>${flag(ms)}</td><td>${flag(us)}</td>
          <td></td></tr>
        <tr><td>lowest &omega; (cm<sup>&minus;1</sup>)</td>
          <td class="mono">${cell((d.dyn||{}).most_neg_cm1,1)}</td>
          <td class="mono">${cell((u.dyn||{}).min_mesh_cm1,1)}</td>
          <td></td></tr>
      </tbody></table>

      <table><thead><tr><th>parameters</th>
        <th><span class="tag mau">MAU</span></th>
        <th><span class="tag ug">UG</span></th><th></th></tr></thead><tbody>
        ${par("m",d.m,u.m,3)}${par("&gamma;",d.gamma,u.gamma,3)}
        ${par("D (eV)",d.D,u.D,4)}${par("&alpha; (1/&Aring;)",d.alpha,u.alpha,4)}
        ${par("r<sub>0</sub> (&Aring;)",d.r0,u.r0,4)}${par("C",d.C,u.C,4)}
        ${par("s<sub>3</sub>",d.s3,u.s3,4)}
        <tr><td>&lambda;<sub>2</sub></td>
          <td class="mono" style="color:var(--ink-3)">0</td>
          <td class="mono">${cell(u.lam2,3)}</td><td></td></tr>
        <tr><td>&lambda;<sub>4</sub></td>
          <td class="mono" style="color:var(--ink-3)">0</td>
          <td class="mono">${cell(u.lam4,3)}</td><td></td></tr>
      </tbody></table>
    </div>
    <div>
      <table><thead><tr><th>mechanical (Hill)</th>
        <th><span class="tag mau">MAU</span></th>
        <th><span class="tag ug">UG</span></th><th></th></tr></thead><tbody>
        ${par("B<sub>H</sub> (GPa)",dm.B_H,um.B_H)}
        ${par("G<sub>H</sub> (GPa)",dm.G_H,um.G_H)}
        ${par("E<sub>H</sub> (GPa)",dm.E_H,um.E_H)}
        ${par("&nu;<sub>H</sub>",dm.nu_H,um.nu_H,3)}
        ${par("Pugh B/G",dm.pugh,um.pugh,3)}
        ${par("A<sup>U</sup>",dm.A_U,um.A_U,3)}
        ${par("Debye &Theta; (K)",dm.debye,um.debye,0)}
        ${par("v<sub>m</sub> (m/s)",dm.v_m,um.v_m,0)}
      </tbody></table>

      <table><thead><tr><th>phonons</th>
        <th><span class="tag mau">MAU</span></th>
        <th><span class="tag ug">UG</span></th><th>exp</th></tr></thead><tbody>
        ${row("&omega;<sub>max</sub> (cm<sup>&minus;1</sup>)",
              (d.ld||{}).maxfreq, (u.ld||{}).maxfreq, null)}
        ${row("S(298) J/mol&middot;K", tvn(d.ld,"S"), tvn(u.ld,"S"),
              d.S298, 2, "exact")}
        ${row("C<sub>v</sub>(298) J/mol&middot;K", tvn(d.ld,"Cv"),
              tvn(u.ld,"Cv"), d.Cp298, 2, "exact")}
        ${(dph||uph)?`<tr><td>measured &omega; error</td>
          <td class="mono">${dph?dph.mae.toFixed(1)+"%":"&mdash;"}</td>
          <td class="mono">${uph?uph.mae.toFixed(1)+"%":"&mdash;"}</td>
          <td></td></tr>`:""}
      </tbody></table>

      <p class="plotnote">The angular factor
      ${gain>=0?`removes ${gain.toFixed(2)} percentage points of elastic error
        here`:`costs ${(-gain).toFixed(2)} percentage points here`}.
      ${us===false?`<strong>But this UG fit is dynamically unstable</strong>
        (lowest frequency ${u.dyn.min_mesh_cm1.toFixed(0)} cm<sup>&minus;1</sup>),
        so the elastic agreement is not usable on its own &mdash; matching every
        C<sub>ij</sub> says nothing about the phonons.`
       :`The fit is dynamically stable on the 8&sup3;&nbsp;&cup;&nbsp;9&sup3;
        union mesh${(u.dyn||{}).n_screened>1?`, and it is the best-scoring of
        ${u.dyn.n_screened} independent searches that is`:""}.`}
      ${uph?` The measured frequencies are out of sample for both columns: the
        objective sees only the q&nbsp;&rarr;&nbsp;0 limit.`:""}</p>
      ${u.comparable===false?`<div class="warn">The two fits used
      <strong>different three-body cutoffs</strong>
      (${u.rcut3_mau.toFixed(2)} vs ${u.rcut3.toFixed(2)} &Aring;), so this
      comparison measures the cutoff as much as the angular term.</div>`:""}
    </div>
  </div>`;
}

/*  The heat-capacity decomposition at 298.15 K, interpolated on the thermo
    grid rather than read at the nearest point, so the numbers on the page are
    the ones the notes quote.  cv: the model; lat: 9 alpha^2 B V T with the
    measured alpha and B; el: gamma*T from refdata_electronic (Kittel Table 2),
    null where the source has no gamma; rest: what is left of the tabulated
    C_p.  */
const ELSRC = "__ELSRC__";
function cvAt(g, T){
  if(!g||!g.thermo) return null;
  const p=[...g.thermo].sort((a,b)=>a.T-b.T);
  for(let i=1;i<p.length;i++) if(p[i-1].T<=T&&T<=p[i].T)
    return p[i-1].Cv+(p[i].Cv-p[i-1].Cv)*(T-p[i-1].T)/(p[i].T-p[i-1].T);
  return null;
}
function cpParts(d){
  const cv=cvAt(d.ld,298.15), lat=cpLat(d);
  const el=(d.gamma_e!=null)?d.gamma_e*298.15/1000:null;
  if(cv==null||!d.Cp298) return null;
  return {cv, lat, el,
          rest:(lat!=null&&el!=null)?d.Cp298-cv-lat-el:null};
}
/*  Library-wide, over every element where all three terms are known -
    counted from DATA, never typed.  */
let CP_LIB=null;
function cpLibrary(){
  if(CP_LIB) return CP_LIB;
  const med=a=>{const s=[...a].sort((x,y)=>x-y), n=s.length;
    return n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2;};
  const rows=Object.values(DATA).map(d=>({d, p:cpParts(d)}))
    .filter(o=>o.p&&o.p.lat!=null&&o.p.el!=null);
  const pct=f=>med(rows.map(o=>100*(f(o.p)-o.d.Cp298)/o.d.Cp298));
  const byRest=[...rows].sort((a,b)=>a.p.rest-b.p.rest);
  CP_LIB={n:rows.length,
    m0:pct(p=>p.cv), m1:pct(p=>p.cv+p.lat), m2:pct(p=>p.cv+p.lat+p.el),
    //  overshoot is only claimed where gamma is large (>= 5 mJ/mol K^2):
    //  beryllium also comes out negative, but because the model's own C_v is
    //  too high, not because of the electronic term
    over:byRest.filter(o=>o.p.rest<-0.75&&o.d.gamma_e>=5).map(o=>o.d.name),
    left:byRest.slice(-4).reverse().map(o=>o.d)};
  return CP_LIB;
}
const nameList=a=>a.length<2?(a[0]||""):a.slice(0,-1).join(", ")+" and "+a[a.length-1];
const ALKALI=["Li","Na","K","Rb","Cs"];

/*  value at the tabulated temperature nearest 298 K, as a number */
function tvn(g,key){
  if(!g||!g.thermo) return null;
  let best=g.thermo[0];
  for(const r of g.thermo) if(Math.abs(r.T-298)<Math.abs(best.T-298)) best=r;
  return best[key];
}

/*  One elastic number out of whichever arm, whatever shape it stores.  The
    root record and UG keep [fit, experiment] pairs; the tapered and re-cut
    arms keep a Cij dict and a scalar B, because they were produced by a
    different driver.  Reading the wrong shape returns undefined and prints an
    em-dash, which looks like "not measured" rather than like a bug - so the
    two shapes are handled in one place instead of at each call.  */
function elVal(r, k){
  if(!r) return null;
  if(k==="Ecoh") return Array.isArray(r.Ecoh)?r.Ecoh[0]
                      :(typeof r.Ecoh==="number"?r.Ecoh:null);
  if(k==="B") return Array.isArray(r.B)?r.B[0]
                   :(typeof r.B==="number"?r.B:null);
  if(r.Cij && typeof r.Cij[k]==="number") return r.Cij[k];
  return Array.isArray(r[k])?r[k][0]:null;
}
/*  Every arm that carries elastic numbers, minus the one already in the "fit"
    column.  The candidates are the point: they were re-cut to reach the same
    targets from a different cutoff, so seeing their C11 beside the published
    one is the most direct statement of what the re-cut did.  */
function elArms(d){
  return [["MAU",d],["UG",d.ug],["tap",d.tap],["tap UG",d.tap_ug],
          ["re-cut",d.rc],["re-cut UG",d.rc_ug]]
    .filter(o=>o[1] && (o[1].Cij || Array.isArray(o[1].C11)))
    .filter(o=>o[1] !== parSet(d))
    .map(o=>[o[0]+(armBad(o[1])?" ✗":""), o[1]]);
}

function erow(label,pair,frozen,mpval,unc,key,arms){
  if(!pair) return "";
  const calc=Math.abs(pair[0]), exp=Math.abs(pair[1]);
  const err=exp?100*(calc-exp)/exp:null;
  const w=err===null?0:Math.min(Math.abs(err),40)*1.4;
  const fz=(frozen===undefined||frozen===null)?"&mdash;":Math.abs(frozen).toFixed(1);
  const mv=(mpval===undefined||mpval===null)?"&mdash;":Math.abs(mpval).toFixed(1);
  const extra=(arms||[]).map(a=>{
    const v=elVal(a[1],key);
    return `<td class="mono" style="color:var(--ink-2)">${
      v===null||v===undefined?"&mdash;":Math.abs(v).toFixed(1)}</td>`;}).join("");
  return `<tr><td>${label}</td>
    <td class="mono"><strong>${calc.toFixed(1)}</strong></td>${extra}
    <td class="mono" style="color:var(--ink-3)">${fz}</td>
    <td class="mono" style="color:var(--phi3)">${mv}</td>
    <td class="mono">${exp.toFixed(unc?(unc>=1?1:unc>=0.1?2:3):1)}${unc?`<span
      style="color:var(--ink-3)"> &plusmn;${
        unc.toFixed(unc>=1?1:unc>=0.1?2:3)}</span>`:""}</td>
    <td class="mono err" style="color:${
      err!==null&&Math.abs(err)>15?"var(--bad)":"var(--ink-2)"}">${
      err===null?"&mdash;":(err>0?"+":"")+err.toFixed(1)+"%"
    }${unc&&err!==null?`<span style="color:var(--ink-3)"> (${
      (Math.abs(calc-exp)/unc).toFixed(1)}&sigma;)</span>`:""
    }<i class="bar" style="width:${w}px"></i></td></tr>`;
}

/*  The parameter block for one potential.

    `u` is the UG record when the UG set is selected, null for MAU.  The two
    differ only by the Legendre factor, so everything below is shared and the
    angular part is added rather than duplicated - a second copy of this
    function is how the UG columns on this page came to be missing the
    thermodynamics and the mechanics in the first place.

    A UG record carries no crystal: the structure, lattice constant and mass are
    properties of the element and live once, on the MAU record.  Both exports
    therefore read the crystal from `d` and the potential from `p`. */
function paramText(d,el,fmtSel,u){
  const p = u || d;
  const ug = !!u;
  const P=[["m",p.m],["gamma",p.gamma],["D (eV)",p.D],["alpha (1/A)",p.alpha],
           ["r0 (A)",p.r0],["C",p.C],["alpha3 (1/A)",p.alpha3],
           ["s3 = alpha3/alpha",p.s3],["rcut2 (A)",p.rcut2],
           ["rcut3 (A)",p.rcut3]];
  if(ug) P.push(["lambda2",p.lam2],["lambda4",p.lam4]);
  if(fmtSel==="json"){
    const o={element:el,potential:ug?"UG":"MAU",structure:d.struct,a0:d.a0};
    if(d.c_over_a) o.c_over_a=d.c_over_a;
    ["m","gamma","D","alpha","r0","C","alpha3","s3","rcut2","rcut3"]
      .forEach(k=>o[k]=p[k]);
    if(ug){o.lam2=p.lam2; o.lam4=p.lam4;}
    return JSON.stringify(o,null,2);
  }
  if(fmtSel==="python"){
    return `# ${el} (${d.struct}) - ${ug?"UG":"MAU"}\n`+
      `from latdyn import Crystal, Potential\n`+
      `cry = Crystal(${JSON.stringify(d.struct)}, ${d.a0}`+
      (d.c_over_a?`, ${d.c_over_a}`:", None")+`, mass=${d.mass})\n`+
      `pot = Potential(m=${p.m}, D=${p.D}, alpha=${p.alpha}, r0=${p.r0},\n`+
      `                gamma=${p.gamma}, C=${p.C}, alpha3=${p.alpha3},\n`+
      `                rcut2=${p.rcut2}, rcut3=${p.rcut3}`+
      (ug?`,\n                lam2=${p.lam2}, lam4=${p.lam4})`
         :`)`)+
      (ug?`\n# lam2/lam4 need the angular tree: angular/latdyn.py, which\n`+
          `# delegates the force constants to angfc.py when they are non-zero.`
        :"");
  }
  let s=`Ugur interatomic potential library - ${el} (${d.struct}) - `
       +`${ug?"UG (Ugur-Guler, angular)":"MAU (modified Akgun-Ugur)"}\n`
       +`phi2 + phi3 form of I. Akgun and G. Ugur, Phys. Rev. B 51, 3458 (1995);`
       +` Nuovo Cimento D 20, 1549 (1998).\nRefitted parameters, not the `
       +`published ones.\n`;
  s+=`  phi2(r)     = D/(m-1) (r0/r)^g [ e^{m a (r0-r)} - m e^{a (r0-r)} ]\n`;
  s+=`  phi3(r1,r2) = C D/(m-1) (r0/x)^g [ e^{m a3 (r0-x)} - m e^{a3 (r0-x)} ]`
     +(ug?` * f(t)\n`:`\n`);
  s+=`              x = r1 + r2,  both legs inside rcut3\n`;
  if(ug){
    s+=`  f(t)        = 1 + lam2 P2(cos t) + lam4 P4(cos t),  t the apex angle\n`;
    s+=`              P2 and P4 average to zero over the sphere, so f only\n`;
    s+=`              redistributes the three-body energy between geometries.\n`;
    s+=`              lam2 = lam4 = 0 gives MAU back exactly.\n`;
  }
  s+="\n";
  P.forEach(([k,v])=>{s+=`  ${k.padEnd(20)} ${Number(v).toPrecision(10)}\n`;});
  s+=`\n  lattice constant     ${d.a0} A`+(d.c_over_a?`,  c/a ${d.c_over_a}`:"")+"\n";
  s+=`  triplets per atom    ${p.ntrip!==undefined?p.ntrip:d.ntrip}\n`;
  s+=`  elastic RMS          ${(ug?u.rms:d.rms).toFixed(2)} %\n`;
  return s;
}

function render(){
  const d=DATA[cur], t=tier(d.rms), g=d.ld||{}, fz=d.frozen||{};
  const mp=d.mp||{}, me=mp.elastic||{}, mpph=mp.phonon;
  const fcbad=(d.fc_check||[]).length>0;
  const P=parSet(d), isUG=(P!==d);
  /*  Mechanics, dynamics and reach all belong to the set on screen.  They
      were read from the root record whatever was selected, so choosing UG or
      a re-cut candidate left the elastic anisotropy, the Poisson range and
      the stability verdict describing MAU while the parameters above them
      described something else - the same "verdict measured from something
      other than what is shown" defect the table cell already fixed once.  */
  const dyn=(P.dyn||d.dyn||{}), mech=(P.mech||d.mech), rch=(P.reach||d.reach);
  /*  The comparison column is the OTHER set: whichever arm is not on screen
      and has mechanics of its own.  It used to be UG unconditionally, which
      put UG beside UG.  */
  const other=(()=>{
    const alt=setsOf(d).filter(o=>o[0]!==PAR_SET
                                  && o[2] && o[2].mech
                                  && (o[0]!=="ug" || d.ug.comparable));
    return alt.length?alt[0]:null;})();
  const umech=other?other[2].mech:null;
  const uname=other?other[1]:"UG";
  /*  Which set is on screen, by key, so the panel can name it instead of
      saying "UG" for anything that is not the root record.  */
  const SETNAME={mau:"MAU",ug:"UG",rc:"re-cut",rc_ug:"re-cut UG"};
  const selName=SETNAME[PAR_SET]||"MAU";
  const isAng = P.lam2!==undefined || P.lam4!==undefined;
  /*  The stability warning has to describe the arm whose parameters are on
      screen.  It read d.dyn whatever was selected, so choosing UG showed the
      MAU record's verdict and the angular arm never got a warning of its own.
      Aluminium and sodium are the cases: unstable on the symmetry path in UG,
      stable in MAU.  */
  const dynSel=(isUG && P && P.dyn) ? P.dyn : (d.dyn||{});
  /*  The bound flags belong to the parameter set actually on screen.  Reading
      them from the root record put MAU's flags under UG's numbers, and hid
      UG's own: eight hard-truncated and seven switched angular records sit at
      the gamma bound and were shown unmarked. */
  const bounds=(P.at_bound||[]);
  $("#panel").innerHTML=`
  <div class="phead">
    <h2 class="disp">${cur}</h2>
    <span class="meta">${d.name} &middot; ${d.struct} &middot; a<sub>0</sub> =
      ${d.a0} &Aring;${d.c_over_a?` &middot; c/a = ${d.c_over_a}`:""}
      &middot; d<sub>nn</sub> = ${fmt(d.dnn,3)} &Aring;</span>
    <span class="rms ${tier(P.rms)}" title="root-mean-square residual of the
      fit against the elastic constants it was given - three numbers for a
      cubic element, five for a hexagonal one. It says nothing about anything
      else on this page.">elastic RMS ${P.rms.toFixed(0)}%</span>
  </div>

  <h3>Fitted parameters${setsOf(d).length>1?` <select id="parset" class="inlinesel">
      ${setsOf(d)
        .filter(o=>o[2]).map(o=>`<option value="${o[0]}"${
          PAR_SET===o[0]?" selected":""}>${o[1]} &mdash; ${
          o[2].rms.toFixed(2)}%${armStab(o[2], o[0])}</option>`)
        .join("")}
    </select>`:""}</h3>
  <div class="pars">
    <div class="par p2"><span class="k">m</span>
      <span class="v mono">${fmt(P.m,3)}</span></div>
    <div class="par p2 ${bounds.includes("gamma")?"flag":""}"${bounds.includes("gamma")?' title="at a search bound, not at an optimum"':""}>
      <span class="k">&gamma;</span>
      <span class="v mono">${fmt(P.gamma,4)}</span></div>
    <div class="par p2"><span class="k">D (eV)</span>
      <span class="v mono">${fmt(P.D)}</span></div>
    <div class="par p2"><span class="k">&alpha; (&Aring;<sup>-1</sup>)</span>
      <span class="v mono">${fmt(P.alpha)}</span></div>
    <div class="par p2"><span class="k">r<sub>0</sub> (&Aring;)</span>
      <span class="v mono">${fmt(P.r0)}</span></div>
    <div class="par p3"><span class="k">C</span>
      <span class="v mono">${fmt(P.C,3)}</span></div>
    <div class="par p3"><span class="k">&alpha;<sub>3</sub> (&Aring;<sup>-1</sup>)</span>
      <span class="v mono">${fmt(P.alpha3)}</span></div>
    <div class="par p3"><span class="k">s<sub>3</sub></span>
      <span class="v mono">${fmt(P.s3,3)}</span></div>
    ${isAng?`<div class="par p3"><span class="k">&lambda;<sub>2</sub></span>
      <span class="v mono">${fmt(P.lam2,3)}</span></div>
    <div class="par p3"><span class="k">&lambda;<sub>4</sub></span>
      <span class="v mono">${fmt(P.lam4,3)}</span></div>`:""}
    <div class="par p3"><span class="k">triplets / atom</span>
      <span class="v mono">${P.ntrip||d.ntrip||"&mdash;"}</span></div>
  </div>
  ${P.rms<0.005?`<p class="note"><strong>An elastic RMS of 0.00 % is real and
  it is narrow.</strong> The residual here is ${P.rms.toExponential(1)} %:
  the fit reproduced the ${d.struct==="hcp"?"five":"three"} elastic constants
  it was handed, to machine precision. ${
    [["MAU",d],["UG",d.ug],["re-cut",d.rc],["re-cut UG",d.rc_ug]]
      .filter(o=>o[1]&&typeof o[1].rms==="number"&&o[1].rms>=0.005)
      .map(o=>`${o[0]} reaches ${o[1].rms.toFixed(1)} %`).join(", ")
    || "Every set for this element does the same"}. It does not mean the
  potential is right about anything it was not fitted to &mdash; the phonons,
  surfaces, stacking faults and thermal expansion further down are all
  predictions, and several of them are not close.</p>`:""}
  ${isUG?`<p class="note">These are the <strong>${selName}</strong>
  parameters &mdash; a separate fit, not MAU with something added. Every one of
  them differs: for this element D is ${fmt(d.D,3)} eV under MAU against
  ${fmt(P.D,3)} here.${isAng?` Setting &lambda;<sub>2</sub> =
  &lambda;<sub>4</sub> = 0 recovers the MAU <em>form</em> exactly, but not
  these values.`:""}${PAR_SET==="rc"&&RCREC(d)?` This is the shell-gap re-cut,
  and for ${d.name} it is the <strong>recommended set for molecular
  dynamics</strong>: it passes every cold and warm screen the rule asks for.`
  :PAR_SET.startsWith("rc")?` The shell-gap re-cut is <strong>not the
  recommended set</strong> ${PAR_SET==="rc_ug"?"in its angular form":"for this element"}:
  it passes every cold screen, and whether it can be used warm is a
  per-element decision.`:""}</p>`:""}
  ${dynSel.stable===false?`<div class="warn">
    <strong>Dynamically unstable.</strong> ${(dynSel.imag_frac*100).toFixed(1)}% of
    modes on the ${String(dynSel.nq).split("+").map(n=>n+"&sup3;").join(" and ")}
    meshes${dynSel.path?" and the standard path":""} are imaginary, down to
    ${(dynSel.most_neg_cm1/CM1_PER_THZ).toFixed(2)} THz${
      dynSel.min_path_cm1!==undefined&&dynSel.min_path_cm1<dynSel.min_mesh_cm1
      ?" &mdash; and the mesh alone does not find it, which is why the "
       +"dispersion above dips below zero where a mesh check said it did not"
      :""}. The fit only ever sees elastic
    constants, which are the q&rarr;0 limit, so nothing in it prevents this:
    Born stability and dynamical stability are separate criteria. Usable for
    elasticity, not for lattice dynamics or molecular dynamics.</div>`:""}
  ${bounds.length?`<div class="warn">Parameter${bounds.length>1?"s":""}
    <strong>${bounds.join(", ")}</strong> ended on the edge of the allowed range,
    so this is a constrained optimum rather than a converged one. The bound on
    m comes from the range of &alpha;/&beta; reported for this potential family
    (about 2 to 11); values far outside it reproduce the elastic constants with
    an implausibly hard core.</div>`:""}
  ${fcbad?`<div class="warn">Force-constant check failed:
    ${d.fc_check[0]}</div>`:""}

  <div class="grid2" style="margin-top:22px">
    <div>
      <h3>Potential</h3>
      <canvas id="plot"></canvas>
      <p class="plotnote"><i class="swatch" style="background:var(--phi2)"></i>
      &phi;<sub>2</sub>(r) &nbsp;&middot;&nbsp;
      <i class="swatch" style="background:var(--phi3)"></i>
      &phi;<sub>3</sub>(r,r) for a symmetric triplet, x = 2r.
      Dotted vertical: the first neighbour distance.
      ${(d.ug&&d.ug.comparable)?`Dashed: the <span class="tag ug">UG</span>
        pair term, a separate fit with its own <i>D</i>, &alpha; and
        <i>r</i><sub>0</sub>. The shaded band is UG's three-body term:
        unlike MAU's it is <strong>not a single curve</strong>, because the
        Legendre factor makes it depend on the apex angle as well as on
        <i>x</i>, so at each <i>x</i> it spans an interval &mdash; here
        h &isin; [${hRange(d.ug)[0].toFixed(2)},
        ${hRange(d.ug)[1].toFixed(2)}] over the full range of &theta;.
        Both potentials share these axes because they are the same quantity in
        the same units; drawing them apart would hide exactly the difference
        worth seeing.`:""}</p>
      ${coreless(d)?`<div class="warn"
        style="border-left-color:var(--bad);margin-top:9px">
        <strong>No repulsive core.</strong> &phi;<sub>2</sub> has no minimum:
        it falls without limit as r &rarr; 0, so nothing stops two atoms
        merging. This is the only fit in the library where that happens, and
        it is why the panel above looks unlike the others &mdash; there is no
        wall to draw. It comes from <i>m</i> = ${d.m.toFixed(2)} sitting on
        its lower bound, where the repulsive exponential decays almost as
        slowly as the attractive one, together with the
        (r<sub>0</sub>/r)<sup>&gamma;</sup> prefactor at
        &gamma; = ${d.gamma.toFixed(2)}. The elastic constants and the phonons
        are unaffected, since both are evaluated near d<sub>nn</sub>; the
        potential is not usable for molecular dynamics.</div>`:""}
    </div>
    <div>
      <h3>Elastic constants (GPa)</h3>
      ${(()=>{const arms=elArms(d);
        /*  "fit" is the set selected above, and every other arm that has
            elastic numbers gets its own column beside it.  The error column
            belongs to the fit column only - the extras are there to be
            compared with it, not each scored separately.  */
        const epair=k=>{const root=d[k];
          if(!Array.isArray(root)) return null;
          const v=elVal(P,k);
          return [(v===null||v===undefined)?root[0]:v, root[1]];};
        return `<table><thead><tr><th></th><th>fit (${selName})</th>${
          arms.map(a=>`<th style="font-weight:400">${a[0]}</th>`).join("")
        }<th>frozen ion</th>
        <th>MP (DFT)</th><th>exp</th><th>error</th></tr></thead><tbody>
        ${erow("E_coh (eV)",epair("Ecoh"),null,null,d.Ecoh_unc,"Ecoh",arms)}
        ${erow("B",epair("B"),null,me.K,null,"B",arms)}
        ${erow("C11",epair("C11"),fz.C11,me.C11,(d.Cij_unc||{}).C11,"C11",arms)}
        ${erow("C12",epair("C12"),fz.C12,me.C12,(d.Cij_unc||{}).C12,"C12",arms)}
        ${d.struct==="hcp"?erow("C13",epair("C13"),fz.C13,me.C13,null,"C13",arms):""}
        ${d.struct==="hcp"?erow("C33",epair("C33"),fz.C33,me.C33,null,"C33",arms):""}
        ${erow("C44",epair("C44"),fz.C44,me.C44,(d.Cij_unc||{}).C44,"C44",arms)}
      </tbody></table>${arms.some(a=>armBad(a[1]))?`<p class="plotnote">
        &#10007; marks an arm with imaginary modes on the mesh or the symmetry
        path: it reaches these elastic constants and is still not a harmonic
        minimum. For a re-cut candidate that is why it failed selection; a
        published arm marked so ships with the same warning in its potential
        file.</p>`:""}`;})()}
      ${mp.mp_id?`<p class="plotnote" style="margin-top:9px">
        Materials Project <strong>${mp.mp_id}</strong>, space group
        ${mp.spacegroup}${mp.matches_structure?"":
          " &mdash; <strong>different structure from ours</strong>"},
        DFT d<sub>nn</sub> = ${mp.a_dft?mp.a_dft.toFixed(4):"?"} &Aring;
        (ours ${d.dnn.toFixed(4)},
        ${((mp.a_dft/d.dnn-1)*100).toFixed(1)}%).</p>`:""}
      ${mp.elastic_rejected?`<div class="warn"><strong>Materials Project's
        elastic tensor for this element is not usable</strong> and is left out
        of the column above: ${mp.elastic_rejected}. It is shown as absent
        rather than plotted, because a broken reference on the page is worse
        than a missing one. The phonon comparison, where MP has one, is a
        separate calculation and is unaffected.</div>`:""}
      ${(me.tensor_source==="ieee_format")?`<p class="plotnote">MP's unrounded
        tensor is not in the standard orientation for this element, so the
        column above uses their rounded one. That costs nothing here, but for a
        very soft metal integer rounding can erase C&prime; entirely.</p>`:""}
      ${mp.e_above_hull>0.005?`<div class="warn">This structure is
        <strong>${(mp.e_above_hull*1000).toFixed(0)} meV/atom above the hull</strong>
        in Materials Project's own 0 K calculation, so their DFT prefers a
        different polymorph. We fit the room-temperature phase, which is the
        right target here, but it does mean their elastic constants and phonons
        describe a structure their functional does not favour &mdash; and a
        structure that is not a minimum can legitimately return a negative C44.
        Treat the MP column for this element with care.</div>`:""}
      ${d.Ecoh_unc?`<p class="plotnote">The cohesive energy, the lattice
      constant and the bulk modulus are <strong>hard constraints</strong>,
      satisfied exactly at every trial point; only the elastic constants are
      scored. Brewer quotes E<sub>coh</sub> for ${d.name} to
      &plusmn;${(100*d.Ecoh_unc/Math.abs(d.Ecoh[1])).toFixed(1)}%, so pinning
      it exactly is more confident than the measurement is &mdash; worth
      remembering wherever the RMS is very small.</p>`:""}
      <p class="plotnote">The relaxed column carries the non-affine
      internal-strain correction. For one atom per primitive cell symmetry makes
      it vanish, so relaxed and frozen agree; for hcp they do not, and only the
      relaxed values are comparable with experiment. Residual pressure
      ${(d.P_resid!==undefined?d.P_resid.toExponential(1):"?")} GPa.</p>
      <div class="warn" style="border-left-color:var(--bad);margin-top:9px">
        <strong>Not for molecular dynamics as published.</strong>
        &phi;<sub>2</sub> is truncated hard and does not vanish at the cutoff:
        for ${d.name} it is
        <strong>${d.md?d.md.step_eV.toFixed(4):"?"} eV</strong> there, so a
        neighbour crossing the sphere changes the energy by that much in one
        step. At a fixed geometry that is consistent, which is why every number
        on this page is sound; in dynamics energy is simply not conserved.
        ${(d.md&&d.md.drift_hard!==undefined)?`Measured: an NVE run at 600 K
          drifts by <strong>${d.md.drift_hard.toFixed(0)} meV/atom/ps</strong>
          with the cutoff truncated, against
          <strong>${Math.abs(d.md.drift_taper).toFixed(2)}</strong> with it
          switched off smoothly &mdash; a factor of
          ${(Math.abs(d.md.drift_hard/d.md.drift_taper)).toFixed(0)}.`:
          `Measured on eight elements spanning the range: the drift runs from
           5 to 7256 meV/atom/ps truncated, and 0.03 to 0.28 switched.`}
        The drift tracks the step rather than anything else &mdash;
        log&ndash;log correlation 0.944 over those eight, palladium stepping by
        0.0002 eV and drifting by 5, chromium by 0.129 and drifting by 7256.
        <br><br>Shifting is not a way out: &phi;<sub>2</sub> carries real
        binding out to r<sub>cut</sub>, and subtracting its value there moves
        the cohesive energy by more than a tenth of itself for 32 of the 38
        elements &mdash; for rhodium by 23.7 eV against a cohesive energy of
        5.75. The fix is to switch the cutoff off over a window and refit; see
        the reachability section for what that costs and what it buys.
        <br><br>The window is fixed at <strong>0.85</strong> of the cutoff for
        every element, and deliberately so. It was scanned &mdash; seven values,
        ten elements, 300 restarts each &mdash; and the response is not
        monotonic: freezing the parameters found at 0.85 and re-evaluating them
        elsewhere <em>without refitting</em> traces the same jagged curve, so it
        is the window and not the search. What moves is which neighbour shells
        fall inside the fade zone. Tuning it per element would fit a shell
        arrangement that only exists at the equilibrium volume at 0&nbsp;K: the
        boundary sits 2.21&nbsp;a<sub>0</sub> against a shell at
        &radic;5&nbsp;a<sub>0</sub>, a gap of <strong>1.17 %</strong> for every
        cubic element here, while thermal displacement at 600 K is several per
        cent of the nearest-neighbour distance. So <code>taper</code> stays a
        truncation constant alongside r<sub>cut2</sub> = 2.6&nbsp;a<sub>0</sub>,
        not an eighth fitted parameter. Nothing jumps when a shell crosses the
        boundary &mdash; S is C&sup2; &mdash; and energy conservation does not
        depend on the value.</div>
    </div>
  </div>

  ${d.ug?ugBlock(d):""}

  ${rch?`<h3>Can this form reproduce this metal?</h3>
  <div class="warn" style="border-left-color:var(--${rch.ok?"good":"bad"})">
    <strong>${rch.ok?"Within reach.":"Out of reach."}</strong>
    Fixing the bulk modulus pins C<sub>11</sub>+2C<sub>12</sub>, so only
    C&prime; = (C<sub>11</sub>&minus;C<sub>12</sub>)/2 and C<sub>44</sub> are
    free, and their ratio has a floor this functional form cannot go under.
    Measured R = C<sub>44</sub>/C&prime; = <strong>${rch.R_exp}</strong>;
    the lowest the &phi;<sub>2</sub>+&phi;<sub>3</sub> form reaches for
    ${d.name}, without the angular factor, is
    <strong>${rch.R_floor}</strong>${rch.from_fit?
      ` &mdash; a bound set by the fitted potential itself, which landed under
        where the search stopped`:""}${rch.ok?
      `, so the anisotropy is attainable and any remaining error is the search,
       not the form.`:
      `. Since ${rch.R_exp} &lt; ${rch.R_floor}, no parameter set with this
       cutoff reproduces the anisotropy however long the fit runs: reproducing
       C<sub>44</sub> forces C&prime; below its physical value, which softens
       the transverse [110] branch and is what the phonons show.`}
    ${(!rch.ok&&d.tap)?`<br><br><strong>&mdash; and the cutoff is why.</strong>
      Truncating &phi;<sub>2</sub> at r<sub>cut</sub> leaves a step there, and
      that step, not the functional form, is what closes the region off.
      ${rch.R_floor_taper!==undefined?`With the cutoff switched off smoothly the
        floor for ${d.name} is not ${rch.R_floor} but
        <strong>${rch.R_floor_taper}</strong>, below the measured
        ${rch.R_exp}.${rch.taper_floor_is_guard?` That number is where the
          search was told to stop rather than where the form does: it sits
          exactly on the shear guard, C<sub>44</sub> at one per cent of
          C<sub>11</sub>, which is why every fcc metal here reports the same
          0.023. The form reaches lower; how much lower is not measured, and a
          crystal below that guard is not a potential anyone could use.`:
          ` That one is measured, not a constraint: the minimum sits at
            C<sub>44</sub>/C<sub>11</sub> near 0.08, well clear of the shear
            guard.`} `:""}Switch
      the cutoff off smoothly over its outer 15 % and refit, and ${d.name}
      comes back at <strong>${d.tap.rms.toFixed(2)} %</strong>
      ${d.tap.R!==undefined&&d.tap.R!==null?`with R = ${d.tap.R.toFixed(2)}
        against a measured ${d.tap.R_exp.toFixed(2)}`:""},
      dynamically stable. All eight metals reported out of reach here do;
      four of them exactly.
      <br><br>It is the discontinuity and not the range: shortening the hard
      cutoff instead, 2.60 &rarr; 2.405 &rarr; 2.210 a<sub>0</sub> at the same
      search budget, moves niobium only 20.4 &rarr; 18.3 &rarr; 16.8 %, where
      the smooth switch reaches 0.00.
      <br><br>The switch is not a free improvement and is not what this library
      ships. Over all 38 elements the median error does not move,
      5.83 &rarr; 5.83 %, and that is a trade rather than no effect: the bcc
      metals gain, 11.45 to 5.41, while the hcp metals lose hard, 7.86 to
      19.81, and the alkalis lose too, 0.00 to 13.02. Cadmium and zinc, which fail on axial anisotropy rather
      than on C<sub>44</sub>/C&prime;, are not rescued either. What it removes
      is this one limitation.
      ${d.tap_ug?`<br><br><strong>And the angular factor is still worth
        something on top of it.</strong> That had to be checked, because
        &lambda;<sub>2</sub> and &lambda;<sub>4</sub> exist to reach exactly
        these eight metals and the switch reaches them without any angular
        term. Both arms switched at the same window, everything else identical:
        the switch alone takes four of the eight to an exact fit, and the
        angular term takes three more.
        ${d.tap_ug.lam_off?`For ${d.name} the search turned the weights
          <em>off</em> &mdash; &lambda; below 0.01 &mdash; so the switch is
          doing all the work here${(d.tap_ug.score*100)>1?`, and
          ${d.tap_ug.score*100 >= d.tap_ug.score_ctrl*100-0.5?
            `neither mechanism reaches this one: it stays at
             <strong>${(100*d.tap_ug.score).toFixed(2)} %</strong>`:
            `it still falls short at ${(100*d.tap_ug.score).toFixed(2)} %`}`:""}.`:
          `For ${d.name} the weights come out
           &lambda;<sub>2</sub> = ${d.tap_ug.lam2.toFixed(3)},
           &lambda;<sub>4</sub> = ${d.tap_ug.lam4.toFixed(3)} and the error
           goes ${(100*d.tap_ug.score_ctrl).toFixed(2)} &rarr;
           <strong>${(100*d.tap_ug.score).toFixed(2)} %</strong>.`}
        The weights are what make this a statement about the form rather than
        the search: every element that gained carries a large &lambda;, and
        every element that did not has &lambda; driven to zero by the search
        itself.`:""}`:""}
    ${rch.ug_below?` The angular factor is not bound by that number, and here it
      is not: the UG fit reaches R = <strong>${rch.R_ug}</strong>, below the
      floor of the form without it. h(cos&theta;) = 1 + &lambda;<sub>2</sub>
      P<sub>2</sub> + &lambda;<sub>4</sub>P<sub>4</sub> weights the triplets by
      bond angle, and C<sub>44</sub> and C&prime; sample different angles, so
      the two stop moving together.`:""}
    ${rch.gap!==null?` Distance from the measured (C&prime;/B, C<sub>44</sub>/B)
      to the reachable set: <strong>${rch.gap}</strong> of the target
      &mdash; the ratio test is necessary, this one is the sharper.`:""}
    ${rch.stale?`<br><br><em>Provisional.</em> This floor was measured with the
      three-body sum truncated at 1.12 d<sub>nn</sub>, and the library is now
      built at 1.50 &mdash; a different potential, carrying 153 triplets per
      atom in fcc instead of 66. The verdict above is therefore the old form's,
      not this one's, and is being re-measured. Two intervening passes are
      already discarded: the first left the search unconstrained, and the
      second, which screened for dynamical stability, still reached
      C<sub>44</sub> = 1.0&times;10<sup>&minus;6</sup> GPa for gold, because
      C<sub>44</sub> &rarr; 0 leaves the crystal marginal rather than unstable
      and every mesh frequency stays real.`:""}
  </div>`:""}

  ${mech?`<h3>Mechanical response of the elastic tensor</h3>
  <div class="grid2">
    <div>
      <table><thead>
        ${umech?`<tr><th></th><th colspan="3"><span class="tag mau">${selName
          }</span></th>
          <th colspan="3"><span class="tag ug">${uname}</span></th></tr>`:""}
        <tr><th></th><th>Voigt</th><th>Reuss</th><th>Hill</th>
        ${umech?`<th>Voigt</th><th>Reuss</th><th>Hill</th>`:""}</tr>
        </thead><tbody>
        <tr><td style="text-align:left">bulk B (GPa)</td>
          <td>${mech.B_V.toFixed(1)}</td><td>${mech.B_R.toFixed(1)}</td>
          <td><strong>${mech.B_H.toFixed(1)}</strong></td>
          ${umech?`<td>${umech.B_V.toFixed(1)}</td><td>${umech.B_R.toFixed(1)}</td>
            <td><strong>${umech.B_H.toFixed(1)}</strong></td>`:""}</tr>
        <tr><td style="text-align:left">shear G (GPa)</td>
          <td>${mech.G_V.toFixed(1)}</td><td>${mech.G_R.toFixed(1)}</td>
          <td><strong>${mech.G_H.toFixed(1)}</strong></td>
          ${umech?`<td>${umech.G_V.toFixed(1)}</td><td>${umech.G_R.toFixed(1)}</td>
            <td><strong>${umech.G_H.toFixed(1)}</strong></td>`:""}</tr>
      </tbody></table>
      <div class="cols3" style="margin-top:12px">
        ${cell3("Young E (Hill)",mech.E_H.toFixed(1)+" GPa",
                umech?umech.E_H.toFixed(1):null,uname)}
        ${cell3("Poisson (Hill)",mech.nu_H.toFixed(3),
                umech?umech.nu_H.toFixed(3):null,uname)}
        ${cell3("anisotropy A<sup>U</sup>",mech.A_U.toFixed(3),
                umech?umech.A_U.toFixed(3):null,uname)}
        ${cell3("Pugh B/G",mech.pugh.toFixed(2),
                umech?umech.pugh.toFixed(2):null,uname)}
        ${cell3("Cauchy C12-C44",mech.cauchy.toFixed(1)+" GPa",
                umech?umech.cauchy.toFixed(1):null,uname)}
        ${mech.debye?cell3("Debye temperature",mech.debye.toFixed(0)+" K",
                umech&&umech.debye?umech.debye.toFixed(0):null,uname):""}
        ${(d.theta_D_exp?cell3("&theta;<sub>D</sub> measured",
                d.theta_D_exp.toFixed(0)+" K","Stewart 1983, at 0 K"):"")}
        ${mech.v_l?cell3("v longitudinal",(mech.v_l/1000).toFixed(2)+" km/s",
                umech&&umech.v_l?(umech.v_l/1000).toFixed(2):null,uname):""}
        ${mech.v_t?cell3("v transverse",(mech.v_t/1000).toFixed(2)+" km/s",
                umech&&umech.v_t?(umech.v_t/1000).toFixed(2):null,uname):""}
      </div>
      <p class="plotnote">A<sup>U</sup> = 5G<sub>V</sub>/G<sub>R</sub> +
      B<sub>V</sub>/B<sub>R</sub> &minus; 6 vanishes only for an isotropic
      solid. B/G above ~1.75 is the usual ductility indicator, and a positive
      Cauchy pressure points the same way. The Debye temperature comes from the
      elastic constants and the density, not from the phonon spectrum &mdash;
      the two are independent estimates.
      ${mech.debye_iso?`<br><br><b>It is averaged over the real slowness
      surface.</b> For each direction the Christoffel problem gives three
      speeds, and the Debye model asks for
      3/v<sub>m</sub><sup>3</sup>&nbsp;=&nbsp;&lang;&sum;<sub>i</sub>
      1/v<sub>i</sub><sup>3</sup>&rang; over the sphere. Until 2026 this page
      built it from the Voigt&ndash;Reuss&ndash;Hill averages instead &mdash;
      one longitudinal and one transverse speed &mdash; which throws the
      anisotropy away <em>before</em> an average that the slow directions
      dominate, and came out high. For ${d.name} that old value was
      ${mech.debye_iso.toFixed(0)}&nbsp;K against
      ${mech.debye.toFixed(0)}&nbsp;K here; across the library the median move
      is &minus;1.5&nbsp;% and lithium's is &minus;8&nbsp;%.
      ${d.theta_D_exp?`<br><br>Do not expect the number beside it to match the
      measured &theta;<sub>D</sub>: that is a T&nbsp;&rarr;&nbsp;0 quantity and
      this one is built from elastic constants at ~300&nbsp;K and a
      room-temperature volume, which is the anchor mismatch described in the
      reference-data notes. Carried to 0&nbsp;K with III/29a's temperature
      coefficients the close-packed metals land on it to within half a per
      cent &mdash; copper 1.004, silver 0.999, gold 1.004, palladium 0.999
      &mdash; and that agreement is what says this averaging is the right
      one. The published comparison behind it is <b>A. Tari</b>,
      <i>The Specific Heat of Matter at Low Temperatures</i> (Imperial College
      Press, 2003), Table 2.4, which puts the calorimetric and elastic
      &theta;<sub>D</sub> side by side at T&nbsp;&rarr;&nbsp;0 for five metals:
      four agree to 0.55&nbsp;% or better and palladium is the outlier at
      1.7&nbsp;%.`:""}`:""}</p>
    </div>
    <div>
      <div class="gen" style="margin-bottom:8px">
        <select id="mprop">
          <option value="E">Young's modulus E</option>
          <option value="beta">linear compressibility</option>
          <option value="G">shear modulus G</option>
          <option value="nu">Poisson's ratio</option>
        </select>
        ${(()=>{const o=marmOpts(d).filter(x=>x[0]!==PAR_SET);
          return o.length>1?`<select id="marm" style="margin-left:6px">${
            o.map(x=>`<option value="${x[0]}">beside ${x[1]}</option>`)
             .join("")}</select>`:"";})()}
      </div>
      <canvas id="polar"></canvas>
      <p class="plotnote">Radius is the value in that direction, so a circle
      means isotropy in the plane.<i class="swatch"
      style="background:var(--phi3);margin-left:8px"></i>dashed: the upper
      envelope over the transverse direction, for the two quantities that need
      one. Dotted circle marks zero where the range crosses it.
      ${(d.mech_planes&&d.mech_planes.d)?`<br><b>The fourth panel is not a
      coordinate plane.</b> Between them x&ndash;y, x&ndash;z and y&ndash;z
      contain [100], [110] and [101] and <b>no member of
      &lang;111&rang;</b>&nbsp;&mdash; which for a cubic crystal is where
      Young's modulus is extremal, so three panels understate the stiffest
      direction badly and with nothing about them looking wrong. The
      (11&#772;0) plane is spanned by [110] and [001] and &lang;111&rang; sits
      inside it at 35.26&deg;.${d.struct==="hcp"?` Hexagonal elasticity is
      transversely isotropic, so for this element it adds no direction the
      x&ndash;z panel did not already have; it is drawn for consistency across
      the page.`:""}`:""}</p>
      <div class="cols3" style="margin-top:10px">
        ${cell3("E range",mech.E_min.toFixed(1)+" &ndash; "
           +mech.E_max.toFixed(1)+" GPa","&times;"+mech.E_aniso.toFixed(2))}
        ${cell3("Poisson range",mech.nu_min.toFixed(3)+" &ndash; "
           +mech.nu_max.toFixed(3))}
        ${cell3("shear range",mech.G_min.toFixed(1)+" &ndash; "
           +mech.G_max.toFixed(1)+" GPa")}
      </div>
      ${mech.nu_min<0?`<p class="plotnote"><strong>Auxetic directions.</strong>
        Poisson's ratio is negative down to ${mech.nu_min.toFixed(3)}, so along
        those the crystal expands sideways when stretched.</p>`:""}
      ${mech.born_stable?"":`<div class="warn">The stiffness matrix has a
        non-positive eigenvalue (${Math.min(...mech.eigenvalues).toFixed(1)}
        GPa), so this tensor is not Born stable and the averages above are not
        meaningful.</div>`}
    </div>
  </div>`:""}

  ${g.std?`<div class="grid2" style="margin-top:26px">
    <div>
      <h3>Phonon dispersion</h3>
      ${(mpph&&mpph.ours)||(d.mc3d&&d.mc3d.ours)
        ||(d.jarvis&&d.jarvis.trusted&&d.jarvis.ours)
        ||(d.ug&&d.ug.comparable)?`<div class="gen"
        style="margin-bottom:8px">
        <select id="dmode">
          ${mpph&&mpph.ours?`<option value="mp">vs Materials Project
            (${mpph.method})</option>`:""}
          ${d.mc3d&&d.mc3d.ours?`<option value="mc3d">vs Materials Cloud MC3D
            (PBEsol)</option>`:""}
          ${d.jarvis&&d.jarvis.trusted&&d.jarvis.ours?`<option value="jarvis">
            vs JARVIS-DFT (NIST)${d.jarvis.checked===false
              ? " – sole reference" : ""}</option>`:""}
          <option value="std">this potential only</option>
          ${d.ug&&d.ug.comparable?`<option value="ug">MAU and UG together
            (this potential only)</option>`:""}
          ${(d.rc&&d.rc.ld&&d.rc.ld.std)?`<option value="recut">MAU and the
            re-cut candidate together (this potential only)${
              armBad(d.rc)?" — REJECTED":""}</option>`:""}
          ${(d.rc_ug&&d.rc_ug.ld&&d.rc_ug.ld.std)?`<option value="recut_ug">MAU
            and the angular re-cut candidate together (this potential only)${
              armBad(d.rc_ug)?" — REJECTED":""}</option>`:""}
        </select>
      </div>`:""}
      <canvas id="disp"></canvas>
      <p class="plotnote"><strong>Which record this curve is.</strong> The
        parameter-set selector offers the HARD-TRUNCATED sets &mdash; its
        &ldquo;MAU&rdquo; is the root record shipped as
        <code>${d.name}.ugur</code> and its &ldquo;UG&rdquo; is
        <code>${d.name}.ugur.ang</code>. The switched-cutoff sets
        (<code>${d.name}_taper.ugur</code>) carry no stored dispersion and are
        reachable only from the download table below. That distinction is not
        cosmetic: measured against neutron and x-ray data over the 33 elements
        that have it, the hard-truncated sets average 9.7&nbsp;% and the switched
        ones 12.6&nbsp;%, and the hard one is closer in 25 of the 33. The switched
        sets are what molecular dynamics has to use, because a hard cut leaves
        the pair energy discontinuous, and that is the arm any
        finite-temperature number quoted elsewhere belongs to.</p>
      ${d.exp_phonon?`<p class="plotnote"><span style="display:inline-block;
        width:9px;height:9px;border:1.4px solid currentColor;border-radius:50%;
        margin-right:5px"></span><strong>Open circles: measured
        frequencies</strong> at
        ${d.exp_phonon.points.map(p=>p.name).join(" and ")}, mean absolute
        error <strong>${d.exp_phonon.mae}%</strong>. Nothing at finite q enters
        the fit, so these are out of sample.<br>Source: ${d.exp_phonon.ref}.
        </p>`:""}
      ${d.model_curve?`<p class="plotnote"><strong>The dashed lines are a
        model, not a measurement.</strong> ${d.model_curve.why}. What the
        source tabulates is the force constants of a ${d.model_curve.kind}
        model${d.model_curve.shells?` &mdash; ${d.model_curve.shells} shells,
        reaching ${d.model_curve.neighbours} neighbours &mdash;`:", "} and
        this is those constants evaluated along the path at
        ${d.model_curve.T_K} K, which
        is why it is drawn as a line and not as the circles used everywhere
        else here for individually measured frequencies. Reconstructing a
        tensor from a published index can go wrong quietly, giving a
        plausible curve rather than an error, so it was checked against
        ${d.model_curve.gate}.${d.model_curve.check?` For the record the
        model returns c<sub>11</sub> ${d.model_curve.check.c11},
        c<sub>44</sub> ${d.model_curve.check.c44} and c&prime;
        ${d.model_curve.check.cp} GPa, against ${d.model_curve.check.src}
        values of ${d.model_curve.check.published_c11},
        ${d.model_curve.check.published_c44} and
        ${d.model_curve.check.published_cp}.`:""}<br>Source:
        ${d.model_curve.ref}.</p>`:""}
      ${d.exp_curve?`<p class="plotnote">${d.exp_curve.points_only
        ?(nMeasured(d.exp_curve)===1
          ?`The circle is <strong>one measured frequency, not a
        dispersion</strong>, drawn at each place the path passes its point of
        the zone. Nothing between is measured and nothing is drawn`
          :`The circles are <strong>${nMeasured(d.exp_curve)} measured frequencies, not a dispersion</strong>.
        This element's was published as a figure with a list of points beside
        it &mdash; the symmetry points, or a handful the authors singled out
        &mdash; so there is no measured curve between them and none is drawn`)
        :d.exp_curve.digitised?`The diamonds along the branches are
        <strong>${Object.values(d.exp_curve.segs).reduce(
        (a,v)=>a+v.length,0)} points read from the published figure</strong>,
        not a table &mdash; see below`
        :`The circles along the branches are the <strong>whole measured
        dispersion</strong>, ${Object.values(d.exp_curve.segs).reduce(
        (a,v)=>a+v.length,0)} points with their quoted errors`}, at
        <strong>${d.exp_curve.T_K} K</strong>. The temperature is
        part of the datum rather than a footnote: frequencies soften measurably
        with it${d.exp_curve.T_K<200?`, and this set is well below room
        temperature &mdash; the full dispersion of this element has not been
        measured at 300 K, which is why the density-functional literature
        compares against ${d.exp_curve.T_K} K as well`:""}.${d.exp_curve.points_only?""
        :d.struct==="hcp"
        ?` The [&zeta;&zeta;0] branch runs past K and on to M; the points
        beyond K belong on M&ndash;K and are drawn there, not folded back
        onto &Gamma;&ndash;K.`
        :d.struct!=="fcc"?""
        :d.exp_curve.segs["U|X"]
        ?` Points on the part of &Sigma; beyond K are drawn on U&ndash;X, the
        same line continued to X by symmetry.`
        :` Points on the part of &Sigma; beyond K, where measured, are not
        drawn: this path does not go there.`}${d.exp_curve.a_meas?`<br><b>Scored at the
        crystal's lattice constant at the measurement temperature</b>,
        a&nbsp;=&nbsp;${d.exp_curve.a_meas.toFixed(4)}&nbsp;&#8491;${
          d.exp_curve.c_meas?`, c&nbsp;=&nbsp;${d.exp_curve.c_meas.toFixed(4)}`
          :""}, not at the fitted a<sub>0</sub> of
        ${d.a0.toFixed(4)}&nbsp;&#8491;. The fit targets a room-temperature
        lattice (a 5&nbsp;K one for rubidium and caesium, 78&nbsp;K for
        lithium) and this curve was taken at ${d.exp_curve.T_K}&nbsp;K, so
        comparing at a<sub>0</sub> would charge the potential for a bookkeeping
        mismatch rather than for physics.${d.exp_curve.a_meas_from?`
        ${d.exp_curve.a_meas_from==="expansion integral"
          ? `This volume is <b>not measured</b>: no second lattice constant is
             published for this crystal at this temperature, so it is carried
             from a<sub>0</sub> by the thermal expansion integral &mdash; the
             Debye heat capacity under the CRC's &alpha;(25&nbsp;&deg;C) with
             the Debye temperature above. Checked against sodium and potassium,
             where Kittel does publish both ends, it runs about 12&nbsp;% high,
             so read it as carrying that much slack.`
          : `This one is a measurement: ${d.exp_curve.a_meas_from}.`}`:""}`:""}
        ${d.exp_curve.digitised?`<br><b>Digitised, not tabulated.</b>
        ${d.exp_curve.digitised}`:""}
        <br>Source: ${d.exp_curve.ref}.</p>`:""}
      <p class="plotnote">${g.std[0].branches.length} branches. From the
      dynamical matrix &mdash; phonons are not fit targets, so this is a
      prediction. The dashed curve is whichever reference the selector shows,
      with our dispersion re-evaluated at that reference's own q-points, so the
      residual is physics rather than interpolation.
      <b>Each reference brings its own symmetry path</b>, so the horizontal
      axis and the labels on it change with the selector. The two curves in
      any one panel are always on the same path as each other; the paths
      across panels are not.
      ${(()=>{const nm=t=>t==="G"?"&Gamma;":t;
        const walk=r=>r&&r.marks?r.marks.map(m=>nm(m[1])).join("&ndash;"):"";
        return [["Materials Project",mpph],["MC3D",d.mc3d],
                ["JARVIS",d.jarvis&&d.jarvis.trusted?d.jarvis:null]]
          .filter(r=>walk(r[1])).map(r=>`<br>${r[0]}: ${walk(r[1])}`).join("");
      })()}
      <br>A reference need not sample the point a label names, so the labels
      sit at fractional positions between samples, and a curve can cross a
      &Gamma; line without having reached zero AT any sample there.
      ${PAR_SET!=="mau"?`<br><b>"This potential only" follows the parameter
      set selected above; the DFT comparisons do not.</b> Their residual was
      computed by re-evaluating our dispersion at the reference's own
      q-points, and that was done for MAU alone &mdash; no such resampling
      exists for the other arms, so those three views stay MAU whatever is
      selected.`:""}
      ${[["Materials Project ("+(mpph&&mpph.method||"")+")",mpph],
         ["Materials Cloud MC3D (PBEsol)",d.mc3d],
         ["JARVIS-DFT (NIST)", d.jarvis&&d.jarvis.trusted?d.jarvis:null]]
        .filter(r=>r[1]&&r[1].stats).map(r=>`<br>
        <i class="swatch" style="background:var(--phi3)"></i>vs <strong>${r[0]}
        </strong>: RMS ${(r[1].stats.rms_cm1/CM1_PER_THZ).toFixed(2)} THz
        (${r[1].stats.rel_pct!==undefined?r[1].stats.rel_pct
          :(100*r[1].stats.rms_cm1/r[1].stats.ref_max).toFixed(1)}% of their
        highest branch), top frequency
        ${(r[1].stats.ours_max/CM1_PER_THZ).toFixed(2)} vs
        ${((r[1].stats.mp_max||r[1].stats.ref_max)/CM1_PER_THZ).toFixed(2)}
        THz.`).join("")}</p>
    </div>
    <div>
      <h3>Phonon thermodynamics</h3>
      <canvas id="thermo"></canvas>
      <p class="plotnote"><i class="swatch" style="background:var(--phi2)"></i>
      S(T) &nbsp;&middot;&nbsp;<i class="swatch"
      style="background:var(--phi3)"></i>C<sub>v</sub>(T), J/(mol&middot;K).
      Dashed: the Dulong-Petit limit 3R.
      ${(d.S298||d.Cp298)?`Rings: the measured S&deg; and C<sub>p</sub>&deg; at
        298.15 K, from the <strong>CRC Handbook of Chemistry and Physics</strong>
        standard thermodynamic tables. Nothing thermal enters the fit, so these
        are out of sample. <strong>The heat-capacity ring is C<sub>p</sub> and
        the curve is C<sub>v</sub></strong>, which are not the same quantity
        &mdash; see the note under the numbers below.`
       :`<strong>No rings:</strong> the measured S&deg; and
        C<sub>p</sub>&deg; are not in the table for this element yet, so there
        is nothing to compare the curves against here.`}</p>
    </div>
  </div>`:""}


  ${g.thermo?`<h4 class="sub">Read off those curves at 298 K</h4>
  <div class="cols3">
    ${cell3("zero point energy",tv(g,"zpe")+" eV")}
    ${cell3("entropy S",tv(g,"S")+" J/(mol K)",d.S298?exact(d.S298):null)}
    ${cell3("heat capacity C_v",tv(g,"Cv")+" J/(mol K)",
       d.Cp298?exact(d.Cp298)+" (C_p)"+(cpLat(d)!=null
         ?"; "+(d.Cp298-cpLat(d)).toFixed(2)+" as C_v":""):null)}
    ${cell3("Helmholtz F",tv(g,"F")+" eV")}
    ${cell3("highest frequency",
       ((g.maxfreq||0)/CM1_PER_THZ).toFixed(2)+" THz ("
       +(g.maxfreq||0).toFixed(0)+" cm-1)")}
  </div>
  <p class="note"><strong>These are lattice-dynamics numbers, not molecular
  dynamics.</strong> They are quasi-harmonic sums over the phonon spectrum
  above, evaluated at 298 K from force constants computed at 0 K, so nothing
  in them depends on a trajectory and nothing in them is affected by the
  finite-temperature behaviour reported further down the page.
  <br><br><strong>The tabulated value is C<sub>p</sub>; the calculation gives
  C<sub>v</sub>.</strong> They differ by the lattice term
  C<sub>p</sub>&nbsp;&minus;&nbsp;C<sub>v</sub>&nbsp;=&nbsp;9&alpha;<sup>2</sup>BV<sub>m</sub>T,
  with &alpha; the linear expansion coefficient.${cpLat(d)!=null?` With this
  element's measured &alpha; and bulk modulus it is
  <strong>${cpLat(d).toFixed(2)}&nbsp;J/(mol&nbsp;K)</strong> at 298&nbsp;K,
  ${(100*cpLat(d)/d.Cp298).toFixed(1)}&nbsp;% of C<sub>p</sub>, and the
  measured value carried to C<sub>v</sub> is
  ${(d.Cp298-cpLat(d)).toFixed(2)}.`:""}
  <br><br><strong>The conduction electrons add a second term, &gamma;T</strong>,
  which no interatomic potential contains.${(()=>{const p=cpParts(d);
    return (p&&p.el!=null)?` For ${d.name} &gamma; =
    ${d.gamma_e}&nbsp;mJ&nbsp;mol<sup>&minus;1</sup>&nbsp;K<sup>&minus;2</sup>
    gives <strong>${p.el.toFixed(2)}&nbsp;J/(mol&nbsp;K)</strong> at
    298&nbsp;K${p.rest!=null?`, and what is left between the curve and the
    ring is <strong>${p.rest.toFixed(2)}&nbsp;J/(mol&nbsp;K)</strong>`:""}.`
    :` The source used has no &gamma; for ${d.name}, so that term is not shown
    and the gap above is not decomposed further.`;})()}
  ${(()=>{const M=cpLibrary(); return `Over the ${M.n} elements where both
    terms are known, the calculated C<sub>v</sub> sits a median
    ${M.m0.toFixed(1)}&nbsp;% from the tabulated C<sub>p</sub>; the lattice term
    takes that to ${M.m1.toFixed(1)}&nbsp;% and the electronic term to
    <strong>${M.m2.toFixed(1)}&nbsp;%</strong>. Two things remain, and both are
    physics rather than the potential. &gamma; is measured at low temperature,
    where it carries the electron&ndash;phonon mass enhancement; that fades
    above the Debye temperature, so &gamma;T at 298&nbsp;K is an upper bound, and
    it overshoots most where &gamma; is large (${nameList(M.over)}). The largest
    remainders are ${nameList(M.left.map(x=>x.name))}.${(()=>{
      const alk=M.left.filter(x=>ALKALI.indexOf(x.sym)>=0).map(x=>x.name),
            rest=M.left.filter(x=>ALKALI.indexOf(x.sym)<0).map(x=>x.name);
      return (alk.length?` The alkali metals among them are hot for their
        melting point at 298&nbsp;K, where anharmonicity and thermal vacancies
        add heat capacity that a harmonic calculation from 0&nbsp;K force
        constants does not have.`:"")
        +(rest.length?` That does not explain ${nameList(rest)}, for which
        298&nbsp;K is far from melting.`:"");})()}`;})()}
  <br><span style="opacity:.75">&gamma;: ${ELSRC}.</span> No thermal quantity
  enters the fit.</p>`:""}

  ${finiteTBlock(d)}

  ${d.elasticT?`<h3>Elastic constants against temperature
    <select id="etq" class="inlinesel">${ET_QS.map(a=>
      `<option value="${a[0]}"${a[0]===ET_Q?" selected":""}>${a[1]}</option>`
      ).join("")}</select></h3>
  <canvas id="elasT"></canvas>
  <p class="plotnote">
    <strong>Every curve is named at its right-hand end</strong>, in its own
    colour &mdash; there are ${Object.keys(d.elasticT).length} of them here and
    no swatch list would carry that. Ours are drawn heavier; the
    published potentials are the thin dotted ones.
    ${Object.values(d.elasticT).some(r=>r.kind==="base")?`
      Those are
      ${Object.values(d.elasticT).filter(r=>r.kind==="base")
         .map(r=>r.label+" &mdash; <code>"+r.file+"</code>").join("; ")},
      shipped with LAMMPS and run through the identical cell, recipe
      and post-processing &mdash; without them these curves would have nothing
      to be measured against.`:""}
    The dashed vertical is the melting point.
    <strong>Molecular dynamics</strong>, unlike everything above it on this
    page: a thermostatted trajectory at each temperature, with the tensor from
    the Born stress-fluctuation method &mdash; the recipe of LAMMPS's own
    <code>examples/ELASTIC_T/BORN_MATRIX</code>.
    ${Object.values(d.elasticT).some(r=>r.pts.some(q=>q.above_melt||!q.born_ok))?`
      <strong>Hollow markers</strong> sit above the melting point or violate the
      Born criteria: a small perfect crystal has nowhere to nucleate from and
      superheats, so those points describe a metastable solid, not the metal.`:""}
    ${Object.values(d.elasticT).some(r=>r.nudge_bad)?`
      <strong style="color:var(--bad)">Dashed red</strong>: this lattice does not
      survive a 1e-5 &Aring; displacement, and thermal motion here is ten
      thousand times that &mdash; the curve is a property of whatever structure
      it relaxes into, not of the one it is labelled with.`:""}
    ${(d.aflow&&d.aflow.usable&&(ET_Q==="B"||ET_Q==="G"))?`
      <strong>Ring at T = 0:</strong> AFLOW's density-functional value
      (B = ${d.aflow.B.toFixed(1)}, G = ${(d.aflow.G||0).toFixed(1)} GPa,
      space group ${d.aflow.spacegroup}). One 0 K point, not a curve, and only
      B and G exist there.`:""}
    ${(d.aflow&&!d.aflow.usable)?`
      <strong>AFLOW has a value for this element but in the wrong phase</strong>
      (${d.aflow.struct||"?"}, space group ${d.aflow.spacegroup}, against this
      library's ${d.struct}), so it is not shown &mdash; an elastic modulus
      belongs to a structure.`:""}
    ${d.elasticT[Object.keys(d.elasticT)[0]].grid==="fixed"?`
      Grid: fixed 0&ndash;1200 K, as in the ruthenium study this reproduces.`:`
      Grid: fractions of this element's own melting point, so the curve spans
      the same physical range for every element.`}</p>`:""}

  ${(d.tap&&d.tap.ground)?`<h3>Is the fitted structure the ground state?</h3>
  <div class="cols3">
    ${["bcc","fcc","hcp"].map(s=>cell3(s+(s===d.tap.ground.want?" (fitted)":""),
       (d.tap.ground.rel[s]>=0?"+":"")+d.tap.ground.rel[s].toFixed(1)+" meV/atom",
       d.tap_ug&&d.tap_ug.ground?
         ((d.tap_ug.ground.rel[s]>=0?"+":"")
          +d.tap_ug.ground.rel[s].toFixed(1)+" (UG)"):null)).join("")}
  </div>
  <p class="note">Each structure built, relaxed under this potential with all
  three cell axes free, and compared per atom at its own relaxed lattice. The
  cell shape is checked afterwards, because a cell allowed to relax freely can
  leave the symmetry it was built in, and then the number in the column would
  not belong to the structure named at the top of it.
  ${d.tap.ground.ok?`<strong>${d.name}'s fitted ${d.tap.ground.want} is the
    lowest of the three</strong> &mdash; four of this library's seventy-six
    records manage that.`:`<strong style="color:var(--bad)">The fitted
    ${d.tap.ground.want} is not the lowest: ${d.tap.ground.lowest} lies
    ${Math.abs(d.tap.ground.rel[d.tap.ground.lowest]).toFixed(1)} meV/atom
    below it.</strong> Seventy-two of this library's seventy-six records are in
    the same position, against one of the nineteen published potentials tested
    the same way.
    ${(d.tap_ug&&d.tap_ug.ground&&d.tap_ug.ground.ok)?`<strong>The angular term
      fixes it:</strong> under UG the fitted ${d.tap.ground.want} is the lowest
      of the three. Four records in the library reach that and all four are
      UG, which is the clearest thing measured so far in the angular term's
      favour.`:""}
    It does not invalidate the elastic tensor &mdash; a
    metastable structure has perfectly well-defined curvatures, and the phonons
    confirm the reference is a local minimum &mdash; but it does mean the
    crystal has somewhere to go, and defect, surface, melting and
    high-temperature properties are not safe here.`}
  ${(d.baseline_ground&&Object.keys(d.baseline_ground).length)?`
    Published potentials for ${d.name} run through the identical test:
    ${Object.entries(d.baseline_ground).map(([f,r])=>
      "<code>"+f+"</code> &rarr; "+r.lowest
      +(r.ok?"":" <strong>(also wrong)</strong>")).join("; ")}.`:""}</p>

  <canvas id="bain"></canvas>
  <p class="plotnote">
    ${armKey(armRecs(d,r=>r.bain&&r.bain.E&&r.bain.E.length))}
    &mdash; energy along the volume-conserving tetragonal strain
    (1+&delta;, 1+&delta;, 1/(1+&delta;)<sup>2</sup>), whose curvature at
    &delta; = 0 is C&prime; = (C<sub>11</sub>&minus;C<sub>12</sub>)/2.
    That curvature is a target of the fit and comes out right; the shape away
    from the origin is not, and is what this shows. The ring marks where the
    well stops being a well.
    ${(d.tap.bain&&d.tap.bain.turn_up!==null&&d.tap.bain.turn_up!==undefined)?`
      For ${d.name} that happens at &delta; = ${d.tap.bain.turn_up.toFixed(3)},
      over a hump of only
      <strong>${(d.tap.bain.barrier_up*1000).toFixed(1)} meV/atom</strong>.
      Across the seven cubic metals measured, the height of that hump orders the
      finite-temperature failures exactly &mdash; Fe 0.7 meV fails at 0.05
      T<sub>m</sub>, V 0.9 at 0.05, Mo 1.0 at 0.10, W 1.3 at 0.08, Nb 1.8 at
      0.60, Ta 3.2 never. A static scan that takes seconds anticipates a
      molecular-dynamics sweep that takes a night on forty cores.`:`
      For ${d.name} the curve rises on both sides throughout the scanned range,
      which is what a sound well looks like: forty-six of the seventy-six
      records reach that.`}</p>`:""}

  ${(d.tap&&d.tap.surface)?`<h3>Surface energy</h3>
  <div style="overflow-x:auto">
  <table><thead><tr><th>facet</th>
    ${armRecs(d,r=>r.surface&&r.surface.gamma).map(o=>
      `<th>${o[2]}</th>`).join("")}
    <th>DFT</th>${d.surface_ref&&d.surface_ref.tyson?"<th>experiment</th>":""}
    ${d.baseline_surface&&Object.keys(d.baseline_surface).length
      ?`<th>published${Object.keys(d.baseline_surface).length>1
        ?" (lowest&ndash;highest)":""}</th>`:""}</tr></thead><tbody>
  ${d.tap.surface.order_want.map(f=>{
    const cells=armRecs(d,r=>r.surface&&r.surface.gamma)
      .map(o=>o[1].surface.gamma[f]);
    const dft=(d.surface_ref&&d.surface_ref.facets)?d.surface_ref.facets[f]:null;
    const bs=Object.values(d.baseline_surface||{})
      .map(s=>s.gamma[f]).filter(v=>v!==undefined&&v!==null);
    return `<tr><td>(${f})</td>
      ${cells.map(v=>`<td>${v!==undefined&&v!==null
        ?v.toFixed(3):"&mdash;"}</td>`).join("")}
      <td>${dft?dft.toFixed(3):"&mdash;"}</td>
      ${d.surface_ref&&d.surface_ref.tyson
        ?`<td>${d.surface_ref.tyson.toFixed(2)}</td>`:""}
      ${bs.length?`<td>${bs.length===1
         ? bs[0].toFixed(3)
         : Math.min(...bs).toFixed(3) + "&ndash;" + Math.max(...bs).toFixed(3)
           + ' <span class="dim">(' + bs.length + ')</span>'}</td>`:
        (d.baseline_surface&&Object.keys(d.baseline_surface).length
         ?"<td>&mdash;</td>":"")}</tr>`;}).join("")}
  </tbody></table></div>
  <p class="note">J/m&sup2;. A slab with vacuum along the normal, the atoms
  relaxed and the cell fixed, against a separately relaxed bulk so that no
  residual pressure is counted as surface energy.
  ${(()=>{const b=d.baseline_surface||{};const k=Object.keys(b);
    if(!k.length) return "";
    return k.length===1
      ? ` The published column is <code>${k[0]}</code>.`
      : ` The published column spans ${k.length} potentials, lowest to
          highest: ${k.map(x=>"<code>"+x+"</code>").join(", ")}; they are
          separate potentials rather than a range of uncertainty.`;})()}
  ${d.surface_ref?`DFT from the Materials Project surface database
    (${d.surface_ref.mp_id}); experiment from Tyson and Miller, Surf. Sci.
    <strong>62</strong>, 267 (1977), one number per element with no facet
    resolution. GGA surface energies run 10&ndash;30 % below experiment, which
    is why both are quoted &mdash; they bracket the answer from either side.`:""}
  ${(d.baseline_surface&&Object.keys(d.baseline_surface).length)?`
    The published potentials went through the identical slab, vacuum, bulk
    reference and relaxation, which is what makes the column beside them a
    measurement of the potential rather than of the recipe.
    <em>They are not all tight:</em> across the library they sit at a median
    0.88 of the density-functional value with the middle half between 0.76 and
    1.01 and the full span 0.48 to 1.97, so a published potential is typically
    a little below DFT and occasionally far from it. That spread is the scale
    this element&rsquo;s numbers should be read against.
    <strong>That is checked from outside, not asserted:</strong> the NIST
    Interatomic Potentials Repository runs its own calculations on the
    potentials it hosts, and over the 74 facets belonging to potentials
    confirmed to be the same file, the median difference from the numbers here
    is 0.033 % and every one of the 74 agrees to within 1 %.
    ${(()=>{const all=Object.values(d.baseline_surface||{});
      const n=all.filter(s=>s.nist&&s.nist.worst_pct!==undefined);
      const bad=all.filter(s=>s.nist&&s.nist.unreliable);
      let out="";
      if(n.length){const w=Math.max(...n.map(s=>s.nist.worst_pct));
        out += ` For ${d.name}, ${n.length} of these potential${
          n.length>1?"s were":" was"} checked that way and agree${
          n.length>1?"":"s"} to within ${
          w<0.1?"0.1":w.toFixed(1)} % on every facet.`;}
      /*  Saying "the largest disagreement is 2828 %" on ruthenium's page read
          as this calculation failing, when what failed was the record it was
          being compared against: NIST gives the same 0.1036 J/m2 for both the
          fcc (111) and the hcp (0001), which is not physics, while its own
          double-hcp entry gives 2.8376 next to the 3.03 here and an
          experimental 3.05.  A number is not worth quoting if its frame
          inverts what it means.                                             */
      if(bad.length) out += ` ${bad.length===all.length?"The":"One"} NIST
        record${bad.length>1?"s":""} for ${d.name} ${bad.length>1?"are":"is"}
        internally inconsistent &mdash; the close-packed surface energy differs
        between crystal structures that differ only in stacking &mdash; so
        ${bad.length>1?"they are":"it is"} not used as a check here.`;
      return out;})()}
    `:""}</p>

  <div class="cols3">
    ${cell3("ours, against DFT",
       d.tap.surface.ratio_dft_median
         ? "&times;" + d.tap.surface.ratio_dft_median.toFixed(2) : "&mdash;",
       d.tap.surface.ratio_exp_median
         ? "&times;" + d.tap.surface.ratio_exp_median.toFixed(2)
           + " against experiment" : null)}
    ${cell3("facet ordering",
       d.tap.surface.order.map(f=>"("+f+")").join(" &lt; "),
       "should be " + d.tap.surface.order_want.map(f=>"("+f+")").join(" < ")
       + (d.tap.surface.order_basis==="rule" ? " (rule)" : ""))}
    ${cell3("anisotropy",
       (100*d.tap.surface.spread).toFixed(1) + " %",
       (d.surface_ref&&d.surface_ref.anisotropy!==undefined
         &&d.surface_ref.anisotropy!==null)
         ? "DFT " + (100*d.surface_ref.anisotropy).toFixed(1) + " %" : null)}
  </div>

  <p class="note">${d.tap.surface.order_ok===null?`The three reference facets
    fall within ${(100*d.tap.surface.order_margin).toFixed(1)} % of one another,
    which is finer than this potential &mdash; or the calculation the reference
    comes from &mdash; can resolve, so the ordering is not scored for this
    element.`:d.tap.surface.order_ok?`The ordering is right: the facets come out
    in the order ${d.tap.surface.order_basis==="rule"
    ?`close packing predicts. That is the target here only because the
    reference does not cover all three facets of this element; where it does,
    it is the reference and not the rule that is scored against`
    :`the reference calculation gives for this element`}.`:`<strong
    style="color:var(--bad)">The ordering is wrong.</strong> ${
    d.tap.surface.order_basis==="rule"
    ?`The reference does not cover all three facets of this element, so the
    target fell back to close packing &mdash; which is a weaker test than the
    rest of this column, because the rule itself is unreliable: the reference
    contradicts it for 17 of the 35 elements it does cover.`
    :`The reference calculation puts these facets in a different order. Close
    packing is the usual rule but not a reliable one &mdash; the reference
    contradicts it for 17 of the 35 elements it covers &mdash; so the target
    here is that element&rsquo;s own calculated ordering, and on the records
    where it is resolved most published potentials reproduce it; this record
    does not.`}`}
  It is the same shortage the vacancy energy exposes, and the usual account of
  it &mdash; that the energy of a bond in this form cannot depend on how many
  other bonds an atom has &mdash; is <em>too strong</em>. A three-body term
  counts neighbour pairs, so it does carry a coordination dependence, and a
  two- plus three-body potential with free radial shapes fitted to
  density-functional energies and forces reproduces coordination-sensitive
  energies to within a few per cent (Xie, Rupp and Hennig, npj Comput. Mater.
  <strong>9</strong>, 162 (2023)).
  So what these numbers measure is this parameterisation
  &mdash; parameters fitted to bulk elastic data, which never saw a surface
  &mdash; rather than the functional form.
  ${uf3Note(DATA.Nb, d.name==="Niobium")}</p>`:""}

  ${(d.tap&&d.tap.stacking)?`<h3>Stacking fault, and the prediction it tests</h3>
  ${(()=>{const hex=d.struct==="hcp";
    return `<p class="note">The ranking above says this metal&rsquo;s preferred
  structure is ${hex?"cubic":"hexagonal"} rather than the
  ${hex?"hexagonal":"face-centred"} one it was fitted to. An intrinsic stacking
  fault <em>is</em> a slab of ${hex?"cubic":"hexagonal"} stacking inside a
  ${hex?"hexagonal":"face-centred"} crystal, so that ranking carries a
  prediction with a sign attached: the fault should cost <strong>less</strong>
  than nothing. This measures it, by shifting half the crystal along the
  ${hex?"basal partial a/3[1&#772;100] on (0001)"
       :"Shockley partial a/6[11&#772;2&#772;] on (111)"} and tilting the cell
  to match, so exactly one fault exists rather than two.</p>`;})()}

  <div style="overflow-x:auto">
  <table><thead><tr><th></th>
    ${armRecs(d,r=>r.stacking).map(o=>`<th>${o[2]}</th>`).join("")}
    ${d.tap.stacking.exp?"<th>experiment</th>":""}
    ${Object.keys(d.baseline_stacking||{}).length
      ?`<th>published${Object.keys(d.baseline_stacking).length>1
        ?" (lowest&ndash;highest)":""}</th>`:""}</tr></thead><tbody>
  ${[["intrinsic fault &gamma;<sub>isf</sub>","isf"],
     ["unstable fault &gamma;<sub>usf</sub>","usf"]].map(row=>{
    const key=row[1];
    const cells=armRecs(d,r=>r.stacking).map(o=>o[1].stacking[key]);
    const bs=Object.values(d.baseline_stacking||{})
      .map(s=>s[key]).filter(v=>v!==undefined&&v!==null);
    return `<tr><td>${row[0]}</td>
      ${cells.map(v=>`<td${key==="isf"&&v<0
        ?' style="color:var(--bad)"':""}>${v!==undefined&&v!==null
        ?v.toFixed(1):"&mdash;"}</td>`).join("")}
      ${d.tap.stacking.exp
        ?`<td>${key==="isf"?d.tap.stacking.exp.toFixed(0):"&mdash;"}</td>`:""}
      ${Object.keys(d.baseline_stacking||{}).length
        ?`<td>${bs.length===0?"&mdash;":bs.length===1
           ? bs[0].toFixed(1)
           : Math.min(...bs).toFixed(1)+"&ndash;"+Math.max(...bs).toFixed(1)
             +' <span class="dim">('+bs.length+')</span>'}</td>`:""}
      </tr>`;}).join("")}
  </tbody></table></div>
  <p class="note">mJ/m&sup2;. The intrinsic fault is the crystal left behind
  after one partial has passed; the unstable fault is the barrier on the way
  to it, and the difference between them is what a fault must climb to heal.
  ${(()=>{const k=Object.keys(d.baseline_stacking||{});
    if(!k.length) return "";
    return ` The published column is ${k.length===1?"":"spanned by "}${
      k.map(x=>"<code>"+x+"</code>").join(", ")}, through the identical
      shift, tilt and relaxation.`;})()}</p>

  <canvas id="gamma"></canvas>
  <p class="plotnote">
    ${armKey(armRecs(d,r=>r.stacking&&r.stacking.gamma&&
                             r.stacking.gamma.length))}
    ${Object.keys(d.baseline_stacking||{}).length?`&nbsp;&middot;&nbsp;
      published potentials in thin grey`:""}
    &mdash; energy against the shift, over one whole period. Two points on it
    have names: the hump at 1/6 is the unstable fault and the value at 1/3
    (dashed) is the intrinsic one. The far side is not a fault at all but the
    arrangement close packing forbids, a layer resting directly on its
    neighbour, and the curve returning to zero at a full period is the check
    that the cell is built right.</p>

  ${(()=>{const g=d.tap.stacking;
    const bs=Object.values(d.baseline_stacking||{}).map(s=>s.isf)
      .filter(v=>v!==undefined&&v!==null);
    const pos=bs.filter(v=>v>0).length;
    const nist=Object.values(d.baseline_stacking||{})
      .filter(s=>s.nist_isf!==undefined&&s.nist_isf!==null);
    if(g.isf>=0) return `<p class="note">The fault costs energy here, so the
      prediction the ranking made is not borne out for ${d.name} &mdash; worth
      recording, because the same argument does hold for most of the library.
      </p>`;
    return `<p class="note"><strong style="color:var(--bad)">The intrinsic
    fault energy is negative.</strong> That is not a large error, it is the
    wrong sign${g.exp?`: the measured value for ${d.name} is
    +${g.exp.toFixed(0)} mJ/m&sup2;`:""}. A negative fault energy means the
    faulted crystal is <em>lower</em> than the perfect one, so a fault that
    forms never heals and the partial dislocations bounding it repel without
    limit. The perfect crystal is still metastable rather than unstable
    &mdash; the ${g.usf.toFixed(0)} mJ/m&sup2; hump has to be climbed first,
    and getting back out costs ${g.back_barrier
      ? g.back_barrier.toFixed(0)
      : (g.usf-g.isf).toFixed(0)} &mdash; so nothing collapses on its own.
    What it rules out is everything past the elastic regime: plasticity,
    dislocation motion, deformation, defect evolution. What it leaves standing
    is the fitted elastic tensor and the phonons, which are curvatures at a
    minimum and do not care that a deeper one exists elsewhere.
    ${bs.length?`<br><br><strong>The machinery is not what produced it.</strong>
      ${pos} of the ${bs.length} published potential${bs.length>1?"s":""} for
      ${d.name} went through the identical code and came back
      ${pos===bs.length?"positive":"mixed"}${
        pos===bs.length&&bs.length>1
          ? ", from "+Math.min(...bs).toFixed(1)+" to "
            +Math.max(...bs).toFixed(1)+" mJ/m&sup2;":""}.`:""}
    ${nist.length?` The same numbers are computed independently by NIST, by a
      different geometry, for ${nist.length} of those file${
        nist.length>1?"s":""}: ${nist.map(s=>s.isf.toFixed(2)+" here against "
        +s.nist_isf.toFixed(2)).join(", ")} mJ/m&sup2;.`:""}
    ${g.predicted?`<br><br><strong>It was predicted before it was
      measured.</strong> If a bond&rsquo;s energy really does not know how many
      neighbours an atom has, the fault should simply cost what its two
      ${d.struct==="hcp"?"cubic":"hexagonal"} layers cost in bulk &mdash;
      ${d.struct==="hcp"
        ?"2(E<sub>fcc</sub>&minus;E<sub>hcp</sub>)/A"
        :"2(E<sub>hcp</sub>&minus;E<sub>fcc</sub>)/A"},
      which for ${d.name} is ${g.predicted.toFixed(1)} against the
      ${g.isf.toFixed(1)} measured here. Those two numbers share no code: one
      relaxes bulk crystals, the other slides half a slab.${g.library
        ?` Across the library the sign agrees ${g.library.sign_ok} times out of
           ${g.library.n} and the ratio has a median of
           ${g.library.ratio_median.toFixed(2)}.`:""}
      The ranking and the fault are not two findings, they are one.`:""}
    <br><br>It is the same shortage as the vacancy, the surface and the ground
    state, seen from a fourth direction. The cause is a fit that was shown one
    bulk structure at one volume and never a sheared plane, rather than the
    functional form itself: the same two- plus three-body expansion, fitted to
    density-functional energies and forces, reproduces coordination-sensitive
    energies to within a few per cent (Xie, Rupp and Hennig, npj Comput. Mater.
    <strong>9</strong>, 162 (2023)).
    ${uf3Note(DATA.Nb, d.name==="Niobium")}</p>`;})()}
  `:""}

  ${(d.tap&&d.tap.expansion&&d.tap.expansion.failed)?`
  <h3>Thermal expansion</h3>
  <p class="note"><strong style="color:var(--bad)">No coefficient: the crystal
  does not survive a barostat.</strong> Under zero applied pressure the cell
  ran away and LAMMPS stopped with <code>${d.tap.expansion.failed}</code>
  ${d.tap.expansion.of?` (${d.tap.expansion.points} of
    ${d.tap.expansion.of} temperatures usable)`:""}. That is the dynamic face
  of the compression escape recorded for this element among the pathologies:
  the fitted lattice is a local minimum on a ledge, and nothing in the
  potential holds the cell once a barostat is free to move it. The absence is
  recorded here rather than left blank, because a missing section and a
  measurement that failed are not the same statement.</p>
  `:""}

  ${(d.tap&&d.tap.expansion&&d.tap.expansion.alpha_1e6!==undefined
     &&!d.tap.expansion.failed)?`
  <h3>Thermal expansion</h3>
  <p class="note">This is the one test on the page that stays <em>inside</em>
  the regime the form describes. The vacancy, the surface, the ground state and
  the stacking fault all ask what happens when an atom&rsquo;s coordination
  changes. Here the atoms keep every neighbour and simply sit further apart, so
  what is measured is the anharmonicity of the same bonds whose curvature was
  fitted &mdash; the third derivative of a curve whose second derivative is a
  target. A barostat at zero pressure, the average cell edge against
  temperature, no fitting and no mode tracking.</p>

  <div class="cols3">
    ${cell3("&alpha; near 300 K",
       d.tap.expansion.alpha_1e6.toFixed(1) + " &times;10<sup>&minus;6</sup>/K",
       d.tap.expansion.alpha_exp_1e6
         ? "experiment " + d.tap.expansion.alpha_exp_1e6.toFixed(1) : null)}
    ${cell3("against experiment",
       d.tap.expansion.ratio ? "&times;" + d.tap.expansion.ratio.toFixed(2)
                             : "&mdash;",
       (()=>{const b=Object.values(d.baseline_expansion||{})
               .map(x=>x.ratio).filter(v=>v);
         return b.length ? "published " + Math.min(...b).toFixed(2) + "&ndash;"
                           + Math.max(...b).toFixed(2) : null;})())}
    ${cell3("Gr&uuml;neisen &gamma;",
       d.tap.expansion.gruneisen!==undefined&&d.tap.expansion.gruneisen!==null
         ? d.tap.expansion.gruneisen.toFixed(2) : "&mdash;",
       "a metal sits between 1 and 3")}
  </div>

  <canvas id="expan"></canvas>
  <p class="plotnote">
    ${armKey(armRecs(d,r=>r.expansion&&r.expansion.T&&
                             r.expansion.T.length>1))}
    ${Object.keys(d.baseline_expansion||{}).length?`&nbsp;&middot;&nbsp;
      published potentials in thin grey`:""}
    ${d.tap.expansion.alpha_exp_1e6?`&nbsp;&middot;&nbsp;
      experiment as the dashed straight line`:""}
    &mdash; each curve is normalised to its own value at the lowest
    temperature, because what is being compared is a slope and not a lattice
    constant.
    ${d.tap.expansion.alpha_exp_1e6?` The experimental line is the linear
      coefficient at 25&nbsp;&deg;C, drawn from the same starting point as the
      computed curves, from David R. Lide, ed., <em>CRC Handbook of Chemistry
      and Physics</em>, Internet Version 2005, CRC Press, Boca Raton, FL, 2005.
      Read it as a slope rather than a precise value: the Handbook's figures
      move between editions, and for the alkalis the spread between
      compilations is several per cent &mdash; for the refractory metals about
      one. It is a different edition from the melting points used elsewhere
      here, which are the 97th (2016).`:""}
    ${d.tap.expansion.window?`<br><br><b>The comparison is converted before it
      is made.</b> The slope above is fitted over the whole temperature grid,
      ${d.tap.expansion.T?d.tap.expansion.T[0].toFixed(0)+"&ndash;"
        +d.tap.expansion.T[d.tap.expansion.T.length-1].toFixed(0)+"&nbsp;K":""},
      and the experimental figure is at 25&nbsp;&deg;C. Those are not the same
      quantity: &alpha; follows the heat capacity, so the raw ratio carries a
      factor &lang;C<sub>V</sub>&rang;<sub>grid</sub>&nbsp;/&nbsp;C<sub>V</sub>(298)
      that belongs to neither the potential nor the experiment. Here that
      factor is <b>${d.tap.expansion.window.toFixed(3)}</b>, so the raw ratio of
      ${d.tap.expansion.ratio_raw!==undefined
        ? "&times;"+d.tap.expansion.ratio_raw.toFixed(2) : "&mdash;"}
      becomes &times;${d.tap.expansion.ratio.toFixed(2)}. It is applied
      identically to every published potential drawn beside ours, and that is
      the reason to believe it: a real conversion cannot tell the two apart,
      and in fact it improves the published sets more than it improves ours.
      For most elements it is within a per cent or two and changes nothing;
      for beryllium, lithium and sodium, whose grids sit worst against their
      own Debye temperature, it is decisive.`:""}</p>

  ${(()=>{const a=d.tap.expansion.alpha_1e6, r=d.tap.expansion.ratio;
    const b=Object.values(d.baseline_expansion||{}).map(x=>x.ratio)
             .filter(v=>v);
    const scale = b.length
      ? ` For ${d.name}, ${b.length} published potential${b.length>1?"s go":
          " goes"} through the identical barostat and land${b.length>1?"":"s"}
          at ${Math.min(...b).toFixed(2)}&ndash;${Math.max(...b).toFixed(2)} of
          experiment; across the library their median is 0.96 against our
          0.68.` : "";
    if(a<0) return `<p class="note"><strong style="color:var(--bad)">This
      record contracts on heating.</strong> A negative expansion coefficient
      for a simple metal is not a small error, and it is not a subtlety of the
      barostat: nine of our 74 records come back negative and <em>none</em> of
      the 51 published potentials does. All nine are body-centred.${scale}
      Thermal expansion is the third derivative of the same energy curve whose
      second derivative was fitted, so nothing constrains it &mdash; getting
      the curvature right at one volume says nothing about how it changes with
      volume.</p>`;
    return `<p class="note">The library sits at a median 0.68 of the measured
      coefficient, systematically low.${scale} That is a milder failure than
      the coordination tests above and it has a different cause: the expansion
      is the third derivative of the same energy curve whose second derivative
      was fitted, and nothing in the fit constrains it.</p>`;})()}
  `:""}

  ${d.tap_nudge?`<h3>An alternative fit that holds its lattice${
      d.tap_nudge.withdrawn?" &mdash; withdrawn for this element":""}</h3>
  <p class="note">${d.name}'s shipped switched record does not survive a
  1e-5 &Aring; displacement. The search was rerun keeping every distinct
  solution it found, and the best one that does survive is carried here beside
  it rather than instead of it &mdash; the trade is the interesting part, and
  choosing silently would hide it.</p>
  ${d.tap_nudge.withdrawn?`<p class="note" style="color:var(--bad)">
    <strong>This candidate does not in fact hold its lattice, and is shown
    only so the record is complete.</strong> Re-measured here against
    ${d.tap_nudge.jiggle.n} independent displacement directions it keeps
    ${Math.min(...d.tap_nudge.jiggle.keep_all).toFixed(0)}&ndash;${
      Math.max(...d.tap_nudge.jiggle.keep_all).toFixed(0)} times the
    displacement in every one of them. The passing verdict originally recorded
    for it belonged to a different solution: the filter read its LAMMPS log
    whether or not the run had started, so a failed invocation returned the
    previous candidate's result, and the job had been resubmitted six times
    with the pools regenerated each time. One element in five was affected.
    The corrected filter is being rerun.</p>`:""}
  <div class="cols3">
    ${cell3("shipped fit, error",d.tap_nudge.rms_best.toFixed(2)+" %",
            "does not hold its lattice")}
    ${cell3("this fit, error",d.tap_nudge.rms.toFixed(2)+" %",
            d.tap_nudge.jiggle.ok
              ?("holds it on all "+d.tap_nudge.jiggle.n
                +" directions, residual &le; "
                +d.tap_nudge.jiggle.keep.toFixed(2))
              :"does not hold it either")}
    ${cell3("position in the pool",d.tap_nudge.rank+
            (d.tap_nudge.rank===1?"st":d.tap_nudge.rank===2?"nd":
             d.tap_nudge.rank===3?"rd":"th")+" best",
            (d.tap_nudge.rms-d.tap_nudge.rms_best>=0?"+":"")
            +(d.tap_nudge.rms-d.tap_nudge.rms_best).toFixed(2)+" points")}
  </div>
  <p class="note">${d.tap_nudge.withdrawn?`No usable alternative has been
    established for ${d.name} yet.`:(d.tap_nudge.rms-d.tap_nudge.rms_best)<5?`The cost is small,
    so for ${d.name} holding the lattice is nearly free and the shipped record
    was needlessly broken.`:`<strong>The cost is the whole fit.</strong> For
    ${d.name} an exact reproduction of the elastic constants and a lattice that
    survives a displacement appear to be mutually exclusive within this
    functional form: the best solution that holds its lattice is the
    ${d.tap_nudge.rank}th, at ${d.tap_nudge.rms.toFixed(2)} %. That is a
    statement about the form, not about the search &mdash; three seeds of three
    hundred restarts each found nothing better.`}
  Parameters: m = ${d.tap_nudge.m.toFixed(3)},
  &gamma; = ${d.tap_nudge.gamma.toFixed(4)},
  D = ${d.tap_nudge.D.toFixed(4)},
  &alpha; = ${d.tap_nudge.alpha.toFixed(4)},
  r&#8320; = ${d.tap_nudge.r0.toFixed(4)},
  C = ${d.tap_nudge.C.toFixed(4)},
  &alpha;&#8323; = ${d.tap_nudge.alpha3.toFixed(4)},
  r<sub>cut2</sub> = ${d.tap_nudge.rcut2.toFixed(4)}&nbsp;&#8491;,
  r<sub>cut3</sub> = ${d.tap_nudge.rcut3.toFixed(4)}&nbsp;&#8491;,
  taper = ${d.tap_nudge.taper?d.tap_nudge.taper.toFixed(2):"off"}
  &mdash; all ten, so this is enough to rebuild it without reading the
  file.</p>`:""}

  ${d.hard_disp?`<h3>An alternative fit that describes the phonons better</h3>
  <p class="note">Carried beside the shipped hard-cutoff record rather than
  instead of it, for the same reason as the section above: the objective cannot
  tell the two apart, so choosing one silently would hide the only interesting
  thing about it. Both reach the elastic targets at
  ${d.hard_disp.rms.toFixed(3)}&nbsp;% residual &mdash; the same targets, the
  same <code>refdata.py</code>, nothing reweighted. What separates them is the
  measured dispersion, which never enters the fit.</p>
  <div class="cols3">
    ${cell3("shipped fit, dispersion",
       d.hard_disp.evidence.dispersion_published.toFixed(2)+" %",
       "of the top measured frequency")}
    ${cell3("this fit, dispersion",
       d.hard_disp.evidence.dispersion_here.toFixed(2)+" %",
       ((d.hard_disp.evidence.dispersion_here
         -d.hard_disp.evidence.dispersion_published)>=0?"+":"")
       +(d.hard_disp.evidence.dispersion_here
         -d.hard_disp.evidence.dispersion_published).toFixed(2)+" points")}
    ${cell3("elastic error", d.hard_disp.rms.toFixed(3)+" %",
       "the shipped record is the same")}
  </div>
  <p class="note"><strong>It also brings two constraints back inside the
  line.</strong> <code>fit.py</code> enforces
  E<sub>3</sub>/E<sub>2</sub>&nbsp;&le;&nbsp;0.30 &mdash; the three-body term
  may correct the pair term but may not cancel it &mdash; and a compression
  guard. Both were tightened for the tapered sets and never applied back to
  the hard-cutoff library, so 22 of its 36 records are above the first and the
  shipped record here fails both. This one is at
  ${d.hard_disp.evidence.E3_over_E2_here.toFixed(3)} against the shipped
  ${d.hard_disp.evidence.E3_over_E2_published.toFixed(3)}, and passes the
  compression guard.
  ${d.hard_disp.evidence.md_pe_minus_lattice_here!==undefined?`Under the MD
    screen both survive, but this one sits
    ${d.hard_disp.evidence.md_pe_minus_lattice_here.toFixed(3)}&nbsp;eV/atom
    above its own static lattice against the shipped record's
    ${d.hard_disp.evidence.md_pe_minus_lattice_published.toFixed(3)}, at
    ${d.hard_disp.evidence.md_T_here}&nbsp;K against
    ${d.hard_disp.evidence.md_T_published}&nbsp;K.`:""}</p>
  <p class="note">What has <em>not</em> been checked, and it is not an
  omission: the vacancy, surface and stacking-fault energies. This is a
  hard-cutoff record, &phi;<sub>2</sub> does not vanish at the cutoff, and a
  relaxation walks straight into the discontinuity &mdash; which is why the
  library carries those tests for the switched arms only.
  Parameters: m = ${d.hard_disp.m.toFixed(3)},
  &gamma; = ${d.hard_disp.gamma.toFixed(4)},
  D = ${d.hard_disp.D.toFixed(4)},
  &alpha; = ${d.hard_disp.alpha.toFixed(4)},
  r&#8320; = ${d.hard_disp.r0.toFixed(4)},
  C = ${d.hard_disp.C.toFixed(4)},
  &alpha;&#8323; = ${d.hard_disp.alpha3.toFixed(4)},
  r<sub>cut2</sub> = ${d.hard_disp.rcut2.toFixed(4)}&nbsp;&#8491;,
  r<sub>cut3</sub> = ${d.hard_disp.rcut3.toFixed(4)}&nbsp;&#8491;,
  taper = off &mdash; all ten, so this is enough to rebuild it without
  reading the file. It also ships as
  <code>${d.sym}_disp.ugur</code>.</p>`:""}

  ${lammpsBlock(d)}

  <h3>Parameters</h3>
  <div class="gen">
    ${setsOf(d).length>1?`<select id="pot">
      ${setsOf(d).map(o=>`<option value="${o[0]}"${
          PAR_SET===o[0]?" selected":""}>${
          o[0]==="ug"?"UG (Ugur-Guler)":o[1]} &mdash; ${
          o[2].rms.toFixed(2)}%${armStab(o[2], o[0])}</option>`)
        .join("")}
    </select>`:""}
    <select id="fmt">
      <option value="text" selected>readable</option>
      <option value="json">JSON</option>
      <option value="python">Python (latdyn)</option>
    </select>
    <button class="primary" id="copy">Copy</button>
    <span class="plotnote" id="copied"></span>
  </div>
  <pre class="mono" id="out"></pre>`;

  const pot=$("#pot"), pset=$("#parset");
  /*  Both selectors drive the same PAR_UG, and switching either re-renders, so
      the chip row at the top and the export at the bottom can never end up
      describing different potentials.  They appear only for the fourteen
      elements that have a UG fit. */
  const setPar=v=>{PAR_SET=v; buildTable(); render();};
  if(pot) pot.onchange=()=>setPar(pot.value);
  if(pset) pset.onchange=()=>setPar(pset.value);
  const upd=()=>{$("#out").textContent=
    /*  the exported text has to be the set on screen, whichever it is - it
        used to be hard-wired to d.ug and would have handed out the angular
        parameters while the panel showed a re-cut candidate  */
    paramText(d,cur,$("#fmt").value, isUG?P:null);};
  $("#fmt").onchange=upd;
  upd();
  $("#copy").onclick=async()=>{
    try{await navigator.clipboard.writeText($("#out").textContent);
      $("#copied").textContent="copied to clipboard";}
    catch(e){$("#copied").textContent="copy blocked - select the text manually";}
    setTimeout(()=>{$("#copied").textContent="";},2200);};
  const mps=$("#mprop");
  if(mps){mps.value=MPROP_SEL;
    mps.onchange=()=>{MPROP_SEL=mps.value;drawPolar(d);};}
  const mrm=$("#marm");
  if(mrm){mrm.value=MARM;
    mrm.onchange=()=>{MARM=mrm.value;drawPolar(d);};}
  const dm=$("#dmode");
  if(dm){dm.value=DISP_MODE;
    dm.onchange=()=>{DISP_MODE=dm.value;drawDisp(d);};}
  const etq=$("#etq");
  if(etq){etq.onchange=()=>{ET_Q=etq.value;drawElasticT(d);};}
  plots(d);
}

/*  ---- the mark, in three dimensions -------------------------------------

    A body-centred cubic cell with one triplet inside it, turned as a single
    rigid object: eight vertices, twelve edges, three atoms and the angle arc,
    all defined in model coordinates, rotated about the vertical axis, tilted
    and projected with a weak perspective.  Near edges come out thicker and
    more opaque than far ones; that is the whole depth cue, there is no
    lighting and none is wanted.

    The geometry is not decorative.  The apex atom IS the body-centre atom of
    the cell and its two neighbours ARE two corners of the top face, so the
    apex angle is arccos(1/3) = 70.53 degrees - a real bcc nearest-neighbour
    triplet, of exactly the kind phi3 sums over.  An fcc cell was tried first
    and looked worse for an honest reason: there the nearest-neighbour triplet
    hangs off a corner, so the figure sat low and to one side of its own cell,
    while in bcc the apex sits on the rotation axis and the thing is centred by
    construction.

    It SWAYS rather than spins.  A full turn passes edge-on twice, where the
    triplet collapses to a line and U lands on G.  The amplitude is +-38
    degrees, which is what fits: swept over the whole sway and both tilt
    extremes, the outermost ink - a letter corner, not the cube - reaches 2.2
    and 61.8 of the 64-unit box.

    Three things stop it being a nuisance: prefers-reduced-motion draws one
    static frame, it stops while the tab is hidden, and it stops once it
    scrolls out of view.  A logo is not worth a wakeup every 16 ms on a
    backgrounded tab.                                                         */
(function(){
  const $m = id => document.getElementById(id);
  const mk = $m("mk"); if(!mk) return;
  //  OFF is how far outside its own atom each letter rides, along the leg.
//  At 7 the letters sat on the atoms and read as labels painted on them;
//  at 10.5 they are clearly beside the geometry rather than part of it.
//  The cube is 22 rather than 26 to buy the letters room: at OFF = 10.5 and
//  S = 26 the outermost ink swept to -1.2 and 65.2 of the 64-unit box, so
//  the letters clipped at the extremes of the sway.  Shrinking the cell
//  costs nothing that matters - the apex angle is set by the directions of
//  J and K and both scale with S, so it is unchanged, and the arc has its
//  own radius - while the stylesheet makes the whole mark larger anyway.
//  Swept over the full sway and both tilts the extremes are now 1.7 and
//  62.3, checked rather than assumed.
  const CX = 32, CY = 35.5, D = 130, S = 22, OFF = 10.5, FS = 14.5;
  const A = [0,0,0], J = [-S/2,S/2,-S/2], K = [S/2,S/2,-S/2];
  const V = [];
  for(const x of [0,1]) for(const y of [0,1]) for(const z of [0,1])
    V.push([S*(x-.5), S*(y-.5), S*(z-.5)]);
  const E = [];
  for(let i=0;i<8;i++) for(let j=i+1;j<8;j++){
    let d=0; for(let c=0;c<3;c++) if(V[i][c]!==V[j][c]) d++;
    if(d===1) E.push([i,j]);
  }
  const norm = v => {const L=Math.hypot(v[0],v[1],v[2]); return v.map(t=>t/L);};
  const uJ = norm(J), uK = norm(K), ARC = [];
  //  the label rides on the bisector of the two legs, just outside the arc:
  //  15.5 against the arc's 10, which a sweep over the whole sway puts 12.5
  //  units clear of the nearest initial at its closest approach
  const BIS = norm(norm(J).map((q,i)=>q+norm(K)[i])).map(q=>15.5*q);
  for(let t=0;t<=10;t++){
    const a=t/10, w=[0,1,2].map(i=>uJ[i]*(1-a)+uK[i]*a);
    const L=Math.hypot(w[0],w[1],w[2]);
    //  radius 10 rather than 8: close enough to the vertex to be the angle
    //  and not a smile, far enough out to be seen at this size
    ARC.push(w.map(q=>10*q/L));
  }

  function project(p, phi, tilt){
    let [x,y,z] = p;
    const cp=Math.cos(phi), sp=Math.sin(phi);
    [x,z] = [x*cp + z*sp, -x*sp + z*cp];
    const ct=Math.cos(tilt), st=Math.sin(tilt);
    [y,z] = [y*ct - z*st, y*st + z*ct];
    const s = D/(D-z);
    return [CX + x*s, CY - y*s, s, z];
  }

  function draw(phi, tilt){
    const P = p => project(p, phi, tilt);
    const vp = V.map(P), Ap = P(A), Jp = P(J), Kp = P(K);
    let z0 = Infinity, z1 = -Infinity;
    for(const v of vp){ if(v[3]<z0) z0=v[3]; if(v[3]>z1) z1=v[3]; }
    E.forEach(([i,j],n)=>{
      const el=$m("mkE"+n); if(!el) return;
      const a=vp[i], b=vp[j], f=((a[3]+b[3])/2 - z0)/(z1-z0 || 1);
      el.setAttribute("x1",a[0].toFixed(2)); el.setAttribute("y1",a[1].toFixed(2));
      el.setAttribute("x2",b[0].toFixed(2)); el.setAttribute("y2",b[1].toFixed(2));
      el.style.strokeWidth=(0.7+0.9*f).toFixed(2);
      //  The cell is scaffolding: it says which lattice the triplet is cut
      //  from and should not compete with the legs, the angle or the
      //  initials.  0.28 to 0.62, with the near edges still clearly ahead of
      //  the far ones so the box keeps reading as a solid rather than a
      //  flat outline.
      el.style.opacity=(0.28+0.34*f).toFixed(2);
    });
    $m("mkArc").setAttribute("points",
      ARC.map(p=>{const q=P(p); return q[0].toFixed(2)+","+q[1].toFixed(2);}).join(" "));
    const leg=(el,e)=>{el.setAttribute("x1",Ap[0].toFixed(2));
      el.setAttribute("y1",Ap[1].toFixed(2));
      el.setAttribute("x2",e[0].toFixed(2)); el.setAttribute("y2",e[1].toFixed(2));
      el.style.strokeWidth=(2.8*e[2]).toFixed(2);};
    leg($m("mkLU"),Jp); leg($m("mkLG"),Kp);
    const atom=(el,p,r)=>{el.setAttribute("cx",p[0].toFixed(2));
      el.setAttribute("cy",p[1].toFixed(2)); el.setAttribute("r",(r*p[2]).toFixed(2));
      el.style.opacity=Math.min(1,0.55+0.45*p[2]).toFixed(3);};
    atom($m("mkA"),Ap,3.9); atom($m("mkJ"),Jp,3.0); atom($m("mkK"),Kp,3.0);
    /*  each letter rides just outside its own atom, along the leg, and stays
        upright: it names the neighbour, it is not painted on it */
    const lab=(el,p)=>{
      const dx=p[0]-Ap[0], dy=p[1]-Ap[1], L=Math.hypot(dx,dy)||1;
      el.setAttribute("x",(p[0]+dx/L*OFF*p[2]).toFixed(2));
      el.setAttribute("y",(p[1]+dy/L*OFF*p[2]).toFixed(2));
      el.style.fontSize=(FS*p[2]).toFixed(1)+"px";
      el.style.opacity=Math.min(1,0.6+0.4*p[2]).toFixed(3);};
    lab($m("mkU"),Jp); lab($m("mkG"),Kp);
    //  upright like the initials, and it fades with depth the same way, so
      //  it reads as sitting in the plane of the angle rather than on the glass
    const Hp = P(BIS), h = $m("mkH");
    h.setAttribute("x",Hp[0].toFixed(2)); h.setAttribute("y",Hp[1].toFixed(2));
    h.style.fontSize=(7*Hp[2]).toFixed(1)+"px";
    h.style.opacity=Math.min(1,0.55+0.45*Hp[2]).toFixed(3);
  }

  const still = matchMedia("(prefers-reduced-motion: reduce)");
  let raf = 0, t0 = 0, visible = true;
  const step = t => {
    if(!t0) t0 = t;
    const u = (t - t0)/1000;
    draw(38*Math.PI/180*Math.sin(u*0.40),
         (22 + 6*Math.sin(u*0.24))*Math.PI/180);
    raf = requestAnimationFrame(step);
  };
  const stop = () => {if(raf){cancelAnimationFrame(raf); raf=0;}};
  const start = () => {
    if(raf || still.matches || !visible || document.hidden) return;
    t0 = 0; raf = requestAnimationFrame(step);
  };
  const sync = () => {still.matches ? (stop(), draw(0,28*Math.PI/180)) : start();};
  still.addEventListener("change", sync);
  document.addEventListener("visibilitychange", () => document.hidden?stop():start());
  if(window.IntersectionObserver)
    new IntersectionObserver(es=>{visible = es[0].isIntersecting;
      visible ? start() : stop();}, {threshold:0}).observe(mk);
  sync();
})();

const redraw=()=>{const d=DATA[cur];
  plots(d);};
render();
addEventListener("resize",redraw);
matchMedia("(prefers-color-scheme:dark)").addEventListener("change",redraw);
new MutationObserver(redraw).observe(document.documentElement,
  {attributes:true,attributeFilter:["data-theme"]});
</script>
"""

#  the counts come from the data, so the page cannot drift out of date: every
#  element refdata knows about was attempted, and the ones missing from
#  library.json are the ones with no solution.  Adding an element to refdata.py
#  moves both numbers on the next build with nothing to edit by hand.
ATTEMPTED = sorted(refdata.ELEMENTS)
FAILED = [e for e in ATTEMPTED if e not in DATA]

out = (HTML
       #  Second guard on the AFLOW withdrawal (see add_elastic_T.py).
       #  The merge was removed there, but a library.json from a backup
       #  or an older run still carries the key, and the page would
       #  republish it silently.  Stripped here as well so that the
       #  licence guarantee does not depend on which file was loaded.
       .replace("__FT__", json.dumps(FT, separators=(",", ":")))
       .replace("__ELSRC__", refdata_electronic.GAMMA_SOURCE)
       .replace("__DATA__", json.dumps(
           {e: {k: v for k, v in r.items() if k != "aflow"}
            for e, r in DATA.items()}, separators=(",", ":")))
       .replace("__NOK__", str(len(DATA)))
       .replace("__NTRIED__", str(len(ATTEMPTED)))
       #  the UG count is the data's, not a number typed into the legend: it
       #  was still saying fourteen while the runs that take it to every
       #  element were on the cluster
       #  The dot used to mean "has a UG fit", which marked all 38 and so
       #  marked nothing.  It marks the re-ranked candidates now, and the
       #  legend has to say which - a marker whose caption describes the
       #  previous meaning is worse than no marker.
       .replace("__NUG__", str(sum(1 for v in DATA.values() if v.get("ug"))))
       #  A visible build stamp.  Twice now a change has been reported as
       #  missing when the file on disk already had it and the browser was
       #  showing an older copy; the stamp settles that in one glance instead
       #  of by reading the page's own JavaScript.
       .replace("__BUILT__", _dt.datetime.now().strftime("%Y-%m-%d %H:%M"))
       #  which tree wrote this file, from the path rather than from a
       #  constant, so a copied script cannot mislabel its own output
       .replace("__TREE__", "screen"
                if "_screen" in HERE else "published")
       .replace("__NSET__", str(sum(
           1 for v in DATA.values()
           if isinstance(v, dict) and "rms" in v
           and sum(1 for k in ("ug", "rc", "rc_ug")
                   if isinstance(v.get(k), dict) and "rms" in v[k]) >= 1))))
path = os.path.join(HERE, "potential.html")
open(path, "w", encoding="utf-8").write(out)
print(f"wrote {path}  ({len(out)/1024:.0f} KB, {len(DATA)} elements)")

#  Parse the page's own JavaScript before claiming it was written.  A syntax
#  error here does not raise anything on this side - the file is just text -
#  and the browser reports it to a console nobody has open, so the page comes
#  up blank or half-drawn and looks like a data problem.  esprima is optional;
#  if it is not installed the check says so rather than passing silently.
try:
    import re as _re
    import esprima as _es
    _js = chr(10).join(m.group(1) for m in
                       _re.finditer(r"<script[^>]*>(.*?)</script>",
                                    out, _re.S))
    _es.parseScript(_js)
    print(f"  javascript parses ({len(_js)/1024:.0f} KB)")
except ImportError:
    print("  javascript NOT checked - pip install esprima to enable")
except Exception as _e:
    raise SystemExit(f"  JAVASCRIPT SYNTAX ERROR in the page just written: {_e}")
#  which elements failed is still worth knowing while building the page - it is
#  simply not something the page argues about any more
if FAILED:
    print(f"  no fit: {', '.join(FAILED)}")
