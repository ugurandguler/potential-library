#!/usr/bin/env python3
"""
Build the standalone viewer.  Everything shown is computed by latdyn.py from the
potential's own derivatives - no external code is involved anywhere in the chain.

    python fit.py && python build_library.py && python make_gui.py

Source is deliberately pure ASCII; Greek letters go in as HTML entities.
"""
import json, math, os

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "library.json")))

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
  <strong>CRC Handbook of Chemistry and Physics</strong> thermodynamic tables. The calculated dispersions drawn for comparison come from the
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
  constants are room temperature (except the five alkali metals, which are 5 K
  values); and III/29a states its own convention plainly &mdash; <i>&ldquo;unless
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
  caesium's 78 K rows are &mdash; which, with their lattice constants already at
  5 K, makes them by accident the two most internally consistent records here.
  Nothing has been changed on this account: changing a target means refitting
  the element, and refitting on this account means refitting the library.

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
</footer>
</div>

<script>
const DATA = __DATA__;
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

const pt=$("#pt");
Object.keys(DATA).sort((a,b)=>DATA[a].pos[0]-DATA[b].pos[0]
  ||DATA[a].pos[1]-DATA[b].pos[1]).forEach(el=>{
  const d=DATA[el], b=document.createElement("button");
  b.className="cell"; b.style.setProperty("--r",d.pos[0]-1);
  b.style.setProperty("--c",d.pos[1]);
  b.setAttribute("aria-pressed",el===cur); b.dataset.el=el;
  const unst=d.dyn&&d.dyn.stable===false;
  /*  Only fourteen metals have been fitted with the angular term, so without a
      mark on the table the comparison is invisible: the page opens on an
      element that has none, and nothing says where to look. */
  const hasUG=!!d.ug;
  b.innerHTML=`<i class="tier ${tier(d.rms)}"></i>
    <span class="sym disp">${el}${unst?'<sup style="color:var(--bad)">&#9888;</sup>':""}</span>
    ${hasUG?'<i class="ugdot" title="UG comparison available"></i>':""}
    <span class="st">${d.struct}</span>`;
  b.title=`${d.name} - ${d.struct}, elastic RMS ${d.rms.toFixed(0)}%`
    +(unst?` - dynamically unstable (${(d.dyn.imag_frac*100).toFixed(1)}% imaginary modes)`:"")
    +(hasUG?` - UG fit available, ${d.ug.rms.toFixed(1)}%`:"");
  b.onclick=()=>{cur=el;
    pt.querySelectorAll(".cell").forEach(c=>
      c.setAttribute("aria-pressed",c.dataset.el===el));
    render();};
  pt.appendChild(b);
});

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
    c.textAlign="center";c.fillText(String(Math.round(t)),X(t),H-9);}
  c.textAlign="left";c.fillText("T (K)",L,H-9);
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
    const col = o.r.nudge_bad ? bad : (ours ? (o.k==="tap"?p2:p3) : ref);
    const pts=o.r.pts.filter(q=>isFinite(q[ET_Q]));
    /*  Line style carries the identity, colour only the warning.  The three
        colours here sit within 1.4 of each other in luminance, so on a
        greyscale print or to a colour-blind reader they are one curve drawn
        three times; the dash pattern is what survives that.  It also matches
        the potential panel, where MAU is already solid and UG dashed.       */
    const dash = o.r.kind==="base" ? [1.5,2.5] : (o.k==="tap_ug" ? [6,3] : []);
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
    lab.push({y:Y(q[ET_Q]), x:X(q.T),
              col:(o.r.nudge_bad?bad:(o.r.kind==="ours"?(o.k==="tap"?p2:p3):ref)),
              txt:txt+(o.r.nudge_bad?" ⚠":"")});
  });
  lab.sort((a,b)=>a.y-b.y);
  for(let i=1;i<lab.length;i++)
    if(lab[i].y-lab[i-1].y<12) lab[i].y=lab[i-1].y+12;
  const over=lab.length?lab[lab.length-1].y-(T+ph):0;
  if(over>0) lab.forEach(l=>l.y-=over);
  c.font="11px system-ui,-apple-system,sans-serif";
  c.textAlign="left";
  lab.forEach(l=>{ c.fillStyle=l.col; c.fillText(l.txt, W-R+6, l.y+4); });

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
  const recs=[["tap",d.tap],["tap_ug",d.tap_ug]]
    .filter(o=>o[1]&&o[1].bain&&o[1].bain.E&&o[1].bain.E.length);
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
    c.textAlign="center";c.fillText(x.toFixed(2),X(x),H-9);}
  c.textAlign="left";c.fillText("tetragonal strain δ",L,H-9);
  c.save();c.translate(13,T+ph/2);c.rotate(-Math.PI/2);c.textAlign="center";
  c.fillText("E − E₀  (meV/atom)",0,0);c.restore();
  frame(c,L,T,pw,ph,ink,[0,1,2,3,4].map(k=>X(-dmax+2*dmax*k/4)),
        [0,1,2,3,4].map(k=>Y(lo+(hi-lo)*k/4)));

  //  the zero line is the whole reading: below it, bcc is not the minimum
  c.save();c.setLineDash([2,3]);c.strokeStyle=ink;
  c.beginPath();c.moveTo(L,Y(0));c.lineTo(W-R,Y(0));c.stroke();c.restore();

  recs.forEach(o=>{
    const b=o[1].bain, col=o[0]==="tap"?p2:p3;
    c.strokeStyle=col;c.lineWidth=2;c.save();
    c.setLineDash(o[0]==="tap_ug"?[6,3]:[]);
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
  const recs=[["tap",d.tap],["tap_ug",d.tap_ug]]
    .filter(o=>o[1]&&o[1].stacking&&o[1].stacking.gamma
                &&o[1].stacking.gamma.length);
  const bases=Object.entries(d.baseline_stacking||{})
    .filter(o=>o[1]&&o[1].gamma&&o[1].gamma.length);
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
    c.textAlign="center";c.fillText(t[1],X(t[0]),H-9);});
  c.textAlign="left";
  c.fillText("shift along [11̄2̄], in periods",L,H-9);
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
    const g=o[1].stacking, col=o[0]==="tap"?p2:p3;
    c.strokeStyle=col;c.lineWidth=2;c.save();
    c.setLineDash(o[0]==="tap_ug"?[6,3]:[]);
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
  const recs=[["tap",d.tap],["tap_ug",d.tap_ug]]
    .filter(o=>o[1]&&o[1].expansion&&o[1].expansion.T&&o[1].expansion.T.length>1);
  const bases=Object.entries(d.baseline_expansion||{})
    .filter(o=>o[1]&&o[1].T&&o[1].T.length>1);
  if(!recs.length&&!bases.length) return;
  const s=setup("#expan",0.52); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=62,R=16,T=12,B=30, pw=W-L-R, ph=H-T-B;

  const series=recs.map(o=>({g:o[1].expansion,col:o[0]==="tap"?p2:p3,
                             dash:o[0]==="tap_ug"?[6,3]:[],wide:2}))
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
    c.textAlign="center";c.fillText(t.toFixed(0),X(t),H-9);}
  c.textAlign="left";c.fillText("T (K)",L,H-9);
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
  const jobs=[["dispersion",drawDisp],["thermo",drawThermo],
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
  const u=ugData(d), hr_=u?hRange(u):null;
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
  const walls=[wallOf(d)].concat(u?[wallOf(u)]:[]).filter(w=>w!==null);
  const r0=walls.length?Math.max(0.45*d.dnn,
              Math.min(0.78*d.dnn, Math.min.apply(null,walls)-0.04*d.dnn))
            :0.78*d.dnn,
        r1=2.05*d.dnn, N=380;
  const xs=[],y2=[],y3=[],u2=[],u3lo=[],u3hi=[];
  for(let i=0;i<N;i++){const r=r0+(r1-r0)*i/(N-1);
    xs.push(r); y2.push(phi2(r,d)); y3.push(phi3(2*r,d));
    if(u){const b=phi3(2*r,u); u2.push(phi2(r,u));
      u3lo.push(Math.min(b*hr_[0],b*hr_[1]));
      u3hi.push(Math.max(b*hr_[0],b*hr_[1]));}}
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
  const panel=(T,sets,label)=>{
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
    frame(c,L,T,pw,ph,ink,[0,1,2,3,4].map(k=>X(r0+(r1-r0)*k/4)),
          [0,1,2,3].map(k=>Y(lo+(hi-lo)*k/3)));
    c.fillStyle=ink;c.textAlign="left";c.fillText(label,L+7,T+13);
  };

  panel(TOP, u?[[y2,p2,2.2],[u2,p2,1.6,[5,3]]]:[[y2,p2,2.2]],
        "φ₂(r)");
  if(u) key(TOP, [["MAU","line",p2],["UG","dash",p2]]);
  panel(TOP+ph+GAP,
        u?[[u3hi,p3,0,null,u3lo],[y3,p3,1.8]]:[[y3,p3,1.8]],
        "φ₃(r,r)");
  if(u) key(TOP+ph+GAP, [["MAU","line",p3],
                         ["UG, all θ","band",p3]]);

  /*  one axis label for both panels: the quantity is the same in each and the
      unit was previously carried inside the panel captions, where it read as
      part of the function name */
  c.save();c.translate(14,TOP+ph+GAP/2);c.rotate(-Math.PI/2);
  c.textAlign="center";c.fillStyle=ink;c.fillText("Energy (eV)",0,0);c.restore();
  c.fillStyle=ink;c.textAlign="center";
  for(let k=0;k<=4;k++){const r=r0+(r1-r0)*k/4;c.fillText(r.toFixed(2),X(r),H-10);}
  c.textAlign="left";c.fillText("r (A)",L,H-10);
  c.textAlign="center";c.fillText("d_nn",X(d.dnn),H-10);
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
let PAR_UG = false;
const parSet = d => (PAR_UG && d.ug) ? d.ug : d;

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
  const g=d.ld||{};
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
  if(DISP_MODE==="ug" && !ugData(d)) DISP_MODE=null;   // element without one
  if(DISP_MODE===null||!(DISP_MODE==="std"||DISP_MODE==="ug"||has(DISP_MODE)))
    DISP_MODE = has("mp") ? "mp"
              : (has("mc3d") ? "mc3d" : (has("jarvis") ? "jarvis" : "std"));
  /*  The DFT references are stored at THEIR q-points and UG at ours, so the
      two cannot share an x axis; "ug" is therefore a view of the standard path
      with the second potential drawn on it, not a reference mode. */
  const withUG = (DISP_MODE === "ug") && ugData(d);
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
    (mpp.marks||[]).forEach(([i,name])=>vline(xs[i],name));
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
    (mpp.marks||[]).forEach(([i,name])=>marker(xs[i],name,Y));
    return;
  }

  /* ---- our curve on the same path, no MP reference available ---- */
  const segs = g.std||[];
  if(!segs.length) return;
  const u=withUG?ugData(d):null, usegs=(u&&u.ld)?u.ld.std:null;
  let hi=0, lo=0;
  const span=set=>set.forEach(p=>p.branches.forEach(b=>b.forEach(v=>{
    if(v===null)return; if(v>hi)hi=v; if(v<lo)lo=v;})));
  span(segs); if(usegs) span(usegs);
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
    x0 += w + (brk[i]?GAP:0);
  });
  if(usegs) legend(c,L+8,T+14,ink,p2,p3);
}

/*  two-entry key, drawn inside the frame so it cannot be clipped */
function legend(c,x,y,ink,p2,p3){
  c.font="11px ui-monospace,Consolas,monospace";c.textAlign="left";
  [["MAU",p2,[]],["UG",p3,[5,3]]].forEach(([t,col,dash],i)=>{
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

function drawPolar(d){
  const pl=d.mech_planes; if(!pl) return;
  const s=setup("#polar",0.36); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const spec=MPROP[MPROP_SEL], sc=spec.scale||1;
  const planes=["xy","xz","yz"];
  /*  Drawn whenever a comparable UG fit exists, not behind the overlay switch:
      the table beside this panel now shows both potentials unconditionally, and
      a plot that disagreed with the table next to it would be a trap.  The
      switch stays for the dispersion and the thermodynamics, where the second
      curve changes how the first is read. */
  const u=(d.ug&&d.ug.comparable)?d.ug:null, upl=u?u.mech_planes:null;
  /* one radial scale for all three panels, so they are comparable - and the
     same scale for both potentials, or the shapes could not be compared */
  let vmax=-Infinity, vmin=Infinity;
  [pl,upl].forEach(src=>{ if(!src) return;
    planes.forEach(p=>spec.keys.forEach(k=>(src[p]||{})[k]&&src[p][k].forEach(v=>{
      const x=v*sc; if(x>vmax)vmax=x; if(x<vmin)vmin=x;})));});
  const floor=Math.min(0,vmin), span=(vmax-floor)||1;
  const pad=26, w=W/3, R=Math.min(w,H)/2-pad;
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
    c.fillText(p[0]+"-"+p[1]+" plane",cx,H-7);
  });
  c.fillStyle=ink;c.textAlign="left";
  c.fillText(`${(vmax).toFixed(vmax<10?2:0)} ${spec.unit}`.trim()+" outer ring"
             +(upl?"   (thin = UG)":""), 6,12);
}

function drawThermo(d){
  const g=d.ld||{}; if(!g.thermo) return;
  const s=setup("#thermo"); if(!s) return;
  const {c,W,H,p2,p3,lin,ink}=s;
  const L=48,R=10,T=10,B=26, pw=W-L-R, ph=H-T-B;
  const th=g.thermo.filter(r=>r.S!==undefined&&r.Cv!==undefined);
  if(!th.length) return;
  const u=ugData(d);          // unconditional, as for the polar panels
  const uth=(u&&u.ld&&u.ld.thermo)
            ?u.ld.thermo.filter(r=>r.S!==undefined&&r.Cv!==undefined):null;
  const Tm=th[th.length-1].T;
  let hi=Math.max(...th.map(r=>Math.max(r.S,r.Cv)), d.S298||0, d.Cp298||0,
                  ...(uth?uth.map(r=>Math.max(r.S,r.Cv)):[0]));
  hi=Math.ceil(hi/10)*10;
  const X=t=>L+t/Tm*pw, Y=v=>T+ph-v/hi*ph;
  c.strokeStyle=lin;c.lineWidth=1;
  c.font="11px ui-monospace,Consolas,monospace";c.fillStyle=ink;
  for(let k=0;k<=4;k++){const v=hi*k/4,y=Y(v);
    c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();
    c.textAlign="right";c.fillText(String(Math.round(v)),L-6,y+4);}
  for(let k=0;k<=4;k++){const t=Tm*k/4;
    c.textAlign="center";c.fillText(String(Math.round(t)),X(t),H-8);}
  c.textAlign="left";c.fillText("T (K)",L,H-8);
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
  if(uth){line(uth,"S",p2,1.3,[5,3]); line(uth,"Cv",p3,1.3,[5,3]);}
  line(th,"S",p2,2,[]); line(th,"Cv",p3,2,[]);
  if(uth){c.fillStyle=ink;c.textAlign="left";
    c.font="11px ui-monospace,Consolas,monospace";
    c.fillText("dashed = UG",L+6,T+12);}
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
function lammpsBlock(d){
  const SETS = [
    ["hard", d,        d.name+".ugur",             "ugur",     "hard truncation"],
    ["tap",  d.tap,    d.name+"_taper.ugur",       "ugur",     "switched"],
    ["ug",   d.ug,     d.name+".ugur.ang",         "ugur/ang", "hard + angular"],
    ["tap_ug", d.tap_ug, d.name+"_taper.ugur.ang", "ugur/ang", "switched + angular"]
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

/*  value at the tabulated temperature nearest 298 K, as a number */
function tvn(g,key){
  if(!g||!g.thermo) return null;
  let best=g.thermo[0];
  for(const r of g.thermo) if(Math.abs(r.T-298)<Math.abs(best.T-298)) best=r;
  return best[key];
}

function erow(label,pair,frozen,mpval,unc){
  if(!pair) return "";
  const calc=Math.abs(pair[0]), exp=Math.abs(pair[1]);
  const err=exp?100*(calc-exp)/exp:null;
  const w=err===null?0:Math.min(Math.abs(err),40)*1.4;
  const fz=(frozen===undefined||frozen===null)?"&mdash;":Math.abs(frozen).toFixed(1);
  const mv=(mpval===undefined||mpval===null)?"&mdash;":Math.abs(mpval).toFixed(1);
  return `<tr><td>${label}</td>
    <td class="mono"><strong>${calc.toFixed(1)}</strong></td>
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
  const dyn=d.dyn||{}, mech=d.mech, rch=d.reach;
  /*  The UG mechanics come from the same producer as the MAU ones, so the two
      are the same quantity computed from two elastic tensors and belong in one
      table rather than in a second one further down the page.  Absent for the
      elements with no angular fit, and the columns simply do not appear. */
  const umech=(d.ug&&d.ug.comparable)?d.ug.mech:null;
  const fcbad=(d.fc_check||[]).length>0;
  const P=parSet(d), isUG=(P!==d);
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
    <span class="rms ${t}">elastic RMS ${d.rms.toFixed(0)}%</span>
  </div>

  <h3>Fitted parameters${d.ug?` <select id="parset" class="inlinesel">
      <option value="mau"${PAR_UG?"":" selected"}>MAU &mdash; ${d.rms.toFixed(2)}%</option>
      <option value="ug"${PAR_UG?" selected":""}>UG &mdash; ${d.ug.rms.toFixed(2)}%</option>
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
    ${isUG?`<div class="par p3"><span class="k">&lambda;<sub>2</sub></span>
      <span class="v mono">${fmt(P.lam2,3)}</span></div>
    <div class="par p3"><span class="k">&lambda;<sub>4</sub></span>
      <span class="v mono">${fmt(P.lam4,3)}</span></div>`:""}
    <div class="par p3"><span class="k">triplets / atom</span>
      <span class="v mono">${P.ntrip||d.ntrip||"&mdash;"}</span></div>
  </div>
  ${isUG?`<p class="note">These are the <strong>UG</strong> parameters &mdash;
  a separate fit, not MAU plus two extra numbers. Every one of them differs:
  for this element D is ${fmt(d.D,3)} eV under MAU against ${fmt(d.ug.D,3)}
  here. Setting &lambda;<sub>2</sub> = &lambda;<sub>4</sub> = 0 recovers the MAU
  <em>form</em> exactly, but not these values.</p>`:""}
  ${dyn.stable===false?`<div class="warn">
    <strong>Dynamically unstable.</strong> ${(dyn.imag_frac*100).toFixed(1)}% of
    modes on the ${String(dyn.nq).split("+").map(n=>n+"&sup3;").join(" and ")}
    meshes are imaginary, down to
    ${(dyn.most_neg_cm1/CM1_PER_THZ).toFixed(2)} THz. The fit only ever sees elastic
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
      <table><thead><tr><th></th><th>fit</th><th>frozen ion</th>
        <th>MP (DFT)</th><th>exp</th><th>error</th></tr></thead><tbody>
        ${erow("E_coh (eV)",d.Ecoh,null,null,d.Ecoh_unc)}
        ${erow("B",d.B,null,me.K)}
        ${erow("C11",d.C11,fz.C11,me.C11,(d.Cij_unc||{}).C11)}
        ${erow("C12",d.C12,fz.C12,me.C12,(d.Cij_unc||{}).C12)}
        ${d.struct==="hcp"?erow("C13",d.C13,fz.C13,me.C13):""}
        ${d.struct==="hcp"?erow("C33",d.C33,fz.C33,me.C33):""}
        ${erow("C44",d.C44,fz.C44,me.C44,(d.Cij_unc||{}).C44)}
      </tbody></table>
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
      ships. Over all 38 elements the median error moves the wrong way,
      5.83 &rarr; 7.67 %: the hcp metals go from 7.86 to 19.81 and the alkalis
      lose heavily. Cadmium and zinc, which fail on axial anisotropy rather
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
        ${umech?`<tr><th></th><th colspan="3"><span class="tag mau">MAU</span></th>
          <th colspan="3"><span class="tag ug">UG</span></th></tr>`:""}
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
                umech?umech.E_H.toFixed(1):null,"UG")}
        ${cell3("Poisson (Hill)",mech.nu_H.toFixed(3),
                umech?umech.nu_H.toFixed(3):null,"UG")}
        ${cell3("anisotropy A<sup>U</sup>",mech.A_U.toFixed(3),
                umech?umech.A_U.toFixed(3):null,"UG")}
        ${cell3("Pugh B/G",mech.pugh.toFixed(2),
                umech?umech.pugh.toFixed(2):null,"UG")}
        ${cell3("Cauchy C12-C44",mech.cauchy.toFixed(1)+" GPa",
                umech?umech.cauchy.toFixed(1):null,"UG")}
        ${mech.debye?cell3("Debye temperature",mech.debye.toFixed(0)+" K",
                umech&&umech.debye?umech.debye.toFixed(0):null,"UG"):""}
        ${mech.v_l?cell3("v longitudinal",(mech.v_l/1000).toFixed(2)+" km/s",
                umech&&umech.v_l?(umech.v_l/1000).toFixed(2):null,"UG"):""}
        ${mech.v_t?cell3("v transverse",(mech.v_t/1000).toFixed(2)+" km/s",
                umech&&umech.v_t?(umech.v_t/1000).toFixed(2):null,"UG"):""}
      </div>
      <p class="plotnote">A<sup>U</sup> = 5G<sub>V</sub>/G<sub>R</sub> +
      B<sub>V</sub>/B<sub>R</sub> &minus; 6 vanishes only for an isotropic
      solid. B/G above ~1.75 is the usual ductility indicator, and a positive
      Cauchy pressure points the same way. The Debye temperature comes from the
      Hill averages and the density, not from the phonon spectrum &mdash; the
      two are independent estimates.</p>
    </div>
    <div>
      <div class="gen" style="margin-bottom:8px">
        <select id="mprop">
          <option value="E">Young's modulus E</option>
          <option value="beta">linear compressibility</option>
          <option value="G">shear modulus G</option>
          <option value="nu">Poisson's ratio</option>
        </select>
      </div>
      <canvas id="polar"></canvas>
      <p class="plotnote">Radius is the value in that direction, so a circle
      means isotropy in the plane.<i class="swatch"
      style="background:var(--phi3);margin-left:8px"></i>dashed: the upper
      envelope over the transverse direction, for the two quantities that need
      one. Dotted circle marks zero where the range crosses it.</p>
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
        </select>
      </div>`:""}
      <canvas id="disp"></canvas>
      ${d.exp_phonon?`<p class="plotnote"><span style="display:inline-block;
        width:9px;height:9px;border:1.4px solid currentColor;border-radius:50%;
        margin-right:5px"></span><strong>Open circles: measured
        frequencies</strong> at
        ${d.exp_phonon.points.map(p=>p.name).join(" and ")}, mean absolute
        error <strong>${d.exp_phonon.mae}%</strong>. Nothing at finite q enters
        the fit, so these are out of sample.<br>Source: ${d.exp_phonon.ref}.
        </p>`:""}
      <p class="plotnote">${g.std[0].branches.length} branches. From the
      dynamical matrix &mdash; phonons are not fit targets, so this is a
      prediction. The dashed curve is whichever reference the selector shows,
      with our dispersion re-evaluated at that reference's own q-points, so the
      residual is physics rather than interpolation.
      ${[["Materials Project ("+(mpph&&mpph.method||"")+")",mpph],
         ["Materials Cloud MC3D (PBEsol)",d.mc3d]]
        .filter(r=>r[1]&&r[1].stats).map(r=>`<br>
        <i class="swatch" style="background:var(--phi3)"></i>vs <strong>${r[0]}
        </strong>: RMS ${(r[1].stats.rms_cm1/CM1_PER_THZ).toFixed(2)} THz
        (${r[1].stats.rel_pct}% of their highest branch), top frequency
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
        are out of sample.`
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
       d.Cp298?exact(d.Cp298)+" (C_p)":null)}
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
  The tabulated experimental value is C<sub>p</sub>; the
  calculation gives C<sub>v</sub>. For metals near 300 K the difference is about
  1-2 J/(mol K), so a slightly lower C<sub>v</sub> is expected. No thermal
  quantity enters the fit.</p>`:""}

  ${d.elasticT?`<h3>Elastic constants against temperature
    <select id="etq" class="inlinesel">${ET_QS.map(a=>
      `<option value="${a[0]}"${a[0]===ET_Q?" selected":""}>${a[1]}</option>`
      ).join("")}</select></h3>
  <canvas id="elasT"></canvas>
  <p class="plotnote">
    <i class="swatch" style="background:var(--phi2)"></i>MAU (solid)
    &nbsp;&middot;&nbsp;
    <i class="swatch" style="background:var(--phi3)"></i>UG (dashed)
    ${Object.values(d.elasticT).some(r=>r.kind==="base")?`
      &nbsp;&middot;&nbsp;<i class="swatch" style="background:var(--ref)"></i>
      published potentials shipped with LAMMPS (dotted)
      (${Object.values(d.elasticT).filter(r=>r.kind==="base")
         .map(r=>r.label+" &mdash; <code>"+r.file+"</code>").join("; ")}),
      run through the identical cell, recipe
      and post-processing &mdash; without them these curves would have nothing
      to be measured against.`:""}
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
    <i class="swatch" style="background:var(--phi2)"></i>MAU (solid)
    ${d.tap_ug&&d.tap_ug.bain?`&nbsp;&middot;&nbsp;
      <i class="swatch" style="background:var(--phi3)"></i>UG (dashed)`:""}
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
  <table><thead><tr><th>facet</th><th>MAU</th>
    ${d.tap_ug&&d.tap_ug.surface?"<th>UG</th>":""}
    <th>DFT</th>${d.surface_ref&&d.surface_ref.tyson?"<th>experiment</th>":""}
    ${d.baseline_surface&&Object.keys(d.baseline_surface).length
      ?`<th>published${Object.keys(d.baseline_surface).length>1
        ?" (lowest&ndash;highest)":""}</th>`:""}</tr></thead><tbody>
  ${d.tap.surface.order_want.map(f=>{
    const ours=d.tap.surface.gamma[f];
    const ug=d.tap_ug&&d.tap_ug.surface?d.tap_ug.surface.gamma[f]:null;
    const dft=(d.surface_ref&&d.surface_ref.facets)?d.surface_ref.facets[f]:null;
    const bs=Object.values(d.baseline_surface||{})
      .map(s=>s.gamma[f]).filter(v=>v!==undefined&&v!==null);
    return `<tr><td>(${f})</td>
      <td>${ours!==undefined&&ours!==null?ours.toFixed(3):"&mdash;"}</td>
      ${d.tap_ug&&d.tap_ug.surface
        ?`<td>${ug!==undefined&&ug!==null?ug.toFixed(3):"&mdash;"}</td>`:""}
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
       "should be " + d.tap.surface.order_want.map(f=>"("+f+")").join(" < "))}
    ${cell3("anisotropy",
       (100*d.tap.surface.spread).toFixed(1) + " %",
       (d.surface_ref&&d.surface_ref.anisotropy!==undefined
         &&d.surface_ref.anisotropy!==null)
         ? "DFT " + (100*d.surface_ref.anisotropy).toFixed(1) + " %" : null)}
  </div>

  <p class="note">${d.tap.surface.order_ok?`The ordering is right: the
    close-packed face is the cheapest, as it is in the metal.`:`<strong
    style="color:var(--bad)">The ordering is wrong.</strong> In a real metal the
    close-packed face is the cheapest, and every published potential tested here
    reproduces that; this record does not.`}
  It is the same shortage the vacancy energy exposes, and the usual account of
  it &mdash; that the energy of a bond in this form cannot depend on how many
  other bonds an atom has &mdash; is <em>too strong</em>. A three-body term
  counts neighbour pairs, so it does carry a coordination dependence, and a
  two- plus three-body potential with free radial shapes fitted to
  density-functional energies and forces reproduces tungsten&rsquo;s surface
  energies to within 4 % (Xie, Rupp and Hennig, npj Comput. Mater. <strong>9
  </strong>, 162 (2023)).
  So what these numbers measure is this parameterisation
  &mdash; parameters fitted to bulk elastic data, which never saw a surface
  &mdash; rather than the functional form.</p>`:""}

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
  <table><thead><tr><th></th><th>MAU</th>
    ${d.tap_ug&&d.tap_ug.stacking?"<th>UG</th>":""}
    ${d.tap.stacking.exp?"<th>experiment</th>":""}
    ${Object.keys(d.baseline_stacking||{}).length
      ?`<th>published${Object.keys(d.baseline_stacking).length>1
        ?" (lowest&ndash;highest)":""}</th>`:""}</tr></thead><tbody>
  ${[["intrinsic fault &gamma;<sub>isf</sub>","isf"],
     ["unstable fault &gamma;<sub>usf</sub>","usf"]].map(row=>{
    const key=row[1];
    const ours=d.tap.stacking[key];
    const ug=d.tap_ug&&d.tap_ug.stacking?d.tap_ug.stacking[key]:null;
    const bs=Object.values(d.baseline_stacking||{})
      .map(s=>s[key]).filter(v=>v!==undefined&&v!==null);
    return `<tr><td>${row[0]}</td>
      <td${key==="isf"&&ours<0?' style="color:var(--bad)"':""}>${
        ours!==undefined&&ours!==null?ours.toFixed(1):"&mdash;"}</td>
      ${d.tap_ug&&d.tap_ug.stacking
        ?`<td>${ug!==undefined&&ug!==null?ug.toFixed(1):"&mdash;"}</td>`:""}
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
    <i class="swatch" style="background:var(--phi2)"></i>MAU (solid)
    ${d.tap_ug&&d.tap_ug.stacking?`&nbsp;&middot;&nbsp;
      <i class="swatch" style="background:var(--phi3)"></i>UG (dashed)`:""}
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
    <strong>9</strong>, 162 (2023)).</p>`;})()}
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
    <i class="swatch" style="background:var(--phi2)"></i>MAU (solid)
    ${d.tap_ug&&d.tap_ug.expansion?`&nbsp;&middot;&nbsp;
      <i class="swatch" style="background:var(--phi3)"></i>UG (dashed)`:""}
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
      here, which are the 97th (2016).`:""}</p>

  ${(()=>{const a=d.tap.expansion.alpha_1e6, r=d.tap.expansion.ratio;
    const b=Object.values(d.baseline_expansion||{}).map(x=>x.ratio)
             .filter(v=>v);
    const scale = b.length
      ? ` For ${d.name}, ${b.length} published potential${b.length>1?"s go":
          " goes"} through the identical barostat and land${b.length>1?"":"s"}
          at ${Math.min(...b).toFixed(2)}&ndash;${Math.max(...b).toFixed(2)} of
          experiment; across the library their median is 0.96 against our
          0.69.` : "";
    if(a<0) return `<p class="note"><strong style="color:var(--bad)">This
      record contracts on heating.</strong> A negative expansion coefficient
      for a simple metal is not a small error, and it is not a subtlety of the
      barostat: nine of our 74 records come back negative and <em>none</em> of
      the 51 published potentials does. All nine are body-centred.${scale}
      Thermal expansion is the third derivative of the same energy curve whose
      second derivative was fitted, so nothing constrains it &mdash; getting
      the curvature right at one volume says nothing about how it changes with
      volume.</p>`;
    return `<p class="note">The library sits at a median 0.69 of the measured
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
  C = ${d.tap_nudge.C.toFixed(4)}.</p>`:""}

  ${lammpsBlock(d)}

  <h3>Parameters</h3>
  <div class="gen">
    ${d.ug?`<select id="pot">
      <option value="mau"${PAR_UG?"":" selected"}>MAU &mdash; ${d.rms.toFixed(2)}%</option>
      <option value="ug"${PAR_UG?" selected":""}>UG (Ugur-Guler) &mdash; ${d.ug.rms.toFixed(2)}%</option>
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
  const setPar=v=>{PAR_UG=(v==="ug"); render();};
  if(pot) pot.onchange=()=>setPar(pot.value);
  if(pset) pset.onchange=()=>setPar(pset.value);
  const upd=()=>{$("#out").textContent=
    paramText(d,cur,$("#fmt").value, isUG?d.ug:null);};
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
       .replace("__DATA__", json.dumps(
           {e: {k: v for k, v in r.items() if k != "aflow"}
            for e, r in DATA.items()}, separators=(",", ":")))
       .replace("__NOK__", str(len(DATA)))
       .replace("__NTRIED__", str(len(ATTEMPTED)))
       #  the UG count is the data's, not a number typed into the legend: it
       #  was still saying fourteen while the runs that take it to every
       #  element were on the cluster
       .replace("__NUG__", str(sum(1 for v in DATA.values() if v.get("ug")))))
path = os.path.join(HERE, "potential.html")
open(path, "w", encoding="utf-8").write(out)
print(f"wrote {path}  ({len(out)/1024:.0f} KB, {len(DATA)} elements)")
#  which elements failed is still worth knowing while building the page - it is
#  simply not something the page argues about any more
if FAILED:
    print(f"  no fit: {', '.join(FAILED)}")
