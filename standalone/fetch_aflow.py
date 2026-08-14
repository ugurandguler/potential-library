#!/usr/bin/env python3
"""
AFLOW's elastic library, and what it is and is not good for here.

    ITS OUTPUT MAY NOT GO INTO THE PUBLISHED PAGE.

AFLOW's terms are "free for scientific, academic and non-commercial purposes.
Any other use is prohibited" (aflowlib.duke.edu), and they require citing
Curtarolo et al., Comput. Mater. Sci. 58, 227 (2012).  The library page is
distributed under CC-BY 4.0, which permits commercial use - so embedding AFLOW
values in it would grant a right nobody here holds.  The three other
comparison sources are compatible and stay: Materials Project and Materials
Cloud MC3D are CC-BY 4.0, JARVIS-DFT is a US Government work.

So this script still runs and is still useful - checking our elastic moduli
against an independent density-functional source is exactly the use AFLOW
permits.  What may not happen is the result being redistributed.  The merge
into library.json was removed in add_elastic_T.py and make_gui.py strips the
key a second time; both carry the reason.  If you re-enable either, you are
changing the licence of the page, and LICENSE.md has to change with it.

The cost of the withdrawal was measured rather than assumed: 22 elements had
an AFLOW record, only 14 were usable - the other 8 are in the wrong phase,
AFLOW's silver being hcp where this library's is fcc - and Materials Project
already supplies elastic data for 33 elements.  What AFLOW alone contributed
was potassium and sodium.

This was looked at once before, early on, while surveying phonon references,
and dismissed in a sentence that never reached the code - so the
conclusion could not be reproduced or checked.  This is that check, written
down.

What it holds: AEL gives Voigt-Reuss-Hill B and G from density functional
theory.  Not the C_ij tensor, not anything at finite temperature.  So it can
contribute a single 0 K reference point, not a curve, and only for two of the
quantities the finite-temperature work reports.

How sparse it is: of thirty-eight elements, a scan of the elemental entries
finds elastic data for a handful, and the handful is erratic - hafnium,
scandium, titanium, barium and strontium have it while copper, aluminium and
magnesium have none.  Several of those are in the wrong phase for this library
(barium as hcp against our bcc, ytterbium as bcc against our fcc), which
matters because an elastic modulus belongs to a structure.

Net: it closes two of the nine elements that have no published potential at
all - hafnium and scandium - and both agree with experiment to better than
1.5 %, which is what makes them worth having.

The AFLUX filter syntax is the trap.  A bare keyword projects the field and
returns nulls; `keyword(*)` and `keyword(1*)` return nothing; `keyword(1*1000)`
answers "DB Fail". The form that works is an open upper bound, `keyword(*max)`,
and without it a scan reports zero everywhere and looks like an empty database.

    python fetch_aflow.py            # scan and write aflow_elastic.json
"""
import json
import os
import subprocess
import sys

import refdata

HERE = os.path.dirname(os.path.abspath(__file__))
API = "http://aflow.org/API/aflux/"
SG = {225: "fcc", 229: "bcc", 194: "hcp", 191: "hcp", 166: "hcp"}


def query(el):
    q = (f"?species({el}),nspecies(1),ael_bulk_modulus_vrh(*100000),"
         f"ael_shear_modulus_vrh,spacegroup_relax,$paging(1,50)")
    #  single quotes: the query contains $paging, and inside double quotes
    #  bash expands it to nothing, which turns every scan into "no data"
    r = subprocess.run(["wsl", "-e", "bash", "-lc",
                        "curl -sL --max-time 90 '" + API + q + "'"],
                       capture_output=True, text=True)
    try:
        return json.loads(r.stdout)
    except Exception:
        return []


def main():
    els = sys.argv[1:] or sorted(refdata.ELEMENTS)
    out = {}
    print(f"{'el':4s}{'ours':>7s}{'AFLOW':>8s}{'B':>8s}{'G':>8s}{'expt B':>9s}"
          f"   verdict")
    print("-" * 56)
    for el in els:
        rows = [x for x in query(el)
                if x.get("ael_bulk_modulus_vrh") is not None]
        if not rows:
            continue
        e = refdata.ELEMENTS[el]
        c = e["Cij"]
        bexp = ((2 * (c["C11"] + c["C12"]) + 4 * c["C13"] + c["C33"]) / 9
                if e["struct"] == "hcp" else (c["C11"] + 2 * c["C12"]) / 3)
        #  an elastic modulus belongs to a structure; a value for the wrong
        #  phase is not a reference for ours
        same = [x for x in rows
                if SG.get(x.get("spacegroup_relax")) == e["struct"]]
        best = same[0] if same else rows[0]
        rec = {"B": best["ael_bulk_modulus_vrh"],
               "G": best.get("ael_shear_modulus_vrh"),
               "spacegroup": best.get("spacegroup_relax"),
               "struct": SG.get(best.get("spacegroup_relax")),
               "usable": bool(same), "B_exp": round(bexp, 1),
               "source": "AFLOW AEL (DFT), aflow.org"}
        out[el] = rec
        print(f"{el:4s}{e['struct']:>7s}{str(rec['struct']):>8s}"
              f"{rec['B']:8.1f}{(rec['G'] or 0):8.1f}{bexp:9.1f}   "
              f"{'kullanilabilir' if same else 'YANLIS FAZ'}")
    json.dump(out, open(os.path.join(HERE, "aflow_elastic.json"), "w"),
              indent=1, sort_keys=True)
    ok = [k for k, v in out.items() if v["usable"]]
    print(f"\n{len(out)} elementte veri, {len(ok)} tanesi dogru fazda: "
          f"{' '.join(ok)}")
    print("-> aflow_elastic.json")


if __name__ == "__main__":
    main()
