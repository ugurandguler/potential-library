#!/usr/bin/env python3
"""
NIST's stacking fault energies for the published potentials we compare against.

This is the external control for the one measurement in the project that comes
back with the wrong SIGN.  Our copper gives an intrinsic fault of -63 mJ/m^2
where the experiment is +45, and the argument that this is the potential rather
than the machinery rests on the published potentials going through the same
code and landing where they should.  They do - five of five come back positive,
Mishin's at 44.4 - but "where they should" was so far only the literature value
for one of them.  NIST computes the same quantity for everything it hosts.

Getting at it took two goes.  The per-property CSV endpoint that serves the
surface energies and the phonon tables has no stacking-fault member: the data
is stored per fault as

    stackingfault.{El}.{prototype}--{plane}sf.{uuid}.json

and the uuid is not derivable from anything we hold.  The numbers are however
rendered into the entry page itself, in a table of class `datatable-sf`, and
the plane each table belongs to is recoverable from the plot filename that
precedes it.  Plain urllib is refused there with a 403 while the CSV endpoint
is not, so the page is fetched with a browser user agent and the `.html`
suffix, which the bare entry path does not accept.

One thing to keep in mind when reading the agreement: NIST relaxes a slab with
a free surface and slices it, where this project shifts a fully periodic cell
and tilts the box to match.  Different geometries for the same quantity, which
makes the comparison worth more than a repeat of our own arithmetic would be.

    python fetch_nist_sf.py
    python fetch_nist_sf.py --only Cu,Al
"""
import io
import json
import os
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/124.0 Safari/537.36"}
ENTRY = "https://www.ctcms.nist.gov/potentials/entry"
TAG = re.compile(r"<[^>]+>")
#  the plot filename carries prototype and plane; the table that follows it
#  carries the numbers
PLOT = re.compile(r"stackingfault\.(\w+)\.([A-Za-z0-9'\-]+?)--(\w+)sf\.")
TABLE = re.compile(r'class="[^"]*datatable-sf[^"]*"[^>]*>(.*?)</table>', re.S)
ROW = re.compile(r"<tr>\s*<td>(.*?)</td>\s*<td>(.*?)</td>\s*</tr>", re.S)


def clean(s):
    return TAG.sub("", s).replace("&#964;", "tau").strip()


def parse(html, el):
    """{prototype: {plane: {label: value}}} from one entry page"""
    out = {}
    #  Each table is attributed by the plot filenames it sits with, because
    #  keying the tables by order alone would put the (100) numbers under
    #  (111) for any potential whose planes come in a different order - the
    #  same trap the surface fetch fell into with prototypes.
    #
    #  The direction matters and is not the obvious one: on the page the table
    #  comes FIRST and its plots follow it, so the plane is in the NEXT
    #  filename, not the previous one.  Reading backwards gave Mishin copper a
    #  (100) table holding 44.39, which is the (111) intrinsic fault - a
    #  number that is right, under a label that is wrong, and (111) missing
    #  altogether.  It showed up only because (100) has no intrinsic fault to
    #  report at all.
    marks = [(m.start(), m.group(2), m.group(3))
             for m in PLOT.finditer(html) if m.group(1) == el]
    for m in TABLE.finditer(html):
        nxt = [x for x in marks if x[0] > m.start()]
        if not nxt:
            continue
        _, proto, plane = nxt[0]
        vals = {}
        for lab, val in ROW.findall(m.group(1)):
            try:
                vals[clean(lab)] = float(clean(val))
            except ValueError:
                pass
        if vals:
            out.setdefault(proto, {})[plane] = vals
    return out


def one(job):
    el, fn, impl, potid = job
    url = f"{ENTRY}/{potid}/{impl}.html"
    try:
        req = urllib.request.Request(url, headers=UA)
        html = urllib.request.urlopen(req, timeout=60).read()
        return el, fn, parse(html.decode("utf-8", "replace"), el)
    except Exception as ex:
        return el, fn, {"error": f"{type(ex).__name__}"}


def main():
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))

    p = os.path.join(HERE, "nist_props.json")
    props = json.load(open(p))

    jobs = []
    for key, r in sorted(props.items()):
        el, fn = key.split("|", 1)
        if only and el not in only:
            continue
        if "impl" not in r or "potid" not in r:
            continue
        jobs.append((el, fn, r["impl"], r["potid"]))
    print(f"{len(jobs)} records to query")

    n = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        for el, fn, sf in ex.map(one, jobs):
            key = f"{el}|{fn}"
            if "error" in sf:
                print(f"  {el:3s} {fn[:34]:34s} {sf['error']}")
                continue
            if not sf:
                continue
            props[key]["stacking"] = sf
            n += 1
            flat = [f"{pr.split('--')[0]}/{pl}"
                    for pr, d in sf.items() for pl in d]
            print(f"  {el:3s} {fn[:34]:34s} {len(flat):2d} faults "
                  f"{','.join(flat[:4])}")

    tmp = p + ".tmp"
    json.dump(props, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, p)
    print(f"\n{n} records carry a stacking fault -> {p}")


if __name__ == "__main__":
    main()
