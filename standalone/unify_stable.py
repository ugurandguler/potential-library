#!/usr/bin/env python3
"""
One meaning for `stable`: no imaginary mode on the mesh AND the symmetry path.

The field used to be written by two rules.  add_taper_overlay.py and
refresh_measured.py sampled the 8^3 and 9^3 Monkhorst-Pack meshes only; the
path screens (add_dynstab.py, add_dynstab_ug.py, add_dynstab_arms.py) sample
those meshes AND the Setyawan-Curtarolo path at a physical -1 cm^-1, and write
the result to `dyn.stable`.  A half-shifted mesh never lands on a symmetry
line, so a mode confined to one is invisible to it however fine it is made -
the iridium mistake, fixed for the candidates in add_candidates.py and never
carried back.  On 2026-09-13, 16 of 204 records said `stable: true` beside
`dyn.stable: false`: tap Cr Cs Fe K Mo Nb Rb W, tap_ug Cr Cs Fe Nb Rb, ug Al Na,
and iridium's top-level record.

The page's red "Dynamically unstable" box and the potential-file headers
already read `dyn`, so they were right; the arm labels read `stable` and said
nothing.  This copies the stronger verdict into `stable` wherever a `dyn` block
exists.  It computes nothing, so it cannot drift from the screens it reads, and
it is idempotent - the producers now write the same thing themselves, and this
remains the check that they did.

A record with no `dyn` block keeps its own `stable`: there is no stronger
verdict to copy.

    python unify_stable.py            # write
    python unify_stable.py --check    # report only, exit 1 if anything differs
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = ("hard", "tap", "tap_ug", "ug", "rc", "rc_ug", "tap_force", "hard_disp",
        "disp", "nudge", "stack", "force")


def records(lib):
    for el in sorted(lib):
        v = lib[el]
        if not isinstance(v, dict):
            continue
        yield el, "top", v
        for a in ARMS:
            r = v.get(a)
            if isinstance(r, dict):
                yield el, a, r


def main():
    check = "--check" in sys.argv
    path = os.path.join(HERE, "library.json")
    lib = json.load(open(path))
    n_dyn, changed = 0, []
    for el, arm, r in records(lib):
        dyn = r.get("dyn")
        if not isinstance(dyn, dict) or "stable" not in dyn:
            continue
        n_dyn += 1
        want = bool(dyn["stable"])
        if r.get("stable") != want:
            changed.append((el, arm, r.get("stable"), want,
                            dyn.get("min_mesh_cm1"), dyn.get("min_path_cm1")))
            r["stable"] = want
    print(f"{n_dyn} records carry a dyn block; {len(changed)} disagreed")
    for el, arm, was, now, mm, mp in changed:
        print(f"   {el:3s} {arm:9s} {str(was):>5s} -> {str(now):5s}"
              f"   mesh {mm}  path {mp}")
    if check:
        return 1 if changed else 0
    if changed:
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(lib, fh, indent=1, sort_keys=True, default=str)
        os.replace(tmp, path)
        print(f"written -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
