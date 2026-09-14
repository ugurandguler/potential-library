#!/usr/bin/env python3
"""
Rebuild the 0 K analytic layer downstream of the fits, in the one order
that works - and NOT the measured layer, which is the larger half.

    python refresh.py                       # rebuild from the current fit.json
    python refresh.py runs/2026-08-01_cluster # merge that run in first

The elastic constants, the mechanical analysis, the phonons, the
thermodynamics and the DFT overlays are recomputed from the merged parameters.
Nothing is carried over except the fetched MP reference data, which is a
property of the element and not of the fit.

**WHAT THIS DOES NOT REBUILD, and it is more than what it does.** The tree
holds twenty-eight producers and this chain runs eleven.  The switched arms
`tap` and `tap_ug` - the two the library actually recommends for molecular
dynamics - the re-cut and force-matched arms, the elastic constants against
temperature, the surface, stacking, ground-state and expansion measurements,
the MD screen, the nudge test and every published-potential baseline come from
the other seventeen, and most of those read cluster output that the published
tree does not carry.

Since `build_library.py` builds its dictionary from scratch and keeps only
`mp`, running this over a populated library DELETES all of that.  Measured:
the published page is 7.0 MB and what this chain regenerates from `fit.json`
alone is 3.2 MB.  So it now refuses when the library holds records it cannot
put back, lists them, and needs `--force` to go ahead anyway.

The order is not arbitrary:

  merge_fits      stability-aware best-of; must run before anything reads
                  fit.json, and re-measures stability itself rather than
                  trusting whatever screen a run happened to use
  build_library   elastic tensor, mechanics.py analysis, dispersion, thermo
  add_dynstab     stability flags on the 8^3 U 9^3 union mesh
  add_mp_overlay  our dispersion re-evaluated at MP's own q-points
  add_mc3d_overlay  the same against Materials Cloud MC3D, which covers
                  Ag Au Cr Mo Nb Ni Pb Pd Rh Ta where MP has nothing
  add_exp_phonon  measured frequencies at X/L (fcc) and H/N (bcc), scored
  add_jarvis_overlay  a third DFT opinion, NIST JARVIS-DFT.  Entries more than
                  25 % away from a reference we already hold are stored but
                  flagged rather than drawn.  **Must come after add_exp_phonon**,
                  which supplies the only reference several elements have: with
                  the two the other way round, iron had nothing to be checked
                  against and its non-spin-polarised curve - 121 cm^-1 against a
                  measured 309 - was drawn as trusted
  add_reachability  C44/C' floor verdict; reads R_floor.json and
                  cprime_region.json, which change only if the form does
  fix_mp_path     MP labels and discontinuities - AFTER add_mp_overlay, which
                  rewrites the phonon record
  add_phonon_curves the measured dispersion curves; add_bvk_curve the
                  reconstructed ones.  Both read refdata_phonon_curves.py and
                  need nothing fetched
  add_dynstab_ug  the same stability screen on the angular arm
  unify_stable    one meaning for `stable` on every record, mesh AND path -
                  AFTER the last stability screen, whose verdict it copies
  add_ld_tap      the switched arm's own 0 K dispersion, which the finite-T
                  panel draws its reference curve from; a clean clone has no
                  `tap` record and the step does nothing there
  recommend_recut which switched set to recommend for MD, per element, by a
                  rule over the warm and cold screens; needs `rc` and `tap`,
                  so it too does nothing on a clean clone
  add_plane_d     the (1 -1 0) polar section.  The three coordinate planes
                  contain no member of <111>, which for a cubic crystal is
                  where the extremum sits
  fix_mark_x      the fractional index of each high-symmetry label - AFTER the
                  overlays, whose q-lists it corrects
  make_gui        the page

WHAT A FRESH CLONE CANNOT REBUILD, and why.  None of these is a fault; each is
data that is not ours to redistribute or is too large to carry, and each step
says so and continues rather than failing:

  the DFT overlays        Materials Project, Materials Cloud MC3D and
                          JARVIS-DFT are fetched data.  Set MP_API_KEY for the
                          first; the other two need their own downloads
  the tapered arm         `tap` and `tap_ug` come from add_taper_overlay.py,
                          which reads the `dense_*.json` search output.  That
                          is gigabytes and is gitignored
  the re-cut candidates   add_candidates.py reads `refit/`, which is not
                          published - and `refit/dyn_candidates.json` is in NO
                          tree at all, so a rebuild also loses the ground-state
                          column.  It says so now rather than skipping quietly.
                          The verdicts already in library.json are unaffected
  the finite-temperature  make_finiteT.py needs the tapered arm above.  The
  panel                   finished table ships as finiteT.json and the page
                          carries it, so nothing is missing from the page -
                          only from a rebuild

The shipped `docs/index.html` is the complete page; a rebuilt one is smaller
and says which comparisons it could not make.

A step that fails stops the run, because a later step reading half-written data
is worse than no rebuild.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

#  Third field: True when the step attaches THIRD-PARTY comparison data that a
#  fresh checkout does not have.  Those datasets are fetched by the fetch_*.py
#  scripts and are not redistributed here, so on a clean clone they are absent -
#  which is not a failure, it is the expected state.  Treating it as one made
#  this script stop at step 3 of 10 for anyone who cloned the repository and
#  followed the README.
#
#  The core steps are NOT optional: they compute from the potential itself, and
#  if one of them fails the page must not be built.
CHAIN = [
    ("build_library.py", "elastic tensor, mechanics, phonons, thermodynamics", False),
    ("add_dynstab.py", "dynamical stability on the union mesh", False),
    ("add_mp_overlay.py", "Materials Project comparison at their q-points", True),
    ("add_mc3d_overlay.py", "Materials Cloud MC3D phonon comparison", True),
    ("add_exp_phonon.py", "measured phonon frequencies, out-of-sample", True),
    ("add_jarvis_overlay.py", "JARVIS-DFT comparison, conflicting entries flagged", True),
    ("add_ug_overlay.py", "UG (angular) results beside MAU, if any exist", True),
    ("add_reachability.py", "is this metal inside the form's reach", False),
    #  Everything below was added after 1.0.0 and was being run by hand, which
    #  meant this script stopped rebuilding what it says it rebuilds: a clean
    #  clone came out with no measured dispersion curves at all.  These five
    #  need nothing but library.json, so they belong in the chain.
    ("add_phonon_curves.py", "measured dispersion curves along the path", False),
    ("add_bvk_curve.py", "the reconstructed curves, from published force "
     "constants", False),
    ("add_dynstab_ug.py", "the stability screen on the angular arm", True),
    ("add_plane_d.py", "the (1 -1 0) polar section, which the three "
     "coordinate planes miss", False),
    ("fix_mp_path.py", "MP high-symmetry labels and path discontinuities", True),
    ("fix_mark_x.py", "the fractional index of each high-symmetry label", True),
    #  Added for 1.3.0.  Both read nothing but library.json, so neither is
    #  optional.  unify_stable must follow every stability screen above;
    #  add_ld_tap must precede the page, or the finite-temperature panel falls
    #  back to drawing no reference curve at all.
    ("unify_stable.py", "one meaning for stable: mesh AND path", False),
    ("add_ld_tap.py", "the switched arm's own 0 K dispersion", False),
    ("recommend_recut.py", "which switched set to recommend for MD", False),
    ("make_gui.py", "potential.html", False),
]


def run(script, *args, optional=False):
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *args],
                       capture_output=True, text=True, cwd=HERE)
    tail = [l for l in p.stdout.strip().splitlines()
            if l.strip() and "Warning" not in l][-1:]
    print(f"    {time.time() - t0:6.1f}s  {tail[0].strip() if tail else 'done'}")
    if not p.returncode:
        return
    #  What may be skipped is "this step's input is not in this checkout" -
    #  either because the script guarded it and stopped itself (SystemExit) or
    #  because it opened a file that is not there (FileNotFoundError).  Some of
    #  these scripts do the first and some the second, and both mean the same
    #  thing.  Anything ELSE that raises is a bug and must still stop the
    #  chain, or a broken overlay would be forgiven by the mechanism that
    #  forgives a missing dataset.
    tb = "Traceback (most recent call last)" in p.stderr
    missing = "FileNotFoundError" in p.stderr
    if optional and (not tb or missing):
        why = (p.stderr.strip() or p.stdout.strip()).splitlines()
        last = why[-1].strip() if why else "no data"
        if missing:
            last = "missing " + last.split("'")[-2].split(os.sep)[-1]                 if "'" in last else last
        print(f"            skipped - {last}")
        return
    print(p.stdout[-2000:])
    print(p.stderr[-2000:])
    raise SystemExit(f"{script} failed ({p.returncode}); stopping so the "
                     f"page is not built from half-written data")



#  Anahtarlar ki bu zincir onlari URETMEZ.  Elle yazilmis bir liste degil:
#  build_library.py'nin uretttigi alanlar ile mevcut kutuphanede duranlar
#  karsilastiriliyor, yani zincire bir adim eklendigi gun liste kendiliginden
#  kisalir.  Sabit bir liste o gun yanlis alarm vermeye baslar.
def would_lose(path):
    """(kayip anahtar -> kac element) - zincirin geri getiremeyecegi her sey

    `build_library.py` her calismada sozlugu sifirdan kurar ve oncekinden
    yalnizca `mp` tasir, cunku geri kalan her sey fit'in bir ozelligi sayilir.
    Bu dogruydu, zincir her seyi ureten tek yer oldugu surece.  Bugun agacta
    yirmi sekiz uretici var ve zincir on birini kosuyor; gerisi kume ciktisi
    ister ve o ciktilar yayimlanan agacta yoktur.  Yani zincir artik bir
    yeniden kurulum degil, KISMI bir yeniden kurulum, ve farki bilmeyen biri
    icin bu sessiz bir silme islemidir.
    """
    import json as _json
    if not os.path.exists(path):
        return {}
    try:
        lib = _json.load(open(path, encoding="utf-8"))
    except ValueError:
        return {}
    #  zincirin yazdigi ust duzey alanlar (build_library + 11 adim)
    #  MEASURED, not guessed from reading the code: this is exactly the key
    #  set a from-scratch chain run produces, taken from one - the public tree
    #  cloned into an empty directory with no library.json, `refresh.py`, and
    #  the keys of what came out.  A hand-written list was wrong on the first
    #  try (it cried wolf over C13, C33 and model_curve, all of which the chain
    #  does make), and a gate that reports what it should not is a gate people
    #  learn to pass.  Re-derive it the same way if the chain gains a step.
    MADE = {"B", "C", "C11", "C12", "C13", "C33", "C44", "Cfull", "Cij_unc",
            "Cp298", "D", "Ecoh", "Ecoh_unc", "P_resid", "S298", "a0",
            "alpha", "alpha3", "at_bound", "c_over_a", "dnn", "dyn", "exp",
            "exp_curve", "exp_phonon", "frozen", "gamma", "ld", "m", "mass",
            "mech", "mech_planes", "model_curve", "ntrip", "r0", "rcut2",
            "rcut3", "reach", "rms", "s3", "struct", "theta_D_exp", "ug",
            #  fetched overlays: build_library keeps `mp`, and the chain's
            #  own steps re-fetch or rewrite these three
            "mp", "jarvis", "mc3d", "min_cm1", "stable", "fc_check"}
    lost = {}
    for el, d in lib.items():
        if not isinstance(d, dict):
            continue
        for k, v in d.items():
            if k in MADE or v is None or v == {} or v == []:
                continue
            lost[k] = lost.get(k, 0) + 1
    return lost


def main(merge_dirs, force=False):
    lost = would_lose(os.path.join(HERE, "library.json"))
    if lost and not force:
        print("DUR.  Mevcut library.json'da bu zincirin GERI GETIREMEYECEGI")
        print("kayitlar var.  build_library.py sozlugu sifirdan kurar ve")
        print("oncekinden yalniz `mp` tasir, yani asagidakiler silinir:")
        print()
        for k, n in sorted(lost.items(), key=lambda x: -x[1]):
            print("    %-22s %d element" % (k, n))
        print()
        print("Bunlari ureten ureticiler bu zincirde YOK - cogu kume ciktisi")
        print("ister ve o ciktilar yayimlanan agacta bulunmaz.  Yayimlanan")
        print("sayfayi bu zincirle yeniden uretemezsiniz; uretilen sayfa")
        print("daha kucuk ve daha az sey soyleyen baska bir sayfadir.")
        print()
        print("Bunun bir de asagi akis sonucu var: kollar gidince")
        print("export_potentials.py 213 yerine 76 dosya yazar.")
        print()
        print("Yalnizca 0 K analitik kismi istiyorsaniz:  refresh.py --force")
        print("Yayimlanan sayfanin tamami zaten docs/index.html'dedir.")
        raise SystemExit(1)
    if merge_dirs:
        print(f"[0/{len(CHAIN)}] merge_fits.py  <- {', '.join(merge_dirs)}")
        run("merge_fits.py", *merge_dirs)
    for i, (script, what, optional) in enumerate(CHAIN, 1):
        print(f"[{i}/{len(CHAIN)}] {script:20s} {what}")
        run(script, optional=optional)
    page = os.path.join(HERE, "potential.html")
    print(f"\n{page}  ({os.path.getsize(page) / 1024:.0f} KB)")


if __name__ == "__main__":
    _a = sys.argv[1:]
    _f = "--force" in _a
    main([x for x in _a if x != "--force"], force=_f)
