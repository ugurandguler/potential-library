# -*- coding: utf-8 -*-
"""Look at the finite-temperature panel, rather than trust that it parses.

The syntax gate says the file is JavaScript; it does not say the curve landed
on the axis, that the dotted stretch is where Gamma is, or that the key does
not sit on top of a branch.  Those are answered by looking.

A probe copy of the page is made per element with a script appended that picks
that element, re-renders, and lifts the panel out to the top-left corner so a
viewport screenshot catches it - the page is one very long column and headless
Chrome photographs the top of it.  The probe copies are scratch; the real page
is never touched.

    python shot_ft.py Cu Nb Ti        # one png per element
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
#  the page beside this script, so the tool travels with whichever tree it
#  sits in rather than pointing at the one it was written in
PAGE = os.path.join(HERE, "potential.html")


def chrome():
    """the first browser on this machine that can take a headless screenshot"""
    for p in (r"C:\Program Files\Google\Chrome\Application\chrome.exe",
              r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
              r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
              r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
              "/usr/bin/google-chrome", "/usr/bin/chromium"):
        if os.path.exists(p):
            return p
    raise SystemExit("no Chrome or Edge found - edit chrome()")

#  cur is the module-level element key and render() redraws from it, so the
#  panel can be reached without clicking anything.  Two frames of delay: one
#  for the initial render the page does for itself, one for ours.
JS = """
<script>
(function(){
  const EL = "%s";
  function lift(){
    const cv = document.getElementById("ftdisp");
    if(!cv){
      //  an EXCLUDED element has no canvas, only the heading and its
      //  sentence - lift those instead, which is the whole point of
      //  looking at one
      const h = [].slice.call(document.querySelectorAll("h3")).filter(
        n => /dispersion at the temperature/i.test(n.textContent))[0];
      if(!h){ document.title = "NO-PANEL"; return; }
      const box = document.createElement("div");
      box.style.cssText = "position:fixed;left:0;top:0;z-index:99999;"
        + "width:820px;padding:24px;background:"
        + getComputedStyle(document.body).backgroundColor;
      document.body.appendChild(box);
      box.appendChild(h); if(box.nextSibling){} 
      const p = document.querySelector("p.plotnote");
      [].slice.call(document.querySelectorAll("p.plotnote")).forEach(function(n){
        if(/outside this study/.test(n.textContent)) box.appendChild(n); });
      document.title = "OK-" + EL; return;
    }
    const w = cv.clientWidth;
    const box = document.createElement("div");
    box.style.cssText = "position:fixed;left:0;top:0;z-index:99999;"
      + "width:" + (w+40) + "px;padding:20px;"
      + "background:" + getComputedStyle(document.body).backgroundColor;
    const take = [cv.previousElementSibling, cv,
                  cv.nextElementSibling, cv.nextElementSibling
                    && cv.nextElementSibling.nextElementSibling];
    document.body.appendChild(box);
    take.forEach(n => { if(n) box.appendChild(n); });
    document.title = "OK-" + EL;
  }
  addEventListener("load", function(){
    setTimeout(function(){
      try { cur = EL; render(); } catch(e) { document.title = "ERR " + e; }
      setTimeout(lift, 400);
    }, 400);
  });
})();
</script>
"""


def main(els):
    exe = chrome()
    src = open(PAGE, encoding="utf-8").read()
    #  the page ends at </script> with no </body> - the browser closes both
    #  for itself - so the probe is appended rather than inserted
    assert src.rstrip().endswith("</script>")
    #  Probe copies and the browser profile go to a scratch directory, never
    #  beside the script: a seven-megabyte probe_Cu.html written into the tree
    #  is the kind of file that gets committed by a wide `git add`.  The PNGs
    #  land in the working directory, because looking at them is the point.
    tmp = tempfile.mkdtemp(prefix="shot_ft_")
    try:
        for el in els:
            probe = os.path.join(tmp, "probe_%s.html" % el)
            with open(probe, "w", encoding="utf-8") as fh:
                fh.write(src + JS % el)
            png = os.path.abspath("ft_%s.png" % el)
            if os.path.exists(png):
                os.remove(png)
            subprocess.run([exe, "--headless=new", "--disable-gpu",
                            "--hide-scrollbars",
                            "--force-device-scale-factor=1",
                            "--user-data-dir=" + os.path.join(tmp, "prof"),
                            "--virtual-time-budget=9000",
                            "--window-size=1200,900",
                            "--screenshot=" + png,
                            "file:///" + probe.replace("\\", "/")],
                           capture_output=True, timeout=180)
            print("%-3s %s" % (el, ("%s  %d bayt" % (png,
                                                     os.path.getsize(png)))
                               if os.path.exists(png) else "PNG YOK"))
            os.remove(probe)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["Cu"])
