#!/usr/bin/env python3
"""
Run refine_R.py over the cubic metals, one process per element.

Sequentially this is about 44 minutes per element at rcut3 = 1.50 - the grid is
7920 points and each of the 249 simplex starts costs several hundred more, and
the 1.50 cutoff carries 91 triplets per atom for bcc against 28 at 1.12.  Twenty
three of those in one process is seventeen hours, and refine_R.py writes its
output only at the end, so a shared-file run also has nothing to show for itself
until it finishes.

One element per process is what refine_R.py is written for: each writes its own
R_floor_<el>.json, so progress is visible and a crash costs one element.

    python run_floors.py            # all cubic metals
    python run_floors.py Nb V       # a subset
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
#  one core left for everything else, and OMP pinned to 1: these are 3x3 and 6x6
#  matrices, so a threaded BLAS only oversubscribes.  On a cluster node set
#  FLOOR_NPROC to the cores allocated - barbun gives 40, which runs all
#  twenty three at once.
NPROC = int(os.environ.get("FLOOR_NPROC", "7"))


def main(els):
    if not els:
        lib = json.load(open(os.path.join(HERE, "library.json")))
        els = sorted(e for e in lib if lib[e]["struct"] != "hcp")
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1")
    print(f"{len(els)} elements, {NPROC} in parallel: {' '.join(els)}\n", flush=True)
    t0 = time.time()
    queue, running, done = list(els), [], []
    while queue or running:
        while queue and len(running) < NPROC:
            el = queue.pop(0)
            log = os.path.join(HERE, f"floor_{el}.log")
            p = subprocess.Popen(
                [sys.executable, "-u", os.path.join(HERE, "refine_R.py"), el],
                cwd=HERE, env=env, stdout=open(log, "w"),
                stderr=subprocess.STDOUT)
            running.append((p, el, time.time()))
        time.sleep(5)
        for rec in list(running):
            p, el, ts = rec
            if p.poll() is not None:
                running.remove(rec)
                done.append(el)
                print(f"[{len(done):2d}/{len(els)}] {el:3s} "
                      f"{(time.time()-ts)/60:5.1f} dk  "
                      f"({(time.time()-t0)/60:.0f} min in total)", flush=True)
    print(f"\nbitti, {(time.time()-t0)/60:.0f} dakika")


if __name__ == "__main__":
    main(sys.argv[1:])
