#!/usr/bin/env python3
"""How often is each classical estimator defined, per interval, for a batch?

Simulation only, no inference. Usage:
    trad_coverage.py <batch-dir> <rb-with-estimators> [complete|firstlast]
"""
import glob, os, re, subprocess, sys, tempfile

batch, rb = sys.argv[1], sys.argv[2]
reporting = sys.argv[3] if len(sys.argv) > 3 else "complete"

# the batch's settings as it was simulated, whichever kind of config it has
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import manifest
settings = manifest.flat(batch)
def get(k, d=None):
    v = settings.get(k)
    return d if v in (None, "none") else float(v)
nint, w, agemax = int(get("NINTERVALS")), get("INTERVAL_WIDTH"), get("AGE_MAX")
breaks = [w * i for i in range(1, nint)]
files = sorted(glob.glob(os.path.join(batch, f"specimens_{reporting}", "*.tsv")))

METHODS = ["fnPerCapitaRates", "fnBoundaryCrosserRates", "fnThreeTimerRates"]
L = ["timeline <- v(" + ", ".join(f"{b:.6f}" for b in breaks) + ")"]
for f in files:
    L.append(f'taxa <- readTaxonData("{os.path.abspath(f)}")')
    for meth in METHODS:
        L.append(f'm <- {meth}(taxa=taxa, timeline=timeline, present=0.0, max_age={agemax:.6f})')
        cells = " + ".join(f'm[{r}][{c}]' for r in (1, 2) for c in range(1, nint + 1))
        # one value per line keeps parsing trivial
        for r in (1, 2):
            for c in range(1, nint + 1):
                L.append(f'print("V\t{meth}\t{r}\t{c}\t" + m[{r}][{c}])')
L.append("q()")

with tempfile.NamedTemporaryFile("w", suffix=".Rev", delete=False) as fh:
    fh.write("\n".join(L) + "\n"); script = fh.name
out = subprocess.run([rb, script], capture_output=True, text=True).stdout

import collections
ok = collections.Counter(); tot = collections.Counter()
for line in out.splitlines():
    p = line.strip().split("\t")
    if len(p) == 5 and p[0] == "V":
        meth, row, col, val = p[1], int(p[2]), int(p[3]), p[4]
        tot[(meth, row, col)] += 1
        if val.lower() not in ("nan", "-nan", "inf", "-inf"):
            ok[(meth, row, col)] += 1

print(f"\n  {os.path.basename(batch)}  ({len(files)} reps, {nint} intervals, reporting={reporting})")
print(f"  {'method':<22}{''.join(f'{i:>5}' for i in range(1, nint + 1))}")
for meth in METHODS:
    pct = [round(100 * ok[(meth,1,c)] / tot[(meth,1,c)]) if tot[(meth,1,c)] else 0
           for c in range(1, nint + 1)]
    print(f"  {meth[2:]:<22}" + "".join(f"{p:>5}" for p in pct))
print("  percent of reps with a defined speciation rate, interval 1 = youngest")
