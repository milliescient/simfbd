#!/usr/bin/env python3
"""Describe an exploratory batch: realised rates, record size, and how often the
classical estimators are defined, overall and by sampling regime.

    explore_report.py <batch-dir> <rb-with-estimators> [complete|firstlast]

max_age is taken as the oldest reported bin edge in each record, which is what an
analyst has. Using the generating origin would hand the classical methods information
the Bayesian model has to estimate.
"""
import glob, os, re, statistics as st, subprocess, sys, tempfile, collections

batch, rb = sys.argv[1], sys.argv[2]
reporting = sys.argv[3] if len(sys.argv) > 3 else "complete"
# the batch's settings as it was simulated, whichever kind of config it has
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import manifest
settings = manifest.flat(batch)
def get(k, d=None):
    v = settings.get(k)
    return d if v in (None, "none") else float(v)
nint, w = int(get("NINTERVALS")), get("INTERVAL_WIDTH")
breaks = [w * i for i in range(1, nint)]

tv = [l.split("\t") for l in open(os.path.join(batch, "true_vals.tsv")).read().strip().split("\n")]
hdr, data = tv[0], tv[1:]
pidx = [i for i, h in enumerate(hdr) if h.startswith("psi")]
lidx = [i for i, h in enumerate(hdr) if h.startswith("lambda")]

recs = []
for n, row in enumerate(data, 1):
    f = os.path.join(batch, f"specimens_{reporting}", f"taxa_{n}.tsv")
    if not os.path.exists(f): continue
    lines = open(f).read().strip().split("\n")[1:]
    if not lines: continue
    cols = [l.split("\t") for l in lines]
    recs.append(dict(
        n=n, file=f,
        psi=sum(float(row[i]) for i in pidx) / len(pidx),
        lam=sum(float(row[i]) for i in lidx) / len(lidx),
        taxa=len({c[0] for c in cols}), occ=len(cols),
        maxage=max(float(c[2]) for c in cols)))

# A record whose oldest bin never reaches the oldest rate shift cannot define that
# interval, and the estimator throws rather than returning NaN, which would abort the
# whole script. Count those separately instead of feeding them in.
spans = [r for r in recs if r["maxage"] > breaks[-1]]
short = len(recs) - len(spans)

METHODS = ["fnPerCapitaRates", "fnBoundaryCrosserRates", "fnThreeTimerRates"]
L = ["timeline <- v(" + ", ".join(f"{b:.6f}" for b in breaks) + ")"]
for r in spans:
    L.append(f'taxa <- readTaxonData("{os.path.abspath(r["file"])}")')
    for meth in METHODS:
        L.append(f'm <- {meth}(taxa=taxa, timeline=timeline, present=0.0, max_age={r["maxage"]+1e-6:.6f})')
        for c in range(1, nint + 1):
            L.append(f'print("V\t{r["n"]}\t{meth}\t{c}\t" + m[1][{c}])')
L.append("q()")
with tempfile.NamedTemporaryFile("w", suffix=".Rev", delete=False) as fh:
    fh.write("\n".join(L) + "\n"); script = fh.name
out = subprocess.run([rb, script], capture_output=True, text=True).stdout

defined = collections.defaultdict(dict)
for line in out.splitlines():
    p = line.strip().split("\t")
    if len(p) == 5 and p[0] == "V":
        rep, meth, col, val = int(p[1]), p[2], int(p[3]), p[4]
        defined[(rep, meth)][col] = val.lower() not in ("nan", "-nan", "inf", "-inf")

interior = list(range(2, nint))          # first and last are structurally undefined for bc/tt
def frac(sub, meth):
    vals = [sum(defined[(r["n"], meth)].get(c, False) for c in interior) / len(interior)
            for r in sub if (r["n"], meth) in defined]
    return 100 * st.mean(vals) if vals else 0.0

spans.sort(key=lambda r: r["psi"])
k = max(1, len(spans) // 3)
strata = [("low psi", spans[:k]), ("mid psi", spans[k:2*k]), ("high psi", spans[2*k:])]

print(f"\n  == {os.path.basename(batch)}  ({len(recs)} reps, {nint} intervals, "
      f"MIN_TAXA={get('MIN_TAXA',0):.0f} MAX_TAXA={get('MAX_TAXA',0):.0f}, {reporting})")
print(f"  {short}/{len(recs)} records do not reach the oldest rate shift and are excluded")
ps = sorted(r["psi"] for r in spans)
print(f"  psi  p10 {ps[len(ps)//10]:.3f}  median {ps[len(ps)//2]:.3f}  p90 {ps[9*len(ps)//10]:.3f}")
print(f"  {'stratum':<10}{'psi':>7}{'taxa':>7}{'occ':>7}   {'per-cap':>8}{'bound-x':>9}{'3-timer':>9}")
for lab, sub in strata:
    if not sub: continue
    print(f"  {lab:<10}{st.median(r['psi'] for r in sub):>7.3f}{st.median(r['taxa'] for r in sub):>7.0f}"
          f"{st.median(r['occ'] for r in sub):>7.0f}   "
          f"{frac(sub,'fnPerCapitaRates'):>8.0f}{frac(sub,'fnBoundaryCrosserRates'):>9.0f}"
          f"{frac(sub,'fnThreeTimerRates'):>9.0f}")
print(f"  {'ALL':<10}{st.median(r['psi'] for r in spans):>7.3f}{st.median(r['taxa'] for r in spans):>7.0f}"
      f"{st.median(r['occ'] for r in spans):>7.0f}   "
      f"{frac(spans,'fnPerCapitaRates'):>8.0f}{frac(spans,'fnBoundaryCrosserRates'):>9.0f}"
      f"{frac(spans,'fnThreeTimerRates'):>9.0f}")
print("  percent of interior intervals with a defined rate")
