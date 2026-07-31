#!/usr/bin/env python3
"""
David's SBC diagnostics (rank histogram, ECDF-difference w/ 95% simultaneous band,
coverage, posterior-median bias) for the simfbd pipeline output. Generalised to N
skyline intervals and any of lambda/mu/psi. Runs without the R `SBC` package
(validated against it conceptually; the R version is in sbc_diagnostics.R).

Usage:
  python sbc_diagnostics.py [DATADIR] [OUTDIR] [mapping]
    DATADIR  default sbc/data   (expects DATADIR/true_vals.tsv: rep, lambda1..ni, mu1..ni, psi1..ni, origin, nsamp)
    OUTDIR   default sbc/output (expects OUTDIR/fbd_<rep>.log)
    mapping  'direct' (default; simfbd youngest-first == rb) | 'reversed' (paleobuddy oldest-first)
"""
import os, sys, numpy as np
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,"data")
OUTDIR  = sys.argv[2] if len(sys.argv)>2 else os.path.join(HERE,"output")
MAPPING = sys.argv[3] if len(sys.argv)>3 else "direct"
RATES=["lambda","mu","psi"]; QUANTS=np.round(np.arange(0.05,0.96,0.05),2)
L=100; BURN=0.10

tv = np.genfromtxt(os.path.join(DATADIR,"true_vals.tsv"), delimiter="\t", names=True)
names = tv.dtype.names
ni = sum(1 for n in names if n.startswith("lambda"))
def tmap(k):  # true interval k (1..ni) -> rb interval index
    return (ni+1-k) if MAPPING=="reversed" else k
print(f"SBC diagnostics: ni={ni} intervals, mapping={MAPPING}, reps in true_vals={tv.shape[0]}")

def thin(s):
    idx=np.linspace(0,len(s)-1,L).round().astype(int); return s[idx]
def field(d,r,k):
    for cand in (f"{r}{k}", f"{r}_{k}_", f"{r}{k}_"):
        if cand in d.dtype.names: return d[cand]
    raise KeyError(f"{r}{k} not in {d.dtype.names}")

cov={f"{r}{k}":np.zeros(len(QUANTS)) for r in RATES for k in range(1,ni+1)}
ranks={f"{r}{k}":[] for r in RATES for k in range(1,ni+1)}
bias={f"{r}{k}":[] for r in RATES for k in range(1,ni+1)}
used=0
for rep in range(1, tv.shape[0]+1):
    p=os.path.join(OUTDIR, f"fbd_{rep}.log")
    if not os.path.exists(p): continue
    d=np.genfromtxt(p, delimiter="\t", names=True)
    if d.ndim==0 or d.shape[0]<20: continue
    s=int(d.shape[0]*BURN); used+=1
    for r in RATES:
        for k in range(1,ni+1):
            samp=field(d,r,tmap(k))[s:]; truev=tv[f"{r}{k}"][rep-1]
            bias[f"{r}{k}"].append((np.median(samp)-truev)/truev)
            ranks[f"{r}{k}"].append(int(np.sum(thin(samp)<truev)))
            for qi,q in enumerate(QUANTS):
                lo=np.quantile(samp,0.5-q/2); hi=np.quantile(samp,0.5+q/2)
                if lo<truev<hi: cov[f"{r}{k}"][qi]+=1

GRID=np.linspace(0,1,201)
def ecdf_diff(u): return np.searchsorted(np.sort(u),GRID,side="right")/len(u) - GRID
def band(n,B=3000):
    rng=np.random.default_rng(7); m=np.empty(B)
    for b in range(B): m[b]=np.max(np.abs(ecdf_diff(rng.integers(0,L+1,n)/L)))
    return np.quantile(m,0.95)
def unifp(rk):
    obs,_=np.histogram(np.array(rk),bins=np.linspace(0,L,11)); exp=np.full(10,len(rk)/10)
    return stats.chi2.sf(np.sum((obs-exp)**2/exp),df=9)

bnd=band(used)
print(f"\nreps used={used}; 95% simultaneous ECDF band = +/-{bnd:.3f}\n")
print(f"{'param':8s} {'cov_err':>8s} {'rankSBC_p':>10s} {'med_bias':>9s}")
nfail=0
for r in RATES:
    for k in range(1,ni+1):
        key=f"{r}{k}"; c=cov[key]/used; e=np.mean(np.abs(c-QUANTS)); p=unifp(ranks[key]); b=np.mean(bias[key])
        fail = p<1e-3; nfail+=fail
        print(f"{key:8s} {e:8.3f} {p:10.1e} {b:+9.3f}{'  FAIL' if fail else ''}")
print(f"\nmean coverage error = {np.mean([np.mean(np.abs(cov[f'{r}{k}']/used-QUANTS)) for r in RATES for k in range(1,ni+1)]):.4f}")
print(f"rank-SBC failures (p<1e-3): {nfail}/{3*ni}")

# ---- plots: rank hist + ECDF-diff (David's plot_rank_hist / plot_ecdf_diff) ----
fig,axes=plt.subplots(2*3, ni, figsize=(3.2*ni, 16), squeeze=False)
for ri,r in enumerate(RATES):
    for k in range(1,ni+1):
        key=f"{r}{k}"; p=unifp(ranks[key])
        ax=axes[ri][k-1]; ax.hist(ranks[key],bins=20,color="tab:blue",edgecolor="white")
        ax.axhline(used/20,color="grey",ls="--",lw=1); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"rank {key} (p={p:.0e})",fontsize=8,color="red" if p<1e-3 else "black")
        u=np.array(ranks[key])/L; curve=ecdf_diff(u); md=np.max(np.abs(curve)); out=md>bnd
        ax2=axes[3+ri][k-1]; ax2.axhspan(-bnd,bnd,color="0.85"); ax2.axhline(0,color="grey",lw=.7)
        ax2.plot(GRID,curve,color="tab:red" if out else "tab:blue"); ax2.set_ylim(-0.15,0.15)
        ax2.set_title(f"ecdf-diff {key} ({'FAIL' if out else 'ok'})",fontsize=8,color="red" if out else "black")
fig.suptitle(f"simfbd SBC diagnostics ({MAPPING} mapping, n={used})",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.98])
out=os.path.join(OUTDIR,"sbc_diagnostics.png"); fig.savefig(out,dpi=110)
print(f"\nsaved {out}")
