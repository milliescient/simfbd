#!/usr/bin/env python3
"""Rank SBC for one analysis: the fraction of posterior draws below the truth, per rate.

Ranks, not point estimates: (est-true)/true divides by a lognormal draw that can be near
zero and reads as a large bias even when the posterior is centred.

Usage:  bin/ranks.py sims/<batch>/runs/<name>
"""
import os
import sys
import numpy as np

ESS_MIN = 100.0
BURN = 0.25


def ess(x):
    """Single-chain effective sample size, Geyer initial positive sequence."""
    n = len(x)
    x = x - x.mean()
    if np.var(x) == 0:
        return 0.0
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conjugate(f))[:n].real
    ac /= ac[0]
    t, s = 1, 0.0
    while t + 1 < n:
        p = ac[t] + ac[t + 1]
        if p <= 0:
            break
        s += p
        t += 2
    return n / (1.0 + 2.0 * s)


def main():
    run = sys.argv[1].rstrip("/")
    batch = os.path.dirname(os.path.dirname(run))
    tv = np.genfromtxt(os.path.join(batch, "true_vals.tsv"), delimiter="\t", names=True)
    ni = sum(1 for c in tv.dtype.names if c.startswith("lambda"))

    # batches with no true_vals_order key predate the fix and are oldest-first
    man = os.path.join(batch, "manifest.tsv")
    youngest = os.path.exists(man) and any(l.startswith("true_vals_order") for l in open(man))
    idx = list(range(1, ni + 1)) if youngest else list(range(ni, 0, -1))
    rates = [f"{r}[{i}]" for r in ("lambda", "mu", "psi") for i in range(1, ni + 1)]
    tcols = [f"{r}{i}" for r in ("lambda", "mu", "psi") for i in idx]

    out = os.path.join(run, "output")
    ranks, kept, seen = [], 0, 0
    for fn in sorted(os.listdir(out)):
        if not fn.endswith(".log") or "_" not in fn:
            continue
        seen += 1
        path = os.path.join(out, fn)
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        try:
            cols = [hdr.index(c) for c in rates]
        except ValueError:
            continue
        try:
            a = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=cols, ndmin=2)
        except Exception:
            continue
        if a.shape[0] < 500:
            continue
        a = a[int(a.shape[0] * BURN):]
        if min(ess(a[:, j]) for j in range(a.shape[1])) < ESS_MIN:
            continue
        rep = int(fn.rsplit("_", 1)[1][:-4])
        ranks.append((a < np.array([tv[c][rep - 1] for c in tcols])).mean(axis=0))
        kept += 1

    print(f"\n{run}  truth {'youngest' if youngest else 'oldest'}-first"
          f"  kept {kept} of {seen}")
    if kept < 2:
        return
    R = np.array(ranks)
    se = 1.0 / np.sqrt(12 * kept)
    print(f"mean rank, 0.5 is calibrated, SE {se:.3f}\n")
    for j, name in enumerate(rates):
        m = R[:, j].mean()
        print(f"{name:10} {m:7.3f}  z {(m - 0.5) / se:+6.2f}")
    print("\nrank below 0.5 means the posterior sits above the truth")


if __name__ == "__main__":
    main()
