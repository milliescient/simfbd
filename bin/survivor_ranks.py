#!/usr/bin/env python3
"""Rank the rates in both survivor arms, paired on the replicate.

Ranks, not point estimates: (est-true)/true divides by a lognormal draw that can be near zero
and reads as a large bias even when the posterior is centred.

A log carries ~140 columns and this needs nine, so it reads only those. Parsing the whole table
costs 15x more and is the entire runtime.
"""
import os
import sys
import numpy as np

BATCH = os.environ.get("BATCH", "sims/survivors")
ROOT = "/research/phyloworks/simfbd"
ESS_MIN = 100.0
BURN = 0.25
RATES = [f"{r}[{i}]" for r in ("lambda", "mu", "psi") for i in (1, 2, 3)]


def ess(x):
    """Single-chain effective sample size, Geyer initial positive sequence."""
    n = len(x)
    x = x - x.mean()
    v = np.var(x)
    if v == 0:
        return 0.0
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conjugate(f))[:n].real
    ac /= ac[0]
    # sum paired autocorrelations while they stay positive
    t, s = 1, 0.0
    while t + 1 < n:
        p = ac[t] + ac[t + 1]
        if p <= 0:
            break
        s += p
        t += 2
    return n / (1.0 + 2.0 * s)


def truth_columns():
    """Batches with no true_vals_order key predate the fix and are oldest-first."""
    man = os.path.join(ROOT, BATCH, "manifest.tsv")
    order = "oldest"
    if os.path.exists(man):
        if any(l.startswith("true_vals_order") for l in open(man)):
            order = "youngest"
    idx = [3, 2, 1] if order == "oldest" else [1, 2, 3]
    return [f"{r}{i}" for r in ("lambda", "mu", "psi") for i in idx], order


def read_arm(dirpath, tv, tcols):
    out = {}
    if not os.path.isdir(dirpath):
        return out
    for fn in sorted(os.listdir(dirpath)):
        if not (fn.startswith("rep_") and fn.endswith(".log")):
            continue
        rep = int(fn[len("rep_"):-len(".log")])
        path = os.path.join(dirpath, fn)
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        try:
            cols = [hdr.index(c) for c in RATES]
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
        true = np.array([tv[c][rep - 1] for c in tcols])
        out[rep] = (a < true).mean(axis=0)
    return out


def main():
    tvpath = os.path.join(ROOT, BATCH, "true_vals.tsv")
    tv = np.genfromtxt(tvpath, delimiter="\t", names=True)
    tcols, order = truth_columns()

    a = read_arm(os.path.join(ROOT, BATCH, "runs/survivors_true/output"), tv, tcols)
    b = read_arm(os.path.join(ROOT, BATCH, "runs/survivors_false/output"), tv, tcols)
    pair = sorted(set(a) & set(b))
    n = len(pair)
    print(f"\n{BATCH}  truth {order}-first  paired {n} (admitted {len(a)}, forbidden {len(b)})")
    if n < 2:
        return
    A = np.array([a[r] for r in pair])
    B = np.array([b[r] for r in pair])
    se = 1.0 / np.sqrt(12 * n)
    print(f"mean rank, 0.5 is calibrated, SE {se:.3f}\n")
    print(f"{'':10} {'admitted':>18} {'forbidden':>18} {'paired diff':>16}")
    for j, name in enumerate(RATES):
        ra, rb = A[:, j].mean(), B[:, j].mean()
        d = A[:, j] - B[:, j]
        sd = d.std(ddof=1) / np.sqrt(n)
        zd = d.mean() / sd if sd > 0 else 0.0
        print(f"{name:10} {ra:8.3f} (z{(ra-0.5)/se:+5.1f}) {rb:8.3f} (z{(rb-0.5)/se:+5.1f})"
              f" {d.mean():+7.3f} (z{zd:+5.1f})")
    print("\nrank below 0.5 means the posterior sits above the truth")


if __name__ == "__main__":
    main()
