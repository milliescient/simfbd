"""Simulate a batch's replicates: its first, or more beside the ones it has.

A batch owns the generating process and every hyperprior. One simulation pass emits all the
reporting models from one record, so the fits under it share a timeline and one
true_vals.tsv, and their comparisons are paired.

Replicates are only ever added. The ones a batch holds are never touched: an extension
simulates what is missing, appends to true_vals.tsv and nums.tsv, and writes its
per-replicate files beside the others. A config that has changed since the earlier replicates
is refused, since adding to them would mix parameters, and so is a simulator that has
changed, unless that is asked for and then the manifest records both.
"""
import os
import subprocess
import sys

from . import config, manifest, tools


def die(msg):
    raise SystemExit(f"error: {msg}")


def warn_if_dirty():
    if tools.dirty():
        sys.stderr.write("warning: simfbd has uncommitted changes, so the commit it records "
                         "will not reproduce this\n")


def replicates(d):
    return int(manifest.flat(d).get("NREPS") or 0)


def extend(d, n, allow_sim_change=False):
    d = os.path.normpath(d)
    cfg = os.path.join(d, "config.yml")
    if not os.path.exists(cfg):
        die(f"no {cfg}")
    if n < 1:
        die("how many must be a positive whole number")
    settings = config.load(cfg, "batch")
    stated = settings.pop("_uses")
    tools.require(stated, cfg)
    cfg_hash = manifest.digest(cfg)
    sim = tools.script("sim.R")
    sim_hash = manifest.digest(sim)

    have = 0
    path, doc = manifest.find(d)
    if path is not None:
        prev = manifest.flat(d).get("config_hash")
        if prev != cfg_hash:
            die(f"{cfg} changed since its replicates were simulated ({prev} -> {cfg_hash}).\n"
                f"       Adding to them or redrawing them would mix parameters. Make a new batch.")
        have = replicates(d)
        was = manifest.flat(d).get("sim_script_hash")
        # same config, different simulator is still a mixed batch, and nothing else would say so
        if have and was != sim_hash and not allow_sim_change:
            die(f"the simulator has changed since replicates up to {have} were simulated "
                f"({was} -> {sim_hash}).\n       Make a new batch, or pass --allow-sim-change "
                f"to add to this one; the manifest records both.")
    last = have + n
    print(f"simulating {d}: reps {have + 1} to {last}, {settings['NINTERVALS']} intervals "
          f"of width {settings['INTERVAL_WIDTH']}", flush=True)
    warn_if_dirty()
    tools.announce()
    env = dict(tools.environment(), **settings, BATCHDIR=d, NREPS=str(last), REPS_HAVE=str(have))
    for stale in ("SIM_TIMEOUT",):
        env.pop(stale, None)
    if subprocess.run([tools.rscript(), str(sim)], env=env).returncode != 0:
        die("simulation failed")
    manifest.write_batch(d, {"first": have + 1, "last": last,
                             "inputs": {"config.yml": cfg_hash, "simfbd": {"sim.R": sim_hash}},
                             "simfbd_commit": tools.version(),
                             "uses": tools.uses(tools.SIMULATES)})
    # the simulator reports what it actually drew with; disagreement means a setting did not reach it
    used = os.path.join(d, "sim_params.tsv")
    if os.path.exists(used):
        want = dict(settings, NREPS=str(last))
        name_of = {key: name for key, name, _ in config.S.BATCH_SETTINGS}
        for line in open(used):
            if "\t" not in line:
                continue
            k, v = line.rstrip("\n").split("\t", 1)
            try:
                same = float(v) == float(want.get(k, "nan"))
            except ValueError:
                same = v == want.get(k)
            if not same:
                die(f"the simulator used {name_of.get(k, k)} = {v}, and the config says {want.get(k)}")
    print(f"  -> {os.path.join(d, 'manifest.yml')}")
    print(f"batch ready: {d}")
