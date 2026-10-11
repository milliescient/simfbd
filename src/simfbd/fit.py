"""Fit a run: every replicate of its batch not yet fitted.

A run chooses only what to read and how to fit it. It inherits the batch's hyperpriors and
has no way to set them, so the generating and the fitting priors cannot drift. Its manifest
records the batch config it ran against, and a config edited after its output exists is
refused rather than mixed.

A replicate whose log is complete is skipped, so a rerun fills gaps, and a run under a batch
that has been extended fits only the new replicates.
"""
import fcntl
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

from . import config, manifest, tools
from .batch import die, replicates, warn_if_dirty

TEMPLATES = {"fbdr": "infer-fbdrp.Rev",      # dnFBDRP: accounts for unsampled lineages
             "bds": "infer-bds.Rev",         # dnBDS: complete lineage sampling, PyRate's assumption
             "fbdsp": "infer-fbdsp.Rev"}     # dnFBDSP: the tree process, the only one with lambda_a

# what a template is handed for a setting the batch left out
UNSET = {"GMRF_SD": "0", "GMRF_SD_MEAN": "-1", "GMRF_SD_SD": "0", "GMRF_PSI_MEAN": "-1",
         "GMRF_PSI_SD": "0", "OU_THETA": "-1", "OU_DREF": "1.5", "OU_PHI_A": "-1", "OU_PHI_B": "-1",
         "TMEAN": "0", "TSD": "0.3", "LANCH": "0", "TANCH": "0", "PANCH": "0", "OU_SD_MEAN": "-1",
         "OU_SD_SD": "0", "OU_PSI_MEAN": "-1", "OU_PSI_SD": "0", "LAMBDA_A": "0", "PRESENT": "0",
         "REPORT_PROB": "1"}
HEAD = ["LMEAN", "LSD", "MMEAN", "MSD", "PMEAN", "PSD", "AGE_MIN", "AGE_MAX"]
TAIL = ["OU_THETA", "OU_DREF", "OU_PHI_A", "OU_PHI_B", "TMEAN", "TSD", "LANCH", "TANCH", "PANCH",
        "OU_SD_MEAN", "OU_SD_SD", "OU_PSI_MEAN", "OU_PSI_SD", "LAMBDA_A", "PRESENT", "REPORT_PROB"]


def lines_in(path):
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def max_first_appearance(taxa):
    """The oldest bin floor across taxa: PyRate's max FA, and the tightest lower bound the
    data puts on the origin. Only used when the origin prior is exponential."""
    oldest = {}
    with open(taxa) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                oldest[f[0]] = max(oldest.get(f[0], 0.0), float(f[1]))
            except (IndexError, ValueError):
                continue
    x = max(oldest.values(), default=0.0)
    return f"{x if x > 0 else 1:.6f}"


def fit(run, force=False, slots=None, slot_dir=None):
    run = os.path.normpath(run)
    batch = os.path.dirname(os.path.dirname(run))
    if os.path.basename(os.path.dirname(run)) != "runs":
        die("fit takes a run, <batch>/runs/<name>")
    bcfg, rcfg = os.path.join(batch, "config.yml"), os.path.join(run, "config.yml")
    for f in (bcfg, rcfg):
        if not os.path.exists(f):
            die(f"no {f}")
    b = config.load(bcfg, "batch")
    r = config.load(rcfg, "run")
    b.pop("_uses")
    tools.require(r.pop("_uses"), rcfg)
    bhash, rhash = manifest.digest(bcfg), manifest.digest(rcfg)
    _, bdoc = manifest.find(batch)
    if bdoc is None or not replicates(batch):
        die(f"{batch} has no replicates yet: simfbd extend {batch} <how many>")
    simulated_from = manifest.flat(batch).get("config_hash")
    if simulated_from != bhash:
        die(f"{bcfg} changed since its replicates were simulated ({simulated_from} -> {bhash}).\n"
            f"       Fitting them under it would describe data that was not drawn that way. Make a new batch.")
    n = replicates(batch)

    # rho is the generating parameter and the fitted one. A batch that sets it wins, and a run
    # that disagrees is refused rather than silently fitting a rho the data was not drawn under.
    if "RHO" in b and float(b["RHO"]) != float(r["RHO"]):
        die(f"{rcfg} fits an extant sampling probability of {r['RHO']}, and {batch} was "
            f"generated under {b['RHO']}.")
    template = tools.script(TEMPLATES[r["MODEL"]])
    complete = "true" if r["INFER"] == "complete" else "false"
    specimens = f"specimens_{r['REPORTING']}"
    if not os.path.isdir(os.path.join(batch, specimens)):
        die(f"the batch has no {specimens}")
    sky = "true" if int(b["NINTERVALS"]) > 1 else "false"

    path, _ = manifest.find(run)
    if path is not None and not force:
        was = manifest.flat(run)
        if was.get("batch_config_hash") != bhash:
            die(f"{run} ran against batch config {was.get('batch_config_hash')}, now {bhash}. "
                f"Its output would mix. Delete output/ or pass --force.")
        if was.get("config_hash") != rhash:
            die(f"{rcfg} changed since its output was written ({was.get('config_hash')} -> {rhash}). "
                f"Delete output/ or pass --force.")

    out, aux = os.path.join(run, "output"), os.path.join(run, "aux")
    for d in (out, aux, os.path.join(run, "results")):
        os.makedirs(d, exist_ok=True)
    print(f"{os.path.basename(run)}: batch={os.path.basename(batch)} reads={r['REPORTING']} "
          f"treated as {r['INFER']} conditioned on {r['COND']} rho={r['RHO']} reps={n}", flush=True)
    warn_if_dirty()
    tools.announce()

    # a finished replicate's log has this many lines; a rerun fills gaps rather than redo work
    full = int(r["GENS"]) // int(r["PRINTGEN"]) + 2
    failures = os.path.join(run, "failures.log")
    open(failures, "w").close()
    absolute = os.path.abspath
    rb, env = tools.rb(), tools.environment()
    slots = int(slots or os.environ.get("SLOTS") or 0)
    slot_dir = slot_dir or os.environ.get("SLOT_DIR") or os.path.join(os.path.dirname(batch), ".slots")
    if slots:
        os.makedirs(slot_dir, exist_ok=True)

    def one(rep):
        log = os.path.join(out, f"rep_{rep}.log")
        if os.path.exists(log) and lines_in(log) == full:
            return
        # one self-contained file per replicate: the settings as Rev variables, then the template
        v = {**UNSET, **b, **r}
        # seeded, so the same replicate fitted twice gives the same chain. Each replicate has
        # its own seed, or their Monte Carlo errors would move together.
        head = [f"seed({int(r['SEED']) + rep})",
                f'rep <- "{rep}"', f'ORIGIN_PRIOR <- "{r["ORIGIN_PRIOR"]}"',
                f"MAXFA <- {max_first_appearance(os.path.join(batch, specimens, f'taxa_{rep}.tsv'))}",
                f"SKYLINE <- {sky}", f"COMPLETE <- {complete}", f'COND <- "{r["COND"]}"',
                f"RHO <- {r['RHO']}", f"GENS <- {r['GENS']}", f"PRINTGEN <- {r['PRINTGEN']}",
                f'BATCHDIR <- "{absolute(batch)}"', f'SPECIMENS <- "{specimens}"',
                f'OUTDIR <- "{absolute(out)}"']
        head += [f"{k} <- {v[k]}" for k in HEAD]
        head += [f"SURVIVORS <- {r['SURVIVORS']}", f"DROP_PRESENT_ONLY <- {r['DROP_PRESENT_ONLY']}"]
        head += [f"{k} <- {v[k]}" for k in ("GMRF_SD", "GMRF_SD_MEAN", "GMRF_SD_SD", "GMRF_PSI_MEAN", "GMRF_PSI_SD")]
        head += [f'BINDIR <- "{tools.SCRIPTS}"']
        head += [f"{k} <- {v[k]}" for k in TAIL]
        head += [f'source("{template}")']
        rev = os.path.join(aux, f"run_{rep}.Rev")
        with open(rev, "w") as fh:
            fh.write("\n".join(head) + "\n")
        # slots cap chains across every fit sharing the directory, so concurrent runs fill the
        # cores between them: each chain holds one slot lock for as long as rb runs
        held = None
        while slots and held is None:
            for s in range(1, slots + 1):
                fh = open(os.path.join(slot_dir, f"slot_{s}"), "w")
                try:
                    fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    held = fh
                    break
                except OSError:
                    fh.close()
            if held is None:
                time.sleep(10)
        try:
            with open(os.path.join(aux, f"rb_{rep}.out"), "w") as so:
                subprocess.run([rb, absolute(rev)], stdin=subprocess.DEVNULL, stdout=so,
                               stderr=subprocess.STDOUT, env=env)
        finally:
            if held is not None:
                held.close()
        if not os.path.exists(log) or lines_in(log) != full:
            why = "no log written"
            for line in open(os.path.join(aux, f"rb_{rep}.out"), errors="replace"):
                if "error" in line.lower():
                    why = line.strip()
                    break
            with open(failures, "a") as fh:
                fh.write(f"{rep} :: {why}\n")

    with ThreadPoolExecutor(max_workers=max(1, int(r["NCORES"]))) as pool:
        list(pool.map(one, range(1, n + 1)))

    rb_path, rb_md5, rb_describe = tools.rb_identity()
    manifest.write_run(run, {"batch_dir": batch, "replicates": n, "specimens": specimens,
                             "inputs": {"config.yml": rhash,
                                        "batch": {"config.yml": bhash,
                                                  specimens: manifest.digest(os.path.join(batch, specimens))},
                                        "simfbd": {template.name: manifest.digest(template)}},
                             "simfbd_commit": tools.version(), "rb_path": rb_path,
                             "rb_md5": rb_md5, "rb_commit": rb_describe})
    done = sum(1 for f in os.listdir(out) if f.startswith("rep_") and f.endswith(".log"))
    print(f"logs: {done}/{n}, {lines_in(failures)} failed")
    print(f"manifest: {os.path.join(run, 'manifest.yml')}")
    print(f"now run: simfbd summarize {run}")
    return lines_in(failures) == 0
