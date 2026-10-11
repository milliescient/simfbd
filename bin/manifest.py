#!/usr/bin/env python3
"""Write and read a batch's or a run's manifest.yml.

    manifest.py write-batch DIR field=value ...   # settings come from the environment
    manifest.py write-run   DIR field=value ...
    manifest.py get  DIR KEY                      # one value, by config key or record field
    manifest.py flat DIR                          # every value, as key<TAB>value lines
    manifest.py status DIR                        # does the directory still match its manifest

A manifest says what a batch or a run is, in words. config.sh keeps its short shell names,
because bin/run.sh sources it, and the SETTINGS tables below are the one place that says what each of
them means: every line of a manifest carries the name a reader wants, and beside it the
config key it came from.

Beside the settings a manifest records what a notebook's lock records, in the same words, so
a batch or a run can be checked with no notebook: `uses` names the code that ran as
name@version, and `inputs` and `outputs` give each path with its hash. The hash is digest()
below, which is the one sloop computes, so a path hashes the same in both.

Scripts read a manifest through `get` and `flat`, which answer by config key. A batch
simulated before manifest.yml existed has a manifest.tsv of key<TAB>value lines instead, and
both commands read that too, so a script never has to know which kind it was handed.
"""
import hashlib
import os
import re
import sys
from pathlib import Path

import yaml

# config key, where it sits in the manifest, what it means
BATCH_SETTINGS = [
    ("NINTERVALS", "timeline.rate_intervals",
     "piecewise-constant rate intervals"),
    ("INTERVAL_WIDTH", "timeline.interval_width_myr",
     "width of each interval; the oldest runs on to the origin"),
    ("PRESENT", "timeline.record_stops_at_myr",
     "the record ends here, not at zero; younger occurrences are discarded"),
    ("AGE_MIN", "origin_age.min_myr", "origin ~ Uniform(min, max)"),
    ("AGE_MAX", "origin_age.max_myr", ""),
    ("ORIGIN_SAMPLED", "origin_age.origin_lineage_must_be_sampled", ""),

    ("LMEAN", "rate_priors.speciation.log_mean",
     "lognormal, on the log scale; one draw per interval"),
    ("LSD", "rate_priors.speciation.log_sd", ""),
    ("MMEAN", "rate_priors.extinction.log_mean", ""),
    ("MSD", "rate_priors.extinction.log_sd", ""),
    ("TMEAN", "rate_priors.turnover.log_mean",
     "stationary prior only: extinction is speciation times turnover"),
    ("TSD", "rate_priors.turnover.log_sd", ""),
    ("PMEAN", "rate_priors.fossil_sampling.log_mean", ""),
    ("PSD", "rate_priors.fossil_sampling.log_sd", ""),
    ("LAMBDA_A", "rate_priors.anagenetic_speciation",
     "above 0, an anagenetic rate is drawn from the speciation prior; 0 is budding only"),

    ("GMRF_SD", "rates_through_time.random_walk.step_sd",
     "sd of the step between adjacent intervals' log rates"),
    ("GMRF_SD_MEAN", "rates_through_time.random_walk.step_sd_prior.log_mean",
     "speciation and extinction: the step sd is drawn per replicate, lognormal"),
    ("GMRF_SD_SD", "rates_through_time.random_walk.step_sd_prior.log_sd", ""),
    ("GMRF_PSI_MEAN", "rates_through_time.random_walk.fossil_sampling_step_sd_prior.log_mean",
     "fossil sampling gets its own step sd"),
    ("GMRF_PSI_SD", "rates_through_time.random_walk.fossil_sampling_step_sd_prior.log_sd", ""),
    ("OU_THETA", "rates_through_time.stationary.mean_reversion_per_myr",
     "fixed; set this or correlation_prior, not both"),
    ("OU_DREF", "rates_through_time.stationary.reference_lag_myr",
     "the lag the correlation is stated at"),
    ("OU_PHI_A", "rates_through_time.stationary.correlation_prior.beta_a",
     "correlation at the reference lag, drawn per replicate from Beta(a, b)"),
    ("OU_PHI_B", "rates_through_time.stationary.correlation_prior.beta_b", ""),
    ("LANCH", "rates_through_time.stationary.level_sd.speciation",
     "sd of each replicate's mean log rate about the prior mean; 0 fixes it there"),
    ("TANCH", "rates_through_time.stationary.level_sd.turnover", ""),
    ("PANCH", "rates_through_time.stationary.level_sd.fossil_sampling", ""),
    ("OU_SD_MEAN", "rates_through_time.stationary.spread_prior.log_mean",
     "speciation and turnover: the sd through time is drawn per replicate, lognormal"),
    ("OU_SD_SD", "rates_through_time.stationary.spread_prior.log_sd", ""),
    ("OU_PSI_MEAN", "rates_through_time.stationary.fossil_sampling_spread_prior.log_mean",
     "fossil sampling gets its own"),
    ("OU_PSI_SD", "rates_through_time.stationary.fossil_sampling_spread_prior.log_sd", ""),

    ("RHO", "record.extant_sampling_probability",
     "chance a lineage alive at the present is seen there"),
    ("REPORT_PROB", "record.interior_reporting",
     "below 1, interior occurrences are thinned at a rate drawn per replicate from Uniform(0, 1)"),
    ("BIN_WIDTH", "record.occurrence_bins.width_myr",
     "occurrence ages are reported to these bins"),
    ("BIN_MAX", "record.occurrence_bins.oldest_edge_myr", ""),

    ("MIN_TAXA", "accepted_replicates.min_taxa",
     "a replicate outside these bounds is redrawn"),
    ("MAX_TAXA", "accepted_replicates.max_taxa", ""),
    ("TAXA_COUNT", "accepted_replicates.taxa_counted",
     "fossil: fossil-sampled taxa only; record: survivors seen at the present too"),
    ("MIN_FOSSIL_TAXA", "accepted_replicates.min_fossil_taxa", ""),

    ("MAX_LINEAGES", "simulation_limits.max_lineages",
     "a draw expected to exceed this is refused unsimulated, which truncates the prior"),
    ("SIM_MAX_DRAWS", "simulation_limits.max_waiting_time_draws",
     "a draw whose simulation needs more waiting times than this is redrawn"),
]

RUN_SETTINGS = [
    ("MODEL", "inference.model", "fbdr: dnFBDRP; bds: dnBDS; fbdsp: dnFBDSP"),
    ("REPORTING", "inference.record_read", "the batch's specimens_<this> directory"),
    ("INFER", "inference.record_treated_as", "complete or incomplete"),
    ("COND", "inference.conditioned_on", ""),
    ("RHO", "inference.extant_sampling_probability", ""),
    ("ORIGIN_PRIOR", "inference.origin_prior", ""),
    ("SURVIVORS", "inference.survivors", ""),
    ("DROP_PRESENT_ONLY", "inference.drop_present_only_taxa", ""),
    ("GENS", "mcmc.generations", ""),
    ("PRINTGEN", "mcmc.sample_every", ""),
]

# what a manifest records beyond settings: the name a script asks for, and where it sits
BATCH_RECORD = {
    "kind": "kind", "batch": "batch", "config_hash": "config_hash",
    "true_vals_order": "true_vals_order", "NREPS": "replicates",
    "screened": "prior_screen.refused", "screened_of": "prior_screen.of",
    "screened_frac": "prior_screen.fraction",
}
RUN_RECORD = {
    "kind": "kind", "run": "run", "batch": "batch",
    "batch_config_hash": "batch_config_hash", "config_hash": "config_hash",
    "true_vals_order": "true_vals_order", "NREPS": "replicates", "specimens": "specimens",
    "infer_script_hash": "toolchain.infer_script_hash",
    "rb_path": "toolchain.rb_path", "rb_md5": "toolchain.rb_md5",
    "rb_describe": "toolchain.rb_describe",
    # the names a manifest.tsv used for these two
    "survivors": "inference.survivors", "drop_present_only": "inference.drop_present_only_taxa",
}

NUMBER = re.compile(r"-?(0|[1-9]\d*)(\.\d+)?([eE][-+]?\d+)?")
WORD = re.compile(r"[A-Za-z][A-Za-z0-9_./-]*")
RESERVED = {"true", "false", "yes", "no", "on", "off", "null", "none", "y", "n"}


def setting(v):
    """A config value as YAML: numbers and true/false bare, anything else as a string."""
    return v if NUMBER.fullmatch(v) or v in ("true", "false") else text(v)


def text(v):
    """A hash, a commit or a path. Quoted unless it is plainly a word, so that a hash made
    only of digits is not read back as a number."""
    v = str(v)
    if WORD.fullmatch(v) and v.lower() not in RESERVED:
        return v
    return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'


def digest(path):
    """A stable hash of a file, or of every file under a directory: SHA-1 over one line per
    file, its relative name and the SHA-1 of its bytes, files in sorted order, first 12 hex
    digits. The same as sloop's digest(), on purpose, and kept in step with it by hand."""
    path = Path(path)
    base = path if path.is_dir() else path.parent
    files = sorted(p for p in path.rglob("*") if p.is_file()) if path.is_dir() else [path]
    h = hashlib.sha1()
    for f in files:
        one = hashlib.sha1()
        with f.open("rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                one.update(block)
        h.update(f"{f.relative_to(base)}\0{one.hexdigest()}\n".encode())
    return h.hexdigest()[:12]


# what a simulation leaves in a batch for something else to read
def batch_outputs(d):
    names = ["true_vals.tsv", "nums.tsv", "times", "true_ranges", "ou_lags.tsv"]
    names += sorted(n for n in os.listdir(d) if n.startswith("specimens_"))
    return [n for n in names if os.path.exists(os.path.join(d, n))]


def used(uses, name):
    """The version of one thing in a `uses` list of name@version."""
    for u in uses or []:
        if str(u).split("@", 1)[0] == name:
            return str(u).split("@", 1)[1]
    return None


def dig(doc, path):
    for part in path.split("."):
        if not isinstance(doc, dict) or part not in doc:
            return None
        doc = doc[part]
    return doc


def shown(v):
    return str(v).lower() if isinstance(v, bool) else str(v)


def settings_lines(settings, value_of, indent=""):
    """The settings that are set, nested by section, each beside the config key it came from."""
    rows, open_path = [], []
    for key, path, meaning in settings:
        v = value_of(key)
        if v is None or v == "":
            continue
        *sections, name = path.split(".")
        common = 0
        while (common < min(len(open_path), len(sections))
               and open_path[common] == sections[common]):
            common += 1
        if common == 0 and rows:
            rows.append(("", ""))
        for depth in range(common, len(sections)):
            rows.append((f"{indent}{'  ' * depth}{sections[depth]}:", ""))
        open_path = sections
        note = f"{key}  {meaning}".rstrip()
        rows.append((f"{indent}{'  ' * len(sections)}{name}: {setting(v)}", note))
    width = max((len(left) for left, note in rows if note), default=0) + 2
    return [f"{left.ljust(width)}# {note}" if note else left for left, note in rows]


def find(d):
    """The manifest in a directory, as (path, parsed YAML or None for a manifest.tsv)."""
    yml, tsv = os.path.join(d, "manifest.yml"), os.path.join(d, "manifest.tsv")
    if os.path.exists(yml):
        with open(yml) as fh:
            return yml, yaml.safe_load(fh) or {}
    if os.path.exists(tsv):
        return tsv, None
    return None, None


def flat(d):
    """Every value by the key a script asks for, from either kind of manifest."""
    path, doc = find(d)
    if path is None:
        return {}
    if doc is None:
        out = {}
        with open(path) as fh:
            for line in fh:
                if "\t" in line:
                    k, v = line.rstrip("\n").split("\t", 1)
                    out[k] = v              # a run repeats RHO; the later line is the fitted one
        if out.get("kind") == "run":
            out["batch_dir"] = os.path.dirname(os.path.dirname(os.path.normpath(d))) or "."
        return out
    out = {}
    run = doc.get("kind") == "run"
    for name, where in (RUN_RECORD if run else BATCH_RECORD).items():
        v = dig(doc, where)
        if v is not None:
            out[name] = shown(v)
    if not run and doc.get("simulated"):
        # what the newest replicates were simulated with, which is what an extension must match
        last = doc["simulated"][-1]
        out["sim_script_hash"] = shown(last["sim_script_hash"])
        out["simfbd_commit"] = shown(used(last.get("uses"), "simfbd") or last.get("simfbd_commit"))
    if run:
        # a run sits at <batch>/runs/<name>, wherever the batch is, so this is never recorded
        out["batch_dir"] = os.path.dirname(os.path.dirname(os.path.normpath(d))) or "."
        tc = doc.get("toolchain") or {}
        out["simfbd_commit"] = shown(used(doc.get("uses"), "simfbd") or tc.get("simfbd_commit"))
        out["rb_commit"] = shown(tc.get("rb_describe") or tc.get("rb_commit"))
    for key, where, _ in BATCH_SETTINGS:
        v = dig(doc, ("simulation." if run else "") + where)
        # a manifest.tsv wrote every key, and readers test for this word
        out[key] = "none" if v is None else shown(v)
    if run:
        for key, where, _ in RUN_SETTINGS:
            v = dig(doc, where)
            if v is not None:
                out[key] = shown(v)
    return out


def fields(args):
    bad = [a for a in args if "=" not in a]
    if bad:
        raise SystemExit(f"manifest.py: expected field=value, got {bad[0]}")
    return dict(a.split("=", 1) for a in args)


def save(d, lines):
    body = "\n".join(lines) + "\n"
    yaml.safe_load(body)                    # a manifest nothing can read back is worse than none
    tmp = os.path.join(d, "manifest.yml.tmp")
    with open(tmp, "w") as fh:
        fh.write(body)
    os.replace(tmp, os.path.join(d, "manifest.yml"))


def write_batch(d, f):
    _, old = find(d)
    first, last = int(f["first"]), int(f["last"])
    simulated = list((old or {}).get("simulated") or [])
    have = int(str(simulated[-1]["replicates"]).split("-")[-1]) if simulated else 0
    if first != have + 1:
        raise SystemExit(f"manifest.py: {d} holds replicates to {have}, "
                         f"so the next cannot start at {first}")
    simulated.append({"replicates": f"{first}-{last}",
                      "sim_script_hash": f["sim_script_hash"],
                      "uses": [f"simfbd@{f['simfbd_commit']}"],
                      "outputs": {n: digest(os.path.join(d, n)) for n in batch_outputs(d)}})
    lines = [
        "# Written by bin/run.sh. Each setting is followed by its config.sh key.",
        "kind: batch",
        f"batch: {text(os.path.basename(os.path.normpath(d)))}",
        f"config_hash: {text(f['config_hash'])}",
        "# interval 1 of each rate in true_vals.tsv is the youngest, as RevBayes numbers them",
        "true_vals_order: youngest_first",
        "",
        f"replicates: {last}",
        "# One entry per extension. outputs are the batch's files as they stood after it, so",
        "# only the last entry can still be checked: simfbd status <batch>",
        "simulated:",
    ]
    for s in simulated:
        uses = s.get("uses") or [f"simfbd@{s['simfbd_commit']}"]
        lines += [f"  - replicates: {text(s['replicates'])}",
                  f"    uses: [{', '.join(text(u) for u in uses)}]",
                  f"    sim_script_hash: {text(s['sim_script_hash'])}"]
        if s.get("outputs"):
            lines.append("    outputs:")
            lines += [f"      {n}: {text(h)}" for n, h in s["outputs"].items()]
    screen = os.path.join(d, "screened.tsv")
    if os.path.exists(screen):
        with open(screen) as fh:
            s = dict(l.rstrip("\n").split("\t", 1) for l in fh if "\t" in l)
        lines += ["# draws refused before simulating by simulation_limits.max_lineages, measured",
                  "prior_screen:",
                  f"  refused: {s['screened']}",
                  f"  of: {s['screened_of']}",
                  f"  fraction: {s['screened_frac']}"]
    lines += ["", "# ---- the generating process, from config.sh ----"]
    lines += settings_lines(BATCH_SETTINGS, os.environ.get)
    save(d, lines)


def write_run(d, f):
    batch_path, batch = find(f["batch_dir"])
    if batch is None:
        raise SystemExit(f"manifest.py: {f['batch_dir']} has no manifest.yml")
    sim = flat(f["batch_dir"])
    lines = [
        "# Written by bin/run.sh. Each setting is followed by its config.sh key.",
        "kind: run",
        f"run: {text(os.path.basename(os.path.normpath(d)))}",
        f"batch: {text(batch['batch'])}",
        f"batch_config_hash: {text(f['batch_config_hash'])}",
        f"config_hash: {text(f['config_hash'])}",
        f"true_vals_order: {text(batch['true_vals_order'])}",
        f"replicates: {f['replicates']}            # in the batch when this last ran",
        f"specimens: {text(f['specimens'])}",
        "",
        "# the code that ran, as name@version; rb's is the commit its binary was built from",
        f"uses: [{text('simfbd@' + f['simfbd_commit'])}, "
        f"{text('rb@' + f['rb_commit'].rsplit('-g', 1)[-1])}]",
        "# read from the batch, and written here; simfbd status <run> checks both",
        "inputs:",
        f"  {f['specimens']}: {text(digest(os.path.join(f['batch_dir'], f['specimens'])))}",
        "outputs:",
        f"  output: {text(digest(os.path.join(d, 'output')))}",
        "",
        "toolchain:",
        f"  infer_script_hash: {text(f['infer_script_hash'])}",
        f"  rb_path: {text(f['rb_path'])}",
        f"  rb_md5: {text(f['rb_md5'])}",
        f"  rb_describe: {text(f['rb_commit'])}",
    ]
    lines += ["", "# ---- the fit, from this run's config.sh ----"]
    lines += settings_lines(RUN_SETTINGS, os.environ.get)
    lines += ["", f"# ---- the generating process, copied from {batch_path} ----",
              "simulation:"]
    lines += settings_lines(BATCH_SETTINGS,
                            lambda k: None if sim.get(k) == "none" else sim.get(k),
                            indent="  ")
    save(d, lines)


def sha1_12(path):
    with open(path, "rb") as fh:
        return hashlib.sha1(fh.read()).hexdigest()[:12]


def status(d):
    """Whether a batch or a run still matches its manifest. Returns the number of problems."""
    path, doc = find(d)
    if path is None:
        raise SystemExit(f"manifest.py: no manifest in {d}")
    if doc is None:
        print(f"{d}: manifest.tsv predates recorded hashes, so there is nothing to check")
        return 0
    rows, bad = [], 0

    def row(state, what, note=""):
        nonlocal bad
        bad += state in ("changed", "missing", "stale")
        rows.append((state, what, note))

    def compare(base, recorded, label=""):
        for name, want in (recorded or {}).items():
            at = os.path.join(base, name)
            if not os.path.exists(at):
                row("missing", label + name)
            else:
                got = digest(at)
                row("current" if got == want else "changed", label + name,
                    "" if got == want else f"{want} -> {got}")

    cfg = os.path.join(d, "config.sh")
    if os.path.exists(cfg):
        got = sha1_12(cfg)
        row("current" if got == doc.get("config_hash") else "changed", "config.sh",
            "" if got == doc.get("config_hash") else f"{doc.get('config_hash')} -> {got}")
    else:
        row("missing", "config.sh")

    if doc.get("kind") == "run":
        batch_dir = os.path.dirname(os.path.dirname(os.path.normpath(d))) or "."
        _, batch = find(batch_dir)
        for u in doc.get("uses") or []:
            row("uses", u)
        compare(batch_dir, doc.get("inputs"), "batch: ")
        compare(d, doc.get("outputs"))
        if not doc.get("inputs") and not doc.get("outputs"):
            row("unrecorded", "inputs and outputs", "written before hashes were recorded")
        if batch:
            if batch.get("config_hash") != doc.get("batch_config_hash"):
                row("changed", "batch: config.sh",
                    f"{doc.get('batch_config_hash')} -> {batch.get('config_hash')}")
            # the edge: what this run read must be what the batch says it wrote
            wrote = (batch.get("simulated") or [{}])[-1].get("outputs") or {}
            for name, read in (doc.get("inputs") or {}).items():
                if name in wrote and wrote[name] != read:
                    row("stale", "batch: " + name,
                        f"the batch has changed since this run was fitted "
                        f"({doc.get('replicates')} -> {batch.get('replicates')} replicates)")
    else:
        last = (doc.get("simulated") or [{}])[-1]
        for u in last.get("uses") or []:
            row("uses", u)
        if last.get("outputs"):
            compare(d, last["outputs"])
        else:
            row("unrecorded", "outputs", "simulated before hashes were recorded")

    width = max(len(s) for s, _, _ in rows)
    print(f"{d}  ({doc.get('kind')}, {doc.get('replicates')} replicates)")
    for state, what, note in rows:
        print(f"  {state.ljust(width)}  {what}" + (f"  ({note})" if note else ""))
    return bad


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    cmd, d, rest = argv[0], argv[1].rstrip("/"), argv[2:]
    if cmd == "write-batch":
        write_batch(d, fields(rest))
    elif cmd == "write-run":
        write_run(d, fields(rest))
    elif cmd == "get":
        v = flat(d).get(rest[0])
        if v is not None:
            print(v)
    elif cmd == "flat":
        for k, v in flat(d).items():
            print(f"{k}\t{v}")
    elif cmd == "status":
        sys.exit(1 if sum(status(x.rstrip("/")) for x in [d] + rest) else 0)
    elif cmd == "digest":
        print(digest(d))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
