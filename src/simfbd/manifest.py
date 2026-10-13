#!/usr/bin/env python3
"""Write, read and check a batch's or a run's manifest.yml.

A config is what was asked for and a manifest is what was done, and a setting is written
once, in the config. The manifest is the record: the hash of the config it was made from, how
many replicates there are, which code ran as name@version under `uses`, and the hash of each
path it read and wrote under `inputs` and `outputs`. Those are the words a notebook's lock
uses, and the hash is digest() below, which is the one sloop computes, so a batch or a run
can be checked with no notebook and a path hashes alike in both.

Scripts read both through flat(), which answers by the name the simulator knows a setting
by: the settings from the config and the record from the manifest. A batch made before
configs were YAML has its settings in its manifest, as a manifest.tsv of key<TAB>value lines
or an earlier manifest.yml, and flat() reads those too, so a script never has to know which
kind it was handed.

    python -m simfbd.manifest get  DIR KEY      # one value
    python -m simfbd.manifest flat DIR          # every value, as key<TAB>value lines
"""
import hashlib
import os
import re
import sys
from pathlib import Path

import yaml

from . import config
from .settings import BATCH_RECORD, BATCH_SETTINGS, RUN_RECORD, RUN_SETTINGS

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
    names = ["true_vals.tsv", "nums.tsv", "times", "true_ranges"]
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
                    out[k] = v              # a run repeats RHO; the later line is the analysis's
        # a manifest.tsv says screened where everything since says filtered
        for k in [k for k in out if k.startswith("screened")]:
            out.setdefault(k.replace("screened", "filtered"), out[k])
        if out.get("kind") == "run":
            out["batch_dir"] = os.path.dirname(os.path.dirname(os.path.normpath(d))) or "."
        return out
    out = {}
    run = doc.get("kind") == "run"
    for name, where in (RUN_RECORD if run else BATCH_RECORD).items():
        v = dig(doc, where)
        if v is None and where.startswith("prior_filter."):
            v = dig(doc, where.replace("prior_filter.", "prior_screen."))
        if v is not None:
            out[name] = shown(v)
    if not run and doc.get("simulated"):
        # what the newest replicates were simulated with, which is what an extension must match
        last = doc["simulated"][-1]
        given = last.get("inputs") or {}
        if given:
            out["config_hash"] = shown(given.get("config.yml"))
        out["sim_script_hash"] = shown((given.get("simfbd") or {}).get("sim.R") or last.get("sim_script_hash"))
        out["simfbd_commit"] = shown(used(last.get("uses"), "simfbd") or last.get("simfbd_commit"))
    if run:
        # a run sits at <batch>/runs/<name>, wherever the batch is, so this is never recorded
        out["batch_dir"] = os.path.dirname(os.path.dirname(os.path.normpath(d))) or "."
        tc = doc.get("toolchain") or {}
        out["simfbd_commit"] = shown(used(doc.get("uses"), "simfbd") or tc.get("simfbd_commit"))
        out["rb_commit"] = shown(tc.get("rb_describe") or tc.get("rb_commit"))
        given = doc.get("inputs") or {}
        if isinstance(given.get("batch"), dict):
            out["config_hash"] = shown(given.get("config.yml"))
            out["batch_config_hash"] = shown(given["batch"].get("config.yml"))
            out["infer_script_hash"] = shown(next(iter((given.get("simfbd") or {}).values()), None))
    batch_dir = out["batch_dir"] if run else d
    if os.path.exists(os.path.join(batch_dir, "config.yml")):
        # the settings are in the configs. A batch's are given as it set them, and a run's
        # with the defaults filled in, which is what a manifest used to hold of each.
        asked = config.load(os.path.join(batch_dir, "config.yml"), "batch", with_defaults=False)
        for key, _, _ in BATCH_SETTINGS:
            out[key] = asked.get(key, "none")       # readers test for this word
        if run:
            out.update(config.load(os.path.join(d, "config.yml"), "run"))
            out.pop("_uses", None)
        return out
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


def save(d, lines):
    body = "\n".join(lines) + "\n"
    yaml.safe_load(body)                    # a manifest nothing can read back is worse than none
    tmp = os.path.join(d, "manifest.yml.tmp")
    with open(tmp, "w") as fh:
        fh.write(body)
    os.replace(tmp, os.path.join(d, "manifest.yml"))


def hash_lines(given, indent):
    """Paths and their hashes as YAML, a group of them under its name."""
    out = []
    for name, v in given.items():
        if isinstance(v, dict):
            out.append(f"{indent}{name}:")
            out += [f"{indent}  {n}: {text(h)}" for n, h in v.items()]
        else:
            out.append(f"{indent}{name}: {text(v)}")
    return out


def write_batch(d, f):
    _, old = find(d)
    first, last = int(f["first"]), int(f["last"])
    simulated = list((old or {}).get("simulated") or [])
    have = int(str(simulated[-1]["replicates"]).split("-")[-1]) if simulated else 0
    if first != have + 1:
        raise SystemExit(f"manifest.py: {d} holds replicates to {have}, "
                         f"so the next cannot start at {first}")
    entry = {"replicates": f"{first}-{last}",
             "uses": list(f.get("uses") or [f"simfbd@{f['simfbd_commit']}"])}
    if "inputs" in f:
        entry["inputs"] = f["inputs"]
    else:
        entry["sim_script_hash"] = f["sim_script_hash"]      # a batch whose config is a config.sh
    entry["outputs"] = {n: digest(os.path.join(d, n)) for n in batch_outputs(d)}
    simulated.append(entry)
    lines = [
        "# Written by simfbd: the record of what was simulated from config.yml. Do not edit it.",
        "kind: batch",
        f"batch: {text(os.path.basename(os.path.normpath(d)))}",
    ]
    if "inputs" not in f:
        lines.append(f"config_hash: {text(f['config_hash'])}")
    lines += [
        "# interval 1 of each rate in true_vals.tsv is the youngest, as RevBayes numbers them",
        "true_vals_order: youngest_first",
        "",
        f"replicates: {last}",
        "# One entry per simfbd extend: what ran, what it read and what the batch held after it.",
        "# Only the last entry can still be checked, since the next one changes the files:",
        "# simfbd status <batch>",
        "simulated:",
    ]
    for s in simulated:
        uses = s.get("uses") or [f"simfbd@{s['simfbd_commit']}"]
        lines += [f"  - replicates: {text(s['replicates'])}",
                  f"    uses: [{', '.join(text(u) for u in uses)}]"]
        if s.get("inputs"):
            lines.append("    inputs:")
            lines += hash_lines(s["inputs"], "      ")
        else:
            lines.append(f"    sim_script_hash: {text(s['sim_script_hash'])}")
        if s.get("outputs"):
            lines.append("    outputs:")
            lines += [f"      {n}: {text(h)}" for n, h in s["outputs"].items()]
    counts = next((os.path.join(d, n) for n in ("filtered.tsv", "screened.tsv")
                   if os.path.exists(os.path.join(d, n))), None)
    if counts:
        with open(counts) as fh:
            s = {k.replace("screened", "filtered"): v
                 for k, v in (l.rstrip("\n").split("\t", 1) for l in fh if "\t" in l)}
        lines += ["# draws refused before simulating by simulation_limits.max_lineages, measured",
                  "prior_filter:",
                  f"  refused: {s['filtered']}",
                  f"  of: {s['filtered_of']}",
                  f"  fraction: {s['filtered_frac']}"]
    save(d, lines)


def write_run(d, f):
    batch_path, batch = find(f["batch_dir"])
    if batch is None:
        raise SystemExit(f"manifest.py: {f['batch_dir']} has no manifest.yml")
    sim = flat(f["batch_dir"])
    lines = [
        "# Written by simfbd: the record of what was analyzed from config.yml. Do not edit it.",
        "kind: run",
        f"run: {text(os.path.basename(os.path.normpath(d)))}",
        f"batch: {text(batch['batch'])}",
    ]
    if "inputs" not in f:
        lines += [f"batch_config_hash: {text(f['batch_config_hash'])}",
                  f"config_hash: {text(f['config_hash'])}"]
    lines += [
        f"true_vals_order: {text(batch['true_vals_order'])}",
        f"replicates: {f['replicates']}            # in the batch when this last ran",
        f"specimens: {text(f['specimens'])}",
        "",
        "# the code that ran, as name@version; rb's is the commit its binary was built from",
        f"uses: [{text('simfbd@' + f['simfbd_commit'])}, "
        f"{text('rb@' + f['rb_commit'].rsplit('-g', 1)[-1])}]",
        "# what it read: its own config, the batch's config and specimens, and the template simfbd",
        "# carries. And what it wrote. simfbd status <run> checks them.",
        "inputs:",
    ]
    if "inputs" in f:
        lines += hash_lines(f["inputs"], "  ")
    else:
        lines.append(f"  {f['specimens']}: {text(digest(os.path.join(f['batch_dir'], f['specimens'])))}")
    lines += [
        "outputs:",
        f"  output: {text(digest(os.path.join(d, 'output')))}",
        "",
        "toolchain:",
    ]
    if "inputs" not in f:
        lines.append(f"  infer_script_hash: {text(f['infer_script_hash'])}")
    lines += [
        f"  rb_path: {text(f['rb_path'])}",
        f"  rb_md5: {text(f['rb_md5'])}",
        f"  rb_describe: {text(f['rb_commit'])}",
    ]
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

    run = doc.get("kind") == "run"
    batch_dir = (os.path.dirname(os.path.dirname(os.path.normpath(d))) or ".") if run else d
    last = (doc.get("simulated") or [{}])[-1]
    given = (doc.get("inputs") if run else last.get("inputs")) or {}
    one_hash = "config.yml" in given          # everything it read is under inputs, hashed alike

    for u in (doc.get("uses") if run else last.get("uses")) or []:
        row("uses", u)
    if one_hash:
        compare(d, {"config.yml": given["config.yml"]})
        compare(batch_dir, given.get("batch"), "batch: ")
        # what simfbd carries is named by its commit. A script that has since changed is not a
        # fault in the batch, so it is said and not counted.
        from . import tools
        for name, want in (given.get("simfbd") or {}).items():
            at = tools.SCRIPTS / name
            got = digest(at) if at.exists() else "gone"
            rows.append(("current" if got == want else "differs", f"simfbd: {name}",
                         "" if got == want else f"this simfbd's is {got}; it ran with {want}"))
    else:
        name = "config.yml" if os.path.exists(os.path.join(d, "config.yml")) else "config.sh"
        cfg = os.path.join(d, name)
        if os.path.exists(cfg):
            got = sha1_12(cfg)
            row("current" if got == doc.get("config_hash") else "changed", name,
                "" if got == doc.get("config_hash") else f"{doc.get('config_hash')} -> {got}")
        else:
            row("missing", name)

    if run:
        _, batch = find(batch_dir)
        read = given.get("batch") if one_hash else given
        if not one_hash:
            compare(batch_dir, given, "batch: ")
        compare(d, doc.get("outputs"))
        if not given and not doc.get("outputs"):
            row("unrecorded", "inputs and outputs", "written before hashes were recorded")
        if batch:
            if not one_hash and batch.get("config_hash") != doc.get("batch_config_hash"):
                row("changed", "batch: its config",
                    f"{doc.get('batch_config_hash')} -> {batch.get('config_hash')}")
            # the edge: what this run read must be what the batch says it wrote
            wrote = (batch.get("simulated") or [{}])[-1].get("outputs") or {}
            for name, got in (read or {}).items():
                if name in wrote and wrote[name] != got:
                    row("stale", "batch: " + name,
                        f"the batch has changed since this run was analyzed "
                        f"({doc.get('replicates')} -> {batch.get('replicates')} replicates)")
    else:
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
    if cmd == "get":
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
