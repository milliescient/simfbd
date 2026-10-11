"""A batch's or a run's config.yml: what was asked for, in the names settings.py gives.

    timeline:
      rate_intervals: 8
      interval_width: 1.5
    origin_age:
      min: 15
      max: 21

A name that is not a setting is refused, with the nearest one that is, so a typo cannot
quietly leave a default in place. A batch does not say how many replicates it has: `simfbd
extend` adds them and the manifest counts them.
"""
import difflib
import hashlib
import os
import subprocess

import yaml

from . import settings as S

KINDS = {"batch": (S.BATCH_SETTINGS, S.BATCH_REQUIRED, {}),
         "run": (S.RUN_SETTINGS, S.RUN_REQUIRED, S.RUN_DEFAULTS)}


def file_hash(path):
    with open(path, "rb") as fh:
        return hashlib.sha1(fh.read()).hexdigest()[:12]


def _flatten(doc, prefix=""):
    for k, v in (doc or {}).items():
        name = f"{prefix}{k}"
        if isinstance(v, dict):
            yield from _flatten(v, name + ".")
        else:
            yield name, v


def _word(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float) and v == int(v) and abs(v) < 1e15:
        return repr(v)                  # 1.0 stays 1.0, as it was written
    return str(v)


def load(path, kind):
    """{internal name: value as the scripts take it}, validated."""
    table, required, defaults = KINDS[kind]
    try:
        doc = yaml.safe_load(open(path)) or {}
    except yaml.YAMLError as e:
        raise SystemExit(f"simfbd: {path} is not valid YAML\n  {e}")
    if not isinstance(doc, dict):
        raise SystemExit(f"simfbd: {path} must be settings under their section names")
    by_name = {name: key for key, name, _ in table}
    out = dict(defaults)
    for name, v in _flatten(doc):
        if name not in by_name:
            other = {n: k for k, n, _ in KINDS["run" if kind == "batch" else "batch"][0]}
            if name in other:
                raise SystemExit(f"simfbd: {path} sets {name}, which belongs to the "
                                 f"{'run' if kind == 'batch' else 'batch'}, not the {kind}")
            near = difflib.get_close_matches(name, by_name, n=1)
            raise SystemExit(f"simfbd: {path} sets {name}, which is not a setting"
                             + (f". Did you mean {near[0]}?" if near else "")
                             + "\n  simfbd config names lists them")
        if v is None:
            continue
        key = by_name[name]
        out[key] = _word(v)
        if key in S.CHOICES and out[key] not in S.CHOICES[key]:
            raise SystemExit(f"simfbd: {path}: {name} is {out[key]}; it is one of "
                             f"{', '.join(S.CHOICES[key])}")
    name_of = {key: name for key, name, _ in table}
    missing = [name_of[k] for k in required if k not in out]
    if missing:
        raise SystemExit(f"simfbd: {path} does not set {', '.join(missing)}")
    return out


def names(kind):
    table, required, defaults = KINDS[kind]
    lines = []
    for key, name, meaning in table:
        mark = "required" if key in required else f"default {defaults[key]}" if key in defaults else ""
        lines.append(f"{name:<68} {mark:<14} {meaning}".rstrip())
    return lines


def from_shell(path):
    """A config.sh as (kind, {internal name: value}, what was left out and why)."""
    script = 'set -a; source "$1" >/dev/null 2>&1; env -0'
    base = subprocess.run(["bash", "-c", "env -0"], capture_output=True).stdout
    got = subprocess.run(["bash", "-c", script, "_", os.path.abspath(path)], capture_output=True).stdout
    before = dict(e.split("=", 1) for e in base.decode().split("\0") if "=" in e)
    after = dict(e.split("=", 1) for e in got.decode().split("\0") if "=" in e)
    new = {k: v for k, v in after.items() if before.get(k) != v and k not in ("_", "SHLVL", "PWD", "OLDPWD")}
    kind = "run" if "MODEL" in new else "batch"
    known = {key for key, _, _ in KINDS[kind][0]}
    notes = []
    if "NREPS" in new:
        notes.append(f"NREPS={new['NREPS']} is not a setting: a batch is given replicates with "
                     f"simfbd extend <batch> {new['NREPS']}")
    if "SIM_TIMEOUT" in new:
        notes.append(f"SIM_TIMEOUT={new['SIM_TIMEOUT']} counted seconds and is gone. Set "
                     f"simulation_limits.max_waiting_time_draws: 100000 refuses about what 20 "
                     f"seconds did, 200000 what 60 did")
    for k in sorted(set(new) - known - {"NREPS", "SIM_TIMEOUT"}):
        notes.append(f"{k}={new[k]} is not a {kind} setting and was left out")
    return kind, {k: v for k, v in new.items() if k in known}, notes


def as_yaml(kind, values):
    from .manifest import settings_lines
    return "\n".join(settings_lines(KINDS[kind][0], values.get)) + "\n"
