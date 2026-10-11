"""The simfbd command.

    simfbd extend <batch> <n>          add n replicates to a batch, or give a new one its first
    simfbd fit <batch>/runs/<name>     fit every replicate of the batch not yet fitted
    simfbd status <batch or run> ...   does it still match its manifest
    simfbd summarize <run> ...         SBC ranks, coverage and diagnostics for a run
    simfbd ranks <run>                 mean rank per parameter
    simfbd config names [batch|run]    every setting a config can hold, and what it means
    simfbd config convert <config.sh>  an old shell config, as the YAML that replaces it
    simfbd script <name> [args]        run one of the other scripts simfbd carries
    simfbd path <name>                 where one of those scripts is, for something that sources it
    simfbd doctor                      the tools this machine would run, and whether they are there
    simfbd --version                   the commit of this code, -dirty if it has uncommitted changes

A batch is a directory holding a config.yml, anywhere on disk, and a run is a directory
under its runs/ holding one too. A batch made before configs were YAML holds a config.sh
instead. It can still be fitted, summarised and checked, by the shell driver it was made
with, which simfbd carries for that and for nothing else.
"""
import argparse
import os
import subprocess
import sys

from . import batch, config, fit, manifest, tools


def old_format(d):
    return os.path.exists(os.path.join(d, "config.sh")) and not os.path.exists(os.path.join(d, "config.yml"))


def shell_driver(script, args, extra=None):
    env = dict(tools.environment(), **(extra or {}))
    raise SystemExit(subprocess.run(["bash", str(tools.script(script)), *args], env=env).returncode)


def cmd_extend(a):
    d = os.path.abspath(a.batch)
    if "/runs/" in d + "/":
        raise SystemExit(f"error: {a.batch} is a run; extend its batch")
    if old_format(d):
        if manifest.find(d)[0] is None:
            raise SystemExit(f"error: {a.batch} has a config.sh and no replicates. A new batch "
                             f"takes a config.yml:\n       simfbd config convert {a.batch}/config.sh "
                             f"> {a.batch}/config.yml")
        shell_driver("extend.sh", [d, str(a.n)],
                     {"ALLOW_SIM_CHANGE": "1"} if a.allow_sim_change else None)
    batch.extend(d if os.path.isabs(a.batch) else a.batch, a.n, a.allow_sim_change)


def cmd_fit(a):
    run = a.run.rstrip("/")
    b = os.path.dirname(os.path.dirname(os.path.abspath(run)))
    if old_format(b):
        shell_driver("run.sh", [os.path.abspath(run)], {"FORCE": "1"} if a.force else None)
    if not fit.fit(run, force=a.force, slots=a.slots, slot_dir=a.slot_dir):
        raise SystemExit(1)


def cmd_status(a):
    raise SystemExit(1 if sum(manifest.status(p.rstrip("/")) for p in a.paths) else 0)


def run_script(name, args):
    path = tools.script(name)
    how = {".R": [tools.rscript()], ".py": [sys.executable], ".sh": ["bash"]}.get(path.suffix)
    if how is None:
        raise SystemExit(f"simfbd: {name} is not something to run; simfbd path {name} says where it is")
    raise SystemExit(subprocess.run([*how, str(path), *args], env=tools.environment()).returncode)


def cmd_config(a):
    if a.what == "names":
        for kind in ([a.arg] if a.arg in ("batch", "run") else ["batch", "run"]):
            print(f"# a {kind}'s config.yml")
            print("\n".join(config.names(kind)))
            print()
        return
    if not a.arg:
        raise SystemExit("usage: simfbd config convert <config.sh>")
    kind, values, notes = config.from_shell(a.arg)
    sys.stdout.write(f"# a {kind}'s config, converted from {os.path.basename(a.arg)}\n")
    sys.stdout.write(config.as_yaml(kind, values))
    for n in notes:
        sys.stderr.write(f"note: {n}\n")


def cmd_doctor(a):
    import shutil
    print(f"simfbd     {tools.version()}   {tools.home()}")
    ok = True
    for label, prog in (("rb", tools.rb()), ("Rscript", tools.rscript())):
        at = shutil.which(prog)
        ok = ok and bool(at)
        print(f"{label:<10} {at or 'NOT FOUND: ' + prog}")
    _, _, describe = tools.rb_identity()
    print(f"           rb was built from {describe}")
    r = subprocess.run([tools.rscript(), "-e",
                        'for (p in c("paleobuddy","posterior","SBC")) cat(sprintf("%-10s %s\\n", p, '
                        'if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))'],
                       capture_output=True, text=True)
    print(r.stdout.rstrip() or r.stderr.strip()[-200:])
    print(f"tools file {tools.TOOLS_FILE}" + ("" if tools.TOOLS_FILE.exists() else "  (none; $RBIN, $RSCRIPT or PATH are used)"))
    raise SystemExit(0 if ok and "NOT INSTALLED" not in r.stdout else 1)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="simfbd", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="store_true")
    ap.add_argument("--home", action="store_true", help="where this code is")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("extend", help="add replicates to a batch")
    p.add_argument("batch")
    p.add_argument("n", type=int)
    p.add_argument("--allow-sim-change", action="store_true",
                   help="add to a batch whose earlier replicates a different simulator made")
    p = sub.add_parser("fit", help="fit a run")
    p.add_argument("run")
    p.add_argument("--force", action="store_true", help="fit although a config has changed since the output was written")
    p.add_argument("--slots", type=int, help="cap chains across every fit sharing --slot-dir")
    p.add_argument("--slot-dir")
    p = sub.add_parser("status", help="does a batch or a run still match its manifest")
    p.add_argument("paths", nargs="+")
    p = sub.add_parser("summarize", help="SBC ranks, coverage and diagnostics")
    p.add_argument("runs", nargs="+")
    p = sub.add_parser("ranks", help="mean rank per parameter")
    p.add_argument("run")
    p = sub.add_parser("config", help="the settings a config can hold, or an old config converted")
    p.add_argument("what", choices=["names", "convert"])
    p.add_argument("arg", nargs="?")
    p = sub.add_parser("script", help="run one of the other scripts simfbd carries")
    p.add_argument("name")
    p.add_argument("args", nargs=argparse.REMAINDER)
    p = sub.add_parser("path", help="where one of those scripts is")
    p.add_argument("name")
    sub.add_parser("doctor", help="the tools this machine would run")

    a = ap.parse_args(argv)
    if a.version:
        print(tools.version())
    elif a.home:
        print(tools.home())
    elif a.cmd == "extend":
        cmd_extend(a)
    elif a.cmd == "fit":
        cmd_fit(a)
    elif a.cmd == "status":
        cmd_status(a)
    elif a.cmd == "summarize":
        run_script("summarize.R", a.runs)
    elif a.cmd == "ranks":
        from . import ranks
        ranks.main([a.run])
    elif a.cmd == "config":
        cmd_config(a)
    elif a.cmd == "script":
        run_script(a.name, a.args)
    elif a.cmd == "path":
        print(tools.script(a.name))
    elif a.cmd == "doctor":
        cmd_doctor(a)
    else:
        ap.print_help()
