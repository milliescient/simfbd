"""The programs simfbd runs, and which version of itself it is.

rb and Rscript are looked for in the environment first ($RBIN, $RSCRIPT), then in
~/.config/simfbd/tools.yml, then on PATH:

    rb: /opt/revbayes/build/rb
    rscript: /opt/R/bin/Rscript

So a machine says once where its builds are, and a single run can still be pointed at
another without editing anything.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE / "scripts"
TOOLS_FILE = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "simfbd" / "tools.yml"


def _configured():
    try:
        return yaml.safe_load(TOOLS_FILE.read_text()) or {}
    except OSError:
        return {}


def _find(var, key, fallback):
    """(the program, which rule supplied it)."""
    if os.environ.get(var):
        return os.environ[var], f"from ${var}"
    if _configured().get(key):
        return str(_configured()[key]), f"from {TOOLS_FILE}"
    return fallback, "from PATH"


def rb():
    return _find("RBIN", "rb", "rb")[0]


def rscript():
    return _find("RSCRIPT", "rscript", "Rscript")[0]


def found():
    """Where each tool was found and why, for a person to read before a long run."""
    return [(name, shutil.which(prog) or f"NOT FOUND: {prog}", rule)
            for name, (prog, rule) in (("rb", _find("RBIN", "rb", "rb")),
                                       ("Rscript", _find("RSCRIPT", "rscript", "Rscript")))]


def announce():
    for name, at, rule in found():
        print(f"  {name:<8} {at}  ({rule})", flush=True)


def environment():
    """What a script simfbd launches is given: the tools, and this interpreter, which is
    the one that has simfbd and its YAML library."""
    return dict(os.environ, RBIN=rb(), RSCRIPT=rscript(), PY=sys.executable)


def script(name):
    path = SCRIPTS / name
    if not path.exists():
        have = ", ".join(sorted(p.name for p in SCRIPTS.iterdir() if p.is_file()))
        raise SystemExit(f"simfbd: no script {name}. There are: {have}")
    return path


def home():
    """The checkout this is running from, or the installed package when there is none."""
    top = _git("rev-parse", "--show-toplevel")
    return top or str(HERE)


def _git(*args):
    try:
        r = subprocess.run(["git", "-C", str(HERE), *args], capture_output=True, text=True)
    except OSError:
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def version():
    """The commit of this code, ending in -dirty when it has uncommitted changes. Run from a
    checkout it is asked of git, and an installed copy carries the commit it was built from,
    so the word is the same either way and a manifest or a notebook can name it."""
    inside = _git("rev-parse", "--is-inside-work-tree") == "true" and _git("ls-files", "--error-unmatch", "tools.py")
    if inside:
        return _git("describe", "--always", "--dirty", "--tags") or "none"
    baked = HERE / "_commit.txt"
    return baked.read_text().strip() if baked.exists() else "none"


def dirty():
    return version().endswith("-dirty")


def which(program):
    return shutil.which(program) or program


# what a batch is simulated with, and what a run is fitted with
SIMULATES = ("simfbd", "R", "paleobuddy")
FITS = ("simfbd", "rb")
KNOWN = ("simfbd", "rb", "R", "paleobuddy")
_seen = {}


def version_of(name):
    """The version of one of the things simfbd runs, as a word a config or a notebook can
    name: a commit for simfbd and rb, a version number for R and paleobuddy."""
    if name in _seen:
        return _seen[name]
    if name == "simfbd":
        v = version()
    elif name == "rb":
        v = rb_identity()[2].rsplit("-g", 1)[-1]
    elif name in ("R", "paleobuddy"):
        r = subprocess.run([rscript(), "-e", 'cat(as.character(getRversion()), "\n"); '
                            'cat(if (requireNamespace("paleobuddy", quietly = TRUE)) '
                            'as.character(packageVersion("paleobuddy")) else "missing", "\n")'],
                           capture_output=True, text=True)
        got = r.stdout.split()
        _seen["R"] = got[0] if len(got) > 0 else "unknown"
        _seen["paleobuddy"] = got[1] if len(got) > 1 else "unknown"
        return _seen[name]
    else:
        raise SystemExit(f"simfbd: {name} is not something simfbd runs; it runs {', '.join(KNOWN)}")
    _seen[name] = v
    return v


def uses(names):
    return [f"{n}@{version_of(n)}" for n in names]


def require(stated, where):
    """Refuse unless everything a config names under `uses` is here at that version.

    A commit may be abbreviated differently in the two places, so two hex strings match when
    one begins the other. Anything else is compared whole, and -dirty never matches.
    """
    import re
    bad = []
    for u in stated:
        name, at, want = str(u).partition("@")
        if name not in KNOWN or not at or not want:
            raise SystemExit(f"error: {where}: uses {u}; write one of {', '.join(KNOWN)} as name@version")
        have = version_of(name)
        hexes = all(re.fullmatch(r"[0-9a-f]{6,40}", v) for v in (have, want))
        if not (have == want or (hexes and (have.startswith(want) or want.startswith(have)))):
            bad.append(f"{name} is at {have} here, and {where} states {name}@{want}")
    if bad:
        raise SystemExit("error: " + "\n       ".join(bad))


def rb_identity():
    """(path, md5 of the binary, the git describe it carries). The binary names what
    actually ran even after its source tree moves on; an md5 identifies a build, not a commit."""
    import hashlib
    import re
    path = which(rb())
    try:
        data = Path(path).read_bytes()
    except OSError:
        return path, "unknown", "unknown"
    md5 = hashlib.md5(data).hexdigest()[:12]
    m = re.search(rb"(?<![\x20-\x7e])([A-Za-z0-9._-]+-[0-9]+-g[0-9a-f]{6,})(?![\x20-\x7e])", data)
    return path, md5, m.group(1).decode() if m else "unknown"
