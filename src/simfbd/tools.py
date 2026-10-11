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


def rb():
    return os.environ.get("RBIN") or _configured().get("rb") or "rb"


def rscript():
    return os.environ.get("RSCRIPT") or _configured().get("rscript") or "Rscript"


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
