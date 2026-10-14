"""The programs simfbdr runs, and which version of itself it is.

rb and Rscript are looked for in the environment first ($RBIN, $RSCRIPT), then in
~/.config/simfbdr/tools.yml, then on PATH:

    rb:
      - /opt/revbayes/build/rb
      - /opt/revbayes/build-older/rb
    rscript: /opt/R/bin/Rscript

So a machine says once where its builds are, and a single run can still be pointed at
another without editing anything.

A machine can list several builds of rb, and each is known by the commit its binary carries.
Which of them analyzes a run is choose_rb()'s to say: the one the run's config states, else the
one that analyzed the run's earlier replicates, else the first listed.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE / "scripts"
TOOLS_FILE = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "simfbdr" / "tools.yml"


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


def _listed():
    """The builds of rb the tools file names, in its order."""
    got = _configured().get("rb") or []
    return [str(b) for b in ([got] if isinstance(got, str) else got)]


def _rb():
    """(the program, which rule supplied it), before a run has had its say."""
    if os.environ.get("RBIN"):
        return os.environ["RBIN"], "from $RBIN"
    if _chosen:
        return _chosen
    if _listed():
        return _listed()[0], f"from {TOOLS_FILE}" + (", the first listed" if len(_listed()) > 1 else "")
    return "rb", "from PATH"


def rb():
    return _rb()[0]


def rb_builds():
    """Every build of rb this machine has, as (path, the commit it was built from): the ones
    the tools file lists, then the one on PATH if it is not among them."""
    paths = [which(b) for b in _listed()]
    on_path = shutil.which("rb")
    if on_path and on_path not in paths:
        paths.append(on_path)
    return [(p, rb_identity(p)[2].rsplit("-g", 1)[-1]) for p in paths]


def same_version(have, want):
    """A commit may be abbreviated differently in two places, so two hex strings match when
    one begins the other. Anything else is compared whole, and -dirty never matches."""
    import re
    hexes = all(re.fullmatch(r"[0-9a-f]{6,40}", v) for v in (have, want))
    return have == want or (hexes and (have.startswith(want) or want.startswith(have)))


_chosen = None


def choose_rb(stated=None, recorded=None, config="", run="", force=False):
    """Settle which build of rb analyzes a run.

    The version its config states under `uses` wins, and is refused when no build here is
    that one. Otherwise the build that analyzed the run's earlier replicates, so one run's logs
    all come from one rb; when that build is gone it is refused, unless forced. Otherwise the
    first listed. $RBIN overrides all of it for one command, and require() still holds it to
    what the config states.
    """
    global _chosen
    _seen.pop("rb", None)
    if os.environ.get("RBIN"):
        return
    for want, why in ((stated, f"{config} states rb@{stated}"),
                      (recorded, f"rb@{recorded} analyzed the earlier replicates of {run}")):
        if not want or want == "unknown":
            continue
        builds = rb_builds()
        for path, commit in builds:
            if same_version(commit, want):
                _chosen = (path, why)
                return
        have = "\n".join(f"         rb@{c}  {p}" for p, c in builds) or "         none"
        if want is stated:
            raise SystemExit(f"error: {why}, and no build of rb here is that one. This machine has\n"
                             f"{have}\n       List the build in {TOOLS_FILE}.")
        if not force:
            raise SystemExit(f"error: {why}, and no build of rb here is that one. "
                             f"This machine has\n{have}\n"
                             f"       Analyzing the rest with another would mix two builds in one run. "
                             f"List the build in {TOOLS_FILE},\n       or pass --force to use the first listed.")


def rscript():
    return _find("RSCRIPT", "rscript", "Rscript")[0]


def found():
    """Where each tool was found and why, for a person to read before a long run."""
    return [(name, shutil.which(prog) or f"NOT FOUND: {prog}", rule)
            for name, (prog, rule) in (("rb", _rb()),
                                       ("Rscript", _find("RSCRIPT", "rscript", "Rscript")))]


def announce():
    for name, at, rule in found():
        print(f"  {name:<8} {at}  ({rule})", flush=True)


def environment():
    """What a script simfbdr launches is given: the tools, and this interpreter, which is
    the one that has simfbdr and its YAML library."""
    return dict(os.environ, RBIN=rb(), RSCRIPT=rscript(), PY=sys.executable)


def script(name):
    path = SCRIPTS / name
    if not path.exists():
        have = ", ".join(sorted(p.name for p in SCRIPTS.iterdir() if p.is_file()))
        raise SystemExit(f"simfbdr: no script {name}. There are: {have}")
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


# what a batch is simulated with, and what a run is analyzed with
SIMULATES = ("simfbdr", "R", "paleobuddy")
FITS = ("simfbdr", "rb")
KNOWN = ("simfbdr", "rb", "R", "paleobuddy")
_seen = {}


def version_of(name):
    """The version of one of the things simfbdr runs, as a word a config or a notebook can
    name: a commit for simfbdr and rb, a version number for R and paleobuddy."""
    if name in _seen:
        return _seen[name]
    if name == "simfbdr":
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
        raise SystemExit(f"simfbdr: {name} is not something simfbdr runs; it runs {', '.join(KNOWN)}")
    _seen[name] = v
    return v


def uses(names):
    return [f"{n}@{version_of(n)}" for n in names]


def require(stated, where):
    """Refuse unless everything a config names under `uses` is here at that version."""
    bad = []
    for u in stated:
        name, at, want = str(u).partition("@")
        if name not in KNOWN or not at or not want:
            raise SystemExit(f"error: {where}: uses {u}; write one of {', '.join(KNOWN)} as name@version")
        have = version_of(name)
        if not same_version(have, want):
            bad.append(f"{name} is at {have} here, and {where} states {name}@{want}")
    if bad:
        raise SystemExit("error: " + "\n       ".join(bad))


_identities = {}


def rb_identity(path=None):
    """(path, md5 of the binary, the git describe it carries), for the rb that would run or
    for the build at a path. The binary names what actually ran even after its source tree
    moves on; an md5 identifies a build, not a commit."""
    import hashlib
    import re
    path = which(path or rb())
    if path not in _identities:
        try:
            data = Path(path).read_bytes()
        except OSError:
            return path, "unknown", "unknown"
        md5 = hashlib.md5(data).hexdigest()[:12]
        m = re.search(rb"(?<![\x20-\x7e])([A-Za-z0-9._-]+-[0-9]+-g[0-9a-f]{6,})(?![\x20-\x7e])", data)
        _identities[path] = (path, md5, m.group(1).decode() if m else "unknown")
    return _identities[path]
