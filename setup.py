"""Only here to write down the commit a copy is built from.

Run from a checkout, simfbdr asks git which commit it is. An installed copy has no checkout to
ask, so the commit is written into the package as it is built, and `simfbdr --version` says
the same word either way.
"""
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


class build_with_commit(build_py):
    def run(self):
        super().run()
        here = Path(__file__).resolve().parent
        r = subprocess.run(["git", "-C", str(here), "describe", "--always", "--dirty", "--tags"],
                           capture_output=True, text=True)
        if r.returncode == 0 and not self.editable_mode:
            out = Path(self.build_lib) / "simfbdr" / "_commit.txt"
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(r.stdout.strip() + "\n")


setup(cmdclass={"build_py": build_with_commit})
