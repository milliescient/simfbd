#!/usr/bin/env python3
"""The manifest module, where the shell driver and the R scripts beside this file look for it.

They call `manifest.py get|flat|write-batch|write-run|status DIR ...` and the Python scripts
here `import manifest`. Both get simfbd.manifest, found from this file's own place in the
package, so it works whether or not the package is installed for the interpreter running it.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from simfbd.manifest import *          # noqa: E402,F401,F403
from simfbd.manifest import main      # noqa: E402

if __name__ == "__main__":
    main(sys.argv[1:])
