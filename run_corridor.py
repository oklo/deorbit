#!/usr/bin/env python3
"""Entry point for the capture-corridor sweep (kept at the repo root, named so it
does NOT collide with the ``deorbit`` package under ``src/``).

The implementation lives at ``src/deorbit/corridor/sweep.py``. Usage is unchanged:

    python3 run_corridor.py --daemon          # local 24/7 daemon (see daemons.json)
    python3 run_corridor.py --increment 1500  # bounded step (cloud routine)

Pure stdlib, no install: it puts ``src/`` on the path and calls
``deorbit.corridor.sweep.run``. ``run`` derives its state dir from ``sys.argv[0]``
(= this file at the repo root), so state stays at ``<repo>/state``.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from deorbit.corridor.sweep import run  # noqa: E402

if __name__ == "__main__":
    run()
