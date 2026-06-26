#!/usr/bin/env python3
"""Entry point for the bounce-off boundary-locating daemon (kept at the repo root so it
does NOT collide with the `deorbit` package, like run_corridor.py).

    python3 run_bounce.py --daemon             # 24/7 local SPH-bisection of the regime boundaries
    python3 run_bounce.py --increment 1800     # bounded step (seconds)
    python3 run_bounce.py --increment 5 --score reduced   # fast self-test of the driver

Implementation: src/deorbit/impact/boundary_daemon.py. LOCAL-only (needs the Mac GPU + Metal
for the SPH score). State -> <repo>/state/bounce/.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from deorbit.impact.boundary_daemon import run  # noqa: E402

if __name__ == "__main__":
    run()
