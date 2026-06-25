"""Invariant: the corridor/site/impact core is PURE STDLIB.

The 24/7 local daemon and the cloud routine run this core with no install and no
third-party packages. This test fails loudly if a refactor ever pulls numpy (or
the rebound/assist high-fidelity stack) into the core import path.
"""
import sys


def test_core_has_no_heavy_deps():
    for mod in ("numpy", "rebound", "assist"):
        sys.modules.pop(mod, None)
    import deorbit.corridor.physics   # noqa: F401
    import deorbit.corridor.sweep     # noqa: F401
    import deorbit.site_selection     # noqa: F401
    import deorbit.impact.grazing     # noqa: F401
    leaked = [m for m in ("numpy", "rebound", "assist") if m in sys.modules]
    assert not leaked, f"stdlib core unexpectedly imported heavy deps: {leaked}"
