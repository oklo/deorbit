#!/usr/bin/env python3
"""Validate the bounce-off phase diagram against real terrestrial IRON impactors (campaign step 5).

Two best-characterized cases bracketing the regimes we map:
  Barringer / Meteor Crater, AZ  -- Canyon Diablo iron, ~50 m, ~12.8 km/s, theta~45 deg
      (Artemieva & Pierazzo 2009/2011; Melosh & Collins 2005). v/v_esc~1.14, theta>>crater onset
      -> diagram predicts CRATER (C). Known outcome: 1.2 km crater, projectile melted/fragmented.
  Campo del Cielo, Argentina     -- IAB iron, very shallow entry ~9 deg, fragments at a few km/s,
      18.5 km crater-chain strewn field. v/v_esc~0.3 -> diagram predicts the intact/dispersed
      grazing corner.

CAVEATS (honest):
 - The SPH setup is SCALE-FREE here (no gravity/atmosphere; strength negligible at these speeds),
   so these test the REGIME prediction in (theta, v/v_esc), not absolute crater size. We run at the
   campaign scale (D=1000 m); a 50 m run would give the same regime. Absolute crater size needs
   gravity + the real target + longer integration (a follow-on).
 - Campo del Cielo's crater CHAIN is an atmospheric-fragmentation product; our vacuum monolith tests
   "what a coherent shallow iron body does," not the strewn-field mechanism.
 - Target here is basalt (the campaign material); the real targets were sandstone (Barringer) and
   loess/sediment (Campo) -- softer. Another reason this is a regime check, not a size match.

Run:  python3 -m deorbit.impact.validate_craters   (-> state/validation/, prints regime vs prediction)
"""
import os
import sys
import json
import math
import subprocess

from . import boundary_daemon as bd

# name -> (theta_deg, v_ms, x_impact_frac for steep IC | None for grazing, predicted regime, note)
CASES = {
    "barringer": (45.0, 12800.0, -0.10, "C",
                  "Canyon Diablo iron ~50 m, ~12.8 km/s, ~45 deg; real: 1.2 km crater, projectile destroyed"),
    "campo_del_cielo": (9.0, 3468.0, None, "I/D",
                        "IAB iron, ~9 deg shallow entry, ~3.5 km/s fragments; real: 18.5 km crater chain (atmospheric breakup)"),
}


def run_case(name, theta, v, x_impact_frac, statedir):
    rundir = os.path.join(statedir, name)
    os.makedirs(rundir, exist_ok=True)
    ic = os.path.join(rundir, "ic.bin")
    dx, N = bd.build_grazing_ic(theta, v, ic, bd.SPH["cppr"], bd.SPH["D"], x_impact_frac=x_impact_frac)
    t_end = 25e3 / max(v * math.cos(math.radians(theta)), 1.0)
    gpu = os.path.join(bd._SPH, "gpu_ic")
    subprocess.run([gpu, ic, f"{dx}", f"{t_end}", "--snap", f"{0.6 * t_end}"], cwd=rundir, check=True,
                   stdout=open(os.path.join(rundir, "run.log"), "w"), stderr=subprocess.STDOUT)
    snaps = sorted(f for f in os.listdir(rundir) if f.startswith("gpu_snap_") and f.endswith(".bin"))
    peak_melt = 0.0
    for line in open(os.path.join(rundir, "run.log")):
        i = line.find("melt=")
        if i >= 0:
            try:
                peak_melt = max(peak_melt, float(line[i + 5:line.index("%", i)]) / 100.0)
            except ValueError:
                pass
    reg, diag = bd.classify_regime(os.path.join(rundir, snaps[-1]), bd.SPH["v_esc"], bd.SPH["imp_mat"], peak_melt, v)
    out = {"case": name, "theta": theta, "v_ms": v, "vr": v / bd.SPH["v_esc"], "regime": reg, **diag}
    json.dump(out, open(os.path.join(rundir, "result.json"), "w"), indent=1)
    for f in os.listdir(rundir):
        if f.endswith(".bin"):
            os.remove(os.path.join(rundir, f))
    return out


def main():
    sd = os.path.join(bd._REPO, "state", "validation")
    print(f"validating against real iron impactors (CPPR={bd.SPH['cppr']}, D={bd.SPH['D']:.0f} m scale)\n")
    for name, (theta, v, xf, pred, note) in CASES.items():
        print(f"[{name}] theta={theta} deg, v={v/1e3:.1f} km/s (v/v_esc={v/bd.SPH['v_esc']:.2f}); predicted {pred}\n  {note}")
        try:
            o = run_case(name, theta, v, xf, sd)
            ok = "MATCH" if o["regime"] in pred else "MISMATCH"
            print(f"  -> SPH regime: {o['regime']}  [{ok} vs predicted {pred}]   "
                  f"v_bulk={o['v_bulk_kms']:.2f} km/s, disp={o['disp_kms']:.2f}, "
                  f"z_mean={o['z_mean_m']:.0f} m, peak_melt={100*(o['peak_melt'] or 0):.0f}%\n")
        except Exception as e:
            print(f"  -> FAILED: {e}\n")


if __name__ == "__main__":
    main()
