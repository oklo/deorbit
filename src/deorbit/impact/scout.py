#!/usr/bin/env python3
"""Reduced-order SCOUT of the grazing-impact bounce-off phase diagram (pure stdlib).

Step 1 of the Paper-§2 campaign (docs/bounce_off_plan.md): cheaply map the regime
TOPOLOGY across the dimensionless (theta, v/v_esc) plane and emit the bracket
intervals that the expensive SPH bisection runs then refine. This is an
order-of-magnitude classifier built on the reduced-order model (`grazing.py`);
it is NOT the answer -- it tells the SPH campaign WHERE to look.

Regimes (for a fiducial strong km IRON body on rock):
  C  conventional crater        -- theta above the ricochet onset
  I  suborbital ricochet, INTACT -- grazing, retained v_out < v_esc, contact pressure < iron melt
  D  suborbital ricochet, DISPERSED -- as I but the leading contact shock-melts/fragments the body
  E  escape (v_out > v_esc)      -- at these speeds the contact always melts -> a dispersed fragment/vapor stream

Key physical couplings (grounded, see grazing.py + docs/bounce_off_plan.md):
  - escape  <=>  retained(theta) * v * cos(theta) > v_esc  (planet gravity via v_esc)
  - survival (intact vs dispersed)  <=>  P_contact = rho_t (v cos theta)^2  <  P_melt,iron (~220 GPa)
So an INTACT escape is essentially impossible for iron (escape needs v cos theta >~ v_esc/retained,
which drives P_contact far past melt) -- escaping material is a melted/fragmented stream (Pierazzo &
Melosh 2000; Schultz & Gault 1990 ring scenario). Cayambe sits in I; Orcus-like escape sits in E.

Run:  python3 -m deorbit.impact.scout            (ASCII map + results/bounce_off_scout.json)
"""
import math
import json
import os

from . import grazing

# escape speeds (m/s)
V_ESC = {"earth": 11186.0, "mars": 5027.0, "moon": 2376.0}
THETA_RIC_DEG = grazing.RICOCHET_ANGLE_DEG     # crater <-> ricochet onset (reduced-order; SPH refines)
P_MELT_IRON = 220e9                            # leading-contact shock-melt threshold

_BODY = grazing.Body()                         # 1 km iron sphere
_TGT = grazing.TARGETS["cayambe"]              # strong volcanic rock


def regime(theta_deg, vr, v_esc=V_ESC["earth"]):
    """Classify one (theta, v/v_esc) point -> 'C'|'I'|'D'|'E'."""
    v = vr * v_esc
    o = grazing.grazing_impact(v, theta_deg, _BODY, _TGT)
    if theta_deg >= THETA_RIC_DEG:
        return "C"
    v_out = o["retained_v_frac"] * v * math.cos(math.radians(theta_deg))
    melts = o["P_contact_GPa"] * 1e9 >= P_MELT_IRON
    if v_out > v_esc:
        return "E"
    return "D" if melts else "I"


def boundaries_along(vary, fixed, lo, hi, v_esc, n=80):
    """Find EVERY regime change along one axis (a line can cross >1 boundary, e.g.
    I->D->E), bisect each to the scout boundary x*, and return a widened SPH search
    bracket per transition with the correct adjacent regimes. vary='theta' or 'vr'."""
    def reg(x):
        return regime(x, fixed, v_esc) if vary == "theta" else regime(fixed, x, v_esc)
    xs = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    out = []
    margin = 0.5 if vary == "theta" else 0.25   # theta boundary (crater<->ricochet) is the most uncertain
    for i in range(1, n):
        ra, rb = reg(xs[i - 1]), reg(xs[i])
        if ra == rb:
            continue
        a, b = xs[i - 1], xs[i]
        for _ in range(40):
            m = 0.5 * (a + b)
            if reg(m) == ra:
                a = m
            else:
                b = m
        xstar = 0.5 * (a + b); w = margin * xstar
        out.append({"vary": vary, "fixed": fixed, "boundary": xstar,
                    "sph_bracket": [max(0.0, xstar - w), xstar + w], "from": ra, "to": rb})
    return out


def run(nt=24, nvr=22):
    thetas = [1.0 + (45.0 - 1.0) * i / (nt - 1) for i in range(nt)]
    vrs = [0.3 + (2.5 - 0.3) * j / (nvr - 1) for j in range(nvr)]
    grid = [[regime(t, vr) for t in thetas] for vr in vrs]   # grid[vr][theta]

    # bisection brackets to hand to the SPH campaign (Earth fiducial); each line may cross several:
    ve = V_ESC["earth"]
    brackets = []
    for vr in (0.5, 0.75, 1.0, 1.5, 2.0):           # vary theta: crater<->ricochet (+ dispersed sub-split)
        brackets += boundaries_along("theta", vr, 1.0, 44.0, ve)
    for th in (2.0, 5.0, 8.0, 12.0):                 # vary v/v_esc: intact<->dispersed<->escape
        brackets += boundaries_along("vr", th, 0.3, 2.5, ve)

    anchors = {
        "cayambe": {"planet": "earth", "v_kms": 7.4, "theta_deg": 2.0,
                    "vr": 7400.0 / V_ESC["earth"], "regime": regime(2.0, 7400.0 / V_ESC["earth"])},
        "orcus":   {"planet": "mars", "v_kms": 10.0, "theta_deg": 30.0,
                    "vr": 10000.0 / V_ESC["mars"], "regime": regime(30.0, 10000.0 / V_ESC["mars"], V_ESC["mars"])},
    }
    return {"thetas_deg": thetas, "vrs": vrs, "grid": grid, "brackets": brackets,
            "anchors": anchors, "theta_ric_deg": THETA_RIC_DEG, "v_esc": V_ESC,
            "note": "order-of-magnitude scout; boundaries are SPH-bisection brackets, not final"}


def ascii_map(res):
    th, vrs, grid = res["thetas_deg"], res["vrs"], res["grid"]
    lines = ["Bounce-off regime scout (Earth fiducial; rows=v/v_esc, cols=theta deg)",
             "  C=crater  I=intact ricochet  D=dispersed ricochet  E=escape"]
    for j in range(len(vrs) - 1, -1, -1):
        lines.append(f"{vrs[j]:4.2f} |" + "".join(grid[j]))
    lines.append("     +" + "-" * len(th))
    lines.append(f"      {th[0]:.0f}{' ' * (len(th) - 5)}{th[-1]:.0f}")
    return "\n".join(lines)


def main():
    res = run()
    print(ascii_map(res))
    print("\nanchors:")
    for k, a in res["anchors"].items():
        print(f"  {k:8s} {a['planet']:5s} v={a['v_kms']:.1f} km/s theta={a['theta_deg']:.0f}deg "
              f"v/v_esc={a['vr']:.2f} -> regime {a['regime']}")
    print("\nSPH bisection brackets to refine (scout boundary -> widened SPH search interval):")
    for b in res["brackets"]:
        lo, hi = b["sph_bracket"]
        if b["vary"] == "theta":
            print(f"  {b['from']}->{b['to']}: theta~{b['boundary']:.1f} deg (search [{lo:.1f},{hi:.1f}]) at v/v_esc={b['fixed']:.2f}")
        else:
            print(f"  {b['from']}->{b['to']}: v/v_esc~{b['boundary']:.2f} (search [{lo:.2f},{hi:.2f}]) at theta={b['fixed']:.0f} deg")
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
                       "results", "bounce_off_scout.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
