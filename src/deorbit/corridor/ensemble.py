"""Atmosphere/drag uncertainty ensemble for the capture corridor (pure stdlib).

The deterministic sweep ([[sweep]]) draws a single capture boundary. This module
asks the harder question the corridor's narrowness demands: how much of that
boundary is real vs. an artifact of uncertain inputs? It Monte-Carlos the
physically-relevant uncertainties and reports the corridor as a PROBABILITY band
plus a sensitivity ranking.

GROUNDING (verified from the deterministic corridor): a beta~5e6 km iron body
sheds its capture energy at periapsis ~2-45 km -- the stratosphere / lower
mesosphere -- NOT the thermosphere. So the factor-of-several thermospheric solar-
cycle density variability is irrelevant here; the body never samples it. The
relevant levers are:
  * the body's drag coefficient Cd (orientation/tumbling), since beta ~ 1/Cd and
    drag energy loss ~ 1/beta -- the most direct lever;
  * the lower/middle-atmosphere density (~+-20%, well constrained: latitude,
    season, weather, sudden stratospheric warmings);
  * Earth co-rotation (prograde vs retrograde x inclination, +-465 m/s eq.),
    which also fixes the inertial-velocity simplification of the base model;
  * scale height / temperature (a few %, expected minor).

Run:  python3 -m deorbit.corridor.ensemble --n 200            (writes results/corridor_ensemble.json)
"""
import math
import json
import os
import random
import argparse

from . import physics

# --- priors (best-practice, grounded; see module docstring) ---
DENS_SIGMA_LN = 0.20            # lower/mid-atmosphere density factor: lognormal, ~+-20%
H_SIGMA = 0.05                  # scale-height (temperature) factor: normal, +-5%
CD_MEDIAN, CD_SIGMA_LN = 1.1, 0.25   # drag coefficient: lognormal (sphere ~1.0; tumbling/irregular ->~2)
CD_CLIP = (0.7, 2.2)

# grid: captures live at periapsis<~45 km; v_inf extended to 6 km/s to resolve the upper tail
V_LO, V_HI = 0.2e3, 6.0e3
P_LO, P_HI = 0.0, 60e3

# Moon-survival filter (apogee-based proxy, calibrated to the high-fidelity REBOUND/ASSIST
# study, see deorbit.highfidelity): a 2-body+drag "capture" only stays bound + comes down to
# Earth if its FIRST APOGEE is below the lunar-ejection scale; below ~130,000 km it aerobrakes
# gently to a grazing impact (the Cayambe scenario), between that and ~Earth's Hill radius the
# Moon pumps it to a steep impact (still bound), and beyond Hill it is ejected. This is a
# phase-averaged proxy for the per-case 3-body result; the boundary is phase-dependent ~1.5-2e6 km.
APO_GENTLE_KM = 130_000.0      # first apogee below this -> gentle aerobrake -> grazing impact
APO_EJECT_KM = 1.5e6           # ~Earth Hill radius; above this the Moon ejects the capture


def _grid(nv, npa):
    vs = [V_LO + (V_HI - V_LO) * i / (nv - 1) for i in range(nv)]
    ps = [P_LO + (P_HI - P_LO) * j / (npa - 1) for j in range(npa)]
    return vs, ps


def sample_params(rng):
    """Draw one atmosphere+body+geometry realization from the priors."""
    f_rho = math.exp(rng.gauss(0.0, DENS_SIGMA_LN))
    f_H = max(0.7, 1.0 + rng.gauss(0.0, H_SIGMA))
    cd = min(CD_CLIP[1], max(CD_CLIP[0], CD_MEDIAN * math.exp(rng.gauss(0.0, CD_SIGMA_LN))))
    sgn = 1.0 if rng.random() < 0.5 else -1.0          # prograde(+) / retrograde(-)
    inc = math.acos(rng.random())                       # isotropic inclination: cos(i) ~ U[0,1]
    omega = sgn * physics.OMEGA_EARTH * math.cos(inc)
    return {"f_rho": f_rho, "f_H": f_H, "cd": cd, "omega": omega}


CLASSES = ("capture", "bound", "gentle")   # capture: eps<0; bound: survives the Moon; gentle: grazing


def corridor_metrics(vs, ps, cd=1.0, f_rho=1.0, f_H=1.0, omega=0.0):
    """Classify the grid for one realization. Each captured cell is further filtered by
    first apogee into 'bound' (survives the Moon, apogee<Hill) and 'gentle' (apogee<130k km,
    aerobrakes to a grazing impact). Returns (masks, n, vmax) dicts keyed by CLASSES."""
    nv, npa = len(vs), len(ps)
    body = physics.Body(cd=cd)
    masks = {c: [[False] * npa for _ in range(nv)] for c in CLASSES}
    n = {c: 0 for c in CLASSES}
    vmax = {c: 0.0 for c in CLASSES}
    for i, v in enumerate(vs):
        seen = {c: False for c in CLASSES}
        for j, p in enumerate(ps):
            r = physics.simulate_passage(v, p, body, f_rho=f_rho, f_H=f_H, omega=omega)
            if r["outcome"] != "capture":
                continue
            apo = r.get("apoapsis_km", float("inf"))
            hits = ["capture"]
            if apo < APO_EJECT_KM:
                hits.append("bound")
            if apo < APO_GENTLE_KM:
                hits.append("gentle")
            for c in hits:
                masks[c][i][j] = True
                n[c] += 1
                seen[c] = True
        for c in CLASSES:
            if seen[c] and v > vmax[c]:
                vmax[c] = v
    return masks, n, {c: vmax[c] / 1e3 for c in CLASSES}


def _pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * len(xs)))]


def run_ensemble(n=200, nv=30, npa=30, seed=0):
    vs, ps = _grid(nv, npa)
    rng = random.Random(seed)

    # nominal (deterministic, non-rotating, Cd=1.0) -- validation anchor
    _, nom_n, nom_vmax = corridor_metrics(vs, ps)

    pmap = {c: [[0] * npa for _ in range(nv)] for c in CLASSES}
    vmaxes = {c: [] for c in CLASSES}
    for _ in range(n):
        pr = sample_params(rng)
        masks, _, vmax = corridor_metrics(vs, ps, pr["cd"], pr["f_rho"], pr["f_H"], pr["omega"])
        for c in CLASSES:
            for i in range(nv):
                for j in range(npa):
                    if masks[c][i][j]:
                        pmap[c][i][j] += 1
            vmaxes[c].append(vmax[c])
    pmap = {c: [[v / n for v in row] for row in pmap[c]] for c in CLASSES}

    # single-axis variance share on the physically-complete (Moon-surviving) corridor's vmax
    def single_axis(axis, m=40):
        r2 = random.Random(seed + 1)
        out = []
        for _ in range(m):
            pr = {"cd": 1.0, "f_rho": 1.0, "f_H": 1.0, "omega": 0.0}
            pr[axis] = sample_params(r2)[axis]
            out.append(corridor_metrics(vs, ps, pr["cd"], pr["f_rho"], pr["f_H"], pr["omega"])[2]["bound"])
        return out

    sens = {}
    for axis in ("cd", "f_rho", "omega", "f_H"):
        vals = single_axis(axis)
        mu = sum(vals) / len(vals)
        sd = (sum((x - mu) ** 2 for x in vals) / len(vals)) ** 0.5
        sens[axis] = {"vmax_std_kms": sd}
    tot_var = sum(s["vmax_std_kms"] ** 2 for s in sens.values()) or 1.0
    for s in sens.values():
        s["variance_share"] = s["vmax_std_kms"] ** 2 / tot_var

    return {
        "priors": {"dens_sigma_ln": DENS_SIGMA_LN, "h_sigma": H_SIGMA,
                   "cd_median": CD_MEDIAN, "cd_sigma_ln": CD_SIGMA_LN, "cd_clip": CD_CLIP},
        "moon_filter": {"apo_gentle_km": APO_GENTLE_KM, "apo_eject_km": APO_EJECT_KM},
        "grid": {"nv": nv, "npa": npa, "v_kms": [V_LO / 1e3, V_HI / 1e3], "peri_km": [P_LO / 1e3, P_HI / 1e3]},
        "n_samples": n,
        "nominal": {"n_capture_cells": nom_n["capture"], "vmax_capture_kms": nom_vmax["capture"],
                    "vmax_bound_kms": nom_vmax["bound"], "vmax_gentle_kms": nom_vmax["gentle"]},
        "envelope_vmax_kms": {c: {"mean": sum(vmaxes[c]) / n, "p05": _pct(vmaxes[c], 0.05),
                                  "p50": _pct(vmaxes[c], 0.5), "p95": _pct(vmaxes[c], 0.95),
                                  "min": min(vmaxes[c]), "max": max(vmaxes[c])} for c in CLASSES},
        "sensitivity_vmax_bound": sens,
        "vs_kms": [v / 1e3 for v in vs],
        "ps_km": [p / 1e3 for p in ps],
        "p_capture": pmap["capture"],
        "p_bound": pmap["bound"],
        "p_gentle": pmap["gentle"],
    }


def ascii_map(res, key, title):
    """Terminal heatmap of a probability grid: rows = periapsis (top=high), cols = v_inf."""
    pcap, vs, ps = res[key], res["vs_kms"], res["ps_km"]
    nv, npa = len(vs), len(ps)
    chars = " .:+#"
    lines = [f"{title}  (N={res['n_samples']}; rows=periapsis km, cols=v_inf km/s)"]
    for j in range(npa - 1, -1, -1):
        row = "".join(chars[min(4, int(pcap[i][j] * 5))] for i in range(nv))
        lines.append(f"{ps[j]:5.1f} |{row}")
    lines.append("      +" + "-" * nv)
    lines.append(f"       {vs[0]:.1f}{' ' * (nv - 6)}{vs[-1]:.1f}")
    lines.append("  legend: ' '<.05  '.'<.25  ':'<.5  '+'<.75  '#'>=.75")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="capture-corridor atmosphere/drag uncertainty ensemble")
    ap.add_argument("--n", type=int, default=200, help="Monte-Carlo samples")
    ap.add_argument("--nv", type=int, default=30)
    ap.add_argument("--npa", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
        "results", "corridor_ensemble.json"))
    a = ap.parse_args()
    res = run_ensemble(a.n, a.nv, a.npa, a.seed)

    print(ascii_map(res, "p_capture", "P(capture)  [2-body+drag]"))
    print()
    print(ascii_map(res, "p_bound", "P(bound capture)  [after Moon-survival filter]"))
    nm = res["nominal"]
    print(f"\nnominal (Cd=1, no rotation): vmax km/s  capture={nm['vmax_capture_kms']:.2f}  "
          f"bound={nm['vmax_bound_kms']:.2f}  gentle={nm['vmax_gentle_kms']:.2f}")
    for c in CLASSES:
        e = res["envelope_vmax_kms"][c]
        print(f"ensemble vmax_{c:7s}: mean {e['mean']:.2f}  p05 {e['p05']:.2f}  p50 {e['p50']:.2f}  "
              f"p95 {e['p95']:.2f}  [min {e['min']:.2f}, max {e['max']:.2f}] km/s")
    print("sensitivity of vmax_bound (variance share):")
    for ax, s in sorted(res["sensitivity_vmax_bound"].items(), key=lambda kv: -kv[1]["variance_share"]):
        print(f"  {ax:6s}  share {s['variance_share']*100:4.0f}%   (std {s['vmax_std_kms']:.3f} km/s)")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
