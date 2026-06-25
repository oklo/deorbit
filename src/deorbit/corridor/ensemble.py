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

# grid: captures live at v_inf<~4.2 km/s and periapsis<~45 km (from the deterministic corridor)
V_LO, V_HI = 0.2e3, 5.0e3
P_LO, P_HI = 0.0, 60e3


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


def corridor_metrics(vs, ps, cd=1.0, f_rho=1.0, f_H=1.0, omega=0.0):
    """Classify the grid for one realization. Returns (capture_mask, n_capture,
    vmax_capture_kms) where capture_mask[i][j] is True for a captured cell."""
    body = physics.Body(cd=cd)
    mask = [[False] * len(ps) for _ in range(len(vs))]
    n_cap = 0
    vmax = 0.0
    for i, v in enumerate(vs):
        captured_here = False
        for j, p in enumerate(ps):
            r = physics.simulate_passage(v, p, body, f_rho=f_rho, f_H=f_H, omega=omega)
            if r["outcome"] == "capture":
                mask[i][j] = True
                n_cap += 1
                captured_here = True
        if captured_here and v > vmax:
            vmax = v
    return mask, n_cap, vmax / 1e3


def run_ensemble(n=200, nv=28, npa=30, seed=0):
    vs, ps = _grid(nv, npa)
    rng = random.Random(seed)

    # nominal (deterministic, non-rotating, Cd=1.0) -- validation anchor
    nom_mask, nom_n, nom_vmax = corridor_metrics(vs, ps)

    pcap = [[0] * npa for _ in range(nv)]
    n_caps, vmaxes = [], []
    for _ in range(n):
        pr = sample_params(rng)
        mask, ncap, vmax = corridor_metrics(vs, ps, pr["cd"], pr["f_rho"], pr["f_H"], pr["omega"])
        for i in range(nv):
            for j in range(npa):
                if mask[i][j]:
                    pcap[i][j] += 1
        n_caps.append(ncap)
        vmaxes.append(vmax)
    pcap = [[c / n for c in row] for row in pcap]

    # single-axis variance share: vary ONE axis (others nominal), measure std(vmax)
    def single_axis(axis, m=40):
        r2 = random.Random(seed + 1)
        out = []
        for _ in range(m):
            pr = {"cd": 1.0, "f_rho": 1.0, "f_H": 1.0, "omega": 0.0}
            s = sample_params(r2)
            pr[axis] = s[axis]
            out.append(corridor_metrics(vs, ps, pr["cd"], pr["f_rho"], pr["f_H"], pr["omega"])[2])
        return out

    sens = {}
    for axis in ("cd", "f_rho", "omega", "f_H"):
        vals = single_axis(axis)
        mu = sum(vals) / len(vals)
        sd = (sum((x - mu) ** 2 for x in vals) / len(vals)) ** 0.5
        sens[axis] = {"vmax_mean_kms": mu, "vmax_std_kms": sd}
    tot_var = sum(s["vmax_std_kms"] ** 2 for s in sens.values()) or 1.0
    for s in sens.values():
        s["variance_share"] = s["vmax_std_kms"] ** 2 / tot_var

    def pct(xs, q):
        xs = sorted(xs)
        return xs[min(len(xs) - 1, int(q * len(xs)))]

    return {
        "priors": {"dens_sigma_ln": DENS_SIGMA_LN, "h_sigma": H_SIGMA,
                   "cd_median": CD_MEDIAN, "cd_sigma_ln": CD_SIGMA_LN, "cd_clip": CD_CLIP},
        "grid": {"nv": nv, "npa": npa, "v_kms": [V_LO / 1e3, V_HI / 1e3], "peri_km": [P_LO / 1e3, P_HI / 1e3]},
        "n_samples": n,
        "nominal": {"n_capture_cells": nom_n, "vmax_capture_kms": nom_vmax},
        "envelope": {
            "vmax_capture_kms": {"mean": sum(vmaxes) / n, "p05": pct(vmaxes, 0.05),
                                 "p50": pct(vmaxes, 0.5), "p95": pct(vmaxes, 0.95),
                                 "min": min(vmaxes), "max": max(vmaxes)},
            "n_capture_cells": {"mean": sum(n_caps) / n, "p05": pct(n_caps, 0.05),
                                "p50": pct(n_caps, 0.5), "p95": pct(n_caps, 0.95)},
        },
        "sensitivity_vmax": sens,
        "vs_kms": [v / 1e3 for v in vs],
        "ps_km": [p / 1e3 for p in ps],
        "p_capture": pcap,
    }


def ascii_map(res):
    """Terminal heatmap of P(capture): rows = periapsis (top=high), cols = v_inf."""
    pcap, vs, ps = res["p_capture"], res["vs_kms"], res["ps_km"]
    nv, npa = len(vs), len(ps)
    chars = " .:+#"
    lines = [f"P(capture) corridor  (N={res['n_samples']}; rows=periapsis km, cols=v_inf km/s)"]
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
    ap.add_argument("--nv", type=int, default=28)
    ap.add_argument("--npa", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
        "results", "corridor_ensemble.json"))
    a = ap.parse_args()
    res = run_ensemble(a.n, a.nv, a.npa, a.seed)

    print(ascii_map(res))
    e = res["envelope"]["vmax_capture_kms"]
    print(f"\nnominal (Cd=1, no rotation): vmax_capture = {res['nominal']['vmax_capture_kms']:.2f} km/s, "
          f"{res['nominal']['n_capture_cells']} capture cells")
    print(f"ensemble vmax_capture: mean {e['mean']:.2f}  p05 {e['p05']:.2f}  p50 {e['p50']:.2f}  "
          f"p95 {e['p95']:.2f}  [min {e['min']:.2f}, max {e['max']:.2f}] km/s")
    print("sensitivity of vmax_capture (variance share):")
    for ax, s in sorted(res["sensitivity_vmax"].items(), key=lambda kv: -kv[1]["variance_share"]):
        print(f"  {ax:6s}  share {s['variance_share']*100:4.0f}%   (std {s['vmax_std_kms']:.3f} km/s)")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
