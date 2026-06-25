#!/usr/bin/env python3
"""Ground-truth checks for the deorbit flight physics (pure stdlib).

Run: python3 validate.py
Distinguishes "the integrator says X" from "X is verified" per CLAUDE.md.
"""
import math
from deorbit.corridor import physics as P


def check_kepler_closes():
    """Drag off: a bound orbit must conserve energy and return to its start
    radius with the same speed after half a period (apsis to apsis)."""
    body = P.Body()
    body.k = 0.0  # drag off
    # start at perigee of a bound ellipse: rp=RE+300km, ra=RE+3000km
    rp, ra = P.RE + 300e3, P.RE + 3000e3
    a = 0.5 * (rp + ra)
    vp = math.sqrt(P.MU * (2.0 / rp - 1.0 / a))
    s = [rp, 0.0, 0.0, vp]
    eps0 = 0.5 * vp * vp - P.MU / rp
    period = 2 * math.pi * math.sqrt(a ** 3 / P.MU)
    t, dt = 0.0, period / 20000.0
    while t < period:
        s = P._rk4(s, dt, 0.0)
        t += dt
    r = math.hypot(s[0], s[1])
    v = math.hypot(s[2], s[3])
    eps1 = 0.5 * v * v - P.MU / r
    drift = abs(eps1 - eps0) / abs(eps0)
    ok = drift < 1e-4 and abs(r - rp) / rp < 2e-3
    print(f"[kepler] energy drift/orbit={drift:.2e}  return radius err={abs(r-rp)/rp:.2e}  -> {'OK' if ok else 'FAIL'}")
    return ok


def check_dv_vs_analytic():
    """Per-pass dv from the integrator vs analytic grazing-column estimate
    dv ~ v * Sigma/beta, Sigma = rho(h_p)*sqrt(2*pi*RE*H), at low v_inf where
    the body skips (so the passage is a clean single graze)."""
    body = P.Body()                         # 1 km iron, Cd=1
    H = 7000.0
    print(f"[dv] beta = {body.beta:.3e} kg/m^2")
    allok = True
    for hp_km in (60, 50, 40):
        hp = hp_km * 1e3
        res = P.simulate_passage(0.5e3, hp, body)
        v = math.sqrt(0.5e3 ** 2 + 2 * P.MU / (P.RE + hp))
        sigma = P.density(hp) * math.sqrt(2 * math.pi * P.RE * H)
        dv_an = v * sigma / body.beta
        ratio = res["dv_ms"] / dv_an if dv_an else float("nan")
        ok = 0.3 < ratio < 3.0              # order-of-magnitude agreement
        allok = allok and ok
        print(f"[dv] hp={hp_km:3d}km  outcome={res['outcome']:7s}  "
              f"dv_sim={res['dv_ms']:7.2f} m/s  dv_analytic={dv_an:7.2f} m/s  "
              f"ratio={ratio:4.2f}  -> {'OK' if ok else 'FAIL'}")
    return allok


def check_regimes():
    """Physical sanity: lowering periapsis at low v_inf must show a skip->capture
    transition (the corridor edge); a vacuum periapsis at/below the surface must
    give a direct impact."""
    body = P.Body()
    print("[regime] v_inf=0.5 km/s, lowering periapsis:")
    seen = set()
    for hp_km in (110, 90, 70, 50, 40, 30, 20, 10, 2):
        res = P.simulate_passage(0.5e3, hp_km * 1e3, body)
        seen.add(res["outcome"])
        extra = ""
        if res["outcome"] == "capture":
            extra = f" apo={res.get('apoapsis_km',0):.0f}km"
        print(f"[regime] hp={hp_km:3d}km -> {res['outcome']:7s} "
              f"min_alt={res['min_alt_km']:5.1f}km q={res['peak_q_MPa']:.2f}MPa{extra}")
    # vacuum periapsis right at the surface -> must hit
    imp = P.simulate_passage(0.5e3, 0.0, body)
    print(f"[regime] hp=  0km -> {imp['outcome']:7s} (direct-impact probe)")
    transition = {"skip", "capture"}.issubset(seen)
    ok = transition and imp["outcome"] == "impact"
    print(f"[regime] outcomes seen: {sorted(seen)}; skip->capture transition="
          f"{transition}; surface-grazer impacts={imp['outcome']=='impact'} -> {'OK' if ok else 'FAIL'}")
    return ok


# --- pytest entry points (the check_* helpers above are the actual physics) ---
def test_kepler_closes():
    assert check_kepler_closes()


def test_dv_vs_analytic():
    assert check_dv_vs_analytic()


def test_regimes():
    assert check_regimes()


if __name__ == "__main__":
    import sys
    results = [check_kepler_closes(), check_dv_vs_analytic(), check_regimes()]
    print("\nALL PASS" if all(results) else "\nSOME CHECKS FAILED")
    sys.exit(0 if all(results) else 1)
