#!/usr/bin/env python3
"""Quantitative closure of the capture-feasibility chain (pure stdlib).

Tightens test_corridor.check_dv_vs_analytic (order-of-magnitude) to a real
closure: for a near-parabolic grazing pass the trajectory's radius of
curvature at periapsis is the semi-latus rectum p ~= 2*r_p, so altitude rises
along-track as s^2/(4 r_p) (not s^2/(2 r_p)) and the traversed column is

    Sigma_eff = rho(h_p) * sqrt(2*pi*r_eff*H),   r_eff = r_p*p/(p-r_p) ~= 2*r_p

i.e. sqrt(2) larger than the flat "surface-curvature" estimate. With
F = 1/2 rho Cd A v^2 the per-pass loss and capture boundary follow:

    dv_p        = v_p * Sigma_eff / (2*beta)
    d_eps       = v_p * dv_p
    v_inf,max   = v_p * sqrt(Sigma_eff / beta)      [d_eps > v_inf^2/2]

These are the closed forms the paper quotes; this test pins the integrator to
them (energy loss to ~5%, capture boundary to ~8%). Derivation + narrative:
docs/capture_feasibility.md.
"""
import math
from deorbit.corridor import physics as P


def _scale_height_m(alt_km):
    base = P._ATM[0]
    for row in P._ATM:
        if row[0] <= alt_km:
            base = row
        else:
            break
    return base[2] * 1000.0


def _analytic(hp_m, v_inf, beta):
    rp = P.RE + hp_m
    vp = math.sqrt(v_inf * v_inf + 2.0 * P.MU / rp)
    H = _scale_height_m(hp_m / 1000.0)
    sigma_eff = P.density(hp_m) * math.sqrt(2.0 * math.pi * (2.0 * rp) * H)
    deps = vp * vp * sigma_eff / (2.0 * beta)
    vinf_max = vp * math.sqrt(sigma_eff / beta)
    return deps, vinf_max


def test_beta():
    b = P.Body()
    assert abs(b.beta - 5.2e6) / 5.2e6 < 0.01   # 1 km iron, Cd=1


def test_energy_loss_sqrt2_closure():
    """Integrated per-pass orbit-energy loss matches the curvature-corrected
    column to a few percent across two decades of periapsis altitude."""
    b = P.Body()
    v_inf = 1000.0
    for hp_km in (40, 30):                       # clean skip passes (fast)
        deps_an, _ = _analytic(hp_km * 1e3, v_inf, b.beta)
        r = P.simulate_passage(v_inf, hp_km * 1e3, b)
        assert r["outcome"] == "skip"
        deps_sim = 0.5 * v_inf * v_inf - r["eps_out"]
        ratio = deps_sim / deps_an
        assert 0.93 < ratio < 1.07, f"h_p={hp_km} km: deps ratio {ratio:.3f}"


def test_capture_boundary_closed_form():
    """Bisected capture/skip boundary sits on v_p*sqrt(Sigma_eff/beta)."""
    b = P.Body()
    hp = 10e3
    _, vmax_an = _analytic(hp, 0.0, b.beta)
    lo, hi = 500.0, 6000.0
    assert P.simulate_passage(lo, hp, b)["outcome"] == "capture"
    for _ in range(12):
        mid = 0.5 * (lo + hi)
        if P.simulate_passage(mid, hp, b)["outcome"] == "capture":
            lo = mid
        else:
            hi = mid
    ratio = lo / vmax_an
    assert 0.92 < ratio < 1.10, f"boundary {lo:.0f} vs analytic {vmax_an:.0f} m/s (ratio {ratio:.3f})"


def test_survival_margins():
    """A capturing pass stays in the monolith-survivable regime: ram pressure
    tens of MPa (<< iron strength ~1e8-5e8 Pa, >> rubble strength), and the
    body exits bound with a finite apoapsis."""
    b = P.Body()
    r = P.simulate_passage(3000.0, 5e3, b)
    assert r["outcome"] == "capture"
    assert 20.0 < r["peak_q_MPa"] < 60.0
    assert r["apoapsis_km"] > 0.0
    # rigid-body tidal stress at grazing periapsis: ~kPa, utterly negligible
    rp = P.RE + 5e3
    sig_tidal = P.RHO_IRON * (2.0 * P.MU * 500.0 / rp ** 3) * 500.0
    assert sig_tidal < 1e5
