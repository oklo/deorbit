"""High-fidelity capture dynamics: REBOUND + ASSIST (JPL DE440 ephemeris) with
Earth J2/J3/J4 + Sun/Moon/planets + GR, and atmospheric drag added by
operator-splitting.

Why this exists: the stdlib `physics.py`/`deorbit.py` corridor sweep is a fast
2-body+drag classifier. It cannot answer whether a "capture" survives the Moon
(first apogees are typically ~2.4e5 km, deep in the lunar zone) or how J2 steers
the long aerobraking evolution. This module does the dynamics correctly:

  * Gravity / solar-system dynamics: ASSIST integrates a massless test particle
    in the field of the DE440 bodies (Sun, planets, Moon, 16 asteroids) with
    Earth harmonics (J2,J3,J4), Sun J2, and GR. Run in ASSIST's GEOCENTRIC mode
    so the particle's coordinates are Earth-centered (altitude is immediate).
  * Drag: ASSIST owns rebound's force hook, so drag is applied by Lie/Strang
    operator-splitting -- integrate gravity for a small dt, then apply the drag
    velocity kick computed from the geocentric altitude and the ATMOSPHERE-
    RELATIVE velocity (Earth co-rotation subtracted).

Units: ASSIST uses AU, AU/day, ICRF equatorial. Drag physics is done in SI and
converted. Atmosphere model and Body (beta) are reused from the validated
stdlib `physics.py` (single source of truth).

See validate_hifi.py for ground-truth checks (Kepler energy closure, drag dv
vs the stdlib model, J2 nodal-regression vs the analytic rate).
"""
import os
import sys
import math

import numpy as np
import rebound
import assist

# reuse the validated stdlib atmosphere + body
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import physics  # noqa: E402

AU = 1.495978707e11        # m
DAY = 86400.0              # s
RE = physics.RE            # m
MU_E = physics.MU          # m^3/s^2
ENTRY_ALT = physics.ENTRY_ALT
OMEGA_E = 7.2921159e-5     # Earth rotation rate, rad/s (axis ~ ICRF z in geo frame)

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
PLANET_EPH = os.path.join(_DATA, "linux_p1550p2650.440")
SMALLBODY_EPH = os.path.join(_DATA, "sb441-n16.bsp")

_EPHEM = None


def ephem():
    global _EPHEM
    if _EPHEM is None:
        _EPHEM = assist.Ephem(PLANET_EPH, SMALLBODY_EPH)
    return _EPHEM


def make_sim(epoch_days=0.0, forces=None):
    """Geocentric ASSIST simulation. epoch_days is days past ephem.jd_ref
    (J2000 = 2451545.0). If `forces` is given (a list) it overrides the default
    force set (use [] for pure point-mass Earth -> exact Kepler)."""
    sim = rebound.Simulation()
    sim.t = epoch_days
    ext = assist.Extras(sim, ephem())
    ext.geocentric = 1
    if forces is not None:
        ext.forces = forces
    return sim, ext


def _drag_kick_si(r_m, v_m, body):
    """Drag acceleration vector (m/s^2) at geocentric position r_m, inertial
    velocity v_m (both numpy length-3, SI). Uses atmosphere-relative velocity."""
    alt = math.sqrt(r_m @ r_m) - RE
    rho = physics.density(alt)
    if rho <= 0.0:
        return np.zeros(3), alt, 0.0
    omega = np.array([0.0, 0.0, OMEGA_E])
    v_rel = v_m - np.cross(omega, r_m)          # atmosphere co-rotates with Earth
    vmag = math.sqrt(v_rel @ v_rel)
    if vmag == 0.0:
        return np.zeros(3), alt, 0.0
    a = -0.5 * rho * vmag * body.cd * body.area / body.mass * v_rel
    q = 0.5 * rho * vmag * vmag
    return a, alt, q


def geo_state_si(sim):
    p = sim.particles[0]
    r = np.array([p.x, p.y, p.z]) * AU
    v = np.array([p.vx, p.vy, p.vz]) * AU / DAY
    return r, v


def set_geo_state_si(sim, r_m, v_m):
    p = sim.particles[0]
    p.x, p.y, p.z = (r_m / AU).tolist()
    p.vx, p.vy, p.vz = (v_m / AU * DAY).tolist()


def initial_state_si(v_inf, peri_alt_m, inclination_deg=0.0, r0=RE + ENTRY_ALT):
    """Inbound geocentric state (SI) at interface r0 for a vacuum orbit with
    hyperbolic-excess speed v_inf and periapsis RE+peri_alt_m, inclined by
    inclination_deg to the equator (node on +x)."""
    rp = RE + peri_alt_m
    vp = math.sqrt(v_inf * v_inf + 2.0 * MU_E / rp)
    h = rp * vp
    v0 = math.sqrt(v_inf * v_inf + 2.0 * MU_E / r0)
    c = max(-1.0, min(1.0, h / (r0 * v0)))
    fpa = math.acos(c)
    vr = -v0 * math.sin(fpa)
    vt = v0 * math.cos(fpa)
    i = math.radians(inclination_deg)
    r_m = np.array([r0, 0.0, 0.0])
    v_m = np.array([vr, vt * math.cos(i), vt * math.sin(i)])
    return r_m, v_m


def simulate_passage(v_inf, peri_alt_m, body, inclination_deg=0.0,
                     epoch_days=0.0, max_steps=2_000_000):
    """Integrate one atmospheric passage with full forces + drag (split).
    Returns outcome dict comparable to physics.simulate_passage."""
    sim, ext = make_sim(epoch_days)
    r_m, v_m = initial_state_si(v_inf, peri_alt_m, inclination_deg)
    sim.add(x=float(r_m[0] / AU), y=float(r_m[1] / AU), z=float(r_m[2] / AU),
            vx=float(v_m[0] / AU * DAY), vy=float(v_m[1] / AU * DAY),
            vz=float(v_m[2] / AU * DAY))
    r0 = RE + ENTRY_ALT
    v_in = math.sqrt(v_m @ v_m)
    min_alt = ENTRY_ALT
    peak_q = 0.0
    passed_peri = False
    outcome = "timeout"
    steps = 0
    while steps < max_steps:
        r, v = geo_state_si(sim)
        rmag = math.sqrt(r @ r)
        alt = rmag - RE
        if alt < min_alt:
            min_alt = alt
        vr = float(r @ v) / rmag
        if vr > 0.0:
            passed_peri = True
        if alt <= 0.0:
            outcome = "impact"
            break
        if passed_peri and rmag >= r0 and vr > 0.0:
            eps = 0.5 * (v @ v) - MU_E / rmag
            outcome = "capture" if eps < 0.0 else "skip"
            break
        dt_s = max(0.02, min(0.5, 200.0 / (abs(vr) + 10.0))) if alt < ENTRY_ALT else 20.0
        sim.integrate(sim.t + dt_s / DAY)
        if alt < ENTRY_ALT:                       # operator-split drag kick
            r2, v2 = geo_state_si(sim)
            a_drag, _, q = _drag_kick_si(r2, v2, body)
            if q > peak_q:
                peak_q = q
            set_geo_state_si(sim, r2, v2 + a_drag * dt_s)
        steps += 1

    r, v = geo_state_si(sim)
    rmag = math.sqrt(r @ r)
    vmag = math.sqrt(v @ v)
    eps = 0.5 * vmag * vmag - MU_E / rmag
    out = {"outcome": outcome, "v_inf": v_inf, "peri_alt_km": peri_alt_m / 1e3,
           "inclination_deg": inclination_deg, "min_alt_km": min_alt / 1e3,
           "peak_q_MPa": peak_q / 1e6, "dv_ms": v_in - vmag, "eps_out": eps,
           "steps": steps}
    if outcome == "capture":
        h = float(np.linalg.norm(np.cross(r, v)))
        a = -MU_E / (2.0 * eps)
        e = math.sqrt(max(0.0, 1.0 + 2.0 * eps * h * h / (MU_E * MU_E)))
        out["apoapsis_km"] = (a * (1.0 + e) - RE) / 1e3
        out["new_peri_km"] = (a * (1.0 - e) - RE) / 1e3
    return out
