#!/usr/bin/env python3
"""Ground-truth checks for the high-fidelity (REBOUND+ASSIST+drag) dynamics.

Run:  .venv/bin/python highfidelity/validate_hifi.py
"""
import math
import numpy as np
import hifi
from hifi import AU, DAY, RE, MU_E


def check_kepler_energy():
    """Full self-consistent force set: a bound orbit closes and specific energy
    is ~conserved over one orbit (residual = real third-body work + integrator
    error). NB geocentric mode is non-inertial, so Earth's monopole rides on the
    'PLANETS' term -- you cannot isolate a 'pure two-body' subset."""
    sim, ext = hifi.make_sim()                      # default (full) forces
    rp, ra = RE + 300e3, RE + 3000e3
    a = 0.5 * (rp + ra)
    vp = math.sqrt(MU_E * (2.0 / rp - 1.0 / a))
    sim.add(x=rp / AU, y=0.0, z=0.0, vx=0.0, vy=vp / AU * DAY, vz=0.0)
    r, v = hifi.geo_state_si(sim)
    eps0 = 0.5 * (v @ v) - MU_E / math.sqrt(r @ r)
    period = 2 * math.pi * math.sqrt(a ** 3 / MU_E)
    sim.integrate(sim.t + period / DAY)
    r, v = hifi.geo_state_si(sim)
    eps1 = 0.5 * (v @ v) - MU_E / math.sqrt(r @ r)
    rmag = math.sqrt(r @ r)
    drift = abs(eps1 - eps0) / abs(eps0)
    ok = drift < 1e-6 and abs(rmag - rp) / rp < 1e-4
    print(f"[kepler] energy drift/orbit={drift:.2e}  return-radius err={abs(rmag-rp)/rp:.2e}  -> {'OK' if ok else 'FAIL'}")
    return ok


def check_drag_vs_stdlib():
    """Per-pass dv from the split-drag hi-fi integrator vs the validated stdlib
    integrator (physics.simulate_passage) at low v_inf grazing passes."""
    import physics
    body = physics.Body()
    print(f"[drag] beta = {body.beta:.3e} kg/m^2")
    allok = True
    for hp_km in (60, 50, 40):
        hp = hp_km * 1e3
        std = physics.simulate_passage(0.5e3, hp, body)
        hi = hifi.simulate_passage(0.5e3, hp, body)
        ratio = hi["dv_ms"] / std["dv_ms"] if std["dv_ms"] else float("nan")
        ok = 0.7 < ratio < 1.4 and hi["outcome"] == std["outcome"]
        allok = allok and ok
        print(f"[drag] hp={hp_km:3d}km  outcome std/hi={std['outcome']}/{hi['outcome']}  "
              f"dv_std={std['dv_ms']:6.2f}  dv_hi={hi['dv_ms']:6.2f} m/s  ratio={ratio:4.2f}  -> {'OK' if ok else 'FAIL'}")
    return allok


def _elements(r, v):
    """Return (inclination, RAAN) from geocentric state vectors (SI)."""
    h = np.cross(r, v)
    n = np.cross([0, 0, 1.0], h)          # node vector
    inc = math.acos(max(-1, min(1, h[2] / np.linalg.norm(h))))
    raan = math.atan2(n[1], n[0])
    return inc, raan


def _node_rate(forces, a, inc_deg, n_orbits=20):
    sim, ext = hifi.make_sim(forces=forces)
    vc = math.sqrt(MU_E / a)
    i = math.radians(inc_deg)
    sim.add(x=a / AU, y=0.0, z=0.0,
            vx=0.0, vy=vc * math.cos(i) / AU * DAY, vz=vc * math.sin(i) / AU * DAY)
    r, v = hifi.geo_state_si(sim)
    _, raan0 = _elements(r, v)
    n_mm = math.sqrt(MU_E / a ** 3)
    T = 2 * math.pi / n_mm
    sim.integrate(sim.t + n_orbits * T / DAY)
    r, v = hifi.geo_state_si(sim)
    _, raan1 = _elements(r, v)
    draan = (raan1 - raan0 + math.pi) % (2 * math.pi) - math.pi
    return draan / (n_orbits * T)


def check_j2_nodal_regression():
    """Isolate the EARTH_HARMONICS effect by DIFFERENCING node rates with vs.
    without it (common frame/lunisolar drift cancels). The difference must match
    the analytic J2 rate  dRAAN/dt = -1.5 n J2 (RE/p)^2 cos i."""
    J2 = 1.08262668e-3
    a = RE + 700e3
    inc_deg = 45.0
    rate_on = _node_rate(["PLANETS", "EARTH_HARMONICS"], a, inc_deg)
    rate_off = _node_rate(["PLANETS"], a, inc_deg)
    rate_sim = rate_on - rate_off
    n_mm = math.sqrt(MU_E / a ** 3)
    rate_an = -1.5 * n_mm * J2 * (RE / a) ** 2 * math.cos(math.radians(inc_deg))
    ratio = rate_sim / rate_an
    ok = 0.85 < ratio < 1.15
    print(f"[J2] nodal regression (J2 isolated) sim={math.degrees(rate_sim)*86400:.4f} deg/day  "
          f"analytic={math.degrees(rate_an)*86400:.4f} deg/day  ratio={ratio:.3f}  -> {'OK' if ok else 'FAIL'}")
    return ok


def check_moon_present():
    """Sanity: the Moon's geocentric distance from the ephemeris is ~3.6-4.1e8 m."""
    e = hifi.ephem()
    earth = e.get_particle("earth", 0.0)
    moon = e.get_particle("moon", 0.0)
    d = math.sqrt((moon.x - earth.x) ** 2 + (moon.y - earth.y) ** 2 + (moon.z - earth.z) ** 2) * AU
    ok = 3.5e8 < d < 4.1e8
    print(f"[moon] geocentric distance at J2000 = {d/1e3:,.0f} km  -> {'OK' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    res = [check_kepler_energy(), check_moon_present(),
           check_j2_nodal_regression(), check_drag_vs_stdlib()]
    print("\nALL PASS" if all(res) else "\nSOME CHECKS FAILED")
