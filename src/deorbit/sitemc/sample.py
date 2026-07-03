"""Arrival sampling + 3D hyperbolic entry-state construction (pure stdlib).

Priors are DELIBERATELY simple and fully recorded per draw so the catalog can
be reweighted post-hoc (docs/site_mc_plan.md §4): v_inf ~ U[0.2, 5] km/s,
arrival asymptote isotropic, B-plane clock angle uniform, vacuum-perigee
altitude ~ U[0, 60 km] (record b so per-encounter flux weights d(b^2) can be
applied later), Cd/f_rho/f_H from the corridor-ensemble priors, epoch uniform
over one 18.6-yr lunar-node cycle.

Geometry: incoming asymptote velocity direction u_hat; B-plane unit B_hat
(perp. to u_hat); angular momentum h_hat = B_hat x u_hat; with e, the
perifocal basis is
    p_hat = (u_hat + sqrt(e^2-1) B_hat)/e
    q_hat = (sqrt(e^2-1) u_hat - B_hat)/e
(from u_hat = (p_hat + sqrt(e^2-1) q_hat)/e as nu -> -nu_inf). The entry state
is evaluated at radius r0 inbound (nu0 = -acos((p/r0-1)/e)). Invariants gated
in tests/test_sitemc_sample.py.
"""
import math
import random

from deorbit.corridor.physics import MU, RE

ENTRY_R0 = RE + 200e3     # spherical entry interface radius (matches 2D)


def entry_state(v_inf, u_hat, theta_B, rp_vac, r0=ENTRY_R0):
    """Inbound state [x,y,z,vx,vy,vz] at radius r0 for the hyperbola defined
    by (v_inf [m/s], asymptote direction u_hat, B-plane angle theta_B,
    vacuum perigee radius rp_vac [m]). Also returns derived (b, i_deg,
    node_deg, h_hat)."""
    ux, uy, uz = u_hat
    # B-plane basis perpendicular to u_hat
    if abs(uz) < 0.999:
        e1 = _norm(_cross(u_hat, (0.0, 0.0, 1.0)))
    else:
        e1 = _norm(_cross(u_hat, (1.0, 0.0, 0.0)))
    e2 = _cross(u_hat, e1)
    ct, st = math.cos(theta_B), math.sin(theta_B)
    B_hat = (ct * e1[0] + st * e2[0], ct * e1[1] + st * e2[1], ct * e1[2] + st * e2[2])

    e = 1.0 + rp_vac * v_inf * v_inf / MU
    b = rp_vac * math.sqrt(1.0 + 2.0 * MU / (rp_vac * v_inf * v_inf))
    h = b * v_inf
    p = h * h / MU
    se = math.sqrt(e * e - 1.0)
    p_hat = tuple((u_hat[i] + se * B_hat[i]) / e for i in range(3))
    q_hat = tuple((se * u_hat[i] - B_hat[i]) / e for i in range(3))

    cnu = max(-1.0, min(1.0, (p / r0 - 1.0) / e))
    nu0 = -math.acos(cnu)                     # inbound branch
    cn, sn = math.cos(nu0), math.sin(nu0)
    pos = tuple(r0 * (cn * p_hat[i] + sn * q_hat[i]) for i in range(3))
    vf = MU / h
    vel = tuple(vf * (-sn * p_hat[i] + (e + cn) * q_hat[i]) for i in range(3))

    h_vec = _cross(pos, vel)
    hn = math.sqrt(sum(c * c for c in h_vec))
    i_deg = math.degrees(math.acos(max(-1.0, min(1.0, h_vec[2] / hn))))
    node_deg = math.degrees(math.atan2(h_vec[0], -h_vec[1])) % 360.0
    return list(pos + vel), {"b_m": b, "e": e, "i_deg": i_deg, "node_deg": node_deg,
                             "h_hat": tuple(c / hn for c in h_vec)}


def draw(seed):
    """One fully-recorded prior draw."""
    rng = random.Random(seed)
    v_inf = 200.0 + 4800.0 * rng.random()
    zc = 2.0 * rng.random() - 1.0             # isotropic asymptote direction
    az = 2.0 * math.pi * rng.random()
    sz = math.sqrt(max(0.0, 1.0 - zc * zc))
    u_hat = (sz * math.cos(az), sz * math.sin(az), zc)
    theta_B = 2.0 * math.pi * rng.random()
    rp_alt = 60e3 * rng.random()
    cd = min(2.2, max(0.7, 1.1 * math.exp(rng.gauss(0.0, 0.25))))
    f_rho = math.exp(rng.gauss(0.0, 0.2))
    f_H = min(1.15, max(0.85, rng.gauss(1.0, 0.05)))
    epoch_jd = 2451545.0 + rng.random() * 6798.0     # one lunar-node cycle
    return {"seed": seed, "v_inf": v_inf, "u_hat": u_hat, "theta_B": theta_B,
            "rp_alt_m": rp_alt, "cd": cd, "f_rho": f_rho, "f_H": f_H,
            "epoch_jd": epoch_jd, "asym_dec_deg": math.degrees(math.asin(zc))}


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)
