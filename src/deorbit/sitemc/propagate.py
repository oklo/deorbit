"""3D geocentric propagator for the site-selection MC (pure stdlib).

Forces: Earth point-mass + J2 + oblate-atmosphere drag (co-rotating wind) +
lunisolar point-mass tides (analytic ephemerides). Every force has a switch
so limit gates can isolate it (J2 off + Moon/Sun off + oblate=False must
reproduce the validated 2D corridor integrator).

State s = [x, y, z, vx, vy, vz] (m, m/s), t = seconds past cfg.jd0.
"""
import math

from deorbit.corridor.physics import RE as RE_SPHERE, density as _rho_profile
from deorbit.sitemc.geodesy import A_EQ, ellipsoid_radius
from deorbit.sitemc import ephem
from deorbit.sitemc.geodesy import density_oblate

MU_E = 3.986004418e14
J2 = 1.08262668e-3
OMEGA_E = 7.2921159e-5
ENTRY_ALT = 200e3         # atmosphere interface (geodetic), same as 2D

# highest terrain ~6.4 km geocentric-relevant; check terrain below this geodetic alt
TERRAIN_ALT = 9.5e3


class Config:
    """Force/model switches + body + epoch. Defaults = full physics."""

    def __init__(self, body_k, jd0=2451545.0, oblate=True, j2=True, moonsun=True,
                 omega=OMEGA_E, f_rho=1.0, f_H=1.0, terrain=None, re_sphere=RE_SPHERE,
                 sph_surface=False):
        self.k = body_k                  # 0.5*Cd*A/m
        self.jd0 = jd0
        self.oblate = oblate
        # sph_surface: impact/altitude REFERENCE surface is the re_sphere ball
        # (full-legacy 2D parity). Default False = real ellipsoid surface even
        # when the DRAG is spherical (oblate=False), so atmosphere-model A/Bs
        # don't smuggle in a phantom impact surface.
        self.sph_surface = sph_surface
        self.j2 = j2
        self.moonsun = moonsun
        self.omega = omega
        self.f_rho = f_rho
        self.f_H = f_H
        self.terrain = terrain           # None or object with rock_radius(latc_rad, lon_e_rad)
        self.re_sphere = re_sphere       # spherical-mode reference radius
        self.cache = ephem.EphemCache(jd0) if moonsun else None


def _accel(t, s, cfg):
    x, y, z, vx, vy, vz = s
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    g = -MU_E / (r2 * r)
    ax, ay, az = g * x, g * y, g * z
    if cfg.j2:
        zr2 = (z * z) / r2
        f = -1.5 * J2 * (MU_E / r2) * (A_EQ * A_EQ / r2) / r
        ax += f * x * (1.0 - 5.0 * zr2)
        ay += f * y * (1.0 - 5.0 * zr2)
        az += f * z * (3.0 - 5.0 * zr2)
    if cfg.moonsun:
        (mx, my, mz), (sx, sy, sz) = cfg.cache.moon_sun(t)
        for (px, py, pz), mu in (((mx, my, mz), ephem.MU_MOON), ((sx, sy, sz), ephem.MU_SUN)):
            dx, dy, dz = px - x, py - y, pz - z
            d3 = (dx * dx + dy * dy + dz * dz) ** 1.5
            p3 = (px * px + py * py + pz * pz) ** 1.5
            ax += mu * (dx / d3 - px / p3)
            ay += mu * (dy / d3 - py / p3)
            az += mu * (dz / d3 - pz / p3)
    # drag (only possibly nonzero below the interface; cheap reject first)
    if r < A_EQ + ENTRY_ALT + 25e3:
        if cfg.oblate:
            rho = density_oblate(x, y, z, cfg.f_rho, cfg.f_H, True)
        else:
            # spherical mode measures altitude from cfg.re_sphere (2D-limit parity)
            rho = _sph_density(r, cfg)
        if rho > 0.0:
            w = cfg.omega
            rvx, rvy, rvz = vx + w * y, vy - w * x, vz
            vrel = math.sqrt(rvx * rvx + rvy * rvy + rvz * rvz)
            if vrel > 0.0:
                fd = cfg.k * rho * vrel
                ax -= fd * rvx
                ay -= fd * rvy
                az -= fd * rvz
    return ax, ay, az


def _sph_density(r, cfg):
    return _rho_profile(r - cfg.re_sphere, cfg.f_rho, cfg.f_H)


def _rk4(t, s, dt, cfg):
    a1 = _accel(t, s, cfg)
    k1 = (s[3], s[4], s[5]) + a1
    s2 = [s[i] + 0.5 * dt * k1[i] for i in range(6)]
    a2 = _accel(t + 0.5 * dt, s2, cfg)
    k2 = (s2[3], s2[4], s2[5]) + a2
    s3 = [s[i] + 0.5 * dt * k2[i] for i in range(6)]
    a3 = _accel(t + 0.5 * dt, s3, cfg)
    k3 = (s3[3], s3[4], s3[5]) + a3
    s4 = [s[i] + dt * k3[i] for i in range(6)]
    a4 = _accel(t + dt, s4, cfg)
    k4 = (s4[3], s4[4], s4[5]) + a4
    return [s[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6)]


def _alt_ref(x, y, z, r, cfg):
    """Reference altitude used for step control / impact test."""
    if cfg.sph_surface:
        return r - cfg.re_sphere
    sl = z / r
    cl = math.sqrt(max(0.0, 1.0 - sl * sl))
    return r - ellipsoid_radius(sl, cl)


def fly(s0, cfg, t0=0.0, max_days=200.0, want_track=False):
    """Integrate from an inbound entry state to an outcome.

    Returns dict: outcome 'impact'|'escape'|'timeout', diagnostics, and the
    impact state if any. Pass bookkeeping: a 'pass' = one excursion below the
    atmosphere interface (entry->exit), so radial-velocity wobbles during a
    single skimming arc do not inflate the count; apogee/eps are recorded at
    atmosphere EXIT (after the full drag dv of the pass). A low pass reached
    geodetic alt < TERRAIN_ALT.
    """
    s = list(s0)
    t = t0
    tmax = t0 + max_days * 86400.0
    n_pass = 0
    n_low_pass = 0
    apo1_km = None
    eps_after_1 = None
    in_atm = False
    pass_min_alt = 1e30
    min_alt_ever = 1e30
    steps = 0
    outcome = "timeout"
    track = [] if want_track else None
    while t < tmax:
        x, y, z, vx, vy, vz = s
        r = math.sqrt(x * x + y * y + z * z)
        v2 = vx * vx + vy * vy + vz * vz
        alt = _alt_ref(x, y, z, r, cfg)
        if alt < pass_min_alt:
            pass_min_alt = alt
        if alt < min_alt_ever:
            min_alt_ever = alt
        # --- impact test ---
        if alt <= 0.0:
            outcome = "impact"
            break
        if cfg.terrain is not None and alt < TERRAIN_ALT:
            sl = z / r
            latc = math.asin(max(-1.0, min(1.0, sl)))
            lon_e = math.atan2(y, x) - ephem.gmst_rad(cfg.jd0 + t / 86400.0)
            if r <= cfg.terrain.rock_radius(latc, lon_e):
                outcome = "impact"
                break
        vr = (x * vx + y * vy + z * vz) / r
        eps = 0.5 * v2 - MU_E / r
        # --- pass bookkeeping (one excursion below the interface = one pass) ---
        if not in_atm and alt < ENTRY_ALT:
            in_atm = True
            pass_min_alt = alt
        elif in_atm and alt >= ENTRY_ALT and vr > 0.0:
            in_atm = False
            n_pass += 1
            if pass_min_alt < TERRAIN_ALT:
                n_low_pass += 1
            if n_pass == 1 and eps < 0.0:
                a = -MU_E / (2.0 * eps)
                hx = y * vz - z * vy
                hy = z * vx - x * vz
                hz = x * vy - y * vx
                h2 = hx * hx + hy * hy + hz * hz
                e2 = max(0.0, 1.0 + 2.0 * eps * h2 / (MU_E * MU_E))
                apo1_km = (a * (1.0 + math.sqrt(e2)) - RE_SPHERE) / 1e3
                eps_after_1 = eps
            pass_min_alt = 1e30
        # --- escape tests ---
        if eps > 0.0 and vr > 0.0 and r > 1e8:
            outcome = "escape"
            break
        if r > 1.5e9:
            outcome = "escape"
            break
        # --- step control (2D corridor rule in-atmosphere; orbital fraction on coast) ---
        if alt < 150e3:
            dt = max(0.02, min(0.5, 200.0 / (abs(vr) + 10.0)))
        else:
            dt = max(1.0, min(600.0, 0.02 * r / math.sqrt(v2)))
        if want_track and alt < TERRAIN_ALT + 5e3:
            track.append((t, x, y, z, alt))
        s = _rk4(t, s, dt, cfg)
        t += dt
        steps += 1

    x, y, z, vx, vy, vz = s
    r = math.sqrt(x * x + y * y + z * z)
    v = math.sqrt(vx * vx + vy * vy + vz * vz)
    vr = (x * vx + y * vy + z * vz) / r
    sl = max(-1.0, min(1.0, z / r))
    res = {
        "outcome": outcome, "t_days": (t - t0) / 86400.0, "steps": steps,
        "n_pass": n_pass, "n_low_pass": n_low_pass, "apo1_km": apo1_km,
        "eps_after_1": eps_after_1, "min_alt_km": min_alt_ever / 1e3,
        "latc_deg": math.degrees(math.asin(sl)),
        "lon_e_deg": math.degrees((math.atan2(y, x)
                                   - ephem.gmst_rad(cfg.jd0 + t / 86400.0)) % (2 * math.pi)),
        "v_end": v, "fpa_deg": math.degrees(math.asin(max(-1.0, min(1.0, vr / v)))) if v > 0 else 0.0,
        "state": s, "t_end": t,
    }
    if want_track:
        res["track"] = track
    return res
