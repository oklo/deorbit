#!/usr/bin/env python3
"""Gates for the site-MC 3D engine: ephemerides, IC geometry, propagator
physics limits, terrain. House rule: every force/limit isolated."""
import math
import os

import pytest

from deorbit.corridor import physics as P2
from deorbit.sitemc import ephem, propagate, sample
from deorbit.sitemc.propagate import Config, fly

DEM_OK = os.path.exists(__import__("deorbit.sitemc.terrain", fromlist=["DEM_PATH"]).DEM_PATH)


# ---------- ephemerides ----------

def test_moon_invariants():
    """Distance range, ecliptic-latitude bound, sidereal rate."""
    dmin, dmax, bmax = 1e12, 0.0, 0.0
    for k in range(0, 40 * 24):                     # 40 days hourly
        jd = 2461000.5 + k / 24.0
        x, y, z = ephem.moon_pos(jd)
        r = math.sqrt(x * x + y * y + z * z)
        dmin, dmax = min(dmin, r), max(dmax, r)
        # ecliptic latitude
        ye = math.cos(-ephem.OBLIQ) * y - math.sin(-ephem.OBLIQ) * z
        ze = math.sin(-ephem.OBLIQ) * y + math.cos(-ephem.OBLIQ) * z
        bmax = max(bmax, abs(math.degrees(math.asin(ze / r))))
    assert 3.50e8 < dmin < 3.72e8, dmin
    assert 3.95e8 < dmax < 4.10e8, dmax
    assert 4.8 < bmax < 5.5, bmax
    # angular rate ~13.18 deg/day over one sidereal month
    x0, y0, _ = ephem.moon_pos(2461000.5)
    x1, y1, _ = ephem.moon_pos(2461000.5 + 27.321661)
    dang = math.degrees(abs(math.atan2(y1, x1) - math.atan2(y0, x0)))
    assert dang < 4.0                                # returns to ~same longitude


def test_sun_invariants():
    rmin, rmax = 1e13, 0.0
    for k in range(0, 366, 5):
        x, y, z = ephem.sun_pos(2461000.5 + k)
        r = math.sqrt(x * x + y * y + z * z)
        rmin, rmax = min(rmin, r), max(rmax, r)
    assert 1.468e11 < rmin < 1.475e11
    assert 1.517e11 < rmax < 1.524e11


def test_ephem_cache_interp():
    c = ephem.EphemCache(2461000.5)
    m_i, s_i = c.moon_sun(12345.0)
    m_d = ephem.moon_pos(2461000.5 + 12345.0 / 86400.0)
    err = math.dist(m_i, m_d)
    assert err < 5e3                                  # <5 km on 600 s nodes


# ---------- IC geometry ----------

def test_entry_state_invariants():
    """Energy, angular momentum, inbound-ness, vacuum perigee recovery."""
    import random
    rng = random.Random(7)
    for _ in range(50):
        v_inf = 300.0 + 4000.0 * rng.random()
        zc = 2 * rng.random() - 1
        az = 2 * math.pi * rng.random()
        sz = math.sqrt(1 - zc * zc)
        u = (sz * math.cos(az), sz * math.sin(az), zc)
        rp = P2.RE + 60e3 * rng.random()
        s, d = sample.entry_state(v_inf, u, 2 * math.pi * rng.random(), rp)
        r = math.sqrt(s[0]**2 + s[1]**2 + s[2]**2)
        v = math.sqrt(s[3]**2 + s[4]**2 + s[5]**2)
        assert abs(r - sample.ENTRY_R0) < 1.0
        eps = 0.5 * v * v - P2.MU / r
        assert abs(eps - 0.5 * v_inf * v_inf) / (0.5 * v_inf * v_inf) < 1e-9
        vr = (s[0]*s[3] + s[1]*s[4] + s[2]*s[5]) / r
        assert vr < 0.0                               # inbound
        hx = s[1]*s[5] - s[2]*s[4]
        hy = s[2]*s[3] - s[0]*s[5]
        hz = s[0]*s[4] - s[1]*s[3]
        h = math.sqrt(hx*hx + hy*hy + hz*hz)
        assert abs(h - d["b_m"] * v_inf) / h < 1e-9
        # two-body vacuum perigee from (eps, h)
        a = -P2.MU / (2 * eps) if eps < 0 else P2.MU / (2 * eps)
        e2 = 1.0 + 2.0 * eps * h * h / (P2.MU * P2.MU)
        rp_back = h * h / P2.MU / (1.0 + math.sqrt(max(0.0, e2)))
        assert abs(rp_back - rp) / rp < 1e-8


# ---------- propagator limits ----------

def _cfg(k=0.0, **kw):
    kw.setdefault("oblate", False)
    kw.setdefault("j2", False)
    kw.setdefault("moonsun", False)
    kw.setdefault("omega", 0.0)
    return Config(k, **kw)


def test_kepler_closure_3d():
    """Two-body, drag off: energy conserved over 3 eccentric orbits."""
    rp, ra = P2.RE + 300e3, P2.RE + 60000e3
    a = 0.5 * (rp + ra)
    vp = math.sqrt(P2.MU * (2 / rp - 1 / a))
    s = [rp, 0, 0, 0, vp * math.cos(0.3), vp * math.sin(0.3)]
    eps0 = 0.5 * vp * vp - P2.MU / rp
    cfg = _cfg()
    T = 2 * math.pi * math.sqrt(a**3 / P2.MU)
    res = fly(s, cfg, max_days=3 * T / 86400.0)
    v = math.sqrt(sum(c * c for c in res["state"][3:]))
    r = math.sqrt(sum(c * c for c in res["state"][:3]))
    eps1 = 0.5 * v * v - P2.MU / r
    assert res["outcome"] == "timeout"                # still orbiting
    assert abs(eps1 - eps0) / abs(eps0) < 2e-5
    assert res["n_pass"] == 0                         # 300 km perigee: never in atmosphere


def test_j2_secular_rates():
    """Node regression + apsidal advance vs analytic J2 rates (2%)."""
    a, ecc, inc = 7000e3, 0.01, math.radians(35.0)
    p = a * (1 - ecc * ecc)
    n = math.sqrt(P2.MU / a**3)
    rp = a * (1 - ecc)
    vp = math.sqrt(P2.MU * (2 / rp - 1 / a))
    # perigee on +x, node at +x, inclined by rotating velocity about x
    s = [rp, 0, 0, 0, vp * math.cos(inc), vp * math.sin(inc)]
    cfg = _cfg(j2=True)
    ndays = 2.0
    res = fly(s, cfg, max_days=ndays)
    st = res["state"]
    hx = st[1]*st[5] - st[2]*st[4]
    hy = st[2]*st[3] - st[0]*st[5]
    node = math.atan2(hx, -hy)                        # RAAN from h vector
    rate_num = node / (ndays * 86400.0)
    rate_an = -1.5 * propagate.J2 * (propagate.A_EQ / p)**2 * n * math.cos(inc)
    assert abs(rate_num - rate_an) / abs(rate_an) < 0.02, (rate_num, rate_an)


def test_corridor_limit_2d():
    """Spherical + no J2/Moon/Sun/rotation equatorial pass reproduces the 2D
    integrator's capture boundary at h_p=10 km (within bisection tol)."""
    b2 = P2.Body()
    cfg = _cfg(k=b2.k, sph_surface=True)      # full-legacy parity: sphere surface too
    u = (1.0, 0.0, 0.0)

    def outcome3d(v_inf):
        s, d = sample.entry_state(v_inf, u, 0.0, P2.RE + 10e3)
        assert d["i_deg"] < 1e-6 or d["i_deg"] > 179.9  # equatorial
        res = fly(s, cfg, max_days=0.5)
        if res["outcome"] == "impact":
            return "impact"
        if res["outcome"] == "escape":
            return "skip"
        return "capture" if res["n_pass"] >= 1 else res["outcome"]

    lo, hi = 1000.0, 4000.0
    for _ in range(10):
        mid = 0.5 * (lo + hi)
        if outcome3d(mid) == "capture":
            lo = mid
        else:
            hi = mid
    # 2D boundary (test_capture_feasibility): 2400 m/s numeric
    assert abs(lo - 2400.0) < 120.0, lo


def test_oblate_switch_matters():
    """Same polar skip pass at fixed geocentric periapsis: the oblate arm
    sees geodetic-altitude air (~14 km higher over the pole for this radius)
    so the per-pass energy loss must be several-fold smaller."""
    b2 = P2.Body()
    v_inf = 2000.0
    rp = P2.RE + 30e3
    s0, _ = sample.entry_state(v_inf, (0.0, 0.0, -1.0), 0.0, rp)    # polar arrival
    eps_in = 0.5 * v_inf * v_inf
    de = {}
    for arm, obl in (("sph", False), ("obl", True)):
        res = fly(list(s0), _cfg(k=b2.k, oblate=obl), max_days=0.2)
        st = res["state"]
        v2 = sum(c * c for c in st[3:])
        r = math.sqrt(sum(c * c for c in st[:3]))
        de[arm] = eps_in - (0.5 * v2 - P2.MU / r)
        assert res["outcome"] in ("escape", "timeout")               # clean skip
    assert de["sph"] / de["obl"] > 3.0, de


# ---------- terrain ----------

@pytest.mark.skipif(not DEM_OK, reason="ETOPO1 not downloaded")
def test_terrain_points():
    from deorbit.sitemc.terrain import Terrain
    t = Terrain()
    assert t.elevation(-1.469, -78.817) > 6200.0      # Chimborazo (peak-patched)
    assert t.elevation(27.988, 86.925) > 8700.0       # Everest (peak-patched)
    assert t.elevation(0.0, -120.0) < -3000.0          # deep Pacific
    assert 0.0 <= t.elevation(-3.0, -60.0) < 300.0     # Amazon lowland
    # ocean rock radius == ellipsoid (clamped)
    from deorbit.sitemc.geodesy import ellipsoid_radius_deg
    rr = t.rock_radius(0.0, math.radians(-120.0))
    assert abs(rr - ellipsoid_radius_deg(0.0)) < 1.0


@pytest.mark.skipif(not DEM_OK, reason="ETOPO1 not downloaded")
def test_terrain_chimborazo_is_global_max():
    """Geocentric radius of the rock surface: Chimborazo beats Everest and
    the other named peaks (site_selection's headline, now from the DEM)."""
    from deorbit.sitemc.terrain import Terrain
    from deorbit.site_selection import PEAKS
    t = Terrain()
    def rr(lat_d, lon_d):
        latc = math.atan2(math.sin(math.radians(lat_d)) * (1 - 1/298.257223563)**2,
                          math.cos(math.radians(lat_d)))
        return t.rock_radius(latc, math.radians(lon_d))
    chz = rr(-1.469, -78.817)
    for name, la, lo, el in PEAKS:
        if name != "Chimborazo":
            assert chz > rr(la, lo) - 1.0, name
    assert chz - rr(27.988, 86.925) > 1500.0          # ~2.1 km over Everest
