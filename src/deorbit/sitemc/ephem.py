"""Low-precision analytic Moon/Sun ephemerides (pure stdlib).

Montenbruck & Gill, "Satellite Orbits" §3.3.2 truncated series: Moon good to
~a few arcmin / ~500 km, Sun to ~0.01 deg. Far more than needed here — the
Moon's tidal effect on the spacecraft cares about the lunar POSITION only at
the 1e-3 level for MC statistics (ephemeris error just re-randomizes lunar
phase across draws). The high-fidelity REBOUND/ASSIST (DE440) layer is the
oracle for cross-checks.

Frame: geocentric equatorial inertial (mean equator/equinox of date ~ J2000;
precession ignored). Units m, s. Time: Julian date (UTC~TT, |dt|<70 s
irrelevant here).
"""
import math

J2000_JD = 2451545.0
MU_MOON = 4.9028e12       # m^3/s^2
MU_SUN = 1.32712440018e20
OBLIQ = math.radians(23.43929111)
_CE, _SE = math.cos(OBLIQ), math.sin(OBLIQ)
ARCSEC = 1.0 / 3600.0     # deg


def _ecl_to_eq(lam, bet, r):
    """Ecliptic lon/lat (rad) + radius -> equatorial xyz."""
    cl, sl = math.cos(lam), math.sin(lam)
    cb, sb = math.cos(bet), math.sin(bet)
    x_ec, y_ec, z_ec = r * cb * cl, r * cb * sl, r * sb
    return (x_ec, _CE * y_ec - _SE * z_ec, _SE * y_ec + _CE * z_ec)


def sun_pos(jd):
    """Geocentric equatorial Sun position (m)."""
    T = (jd - J2000_JD) / 36525.0
    M = math.radians((357.5256 + 35999.049 * T) % 360.0)
    lam = math.radians(282.9400) + M + math.radians(
        6892.0 * ARCSEC * math.sin(M) + 72.0 * ARCSEC * math.sin(2 * M))
    r = (149.619 - 2.499 * math.cos(M) - 0.021 * math.cos(2 * M)) * 1e9
    return _ecl_to_eq(lam, 0.0, r)


def moon_pos(jd):
    """Geocentric equatorial Moon position (m)."""
    T = (jd - J2000_JD) / 36525.0
    d2r = math.radians
    L0 = (218.31617 + 481267.88088 * T) % 360.0
    l = d2r((134.96292 + 477198.86753 * T) % 360.0)
    lp = d2r((357.52543 + 35999.04944 * T) % 360.0)
    F = d2r((93.27283 + 483202.01873 * T) % 360.0)
    D = d2r((297.85027 + 445267.11135 * T) % 360.0)
    dlam = (22640.0 * math.sin(l) + 769.0 * math.sin(2 * l)
            - 4586.0 * math.sin(l - 2 * D) + 2370.0 * math.sin(2 * D)
            - 668.0 * math.sin(lp) - 412.0 * math.sin(2 * F)
            - 212.0 * math.sin(2 * l - 2 * D) - 206.0 * math.sin(l + lp - 2 * D)
            + 192.0 * math.sin(l + 2 * D) - 165.0 * math.sin(lp - 2 * D)
            + 148.0 * math.sin(l - lp) - 125.0 * math.sin(D)
            - 110.0 * math.sin(l + lp) - 55.0 * math.sin(2 * F - 2 * D)) * ARCSEC
    lam = d2r((L0 + dlam) % 360.0)
    s = F + d2r(dlam + 412.0 * ARCSEC * math.sin(2 * F) + 541.0 * ARCSEC * math.sin(lp))
    bet = d2r((18520.0 * math.sin(s) - 526.0 * math.sin(F - 2 * D)
               + 44.0 * math.sin(l + F - 2 * D) - 31.0 * math.sin(-l + F - 2 * D)
               - 25.0 * math.sin(-2 * l + F) - 23.0 * math.sin(lp + F - 2 * D)
               + 21.0 * math.sin(-l + F) + 11.0 * math.sin(-lp + F - 2 * D)) * ARCSEC)
    r = (385000.0 - 20905.0 * math.cos(l) - 3699.0 * math.cos(2 * D - l)
         - 2956.0 * math.cos(2 * D) - 570.0 * math.cos(2 * l)
         + 246.0 * math.cos(2 * l - 2 * D) - 205.0 * math.cos(lp - 2 * D)
         - 171.0 * math.cos(l + 2 * D) - 152.0 * math.cos(l + lp - 2 * D)) * 1e3
    return _ecl_to_eq(lam, bet, r)


def gmst_rad(jd):
    """Greenwich mean sidereal time (rad)."""
    d = jd - J2000_JD
    return math.radians((280.46061837 + 360.98564736629 * d) % 360.0)


class EphemCache:
    """Piecewise-linear Moon/Sun positions on a lazy 600 s grid keyed to an
    epoch JD. Linear-interp error ~120 m for the Moon (a_moon*dt^2/8) —
    nothing. Cuts ~60 trig calls per force evaluation to ~12 mul/add."""
    DT = 600.0

    def __init__(self, jd0):
        self.jd0 = jd0
        self._nodes = {}

    def _node(self, i):
        n = self._nodes.get(i)
        if n is None:
            jd = self.jd0 + i * self.DT / 86400.0
            n = (moon_pos(jd), sun_pos(jd))
            self._nodes[i] = n
        return n

    def moon_sun(self, t):
        """Positions at t seconds past epoch."""
        fi = t / self.DT
        i = math.floor(fi)
        w = fi - i
        (m0, s0), (m1, s1) = self._node(i), self._node(i + 1)
        moon = tuple(a + w * (b - a) for a, b in zip(m0, m1))
        sun = tuple(a + w * (b - a) for a, b in zip(s0, s1))
        return moon, sun
