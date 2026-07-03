#!/usr/bin/env python3
"""Gates for the site-MC geodesy/oblate-atmosphere module (pure stdlib)."""
import math

from deorbit.corridor.physics import density as rho_profile, RE as RE_SPHERE
from deorbit.sitemc import geodesy as G


def test_wgs84_radii():
    assert abs(G.ellipsoid_radius_deg(0.0) - 6378137.0) < 1.0
    assert abs(G.ellipsoid_radius_deg(90.0) - 6356752.3) < 1.0
    assert abs(G.ellipsoid_radius_deg(-90.0) - 6356752.3) < 1.0
    bulge = G.ellipsoid_radius_deg(0.0) - G.ellipsoid_radius_deg(90.0)
    assert abs(bulge - 21384.7) < 1.0          # the 21.4 km lever arm


def test_site_selection_consistency():
    """Same ellipsoid as the site-selection module (one geoid in the repo).
    site_selection takes GEODETIC latitude; geodesy takes GEOCENTRIC (what a
    state vector gives) — convert, then the radii must agree to sub-meter."""
    from deorbit import site_selection as S
    for lat_d in (0.0, -1.469, 23.5, 45.0, 90.0):   # incl. Chimborazo's latitude
        p = math.radians(lat_d)
        lat_c = math.degrees(math.atan2((G.B_POL / G.A_EQ) ** 2 * math.sin(p), math.cos(p)))
        assert abs(S.ellipsoid_geocentric_radius(lat_d) - G.ellipsoid_radius_deg(lat_c)) < 0.5


def test_oblate_density_contrast():
    """At FIXED geocentric radius (equatorial surface + 8 km), the pole sits
    ~29 km up: density must be ~20x lower there (3 scale heights)."""
    r = G.A_EQ + 8000.0
    rho_eq = G.density_oblate(r, 0.0, 0.0)
    rho_pole = G.density_oblate(0.0, 0.0, r)
    assert abs(rho_eq / rho_profile(8000.0) - 1.0) < 1e-12   # equator == profile at 8 km
    ratio = rho_eq / rho_pole
    assert 15.0 < ratio < 30.0, f"equator/pole density ratio {ratio:.1f}"


def test_spherical_fallback_exact():
    """oblate=False reproduces the legacy spherical model bit-for-bit."""
    for r, zfrac in ((RE_SPHERE + 8e3, 0.3), (RE_SPHERE + 50e3, 0.9), (RE_SPHERE + 120e3, 0.0)):
        z = r * zfrac
        x = math.sqrt(r * r - z * z)
        assert G.density_oblate(x, 0.0, z, oblate=False) == rho_profile(r - RE_SPHERE)


def test_geodetic_alt_roundtrip():
    """A point built at altitude h above the ellipsoid at latitude lat must
    return ~h (radial-vs-normal error << scale height)."""
    for lat in (0.0, -30.0, 60.0):
        for h in (0.0, 8e3, 100e3):
            la = math.radians(lat)
            r = G.ellipsoid_radius_deg(lat) + h
            x, z = r * math.cos(la), r * math.sin(la)
            assert abs(G.geodetic_alt(x, 0.0, z) - h) < 25.0
