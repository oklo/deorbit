"""WGS84 ellipsoid geometry + oblate-atmosphere density (pure stdlib).

The site-selection Monte Carlo's load-bearing physics upgrade
(docs/site_mc_plan.md): atmospheric drag must be evaluated at GEODETIC
altitude above the ellipsoid, not r - R_sphere. The 21.4 km equator-pole
radius difference is ~3 scale heights, i.e. a ~20x density contrast at fixed
geocentric radius — first-order for where an inclined decaying orbit dies.

Altitude approximation: h = r - R_ell(phi_c), the radial height above the
ellipsoid surface point at the same geocentric latitude. This differs from
the true normal (geodetic) altitude by <~20 m for h < 200 km (the surface
normal tilts <=0.19 deg from radial) — negligible against a 7 km scale
height. The J2 used by the integrator must correspond to the same
flattening (consistency gate in tests/test_sitemc_geodesy.py).
"""
import math

from deorbit.corridor.physics import density as _density_profile, RE as RE_SPHERE

# WGS84
A_EQ = 6378137.0                 # equatorial radius, m
F = 1.0 / 298.257223563          # flattening
B_POL = A_EQ * (1.0 - F)         # polar radius, m


def ellipsoid_radius(sin_latc, cos_latc):
    """Geocentric radius of the WGS84 surface at geocentric latitude
    (given as sin/cos so the integrator never calls atan2/asin)."""
    return (A_EQ * B_POL) / math.sqrt((A_EQ * sin_latc) ** 2 + (B_POL * cos_latc) ** 2)


def ellipsoid_radius_deg(latc_deg):
    la = math.radians(latc_deg)
    return ellipsoid_radius(math.sin(la), math.cos(la))


def geodetic_alt(x, y, z):
    """Radial altitude above the WGS84 ellipsoid for a geocentric ECEF/ECI
    position (m). (Latitude here is geocentric; see module docstring.)"""
    r = math.sqrt(x * x + y * y + z * z)
    sin_latc = z / r
    cos_latc = math.sqrt(max(0.0, 1.0 - sin_latc * sin_latc))
    return r - ellipsoid_radius(sin_latc, cos_latc)


def density_oblate(x, y, z, f_rho=1.0, f_H=1.0, oblate=True):
    """Atmospheric density at a 3D position: the corridor's Vallado profile
    evaluated at geodetic altitude (oblate=True) or at r - RE_sphere
    (oblate=False, the legacy spherical model for A/B isolation)."""
    if oblate:
        h = geodetic_alt(x, y, z)
    else:
        h = math.sqrt(x * x + y * y + z * z) - RE_SPHERE
    return _density_profile(h, f_rho, f_H)
