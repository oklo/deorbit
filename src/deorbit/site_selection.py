#!/usr/bin/env python3
"""Impact-site selection for an ultra-grazing de-orbiting body (pure stdlib).

Physics (established earlier in this project): a body that decays *slowly* and
grazes nearly tangent to the figure of the Earth impacts at the point of
MAXIMUM GEOCENTRIC RADIUS that its ground track can reach. The trajectory is a
conic about Earth's CENTER, so the selecting metric is geocentric radius -- and
there the ~21 km equatorial bulge dwarfs the ~9 km of topography. As the orbit
slowly decays its periapsis radius descends while nodal precession sweeps the
ground track across longitudes, so the body samples many great circles within
its latitude band (|lat| <= inclination) and homes in on the GLOBAL geocentric
maximum in that band.

This module:
  1. Geocentric surface radius = exact WGS84 ellipsoid + topography.
  2. Validates the pipeline against the known fact that Chimborazo is the point
     farthest from Earth's center (beating Everest by ~2 km).
  3. Maps an orbit plane (inclination, node longitude) -> the max-geocentric-
     radius point its ground track reaches = the predicted impact site, and the
     nodal-precession scan -> the global optimum in the latitude band.

Topography is a curated table of the highest-geocentric-radius summits (a small,
known set). A full DEM is the stress-test upgrade.
"""
import math

# WGS84
A_EQ = 6378137.0
F = 1.0 / 298.257223563
B_POL = A_EQ * (1.0 - F)

# (name, geodetic lat deg, lon deg [E+], elevation m)
PEAKS = [
    ("Chimborazo",      -1.469,  -78.817, 6263),
    ("Cayambe",          0.029,  -77.986, 5790),
    ("Cotopaxi",        -0.684,  -78.438, 5897),
    ("Huascaran Sur",   -9.121,  -77.604, 6768),
    ("Yerupaja",       -10.270,  -76.910, 6635),
    ("Kilimanjaro",     -3.066,   37.355, 5895),
    ("Mount Kenya",     -0.151,   37.307, 5199),
    ("Cerro Bonete",   -27.967,  -68.758, 6759),
    ("Aconcagua",      -32.653,  -70.011, 6961),
    ("Everest",         27.988,   86.925, 8849),
    ("Denali",          63.069, -151.007, 6190),
    ("Pico Bolivar",     8.541,  -71.046, 4978),
]


def ellipsoid_geocentric_radius(lat_deg):
    """Geocentric radius (m) of the WGS84 ellipsoid surface at geodetic lat."""
    p = math.radians(lat_deg)
    c, s = math.cos(p), math.sin(p)
    num = (A_EQ ** 2 * c) ** 2 + (B_POL ** 2 * s) ** 2
    den = (A_EQ * c) ** 2 + (B_POL * s) ** 2
    return math.sqrt(num / den)


def ang_dist_deg(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def topo(lat_deg, lon_deg, width_deg=1.0):
    """Topographic height (m) from the peak table, each peak a ~width_deg bump
    (mountain-range scale ~111 km). Crude stand-in for a DEM."""
    h = 0.0
    for _, la, lo, el in PEAKS:
        d = ang_dist_deg(lat_deg, lon_deg, la, lo)
        h = max(h, el * math.exp(-(d / width_deg) ** 2))
    return h


def surface_geocentric_radius(lat_deg, lon_deg):
    return ellipsoid_geocentric_radius(lat_deg) + topo(lat_deg, lon_deg)


def peak_geocentric_radius(lat_deg, elev_m):
    return ellipsoid_geocentric_radius(lat_deg) + elev_m


def rank_peaks():
    rows = [(n, peak_geocentric_radius(la, el), la) for n, la, lo, el in PEAKS]
    return sorted(rows, key=lambda r: -r[1])


def ground_track(inc_deg, node_lon_deg, n=2880):
    """Points (lat, lon) along the great circle with given inclination and
    Earth-fixed ascending-node longitude."""
    i = math.radians(inc_deg)
    pts = []
    for k in range(n):
        u = 2 * math.pi * k / n                    # argument of latitude
        lat = math.degrees(math.asin(math.sin(i) * math.sin(u)))
        dlon = math.degrees(math.atan2(math.cos(i) * math.sin(u), math.cos(u)))
        lon = ((node_lon_deg + dlon + 180.0) % 360.0) - 180.0
        pts.append((lat, lon))
    return pts


def impact_site(inc_deg, node_lon_deg):
    """Max-geocentric-radius point the ground track reaches = predicted impact
    site for a slowly-decaying grazing body."""
    best = None
    for lat, lon in ground_track(inc_deg, node_lon_deg):
        R = surface_geocentric_radius(lat, lon)
        if best is None or R > best[2]:
            best = (lat, lon, R)
    lat, lon, R = best
    near = min(PEAKS, key=lambda p: ang_dist_deg(lat, lon, p[1], p[2]))
    return dict(lat=lat, lon=lon, R_km=R / 1e3, nearest=near[0],
                nearest_dist_deg=ang_dist_deg(lat, lon, near[1], near[2]))


def precession_scan(inc_deg, n_nodes=360):
    """Slow decay + nodal precession sweeps the node longitude over all values;
    the body homes in on the global geocentric optimum in its latitude band."""
    best = None
    for j in range(n_nodes):
        s = impact_site(inc_deg, -180.0 + 360.0 * j / n_nodes)
        if best is None or s["R_km"] > best["R_km"]:
            best = s
    return best


def bulge_reach(fpa_deg, dR_km):
    """A grazing trajectory at flight-path angle fpa descends dR over an along-
    track distance dR/tan(fpa); return that distance (km) and the latitude span."""
    along_km = dR_km / math.tan(math.radians(fpa_deg))
    dlat_deg = along_km / 111.32
    return along_km, dlat_deg


if __name__ == "__main__":
    print("=== 1. Geocentric-radius ranking of high summits (validates pipeline) ===")
    print(f"  WGS84 equatorial bulge: R_eq - R_pole = "
          f"{(ellipsoid_geocentric_radius(0) - ellipsoid_geocentric_radius(90))/1e3:.1f} km"
          f"  (vs topography <~ 9 km)")
    for n, R, la in rank_peaks():
        print(f"    {n:14s} lat={la:7.2f}  geocentric R = {R/1e3:10.3f} km")
    top = rank_peaks()[0]
    print(f"  -> global max: {top[0]} ({top[1]/1e3:.3f} km)  "
          f"{'[CHIMBORAZO as expected]' if top[0]=='Chimborazo' else '[!!]'}")

    print("\n=== 2. Bulge dominance / grazing steering reach ===")
    for lat in (1, 2, 5, 10, 30):
        dR = (ellipsoid_geocentric_radius(0) - ellipsoid_geocentric_radius(lat)) / 1e3
        print(f"    bulge drop equator->{lat:2d} deg lat = {dR:6.3f} km")
    for fpa in (1.0, 2.0, 4.0):
        a, dl = bulge_reach(fpa, 6.0)   # 6 km ~ Andes topographic relief
        print(f"    fpa={fpa:.0f} deg: a 6 km surface rise is reached over "
              f"{a:6.0f} km along-track (~{dl:.1f} deg lat) -> body is drawn to it")

    print("\n=== 3. Trajectory -> impact site ===")
    print("  near-equatorial track (i=5 deg) with node over the Andes (-78 deg):")
    s = impact_site(5.0, -78.0)
    print(f"    impact lat={s['lat']:.2f}, lon={s['lon']:.2f}, R={s['R_km']:.3f} km "
          f"-> nearest summit: {s['nearest']} ({s['nearest_dist_deg']:.2f} deg)")
    print("  nodal-precession scan (all node longitudes) per inclination:")
    for inc in (1.0, 2.0, 5.0, 10.0, 23.5):
        b = precession_scan(inc)
        print(f"    i={inc:5.1f} deg -> global site lat={b['lat']:6.2f} lon={b['lon']:8.2f} "
              f"R={b['R_km']:.3f} km  nearest: {b['nearest']}")
