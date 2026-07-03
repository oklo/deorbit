"""ETOPO1 terrain for the site MC (pure stdlib, mmap).

Data: data/dem/etopo1_ice_g_i2.bin — grid-registered global relief, int16 LSB,
21601 x 10801, cell 1/60 deg, row 0 = lat +90, col 0 = lon -180 (see .hdr;
sha256 in data/dem/README). Ice-surface elevations, bathymetry negative.

Known systematic: a 1-arc-min grid shaves sharp summits (file max 8271 m vs
Everest 8849). Fix: the curated site_selection.PEAKS are patched in as cones
(peak elev at the summit, sloping at PEAK_SLOPE back to the DEM), so summit
CATCH cross-sections are not silently deflated. PEAK_SLOPE is a systematic to
vary (0.10/0.15/0.20).

rock_radius(latc, lon_e) returns the geocentric radius of the solid/liquid
surface (ocean surface = ellipsoid: bathymetry clamped to 0) — the impact
catch test is trajectory r <= rock_radius.
"""
import math
import mmap
import os
import struct

from deorbit.sitemc.geodesy import A_EQ, B_POL, ellipsoid_radius

NC, NR = 21601, 10801
CELL = 1.0 / 60.0
PEAK_SLOPE = 0.15
_E2FAC = (B_POL / A_EQ) ** 2          # geocentric <-> geodetic tan factor
DEM_PATH = os.environ.get(
    "DEORBIT_DEM_ETOPO",
    os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))), "data", "dem", "etopo1_ice_g_i2.bin"))


class Terrain:
    def __init__(self, path=DEM_PATH, peaks=True):
        f = open(path, "rb")
        self.mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        assert len(self.mm) == NC * NR * 2, "ETOPO1 size mismatch"
        self.peaks = []
        if peaks:
            from deorbit.site_selection import PEAKS
            for name, lat, lon, elev in PEAKS:
                self.peaks.append((name, lat, lon, float(elev),
                                   math.cos(math.radians(lat))))

    def elevation(self, lat_deg, lon_deg):
        """Bilinear DEM elevation (m) at geodetic lat/lon; peaks patched."""
        rr = (90.0 - lat_deg) / CELL
        cc = ((lon_deg + 180.0) % 360.0) / CELL
        r0 = min(NR - 2, max(0, int(rr)))
        c0 = int(cc) % (NC - 1)
        fr, fc = rr - r0, cc - int(cc)
        c1 = (c0 + 1) % (NC - 1)
        u = struct.unpack_from
        mm = self.mm
        v00 = u("<h", mm, 2 * (r0 * NC + c0))[0]
        v01 = u("<h", mm, 2 * (r0 * NC + c1))[0]
        v10 = u("<h", mm, 2 * ((r0 + 1) * NC + c0))[0]
        v11 = u("<h", mm, 2 * ((r0 + 1) * NC + c1))[0]
        e = (v00 * (1 - fr) * (1 - fc) + v01 * (1 - fr) * fc
             + v10 * fr * (1 - fc) + v11 * fr * fc)
        # summit-cone patch
        for name, plat, plon, pelev, pcos in self.peaks:
            if abs(lat_deg - plat) < 0.12:
                dlon = abs((lon_deg - plon + 180.0) % 360.0 - 180.0)
                if dlon < 0.12:
                    d_m = 111.32e3 * math.hypot(lat_deg - plat, dlon * pcos)
                    e = max(e, pelev - PEAK_SLOPE * d_m)
        return e

    def rock_radius(self, latc_rad, lon_e_rad):
        """Geocentric radius (m) of the surface under geocentric lat / Earth-
        fixed lon. Ocean surface = ellipsoid (elev clamped to >= 0)."""
        sl, cl = math.sin(latc_rad), math.cos(latc_rad)
        lat_d = math.degrees(math.atan2(sl, _E2FAC * cl))     # geodetic
        e = self.elevation(lat_d, math.degrees(lon_e_rad))
        return ellipsoid_radius(sl, cl) + max(0.0, e)

    def area_elev_hist(self, lat_lo, lat_hi, edges, stride=4):
        """cos(lat)-area-weighted elevation histogram over a latitude band
        (the eta-statistic null). Returns list of weights per bin defined by
        `edges` (ocean = elev <= 0 goes to bin 0 edge handling by caller)."""
        w = [0.0] * (len(edges) + 1)
        r_lo = max(0, int((90.0 - lat_hi) / CELL))
        r_hi = min(NR - 1, int((90.0 - lat_lo) / CELL))
        u = struct.unpack_from
        mm = self.mm
        for r in range(r_lo, r_hi + 1, stride):
            cw = math.cos(math.radians(90.0 - r * CELL))
            for c in range(0, NC - 1, stride):
                v = u("<h", mm, 2 * (r * NC + c))[0]
                i = 0
                for k, edge in enumerate(edges):
                    if v > edge:
                        i = k + 1
                w[i] += cw
        return w
