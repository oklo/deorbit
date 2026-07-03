# Site-MC findings log

## Phase B (2026-07-03): oblate-vs-spherical A/B + the discrimination race

Setup: 30 matched draws per inclination (v_inf 0.2-0.8 km/s, rp_alt 2-12 km,
random lunar epoch), i = 5/25/45/65 deg, both arms with the REAL
ellipsoid+DEM impact surface; only the ATMOSPHERE'S latitude dependence
toggles (spherical arm referenced to A_EQ so an equatorial orbit is
identical in both arms). Scratch: scratchpad phaseb_ab.py.

| i | arm | deorbit impacts | med abs(lat) | q90 abs(lat) | grazing | med low passes |
|---|---|---|---|---|---|---|
| 5  | oblate | 25 | 0.83 | 4.1  | 25/25 | 2 |
| 5  | sph    | 25 | 0.79 | 3.8  | 25/25 | 2 |
| 25 | oblate | 25 | 3.24 | 19.9 | 25/25 | 2 |
| 25 | sph    | 25 | 5.43 | 13.9 | 25/25 | 2 |
| 45 | oblate | 25 | 6.02 | 25.6 | 25/25 | 2 |
| 45 | sph    | 26 | 7.27 | 21.8 | 26/26 | 2 |
| 65 | oblate | 25 | 6.59 | 18.4 | 25/25 | 2 |
| 65 | sph    | 25 | 9.37 | 28.7 | 25/25 | 3 |

**Finding B1 — impact latitude concentrates equatorward FAR inside the
inclination band, and mostly NOT because of the oblate atmosphere.** Even
with spherical air, i=65 deg impacts have median |lat| ~ 9 deg (band edge
65). Mechanism: under J2 the near-circular decaying orbit's radius varies
with latitude by only ~J2-scale (~7-9 km, orbit "partially rides" the
bulge), while the solid surface drops 21.4 km equator->pole — so closest
approach, peak drag, and death all migrate to the low-latitude part of the
track. This REFINES the original geocentric-radius site-selection argument:
J2 shrinks the effective bulge lever from 21.4 km to ~12-14 km — still well
above the 9 km of topography, so the equatorward steering survives, and it
operates through the ORBIT SHAPE, not only through which terrain is tallest.

**Finding B2 — the oblate atmosphere is a secondary but real tightener**
(median |lat| shrinks a further ~30-40% at i>=25: 5.4->3.2 at i=25,
9.4->6.6 at i=65). The spherical-atmosphere corridor results are therefore
not wrong about capture ENERGETICS, but any site statistics computed with a
spherical atmosphere would be biased poleward.

**Finding B3 — the discrimination race is fast, as estimated.** Median LOW
(terrain-window) passes = 2 across all inclinations: the body does NOT
slowly comb the latitude band for the tallest summit; it dies within ~2
window passes + one sub-orbital plunge. Longitude sampling is a handful of
long arcs, so summit-scale cross-section enhancement (if any) must come
from the plunge geometry (~90-350 km terminal footprint), not from
band-max search. The original "precession samples the band before decay"
assumption is quantitatively FALSE; equatorward steering happens anyway,
via B1/B2.

**B4 — terminal state is robustly ultra-oblique:** 201/201 deorbit impacts
across all Phase B arms have |fpa| < 5 deg (typ. 0.7-1.3 deg at ~7.3-7.7
km/s) — matches the hi-fi reference (7.39 km/s, -1.16 deg).

Engine bug found by Phase B and fixed before any statistics: the pass
counter triggered on vr wobbles within a single skimming arc and apo1 was
measured mid-pass (159,000 "km" for a 3-day cascade); passes are now
interface-excursions with apogee at atmosphere EXIT. Also fixed: the
spherical A/B arm originally impacted a phantom A_EQ sphere (21 km above
the real geoid at the poles) — the impact surface is now ALWAYS the
ellipsoid+DEM (cfg.sph_surface exists only for the 2D-parity gate).

## Phase C pilot (300 seeds, full prior box)

escape 273 / direct impact 20 / deorbit impact 4 / timeout 3 (200-day cap).
Deorbit-impact rate ~1.3% of arrivals in the (v_inf<=5 km/s, rp<=60 km,
isotropic) box; all 4 grazing; 3/4 with apo1 < 130,000 km. Throughput ~65
seeds/s on 7 workers (pure Python) -> 1e6 seeds ~ 4-5 h.
