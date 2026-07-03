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

## Early-catalog scorecard (2026-07-03 ~14:15, 42.9k seeds, 1,353 deorbit impacts)

Pre-registered predictions (conversation, before peeking) vs data:

| prediction | verdict |
|---|---|
| Latitude law dominates | DIRECTION YES (med abs(lat) 13.7 deg vs 30 isotropic-null) but SOFTER than predicted (predicted 50%@10/90%@25-30; got 13.7/34.4). Cause: Phase B sampled only low-v_inf deep captures at i<=65; the isotropic prior's high-i captures widen the belt. |
| Longitude ~uniform; Pacific plurality; land ~ area share | YES: 30-deg bins 94-146 (flat); Pacific-ish 42%; land 28.1% (pred 25-30); eta(ocean)=0.97, eta(0-1km)=1.00. |
| eta rises with elevation; 3-5km ~2-6; >5km CI wide but >1 | YES: 1.00 / 1.14 / 3.12 [2.08-4.04] / 4.19 [1.40-9.07]. Range-scale enhancement is REAL (>5km CI excludes 1). 1-3km came in LOW (1.14 vs pred 1.5-2.5). |
| Andes dominate highlands; Himalaya excluded by latitude; New Guinea sleeper | HALF-WRONG, interestingly: elev>2km hits split Andes 16-18 vs Tibet/Himalaya 16-17 — and the Tibetan hits are GRAZING-channel (16/17), delivered by high-inclination captures reaching 28-35N. TIBET IS THE SECOND HIGHLAND ATTRACTOR (huge 4-5 km AREA beats wall geometry at ensemble level). New Guinea: only 2 (flop so far). Isolated cones (E Africa 1) underperform as predicted. |
| Trophy summits: no summit-cone enhancement (Greg) | HOLDING: 0 impacts near Chimborazo/Cayambe in 1,353; named-peak <150 km proximity total 7 (~0.5%), led by Huascaran 3 (closest 36.6 km at 5,126 m elev — the latitude-sweet-spot call). |
| ~100% grazing; apo1<130k ~75% | 90.5% grazing (over-claimed); the 9.5% steep have apo1>130k in 130/130 cases — an EXACT, independent reproduction of the hi-fi lunar-pumping boundary by the stdlib engine. apo1<130k = 72%. |
| Slight retrograde excess (co-rotation) | +1.8 pp (51.6 vs 49.8% in escapes), right direction, ~1.4 sigma — pending full N. |

Emerging map narrative: an equatorial-belt-dominant, ocean-plurality map with
TWO highland attractors — the Andes (wall geometry, low-lat channel) and
Tibet (area, high-i grazing channel) — real order-3-4 range-scale elevation
enhancement, and no summit-cone enhancement. Timeouts 2.6% (200-day cap) to
revisit at full N.

## Phase E (2026-07-03): arrival prior built — THE MAP IS PRIOR-ROBUST

Machinery (src/deorbit/sitemc/arrival_prior.py, gated): empirical torus-MC
encounter kinematics (no Opik closed forms to get wrong; every hit checks
against Tisserand, median err 0.9%), pluggable population density (flat vs
Granvik-shaped tilt), KDE over local (U_r,U_t,U_z), per-seed weights
w ~ f(U)*U^3*db^2/drp with u_hat regenerated deterministically from the seed.

RESULT: the realistic arrival prior BARELY MOVES the site statistics.
ESS 84-89% (weights bulk within ~3x); weighted |lat| q50 = 14.3 (tilt) /
14.6 (flat) vs 13.7 isotropic; elev>2km highland share Andes ~28-29% vs
Tibet/Himalaya ~34-35% under BOTH densities — Tibet SURVIVES the realistic
prior. Physics: Tisserand pins sqrt(p)cos(i)~1 for U<5 km/s but leaves the
U-direction broad (radial/tangential/vertical branches all populated), so
low-U radiants are NOT an ecliptic-thin cone and high-geocentric-inclination
captures remain fully supported. Reweighted flashlight map (rejection,
acceptance 32%): figures/sitemc_flashlight_rw.png — visually the isotropic
beam. => The impact-site map is DYNAMICS-dominated, not prior-dominated;
prior uncertainty (Granvik vs NEOMOD vs M-source conditioning) is a
second-order systematic for the map shape (it matters for the RATE only).
Caveats: analytic tilt is Granvik-shaped placeholder (NEOMOD grid ingestion
open); KDE h=0.4 km/s coarse below v_inf~0.5; Earth eccentricity ignored.

## FULL-CATALOG HARVEST (2026-07-03 evening; 1e6 seeds, 32,734 deorbit impacts)

Outcomes: escape 898.0k / direct impact 43.3k / deorbit impact 32.7k /
timeout 26.0k (200-day cap; censoring caveat, 2.6%).

**Final eta ladder (latitude-controlled, bootstrap 95% CI):**
| class | P_impact | P_area(lat) | eta | CI |
|---|---|---|---|---|
| ocean | 0.7145 | 0.7417 | 0.96 | [0.96, 0.97] |
| 0-1 km | 0.2045 | 0.2070 | 0.99 | [0.97, 1.01] |
| 1-3 km | 0.0627 | 0.0431 | 1.45 | [1.38, 1.52] |
| 3-5 km | 0.0151 | 0.0070 | 2.15 | [1.97, 2.33] |
| >5 km | 0.0032 | 0.0012 | 2.75 | [2.25, 3.28] |

|lat| quartiles 7.5 / 15.2 / 25.1 / 36.0 deg. Grazing 91.7%; apo1<130k
75.3%; land fraction 28.5% (vs 25.8% area null -> overall land
enhancement only ~1.11x). Weighted (realistic prior, ESS 88.7%):
quartiles 7.2/14.8/24.3/35.0 — unchanged; elev>2km share Andes 34.3%
vs Tibet/Himalaya 29.2% (isotropic ~even) — the two-attractor structure
holds under the realistic prior with a mild Andes lead.

**Named-summit census (impacts within 150 km): CAYAMBE LEADS with 26**,
then Cerro Bonete 17, Everest 17, Pico Bolivar 15, Huascaran 14,
Yerupaja 14, Kilimanjaro 13, Cotopaxi 12, Mt Kenya 9, Chimborazo 7,
Aconcagua 7, Denali 1. Cayambe's lead is LATITUDE-density driven (it
sits at 0.03 deg N, the peak of the impact-latitude distribution), not
height-driven — Chimborazo (taller, 1.5 deg S) draws 4x fewer.
P(within 150 km of Cayambe) = 8e-4 per deorbit impact.

**Pre-registered hypothesis verdict:** split, leaning Greg. At the
summit-cone scale there is NO enhancement beyond the range-scale
factor (trophy peaks collect ~area x eta(elev-class) x latitude-density
— no extra cross-section). The range-scale elevation steering is REAL
but MODEST: eta rises monotonically to 2.15 (3-5 km) and 2.75 (>5 km),
at the LOW end of the counter-prediction's 3-4x. Ocean is significantly
but barely suppressed (0.96). The dominant selection remains LATITUDE
(bulge, via J2 orbit-shape + oblate air), not elevation.

Showcase mechanism case (seed 402938, figures/sitemc_traj*): Pacific
death-run -> caught by the west-facing Cordillera Occidental rampart at
3,159 m (0.53 N, 78.48 W, Pinan paramo NW of Cotacachi), where terrain
rose ~2,000 m under the final 30 km of a -0.84 deg glide (~4 deg wall
vs ~1 deg glider) — the wall does the catching, not the trajectory.
