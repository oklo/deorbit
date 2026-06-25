# Findings — high-fidelity capture & aerobraking (REBOUND + ASSIST + drag)

All results below are from the full-fidelity engine (JPL DE440 ephemeris,
Earth J2/J3/J4, Sun/Moon/planets, GR, geocentric frame) with atmospheric drag
by operator-splitting, validated in `validate_hifi.py`. Body: 1 km iron
monolith, β ≈ 5.2×10⁶ kg/m². Verified, not diagnostic-only.

## 1. The Moon sets a sharp upper limit on the post-capture first apogee

Survival-map sweep (`aerobrake_sweep.py`), v∞ = 0.2 km/s, inc = 5°, over
periapsis (→ first apogee) and 4 Moon phases:

| first apogee | regime | impact flight-path angle |
|---|---|---|
| ≲ ~130,000 km (≲ ⅓ lunar distance) | **grazing impact** (target) | −0.2° to −3.7° |
| ~140,000 – 600,000 km | **steep impact** (Moon-pumped) | −6° → −68° (steepens with apogee) |
| ≳ ~2,000,000 km (> Earth Hill) | **ejected** | — |

- Boundary is **below** the Moon's orbit (384,400 km): periapsis pumping grows
  steeply with apogee, so the Moon destabilizes captures even when apogee never
  reaches lunar distance.
- **Phase-dependent**: usable boundary ~100,000–137,000 km depending on Moon phase.
- Deep-apogee cases (apo ~1.1×10⁶ km) impact at fpa −55° to −69° **and latitude
  17–38°** — the Moon completely reorients a 5° capture.

## 2. Pure-coast verification of lunar periapsis pumping

A barely-bound orbit (apo ~604,000 km, peri alt +42 km), drag OFF, one 19.7-day
orbit: osculating periapsis goes **+42 km → −2582 km**. One-shot integrate and
the stepped coast agree to the meter → real lunar dynamics, not a numerical
artifact. Periapsis only needs to drop ~42 km to impact; the Moon moves it ~2600.

## 3. The optimistic low-apogee case aerobrakes cleanly to a grazing impact

v∞ = 0.3 km/s, peri 8 km, inc 5° → first apogee 100,926 km (inside lunar orbit):
captures in one pass (Δv 327 m/s, q 24 MPa), then **66 passes / 26 days** to
impact at **7.39 km/s, fpa −1.16°, latitude 4.95°**. Inclination preserved
(5.00→4.96°; drag doesn't damp i); RAAN regresses 0→−16.8° (J2). No lunar
encounter (closest 355,000 km). A non-obvious self-limiting behavior appears:
the deep first pass over-circularizes, lifting periapsis out of the dense air,
and J2/lunar precession later walks it back down.

## Implication

The clean capture → aerobrake → **ultra-oblique grazing impact** chain works,
but only in a narrow, demanding corner: **low v∞** (slow encounter) AND a
**deep capture pass** that drops the first apogee below ~⅓ of the lunar distance.
Above that, the Moon pumps periapsis into a steep impact or ejects the body. This
quantifies the qualitative worry that "the Moon is creating an issue here."

## 4. Site selection: the bulge steers grazing impacts to the equatorial Andes

(`../site_selection.py`, pure stdlib: exact WGS84 ellipsoid + a curated table of
the highest-geocentric-radius summits.)

- **Pipeline validation:** ranking the world's high summits by geocentric radius
  reproduces the known fact — **Chimborazo is the global maximum** (6384.39 km),
  beating Everest (6382.31 km) by ~2.1 km; Huascarán essentially tied; Cayambe
  top-tier (6383.93 km). All leaders are equatorial-ish high mountains.
- **Bulge dominance:** equator→pole geocentric drop is 21.4 km vs ≤9 km topo. A
  grazing body (fpa 1–4°) descends slowly enough that a 6 km surface rise is
  reached over 86–344 km along-track — so it is drawn to nearby high terrain.
- **Result (slow decay + J2 precession of ω, Ω samples the |lat|≤i band):** the
  body impacts at the global geocentric-radius maximum within its latitude band:
    - i ≈ 1°  → Cotopaxi/Cayambe (lat ~0, the equator crossing over the Andes)
    - i ≈ 2–5° → **Chimborazo** (the global optimum)
    - i ≈ 10–23° → Huascarán Sur (tied geocentric radius at higher lat)
  A near-equatorial track with its node over the Andes impacts at **Cayambe**
  (lat 0.00, lon −78.0). So the equatorial bulge + Andes topography make
  **Chimborazo/Cayambe the bullseye** — validating the target choice as a
  geometric consequence, not an arbitrary pick.
- **Caveat to stress-test:** assumes nodal/apsidal precession samples the band
  faster than periapsis decays (true for slow aerobraking, to be confirmed); and
  topography is a peak-table stand-in, not a full DEM.

## Open / next

- Distinguish *captured-then-aerobraked* grazing impacts (many passes, ~7.4 km/s)
  from *direct* grazing impactors (1 pass, ~10–11 km/s) — both pass the fpa<5°
  test but are different histories.
- Refine the boundary (finer periapsis near 10–14 km), add v∞ and inclination
  dimensions, resolve the few budget-limited cells.
- Site selection: feed grazing-impact lat/lon + final great circle into the
  geocentric-radius (Chimborazo/Cayambe) intercept.
