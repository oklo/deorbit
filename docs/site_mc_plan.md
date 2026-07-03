# Terminal-impact-site Monte Carlo — design (2026-07-03)

**Question:** For the km-iron deorbit scenario, what is the probability
distribution of the terminal impact site — and do high-drama targets
(Chimborazo, Cayambe, Huascarán, Kilimanjaro, Mt Kenya) get a target
cross-section significantly ABOVE their geometric area?

**Pre-registered hypotheses.** Greg: summit-class targets provide no
significant enhancement (η ≈ 1, impacts land ~area-weighted within the
dynamically selected latitude band). Counter (the steering argument): the
equatorial bulge + slow grazing descent give high equatorial terrain a large
enhancement, at least at the mountain-RANGE scale (the Andes as a meridional
wall), possibly less at the summit-cone scale.

**Headline statistic:** η(elevation class | lat band) =
P_impact(elev) / P_area(elev), with bootstrap CIs, plus a
P(ocean / lowland / >3 km / >5 km / named-summit) table and impact lat/lon
density maps. Chaos handled at the distribution level, never per-trajectory.

## Architecture: one 3D pipeline, priors decoupled by reweighting

Each Monte Carlo draw flies ONE continuous 3D trajectory from arrival
asymptote to ground: capture pass -> aerobraking cascade -> terminal plunge ->
DEM intercept. No seams. The astronomical prior is NEVER baked into the
dynamics: every sampled quantity (v_inf, B-plane point, arrival direction,
epoch/GMST, Cd, f_rho, f_H) is stored per row in the impact catalog, so priors
are applied/upgraded by post-hoc reweighting at zero CPU cost.

### Components

1. **Physics upgrades** (the load-bearing pair):
   - **Oblate atmosphere**: density evaluated at GEODETIC altitude above the
     WGS84 ellipsoid, not r − R_sphere. The 21.4 km equator-pole radius
     difference is ~3 scale heights = ~20x density contrast at fixed
     geocentric radius — first-order for WHERE an inclined orbit dies.
     Must use the same flattening as the J2 term (consistency gate).
   - **3D co-rotating wind**: v_atm = omega x r (the 2D code has the planar
     version).
   - Deliberate exclusions (all justified, state in paper): ablation
     (dm/m <= 3e-4/pass), thermospheric variability (periapses at 2–45 km),
     winds (~30 m/s vs 8 km/s), lower-atmosphere latitude structure (~±10%,
     folded into the f_rho ensemble knob), solar radiation pressure.

2. **Integrator** (`src/deorbit/sitemc/`, pure stdlib like the corridor
   layer — daemon/cloud-ready with zero setup): geocentric 3D, point-mass +
   J2 + oblate-atmosphere drag + analytic Moon AND Sun (Meeus truncated
   series; arcminute accuracy is orders beyond need). RK4 with the corridor's
   adaptive-dt pattern (fine in atmosphere, coarse on coast arcs).
   The planar corridor code is a validation limit; the REBOUND+ASSIST layer
   is the oracle (its drag hook is ours — oblate density drops in for
   cross-checks).

3. **Terrain**: global DEM at 1 arc-min (ETOPO int16 flat binary; band
   |lat|<=35° only, ~180 MB; data/fetch.sh + checksum per repo convention;
   if only netCDF4/HDF5 distributions are available, a one-time venv
   conversion tool writes the flat binary the stdlib code mmaps).
   Catch test: bilinear terrain geocentric radius along track wherever
   geodetic altitude < ~10 km; impact when r < r_terrain. Known systematic:
   1' grids shave sharp summits — patch in the curated summit table
   (site_selection.py PEAKS) as local overrides.

4. **Priors** (first pass, upgradeable by reweighting): v_inf uniform-ish over
   the capturable tail [0.1, 5] km/s (report conditionals); arrival radiant in
   a near-ecliptic cone; B-plane point uniform over the corridor annulus;
   Cd/f_rho/f_H from the existing ensemble priors; epoch/GMST uniform.
   Phase E upgrade: Opik/Granvik-weighted arrival distribution -> reweight the
   stored catalog. Free byproduct: the capture-INCLINATION distribution (the
   open question from the Cayambe discussion).

5. **Statistics & diagnostics**: besides eta and the P-table — terminal-phase
   "discrimination race" numbers (window pass count, low-arc lengths,
   death-longitude vs uniform), conditionals on inclination class, and the
   isolating A/B: spherical vs oblate atmosphere on identical seeds.
   Cheap sample multiplier: re-perturb each cascade's terminal phase
   (~10x effective draws; the endgame is chaotic, the cascade is not).

### Validation gates (pytest, house style)

- Drag-off: Kepler closure + analytic J2 secular rates (dOmega/dt, domega/dt)
  in 3D.
- Spherical-equatorial limit reproduces the 2D corridor capture boundary
  (the sqrt(2)-closure test carries over).
- Oblate-off flag reproduces the spherical integrator bit-for-bit intent.
- Same-IC vs the hi-fi oracle: pass-by-pass apogee ladder + node regression
  until chaotic divergence; then ensemble-level agreement on a seed set.
- DEM pipeline reproduces the site_selection.py summit ranking (Chimborazo
  global max 6384.4 km, etc.).
- dt-convergence of impact DISTRIBUTIONS (not trajectories) on a fixed
  seed set.

### Compute budget (laptop economics)

~30–60 s/trajectory pure Python (~2e5 in-atmosphere RK4 steps + coast).
10^4 trajectories ≈ half a day on 14 cores; 10^5 ≈ 3–5 days — resolves
P ~ 1e-3 directly; the eta curve needs far less. Idempotent checkpointed
sweep (foreground QoS + caffeinate — the euler lessons), registered on the
board. Scheduling: euler AF ladder owns the cores through ~2026-07-05;
Phases A–B cost ~no CPU and proceed immediately; the MC grind slots behind
the AF sweep (or 7/7 core split if both should move).

### Phases

- **A — build + gates** (1–2 sessions, no CPU): ellipsoid/oblate-atmosphere
  module, 3D integrator, Moon/Sun series, DEM fetch + terrain module.
- **B — cheap physics A/Bs** (hours): oblate-vs-spherical reference cascade
  (how strongly does drag death concentrate equatorward?); terminal
  window-pass statistics. May largely settle the summit question early.
- **C — MC grind** (background, 1–2 weeks wall).
- **D — analysis + figures.**
- **E — arrival-prior upgrade** (reweight only, zero re-simulation).

### Context (what led here — see docs/capture_feasibility.md §6–8)

Established: equatorward bias (bulge 21.4 km >> topo 9 km) is robust; drag
does not damp inclination; the final plunge is sub-orbital. NOT established:
that the discrimination window (summit radius -> sea level, ~6 km) is crossed
slowly enough for band-max steering — quick estimate says only ~2–5 window
passes (per-orbit loss grows as sqrt(r_p/e) as the orbit rounds), each a
~±2000 km arc at ~30° longitude spacing, then a final descent whose last 6 km
of altitude covers only ~90–350 km of track (fpa 1–4°). Also open: arrival
inclinations for low-v_inf encounters are typically ~10–25° (ecliptic tilt),
not the i≈1° Cayambe corner. The current spherical-atmosphere drag is the
identified modeling gap; the Andes-as-a-wall downstream of the Pacific
low-drag stretch is the identified structural bias in favor of Andean
impacts among land outcomes.
