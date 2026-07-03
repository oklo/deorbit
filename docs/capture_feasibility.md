# Capture feasibility — the calculation chain, revisited (2026-07-03)

**Claim being established:** deorbital capture of a km-scale *monolithic* M-type
(metallic) asteroid by Earth's atmosphere is physically possible with current
solar-system and terrestrial parameters — no fine-tuning beyond a low-probability
(but observed-to-exist) encounter geometry.

Every number below is re-derived analytically and pinned to the integrator
(`src/deorbit/corridor/physics.py`) by `tests/test_capture_feasibility.py`
(energy loss to ~5%, capture boundary to ~8%); the corridor map + uncertainty
ensemble are `results/corridor_latest.json` / `results/corridor_ensemble.json`.

## 1. The object

1 km diameter iron sphere (ρ = 7800 kg/m³): m = 4.08×10¹² kg, A = 7.85×10⁵ m²,

    beta = m/(Cd·A) = 5.2×10⁶ kg/m²   (Cd = 1)

This enormous ballistic coefficient is the whole game: the body barely notices
the atmosphere on any single pass, which is exactly why a *narrow* corridor of
grazing passes can shed just enough energy to capture without destroying it.

## 2. Why aerocapture is the only channel

- **Pure gravitational capture is impossible** — two-body energy is conserved;
  a hyperbolic arrival leaves hyperbolic.
- **Binary-exchange capture** (Agnor & Hamilton 2006, Triton) fails by 3–4
  orders of magnitude at km scale: the binary internal velocity
  v_bin ~ √(Gm/a_bin) is sub-m/s against km/s encounter excesses, and even a
  successful exchange parks the body in a high, non-decaying orbit.
- **Tidal-disruption capture** (the Tomkins et al. 2024 Ordovician-ring
  mechanism) requires a rubble pile — see §5: at grazing periapsis a rubble
  pile disrupts, a monolith does not. The monolith channel and the ring channel
  are *mutually exclusive outcomes of the same flyby geometry* (a bifurcation
  on internal strength).
- **Aerocapture** is what remains, and it works — quantitatively, below.

## 3. The grazing-pass drag budget (closed form)

Atmosphere: piecewise-exponential (Vallado), density ρ(h), scale height H(h).
A hyperbolic pass with vacuum periapsis altitude h_p has periapsis speed
v_p = √(v_inf² + 2μ/r_p) ≈ 11.2 km/s (escape speed) for v_inf ≪ v_esc.

**Traversed column.** Near periapsis the altitude rises along-track as
h(s) = h_p + (s²/2)(1/r_p − 1/p) where the trajectory's radius of curvature at
periapsis is the semi-latus rectum p. For a near-parabolic graze p ≈ 2 r_p, so
h(s) ≈ h_p + s²/(4 r_p) — the path hugs the shell **√2 longer** than the naive
surface-curvature estimate:

    Sigma_eff = rho(h_p) · sqrt(2π · r_eff · H),   r_eff = r_p·p/(p−r_p) ≈ 2 r_p

(The integrator confirms this factor: integrated per-pass energy loss sits at
1.38–1.42× the naive flat-column estimate across h_p = 5–40 km, i.e. within a
few % of √2 = 1.414. Sanity limits: e→0 gives r_eff→∞, a circular orbit never
leaves the shell; e≫1 gives r_eff→r_p, the straight-line flyby.)

**Per-pass loss and capture condition.** With F = ½ρ C_d A v²,

    dv_p    = v_p · Sigma_eff / (2·beta)
    d_eps   = v_p · dv_p                       (specific orbital energy)
    capture ⇔ d_eps > v_inf²/2
    ⇒  v_inf,max(h_p) = v_p · sqrt(Sigma_eff / beta)

**Numbers** (analytic | integrator-bisected boundary, km/s):

| h_p (km) | Σ_eff (kg/m²) | v_inf,max analytic | v_inf,max numeric |
|---|---|---|---|
| 30 | 1.3×10⁴ | 0.56 | 0.56 |
| 20 | 5.9×10⁴ | 1.19 | 1.18 |
| 10 | 2.4×10⁵ | 2.38 | 2.40 |
| 5  | 4.7×10⁵ | 3.35 | 3.42 |
| 2  | 7.1×10⁵ | 4.13 | 4.26 |

(The residual few-% excess at the deepest rows is the drag-lowered actual
periapsis dipping below the vacuum h_p.) The full sweep (60² grid, current
`results/status_latest.txt`) gives **capture out to v_inf = 4.19 km/s** with
periapsis windows up to **~43 km wide** at low v_inf — a corridor, not a razor
edge. The 160-sample uncertainty ensemble puts vmax_capture at 4.0 nominal,
p05–p95 = 3.2–5.6 km/s; the variance is dominated by C_d (~56%) and
lower-atmosphere density (~40%). Because capture periapses are at **2–45 km
altitude** (stratosphere/mesosphere), thermospheric solar-cycle variability is
irrelevant — the corridor is set by the body's drag coefficient and the
well-constrained lower atmosphere.

## 4. Example capture (integrator)

v_inf = 3 km/s, h_p = 5 km: **capture**, Δv = 521 m/s in one 47 s graze
(min altitude 4.9 km), exiting on a bound orbit with apogee 292,600 km.
(This particular apogee exceeds the Moon filter — see §6; gentle-chain
captures need v_inf ≲ 3.2 km/s *and* a deep pass.)

## 5. Survival at periapsis — where "monolithic M-type" is load-bearing

| stress | magnitude | monolithic iron | rubble pile |
|---|---|---|---|
| peak ram pressure ½ρv² | ~40 MPa (h_p = 5 km) | survives (strength ~100–500 MPa, Canyon Diablo) | shreds (strength ~kPa–MPa) |
| tidal (rigid) | ~6 kPa across 1 km | negligible | comparable to self-gravity central pressure (~2 kPa) → stripped |
| fluid Roche limit | 2.17 R_E for ρ = 7800 | irrelevant (strength-dominated) | grazing r_p = 1.001 R_E is **inside** it → disruption |
| ablation (Sutton–Graves, all-heat-to-melt bound) | ≤ ~5 cm skin per pass | Δm/m ≤ 3×10⁻⁴, negligible (and mass loss *lowers* β, aiding capture) | — |

So the same grazing encounter bifurcates on internal strength: **rubble pile →
tidal disruption → ring** (the Ordovician/Tomkins channel); **monolith →
survives intact → aerocapture cascade** (this paper's channel).

## 6. After the first pass: the Moon filter and the aerobraking cascade

(High-fidelity REBOUND+ASSIST layer, `docs/highfidelity_findings.md`.)
Capture alone is not the end state — the first apogee decides:

- first apogee ≲ 130,000 km (⅓ lunar distance): **gentle** — the orbit
  aerobrakes pass-by-pass to an ultra-oblique grazing terminal impact
  (nominal case: 66 passes / 26 days → 7.39 km/s at flight-path angle −1.16°).
  Requires v_inf ≲ 3.2 km/s (nominal; ensemble p05–p95 = 2.2–5.0).
- ~140,000–600,000 km: the Moon pumps periapsis → steep impact.
- ≳ 2×10⁶ km: ejected. (Marginal shallow captures die here — quantified in the
  ensemble's "bound" filter.)

Drag does not damp inclination; J2 precesses the node/argument, so the
terminal grazing impact samples the |lat| ≤ i band and lands at the maximum
*geocentric radius* in that band (equatorial bulge ≫ topography) —
the site-selection result (Chimborazo/Cayambe for near-equatorial captures).

## 7. Does the current solar system supply such encounters?

Two observational existence proofs, no modeling required:

- **km-scale metallic NEAs exist now:** (6178) 1986 DA, radar-confirmed
  metal-rich, D ≈ 2.4–2.9 km (+ candidate 2016 ED85) — Sanchez et al. 2021,
  PSJ 2, 205 (doi:10.3847/PSJ/ac235f).
- **near-zero v_inf encounter geometries occur:** the observed temporarily
  captured minimoons 2006 RH120 and 2020 CD3 (meter-scale, captured with *no*
  atmospheric assist) demonstrate the v_inf ≈ 0 tail of the encounter
  distribution is populated (cf. Granvik, Vaubaillon & Jedicke 2012).

Feasibility therefore does not hinge on a rate estimate: the corridor exists
(§3), the survivor class exists (§5, §7a), and the encounter geometry class
exists (§7b). **TODO (rate, for the paper's discussion):** fold the Granvik et
al. 2018 NEO model velocity distribution against the corridor cross-section
(v_inf ≤ 4 km/s, ~40 km periapsis window) × the M-type fraction for an
events-per-Gyr estimate.

## 8. Caveats / sensitivity (state these in the paper)

1. **C_d and shape/tumbling** dominate the corridor-width uncertainty (56% of
   ensemble variance); a factor-2 in C_d moves vmax by ~0.6 km/s.
2. Phase-1a is planar two-body + drag; the three-body (Moon) and J2 effects
   are layered on via the high-fidelity runs, which *shrink* the usable
   corridor (bound/gentle filters) — the feasibility claim uses the filtered
   numbers, not the raw capture map.
3. Current-epoch parameters. For the Ordovician application: the Moon was a
   few % closer (apogee filter scales with lunar distance), day shorter
   (co-rotation term, ~3% of variance), paleo-atmosphere density uncertain —
   all second-order against the C_d spread.
4. Single-pass energetics assume no fragmentation; §5 is the argument that a
   monolith stays whole (40 MPa ≪ strength). A cracked-but-coherent body
   (Campo del Cielo class) is the interesting marginal case.
