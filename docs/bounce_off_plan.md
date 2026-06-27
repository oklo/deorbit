# Plan of record: the grazing-impact "bounce-off" phase diagram (Paper §2)

Status: PLAN (not yet executed). Owner: Greg + Claude. Codes: `hydrocodes/sph/` (SPH,
primary), `hydrocodes/euler/` (Eulerian cross-check), `src/deorbit/impact/grazing.py`
(reduced-order scout). This document is the plan a boundary-locating campaign (and a
future session) executes against.

## 1. The question and the gap

For an impactor striking a planet at angle θ from the horizontal, delineate the boundary
between three outcomes:

1. **conventional crater** — the projectile couples its energy into the target;
2. **scar + suborbital trajectory** — the projectile ploughs/ricochets and is launched
   downrange below escape speed (remaining intact, or dispersing into a fragment train);
3. **escape** — projectile material leaves the planet (v_out > v_esc), e.g. the Orcus-Patera
   "impact and leave" picture.

The literature has the *pieces* but not the *map*. Crater-morphology studies (Elbeshausen,
Wünnemann & Collins 2009/2013; Davison & Collins 2022; Suzuki et al. 2021) give circular→
elliptical→scar thresholds; projectile-fate studies (Gault & Wedekind 1978; Schultz & Gault
1990; **Pierazzo & Melosh 2000**, the closest ancestor; Daly & Schultz 2018) give ricochet/
escape behavior — but no systematic 3-regime phase diagram exists across size/type/planet/
velocity/angle, and especially **not for strong km-scale iron at ultra-grazing angles**.
That gap is the §2 deliverable.

## 2. The dimensionless space

The raw 5 parameters collapse onto impact-scaling π-groups (Holsapple 1993; Housen &
Holsapple 2007). The regimes are organized by a PRIMARY 2-D plane plus secondary cuts.

Raw → π map:

| raw | enters as |
|---|---|
| impact angle θ | **θ** (primary; ricochet onset) |
| velocity v, planet mass/radius | **v/v_esc** (v_esc=√(2 g R_pl)); also π₂ and the survival group |
| impactor radius R_p | gravity-scaled size **π₂ = g R_p / v²**; size-scaling of fragmentation |
| impactor type (Fe/stony) | density/impedance ratio **ρ_p/ρ_t**; projectile melt/strength threshold |
| target | strength **π₃ = Y_t/(ρ_t v²)**, ρ_t |

**Primary phase diagram: (θ, v/v_esc)** for a fiducial *strong km iron on rock*. Two boundaries:
- **crater ↔ scar/ricochet**: near θ_c ≈ 30° (experiments), modulated by v and size;
- **suborbital ↔ escape**: where retained downrange speed `ε(θ)·v·cosθ > v_esc`.

**Secondary cuts** (sensitivity + generalization):
- projectile-survival group `Π_s = ρ_t (v sinθ)² / σ_melt,p` — intact vs. fragment/melt, set by
  the NORMAL velocity component (the grazing decoupling: v_n=v sinθ craters, v_t=v cosθ carries
  the body downrange);
- gravity-scaled size **π₂** (hold-down on ejecta + crater regime);
- density/impedance ratio **ρ_p/ρ_t** (iron vs. stony).

**Anchor cases (the framing test — both must land sensibly):**
- **Cayambe** (Earth; v≈7.4 km/s, θ≈2°): v/v_esc = 7.4/11.2 = **0.66 < 1** → escape impossible;
  v_n = v sinθ ≈ 260 m/s → peak normal pressure ρ_t v_n² ≈ 180 MPa ≪ iron melt → survives →
  **regime 2 (suborbital ricochet → reimpact)**. The narrative case.
- **Orcus** (Mars; v≈10 km/s, v_esc≈5 → v/v_esc≈2) but a basin formed → steeper θ, **regime 1/
  escape-adjacent**. (See [[deorbit-orcus-patera]].)

Falsifiable prediction the diagram encodes: the **escape regime opens only at v/v_esc ≳ 1**
(interplanetary-speed impactors), so a sub-escape body like Cayambe *cannot* leave — only
ricochet suborbitally.

## 3. Per-run regime classifier

From each SPH run (most already computed for the Cayambe run):
- largest connected projectile-material clump → **retained mass fraction**;
- its **bulk v_out + flight-path angle** and velocity **dispersion** (intact vs. dispersing);
- **melt/vapor fraction** (peak-pressure based; cross-check with Euler);
- crater/scar geometry (depth, length, ellipticity — precise shape deferred to Euler).

Classification:
- **(1) crater** — retained fraction low / v_out small;
- **(2) suborbital** — retained, v_out < v_esc (sub-split intact vs. dispersed via the measured
  dispersion, then handed to analytic ballistic continuation for the reimpact footprint);
- **(3) escape** — v_out > v_esc.

## 4. Boundary-locating campaign (NOT grid-filling)

1. **Scout** the π-space with `src/deorbit/impact/grazing.py` (reduced-order) → approximate
   regime map. Cheap, analytic, sets the bisection brackets.
2. **Bisect to the boundaries with SPH**: along ~5 lines of constant v/v_esc, bisect on θ to
   find each transition to a few % (~6 runs/boundary). ≈ 5 × 2 × 6 ≈ **60 SPH runs** for the
   primary diagram.
3. **Resolution convergence** (CPPR 6→8→10) at ~3 boundary points — non-negotiable for the
   thresholds (the π-scaling gate showed laptop CPPR≈4 is penetration-dominated/under-resolved).
4. **Two-code cross-check** at ~4 overlapping points: Euler (with the material tracer, §6) vs SPH.
5. **Validate** (credibility anchors, do FIRST): reproduce the Gault & Wedekind ~30° ricochet
   onset + retained-velocity fraction, and a Pierazzo & Melosh 2000 point (10 km dunite, 20
   km/s, 30°/15° → escaping fraction).

**Division of labor (the two-code strategy):**
- **SPH owns** the fate question — coherence, fragmentation, retained mass/velocity, the grazing
  ricochet regime (Lagrangian, tracks projectile material natively).
- **Euler owns** final crater morphology + peak-pressure/melt, and provides a *diffusive*
  material-partition cross-check via the tracer scalar (§6).
- **Analytic ballistic** continuation determines escape/reimpact once fragments clear the
  violent phase (don't carry them to escape on the grid).

## 5. Compute and execution model

Thresholds need CPPR≳8–10 → millions of particles → the SPH **GPU** port (~35–70 min per
grazing traverse at ~5.5M). ≈60 runs + convergence ≈ **~50 GPU-hours** — a multi-day
laptop-GPU campaign. This maps onto the project's `engine.py` anytime-search pattern: a
**boundary-locating daemon** that runs bounded SPH increments, bisects each transition, and
checkpoints the phase-diagram state (regime per (θ, v/v_esc) node + the bracket intervals).
The π-reduction + bisection is what makes this feasible instead of astronomical (a raw 5-D
grid would be ~10⁴–10⁵ runs).

## 6. Prerequisite capability: the Euler material tracer (a milestone)

To let Euler answer the *fate* question (genuine two-code cross-check, and to color the
projectile in the Orcus 3-way figure), the Euler code needs a passive material tracer.
Scoped as a standalone milestone — see `hydrocodes/euler/README.md` (M-tag) and §"Euler
tracer milestone" below. Honest caveat: Eulerian advection is diffusive, so the tracer gives
the material **partition + bulk velocity** well but **not** fine coherence/fragmentation —
which is precisely why SPH remains primary for fate.

## Campaign status

- **Step 1 (reduced-order scout) — DONE.** `python3 -m deorbit.impact.scout`
  (`src/deorbit/impact/scout.py`, pure stdlib; result `results/bounce_off_scout.json`)
  classifies the (θ, v/v_esc) plane into **C** crater / **I** intact ricochet / **D** dispersed
  ricochet / **E** escape, sub-splitting ricochet by the contact-pressure-vs-iron-melt
  criterion. Topology confirmed and both anchors land correctly: **Cayambe → I**, **Orcus → C**.
  Key boundaries (Earth fiducial, to be refined by SPH): crater onset θ≈15° (search [8,22]°);
  intact→dispersed at v/v_esc≈0.78 (iron-melt threshold, v·cosθ≈8.7 km/s); dispersed→escape at
  v/v_esc≈1.44 (θ=2°), rising with θ. The scout emits a widened SPH search bracket per transition.
- **Step 2 (SPH bisection daemon) — BUILT + RUNNING.** `run_bounce.py`
  (`src/deorbit/impact/boundary_daemon.py`, pure stdlib) bisects each scout bracket with real
  GPU-SPH grazing runs: a generic iron-on-basalt grazing IC → `gpu_ic` (CPPR=8) → a
  projectile-clump regime classifier (bulk v_out vs v_esc; dispersion for intact-vs-dispersed).
  Anytime/checkpointed (`state/bounce/`), `--daemon`/`--increment`, registered for 24/7
  launchd supervision (`deorbit-bounce`) + on the dashboard. Validated end-to-end: the driver
  (bisection/checkpoint/resume) and the SPH score on a Cayambe-like point (θ=2°, v/v_esc=0.66)
  → **I**, v_bulk 7.19 km/s, dispersion 0.78 km/s, projectile climbing off the surface — a
  textbook intact ricochet, matching the scout + the narrative. (Classifier thresholds for the
  D/E boundaries are seeded by this anchor and refine as the campaign's runs accumulate.)
- **Step 2 — CONVERGED (69 SPH runs, 13/13 boundaries; `results/bounce_off_phase_diagram.json`,
  regenerate with `run_bounce.py --report`).** The (θ, v/v_esc) phase diagram for km iron on basalt:
  - **crater onset** (ricochet/dispersed → retained): θ_crit = **36.2 ± 0.6°** @ v/v_esc=1.0 —
    consistent with Gault & Wedekind (1978) ~30° grazing ricochet onset;
  - **intact → dispersed** (shock-melt onset): v/v_esc ≈ 0.96 at θ=2–5°, dropping to ≈0.60 at θ=12°
    (steeper ⇒ more normal-component shock ⇒ melts/disperses at lower speed);
  - **dispersed → escape**: v/v_esc ≈ 1.09 at θ=2° rising to ≈1.43 at θ=5° (escape harder at
    steeper angles); at v/v_esc=1.5–2.0 escape only below θ≈4–8°.
  Both anchors land correctly: **Cayambe → I**, the **Orcus-like steep/fast corner → C/E**.
  (Daemon continues: 4 crater brackets added at v/v_esc=0.5/0.75/1.5/2.0 to confirm the onset's
  velocity dependence, since it currently rests on the v/v_esc=1.0 line.)
- **Steps 3–5 — NEXT:** CPPR 8→10 convergence at 2–3 boundary points (confirm thresholds are
  resolution-independent), SPH↔Euler+tracer cross-check, experiment + Pierazzo&Melosh validation.

## 7. Deliverable + novelty

A dimensionless **(θ, v/v_esc) phase diagram** — crater / suborbital-ricochet / escape — for
strong km-scale iron at ultra-grazing angles, with retained-mass and v_out contours, secondary
cuts in Π_s / ρ_p/ρ_t / π₂, validated against the ricochet experiments + Pierazzo & Melosh +
the SPH↔Euler cross-check, with Cayambe and Orcus as anchor points. No such diagram exists in
the literature for this regime → a genuinely novel, referee-defensible §2.
