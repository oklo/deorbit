# PICKUP — Euler Phase 3 (AF calibration to the depth–diameter curve)

Hand-off for a fresh context. Phase 2b is COMPLETE (cylindrical elastic-plastic strength + the full AF
model — block activation, shock seeding, fluidized viscosity, vib transport — all gated, GPU==CPU; see
PICKUP-phase2b.md). Phase 3 = calibrate the AF parameters (TDEC, ETA_AF) against the **Wünnemann & Ivanov
2003** depth–diameter curve (Planet. Space Sci. 51:831), then apply the fixed params to the 3D oblique
Chicxulub run (Collins et al. 2020).

## What's DONE (committed on master)
- **`crater` driver (CPU)** — commit eae390a. Axisymmetric vertical impact: a basalt impactor sphere into a
  basalt half-space under Mars gravity, composing route-1 WB substrate + vacuum-aware free surface + strength
  + Grady-Kipp damage + the full AF model, all in axisymmetric r-z. Args:
  `crater [a_m] [U_ms] [TDEC] [ETA_AF] [g] [cppr] [tend(<=0=auto)] [Rfac] [Zfac] [profile_path]`.
  - Runs STABLY across sizes; stays finite (max|v| settles to ~50 m/s post-impact). Transient depth matches
    gravity-regime π-scaling order-of-magnitude (a=500m → ~2.9 km transient).
  - Dumps the final surface profile z_s(r)−z0 to a file for ROBUST OFFLINE depth-diameter measurement.
    The in-loop `RESULT a D_app depth transient d/D` line is PRELIMINARY (the apparent-diameter datum-crossing
    is fragile for broadly-disturbed large craters — see open problems).
- **Sweep harness `sweep_crater.py`** — runs `crater` over sizes × dimensionless AF params (Tfrac, Efrac),
  dumping every profile to state/profiles/ and a summary row to state/results.csv. Idempotent/resumable
  (skips rows already done). `--increment N` (cloud), `--daemon` (local grind), `--list`.
  AF scaling (impactor-scaled, W&I-style, size-INDEPENDENT so one pair calibrates the whole curve):
  `TDEC = Tfrac·a/c_s`, `ETA_AF = Efrac·ρ·c_s·a` (c_s≈3000, ρ=2700). Grid: SIZES 300..3000 m × {off + 8 AF pairs}.

## UPDATE — friction/ROCK yield is now IN (commit f0ac4aa); fixes collapse-to-flat, but tuning remains
The pressure-dependent ROCK yield (intact Lundborg + cohesionless damaged friction, AF-coupled) is implemented
+ gated (`friction` gate, machine-precision; ROCK flag, off by default so all prior gates unchanged). It FIXES
the collapse-to-flat blocker below: a=300 m AF-off went from 0 (flat, cohesion-only) to a HELD crater with ROCK.
BUT it is not yet a settled simple crater: ROCK a=300 m AF-off creeps 0.45 km @76 s → 0.25 km @150 s
(max|v| 24.6→14.3 m/s, still slowing) and is OVER-RELAXED (d/D≈0.04 vs the ~0.2 of a real simple crater).
Root: pure cohesionless damaged friction Y_d=μ_d·P is marginal at the low confining pressure of a small/shallow
crater. LIKELY FIX (next): add a small DAMAGED COHESION Y_d0 (iSALE keeps residual breccia cohesion: Y_d=Y_d0+μ_d·P),
plus the AF combination + resolution convergence, then measure the SETTLED state. These are calibration knobs →
hand to the sweep daemon (free compute), not interactive 15–30 min runs.

## (historical) the blocker this replaced — missing friction (Drucker-Prager) strength model
The driver excavates a PHYSICAL transient crater (a=300 m, cppr=6: d/D=0.22 at t=15 s — a textbook simple
bowl). But the crater then **slowly flows back to flat** under gravity (a=300 m AF-off: depth 1.40 km @15 s
→ 1.10 @30 s → 0.15 @50 s → 0.00 @76 s). It is NOT ringing — it monotonically relaxes.
ROOT CAUSE: the strength model is cohesion-only von Mises with `Y=(1−D)(1−af)·Y0`. The impact damages the
rock (D→1), so Y→0, and **damaged rock becomes a strengthless liquid that relaxes flat under gravity**. AF
only weakens strength further, so it cannot hold the crater either (it accelerates collapse). High AF
viscosity merely *slows* the collapse (a=1500 m, ETA=1e9 freezes it mid-collapse at 3.4 km — not a settled
crater).
This is exactly why production cratering codes use a **pressure-dependent (friction) yield with a damaged
branch** (iSALE "ROCK" model, Collins, Melosh & Ivanov 2004; Lundborg yield): damaged breccia retains
FRICTION strength `Y_d ≈ μ_d·P` (cohesionless but pressure-dependent), which holds a simple bowl and sets
the realistic d/D. WITHOUT it, no simple crater is held and the AF depth–diameter calibration is meaningless.

## NEXT STEP (the unblock) — add a ROCK-style friction yield, THEN calibrate AF
Replace the cohesion-only cap in `vonmises` with a pressure- and damage-dependent yield (own gate first):
  - Intact:  `Y_i(P) = Y0 + μ_i·P/(1 + μ_i·P/(Y_m − Y0))`   (Lundborg; Y0=cohesion, μ_i≈1.2, Y_m≈2.5 GPa)
  - Damaged: `Y_d(P) = min(μ_d·P, Y_m)`                      (cohesionless friction, μ_d≈0.6)
  - `Y = [(1−D)·Y_i + D·Y_d]·(1−af)`   (AF fluidizes the friction-supported state — the W&I picture)
  Needs new Material params (μ_i, μ_d, Y_m, and a sensible cohesion Y0 — the current basalt Y=350 MPa is too
  high for cohesion; iSALE basalt ~ Y0~1e7, μ_i~1.2, μ_d~0.6, Y_m~2.5e9). GATE: a confined-shear test —
  the yield surface must follow Y_i(P)/Y_d(P) vs pressure (an analytic oracle), GPU==CPU. THEN re-run `crater`:
  a small crater must HOLD a settled bowl (d/D~0.2) at long tend; only then does the AF sweep mean anything.

## OPEN PROBLEMS / the Phase 3 grind (after the friction unblock)
1. **Friction/ROCK yield — DONE (commit f0ac4aa), gated.** Remaining tuning: add a small damaged cohesion
   Y_d0 (the pure-friction crater still creeps/over-relaxes, d/D≈0.04 — see UPDATE at top); pick Y0/μ_i/μ_d/Y_m/Y_d0
   from the literature; combine with AF; converge.
2. **Measure the SETTLED crater** (run to settling / max|v|→0), not a fixed tend — the current in-loop d/D and
   the datum-crossing D_app are unreliable for the still-relaxing surface. Do it OFFLINE from the dumped
   state/profiles/* (rim crest + floor + apparent diameter, robust edge logic).
3. **Resolution convergence** (cppr 8, 12 on a couple of sizes) before trusting absolute d/D.
4. **Digitize the W&I 2003 depth–diameter curve** (and Mars d–D data) as the oracle; tune (Tfrac, Efrac) to
   match, esp. the ~7 km Mars simple→complex transition. The sweep harness `sweep_crater.py` is ready to grind
   this once the friction model makes craters physical (do NOT register the 24/7 daemon before then — it would
   only grind collapse-to-flat).
5. **GPU port of `crater`** for the production sweep / 3D.
6. Then → 3D oblique Chicxulub with the FIXED calibrated AF params (60°, 17 km granite, 12 km/s; Tillotson
   basalt crust proxy — ANEOS granite/dunite still a flagged EOS gap).

## Run / grind
- Build: `clang++ -std=c++17 -O2 hydro_cpu.cpp -o hydro_cpu`
- One crater: `./hydro_cpu crater 500 12000 30 1e8 3.71 6 -1 18 22 /tmp/p.txt`
- Sweep: `python3 sweep_crater.py --list | --increment N | --daemon`. Results: state/results.csv; profiles: state/profiles/.
- The local supervisor (~/investigations/daemons.json) keeps the `--daemon` grind alive 24/7 (entry `deorbit-crater`,
  cwd hydrocodes/euler, log state/sweep.out). **LIVE as of 2026-06-28** — grinding the first 42-run pass
  (Y_d0 in {1,5,10} MPa x AF off/on x 7 sizes, ROCK friction on). Monitor: `cat state/results.csv`,
  `tail state/sweep.out`, `python3 sweep_crater.py --list`. Each run ~15-30 min (cppr=6), so a full pass ~12-20 h.
  After the pass: analyse profiles OFFLINE for settled d/D-vs-D, pick the (Y_d0, AF) that holds simple craters at
  d/D~0.2 and gives the complex transition, vs Wuennemann & Ivanov 2003.
