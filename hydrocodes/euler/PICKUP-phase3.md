> ## BINGHAM LADDER VERDICT (2026-07-05; 10/10 done, state/bingham/, harvested dD.csv)
> **Mid sizes: the FIRST genuine slump.** a=1000: d/D 0.146 (ratio 2.45-2.48 vs
> friction 2.85, flat-ladder 3.3); a=2000: 0.123-0.143 (ratio 3.1-3.3 vs 4.0/4.3+).
> Mechanism visible in traces: depth holds ~transient while D_app WIDENS ~70%
> (21->34 km at a=1000) -- outward slump, truncated when vib decays. Y_B 1->10 MPa
> is nearly IRRELEVANT (outcomes identical) => the arrest is still TDEC-limited.
> a=300 guard: ratio 1.48->1.66, D widened into the 'complex' band (mild cost).
> **a=3000 REINTERPRETED: there was never runaway excavation.** Traces (both the
> Bingham rungs AND yesterday's no-floor TDEC=1000/2000 runs) show depth PINNED at
> 15.90 km for ~4000 s, then an ABRUPT jump to 40.5-41.7 km at t~4000-4200 s.
> Signature of a late BREAKTHROUGH artifact, not physics: af_max stays ~0.97-1.0
> the whole run (immortal fully-fluidized pocket) while af_mean decays; Vexc grows
> monotonically (1.2e7->2.4e7, 'crater volume' doubling with no depth change);
> voidzero refills ~5.3e8 events (~1600 cells/step = a standing evacuated region
> churning below the floor). HYPOTHESIS: cavitation<->shock-gate feedback -- void
> cells cycling (reset to ambient, recompressed) spike dP > P_ACT each cycle,
> perpetually re-seeding vib -> an immortal fluidized subsurface pocket that
> churns until it connects to the surface datum (depth jumps to the void column).
> The short-TDEC flat ladder never saw it only because its windows ended at 1869 s.
> LIKELY FIX: refractory seeding gate (no vib re-seed in cells that were void
> within N steps, or require rho near reference at seed time). DISCUSS WITH GREG
> before surgery. Pre-breakthrough state (15.9 km held with floor active) is
> arguably the physical answer at a=3000 (ratio ~4.3, matching the a=2000 trend).

> ## AF TIMESCALE TEST VERDICT (2026-07-04 runs, harvested 2026-07-05): the middle
> regime DOES NOT EXIST in this AF implementation. TDEC=1000/2000 s at a=3000
> (Efrac 0.01) produced RUNAWAY EXCAVATION, not late slump: transient ~42 km deep
> (vs 15.9 at TDEC<=200), final depth 38 km, D_app ~110 km (crater FILLED the
> Rfac=30 domain -> numbers qualitative), d/D 0.34-0.36, ratio vs Garvin ~10x,
> still max|v|~80 m/s at cutoff. MECHANISM: while vib is active the (1-af)
> strength cut is ~total and Efrac=0.01 viscosity is far too low to arrest the
> flow -> the excavation cavity keeps growing quasi-hydrodynamically for ~TDEC;
> when vib decays, pressure-dependent friction (Y_d ~ 0.6*rho*g*z ~ 2e8 Pa at
> 38 km) locks ANY cavity instantly. Short TDEC freezes a too-deep transient;
> long TDEC digs a monster. The W&I-intended behavior (normal transient, then
> gentle viscosity-limited slump) requires the fluidized state to carry a
> SUBSTANTIAL effective resistance during flow. SURGERY OPTIONS (Greg to pick):
>   (A) cap af (Y_eff=(1-af)Y with af<=~0.8) — one-line, keeps 20% strength while
>       fluidized; (B) Bingham floor: Y_eff=max(Y_B,(1-af)Y), Y_B~1-5 MPa;
>   (C) raise ETA 10-30x with viscous substepping (closest to W&I block model,
>       biggest change). (A)/(B) are cheap and can be laddered same-day.
> Test artifacts: state/aftest/ (profiles, logs with af-state stream, dD.csv).

> ## AF LADDER VERDICT (2026-07-04 morning; 36/36 done, harvested to state/dD_af.csv)
> **The (Tfrac,Efrac) grid is FLAT and the late slump NEVER APPEARS.** Depths are
> essentially identical across the whole ladder (a300: 1.80-1.86 km; a1000: 5.40;
> a2000: 10.0-10.4; a3000: 14.4-15.0), depth ~= 94% of transient everywhere, and
> every rung is WORSE vs Garvin than friction-only (ratios 1.9-2.0 simple,
> 4.5-5.8 at a=3000, vs off_yd5M 1.5/4.0). Efrac has NO effect at any size.
> DIAGNOSIS: TDEC = Tfrac*a/c_s spans 6-200 s across the grid, but the slump
> timescale is ~2*pi*sqrt(D/g) ~ 700+ s for the big craters -- the fluidization
> decays before collapse begins at EVERY grid point; all AF can do in this
> parameterization is damp the rebound (deepening the transient, in proportion
> to TDEC). The impactor-acoustic scaling cannot reach the slump regime.
> PROPOSED NEXT (cheap, decisive; awaiting Greg): one pair at a=3000 with TDEC
> ~1000-2000 s (i.e., keyed to the crater gravitational timescale, not the
> impactor), Efrac 0.01 -> if slump appears, re-key TDEC to the transient-crater
> scale and re-ladder; if not, the AF block model itself needs revisiting.
> E22 cppr20: ~38 h wall, output file still empty, duty ~26% (decision today).

> ## RESUME HERE (2026-07-03 evening — vamb sweep HARVESTED; AF ladder RUNNING)
> **CLEAN (vamb=0.27) GARVIN VERDICT (42/42 settled, analyze_crater -> mars_dD_oracle):**
> the artifact-free baseline UNDER-COLLAPSES EVERYWHERE, and the deficit GROWS with size:
> simple branch off-arms ratio 1.46-1.62 (yd1M-yd10M); complex branch off_yd5M 1.9 (D=9.5)
> -> 4.0 (D=69); off_yd1M 2.0 -> 3.9; AF-on Tfrac=30 is WORSE everywhere (2.1 -> 5.8,
> the known rebound-damping over-deepening). The 2026-07-01 published-quality claim
> "complex branch falls too steeply" is FORMALLY DEAD (it was the voidzero refill).
> Y_d0 ordering is now SECONDARY (1-10 MPa arms differ ~20-30%); the missing physics is
> LATE-STAGE COLLAPSE -> the AF (Tfrac 60/100/200 x Efrac 0.003/0.01/0.03 @ yd5M,
> a=300/1000/2000/3000) ladder is THE calibration lever and is RUNNING (--jobs 7,
> launched 2026-07-03 ~17:50, nohup+caffeinate, log state/sweep_af.log; 36 jobs,
> ETA ~overnight; 7/7 core split with the sitemc MC per Greg). NOTE the required
> collapse enhancement grows with D -- watch whether the impactor-scaled AF
> (TDEC ~ a, ETA ~ a) reproduces that trend or needs stronger size-scaling.
> **WHEN AF DONE:** same harvest (analyze_crater.py state/profiles/af5M_T*_E*.txt
> --csv state/dD_af.csv; mars_dD_oracle.py) -> pick (Tfrac,Efrac) that flattens the
> ratio-vs-D trend toward 1 WITHOUT wrecking the a=300 simple guard.
> E22 cppr20: decision point 2026-07-04 morning (PICKUP-pierazzo.md).

# PICKUP — Euler Phase 3: AF + strength calibration to the depth–diameter curve

> ## RESUME HERE (2026-07-03 ~midday, master==origin/master @ 8360397; all code+docs committed)
> **TWO RUNS IN FLIGHT (both nohup+caffeinate detached; verify with `ps aux | grep hydro` FIRST):**
>   1. **Clean vamb=0.27 recalibration sweep** (Greg-approved): 42 jobs, 10 parallel CPU workers,
>      started ~10:00, ~1-2 h/job -> ETA late afternoon. Progress: `python3 sweep_crater.py --list`
>      + `tail state/sweep_vamb.log`. Rows -> state/results.csv (legacy sweep archived as
>      results_legacy_vamb0.csv / profiles_legacy_vamb0/ / dD_settled_legacy_vamb0.csv).
>      **WHEN DONE:** `python3 analyze_crater.py state/profiles/*.txt --csv state/dD_settled.csv`
>      then `python3 mars_dD_oracle.py state/dD_settled.csv` = the CLEAN Garvin comparison.
>      EARLY READ (7 rows in): friction-only craters now BARELY COLLAPSE (off_yd1M depth==transient,
>      12 km deep at a=3000) -> the complex branch needs REAL collapse physics = the AF (Tfrac,Efrac)
>      sweep is now the load-bearing calibration step. Watch: af_yd1M a750 came out DEEPER than
>      off (4.35 vs 3.30 km) = the known Tfrac=30 rebound-damping over-deepening; AF tuning likely
>      needs LONGER decay times (late slump, not rebound damping).
>   2. **E22 20 km/s cppr20 oblique rung** (state/pz45_rcfl100/, out_u20000_c20.txt empty = running):
>      healthy but slow (dt pathology partially survives rcfl=100 at this corner; ~21-26% CPU duty,
>      ~180k steps by 2026-07-03 10:00). Greg: DO NOT kill midday; decision point = 2026-07-04
>      morning — if still running, kill + either stand on cppr16 (E22 already PASSES both velocities,
>      slopes converged; see PICKUP-pierazzo.md verdict table) or spot-check rcfl=500 vs the
>      cppr12/16 slope trend.
> Context for both: the 2026-07-03 voidzero block below + PICKUP-pierazzo.md + memory
> [[deorbit-euler-code]]. Gates all green (39-40 PASS) @ 8360397. After the Garvin comparison:
> design the AF sweep grid (SETTINGS in sweep_crater.py; Tfrac ladder e.g. 60/100/200 x Efrac
> 0.003/0.01/0.03 at Y_d0 5MPa on 3-4 anchor sizes) and grind it the same way.

> ## SESSION UPDATE 2026-07-03 — VOIDZERO REFILL ARTIFACT CONFIRMED: the complex-branch
> ## over-collapse was largely NUMERICAL. RECALIBRATION REQUIRED before trusting the d-D curve.
> The 2026-07-02 prater work exposed that `void_cells` (CPU, hydro_cpu.cpp:230) / `voidzero` (GPU)
> reset evacuated (rho<100) cells to REF_R0 = the LITHOSTATIC reference — i.e. a below-surface void
> cell gets REFILLED with solid rock. Instrumented (VZ_N/VZ_DM counters) + added crater arg14
> `vamb` (VOID_AMB: >0 = reset voids to that AMBIENT density instead; default 0 = legacy).
> **A/B (off_yd5M config: U=12, g=3.71, cppr5, ROCK, Y_d0=5MPa, AF off; analyze_crater.py numbers,
> same analyzer as the sweep — legacy arms REPRODUCE the sweep exactly):**
> | run | d/D legacy | d/D ambient(0.27) | depth L->A (km) | legacy refills |
> |---|---|---|---|---|
> | a=300 | 0.186 (==sweep) | 0.222 | 1.08 -> 1.32 | 437 |
> | a=3000 | 0.018 (==sweep) | 0.165 | 0.80 -> 11.40 | 627 |
> MECHANISM at a=3000: transient depth ~12 km in BOTH arms; under legacy the floor comes back UP
> to 0.8 km — the 627 refill events inject rock of the same ORDER as the whole crater bowl mass
> (each event converts a near-empty 600 m cell at r~10-20 km to rho 2700; bowl ~1.6e16 kg) exactly
> at the floor/walls during collapse. Ambient arms also SETTLE (legacy arms hit the cap) — the
> refill was fighting the settling detector too. Simple branch mildly affected (+19% depth at a=300).
> **CONSEQUENCES:** (1) the published-quality claim "complex branch falls too steeply" was the
> artifact, not physics; (2) with ambient voids the complex branch now UNDER-collapses vs Garvin
> (11.4 km deep at D_app~69 km vs Garvin ~2.9 km) -> Y_d0/AF recalibration needed (expected — the
> old calibration was tuned ON TOP of the artifact); (3) the GPU `voidzero` needs the mirrored
> vamb option before any GPU crater runs (CPU-only switch so far).
> **NEXT (supersedes step 1 below):** re-run the sweep with vamb=0.27 (edit sweep_crater.py to pass
> arg14) -> re-do the Garvin comparison -> THEN the AF (Tfrac,Efrac) sweep on the clean base.

> ## RESUME HERE (as of 2026-07-01, master==origin/master @ 6bea00e)
> The AF crux is FIXED and the d-D pipeline WORKS end-to-end. We have the FIRST COMPLETE simple+complex
> d/D-vs-D curve (sweep 42/42, sponge fix) vs the Garvin-2003 Mars oracle. Read the two "SESSION UPDATE"
> blocks below (2026-06-30 and 2026-07-01) for the numbers.
> **State of the science:** the model reproduces the d/D DECREASE with size; Y_d0=5-10MPa BRACKETS the Garvin
> LEVEL; the COMPLEX branch falls TOO STEEPLY (big craters over-collapse) BUT the big-a craters are
> NOT-SETTLED at cap=6*t_auto, so that steep falloff is the least-trusted part.
> **IMMEDIATE NEXT STEPS (ordered):**
>   1. **Bump the big-a settling cap** (cheap, removes the biggest caveat): the 12km transients need >6*t_auto.
>      Raise the crater cap (hydro_cpu.cpp/hydro_gpu.cpp `6.0*t_auto`) for big a, or make it size-aware, and
>      re-measure the complex branch (a>=1500). Check whether the complex falloff flattens toward Garvin.
>   2. **AF (Tfrac,Efrac) sweep (task #8)** around Y_d0~7MPa: tune AF to lift the complex branch (AF's intended
>      role) and match the LEVEL. Edit sweep_crater.py SETTINGS; run `caffeinate -i -s python3 sweep_crater.py
>      --jobs 14` (FP64 CPU oracle -- GPU collapse is FP-chaotic for small a; big-a GPU==CPU to ~5%).
>   3. **cppr 8/12 convergence + repeat-run averaging** on 2-3 anchor sizes (single-run d/D ~20% chaotic scatter).
>   4. Then the 3D oblique Chicxulub run (the original goal; Collins et al. 2020).
> **Run/measure loop:** `./hydro_cpu crater a U TDEC ETA g cppr -1 30 22 prof.txt 1 Yd0` (settling; profile out)
> -> `python3 analyze_crater.py <profiles> --csv state/dD_settled.csv` -> `python3 mars_dD_oracle.py state/dD_settled.csv`.
> Gates: `./gates.sh quick` (34 PASS, GPU==CPU). state/ is gitignored (profiles+CSV are local).

## SCOPED TASK 2026-07-01: Pierazzo-2008 ALUMINUM peak-pressure parity (2D iSALE head-to-head)
Decision (Greg): the "2D iSALE replication" head-to-head benchmark = **Pierazzo et al. 2008** (MAPS 43:1917,
"Validation of numerical codes for impact and explosion cratering: strengthless and metal targets"), the
**aluminum-exact** case. **Collins 2020 is a 3D-OBLIQUE APPLICATION, NOT a 2D parity test** — keep it as the
endgame, don't conflate it with code validation.

FINDING (ground-truthed vs the paper + M5 gate): Pierazzo's ALUMINUM benchmark is **PEAK-SHOCK-PRESSURE vs
distance** (documented inter-code envelope ~10-20%), NOT crater growth — crater growth is the strengthless
DUNITE case. So the aluminum-exact head-to-head = the peak-pressure decay, and it is **~90% already done** (M5
gate, README): 1D Al-on-Al reproduces the analytic Tillotson Hugoniot to 0.0%; 3D Al-sphere isobaric core =
93% of Hugoniot. **OPEN CAVEAT** = the near-field decay exponent n~1.3 over r/a 1-6 is UNCONVERGED (canonical
far-field ~2); domain too small + cppr not converged.

WIN CONDITION: aluminum P_peak(r/a) along the impact axis lands **within the Pierazzo ~10-20% inter-code
envelope** across the documented r/a range, at converged resolution, with the fitted **far-field exponent -> ~2**.
That is a defensible "we reproduce iSALE's validated shock behavior" claim (peak P is the credibility-critical
quantity Pierazzo stresses). Zero EOS ambiguity (Al Tillotson is exact; `Material::aluminum()` in eos.hpp,
`AL` in hydro_gpu.cpp).

CONCRETE STEPS:
  1. **Get the paper's exact Al setup + published curve.** Al sphere into Al half-space, vertical. Extract the
     impactor diameter, velocity(ies), and DIGITIZE the cross-code peak-pressure-vs-distance figure (iSALE line
     + the inter-code envelope). Sources: MAPS 43:1917 (Wiley doi 10.1111/j.1945-5100.2008.tb00653.x; ADS
     2008M&PS...43.1917P; Arizona repo repository.arizona.edu/handle/10150/656500) and the MDPI-2021 follow-up
     "Benchmarking Numerical Methods..." (mdpi.com/2076-3417/11/6/2504) which re-plots these.
  2. **Far-field domain + cppr convergence** on the existing `pierazzo` mode (hydro_gpu.cpp: nx=ny=80,nz=100,
     dx=0.1, a=1m=10cppr, U=10km/s; diagnostic writes pierazzo_decay.txt + fits exponent n). Enlarge the domain
     so r/a reaches the far field (n~2 regime); run cppr 10/15/20 and show P_core + n(r/a) converge.
  3. **Compare** digitized iSALE/envelope vs our P(r/a): per-point ratio + fitted far-field exponent. PASS if
     inside the envelope with n->~2.
  4. **Gate it** (extend M5 / add a `pierazzo` convergence check to gates.sh so it stays green).

SCOPE LIMIT (state it in any writeup): this validates EOS + shock capture, NOT gravity/excavation/collapse.
The crater-FORMATION head-to-head needs either the DUNITE strengthless case (no ANEOS dunite here -> dense
Tillotson PROXY, an EOS mismatch that weakens a strict "replication" claim) or stays on the W&I-2003 d-D
calibration already underway (RESUME banner below). **Phase B (optional, only if Al pressure parity passes):**
dunite-proxy transient-crater comparison (transient, strengthless -> dodges the settling/collapse chaos).

---
Self-contained hand-off for a fresh context. **Repo:** github.com/oklo/deorbit, code in
`deorbit/hydrocodes/euler/` (CPU oracle `hydro_cpu.cpp`, GPU `hydro.metal` + `hydro_gpu.cpp`,
shared EOS `../common/eos.hpp`). master == origin/master as of commit 54ddc37.

## The goal
Calibrate the acoustic-fluidization (AF) **and** rock-strength parameters so the Eulerian code reproduces
the **Wünnemann & Ivanov 2003** crater depth–diameter dependence (Planet. Space Sci. 51:831), then apply the
*fixed* parameters to the 3D oblique **Chicxulub** peak-ring impact (Collins et al. 2020, Nat. Commun. 11:1480,
PDF in ~/Downloads/s41467-020-15269-x.pdf). Calibrating against depth–diameter (not Chicxulub) keeps Chicxulub
a genuine prediction. Vertical craters are axisymmetric → run cheap in 2D (r,z); that's the `crater` mode.

Phase 2b (cylindrical elastic-plastic strength + the full AF model — block activation, shock seeding,
fluidized viscosity, vib transport) is COMPLETE, all gated, GPU==CPU (see PICKUP-phase2b.md).

## SESSION UPDATE 2026-06-28 (run-to-settling DONE; AF creep crux FOUND)
Steps 1-3 of the prior NEXT-STEPS are DONE (not yet committed at time of writing):
- **Run-to-settling** in the `crater` driver (CPU+GPU). Criterion = WINDOWED-MEAN DRIFT of the
  r-weighted excavated volume Vexc: settled when two consecutive means over a window Wwin=2*t_auto
  (>= one oscillation period) agree within tolM=2%; cap=10*t_auto. (A naive max|v|<2 floor and a
  naive 0.5% instantaneous Vexc test were both tried and REJECTED: the settled state has residual
  max|v|~5 m/s sloshing + Vexc has ~few-% density-threshold jitter that never damps — only a
  window-mean is robust.) Progress + AF state stream to stderr (sweep DEVNULLs it).
- **GPU `crater` port** complete (hydro_gpu.cpp): IC + impactor-after-WB-freeze + floor-pin + surface
  diagnostic + CLI + settling, reusing the gated kernels. All 19 GPU + 16 CPU gates still PASS, GPU==CPU.
- **Offline analyzer** `analyze_crater.py` (median-despike + datum-crossing d/D). `sweep_crater.py`
  gained `--gpu`/`--bin PATH` (use --jobs 1 with --gpu; single shared device).

**THE CRUX FINDING (ground-truthed, a=300, Y_d0=10MPa): AF makes craters creep to FLAT.**
- AF ON  (TDEC=20,ETA=2.4e7): transient 3.15km -> floor crept to **0.24km** (d/D garbage), NOT-SETTLED@cap.
- AF OFF (TDEC=0,ETA=0):       transient 1.59km -> floor **HELD at ~1.5km** stably for 700s.
- This **overturns the prior "ROCK friction FIXES it; the crater holds"** finding — that rested on
  UNSETTLED runs (stopped t<320). Run long, AF-on craters creep toward flat and never equilibrate.
- MECHANISM (ground-truthed with af + damage instruments; the D->1 hypothesis was REFUTED):
  * af MEAN decays (0.44->0.007 by t=236) BUT af MAX stays ~0.85-0.99 for a long time -- a SUBSET of
    critical crater-wall cells stay persistently fluidized (continuously re-seeded by the ongoing
    collapse flow), even as the bulk re-solidifies.
  * Damage is NOT the discriminator: AF-OFF is MORE damaged (D mean 0.89, 86% of cells D>0.95) yet
    HOLDS; AF-ON is LESS damaged (D 0.18->0.76) yet CREEPS. So the damaged-friction branch Y_d alone
    holds the crater fine -- the creep is the AF weakening of the still-fluidized wall cells (the
    (1-af) strength cut + the af-scaled viscous flow), NOT a strength/damage deficit.
  * => the lever is the AF ACTIVATION/DECAY (why do wall cells keep re-fluidizing under slow post-
    collapse flow?), not Y_d/mu_d. Lit: W&I AF relaxes to a STABLE final crater, so this persistent
    re-fluidization is most likely a model/param artifact (activation should key on a SHOCK / dP/dt,
    not slow shear) -- the prime suspect to investigate/fix next.
- SECONDARY: the in-loop & offline **D_app are unreliable** — the relaxing rim spreads to the domain
  edge (Rfac=18 too small), contaminating the surf(nx-1) datum. The robust depth metric is
  zsurf-surf (vs the FIXED original surface), which the progress log prints; the RESULT depth (z0-based) is not.

## STATUS — what's DONE (all committed on master)
- **`crater` driver (CPU)** — commit eae390a. Axisymmetric vertical impact: a basalt impactor sphere into a
  basalt half-space under Mars gravity, composing route-1 WB substrate + vacuum-aware free surface + strength
  + Grady-Kipp damage + the full AF model in r-z. Runs STABLY across sizes; dumps the final surface profile
  z_s(r)-z0 for offline measurement.
  CLI: `./hydro_cpu crater [a_m=500] [U=12000] [TDEC=0] [ETA_AF=0] [g=3.71] [cppr=5] [tend(<=0=auto)] [Rfac=18] [Zfac=22] [profile] [rock=1] [Y_d0=0]`
- **ROCK pressure-dependent friction yield (CPU commit f0ac4aa + Y_d0 6111263; GPU commit a4ffdf2; GPU==CPU).**
  Replaces cohesion-only von Mises. Intact `Y_i(P)=Y0+mu_i*P/(1+mu_i*P/(Y_m-Y0))` (Lundborg) + damaged
  `Y_d(P)=min(Y_d0+mu_d*P, Y_m)` (residual cohesion + friction), `Y=[(1-D)Y_i+D Y_d](1-af)`. Params (globals
  in hydro_cpu.cpp; RockP struct on GPU): Y0=1e7, mu_i=1.2, mu_d=0.6, Y_d0=0 (default), Y_m=2.5e9. ROCK flag
  OFF by default → all M3-M5 strength gates bit-unchanged; `crater` defaults ROCK on (arg 12=0 for A/B).
  GATE `friction` (CPU+GPU, in gates.sh): capped stress follows Y_i(P)/Y_d(P) across P=2.7e7..1.4e9, D=0 & D=1,
  to machine precision (CPU) / FP32 (GPU). **GPU ROCK is the only GPU crater piece done so far.**
- **Sweep harness `sweep_crater.py`** — parallel worker pool (`--jobs K`, commit 8f15536). Grids impactor size
  x dimensionless AF (TDEC=Tfrac*a/c_s, ETA=Efrac*rho*c_s*a, size-independent) x Y_d0. Dumps each profile +
  a state/results.csv row; idempotent/resumable; `--list`/`--increment N`/`--daemon`/`--jobs K`.

## SESSION UPDATE 2026-06-30 (first SETTLED sweep + Garvin calibration)
- The 42-run settled sweep RAN (CPU FP64 oracle, cppr=5, cap 6*t_auto; sweep_crater.py --jobs 14 under
  `caffeinate` -- NOTE the throttle was the laptop SLEEPING when unattended, not QoS; caffeinate fixes it).
  Y_d0 {1,5,10 MPa} x AF{off, on Tfrac30/Efrac0.01} x 7 sizes (a=300..3000). Profiles in state/profiles/.
- Calibration oracle `mars_dD_oracle.py` = Garvin et al. 2003 MOLA Mars fresh-crater fits (simple
  d=0.21D^0.81, complex d=0.36D^0.49 km; transition ~7km). `analyze_crater.py ... --csv state/dD_settled.csv`
  then `mars_dD_oracle.py state/dD_settled.csv` prints sim-vs-oracle (ratio = sim d/D ÷ Garvin).
- **FIRST-PASS RESULT (reliable subset only — see blocker): AF-OFF, Y_d0≈5 MPa best tracks Garvin.**
  off_yd5M ratios across the measurable range: 1.03 (D=5.8), 1.14 (8.9), 1.18 (13.1), 0.96 (17.4 km) — and
  d/D DECREASES with D correctly (0.155→0.135→0.115→0.080). Y_d0=1MPa too weak/shallow; Y_d0=10MPa a bit
  deep; **AF-ON (Tfrac=30) systematically OVER-deepens (ratio 1.8-2.1)** — short TDEC + shock-gating lets it
  dig a deeper transient then freeze. AF param space (Tfrac,Efrac) is essentially UNtuned (one slice only).
- **NEW BLOCKER (complex branch): a>=1500 are UNMEASURABLE — D_app pins to the domain edge** (89.7/119.6/
  179.4 km = ~2*r_max). The outer-boundary far-field SAG scales with crater size (sag -3.6km for a=3000),
  AND the big crater fills the domain, so the mid-field-ring datum fails. Only a<=1000 (D<=~17km, simple→just
  past transition) give trustworthy d/D. The complex branch + the ~7km transition shape are NOT yet constrained.
- NEXT: (1) make big craters measurable — proper outer radial BC (reflective/WB so the far field stays flat),
  and/or scale Rfac with a; (2) an AF (Tfrac,Efrac) sweep to see if AF-on can match (and to get complex
  collapse right); (3) resolution convergence (cppr 8/12) + repeat-run scatter (collapse is FP/chaotic) on a
  few anchor sizes before trusting absolute d/D. Single-run d/D here is one chaotic realization.

## SESSION UPDATE 2026-07-01 (#7 DONE: far-field sponge -> big craters measurable)
- Root cause of the big-a d/D failure: the axisym FREE SURFACE SINKS over long runs -- a BOUNDARY-edge
  artifact (outer ~15%) that accumulates in time and propagates inward. a=3000 sank domain-wide to -6.6km
  (center looked like a PEAK). Isolation (small impactor, Rfac=90, 1024s): the INTERIOR is well-balanced
  (max|v|~5.5, flat far field 0..-60m) -- NOT a bulk WB failure, so NO radial-WB core surgery needed.
- FIX = a far-field SPONGE (commits e7eb6ca CPU + a2c1a34 GPU, crater step loop, crater-mode-only, all gates
  PASS GPU==CPU): each step relax the outer 15% radial zone toward the frozen WB reference (ramp sp=xi^2,
  0->1 to the edge), absorbing the edge sink + outgoing waves. a=3000: far field now FLAT, D_app 179km(edge
  garbage)->41km, d/D=0.043 (CPU) / 0.045 (GPU) ~= Garvin complex (~0.05). a=300 unharmed (d/D~0.186).
  NOTE: big-crater collapse is LESS FP-chaotic than small-a (GPU==CPU to ~5% at a=3000 vs ~3x at a=300).
- FULL 42-run sponged sweep DONE (state/dD_settled.csv; view: `python3 mars_dD_oracle.py state/dD_settled.csv`).
  **First complete simple+complex d-D CURVE (D=5..45km).** RESULT vs Garvin (sim d/D, ratio=sim/Garvin):
    off_yd5M : D5.8=0.186(1.2) D9.5=0.168(1.5) D13.4=0.157(1.6) D17.8=0.124(1.5) D25.5=0.082(1.2) D31.6=0.051(0.8) D43.8=0.018(0.35)
    off_yd10M: D5.0=0.229(1.5) D8.7=0.195(1.6) D12.8=0.176(1.8) D16.6=0.157(1.8) D23.7=0.127(1.8) D30=0.093(1.5) D41.4=0.043(0.8)
  FINDINGS: (1) the model REPRODUCES the d/D DECREASE with size (simple->complex) -- invisible before #7.
  (2) LEVEL: Y_d0=5-10MPa BRACKETS Garvin (yd5 in the fresh band at small D; yd10 ~uniformly 1.6x too deep).
  (3) SLOPE: the COMPLEX branch falls TOO STEEPLY -- big craters OVER-collapse (off_yd5M 0.018 at D=44 vs
  Garvin 0.054). AF-ON over-deepens small D and is erratic at big D (untuned Tfrac=30). CAVEATS on the complex
  branch: big-a are NOT-SETTLED@cap (cap 6*t_auto too short for the 12km transients -> the steep falloff is the
  LEAST trustworthy part), + ~20% single-run chaotic scatter (a=300: 0.155 then 0.186 on reruns). The 2 weakest
  biggest configs (af/off_yd1M_a3000) still fail to measure (crater too shallow/spread -> D_app=domain edge).

## KEY FINDINGS (the science — read these before re-running anything)
1. **Cohesion-only strength can't hold craters** → impact damage drives Y→0 → damaged rock flows flat under
   gravity (a=300m crater relaxed 1.4km→0). ROCK friction FIXES this: the crater now holds.
2. **The damaged-cohesion knob Y_d0 works** (monotonic, the calibration lever): in the first (unsettled) pass,
   Y_d0=5 MPa craters were ~1.5-2x deeper than 1 MPa across all sizes. More cohesion → more support → deeper.
3. **THE REAL BLOCKER — craters are NOT SETTLED at the auto-tend `2*sqrt(Rfac*a/g)`.** Profiles show
   concentric collapse ripples (e.g. af_yd1M a=300: surface oscillates -1800/+350/-1750/+350… with radius)
   + a central-axis artifact. So the in-loop d/D and the datum-crossing D_app are GARBAGE (values >1 seen).
   The DEPTH column is the only trustworthy output right now. Settling needs ~several oscillation periods
   (period ~2π√(D/g) ~ hundreds of s) i.e. **~10x longer runs** than the auto-tend.
4. **AF viscosity has an upper bound:** ETA=Efrac~0.01 (≈2.4e7 for a=300) settles fine (max|v|~6 m/s);
   10x (2.4e8) puts the explicit solver in a tiny-timestep regime and never settles (killed at 67 min). Keep
   Efrac small.
5. **CPU speed was a launchd QoS artifact, now solved:** the supervised daemon ran its children at background
   QoS = 32% CPU/core. Running the sweep FOREGROUND (`sweep_crater.py --jobs 14`, normal QoS) = full speed
   (a=3000 in ~5 min, the 42-run 2D pass in ~30 min). The `deorbit-crater` entry in
   ~/investigations/daemons.json was REMOVED — re-add ONLY with a QoS fix, else it throttles. The unsettled
   first pass is archived at state/results_unsettled.csv (depth signal only; d/D unusable).

## NEXT STEPS (ordered — REVISED; AF crux DIAGNOSED + FIXED)
0. **DONE — AF direction chosen (Greg) = shock-gate activation, and IMPLEMENTED.** Root cause (ground-truthed):
   update_af seeded vib on ANY new per-cell pressure peak; with Eulerian cells the slow collapse flow kept
   re-seeding former near-surface cells -> persistent wall fluidization (af max stuck ~0.85) -> creep to flat.
   FIX (commit 1691f1b CPU; GPU ported, gates PASS GPU==CPU): new global/scalar P_ACT/pact -- only a
   SHOCK-sized new-peak rise (dP>P_ACT in a step) seeds vib; default 0 (= legacy, gates bit-identical);
   crater sets P_ACT=1e8 with AF on. VALIDATED a=300 AF-on: af MAX now decays 0.99->0.05 (re-solidifies, was
   stuck ~0.85) and depth HOLDS ~0.57km (was creeping to 0.24). AF now relaxes to a STABLE crater (W&I-intended).
1. **CHAOTIC-COLLAPSE / GPU!=CPU caveat (found porting the fix):** the shock-gated crater HOLDS on both codes
   (no creep) but at DIFFERENT depths -- CPU(FP64) 0.57km vs GPU(FP32) 1.77km for the same a=300 config. The
   sharp dP>P_ACT seeding threshold AMPLIFIES FP divergence (a cell near threshold seeds-or-not differently ->
   different amount fluidizes -> different collapse). Even pre-fix the crept GPU/CPU craters differed ~2-3x, so
   the collapse is intrinsically FP-sensitive/chaotic. IMPLICATION: a single-run d/D is NOT a robust number yet.
   For calibration treat the FP64 CPU as the ORACLE; converge resolution (cppr 8,12); consider SOFTENING the
   threshold (smooth ramp seed ~ smoothstep(dP/P_ACT) instead of a hard gate) and/or more collapse damping to
   tame the sensitivity. (The gates remain GPU==CPU -- this is the chaotic integration, not a kernel mismatch.)
2. **Datum/domain — DONE via measurement workaround (commit 58acb88); d/D now WORKS.** Findings: Rfac 18->30,
   fixed-zsurf datum, and DEPTH-based settling (Vexc never settles -- conflates the stable bowl with slow
   lateral spreading; the bowl depth is the physical signal, settles ~t=886 for a=300). The remaining
   contamination is an OUTER-BOUNDARY SAG (far-field surface droops, present early ~-60m, grows to ~-900m at
   the r=9km edge by t=985 -- worst AT the boundary => outer radial BC + long-time WB drift). WORKED AROUND in
   analyze_crater.py: datum = median surface over the mid-field ring [0.45,0.80]*r_max (excludes inner crater +
   outer-boundary zone), measure d/D below that plain. RESULT: D_app consistent ~6.7km; settled a=300
   (shock-gated AF, Y_d0=10MPa) = depth 1.26km, D_app 6.72km, d/D=0.188 (sensible at the ~7km transition).
   OPTIONAL cleaner fix: make the outer radial face reflective/WB so the far field stays flat (touches the
   flux sweep + re-gate + GPU).
3. **Re-run a SETTLED sweep** (GPU `--gpu --jobs 1`, or CPU oracle) over Y_d0 x AF x size; analyze with the
   local-datum analyzer -> d/D-vs-D table. Use the FP64 CPU as the oracle (GPU collapse is FP-sensitive, item 1).
   [DONE 2026-06-30/07-01 — see the SESSION UPDATE blocks above; stale duplicate steps removed 2026-07-02.]
4. **Calibrate:** digitize the W&I 2003 depth–diameter curve (and Mars d–D data) as the oracle; tune
   (Y_d0, Tfrac, Efrac) so model d/D-vs-D matches, esp. the ~7 km Mars simple→complex transition. Converge
   resolution (cppr 8, 12 on a couple of sizes) before trusting absolute numbers.
6. **Then → 3D oblique Chicxulub** with the FIXED calibrated params (60°, 17 km granite, 12 km/s, NE→SW;
   500 m cells; T=300 s; diagnostics: transient→collapse→peak ring, melt>60 GPa, peak-shock tracer). NOTE the
   EOS gap: we use Tillotson basalt (crust proxy); no ANEOS granite/dunite yet (a dunite-like denser Tillotson
   mantle is the minimal add).

## Build / run / test
- Build CPU: `clang++ -std=c++17 -O2 hydro_cpu.cpp -o hydro_cpu`
- Build GPU: `xcrun -sdk macosx metal -O3 -ffast-math -fmodules-cache-path=.clang-module-cache -c hydro.metal -o hydro.air && xcrun -sdk macosx metallib hydro.air -o hydro.metallib && clang++ -std=c++17 -O2 -I../common/metal-cpp hydro_gpu.cpp -framework Metal -framework Foundation -framework QuartzCore -o hydro_gpu`
- All gates: `./gates.sh` (or `./gates.sh quick` to skip the slow CPU sedov). Every gate prints PASS/CHECK;
  current state = all green (33 PASS in quick). Gates relevant here: `friction` (CPU+GPU).
- One crater: `./hydro_cpu crater 500 12000 30 1e8 3.71 6 -1 18 22 /tmp/p.txt 1 5e6`
- Sweep (FOREGROUND for speed): `python3 sweep_crater.py --jobs 14` (or --list / --increment N).
  Results → state/results.csv; profiles → state/profiles/. state/ is gitignored (local).

## Code map
- `hydro_cpu.cpp`: `crater` mode (the driver + surface diagnostic), `friction` gate, `vonmises` (ROCK branch),
  ROCK globals (ROCK, Y_I0, MU_I, MU_D, Y_D0, Y_M) near the top with GZ/TDEC/etc. step_rk2/Lop unchanged in form.
- `hydro.metal`: `vonmises` kernel (ROCK branch, RockP struct), `strength`/`lop`/`update_af`/`grow_damage`/
  `voidzero` (all the crater machinery). SParams + RockP structs near the top (Metal 31-buffer cap → packed scalars).
- `hydro_gpu.cpp`: the host step loop (the `while(t<tend)` block ~line 195-220) is what a GPU `crater` mode
  reuses; `friction`/`lame`/`af_visc` are pre-loop single-eval gate blocks (the pattern for adding `crater`).
- `sweep_crater.py`: SETTINGS grid (Y_d0 x AF) + SIZES; launch/reap pool; BIN points at hydro_cpu (repoint to
  hydro_gpu after the GPU port).

## Conventions
- Commit as `user.name='oklo' user.email='oklo@mac.com'`; end messages with
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. `git fetch && rebase origin/master` before push
  (a daily cloud routine commits there); keep history linear, no force-push.
- CPU oracle first → gate → GPU port (GPU==CPU), the established rhythm. Ground-truth every claim with a gate.
- euler/ only tracks `*.cpp *.metal *.sh *.md *.py` (gitignore keeps the ~130GB run outputs + binaries out).
  Leave `docs/codex_gpt55_review.md` untracked. Don't commit state/.
- The box is 14 cores (10 perf + 4 efficiency), dedicated to this research → run heavy work FOREGROUND
  (normal QoS), `--jobs 14`. Avoid the launchd supervisor for compute-heavy crater runs (it throttles to 32%).
