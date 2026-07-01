# PICKUP — Euler Phase 3: AF + strength calibration to the depth–diameter curve

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
- Re-running the FULL 42-run sweep WITH the sponge (state/dD_settled.csv) for the complete simple+complex
  d-D curve. CAVEAT still: single-run d/D has ~20% chaotic scatter (a=300 gave 0.155 then 0.186 on reruns) ->
  want repeat-run averaging on anchor sizes; and big-a are NOT-SETTLED@cap (max|v|~74, 12km transient still
  relaxing) so their d/D is a lower bound on settling-completeness.

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
2. **Fix the datum/domain:** measure depth against the FIXED original surface zsurf (already in the progress
   log), not surf(nx-1); and/or enlarge Rfac (18 -> ~30) so the relaxing rim doesn't reach the boundary.
   Make the RESULT depth + the offline analyzer use the robust datum so d/D becomes meaningful.
3. **Re-run a SETTLED sweep** (GPU, `sweep_crater.py --gpu --jobs 1`) over Y_d0 x AF x size — only AFTER the
   AF direction is set and the datum is fixed, else the d/D column stays garbage.
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
