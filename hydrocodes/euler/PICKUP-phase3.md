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

## NEXT STEPS (ordered — this is the remaining Phase 3 work)
1. **Add run-to-settling to the `crater` driver (CPU first).** Replace the fixed tend: run until max|v| over
   dense (rho>1350) cells drops below a small threshold (e.g. a few m/s, or ~0.01*sqrt(g*a)) sustained, capped
   at a generous max tend (~10x the current auto-tend). Report the settled crater. (Optional: also report the
   time-to-settle.) This makes profiles smooth/final so d/D is meaningful.
2. **Port the `crater` driver to GPU** (`hydro_gpu.cpp`). The step LOOP already exists on GPU (lop/strength/
   vonmises/update_af/grow_damage/voidzero, with wb + axisym + rvac + ROCK all present). Remaining: the
   impactor+WB-substrate+ambient IC, the deep-floor pin each step (mirror the CPU substrate floor-pin), the
   surface-profile diagnostic + RESULT print, the CLI parse, and run-to-settling. Validate GPU crater ≈ CPU
   crater on one config. GPU makes the long settled runs (step 1) affordable.
3. **Robust OFFLINE d/D measurement** from the dumped profiles (write a small analyzer): floor depth, apparent
   diameter at the datum (robust edge logic + smoothing, handle the central-axis cell), rim crest. Apply to
   the settled profiles → a clean d/D-vs-D table per (Y_d0, AF).
4. **Re-run a SETTLED sweep** (GPU, `sweep_crater.py` repointed to hydro_gpu) over Y_d0 x AF x size.
5. **Calibrate:** digitize the W&I 2003 depth–diameter curve (and Mars d–D data) as the oracle; tune
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
