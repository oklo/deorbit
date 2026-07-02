# PICKUP — Pierazzo-2008 aluminum peak-pressure benchmark scale-up

> ## RESUME HERE (2026-07-02 ~15:00, master==origin/master @ 2a6ea95; all work committed+pushed)
> **The head-to-head is DONE and PASSED** — read the SESSION UPDATE block below for the verdict table.
> TWO RUNS STILL IN FLIGHT (both nohup+caffeinate detached, survive session/clear):
>   1. **Batch-3 half24** (3D boundary confirm, ETA ~17:30): `tail -5 state/pierazzo/ladder.log`.
>      Expected: n_far == 1.537 (same as half10/half16). If so: one-line confirm in this file + in
>      [[deorbit-euler-code]]; the 3D 10 km/s study is then fully closed (internal study, not the benchmark).
>      If it DIFFERS: the deep slope is half-sensitive after all -> re-open the boundary question (unlikely;
>      pierazzo2d reach 12 vs 16 was bit-identical).
>   2. **GPU 20 km/s cross-check: DONE + PASS (14:59).** GPU P_cc=380.5 GPa n=2.2796 == CPU cppr20
>      (380.5 / 2.280) to 4 sig figs. E16 GPU claim complete at BOTH velocities.
> **THEN the next Methods anchor (Greg-approved ladder): the 45° OBLIQUE Al benchmark** — 3D, U=5+20 km/s,
> sphere velocity 45° from vertical; measure peak-P along (a) straight down from the impact point and
> (b) 45° downrange starting ~1 projectile radius downrange (paper Fig. 2a geometry: d=sqrt((x-R_pr)^2+z^2)).
> Oracle values in pierazzo2008_oracle.md (Table 1 lower block; iSALE3D: n=1.53±0.02 (0°) / 2.304±0.009
> (45°) at 20 km/s; check the 5 km/s column against the PDF pp. 1923 when implementing — the 45° block
> transcription there is the least certain part of the oracle file). Also expect the max peak-P ~1 R_pr
> DOWNRANGE of the impact point (Fig. 2) — a qualitative check our 3D mode can do. Budget note: 3D at
> cppr 20 was ~1 h (GPU, 9M cells) at 10 km/s; 20 km/s costs ~2-3x steps. After that: Prater-1970 Al-alloy
> strength validation (Table 4 data in the oracle file) -> Chicxulub 60° (Collins 2020) capstone.
> Paper Methods skeleton to keep updated: docs/validation_matrix.md (E16 done; E22/E23/E24 next).

> ## SESSION UPDATE 2026-07-02 (PAPER IN HAND — course correction + `pierazzo2d` built)
> Greg supplied PierazzoEtAl2008.pdf. **The paper's benchmark velocities are 5 and 20 km/s — there is NO
> 10 km/s case**, so the 3D ladder below is an INTERNAL convergence/boundary study, not the head-to-head.
> The quantitative oracle (Table 1 per-code values, CC windows, fit range = d/a 4-10) is extracted to
> **`pierazzo2008_oracle.md`** — no figure digitization needed for pass/fail (Fig. 1 overlay = optional bonus).
> **BUILT (committed on this branch): `pierazzo2d`** — the paper-exact AXISYMMETRIC benchmark mode, CPU
> oracle + GPU port + `gates.sh` gate (both codes, default cppr=12 U=5 km/s). It measures the paper's own
> metrics (P_cc contact/compression mean, decay slope n over d/a 4-10, P at 4 diameters) and gates against
> the PUBLISHED inter-code range. **FIRST RESULTS: 5 km/s PASSES the envelope** (P_cc 39-42 GPa vs 40.4±6.2;
> n 1.26-1.34 vs range 1.13-1.41, converging downward toward iSALE's 1.13-1.2 with cppr; GPU==CPU to 4
> significant digits, ~2 s/run on GPU). 20 km/s: P_cc 358 GPa in range [335,411] at cppr=8; n ladder in
> flight (`state/pz2d/*/out.txt`). 2D axisym cost is trivial (16 s CPU @ cppr12) vs hours in 3D — as planned.
>
> **LADDER VERDICT (14:56, all runs done — BOTH PAPER CASES PASS, CONVERGED ON THE MULTI-CODE MEANS):**
> | U | cppr 8 | 12 | 16 | 20 | 28 | paper mean | code range |
> |---|---|---|---|---|---|---|---|
> | 5 km/s: P_cc GPa | 37.5 | 39.3 | 39.3 | 41.7 | **41.0** | 40.4±6.2 | 28.4-48.0 |
> | 5 km/s: n | 1.444 | 1.341 | 1.261 | 1.263 | **1.243** | 1.2±0.1 | 1.13-1.41 |
> | 20 km/s: P_cc GPa | 357.6 | 366.3 | 375.7 | 380.5 | **383.2** | 379±26 | 334.7-411.1 |
> | 20 km/s: n | 2.144 | 2.202 | 2.252 | 2.280 | **2.316** | 2.3±0.1 | 2.27-2.53 |
> P(4 diam) 5 km/s -> 3.04 GPa (paper 3.2±0.5). Under-resolved runs FAIL exactly per the paper's own
> resolution study (peak P underestimated, slope flattened; cppr>=20 required = the paper's guidance) —
> we reproduce their Table-2 resolution behaviour too. Domain check reach 16 vs 12 @ cppr20: bit-identical.
> The old M5 caveat is CLOSED: the velocity-dependent decay slope (1.2 @ 5 km/s vs 2.3 @ 20 km/s) is
> REPRODUCED; "n~1.3 @ 10 km/s" was simply intermediate. Cost note: 20 km/s cppr28 = 132k steps (~3.5 h
> CPU; ambient tiny-dt pathology) — cppr20 (100k steps) is the practical paper-grade resolution.
> **BATCH-3 half16 VERDICT (12:13): n_far[8-17]=1.537 = IDENTICAL to half=10** -> the 3D deep-slope
> flattening past r/a~10 is NOT lateral-boundary contamination (domain-independent, likely physical
> transition structure at 10 km/s). half24 confirms ~17:30. GPU 20 km/s cppr20 cross-check launched
> (state/pz2d/gpu_u20000_c20/, detached).
> NOTE: 20 km/s runs hit the known ambient tiny-dt pathology (~10x steps of 5 km/s). Batch 3 (3D half=16/24
> boundary check) was found DEAD (died with the prior session ~07:40); RELAUNCHED 09:26 detached
> (nohup+caffeinate, stderr kept) — half16 ~11:50, half24 ~17:15. Validation-matrix skeleton (paper Methods)
> = `docs/validation_matrix.md`. NEXT: read ladder + batch-3 verdicts -> commit -> 45° oblique benchmark
> (Table 1 lower block) as the 3D-oblique anchor -> Chicxulub 60° capstone (Greg-approved for Methods).

> ## RESUME HERE (2026-07-02 morning, superseded by the update above; 3D-ladder context below still valid)
> Scaling up the `pierazzo` 3D Al-on-Al peak-shock-pressure benchmark (domain + resolution) to close the
> **M5 caveat** (near-field n~1.3 unconverged; far-field slope unseen because the old grid reached only r/a~7).
> The parent scoped task + rationale is in **PICKUP-phase3.md → "SCOPED TASK 2026-07-01"** (read it first);
> memory [[deorbit-euler-code]] has the multi-session history. Repo: github.com/oklo/deorbit, code in
> `hydrocodes/euler/`. Runs write to gitignored `state/pierazzo/`; log = `state/pierazzo/ladder.log`.
>
> **HEADLINE RESULT (batch 1 done, all PASS):** the scale-up to r/a=20 resolved the decay STRUCTURE the old
> r/a<=7 benchmark could not see. Local decay exponent n (P ~ (r/a)^-n), converged across cppr 8/10/12:
> | r/a band | 1-3 | 3-6 | **6-10** | 10-15 | 15-20 |
> |---|---|---|---|---|---|
> | n | 0.85 | ~1.6 | **1.75** | ~1.5 | ~1.42 |
> core P converged to **93% of the analytic Hugoniot** (the documented ~7% deficit; plateaus by cppr 10).
> INTERPRETATION (literature-grounded): impact peak-P decay is physically a WEAK-decay near field (isobaric
> core, out to r/a~7-15) -> STRONG far-field power law. The old "n~1.3" sat ENTIRELY in the weak near field --
> NOT a bug. Our extended domain now reaches the far-field onset (n climbs to ~1.75, well converged).
>
> **OPEN QUESTION being tested RIGHT NOW (batch 3 running):** the apparent RE-flattening past r/a~10
> (1.75 -> 1.42) is SUSPECT -- the lateral half-width was half=10 (CLOSER than the r/a=20 measurement depth).
> Batch 3 = cppr10 reach20 at half=16 and half=24 to see if the deep slope steepens with half (=> boundary
> contamination, true far-field ~1.75+) or holds (=> physical flattening). Lateral BC is zero-gradient/outflow
> (hydro.metal lines 156/250) -- so contamination should be modest but not necessarily zero.
>
> **IMMEDIATE NEXT STEPS (ordered):**
>   1. **Read batch-3 verdict:** `tail -30 state/pierazzo/ladder.log`; then re-run the slope-bin table (script
>      below) including decay_half16.txt / decay_half24.txt. Compare the 10-15 & 15-20 bins vs half=10. That
>      settles the far-field slope. If half<measurement-depth is the problem, the RULE is half >~ reach.
>   2. **Digitize the PUBLISHED Pierazzo iSALE Al P(r/a) curve** for the strict "within 10-20% inter-code
>      envelope" pass/fail -- the ONE missing piece for the actual head-to-head. Paper: Pierazzo et al. 2008,
>      MAPS 43:1917 (Wiley doi 10.1111/j.1945-5100.2008.tb00653.x; ADS 2008M&PS...43.1917P; open copy maybe at
>      repository.arizona.edu/handle/10150/656500). WebFetch got 403 on MDPI; may need the PDF from Greg.
>   3. **Commit** the harness + a short findings writeup (see "COMMIT CHECKLIST").
>   4. **STRATEGIC (recommended):** port `pierazzo` to AXISYMMETRIC (2D r,z) -- iSALE itself is 2D-axisym for
>      this benchmark, and 2D makes cppr 20-40 + half>=reach TRIVIALLY affordable (3D here is ~cppr^5-6 cost).
>      The axisym solver is already validated (sedov_axi = on-axis point blast vs analytic Sedov). This is the
>      clean way to a converged, boundary-clean far-field slope AND is truer to "2D iSALE replication."

## WHAT WAS BUILT (uncommitted; `git status` in hydrocodes/euler)
- **`hydro_gpu.cpp` — parameterized the `pierazzo` mode** (was hardcoded cppr=10, r/a<=7). New optional CLI:
  `./hydro_gpu pierazzo [cppr=10] [reach=7] [half=4] [U=10000] [tend=-1auto] [rcfl=0]`
  - `reach` = max on-axis depth in r/a units; `half` = lateral halfwidth in a units (a=1 m); grid built as
    dx=1/cppr, nz=(reach+3)*cppr, nx=ny=2*half*cppr; auto tend = reach/5833*1.35 (front ~5.8 km/s + margin).
  - **Defaults reproduce the M5 gate EXACTLY** (80x80x100, core 92.8%, n=1.262, PASS) -> gates.sh untouched.
  - Diagnostic now also fits a **far-field slope n_far** over [max(6,0.4*reach), 0.85*reach] and prints a
    machine-readable **PZSUMMARY** line (cppr reach half ncells corepct nnear nfar tend steps).
  - **`rcfl` arg is INERT/DO-NOT-USE for pierazzo:** a CFL density floor was TESTED and it DESTABILIZES the run
    (NaN after ~21 steps) -- unlike `crater`, pierazzo has no voidzero, so the live low-density ambient cells
    NEED the small dt. That is why 3D cost is high (~cppr^5-6): the timestep is throttled by those cells and
    their spurious wavespeed grows with resolution (implied cmax ~1.4e6 m/s, ~100x physical). Left as-is;
    the real fix is axisym (step 4) or adding voidzero to pierazzo (a physics change needing M5 re-validation).
- **`run_pierazzo_ladder.sh`** — the runner (deadline-guarded, logs to state/pierazzo/ladder.log). Batch 1 used
  it; batches 2/3 were inline loops appending to the same log.
- Rebuild: `clang++ -std=c++17 -O2 -I../common/metal-cpp hydro_gpu.cpp -framework Metal -framework Foundation -framework QuartzCore -o hydro_gpu`
  (metallib already built; only rebuild it if hydro.metal changed).

## RUNS DONE / IN FLIGHT (all in state/pierazzo/)
- Batch 1 (DONE, 4.4 h, all PASS): conv_cppr8/10/12 @ reach20 half10 + deep_r30 @ cppr10 reach30 half12.
  Decay files: decay_conv_cppr8.txt, _cppr10, _cppr12, decay_deep_r30.txt.
  Numbers: core% 89.2/92.8/92.6/92.8; n_far[8-17] 1.51/1.54/1.58; deep_r30 n_far[12-25.5]=1.66 (BUT its far
  end r/a>~20 is un-arrived noise floor -> ignore; reach20 already fills to r/a~20, so reach30 added no range).
- Batch 2 (KILLED): cppr14/16 @ half10 -- abandoned (half10 => same deep contamination; near field already
  converged). Its partial decay files were removed.
- **Batch 3 (RUNNING): cppr10 reach20 @ half=16 then half=24.** half16 started 07:40, ~2.4 h; half24 ~5.4 h.
  Will write decay_half16.txt / decay_half24.txt + PZSUMMARY lines to ladder.log. If a fresh session finds
  them missing (process died on context-clear/sleep), RE-RUN:
  `cd state/pierazzo && for h in 16 24; do ../../hydro_gpu pierazzo 10 20 $h 10000 >rung_half$h.out 2>/dev/null; mv pierazzo_decay.txt decay_half$h.txt; grep PZSUMMARY rung_half$h.out; done`
  (wrap the whole thing in `caffeinate -i -m -s` -- see GOTCHA).

## SLOPE-BIN ANALYSIS SCRIPT (stdlib only -- no numpy on this box)
Run from `state/pierazzo/`. Add ('half16','decay_half16.txt') etc. to `files` once batch 3 lands.
```python
python3 -c "
import math
bins=[(1.2,3),(3,6),(6,10),(10,15),(15,20)]
files=[('cppr8','decay_conv_cppr8.txt'),('cppr10','decay_conv_cppr10.txt'),('cppr12','decay_conv_cppr12.txt')]
def slope(xs,ys):
    n=len(xs);sx=sum(xs);sy=sum(ys);sxx=sum(x*x for x in xs);sxy=sum(x*y for x,y in zip(xs,ys))
    return (n*sxy-sx*sy)/(n*sxx-sx*sx)
print(f'{\"run\":11s} '+' '.join(f'{lo}-{hi}'.rjust(8) for lo,hi in bins))
for name,f in files:
    ra=[];P=[]
    for line in open(f): a,b=line.split();ra.append(float(a));P.append(float(b))
    core=max(p for r,p in zip(ra,P) if 0.5<=r<1.5); row=[]
    for lo,hi in bins:
        xs=[math.log(r) for r,p in zip(ra,P) if lo<=r<=hi and p>0.005*core]
        ys=[math.log(p) for r,p in zip(ra,P) if lo<=r<=hi and p>0.005*core]
        row.append(f'{-slope(xs,ys):8.2f}' if len(xs)>2 else '   --   ')
    print(f'{name:11s} '+' '.join(row))
"
```

## GOTCHA (cost me ~6.5 h overnight)
Batches 2/3 were launched WITHOUT `caffeinate` -> the laptop SLEPT ~01:00-07:40 and the run stalled the whole
time (matches the known "unattended slowdown = laptop SLEEPING" note in [[deorbit-euler-code]]). ALWAYS wrap
long runs: `caffeinate -i -m -s ./whatever` (or `nohup caffeinate -i -m -s -t <sec> &` as a standalone keep-awake).
A protective `caffeinate -t 32400` (~9 h) was started at pickup time to cover batch 3.

## COMMIT CHECKLIST (when the picture is complete)
- Commit `hydro_gpu.cpp` (parameterized pierazzo + far-field diagnostic + inert rcfl arg) and
  `run_pierazzo_ladder.sh` and this `PICKUP-pierazzo.md`. NOT `state/` (gitignored). Confirm `./gates.sh quick`
  still green (esp. `[pierazzo] ... PASS`, unchanged).
- Add a dated "SESSION UPDATE" block to PICKUP-phase3.md with the headline numbers + the boundary verdict.
- Commit as user.name='oklo' user.email='oklo@mac.com'; end msg with
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`; `git fetch && rebase origin/master` before push.
- Update memory [[deorbit-euler-code]] with the result (near->far decay structure resolved; far-field-onset
  n~1.75; 93% Hugoniot; the old n~1.3 = near-field-only artifact EXPLAINED).
