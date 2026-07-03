# Validation matrix — deorbit codes (paper Methods skeleton)

Status: DRAFT 2026-07-02. One row per verification/validation item, for the peer-reviewable
methods section. Terminology follows Pierazzo et al. 2008: **verification** = solves the
equations right (analytic/manufactured oracles); **validation** = solves the right equations
(experiments); **benchmark** = code-to-code comparison; **calibration** = parameter tuning
(NOT validation — list separately, a referee will check we don't conflate them).

Legend: ✅ done+gated · 🟡 done, caveat noted · 🔶 in flight · ⬜ planned · ❌ known gap

## 1. Euler code (GPU/Metal + FP64 CPU oracle, `hydrocodes/euler/`)

Every row is a scripted gate in `gates.sh` (prints PASS/CHECK; run `./gates.sh quick`).
All gates currently PASS on both codes, GPU==CPU.

| # | Test (gate) | Oracle | Result | Status |
|---|---|---|---|---|
| E1 | Sod shock tube (`sod`) | exact Riemann | L1(rho) ~0.004 (MUSCL) | ✅ |
| E2 | Sedov 3D point blast (`sedov`) | analytic self-similar | shock radius <1% (FP64), 6.1% @64^3 note | ✅ |
| E3 | Sedov axisymmetric (`sedov_axi`) | same analytic (on-axis blast = spherical) | R err 4.2% — validates cylindrical geometry | ✅ |
| E4 | Tillotson free surface static (`surface`) | v=0 | max\|v\|=0, no gravity | ✅ |
| E5 | Basalt shock GPU==CPU (`bshock`) | FP64 CPU | GPU-CPU 3.3e-7 | ✅ |
| E6 | Elastic shear wave (`shear`) | c_s=sqrt(G/rho) | 1.3% | ✅ |
| E7 | Von Mises yield cap (`yield`) | sqrt(3 J2)<=Y | exact cap | ✅ |
| E8 | ROCK friction yield (`friction`) | Lundborg Y_i(P) + damaged Y_d(P) analytic | machine precision (CPU) / FP32 (GPU) | ✅ |
| E9 | Free fall (`freefall`) | v=-g t | err <1e-3 | ✅ |
| E10 | Hydrostatic atmosphere (`atmos`) | static column | interior max\|v\|<0.02 c_s | ✅ |
| E11 | Gravity-loaded substrate w/ free surface (`substrate`) | stays static, no damping | max\|v\|~1.4 m/s, 200 s (Audusse WB) | ✅ |
| E12 | Vacuum Riemann flux (`vacuum`) | Toro rarefaction-into-vacuum fan | L1 rho 0.0012, u 0.0027 | ✅ |
| E13 | Grady-Kipp tensile damage (`tensile`) | D->1, strength->0 | as expected | ✅ |
| E14 | Material tracer (`tracer`) | conservation + advection | sum 1e-7, centroid 1e-7 | ✅ |
| E15 | 1D Al-on-Al planar impact (`alimpact`) | analytic Tillotson Hugoniot (R-H bisection) | peak P err 0.0% | ✅ |
| E16 | **Pierazzo-2008 Al benchmark, axisym (`pierazzo2d`)** | published inter-code envelope, Table 1 (see `hydrocodes/euler/pierazzo2008_oracle.md`) | **BOTH cases PASS, converged @ cppr 20-28**: 5 km/s P_cc 41.0 GPa (paper 40.4±6.2), n 1.243 (1.2±0.1), P(4D) 3.0 GPa (3.2±0.5); 20 km/s P_cc 383 GPa (379±26), n 2.316 (2.3±0.1). Velocity dependence of the decay slope reproduced; under-resolution fails per the paper's own Table-2 study (cppr>=20 = their guidance); domain-doubling bit-identical; GPU==CPU 4 sig figs | ✅ |
| E17 | Pierazzo 3D decay structure (`pierazzo`, GPU) | Hugoniot core + power-law decay | core 93% Hugoniot; near-field n~1.26 -> far-field onset ~1.75; deep-slope flattening domain-INDEPENDENT (n_far=1.537 identical across half 10/16/24 — closed 2026-07-02) | ✅ internal study (no 10 km/s paper case) |
| E18 | AF activation/decay (`af_activate`) | shock seeds vib behind front only; af decays | exact behaviour, GPU==CPU | ✅ |
| E19 | Cylindrical strength — hoop strain + Lame equilibrium (`lame`) | analytic thick-cylinder deviatoric field | 1e-15 (CPU) / 2.6e-6 (GPU); equilibrium residual 0.4% | ✅ |
| E20 | Cylindrical AF viscosity (`af_visc`) | (lap v)_r = 3A for u=Ar^2 (vs 2A Cartesian) | 4e-13 (CPU) / 4.3e-4 (GPU) | ✅ |
| E21 | vib advection (`vib_advect`) | translation at v0, conservation | exact | ✅ |
| E22 | **45° oblique Al benchmark (3D, `pierazzo ... ang=45`)** | Pierazzo-2008 Table 1 lower block (n 0°/45° over d/a 4-10; peak-P max ~1 R_pr downrange), verified vs PDF p. 1923 | **BOTH velocities PASS the published envelope at paper-grade resolution**: 5 km/s cppr20 n0=1.018/n45=1.182 (means 1.1±0.1/1.1±0.1, ranges [0.90,1.226]/[0.87,1.318]); 20 km/s cppr16 n0=1.645/n45=2.264 (means 1.5±0.3/2.1±0.6, ranges [1.06,1.95]/[1.11,2.57]); converging onto the means with cppr; peak-P max 0.7-1.7 R_pr downrange and the paper's slope ordering (n0_oblique < n45 < n_vertical at 20 km/s) reproduced; two free-surface treatments agree ≲4% (legacy ambient-live vs rcfl voidzero-to-ambient — the latter ~200x cheaper and required at 20 km/s cppr≥16 where legacy NaNs); 20 km/s cppr20 rung in flight (dt pathology partially survives rcfl at that corner) | ✅ (cppr20 20 km/s pending) |
| E23 | **Al validation #2 (`prater`): 6.35 mm Al sphere -> Al 6061-T6 @ 7 km/s** | EXPERIMENT (Prater 1970; Table 4 R(t),d(t) transcribed in `prater1970_oracle.md`) | von Mises Y=414 MPa (iSALE's value, untuned): growth phase (3-15 us) INSIDE the inter-code band (cppr20: R -8.3%, D -1.4% vs codes' R -5..-13%); plateau low (R -22%, D -10%; codes +4..+12% on D) — no-rim-uplift/wall-extrusion limitation of the strengthless vacuum-face flux, documented; converging upward with cppr; CPU==GPU (c10: D identical, R within 2.5%); gates in gates.sh (cppr6 smoke quick, cppr10 envelope full) | 🟡 |
| E24 | **Chicxulub 60° peak-ring reproduction (3D oblique)** | Collins et al. 2020 (published iSALE3D application) | — the Methods capstone (Greg-approved) | ⬜ |

**Calibration (separate from validation):** Mars d/D-vs-D curve vs Garvin et al. 2003 MOLA fits
(`mars_dD_oracle.py`): first complete simple+complex curve done (42-run settled sweep, sponge fix);
Y_d0=5-10 MPa brackets the level; d/D decrease with D reproduced. **INVALIDATED IN PART 2026-07-03:
the voidzero below-surface refill artifact (void cells reset to the lithostatic reference) was
driving the complex-branch over-collapse — A/B at a=3000 (off_yd5M): d/D 0.018 (legacy, ==sweep)
vs 0.165 (ambient voids), the refills injecting ~a bowl's worth of rock at the floor during
collapse. Simple branch mildly affected (a=300: 0.186 vs 0.222). The sweep must be re-run with
`vamb=0.27` (crater arg14) and Y_d0/AF recalibrated on the clean base (details:
hydrocodes/euler/PICKUP-phase3.md 2026-07-03 block).** Other caveats stand: ~20% single-run
chaotic scatter (needs repeat-run averaging); cppr 8/12 convergence pending. W&I-2003 comparison
= the intended calibration narrative for AF params (Tfrac, Efrac).

**Known Euler gaps (state in the paper):** Eulerian tracer advection is diffusive (partition yes,
fine fragmentation no — that stays SPH's job); no ANEOS (Tillotson only: basalt/Al proxies, no
granite/dunite); crater collapse is chaotic/FP-sensitive (FP64 CPU = oracle; single-run d/D ±20%);
`collapse` M6 demo gate still CHECK (physically explained; not in the pass/fail set).

## 2. SPH code (C++ CPU + Metal GPU, `hydrocodes/sph/`)

| # | Test | Oracle | Result | Status |
|---|---|---|---|---|
| S1 | Sod shock tube | exact Riemann | rho/P RMS <5%; v ~10% = standard-SPH contact artifact (method-intrinsic, documented) | 🟡 |
| S2 | Tillotson EOS unit tests (`eos_test`) | P(rho0,0)=0; c0=sqrt(A/rho0) vs known | water 1485 vs 1481 m/s, iron 4051, basalt 3145 | ✅ |
| S3 | Strength unit tests (`strength_test`) | dS_xy/dt=2G eps-dot; yield cap; elastic wave | 0.2%; exact cap; c_L 4250 vs 4593 (7%) | 🟡 |
| S4 | Crater pi-scaling (`pi_scaling`) | coupling exponent mu ~0.55 (strength regime) | mu=0.58 for cratering efficiency; BUT CPPR~4 => absolute dims resolution-limited | 🟡 needs CPPR>=8-10 rerun |
| S5 | GPU port parity | CPU SPH oracle | validated, ~70x faster (bandwidth-bound) | ✅ |
| S6 | Regime validation: Barringer (45°, v/v_esc 1.14) | real crater (Canyon Diablo) | classifier -> C (crater, destroyed projectile) = reality | ✅ regime-level |
| S7 | Regime validation: Campo del Cielo (9°) | real strewn field | classifier -> I (ricochet, 89% v retained); real fragmentation was atmospheric (unmodeled) — regime confirmed | 🟡 |
| S8 | Boundary-point convergence CPPR 8->10 | self-convergence | — phase-diagram boundaries single-resolution so far | ⬜ |
| S9 | Oblique-impact literature anchors | Gault & Wedekind 1978 ricochet onset ~30°; Pierazzo & Melosh 2000 | crater onset theta_crit=36.2±0.6° consistent; make the comparison explicit+quantitative | 🟡 |
| S10 | Euler cross-check on ricochet fate (M-tag partition/bulk-v) | independent method | — planned on 2-3 phase-diagram anchors | ⬜ |

**Known SPH gaps:** standard-SPH velocity contact artifact (S1); no surface tension / liquid
fragmentation (iron-blob lacework is numerical — documented honestly in the Cayambe analysis);
tensile instability at free surfaces (why the Euler code exists for late-stage morphology).

## 3. Orbital layers (`src/deorbit/`)

| # | Test | Oracle | Result | Status |
|---|---|---|---|---|
| O1 | Corridor integrator (pytest, CI) | Kepler closure; per-pass dv vs analytic column | passes, CI-enforced on every push | ✅ |
| O2 | Hi-fi REBOUND+ASSIST (`validate_hifi.py`) | Kepler 6.6e-8/orbit; J2 nodal regression vs analytic 1.2%; drag dv vs stdlib ~10%; Moon at 402,449 km | passes | ✅ |
| O3 | Uncertainty ensemble (N=160) | MC over Cd/atmosphere/rotation | corridor narrow-but-robust; vmax p05-p95 quantified | ✅ |
| O4 | Site selection | WGS84 exact; Chimborazo 6384.39 km global max reproduced | known-fact check passes | ✅ |

## 4. Cross-cutting to-dos for the Methods section
1. E16 finish (convergence ladder both velocities) + overlay figure vs digitized paper Fig. 1.
2. E17 boundary verdict -> fold the 3D far-field-structure study into supporting material.
3. S4 rerun at CPPR>=8-10 (cloud/overnight GPU) — the referee-obvious SPH weakness.
4. Repeat-run averaging + settling-cap bump for the d/D calibration (before AF tuning).
5. E22 45° benchmark -> E24 Chicxulub 60° (Collins 2020) as the 3D-oblique ladder.
6. Write the "honest limitations" paragraph from the gap lists above (referees reward this).
