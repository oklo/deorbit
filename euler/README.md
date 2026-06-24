# euler — GPU Eulerian shock-physics code ("mini-iSALE")

A Metal/GPU Eulerian multi-material shock-physics solver with strength, damage, and
acoustic fluidization — built as an **independent cross-check** of the SPH solver
(`../gpu/`) for impact cratering, where SPH is weak (free-surface noise, tensile
instability, no late-stage collapse). Motivated by the Orcus Patera oblique-impact
study: SPH gives the early dynamics + projectile fate, but the final crater
morphology needs a grid code. We don't have iSALE-3D access, so we build a
principled, **validation-gated** equivalent.

## Why Eulerian + GPU
- Eulerian regular grid: no tensile instability, clean free surface, handles ejecta
  + large deformation — exactly SPH's weak spots.
- Regular grid is *more* GPU-friendly than SPH (coalesced memory, no neighbour
  search) → iSALE's hours-to-days CPU runs become minutes-to-hours.

## Architecture
Operator-split Eulerian finite-volume (the iSALE/CTH lineage):
- hydro: dimensionally-split MUSCL-Hancock reconstruction + HLLC Riemann solver
- multi-material: VOF volume fractions + vacuum cells (free surface)
- strength: elastic-plastic (Jaumann rate, von Mises), advected across the grid
- damage: Grady-Kipp (reused from the SPH)
- acoustic fluidization: Wünnemann-Ivanov block model (the key collapse ingredient)
- EOS: Tillotson (reused) / ideal gas (for the hydro gates)
- gravity body force
Reuses the SPH port's Metal host framework + EOS/strength/damage + validation gates.

## VALIDATION IS THE POINT
A home-built hydrocode is only a credible cross-check if it passes the field's
standard benchmarks. Every milestone is gated against an oracle (analytic or
published). **We do not advance on an unvalidated stage, and we stop honestly if
M3 or M5 won't validate.**

## Milestones (gated)
- [x] **M1a Euler hydro (ideal gas), CPU reference** — gates PASSED: Sod MUSCL L1
      rho=0.004/p=0.003/u=0.008 vs exact Riemann; Sedov 3-D shock radius 0.8% vs
      analytic. Unsplit MUSCL(minmod)+HLLC+SSP-RK2. (peak compression resolution-limited.)
- [x] **M1b Metal GPU port** — DONE: FP32 GPU matches FP64 CPU to printed precision
      on both gates (Sod L1 0.0042 identical; Sedov R 0.776 / 0.8% / comp 3.12 identical).
      Per-cell divergence (MUSCL+HLLC, 13-pt stencil) + RK2; regular grid = clean GPU map.
- [x] **M2 Tillotson EOS + free surface — DONE (CPU + GPU).** Pluggable EOS
      (ideal | Tillotson via ../sph/eos.hpp & MSL port of sph_force.metal's
      till/ssound); reconstruct internal energy; P>=0 fluid floor + low-density
      ambient. GATES: Sod regression L1 0.0044 (CPU==GPU); basalt free surface
      max|v|=0.000 (static, CPU==GPU); basalt shock GPU-vs-CPU rel 3.3e-7.
      (Impactor color tag deferred to M7 -- passive scalar, not a physics gate.)
- [x] **M3 strength (elastic-plastic, Jaumann, von Mises) — DONE (CPU + GPU).**
      Deviatoric stress S (6 comp): Jaumann rate dS=2G*edev+(SW-WS), v.grad advection,
      div S -> momentum + div(S.v) -> energy, von Mises radial return; gated on G>0.
      GATES (CPU==GPU): elastic shear wave c_s 2863 vs sqrt(G/rho) 2900 (1.3%); von
      Mises cap sqrt(3J2)->Y exactly; Sod/surface regress (strength no-op for ideal gas).
      THE FLAGGED RISK GATE PASSES on both paths. Metal kernels: strength/vonmises/rk1s/rk2s.
- [x] **M4 damage (Grady-Kipp) + gravity — DONE (CPU + GPU).** Gravity: body-force
      source; free-fall v=-g*t exact, hydrostatic atmosphere interior well-balanced
      (5e-4 of cs). Damage: D scalar advected, grows via crack ODE d^(1/3)+=(cg/Rs)dt
      when tensile strain sigmax/Emod > Weibull activation (1/(wk*dx^3))^(1/wm); degrades
      Y->(1-D)Y. GATES (CPU==GPU): tensile stretch -> D->1, shear strength -> 0;
      free-fall + hydrostatic; all prior gates regress. (wk=1e61 overflows FP32 ->
      eps_act precomputed on host.)
- [x] **M5 Pierazzo Al-on-Al shock-physics validation (CREDIBILITY GATE).** Validated on
      the credibility-critical quantity, the PEAK SHOCK PRESSURE: (M5a) 1D Al-on-Al planar
      impact U=10 km/s reproduces the analytic Tillotson Hugoniot EXACTLY (1.638e11 vs
      1.637e11 Pa, 0.0%); (M5b) 3D Al-sphere impact (10 cppr) isobaric core = 93% of the
      Hugoniot (7% deficit = 3D geometry + resolution, per Pierazzo), clean monotonic peak-
      pressure decay P(r/a) (per-cell Pmax tracker). CAVEAT: near-field decay exponent ~1.3
      over r/a 1-6 (the transition regime is shallower than the canonical far-field ~2; a
      far-field domain + higher cppr convergence study would pin it -- a refinement, not done
      here per the option-1 physics-based plan). The shock-capture + EOS give correct impact
      pressures -> the cross-check is credible. Metal: pmax_update kernel; mode 'pierazzo'.
- [~] **M6 acoustic fluidization (block model) — PARTIAL (honest checkpoint).** AF
      mechanism IMPLEMENTED + regression-safe: fluidization field `af` reduces shear
      strength Y->(1-D)(1-af)Y; Maxwell-style decay (TDEC); Newtonian viscosity
      ETA_AF (eta*grad^2 v) to damp fluidized flow; void-CFL exclusion (RHO_CFL) so
      near-vacuum ambient cells don't tank dt. BLOCKED on the collapse DEMO by a
      separate numerics issue: the sharp basalt|ambient free surface is NOT hydrostatically
      well-balanced under gravity (surface cells get a spurious net force -> rarefy; M2
      surface gate passed only because it had no gravity). Needs a well-balanced interface
      reconstruction (or fixed-substrate treatment) -- this ALSO gates M7's loaded substrate.
      Next focused work: free-surface + gravity well-balancing, then the collapse demo + calibration.
- [ ] **M7 Orcus cross-check** — cross-validate vs SPH (early dynamics), then the
      3-way figure: MOLA | SPH | euler.

## Status
**M5 DONE** (credibility gate, peak shock pressure). **M6 PARTIAL**: acoustic-fluidization
mechanism implemented (strength reduction + decay + viscosity + void-CFL), regression-safe,
but the collapse demo is blocked by free-surface + gravity well-balancing (a sharp material
interface under gravity isn't hydrostatically balanced -> surface cells rarefy). This same
issue gates M7's loaded substrate. **Next focused work: well-balanced free-surface gravity**
(then M6 collapse demo + M7 Orcus 3-way). 5 of 7 milestones validated; M6/M7 need the
well-balancing fix. SPH Orcus run done; renders committed.
