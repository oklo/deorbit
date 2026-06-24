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
- [~] **M4 damage (Grady-Kipp) + gravity.** Gravity DONE (CPU+GPU): body-force source;
      free-fall v=-g*t exact, hydrostatic atmosphere interior well-balanced (5e-4 of cs).
      Damage CPU DONE: D scalar advected, grows via crack ODE d^(1/3)+=(cg/Rs)dt when
      tensile strain sigmax/Emod > Weibull activation (1/(wk*dx^3))^(1/wm); degrades
      Y->(1-D)Y. GATE: tensile stretch -> D->1, shear strength -> 0. Next: M4 damage GPU.
- [ ] **M5 Pierazzo et al. (2008) Al-on-Al impact benchmark** — THE CREDIBILITY GATE:
      match the published code-comparison results. Pass/fail for the whole effort.
- [ ] **M6 acoustic fluidization (block model)** — gate: reproduce a known complex
      crater's depth/rim (calibration-dependent by nature).
- [ ] **M7 Orcus cross-check** — cross-validate vs SPH (early dynamics), then the
      3-way figure: MOLA | SPH | euler.

## Status
**M3 DONE (CPU + GPU)** — elastic-plastic strength validated on both paths
(shear-wave speed 1.3%, von Mises cap exact, ideal-gas no-op regression). The
flagged risk gate passes. Next: **M4** (Grady-Kipp damage + gravity, mostly reused
from SPH), then **M5 Pierazzo** -- the pass/fail credibility gate. SPH runs in parallel.
