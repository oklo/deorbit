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
- [~] **M3 strength (elastic-plastic, Jaumann, von Mises).** M3a (CPU) DONE: deviatoric
      stress S (6 comp), Jaumann rate dS=2G*edev+(SW-WS), v.grad advection, div S ->
      momentum + div(S.v) -> energy, von Mises radial return. GATES: elastic shear
      wave c_s 2863 vs sqrt(G/rho) 2900 (1.2%); von Mises cap sqrt(3J2)->Y exactly;
      Sod/surface regress (strength no-op when G=0). Next: M3b Metal port.
- [ ] **M4 damage (Grady-Kipp) + gravity** — gate: reuse SPH gates; lithostatic balance.
- [ ] **M5 Pierazzo et al. (2008) Al-on-Al impact benchmark** — THE CREDIBILITY GATE:
      match the published code-comparison results. Pass/fail for the whole effort.
- [ ] **M6 acoustic fluidization (block model)** — gate: reproduce a known complex
      crater's depth/rim (calibration-dependent by nature).
- [ ] **M7 Orcus cross-check** — cross-validate vs SPH (early dynamics), then the
      3-way figure: MOLA | SPH | euler.

## Status
**M3a DONE (CPU)** — elastic-plastic strength validated: shear-wave speed 1.2%,
von Mises cap exact. The flagged risk gate passes. Next: **M3b** Metal port of the
strength kernels (GPU==CPU). Then M4 (damage + gravity), M5 (Pierazzo credibility
gate). SPH (../gpu/) runs in parallel.
