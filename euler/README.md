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
- [ ] **M1 Euler hydro (ideal gas)** — gate: Sod shock tube + Sedov blast (analytic).
      M1a CPU reference (scheme correctness) -> M1b Metal port (GPU vs CPU ~1e-6).
- [ ] **M2 Tillotson EOS + multi-material VOF + vacuum/free surface** — gate: 2-material
      shock; clean static free surface (no spurious flow).
- [ ] **M3 strength (elastic-plastic, Jaumann, von Mises)** — THE HARD PART (advecting
      the stress tensor through the remap). gate: elastic-plastic benchmark.
- [ ] **M4 damage (Grady-Kipp) + gravity** — gate: reuse SPH gates; lithostatic balance.
- [ ] **M5 Pierazzo et al. (2008) Al-on-Al impact benchmark** — THE CREDIBILITY GATE:
      match the published code-comparison results. Pass/fail for the whole effort.
- [ ] **M6 acoustic fluidization (block model)** — gate: reproduce a known complex
      crater's depth/rim (calibration-dependent by nature).
- [ ] **M7 Orcus cross-check** — cross-validate vs SPH (early dynamics), then the
      3-way figure: MOLA | SPH | euler.

## Status
M1a: **1st-order HLLC Sod PASS** vs exact Riemann (L1 rho=0.012, p=0.010, u=0.019).
Next: MUSCL (2nd order, sharper) + 3D Sedov blast, then M1b Metal port.
SPH (../gpu/) runs in parallel as the early-dynamics oracle.
