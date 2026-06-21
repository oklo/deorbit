# 3D SPH hydrocode (Phase 2b)

A from-scratch 3D smoothed-particle-hydrodynamics solver in C++ for the
ultra-oblique terminal impact (Cayambe / deep ocean / Amazon). SPH is the
standard tool for oblique and marine impacts: naturally 3D, handles large
deformation, free surfaces, and ricochet. Built incrementally behind hard
validation gates.

## Status

**Gate 1 — Sod shock tube (DONE, core validated).** 3D slab, ideal gas,
compared to the exact Riemann solution (`analyze_sod.py`), with **adaptive
smoothing length** (h_i = η(m_i/ρ_i)^{1/3}, iterated each step):
- shock-front speed correct (within ~3 bins; SPH shocks smear over ~2-3 h);
- density & pressure RMS < 5% in the smooth regions (3.9% / 4.5%);
- all wave structure (rarefaction, contact, shock) reproduced.
- velocity RMS ~9%: **method-intrinsic** to standard SPH + artificial viscosity.
  NB adaptive-h only marginally improved this over constant-h (it mainly helps
  density); the velocity/shock smearing is the documented accuracy of the method
  at this resolution, not a bug. grad-h correction terms or a Riemann-SPH /
  modern viscosity switch would tighten it — diminishing returns for our goals
  (impact gouge/melt/ricochet need ~10%). Adaptive-h is kept because it is
  *required* for the large density contrasts in impacts, where constant-h fails.

Perf note: the adaptive-h Sod (≈29k particles, 3 ρ↔h iterations/step) took ~7 min.
Before large impact runs this needs optimization (persistent neighbour lists,
fewer h-iterations, OpenMP).

**Gate 1b — Tillotson EOS (DONE, `eos_test.cpp`).** Iron / basalt / water:
P(ρ₀,0)=0 exactly, and the small-strain sound speed c₀=√(A/ρ₀) matches known
bulk sound speeds (iron 4051, basalt 3145, **water 1485 vs. real 1481 m/s**);
10%-compression pressures sensible (iron 13.9 GPa, basalt 2.9, water 0.32).

**Optimization (DONE).** Density + force loops parallelised with `std::thread`;
round-based cell rebuilds. Sod (≈29k particles) went **~407 s → 37 s (~11×)** on
14 cores, identical accuracy.

**Gate 1c — strength model (DONE, `strength_test.cpp`).** Elastic-perfectly-
plastic deviatoric stress (Jaumann rate, von Mises yield):
- elastic shear: dS_xy/dt = 2G·ε̇_xy to **0.2%**;
- von Mises yield: radial return caps the von Mises stress at Y exactly;
- elastic **longitudinal wave speed** 4250 vs. c_L=√((K+4G/3)/ρ)=4593 m/s (7.5%
  low — SPH dispersion), validating the stress-divergence coupling in momentum.

## What's implemented (`sph.hpp`, `eos.hpp`)

- cubic-spline kernel (3D); per-dimension cell-list neighbours, periodic in y,z;
- **adaptive smoothing length** (gather density with h_i; h̄-symmetrized,
  momentum-conserving forces);
- Monaghan artificial viscosity + Price artificial conductivity;
- KDK leapfrog with trapezoidal energy update;
- **multi-material EOS**: ideal gas + **Tillotson** (iron/basalt/water);
- **material strength**: elastic-perfectly-plastic deviatoric stress (Jaumann
  rate, von Mises radial-return yield), with the elastic signal speed in AV/CFL;
- **`std::thread` parallelism** over particles (density + forces).

**Gate 2 — crater π-scaling (QUALIFIED PASS, `pi_scaling.cpp`, see
`pi_scaling_results.txt`).** Vertical iron→basalt impacts, U=0.6–2.2 km/s,
CPPR≈4. The cratering-efficiency (volume) velocity scaling gives coupling
exponent **μ ≈ 0.58**, matching the canonical strength-regime value (~0.55) for
competent rock → the energy/momentum coupling is correct. Caveat: at CPPR≈4 the
crater is penetration-dominated (deep+narrow), so absolute dimensions/shape are
resolution-limited; a fully-resolved validation needs CPPR ≥ 8–10 (a cloud/HPC
run, not a laptop sweep).

The physics is complete (EOS + strength) and the scaling is validated. Roadmap:

1. **Oblique runs**: the Phase-1 terminal conditions (v≈7.4 km/s, γ≈1–4°) into
   the three targets; compare gouge/ricochet/melt to the Phase-2a estimates.
   (Same resolution caveat → laptop gives qualitative results; quantitative
   dimensions want a higher-resolution cloud run.)

## Build & validate

```bash
cd sph && make && ./sod && python3 analyze_sod.py
```
