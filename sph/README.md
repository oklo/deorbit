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

## What's implemented (`sph.hpp`)

- cubic-spline kernel (3D); per-dimension cell-list neighbours, periodic in y,z;
- **adaptive smoothing length** (gather density with h_i; h̄-symmetrized,
  momentum-conserving forces);
- Monaghan artificial viscosity + Price artificial conductivity;
- KDK leapfrog with trapezoidal energy update;
- pluggable EOS (`IdealGas` now; Tillotson next).

## Roadmap to the impact runs

1. **Tillotson EOS** for iron / basalt / water (+ a von Mises strength model).
2. **Gate 2 — vertical impact π-scaling**: reproduce crater-size scaling for a
   vertical impact before trusting oblique runs.
3. **Oblique runs**: the Phase-1 terminal conditions (v≈7.4 km/s, γ≈1–4°) into
   the three targets; compare gouge/ricochet/melt to the Phase-2a estimates.

## Build & validate

```bash
cd sph && make && ./sod && python3 analyze_sod.py
```
