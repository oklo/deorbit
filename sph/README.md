# 3D SPH hydrocode (Phase 2b)

A from-scratch 3D smoothed-particle-hydrodynamics solver in C++ for the
ultra-oblique terminal impact (Cayambe / deep ocean / Amazon). SPH is the
standard tool for oblique and marine impacts: naturally 3D, handles large
deformation, free surfaces, and ricochet. Built incrementally behind hard
validation gates.

## Status

**Gate 1 — Sod shock tube (DONE, core validated).** 3D slab, ideal gas,
compared to the exact Riemann solution (`analyze_sod.py`):
- shock-front speed correct (within 2 bins of exact);
- density & pressure RMS < 5% in the smooth regions;
- all wave structure (rarefaction, contact, shock) reproduced.
- velocity RMS ~10%: the well-known **standard-SPH contact-discontinuity
  artifact** (left-of-contact over-density biases the star-region velocity).
  Density-summation SPH cannot hold a sharp contact; this does NOT converge with
  resolution. Artificial thermal conductivity (Price 2008, implemented) tightens
  density; the cure for the velocity is **adaptive smoothing length** — the next
  step, and required for impacts regardless (huge density contrasts).

## What's implemented (`sph.hpp`)

- cubic-spline kernel (3D); cell-list neighbours, periodic in y,z;
- Monaghan artificial viscosity + Price artificial conductivity;
- KDK leapfrog with trapezoidal energy update;
- pluggable EOS (`IdealGas` now; Tillotson next).

## Roadmap to the impact runs

1. **Adaptive h** (h_i ∝ (m/ρ_i)^{1/3}) — fixes the contact artifact, handles
   impact density contrasts.
2. **Tillotson EOS** for iron / basalt / water (+ a von Mises strength model).
3. **Gate 2 — vertical impact π-scaling**: reproduce crater-size scaling for a
   vertical impact before trusting oblique runs.
4. **Oblique runs**: the Phase-1 terminal conditions (v≈7.4 km/s, γ≈1–4°) into
   the three targets; compare gouge/ricochet/melt to the Phase-2a estimates.

## Build & validate

```bash
cd sph && make && ./sod && python3 analyze_sod.py
```
