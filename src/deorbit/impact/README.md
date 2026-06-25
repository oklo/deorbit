# Phase 2 — terminal ultra-oblique impact

The de-orbiting km iron monolith arrives **grazing** (flight-path angle ~1–4°,
v ~ 7.4–11 km/s). This is *not* vertical cratering: the normal velocity is small
(subsonic in rock) while the along-track component is hypersonic, so the body
**ploughs/ricochets** and shock-processes its leading contact — the
decapitation / long-gouge regime (Schultz & Gault 1990).

## 2a — reduced-order model (`grazing_impact.py`, stdlib) — DONE

Order-of-magnitude scaling to bracket outcomes and define hydrocode inputs.
Nominal aerobraked case (v=7.39 km/s, γ=1.16°, 1 km iron, KE ≈ **27 Gt TNT**):

- normal v ≈ 150 m/s (Mach ~0.03–0.1) vs along-track 7.4 km/s → grazing plough.
- leading-contact pressure 55–160 GPa (≫ strength → fluidized contact). At this
  speed it's **below** iron's ~220 GPa shock-melt threshold → the iron largely
  **survives and ricochets** (~74% velocity retained), leaving a melt-lined gouge.
- bulldozer e-fold gouge length ~3.6 km (rock) / 5.2 km (sediment) / 10.5 km (water).
- **Ocean**: punches through the 4 km water column at ~5 km/s to the seafloor →
  large water cavity + tsunami.
- melt is **shock-coupling-limited, not energy-limited** (deposited energy ≈ 14× the
  energy to melt the whole body). Higher-energy variant (11 km/s, 4°, 60 Gt):
  contact reaches ~350 GPa at Cayambe → melts iron at the contact.

These are scaling estimates, grounded in Schultz–Gault ricochet behavior and
Newtonian penetration; validate against the Rio Cuarto elongated-crater field.
NOT a substitute for a hydrocode.

## 2b — 3D hydrocode (the detailed gouge / ricochet / ejecta / tsunami) — TODO

The grazing geometry is inherently 3D, so **iSALE2D (axisymmetric) cannot do it**.
Options, in order of recommendation:

1. **Write a 3D SPH solver** (Tillotson/ANEOS EOS for iron + basalt + water, with
   a strength model) — SPH is the standard tool for oblique/marine impacts, is
   naturally 3D + large-deformation + free-surface, and fits this project's
   "write the code" pattern (numpy/CPU). Largest effort, fullest control.
2. **iSALE3D** — the canonical impact hydrocode, but **gated** (not in the public
   Dellen release; needs a developer account from Imperial) and a hard Fortran
   build. Would require Greg's access negotiation.
3. **GPU SPH (miluphcuda)** — impact-oriented, open, but needs CUDA (not available
   on this macOS/arm64 machine).

Recommended: option 1, built incrementally and validated against vertical-impact
crater scaling (π-scaling) before trusting the oblique runs.
