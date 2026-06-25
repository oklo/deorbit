# PICKUP — Route 1: well-balanced free surface under gravity (euler code)

Handoff for a fresh session. The euler grid code (`deorbit/euler/`) is a validation-gated
cross-check of the SPH for impact cratering. **M1–M5 are DONE and gated** (see
`euler/README.md` + the `deorbit-euler-code` memory). This file is ONE focused task that
unblocks M6 (crater collapse) and M7 (the Orcus run).

## The one blocker
A gravity-loaded basalt substrate with a free surface (basalt over a low-density ambient,
ρ≈0.27) is **not hydrostatically well-balanced**, so the free surface erodes. Diagnosed
precisely:
- The minmod limiter **under-resolves the hydrostatic pressure gradient in the surface cell**
  (it flattens the slope at the sharp basalt/ambient interface → the bottom-face reconstructed
  P comes out ~¾ of lithostatic → a residual net downward force → the surface erodes mass
  top-down, every step).
- Bulk imbalance is tiny (~370 m drift over 200 s), so the problem is **specifically the
  interface/surface**, not the bulk.
- **Route 2 (fixed floor + relaxation damping + void-CFL) was tried and is INSUFFICIENT**:
  damping bounds the *velocity* (max|v|→0) but mass erosion is a *flux* effect, so the
  substrate still rarefies catastrophically (basalt cells 2400 → 120). Don't redo route 2.

## The fix: Route 1 — hydrostatic reconstruction (Audusse et al. 2004)
Reconstruct the pressure **deviation from a hydrostatic reference**, so a hydrostatic state
has identically zero net flux (gravity source ↔ pressure-gradient cancel discretely) and the
limiter can't break the balance.

**Suggested pragmatic implementation** (in the z-sweep of `Lop`, both `hydro_cpu.cpp` and
`hydro.metal`):
1. Build a fixed hydrostatic reference `P0(z)` for the substrate (the lithostatic profile from
   the IC; for basalt `P0 = ρ0 g (z_surf − z)`, ambient `P0=0`). Pass it in (a per-k array, or
   recompute from a stored reference density).
2. In the z-direction reconstruction, reconstruct `δP = P − P0(z)` (and the corresponding state)
   instead of `P`. For a near-hydrostatic column δP≈0 → smooth → the limiter doesn't flatten the
   balance. The shock/crater perturbation lives in δP (large), so the limiter still works there.
3. At the z-face, `P_face = P0(z_face) + δP_reconstructed`. Discretize the gravity source
   consistently with the reference (the well-balanced source term) so a hydrostatic state has
   exactly zero net `(flux divergence + gravity)`.
4. The sharp basalt/ambient interface is the analogue of Audusse's discontinuous bathymetry —
   handle the P0 jump at the interface with the standard hydrostatic-reconstruction interface
   values (min/clip of the reconstructed depths/pressures from each side).

Only the z-direction (gravity) needs this; x,y sweeps are unchanged. Keep everything else
(MUSCL+HLLC+RK2, strength, damage, AF) as-is.

## Validation gates (already coded — must flip CHECK → PASS)
- `./hydro_cpu substrate` — relaxed gravity-loaded basalt substrate must **settle stable, no
  rarefaction** (currently: basalt 2400→120 = FAIL). Target: basalt cell count stays ~constant,
  max|v| small. THIS IS THE ROUTE-1 GATE.
- `./hydro_cpu collapse` — AF-on basalt step should slump (<0.5·h0), AF-off should hold
  (>0.7·h0). Needs the stable substrate first; this is the M6 acoustic-fluidization demo.
- **Regression (must stay green): `for m in sod sedov surface shear yield tensile freefall
  atmos alimpact; do ./hydro_cpu $m; done`** and the GPU equivalents — every gate prints
  PASS/CHECK. The route-1 change must keep all M1–M5 gates PASS.

## Then
- M6: tune AF (TDEC, ETA_AF) so the fluidized transient crater collapses to the right
  depth/diameter (calibration-dependent).
- M7: build the Orcus IC (basalt half-space + ~50–70 km basalt impactor, ~10 km/s, ~7–10°,
  Mars g=3.71) on the stable substrate; run; produce the 3-way figure MOLA | SPH | euler.
  (impactor color tag: add a passive advected scalar to distinguish projectile material.)

## Build / run
```
cd euler
clang++ -std=c++17 -O2 hydro_cpu.cpp -o hydro_cpu                                  # CPU oracle
xcrun -sdk macosx metal -O3 -ffast-math -c hydro.metal -o hydro.air && \
  xcrun -sdk macosx metallib hydro.air -o hydro.metallib                          # GPU lib
clang++ -std=c++17 -O2 -I../gpu/metal-cpp hydro_gpu.cpp \
  -framework Metal -framework Foundation -framework QuartzCore -o hydro_gpu        # GPU host
```
Develop route 1 on the CPU (`hydro_cpu.cpp`) first (it's the FP64 oracle), gate with
`substrate`, then port to `hydro.metal`/`hydro_gpu.cpp` (GPU==CPU), then M6/M7.

## Key files / globals
- `hydro_cpu.cpp` (oracle), `hydro.metal` + `hydro_gpu.cpp` (FP32 GPU, metal-cpp in `../gpu/`).
- Globals: `GZ` (gravity), `TDEC`/`ETA_AF` (acoustic fluidization), `RHO_CFL` (void-CFL),
  `DAMP` (relaxation). Grid fields incl. `af` (fluidization), `D` (damage), `S*` (stress).
- EOS shared from `../sph/eos.hpp` (Material::basalt/aluminum/ideal).
- Commit as user.name='oklo' user.email='oklo@mac.com'; end messages with the Co-Authored-By line.
