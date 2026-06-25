# PICKUP — Route 1 DONE; next is M6 (collapse) + M7 (Orcus)

**Route 1 (well-balanced free surface under gravity) is COMPLETE and validated on CPU+GPU.**
The gravity-loaded basalt substrate with a free surface stays hydrostatically stable with no
damping. See `euler/README.md` (the Route-1 milestone) for the full description and the
`deorbit-euler-code` memory. This file now just hands off the REMAINING work.

## What route 1 delivered (so you don't redo it)
The `substrate` gate PASSES on both paths (basalt cells 2400->2400, max|v|~1.4 m/s, 200 s, no
damping). All M1-M5 gates stay green. Build/run below. The solution (z-direction only):
1. Audusse hydrostatic reconstruction: reconstruct the pressure DEVIATION from a frozen
   reference P0(z), inject linearly into a pressure-aware HLLC (`hllc_p`/`faceflux_wb`);
   reference face P = EOS(avg ref density) (=> lithostatic in bulk, ~0 at the free surface).
   Well-balanced gravity source cancels the reference flux divergence.
2. void cells (rho<RHO_CFL) = passive vacuum: reset to reference each step (mom/rho/E/stress);
   excluded from the CFL.
3. void-aware strength: centre velocity for void neighbours (traction-free surface), no strength
   in void; predictor void-cleaned before the strength read. (This was the key FP32 fix.)
4. deep far-field floor pinned to the reference.
Globals/flags: REF_R0/REF_P0 (CPU) & bRR0/bRP0+wb (GPU); RHO_CFL; DAMP/dampf (=1, off).

## Remaining: M6 collapse (acoustic fluidization demo)
`./hydro_cpu collapse` and `./hydro_gpu` are still BLOCKED — but on a DIFFERENT issue than
route 1. The collapse test uses a basalt STEP; its vertical CLIFF is a free surface in the
**x-direction**, and route 1 well-balances only **z** (gravity). The cliff/vacuum interface
explodes (max|v|~1e9 in ~50 steps). To unblock M6:
- Generalize the free-surface handling to x,y: the WB hydrostatic reconstruction is z-only and
  correct as-is (gravity is in z), but the **void/free-surface robustness** (void-cell reset +
  void-aware strength) should already help the cliff — verify whether the cliff still explodes
  with the current void-aware code (the collapse mode predates some of it; re-test first).
- If it still explodes, the cliff needs proper vacuum-Riemann handling in the x-sweep (the
  current x-flux has no special free-surface treatment) and/or the void-clean applied so the
  cliff ambient never drives runaway flux.
- Then tune AF (TDEC, ETA_AF) so AF-on slumps (<0.5*h0) and AF-off holds (>0.7*h0).

## Then M7 — Orcus cross-check (now unblocked for a FLAT start)
The flat loaded substrate is stable, so M7 can start. Build the Orcus IC (basalt half-space +
~50-70 km basalt impactor, ~10 km/s, ~7-10 deg, Mars g=3.71) on the substrate; run on GPU;
produce the 3-way figure MOLA | SPH | euler. Add a passive advected scalar to tag projectile
material. (Crater-collapse morphology needs M6's AF tuning.)

## Build / run
```
cd euler
clang++ -std=c++17 -O2 hydro_cpu.cpp -o hydro_cpu                                  # CPU oracle
xcrun -sdk macosx metal -O3 -ffast-math -c hydro.metal -o hydro.air && \
  xcrun -sdk macosx metallib hydro.air -o hydro.metallib                          # GPU lib
clang++ -std=c++17 -O2 -I../gpu/metal-cpp hydro_gpu.cpp \
  -framework Metal -framework Foundation -framework QuartzCore -o hydro_gpu        # GPU host
# gates: ./hydro_cpu substrate ; ./hydro_gpu substrate ; regression: sod surface shear yield
#        tensile freefall atmos alimpact (CPU) / + pierazzo (GPU)
```
Commit as user.name='oklo' user.email='oklo@mac.com'; end messages with the Co-Authored-By line.
