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

## Vacuum-aware Riemann flux — ALSO DONE (gate `vacuum`, CPU+GPU)
Material expanding into (near-)vacuum is handled by the exact rarefaction-into-vacuum flux sampled
at the face (effective gamma g=rho*c^2/p; exact for ideal gas, stiff-limit for Tillotson). Gated
vs Toro's analytic centred rarefaction (L1 rho=0.0012, u=0.0027, positive; GPU==CPU). Off by
default (RHO_VAC=0); on for collapse. Files: vac_flux + RHO_VAC (CPU); vac_flux + rvac threaded
through hllc/faceflux/lop (GPU).

## Remaining: M6 collapse — now PHYSICS calibration, not numerics
The collapse cliff (vertical x-direction free surface) is now NUMERICALLY STABLE — the void-aware
strength + vacuum-aware flux brought max|v| from ~1e9 down to ~8-20 m/s, no blowup. The `collapse`
gate still CHECKs, but for a PHYSICS reason: a 4 km basalt step on Mars has lithostatic shear
~40 MPa << Y=350 MPa, so it holds elastically (AF-off correctly holds). AF-on fails to slump
because ETA_AF=1e9 Pa.s viscously freezes the fluidized flow (eta*grad v ~ 125 MPa resists it).
NEXT (M6 tuning, CPU oracle first): lower ETA_AF (try 1e6-1e8), set a finite TDEC, until AF-on
slumps (<0.5*h0) while AF-off holds (>0.7*h0). Then add a `collapse` mode to the GPU host and
re-gate GPU==CPU. (Note: the WB z-path `hllc_p` does NOT yet have vacuum-awareness; for M7 ejecta
flying up into vacuum, consider threading RHO_VAC into faceflux_wb too.)

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
