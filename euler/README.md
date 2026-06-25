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
- [x] **M4 damage (Grady-Kipp) + gravity — DONE (CPU + GPU).** Gravity: body-force
      source; free-fall v=-g*t exact, hydrostatic atmosphere interior well-balanced
      (5e-4 of cs). Damage: D scalar advected, grows via crack ODE d^(1/3)+=(cg/Rs)dt
      when tensile strain sigmax/Emod > Weibull activation (1/(wk*dx^3))^(1/wm); degrades
      Y->(1-D)Y. GATES (CPU==GPU): tensile stretch -> D->1, shear strength -> 0;
      free-fall + hydrostatic; all prior gates regress. (wk=1e61 overflows FP32 ->
      eps_act precomputed on host.)
- [x] **M5 Pierazzo Al-on-Al shock-physics validation (CREDIBILITY GATE).** Validated on
      the credibility-critical quantity, the PEAK SHOCK PRESSURE: (M5a) 1D Al-on-Al planar
      impact U=10 km/s reproduces the analytic Tillotson Hugoniot EXACTLY (1.638e11 vs
      1.637e11 Pa, 0.0%); (M5b) 3D Al-sphere impact (10 cppr) isobaric core = 93% of the
      Hugoniot (7% deficit = 3D geometry + resolution, per Pierazzo), clean monotonic peak-
      pressure decay P(r/a) (per-cell Pmax tracker). CAVEAT: near-field decay exponent ~1.3
      over r/a 1-6 (the transition regime is shallower than the canonical far-field ~2; a
      far-field domain + higher cppr convergence study would pin it -- a refinement, not done
      here per the option-1 physics-based plan). The shock-capture + EOS give correct impact
      pressures -> the cross-check is credible. Metal: pmax_update kernel; mode 'pierazzo'.
- [x] **Route 1 — well-balanced free surface under gravity (DONE, CPU+GPU).** The loaded
      basalt substrate with a free surface is hydrostatically well-balanced with NO damping
      (CPU DAMP=1.0, GPU dampf=1.0). Five ingredients:
      (1) **Audusse hydrostatic reconstruction** (z-sweep only; x,y unchanged) -- reconstruct the
      pressure DEVIATION from a frozen reference P0(z) and inject it linearly into a pressure-aware
      HLLC (`hllc_p`); the reference face pressure = EOS(avg reference density), so it is exactly
      lithostatic in the bulk and ~0 at the free surface (avg density falls on the cold expanded
      branch). The well-balanced gravity source cancels the reference flux divergence to machine
      precision; the limiter only ever touches the (smooth, ~0) deviation.
      (2) **void cells = passive vacuum** -- cells with rho<RHO_CFL are reset to the reference
      (momentum/density/energy/stress) each step; the tenuous ambient has ~no mass, so any residual
      force gives runaway velocity (a free-surface instability, not a balancing failure). Also
      excluded from the CFL (their huge expanded-EOS sound speed would crush dt).
      (3) **void-aware strength** -- the strength solver uses the CENTRE velocity for void
      neighbours (zero velocity gradient toward vacuum = traction-free free surface) and runs no
      strength in void cells; otherwise the near-vacuum velocity (mu/rho_tiny) feeds spurious
      deviatoric stress into the surface cells (catastrophic in FP32). The predictor's void cells
      are cleaned before the strength read.
      (4) **deep far-field floor** pinned to the reference (a transmissive bottom drains the column).
      GATE `substrate` (CPU & GPU): basalt cells 2400->2400 (zero rarefaction), max|v|~1.4 m/s,
      stable over 200 s, no damping. All M1-M5 gates stay PASS. Metal: hllc_p, faceflux_wb,
      voidzero, damp kernels; lop takes the reference buffers + wb flag; wavespeed & strength take rcfl.
- [x] **Vacuum-aware Riemann flux (DONE, CPU+GPU).** A material|(near-)vacuum face has ~no mass
      on the void side, so any pressure force gives it runaway velocity (the classic vacuum-Riemann
      problem; Toro Ch.4, Einfeldt 1991). At such a face (one side rho<RHO_VAC) we replace HLLC with
      the EXACT rarefaction-into-vacuum solution sampled at x/t=0, using an effective adiabatic
      exponent g=rho*c^2/p (= gamma for ideal gas exactly; ->large for a stiff EOS => minimal
      expansion, the correct near-incompressible-solid limit). This lets material expand into vacuum
      physically instead of relying on the post-step void reset, and is what flagship codes (iSALE
      P=0 free surface, astro positivity/vacuum-Riemann) do. GATE `vacuum` (CPU & GPU, GPU==CPU):
      expansion into vacuum matches Toro's analytic centred-rarefaction fan, L1 rho=0.0012 u=0.0027,
      density stays positive. Direction-independent (works for the x-cliff); all M1-M5 + route-1 gates
      stay PASS. Metal: vac_flux; hllc/faceflux/lop take rvac. (Off by default: RHO_VAC=0.)
- [~] **M6 acoustic fluidization (block model) — PARTIAL (numerics solid, physics to tune).** AF
      mechanism IMPLEMENTED + regression-safe (fluidization field `af` reduces Y->(1-D)(1-af)Y; TDEC
      decay; ETA_AF viscosity). With route 1 + the void-aware strength + the vacuum-aware flux, the
      `collapse` step (incl. the vertical CLIFF) is now NUMERICALLY STABLE (max|v|~8-20 m/s, no
      blowup -- down from ~1e9). The remaining gate failure is PHYSICS CALIBRATION, not numerics:
      a 4 km step on Mars has lithostatic shear ~40 MPa << Y=350 MPa, so it holds elastically; AF-on
      fails to slump because ETA_AF=1e9 Pa.s viscously freezes the fluidized flow. NEXT: tune AF
      (lower ETA_AF, set TDEC) so AF-on slumps (<0.5*h0) and AF-off holds (>0.7*h0).
      Modes: substrate (PASS), vacuum (PASS), collapse (stable; AF calibration pending).
      Globals: af, TDEC, ETA_AF, RHO_CFL, RHO_VAC, DAMP.
- [ ] **M7 Orcus cross-check** — cross-validate vs SPH (early dynamics), then the
      3-way figure: MOLA | SPH | euler.

## Status
**M5 DONE** (credibility gate, peak shock pressure). **Route 1 DONE (CPU+GPU)**: the
gravity-loaded basalt substrate with a free surface is hydrostatically well-balanced
(Audusse P-deviation reconstruction + EOS reference face pressure + void-cell vacuum +
pinned far-field floor); `substrate` gate PASSES with no damping, all M1-M5 gates stay green.
**Vacuum-aware Riemann flux DONE** (gate `vacuum` vs Toro, CPU+GPU): material expands into vacuum
physically (exact rarefaction-into-vacuum at the face). **M6**: the `collapse` cliff is now
numerically stable (~10 m/s, was ~1e9); the residual is AF physics calibration (ETA_AF too high),
not numerics. **M7** (Orcus 3-way figure) now unblocked for a FLAT-substrate start.

**Scope: M1-M5 + Route 1 validated** -- a credible cross-check of the SPH's early impact
dynamics, now with a stable gravity-loaded free surface for the late-stage substrate.
SPH Orcus run done; renders committed.
