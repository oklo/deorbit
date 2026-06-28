# PICKUP — Euler AF calibration, Phase 2b (cylindrical strength) and the road to Chicxulub

Hand-off for a fresh context. Goal of this whole thread: **calibrate acoustic fluidization (AF)
in the euler (finite-difference / Eulerian) code correctly, so we can reproduce the 60° Chicxulub
peak-ring impact of Collins et al. 2020** (Nat. Commun. 11:1480, PDF in ~/Downloads/s41467-020-15269-x.pdf).
AF is the gating physics — it drives the transient-crater collapse that forms the peak ring; without
a calibrated AF model there is no peak ring.

## Why we're doing it this way (the no-shortcuts plan)
Mirror the literature workflow: build the AF block model → calibrate (TDEC, viscosity) on **cheap
vertical craters** vs the depth–diameter curve → apply the *fixed* parameters to the 3D oblique
Chicxulub run. Calibrating against depth–diameter (not against Chicxulub) keeps Chicxulub a genuine
prediction. Vertical craters are axisymmetric, so we run them in a cheap 2D (r,z) slab — which is why
Phase 2 added cylindrical geometry.

## Status (what's DONE and validated, all on master)
- **Phase 1 — shock-activated AF block model (CPU+GPU, GPU==CPU).** Commits c77ed6f (CPU), d9dea58 (GPU).
  - State: `vib` (vibrational velocity) field + `Pmax` (shock-arrival detector). ACTIVATION: a new
    pressure peak (shock arrival) seeds `vib = C_ACT * post-shock speed`. DECAY: `vib *= exp(-dt/TDEC)`.
    RHEOLOGY: `af = p_vib/(p_vib+p_ov)`, `p_vib = rho*c_s*vib`, `p_ov = max(P, P_COH)`. `af` is derived
    each step in `update_af()` (CPU) / the `update_af` kernel (GPU) BEFORE strength reads it.
  - `af` degrades shear strength via `vonmises`: `Y = (1-D)(1-af)*Y0`.
  - GATE `af_activate` (CPU+GPU): a planar basalt shock seeds vib=457.4 m/s behind the front, ~0 ahead;
    fluidized rock (af 1.0) re-solidifies (af 0.03) as vib decays. GPU==CPU exact.
- **Phase 2a — axisymmetric cylindrical (r,z) HYDRO (CPU+GPU, GPU==CPU).** Commits 4a47616 (CPU), e39dad3 (GPU).
  - `AXISYM` flag. x→r (axis at r=0), z→axial, y unused (ny=1). Finite-volume area-weighting: r-sweep
    fluxes weighted by face radius (axis face r=0 → zero area = automatic reflective axis BC, no ghosts)
    + cylindrical pressure source `+p/r_c` on the radial momentum. AXISYM (r) and route-1 WB-gravity (z)
    compose.
  - GATE `sedov_axi` (CPU+GPU): on-axis point blast = 3D spherical Sedov, R=0.750 vs analytic 0.783 (4.2%).

**Key point:** Phase 2a did HYDRO only. Strength is still CARTESIAN. So a vertical crater with strength +
gravity + AF will be WRONG until Phase 2b makes the strength cylindrical.

## PHASE 2b — what to build (cylindrical elastic-plastic strength + AF transport)

Axisymmetric (r,z), no swirl (v_θ=0, ∂/∂θ=0). Map **x→r, y→θ, z→z**. The deviatoric stress has 4
independent components: S_rr↔Sxx, S_θθ↔Syy, S_zz↔Szz, S_rz↔Sxz; and S_rθ↔Sxy=0, S_zθ↔Syz=0 (no swirl).
So the existing 6-component machinery can REPRESENT the cylindrical stress — the changes are geometric.

### 1. Hoop strain rate (the crux)
The θθ strain rate is **geometric**: `ε_θθ = u/r = v_x/r` (radial motion stretches the hoop), NOT
`∂v_y/∂y` (which is 0 for ny=1). In the strength routine, when AXISYM, set the θθ (=yy) strain rate to
`vx(c)/r_c` instead of the y-gradient `Lyy`. This also changes the trace `tr = (ε_rr+ε_zz+ε_θθ)/3`,
hence the deviatoric strain for ALL components.

### 2. Geometric stress-divergence source terms (added to momentum)
Full cylindrical divergence of the symmetric deviatoric stress (axisym, no swirl):
- `(∇·S)_r = ∂S_rr/∂r + ∂S_rz/∂z + (S_rr − S_θθ)/r`   → the existing Cartesian `∂Sxx/∂x + ∂Sxz/∂z`
  PLUS a geometric source `(Sxx − Syy)/r_c` added to d.mu (radial momentum).
- `(∇·S)_z = ∂S_rz/∂r + ∂S_zz/∂z + S_rz/r`             → existing Cartesian `∂Sxz/∂x + ∂Szz/∂z`
  PLUS a geometric source `Sxz/r_c` added to d.mw (axial momentum).
- `(∇·S)_θ = 0`.
(Reason: (1/r)∂(r S_rr)/∂r = ∂S_rr/∂r + S_rr/r, and the −S_θθ/r comes from the curvilinear basis.)

### 3. Energy strength-work cylindrical correction
The deviatoric work term in the energy equation, currently `d.E += div(S·v)` (Cartesian central diffs),
needs the cylindrical geometric term: `(∇·(S·v))_cyl = ∂(S·v)_r/∂r + ∂(S·v)_z/∂z + (S·v)_r/r`. Add the
`(S·v)_r/r_c = (Sxx*vx + Sxz*vz)/r_c` geometric term to d.E when AXISYM. (Alternatively use the
deviatoric work form S:ε_dev, which is frame-independent — but match whatever the code already does.)

### 4. Jaumann rate
In axisym the only spin is `W_rz = ½(∂u/∂z − ∂w/∂r)`; W_rθ=W_zθ=0. θ is a principal direction
(S_rθ=S_zθ=0), so the Jaumann correction for S_θθ (=Syy) should be ≈0. The existing 6-component Jaumann
should reduce correctly IF Sxy=Syz=0 are maintained (no swirl) and Lyy is replaced by the hoop term.
VERIFY this — a stray nonzero W from the y-row could pollute S_θθ.

### 5. AF Newtonian viscosity in cylindrical (+ finish the GPU AF viscosity, deferred from Phase 1)
The AF viscosity is `η·∇²v`. The cylindrical vector Laplacian of the radial velocity has extra terms:
`(∇²v)_r = ∇²u − u/r²` (and the scalar Laplacian itself has the +1/r ∂/∂r term). The CPU has the AF
viscosity (Cartesian Laplacian) in `Lop`; the GPU does NOT (it was deferred in Phase 1). In 2b: (a) add
the cylindrical correction to the CPU AF viscosity, (b) port the AF viscosity to the GPU strength kernel
(needs `af` + `ETA_AF`; the strength kernel already has the velocity stencil).

### 6. vib advection (deferred from Phase 1)
`vib` is a material property and must advect with the flow. Reuse the rc (tracer) advection machinery
— advect `rho*vib` like `rc`, with the cylindrical area-weighting now in place. Add seed/decay in
`update_af` operating on the advected field. (Without this, the AF pattern is frozen in the Eulerian
grid while material flows km during collapse — a real error for quantitative calibration.)

## Code locations
- CPU strength: `hydro_cpu.cpp`, the strength block inside `Lop` (deviatoric stress: Jaumann, the
  `exx/eyy/ezz/tr` strain rates, the AF viscosity `if(g.af[c]>0&&ETA_AF>0)`, and the
  `d.mu += dSxx_x+dSxy_y+dSxz_z` stress-divergence-to-momentum). `vonmises` applies the (1-D)(1-af) cap.
  `update_af` is just above `step_rk2`.
- GPU strength: `hydro.metal`, `kernel void strength` (has the VX/VY/VZ stencil, writes dS + div-S to
  momentum). `kernel void vonmises` (reads `af`). `kernel void update_af`. Host wiring +
  buffers/dispatch in `hydro_gpu.cpp` (search `Pstr`, `Pvm`, `Pupaf`, `axisym`, `baf`, `bvib`).
- The axisym hydro pattern to copy: `hydro_cpu.cpp` the `bool axi=(AXISYM&&dir==0)` block in the sweep;
  GPU `hydro.metal` lop `if(axisym!=0 && dir==0)` area-weighting + `+p/r_c` source.

## Validation gates for Phase 2b
1. **`lame` (recommended primary)** — thick-walled cylinder under internal pressure, elastostatic Lamé
   solution: `S_rr(r) = A − B/r²`, `S_θθ(r) = A + B/r²` with A,B from the BCs. Run to elastic
   equilibrium (low velocity), compare the radial profiles of S_rr and S_θθ. This directly validates the
   hoop stress + geometric stress divergence. (Pure elastic, no AF.) Tolerance ~few %.
2. **Cartesian regression** — `AXISYM=false` must keep ALL existing gates bit-identical (the geometric
   terms are guarded by AXISYM). Run `./gates.sh` (CPU+GPU).
3. **GPU==CPU** on `lame` and on a small axisym strength run.
4. (Stretch) a vertical impact crater with gravity+strength+AF that runs STABLY and produces a transient
   crater → collapse — qualitative, sets up Phase 3.

## Pitfalls / watch-list
- Divide-by-r at the axis: r_c = (i+0.5)*dx ≠ 0 for cell centers, so the geometric SOURCES are finite.
  But guard anyway.
- Keep Sxy=Syz=0 in axisym (no swirl). If the Jaumann/advection leaks nonzero into them, the cylindrical
  mapping breaks.
- The strength block uses CENTRAL differences (not the FV flux form), so geometric terms are added as
  SOURCES (not area-weighting). Don't double-count with the hydro area-weighting (that's only in the
  pressure/advective flux, a different operator).
- FP32 on GPU: expect tiny reorder noise (as seen in the substrate gate, 1.448→2.845 m/s, still <gate).
  Judge by the gate tolerance, not bit-equality, for sensitive near-equilibrium runs.
- vib advection + cylindrical area-weighting must use the SAME face-radius weighting as the hydro fluxes.

## After Phase 2b → Phase 3 (calibration) → Chicxulub
- **Phase 3:** run vertical impacts at a few sizes in axisym; tune (TDEC, ETA_AF) — in dimensionless,
  impactor-scaled form — to match the **Wünnemann & Ivanov 2003** depth–diameter curve (the published
  numerical benchmark; "Numerical modelling of the impact crater depth-diameter dependence in an
  acoustically fluidized target", Planet. Space Sci. 51:831). Validate on an independent peak-ring crater.
- **Chicxulub 60° target (Collins et al. 2020, Fig. 2):** 17 km granite impactor, ρ=2630 kg/m³,
  **12 km/s, 60° to horizontal**, NE→SW. Two-layer target: granite crust (33 km) over dunite mantle,
  ANEOS EOS. 500 m cells (~17 cells/impactor radius), run to T=300 s. Diagnostics: transient crater →
  collapse → peak ring; peak-ring/mantle-uplift offsets vs crater centre (their Fig. 5); melt (>60 GPa)
  + peak-shock-pressure tracer fields. **We use Tillotson basalt (crust proxy); no ANEOS granite/dunite
  yet** — a separate EOS gap to flag (a dunite-like denser Tillotson mantle would be the minimal add).
  This is a ~200 km, ~5-minute, 3D, gravity+AF run — the big one; scope compute carefully.

## Build / test
- Build CPU: `clang++ -std=c++17 -O2 hydro_cpu.cpp -o hydro_cpu`
- Build GPU: rebuild metallib then host:
  `xcrun -sdk macosx metal -O3 -ffast-math -fmodules-cache-path=.clang-module-cache -c hydro.metal -o hydro.air && xcrun -sdk macosx metallib hydro.air -o hydro.metallib`
  `clang++ -std=c++17 -O2 -I../common/metal-cpp hydro_gpu.cpp -framework Metal -framework Foundation -framework QuartzCore -o hydro_gpu`
- All gates: `./gates.sh` (or `./gates.sh quick` to skip the slow CPU sedov). Add `lame` to the CPU+GPU
  lists in gates.sh when built.

## Conventions (project)
- Commit as `user.name='oklo' user.email='oklo@mac.com'`; end messages with
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
- The remote `origin/master` gets a daily cloud-routine corridor-progress commit — `git fetch && rebase
  origin/master` before pushing (keep history linear; no force-push).
- CPU oracle first → gate → GPU port (GPU==CPU), the established rhythm.
- Don't delete the ~130 GB of run outputs; leave `docs/codex_gpt55_review.md` untracked.
