# GPU (Metal) port of the SPH solver

Porting the neighbour loops to the M4 Max's 32-core GPU. The CPU code (`../sph/`)
+ its validation gates are the **reference oracle**: every GPU kernel is diffed
against the FP64 CPU result before it is trusted.

## Setup (one-time)
1. **Metal Toolchain** (for `.metal` -> `.metallib`):
   `xcodebuild -downloadComponent MetalToolchain`
2. **metal-cpp** (C++ host bindings, header-only, redistributable):
   download from https://developer.apple.com/metal/cpp/ and unzip to `gpu/metal-cpp/`
   (we used `metal-cpp_macOS15_iOS18` — forward-compatible).

## Build & run
```
./build.sh          # compiles kernels (.metallib) + hosts
./build.sh host     # hosts only (works before the toolchain is installed)
./hello             # toolchain end-to-end check
./density_test      # FP32-GPU vs FP64-CPU density gate + throughput
```

## Status
- **Milestone 1 (density) — PASSED.** FP32 vs FP64 RMS rel err **1.3e-7** (gate 1e-3);
  brute-force GPU **156x** a single CPU thread. Precision is a non-issue for the
  kernel math; positions over the ~10 km domain will use relative coords.
- **Milestone 2 (cell list) — PASSED.** Counting-sort cell list (hash + atomic
  count + CPU exclusive scan over the shared buffer + scatter) -> O(N) neighbours.
  Cell-list density matches FP64 to **1.3e-7**; **4.1M particles in 40 ms
  (~100 M-particles/s)**, interior rho ~ rho0. The GPU sort risk is retired.
- **Milestone 3a (Tillotson EOS) — PASSED.** FP32 P (|P|>10 MPa) vs FP64 **7e-6**,
  signal speed csig (AV+dt) **1.1e-4**. NB the finite-difference sound speed is
  FP32-broken (catastrophic cancellation when the density term dwarfs the energy
  term) -> switched to **ANALYTIC derivatives** (all 3 branches), which are exact
  + faster. TODO: fold the analytic sound speed back into the CPU eos.hpp so both
  use one method (re-pass strength/pi gates; Sod uses ideal-gas cs, unaffected).
  DONE: CPU eos.hpp now uses the analytic sound speed (gates re-passed).
- **Milestone 3b (strain/stress + forces) — PASSED.** The full per-step force
  evaluation (Jaumann dS/dt; momentum = -gradP + artificial viscosity +
  deviatoric stress + Monaghan artificial stress; energy incl. deviatoric power;
  damage (1-D) scaling) matches the FP64 CPU oracle in FP32: dS/dt 7e-7, acc
  1.5e-6, dudt 5e-7. Validated via cpu_dump (a sheared multi-material state run
  forward so rho/h/S/D are realistic) + force_test.

## Precision plan (FP32)
- Apple GPUs have no fast FP64 -> kernels are FP32/mixed.
- Positions: domain-relative coords so FP32's ~7 digits sit where the action is.
- Density/energy sums: compensated (Kahan) summation; the energy-conservation
  gate catches drift.

## Roadmap
1. ✅ density (brute-force) — toolchain + FP32 + neighbour kernel
2. ✅ GPU cell-list (counting sort) — O(N) neighbours, 4M particles in 40 ms
3. forces + strain/stress + Tillotson EOS + Benz-Asphaug damage kernels
   - ✅ 3a Tillotson EOS (analytic sound speed)
   - ✅ 3b strain/stress (Jaumann) + forces (P + AV + deviatoric + artificial
        stress + damage scaling): acc 1.5e-6, dudt 5e-7, dS/dt 7e-7 vs CPU
   - ✅ 3c Benz-Asphaug damage growth (eig_max + grow_damage; in the step)
4. ✅ FULL GPU KDK STEP (kick/drift + adaptive-h density loop + EOS + sigmax +
      strain/stress + forces + grow-damage + finish[kick2+trapezoidal u/S+yield]).
      One GPU step vs CPU one-step oracle: pos 1e-7, vel 1e-7, u 8e-8, S 3e-7.
5. ✅ multi-step energy + benchmark (gpu_run.cpp + cpu_bench.cpp).
6. ✅ PRODUCTION DRIVER + full-length FP32 check (gpu_ic.cpp). Loads the real IC
      (run_ic format, incl frozen flags), adaptive CFL dt, full-res snapshots
      (gpu_snap_*), frac (host-seeded Weibull eps_act). On cayambe80.bin (291k,
      dx=80) to t=0.3: GPU matches CPU run_ic2 EXACTLY (x,z,|v|,vz,melt,umax all
      identical to printed precision, SAME 283 steps) -> FP32 position-drift
      concern RETIRED. THE GPU WORKFLOW IS TURNKEY:
        ./gpu_ic ic.bin dx t_end [walltime] [--frac]
   SPEEDUP (clean, 5.5M dx=30 vs CPU ~20 s/step): ~11x baseline -> ~54x after
     optimization (370 ms/step); full arc in item 7. t=2.0 run ~43 min (was ~40 h
     CPU). All opts validated EXACT on cayambe80. (NB the early "~19x" at dx=80
     was contention-inflated; trust the 5.5M numbers.)
7. ✅ OPTIMIZATION (profile-driven; gpu_ic has per-phase timers). Profile showed
   neighbour kernels are 99% (sync/PSO/reductions <1%, dropped). All exact-validated
   on cayambe80 (identical iron metrics to CPU):
   (a) PSO caching; (b) density 4->2 iters (h carried -> converges in 2);
   (c) FUSED strain_stress+forces into one neighbour pass (strain_forces kernel);
   (d) sqrt-avoidance (r^2 cutoff test; tiny -> confirmed memory-bound);
   (e) PARTICLE REORDERING by cell each step (gather_w/gather_u8/set_iota: gather
       carried state into cell-sorted order, sorted=identity -> coalesced
       neighbour reads). This was the big one: 2.55x.
   Cumulative 5.5M: 1809 -> 370 ms/step = ~11x -> ~54x vs CPU. A t=2.0 frac run
   is ~43 min (was ~40 h CPU). NB: diagnostics (energy/report/snapshot) must read
   the REORDERED GPU bmass/bmat, not the host arrays (that bug looked like broken
   physics but was just mismatched metrics).
   Further (optional): threadgroup-memory tiling (~1.3x more); ideal-gas GPU Sod.
