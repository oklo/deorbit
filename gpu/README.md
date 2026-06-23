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
   SPEEDUP (clean, production scale): ~10-11x. Measured at 5.5M (dx=30): GPU
   ~1.88 s/step vs CPU ~20 s/step -> the ~40 h frac run becomes ~3.5-4 h. NB the
   earlier "~19x" (dx=80, 24 vs 462 s) was CONTENTION-INFLATED -- that CPU
   validation shared cores with the then-running CPU frac sim, so the CPU side
   was ~2x slow; the clean uncontended figure is ~10-11x. (Synthetic gpu_run vs
   cpu_bench numbers were also contended + N-dependent; trust the 5.5M figure.)
   Remaining (optional, for ~30-50x): GPU prefix-scan + max-reduction (adaptive
   cells) + batched encoders; ideal-gas EOS for a direct GPU Sod.
