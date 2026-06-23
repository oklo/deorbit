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
   - ⬜ 3b forces (pressure + AV + deviatoric + artificial stress) + strain/stress
   - ⬜ 3c Benz-Asphaug damage
4. KDK integrator + adaptive-h iteration (on-GPU step loop)
5. validate full pipeline vs the CPU gates (Sod, strength, pi-scaling, energy)
6. benchmark; then high-resolution / multi-case production runs
