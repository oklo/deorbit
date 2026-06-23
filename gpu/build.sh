#!/bin/bash
# Build the Metal GPU SPH pieces. The .metal -> .metallib steps need the Metal
# Toolchain component (xcodebuild -downloadComponent MetalToolchain); the C++
# hosts only need Metal.framework (already in the OS) + the vendored metal-cpp.
set -e
cd "$(dirname "$0")"
SDK=macosx
FRAMEWORKS="-framework Metal -framework Foundation -framework QuartzCore"

build_kernel() {  # $1 = name (without .metal)
    xcrun -sdk $SDK metal -O3 -ffast-math -c "$1.metal" -o "$1.air"
    xcrun -sdk $SDK metallib "$1.air" -o "$1.metallib"
}
build_host() {    # $1 = name (without .cpp)
    clang++ -std=c++17 -O2 -Imetal-cpp "$1.cpp" $FRAMEWORKS -o "$1"
}

MODE="${1:-all}"
if [ "$MODE" = "host" ]; then         # host-only: works before the toolchain lands
    build_host hello
    build_host density_test
    echo "hosts built (kernels still need the Metal toolchain)"
else
    build_kernel hello;        build_host hello
    build_kernel sph_density;  build_host density_test
    build_kernel sph_cells;    build_host cells_test
    build_kernel sph_eos;      build_host eos_gpu_test
    build_kernel sph_force;    build_host force_test;  build_host step_test;  build_host gpu_run;  build_host gpu_ic
    clang++ -std=c++17 -O2 -pthread cpu_dump.cpp  -o cpu_dump     # CPU oracle dumper
    clang++ -std=c++17 -O2 -pthread cpu_bench.cpp -o cpu_bench    # CPU benchmark
    echo "built. run: ./cpu_dump && ./force_test  (also hello/density_test/cells_test/eos_gpu_test)"
fi
