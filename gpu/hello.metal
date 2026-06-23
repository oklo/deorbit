#include <metal_stdlib>
using namespace metal;

// Trivial kernel: verifies the toolchain end-to-end (compile -> metallib ->
// dispatch -> unified-memory readback) before we trust anything fancier.
kernel void vadd(device const float* a [[buffer(0)]],
                 device const float* b [[buffer(1)]],
                 device float*       c [[buffer(2)]],
                 uint i [[thread_position_in_grid]]) {
    c[i] = a[i] + b[i];
}
