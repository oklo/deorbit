#include <metal_stdlib>
using namespace metal;

// Cubic-spline kernel, FP32 — must match the CPU kW in sph.hpp exactly:
//   q=r/h, s=1/(pi h^3); q<1: s(1-1.5q^2+0.75q^3); q<2: s*0.25(2-q)^3; else 0.
inline float kW(float r, float h) {
    float q = r / h, s = 1.0f / (M_PI_F * h * h * h);
    if (q < 1.0f) return s * (1.0f - 1.5f * q * q + 0.75f * q * q * q);
    if (q < 2.0f) { float t = 2.0f - q; return s * 0.25f * t * t * t; }
    return 0.0f;
}

// Brute-force O(N^2) density (correctness milestone): rho_i = sum_j m_j W(|ri-rj|, h).
// packed_float3 -> tight 12-byte stride so the host can pass a flat [x,y,z,...] array.
kernel void density_bruteforce(device const packed_float3* pos  [[buffer(0)]],
                               device const float*         mass [[buffer(1)]],
                               device float*               rho  [[buffer(2)]],
                               constant uint&  n [[buffer(3)]],
                               constant float& h [[buffer(4)]],
                               uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    float3 ri = float3(pos[i]);
    float twoh = 2.0f * h;
    float s = 0.0f;
    for (uint j = 0; j < n; j++) {
        float3 d = ri - float3(pos[j]);
        float r = length(d);
        if (r < twoh) s += mass[j] * kW(r, h);
    }
    rho[i] = s;
}
