#include <metal_stdlib>
using namespace metal;

inline float kW(float r, float h) {                   // cubic spline, matches sph.hpp
    float q = r / h, s = 1.0f / (M_PI_F * h * h * h);
    if (q < 1.0f) return s * (1.0f - 1.5f * q * q + 0.75f * q * q * q);
    if (q < 2.0f) { float t = 2.0f - q; return s * 0.25f * t * t * t; }
    return 0.0f;
}
inline int cell_of(float3 p, float lox, float loy, float loz, float cw,
                   int ncx, int ncy, int ncz, thread int& bx, thread int& by, thread int& bz) {
    bx = clamp(int((p.x - lox) / cw), 0, ncx - 1);
    by = clamp(int((p.y - loy) / cw), 0, ncy - 1);
    bz = clamp(int((p.z - loz) / cw), 0, ncz - 1);
    return (bx * ncy + by) * ncz + bz;
}

// counting-sort cell list, step 1: hash each particle to a cell, count per cell
// (atomic), and record the particle's slot-within-cell for the scatter.
kernel void cell_hash(device const packed_float3* pos [[buffer(0)]],
                      device int*          hash       [[buffer(1)]],
                      device atomic_uint*  cell_count [[buffer(2)]],
                      device uint*         slot       [[buffer(3)]],
                      constant float& lox [[buffer(4)]], constant float& loy [[buffer(5)]],
                      constant float& loz [[buffer(6)]], constant float& cw  [[buffer(7)]],
                      constant int& ncx [[buffer(8)]], constant int& ncy [[buffer(9)]],
                      constant int& ncz [[buffer(10)]], constant uint& n [[buffer(11)]],
                      uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    int bx, by, bz;
    int c = cell_of(float3(pos[i]), lox, loy, loz, cw, ncx, ncy, ncz, bx, by, bz);
    hash[i] = c;
    slot[i] = atomic_fetch_add_explicit(&cell_count[c], 1u, memory_order_relaxed);
}

// step 3 (after the CPU exclusive-scans cell_count -> cell_start): scatter each
// particle index into its sorted slot.
kernel void cell_scatter(device const int*  hash       [[buffer(0)]],
                         device const uint* slot       [[buffer(1)]],
                         device const uint* cell_start [[buffer(2)]],
                         device uint*       sorted     [[buffer(3)]],
                         constant uint& n [[buffer(4)]],
                         uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    sorted[cell_start[hash[i]] + slot[i]] = i;
}

// step 4: O(N) density via the 27-cell stencil (cell width = 2h covers the support).
kernel void density_cells(device const packed_float3* pos [[buffer(0)]],
                          device const float* mass        [[buffer(1)]],
                          device const uint*  sorted      [[buffer(2)]],
                          device const uint*  cell_start  [[buffer(3)]],
                          device const uint*  cell_count  [[buffer(4)]],
                          device float*       rho         [[buffer(5)]],
                          constant float& lox [[buffer(6)]], constant float& loy [[buffer(7)]],
                          constant float& loz [[buffer(8)]], constant float& cw  [[buffer(9)]],
                          constant int& ncx [[buffer(10)]], constant int& ncy [[buffer(11)]],
                          constant int& ncz [[buffer(12)]], constant float& h [[buffer(13)]],
                          constant uint& n [[buffer(14)]],
                          uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    float3 ri = float3(pos[i]);
    int bx, by, bz;
    cell_of(ri, lox, loy, loz, cw, ncx, ncy, ncz, bx, by, bz);
    float twoh = 2.0f * h, s = 0.0f;
    for (int dx = -1; dx <= 1; dx++) {
        int cx = bx + dx; if (cx < 0 || cx >= ncx) continue;
        for (int dy = -1; dy <= 1; dy++) {
            int cy = by + dy; if (cy < 0 || cy >= ncy) continue;
            for (int dz = -1; dz <= 1; dz++) {
                int cz = bz + dz; if (cz < 0 || cz >= ncz) continue;
                int c = (cx * ncy + cy) * ncz + cz;
                uint start = cell_start[c], cnt = cell_count[c];
                for (uint k = 0; k < cnt; k++) {
                    uint j = sorted[start + k];
                    float3 d = ri - float3(pos[j]);
                    float r = length(d);
                    if (r < twoh) s += mass[j] * kW(r, h);
                }
            }
        }
    }
    rho[i] = s;
}
