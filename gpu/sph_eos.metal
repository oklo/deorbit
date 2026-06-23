#include <metal_stdlib>
using namespace metal;

// Tillotson EOS in FP32 — must match eos.hpp::till / pressure / sound_speed.
// 11 floats, no vec members -> identical layout on host and device.
struct GMat { float rho0, A, B, a, b, alpha, beta, u0, uiv, ucv, G, Emod, Y; };

inline float till(float rho, float u, constant GMat& m) {
    if (rho <= 0.0f) return 0.0f;
    float eta = rho / m.rho0, mu = eta - 1.0f;
    float w0 = u / (m.u0 * eta * eta) + 1.0f;
    float Pc = (m.a + m.b / w0) * rho * u + m.A * mu + m.B * mu * mu;     // condensed
    if (rho >= m.rho0 || u <= m.uiv) return Pc;
    float z = m.rho0 / rho - 1.0f;
    float Pe = m.a * rho * u + (m.b * rho * u / w0 + m.A * mu * exp(-m.beta * z)) * exp(-m.alpha * z * z);
    if (u >= m.ucv) return Pe;                                            // expanded hot
    return ((u - m.uiv) * Pe + (m.ucv - u) * Pc) / (m.ucv - m.uiv);       // intermediate
}
// ANALYTIC derivatives (FP32-safe: the finite-difference dP/du suffers
// catastrophic cancellation when the density term dwarfs the energy term).
// Condensed branch P and (dP/drho, dP/du):
inline void cond(float rho, float u, constant GMat& m, thread float& P, thread float& dPdr, thread float& dPdu) {
    float eta = rho / m.rho0, mu = eta - 1.0f;
    float w = 1.0f + u / (m.u0 * eta * eta);
    float aob = m.a + m.b / w;
    P = aob * rho * u + m.A * mu + m.B * mu * mu;
    dPdu = aob * rho - m.b * rho * u / (w * w * m.u0 * eta * eta);
    float dwdr = -2.0f * u / (m.u0 * eta * eta * eta * m.rho0);
    dPdr = aob * u + rho * u * (-m.b / (w * w)) * dwdr + m.A / m.rho0 + 2.0f * m.B * mu / m.rho0;
}
// Expanded-hot branch:
inline void expd(float rho, float u, constant GMat& m, thread float& P, thread float& dPdr, thread float& dPdu) {
    float eta = rho / m.rho0, mu = eta - 1.0f;
    float w = 1.0f + u / (m.u0 * eta * eta);
    float z = m.rho0 / rho - 1.0f, dzdr = -m.rho0 / (rho * rho);
    float E1 = exp(-m.alpha * z * z), E2 = exp(-m.beta * z);
    float G = m.b * rho * u / w + m.A * mu * E2;
    P = m.a * rho * u + G * E1;
    float dwdr = -2.0f * u / (m.u0 * eta * eta * eta * m.rho0);
    float dGu = m.b * rho / w - m.b * rho * u / (w * w * m.u0 * eta * eta);
    dPdu = m.a * rho + dGu * E1;
    float dbruw_dr = m.b * u / w + m.b * rho * u * (-1.0f / (w * w)) * dwdr;
    float dGr = dbruw_dr + m.A * E2 / m.rho0 + m.A * mu * (E2 * (-m.beta) * dzdr);
    float dE1dr = E1 * (-m.alpha * 2.0f * z) * dzdr;
    dPdr = m.a * u + dGr * E1 + G * dE1dr;
}
inline float ssound(float rho, float u, constant GMat& m) {
    if (rho <= 0.0f) return 1e-3f;
    float P, dPdr, dPdu;
    if (rho >= m.rho0 || u <= m.uiv) { cond(rho, u, m, P, dPdr, dPdu); }
    else if (u >= m.ucv) { expd(rho, u, m, P, dPdr, dPdu); }
    else {                                              // intermediate: blend both
        float Pc, dPcr, dPcu, Pe, dPer, dPeu;
        cond(rho, u, m, Pc, dPcr, dPcu); expd(rho, u, m, Pe, dPer, dPeu);
        float inv = 1.0f / (m.ucv - m.uiv);
        P    = ((u - m.uiv) * Pe + (m.ucv - u) * Pc) * inv;
        dPdu = ((Pe + (u - m.uiv) * dPeu) + (-Pc + (m.ucv - u) * dPcu)) * inv;
        dPdr = ((u - m.uiv) * dPer + (m.ucv - u) * dPcr) * inv;
    }
    return sqrt(max(dPdr + P / (rho * rho) * dPdu, 1e-6f));
}

kernel void eos(device const float* rho  [[buffer(0)]],
                device const float* u    [[buffer(1)]],
                device const int*   mat  [[buffer(2)]],
                device float* P     [[buffer(3)]],
                device float* cs    [[buffer(4)]],
                device float* csig  [[buffer(5)]],
                constant GMat* mats [[buffer(6)]],
                constant uint& n    [[buffer(7)]],
                uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    constant GMat& m = mats[mat[i]];
    float p = till(rho[i], u[i], m);
    float c = ssound(rho[i], u[i], m);
    P[i] = p; cs[i] = c;
    csig[i] = sqrt(c * c + (4.0f / 3.0f) * m.G / max(rho[i], 1e-30f));
}
