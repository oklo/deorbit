// Density milestone / go-no-go: build a jittered lattice of particles, compute
// SPH density on the GPU (FP32) and on the CPU (FP64, the oracle), and report the
// relative error + the GPU vs CPU throughput. Retires the three big risks at
// once: toolchain, FP32 fidelity, a working neighbour kernel.
#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <Foundation/Foundation.hpp>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

static double kW64(double r, double h) {              // FP64 reference, matches sph.hpp
    double q = r / h, s = 1.0 / (M_PI * h * h * h);
    if (q < 1.0) return s * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    if (q < 2.0) { double t = 2.0 - q; return s * 0.25 * t * t * t; }
    return 0.0;
}
static double now() { return std::chrono::duration<double>(
    std::chrono::steady_clock::now().time_since_epoch()).count(); }

int main() {
    // jittered cubic lattice (interior particles -> rho ~ rho0, a real density test)
    const int L = 40; const int n = L * L * L;        // 64,000 particles
    const double dx = 1.0, h = 1.3 * dx, rho0 = 1000.0, m = rho0 * dx * dx * dx;
    std::vector<float> pos(3 * n), mass(n, (float)m);
    std::mt19937 rng(1); std::uniform_real_distribution<double> J(-0.2 * dx, 0.2 * dx);
    for (int a = 0, k = 0; a < L; a++) for (int b = 0; b < L; b++) for (int c = 0; c < L; c++, k++) {
        pos[3 * k] = a * dx + J(rng); pos[3 * k + 1] = b * dx + J(rng); pos[3 * k + 2] = c * dx + J(rng);
    }

    // ---- CPU FP64 reference ----
    std::vector<double> rho_cpu(n);
    double t0 = now();
    for (int i = 0; i < n; i++) {
        double rx = pos[3*i], ry = pos[3*i+1], rz = pos[3*i+2], s = 0;
        for (int j = 0; j < n; j++) {
            double d0 = rx - pos[3*j], d1 = ry - pos[3*j+1], d2 = rz - pos[3*j+2];
            double r = std::sqrt(d0*d0 + d1*d1 + d2*d2);
            if (r < 2.0 * h) s += mass[j] * kW64(r, h);
        }
        rho_cpu[i] = s;
    }
    double cpu_s = now() - t0;

    // ---- GPU FP32 ----
    MTL::Device* dev = MTL::CreateSystemDefaultDevice();
    printf("GPU: %s\n", dev->name()->utf8String());
    NS::Error* err = nullptr;
    auto lib = dev->newLibrary(NS::String::string("sph_density.metallib", NS::UTF8StringEncoding), &err);
    if (!lib) { printf("lib load failed: %s\n", err ? err->localizedDescription()->utf8String() : "?"); return 1; }
    auto fn = lib->newFunction(NS::String::string("density_bruteforce", NS::UTF8StringEncoding));
    auto pso = dev->newComputePipelineState(fn, &err);
    auto queue = dev->newCommandQueue();
    auto bpos = dev->newBuffer(pos.data(),  3 * n * sizeof(float), MTL::ResourceStorageModeShared);
    auto bmas = dev->newBuffer(mass.data(),     n * sizeof(float), MTL::ResourceStorageModeShared);
    auto brho = dev->newBuffer(            n * sizeof(float), MTL::ResourceStorageModeShared);
    uint32_t nn = n; float hh = (float)h;

    double tg = now();
    auto cmd = queue->commandBuffer();
    auto enc = cmd->computeCommandEncoder();
    enc->setComputePipelineState(pso);
    enc->setBuffer(bpos, 0, 0); enc->setBuffer(bmas, 0, 1); enc->setBuffer(brho, 0, 2);
    enc->setBytes(&nn, sizeof(nn), 3); enc->setBytes(&hh, sizeof(hh), 4);
    enc->dispatchThreads(MTL::Size(n, 1, 1), MTL::Size(256, 1, 1));
    enc->endEncoding(); cmd->commit(); cmd->waitUntilCompleted();
    double gpu_s = now() - tg;

    // ---- compare (interior particles only; surfaces have legitimately low rho) ----
    float* rho_gpu = (float*)brho->contents();
    double maxrel = 0, sumsq = 0; int cnt = 0;
    for (int i = 0; i < n; i++) {
        if (rho_cpu[i] < 0.5 * rho0) continue;        // skip edges
        double rel = std::fabs(rho_gpu[i] - rho_cpu[i]) / rho_cpu[i];
        maxrel = std::max(maxrel, rel); sumsq += rel * rel; cnt++;
    }
    double rms = std::sqrt(sumsq / cnt);
    printf("particles: %d   interior compared: %d\n", n, cnt);
    printf("FP32-vs-FP64 density:  RMS rel err = %.2e   max rel err = %.2e\n", rms, maxrel);
    printf("timing:  CPU %.3f s   GPU %.4f s   -> %.0fx (brute-force O(N^2); cell-list later)\n",
           cpu_s, gpu_s, cpu_s / gpu_s);
    printf("%s\n", (rms < 1e-3) ? "DENSITY GATE: PASS (FP32 within 1e-3 of FP64)"
                                : "DENSITY GATE: investigate (FP32 error > 1e-3)");
    return 0;
}
