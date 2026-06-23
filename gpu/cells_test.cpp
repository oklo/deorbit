// Milestone 2: GPU counting-sort cell list -> O(N) neighbours. Validates the
// cell-list density against the FP64 CPU oracle (must match the brute-force
// result), then scales to millions of particles that brute-force can't touch.
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
#include <cstring>

static double kW64(double r, double h) {
    double q = r / h, s = 1.0 / (M_PI * h * h * h);
    if (q < 1.0) return s * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    if (q < 2.0) { double t = 2.0 - q; return s * 0.25 * t * t * t; }
    return 0.0;
}
static double now() { return std::chrono::duration<double>(
    std::chrono::steady_clock::now().time_since_epoch()).count(); }

static MTL::Device* gDev; static MTL::CommandQueue* gQ;
static MTL::ComputePipelineState *psoHash, *psoScatter, *psoDens;

static MTL::ComputePipelineState* mkpso(MTL::Library* lib, const char* name) {
    NS::Error* e = nullptr;
    auto fn = lib->newFunction(NS::String::string(name, NS::UTF8StringEncoding));
    auto p = gDev->newComputePipelineState(fn, &e);
    if (!p) { printf("pipeline %s failed: %s\n", name, e ? e->localizedDescription()->utf8String() : "?"); exit(1); }
    return p;
}
static MTL::Buffer* buf(size_t bytes, const void* src = nullptr) {
    return src ? gDev->newBuffer(src, bytes, MTL::ResourceStorageModeShared)
               : gDev->newBuffer(bytes, MTL::ResourceStorageModeShared);
}
static void dispatch1D(MTL::ComputePipelineState* pso, uint32_t n,
                       std::vector<MTL::Buffer*> bufs, std::vector<std::pair<const void*, size_t>> bytes) {
    auto cb = gQ->commandBuffer(); auto enc = cb->computeCommandEncoder();
    enc->setComputePipelineState(pso);
    for (size_t i = 0; i < bufs.size(); i++) enc->setBuffer(bufs[i], 0, i);
    for (size_t i = 0; i < bytes.size(); i++) enc->setBytes(bytes[i].first, bytes[i].second, bufs.size() + i);
    enc->dispatchThreads(MTL::Size(n, 1, 1), MTL::Size(256, 1, 1));
    enc->endEncoding(); cb->commit(); cb->waitUntilCompleted();
}

// GPU cell-list density into rho_out; returns GPU wall time (s)
static double density_cells_gpu(const std::vector<float>& pos, const std::vector<float>& mass,
                                int n, float h, std::vector<float>& rho_out) {
    float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
    for (int i = 0; i < n; i++) for (int a = 0; a < 3; a++) {
        lo[a] = std::min(lo[a], pos[3 * i + a]); hi[a] = std::max(hi[a], pos[3 * i + a]); }
    float cw = 2.0f * h;
    int nc[3]; for (int a = 0; a < 3; a++) nc[a] = std::max(1, (int)((hi[a] - lo[a]) / cw) + 1);
    int Ncells = nc[0] * nc[1] * nc[2];

    auto bpos = buf(3 * n * sizeof(float), pos.data());
    auto bmass = buf(n * sizeof(float), mass.data());
    auto bhash = buf(n * sizeof(int));
    auto bcount = buf(Ncells * sizeof(uint32_t));
    auto bslot = buf(n * sizeof(uint32_t));
    auto bstart = buf(Ncells * sizeof(uint32_t));
    auto bsorted = buf(n * sizeof(uint32_t));
    auto brho = buf(n * sizeof(float));
    memset(bcount->contents(), 0, Ncells * sizeof(uint32_t));
    uint32_t nn = n;

    double t0 = now();
    dispatch1D(psoHash, n, {bpos, bhash, bcount, bslot},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
    // exclusive prefix scan on the shared count buffer (Ncells << N, trivial on CPU)
    uint32_t* cnt = (uint32_t*)bcount->contents(); uint32_t* start = (uint32_t*)bstart->contents();
    uint32_t acc = 0; for (int c = 0; c < Ncells; c++) { start[c] = acc; acc += cnt[c]; }
    dispatch1D(psoScatter, n, {bhash, bslot, bstart, bsorted}, {{&nn,4}});
    dispatch1D(psoDens, n, {bpos, bmass, bsorted, bstart, bcount, brho},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&h,4},{&nn,4}});
    double dt = now() - t0;

    rho_out.resize(n);
    memcpy(rho_out.data(), brho->contents(), n * sizeof(float));
    printf("   [cells: %dx%dx%d=%d cells, %.1f particles/cell]\n", nc[0],nc[1],nc[2],Ncells,(double)n/Ncells);
    return dt;
}

static int make_lattice(int L, std::vector<float>& pos, std::vector<float>& mass, double dx, double rho0) {
    int n = L * L * L; pos.resize(3 * n); mass.assign(n, (float)(rho0 * dx * dx * dx));
    std::mt19937 rng(1); std::uniform_real_distribution<double> J(-0.2 * dx, 0.2 * dx);
    for (int a = 0, k = 0; a < L; a++) for (int b = 0; b < L; b++) for (int c = 0; c < L; c++, k++) {
        pos[3*k]=a*dx+J(rng); pos[3*k+1]=b*dx+J(rng); pos[3*k+2]=c*dx+J(rng); }
    return n;
}

int main() {
    gDev = MTL::CreateSystemDefaultDevice(); gQ = gDev->newCommandQueue();
    printf("GPU: %s\n", gDev->name()->utf8String());
    NS::Error* err = nullptr;
    auto lib = gDev->newLibrary(NS::String::string("sph_cells.metallib", NS::UTF8StringEncoding), &err);
    if (!lib) { printf("lib failed: %s\n", err ? err->localizedDescription()->utf8String() : "?"); return 1; }
    psoHash = mkpso(lib, "cell_hash"); psoScatter = mkpso(lib, "cell_scatter"); psoDens = mkpso(lib, "density_cells");

    // ---- validation vs FP64 oracle (N=64k) ----
    const double dx = 1.0, h = 1.3 * dx, rho0 = 1000.0;
    std::vector<float> pos, mass;
    int n = make_lattice(40, pos, mass, dx, rho0);
    std::vector<double> rho_cpu(n);
    double t0 = now();
    for (int i = 0; i < n; i++) { double rx=pos[3*i],ry=pos[3*i+1],rz=pos[3*i+2],s=0;
        for (int j = 0; j < n; j++) { double d0=rx-pos[3*j],d1=ry-pos[3*j+1],d2=rz-pos[3*j+2];
            double r=std::sqrt(d0*d0+d1*d1+d2*d2); if(r<2*h) s+=mass[j]*kW64(r,h);} rho_cpu[i]=s; }
    double cpu_s = now() - t0;
    std::vector<float> rho_gpu;
    double g64 = density_cells_gpu(pos, mass, n, h, rho_gpu);
    double maxrel=0,sumsq=0; int cnt=0;
    for (int i=0;i<n;i++){ if(rho_cpu[i]<0.5*rho0) continue;
        double rel=std::fabs(rho_gpu[i]-rho_cpu[i])/rho_cpu[i]; maxrel=std::max(maxrel,rel); sumsq+=rel*rel; cnt++; }
    printf("N=%d  cell-list vs FP64 CPU: RMS %.2e  max %.2e   (CPU %.2fs, GPU %.4fs, %.0fx)\n",
           n, std::sqrt(sumsq/cnt), maxrel, cpu_s, g64, cpu_s/g64);
    printf("%s\n", std::sqrt(sumsq/cnt)<1e-3 ? "CELL-LIST GATE: PASS" : "CELL-LIST GATE: investigate");

    // ---- scaling to millions (cell-list only; brute-force can't) ----
    std::vector<float> pos2, mass2;
    int n2 = make_lattice(160, pos2, mass2, dx, rho0);     // 4.096M particles
    std::vector<float> rho2;
    double g2 = density_cells_gpu(pos2, mass2, n2, h, rho2);
    // spot check a geometrically-central (true interior) particle: rho ~ rho0
    int L2 = 160, mid = ((L2/2) * L2 + L2/2) * L2 + L2/2;
    printf("SCALING: N=%d  GPU %.3fs  (%.0f Mparticles/s)  interior rho=%.0f (rho0=%.0f)\n",
           n2, g2, n2/g2/1e6, rho2[mid], rho0);
    return 0;
}
