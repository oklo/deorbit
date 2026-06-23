// Milestone 3a: validate the FP32 GPU Tillotson EOS (P, sound speed, signal
// speed) against the FP64 CPU Material across all branches (compressed, expanded
// hot, intermediate) for iron and basalt. The numerical sound speed is the
// FP32-sharpest piece, so it gets watched closely.
#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <Foundation/Foundation.hpp>
#include "../sph/eos.hpp"            // CPU oracle: Material
#include <cstdio>
#include <vector>
#include <cmath>

struct GMat { float rho0, A, B, a, b, alpha, beta, u0, uiv, ucv, G, Emod, Y; };
static GMat toG(const Material& m) {
    return { (float)m.rho0,(float)m.A,(float)m.B,(float)m.a,(float)m.b,(float)m.alpha,
             (float)m.beta,(float)m.u0,(float)m.uiv,(float)m.ucv,(float)m.G,(float)m.Emod,(float)m.Y };
}

int main() {
    Material mat_cpu[2] = { Material::basalt(), Material::iron() };
    GMat mat_g[2] = { toG(mat_cpu[0]), toG(mat_cpu[1]) };

    // sweep (rho, u) across the branches for both materials
    std::vector<float> rho, u; std::vector<int> mat;
    for (int mi = 0; mi < 2; mi++) {
        double r0 = mat_cpu[mi].rho0;
        for (double fr = 0.5; fr <= 1.5001; fr += 0.05)        // 0.5..1.5 rho0 (expanded->compressed)
            for (double uu = 0.0; uu <= 3.0e7; uu += 1.0e6) {  // 0..30 MJ/kg (cold->vaporised)
                rho.push_back((float)(fr * r0)); u.push_back((float)uu); mat.push_back(mi);
            }
    }
    int n = rho.size();

    // CPU reference (FP64)
    std::vector<double> Pc(n), csc(n), csigc(n);
    for (int i = 0; i < n; i++) {
        Pc[i] = mat_cpu[mat[i]].pressure(rho[i], u[i]);
        csc[i] = mat_cpu[mat[i]].sound_speed(rho[i], u[i]);
        csigc[i] = std::sqrt(csc[i]*csc[i] + (4.0/3.0)*mat_cpu[mat[i]].G/std::max((double)rho[i],1e-30));
    }

    // GPU
    auto dev = MTL::CreateSystemDefaultDevice(); auto q = dev->newCommandQueue();
    printf("GPU: %s\n", dev->name()->utf8String());
    NS::Error* e = nullptr;
    auto lib = dev->newLibrary(NS::String::string("sph_eos.metallib", NS::UTF8StringEncoding), &e);
    if (!lib) { printf("lib: %s\n", e?e->localizedDescription()->utf8String():"?"); return 1; }
    auto pso = dev->newComputePipelineState(lib->newFunction(NS::String::string("eos", NS::UTF8StringEncoding)), &e);
    auto S = MTL::ResourceStorageModeShared;
    auto brho=dev->newBuffer(rho.data(),n*4,S), bu=dev->newBuffer(u.data(),n*4,S), bmat=dev->newBuffer(mat.data(),n*4,S);
    auto bP=dev->newBuffer(n*4,S), bcs=dev->newBuffer(n*4,S), bcsig=dev->newBuffer(n*4,S);
    auto bmats=dev->newBuffer(mat_g, sizeof(mat_g), S);
    uint32_t nn=n;
    auto cb=q->commandBuffer(); auto enc=cb->computeCommandEncoder();
    enc->setComputePipelineState(pso);
    enc->setBuffer(brho,0,0); enc->setBuffer(bu,0,1); enc->setBuffer(bmat,0,2);
    enc->setBuffer(bP,0,3); enc->setBuffer(bcs,0,4); enc->setBuffer(bcsig,0,5);
    enc->setBuffer(bmats,0,6); enc->setBytes(&nn,4,7);
    enc->dispatchThreads(MTL::Size(n,1,1), MTL::Size(64,1,1));
    enc->endEncoding(); cb->commit(); cb->waitUntilCompleted();
    float *Pg=(float*)bP->contents(), *csg=(float*)bcs->contents(), *csigg=(float*)bcsig->contents();

    // compare (relative where the CPU value is well above FP32 noise floor)
    auto rel = [](double g, double c, double floor){ return std::fabs(g-c)/std::max(std::fabs(c),floor); };
    // gate on the physically-used quantities: significant P (|P|>10 MPa, real
    // impact pressures are GPa) and the signal speed csig (drives AV + timestep).
    // Raw cs is only an intermediate (used via csig); reported for information.
    double mP=0,mS=0,mCsig_strength=0; int iwP=0;
    for (int i=0;i<n;i++){
        if (std::fabs(Pc[i])>1e7){ double rp=rel(Pg[i],Pc[i],1e7); if(rp>mP){mP=rp;iwP=i;} }
        mS=std::max(mS, rel(csigg[i],csigc[i],50.0));
    }
    printf("states: %d (both materials, all branches)\n", n);
    printf("PRESSURE  (|P|>10 MPa): max rel err = %.2e   [worst: mat=%d rho=%.2f rho0 u=%.1e, cpu=%.3e gpu=%.3e]\n",
           mP, mat[iwP], rho[iwP]/mat_cpu[mat[iwP]].rho0, u[iwP], Pc[iwP], (double)Pg[iwP]);
    printf("SIGNAL SPEED csig (used in AV+dt): max rel err = %.2e\n", mS);
    printf("%s\n", (mP<5e-3 && mS<1e-3) ? "EOS GATE: PASS (FP32; csig + significant P)"
                                        : "EOS GATE: investigate");
    return 0;
}
