// Milestone 3b-i: validate the GPU strain/stress kernel (Jaumann dS/dt) against
// the CPU dump (force_ref.bin from cpu_dump). Builds the cell list, runs
// strain_stress, diffs dS/dt vs the FP64 CPU oracle.
#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <Foundation/Foundation.hpp>
#include "../common/eos.hpp"
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>

struct GMat { float rho0,A,B,a,b,alpha,beta,u0,uiv,ucv,G,Emod,Y; };
static GMat toG(const Material& m){ return {(float)m.rho0,(float)m.A,(float)m.B,(float)m.a,(float)m.b,
    (float)m.alpha,(float)m.beta,(float)m.u0,(float)m.uiv,(float)m.ucv,(float)m.G,(float)m.Emod,(float)m.Y}; }

static MTL::Device* gDev; static MTL::CommandQueue* gQ; static MTL::Library* gLib;
static MTL::Buffer* buf(size_t b, const void* s=nullptr){ return s ? gDev->newBuffer(s,b,MTL::ResourceStorageModeShared)
                                                                   : gDev->newBuffer(b,MTL::ResourceStorageModeShared); }
static MTL::ComputePipelineState* pso(const char* nm){ NS::Error* e=nullptr;
    auto p=gDev->newComputePipelineState(gLib->newFunction(NS::String::string(nm,NS::UTF8StringEncoding)),&e);
    if(!p){printf("pso %s: %s\n",nm,e?e->localizedDescription()->utf8String():"?");exit(1);} return p; }
static void run(MTL::ComputePipelineState* p, uint32_t n, std::vector<MTL::Buffer*> bs,
                std::vector<std::pair<const void*,size_t>> by){
    auto cb=gQ->commandBuffer(); auto e=cb->computeCommandEncoder(); e->setComputePipelineState(p);
    for(size_t i=0;i<bs.size();i++) e->setBuffer(bs[i],0,i);
    for(size_t i=0;i<by.size();i++) e->setBytes(by[i].first,by[i].second,bs.size()+i);
    e->dispatchThreads(MTL::Size(n,1,1),MTL::Size(256,1,1)); e->endEncoding(); cb->commit(); cb->waitUntilCompleted();
}

int main(){
    FILE* f=std::fopen("force_ref.bin","rb");
    if(!f){printf("run ./cpu_dump first\n");return 1;}
    int32_t n; std::fread(&n,4,1,f); double hdr[6]; std::fread(hdr,8,6,f);
    std::vector<double> rec(29*(size_t)n); std::fread(rec.data(),8,29*(size_t)n,f); std::fclose(f);

    std::vector<float> pos(3*n),vel(3*n),rho(n),hh(n),mass(n),Sin(6*n),dSdt_cpu(6*n),
                       u(n),D(n),acc_cpu(3*n),dudt_cpu(n);
    std::vector<int> mat(n);
    for(int i=0;i<n;i++){ double* r=&rec[29*i];
        for(int k=0;k<3;k++){pos[3*i+k]=r[k]; vel[3*i+k]=r[3+k]; acc_cpu[3*i+k]=r[19+k];}
        u[i]=r[7]; mass[i]=r[6]; rho[i]=r[8]; hh[i]=r[9]; mat[i]=(int)r[10]; D[i]=r[17]; dudt_cpu[i]=r[22];
        for(int k=0;k<6;k++){ Sin[6*i+k]=r[11+k]; dSdt_cpu[6*i+k]=r[23+k]; } }
    float eta=(float)hdr[0],alpha=(float)hdr[1],beta=(float)hdr[2],eps=(float)hdr[3],eps_as=(float)hdr[5];

    gDev=MTL::CreateSystemDefaultDevice(); gQ=gDev->newCommandQueue();
    printf("GPU: %s   n=%d\n", gDev->name()->utf8String(), n);
    NS::Error* e=nullptr; gLib=gDev->newLibrary(NS::String::string("sph_force.metallib",NS::UTF8StringEncoding),&e);
    if(!gLib){printf("lib: %s\n",e?e->localizedDescription()->utf8String():"?");return 1;}
    GMat mats[2]={toG(Material::basalt()),toG(Material::iron())};

    // cell list (cell size 2*max hh)
    float lo[3]={1e30f,1e30f,1e30f},hi[3]={-1e30f,-1e30f,-1e30f};
    for(int i=0;i<n;i++)for(int a=0;a<3;a++){lo[a]=std::min(lo[a],pos[3*i+a]);hi[a]=std::max(hi[a],pos[3*i+a]);}
    float hmax=0; for(int i=0;i<n;i++) hmax=std::max(hmax,hh[i]);
    float cw=2.0f*hmax; int nc[3]; for(int a=0;a<3;a++)nc[a]=std::max(1,(int)((hi[a]-lo[a])/cw)+1);
    int Nc=nc[0]*nc[1]*nc[2]; uint32_t nn=n;
    auto bpos=buf(3*n*4,pos.data()), bvel=buf(3*n*4,vel.data()), brho=buf(n*4,rho.data()),
         bhh=buf(n*4,hh.data()), bmass=buf(n*4,mass.data()), bS=buf(6*n*4,Sin.data()), bmat=buf(n*4,mat.data());
    auto bhash=buf(n*4), bcount=buf(Nc*4), bslot=buf(n*4), bstart=buf(Nc*4), bsorted=buf(n*4);
    auto bmats=buf(sizeof(mats),mats), bdS=buf(6*n*4), bepsw=buf(n*4);
    memset(bcount->contents(),0,Nc*4);
    run(pso("cell_hash"), n, {bpos,bhash,bcount,bslot},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
    uint32_t* cnt=(uint32_t*)bcount->contents(); uint32_t* st=(uint32_t*)bstart->contents();
    uint32_t acc=0; for(int c=0;c<Nc;c++){st[c]=acc;acc+=cnt[c];}
    run(pso("cell_scatter"), n, {bhash,bslot,bstart,bsorted}, {{&nn,4}});
    run(pso("strain_stress"), n, {bpos,bvel,brho,bhh,bmass,bS,bmat,bsorted,bstart,bcount,bmats,bdS,bepsw},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});

    // diff dS/dt (relative to the dataset scale, since individual components vary in sign/size)
    float* dSg=(float*)bdS->contents();
    double scale=0; for(size_t k=0;k<6*(size_t)n;k++) scale=std::max(scale,(double)std::fabs(dSdt_cpu[k]));
    double maxrel=0,sumsq=0; int cnt2=0;
    for(size_t k=0;k<6*(size_t)n;k++){ double rel=std::fabs(dSg[k]-dSdt_cpu[k])/scale; maxrel=std::max(maxrel,rel); sumsq+=rel*rel; cnt2++; }
    printf("dS/dt scale=%.3e   RMS rel err=%.2e   max rel err=%.2e   %s\n", scale, std::sqrt(sumsq/cnt2), maxrel,
           maxrel<2e-3 ? "STRAIN/STRESS GATE: PASS" : "STRAIN/STRESS GATE: investigate");

    // ---- forces: eos_full -> forces, diff acc + dudt ----
    auto bu=buf(n*4,u.data()), bD=buf(n*4,D.data());
    auto bPf=buf(n*4), bcsig=buf(n*4), bdscale=buf(n*4), bRart=buf(n*4), baccb=buf(3*n*4), bdudt=buf(n*4);
    int32_t damon=1;
    run(pso("eos_full"), n, {brho,bu,bmat,bD,bPf,bcsig,bdscale,bRart,bmats},
        {{&eps_as,4},{&damon,4},{&nn,4}});
    run(pso("forces"), n, {bpos,bvel,brho,bhh,bmass,bS,bPf,bcsig,bdscale,bRart,bepsw,bsorted,bstart,bcount,baccb,bdudt},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},
         {&eta,4},{&alpha,4},{&beta,4},{&eps,4},{&eps_as,4},{&nn,4}});
    float* accg=(float*)baccb->contents(); float* dudtg=(float*)bdudt->contents();
    double ascale=0,dscale_=0;
    for(int i=0;i<3*n;i++) ascale=std::max(ascale,(double)std::fabs(acc_cpu[i]));
    for(int i=0;i<n;i++) dscale_=std::max(dscale_,(double)std::fabs(dudt_cpu[i]));
    double aM=0,aS=0,dM=0,dS_=0;
    for(int i=0;i<3*n;i++){ double r=std::fabs(accg[i]-acc_cpu[i])/ascale; aM=std::max(aM,r); aS+=r*r; }
    for(int i=0;i<n;i++){ double r=std::fabs(dudtg[i]-dudt_cpu[i])/dscale_; dM=std::max(dM,r); dS_+=r*r; }
    printf("acc   scale=%.3e   RMS rel err=%.2e   max rel err=%.2e\n", ascale, std::sqrt(aS/(3*n)), aM);
    printf("dudt  scale=%.3e   RMS rel err=%.2e   max rel err=%.2e\n", dscale_, std::sqrt(dS_/n), dM);
    printf("%s\n", (aM<3e-3 && dM<3e-3) ? "FORCES GATE: PASS (FP32)" : "FORCES GATE: investigate");
    return 0;
}
