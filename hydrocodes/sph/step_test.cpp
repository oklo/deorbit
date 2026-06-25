// Milestone 4: validate a full GPU KDK step (kick/drift, adaptive-h density,
// EOS, sigmax, strain/stress, forces, grow-damage, finish=kick2+trapezoidal
// u/S+yield) against the CPU one-step oracle from cpu_dump (force_ref.bin:
// pre-state + acc/dudt/dSdt, then dt + post-step pos/vel/u/S).
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
static GMat toG(const Material& m){ return {(float)m.rho0,(float)m.A,(float)m.B,(float)m.a,(float)m.b,(float)m.alpha,
    (float)m.beta,(float)m.u0,(float)m.uiv,(float)m.ucv,(float)m.G,(float)m.Emod,(float)m.Y}; }

static MTL::Device* D_; static MTL::CommandQueue* Q_; static MTL::Library* L_;
static MTL::Buffer* buf(size_t b,const void* s=nullptr){ return s?D_->newBuffer(s,b,MTL::ResourceStorageModeShared):D_->newBuffer(b,MTL::ResourceStorageModeShared); }
static MTL::ComputePipelineState* pso(const char* nm){ NS::Error* e=nullptr;
    auto p=D_->newComputePipelineState(L_->newFunction(NS::String::string(nm,NS::UTF8StringEncoding)),&e);
    if(!p){printf("pso %s: %s\n",nm,e?e->localizedDescription()->utf8String():"?");exit(1);} return p; }
static void run(MTL::ComputePipelineState* p,uint32_t n,std::vector<MTL::Buffer*> bs,std::vector<std::pair<const void*,size_t>> by){
    auto cb=Q_->commandBuffer(); auto e=cb->computeCommandEncoder(); e->setComputePipelineState(p);
    for(size_t i=0;i<bs.size();i++) e->setBuffer(bs[i],0,i);
    for(size_t i=0;i<by.size();i++) e->setBytes(by[i].first,by[i].second,bs.size()+i);
    e->dispatchThreads(MTL::Size(n,1,1),MTL::Size(256,1,1)); e->endEncoding(); cb->commit(); cb->waitUntilCompleted();
}

int main(){
    FILE* f=std::fopen("force_ref.bin","rb"); if(!f){printf("run ./cpu_dump first\n");return 1;}
    int32_t n; std::fread(&n,4,1,f); double hdr[6]; std::fread(hdr,8,6,f);
    std::vector<double> rec(29*(size_t)n); std::fread(rec.data(),8,29*(size_t)n,f);
    double dtv; std::fread(&dtv,8,1,f);
    std::vector<double> post(13*(size_t)n); std::fread(post.data(),8,13*(size_t)n,f); std::fclose(f);
    float eta=(float)hdr[0],alpha=(float)hdr[1],beta=(float)hdr[2],eps=(float)hdr[3],eps_as=(float)hdr[5];
    float h_init=eta*30.0f; float dt=(float)dtv; uint32_t nn=n; int damon=1;

    std::vector<float> pos(3*n),vel(3*n),u(n),Sv(6*n),rho(n),hh(n),Dd(n),mass(n),epsact(n),
                       accold(3*n),dudtold(n),dSdtold(6*n); std::vector<int> mat(n); std::vector<unsigned char> fixed(n,0);
    for(int i=0;i<n;i++){ double* r=&rec[29*i];
        for(int k=0;k<3;k++){pos[3*i+k]=r[k];vel[3*i+k]=r[3+k];accold[3*i+k]=r[19+k];}
        mass[i]=r[6];u[i]=r[7];rho[i]=r[8];hh[i]=r[9];mat[i]=(int)r[10];
        for(int k=0;k<6;k++){Sv[6*i+k]=r[11+k];dSdtold[6*i+k]=r[23+k];}
        Dd[i]=r[17];epsact[i]=r[18];dudtold[i]=r[22]; }

    D_=MTL::CreateSystemDefaultDevice(); Q_=D_->newCommandQueue();
    printf("GPU: %s  n=%d  dt=%.3e\n",D_->name()->utf8String(),n,dt);
    NS::Error* e=nullptr; L_=D_->newLibrary(NS::String::string("sph_force.metallib",NS::UTF8StringEncoding),&e);
    if(!L_){printf("lib: %s\n",e?e->localizedDescription()->utf8String():"?");return 1;}
    GMat mats[2]={toG(Material::basalt()),toG(Material::iron())};
    auto bmats=buf(sizeof(mats),mats);
    auto bpos=buf(3*n*4,pos.data()),bvel=buf(3*n*4,vel.data()),bu=buf(n*4,u.data()),bS=buf(6*n*4,Sv.data()),
         brho=buf(n*4,rho.data()),bhh=buf(n*4,hh.data()),bD=buf(n*4,Dd.data());
    auto bmass=buf(n*4,mass.data()),bmat=buf(n*4,mat.data()),beps=buf(n*4,epsact.data()),bfix=buf(n,fixed.data());
    auto baccold=buf(3*n*4,accold.data()),bduold=buf(n*4,dudtold.data()),bdSold=buf(6*n*4,dSdtold.data());
    auto bacc=buf(3*n*4),bdudt=buf(n*4),bdSdt=buf(6*n*4),bepsw=buf(n*4),bsig=buf(n*4);
    auto bPf=buf(n*4),bcsig=buf(n*4),bdsc=buf(n*4),bRart=buf(n*4);
    auto bhash=buf(n*4),bslot=buf(n*4),bsorted=buf(n*4);

    float lo[3],cw; int nc[3];
    auto rebuild=[&](){
        float hi[3]={-1e30f,-1e30f,-1e30f}; for(int a=0;a<3;a++) lo[a]=1e30f;
        float* P=(float*)bpos->contents(); float* H=(float*)bhh->contents(); float hmax=0;
        for(int i=0;i<n;i++){ for(int a=0;a<3;a++){lo[a]=std::min(lo[a],P[3*i+a]);hi[a]=std::max(hi[a],P[3*i+a]);} hmax=std::max(hmax,H[i]); }
        cw=2.0f*hmax; for(int a=0;a<3;a++) nc[a]=std::max(1,(int)((hi[a]-lo[a])/cw)+1);
        int Nc=nc[0]*nc[1]*nc[2];
        static MTL::Buffer* bcount=nullptr; static MTL::Buffer* bstart=nullptr; static int Ncprev=0;
        if(Nc!=Ncprev){ bcount=buf(Nc*4); bstart=buf(Nc*4); Ncprev=Nc; }
        memset(bcount->contents(),0,Nc*4);
        run(pso("cell_hash"),n,{bpos,bhash,bcount,bslot},{{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
        uint32_t* cnt=(uint32_t*)bcount->contents(); uint32_t* st=(uint32_t*)bstart->contents(); uint32_t a=0;
        for(int c=0;c<Nc;c++){st[c]=a;a+=cnt[c];}
        run(pso("cell_scatter"),n,{bhash,bslot,bstart,bsorted},{{&nn,4}});
        return std::make_pair(bcount,bstart);
    };

    // ---- KDK step ----
    float halfdt=0.5f*dt;
    run(pso("kick"), n, {bvel,baccold,bfix}, {{&halfdt,4},{&nn,4}});
    run(pso("drift"), n, {bpos,bvel,bfix}, {{&dt,4},{&nn,4}});
    MTL::Buffer *bcount,*bstart;
    for(int round=0;round<4;round++){ auto cs=rebuild(); bcount=cs.first; bstart=cs.second;
        float hmax0=0; { float* H=(float*)bhh->contents(); for(int i=0;i<n;i++) hmax0=std::max(hmax0,H[i]); }
        for(int it=0;it<2;it++) run(pso("density_adaptive"),n,{bpos,bmass,bsorted,bstart,bcount,brho,bhh},
            {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&eta,4},{&h_init,4},{&nn,4}});
        float hm=0; { float* H=(float*)bhh->contents(); for(int i=0;i<n;i++) hm=std::max(hm,H[i]); }
        if(hm<=hmax0*1.02f) break;
    }
    auto cs=rebuild(); bcount=cs.first; bstart=cs.second;
    run(pso("eos_full"),n,{brho,bu,bmat,bD,bPf,bcsig,bdsc,bRart,bmats},{{&eps_as,4},{&damon,4},{&nn,4}});
    run(pso("compute_sigmax"),n,{brho,bu,bmat,bS,bsig,bmats},{{&nn,4}});
    run(pso("strain_stress"),n,{bpos,bvel,brho,bhh,bmass,bS,bmat,bsorted,bstart,bcount,bmats,bdSdt,bepsw},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
    run(pso("forces"),n,{bpos,bvel,brho,bhh,bmass,bS,bPf,bcsig,bdsc,bRart,bepsw,bsorted,bstart,bcount,bacc,bdudt},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&eta,4},{&alpha,4},{&beta,4},{&eps,4},{&eps_as,4},{&nn,4}});
    run(pso("grow_damage"),n,{bsig,bD,beps,bmat,brho,bmats},{{&h_init,4},{&eta,4},{&dt,4},{&nn,4}});
    run(pso("finish"),n,{bvel,bacc,bu,bS,bdudt,bduold,bdSdt,bdSold,bfix,bmat,bmats},{{&dt,4},{&nn,4}});

    // ---- compare to CPU post-step ----
    float* Pg=(float*)bpos->contents(); float* Vg=(float*)bvel->contents(); float* Ug=(float*)bu->contents(); float* Sg=(float*)bS->contents();
    auto cmp=[&](const char* nm, float* g, int comp, int off){
        double sc=0,mx=0,ss=0; int c=0;
        for(int i=0;i<n;i++) for(int k=0;k<comp;k++) sc=std::max(sc,std::fabs(post[13*i+off+k]));
        for(int i=0;i<n;i++) for(int k=0;k<comp;k++){ double rel=std::fabs(g[comp*i+k]-post[13*i+off+k])/sc; mx=std::max(mx,rel); ss+=rel*rel; c++; }
        printf("  %-5s scale=%.3e  RMS=%.2e  max=%.2e\n", nm, sc, std::sqrt(ss/c), mx); return mx; };
    printf("GPU step vs CPU post-step:\n");
    double m1=cmp("pos",Pg,3,0), m2=cmp("vel",Vg,3,3), m3=cmp("u",Ug,1,6), m4=cmp("S",Sg,6,7);
    printf("%s\n",(m1<2e-3&&m2<2e-3&&m3<2e-3&&m4<2e-3)?"STEP GATE: PASS (FP32)":"STEP GATE: investigate");
    return 0;
}
