// Milestone 5: GPU multi-step run -> benchmark (ms/step, particles/s) + energy
// conservation. Fixed uniform cell grid (cell = 4*h_init >= 2*max clamped h) so
// there's no per-step O(N) CPU readback; only the small cell-count scan is on CPU.
//   ./gpu_run <L> <steps>     (N = L^3 particles)
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
#include <chrono>
#include <random>

struct GMat { float rho0,A,B,a,b,alpha,beta,u0,uiv,ucv,G,Emod,Y; };
static GMat toG(const Material& m){ return {(float)m.rho0,(float)m.A,(float)m.B,(float)m.a,(float)m.b,(float)m.alpha,
    (float)m.beta,(float)m.u0,(float)m.uiv,(float)m.ucv,(float)m.G,(float)m.Emod,(float)m.Y}; }
static MTL::Device* D_; static MTL::CommandQueue* Q_; static MTL::Library* L_;
static MTL::Buffer* buf(size_t b,const void* s=nullptr){ return s?D_->newBuffer(s,b,MTL::ResourceStorageModeShared):D_->newBuffer(b,MTL::ResourceStorageModeShared); }
static MTL::ComputePipelineState* pso(const char* nm){ NS::Error* e=nullptr;
    auto p=D_->newComputePipelineState(L_->newFunction(NS::String::string(nm,NS::UTF8StringEncoding)),&e);
    if(!p){printf("pso %s\n",nm);exit(1);} return p; }
static void run(MTL::ComputePipelineState* p,uint32_t n,std::vector<MTL::Buffer*> bs,std::vector<std::pair<const void*,size_t>> by){
    auto cb=Q_->commandBuffer(); auto e=cb->computeCommandEncoder(); e->setComputePipelineState(p);
    for(size_t i=0;i<bs.size();i++) e->setBuffer(bs[i],0,i);
    for(size_t i=0;i<by.size();i++) e->setBytes(by[i].first,by[i].second,bs.size()+i);
    e->dispatchThreads(MTL::Size(n,1,1),MTL::Size(256,1,1)); e->endEncoding(); cb->commit(); cb->waitUntilCompleted();
}
static double now(){ return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc,char** argv){
    int Lc = argc>1?atoi(argv[1]):100; int steps = argc>2?atoi(argv[2]):20;
    const double dx=30.0, h_initd=1.3*dx, R=0.2*Lc*dx;
    int n=Lc*Lc*Lc; double lo0=-Lc*dx/2.0;
    std::vector<float> pos(3*n),vel(3*n),u(n),mass(n); std::vector<int> mat(n); std::vector<unsigned char> fixed(n,0);
    std::mt19937 rng(1); std::uniform_real_distribution<double> J(-0.2*dx,0.2*dx);
    for(int a=0,k=0;a<Lc;a++)for(int b=0;b<Lc;b++)for(int c=0;c<Lc;c++,k++){
        double x=lo0+(a+0.5)*dx+J(rng),y=lo0+(b+0.5)*dx+J(rng),z=lo0+(c+0.5)*dx+J(rng);
        bool iron=(x*x+y*y+z*z)<R*R; mat[k]=iron?1:0;
        pos[3*k]=x;pos[3*k+1]=y;pos[3*k+2]=z; vel[3*k]=80.0*(2*b/(double)Lc-1); vel[3*k+1]=0; vel[3*k+2]=-40.0*(2*c/(double)Lc-1);
        u[k]=2e5; mass[k]=(iron?7800.0:2700.0)*dx*dx*dx; }

    D_=MTL::CreateSystemDefaultDevice(); Q_=D_->newCommandQueue();
    NS::Error* e=nullptr; L_=D_->newLibrary(NS::String::string("sph_force.metallib",NS::UTF8StringEncoding),&e);
    if(!L_){printf("lib\n");return 1;}
    GMat mats[2]={toG(Material::basalt()),toG(Material::iron())};
    float h_init=(float)h_initd, eta=1.3f, alpha=1.0f, beta=2.0f, eps=0.01f, eps_as=0.3f; int damon=1; uint32_t nn=n;
    // fixed grid: cell = 4*h_init (>= 2*max clamped hh = 2*2*h_init)
    float cw=4.0f*h_init, lo[3]={(float)(lo0-dx),(float)(lo0-dx),(float)(lo0-dx)};
    float ext=(float)(Lc*dx+2*dx); int nc[3]={(int)(ext/cw)+1,(int)(ext/cw)+1,(int)(ext/cw)+1};
    int Nc=nc[0]*nc[1]*nc[2];
    std::vector<float> hh(n,h_init);
    auto bpos=buf(3*n*4,pos.data()),bvel=buf(3*n*4,vel.data()),bu=buf(n*4,u.data()),
         bhh=buf(n*4,hh.data()),bmass=buf(n*4,mass.data());
    auto bS=buf(6*n*4); memset(bS->contents(),0,6*n*4);
    auto bD=buf(n*4); memset(bD->contents(),0,n*4);
    auto beps=buf(n*4); { float* p=(float*)beps->contents(); for(int i=0;i<n;i++)p[i]=1e30f; }
    auto bmat=buf(n*4,mat.data()),bfix=buf(n,fixed.data()),bmats=buf(sizeof(mats),mats);
    auto brho=buf(n*4),bacc=buf(3*n*4),bdudt=buf(n*4),bdSdt=buf(6*n*4),bepsw=buf(n*4),bsig=buf(n*4);
    auto bduold=buf(n*4),bdSold=buf(6*n*4),bPf=buf(n*4),bcsig=buf(n*4),bdsc=buf(n*4),bRart=buf(n*4);
    auto bhash=buf(n*4),bcount=buf(Nc*4),bslot=buf(n*4),bstart=buf(Nc*4),bsorted=buf(n*4);

    auto build=[&](){ memset(bcount->contents(),0,Nc*4);
        run(pso("cell_hash"),n,{bpos,bhash,bcount,bslot},{{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
        uint32_t* cnt=(uint32_t*)bcount->contents(); uint32_t* st=(uint32_t*)bstart->contents(); uint32_t a=0;
        for(int c=0;c<Nc;c++){st[c]=a;a+=cnt[c];}
        run(pso("cell_scatter"),n,{bhash,bslot,bstart,bsorted},{{&nn,4}}); };
    auto density=[&](){ for(int it=0;it<4;it++) run(pso("density_adaptive"),n,{bpos,bmass,bsorted,bstart,bcount,brho,bhh},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&eta,4},{&h_init,4},{&nn,4}}); };
    auto force_eval=[&](){
        run(pso("eos_full"),n,{brho,bu,bmat,bD,bPf,bcsig,bdsc,bRart,bmats},{{&eps_as,4},{&damon,4},{&nn,4}});
        run(pso("compute_sigmax"),n,{brho,bu,bmat,bS,bsig,bmats},{{&nn,4}});
        run(pso("strain_stress"),n,{bpos,bvel,brho,bhh,bmass,bS,bmat,bsorted,bstart,bcount,bmats,bdSdt,bepsw},
            {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
        run(pso("forces"),n,{bpos,bvel,brho,bhh,bmass,bS,bPf,bcsig,bdsc,bRart,bepsw,bsorted,bstart,bcount,bacc,bdudt},
            {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&eta,4},{&alpha,4},{&beta,4},{&eps,4},{&eps_as,4},{&nn,4}}); };
    auto energy=[&](){ float* V=(float*)bvel->contents(); float* U=(float*)bu->contents(); float* M=(float*)bmass->contents();
        double KE=0,IE=0; for(int i=0;i<n;i++){ KE+=0.5*M[i]*(V[3*i]*V[3*i]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]); IE+=M[i]*U[i]; } return KE+IE; };

    build(); density(); build(); force_eval();            // init
    double E0=energy(); float dt=2.0e-4f, halfdt=0.5f*dt;
    double t0=now();
    for(int s=0;s<steps;s++){
        auto bl=Q_->commandBuffer(); auto be=bl->blitCommandEncoder();
        be->copyFromBuffer(bdudt,0,bduold,0,n*4); be->copyFromBuffer(bdSdt,0,bdSold,0,6*n*4);
        be->endEncoding(); bl->commit(); bl->waitUntilCompleted();
        run(pso("kick"),n,{bvel,bacc,bfix},{{&halfdt,4},{&nn,4}});
        run(pso("drift"),n,{bpos,bvel,bfix},{{&dt,4},{&nn,4}});
        build(); density(); force_eval();
        run(pso("grow_damage"),n,{bsig,bD,beps,bmat,brho,bmats},{{&h_init,4},{&eta,4},{&dt,4},{&nn,4}});
        run(pso("finish"),n,{bvel,bacc,bu,bS,bdudt,bduold,bdSdt,bdSold,bfix,bmat,bmats},{{&dt,4},{&nn,4}});
    }
    double el=now()-t0, E1=energy();
    float* V=(float*)bvel->contents(); double vmax=0; for(int i=0;i<n;i++) vmax=std::max(vmax,(double)std::hypot(std::hypot(V[3*i],V[3*i+1]),V[3*i+2]));
    printf("GPU n=%d (%dx^3)  %d steps  %.1f ms/step  %.0f Mparticle-steps/s   E drift=%.2e  vmax=%.0f\n",
           n, Lc, steps, 1000*el/steps, n/(el/steps)/1e6, (E1-E0)/E0, vmax);
    return 0;
}
