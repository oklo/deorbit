// Production GPU SPH driver: loads a real IC (run_ic format, incl frozen flags),
// adaptive CFL timestep, periodic full-res snapshots, walltime cap. Mirrors the
// CPU run_ic. Snapshots are "gpu_snap_NNN.bin" (10 col, like run_ic) so they
// don't collide with a CPU run's frac_snap_*.
//   ./gpu_ic <ic.bin> <dx> <t_end> [walltime] [--frac]
#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <Foundation/Foundation.hpp>
#include "../common/eos.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cstring>
#include <chrono>
#include <random>
#include <string>

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
    bool frac=false; float gz=0.0f; float snap_dt=0.05f; float damp=1.0f; std::vector<char*> av;
    for(int i=0;i<argc;i++){ std::string s=argv[i]; if(s=="--frac") frac=true; else if(s=="--g"&&i+1<argc){ gz=atof(argv[++i]); } else if(s=="--snap"&&i+1<argc){ snap_dt=atof(argv[++i]); } else if(s=="--damp"&&i+1<argc){ damp=atof(argv[++i]); } else av.push_back(argv[i]); }
    if(av.size()<4){ printf("usage: gpu_ic ic.bin dx t_end [walltime] [--frac] [--g g_mars]\n"); return 1; }
    const char* icf=av[1]; float dx=atof(av[2]); float t_end=atof(av[3]); double walltime=av.size()>4?atof(av[4]):1e30;

    // ---- load IC (run_ic format) ----
    FILE* f=fopen(icf,"rb"); if(!f){printf("cannot open %s\n",icf);return 1;}
    double hdr[7]; fread(hdr,8,7,f); long n=llround(hdr[6]);
    std::vector<float> pos(3*n),vel(3*n),mass(n),u(n); std::vector<int> mat(n); std::vector<unsigned char> fixed(n);
    int n_iron=0;
    for(long i=0;i<n;i++){ double r[10]; fread(r,8,10,f);
        pos[3*i]=r[0];pos[3*i+1]=r[1];pos[3*i+2]=r[2]; vel[3*i]=r[3];vel[3*i+1]=r[4];vel[3*i+2]=r[5];
        mass[i]=r[6]; u[i]=r[7]; mat[i]=(int)r[8]; fixed[i]=(r[9]!=0.0); if(mat[i]==1)n_iron++; }
    fclose(f);
    printf("loaded %s: n=%ld (%d iron)  dx=%.0f  frac=%d  g=%.2f\n", icf, n, n_iron, dx, frac, gz);

    float h_init=1.3f*dx, eta=1.3f, alpha=1.0f, beta=2.0f, eps=0.01f, cfl=0.2f;
    float eps_as = frac?0.3f:0.0f; int damon = frac?1:0;
    // eps_act: seed on host matching sph.hpp::seed_damage (basalt wk=1e61,wm=16; iron wk=0)
    std::vector<float> epsact(n,1e30f);
    if(frac){ std::mt19937 rng(12345); std::uniform_real_distribution<double> U(1e-9,1.0); double V=(double)dx*dx*dx;
        for(long i=0;i<n;i++){ if(mat[i]!=1) epsact[i]=(float)std::pow(U(rng)/(1e61*V),1.0/16.0); } }  // mat 0,2 = basalt-type

    D_=MTL::CreateSystemDefaultDevice(); Q_=D_->newCommandQueue();
    NS::Error* e=nullptr;
    std::string exe=argv[0]; size_t sl=exe.find_last_of('/');   // resolve the metallib relative to the executable, not the CWD (so gpu_ic runs from any dir)
    std::string libp=(sl==std::string::npos?std::string("."):exe.substr(0,sl))+"/sph_force.metallib";
    L_=D_->newLibrary(NS::String::string(libp.c_str(),NS::UTF8StringEncoding),&e);
    if(!L_){fprintf(stderr,"metallib load failed: %s\n  path: %s\n",e?e->localizedDescription()->utf8String():"(no NS::Error)",libp.c_str());return 1;}
    GMat mats[3]={toG(Material::basalt()),toG(Material::iron()),toG(Material::basalt())}; uint32_t nn=n;  // mat2 = basalt impactor (trackable)
    // fixed uniform grid: cell = 4*h_init (>= 2*max clamped hh)
    float cw=4.0f*h_init, lo[3]={(float)hdr[0],(float)hdr[2],(float)hdr[4]};
    int nc[3]={ (int)(((float)hdr[1]-lo[0])/cw)+1, (int)(((float)hdr[3]-lo[1])/cw)+1, (int)(((float)hdr[5]-lo[2])/cw)+1 };
    int Nc=nc[0]*nc[1]*nc[2];
    std::vector<float> hh(n,h_init);
    auto bpos=buf(3*n*4,pos.data()),bvel=buf(3*n*4,vel.data()),bu=buf(n*4,u.data()),bhh=buf(n*4,hh.data()),bmass=buf(n*4,mass.data());
    auto bS=buf(6*n*4); memset(bS->contents(),0,6*n*4);
    auto bD=buf(n*4); memset(bD->contents(),0,n*4);
    auto beps=buf(n*4,epsact.data()), bmat=buf(n*4,mat.data()), bfix=buf(n,fixed.data()), bmats=buf(sizeof(mats),mats);
    auto brho=buf(n*4),bacc=buf(3*n*4),bdudt=buf(n*4),bdSdt=buf(6*n*4),bepsw=buf(n*4),bsig=buf(n*4);
    auto bduold=buf(n*4),bdSold=buf(6*n*4),bPf=buf(n*4),bcsig=buf(n*4),bdsc=buf(n*4),bRart=buf(n*4);
    auto bhash=buf(n*4),bcount=buf(Nc*4),bslot=buf(n*4),bstart=buf(Nc*4),bsorted=buf(n*4);

    // cache pipeline states ONCE (were being recreated every step -> wasteful)
    auto PH=pso("cell_hash"),PSC=pso("cell_scatter"),PD=pso("density_adaptive"),PE=pso("eos_full"),
         PSG=pso("compute_sigmax"),PSF=pso("strain_forces"),PK=pso("kick"),
         PDR=pso("drift"),PG=pso("grow_damage"),PFN=pso("finish");
    auto PGW=pso("gather_w"),PGU8=pso("gather_u8"),PIOTA=pso("set_iota");
    double T_cells=0,T_dens=0,T_force=0,T_small=0,T_dt=0,T_blit=0,T_reord=0;   // phase timers
    // double buffers for particle reordering (coalesced neighbour reads)
    auto bpos_b=buf(3*n*4),bvel_b=buf(3*n*4),bS_b=buf(6*n*4),bdSdt_b=buf(6*n*4),
         bu_b=buf(n*4),bhh_b=buf(n*4),bmass_b=buf(n*4),bmat_b=buf(n*4),bD_b=buf(n*4),beps_b=buf(n*4),bdudt_b=buf(n*4),bfix_b=buf(n);
    auto gw=[&](MTL::Buffer*& src,MTL::Buffer*& dst,uint32_t W){ uint32_t w=W; run(PGW,n,{src,dst,bsorted},{{&w,4},{&nn,4}}); std::swap(src,dst); };
    auto reorder=[&](){     // gather carried state into cell-sorted order, then sorted=identity
        gw(bpos,bpos_b,3); gw(bvel,bvel_b,3); gw(bS,bS_b,6); gw(bdSdt,bdSdt_b,6);
        gw(bu,bu_b,1); gw(bhh,bhh_b,1); gw(bmass,bmass_b,1); gw(bmat,bmat_b,1); gw(bD,bD_b,1); gw(beps,beps_b,1); gw(bdudt,bdudt_b,1);
        run(PGU8,n,{bfix,bfix_b,bsorted},{{&nn,4}}); std::swap(bfix,bfix_b);
        run(PIOTA,n,{bsorted},{{&nn,4}}); };

    auto build=[&](){ memset(bcount->contents(),0,Nc*4);
        run(PH,n,{bpos,bhash,bcount,bslot},{{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&nn,4}});
        uint32_t* c=(uint32_t*)bcount->contents(); uint32_t* s=(uint32_t*)bstart->contents(); uint32_t a=0; for(int i=0;i<Nc;i++){s[i]=a;a+=c[i];}
        run(PSC,n,{bhash,bslot,bstart,bsorted},{{&nn,4}}); };
    auto density=[&](){ for(int it=0;it<1;it++) run(PD,n,{bpos,bmass,bsorted,bstart,bcount,brho,bhh},
        {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&eta,4},{&h_init,4},{&nn,4}}); };
    auto force_eval=[&](){
        run(PE,n,{brho,bu,bmat,bD,bPf,bcsig,bdsc,bRart,bmats},{{&eps_as,4},{&damon,4},{&nn,4}});
        run(PSG,n,{brho,bu,bmat,bS,bsig,bmats},{{&nn,4}});
        run(PSF,n,{bpos,bvel,brho,bhh,bmass,bS,bmat,bPf,bcsig,bdsc,bRart,bsorted,bstart,bcount,bmats,bdSdt,bacc,bdudt},
            {{&lo[0],4},{&lo[1],4},{&lo[2],4},{&cw,4},{&nc[0],4},{&nc[1],4},{&nc[2],4},{&eta,4},{&alpha,4},{&beta,4},{&eps,4},{&eps_as,4},{&nn,4}}); };
    auto adapt_dt=[&](){ float* H=(float*)bhh->contents(); float* C=(float*)bcsig->contents(); double dt=1e30;
        for(long i=0;i<n;i++){ double d=cfl*H[i]/(C[i]*(1.0+alpha)+1e-30); if(d<dt) dt=d; } return (float)dt; };
    auto energy=[&](){ float* V=(float*)bvel->contents(); float* Uu=(float*)bu->contents(); float* M=(float*)bmass->contents(); double KE=0,IE=0;
        for(long i=0;i<n;i++){ KE+=0.5*M[i]*(V[3*i]*V[3*i]+V[3*i+1]*V[3*i+1]+V[3*i+2]*V[3*i+2]); IE+=M[i]*Uu[i]; } return KE+IE; };
    auto report=[&](float t,int step,double wsec){ float* P=(float*)bpos->contents(); float* Vp=(float*)bvel->contents(); float* Uu=(float*)bu->contents(); int* MT=(int*)bmat->contents();
        double cx=0,cz=0,vz=0,sp=0; int ni=0,nm=0; double um=0;
        for(long i=0;i<n;i++){ if(MT[i]<1)continue; cx+=P[3*i];cz+=P[3*i+2];vz+=Vp[3*i+2];   // mat>=1 = projectile (iron@cayambe / mat2@orcus)
            sp+=std::sqrt(Vp[3*i]*Vp[3*i]+Vp[3*i+1]*Vp[3*i+1]+Vp[3*i+2]*Vp[3*i+2]); if(Uu[i]>1e6)nm++; um=std::max(um,(double)Uu[i]); ni++; }
        printf("  t=%.4f step=%d proj: x=%.0f z=%.0f |v|=%.0f vz=%.0f melt=%.0f%% umax=%.1e [%.0fs]\n",
            t,step,cx/ni,cz/ni,sp/ni,vz/ni,100.0*nm/ni,um,wsec); fflush(stdout); };
    auto snapshot=[&](int idx){ char fn[64]; snprintf(fn,64,"gpu_snap_%03d.bin",idx); FILE* o=fopen(fn,"wb");
        long N=n; fwrite(&N,8,1,o); float* P=(float*)bpos->contents(); float* Vp=(float*)bvel->contents();
        float* Uu=(float*)bu->contents(); float* Rho=(float*)brho->contents(); float* H=(float*)bhh->contents(); int* MT=(int*)bmat->contents();
        std::vector<double> rec(10*n); for(long i=0;i<n;i++){ double* r=&rec[10*i];
            r[0]=P[3*i];r[1]=P[3*i+1];r[2]=P[3*i+2]; r[3]=Vp[3*i];r[4]=Vp[3*i+1];r[5]=Vp[3*i+2];
            r[6]=Uu[i];r[7]=Rho[i];r[8]=H[i];r[9]=MT[i]; } fwrite(rec.data(),8,10*n,o); fclose(o); };

    build(); density(); build(); force_eval();
    double E0=energy(); printf("init done, E0=%.4e\n",E0);
    double t=0,wall0=now(); int step=0,isnap=0; float next_snap=snap_dt;
    while(t<t_end){
        double a;
        a=now(); float dt=adapt_dt(); T_dt+=now()-a; if(t+dt>t_end) dt=t_end-t; float halfdt=0.5f*dt;
        a=now(); run(PK,n,{bvel,bacc,bfix},{{&halfdt,4},{&gz,4},{&nn,4}}); run(PDR,n,{bpos,bvel,bfix},{{&dt,4},{&nn,4}}); T_small+=now()-a;
        a=now(); build(); T_cells+=now()-a;
        a=now(); reorder(); T_reord+=now()-a;
        a=now(); auto bl=Q_->commandBuffer(); auto be=bl->blitCommandEncoder();
        be->copyFromBuffer(bdudt,0,bduold,0,n*4); be->copyFromBuffer(bdSdt,0,bdSold,0,6*n*4); be->endEncoding(); bl->commit(); bl->waitUntilCompleted(); T_blit+=now()-a;
        a=now(); density(); T_dens+=now()-a;
        a=now(); force_eval(); T_force+=now()-a;
        a=now(); run(PG,n,{bsig,bD,beps,bmat,brho,bmats},{{&h_init,4},{&eta,4},{&dt,4},{&nn,4}});
        run(PFN,n,{bvel,bacc,bu,bS,bdudt,bduold,bdSdt,bdSold,bfix,bmat,bmats},{{&dt,4},{&gz,4},{&damp,4},{&nn,4}}); T_small+=now()-a;
        t+=dt; step++;
        if(step%25==0||t>=t_end) report(t,step,now()-wall0);
        if(t>=next_snap){ snapshot(isnap++); next_snap+=snap_dt; }
        if(now()-wall0>walltime){ printf("walltime hit\n"); break; }
    }
    double E1=energy();
    printf("DONE: %d steps to t=%.4f  %.1f ms/step  E drift=%.2e\n", step, t, 1000*(now()-wall0)/step, (E1-E0)/E0);
    printf("PROFILE (ms/step): cells=%.1f reorder=%.1f density=%.1f force_eval=%.1f small=%.1f dt_read=%.1f blit=%.1f\n",
        1000*T_cells/step,1000*T_reord/step,1000*T_dens/step,1000*T_force/step,1000*T_small/step,1000*T_dt/step,1000*T_blit/step);
    return 0;
}
