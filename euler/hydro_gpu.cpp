// euler M1b: Metal GPU hydro host. Runs Sod + Sedov, validates vs analytic
// (same gates as the CPU oracle hydro_cpu.cpp). FP32.
//   ./hydro_gpu sod | sedov
#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <Foundation/Foundation.hpp>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
static float GAM=1.4f;
static MTL::Device* D_; static MTL::CommandQueue* Q_; static MTL::Library* L_;
static MTL::Buffer* buf(size_t b){ return D_->newBuffer(b, MTL::ResourceStorageModeShared); }
static MTL::ComputePipelineState* pso(const char* nm){ NS::Error* e=nullptr;
    auto p=D_->newComputePipelineState(L_->newFunction(NS::String::string(nm,NS::UTF8StringEncoding)),&e);
    if(!p){printf("pso %s failed\n",nm);exit(1);} return p; }
static void run(MTL::ComputePipelineState* p,uint32_t n,vector<MTL::Buffer*> bs,vector<pair<const void*,size_t>> by){
    auto cb=Q_->commandBuffer(); auto e=cb->computeCommandEncoder(); e->setComputePipelineState(p);
    for(size_t i=0;i<bs.size();i++) e->setBuffer(bs[i],0,i);
    for(size_t i=0;i<by.size();i++) e->setBytes(by[i].first,by[i].second,bs.size()+i);
    e->dispatchThreads(MTL::Size(n,1,1),MTL::Size(256,1,1)); e->endEncoding(); cb->commit(); cb->waitUntilCompleted();
}
// exact Sod (Toro)
static void exact_sod(double S,double&r,double&u,double&p){
    double rL=1,uL=0,pL=1,rR=0.125,uR=0,pR=0.1,g=GAM; double cL=sqrt(g*pL/rL),cR=sqrt(g*pR/rR);
    double G1=(g-1)/(2*g),G3=2*g/(g-1),G4=2/(g-1),G5=2/(g+1),G6=(g-1)/(g+1),G2=(g+1)/(2*g);
    auto fk=[&](double P,double rk,double pk,double ck){if(P>pk){double A=G5/rk,B=G6*pk;return(P-pk)*sqrt(A/(P+B));}return G4*ck*(pow(P/pk,G1)-1);};
    auto fp=[&](double P,double rk,double pk,double ck){if(P>pk){double A=G5/rk,B=G6*pk;return sqrt(A/(B+P))*(1-(P-pk)/(2*(B+P)));}return pow(P/pk,-G2)/(rk*ck);};
    double P=0.5*(pL+pR);for(int it=0;it<100;it++){double f=fk(P,rL,pL,cL)+fk(P,rR,pR,cR)+(uR-uL);double d=fp(P,rL,pL,cL)+fp(P,rR,pR,cR);double Pn=fabs(P-f/d);if(fabs(Pn-P)<1e-12){P=Pn;break;}P=Pn;}
    double us=0.5*(uL+uR)+0.5*(fk(P,rR,pR,cR)-fk(P,rL,pL,cL));
    if(S<=us){if(P>pL){double SL=uL-cL*sqrt(G2*P/pL+G1);if(S<SL){r=rL;u=uL;p=pL;}else{r=rL*((P/pL+G6)/(G6*P/pL+1));u=us;p=P;}}
        else{double SHL=uL-cL,STL=us-cL*pow(P/pL,G1);if(S<SHL){r=rL;u=uL;p=pL;}else if(S>STL){r=rL*pow(P/pL,1/g);u=us;p=P;}else{u=G5*(cL+(g-1)/2*uL+S);double c=G5*(cL+(g-1)/2*(uL-S));r=rL*pow(c/cL,G4);p=pL*pow(c/cL,G3);}}}
    else{if(P>pR){double SR=uR+cR*sqrt(G2*P/pR+G1);if(S>SR){r=rR;u=uR;p=pR;}else{r=rR*((P/pR+G6)/(G6*P/pR+1));u=us;p=P;}}
        else{double SHR=uR+cR,STR=us+cR*pow(P/pR,G1);if(S>SHR){r=rR;u=uR;p=pR;}else if(S<STR){r=rR*pow(P/pR,1/g);u=us;p=P;}else{u=G5*(-cR+(g-1)/2*uR+S);double c=G5*(cR-(g-1)/2*(uR-S));r=rR*pow(c/cR,G4);p=pR*pow(c/cR,G3);}}}
}
int main(int argc,char**argv){
    string mode=argc>1?argv[1]:"sod";
    int nx,ny,nz; double Ldom,tend,CFL; double rho0=1,E0=1,p0=1e-4;
    if(mode=="sod"){nx=200;ny=nz=1;Ldom=1.0;tend=0.2;CFL=0.4;}
    else {nx=ny=nz=64;Ldom=2.0;tend=0.5;CFL=0.3;}
    double dx=Ldom/nx; uint32_t n=nx*ny*nz; float invdx=1.0f/dx;
    D_=MTL::CreateSystemDefaultDevice(); Q_=D_->newCommandQueue(); NS::Error* e=nullptr;
    L_=D_->newLibrary(NS::String::string("hydro.metallib",NS::UTF8StringEncoding),&e);
    if(!L_){printf("metallib load failed\n");return 1;}
    auto br=buf(n*4),bmu=buf(n*4),bmv=buf(n*4),bmw=buf(n*4),bE=buf(n*4);
    auto br1=buf(n*4),bmu1=buf(n*4),bmv1=buf(n*4),bmw1=buf(n*4),bE1=buf(n*4);
    auto bdr=buf(n*4),bdmu=buf(n*4),bdmv=buf(n*4),bdmw=buf(n*4),bdE=buf(n*4),bsp=buf(n*4);
    float*r=(float*)br->contents(),*mu=(float*)bmu->contents(),*mv=(float*)bmv->contents(),*mw=(float*)bmw->contents(),*E=(float*)bE->contents();
    // IC
    if(mode=="sod"){ for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5?1.0f:0.125f,pp=x<0.5?1.0f:0.1f;
        r[i]=rr;mu[i]=mv[i]=mw[i]=0;E[i]=pp/(GAM-1);} }
    else { for(uint32_t c=0;c<n;c++){r[c]=rho0;mu[c]=mv[c]=mw[c]=0;E[c]=p0/(GAM-1);}
        int ic=nx/2; int cc=(ic*ny+ic)*nz+ic; E[cc]+=E0/(dx*dx*dx); }
    auto Plop=pso("lop"),Pws=pso("wavespeed"),Prk1=pso("rk1"),Prk2=pso("rk2");
    int NX=nx,NY=ny,NZ=nz;
    double t=0; int step=0;
    while(t<tend){
        run(Pws,n,{br,bmu,bmv,bmw,bE,bsp},{{&GAM,4},{&n,4}});
        float*sp=(float*)bsp->contents(); double smax=1e-30; for(uint32_t c=0;c<n;c++) smax=max(smax,(double)sp[c]);
        float dt=CFL*dx/smax; if(t+dt>tend)dt=tend-t;
        run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE},{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&GAM,4},{&n,4}});
        run(Prk1,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,br1,bmu1,bmv1,bmw1,bE1},{{&dt,4},{&n,4}});
        run(Plop,n,{br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&GAM,4},{&n,4}});
        run(Prk2,n,{br,bmu,bmv,bmw,bE,br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},{{&dt,4},{&n,4}});
        t+=dt; step++;
    }
    auto pres=[&](uint32_t c){ float ke=0.5f*(mu[c]*mu[c]+mv[c]*mv[c]+mw[c]*mw[c])/r[c]; return (GAM-1)*(E[c]-ke); };
    if(mode=="sod"){
        double l1r=0,l1p=0,l1u=0; for(int i=0;i<nx;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);
            l1r+=fabs(r[i]-re);l1p+=fabs(pres(i)-pe);l1u+=fabs(mu[i]/r[i]-ue);} l1r/=nx;l1p/=nx;l1u/=nx;
        // dump rho profile for GPU-vs-CPU diff
        FILE*o=fopen("sod_gpu.txt","w"); for(int i=0;i<nx;i++)fprintf(o,"%.8f %.8f\n",(i+0.5)*dx,r[i]); fclose(o);
        printf("GPU Sod steps=%d  L1 vs EXACT: rho=%.4f p=%.4f u=%.4f\n",step,l1r,l1p,l1u);
        printf("GATE (match CPU ~0.004): %s\n",(l1r<0.006&&l1p<0.006&&l1u<0.008)?"PASS":"CHECK");
    } else {
        double rmax=0,Rsh=0,xc=(nx/2+0.5)*dx;
        for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){uint32_t c=(i*ny+j)*nz+k;
            if(r[c]>rmax){rmax=r[c];double X=(i+0.5)*dx-xc,Y=(j+0.5)*dx-xc,Z=(k+0.5)*dx-xc;Rsh=sqrt(X*X+Y*Y+Z*Z);}}
        double Ranal=pow(E0*t*t/(0.851*rho0),0.2);
        printf("GPU Sedov steps=%d t=%.3f  R_num=%.4f R_analytic=%.4f (err %.1f%%)  peakcomp=%.2f\n",step,t,Rsh,Ranal,100*fabs(Rsh-Ranal)/Ranal,rmax/rho0);
        printf("GATE (R<10%%): %s\n",(fabs(Rsh-Ranal)/Ranal<0.10)?"PASS":"CHECK");
    }
    return 0;
}
