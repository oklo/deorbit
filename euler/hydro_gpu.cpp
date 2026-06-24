// euler M2b: Metal GPU hydro host with pluggable EOS (ideal | Tillotson) + free surface.
//   ./hydro_gpu sod|sedov   (ideal-gas gates)   |   surface|bshock (basalt)
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
struct GMat { float rho0,A,B,a,b,alpha,beta,u0,uiv,ucv; };
static GMat BASALT={2700,2.67e10f,2.67e10f,0.5f,1.5f,5,5,4.87e8f,4.72e6f,1.82e7f};
static float GAM=1.4f;
static MTL::Device* D_; static MTL::CommandQueue* Q_; static MTL::Library* L_;
static MTL::Buffer* buf(size_t b){ return D_->newBuffer(b, MTL::ResourceStorageModeShared); }
static MTL::ComputePipelineState* pso(const char* nm){ NS::Error* e=nullptr;
    auto p=D_->newComputePipelineState(L_->newFunction(NS::String::string(nm,NS::UTF8StringEncoding)),&e);
    if(!p){printf("pso %s\n",nm);exit(1);} return p; }
static void run(MTL::ComputePipelineState* p,uint32_t n,vector<MTL::Buffer*> bs,vector<pair<const void*,size_t>> by){
    auto cb=Q_->commandBuffer(); auto e=cb->computeCommandEncoder(); e->setComputePipelineState(p);
    for(size_t i=0;i<bs.size();i++) e->setBuffer(bs[i],0,i);
    for(size_t i=0;i<by.size();i++) e->setBytes(by[i].first,by[i].second,bs.size()+i);
    e->dispatchThreads(MTL::Size(n,1,1),MTL::Size(256,1,1)); e->endEncoding(); cb->commit(); cb->waitUntilCompleted();
}
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
    int nx,ny,nz,emode; double Ldom,tend,CFL;
    if(mode=="sod"){nx=200;ny=nz=1;Ldom=1.0;tend=0.2;CFL=0.4;emode=0;}
    else if(mode=="sedov"){nx=ny=nz=64;Ldom=2.0;tend=0.5;CFL=0.3;emode=0;}
    else if(mode=="surface"){nx=1;ny=1;nz=64;Ldom=200e3;tend=2.0;CFL=0.4;emode=1;}
    else {nx=400;ny=nz=1;Ldom=400e3;tend=4.0;CFL=0.4;emode=1;}   // bshock
    double dx=Ldom/(mode=="surface"?nz:nx); uint32_t n=nx*ny*nz; float invdx=1.0f/dx; float gam=GAM;
    D_=MTL::CreateSystemDefaultDevice(); Q_=D_->newCommandQueue(); NS::Error* e=nullptr;
    L_=D_->newLibrary(NS::String::string("hydro.metallib",NS::UTF8StringEncoding),&e);
    if(!L_){printf("metallib load failed\n");return 1;}
    auto br=buf(n*4),bmu=buf(n*4),bmv=buf(n*4),bmw=buf(n*4),bE=buf(n*4);
    auto br1=buf(n*4),bmu1=buf(n*4),bmv1=buf(n*4),bmw1=buf(n*4),bE1=buf(n*4);
    auto bdr=buf(n*4),bdmu=buf(n*4),bdmv=buf(n*4),bdmw=buf(n*4),bdE=buf(n*4),bsp=buf(n*4);
    float*r=(float*)br->contents(),*mu=(float*)bmu->contents(),*mv=(float*)bmv->contents(),*mw=(float*)bmw->contents(),*E=(float*)bE->contents();
    if(mode=="sod"){for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5?1.0f:0.125f,pp=x<0.5?1.0f:0.1f;r[i]=rr;mu[i]=mv[i]=mw[i]=0;E[i]=pp/(GAM-1);}}
    else if(mode=="sedov"){for(uint32_t c=0;c<n;c++){r[c]=1;mu[c]=mv[c]=mw[c]=0;E[c]=1e-4/(GAM-1);}int ic=nx/2;E[(ic*ny+ic)*nz+ic]+=1.0/(dx*dx*dx);}
    else if(mode=="surface"){for(int k=0;k<nz;k++){double z=(k+0.5)*dx;float rr=z<0.5*Ldom?2700.0f:0.27f;r[k]=rr;mu[k]=mv[k]=mw[k]=0;E[k]=0;}}
    else {for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5*Ldom?3000.0f:2700.0f,ee=x<0.5*Ldom?1e6f:0.0f;r[i]=rr;mu[i]=mv[i]=mw[i]=0;E[i]=rr*ee;}}
    auto Plop=pso("lop"),Pws=pso("wavespeed"),Prk1=pso("rk1"),Prk2=pso("rk2");
    int NX=nx,NY=ny,NZ=nz;
    double t=0;int step=0;
    while(t<tend){
        run(Pws,n,{br,bmu,bmv,bmw,bE,bsp},{{&emode,4},{&gam,4},{&BASALT,sizeof(GMat)},{&n,4}});
        float*sp=(float*)bsp->contents();double smax=1e-30;for(uint32_t c=0;c<n;c++)smax=max(smax,(double)sp[c]);
        float dt=CFL*dx/smax;if(t+dt>tend)dt=tend-t;
        auto lopargs=vector<pair<const void*,size_t>>{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&emode,4},{&gam,4},{&BASALT,sizeof(GMat)},{&n,4}};
        run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE},lopargs);
        run(Prk1,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,br1,bmu1,bmv1,bmw1,bE1},{{&dt,4},{&n,4}});
        run(Plop,n,{br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},lopargs);
        run(Prk2,n,{br,bmu,bmv,bmw,bE,br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},{{&dt,4},{&n,4}});
        t+=dt;step++;
    }
    if(mode=="sod"){
        auto pres=[&](int c){float ke=0.5f*(mu[c]*mu[c]+mv[c]*mv[c]+mw[c]*mw[c])/r[c];return (GAM-1)*(E[c]-ke);};
        double l1r=0,l1p=0,l1u=0;for(int i=0;i<nx;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);l1r+=fabs(r[i]-re);l1p+=fabs(pres(i)-pe);l1u+=fabs(mu[i]/r[i]-ue);}
        l1r/=nx;l1p/=nx;l1u/=nx; printf("GPU Sod steps=%d L1 vs EXACT: rho=%.4f p=%.4f u=%.4f  GATE(==CPU): %s\n",step,l1r,l1p,l1u,(l1r<0.007)?"PASS":"CHECK");
    } else if(mode=="surface"){
        double vmax=0;for(int k=0;k<nz;k++)vmax=max(vmax,(double)fabs(mw[k]/r[k]));
        printf("GPU free surface steps=%d t=%.2f max|v|=%.3f m/s  GATE(static): %s\n",step,t,vmax,vmax<5.0?"PASS":"CHECK");
    } else if(mode=="bshock"){
        FILE*o=fopen("bshock_gpu.txt","w");for(int i=0;i<nx;i++)fprintf(o,"%.8e\n",r[i]);fclose(o);
        printf("GPU bshock steps=%d t=%.2f  wrote bshock_gpu.txt (diff vs CPU)\n",step,t);
    } else printf("GPU sedov done steps=%d\n",step);
    return 0;
}
