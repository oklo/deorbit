// euler M3b: Metal GPU hydro + EOS + free surface + STRENGTH (elastic-plastic).
//   ./hydro_gpu sod|sedov | surface|bshock | shear|yield
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
struct GMat { float rho0,A,B,a,b,alpha,beta,u0,uiv,ucv,G,Y,Emod; };
static GMat BASALT={2700,2.67e10f,2.67e10f,0.5f,1.5f,5,5,4.87e8f,4.72e6f,1.82e7f,2.27e10f,3.5e8f,5.306e10f};
static float GAM=1.4f;
static MTL::Device* D_; static MTL::CommandQueue* Q_; static MTL::Library* L_;
static MTL::Buffer* buf(size_t b){ auto p=D_->newBuffer(b, MTL::ResourceStorageModeShared); memset(p->contents(),0,b); return p; }
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
    double rL=1,uL=0,pL=1,rR=0.125,uR=0,pR=0.1,g=GAM;double cL=sqrt(g*pL/rL),cR=sqrt(g*pR/rR);
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
    string mode=argc>1?argv[1]:"sod"; double rho0=2700,G=BASALT.G,Y=BASALT.Y,csb=sqrt(G/rho0);
    int nx,ny,nz,emode; double Ldom,tend,CFL; float gz=0.0f;
    if(mode=="sod"){nx=200;ny=nz=1;Ldom=1.0;tend=0.2;CFL=0.4;emode=0;}
    else if(mode=="sedov"){nx=ny=nz=64;Ldom=2.0;tend=0.5;CFL=0.3;emode=0;}
    else if(mode=="surface"){nx=1;ny=1;nz=64;Ldom=200e3;tend=2.0;CFL=0.4;emode=1;}
    else if(mode=="bshock"){nx=400;ny=nz=1;Ldom=400e3;tend=4.0;CFL=0.4;emode=1;}
    else if(mode=="shear"){nx=800;ny=nz=1;Ldom=400e3;tend=20e3/csb;CFL=0.3;emode=1;}
    else if(mode=="freefall"){nx=ny=1;nz=50;Ldom=500.0;tend=10.0;CFL=0.4;emode=0;gz=9.8f;}
    else if(mode=="atmos"){nx=ny=1;nz=200;double H=1e5/9.8;Ldom=4*H;tend=0.05*Ldom/sqrt(GAM*1e5);CFL=0.4;emode=0;gz=9.8f;}
    else if(mode=="tensile"){nx=100;ny=nz=1;Ldom=10e3;tend=0.1;CFL=0.3;emode=1;}
    else {nx=400;ny=nz=1;Ldom=400e3;tend=15e3/csb;CFL=0.3;emode=1;}   // yield
    double dx=Ldom/((nx==1&&ny==1)?nz:nx); uint32_t n=nx*ny*nz; float invdx=1.0f/dx; float gam=GAM; bool strn=(emode==1);
    float epsact=(float)pow(1.0/(1e61*dx*dx*dx),1.0/16.0);   // Weibull weakest-flaw activation strain (host: wk=1e61 overflows FP32)
    D_=MTL::CreateSystemDefaultDevice(); Q_=D_->newCommandQueue(); NS::Error* e=nullptr;
    L_=D_->newLibrary(NS::String::string("hydro.metallib",NS::UTF8StringEncoding),&e);
    if(!L_){printf("metallib load failed\n");return 1;}
    auto br=buf(n*4),bmu=buf(n*4),bmv=buf(n*4),bmw=buf(n*4),bE=buf(n*4);
    auto br1=buf(n*4),bmu1=buf(n*4),bmv1=buf(n*4),bmw1=buf(n*4),bE1=buf(n*4);
    auto bdr=buf(n*4),bdmu=buf(n*4),bdmv=buf(n*4),bdmw=buf(n*4),bdE=buf(n*4),bsp=buf(n*4);
    auto bxx=buf(n*4),byy=buf(n*4),bzz=buf(n*4),bxy=buf(n*4),bxz=buf(n*4),byz=buf(n*4);
    auto dxx=buf(n*4),dyy=buf(n*4),dzz=buf(n*4),dxy=buf(n*4),dxz=buf(n*4),dyz=buf(n*4);
    auto sxx1=buf(n*4),syy1=buf(n*4),szz1=buf(n*4),sxy1=buf(n*4),sxz1=buf(n*4),syz1=buf(n*4);
    auto bD=buf(n*4),bdD=buf(n*4),bD1=buf(n*4);   // damage, rate, RK temp
    float*r=(float*)br->contents(),*mu=(float*)bmu->contents(),*mv=(float*)bmv->contents(),*mw=(float*)bmw->contents(),*E=(float*)bE->contents();
    float*pxy=(float*)bxy->contents(); // for shear/yield diagnostics
    if(mode=="sod"){for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5?1.0f:0.125f,pp=x<0.5?1.0f:0.1f;r[i]=rr;E[i]=pp/(GAM-1);}}
    else if(mode=="sedov"){for(uint32_t c=0;c<n;c++){r[c]=1;E[c]=1e-4/(GAM-1);}int ic=nx/2;E[(ic*ny+ic)*nz+ic]+=1.0/(dx*dx*dx);}
    else if(mode=="surface"){for(int k=0;k<nz;k++){double z=(k+0.5)*dx;r[k]=z<0.5*Ldom?2700.0f:0.27f;E[k]=0;}}
    else if(mode=="bshock"){for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5*Ldom?3000.0f:2700.0f,ee=x<0.5*Ldom?1e6f:0.0f;r[i]=rr;E[i]=rr*ee;}}
    else if(mode=="freefall"){for(int k=0;k<nz;k++){r[k]=1.0f;E[k]=1e5/(GAM-1);}}
    else if(mode=="atmos"){double H=1e5/9.8;for(int k=0;k<nz;k++){double z=(k+0.5)*dx;r[k]=exp(-z/H);E[k]=(1e5*exp(-z/H))/(GAM-1);}}
    else if(mode=="tensile"){double rate=1e-2;for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=rho0;double vx=rate*(x-0.5*Ldom);mu[i]=rho0*vx;E[i]=0.5*rho0*vx*vx;}}
    else {double x0=0.3*Ldom,wid=8e3,A=(mode=="shear"?1.0:2000.0);for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=rho0;double vy=A*exp(-((x-x0)/wid)*((x-x0)/wid));mv[i]=rho0*vy;E[i]=0.5*rho0*vy*vy;}}
    auto Plop=pso("lop"),Pws=pso("wavespeed"),Prk1=pso("rk1"),Prk2=pso("rk2"),Pstr=pso("strength"),Pvm=pso("vonmises"),Prk1s=pso("rk1s"),Prk2s=pso("rk2s"),Pgd=pso("grow_damage");
    int NX=nx,NY=ny,NZ=nz;
    auto lopA=vector<pair<const void*,size_t>>{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&emode,4},{&gam,4},{&BASALT,sizeof(GMat)},{&n,4},{&gz,4}};
    auto strA=vector<pair<const void*,size_t>>{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&BASALT,sizeof(GMat)},{&n,4}};
    auto vmA=vector<pair<const void*,size_t>>{{&BASALT,sizeof(GMat)},{&n,4}};
    double t=0;int step=0;
    while(t<tend){
        run(Pws,n,{br,bmu,bmv,bmw,bE,bsp},{{&emode,4},{&gam,4},{&BASALT,sizeof(GMat)},{&n,4}});
        float*sp=(float*)bsp->contents();double smax=1e-30;for(uint32_t c=0;c<n;c++)smax=max(smax,(double)sp[c]);
        float dt=CFL*dx/smax;if(t+dt>tend)dt=tend-t; auto dtA=vector<pair<const void*,size_t>>{{&dt,4},{&n,4}};
        auto gdA=vector<pair<const void*,size_t>>{{&BASALT,sizeof(GMat)},{&epsact,4},{&invdx,4},{&dt,4},{&n,4}};
        run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE},lopA);
        if(strn) run(Pstr,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bxy,bxz,byz,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD,bdD},strA);
        run(Prk1,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,br1,bmu1,bmv1,bmw1,bE1},dtA);
        if(strn){ run(Prk1s,n,{bxx,byy,bzz,bxy,bxz,byz,dxx,dyy,dzz,dxy,dxz,dyz,sxx1,syy1,szz1,sxy1,sxz1,syz1,bD,bdD,bD1},dtA);
                  run(Pvm,n,{sxx1,syy1,szz1,sxy1,sxz1,syz1,bD1},vmA); }
        run(Plop,n,{br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},lopA);
        if(strn) run(Pstr,n,{br1,bmu1,bmv1,bmw1,bE1,sxx1,syy1,szz1,sxy1,sxz1,syz1,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD1,bdD},strA);
        run(Prk2,n,{br,bmu,bmv,bmw,bE,br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},dtA);
        if(strn){ run(Prk2s,n,{bxx,byy,bzz,bxy,bxz,byz,sxx1,syy1,szz1,sxy1,sxz1,syz1,dxx,dyy,dzz,dxy,dxz,dyz,bD,bD1,bdD},dtA);
                  run(Pgd,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bD},gdA);
                  run(Pvm,n,{bxx,byy,bzz,bxy,bxz,byz,bD},vmA); }
        t+=dt;step++;
    }
    if(mode=="sod"){
        auto pres=[&](int c){float ke=0.5f*(mu[c]*mu[c]+mv[c]*mv[c]+mw[c]*mw[c])/r[c];return (GAM-1)*(E[c]-ke);};
        double l1r=0,l1p=0,l1u=0;for(int i=0;i<nx;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);l1r+=fabs(r[i]-re);l1p+=fabs(pres(i)-pe);l1u+=fabs(mu[i]/r[i]-ue);}
        printf("GPU Sod L1 vs EXACT: rho=%.4f p=%.4f u=%.4f  GATE(==CPU): %s\n",l1r/nx,l1p/nx,l1u/nx,(l1r/nx<0.007)?"PASS":"CHECK");
    } else if(mode=="surface"){
        double vmax=0;for(int k=0;k<nz;k++)vmax=max(vmax,(double)fabs(mw[k]/r[k]));
        printf("GPU free surface max|v|=%.3f m/s  GATE(static): %s\n",vmax,vmax<5.0?"PASS":"CHECK");
    } else if(mode=="shear"){
        double x0=0.3*Ldom; int ip0=(int)(x0/dx); double pk=0,xpk=0;for(int i=ip0+3;i<nx-2;i++){double vy=mv[i]/r[i];if(vy>pk){pk=vy;xpk=(i+0.5)*dx;}}
        double cs_num=(xpk-x0)/t;
        printf("GPU shear wave: c_s_num=%.0f vs sqrt(G/rho)=%.0f (err %.1f%%) amp=%.3f  GATE: %s\n",cs_num,csb,100*fabs(cs_num-csb)/csb,pk,(fabs(cs_num-csb)/csb<0.03)?"PASS":"CHECK");
    } else if(mode=="freefall"){
        double vexp=-9.8*t,emax=0;for(int k=0;k<nz;k++)emax=max(emax,(double)fabs(mw[k]/r[k]-vexp));
        printf("GPU free-fall: v_z=%.4f expected=%.4f max err=%.3e  GATE: %s\n",mw[nz/2]/r[nz/2],vexp,emax,(emax<1e-3*fabs(vexp))?"PASS":"CHECK");
    } else if(mode=="atmos"){
        double cs=sqrt(GAM*1e5),vmax=0;for(int k=80;k<120;k++)vmax=max(vmax,(double)fabs(mw[k]/r[k]));
        printf("GPU hydrostatic atmos: deep max|v|=%.3f cs=%.0f ratio=%.4f  GATE: %s\n",vmax,cs,vmax/cs,(vmax<0.02*cs)?"PASS":"CHECK");
    } else if(mode=="yield"){
        float*Sxx=(float*)bxx->contents(),*Syy=(float*)byy->contents(),*Szz=(float*)bzz->contents(),*Sxy=pxy,*Sxz=(float*)bxz->contents(),*Syz=(float*)byz->contents();
        double smax=0;for(int i=0;i<nx;i++){double J2=0.5*(Sxx[i]*Sxx[i]+Syy[i]*Syy[i]+Szz[i]*Szz[i])+Sxy[i]*Sxy[i]+Sxz[i]*Sxz[i]+Syz[i]*Syz[i];smax=max(smax,sqrt(3*J2));}
        printf("GPU yield: max sqrt(3J2)=%.4e Y=%.4e ratio=%.4f  GATE: %s\n",smax,Y,smax/Y,(smax<=1.01*Y&&smax>0.9*Y)?"PASS":"CHECK");
    } else if(mode=="tensile"){
        float*Dd=(float*)bD->contents(),*Sxx=(float*)bxx->contents();
        double Dmax=0,Sxxmax=0;for(int i=20;i<nx-20;i++){Dmax=max(Dmax,(double)Dd[i]);Sxxmax=max(Sxxmax,(double)fabs(Sxx[i]));}
        printf("GPU tensile damage: max D=%.4f max|Sxx|=%.3e (Y=%.3e ratio %.4f)  GATE: %s\n",Dmax,Sxxmax,Y,Sxxmax/Y,(Dmax>0.9&&Sxxmax<0.5*Y)?"PASS":"CHECK");
    } else if(mode=="bshock"){
        FILE*o=fopen("bshock_gpu.txt","w");for(int i=0;i<nx;i++)fprintf(o,"%.8e\n",r[i]);fclose(o);printf("GPU bshock wrote bshock_gpu.txt\n");
    } else printf("GPU %s done steps=%d\n",mode.c_str(),step);
    return 0;
}
