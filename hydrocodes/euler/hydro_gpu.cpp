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
static GMat AL={2700,7.52e10f,6.5e10f,0.5f,1.63f,5,5,5.0e6f,3.0e6f,1.39e7f,0,0,0};   // aluminum, hydro (Pierazzo)
static float GAM=1.4f;
// host Tillotson pressure (mirrors hydro.metal till) -- for the route-1 hydrostatic reference P0
static float tillP(float rho,float u,const GMat&m){
    if(rho<=0.0f) return 0.0f;
    float eta=rho/m.rho0, mu=eta-1.0f, w0=u/(m.u0*eta*eta)+1.0f;
    float Pc=(m.a+m.b/w0)*rho*u+m.A*mu+m.B*mu*mu;
    if(rho>=m.rho0||u<=m.uiv) return Pc;
    float z=m.rho0/rho-1.0f;
    float Pe=m.a*rho*u+(m.b*rho*u/w0+m.A*mu*expf(-m.beta*z))*expf(-m.alpha*z*z);
    if(u>=m.ucv) return Pe;
    return ((u-m.uiv)*Pe+(m.ucv-u)*Pc)/(m.ucv-m.uiv);
}
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
    int nx,ny,nz,emode; double Ldom,tend,CFL; float gz=0.0f, rcfl=0.0f, rvac=0.0f; int wb=0;
    float cact=0.0f, tdecf=0.0f, pcoh=1.0e6f;   // shock-activated AF: activation coupling, vibration decay time, cohesion floor (0 = AF off)
    if(mode=="sod"){nx=200;ny=nz=1;Ldom=1.0;tend=0.2;CFL=0.4;emode=0;}
    else if(mode=="substrate"){nx=60;ny=1;nz=60;Ldom=30000.0;tend=200.0;CFL=0.4;emode=1;gz=3.71f;wb=1;rcfl=100.0f;}  // route 1: well-balanced gravity-loaded substrate
    else if(mode=="tracer"){nx=200;ny=nz=1;Ldom=1.0;tend=0.15;CFL=0.4;emode=0;}  // M-tag: passive material tracer advection
    else if(mode=="vacuum"){nx=400;ny=nz=1;Ldom=1.0;tend=0.10;CFL=0.4;emode=0;rvac=1e-3f;}  // Toro expansion-into-vacuum
    else if(mode=="sedov"){nx=ny=nz=64;Ldom=2.0;tend=0.5;CFL=0.3;emode=0;}
    else if(mode=="surface"){nx=1;ny=1;nz=64;Ldom=200e3;tend=2.0;CFL=0.4;emode=1;}
    else if(mode=="bshock"){nx=400;ny=nz=1;Ldom=400e3;tend=4.0;CFL=0.4;emode=1;}
    else if(mode=="shear"){nx=800;ny=nz=1;Ldom=400e3;tend=20e3/csb;CFL=0.3;emode=1;}
    else if(mode=="freefall"){nx=ny=1;nz=50;Ldom=500.0;tend=10.0;CFL=0.4;emode=0;gz=9.8f;}
    else if(mode=="atmos"){nx=ny=1;nz=200;double H=1e5/9.8;Ldom=4*H;tend=0.05*Ldom/sqrt(GAM*1e5);CFL=0.4;emode=0;gz=9.8f;}
    else if(mode=="tensile"){nx=100;ny=nz=1;Ldom=10e3;tend=0.1;CFL=0.3;emode=1;}
    else if(mode=="pierazzo"){nx=ny=80;nz=100;Ldom=8.0;tend=1.2e-3;CFL=0.3;emode=1;}   // 3D Al sphere impact (dx=0.1m, a=1m=10cppr)
    else if(mode=="yield"){nx=400;ny=nz=1;Ldom=400e3;tend=15e3/csb;CFL=0.3;emode=1;}
    else if(mode=="af_activate"){nx=200;ny=nz=1;Ldom=100e3;tend=3.0;CFL=0.4;emode=1;cact=0.5f;tdecf=10.0f;}  // shock-activated AF unit test (matches CPU Run A)
    else { fprintf(stderr,"unknown mode '%s' (modes: sod sedov surface bshock shear yield freefall atmos tensile pierazzo vacuum substrate tracer af_activate)\n",mode.c_str()); return 2; }   // fatal: never silently validate the wrong physics
    double dx=Ldom/((nx==1&&ny==1)?nz:nx); uint32_t n=nx*ny*nz; float invdx=1.0f/dx; float gam=GAM;
    GMat MG=(mode=="pierazzo")?AL:BASALT; bool strn=(emode==1 && MG.G>0);   // Al = pure hydro (no strength/damage)
    float epsact=(float)pow(1.0/(1e61*dx*dx*dx),1.0/16.0);   // Weibull weakest-flaw activation strain (host: wk=1e61 overflows FP32)
    D_=MTL::CreateSystemDefaultDevice();
    if(!D_){fprintf(stderr,"no Metal device (this gate needs a Mac GPU)\n");return 1;}
    Q_=D_->newCommandQueue(); NS::Error* e=nullptr;
    string exe=argv[0]; size_t sl=exe.find_last_of('/'); string libp=(sl==string::npos?string("."):exe.substr(0,sl))+"/hydro.metallib";  // resolve relative to the executable, not the CWD
    L_=D_->newLibrary(NS::String::string(libp.c_str(),NS::UTF8StringEncoding),&e);
    if(!L_){ fprintf(stderr,"metallib load failed: %s\n  path: %s\n",e?e->localizedDescription()->utf8String():"(no NS::Error)",libp.c_str()); return 1; }
    auto br=buf(n*4),bmu=buf(n*4),bmv=buf(n*4),bmw=buf(n*4),bE=buf(n*4);
    auto br1=buf(n*4),bmu1=buf(n*4),bmv1=buf(n*4),bmw1=buf(n*4),bE1=buf(n*4);
    auto bdr=buf(n*4),bdmu=buf(n*4),bdmv=buf(n*4),bdmw=buf(n*4),bdE=buf(n*4),bsp=buf(n*4);
    auto bxx=buf(n*4),byy=buf(n*4),bzz=buf(n*4),bxy=buf(n*4),bxz=buf(n*4),byz=buf(n*4);
    auto dxx=buf(n*4),dyy=buf(n*4),dzz=buf(n*4),dxy=buf(n*4),dxz=buf(n*4),dyz=buf(n*4);
    auto sxx1=buf(n*4),syy1=buf(n*4),szz1=buf(n*4),sxy1=buf(n*4),sxz1=buf(n*4),syz1=buf(n*4);
    auto bD=buf(n*4),bdD=buf(n*4),bD1=buf(n*4);   // damage, rate, RK temp
    auto bPmax=buf(n*4);   // peak pressure each cell experiences (M5b Pierazzo decay)
    auto bRR0=buf(n*4),bRP0=buf(n*4);   // route 1: frozen hydrostatic reference (density, pressure)
    auto brc=buf(n*4),brc1=buf(n*4),bdrc=buf(n*4);   // M-tag: passive material tracer rc=rho*c (+ RK temp, derivative)
    auto bvib=buf(n*4),baf=buf(n*4);   // AF: vibrational velocity + fluidization fraction (bPmax reused as the shock-arrival detector)
    float*r=(float*)br->contents(),*mu=(float*)bmu->contents(),*mv=(float*)bmv->contents(),*mw=(float*)bmw->contents(),*E=(float*)bE->contents();
    float*pxy=(float*)bxy->contents(); // for shear/yield diagnostics
    if(mode=="sod"){for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5?1.0f:0.125f,pp=x<0.5?1.0f:0.1f;r[i]=rr;E[i]=pp/(GAM-1);}}
    else if(mode=="sedov"){for(uint32_t c=0;c<n;c++){r[c]=1;E[c]=1e-4/(GAM-1);}int ic=nx/2;E[(ic*ny+ic)*nz+ic]+=1.0/(dx*dx*dx);}
    else if(mode=="surface"){for(int k=0;k<nz;k++){double z=(k+0.5)*dx;r[k]=z<0.5*Ldom?2700.0f:0.27f;E[k]=0;}}
    else if(mode=="substrate"){float A=2.67e10f;for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){double z=(k+0.5)*dx;int c=i*nz+k;if(z<20e3)r[c]=2700.0f*(1.0f+2700.0f*gz*(20e3-z)/A);else r[c]=0.27f;E[c]=0;}}
    else if(mode=="tracer"){float*RC=(float*)brc->contents();double v0=2.0;for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=1.0f;mu[i]=1.0f*v0;E[i]=1.0/(GAM-1)+0.5*1.0*v0*v0;RC[i]=(x>0.3&&x<0.5)?1.0f:0.0f;}}
    else if(mode=="vacuum"){for(int i=0;i<nx;i++){double x=(i+0.5)*dx;if(x<0.5){r[i]=1.0f;E[i]=1.0f/(GAM-1);}else{r[i]=1e-6f;E[i]=1e-6f/(GAM-1);}}}  // material at rest | vacuum
    else if(mode=="bshock"){for(int i=0;i<nx;i++){double x=(i+0.5)*dx;float rr=x<0.5*Ldom?3000.0f:2700.0f,ee=x<0.5*Ldom?1e6f:0.0f;r[i]=rr;E[i]=rr*ee;}}
    else if(mode=="freefall"){for(int k=0;k<nz;k++){r[k]=1.0f;E[k]=1e5/(GAM-1);}}
    else if(mode=="atmos"){double H=1e5/9.8;for(int k=0;k<nz;k++){double z=(k+0.5)*dx;r[k]=exp(-z/H);E[k]=(1e5*exp(-z/H))/(GAM-1);}}
    else if(mode=="tensile"){double rate=1e-2;for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=rho0;double vx=rate*(x-0.5*Ldom);mu[i]=rho0*vx;E[i]=0.5*rho0*vx*vx;}}
    else if(mode=="pierazzo"){double a=1.0,U=10000.0,zsurf=7.0,cx=0.5*Ldom,cy=0.5*Ldom,cz=zsurf+a;   // Al sphere tangent to surface, moving down
        for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){uint32_t c=(i*ny+j)*nz+k;
            double x=(i+0.5)*dx,y=(j+0.5)*dx,z=(k+0.5)*dx; double dr=sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz));
            if(dr<a){r[c]=rho0;mw[c]=-rho0*U;E[c]=0.5*rho0*U*U;}            // projectile
            else if(z<zsurf){r[c]=rho0;E[c]=0;}                            // Al half-space
            else {r[c]=0.27f;E[c]=0;}                                      // low-density ambient
        }}
    else if(mode=="af_activate"){double V=2000.0;for(int i=0;i<nx;i++){r[i]=rho0;E[i]=0;if(i<nx/2)mu[i]=rho0*V;}}  // left half -> right: a planar basalt shock
    else {double x0=0.3*Ldom,wid=8e3,A=(mode=="shear"?1.0:2000.0);for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=rho0;double vy=A*exp(-((x-x0)/wid)*((x-x0)/wid));mv[i]=rho0*vy;E[i]=0.5*rho0*vy*vy;}}
    auto Plop=pso("lop"),Pws=pso("wavespeed"),Prk1=pso("rk1"),Prk2=pso("rk2"),Pstr=pso("strength"),Pvm=pso("vonmises"),Prk1s=pso("rk1s"),Prk2s=pso("rk2s"),Pgd=pso("grow_damage"),Ppm=pso("pmax_update"),Pvoid=pso("voidzero"),Pdamp=pso("damp"),Prk1c=pso("rk1c"),Prk2c=pso("rk2c"),Pupaf=pso("update_af");
    float dampf=1.0f;   // route 1: relaxation damping (1 = off; the void/strength fix makes the substrate stable without it)
    int NX=nx,NY=ny,NZ=nz;
    float*RR0=(float*)bRR0->contents(),*RP0=(float*)bRP0->contents();
    if(wb){ for(uint32_t c=0;c<n;c++){ RR0[c]=r[c]; float p=tillP(r[c],0.0f,MG); RP0[c]=p<0?0:p; } }   // route 1: freeze IC as the hydrostatic reference
    auto lopA=vector<pair<const void*,size_t>>{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&emode,4},{&gam,4},{&MG,sizeof(GMat)},{&n,4},{&gz,4},{&wb,4},{&rvac,4}};
    auto strA=vector<pair<const void*,size_t>>{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&MG,sizeof(GMat)},{&n,4},{&rcfl,4}};
    auto vmA=vector<pair<const void*,size_t>>{{&MG,sizeof(GMat)},{&n,4}};
    double rc0=0, cen0=0;
    if(mode=="tracer"){float*RC=(float*)brc->contents();for(int i=0;i<nx;i++){rc0+=RC[i];cen0+=RC[i]*(i+0.5)*dx;} cen0/=rc0;}
    double t=0;int step=0;
    while(t<tend){
        run(Pws,n,{br,bmu,bmv,bmw,bE,bsp},{{&emode,4},{&gam,4},{&MG,sizeof(GMat)},{&n,4},{&rcfl,4}});
        float*sp=(float*)bsp->contents();double smax=1e-30;for(uint32_t c=0;c<n;c++)smax=max(smax,(double)sp[c]);
        float dt=CFL*dx/smax;if(t+dt>tend)dt=tend-t; auto dtA=vector<pair<const void*,size_t>>{{&dt,4},{&n,4}};
        auto gdA=vector<pair<const void*,size_t>>{{&MG,sizeof(GMat)},{&epsact,4},{&invdx,4},{&dt,4},{&n,4}};
        if(cact>0.0f) run(Pupaf,n,{br,bmu,bmv,bmw,bE,bvib,bPmax,baf},   // shock-activated AF: refresh vib + derive af BEFORE strength reads it (once per step, on the start-of-step state)
            {{&dt,4},{&cact,4},{&tdecf,4},{&pcoh,4},{&emode,4},{&gam,4},{&MG,sizeof(GMat)},{&n,4}});
        run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,bRR0,bRP0,brc,bdrc},lopA);
        if(strn) run(Pstr,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bxy,bxz,byz,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD,bdD},strA);
        run(Prk1,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,br1,bmu1,bmv1,bmw1,bE1},dtA);
        run(Prk1c,n,{brc,bdrc,brc1},dtA);   // tracer predictor: rc1 = rc + dt*drc
        if(strn){ run(Prk1s,n,{bxx,byy,bzz,bxy,bxz,byz,dxx,dyy,dzz,dxy,dxz,dyz,sxx1,syy1,szz1,sxy1,sxz1,syz1,bD,bdD,bD1},dtA);
                  run(Pvm,n,{sxx1,syy1,szz1,sxy1,sxz1,syz1,bD1,baf},vmA); }
        if(rcfl>0) run(Pvoid,n,{bmu1,bmv1,bmw1,br1,bE1,sxx1,syy1,szz1,sxy1,sxz1,syz1,bRR0},{{&rcfl,4},{&n,4}});   // clean predictor void cells before strength reads near-vacuum velocity/stress
        run(Plop,n,{br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE,bRR0,bRP0,brc1,bdrc},lopA);   // predictor reads rc1
        if(strn) run(Pstr,n,{br1,bmu1,bmv1,bmw1,bE1,sxx1,syy1,szz1,sxy1,sxz1,syz1,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD1,bdD},strA);
        run(Prk2,n,{br,bmu,bmv,bmw,bE,br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},dtA);
        run(Prk2c,n,{brc,brc1,bdrc},dtA);   // tracer corrector: rc = 0.5*(rc + rc1 + dt*drc)
        if(strn){ run(Prk2s,n,{bxx,byy,bzz,bxy,bxz,byz,sxx1,syy1,szz1,sxy1,sxz1,syz1,dxx,dyy,dzz,dxy,dxz,dyz,bD,bD1,bdD},dtA);
                  run(Pgd,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bD},gdA);
                  run(Pvm,n,{bxx,byy,bzz,bxy,bxz,byz,bD,baf},vmA); }
        if(mode=="pierazzo") run(Ppm,n,{br,bmu,bmv,bmw,bE,bPmax},{{&MG,sizeof(GMat)},{&n,4}});
        if(dampf<1.0f) run(Pdamp,n,{bmu,bmv,bmw},{{&dampf,4},{&n,4}});   // route 1: quench FP32-noise velocities
        if(rcfl>0) run(Pvoid,n,{bmu,bmv,bmw,br,bE,bxx,byy,bzz,bxy,bxz,byz,bRR0},{{&rcfl,4},{&n,4}});   // route 1: void cells = passive vacuum (reset to reference)
        if(mode=="substrate"){ for(int i=0;i<nx;i++)for(int kk=0;kk<2;kk++){int c=i*nz+kk;r[c]=RR0[c];mu[c]=mv[c]=mw[c]=0;E[c]=0;} }  // pin deep far-field floor to reference
        t+=dt;step++;
    }
    if(mode=="sod"){
        auto pres=[&](int c){float ke=0.5f*(mu[c]*mu[c]+mv[c]*mv[c]+mw[c]*mw[c])/r[c];return (GAM-1)*(E[c]-ke);};
        double l1r=0,l1p=0,l1u=0;for(int i=0;i<nx;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);l1r+=fabs(r[i]-re);l1p+=fabs(pres(i)-pe);l1u+=fabs(mu[i]/r[i]-ue);}
        printf("GPU Sod L1 vs EXACT: rho=%.4f p=%.4f u=%.4f  GATE(==CPU): %s\n",l1r/nx,l1p/nx,l1u/nx,(l1r/nx<0.007)?"PASS":"CHECK");
    } else if(mode=="af_activate"){   // shock-activated AF (matches CPU Run A): shock seeds vib behind the front, none ahead
        float*Vb=(float*)bvib->contents(); double vib_b=0,vib_a=0;
        for(int i=101;i<=120;i++) vib_b=max(vib_b,(double)Vb[i]);
        for(int i=180;i<nx;i++)   vib_a=max(vib_a,(double)Vb[i]);
        printf("GPU AF activate (shock seeding): vib behind front=%.1f m/s, ahead=%.1e m/s  GATE(==CPU, vib_b>1 & vib_a<1e-3): %s\n",
               vib_b,vib_a,(vib_b>1.0&&vib_a<1e-3)?"PASS":"CHECK");
    } else if(mode=="tracer"){
        float*RC=(float*)brc->contents(); double rc1=0,cen1=0,cmin=1e30,cmax=-1e30,v0=2.0;
        for(int i=0;i<nx;i++){double cc=RC[i]/r[i];rc1+=RC[i];cen1+=RC[i]*(i+0.5)*dx;cmin=min(cmin,cc);cmax=max(cmax,cc);}
        cen1/=rc1; double vcen=(cen1-cen0)/t, merr=fabs(rc1-rc0)/rc0, verr=fabs(vcen-v0)/v0;
        printf("GPU tracer: sum(rc) %.6e->%.6e (err %.2e); centroid v=%.4f vs v0=%.1f (err %.2e); c in [%.2e, %.4f]  GATE(==CPU, mass<1e-6 & v<1%% & c in[0,1]): %s\n",
               rc0,rc1,merr,vcen,v0,verr,cmin,cmax,(merr<1e-6&&verr<0.01&&cmin>=-1e-9&&cmax<=1.0+1e-9)?"PASS":"CHECK");
    } else if(mode=="sedov"){
        double E0=1.0,rho0=1.0,cen=0.5*nx*dx,rpk=0,rhomax=0;   // shock radius = radius of peak density
        for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){double rr=r[(i*ny+j)*nz+k];
            if(rr>rhomax){rhomax=rr;double xx=(i+0.5)*dx-cen,yy=(j+0.5)*dx-cen,zz=(k+0.5)*dx-cen;rpk=sqrt(xx*xx+yy*yy+zz*zz);}}
        double Ran=pow(E0/(0.851*rho0),0.2)*pow(t,0.4);   // Sedov 3D, gamma=1.4 (alpha=0.851)
        printf("GPU Sedov 3D: shock R_num=%.3f analytic=%.3f (err %.1f%%) compression=%.2f  GATE(<10%%, resolution-limited): %s\n",rpk,Ran,100*fabs(rpk-Ran)/Ran,rhomax/rho0,(fabs(rpk-Ran)/Ran<0.10)?"PASS":"CHECK");
    } else if(mode=="surface"){
        double vmax=0;for(int k=0;k<nz;k++)vmax=max(vmax,(double)fabs(mw[k]/r[k]));
        printf("GPU free surface max|v|=%.3f m/s  GATE(static): %s\n",vmax,vmax<5.0?"PASS":"CHECK");
    } else if(mode=="vacuum"){
        double cL=sqrt(GAM*1.0/1.0),l1r=0,l1u=0,rmin=1e30;int np=0;   // vs Toro left-rarefaction-into-vacuum (u_L=0)
        for(int i=0;i<nx;i++){double x=(i+0.5)*dx,S=(x-0.5)/t,ra,ua;
            if(S<=-cL){ra=1.0;ua=0;} else if(S<2*cL/(GAM-1)){double u=2.0/(GAM+1)*(cL+S),cc=2.0/(GAM+1)*(cL-0.5*(GAM-1)*S);ra=pow(cc/cL,2.0/(GAM-1));ua=u;} else {ra=0;ua=0;}
            rmin=min(rmin,(double)r[i]); if(ra>0.02){l1r+=fabs(r[i]-ra);l1u+=fabs(mu[i]/r[i]-ua);np++;} }
        l1r/=np;l1u/=np;
        printf("GPU vacuum expansion (Toro): L1 rho=%.4f u=%.4f (npts=%d) min rho=%.2e  GATE(==CPU,L1 rho<0.03&u<0.05&pos): %s\n",l1r,l1u,np,rmin,(l1r<0.03&&l1u<0.05&&rmin>0)?"PASS":"CHECK");
    } else if(mode=="substrate"){
        int nb0=0,nb=0;double vmx=0;for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){double z=(k+0.5)*dx;if(z<20e3)nb0++;}
        for(uint32_t c=0;c<n;c++){if(r[c]>1350){double v=sqrt(mu[c]*mu[c]+mv[c]*mv[c]+mw[c]*mw[c])/r[c];vmx=max(vmx,v);nb++;}}
        printf("GPU substrate (route 1): max|v|=%.3f m/s, basalt cells %d->%d, steps %d  GATE(==CPU,<5 & nb>0.97): %s\n",vmx,nb0,nb,step,(vmx<5.0&&nb>0.97*nb0)?"PASS":"CHECK");
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
    } else if(mode=="pierazzo"){
        float*Pm=(float*)bPmax->contents(); double a=1.0,zsurf=7.0; int ci=nx/2,cj=ny/2,ksurf=(int)(zsurf/dx);
        double Pcore=0; for(int k=ksurf;k>=ksurf-(int)(a/dx) && k>=0;k--){Pcore=max(Pcore,(double)Pm[(ci*ny+cj)*nz+k]);}
        double Sx=0,Sy=0,Sx2=0,Sxy=0; int np=0; FILE*o=fopen("pierazzo_decay.txt","w");
        for(int k=ksurf-1;k>=0;k--){double depth=zsurf-(k+0.5)*dx,ra=depth/a; float P=Pm[(ci*ny+cj)*nz+k];
            if(ra>=0.5&&P>0)fprintf(o,"%.4f %.6e\n",ra,P);
            if(ra>=1.2&&ra<=6.0&&P>0.05*Pcore){double lx=log(ra),ly=log(P);Sx+=lx;Sy+=ly;Sx2+=lx*lx;Sxy+=lx*ly;np++;}}  // developed region only
        fclose(o); double nexp=-(np*Sxy-Sx*Sy)/(np*Sx2-Sx*Sx),Phug=1.637e11;
        printf("Pierazzo Al-sphere (U=10km/s,10cppr): isobaric-core P=%.3e Pa (Hugoniot %.3e, %.0f%%); decay P~(r/a)^-n, n=%.2f (npts=%d)\n",Pcore,Phug,100*Pcore/Phug,nexp,np);
        printf("GATE (core P >0.7*Hugoniot & decay n in [1,3]): %s\n",(Pcore>0.7*Phug&&nexp>1.0&&nexp<3.0)?"PASS":"CHECK");
    } else printf("GPU %s done steps=%d\n",mode.c_str(),step);
    return 0;
}
