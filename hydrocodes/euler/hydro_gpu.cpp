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
struct SParams { int nx,ny,nz; float invdx; uint32_t n; float rcfl; int axisym; float eta_af; };   // mirrors hydro.metal SParams (strength kernel scalars)
struct RockP { int rock; float Yi0, mui, mud, yd0, Ym; };   // mirrors hydro.metal RockP (pressure-dependent yield; rock=0 -> cohesion-only)
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
    float etaaf=0.0f;   // AF Newtonian viscosity coefficient (Pa.s); 0 = AF viscosity off
    int axisym=0;   // cylindrical (r,z) axisymmetric geometry (x->r, axis at r=0); 0 = Cartesian
    // Phase-3 crater CLI (mirrors hydro_cpu crater): a U TDEC ETA g cppr tend Rfac Zfac profile rock Yd0
    double cr_a=500,cr_U=12000,cr_TDEC=0,cr_ETA=0,cr_g=3.71,cr_cppr=5,cr_targ=-1,cr_Rfac=18,cr_Zfac=22,cr_Yd0=0; int cr_rock=1; double cr_zsurf=0; string cr_prof="crater_profile.txt";
    if(mode=="crater"){
        if(argc>2)cr_a=atof(argv[2]); if(argc>3)cr_U=atof(argv[3]); if(argc>4)cr_TDEC=atof(argv[4]); if(argc>5)cr_ETA=atof(argv[5]);
        if(argc>6)cr_g=atof(argv[6]); if(argc>7)cr_cppr=atof(argv[7]); if(argc>8)cr_targ=atof(argv[8]); if(argc>9)cr_Rfac=atof(argv[9]);
        if(argc>10)cr_Zfac=atof(argv[10]); if(argc>11)cr_prof=argv[11]; if(argc>12)cr_rock=atoi(argv[12]); if(argc>13)cr_Yd0=atof(argv[13]);
    }
    double cr_tauto=2.0*sqrt((cr_Rfac*cr_a)/cr_g);   // gravity formation/collapse timescale (settling cap = 10x)
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
    else if(mode=="sedov_axi"){nx=64;ny=1;nz=128;Ldom=2.0;tend=0.5;CFL=0.3;emode=0;axisym=1;}  // axisymmetric (r,z) on-axis point blast = 3D spherical Sedov
    else if(mode=="lame"){nx=200;ny=nz=1;Ldom=200*500.0;tend=0.0;CFL=0.4;emode=1;axisym=1;}  // Phase-2b: cylindrical strength gate (single-eval hoop strain + Lame equilibrium); tend=0 -> handled before the main loop
    else if(mode=="af_visc"){nx=120;ny=nz=1;Ldom=120*500.0;tend=0.0;CFL=0.4;emode=1;axisym=1;etaaf=1e9f;}  // Phase-2b item 5: cylindrical AF viscosity (single-eval af=1 vs af=0); pre-loop block
    else if(mode=="vib_advect"){nx=200;ny=nz=1;Ldom=200*500.0;tend=10.0;CFL=0.4;emode=1;}  // Phase-2b item 6: vib advects with the flow (uniform v0 -> bump translates); cact=etaaf=0 -> pure advection
    else if(mode=="friction"){nx=8;ny=nz=1;Ldom=8.0;tend=0.0;CFL=0.4;emode=1;}  // Phase-3: ROCK pressure-dependent yield gate (single vonmises eval; pre-loop block)
    else if(mode=="crater"){ double dxc=cr_a/cr_cppr; nx=(int)(cr_Rfac*cr_a/dxc+0.5); ny=1; nz=(int)(cr_Zfac*cr_a/dxc+0.5);   // Phase-3: axisymmetric vertical impact crater (mirrors hydro_cpu crater)
        Ldom=nx*dxc; CFL=0.4; emode=1; axisym=1; wb=1; gz=(float)cr_g; rcfl=100.0f; rvac=100.0f;
        cact=(cr_TDEC>0?0.5f:0.0f); tdecf=(float)cr_TDEC; etaaf=(float)cr_ETA; pcoh=1.0e6f;
        tend=(cr_targ>0?cr_targ:10.0*cr_tauto); }   // explicit tend>0 -> fixed run; else run-to-settling capped at 10x t_auto
    else { fprintf(stderr,"unknown mode '%s' (modes: sod sedov surface bshock shear yield freefall atmos tensile pierazzo vacuum substrate tracer af_activate sedov_axi lame af_visc vib_advect friction crater)\n",mode.c_str()); return 2; }   // fatal: never silently validate the wrong physics
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
    auto bvib1=buf(n*4),bdvib=buf(n*4);   // AF: vib RK predictor temp + advection rate (item 6)
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
    else if(mode=="vib_advect"){double v0=2000.0,x0=0.3*Ldom,wid=8e3;float*Vb=(float*)bvib->contents();for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=rho0;mu[i]=rho0*v0;E[i]=0.5*rho0*v0*v0;Vb[i]=exp(-((x-x0)/wid)*((x-x0)/wid));}}  // uniform v0 flow + a vib bump
    else if(mode=="sedov_axi"){for(uint32_t c=0;c<n;c++){r[c]=1.0f;E[c]=1e-4/(GAM-1);} int kc=nz/2; E[kc]+=1.0/(3.141592653589793*dx*dx*dx);}  // on-axis point blast (cell i=0,k=kc -> linear index kc); 3D volume = pi*dx^3
    else if(mode=="crater"){ double Abulk=BASALT.A; double Zdom=nz*dx, above=5.0*cr_a; cr_zsurf=Zdom-above;   // lithostatic basalt below the free surface, low-density ambient above (impactor added AFTER the WB reference is frozen)
        for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){ double z=(k+0.5)*dx; uint32_t c=i*nz+k;
            if(z<cr_zsurf) r[c]=(float)(rho0*(1.0+rho0*gz*(cr_zsurf-z)/Abulk)); else r[c]=0.27f; E[c]=0; } }
    else {double x0=0.3*Ldom,wid=8e3,A=(mode=="shear"?1.0:2000.0);for(int i=0;i<nx;i++){double x=(i+0.5)*dx;r[i]=rho0;double vy=A*exp(-((x-x0)/wid)*((x-x0)/wid));mv[i]=rho0*vy;E[i]=0.5*rho0*vy*vy;}}
    auto Plop=pso("lop"),Pws=pso("wavespeed"),Prk1=pso("rk1"),Prk2=pso("rk2"),Pstr=pso("strength"),Pvm=pso("vonmises"),Prk1s=pso("rk1s"),Prk2s=pso("rk2s"),Pgd=pso("grow_damage"),Ppm=pso("pmax_update"),Pvoid=pso("voidzero"),Pdamp=pso("damp"),Prk1c=pso("rk1c"),Prk2c=pso("rk2c"),Pupaf=pso("update_af");
    float dampf=1.0f;   // route 1: relaxation damping (1 = off; the void/strength fix makes the substrate stable without it)
    int NX=nx,NY=ny,NZ=nz;
    float*RR0=(float*)bRR0->contents(),*RP0=(float*)bRP0->contents();
    if(wb){ for(uint32_t c=0;c<n;c++){ RR0[c]=r[c]; float p=tillP(r[c],0.0f,MG); RP0[c]=p<0?0:p; } }   // route 1: freeze IC as the hydrostatic reference
    if(mode=="crater"){ double zc=cr_zsurf+cr_a;   // impactor sphere on the axis, tangent to the surface, moving down (added after the WB reference so the reference is impactor-free)
        for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){ double rr=(i+0.5)*dx, z=(k+0.5)*dx; uint32_t c=i*nz+k;
            if(sqrt(rr*rr+(z-zc)*(z-zc))<cr_a){ r[c]=(float)rho0; mw[c]=(float)(-rho0*cr_U); E[c]=(float)(0.5*rho0*cr_U*cr_U); } } }
    auto lopA=vector<pair<const void*,size_t>>{{&NX,4},{&NY,4},{&NZ,4},{&invdx,4},{&emode,4},{&gam,4},{&MG,sizeof(GMat)},{&n,4},{&gz,4},{&wb,4},{&rvac,4},{&axisym,4}};
    SParams sp{NX,NY,NZ,invdx,n,rcfl,axisym,etaaf};   // strength scalars packed into one buffer (Metal bind-point cap)
    auto strA=vector<pair<const void*,size_t>>{{&MG,sizeof(GMat)},{&sp,sizeof(SParams)}};
    RockP rp{0,1.0e7f,1.2f,0.6f,0.0f,2.5e9f};   // default rock=0 (cohesion-only); crater/friction modes enable + set
    if(mode=="crater" && cr_rock) rp = RockP{1, 1.0e7f, 1.2f, 0.6f, (float)cr_Yd0, 2.5e9f};   // ROCK pressure-dependent yield (mirrors CPU globals Y0/mu_i/mu_d/Y_d0/Y_m)
    auto vmA=vector<pair<const void*,size_t>>{{&MG,sizeof(GMat)},{&n,4},{&rp,sizeof(RockP)}};
    if(mode=="lame"){   // Phase-2b GPU: cylindrical strength gate -- two single-evaluation sub-checks, mirrors the CPU oracle (GPU==CPU)
        float Gm=MG.G; float dxf=(float)dx; bool A=false,B=false;
        auto evalLopStr=[&](){   // one lop + one strength evaluation (no time integration); reads dS / dmu back
            run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,bRR0,bRP0,brc,bdrc},lopA);
            run(Pstr,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bxy,bxz,byz,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD,bdD,baf,bvib,bdvib},strA); };
        // (A) HOOP STRAIN RATE: u(r)=edot*r, zero stress -> dS driven by e_thth=u/r (geometric)
        float edot=1e-6f;
        for(uint32_t c=0;c<n;c++){ r[c]=MG.rho0; E[c]=0; mu[c]=MG.rho0*edot*(((float)c+0.5f)*dxf); mv[c]=mw[c]=0; }
        for(auto b:{bxx,byy,bzz,bxy,bxz,byz,bdmu,bdmv,bdmw,bdE}) memset(b->contents(),0,n*4);
        evalLopStr();
        float*DSxx=(float*)dxx->contents(),*DSyy=(float*)dyy->contents(),*DSzz=(float*)dzz->contents();
        int c=100; float trc=2.0f*edot/3.0f, an_rr=2*Gm*(edot-trc), an_th=2*Gm*(edot-trc), an_zz=2*Gm*(0.0f-trc);
        float erra=fabs(DSxx[c]-an_rr)+fabs(DSyy[c]-an_th)+fabs(DSzz[c]-an_zz), sca=fabs(an_rr)+fabs(an_zz)+1e-30f;
        A=(erra/sca<1e-3);
        printf("GPU Lame (A hoop): dS_rr=%.4e/%.4e  dS_thth=%.4e/%.4e  dS_zz=%.4e/%.4e (num/analytic, rel err %.2e)\n",
            DSxx[c],an_rr,DSyy[c],an_th,DSzz[c],an_zz,erra/sca);
        // (B) LAME STATIC EQUILIBRIUM: S_rr=-Bc/r^2, S_thth=+Bc/r^2, S_zz=0 at rest, P=0 -> radial-momentum residual ~ truncation
        int ia=20,ib=180; float ia_r=(ia+0.5f)*dxf, Bc=2.0e6f*ia_r*ia_r;
        float*Sx=(float*)bxx->contents(),*Sy=(float*)byy->contents(),*Sz=(float*)bzz->contents();
        for(uint32_t c2=0;c2<n;c2++){ r[c2]=MG.rho0; E[c2]=0; mu[c2]=mv[c2]=mw[c2]=0; float rr=((float)c2+0.5f)*dxf;
            if((int)c2>=ia&&(int)c2<=ib){ Sx[c2]=-Bc/(rr*rr); Sy[c2]=Bc/(rr*rr); Sz[c2]=0; } else { Sx[c2]=Sy[c2]=Sz[c2]=0; } }
        for(auto b:{bxy,bxz,byz,bdmu,bdmv,bdmw,bdE}) memset(b->contents(),0,n*4);
        evalLopStr();
        float*DMU=(float*)bdmu->contents();
        double resmax=0,termmax=0,baremax=0;
        for(int i=ia+2;i<=ib-2;i++){ double rr=(i+0.5)*dxf;
            resmax=max(resmax,(double)fabs(DMU[i]));
            termmax=max(termmax,(double)fabs((Sx[i]-Sy[i])/rr));
            baremax=max(baremax,(double)fabs(DMU[i]-(Sx[i]-Sy[i])/rr)); }
        B=(resmax<1e-2*termmax);
        printf("GPU Lame (B equilib): max|d.mu|=%.3e  geom-term scale=%.3e (ratio %.2e)  vs no-geom residual=%.3e\n",resmax,termmax,resmax/termmax,baremax);
        printf("GATE(==CPU, cylindrical strength: hoop strain rate exact & Lame equilibrium balances): %s\n",(A&&B)?"PASS":"CHECK");
        return 0;
    }
    if(mode=="af_visc"){   // Phase-2b item 5 GPU: cylindrical AF viscosity. Difference af=1 vs af=0 isolates the af-dependent viscous term (hydro/pressure cancel). GPU==CPU.
        float dxf=(float)dx, A=1e-9f; int ic=nx/2;
        for(uint32_t c=0;c<n;c++){ r[c]=MG.rho0; E[c]=0; mv[c]=mw[c]=0; float rr=((float)c+0.5f)*dxf; mu[c]=MG.rho0*A*rr*rr; }  // u(r)=A*r^2
        for(auto b:{bxx,byy,bzz,bxy,bxz,byz,bvib}) memset(b->contents(),0,n*4);
        float*AF=(float*)baf->contents();
        auto eval=[&](float afv)->double{   // one lop+strength eval with af=afv everywhere; return d.mu at the interior cell
            for(uint32_t c=0;c<n;c++) AF[c]=afv;
            for(auto b:{bdmu,bdmv,bdmw,bdE}) memset(b->contents(),0,n*4);
            run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,bRR0,bRP0,brc,bdrc},lopA);
            run(Pstr,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bxy,bxz,byz,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD,bdD,baf,bvib,bdvib},strA);
            return (double)((float*)bdmu->contents())[ic]; };
        double dmu1=eval(1.0f), dmu0=eval(0.0f), visc=dmu1-dmu0;
        double an_cyl=(double)etaaf*3.0*A, an_cart=(double)etaaf*2.0*A;
        bool P=(fabs(visc-an_cyl)<1e-3*fabs(an_cyl));
        printf("GPU AF visc (cylindrical): isolated viscous d.mu=%.6e  analytic cyl 3A=%.6e (Cartesian-only 2A=%.6e)  rel err %.2e\n",visc,an_cyl,an_cart,fabs(visc-an_cyl)/fabs(an_cyl));
        printf("GATE(==CPU, cylindrical AF viscosity (lap v)_r = 3A not 2A): %s\n",P?"PASS":"CHECK");
        return 0;
    }
    if(mode=="friction"){   // Phase-3 GPU: ROCK pressure-dependent yield. 4 pressures x {intact,damaged}; capped stress must follow Y_i(P)/Y_d(P). GPU==CPU.
        rp = RockP{1, 1.0e7f, 1.2f, 0.6f, 5.0e6f, 2.5e9f};   // rock on; Y_d0=5e6 to exercise the damaged-cohesion term
        float muset[4]={0.001f,0.005f,0.02f,0.05f};
        float*Sxy=(float*)bxy->contents(); float*Dd=(float*)bD->contents();
        for(int j=0;j<4;j++)for(int dam=0;dam<2;dam++){ int c=j*2+dam;
            r[c]=BASALT.rho0*(1.0f+muset[j]); E[c]=0; mu[c]=mv[c]=mw[c]=0; Dd[c]=(float)dam; Sxy[c]=5.0e9f; }
        for(auto b:{bxx,bzz,bxz,byz,byy,baf}) memset(b->contents(),0,n*4);   // only Sxy nonzero; af=0
        run(Pvm,n,{bxx,byy,bzz,bxy,bxz,byz,bD,baf,br,bmu,bmv,bmw,bE},vmA);
        auto tillP=[&](float rho){ float eta=rho/BASALT.rho0,m=eta-1.0f; return BASALT.A*m+BASALT.B*m*m; };  // e=0 condensed
        bool pass=true;
        printf("GPU ROCK yield: Y_I0=%.2e MU_I=%.2f MU_D=%.2f Y_D0=%.2e Y_M=%.2e\n",rp.Yi0,rp.mui,rp.mud,rp.yd0,rp.Ym);
        for(int j=0;j<4;j++)for(int dam=0;dam<2;dam++){ int c=j*2+dam;
            double P=tillP(BASALT.rho0*(1.0f+muset[j]));
            double vm=sqrt(3.0)*fabs(Sxy[c]);
            double Yi=rp.Yi0+rp.mui*P/(1.0+rp.mui*P/(rp.Ym-rp.Yi0)), Yd=fmin(rp.yd0+rp.mud*P,rp.Ym), Yan=(dam?Yd:Yi);
            double err=fabs(vm-Yan)/Yan; if(err>1e-3) pass=false;
            printf("  P=%.3e D=%d: capped vm=%.4e  analytic Y=%.4e (%s, err %.1e)\n",P,dam,vm,Yan,dam?"damaged friction":"intact Lundborg",err); }
        printf("GATE(==CPU, ROCK yield: capped stress follows Y_i(P)/Y_d(P)): %s\n",pass?"PASS":"CHECK");
        return 0;
    }
    double rc0=0, cen0=0;
    if(mode=="tracer"){float*RC=(float*)brc->contents();for(int i=0;i<nx;i++){rc0+=RC[i];cen0+=RC[i]*(i+0.5)*dx;} cen0/=rc0;}
    double vs0=0,vcen0=0;
    if(mode=="vib_advect"){float*Vb=(float*)bvib->contents();for(int i=0;i<nx;i++){vs0+=Vb[i];vcen0+=Vb[i]*(i+0.5)*dx;} vcen0/=vs0;}
    // crater run-to-settling state + host-side surface diagnostics (read live from the shared buffers)
    auto crsurf=[&](int i)->double{ for(int k=nz-1;k>=0;k--){ if(r[i*nz+k]>1350.0f) return (k+0.5)*dx; } return 0.0; };
    auto crvexc=[&](){ double V=0; for(int i=0;i<nx;i++){ double d=cr_zsurf-crsurf(i); if(d>0) V+=d*(i+0.5); } return V; };   // r-weighted excavated cross-section ~ crater volume (continuous settling signal)
    auto crvmax=[&](){ double v=0; for(uint32_t c=0;c<n;c++) if(r[c]>1350.0f){ double vv=sqrt((double)mu[c]*mu[c]+(double)mv[c]*mv[c]+(double)mw[c]*mw[c])/r[c]; v=max(v,vv);} return v; };
    double cr_tolM=0.02, cr_Wwin=2.0*cr_tauto, cr_tsettled=-1.0, cr_nextlog=cr_tauto, cr_dmaxT=0.0;   // window-mean drift settling (mirrors CPU)
    double cr_sumV=0.0, cr_meanPrev=-1.0, cr_winStart=-1.0; long cr_cntV=0;
    double t=0;int step=0;
    while(t<tend){
        run(Pws,n,{br,bmu,bmv,bmw,bE,bsp},{{&emode,4},{&gam,4},{&MG,sizeof(GMat)},{&n,4},{&rcfl,4}});
        float*sp=(float*)bsp->contents();double smax=1e-30;for(uint32_t c=0;c<n;c++)smax=max(smax,(double)sp[c]);
        float dt=CFL*dx/smax;if(t+dt>tend)dt=tend-t; auto dtA=vector<pair<const void*,size_t>>{{&dt,4},{&n,4}};
        auto gdA=vector<pair<const void*,size_t>>{{&MG,sizeof(GMat)},{&epsact,4},{&invdx,4},{&dt,4},{&n,4}};
        if(cact>0.0f) run(Pupaf,n,{br,bmu,bmv,bmw,bE,bvib,bPmax,baf},   // shock-activated AF: refresh vib + derive af BEFORE strength reads it (once per step, on the start-of-step state)
            {{&dt,4},{&cact,4},{&tdecf,4},{&pcoh,4},{&emode,4},{&gam,4},{&MG,sizeof(GMat)},{&n,4}});
        run(Plop,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,bRR0,bRP0,brc,bdrc},lopA);
        if(strn) run(Pstr,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bxy,bxz,byz,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD,bdD,baf,bvib,bdvib},strA);
        run(Prk1,n,{br,bmu,bmv,bmw,bE,bdr,bdmu,bdmv,bdmw,bdE,br1,bmu1,bmv1,bmw1,bE1},dtA);
        run(Prk1c,n,{brc,bdrc,brc1},dtA);   // tracer predictor: rc1 = rc + dt*drc
        if(strn){ run(Prk1c,n,{bvib,bdvib,bvib1},dtA);   // vib predictor: vib1 = vib + dt*dvib (reuses the generic scalar-RK kernel)
                  run(Prk1s,n,{bxx,byy,bzz,bxy,bxz,byz,dxx,dyy,dzz,dxy,dxz,dyz,sxx1,syy1,szz1,sxy1,sxz1,syz1,bD,bdD,bD1},dtA);
                  run(Pvm,n,{sxx1,syy1,szz1,sxy1,sxz1,syz1,bD1,baf,br1,bmu1,bmv1,bmw1,bE1},vmA); }   // predictor state for ROCK pressure
        if(rcfl>0) run(Pvoid,n,{bmu1,bmv1,bmw1,br1,bE1,sxx1,syy1,szz1,sxy1,sxz1,syz1,bRR0},{{&rcfl,4},{&n,4}});   // clean predictor void cells before strength reads near-vacuum velocity/stress
        run(Plop,n,{br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE,bRR0,bRP0,brc1,bdrc},lopA);   // predictor reads rc1
        if(strn) run(Pstr,n,{br1,bmu1,bmv1,bmw1,bE1,sxx1,syy1,szz1,sxy1,sxz1,syz1,bdmu,bdmv,bdmw,bdE,dxx,dyy,dzz,dxy,dxz,dyz,bD1,bdD,baf,bvib1,bdvib},strA);   // corrector reads predictor state (vib1; af held through the step)
        run(Prk2,n,{br,bmu,bmv,bmw,bE,br1,bmu1,bmv1,bmw1,bE1,bdr,bdmu,bdmv,bdmw,bdE},dtA);
        run(Prk2c,n,{brc,brc1,bdrc},dtA);   // tracer corrector: rc = 0.5*(rc + rc1 + dt*drc)
        if(strn){ run(Prk2c,n,{bvib,bvib1,bdvib},dtA);   // vib corrector: vib = 0.5*(vib + vib1 + dt*dvib)
                  run(Prk2s,n,{bxx,byy,bzz,bxy,bxz,byz,sxx1,syy1,szz1,sxy1,sxz1,syz1,dxx,dyy,dzz,dxy,dxz,dyz,bD,bD1,bdD},dtA);
                  run(Pgd,n,{br,bmu,bmv,bmw,bE,bxx,byy,bzz,bD},gdA);
                  run(Pvm,n,{bxx,byy,bzz,bxy,bxz,byz,bD,baf,br,bmu,bmv,bmw,bE},vmA); }   // corrector state for ROCK pressure
        if(mode=="pierazzo") run(Ppm,n,{br,bmu,bmv,bmw,bE,bPmax},{{&MG,sizeof(GMat)},{&n,4}});
        if(dampf<1.0f) run(Pdamp,n,{bmu,bmv,bmw},{{&dampf,4},{&n,4}});   // route 1: quench FP32-noise velocities
        if(rcfl>0) run(Pvoid,n,{bmu,bmv,bmw,br,bE,bxx,byy,bzz,bxy,bxz,byz,bRR0},{{&rcfl,4},{&n,4}});   // route 1: void cells = passive vacuum (reset to reference)
        if(mode=="substrate"||mode=="crater"){ for(int i=0;i<nx;i++)for(int kk=0;kk<2;kk++){int c=i*nz+kk;r[c]=RR0[c];mu[c]=mv[c]=mw[c]=0;E[c]=0;} }  // pin deep far-field floor to reference
        t+=dt;step++;
        if(mode=="crater"){
            if(step%20==0){ double dn=0; for(int i=0;i<nx;i++) dn=max(dn,cr_zsurf-crsurf(i)); cr_dmaxT=max(cr_dmaxT,dn); }   // track transient excavation depth
            if(t>=cr_nextlog){ double dn=0; for(int i=0;i<nx;i++) dn=max(dn,cr_zsurf-crsurf(i));   // progress to stderr (mirrors CPU)
                fprintf(stderr,"  [GPU crater a=%.0f] t=%.1f/%.0fs steps=%d max|v|=%.2f Vexc=%.4e depth=%.2fkm meanPrev=%.4e\n",cr_a,t,tend,step,crvmax(),crvexc(),dn/1e3,cr_meanPrev); cr_nextlog+=20.0; }
            if(cr_targ<=0 && t>cr_tauto && step%5==0){ double V=crvexc();   // run-to-settling: compare successive window-means of the excavated volume
                if(cr_winStart<0)cr_winStart=t; cr_sumV+=V; cr_cntV++;
                if(t-cr_winStart>=cr_Wwin){ double meanNow=cr_sumV/cr_cntV;
                    if(cr_meanPrev>0 && fabs(meanNow-cr_meanPrev)<cr_tolM*meanNow){ cr_tsettled=t; break; }
                    cr_meanPrev=meanNow; cr_sumV=0; cr_cntV=0; cr_winStart=t; } }
        }
    }
    if(mode=="sod"){
        auto pres=[&](int c){float ke=0.5f*(mu[c]*mu[c]+mv[c]*mv[c]+mw[c]*mw[c])/r[c];return (GAM-1)*(E[c]-ke);};
        double l1r=0,l1p=0,l1u=0;for(int i=0;i<nx;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);l1r+=fabs(r[i]-re);l1p+=fabs(pres(i)-pe);l1u+=fabs(mu[i]/r[i]-ue);}
        printf("GPU Sod L1 vs EXACT: rho=%.4f p=%.4f u=%.4f  GATE(==CPU): %s\n",l1r/nx,l1p/nx,l1u/nx,(l1r/nx<0.007)?"PASS":"CHECK");
    } else if(mode=="crater"){   // Phase-3 GPU crater diagnostic (mirrors hydro_cpu): floor depth, apparent diameter, rim, d/D + profile dump
        double vmx=0,pmx=0; bool fin=true;
        for(uint32_t c=0;c<n;c++){ if(!isfinite(r[c])||!isfinite(E[c])) fin=false;
            if(r[c]>1350.0f){ double v=sqrt((double)mu[c]*mu[c]+(double)mv[c]*mv[c]+(double)mw[c]*mw[c])/r[c]; vmx=max(vmx,v);
                double u=(E[c]-0.5*((double)mu[c]*mu[c]+(double)mv[c]*mv[c]+(double)mw[c]*mw[c])/r[c])/r[c]; pmx=max(pmx,(double)tillP(r[c],(float)u,MG)); } }
        double z0=crsurf(nx-1);
        double zfloor=z0; int ifloor=0; for(int i=0;i<nx;i++){double zs=crsurf(i); if(zs<zfloor){zfloor=zs;ifloor=i;}}
        int iD=nx-1; for(int i=ifloor;i<nx;i++){ if(crsurf(i)>=z0-0.25*dx){ iD=i; break; } }
        double rD=(iD+0.5)*dx, Dapp=2.0*rD, dapp=z0-zfloor, dD=(Dapp>0?dapp/Dapp:0.0);
        double zrim=z0,rrim=0; for(int i=0;i<nx;i++){double zs=crsurf(i); if(zs>zrim){zrim=zs;rrim=(i+0.5)*dx;}}
        printf("GPU CRATER a=%.0fm U=%.0f g=%.2f TDEC=%.3g ETA=%.3g | grid %dx%d dx=%.0f zsurf=%.1fkm tend=%.1fs steps=%d\n",cr_a,cr_U,cr_g,cr_TDEC,cr_ETA,nx,nz,dx,cr_zsurf/1e3,tend,step);
        printf("  settling: %s at t=%.1fs (t_auto=%.1fs cap=%.1fs tolM=%.1f%% Wwin=%.1fs final max|v|=%.1f m/s)\n",(cr_tsettled>0?"SETTLED":(cr_targ>0?"FIXED-RUN":"NOT-SETTLED@cap")),(cr_tsettled>0?cr_tsettled:t),cr_tauto,tend,cr_tolM*100,cr_Wwin,vmx);
        printf("  D_app=%.2f km  depth=%.2f km  transient depth=%.2f km  rim uplift=%.2f km @ r=%.2f km  d/D=%.3f\n",Dapp/1e3,dapp/1e3,cr_dmaxT/1e3,(zrim-z0)/1e3,rrim/1e3,dD);
        printf("  stability: max|v|=%.1f m/s  maxP=%.2e Pa  %s\n",vmx,pmx,fin?"FINITE":"NONFINITE-BLOWUP");
        printf("RESULT %.1f %.4f %.4f %.4f %.4f\n",cr_a,Dapp,dapp,cr_dmaxT,dD);
        FILE*pf=fopen(cr_prof.c_str(),"w"); if(pf){ fprintf(pf,"# r_m  zsurf-z0_m   (GPU a=%.0f U=%.0f g=%.2f TDEC=%.3g ETA=%.3g dx=%.0f z0=%.1f)\n",cr_a,cr_U,cr_g,cr_TDEC,cr_ETA,dx,z0);
            for(int i=0;i<nx;i++) fprintf(pf,"%.1f %.3f\n",(i+0.5)*dx, crsurf(i)-z0); fclose(pf); }
    } else if(mode=="sedov_axi"){   // axisymmetric geometry check: on-axis point blast = 3D spherical Sedov (same analytic as the 3D gate)
        double rpk=0,rhomax=0,zc=(nz/2+0.5)*dx;
        for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){int c=i*nz+k;double rr=r[c];
            if(rr>rhomax){rhomax=rr;double rad=(i+0.5)*dx,zz=(k+0.5)*dx-zc;rpk=sqrt(rad*rad+zz*zz);}}
        double Ran=pow(1.0/0.851,0.2)*pow(tend,0.4);
        printf("GPU Sedov axisym (on-axis = 3D spherical): R_num=%.3f analytic=%.3f (err %.1f%%), compression=%.2f  GATE(==CPU,<10%%): %s\n",
               rpk,Ran,100*fabs(rpk-Ran)/Ran,rhomax,(fabs(rpk-Ran)/Ran<0.10)?"PASS":"CHECK");
    } else if(mode=="af_activate"){   // shock-activated AF (matches CPU Run A): shock seeds vib behind the front, none ahead
        float*Vb=(float*)bvib->contents(); double vib_b=0,vib_a=0;
        for(int i=101;i<=120;i++) vib_b=max(vib_b,(double)Vb[i]);
        for(int i=180;i<nx;i++)   vib_a=max(vib_a,(double)Vb[i]);
        printf("GPU AF activate (shock seeding): vib behind front=%.1f m/s, ahead=%.1e m/s  GATE(==CPU, vib_b>1 & vib_a<1e-3): %s\n",
               vib_b,vib_a,(vib_b>1.0&&vib_a<1e-3)?"PASS":"CHECK");
    } else if(mode=="vib_advect"){   // Phase-2b item 6: vib advects with the flow -> bump translates at v0, sum conserved
        float*Vb=(float*)bvib->contents(); double s1=0,cen1=0,vmax=0,v0=2000.0;
        for(int i=0;i<nx;i++){s1+=Vb[i];cen1+=Vb[i]*(i+0.5)*dx;vmax=max(vmax,(double)Vb[i]);}
        cen1/=s1; double vcen=(cen1-vcen0)/t,merr=fabs(s1-vs0)/vs0,verr=fabs(vcen-v0)/v0;
        printf("GPU vib advect: sum(vib) %.4f->%.4f (err %.2e); centroid v=%.1f vs v0=%.0f (err %.2e); peak=%.3f  GATE(==CPU, sum<1%% & v<2%% & peak>0): %s\n",
               vs0,s1,merr,vcen,v0,verr,vmax,(merr<0.01&&verr<0.02&&vmax>0)?"PASS":"CHECK");
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
