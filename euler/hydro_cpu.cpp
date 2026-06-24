// euler M2: CPU reference -- 3D FV hydro (unsplit MUSCL+HLLC+RK2) with a PLUGGABLE
// EOS (ideal gas or Tillotson, via ../sph/eos.hpp), P>=0 floor + low-density ambient
// for a clean free surface. Reconstructs internal energy e (general EOS).
//   ./hydro_cpu sod | sedov         -> ideal-gas gates (regression)
//   ./hydro_cpu surface             -> basalt block + ambient, no gravity, must stay static
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include "../sph/eos.hpp"
using namespace std;
static Material MAT = Material::ideal(1.4);
static double GAM = 1.4;
static inline void eos_pc(double rho,double e,double&p,double&cs){   // (rho, specific internal energy) -> p, sound speed
    p = MAT.pressure(rho, e); if(p<0) p=0;                          // fluid: no tension (strength added in M3)
    cs = MAT.sound_speed(rho, e);
}
struct F5 { double r, mn, mt1, mt2, E; };
// HLLC with general EOS: prim L/R = (rho, vn, t1, t2, e). returns flux (r,mn,mt1,mt2,E).
static F5 hllc(double rL,double vnL,double t1L,double t2L,double eL,
               double rR,double vnR,double t1R,double t2R,double eR){
    double pL,cL,pR,cR; eos_pc(rL,eL,pL,cL); eos_pc(rR,eR,pR,cR);
    double SL=min(vnL-cL,vnR-cR), SR=max(vnL+cL,vnR+cR);
    double Ss=(pR-pL+rL*vnL*(SL-vnL)-rR*vnR*(SR-vnR))/(rL*(SL-vnL)-rR*(SR-vnR));
    double EL=rL*eL+0.5*rL*(vnL*vnL+t1L*t1L+t2L*t2L), ER=rR*eR+0.5*rR*(vnR*vnR+t1R*t1R+t2R*t2R);
    auto Fp=[&](double r,double vn,double t1,double t2,double p,double E){ return F5{r*vn, r*vn*vn+p, r*vn*t1, r*vn*t2, (E+p)*vn}; };
    if(SL>=0) return Fp(rL,vnL,t1L,t2L,pL,EL);
    if(SR<=0) return Fp(rR,vnR,t1R,t2R,pR,ER);
    auto star=[&](double r,double vn,double t1,double t2,double p,double E,double S,F5 Fk){
        double f=r*(S-vn)/(S-Ss);
        F5 Us{f, f*Ss, f*t1, f*t2, f*(E/r+(Ss-vn)*(Ss+p/(r*(S-vn))))};
        F5 U{r, r*vn, r*t1, r*t2, E};
        return F5{Fk.r+S*(Us.r-U.r),Fk.mn+S*(Us.mn-U.mn),Fk.mt1+S*(Us.mt1-U.mt1),Fk.mt2+S*(Us.mt2-U.mt2),Fk.E+S*(Us.E-U.E)}; };
    if(Ss>=0) return star(rL,vnL,t1L,t2L,pL,EL,SL,Fp(rL,vnL,t1L,t2L,pL,EL));
    return star(rR,vnR,t1R,t2R,pR,ER,SR,Fp(rR,vnR,t1R,t2R,pR,ER));
}
static inline double mm(double a,double b){ return (a*b<=0)?0.0:(fabs(a)<fabs(b)?a:b); }
struct Grid { int nx,ny,nz; double dx; vector<double> r,mu,mv,mw,E;
    int idx(int i,int j,int k)const{ return (i*ny+j)*nz+k; }
    Grid(int X,int Y,int Z,double d):nx(X),ny(Y),nz(Z),dx(d),r(X*Y*Z),mu(X*Y*Z),mv(X*Y*Z),mw(X*Y*Z),E(X*Y*Z){} };
static inline double eint(const Grid&g,int c){ double ke=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; return (g.E[c]-ke)/g.r[c]; }
static inline double Pcell(const Grid&g,int c){ double p,cs; eos_pc(g.r[c],eint(g,c),p,cs); return p; }

static void Lop(const Grid&g, vector<double>&dr,vector<double>&dmu,vector<double>&dmv,vector<double>&dmw,vector<double>&dE){
    fill(dr.begin(),dr.end(),0.0);fill(dmu.begin(),dmu.end(),0.0);fill(dmv.begin(),dmv.end(),0.0);fill(dmw.begin(),dmw.end(),0.0);fill(dE.begin(),dE.end(),0.0);
    double invdx=1.0/g.dx;
    auto prim=[&](int c,double*W){ W[0]=g.r[c]; W[1]=g.mu[c]/W[0]; W[2]=g.mv[c]/W[0]; W[3]=g.mw[c]/W[0]; W[4]=eint(g,c); };
    auto sweep=[&](int d){
        int nL=(d==0?g.nx:(d==1?g.ny:g.nz)), nA=(d==0?g.ny:g.nx), nB=(d==2?g.ny:g.nz);
        auto CELL=[&](int p,int a,int b){ return d==0?g.idx(p,a,b):(d==1?g.idx(a,p,b):g.idx(a,b,p)); };
        int in=(d==0?1:(d==1?2:3)), it1=(d==0?2:1), it2=(d==2?2:3);
        vector<double> Fr(nL+1),Fmn(nL+1),Ft1(nL+1),Ft2(nL+1),FE(nL+1);
        for(int a=0;a<nA;a++)for(int b=0;b<nB;b++){
            for(int f=0;f<=nL;f++){
                int iL=max(0,f-1),iR=min(nL-1,f),iLL=max(0,f-2),iRR=min(nL-1,f+1);
                double A[5],B[5],AA[5],BB[5],L[5],R[5];
                prim(CELL(iL,a,b),A);prim(CELL(iR,a,b),B);prim(CELL(iLL,a,b),AA);prim(CELL(iRR,a,b),BB);
                for(int q=0;q<5;q++){ L[q]=A[q]+0.5*mm(A[q]-AA[q],B[q]-A[q]); R[q]=B[q]-0.5*mm(B[q]-A[q],BB[q]-B[q]); }
                if(L[0]<1e-9)L[0]=A[0]; if(R[0]<1e-9)R[0]=B[0]; if(L[4]<0)L[4]=A[4]; if(R[4]<0)R[4]=B[4];
                F5 fl=hllc(L[0],L[in],L[it1],L[it2],L[4], R[0],R[in],R[it1],R[it2],R[4]);
                Fr[f]=fl.r;Fmn[f]=fl.mn;Ft1[f]=fl.mt1;Ft2[f]=fl.mt2;FE[f]=fl.E;
            }
            for(int p=0;p<nL;p++){ int c=CELL(p,a,b);
                dr[c]+=-(Fr[p+1]-Fr[p])*invdx; dE[c]+=-(FE[p+1]-FE[p])*invdx;
                double dmn=-(Fmn[p+1]-Fmn[p])*invdx,dt1=-(Ft1[p+1]-Ft1[p])*invdx,dt2=-(Ft2[p+1]-Ft2[p])*invdx;
                if(d==0){dmu[c]+=dmn;dmv[c]+=dt1;dmw[c]+=dt2;} else if(d==1){dmv[c]+=dmn;dmu[c]+=dt1;dmw[c]+=dt2;} else {dmw[c]+=dmn;dmu[c]+=dt1;dmv[c]+=dt2;}
            }
        }
    };
    sweep(0);sweep(1);sweep(2);
}
static double maxspeed(const Grid&g){ double s=1e-30; int n=g.r.size();
    for(int c=0;c<n;c++){ double p,cs; eos_pc(g.r[c],eint(g,c),p,cs); double v=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; s=max(s,v+cs);} return s; }
static void step_rk2(Grid&g,double dt){ int n=g.r.size(); vector<double> dr(n),dmu(n),dmv(n),dmw(n),dE(n); Grid g1=g;
    Lop(g,dr,dmu,dmv,dmw,dE);
    for(int c=0;c<n;c++){g1.r[c]=g.r[c]+dt*dr[c];g1.mu[c]=g.mu[c]+dt*dmu[c];g1.mv[c]=g.mv[c]+dt*dmv[c];g1.mw[c]=g.mw[c]+dt*dmw[c];g1.E[c]=g.E[c]+dt*dE[c];}
    Lop(g1,dr,dmu,dmv,dmw,dE);
    for(int c=0;c<n;c++){g.r[c]=0.5*(g.r[c]+g1.r[c]+dt*dr[c]);g.mu[c]=0.5*(g.mu[c]+g1.mu[c]+dt*dmu[c]);g.mv[c]=0.5*(g.mv[c]+g1.mv[c]+dt*dmv[c]);g.mw[c]=0.5*(g.mw[c]+g1.mw[c]+dt*dmw[c]);g.E[c]=0.5*(g.E[c]+g1.E[c]+dt*dE[c]);}
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
    if(mode=="sod"){
        MAT=Material::ideal(1.4); int N=200; double dx=1.0/N,tend=0.2,CFL=0.4; Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;double r=x<0.5?1:0.125,p=x<0.5?1:0.1;int c=g.idx(i,0,0);g.r[c]=r;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=p/(GAM-1);}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double l1r=0,l1p=0,l1u=0;for(int i=0;i<N;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);int c=g.idx(i,0,0);l1r+=fabs(g.r[c]-re);l1p+=fabs(Pcell(g,c)-pe);l1u+=fabs(g.mu[c]/g.r[c]-ue);}
        l1r/=N;l1p/=N;l1u/=N;
        printf("Sod (general-EOS) N=%d steps=%d  L1 vs EXACT: rho=%.4f p=%.4f u=%.4f\n",N,s,l1r,l1p,l1u);
        printf("GATE (regression, L1<0.007): %s\n",(l1r<0.007&&l1p<0.007&&l1u<0.009)?"PASS":"CHECK");
    } else if(mode=="sedov"){
        MAT=Material::ideal(1.4); int N=64;double L=2.0,dx=L/N,CFL=0.3,tend=0.5,rho0=1,E0=1,p0=1e-4; Grid g(N,N,N,dx);
        for(int c=0;c<N*N*N;c++){g.r[c]=rho0;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=p0/(GAM-1);} int ic=N/2;g.E[g.idx(ic,ic,ic)]+=E0/(dx*dx*dx);
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double rmax=0,Rsh=0,xc=(ic+0.5)*dx;for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++){int c=g.idx(i,j,k);if(g.r[c]>rmax){rmax=g.r[c];double X=(i+0.5)*dx-xc,Y=(j+0.5)*dx-xc,Z=(k+0.5)*dx-xc;Rsh=sqrt(X*X+Y*Y+Z*Z);}}
        double Ranal=pow(E0*t*t/(0.851*rho0),0.2);
        printf("Sedov N=%d steps=%d  R_num=%.4f R_anal=%.4f (err %.1f%%)\n",N,s,Rsh,Ranal,100*fabs(Rsh-Ranal)/Ranal);
        printf("GATE (R<10%%): %s\n",(fabs(Rsh-Ranal)/Ranal<0.10)?"PASS":"CHECK");
    } else if(mode=="bshock"){ // basalt shock (no analytic) -- CPU vs GPU dynamic-EOS cross-check
        MAT=Material::basalt(); int N=400;double L=400e3,dx=L/N,tend=4.0,CFL=0.4; Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;double rr=x<0.5*L?3000:2700,ee=x<0.5*L?1e6:0;int c=g.idx(i,0,0);g.r[c]=rr;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=rr*ee;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        FILE*o=fopen("bshock_cpu.txt","w");for(int i=0;i<N;i++)fprintf(o,"%.8e\n",g.r[g.idx(i,0,0)]);fclose(o);
        printf("CPU bshock steps=%d t=%.2f  wrote bshock_cpu.txt\n",s,t);
    } else { // surface: basalt block + low-density ambient, NO gravity -> must stay static
        MAT=Material::basalt(); int N=64;double L=200e3,dx=L/N,tend=2.0,CFL=0.4; double rho0=2700,ramb=0.27;
        Grid g(1,1,N,dx);   // 1D column in z
        for(int k=0;k<N;k++){double z=(k+0.5)*dx; double rr=z<0.5*L?rho0:ramb; int c=g.idx(0,0,k); g.r[c]=rr;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=0;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double vmax=0; for(int k=0;k<N;k++){int c=g.idx(0,0,k); vmax=max(vmax,fabs(g.mw[c]/g.r[c]));}
        printf("Free surface (basalt|ambient, no gravity) N=%d steps=%d t=%.2fs  max|v|=%.3f m/s\n",N,s,t,vmax);
        printf("GATE (static, max|v|<5 m/s vs 3145 m/s sound speed): %s\n", vmax<5.0?"PASS":"CHECK");
    }
    return 0;
}
