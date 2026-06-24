// euler M3: CPU reference -- adds elastic-plastic STRENGTH to the hydro.
// Deviatoric stress S (6 comp): Jaumann rate dS=2G*edev + (S.W - W.S), advected
// (v.grad S), coupled to momentum (div S) + energy (div(S.v)), von Mises radial
// return. Strength is active only when G>0 (ideal gas -> no-op -> Sod/Sedov regress).
//   ./hydro_cpu sod|sedov|surface|bshock   |   shear  (elastic shear-wave speed gate)
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
static double GZ = 0.0;   // gravity (accel toward -z), m/s^2
static double TDEC = 0.0;  // acoustic-fluidization decay time (s); 0 = AF off
static double ETA_AF = 0.0; // acoustic-fluidization viscosity (Pa.s); damps fluidized flow (Wunnemann-Ivanov)
static double RHO_CFL = 0.0; // exclude near-void cells (rho<RHO_CFL) from the CFL (their hi-v dynamics are negligible)
static inline void eos_pc(double rho,double e,double&p,double&cs){ p=MAT.pressure(rho,e); if(p<0)p=0; cs=MAT.sound_speed(rho,e); }
struct F5 { double r, mn, mt1, mt2, E; };
static F5 hllc(double rL,double vnL,double t1L,double t2L,double eL,double rR,double vnR,double t1R,double t2R,double eR){
    double pL,cL,pR,cR; eos_pc(rL,eL,pL,cL); eos_pc(rR,eR,pR,cR);
    double SL=min(vnL-cL,vnR-cR), SR=max(vnL+cL,vnR+cR);
    double Ss=(pR-pL+rL*vnL*(SL-vnL)-rR*vnR*(SR-vnR))/(rL*(SL-vnL)-rR*(SR-vnR));
    double EL=rL*eL+0.5*rL*(vnL*vnL+t1L*t1L+t2L*t2L), ER=rR*eR+0.5*rR*(vnR*vnR+t1R*t1R+t2R*t2R);
    auto Fp=[&](double r,double vn,double t1,double t2,double p,double E){ return F5{r*vn,r*vn*vn+p,r*vn*t1,r*vn*t2,(E+p)*vn}; };
    if(SL>=0) return Fp(rL,vnL,t1L,t2L,pL,EL);
    if(SR<=0) return Fp(rR,vnR,t1R,t2R,pR,ER);
    auto star=[&](double r,double vn,double t1,double t2,double p,double E,double S,F5 Fk){ double f=r*(S-vn)/(S-Ss);
        F5 Us{f,f*Ss,f*t1,f*t2,f*(E/r+(Ss-vn)*(Ss+p/(r*(S-vn))))}; F5 U{r,r*vn,r*t1,r*t2,E};
        return F5{Fk.r+S*(Us.r-U.r),Fk.mn+S*(Us.mn-U.mn),Fk.mt1+S*(Us.mt1-U.mt1),Fk.mt2+S*(Us.mt2-U.mt2),Fk.E+S*(Us.E-U.E)}; };
    if(Ss>=0) return star(rL,vnL,t1L,t2L,pL,EL,SL,Fp(rL,vnL,t1L,t2L,pL,EL));
    return star(rR,vnR,t1R,t2R,pR,ER,SR,Fp(rR,vnR,t1R,t2R,pR,ER));
}
static inline double mm(double a,double b){ return (a*b<=0)?0.0:(fabs(a)<fabs(b)?a:b); }
struct Grid { int nx,ny,nz; double dx; vector<double> r,mu,mv,mw,E, Sxx,Syy,Szz,Sxy,Sxz,Syz, D, af;
    int idx(int i,int j,int k)const{ return (i*ny+j)*nz+k; }
    Grid(int X,int Y,int Z,double d):nx(X),ny(Y),nz(Z),dx(d){ int n=X*Y*Z; r.assign(n,0);mu.assign(n,0);mv.assign(n,0);mw.assign(n,0);E.assign(n,0);
        Sxx.assign(n,0);Syy.assign(n,0);Szz.assign(n,0);Sxy.assign(n,0);Sxz.assign(n,0);Syz.assign(n,0);D.assign(n,0);af.assign(n,0);} };
static inline double eint(const Grid&g,int c){ double ke=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; return (g.E[c]-ke)/g.r[c]; }
static inline double Pcell(const Grid&g,int c){ double p,cs; eos_pc(g.r[c],eint(g,c),p,cs); return p; }

struct DU { vector<double> r,mu,mv,mw,E,Sxx,Syy,Szz,Sxy,Sxz,Syz,D; DU(int n){r.assign(n,0);mu.assign(n,0);mv.assign(n,0);mw.assign(n,0);E.assign(n,0);Sxx.assign(n,0);Syy.assign(n,0);Szz.assign(n,0);Sxy.assign(n,0);Sxz.assign(n,0);Syz.assign(n,0);D.assign(n,0);} };
static void Lop(const Grid&g, DU&d){
    int n=g.nx*g.ny*g.nz; double invdx=1.0/g.dx; double G=MAT.G;
    for(int q=0;q<n;q++){d.r[q]=d.mu[q]=d.mv[q]=d.mw[q]=d.E[q]=0;d.Sxx[q]=d.Syy[q]=d.Szz[q]=d.Sxy[q]=d.Sxz[q]=d.Syz[q]=0;d.D[q]=0;}
    auto prim=[&](int c,double*W){ W[0]=g.r[c]; W[1]=g.mu[c]/W[0]; W[2]=g.mv[c]/W[0]; W[3]=g.mw[c]/W[0]; W[4]=eint(g,c); };
    // ---- hydro: unsplit MUSCL+HLLC (P only) ----
    auto sweep=[&](int dir){
        int nL=(dir==0?g.nx:(dir==1?g.ny:g.nz)),nA=(dir==0?g.ny:g.nx),nB=(dir==2?g.ny:g.nz);
        auto CELL=[&](int p,int a,int b){ return dir==0?g.idx(p,a,b):(dir==1?g.idx(a,p,b):g.idx(a,b,p)); };
        int in=(dir==0?1:(dir==1?2:3)),it1=(dir==0?2:1),it2=(dir==2?2:3);
        vector<double> Fr(nL+1),Fmn(nL+1),Ft1(nL+1),Ft2(nL+1),FE(nL+1);
        for(int a=0;a<nA;a++)for(int b=0;b<nB;b++){
            for(int f=0;f<=nL;f++){ int iL=max(0,f-1),iR=min(nL-1,f),iLL=max(0,f-2),iRR=min(nL-1,f+1);
                double A[5],B[5],AA[5],BB[5],L[5],R[5]; prim(CELL(iL,a,b),A);prim(CELL(iR,a,b),B);prim(CELL(iLL,a,b),AA);prim(CELL(iRR,a,b),BB);
                for(int q=0;q<5;q++){L[q]=A[q]+0.5*mm(A[q]-AA[q],B[q]-A[q]);R[q]=B[q]-0.5*mm(B[q]-A[q],BB[q]-B[q]);}
                if(L[0]<1e-9)L[0]=A[0];if(R[0]<1e-9)R[0]=B[0];if(L[4]<0)L[4]=A[4];if(R[4]<0)R[4]=B[4];
                F5 fl=hllc(L[0],L[in],L[it1],L[it2],L[4],R[0],R[in],R[it1],R[it2],R[4]); Fr[f]=fl.r;Fmn[f]=fl.mn;Ft1[f]=fl.mt1;Ft2[f]=fl.mt2;FE[f]=fl.E; }
            for(int p=0;p<nL;p++){ int c=CELL(p,a,b); d.r[c]+=-(Fr[p+1]-Fr[p])*invdx; d.E[c]+=-(FE[p+1]-FE[p])*invdx;
                double dmn=-(Fmn[p+1]-Fmn[p])*invdx,dt1=-(Ft1[p+1]-Ft1[p])*invdx,dt2=-(Ft2[p+1]-Ft2[p])*invdx;
                if(dir==0){d.mu[c]+=dmn;d.mv[c]+=dt1;d.mw[c]+=dt2;}else if(dir==1){d.mv[c]+=dmn;d.mu[c]+=dt1;d.mw[c]+=dt2;}else{d.mw[c]+=dmn;d.mu[c]+=dt1;d.mv[c]+=dt2;}
            }
        }
    };
    sweep(0);sweep(1);sweep(2);
    if(GZ!=0.0){ for(int c=0;c<n;c++){ d.mw[c]+=-GZ*g.r[c]; d.E[c]+=-GZ*g.mw[c]; } }   // gravity body force
    if(G<=0) return;                      // no strength (ideal gas) -> done
    // ---- strength sources (per cell, central diffs; clamp at boundaries) ----
    auto vx=[&](int c){return g.mu[c]/g.r[c];}; auto vy=[&](int c){return g.mv[c]/g.r[c];}; auto vz=[&](int c){return g.mw[c]/g.r[c];};
    for(int i=0;i<g.nx;i++)for(int j=0;j<g.ny;j++)for(int k=0;k<g.nz;k++){ int c=g.idx(i,j,k);
        int xm=g.idx(max(0,i-1),j,k),xp=g.idx(min(g.nx-1,i+1),j,k),ym=g.idx(i,max(0,j-1),k),yp=g.idx(i,min(g.ny-1,j+1),k),zm=g.idx(i,j,max(0,k-1)),zp=g.idx(i,j,min(g.nz-1,k+1));
        double hx=( (i>0&&i<g.nx-1)?2.0:1.0)*g.dx, hy=((j>0&&j<g.ny-1)?2.0:1.0)*g.dx, hz=((k>0&&k<g.nz-1)?2.0:1.0)*g.dx;
        // velocity gradient L[a][b]=d v_a / d x_b
        double Lxx=(vx(xp)-vx(xm))/hx, Lxy=(vx(yp)-vx(ym))/hy, Lxz=(vx(zp)-vx(zm))/hz;
        double Lyx=(vy(xp)-vy(xm))/hx, Lyy=(vy(yp)-vy(ym))/hy, Lyz=(vy(zp)-vy(zm))/hz;
        double Lzx=(vz(xp)-vz(xm))/hx, Lzy=(vz(yp)-vz(ym))/hy, Lzz=(vz(zp)-vz(zm))/hz;
        double exx=Lxx,eyy=Lyy,ezz=Lzz,exy=0.5*(Lxy+Lyx),exz=0.5*(Lxz+Lzx),eyz=0.5*(Lyz+Lzy),tr=(exx+eyy+ezz)/3.0;
        double Wxy=0.5*(Lxy-Lyx),Wxz=0.5*(Lxz-Lzx),Wyz=0.5*(Lyz-Lzy);   // spin (antisym)
        double sxx=g.Sxx[c],syy=g.Syy[c],szz=g.Szz[c],sxy=g.Sxy[c],sxz=g.Sxz[c],syz=g.Syz[c];
        // Jaumann: dS = 2G edev + (S.W - W.S);  W=[[0,Wxy,Wxz],[-Wxy,0,Wyz],[-Wxz,-Wyz,0]]
        double JSxx=2.0*(sxy*(-Wxy)+sxz*(-Wxz));          // (S.W-W.S)_xx = 2(Sxy*Wyx + Sxz*Wzx)=2(-Sxy Wxy - Sxz Wxz)
        double JSyy=2.0*(sxy*Wxy + syz*(-Wyz));
        double JSzz=2.0*(sxz*Wxz + syz*Wyz);
        double JSxy=sxx*Wxy - syy*Wxy + sxz*(-Wyz) - syz*(-Wxz);  // = (Sxx-Syy)Wxy - Sxz Wyz + Syz Wxz... (recompute below carefully)
        // careful (S.W - W.S) for symmetric S, antisym W: M=S.W-W.S, M symmetric.
        double Mxx=2*(sxy*Wxy*0); // placeholder; compute via matrices
        double Sm[3][3]={{sxx,sxy,sxz},{sxy,syy,syz},{sxz,syz,szz}}, Wm[3][3]={{0,Wxy,Wxz},{-Wxy,0,Wyz},{-Wxz,-Wyz,0}}, Jm[3][3];
        for(int a=0;a<3;a++)for(int bb=0;bb<3;bb++){double s=0;for(int cc=0;cc<3;cc++)s+=Sm[a][cc]*Wm[cc][bb]-Wm[a][cc]*Sm[cc][bb];Jm[a][bb]=s;}
        // advection v.grad S (upwind)
        double u=vx(c),v=vy(c),w=vz(c);
        auto adv=[&](const vector<double>&S){ double ax=u>0?(S[c]-S[xm])/((i>0)?g.dx:1e30):(S[xp]-S[c])/((i<g.nx-1)?g.dx:1e30);
            double ay=v>0?(S[c]-S[ym])/((j>0)?g.dx:1e30):(S[yp]-S[c])/((j<g.ny-1)?g.dx:1e30);
            double az=w>0?(S[c]-S[zm])/((k>0)?g.dx:1e30):(S[zp]-S[c])/((k<g.nz-1)?g.dx:1e30);
            return u*ax+v*ay+w*az; };
        d.Sxx[c]+=2*G*(exx-tr)+Jm[0][0]-adv(g.Sxx);
        d.Syy[c]+=2*G*(eyy-tr)+Jm[1][1]-adv(g.Syy);
        d.Szz[c]+=2*G*(ezz-tr)+Jm[2][2]-adv(g.Szz);
        d.Sxy[c]+=2*G*exy+Jm[0][1]-adv(g.Sxy);
        d.Sxz[c]+=2*G*exz+Jm[0][2]-adv(g.Sxz);
        d.Syz[c]+=2*G*eyz+Jm[1][2]-adv(g.Syz);
        d.D[c]+=-adv(g.D);   // damage advects with the flow (Grady-Kipp growth applied post-step)
        // deviatoric stress -> momentum (div S) + energy (div(S.v)), central diffs
        double dSxx_x=(g.Sxx[xp]-g.Sxx[xm])/hx, dSxy_y=(g.Sxy[yp]-g.Sxy[ym])/hy, dSxz_z=(g.Sxz[zp]-g.Sxz[zm])/hz;
        double dSxy_x=(g.Sxy[xp]-g.Sxy[xm])/hx, dSyy_y=(g.Syy[yp]-g.Syy[ym])/hy, dSyz_z=(g.Syz[zp]-g.Syz[zm])/hz;
        double dSxz_x=(g.Sxz[xp]-g.Sxz[xm])/hx, dSyz_y=(g.Syz[yp]-g.Syz[ym])/hy, dSzz_z=(g.Szz[zp]-g.Szz[zm])/hz;
        d.mu[c]+=dSxx_x+dSxy_y+dSxz_z; d.mv[c]+=dSxy_x+dSyy_y+dSyz_z; d.mw[c]+=dSxz_x+dSyz_y+dSzz_z;
        if(g.af[c]>0&&ETA_AF>0){ double eta=g.af[c]*ETA_AF, id2=1.0/(g.dx*g.dx);   // AF Newtonian viscosity (damps fluidized flow)
            d.mu[c]+=eta*(vx(xp)+vx(xm)+vx(yp)+vx(ym)+vx(zp)+vx(zm)-6*vx(c))*id2;
            d.mv[c]+=eta*(vy(xp)+vy(xm)+vy(yp)+vy(ym)+vy(zp)+vy(zm)-6*vy(c))*id2;
            d.mw[c]+=eta*(vz(xp)+vz(xm)+vz(yp)+vz(ym)+vz(zp)+vz(zm)-6*vz(c))*id2; }
        // energy: div(S.v). (S.v)_x = Sxx u+Sxy v+Sxz w, etc.
        auto Svx=[&](int cc){return g.Sxx[cc]*vx(cc)+g.Sxy[cc]*vy(cc)+g.Sxz[cc]*vz(cc);};
        auto Svy=[&](int cc){return g.Sxy[cc]*vx(cc)+g.Syy[cc]*vy(cc)+g.Syz[cc]*vz(cc);};
        auto Svz=[&](int cc){return g.Sxz[cc]*vx(cc)+g.Syz[cc]*vy(cc)+g.Szz[cc]*vz(cc);};
        d.E[c]+=(Svx(xp)-Svx(xm))/hx+(Svy(yp)-Svy(ym))/hy+(Svz(zp)-Svz(zm))/hz;
    }
}
static void vonmises(Grid&g){ double Y0=MAT.Y; if(Y0<=0)return; int n=g.nx*g.ny*g.nz;
    for(int c=0;c<n;c++){ double Y=(1.0-g.D[c])*(1.0-g.af[c])*Y0;   // damage + acoustic fluidization degrade shear strength
        double sxx=g.Sxx[c],syy=g.Syy[c],szz=g.Szz[c],sxy=g.Sxy[c],sxz=g.Sxz[c],syz=g.Syz[c];
        double J2=0.5*(sxx*sxx+syy*syy+szz*szz)+sxy*sxy+sxz*sxz+syz*syz; double vm=sqrt(3.0*J2);
        if(vm>Y){double f=(vm>0?Y/vm:0.0); g.Sxx[c]*=f;g.Syy[c]*=f;g.Szz[c]*=f;g.Sxy[c]*=f;g.Sxz[c]*=f;g.Syz[c]*=f;} } }
static void grow_damage(Grid&g,double dt){ double Em=MAT.Emod,wk=MAT.wk,wm=MAT.wm; if(Em<=0||wk<=0)return;
    double dx=g.dx,eps_act=pow(1.0/(wk*dx*dx*dx),1.0/wm); int n=g.nx*g.ny*g.nz;   // weakest-flaw activation strain
    for(int c=0;c<n;c++){ if(g.D[c]>=1.0)continue; double P=Pcell(g,c);
        double sigmax=max(-P+g.Sxx[c],max(-P+g.Syy[c],-P+g.Szz[c]));   // max tensile principal stress (diagonal proxy)
        if(sigmax>0.0 && sigmax/Em>eps_act){ double cg=0.4*sqrt(Em/max(g.r[c],1e-30)),Rs=0.5*dx;
            double d13=pow(g.D[c],1.0/3.0)+(cg/Rs)*dt; g.D[c]=min(1.0,d13*d13*d13); } } }
static double maxspeed(const Grid&g){ double s=1e-30; int n=g.r.size(); double G=MAT.G;
    for(int c=0;c<n;c++){ if(g.r[c]<RHO_CFL) continue; double p,cs; eos_pc(g.r[c],eint(g,c),p,cs); double cel=sqrt(cs*cs+ (G>0?(4.0/3.0)*G/g.r[c]:0.0));
        double v=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; s=max(s,v+cel);} return s; }
static void step_rk2(Grid&g,double dt){ int n=g.r.size(); DU d(n); Grid g1=g;
    Lop(g,d);
    for(int c=0;c<n;c++){g1.r[c]=g.r[c]+dt*d.r[c];g1.mu[c]=g.mu[c]+dt*d.mu[c];g1.mv[c]=g.mv[c]+dt*d.mv[c];g1.mw[c]=g.mw[c]+dt*d.mw[c];g1.E[c]=g.E[c]+dt*d.E[c];
        g1.Sxx[c]=g.Sxx[c]+dt*d.Sxx[c];g1.Syy[c]=g.Syy[c]+dt*d.Syy[c];g1.Szz[c]=g.Szz[c]+dt*d.Szz[c];g1.Sxy[c]=g.Sxy[c]+dt*d.Sxy[c];g1.Sxz[c]=g.Sxz[c]+dt*d.Sxz[c];g1.Syz[c]=g.Syz[c]+dt*d.Syz[c];g1.D[c]=g.D[c]+dt*d.D[c];}
    vonmises(g1); Lop(g1,d);
    for(int c=0;c<n;c++){g.r[c]=0.5*(g.r[c]+g1.r[c]+dt*d.r[c]);g.mu[c]=0.5*(g.mu[c]+g1.mu[c]+dt*d.mu[c]);g.mv[c]=0.5*(g.mv[c]+g1.mv[c]+dt*d.mv[c]);g.mw[c]=0.5*(g.mw[c]+g1.mw[c]+dt*d.mw[c]);g.E[c]=0.5*(g.E[c]+g1.E[c]+dt*d.E[c]);
        g.Sxx[c]=0.5*(g.Sxx[c]+g1.Sxx[c]+dt*d.Sxx[c]);g.Syy[c]=0.5*(g.Syy[c]+g1.Syy[c]+dt*d.Syy[c]);g.Szz[c]=0.5*(g.Szz[c]+g1.Szz[c]+dt*d.Szz[c]);g.Sxy[c]=0.5*(g.Sxy[c]+g1.Sxy[c]+dt*d.Sxy[c]);g.Sxz[c]=0.5*(g.Sxz[c]+g1.Sxz[c]+dt*d.Sxz[c]);g.Syz[c]=0.5*(g.Syz[c]+g1.Syz[c]+dt*d.Syz[c]);g.D[c]=0.5*(g.D[c]+g1.D[c]+dt*d.D[c]);}
    grow_damage(g,dt); vonmises(g);
    if(TDEC>0){ double f=exp(-dt/TDEC); int N=g.r.size(); for(int c=0;c<N;c++) g.af[c]*=f; }   // AF vibrations decay
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
// analytic Tillotson Hugoniot: shock into material at rest, particle velocity up behind shock.
// solve Rankine-Hugoniot (mass+momentum+energy) coupled to P=till(rho,u) by bisection on shock speed Us.
static double hugoniot_P(double up){
    double rho0=MAT.rho0, lo=up*1.0001, hi=up*30.0;
    for(int it=0;it<200;it++){ double Us=0.5*(lo+hi);
        double rho=rho0*Us/(Us-up), Pmom=rho0*Us*up, uint=0.5*Pmom*(1.0/rho0-1.0/rho);
        double Peos=MAT.pressure(rho,uint); if(Pmom-Peos<0) lo=Us; else hi=Us; }
    double Us=0.5*(lo+hi); return rho0*Us*up;
}
int main(int argc,char**argv){
    string mode=argc>1?argv[1]:"sod";
    if(mode=="sod"){
        MAT=Material::ideal(1.4);int N=200;double dx=1.0/N,tend=0.2,CFL=0.4;Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;double r=x<0.5?1:0.125,p=x<0.5?1:0.1;int c=g.idx(i,0,0);g.r[c]=r;g.E[c]=p/(GAM-1);}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double l1r=0,l1p=0,l1u=0;for(int i=0;i<N;i++){double x=(i+0.5)*dx,re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);int c=g.idx(i,0,0);l1r+=fabs(g.r[c]-re);l1p+=fabs(Pcell(g,c)-pe);l1u+=fabs(g.mu[c]/g.r[c]-ue);}
        printf("Sod (M3 regression) L1: rho=%.4f p=%.4f u=%.4f  GATE: %s\n",l1r/N,l1p/N,l1u/N,(l1r/N<0.007)?"PASS":"CHECK");
    } else if(mode=="shear"){
        // small-amplitude elastic shear pulse v_y(x); must propagate at c_s=sqrt(G/rho0)
        MAT=Material::basalt(); int N=800;double L=400e3,dx=L/N,CFL=0.3; double rho0=2700,G=MAT.G,cs=sqrt(G/rho0);
        Grid g(N,1,1,dx); double x0=0.3*L,wid=8e3,A=1.0;
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;int c=g.idx(i,0,0);g.r[c]=rho0;double vy=A*exp(-((x-x0)/wid)*((x-x0)/wid));g.mv[c]=rho0*vy;g.E[c]=0.5*rho0*vy*vy;}
        double tend=20e3/cs;   // travel ~20 km
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        // right-going pulse peak position (search x>x0)
        int ip0=(int)(x0/dx);double pk=0,xpk=0;for(int i=ip0+3;i<N-2;i++){int c=g.idx(i,0,0);double vy=g.mv[c]/g.r[c];if(vy>pk){pk=vy;xpk=(i+0.5)*dx;}}
        double cs_num=(xpk-x0)/t;
        FILE*o=fopen("shear_cpu.txt","w");for(int i=0;i<N;i++)fprintf(o,"%.2f %.6e\n",(i+0.5)*dx,g.mv[g.idx(i,0,0)]/g.r[g.idx(i,0,0)]);fclose(o);
        printf("Shear wave: c_s_num=%.0f m/s  analytic sqrt(G/rho)=%.0f  (err %.1f%%)  peak amp=%.3f\n",cs_num,cs,100*fabs(cs_num-cs)/cs,pk);
        printf("GATE (elastic shear speed, <3%%): %s\n",(fabs(cs_num-cs)/cs<0.03)?"PASS":"CHECK");
    } else if(mode=="yield"){
        // large-amplitude shear -> deviatoric stress must cap at the von Mises surface Y/sqrt(3)
        MAT=Material::basalt(); int N=400;double L=400e3,dx=L/N,CFL=0.3; double rho0=2700,G=MAT.G,Y=MAT.Y,cs=sqrt(G/rho0);
        Grid g(N,1,1,dx); double x0=0.3*L,wid=8e3,A=2000.0;   // A big -> stress drives past yield
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;int c=g.idx(i,0,0);g.r[c]=rho0;double vy=A*exp(-((x-x0)/wid)*((x-x0)/wid));g.mv[c]=rho0*vy;g.E[c]=0.5*rho0*vy*vy;}
        double tend=15e3/cs;double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double smax=0;for(int i=0;i<N;i++){int c=g.idx(i,0,0);double J2=0.5*(g.Sxx[c]*g.Sxx[c]+g.Syy[c]*g.Syy[c]+g.Szz[c]*g.Szz[c])+g.Sxy[c]*g.Sxy[c]+g.Sxz[c]*g.Sxz[c]+g.Syz[c]*g.Syz[c];smax=max(smax,sqrt(3*J2));}
        printf("Yield: max sqrt(3 J2)=%.4e  Y=%.4e  ratio=%.4f\n",smax,Y,smax/Y);
        printf("GATE (von Mises cap, max<=1.01*Y & yield reached >0.9*Y): %s\n",(smax<=1.01*Y&&smax>0.9*Y)?"PASS":"CHECK");
    } else if(mode=="bshock"){
        MAT=Material::basalt();int N=400;double L=400e3,dx=L/N,tend=4.0,CFL=0.4;Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;double rr=x<0.5*L?3000:2700,ee=x<0.5*L?1e6:0;int c=g.idx(i,0,0);g.r[c]=rr;g.E[c]=rr*ee;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        FILE*o=fopen("bshock_cpu.txt","w");for(int i=0;i<N;i++)fprintf(o,"%.8e\n",g.r[g.idx(i,0,0)]);fclose(o);printf("CPU bshock done\n");
    } else if(mode=="freefall"){ // uniform medium under gravity -> v_z = -g*t exactly
        MAT=Material::ideal(1.4); GZ=9.8; int N=50;double dx=10.0,tend=10.0,CFL=0.4;Grid g(1,1,N,dx);
        for(int k=0;k<N;k++){int c=g.idx(0,0,k);g.r[c]=1.0;g.E[c]=1e5/(GAM-1);}
        double t=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;}
        double vexp=-GZ*t,emax=0;for(int k=0;k<N;k++){int c=g.idx(0,0,k);emax=max(emax,fabs(g.mw[c]/g.r[c]-vexp));}
        printf("Free-fall: v_z=%.4f expected -g*t=%.4f  max err=%.3e\n",g.mw[g.idx(0,0,N/2)]/g.r[g.idx(0,0,N/2)],vexp,emax);
        printf("GATE (v=-g*t, err<1e-3*|v|): %s\n",(emax<1e-3*fabs(vexp))?"PASS":"CHECK");
    } else if(mode=="atmos"){ // ideal-gas isothermal hydrostatic atmosphere -> must stay ~static
        MAT=Material::ideal(1.4); GZ=9.8; double P0=1e5,rho0=1.0,H=P0/(rho0*GZ),cs=sqrt(GAM*P0/rho0);
        int N=200;double L=4*H,dx=L/N,CFL=0.4; Grid g(1,1,N,dx);
        for(int k=0;k<N;k++){double z=(k+0.5)*dx;int c=g.idx(0,0,k);g.r[c]=rho0*exp(-z/H);g.E[c]=(P0*exp(-z/H))/(GAM-1);}
        double tend=0.05*L/cs,t=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;}
        double vmax=0;for(int k=80;k<120;k++){int c=g.idx(0,0,k);vmax=max(vmax,fabs(g.mw[c]/g.r[c]));}  // deep interior, pre-boundary-contamination
        printf("Hydrostatic atmosphere (interior balance): deep max|v|=%.3f m/s  cs=%.0f  ratio=%.4f  t=%.2fs\n",vmax,cs,vmax/cs,t);
        printf("GATE (interior well-balanced, max|v|<0.02*cs): %s\n",(vmax<0.02*cs)?"PASS":"CHECK");
    } else if(mode=="tensile"){ // stretch a basalt block -> tensile damage grows -> shear strength degrades
        MAT=Material::basalt(); int N=100;double L=10e3,dx=L/N,tend=0.1,CFL=0.3; double rho0=2700,rate=1e-2,Y=MAT.Y;
        Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;int c=g.idx(i,0,0);g.r[c]=rho0;double vx=rate*(x-0.5*L);g.mu[c]=rho0*vx;g.E[c]=0.5*rho0*vx*vx;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double Dmax=0,Sxxmax=0;for(int i=20;i<N-20;i++){int c=g.idx(i,0,0);Dmax=max(Dmax,g.D[c]);Sxxmax=max(Sxxmax,fabs(g.Sxx[c]));}
        printf("Tensile damage: max D=%.4f  max|Sxx|=%.3e (Y=%.3e, ratio %.4f)\n",Dmax,Sxxmax,Y,Sxxmax/Y);
        printf("GATE (D>0.9 grows & strength degrades Sxx<0.5Y): %s\n",(Dmax>0.9&&Sxxmax<0.5*Y)?"PASS":"CHECK");
    } else if(mode=="alimpact"){ // Pierazzo-style: 1D Al-on-Al planar impact -> peak shock P vs analytic Hugoniot
        MAT=Material::aluminum(); int N=400;double L=2.0,dx=L/N,CFL=0.3,tend=3e-5; double rho0=2700,up=5000;  // U=10 km/s symmetric (small domain: Hugoniot is scale-free)
        Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;int c=g.idx(i,0,0);g.r[c]=rho0;double vx=(x<0.5*L?up:-up);g.mu[c]=rho0*vx;g.E[c]=0.5*rho0*vx*vx;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double Pmax=0;for(int i=0;i<N;i++)Pmax=max(Pmax,Pcell(g,g.idx(i,0,0)));
        double Phug=hugoniot_P(up);
        printf("Al-on-Al planar (U=%.0f km/s, up=%.0f m/s): peak P=%.3e Pa  Tillotson Hugoniot=%.3e (err %.1f%%)\n",2*up/1000.0,up,Pmax,Phug,100*fabs(Pmax-Phug)/Phug);
        printf("GATE (shock P matches analytic Hugoniot, <5%%): %s\n",(fabs(Pmax-Phug)/Phug<0.05)?"PASS":"CHECK");
    } else if(mode=="collapse"){ // acoustic-fluidization demo: a basalt step slumps if fluidized, holds if not
        MAT=Material::basalt(); GZ=3.71; double rho0=2700,A=2.67e10; int NXc=80,NZc=40;double dx=500.0,CFL=0.4;
        double zlo=10e3,h0=4e3,zhi=zlo+h0,tend=8.0; double hfin[2]; const char*lbl[2]={"AF-on ","AF-off"};
        for(int run=0;run<2;run++){
            TDEC=(run==0?1e6:0.0); ETA_AF=(run==0?1e9:0.0); RHO_CFL=100.0; Grid g(NXc,1,NZc,dx);   // exclude only the true ambient (0.27)
            for(int i=0;i<NXc;i++)for(int k=0;k<NZc;k++){double x=(i+0.5)*dx,z=(k+0.5)*dx;int c=g.idx(i,0,k);double zs=(x<0.5*NXc*dx?zhi:zlo);
                if(z<zs){g.r[c]=rho0*(1.0+rho0*GZ*(zs-z)/A); if(run==0)g.af[c]=1.0;} else g.r[c]=0.27; g.E[c]=0;}   // cold lithostatic step
            double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);
                for(int i=0;i<NXc;i++)for(int k=0;k<2;k++){double zs=(((i+0.5)*dx)<0.5*NXc*dx?zhi:zlo),z=(k+0.5)*dx;int c=g.idx(i,0,k);
                    g.r[c]=rho0*(1.0+rho0*GZ*(zs-z)/A);g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=0;g.Sxx[c]=g.Syy[c]=g.Szz[c]=g.Sxy[c]=g.Sxz[c]=g.Syz[c]=0;}  // held floor
                t+=dt;s++;}
            double smx=0,smn=1e30,vmx=0;for(int i=0;i<NXc;i++)for(int k=0;k<NZc;k++){int c=g.idx(i,0,k);if(g.r[c]>1350)vmx=max(vmx,sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]);}
            for(int i=0;i<NXc;i++){double sf=0;for(int k=NZc-1;k>=0;k--){if(g.r[g.idx(i,0,k)]>1350){sf=(k+0.5)*dx;break;}}smx=max(smx,sf);smn=min(smn,sf);}
            hfin[run]=smx-smn; printf("  %s: step height=%.2f km, max|v|=%.1f m/s (steps %d, t=%.1f)\n",lbl[run],hfin[run]/1e3,vmx,s,t);
        }
        printf("AF collapse: h0=%.1f km -> AF-on %.2f km, AF-off %.2f km\n",h0/1e3,hfin[0]/1e3,hfin[1]/1e3);
        printf("GATE (AF slumps step <0.5*h0 & no-AF holds >0.7*h0): %s\n",(hfin[0]<0.5*h0&&hfin[1]>0.7*h0)?"PASS":"CHECK");
    } else { // surface
        MAT=Material::basalt();int N=64;double L=200e3,dx=L/N,tend=2.0,CFL=0.4;Grid g(1,1,N,dx);
        for(int k=0;k<N;k++){double z=(k+0.5)*dx;double rr=z<0.5*L?2700:0.27;int c=g.idx(0,0,k);g.r[c]=rr;g.E[c]=0;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double vmax=0;for(int k=0;k<N;k++){int c=g.idx(0,0,k);vmax=max(vmax,fabs(g.mw[c]/g.r[c]));}
        printf("Free surface max|v|=%.3f m/s  GATE: %s\n",vmax,vmax<5.0?"PASS":"CHECK");
    }
    return 0;
}
