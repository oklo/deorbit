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
#include "../common/eos.hpp"
using namespace std;
static Material MAT = Material::ideal(1.4);
static double GAM = 1.4;
static double GZ = 0.0;   // gravity (accel toward -z), m/s^2
static double TDEC = 0.0;  // acoustic-fluidization vibration decay time (s); 0 = AF off
static double ETA_AF = 0.0; // acoustic-fluidization viscosity (Pa.s); damps fluidized flow (Wunnemann-Ivanov)
static double C_ACT = 0.0;  // AF shock-activation coupling: vib seeded at shock = C_ACT*|post-shock v|; 0 = shock activation off (af stays at IC)
static double P_COH = 1.0e6; // cohesion floor (Pa) for the overburden in the fluidization ratio (avoids /0 at the free surface)
static double P_ACT = 0.0;  // AF activation: only a SHOCK-sized pressure jump (new-peak rise > P_ACT in a step) seeds vib; 0 = seed on any new peak (legacy). Stops slow collapse flow from re-fluidizing Eulerian wall cells.
static double RHO_CFL = 0.0; // exclude near-void cells (rho<RHO_CFL) from the CFL (their hi-v dynamics are negligible)
static double DAMP = 1.0;    // relaxation velocity damping per step (route 2: settle a gravity-loaded substrate); 1 = off
static double RHO_VAC = 0.0; // vacuum-aware Riemann flux: a face side with rho<RHO_VAC is treated as vacuum (free surface); 0 = off
static bool AXISYM = false; // cylindrical (r,z) axisymmetric geometry: x->r (radial, axis at r=0), z->axial, y unused (ny=1)
// ROCK pressure-dependent yield (Lundborg intact + cohesionless damaged friction; iSALE, Collins/Melosh/Ivanov 2004).
// false => cohesion-only von Mises Y=(1-D)(1-af)*MAT.Y (back-compat; all M3-M5 strength gates keep ROCK off).
static bool ROCK = false;
static double Y_I0 = 1.0e7;  // intact cohesion (Pa)
static double MU_I = 1.2;    // intact coefficient of internal friction
static double MU_D = 0.6;    // damaged (comminuted) friction coefficient
static double Y_D0 = 0.0;    // damaged residual cohesion (Pa); 0 = cohesionless friction (iSALE breccia ~ a few MPa)
static double Y_M  = 2.5e9;  // limiting (von Mises) strength at high pressure (Pa)
static double Y_BINGHAM = 0.0; // AF Bingham floor (Pa): fluidized cells keep this flow resistance (option B, 2026-07-05); 0 = legacy (1-af) cut
static bool VISC_IMP = false;  // Option C (2026-07-07): implicit AF viscous update -- eta_eff unbounded by the explicit dt limit; replaces the explicit term in Lop when on
static vector<double> REF_R0, REF_P0; // route 1 (Audusse): frozen hydrostatic reference (density,pressure) per cell; empty => WB off
static inline void eos_pc(double rho,double e,double&p,double&cs){ p=MAT.pressure(rho,e); if(p<0)p=0; cs=MAT.sound_speed(rho,e); }
struct F5 { double r, mn, mt1, mt2, E; };
// vacuum-aware flux: a (near-)vacuum side has ~no mass -> any pressure force gives it runaway velocity.
// Replace HLLC there with the exact rarefaction-into-vacuum solution sampled at the cell face (x/t=0).
// Effective adiabatic exponent g=rho*c^2/p (=gamma for ideal gas, exact; ->large = stiff -> minimal expansion).
static F5 vac_flux(double rM,double vnM,double t1M,double t2M,double eM,int side){ // side=+1: material is LEFT, vacuum right; -1: material RIGHT
    double pM,cM; eos_pc(rM,eM,pM,cM);
    double g=rM*cM*cM/max(pM,1e-30); if(g<1.01)g=1.01; if(g>1e5)g=1e5;
    double s=(double)side, vn=s*vnM;                 // map to "material expanding to the right" frame
    if(vn-cM>=0){ double E=rM*eM+0.5*rM*(vnM*vnM+t1M*t1M+t2M*t2M); return F5{rM*vnM,rM*vnM*vnM+pM,rM*vnM*t1M,rM*vnM*t2M,(E+pM)*vnM}; } // supersonic outflow: whole state
    if(vn+2*cM/(g-1)<=0) return F5{0,0,0,0,0};        // face already inside vacuum
    double u0=2.0/(g+1)*(cM+0.5*(g-1)*vn), c0=u0;     // sample the rarefaction fan at x/t=0 (u=c there)
    double r0=rM*pow(c0/cM,2.0/(g-1)), p0=pM*pow(r0/rM,g), e0=p0/((g-1)*r0), un0=s*u0;  // map back
    double E0=r0*e0+0.5*r0*(un0*un0+t1M*t1M+t2M*t2M);
    return F5{r0*un0, r0*un0*un0+p0, r0*un0*t1M, r0*un0*t2M, (E0+p0)*un0};
}
static F5 hllc(double rL,double vnL,double t1L,double t2L,double eL,double rR,double vnR,double t1R,double t2R,double eR){
    if(RHO_VAC>0){ bool vL=(rL<RHO_VAC), vR=(rR<RHO_VAC);
        if(vL&&vR) return F5{0,0,0,0,0};                                      // vacuum | vacuum
        if(vR) return vac_flux(rL,vnL,t1L,t2L,eL,+1);                          // material L | vacuum R
        if(vL) return vac_flux(rR,vnR,t1R,t2R,eR,-1);                          // vacuum L | material R
    }
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
// pressure-aware HLLC: caller supplies the (well-balanced, deviation-reconstructed) face pressures
// pL,pR directly; sound speeds come from (rho,e). Used only in the z-sweep under route-1 WB.
static F5 hllc_p(double rL,double vnL,double t1L,double t2L,double eL,double pL,double rR,double vnR,double t1R,double t2R,double eR,double pR){
    double cL=MAT.sound_speed(rL,eL), cR=MAT.sound_speed(rR,eR); if(pL<0)pL=0; if(pR<0)pR=0;
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
struct Grid { int nx,ny,nz; double dx; vector<double> r,mu,mv,mw,E, Sxx,Syy,Szz,Sxy,Sxz,Syz, D, af, rc, vib, Pmax;  // rc=rho*c tracer; vib=AF vibrational velocity; Pmax=peak pressure (shock-arrival detector)
    int idx(int i,int j,int k)const{ return (i*ny+j)*nz+k; }
    Grid(int X,int Y,int Z,double d):nx(X),ny(Y),nz(Z),dx(d){ int n=X*Y*Z; r.assign(n,0);mu.assign(n,0);mv.assign(n,0);mw.assign(n,0);E.assign(n,0);
        Sxx.assign(n,0);Syy.assign(n,0);Szz.assign(n,0);Sxy.assign(n,0);Sxz.assign(n,0);Syz.assign(n,0);D.assign(n,0);af.assign(n,0);rc.assign(n,0);vib.assign(n,0);Pmax.assign(n,0);} };
static inline double eint(const Grid&g,int c){ double ke=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; return (g.E[c]-ke)/g.r[c]; }
static inline double Pcell(const Grid&g,int c){ double p,cs; eos_pc(g.r[c],eint(g,c),p,cs); return p; }

struct DU { vector<double> r,mu,mv,mw,E,Sxx,Syy,Szz,Sxy,Sxz,Syz,D,rc,vib; DU(int n){r.assign(n,0);mu.assign(n,0);mv.assign(n,0);mw.assign(n,0);E.assign(n,0);Sxx.assign(n,0);Syy.assign(n,0);Szz.assign(n,0);Sxy.assign(n,0);Sxz.assign(n,0);Syz.assign(n,0);D.assign(n,0);rc.assign(n,0);vib.assign(n,0);} };
static void Lop(const Grid&g, DU&d){
    int n=g.nx*g.ny*g.nz; double invdx=1.0/g.dx; double G=MAT.G;
    for(int q=0;q<n;q++){d.r[q]=d.mu[q]=d.mv[q]=d.mw[q]=d.E[q]=0;d.Sxx[q]=d.Syy[q]=d.Szz[q]=d.Sxy[q]=d.Sxz[q]=d.Syz[q]=0;d.D[q]=0;d.rc[q]=0;d.vib[q]=0;}
    auto prim=[&](int c,double*W){ W[0]=g.r[c]; W[1]=g.mu[c]/W[0]; W[2]=g.mv[c]/W[0]; W[3]=g.mw[c]/W[0]; W[4]=eint(g,c); };
    bool WB = !REF_P0.empty();   // route 1: well-balanced hydrostatic reconstruction (z-sweep only)
    // ---- hydro: unsplit MUSCL+HLLC (P only; z-dir reconstructs P-deviation when WB) ----
    auto sweep=[&](int dir,bool wb){
        int nL=(dir==0?g.nx:(dir==1?g.ny:g.nz)),nA=(dir==0?g.ny:g.nx),nB=(dir==2?g.ny:g.nz);
        auto CELL=[&](int p,int a,int b){ return dir==0?g.idx(p,a,b):(dir==1?g.idx(a,p,b):g.idx(a,b,p)); };
        int in=(dir==0?1:(dir==1?2:3)),it1=(dir==0?2:1),it2=(dir==2?2:3);
        vector<double> Fr(nL+1),Fmn(nL+1),Ft1(nL+1),Ft2(nL+1),FE(nL+1),Frc(nL+1);
        for(int a=0;a<nA;a++)for(int b=0;b<nB;b++){
            for(int f=0;f<=nL;f++){ int iL=max(0,f-1),iR=min(nL-1,f),iLL=max(0,f-2),iRR=min(nL-1,f+1);
                double A[5],B[5],AA[5],BB[5],L[5],R[5]; prim(CELL(iL,a,b),A);prim(CELL(iR,a,b),B);prim(CELL(iLL,a,b),AA);prim(CELL(iRR,a,b),BB);
                for(int q=0;q<5;q++){L[q]=A[q]+0.5*mm(A[q]-AA[q],B[q]-A[q]);R[q]=B[q]-0.5*mm(B[q]-A[q],BB[q]-B[q]);}
                if(L[0]<1e-9)L[0]=A[0];if(R[0]<1e-9)R[0]=B[0];if(L[4]<0)L[4]=A[4];if(R[4]<0)R[4]=B[4];
                F5 fl;
                if(wb){ // reconstruct the pressure DEVIATION from the reference; inject P linearly (robust); reference face P from the EOS of the averaged reference density (=> ~0 at the free surface)
                    int cLL=CELL(iLL,a,b),cL=CELL(iL,a,b),cR=CELL(iR,a,b),cRR=CELL(iRR,a,b);
                    double dLL=Pcell(g,cLL)-REF_P0[cLL],dL=Pcell(g,cL)-REF_P0[cL],dR=Pcell(g,cR)-REF_P0[cR],dRR=Pcell(g,cRR)-REF_P0[cRR];
                    double dLf=dL+0.5*mm(dL-dLL,dR-dL), dRf=dR-0.5*mm(dR-dL,dRR-dR);
                    double p0f=MAT.pressure(0.5*(REF_R0[cL]+REF_R0[cR]),0.0); if(p0f<0)p0f=0;
                    fl=hllc_p(L[0],L[in],L[it1],L[it2],L[4],p0f+dLf, R[0],R[in],R[it1],R[it2],R[4],p0f+dRf);
                } else fl=hllc(L[0],L[in],L[it1],L[it2],L[4],R[0],R[in],R[it1],R[it2],R[4]);
                Fr[f]=fl.r;Fmn[f]=fl.mn;Ft1[f]=fl.mt1;Ft2[f]=fl.mt2;FE[f]=fl.E;
                int cl=CELL(iL,a,b),cr=CELL(iR,a,b);   // passive material tracer: species flux = mass flux * upwind c
                Frc[f]=fl.r*(fl.r>=0.0 ? g.rc[cl]/g.r[cl] : g.rc[cr]/g.r[cr]); }
            bool axi=(AXISYM && dir==0);   // cylindrical r-sweep: area-weight fluxes by face radius (axis face r=0 -> zero area = automatic axis BC)
            for(int p=0;p<nL;p++){ int c=CELL(p,a,b);
                double w0,w1,inv;
                if(axi){ w0=p*g.dx; w1=(p+1)*g.dx; inv=1.0/(((p+0.5)*g.dx)*g.dx); } else { w0=w1=1.0; inv=invdx; }
                d.r[c]+=-(w1*Fr[p+1]-w0*Fr[p])*inv; d.E[c]+=-(w1*FE[p+1]-w0*FE[p])*inv; d.rc[c]+=-(w1*Frc[p+1]-w0*Frc[p])*inv;
                double dmn=-(w1*Fmn[p+1]-w0*Fmn[p])*inv,dt1=-(w1*Ft1[p+1]-w0*Ft1[p])*inv,dt2=-(w1*Ft2[p+1]-w0*Ft2[p])*inv;
                if(axi) dmn+=Pcell(g,c)/((p+0.5)*g.dx);   // cylindrical pressure geometric source for the radial momentum (corrects the area-weighted p flux)
                if(dir==0){d.mu[c]+=dmn;d.mv[c]+=dt1;d.mw[c]+=dt2;}else if(dir==1){d.mv[c]+=dmn;d.mu[c]+=dt1;d.mw[c]+=dt2;}else{d.mw[c]+=dmn;d.mu[c]+=dt1;d.mv[c]+=dt2;}
            }
        }
    };
    sweep(0,false);sweep(1,false);sweep(2,WB);
    if(GZ!=0.0){
        if(!WB){ for(int c=0;c<n;c++){ d.mw[c]+=-GZ*g.r[c]; d.E[c]+=-GZ*g.mw[c]; } }   // plain gravity body force
        else { for(int i=0;i<g.nx;i++)for(int j=0;j<g.ny;j++)for(int k=0;k<g.nz;k++){ int c=g.idx(i,j,k);
            double r0t=0.5*(REF_R0[c]+REF_R0[g.idx(i,j,min(g.nz-1,k+1))]), r0b=0.5*(REF_R0[g.idx(i,j,max(0,k-1))]+REF_R0[c]);
            double p0t=MAT.pressure(r0t,0.0); if(p0t<0)p0t=0; double p0b=MAT.pressure(r0b,0.0); if(p0b<0)p0b=0;  // EOS reference face P (matches the flux)
            d.mw[c]+=-GZ*(g.r[c]-REF_R0[c])+(p0t-p0b)*invdx;   // WB source: cancels the reference flux divergence to machine precision
            d.E[c]+=-GZ*g.mw[c]; } }
    }
    if(G<=0) return;                      // no strength (ideal gas) -> done
    // ---- strength sources (per cell, central diffs; clamp at boundaries) ----
    auto vx=[&](int c){return g.mu[c]/g.r[c];}; auto vy=[&](int c){return g.mv[c]/g.r[c];}; auto vz=[&](int c){return g.mw[c]/g.r[c];};
    for(int i=0;i<g.nx;i++)for(int j=0;j<g.ny;j++)for(int k=0;k<g.nz;k++){ int c=g.idx(i,j,k);
        if(RHO_CFL>0&&g.r[c]<RHO_CFL) continue;   // no strength in vacuum cells
        int xm=g.idx(max(0,i-1),j,k),xp=g.idx(min(g.nx-1,i+1),j,k),ym=g.idx(i,max(0,j-1),k),yp=g.idx(i,min(g.ny-1,j+1),k),zm=g.idx(i,j,max(0,k-1)),zp=g.idx(i,j,min(g.nz-1,k+1));
        double hx=( (i>0&&i<g.nx-1)?2.0:1.0)*g.dx, hy=((j>0&&j<g.ny-1)?2.0:1.0)*g.dx, hz=((k>0&&k<g.nz-1)?2.0:1.0)*g.dx;
        // void neighbours -> centre velocity (zero gradient toward vacuum = traction-free free surface)
        auto vxw=[&](int cc){return (RHO_CFL>0&&g.r[cc]<RHO_CFL)?g.mu[c]/g.r[c]:g.mu[cc]/g.r[cc];};
        auto vyw=[&](int cc){return (RHO_CFL>0&&g.r[cc]<RHO_CFL)?g.mv[c]/g.r[c]:g.mv[cc]/g.r[cc];};
        auto vzw=[&](int cc){return (RHO_CFL>0&&g.r[cc]<RHO_CFL)?g.mw[c]/g.r[c]:g.mw[cc]/g.r[cc];};
        // velocity gradient L[a][b]=d v_a / d x_b
        double Lxx=(vxw(xp)-vxw(xm))/hx, Lxy=(vxw(yp)-vxw(ym))/hy, Lxz=(vxw(zp)-vxw(zm))/hz;
        double Lyx=(vyw(xp)-vyw(xm))/hx, Lyy=(vyw(yp)-vyw(ym))/hy, Lyz=(vyw(zp)-vyw(zm))/hz;
        double Lzx=(vzw(xp)-vzw(xm))/hx, Lzy=(vzw(yp)-vzw(ym))/hy, Lzz=(vzw(zp)-vzw(zm))/hz;
        double rcyl=(i+0.5)*g.dx;   // cylindrical radius of cell centre (AXISYM geometric terms; r_c>0 for cell centres)
        // hoop (theta-theta) strain rate is GEOMETRIC in cylindrical coords: e_thth = u/r = v_x/r,
        // NOT dv_y/dy (=Lyy=0 for ny=1). This also feeds the trace, hence the deviatoric strain of ALL components.
        double eyy=AXISYM?(vx(c)/rcyl):Lyy;
        double exx=Lxx,ezz=Lzz,exy=0.5*(Lxy+Lyx),exz=0.5*(Lxz+Lzx),eyz=0.5*(Lyz+Lzy),tr=(exx+eyy+ezz)/3.0;
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
        d.vib[c]+=-adv(g.vib);   // AF vibrational velocity advects with the flow (intensive material state, like S; seed/decay are operator-split in update_af). Scalar -> no cylindrical geometric term.
        // deviatoric stress -> momentum (div S) + energy (div(S.v)), central diffs
        double dSxx_x=(g.Sxx[xp]-g.Sxx[xm])/hx, dSxy_y=(g.Sxy[yp]-g.Sxy[ym])/hy, dSxz_z=(g.Sxz[zp]-g.Sxz[zm])/hz;
        double dSxy_x=(g.Sxy[xp]-g.Sxy[xm])/hx, dSyy_y=(g.Syy[yp]-g.Syy[ym])/hy, dSyz_z=(g.Syz[zp]-g.Syz[zm])/hz;
        double dSxz_x=(g.Sxz[xp]-g.Sxz[xm])/hx, dSyz_y=(g.Syz[yp]-g.Syz[ym])/hy, dSzz_z=(g.Szz[zp]-g.Szz[zm])/hz;
        d.mu[c]+=dSxx_x+dSxy_y+dSxz_z; d.mv[c]+=dSxy_x+dSyy_y+dSyz_z; d.mw[c]+=dSxz_x+dSyz_y+dSzz_z;
        if(AXISYM){   // geometric stress-divergence sources (curvilinear basis): (div S)_r += (S_rr-S_thth)/r, (div S)_z += S_rz/r
            d.mu[c]+=(g.Sxx[c]-g.Syy[c])/rcyl;   // radial: the -S_thth/r curvature term plus the +S_rr/r from (1/r)d(r S_rr)/dr
            d.mw[c]+=g.Sxz[c]/rcyl; }            // axial: the +S_rz/r from (1/r)d(r S_rz)/dr
        if(g.af[c]>0&&ETA_AF>0&&!VISC_IMP){ double eta=g.af[c]*ETA_AF, id2=1.0/(g.dx*g.dx);   // AF Newtonian viscosity (damps fluidized flow); Cartesian vector Laplacian. VISC_IMP -> handled implicitly in implicit_visc (no double counting)
            d.mu[c]+=eta*(vx(xp)+vx(xm)+vx(yp)+vx(ym)+vx(zp)+vx(zm)-6*vx(c))*id2;
            d.mv[c]+=eta*(vy(xp)+vy(xm)+vy(yp)+vy(ym)+vy(zp)+vy(zm)-6*vy(c))*id2;
            d.mw[c]+=eta*(vz(xp)+vz(xm)+vz(yp)+vz(ym)+vz(zp)+vz(zm)-6*vz(c))*id2;
            if(AXISYM){   // cylindrical vector-Laplacian corrections (add to the Cartesian d^2/dr^2+d^2/dz^2 above):
                d.mu[c]+=eta*(Lxx/rcyl - vx(c)/(rcyl*rcyl));   // (lap v)_r = lap_scalar(u) - u/r^2, lap_scalar adds (1/r)du/dr
                d.mw[c]+=eta*(Lzx/rcyl); } }                   // (lap v)_z = lap_scalar(w), adds (1/r)dw/dr (no -w/r^2: z is Cartesian-like)
        // energy: div(S.v). (S.v)_x = Sxx u+Sxy v+Sxz w, etc.
        auto Svx=[&](int cc){return g.Sxx[cc]*vx(cc)+g.Sxy[cc]*vy(cc)+g.Sxz[cc]*vz(cc);};
        auto Svy=[&](int cc){return g.Sxy[cc]*vx(cc)+g.Syy[cc]*vy(cc)+g.Syz[cc]*vz(cc);};
        auto Svz=[&](int cc){return g.Sxz[cc]*vx(cc)+g.Syz[cc]*vy(cc)+g.Szz[cc]*vz(cc);};
        d.E[c]+=(Svx(xp)-Svx(xm))/hx+(Svy(yp)-Svy(ym))/hy+(Svz(zp)-Svz(zm))/hz;
        if(AXISYM) d.E[c]+=Svx(c)/rcyl;   // cylindrical geometric term of div(S.v): + (S.v)_r/r
    }
}
static void vonmises(Grid&g){ double Y0=MAT.Y; if(Y0<=0&&!ROCK)return; int n=g.nx*g.ny*g.nz;
    for(int c=0;c<n;c++){ double Y;
        if(ROCK){ double P=max(Pcell(g,c),0.0);   // pressure-dependent yield (friction only under compression)
            double Yi=Y_I0 + MU_I*P/(1.0 + MU_I*P/max(Y_M-Y_I0,1.0));   // intact (Lundborg): cohesion + saturating friction
            double Yd=min(Y_D0 + MU_D*P, Y_M);                           // damaged: residual cohesion + friction, capped at the limit
            double Yb=(1.0-g.D[c])*Yi + g.D[c]*Yd;                       // unfluidized (base) strength
            Y=Yb*(1.0-g.af[c]);                                          // damage blends intact->friction; AF fluidizes
            // Bingham floor (AF surgery option B): a fluidized cell retains a finite flow
            // resistance Y_BINGHAM -- never exceeding its unfluidized strength -- so that
            // long-TDEC fluidization slumps instead of digging without limit. 0 = legacy.
            if(Y_BINGHAM>0 && g.af[c]>1e-3){ double Yf=min(Y_BINGHAM,Yb); if(Y<Yf) Y=Yf; } }
        else Y=(1.0-g.D[c])*(1.0-g.af[c])*Y0;                            // cohesion-only von Mises (back-compat)
        double sxx=g.Sxx[c],syy=g.Syy[c],szz=g.Szz[c],sxy=g.Sxy[c],sxz=g.Sxz[c],syz=g.Syz[c];
        double J2=0.5*(sxx*sxx+syy*syy+szz*szz)+sxy*sxy+sxz*sxz+syz*syz; double vm=sqrt(3.0*J2);
        if(vm>Y){double f=(vm>0?Y/vm:0.0); g.Sxx[c]*=f;g.Syy[c]*=f;g.Szz[c]*=f;g.Sxy[c]*=f;g.Sxz[c]*=f;g.Syz[c]*=f;} } }
static void grow_damage(Grid&g,double dt){ double Em=MAT.Emod,wk=MAT.wk,wm=MAT.wm; if(Em<=0||wk<=0)return;
    double dx=g.dx,eps_act=pow(1.0/(wk*dx*dx*dx),1.0/wm); int n=g.nx*g.ny*g.nz;   // weakest-flaw activation strain
    for(int c=0;c<n;c++){ if(g.D[c]>=1.0)continue; double P=Pcell(g,c);
        double sigmax=max(-P+g.Sxx[c],max(-P+g.Syy[c],-P+g.Szz[c]));   // max tensile principal stress (diagonal proxy)
        if(sigmax>0.0 && sigmax/Em>eps_act){ double cg=0.4*sqrt(Em/max(g.r[c],1e-30)),Rs=0.5*dx;
            double d13=pow(g.D[c],1.0/3.0)+(cg/Rs)*dt; g.D[c]=min(1.0,d13*d13*d13); } } }
static void set_ref(const Grid&g){ int n=g.nx*g.ny*g.nz; REF_R0.assign(n,0); REF_P0.assign(n,0);   // freeze IC as the hydrostatic reference (route 1)
    for(int c=0;c<n;c++){ REF_R0[c]=g.r[c]; REF_P0[c]=Pcell(g,c); } }
// route 1: void cells = passive vacuum (reset to reference). Also kills the near-vacuum velocity (mu/rho_tiny) that
// would otherwise feed spurious deviatoric stress into the strength solver at the free surface.
static double VOID_AMB = 0.0;  // >0: void cells reset to this AMBIENT density instead of REF_R0 (a lithostatic reference REFILLS evacuated below-surface cells with solid rock -- the 2026-07-02 prater finding); 0 = legacy
static long VZ_N=0; static double VZ_DM=0;  // instrumentation: count of SOLID-reference void resets (REF_R0>half-density = below-surface refill under legacy) + the rho they inject(ed)
static long SEED_N=0;          // diagnostic (2026-07-06): cumulative vib seeding events -- discriminates cavitation-reseed (ongoing seeds) from af-ratio persistence (seeds stop after excavation)
static double SNAP_DT=0.0; static const char*SNAP_DIR=0;   // env CRATER_SNAP_DT/CRATER_SNAP_DIR: periodic full-field snapshots (diagnostic only)
static void void_cells(Grid&g){ if(RHO_CFL<=0)return; bool wb=!REF_R0.empty(); int N=g.r.size();
    for(int c=0;c<N;c++) if(g.r[c]<RHO_CFL){ g.mu[c]=g.mv[c]=g.mw[c]=0; if(wb){
        if(REF_R0[c]>1350.0){VZ_N++;VZ_DM+=REF_R0[c]-g.r[c];}   // a below-surface (solid-reference) void: legacy refills it with rock
        g.r[c]=(VOID_AMB>0?VOID_AMB:REF_R0[c]);g.E[c]=0;g.Sxx[c]=g.Syy[c]=g.Szz[c]=g.Sxy[c]=g.Sxz[c]=g.Syz[c]=0;g.rc[c]=0;} } }
// ---- Option C (2026-07-07): implicit AF viscous update (operator-split, backward Euler) ----
// Fluidized cells obey d(rho v)/dt = div(eta grad v), eta = af*ETA_AF, solved IMPLICITLY so
// eta is NOT capped by the explicit viscous dt limit (~2.4e8 Pa s at crater scales, found
// 2026-06-28); the block-model physics needs 1e9-1e10 Pa s (Wuennemann & Ivanov 2003).
// Finite-volume, variable coefficient, face-HARMONIC eta: a zero-eta side (void, or
// unfluidized af<=1e-3) gives zero face flux, so those cells reduce to the identity and stay
// in the system -- no masking special cases. AXISYM: fluxes are face-area weighted (radial
// face area ~ r_face; the r=0 axis face has zero area = the regularity BC) and the r/hoop
// components carry the -v/r^2 curvature term as a POSITIVE diagonal, so v_r = A*r is an
// exact discrete null mode (gate visc_diff B). The volume-weighted system is SPD;
// Jacobi-preconditioned CG to 1e-8 relative. ENERGY: the KE removed by the smoothing is
// returned as heat, deposited per cell proportional to the local face dissipation
// eta*(dv)^2 (positive definite) and rescaled so total energy is conserved EXACTLY
// (the explicit Lop term never heated; this path does -- gate visc_energy).
static long VISC_IT=0;   // instrumentation: cumulative CG iterations
static void implicit_visc(Grid&g,double dt){
    if(!VISC_IMP||ETA_AF<=0) return;
    const int nx=g.nx,ny=g.ny,nz=g.nz,n=nx*ny*nz; const double dx=g.dx;
    static vector<double> eta,Vc,cfx,cfy,cfz,dsum,dax,bb,rr,zz,pp,qq,hr,u1,v1,w1;
    eta.assign(n,0.0); bool any=false; double rmin=max(RHO_CFL,1e-12);
    for(int c=0;c<n;c++) if(g.r[c]>rmin&&g.af[c]>1e-3){ eta[c]=g.af[c]*ETA_AF; any=true; }
    if(!any) return;
    Vc.assign(n,0);cfx.assign(n,0);cfy.assign(n,0);cfz.assign(n,0);dsum.assign(n,0);dax.assign(n,0);
    auto harm=[](double a,double b){ return (a>0&&b>0)?2.0*a*b/(a+b):0.0; };
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){ int c=g.idx(i,j,k);
        double rc=(i+0.5)*dx, wV=(AXISYM?rc:dx); Vc[c]=wV*dx*dx; dsum[c]=g.r[c]*Vc[c];
        if(i<nx-1)      cfx[c]=dt*((AXISYM?(i+1)*dx:dx))*harm(eta[c],eta[g.idx(i+1,j,k)]);   // dt*A_f*eta_f/dx, radial face at r=(i+1)dx
        if(ny>1&&j<ny-1)cfy[c]=dt*dx*harm(eta[c],eta[g.idx(i,j+1,k)]);
        if(k<nz-1)      cfz[c]=dt*wV*harm(eta[c],eta[g.idx(i,j,k+1)]);
        if(AXISYM)      dax[c]=dt*eta[c]*Vc[c]/(rc*rc);   // -v/r^2 curvature term (r and hoop components)
    }
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){ int c=g.idx(i,j,k);
        if(i<nx-1){ dsum[c]+=cfx[c]; dsum[g.idx(i+1,j,k)]+=cfx[c]; }
        if(ny>1&&j<ny-1){ dsum[c]+=cfy[c]; dsum[g.idx(i,j+1,k)]+=cfy[c]; }
        if(k<nz-1){ dsum[c]+=cfz[c]; dsum[g.idx(i,j,k+1)]+=cfz[c]; }
    }
    auto Amul=[&](const vector<double>&x,vector<double>&y,bool wdax){
        for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){ int c=g.idx(i,j,k);
            double s=(dsum[c]+(wdax?dax[c]:0.0))*x[c];
            if(i<nx-1) s-=cfx[c]*x[g.idx(i+1,j,k)];
            if(i>0){ int cl=g.idx(i-1,j,k); s-=cfx[cl]*x[cl]; }
            if(ny>1){ if(j<ny-1) s-=cfy[c]*x[g.idx(i,j+1,k)];
                      if(j>0){ int cl=g.idx(i,j-1,k); s-=cfy[cl]*x[cl]; } }
            if(k<nz-1) s-=cfz[c]*x[g.idx(i,j,k+1)];
            if(k>0){ int cl=g.idx(i,j,k-1); s-=cfz[cl]*x[cl]; }
            y[c]=s; } };
    auto cg=[&](vector<double>&x,const vector<double>&b,bool wdax){
        rr.assign(n,0);zz.assign(n,0);pp.assign(n,0);qq.assign(n,0);
        double bn=0; for(int c=0;c<n;c++) bn+=b[c]*b[c];
        if(bn<=0){ for(int c=0;c<n;c++) x[c]=0; return; }
        Amul(x,qq,wdax);
        double rn=0; for(int c=0;c<n;c++){ rr[c]=b[c]-qq[c]; rn+=rr[c]*rr[c]; }
        double tol2=1e-16*bn, rz=0;
        for(int c=0;c<n;c++){ zz[c]=rr[c]/(dsum[c]+(wdax?dax[c]:0.0)); pp[c]=zz[c]; rz+=rr[c]*zz[c]; }
        for(int it=0; it<800&&rn>tol2; it++){
            Amul(pp,qq,wdax);
            double pq=0; for(int c=0;c<n;c++) pq+=pp[c]*qq[c];
            if(pq<=0) break;
            double al=rz/pq; rn=0;
            for(int c=0;c<n;c++){ x[c]+=al*pp[c]; rr[c]-=al*qq[c]; rn+=rr[c]*rr[c]; }
            VISC_IT++;
            if(rn<=tol2) break;
            double rz2=0; for(int c=0;c<n;c++){ zz[c]=rr[c]/(dsum[c]+(wdax?dax[c]:0.0)); rz2+=rr[c]*zz[c]; }
            double be=rz2/rz; rz=rz2;
            for(int c=0;c<n;c++) pp[c]=zz[c]+be*pp[c];
        } };
    u1.assign(n,0);v1.assign(n,0);w1.assign(n,0);bb.assign(n,0);
    for(int c=0;c<n;c++){ u1[c]=g.mu[c]/g.r[c]; bb[c]=Vc[c]*g.mu[c]; } cg(u1,bb,AXISYM);
    for(int c=0;c<n;c++){ v1[c]=g.mv[c]/g.r[c]; bb[c]=Vc[c]*g.mv[c]; } cg(v1,bb,AXISYM);   // hoop component carries the same -v/r^2
    for(int c=0;c<n;c++){ w1[c]=g.mw[c]/g.r[c]; bb[c]=Vc[c]*g.mw[c]; } cg(w1,bb,false);
    hr.assign(n,0); double totold=0,totnew=0,htot=0;
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){ int c=g.idx(i,j,k);
        double keo=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c];
        double ken=0.5*g.r[c]*(u1[c]*u1[c]+v1[c]*v1[c]+w1[c]*w1[c]);
        totold+=keo*Vc[c]; totnew+=ken*Vc[c];
        auto fh=[&](int cn,double cf){ double du=u1[c]-u1[cn],dv=v1[c]-v1[cn],dw=w1[c]-w1[cn];
            double h=cf*(du*du+dv*dv+dw*dw); hr[c]+=0.5*h; hr[cn]+=0.5*h; htot+=h; };
        if(i<nx-1&&cfx[c]>0) fh(g.idx(i+1,j,k),cfx[c]);
        if(ny>1&&j<ny-1&&cfy[c]>0) fh(g.idx(i,j+1,k),cfy[c]);
        if(k<nz-1&&cfz[c]>0) fh(g.idx(i,j,k+1),cfz[c]);
        if(AXISYM&&dax[c]>0){ double h=dax[c]*(u1[c]*u1[c]+v1[c]*v1[c]); hr[c]+=h; htot+=h; }
    }
    double loss=totold-totnew, sc=(loss>0&&htot>0)?loss/htot:0.0;
    for(int c=0;c<n;c++){
        double keo=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c];
        g.mu[c]=g.r[c]*u1[c]; g.mv[c]=g.r[c]*v1[c]; g.mw[c]=g.r[c]*w1[c];
        double ken=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c];
        g.E[c]+=(ken-keo)+sc*hr[c]/Vc[c];   // internal energy rises only by the (rescaled, positive) local dissipation
    }
}
static double maxspeed(const Grid&g){ double s=1e-30; int n=g.r.size(); double G=MAT.G;
    for(int c=0;c<n;c++){ if(g.r[c]<RHO_CFL) continue; double p,cs; eos_pc(g.r[c],eint(g,c),p,cs); double cel=sqrt(cs*cs+ (G>0?(4.0/3.0)*G/g.r[c]:0.0));
        double v=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; s=max(s,v+cel);} return s; }
// Shock-activated acoustic fluidization (Wuennemann-Ivanov block model, iSALE-practice form).
// The passing shock seeds a vibrational velocity vib ~ post-shock particle speed; vib decays with
// TDEC; the rock is fluidized where the acoustic vibrational pressure p_vib=rho*c_s*vib exceeds the
// overburden p_ov=max(P,P_COH). The resulting af in [0,1) degrades shear strength (vonmises) and
// drives viscous flow (Lop). (Phase 1: per-cell activation/decay; advection of vib is Phase 2.)
static void update_af(Grid&g,double dt){
    if(C_ACT<=0||TDEC<=0) return;            // shock-activated AF off -> af stays at its IC value
    int n=g.nx*g.ny*g.nz; double dec=exp(-dt/TDEC);
    for(int c=0;c<n;c++){
        double P=Pcell(g,c);
        if(P>g.Pmax[c]){                     // new pressure peak: track it; seed vib only if the rise is SHOCK-sized (> P_ACT)
            double dP=P-g.Pmax[c]; g.Pmax[c]=P;
            if(dP>P_ACT){ double sp=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/max(g.r[c],1e-30);
                g.vib[c]=max(g.vib[c],C_ACT*sp); SEED_N++; } }
        g.vib[c]*=dec;                        // acoustic vibrations decay
        double cs=MAT.sound_speed(g.r[c],eint(g,c)), pvib=g.r[c]*cs*g.vib[c], pov=max(P,P_COH);
        g.af[c]=pvib/(pvib+pov);              // fluidization ratio in [0,1)
    }
}
static void step_rk2(Grid&g,double dt){ int n=g.r.size(); DU d(n);
    update_af(g,dt);   // shock-activated AF: refresh vib + derive af BEFORE strength/viscosity read it
    Grid g1=g;
    Lop(g,d);
    for(int c=0;c<n;c++){g1.r[c]=g.r[c]+dt*d.r[c];g1.mu[c]=g.mu[c]+dt*d.mu[c];g1.mv[c]=g.mv[c]+dt*d.mv[c];g1.mw[c]=g.mw[c]+dt*d.mw[c];g1.E[c]=g.E[c]+dt*d.E[c];
        g1.Sxx[c]=g.Sxx[c]+dt*d.Sxx[c];g1.Syy[c]=g.Syy[c]+dt*d.Syy[c];g1.Szz[c]=g.Szz[c]+dt*d.Szz[c];g1.Sxy[c]=g.Sxy[c]+dt*d.Sxy[c];g1.Sxz[c]=g.Sxz[c]+dt*d.Sxz[c];g1.Syz[c]=g.Syz[c]+dt*d.Syz[c];g1.D[c]=g.D[c]+dt*d.D[c];g1.rc[c]=g.rc[c]+dt*d.rc[c];g1.vib[c]=g.vib[c]+dt*d.vib[c];}
    void_cells(g1); vonmises(g1); Lop(g1,d);   // clean the predictor's void cells so strength sees no near-vacuum velocity
    for(int c=0;c<n;c++){g.r[c]=0.5*(g.r[c]+g1.r[c]+dt*d.r[c]);g.mu[c]=0.5*(g.mu[c]+g1.mu[c]+dt*d.mu[c]);g.mv[c]=0.5*(g.mv[c]+g1.mv[c]+dt*d.mv[c]);g.mw[c]=0.5*(g.mw[c]+g1.mw[c]+dt*d.mw[c]);g.E[c]=0.5*(g.E[c]+g1.E[c]+dt*d.E[c]);
        g.Sxx[c]=0.5*(g.Sxx[c]+g1.Sxx[c]+dt*d.Sxx[c]);g.Syy[c]=0.5*(g.Syy[c]+g1.Syy[c]+dt*d.Syy[c]);g.Szz[c]=0.5*(g.Szz[c]+g1.Szz[c]+dt*d.Szz[c]);g.Sxy[c]=0.5*(g.Sxy[c]+g1.Sxy[c]+dt*d.Sxy[c]);g.Sxz[c]=0.5*(g.Sxz[c]+g1.Sxz[c]+dt*d.Sxz[c]);g.Syz[c]=0.5*(g.Syz[c]+g1.Syz[c]+dt*d.Syz[c]);g.D[c]=0.5*(g.D[c]+g1.D[c]+dt*d.D[c]);g.rc[c]=0.5*(g.rc[c]+g1.rc[c]+dt*d.rc[c]);g.vib[c]=0.5*(g.vib[c]+g1.vib[c]+dt*d.vib[c]);}
    implicit_visc(g,dt);   // Option C: operator-split implicit AF viscosity (no-op unless VISC_IMP)
    grow_damage(g,dt); vonmises(g);
    // (AF vibration update moved to update_af() at the start of the step; af is now derived from vib, not decayed here)
    if(DAMP<1.0){ int N=g.r.size(); for(int c=0;c<N;c++){ g.mu[c]*=DAMP; g.mv[c]*=DAMP; g.mw[c]*=DAMP; } }   // relaxation damping
    void_cells(g);   // void cells = passive vacuum (reset to reference)
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
    } else if(mode=="tracer"){ // M-tag: a c=1 material blob advects with a uniform flow -- must conserve, move at v0, stay in [0,1]
        MAT=Material::ideal(1.4); GAM=1.4; int N=200;double dx=1.0/N,CFL=0.4; double rho0=1.0,p0=1.0,v0=2.0,tend=0.15; Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;int c=g.idx(i,0,0);g.r[c]=rho0;g.mu[c]=rho0*v0;g.E[c]=p0/(GAM-1)+0.5*rho0*v0*v0;
            g.rc[c]=rho0*((x>0.3&&x<0.5)?1.0:0.0);}
        double rc0=0,cen0=0;for(int i=0;i<N;i++){int c=g.idx(i,0,0);rc0+=g.rc[c];cen0+=g.rc[c]*(i+0.5)*dx;} cen0/=rc0;
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double rc1=0,cen1=0,cmin=1e30,cmax=-1e30;
        for(int i=0;i<N;i++){int c=g.idx(i,0,0);double cc=g.rc[c]/g.r[c];rc1+=g.rc[c];cen1+=g.rc[c]*(i+0.5)*dx;cmin=min(cmin,cc);cmax=max(cmax,cc);}
        cen1/=rc1; double vcen=(cen1-cen0)/t, merr=fabs(rc1-rc0)/rc0, verr=fabs(vcen-v0)/v0;
        printf("Tracer: sum(rc) %.6e->%.6e (err %.2e); centroid v=%.4f vs v0=%.1f (err %.2e); c in [%.2e, %.4f]\n",rc0,rc1,merr,vcen,v0,verr,cmin,cmax);
        printf("GATE (sum(rc) conserved <1e-6 & centroid v within 1%% & c in [0,1]): %s\n",(merr<1e-6&&verr<0.01&&cmin>=-1e-9&&cmax<=1.0+1e-9)?"PASS":"CHECK");
    } else if(mode=="sedov"){ // 3D Sedov-Taylor point blast: shock radius vs analytic self-similar solution
        MAT=Material::ideal(1.4); GAM=1.4; int N=64;double L=2.0,dx=L/N,CFL=0.3,tend=0.5; double E0=1.0,rho0=1.0; Grid g(N,N,N,dx);
        for(int c=0;c<N*N*N;c++){g.r[c]=rho0;g.E[c]=1e-4/(GAM-1);}
        int ic=N/2; g.E[g.idx(ic,ic,ic)]+=E0/(dx*dx*dx);   // point energy in the central cell
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double rpk=0,rhomax=0,cen=(ic+0.5)*dx;   // shock radius = radius of peak density; compression = max density
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++){int c=g.idx(i,j,k);double rr=g.r[c];
            if(rr>rhomax){rhomax=rr;double xx=(i+0.5)*dx-cen,yy=(j+0.5)*dx-cen,zz=(k+0.5)*dx-cen;rpk=sqrt(xx*xx+yy*yy+zz*zz);}}
        double alpha=0.851;   // Sedov dimensionless energy constant, 3D gamma=1.4
        double Ran=pow(E0/(alpha*rho0),0.2)*pow(tend,0.4);
        printf("Sedov 3D: shock R_num=%.3f analytic=%.3f (err %.1f%%), compression=%.2f (strong-shock %.1f, resolution-limited)\n",rpk,Ran,100*fabs(rpk-Ran)/Ran,rhomax/rho0,(GAM+1)/(GAM-1));
        printf("GATE (shock radius within 10%% of analytic Sedov, 64^3 resolution-limited): %s\n",(fabs(rpk-Ran)/Ran<0.10)?"PASS":"CHECK");
    } else if(mode=="sedov_axi"){ // axisymmetric (r,z) point blast ON THE AXIS = a 3D spherical blast; validates the cylindrical geometry vs the same analytic Sedov solution as the 3D gate
        MAT=Material::ideal(1.4); GAM=1.4; AXISYM=true; double PI=3.141592653589793;
        int Nr=64,Nz=128; double dx=2.0/64,CFL=0.3,tend=0.5, E0=1.0,rho0=1.0; Grid g(Nr,1,Nz,dx);
        for(int c=0;c<Nr*Nz;c++){g.r[c]=rho0;g.E[c]=1e-4/(GAM-1);}
        int kc=Nz/2; g.E[g.idx(0,0,kc)]+=E0/(PI*dx*dx*dx);   // point energy in the on-axis cell (its 3D volume = 2*pi*(dx/2)*dx*dx = pi*dx^3)
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double rpk=0,rhomax=0,zc=(kc+0.5)*dx;
        for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){int c=g.idx(i,0,k);double rr=g.r[c];
            if(rr>rhomax){rhomax=rr;double rad=(i+0.5)*dx,zz=(k+0.5)*dx-zc;rpk=sqrt(rad*rad+zz*zz);}}
        double alpha=0.851, Ran=pow(E0/(alpha*rho0),0.2)*pow(tend,0.4);
        printf("Sedov axisym (on-axis point blast = 3D spherical): shock R_num=%.3f analytic=%.3f (err %.1f%%), compression=%.2f\n",rpk,Ran,100*fabs(rpk-Ran)/Ran,rhomax/rho0);
        printf("GATE (cylindrical geometry: shock radius within 10%% of spherical Sedov): %s\n",(fabs(rpk-Ran)/Ran<0.10)?"PASS":"CHECK");
    } else if(mode=="lame"){ // Phase-2b: validate cylindrical elastic STRENGTH -- hoop strain rate + geometric stress divergence -- vs the analytic Lame thick-cylinder solution
        MAT=Material::basalt(); AXISYM=true; double G=MAT.G, dx=500.0; bool A=false,B=false;
        // (A) HOOP STRAIN RATE: a uniform radial expansion u(r)=edot*r must drive dS via e_thth=u/r (geometric), not dv_y/dy(=0).
        //     From zero stress: e_rr=e_thth=edot, e_zz=0, tr=2edot/3 -> dS_rr=dS_thth=2G(edot-tr), dS_zz=2G(-tr); Jaumann+advection vanish (S=0).
        {   int Nr=40; Grid g(Nr,1,1,dx); double edot=1e-6;   // tiny rate -> stress stays << yield over one Lop
            for(int i=0;i<Nr;i++){int c=g.idx(i,0,0); g.r[c]=MAT.rho0; g.E[c]=0; g.mu[c]=MAT.rho0*edot*((i+0.5)*dx);}
            DU d(Nr); Lop(g,d);
            int c=g.idx(Nr/2,0,0); double trc=2.0*edot/3.0;
            double an_rr=2*G*(edot-trc), an_th=2*G*(edot-trc), an_zz=2*G*(0.0-trc);
            double err=fabs(d.Sxx[c]-an_rr)+fabs(d.Syy[c]-an_th)+fabs(d.Szz[c]-an_zz), scale=fabs(an_rr)+fabs(an_zz)+1e-30;
            A=(err/scale<1e-3);
            printf("Lame (A hoop): dS_rr=%.4e/%.4e  dS_thth=%.4e/%.4e  dS_zz=%.4e/%.4e (num/analytic, rel err %.2e)\n",
                d.Sxx[c],an_rr,d.Syy[c],an_th,d.Szz[c],an_zz,err/scale); }
        // (B) LAME STATIC EQUILIBRIUM: prescribe the deviatoric stress of a thick cylinder under internal pressure,
        //     S_rr=-Bc/r^2, S_thth=+Bc/r^2, S_zz=0 (traceless), at rest with uniform P=0. The continuum radial balance
        //     dS_rr/dr + (S_rr-S_thth)/r = 0 is exact, so the discrete radial-momentum residual d.mu must be truncation-small.
        {   int Nr=200; int ia=20,ib=180; double ia_r=(ia+0.5)*dx, Bc=2.0e6*ia_r*ia_r;   // ~2 MPa peak at inner wall, << Y
            Grid g(Nr,1,1,dx);
            for(int i=0;i<Nr;i++){int c=g.idx(i,0,0); g.r[c]=MAT.rho0; g.E[c]=0; double rr=(i+0.5)*dx;
                if(i>=ia&&i<=ib){ g.Sxx[c]=-Bc/(rr*rr); g.Syy[c]=Bc/(rr*rr); g.Szz[c]=0; } }
            DU d(Nr); Lop(g,d);
            double resmax=0,termmax=0,baremax=0;
            for(int i=ia+2;i<=ib-2;i++){int c=g.idx(i,0,0); double rr=(i+0.5)*dx;
                resmax=max(resmax,fabs(d.mu[c]));                       // full residual (should cancel to truncation)
                termmax=max(termmax,fabs((g.Sxx[c]-g.Syy[c])/rr));      // scale of the geometric source being cancelled
                baremax=max(baremax,fabs(d.mu[c]-(g.Sxx[c]-g.Syy[c])/rr)); } // residual WITHOUT the geometric term (~the bare dS_rr/dr)
            B=(resmax<1e-2*termmax);
            printf("Lame (B equilib): max|d.mu|=%.3e  geom-term scale=%.3e (ratio %.2e)  vs no-geom residual=%.3e\n",resmax,termmax,resmax/termmax,baremax); }
        printf("GATE (cylindrical strength: hoop strain rate exact & Lame equilibrium balances): %s\n",(A&&B)?"PASS":"CHECK");
    } else if(mode=="af_visc"){ // Phase-2b item 5: cylindrical AF Newtonian viscosity. Difference af=1 vs af=0 to ISOLATE the af-dependent viscous term (hydro/pressure are af-independent -> cancel).
        MAT=Material::basalt(); AXISYM=true; double dx=500.0; ETA_AF=1e9;
        int Nr=120; double A=1e-9;   // u(r)=A*r^2 -> discrete cylindrical (lap v)_r = 3A exactly (Cartesian-only would give 2A)
        auto setup=[&](Grid&g){ for(int i=0;i<Nr;i++){int c=g.idx(i,0,0); g.r[c]=MAT.rho0; g.E[c]=0; double r=(i+0.5)*dx; g.mu[c]=MAT.rho0*A*r*r;} };
        Grid g1(Nr,1,1,dx); setup(g1); for(int i=0;i<Nr;i++) g1.af[g1.idx(i,0,0)]=1.0;   // fluidized
        Grid g0(Nr,1,1,dx); setup(g0);                                                    // af=0 (default)
        DU d1(Nr),d0(Nr); Lop(g1,d1); Lop(g0,d0);
        int c=Nr/2; double visc=d1.mu[c]-d0.mu[c], an_cyl=ETA_AF*3.0*A, an_cart=ETA_AF*2.0*A;
        bool P=(fabs(visc-an_cyl)<1e-3*fabs(an_cyl));
        printf("AF visc (cylindrical): isolated viscous d.mu=%.6e  analytic cyl 3A=%.6e (Cartesian-only 2A=%.6e)  rel err %.2e\n",visc,an_cyl,an_cart,fabs(visc-an_cyl)/fabs(an_cyl));
        printf("GATE (cylindrical AF viscosity (lap v)_r = 3A not 2A): %s\n",P?"PASS":"CHECK");
    } else if(mode=="vib_advect"){ // Phase-2b item 6: vib (AF vibrational velocity) advects with the flow. Uniform v0 -> a vib bump translates at v0; sum(vib) conserved (conservative for uniform v).
        MAT=Material::basalt(); double dx=500.0,CFL=0.4; int N=200; double rho0=MAT.rho0,v0=2000.0;
        Grid g(N,1,1,dx); double x0=0.3*N*dx,wid=8e3;
        for(int i=0;i<N;i++){int c=g.idx(i,0,0); g.r[c]=rho0; g.mu[c]=rho0*v0; g.E[c]=0.5*rho0*v0*v0; double x=(i+0.5)*dx; g.vib[c]=exp(-((x-x0)/wid)*((x-x0)/wid));}
        // C_ACT=0,TDEC=0 -> update_af no-op -> vib is purely advected by the strength solver (no seed/decay)
        double s0=0,cen0=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); s0+=g.vib[c]; cen0+=g.vib[c]*(i+0.5)*dx;} cen0/=s0;
        double tend=20e3/v0,t=0; while(t<tend){double dt=CFL*dx/maxspeed(g); if(t+dt>tend)dt=tend-t; step_rk2(g,dt); t+=dt;}
        double s1=0,cen1=0,vmax=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); s1+=g.vib[c]; cen1+=g.vib[c]*(i+0.5)*dx; vmax=max(vmax,g.vib[c]);}
        cen1/=s1; double vcen=(cen1-cen0)/t, merr=fabs(s1-s0)/s0, verr=fabs(vcen-v0)/v0;
        printf("vib advect: sum(vib) %.4f->%.4f (err %.2e); centroid v=%.1f vs v0=%.0f (err %.2e); peak vib=%.3f\n",s0,s1,merr,vcen,v0,verr,vmax);
        printf("GATE (vib advects: sum conserved <1%% & centroid v within 2%% & peak positive): %s\n",(merr<0.01&&verr<0.02&&vmax>0)?"PASS":"CHECK");
    } else if(mode=="friction"){ // Phase 3: pressure-dependent ROCK yield (Lundborg intact + cohesionless damaged friction). Validate the capped stress vs the analytic yield surface.
        MAT=Material::basalt(); ROCK=true; Y_D0=5.0e6; bool pass=true;   // nonzero damaged cohesion to exercise the Y_d0 term
        printf("ROCK yield: Y_I0=%.2e MU_I=%.2f MU_D=%.2f Y_D0=%.2e Y_M=%.2e\n",Y_I0,MU_I,MU_D,Y_D0,Y_M);
        for(double mu : {0.001,0.005,0.02,0.05}){            // compression -> confining pressure P (condensed Tillotson, e=0)
            for(int dam=0;dam<2;dam++){
                Grid g(1,1,1,500.0); int c=0; g.r[c]=MAT.rho0*(1.0+mu); g.E[c]=0; g.D[c]=(dam?1.0:0.0); g.af[c]=0;
                double P=Pcell(g,c); g.Sxy[c]=5.0e9;          // pure shear well above yield -> vonmises returns it to the surface
                vonmises(g);
                double J2=0.5*(g.Sxx[c]*g.Sxx[c]+g.Syy[c]*g.Syy[c]+g.Szz[c]*g.Szz[c])+g.Sxy[c]*g.Sxy[c]+g.Sxz[c]*g.Sxz[c]+g.Syz[c]*g.Syz[c];
                double vm=sqrt(3.0*J2);
                double Yi=Y_I0+MU_I*P/(1.0+MU_I*P/(Y_M-Y_I0)), Yd=min(Y_D0+MU_D*P,Y_M), Yan=(dam?Yd:Yi);
                double err=fabs(vm-Yan)/Yan; if(err>1e-6) pass=false;
                printf("  P=%.3e D=%d: capped vm=%.4e  analytic Y=%.4e (%s, err %.1e)\n",P,dam,vm,Yan,dam?"damaged friction":"intact Lundborg",err);
            }
        }
        printf("GATE (ROCK yield: capped stress follows Y_i(P) intact / Y_d(P) damaged friction): %s\n",pass?"PASS":"CHECK");
    } else if(mode=="crater"){ // Phase 3: axisymmetric vertical impact crater (basalt impactor -> basalt half-space under gravity, strength, AF, free surface). Composes M1-Phase2b.
        // args: crater [a_m=500] [U_ms=12000] [TDEC=0] [ETA_AF=0] [g=3.71] [cppr=5] [tend_s=auto] [Rfac=18] [Zfac=22]
        MAT=Material::basalt(); double rho0=MAT.rho0, Abulk=MAT.A;
        double a   = argc>2?atof(argv[2]):500.0;     // impactor radius (m)
        double U   = argc>3?atof(argv[3]):12000.0;   // vertical impact speed (m/s)
        TDEC       = argc>4?atof(argv[4]):0.0;        // AF decay time (s); 0 = AF off
        ETA_AF     = argc>5?atof(argv[5]):0.0;        // AF viscosity (Pa.s)
        GZ         = argc>6?atof(argv[6]):3.71;       // gravity (Mars)
        double cppr= argc>7?atof(argv[7]):5.0;        // cells per impactor radius
        double Rfac= argc>9?atof(argv[9]):30.0, Zfac=argc>10?atof(argv[10]):22.0;   // domain size in impactor radii (Rfac 30: keep undisturbed far-field so the datum isn't contaminated by the spreading rim)
        double dx=a/cppr; AXISYM=true; RHO_CFL=100.0; RHO_VAC=100.0; double CFL=0.4;
        C_ACT=(TDEC>0?0.5:0.0); P_COH=1.0e6; P_ACT=(TDEC>0?1.0e8:0.0);   // shock-gate AF activation (>100 MPa jump): the impact shock fluidizes, slow collapse flow does NOT re-fluidize Eulerian wall cells
        ROCK=(argc>12?atoi(argv[12]):1);   // pressure-dependent friction yield (default on); arg 0 = old cohesion-only von Mises for A/B
        if(argc>13) Y_D0=atof(argv[13]);   // damaged residual cohesion (Pa); the key knob to arrest small-crater creep
        if(argc>14) VOID_AMB=atof(argv[14]);   // arg14>0 (e.g. 0.27): void cells reset to AMBIENT, not the lithostatic reference (A/B test of the below-surface refill artifact); default 0 = legacy
        if(argc>15) Y_BINGHAM=atof(argv[15]);  // arg15>0 (Pa): AF Bingham floor -- fluidized cells retain this flow resistance; default 0 = legacy
        if(argc>16) VISC_IMP=atoi(argv[16])!=0;   // arg16=1: Option C implicit viscous update (eta unbounded by the explicit dt limit)
        int Nr=(int)(Rfac*a/dx+0.5), Nz=(int)(Zfac*a/dx+0.5);
        double Zdom=Nz*dx, above=5.0*a, zsurf=Zdom-above;   // free surface; 'above' = room for impactor + ejecta
        Grid g(Nr,1,Nz,dx);
        for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){double z=(k+0.5)*dx;int c=g.idx(i,0,k);
            if(z<zsurf){g.r[c]=rho0*(1.0+rho0*GZ*(zsurf-z)/Abulk);} else g.r[c]=0.27; g.E[c]=0;}   // lithostatic basalt below, ambient above
        set_ref(g);   // route 1: freeze the pre-impact lithostatic+ambient column as the WB reference
        double zc=zsurf+a;   // impactor sphere on the axis, tangent to the surface, moving down
        for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){double r=(i+0.5)*dx,z=(k+0.5)*dx;int c=g.idx(i,0,k);
            if(sqrt(r*r+(z-zc)*(z-zc))<a){g.r[c]=rho0;g.mw[c]=-rho0*U;g.E[c]=0.5*rho0*U*U;}}
        double targ=argc>8?atof(argv[8]):-1.0;   // explicit tend>0 -> fixed run (old behaviour); <=0 -> run-to-settling
        double t_auto=2.0*sqrt((Rfac*a)/GZ);      // gravity formation/collapse timescale
        bool fixed=(targ>0); double tend=fixed?targ:6.0*t_auto;   // settling cap = 6x t_auto (depth ~converged by ~5x; keeps run-to-settling affordable)
        double tolM=0.02;                         // settled when consecutive Vexc WINDOW-MEANS agree within 2% (robust to residual sloshing + density-threshold jitter)
        double Wwin=2.0*t_auto;                    // averaging window >= one oscillation period (2pi*sqrt(D/g)) so fast jitter/sloshing averages out
        // Column-mass surface metric (Option C hygiene, 2026-07-07): depth = per-column MASS DEFICIT vs the frozen
        // pre-impact reference, / rho0. Threshold-free: a smoothly dilating sub-floor column reads as gradually
        // increasing depth, never a jump -- the a=3000 "40 km breakthrough" was the old first-rho>1350 surface
        // crossing such a column (2026-07-06 diagnostic). Negative = uplift (rim) / mass overhead (ejecta).
        auto dcol=[&](int i){ double m=0; for(int k=0;k<Nz;k++){ int c=g.idx(i,0,k); m+=REF_R0[c]-g.r[c]; } return m*dx/rho0; };
        auto vmax_dense=[&](){ double v=0; for(uint32_t c=0;c<g.r.size();c++) if(g.r[c]>1350.0){ double vv=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; v=max(v,vv);} return v; };
        auto vexc=[&](){ double V=0; for(int i=0;i<Nr;i++){ double d=dcol(i); if(d>0) V+=d*(i+0.5); } return V; };   // r-weighted excavated cross-section ~ crater volume (logged for info)
        auto dfloor=[&](){ double V=0; for(int i=0;i<Nr;i++){ double d=dcol(i); if(d>V)V=d; } return V; };   // max excavation depth (bowl floor near the axis): the SETTLING signal -- insensitive to peripheral lateral spreading + far-field/boundary substrate sag (which corrupt Vexc & the datum over long runs)
        auto afstat=[&](double&amax){ double s=0;long c2=0;amax=0; for(uint32_t c=0;c<g.r.size();c++) if(g.r[c]>1350.0){ s+=g.af[c]; amax=max(amax,(double)g.af[c]); c2++; } return c2?s/c2:0.0; };   // AF diagnostic: is the rock still fluidized late in the run?
        auto Dstat=[&](double&dfrac){ double s=0;long c2=0,nhi=0; for(uint32_t c=0;c<g.r.size();c++) if(g.r[c]>1350.0){ s+=g.D[c]; if(g.D[c]>0.95)nhi++; c2++; } dfrac=c2?(double)nhi/c2:0.0; return c2?s/c2:0.0; };   // damage diagnostic: mean D + fraction fully-damaged (D>0.95) -> does fluidized collapse drive D->1?
        { const char*sd=getenv("CRATER_SNAP_DT"); if(sd){ SNAP_DT=atof(sd); SNAP_DIR=getenv("CRATER_SNAP_DIR"); } }
        // Lagrangian tracer lattice (env CRATER_TRACE=1; diagnostic, Collins et al. 2020 Fig.2 style):
        // one tracer per cell in the pre-impact target, advected by the bilinearly interpolated
        // velocity, each recording its own peak pressure. Written as tracers_XXXXXXX.txt next to
        // the field snapshots (cols: id i0 k0 r_m z_m Pmax).
        bool TRACE = getenv("CRATER_TRACE") && SNAP_DT>0 && SNAP_DIR;
        vector<double> tx,tz,tpm,tvr,tvz; vector<int> ti0,tk0;
        int ksurf=(int)(zsurf/dx);
        if(TRACE){ for(int i=0;i<Nr;i++)for(int k=0;k<ksurf;k++){
                tx.push_back((i+0.5)*dx); tz.push_back((k+0.5)*dx); tpm.push_back(0.0);
                tvr.push_back(0.0); tvz.push_back(0.0);
                ti0.push_back(i); tk0.push_back(k); } }
        auto interp=[&](const vector<double>&F,bool byrho,double x,double z)->double{   // bilinear at (x,z), cell-centered
            double fi=x/dx-0.5, fk=z/dx-0.5;
            int i0=(int)floor(fi), k0=(int)floor(fk); double ai=fi-i0, ak=fk-k0;
            i0=max(0,min(Nr-2,i0)); k0=max(0,min(Nz-2,k0));
            auto val=[&](int i,int k){ int c=g.idx(i,0,k); return byrho? F[c]/max(g.r[c],1e-30) : F[c]; };
            return (1-ai)*(1-ak)*val(i0,k0)+ai*(1-ak)*val(i0+1,k0)+(1-ai)*ak*val(i0,k0+1)+ai*ak*val(i0+1,k0+1); };
        double dmaxT=0; double t=0;int s=0; double t_settled=-1.0,next_log=t_auto; double snap_next=SNAP_DT;
        double sumV=0.0,meanPrev=-1.0,winStart=-1.0; long cntV=0;   // non-overlapping window-mean drift detector
        while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);
            for(int i=0;i<Nr;i++)for(int k=0;k<2;k++){int c=g.idx(i,0,k);g.r[c]=REF_R0[c];g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=0;}  // pin deep far-field floor
            { int isp=(int)(0.85*Nr);   // far-field radial SPONGE: relax the outer ~15% toward the WB reference each step -> absorbs the slow boundary-edge sink + outgoing waves (the crater sits in the flat interior, r<~0.5 r_max)
              for(int i=isp;i<Nr;i++){ double xi=(double)(i-isp)/max(1,Nr-1-isp), sp=xi*xi;   // ramp 0 (inner edge of sponge) -> 1 (domain edge)
                for(int k=0;k<Nz;k++){ int c=g.idx(i,0,k); g.r[c]=(1-sp)*g.r[c]+sp*REF_R0[c]; g.mu[c]*=(1-sp); g.mv[c]*=(1-sp); g.mw[c]*=(1-sp); g.E[c]*=(1-sp); } } }
            if(s%20==0){double dnow=dfloor();dmaxT=max(dmaxT,dnow);}   // track transient excavation
            if(TRACE){ for(size_t j=0;j<tx.size();j++){
                    int cj=g.idx(min((int)(tx[j]/dx),Nr-1),0,min((int)(tz[j]/dx),Nz-1));
                    if(g.r[cj]>300.0){ tvr[j]=interp(g.mu,true,tx[j],tz[j]); tvz[j]=interp(g.mw,true,tx[j],tz[j]); }
                    else tvz[j]-=GZ*dt;                                     // in void: ballistic flight
                    double vr=tvr[j], vz=tvz[j];
                    tx[j]+=vr*dt; tz[j]+=vz*dt;
                    if(tx[j]<0) tx[j]=-tx[j];                               // axisym reflection
                    tx[j]=min(tx[j],(Nr-0.51)*dx); tz[j]=max(0.51*dx,min(tz[j],(Nz-0.51)*dx));
                    int ci=(int)(tx[j]/dx), ck=(int)(tz[j]/dx);
                    double Pt=Pcell(g,g.idx(min(ci,Nr-1),0,min(ck,Nz-1)));
                    if(Pt>tpm[j]) tpm[j]=Pt; } }
            t+=dt;s++;
            if(t>=next_log){double vn=vmax_dense(),dn=dfloor(),amax,dfrac;double amean=afstat(amax),Dmean=Dstat(dfrac);   // progress to stderr (sweep DEVNULLs stderr); watch settling + AF + damage state
                fprintf(stderr,"  [crater a=%.0f] t=%.1f/%.0fs steps=%d max|v|=%.2f Vexc=%.4e depth=%.2fkm af=%.3f/%.3f D(mean,frac>.95)=%.3f,%.2f seeds=%ld voidrst=%ld\n",a,t,tend,s,vn,vexc(),dn/1e3,amean,amax,Dmean,dfrac,SEED_N,VZ_N); next_log+=20.0; }
            if(SNAP_DT>0 && SNAP_DIR && t>=snap_next){   // diagnostic field snapshot (axisym r-z slice)
                char fn[512]; snprintf(fn,sizeof fn,"%s/snap_%07d.txt",SNAP_DIR,(int)(t+0.5));
                FILE*fs=fopen(fn,"w");
                if(fs){ fprintf(fs,"# t=%.2f nx=%d nz=%d dx=%.1f cols: i k rho P vib af D Pmax vr vz eint\n",t,Nr,Nz,dx);
                    for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){ int c=g.idx(i,0,k); double ir=1.0/max(g.r[c],1e-30);
                        fprintf(fs,"%d %d %.4e %.4e %.4e %.4e %.3f %.4e %.3e %.3e %.3e\n",i,k,g.r[c],Pcell(g,c),g.vib[c],g.af[c],g.D[c],g.Pmax[c],g.mu[c]*ir,g.mw[c]*ir,eint(g,c)); }
                    fclose(fs); }
                if(TRACE){ char ft[512]; snprintf(ft,sizeof ft,"%s/tracers_%07d.txt",SNAP_DIR,(int)(t+0.5));
                    FILE*f2=fopen(ft,"w");
                    if(f2){ fprintf(f2,"# t=%.2f n=%zu dx=%.1f ksurf=%d cols: id i0 k0 r_m z_m Pmax\n",t,tx.size(),dx,ksurf);
                        for(size_t j=0;j<tx.size();j++) fprintf(f2,"%zu %d %d %.2f %.2f %.4e\n",j,ti0[j],tk0[j],tx[j],tz[j],tpm[j]);
                        fclose(f2); } }
                snap_next+=SNAP_DT; }
            if(!fixed && t>t_auto && s%5==0){double V=dfloor();   // run-to-settling: compare successive window-means of the BOWL FLOOR DEPTH (settles fast ~t_auto; stops before substrate sag accrues)
                if(winStart<0)winStart=t; sumV+=V; cntV++;
                if(t-winStart>=Wwin){ double meanNow=sumV/cntV;
                    if(meanPrev>0 && fabs(meanNow-meanPrev)<tolM*meanNow){ t_settled=t; break; }   // two consecutive window-means agree -> settled
                    meanPrev=meanNow; sumV=0; cntV=0; winStart=t; } }
        }
        double vmx=0,pmx=0; for(uint32_t c=0;c<g.r.size();c++){if(g.r[c]>1350){double v=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c];vmx=max(vmx,v);pmx=max(pmx,Pcell(g,c));}}
        // final geometry from the column-mass metric (datum = mass deficit 0, exact for undisturbed columns)
        double dapp=0; int ifloor=0; for(int i=0;i<Nr;i++){double d=dcol(i); if(d>dapp){dapp=d;ifloor=i;}}   // deepest point
        int iD=Nr-1; for(int i=ifloor;i<Nr;i++){ if(dcol(i)<=0.25*dx){ iD=i; break; } }   // apparent crater radius: mass surface back to the datum, scanning out from the floor
        double rD=(iD+0.5)*dx, Dapp=2.0*rD, dD=(Dapp>0?dapp/Dapp:0.0);
        double zrim=0,rrim=0; for(int i=0;i<Nr;i++){double u=-dcol(i); if(u>zrim){zrim=u;rrim=(i+0.5)*dx;}}   // rim crest (mass uplift above datum)
        bool finite=isfinite(vmx)&&isfinite(pmx)&&isfinite(dapp);
        printf("CRATER a=%.0fm U=%.0f g=%.2f TDEC=%.3g ETA=%.3g%s | grid %dx%d dx=%.0f zsurf=%.1fkm tend=%.1fs steps=%d\n",a,U,GZ,TDEC,ETA_AF,VISC_IMP?" IMPLICIT":"",Nr,Nz,dx,zsurf/1e3,tend,s);
        printf("  settling: %s at t=%.1fs (t_auto=%.1fs cap=%.1fs tolM=%.1f%% Wwin=%.1fs final max|v|=%.1f m/s)\n",(t_settled>0?"SETTLED":(fixed?"FIXED-RUN":"NOT-SETTLED@cap")),(t_settled>0?t_settled:t),t_auto,tend,tolM*100,Wwin,vmx);
        printf("  D_app=%.2f km  depth(below datum)=%.2f km  transient depth=%.2f km  rim uplift=%.2f km @ r=%.2f km  d/D=%.3f  (column-mass metric)\n",Dapp/1e3,dapp/1e3,dmaxT/1e3,zrim/1e3,rrim/1e3,dD);
        printf("  stability: max|v|=%.1f m/s  maxP=%.2e Pa  %s\n",vmx,pmx,finite?"FINITE":"NONFINITE-BLOWUP");
        { double mcell=VZ_DM*2*3.14159265358979*0.5*dx*dx*dx;   // injected mass (kg, full 2pi; r~0.5dx underestimates outer cells -- order-of-magnitude only)
          double mtr=2700.0*4.18879*a*a*a;   // impactor mass scale for context
          printf("  voidzero: below-surface refills=%ld  injected rho-sum=%.3e (order ~%.2e kg vs impactor %.2e kg)  mode=%s\n",VZ_N,VZ_DM,mcell,mtr,VOID_AMB>0?"AMBIENT":"LEGACY-REF"); }
        printf("RESULT %.1f %.4f %.4f %.4f %.4f\n",a,Dapp,dapp,dmaxT,dD);   // machine-readable (preliminary, in-loop): a D_app depth transient d/D (metres)
        const char*prof=argc>11?argv[11]:"crater_profile.txt";   // dump the final surface profile z_s(r)-z0 for ROBUST offline depth-diameter measurement
        FILE*pf=fopen(prof,"w"); if(pf){ fprintf(pf,"# r_m  -dcol_m (column-mass surface vs datum)   (a=%.0f U=%.0f g=%.2f TDEC=%.3g ETA=%.3g%s dx=%.0f zsurf=%.1f)\n",a,U,GZ,TDEC,ETA_AF,VISC_IMP?" IMP":"",dx,zsurf);
            for(int i=0;i<Nr;i++) fprintf(pf,"%.1f %.3f\n",(i+0.5)*dx, -dcol(i)); fclose(pf); }
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
    } else if(mode=="vacuum"){ // Toro expansion-into-vacuum: material at rest (left) expands into vacuum (right) -> centred rarefaction
        MAT=Material::ideal(1.4); GAM=1.4; RHO_VAC=1e-3;
        int N=400;double dx=1.0/N,CFL=0.4,tend=0.10; double rL=1.0,pL=1.0,cL=sqrt(GAM*pL/rL); Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;int c=g.idx(i,0,0); if(x<0.5){g.r[c]=rL;g.E[c]=pL/(GAM-1);} else {g.r[c]=1e-6;g.E[c]=1e-6/(GAM-1);}}
        double t=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;}
        double l1r=0,l1u=0,rmin=1e30;int np=0;   // compare vs analytic left-rarefaction-into-vacuum (u_L=0)
        for(int i=0;i<N;i++){int c=g.idx(i,0,0);double x=(i+0.5)*dx,S=(x-0.5)/t,ra,ua;
            if(S<=-cL){ra=rL;ua=0;} else if(S<2*cL/(GAM-1)){double u=2.0/(GAM+1)*(cL+S),cc=2.0/(GAM+1)*(cL-0.5*(GAM-1)*S);ra=rL*pow(cc/cL,2.0/(GAM-1));ua=u;} else {ra=0;ua=0;}
            rmin=min(rmin,g.r[c]); if(ra>0.02){l1r+=fabs(g.r[c]-ra);l1u+=fabs(g.mu[c]/g.r[c]-ua);np++;} }
        l1r/=np;l1u/=np;
        printf("Vacuum expansion (Toro): L1 rho=%.4f u=%.4f (npts=%d), min rho=%.2e\n",l1r,l1u,np,rmin);
        printf("GATE (matches Toro fan, L1 rho<0.03 & u<0.05 & positive): %s\n",(l1r<0.03&&l1u<0.05&&rmin>0)?"PASS":"CHECK");
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
    } else if(mode=="pierazzo2d"){ // Pierazzo et al. 2008 Al-on-Al peak-shock-pressure benchmark, AXISYMMETRIC (r,z).
        // Paper cases (oracle numbers in pierazzo2008_oracle.md): 1 km-DIAMETER Al sphere, vertical, U=5 or
        // 20 km/s, strengthless. No gravity -> scale-free -> a=1 m here; paper distance (km = proj DIAMETERS)
        // maps to d/a = 2*d_km. Metrics vs Table 1: P_cc = mean peak P over depth 0-0.6a (5 km/s) / 0-1.2a
        // (20 km/s) below the impact point; decay slope n over d/a in [4,10] (the paper's 2-5 diameters).
        // Caveat: the paper records LAGRANGIAN tracer peaks; we record the EULERIAN per-cell peak. Equivalent
        // in the fit window (tracer displacement <~1% of distance there, per the paper) but softer in the CC zone.
        // args: pierazzo2d [cppr=12] [reach=12] [U=5000] [tend=-1(auto)]   (radial halfwidth = reach: half>=reach by construction)
        MAT=Material::aluminum(); AXISYM=true; double CFL=0.3;
        double cppr=argc>2?atof(argv[2]):12.0, reach=argc>3?atof(argv[3]):12.0, U=argc>4?atof(argv[4]):5000.0;
        double a=1.0, dx=a/cppr, rho0=MAT.rho0, zsurf=reach*a, zc=zsurf+a;
        int Nr=(int)(reach*a/dx+0.5), Nz=(int)((reach+3.0)*a/dx+0.5);
        Grid g(Nr,1,Nz,dx);
        for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){double r=(i+0.5)*dx,z=(k+0.5)*dx;int c=g.idx(i,0,k);
            if(sqrt(r*r+(z-zc)*(z-zc))<a){g.r[c]=rho0;g.mw[c]=-rho0*U;g.E[c]=0.5*rho0*U*U;}   // projectile, tangent to the surface, moving down
            else if(z<zsurf){g.r[c]=rho0;g.E[c]=0;}                                            // Al half-space
            else {g.r[c]=0.27;g.E[c]=0;}}                                                      // low-density ambient (M2 free-surface treatment; matches the 3D pierazzo mode)
        double targ=argc>5?atof(argv[5]):-1.0, c0=sqrt(MAT.A/rho0);
        double tend=(targ>0?targ:1.3*reach*a/c0);   // the front never falls below the bulk speed c0 -> 1.3x covers the decelerating tail
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);
            for(uint32_t c=0;c<g.r.size();c++){double P=Pcell(g,c);if(P>g.Pmax[c])g.Pmax[c]=P;}   // Eulerian peak-P tracking (update_af is a no-op here: AF off)
            t+=dt;s++;}
        double ccd=(U>=12500.0?1.2:0.6);   // contact/compression averaging depth in a (paper Table 1: 0-300 m / 0-600 m of the 500 m radius)
        double Pcore=0; for(int k=0;k<Nz;k++){double d=zsurf-(k+0.5)*dx;if(d>0&&d<=1.5*a)Pcore=max(Pcore,g.Pmax[g.idx(0,0,k)]);}
        double Pcc=0;int ncc=0; double Sx=0,Sy=0,Sx2=0,Sy2=0,Sxy=0;int np=0; double P4D=0,best=1e30;
        FILE*o=fopen("pierazzo2d_decay.txt","w");
        if(o)fprintf(o,"# d_over_a  Ppeak_Pa   (U=%.0f cppr=%g reach=%g axisym CPU)\n",U,cppr,reach);
        for(int k=Nz-1;k>=0;k--){double d=zsurf-(k+0.5)*dx; if(d<=0)continue; double ra=d/a,P=g.Pmax[g.idx(0,0,k)];
            if(o&&P>0)fprintf(o,"%.4f %.6e\n",ra,P);
            if(ra<=ccd){Pcc+=P;ncc++;}
            if(ra>=4.0&&ra<=10.0&&P>0.005*Pcore){double lx=log(ra),ly=log(P);Sx+=lx;Sy+=ly;Sx2+=lx*lx;Sy2+=ly*ly;Sxy+=lx*ly;np++;}
            if(fabs(ra-8.0)<best){best=fabs(ra-8.0);P4D=P;}}   // spot anchor: 4 projectile diameters (paper: 3.2+-0.5 GPa at 5 km/s)
        if(o)fclose(o);
        Pcc=(ncc?Pcc/ncc:0);
        double nfit=np>2?-(np*Sxy-Sx*Sy)/(np*Sx2-Sx*Sx):0;
        double Rcor=np>2?fabs((np*Sxy-Sx*Sy)/sqrt(max((np*Sx2-Sx*Sx)*(np*Sy2-Sy*Sy),1e-300))):0;
        double Phug=hugoniot_P(0.5*U);   // planar impedance-match reference (up=U/2, same material both sides)
        printf("Pierazzo-2008 Al axisym (U=%.0f km/s cppr=%g reach=%g | %dx%d cells, %d steps, tend=%.3e):\n",U/1000,cppr,reach,Nr,Nz,s,tend);
        printf("  P_cc(0-%.1fa)=%.1f GPa  n[d/a 4-10]=%.3f (R=%.4f np=%d)  P(4 diam)=%.2f GPa  [planar Hugoniot %.1f GPa]\n",ccd,Pcc/1e9,nfit,Rcor,np,P4D/1e9,Phug/1e9);
        printf("PZ2D U=%.0f cppr=%g reach=%g Pcc=%.4e n=%.4f R=%.4f P4D=%.4e steps=%d\n",U,cppr,reach,Pcc,nfit,Rcor,P4D,s);
        if(fabs(U-5000)<1||fabs(U-20000)<1){ bool five=(fabs(U-5000)<1);   // gate vs the PUBLISHED inter-code range (Table 1; the paper's own pass criterion is "within the code spread")
            double plo=five?28.4e9:334.7e9, phi=five?48.0e9:411.1e9, nlo=five?1.13:2.27, nhi=five?1.41:2.53;
            printf("  paper Table 1 range: P_cc [%.1f,%.1f] GPa (mean %s)  n [%.2f,%.2f] (mean %s)\n",plo/1e9,phi/1e9,five?"40.4+-6.2":"379+-26",nlo,nhi,five?"1.2+-0.1":"2.3+-0.1");
            printf("GATE (Pierazzo-2008 inter-code envelope: P_cc & n inside the published code range): %s\n",
                (Pcc>=plo&&Pcc<=phi&&nfit>=nlo&&nfit<=nhi)?"PASS":"CHECK");
        }
    } else if(mode=="prater"){ // Prater 1970 Al-on-Al crater growth (Pierazzo-2008 validation #2; oracle = prater1970_oracle.md).
        // 6.35 mm Al 2017-T4 sphere -> quasi-infinite Al 6061-T6 at ~7 km/s, vertical -> axisym; flash-X-ray R(t), D(t).
        // Strength = von Mises Y=414 MPa (iSALE Table 6 choice), G=27.6 GPa (AUTODYN SG G0); no damage (ductile, wk=0),
        // no gravity (70 us). Projectile modeled as the target alloy (2017-T4 rho 2.79 vs 2.70 -- documented approximation).
        // Expected von Mises signature (paper Fig. 14): radius within a few %, depth OVERestimated up to ~20%.
        // args: prater [cppr=10] [U=7000] [Y=414e6] [G=27.6e9] [Rfac=8] [tend=70e-6] [rcfl=1350]
        MAT=Material::aluminum(); AXISYM=true; RHO_VAC=100.0; double CFL=0.3;
        double cppr=argc>2?atof(argv[2]):10.0, U=argc>3?atof(argv[3]):7000.0;
        MAT.Y=argc>4?atof(argv[4]):414e6; MAT.G=argc>5?atof(argv[5]):27.6e9;
        double Rfac=argc>6?atof(argv[6]):8.0, tend=argc>7?atof(argv[7]):70e-6;
        RHO_CFL=argc>8?atof(argv[8]):100.0;   // 100 = crater-mode-validated. NOTE: raising it to 1350 (half-density) cuts steps ~50x BUT corrupts the solution (CFL-violating skirt cells refill the cavity) -- tested 2026-07-02, do not
        double a=3.175e-3, dx=a/cppr, rho0=MAT.rho0;
        int Nr=(int)(Rfac*a/dx+0.5), Nz=(int)((Rfac+3.0)*a/dx+0.5);
        double zsurf=Rfac*a, zc=zsurf+a; int ksurf=(int)(zsurf/dx+0.5);
        Grid g(Nr,1,Nz,dx);
        for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){double r=(i+0.5)*dx,z=(k+0.5)*dx;int c=g.idx(i,0,k);
            if(sqrt(r*r+(z-zc)*(z-zc))<a){g.r[c]=rho0;g.mw[c]=-rho0*U;g.E[c]=0.5*rho0*U*U;}   // projectile, tangent to the surface
            else if(z<zsurf){g.r[c]=rho0;g.E[c]=0;}                                            // Al target half-space
            else {g.r[c]=0.27;g.E[c]=0;}}                                                      // low-density ambient (vacuum-flagged by RHO_VAC)
        double cut=0.5*rho0;   // crater boundary = density cutoff (paper codes used the same idea; Table 5 note)
        auto sctg=[&](int i){ int k=0; while(k<Nz && g.r[g.idx(i,0,k)]>=cut) k++; return k*dx; };   // ground level = top of the dense column CONTIGUOUS from the bottom (flying ejecta/jet debris ignored)
        auto RD=[&](double&R,double&D){ double zfl=zsurf; int ifl=0; for(int i=0;i<Nr;i++){ double zs=sctg(i); if(zs<zfl){zfl=zs;ifl=i;} }
            D=zsurf-sctg(0); R=Nr*dx;   // depth on the axis (X-ray convention); radius = profile recovers to the original surface plane, scanning out from the floor
            for(int i=ifl;i<Nr;i++){ if(sctg(i)>=zsurf-0.15*a){ R=(i+0.5)*dx; break; } } };   // 0.15a tolerance dodges the sub-cell free-surface dip that flickers a single fixed row
        FILE*o=fopen("prater_growth.txt","w");
        if(o)fprintf(o,"# t_us R_cm D_cm (U=%.0f Y=%.3e G=%.3e cppr=%g CPU; exp oracle = prater1970_oracle.md Table 4)\n",U,MAT.Y,MAT.G,cppr);
        double t=0,tlog=0,Rsum=0,Dsum=0,Rmx=0,Dmx=0;int s=0,nfin=0;
        while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;
            for(uint32_t c=0;c<g.r.size();c++){ if(g.r[c]<RHO_CFL){ g.r[c]=0.27;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=0;   // void cleanup to AMBIENT (mirrors GPU voidzero): without it, sub-100 junk accretes onto the crater floor and bridges the axis (tested: D -> 0 at ~28 us)
                g.Sxx[c]=g.Syy[c]=g.Szz[c]=g.Sxy[c]=g.Sxz[c]=g.Syz[c]=0; } }
            if(t>=tlog){double R,D;RD(R,D);if(o)fprintf(o,"%.3f %.4f %.4f\n",t*1e6,R*100,D*100);tlog+=0.5e-6;
                Rmx=max(Rmx,R);Dmx=max(Dmx,D);
                if(t>=tend-10e-6){Rsum+=R;Dsum+=D;nfin++;}}}   // final values = mean over the last 10 us (experiment plateaus by ~25 us)
        if(o)fclose(o);
        FILE*pf=fopen("prater_profile.txt","w");   // final crater profile (paper-figure material; mirrors the X-ray cross-sections)
        if(pf){ fprintf(pf,"# r_cm  (surf-z0)_cm  (CPU U=%.0f Y=%.3e cppr=%g t=%.0fus)\n",U,MAT.Y,cppr,tend*1e6);
            for(int i=0;i<Nr;i++) fprintf(pf,"%.4f %.4f\n",(i+0.5)*dx*100,(sctg(i)-zsurf)*100); fclose(pf); }
        double Rf=nfin?Rsum/nfin:0, Df=nfin?Dsum/nfin:0, Rexp=1.31e-2, Dexp=1.31e-2;   // 6061-T6 plateau values (Table 4)
        printf("Prater-1970 Al crater growth (U=%.0f km/s cppr=%g Y=%.0f MPa | %dx%d cells, %d steps, tend=%.0f us):\n",U/1000,cppr,MAT.Y/1e6,Nr,Nz,s,tend*1e6);
        printf("  final R=%.3f cm (exp %.2f, err %+.1f%%)  D=%.3f cm (exp %.2f, err %+.1f%%)  transient max R=%.3f D=%.3f\n",Rf*100,Rexp*100,100*(Rf-Rexp)/Rexp,Df*100,Dexp*100,100*(Df-Dexp)/Dexp,Rmx*100,Dmx*100);
        printf("PRSUMMARY U=%.0f cppr=%g Y=%.4e G=%.4e Rf=%.4e Df=%.4e Rmx=%.4e Dmx=%.4e steps=%d\n",U,cppr,MAT.Y,MAT.G,Rf,Df,Rmx,Dmx,s);
        if(fabs(U-7000)<1&&fabs(MAT.Y-414e6)<1){   // paper-exact case. Growth phase tracks experiment inside the inter-code band; the plateau sits low (no rim uplift + slow wall
            // extrusion: the vacuum-face Riemann flux carries no deviatoric strength) -- documented limitation, same league as iSALE-von-Mises' +20% depth overshoot (Fig. 14).
            if(cppr>=10) printf("GATE (Prater-1970 6061-T6 von Mises: Dmax within 10%% & final D in [-20%%,+25%%] & final R in [-30%%,+10%%]): %s\n",
                (fabs(Dmx-Dexp)<=0.10*Dexp && Df-Dexp>=-0.20*Dexp && Df-Dexp<=0.25*Dexp && Rf-Rexp>=-0.30*Rexp && Rf-Rexp<=0.10*Rexp)?"PASS":"CHECK");
            else printf("GATE (prater smoke, under-resolved cppr<10: crater forms, Dmax>0.7*exp & final R>0.35*exp): %s\n",
                (Dmx>0.7*Dexp && Rf>0.35*Rexp)?"PASS":"CHECK");
        }
    } else if(mode=="substrate"){ // route 1: well-balanced gravity-loaded substrate must stay stable (no rarefaction), no damping
        MAT=Material::basalt(); GZ=3.71; double rho0=2700,A=2.67e10; int NXc=60,NZc=60;double dx=500.0,CFL=0.4;
        double zsurf=20e3,tend=200.0; DAMP=1.0; RHO_CFL=100.0;
        Grid g(NXc,1,NZc,dx); int nb0=0;
        for(int i=0;i<NXc;i++)for(int k=0;k<NZc;k++){double z=(k+0.5)*dx;int c=g.idx(i,0,k);
            if(z<zsurf){g.r[c]=rho0*(1.0+rho0*GZ*(zsurf-z)/A);nb0++;} else g.r[c]=0.27; g.E[c]=0;}
        set_ref(g);   // route 1: freeze the lithostatic IC as the WB reference
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);
            for(int i=0;i<NXc;i++)for(int k=0;k<2;k++){int c=g.idx(i,0,k);g.r[c]=REF_R0[c];g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=0;}  // deep far-field floor pinned to the undisturbed reference
            t+=dt;s++;}
        double vmx=0;int nb=0;for(int i=0;i<NXc;i++)for(int k=0;k<NZc;k++){int c=g.idx(i,0,k);if(g.r[c]>1350){vmx=max(vmx,sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]);nb++;}}
        printf("Substrate relax (DAMP=%.2f): max|v|=%.3f m/s, basalt cells %d->%d, steps %d\n",DAMP,vmx,nb0,nb,s);
        printf("GATE (settled max|v|<5 m/s & no rarefaction nb>0.97*nb0): %s\n",(vmx<5.0&&nb>0.97*nb0)?"PASS":"CHECK");
    } else if(mode=="af_activate"){ // Phase-1 unit test for shock-activated AF (Wuennemann-Ivanov):
        // (A) a shock seeds vibrations in the rock it passes through, and not ahead of the front;
        // (B) seeded vibrations at relaxed (low) pressure fluidize the rock, which then re-solidifies as they decay.
        MAT=Material::basalt(); double CFL=0.4;
        // (A) activation + spatial selectivity: a planar basalt shock (left half -> right) seeds vib behind the front.
        int N=200; double dx=500.0, V=2000.0, tA=3.0; TDEC=10.0; C_ACT=0.5; ETA_AF=0.0;
        Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){int c=g.idx(i,0,0); g.r[c]=MAT.rho0; g.E[c]=0; if(i<N/2) g.mu[c]=MAT.rho0*V;}
        double t=0; while(t<tA){double dt=CFL*dx/maxspeed(g); if(t+dt>tA)dt=tA-t; step_rk2(g,dt); t+=dt;}
        double vib_b=0,vib_a=0; for(int i=101;i<=120;i++)vib_b=max(vib_b,g.vib[g.idx(i,0,0)]);
        for(int i=180;i<N;i++)vib_a=max(vib_a,g.vib[g.idx(i,0,0)]);
        bool A=(vib_b>1.0 && vib_a<1e-3);
        printf("AF activate (A: shock seeding): vib behind front=%.1f m/s, ahead=%.1e m/s\n",vib_b,vib_a);
        // (B) fluidization + re-solidification: seeded vibrations at P~0 give af~1, then decay back to ~0.
        int M=50; double tB=5.0, vib0=80.0; TDEC=0.5; C_ACT=0.5;
        Grid h(M,1,1,dx);
        for(int i=0;i<M;i++){int c=h.idx(i,0,0); h.r[c]=MAT.rho0; h.E[c]=0; h.vib[c]=vib0;}
        update_af(h,0.0); double af0=h.af[h.idx(M/2,0,0)];
        double tt=0; while(tt<tB){double dt=CFL*dx/maxspeed(h); if(tt+dt>tB)dt=tB-tt; step_rk2(h,dt); tt+=dt;}
        double af1=h.af[h.idx(M/2,0,0)];
        bool B=(af0>0.9 && af1<0.05);
        printf("AF activate (B: fluidize->decay): af %.2f -> %.3f over %.0fs (TDEC=0.5)\n",af0,af1,tB);
        printf("GATE (shock seeds vib behind front not ahead; fluidized rock re-solidifies as vib decays): %s\n",(A&&B)?"PASS":"CHECK");
    } else if(mode=="collapse"){ // acoustic-fluidization demo: a basalt step slumps if fluidized, holds if not
        MAT=Material::basalt(); GZ=3.71; double rho0=2700,A=2.67e10; int NXc=80,NZc=40;double dx=500.0,CFL=0.4;
        double zlo=10e3,h0=4e3,zhi=zlo+h0,tend=8.0; double hfin[2]; const char*lbl[2]={"AF-on ","AF-off"};
        for(int run=0;run<2;run++){
            TDEC=(run==0?1e6:0.0); ETA_AF=(run==0?1e9:0.0); RHO_CFL=100.0; RHO_VAC=100.0; Grid g(NXc,1,NZc,dx);   // ambient(0.27): void-CFL + vacuum-aware free surface (the vertical cliff)
            for(int i=0;i<NXc;i++)for(int k=0;k<NZc;k++){double x=(i+0.5)*dx,z=(k+0.5)*dx;int c=g.idx(i,0,k);double zs=(x<0.5*NXc*dx?zhi:zlo);
                if(z<zs){g.r[c]=rho0*(1.0+rho0*GZ*(zs-z)/A); if(run==0)g.af[c]=1.0;} else g.r[c]=0.27; g.E[c]=0;}   // cold lithostatic step
            set_ref(g);   // route 1: per-column lithostatic WB reference (vertical balance; the horizontal step drives the slump)
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
    } else if(mode=="surface"){
        MAT=Material::basalt();int N=64;double L=200e3,dx=L/N,tend=2.0,CFL=0.4;Grid g(1,1,N,dx);
        for(int k=0;k<N;k++){double z=(k+0.5)*dx;double rr=z<0.5*L?2700:0.27;int c=g.idx(0,0,k);g.r[c]=rr;g.E[c]=0;}
        double t=0;int s=0;while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double vmax=0;for(int k=0;k<N;k++){int c=g.idx(0,0,k);vmax=max(vmax,fabs(g.mw[c]/g.r[c]));}
        printf("Free surface max|v|=%.3f m/s  GATE: %s\n",vmax,vmax<5.0?"PASS":"CHECK");
    } else if(mode=="visc_diff"){ // Option C gate 3.4b: the implicit viscous solve is a correct diffusion operator at dt >> the explicit stability limit.
        // (A) Cartesian: transverse cosine modes v_y,v_z on a fluidized 1D slab. cos(k(i+.5)dx), k=m*pi/L is the
        //     discrete Neumann eigenmode: each backward-Euler step must scale it by exactly f=1/(1+nu*lam*dt),
        //     lam=(2-2cos(k dx))/dx^2. The implied nu_eff=(1/f-1)/(lam*dt) must equal nu to 1e-6 (solver exact),
        //     and the cumulative decay must track exp(-nu k^2 t) within 10% (k dx small; BE time bias ~x^2/2).
        // (B) axisym: v_r=A*r is an exact DISCRETE null mode of the cylindrical operator (the area-weighted flux
        //     divergence cancels the -v/r^2 curvature diagonal) -> near-axis cells untouched; v_z=const invariant.
        MAT=Material::ideal(1.4); VISC_IMP=true; bool okA;
        double worstnu,ampA,exact,zerr;
        { int N=256; double dx=1.0,rho0=1.0,nu=1.0; ETA_AF=nu*rho0; AXISYM=false;
          int m=4; double L=N*dx,kk=m*M_PI/L,lam=(2.0-2.0*cos(kk*dx))/(dx*dx);
          double dtv=20.0*dx*dx/(2.0*nu);   // 20x the 1D explicit stability limit
          Grid g(N,1,1,dx); double A0=1.0; int NS=12;
          for(int i=0;i<N;i++){int c=g.idx(i,0,0); g.r[c]=rho0; g.af[c]=1.0; double vy=A0*cos(kk*(i+0.5)*dx);
              g.mv[c]=rho0*vy; g.mw[c]=0.5*rho0*vy; g.E[c]=1e5+0.5*(g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/rho0; }
          double fex=1.0/(1.0+nu*lam*dtv), amp=A0; worstnu=0;
          for(int s2=0;s2<NS;s2++){ implicit_visc(g,dtv);
              double num=0,den=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); double ph=cos(kk*(i+0.5)*dx); num+=(g.mv[c]/g.r[c])*ph; den+=ph*ph; }
              double a1=num/den, f=a1/amp; amp=a1;
              worstnu=max(worstnu,fabs((1.0/f-1.0)/(lam*dtv)-nu)/nu); }
          double numz=0,den=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); double ph=cos(kk*(i+0.5)*dx); numz+=(g.mw[c]/g.r[c])*ph; den+=ph*ph; }
          double ampz=numz/den, expz=0.5*A0*pow(fex,NS);
          ampA=amp; exact=A0*exp(-nu*kk*kk*NS*dtv); zerr=fabs(ampz-expz)/expz;
          okA=(worstnu<1e-6 && fabs(amp-exact)/exact<0.10 && zerr<1e-6);
          printf("visc_diff A: nu_eff err(max)=%.2e; vy amp %.4e vs exp(-nu k2 t) %.4e (%+.1f%%); vz BE-factor err=%.2e\n",
              worstnu,ampA,exact,100*(ampA-exact)/exact,zerr); }
        bool okB;
        { int Nr=64,Nz=8; double dx=1.0,rho0=1.0,nu=1.0; ETA_AF=nu*rho0; AXISYM=true;
          Grid g(Nr,1,Nz,dx); double A0=1e-2,W0=3.0;
          for(int i=0;i<Nr;i++)for(int k=0;k<Nz;k++){int c=g.idx(i,0,k); g.r[c]=rho0; g.af[c]=1.0; double r=(i+0.5)*dx;
              g.mu[c]=rho0*A0*r; g.mw[c]=rho0*W0; g.E[c]=1e5+0.5*(g.mu[c]*g.mu[c]+g.mw[c]*g.mw[c])/rho0; }
          double dtv=10.0*dx*dx/nu;   // ~40x the 2D explicit limit; outer-wall Neumann defect decays e^(-di/3.2)
          implicit_visc(g,dtv);
          double eru=0,erw=0;
          for(int i=0;i<20;i++)for(int k=0;k<Nz;k++){int c=g.idx(i,0,k); double r=(i+0.5)*dx;
              eru=max(eru,fabs(g.mu[c]/g.r[c]-A0*r)/(A0*Nr*dx)); erw=max(erw,fabs(g.mw[c]/g.r[c]-W0)/W0); }
          okB=(eru<1e-5&&erw<1e-10);
          printf("visc_diff B (axisym): null-mode v_r=A*r residual=%.2e (tol 1e-5, i<20); v_z=const err=%.2e\n",eru,erw); }
        AXISYM=false; VISC_IMP=false; ETA_AF=0;
        printf("GATE (implicit viscous solve: exact BE diffusion at 20-40x explicit dt & cylindrical null mode): %s\n",(okA&&okB)?"PASS":"CHECK");
    } else if(mode=="visc_agree"){ // Option C gate 3.4c: explicit and implicit viscous paths agree at small eta (both stable), full step_rk2.
        // Fluidized basalt slab (af=1 -> Y=0: strength inert), transverse cosine mode: hydro is inert (no normal flow,
        // uniform P), so the two runs differ ONLY in the viscous scheme (explicit-in-RK2 vs backward Euler). At
        // nu*lam*dt ~ 6e-4 the schemes differ by O(x^2/2) per step -> fields must agree to ~1e-4, gate at 1e-3.
        MAT=Material::basalt(); double dx=500.0,CFL=0.4; int N=256; double rho0=MAT.rho0;
        double eta=1e8; int m=24; double L=N*dx,kk=m*M_PI/L; double A0=1.0,T=25.0;
        vector<double> vy1(N),vy2(N); double amp1=0,amp2=0;
        for(int run=0;run<2;run++){
            VISC_IMP=(run==1); ETA_AF=eta; TDEC=0;C_ACT=0;
            Grid g(N,1,1,dx);
            for(int i=0;i<N;i++){int c=g.idx(i,0,0); g.r[c]=rho0; g.af[c]=1.0; double vy=A0*cos(kk*(i+0.5)*dx);
                g.mv[c]=rho0*vy; g.E[c]=0.5*rho0*vy*vy; }
            double t=0; while(t<T){ double dt=CFL*dx/maxspeed(g); if(t+dt>T)dt=T-t; step_rk2(g,dt); t+=dt; }
            double num=0,den=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); double vv=g.mv[c]/g.r[c];
                (run?vy2:vy1)[i]=vv; double ph=cos(kk*(i+0.5)*dx); num+=vv*ph; den+=ph*ph; }
            if(run)amp2=num/den; else amp1=num/den; }
        VISC_IMP=false; ETA_AF=0;
        double dmax=0; for(int i=0;i<N;i++) dmax=max(dmax,fabs(vy1[i]-vy2[i]));
        double dec=1.0-amp1/A0;
        printf("visc_agree: explicit amp=%.6f implicit amp=%.6f (decay %.1f%%), max field diff=%.3e (A0=%g)\n",amp1,amp2,100*dec,dmax,A0);
        printf("GATE (explicit vs implicit at small eta: fields agree <1e-3*A0 & decay >1%%): %s\n",(dmax<1e-3*A0&&dec>0.01)?"PASS":"CHECK");
    } else if(mode=="visc_energy"){ // Option C gate 3.4e: viscous-dominated window -- total energy conserved, dissipation heats.
        // nu*dt/dx^2 ~ 30 (impossible explicitly): the mode's KE is wiped out in ~10 steps; ALL of it must reappear
        // as internal energy (the explicit path silently deleted it). Gate: |dE_total| < 1e-4*KE0 and IE gain
        // matches KE loss to 1e-4, with >50% of KE0 actually dissipated in the window.
        MAT=Material::basalt(); double dx=500.0,CFL=0.4; int N=256; double rho0=MAT.rho0;
        VISC_IMP=true; ETA_AF=5e11; TDEC=0;C_ACT=0;
        int m=24; double L=N*dx,kk=m*M_PI/L; double A0=2.0,T=5.0;
        Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){int c=g.idx(i,0,0); g.r[c]=rho0; g.af[c]=1.0; double vy=A0*cos(kk*(i+0.5)*dx);
            g.mv[c]=rho0*vy; g.E[c]=0.5*rho0*vy*vy; }
        double E0=0,KE0=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); E0+=g.E[c]; KE0+=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; }
        double t=0; while(t<T){ double dt=CFL*dx/maxspeed(g); if(t+dt>T)dt=T-t; step_rk2(g,dt); t+=dt; }
        double E1=0,KE1=0,IE1=0; for(int i=0;i<N;i++){int c=g.idx(i,0,0); E1+=g.E[c];
            double ke=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; KE1+=ke; IE1+=g.E[c]-ke; }
        VISC_IMP=false; ETA_AF=0;
        double drift=fabs(E1-E0)/KE0, keloss=(KE0-KE1)/KE0, heat=IE1/KE0;   // IE0=0 by construction
        printf("visc_energy: total E drift=%.2e of KE0; KE lost %.4f -> IE gained %.4f (of KE0); CG iters=%ld\n",drift,keloss,heat,VISC_IT);
        printf("GATE (viscous energy: |dE_tot|<1e-4*KE0 & heating==KE loss to 1e-4 & >50%% dissipated): %s\n",(drift<1e-4&&fabs(heat-keloss)<1e-4&&keloss>0.5)?"PASS":"CHECK");
    } else if(mode=="bingham_slope"){ // Option C gate 3.4d: Bingham arrest. A fluidized triangular HEAP on an intact flat
        // substrate (2D x-z). Driving basal shear tau ~ rho*g*h*|dh/dx|; the von Mises cap Y_B yields in shear at
        // Y_B/sqrt(3). Analytic arrest: a Bingham heap spreads until rho*g*h*|dh/dx| <= tau_y everywhere (parabolic
        // final flank h^2 = 2*(tau_y/rho g)*d).
        // Three-way rate contrast (robust to the two known numerical floors: free-surface cell churn ~0.5 m/s and
        // slow vacuum-face sloughing; NB the heap is NOT shallow (H/W~0.3) so 2D stress redistribution legitimately
        // arrests it a factor of a few above the naive infinite-slope envelope -- the envelope ratio is only bounded,
        // not pinned to 1):
        //   (A) tau_y = 2.5*tau_max -> STATIC: peak holds (>=92% of H), only the sloughing floor moves.
        //   (B) tau_y ~ 0.2*tau_max -> FLOWS (peak drops >25%) then SELF-ARRESTS: late-window peak-decay rate <15%
        //       of the early rate, af still 1 throughout -- TDEC=0, so arrest is by STRESS BALANCE alone (Option C
        //       acceptance signature), and the flank stays within a resolution-tolerant band of the yield envelope.
        //   (C) Y_B=0 (Newtonian control, same eta) -> must KEEP collapsing: proves B's arrest is the yield floor,
        //       not viscosity or numerics.
        MAT=Material::basalt(); GZ=3.71; ROCK=true; RHO_CFL=100.0; RHO_VAC=100.0; VISC_IMP=true; ETA_AF=5e7; TDEC=0;C_ACT=0;
        VOID_AMB=0.27;   // void resets go to AMBIENT -- a lithostatic reference would REFILL evacuated cells as the heap subsides (the 2026-07-02 voidzero finding)
        double rho0=MAT.rho0,Ab=MAT.A; int NX=112,NZc=20; double dx=250.0,CFL=0.4;
        double z0=1500.0, H=2000.0, W=6000.0, xc=0.5*NX*dx;
        double taumax=rho0*GZ*H*(H/W);   // ~ rho g h dh/dx at the peak flank
        double res[3][4]; // per case: hmax/H at end, early peak-decay (m), late peak-decay (m), flank stress ratio
        vector<double> hsamp;
        for(int cs=0;cs<3;cs++){
            Y_BINGHAM=(cs==0? 2.5*sqrt(3.0)*taumax : (cs==1? 0.2*sqrt(3.0)*taumax : 0.0));
            double tauy=(Y_BINGHAM>0?Y_BINGHAM/sqrt(3.0):1e30);
            Grid g(NX,1,NZc,dx);
            auto hs0=[&](int i){ double xx=fabs((i+0.5)*dx-xc); return xx<W? H*(1.0-xx/W):0.0; };
            for(int i=0;i<NX;i++)for(int k=0;k<NZc;k++){ double z=(k+0.5)*dx; int c=g.idx(i,0,k);
                double zs=z0+hs0(i);
                if(z<zs){ g.r[c]=rho0*(1.0+rho0*GZ*(zs-z)/Ab); } else g.r[c]=0.27; g.E[c]=0;
                g.af[c]=(z>z0)?1.0:0.0; }   // heap + everything above the substrate plane is fluidized; substrate intact
            set_ref(g);
            auto surfM=[&](int i){ double m2=0; for(int k=0;k<NZc;k++){ double r=g.r[g.idx(i,0,k)]; if(r>100.0) m2+=r; } return m2*dx/rho0; };   // column-mass surface (threshold-robust)
            vector<double> h0(NX); for(int i=0;i<NX;i++) h0[i]=surfM(i);
            double Tc=(cs==0?400.0:(cs==1?1200.0:600.0)), t=0, nlog=100.0;
            hsamp.assign(1,H);   // hmax samples every 100 s
            while(t<Tc){ double dt=CFL*dx/maxspeed(g); if(t+dt>Tc)dt=Tc-t; step_rk2(g,dt);
                for(int i=0;i<NX;i++)for(int k=0;k<2;k++){int c=g.idx(i,0,k); g.r[c]=REF_R0[c]; g.mu[c]=g.mv[c]=g.mw[c]=0; g.E[c]=0;}  // pinned bedrock floor
                t+=dt;
                if(t>=nlog){ double hm=0; for(int i=2;i<NX-2;i++) hm=max(hm,surfM(i)-z0); hsamp.push_back(hm); nlog+=100.0; } }
            double hmax=0,ratio=0;
            for(int i=2;i<NX-2;i++){ double h=surfM(i)-z0; hmax=max(hmax,h);
                double hp=surfM(i+1)-z0, hm2=surfM(i-1)-z0;
                if(h>500.0) ratio=max(ratio, rho0*GZ*h*fabs(hp-hm2)/(2*dx) / tauy); }   // steepest-flank stress vs yield envelope
            int ns=hsamp.size();
            double early=hsamp[0]-hsamp[min(4,ns-1)];                    // peak decay over the first 400 s
            double late =hsamp[max(0,ns-5)]-hsamp[ns-1];                 // ... and the last 400 s
            res[cs][0]=hmax/H; res[cs][1]=early; res[cs][2]=late; res[cs][3]=ratio;
            printf("bingham_slope %s: Y_B=%.3e  hmax/H=%.3f  peak decay early/late=%.0f/%.0f m  flank tau/tau_y=%.2f\n",
                cs==0?"A(static) ":(cs==1?"B(arrest) ":"C(control)"),Y_BINGHAM,res[cs][0],early,late,ratio); }
        VISC_IMP=false; ETA_AF=0; Y_BINGHAM=0; ROCK=false; GZ=0; VOID_AMB=0;
        bool A=(res[0][0]>0.92);                                            // sub-yield: peak holds (sloughing floor only)
        bool B=(res[1][0]<0.75 && res[1][2]<0.15*res[1][1] && res[1][3]<8.0);// flowed, then late rate <15% of early = ARREST, flank bounded
        bool C=(res[2][2]>3.0*res[1][2] && res[2][0]<res[1][0]);             // Newtonian control keeps collapsing
        printf("GATE (Bingham arrest: static below yield; flow-then-self-arrest above it; Y_B=0 control keeps flowing): %s\n",(A&&B&&C)?"PASS":"CHECK");
    } else {
        fprintf(stderr,"unknown mode '%s' (modes: sod sedov shear yield vacuum bshock freefall atmos tensile alimpact pierazzo2d prater substrate collapse surface tracer af_activate sedov_axi lame af_visc vib_advect friction crater visc_diff visc_agree visc_energy bingham_slope)\n",mode.c_str());
        return 2;   // fatal: never silently validate the wrong physics
    }
    return 0;
}
