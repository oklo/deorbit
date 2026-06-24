// euler M1a: CPU reference -- 3D finite-volume compressible hydro (ideal gas),
// dimensionally split, HLLC Riemann flux. First-order (validate the flux/update
// vs the EXACT Sod solution); MUSCL reconstruction added next.
//   ./hydro_cpu sod    -> Sod shock tube, prints L1 error vs exact Riemann solution
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
static double GAM = 1.4;

struct V5 { double r, mu, mv, mw, E; };               // conserved (mu=rho*u etc.)
static inline double pres(const V5&c){ double ke=0.5*(c.mu*c.mu+c.mv*c.mv+c.mw*c.mw)/c.r; return (GAM-1.0)*(c.E-ke); }

// HLLC flux in direction d (0=x,1=y,2=z): returns conserved flux. uses primitive L/R.
// prim: r, vn (normal vel), vt1, vt2, p.
static V5 hllc(double rL,double vnL,double t1L,double t2L,double pL,
               double rR,double vnR,double t1R,double t2R,double pR){
    double cL=sqrt(GAM*pL/rL), cR=sqrt(GAM*pR/rR);
    double SL=min(vnL-cL, vnR-cR), SR=max(vnL+cL, vnR+cR);
    // contact speed (Toro 10.37)
    double Sstar=(pR-pL + rL*vnL*(SL-vnL) - rR*vnR*(SR-vnR)) / (rL*(SL-vnL) - rR*(SR-vnR));
    auto F=[&](double r,double vn,double t1,double t2,double p){   // physical flux (normal dir)
        double E=p/(GAM-1.0)+0.5*r*(vn*vn+t1*t1+t2*t2);
        V5 f; f.r=r*vn; f.mu=r*vn*vn+p; f.mv=r*vn*t1; f.mw=r*vn*t2; f.E=(E+p)*vn; return f;
    };
    auto cons=[&](double r,double vn,double t1,double t2,double p){
        V5 c; c.r=r; c.mu=r*vn; c.mv=r*t1; c.mw=r*t2; c.E=p/(GAM-1.0)+0.5*r*(vn*vn+t1*t1+t2*t2); return c; };
    if(SL>=0) return F(rL,vnL,t1L,t2L,pL);
    if(SR<=0) return F(rR,vnR,t1R,t2R,pR);
    // star states (Toro 10.39)
    auto star=[&](double r,double vn,double t1,double t2,double p,double S){
        double f=r*(S-vn)/(S-Sstar); V5 u;
        u.r=f; u.mu=f*Sstar; u.mv=f*t1; u.mw=f*t2;
        double E=p/(GAM-1.0)+0.5*r*(vn*vn+t1*t1+t2*t2);
        u.E=f*(E/r + (Sstar-vn)*(Sstar + p/(r*(S-vn)))); return u; };
    if(Sstar>=0){ V5 FL=F(rL,vnL,t1L,t2L,pL), UsL=star(rL,vnL,t1L,t2L,pL,SL), UL=cons(rL,vnL,t1L,t2L,pL);
        V5 r; r.r=FL.r+SL*(UsL.r-UL.r); r.mu=FL.mu+SL*(UsL.mu-UL.mu); r.mv=FL.mv+SL*(UsL.mv-UL.mv);
        r.mw=FL.mw+SL*(UsL.mw-UL.mw); r.E=FL.E+SL*(UsL.E-UL.E); return r; }
    V5 FR=F(rR,vnR,t1R,t2R,pR), UsR=star(rR,vnR,t1R,t2R,pR,SR), UR=cons(rR,vnR,t1R,t2R,pR);
    V5 r; r.r=FR.r+SR*(UsR.r-UR.r); r.mu=FR.mu+SR*(UsR.mu-UR.mu); r.mv=FR.mv+SR*(UsR.mv-UR.mv);
    r.mw=FR.mw+SR*(UsR.mw-UR.mw); r.E=FR.E+SR*(UsR.E-UR.E); return r;
}

// ---- exact Sod (Toro exact Riemann solver, ideal gas), sample at x/t = S ----
static void exact_sod(double xs, double t, double&r,double&u,double&p){
    double rL=1,uL=0,pL=1, rR=0.125,uR=0,pR=0.1;
    double cL=sqrt(GAM*pL/rL), cR=sqrt(GAM*pR/rR);
    double g=GAM, G1=(g-1)/(2*g), G2=(g+1)/(2*g), G3=2*g/(g-1), G4=2/(g-1), G5=2/(g+1), G6=(g-1)/(g+1);
    auto fk=[&](double P,double rk,double pk,double ck){
        if(P>pk){ double A=G5/rk,B=G6*pk; return (P-pk)*sqrt(A/(P+B)); }
        else return G4*ck*(pow(P/pk,G1)-1); };
    auto fkp=[&](double P,double rk,double pk,double ck){
        if(P>pk){ double A=G5/rk,B=G6*pk; return sqrt(A/(B+P))*(1-(P-pk)/(2*(B+P))); }
        else return 1.0/(rk*ck)*pow(P/pk,-G2); };
    double P=0.5*(pL+pR);
    for(int it=0;it<100;it++){ double f=fk(P,rL,pL,cL)+fk(P,rR,pR,cR)+(uR-uL);
        double fp=fkp(P,rL,pL,cL)+fkp(P,rR,pR,cR); double Pn=P-f/fp; if(fabs(Pn-P)<1e-12){P=Pn;break;} P=fabs(Pn);}
    double us=0.5*(uL+uR)+0.5*(fk(P,rR,pR,cR)-fk(P,rL,pL,cL));
    double S=xs/t;
    if(S<=us){ // left of contact
        if(P>pL){ double SL=uL-cL*sqrt(G2*P/pL+G1); if(S<SL){r=rL;u=uL;p=pL;} else {r=rL*((P/pL+G6)/(G6*P/pL+1));u=us;p=P;} }
        else { double SHL=uL-cL, STL=us-cL*pow(P/pL,G1);
            if(S<SHL){r=rL;u=uL;p=pL;} else if(S>STL){r=rL*pow(P/pL,1/g);u=us;p=P;}
            else { u=G5*(cL+(g-1)/2*uL+S); double c=G5*(cL+(g-1)/2*(uL-S)); r=rL*pow(c/cL,G4); p=pL*pow(c/cL,G3);} }
    } else { // right of contact
        if(P>pR){ double SR=uR+cR*sqrt(G2*P/pR+G1); if(S>SR){r=rR;u=uR;p=pR;} else {r=rR*((P/pR+G6)/(G6*P/pR+1));u=us;p=P;} }
        else { double SHR=uR+cR, STR=us+cR*pow(P/pR,G1);
            if(S>SHR){r=rR;u=uR;p=pR;} else if(S<STR){r=rR*pow(P/pR,1/g);u=us;p=P;}
            else { u=G5*(-cR+(g-1)/2*uR+S); double c=G5*(cR-(g-1)/2*(uR-S)); r=rR*pow(c/cR,G4); p=pR*pow(c/cR,G3);} }
    }
}

int main(int argc,char**argv){
    string mode=argc>1?argv[1]:"sod";
    int N=200; double L=1.0, dx=L/N, tend=0.2, CFL=0.4;
    vector<V5> U(N), Un(N);
    for(int i=0;i<N;i++){ double x=(i+0.5)*dx; double r=x<0.5?1.0:0.125, p=x<0.5?1.0:0.1;
        U[i]={r,0,0,0,p/(GAM-1.0)}; }
    double t=0; int step=0;
    while(t<tend){
        double smax=1e-30; for(auto&c:U){ double p=pres(c), cs=sqrt(GAM*p/c.r); smax=max(smax, fabs(c.mu/c.r)+cs); }
        double dt=CFL*dx/smax; if(t+dt>tend)dt=tend-t;
        // x-sweep, first order, transmissive BCs
        vector<V5> Fl(N+1);
        for(int f=0;f<=N;f++){ int iL=max(0,f-1), iR=min(N-1,f);
            V5 a=U[iL], b=U[iR];
            Fl[f]=hllc(a.r,a.mu/a.r,a.mv/a.r,a.mw/a.r,pres(a), b.r,b.mu/b.r,b.mv/b.r,b.mw/b.r,pres(b)); }
        for(int i=0;i<N;i++){ double k=dt/dx;
            Un[i].r =U[i].r -k*(Fl[i+1].r -Fl[i].r);
            Un[i].mu=U[i].mu-k*(Fl[i+1].mu-Fl[i].mu);
            Un[i].mv=U[i].mv-k*(Fl[i+1].mv-Fl[i].mv);
            Un[i].mw=U[i].mw-k*(Fl[i+1].mw-Fl[i].mw);
            Un[i].E =U[i].E -k*(Fl[i+1].E -Fl[i].E); }
        U=Un; t+=dt; step++;
    }
    // compare to exact
    double l1r=0,l1p=0,l1u=0;
    for(int i=0;i<N;i++){ double x=(i+0.5)*dx; double re,ue,pe; exact_sod(x-0.5,tend,re,ue,pe);
        l1r+=fabs(U[i].r-re); l1p+=fabs(pres(U[i])-pe); l1u+=fabs(U[i].mu/U[i].r-ue); }
    l1r/=N; l1p/=N; l1u/=N;
    printf("Sod N=%d steps=%d t=%.3f  L1 errors vs EXACT: rho=%.4f  p=%.4f  u=%.4f\n",N,step,t,l1r,l1p,l1u);
    printf("GATE (1st-order HLLC, expect L1<~0.02): %s\n", (l1r<0.02&&l1p<0.02&&l1u<0.02)?"PASS":"check");
    return 0;
}
