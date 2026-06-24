// euler M1: CPU reference -- 3D finite-volume compressible hydro (ideal gas).
// Unsplit MUSCL (minmod) + HLLC + SSP-RK2 (2nd order space & time).
//   ./hydro_cpu sod    -> Sod tube vs EXACT Riemann (L1 errors)
//   ./hydro_cpu sedov  -> 3D point blast vs analytic Sedov-Taylor (shock radius)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
static double GAM = 1.4;

// ---------- HLLC flux (normal vn, transverse t1,t2); returns (f_r,f_mn,f_mt1,f_mt2,f_E) ----------
struct F5 { double r, mn, mt1, mt2, E; };
static F5 hllc(double rL,double vnL,double t1L,double t2L,double pL,
               double rR,double vnR,double t1R,double t2R,double pR){
    double cL=sqrt(GAM*pL/rL), cR=sqrt(GAM*pR/rR);
    double SL=min(vnL-cL,vnR-cR), SR=max(vnL+cL,vnR+cR);
    double Ss=(pR-pL+rL*vnL*(SL-vnL)-rR*vnR*(SR-vnR))/(rL*(SL-vnL)-rR*(SR-vnR));
    auto F=[&](double r,double vn,double t1,double t2,double p){
        double E=p/(GAM-1.0)+0.5*r*(vn*vn+t1*t1+t2*t2);
        return F5{r*vn, r*vn*vn+p, r*vn*t1, r*vn*t2, (E+p)*vn}; };
    if(SL>=0) return F(rL,vnL,t1L,t2L,pL);
    if(SR<=0) return F(rR,vnR,t1R,t2R,pR);
    auto star=[&](double r,double vn,double t1,double t2,double p,double S,F5 Fk){
        double f=r*(S-vn)/(S-Ss); double E=p/(GAM-1.0)+0.5*r*(vn*vn+t1*t1+t2*t2);
        F5 Us{f, f*Ss, f*t1, f*t2, f*(E/r+(Ss-vn)*(Ss+p/(r*(S-vn))))};
        F5 U{r, r*vn, r*t1, r*t2, E};
        return F5{Fk.r+S*(Us.r-U.r), Fk.mn+S*(Us.mn-U.mn), Fk.mt1+S*(Us.mt1-U.mt1),
                  Fk.mt2+S*(Us.mt2-U.mt2), Fk.E+S*(Us.E-U.E)}; };
    if(Ss>=0) return star(rL,vnL,t1L,t2L,pL,SL,F(rL,vnL,t1L,t2L,pL));
    return star(rR,vnR,t1R,t2R,pR,SR,F(rR,vnR,t1R,t2R,pR));
}
static inline double minmod(double a,double b){ return (a*b<=0)?0.0:(fabs(a)<fabs(b)?a:b); }

// ---------- 3D solver ----------
struct Grid {
    int nx,ny,nz; double dx;
    vector<double> r,mu,mv,mw,E;
    int idx(int i,int j,int k)const{ return (i*ny+j)*nz+k; }
    Grid(int X,int Y,int Z,double d):nx(X),ny(Y),nz(Z),dx(d),
        r(X*Y*Z),mu(X*Y*Z),mv(X*Y*Z),mw(X*Y*Z),E(X*Y*Z){}
};
static inline double P(const Grid&g,int c){ double ke=0.5*(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; return (GAM-1.0)*(g.E[c]-ke); }

// L(U) -> dU/dt (5 arrays), unsplit MUSCL+HLLC, transmissive BC by index clamp.
static void Lop(const Grid&g, vector<double>&dr,vector<double>&dmu,vector<double>&dmv,vector<double>&dmw,vector<double>&dE){
    fill(dr.begin(),dr.end(),0.0);fill(dmu.begin(),dmu.end(),0.0);
    fill(dmv.begin(),dmv.end(),0.0);fill(dmw.begin(),dmw.end(),0.0);fill(dE.begin(),dE.end(),0.0);
    double invdx=1.0/g.dx;
    auto prim=[&](int c,double*W){ W[0]=g.r[c]; W[1]=g.mu[c]/W[0]; W[2]=g.mv[c]/W[0]; W[3]=g.mw[c]/W[0]; W[4]=P(g,c); };
    // per-direction, per-line sweep; faces 0..nL (incl. transmissive boundaries via clamp)
    auto sweep=[&](int d){
        int nL=(d==0?g.nx:(d==1?g.ny:g.nz));
        int nA=(d==0?g.ny:g.nx), nB=(d==2?g.ny:g.nz);
        auto CELL=[&](int p,int a,int b){ return d==0?g.idx(p,a,b):(d==1?g.idx(a,p,b):g.idx(a,b,p)); };
        int in=(d==0?1:(d==1?2:3)), it1=(d==0?2:1), it2=(d==2?2:3);
        vector<double> Fr(nL+1),Fmn(nL+1),Ft1(nL+1),Ft2(nL+1),FE(nL+1);
        for(int a=0;a<nA;a++)for(int b=0;b<nB;b++){
            for(int f=0;f<=nL;f++){
                int iL=max(0,f-1),iR=min(nL-1,f),iLL=max(0,f-2),iRR=min(nL-1,f+1);
                double A[5],B[5],AA[5],BB[5],L[5],R[5];
                prim(CELL(iL,a,b),A); prim(CELL(iR,a,b),B); prim(CELL(iLL,a,b),AA); prim(CELL(iRR,a,b),BB);
                for(int q=0;q<5;q++){ L[q]=A[q]+0.5*minmod(A[q]-AA[q],B[q]-A[q]); R[q]=B[q]-0.5*minmod(B[q]-A[q],BB[q]-B[q]); }
                if(L[0]<1e-9)L[0]=A[0]; if(R[0]<1e-9)R[0]=B[0]; if(L[4]<1e-9)L[4]=A[4]; if(R[4]<1e-9)R[4]=B[4];
                F5 fl=hllc(L[0],L[in],L[it1],L[it2],L[4], R[0],R[in],R[it1],R[it2],R[4]);
                Fr[f]=fl.r;Fmn[f]=fl.mn;Ft1[f]=fl.mt1;Ft2[f]=fl.mt2;FE[f]=fl.E;
            }
            for(int p=0;p<nL;p++){ int c=CELL(p,a,b);
                dr[c]+=-(Fr[p+1]-Fr[p])*invdx; dE[c]+=-(FE[p+1]-FE[p])*invdx;
                double dmn=-(Fmn[p+1]-Fmn[p])*invdx, dt1=-(Ft1[p+1]-Ft1[p])*invdx, dt2=-(Ft2[p+1]-Ft2[p])*invdx;
                if(d==0){dmu[c]+=dmn;dmv[c]+=dt1;dmw[c]+=dt2;}
                else if(d==1){dmv[c]+=dmn;dmu[c]+=dt1;dmw[c]+=dt2;}
                else {dmw[c]+=dmn;dmu[c]+=dt1;dmv[c]+=dt2;}
            }
        }
    };
    sweep(0); sweep(1); sweep(2);
}
static double maxspeed(const Grid&g){ double s=1e-30; int n=g.r.size();
    for(int c=0;c<n;c++){ double p=P(g,c); double cs=sqrt(GAM*max(p,1e-12)/g.r[c]);
        double v=sqrt(g.mu[c]*g.mu[c]+g.mv[c]*g.mv[c]+g.mw[c]*g.mw[c])/g.r[c]; s=max(s,v+cs);} return s; }

static void step_rk2(Grid&g,double dt){
    int n=g.r.size();
    vector<double> dr(n),dmu(n),dmv(n),dmw(n),dE(n);
    Grid g1=g;
    Lop(g,dr,dmu,dmv,dmw,dE);
    for(int c=0;c<n;c++){ g1.r[c]=g.r[c]+dt*dr[c]; g1.mu[c]=g.mu[c]+dt*dmu[c]; g1.mv[c]=g.mv[c]+dt*dmv[c];
        g1.mw[c]=g.mw[c]+dt*dmw[c]; g1.E[c]=g.E[c]+dt*dE[c]; }
    Lop(g1,dr,dmu,dmv,dmw,dE);
    for(int c=0;c<n;c++){ g.r[c]=0.5*(g.r[c]+g1.r[c]+dt*dr[c]); g.mu[c]=0.5*(g.mu[c]+g1.mu[c]+dt*dmu[c]);
        g.mv[c]=0.5*(g.mv[c]+g1.mv[c]+dt*dmv[c]); g.mw[c]=0.5*(g.mw[c]+g1.mw[c]+dt*dmw[c]); g.E[c]=0.5*(g.E[c]+g1.E[c]+dt*dE[c]); }
}

// exact Sod (Toro), sample at x/t
static void exact_sod(double S,double&r,double&u,double&p){
    double rL=1,uL=0,pL=1, rR=0.125,uR=0,pR=0.1, g=GAM;
    double cL=sqrt(g*pL/rL), cR=sqrt(g*pR/rR);
    double G1=(g-1)/(2*g),G3=2*g/(g-1),G4=2/(g-1),G5=2/(g+1),G6=(g-1)/(g+1),G2=(g+1)/(2*g);
    auto fk=[&](double Pp,double rk,double pk,double ck){ if(Pp>pk){double A=G5/rk,B=G6*pk;return(Pp-pk)*sqrt(A/(Pp+B));}return G4*ck*(pow(Pp/pk,G1)-1);};
    auto fkp=[&](double Pp,double rk,double pk,double ck){ if(Pp>pk){double A=G5/rk,B=G6*pk;return sqrt(A/(B+Pp))*(1-(Pp-pk)/(2*(B+Pp)));}return pow(Pp/pk,-G2)/(rk*ck);};
    double Pp=0.5*(pL+pR); for(int it=0;it<100;it++){double f=fk(Pp,rL,pL,cL)+fk(Pp,rR,pR,cR)+(uR-uL);double fp=fkp(Pp,rL,pL,cL)+fkp(Pp,rR,pR,cR);double Pn=fabs(Pp-f/fp);if(fabs(Pn-Pp)<1e-12){Pp=Pn;break;}Pp=Pn;}
    double us=0.5*(uL+uR)+0.5*(fk(Pp,rR,pR,cR)-fk(Pp,rL,pL,cL));
    if(S<=us){ if(Pp>pL){double SL=uL-cL*sqrt(G2*Pp/pL+G1);if(S<SL){r=rL;u=uL;p=pL;}else{r=rL*((Pp/pL+G6)/(G6*Pp/pL+1));u=us;p=Pp;}}
        else{double SHL=uL-cL,STL=us-cL*pow(Pp/pL,G1);if(S<SHL){r=rL;u=uL;p=pL;}else if(S>STL){r=rL*pow(Pp/pL,1/g);u=us;p=Pp;}else{u=G5*(cL+(g-1)/2*uL+S);double c=G5*(cL+(g-1)/2*(uL-S));r=rL*pow(c/cL,G4);p=pL*pow(c/cL,G3);}}}
    else{ if(Pp>pR){double SR=uR+cR*sqrt(G2*Pp/pR+G1);if(S>SR){r=rR;u=uR;p=pR;}else{r=rR*((Pp/pR+G6)/(G6*Pp/pR+1));u=us;p=Pp;}}
        else{double SHR=uR+cR,STR=us+cR*pow(Pp/pR,G1);if(S>SHR){r=rR;u=uR;p=pR;}else if(S<STR){r=rR*pow(Pp/pR,1/g);u=us;p=Pp;}else{u=G5*(-cR+(g-1)/2*uR+S);double c=G5*(cR-(g-1)/2*(uR-S));r=rR*pow(c/cR,G4);p=pR*pow(c/cR,G3);}}}
}

int main(int argc,char**argv){
    string mode=argc>1?argv[1]:"sod";
    if(mode=="sod"){
        int N=200; double L=1.0,dx=L/N,tend=0.2,CFL=0.4;
        Grid g(N,1,1,dx);
        for(int i=0;i<N;i++){double x=(i+0.5)*dx;double r=x<0.5?1:0.125,p=x<0.5?1:0.1;int c=g.idx(i,0,0);g.r[c]=r;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=p/(GAM-1);}
        double t=0;int s=0; while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        double l1r=0,l1p=0,l1u=0; for(int i=0;i<N;i++){double x=(i+0.5)*dx;double re,ue,pe;exact_sod((x-0.5)/tend,re,ue,pe);int c=g.idx(i,0,0);l1r+=fabs(g.r[c]-re);l1p+=fabs(P(g,c)-pe);l1u+=fabs(g.mu[c]/g.r[c]-ue);}
        l1r/=N;l1p/=N;l1u/=N;
        printf("Sod MUSCL N=%d steps=%d  L1 vs EXACT: rho=%.4f p=%.4f u=%.4f\n",N,s,l1r,l1p,l1u);
        printf("GATE (2nd-order, expect L1<~0.006): %s\n",(l1r<0.006&&l1p<0.006&&l1u<0.008)?"PASS":"CHECK");
    } else { // sedov
        int N=64; double L=2.0,dx=L/N,CFL=0.3,tend=0.5; double rho0=1.0,E0=1.0,p0=1e-4;
        Grid g(N,N,N,dx);
        for(int c=0;c<N*N*N;c++){g.r[c]=rho0;g.mu[c]=g.mv[c]=g.mw[c]=0;g.E[c]=p0/(GAM-1);}
        // deposit E0 in the central cell
        int ic=N/2; int cc=g.idx(ic,ic,ic); double Vc=dx*dx*dx; g.E[cc]+=E0/Vc;
        double t=0;int s=0; while(t<tend){double dt=CFL*dx/maxspeed(g);if(t+dt>tend)dt=tend-t;step_rk2(g,dt);t+=dt;s++;}
        // numerical shock radius = radius of peak density; peak compression
        double rmax=0,Rsh=0; double xc=(ic+0.5)*dx;
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++){int c=g.idx(i,j,k);
            if(g.r[c]>rmax){rmax=g.r[c]; double X=(i+0.5)*dx-xc,Y=(j+0.5)*dx-xc,Z=(k+0.5)*dx-xc; Rsh=sqrt(X*X+Y*Y+Z*Z);}}
        double alpha=0.851; double Ranal=pow(E0*t*t/(alpha*rho0),0.2);  // Sedov-Taylor, gamma=1.4, 3D
        double comp=(GAM+1)/(GAM-1);
        printf("Sedov N=%d steps=%d t=%.3f  R_num=%.4f  R_analytic=%.4f  (err %.1f%%)\n",N,s,t,Rsh,Ranal,100*fabs(Rsh-Ranal)/Ranal);
        printf("  peak compression num=%.2f  strong-shock limit=%.2f\n",rmax/rho0,comp);
        printf("GATE (R within ~10%%): %s\n",(fabs(Rsh-Ranal)/Ranal<0.10)?"PASS":"CHECK");
    }
    return 0;
}
