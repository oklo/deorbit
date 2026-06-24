#include <metal_stdlib>
using namespace metal;
// euler M1b: GPU hydro kernels (ideal gas). FP32. Mirrors hydro_cpu.cpp (the oracle).
// struct holds 5 conserved (r,mu,mv,mw,E) OR a flux (r,mn,mt1,mt2,E).
struct C5 { float r, mu, mv, mw, E; };
inline float presC(C5 c,float gam){ float ke=0.5f*(c.mu*c.mu+c.mv*c.mv+c.mw*c.mw)/c.r; return (gam-1.0f)*(c.E-ke); }
inline float mmod(float a,float b){ return (a*b<=0.0f)?0.0f:(fabs(a)<fabs(b)?a:b); }

// HLLC: prim L/R (normal vn, transverse t1,t2) -> flux (r,mn,mt1,mt2,E) packed in C5
inline C5 hllc(float rL,float vnL,float t1L,float t2L,float pL,
               float rR,float vnR,float t1R,float t2R,float pR,float gam){
    float cL=sqrt(gam*pL/rL), cR=sqrt(gam*pR/rR);
    float SL=min(vnL-cL,vnR-cR), SR=max(vnL+cL,vnR+cR);
    float Ss=(pR-pL+rL*vnL*(SL-vnL)-rR*vnR*(SR-vnR))/(rL*(SL-vnL)-rR*(SR-vnR));
    float EL=pL/(gam-1.0f)+0.5f*rL*(vnL*vnL+t1L*t1L+t2L*t2L);
    float ER=pR/(gam-1.0f)+0.5f*rR*(vnR*vnR+t1R*t1R+t2R*t2R);
    C5 FL={rL*vnL, rL*vnL*vnL+pL, rL*vnL*t1L, rL*vnL*t2L, (EL+pL)*vnL};
    C5 FR={rR*vnR, rR*vnR*vnR+pR, rR*vnR*t1R, rR*vnR*t2R, (ER+pR)*vnR};
    if(SL>=0.0f) return FL;
    if(SR<=0.0f) return FR;
    if(Ss>=0.0f){ float f=rL*(SL-vnL)/(SL-Ss);
        C5 U={rL,rL*vnL,rL*t1L,rL*t2L,EL};
        C5 Us={f,f*Ss,f*t1L,f*t2L,f*(EL/rL+(Ss-vnL)*(Ss+pL/(rL*(SL-vnL))))};
        return (C5){FL.r+SL*(Us.r-U.r),FL.mu+SL*(Us.mu-U.mu),FL.mv+SL*(Us.mv-U.mv),FL.mw+SL*(Us.mw-U.mw),FL.E+SL*(Us.E-U.E)}; }
    float f=rR*(SR-vnR)/(SR-Ss);
    C5 U={rR,rR*vnR,rR*t1R,rR*t2R,ER};
    C5 Us={f,f*Ss,f*t1R,f*t2R,f*(ER/rR+(Ss-vnR)*(Ss+pR/(rR*(SR-vnR))))};
    return (C5){FR.r+SR*(Us.r-U.r),FR.mu+SR*(Us.mu-U.mu),FR.mv+SR*(Us.mv-U.mv),FR.mw+SR*(Us.mw-U.mw),FR.E+SR*(Us.E-U.E)};
}
// flux at face between conserved cells b(=cL) and cc(=cR); a=cLL, d=cRR. returns conserved flux.
inline C5 faceflux(C5 a,C5 b,C5 cc,C5 d,int dir,float gam){
    float Pa[5]={a.r,a.mu/a.r,a.mv/a.r,a.mw/a.r,presC(a,gam)};
    float Pb[5]={b.r,b.mu/b.r,b.mv/b.r,b.mw/b.r,presC(b,gam)};
    float Pc[5]={cc.r,cc.mu/cc.r,cc.mv/cc.r,cc.mw/cc.r,presC(cc,gam)};
    float Pd[5]={d.r,d.mu/d.r,d.mv/d.r,d.mw/d.r,presC(d,gam)};
    float L[5],R[5];
    for(int q=0;q<5;q++){ L[q]=Pb[q]+0.5f*mmod(Pb[q]-Pa[q],Pc[q]-Pb[q]); R[q]=Pc[q]-0.5f*mmod(Pc[q]-Pb[q],Pd[q]-Pc[q]); }
    if(L[0]<1e-9f)L[0]=Pb[0]; if(R[0]<1e-9f)R[0]=Pc[0]; if(L[4]<1e-9f)L[4]=Pb[4]; if(R[4]<1e-9f)R[4]=Pc[4];
    int in=(dir==0?1:(dir==1?2:3)), it1=(dir==0?2:1), it2=(dir==2?2:3);
    C5 fl=hllc(L[0],L[in],L[it1],L[it2],L[4], R[0],R[in],R[it1],R[it2],R[4],gam);
    C5 f; f.r=fl.r; f.E=fl.E;
    if(dir==0){f.mu=fl.mu;f.mv=fl.mv;f.mw=fl.mw;}
    else if(dir==1){f.mv=fl.mu;f.mu=fl.mv;f.mw=fl.mw;}
    else {f.mw=fl.mu;f.mu=fl.mv;f.mv=fl.mw;}
    return f;
}
inline C5 rd(device const float*r,device const float*mu,device const float*mv,device const float*mw,device const float*E,int c){
    return (C5){r[c],mu[c],mv[c],mw[c],E[c]}; }

kernel void lop(device const float*r[[buffer(0)]],device const float*mu[[buffer(1)]],device const float*mv[[buffer(2)]],
                device const float*mw[[buffer(3)]],device const float*E[[buffer(4)]],
                device float*dr[[buffer(5)]],device float*dmu[[buffer(6)]],device float*dmv[[buffer(7)]],
                device float*dmw[[buffer(8)]],device float*dE[[buffer(9)]],
                constant int&nx[[buffer(10)]],constant int&ny[[buffer(11)]],constant int&nz[[buffer(12)]],
                constant float&invdx[[buffer(13)]],constant float&gam[[buffer(14)]],constant uint&n[[buffer(15)]],
                uint gid[[thread_position_in_grid]]){
    if(gid>=n) return;
    int k=gid%nz, j=(gid/nz)%ny, i=gid/(ny*nz);
    C5 acc={0,0,0,0,0};
    for(int dir=0;dir<3;dir++){
        int m2,m1,p0,p1,p2;
        if(dir==0){int a=max(0,i-2),b=max(0,i-1),d=min(nx-1,i+1),e=min(nx-1,i+2);
            m2=(a*ny+j)*nz+k;m1=(b*ny+j)*nz+k;p0=(i*ny+j)*nz+k;p1=(d*ny+j)*nz+k;p2=(e*ny+j)*nz+k;}
        else if(dir==1){int a=max(0,j-2),b=max(0,j-1),d=min(ny-1,j+1),e=min(ny-1,j+2);
            m2=(i*ny+a)*nz+k;m1=(i*ny+b)*nz+k;p0=(i*ny+j)*nz+k;p1=(i*ny+d)*nz+k;p2=(i*ny+e)*nz+k;}
        else{int a=max(0,k-2),b=max(0,k-1),d=min(nz-1,k+1),e=min(nz-1,k+2);
            m2=(i*ny+j)*nz+a;m1=(i*ny+j)*nz+b;p0=(i*ny+j)*nz+k;p1=(i*ny+j)*nz+d;p2=(i*ny+j)*nz+e;}
        C5 cm2=rd(r,mu,mv,mw,E,m2),cm1=rd(r,mu,mv,mw,E,m1),c0=rd(r,mu,mv,mw,E,p0),cp1=rd(r,mu,mv,mw,E,p1),cp2=rd(r,mu,mv,mw,E,p2);
        C5 FR=faceflux(cm1,c0,cp1,cp2,dir,gam), FL=faceflux(cm2,cm1,c0,cp1,dir,gam);
        acc.r-=(FR.r-FL.r); acc.mu-=(FR.mu-FL.mu); acc.mv-=(FR.mv-FL.mv); acc.mw-=(FR.mw-FL.mw); acc.E-=(FR.E-FL.E);
    }
    dr[gid]=acc.r*invdx; dmu[gid]=acc.mu*invdx; dmv[gid]=acc.mv*invdx; dmw[gid]=acc.mw*invdx; dE[gid]=acc.E*invdx;
}
kernel void wavespeed(device const float*r[[buffer(0)]],device const float*mu[[buffer(1)]],device const float*mv[[buffer(2)]],
                device const float*mw[[buffer(3)]],device const float*E[[buffer(4)]],device float*s[[buffer(5)]],
                constant float&gam[[buffer(6)]],constant uint&n[[buffer(7)]],uint gid[[thread_position_in_grid]]){
    if(gid>=n)return; C5 c=rd(r,mu,mv,mw,E,gid); float p=presC(c,gam); float cs=sqrt(gam*max(p,1e-12f)/c.r);
    s[gid]=sqrt(c.mu*c.mu+c.mv*c.mv+c.mw*c.mw)/c.r+cs;
}
kernel void rk1(device const float*r[[buffer(0)]],device const float*mu[[buffer(1)]],device const float*mv[[buffer(2)]],
                device const float*mw[[buffer(3)]],device const float*E[[buffer(4)]],
                device const float*dr[[buffer(5)]],device const float*dmu[[buffer(6)]],device const float*dmv[[buffer(7)]],
                device const float*dmw[[buffer(8)]],device const float*dE[[buffer(9)]],
                device float*r1[[buffer(10)]],device float*mu1[[buffer(11)]],device float*mv1[[buffer(12)]],
                device float*mw1[[buffer(13)]],device float*E1[[buffer(14)]],
                constant float&dt[[buffer(15)]],constant uint&n[[buffer(16)]],uint g[[thread_position_in_grid]]){
    if(g>=n)return; r1[g]=r[g]+dt*dr[g]; mu1[g]=mu[g]+dt*dmu[g]; mv1[g]=mv[g]+dt*dmv[g]; mw1[g]=mw[g]+dt*dmw[g]; E1[g]=E[g]+dt*dE[g];
}
kernel void rk2(device float*r[[buffer(0)]],device float*mu[[buffer(1)]],device float*mv[[buffer(2)]],
                device float*mw[[buffer(3)]],device float*E[[buffer(4)]],
                device const float*r1[[buffer(5)]],device const float*mu1[[buffer(6)]],device const float*mv1[[buffer(7)]],
                device const float*mw1[[buffer(8)]],device const float*E1[[buffer(9)]],
                device const float*dr[[buffer(10)]],device const float*dmu[[buffer(11)]],device const float*dmv[[buffer(12)]],
                device const float*dmw[[buffer(13)]],device const float*dE[[buffer(14)]],
                constant float&dt[[buffer(15)]],constant uint&n[[buffer(16)]],uint g[[thread_position_in_grid]]){
    if(g>=n)return;
    r[g]=0.5f*(r[g]+r1[g]+dt*dr[g]); mu[g]=0.5f*(mu[g]+mu1[g]+dt*dmu[g]); mv[g]=0.5f*(mv[g]+mv1[g]+dt*dmv[g]);
    mw[g]=0.5f*(mw[g]+mw1[g]+dt*dmw[g]); E[g]=0.5f*(E[g]+E1[g]+dt*dE[g]);
}
