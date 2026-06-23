#include <metal_stdlib>
using namespace metal;

struct GMat { float rho0, A, B, a, b, alpha, beta, u0, uiv, ucv, G, Emod, Y; };

inline float kW(float r, float h) {
    float q = r / h, s = 1.0f / (M_PI_F * h * h * h);
    if (q < 1.0f) return s * (1.0f - 1.5f * q * q + 0.75f * q * q * q);
    if (q < 2.0f) { float t = 2.0f - q; return s * 0.25f * t * t * t; }
    return 0.0f;
}
inline float kdWdr(float r, float h) {
    float q = r / h, s = 1.0f / (M_PI_F * h * h * h);
    if (q < 1.0f) return s * (-3.0f * q + 2.25f * q * q) / h;
    if (q < 2.0f) { float t = 2.0f - q; return s * (-0.75f * t * t) / h; }
    return 0.0f;
}
inline int cell_of(float3 p, float lox, float loy, float loz, float cw,
                   int ncx, int ncy, int ncz, thread int& bx, thread int& by, thread int& bz) {
    bx = clamp(int((p.x - lox) / cw), 0, ncx - 1);
    by = clamp(int((p.y - loy) / cw), 0, ncy - 1);
    bz = clamp(int((p.z - loz) / cw), 0, ncz - 1);
    return (bx * ncy + by) * ncz + bz;
}

// ---- counting-sort cell list (same as sph_cells) ----
kernel void cell_hash(device const packed_float3* pos [[buffer(0)]],
                      device int* hash [[buffer(1)]], device atomic_uint* cell_count [[buffer(2)]],
                      device uint* slot [[buffer(3)]],
                      constant float& lox [[buffer(4)]], constant float& loy [[buffer(5)]],
                      constant float& loz [[buffer(6)]], constant float& cw [[buffer(7)]],
                      constant int& ncx [[buffer(8)]], constant int& ncy [[buffer(9)]],
                      constant int& ncz [[buffer(10)]], constant uint& n [[buffer(11)]],
                      uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    int bx, by, bz;
    int c = cell_of(float3(pos[i]), lox, loy, loz, cw, ncx, ncy, ncz, bx, by, bz);
    hash[i] = c;
    slot[i] = atomic_fetch_add_explicit(&cell_count[c], 1u, memory_order_relaxed);
}
kernel void cell_scatter(device const int* hash [[buffer(0)]], device const uint* slot [[buffer(1)]],
                         device const uint* cell_start [[buffer(2)]], device uint* sorted [[buffer(3)]],
                         constant uint& n [[buffer(4)]], uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    sorted[cell_start[hash[i]] + slot[i]] = i;
}

// ---- strain rate -> spin -> Jaumann dS/dt + deviatoric power (matches CPU) ----
// S stored interleaved: Sin[6*i + {xx,xy,xz,yy,yz,zz}].
kernel void strain_stress(device const packed_float3* pos [[buffer(0)]],
                          device const packed_float3* vel [[buffer(1)]],
                          device const float* rho [[buffer(2)]], device const float* hh [[buffer(3)]],
                          device const float* mass [[buffer(4)]], device const float* Sin [[buffer(5)]],
                          device const int* mat [[buffer(6)]],
                          device const uint* sorted [[buffer(7)]], device const uint* cell_start [[buffer(8)]],
                          device const uint* cell_count [[buffer(9)]], constant GMat* mats [[buffer(10)]],
                          device float* dSdt [[buffer(11)]], device float* epsw [[buffer(12)]],
                          constant float& lox [[buffer(13)]], constant float& loy [[buffer(14)]],
                          constant float& loz [[buffer(15)]], constant float& cw [[buffer(16)]],
                          constant int& ncx [[buffer(17)]], constant int& ncy [[buffer(18)]],
                          constant int& ncz [[buffer(19)]], constant uint& n [[buffer(20)]],
                          uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    float G = mats[mat[i]].G;
    if (G <= 0.0f) { for (int k = 0; k < 6; k++) dSdt[6*i+k] = 0.0f; epsw[i] = 0.0f; return; }
    float3 ri = float3(pos[i]), vi = float3(vel[i]);
    int bx, by, bz; cell_of(ri, lox, loy, loz, cw, ncx, ncy, ncz, bx, by, bz);
    float L[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int dx = -1; dx <= 1; dx++) { int cx = bx+dx; if (cx<0||cx>=ncx) continue;
      for (int dy = -1; dy <= 1; dy++) { int cy = by+dy; if (cy<0||cy>=ncy) continue;
        for (int dz = -1; dz <= 1; dz++) { int cz = bz+dz; if (cz<0||cz>=ncz) continue;
          int c = (cx*ncy+cy)*ncz+cz; uint st = cell_start[c], cn = cell_count[c];
          for (uint q = 0; q < cn; q++) { uint j = sorted[st+q]; if (j == i) continue;
            float3 rij = ri - float3(pos[j]); float r = length(rij);
            float hbar = 0.5f*(hh[i]+hh[j]);
            if (r >= 2.0f*hbar || r == 0.0f) continue;
            float3 gW = rij * (kdWdr(r, hbar) / r);
            float Vj = mass[j] / rho[j];
            float3 dv = float3(vel[j]) - vi;
            float d[3] = {dv.x, dv.y, dv.z}, g[3] = {gW.x, gW.y, gW.z};
            for (int a=0;a<3;a++) for (int b=0;b<3;b++) L[a][b] += Vj*d[a]*g[b];
          } } } }
    float exx=L[0][0],eyy=L[1][1],ezz=L[2][2];
    float exy=0.5f*(L[0][1]+L[1][0]),exz=0.5f*(L[0][2]+L[2][0]),eyz=0.5f*(L[1][2]+L[2][1]);
    float Rxy=0.5f*(L[0][1]-L[1][0]),Rxz=0.5f*(L[0][2]-L[2][0]),Ryz=0.5f*(L[1][2]-L[2][1]);
    float tr=(exx+eyy+ezz)/3.0f;
    float sxx=Sin[6*i],sxy=Sin[6*i+1],sxz=Sin[6*i+2],syy=Sin[6*i+3],syz=Sin[6*i+4],szz=Sin[6*i+5];
    float Rm[3][3]={{0,Rxy,Rxz},{-Rxy,0,Ryz},{-Rxz,-Ryz,0}};
    float Sm[3][3]={{sxx,sxy,sxz},{sxy,syy,syz},{sxz,syz,szz}};
    float J[3][3];
    for (int a=0;a<3;a++) for (int b=0;b<3;b++){ float rs=0,sr=0;
        for (int c=0;c<3;c++){ rs+=Rm[a][c]*Sm[c][b]; sr+=Sm[a][c]*Rm[c][b]; } J[a][b]=rs-sr; }
    dSdt[6*i+0]=2.0f*G*(exx-tr)+J[0][0];
    dSdt[6*i+3]=2.0f*G*(eyy-tr)+J[1][1];
    dSdt[6*i+5]=2.0f*G*(ezz-tr)+J[2][2];
    dSdt[6*i+1]=2.0f*G*exy+J[0][1];
    dSdt[6*i+2]=2.0f*G*exz+J[0][2];
    dSdt[6*i+4]=2.0f*G*eyz+J[1][2];
    epsw[i]=(sxx*exx+syy*eyy+szz*ezz+2.0f*(sxy*exy+sxz*exz+syz*eyz))/rho[i];
}

// ---- Tillotson EOS (analytic sound speed), same as sph_eos.metal ----
inline float till(float rho, float u, constant GMat& m) {
    if (rho <= 0.0f) return 0.0f;
    float eta = rho / m.rho0, mu = eta - 1.0f;
    float w0 = u / (m.u0 * eta * eta) + 1.0f;
    float Pc = (m.a + m.b / w0) * rho * u + m.A * mu + m.B * mu * mu;
    if (rho >= m.rho0 || u <= m.uiv) return Pc;
    float z = m.rho0 / rho - 1.0f;
    float Pe = m.a * rho * u + (m.b * rho * u / w0 + m.A * mu * exp(-m.beta * z)) * exp(-m.alpha * z * z);
    if (u >= m.ucv) return Pe;
    return ((u - m.uiv) * Pe + (m.ucv - u) * Pc) / (m.ucv - m.uiv);
}
inline void cond(float rho, float u, constant GMat& m, thread float& P, thread float& dPdr, thread float& dPdu) {
    float eta = rho/m.rho0, mu = eta-1.0f, w = 1.0f + u/(m.u0*eta*eta), aob = m.a + m.b/w;
    P = aob*rho*u + m.A*mu + m.B*mu*mu;
    dPdu = aob*rho - m.b*rho*u/(w*w*m.u0*eta*eta);
    float dwdr = -2.0f*u/(m.u0*eta*eta*eta*m.rho0);
    dPdr = aob*u + rho*u*(-m.b/(w*w))*dwdr + m.A/m.rho0 + 2.0f*m.B*mu/m.rho0;
}
inline void expd(float rho, float u, constant GMat& m, thread float& P, thread float& dPdr, thread float& dPdu) {
    float eta=rho/m.rho0, mu=eta-1.0f, w=1.0f+u/(m.u0*eta*eta);
    float z=m.rho0/rho-1.0f, dzdr=-m.rho0/(rho*rho);
    float E1=exp(-m.alpha*z*z), E2=exp(-m.beta*z), G=m.b*rho*u/w+m.A*mu*E2;
    P=m.a*rho*u+G*E1;
    float dwdr=-2.0f*u/(m.u0*eta*eta*eta*m.rho0);
    float dGu=m.b*rho/w - m.b*rho*u/(w*w*m.u0*eta*eta);
    dPdu=m.a*rho+dGu*E1;
    float dbruw=m.b*u/w + m.b*rho*u*(-1.0f/(w*w))*dwdr;
    float dGr=dbruw + m.A*E2/m.rho0 + m.A*mu*(E2*(-m.beta)*dzdr);
    dPdr=m.a*u + dGr*E1 + G*(E1*(-m.alpha*2.0f*z)*dzdr);
}
inline float ssound(float rho, float u, constant GMat& m) {
    if (rho<=0.0f) return 1e-3f;
    float P,dPdr,dPdu;
    if (rho>=m.rho0 || u<=m.uiv) cond(rho,u,m,P,dPdr,dPdu);
    else if (u>=m.ucv) expd(rho,u,m,P,dPdr,dPdu);
    else { float Pc,dPcr,dPcu,Pe,dPer,dPeu; cond(rho,u,m,Pc,dPcr,dPcu); expd(rho,u,m,Pe,dPer,dPeu);
        float inv=1.0f/(m.ucv-m.uiv);
        P=((u-m.uiv)*Pe+(m.ucv-u)*Pc)*inv; dPdu=((Pe+(u-m.uiv)*dPeu)+(-Pc+(m.ucv-u)*dPcu))*inv;
        dPdr=((u-m.uiv)*dPer+(m.ucv-u)*dPcr)*inv; }
    return sqrt(max(dPdr + P/(rho*rho)*dPdu, 1e-6f));
}
inline float3 Sdot(float xx,float xy,float xz,float yy,float yz,float zz, float3 g){
    return float3(xx*g.x+xy*g.y+xz*g.z, xy*g.x+yy*g.y+yz*g.z, xz*g.x+yz*g.y+zz*g.z);
}

// per-particle EOS + damage scaling (Pf, csig, dscale, Rart)
kernel void eos_full(device const float* rho [[buffer(0)]], device const float* u [[buffer(1)]],
                     device const int* mat [[buffer(2)]], device const float* D [[buffer(3)]],
                     device float* Pf [[buffer(4)]], device float* csig [[buffer(5)]],
                     device float* dscale [[buffer(6)]], device float* Rart [[buffer(7)]],
                     constant GMat* mats [[buffer(8)]], constant float& eps_as [[buffer(9)]],
                     constant int& damage_on [[buffer(10)]], constant uint& n [[buffer(11)]],
                     uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    constant GMat& m = mats[mat[i]];
    float P = till(rho[i], u[i], m), cs = ssound(rho[i], u[i], m);
    csig[i] = sqrt(cs*cs + (4.0f/3.0f)*m.G/max(rho[i],1e-30f));
    float ds = 1.0f, pf = P;
    if (damage_on) { ds = 1.0f - D[i]; if (P < 0.0f) pf = ds * P; }
    dscale[i] = ds; Pf[i] = pf;
    Rart[i] = (eps_as > 0.0f && pf < 0.0f) ? eps_as*(-pf)/(rho[i]*rho[i]) : 0.0f;
}

// momentum + energy: -grad P (damaged) + artificial viscosity + deviatoric stress
// (damage-scaled) + Monaghan artificial stress. Conductivity off (alpha_u=0).
kernel void forces(device const packed_float3* pos [[buffer(0)]], device const packed_float3* vel [[buffer(1)]],
                   device const float* rho [[buffer(2)]], device const float* hh [[buffer(3)]],
                   device const float* mass [[buffer(4)]], device const float* Sin [[buffer(5)]],
                   device const float* Pf [[buffer(6)]], device const float* csig [[buffer(7)]],
                   device const float* dscale [[buffer(8)]], device const float* Rart [[buffer(9)]],
                   device const float* epsw [[buffer(10)]], device const uint* sorted [[buffer(11)]],
                   device const uint* cell_start [[buffer(12)]], device const uint* cell_count [[buffer(13)]],
                   device packed_float3* acc [[buffer(14)]], device float* dudt [[buffer(15)]],
                   constant float& lox [[buffer(16)]], constant float& loy [[buffer(17)]],
                   constant float& loz [[buffer(18)]], constant float& cw [[buffer(19)]],
                   constant int& ncx [[buffer(20)]], constant int& ncy [[buffer(21)]],
                   constant int& ncz [[buffer(22)]], constant float& eta [[buffer(23)]],
                   constant float& alpha [[buffer(24)]], constant float& beta [[buffer(25)]],
                   constant float& eps [[buffer(26)]], constant float& eps_as [[buffer(27)]],
                   constant uint& n [[buffer(28)]], uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    float3 ri = float3(pos[i]), vi = float3(vel[i]);
    float iri2 = 1.0f/(rho[i]*rho[i]), Pi_over = Pf[i]*iri2, di = dscale[i];
    float sxx=Sin[6*i],sxy=Sin[6*i+1],sxz=Sin[6*i+2],syy=Sin[6*i+3],syz=Sin[6*i+4],szz=Sin[6*i+5];
    int bx,by,bz; cell_of(ri,lox,loy,loz,cw,ncx,ncy,ncz,bx,by,bz);
    float3 a = float3(0.0f); float du = 0.0f;
    for (int dx=-1;dx<=1;dx++){ int cx=bx+dx; if(cx<0||cx>=ncx) continue;
      for (int dy=-1;dy<=1;dy++){ int cy=by+dy; if(cy<0||cy>=ncy) continue;
        for (int dz=-1;dz<=1;dz++){ int cz=bz+dz; if(cz<0||cz>=ncz) continue;
          int c=(cx*ncy+cy)*ncz+cz; uint s0=cell_start[c], cn=cell_count[c];
          for (uint q=0;q<cn;q++){ uint j=sorted[s0+q]; if(j==i) continue;
            float3 rij = ri - float3(pos[j]); float r = length(rij);
            float hbar = 0.5f*(hh[i]+hh[j]); if (r>=2.0f*hbar || r==0.0f) continue;
            float3 gradW = rij*(kdWdr(r,hbar)/r);
            float3 vij = vi - float3(vel[j]); float rbar = 0.5f*(rho[i]+rho[j]);
            float dvr = dot(vij,rij), Pij = 0.0f;
            if (dvr < 0.0f){ float mu=hbar*dvr/(r*r+eps*hbar*hbar); float cbar=0.5f*(csig[i]+csig[j]);
                Pij = (-alpha*cbar*mu + beta*mu*mu)/rbar; }
            float jrj2 = 1.0f/(rho[j]*rho[j]); float term = Pi_over + Pf[j]*jrj2 + Pij;
            du += 0.5f*mass[j]*term*dot(vij,gradW);
            float mterm = term;
            if (eps_as > 0.0f){ float wdp=kW(hbar/eta,hbar); float fr=(wdp>0.0f? kW(r,hbar)/wdp:0.0f);
                float fn=fr*fr; fn*=fn; mterm += (Rart[i]+Rart[j])*fn; }
            a -= gradW*(mass[j]*mterm);
            float3 si = Sdot(sxx,sxy,sxz,syy,syz,szz, gradW);
            float3 sj = Sdot(Sin[6*j],Sin[6*j+1],Sin[6*j+2],Sin[6*j+3],Sin[6*j+4],Sin[6*j+5], gradW);
            a += (si*(iri2*di) + sj*(jrj2*dscale[j])) * mass[j];
          } } } }
    acc[i] = packed_float3(a); dudt[i] = du + epsw[i]*di;
}

// ---- step kernels (milestone 4) ----
inline float eig_max(float sxx,float sxy,float sxz,float syy,float syz,float szz){  // Smith trig, == CPU
    float p1 = sxy*sxy + sxz*sxz + syz*syz, q = (sxx+syy+szz)/3.0f;
    if (p1 < 1e-30f) return max(max(sxx,syy),szz);
    float p2 = (sxx-q)*(sxx-q)+(syy-q)*(syy-q)+(szz-q)*(szz-q)+2.0f*p1, p = sqrt(p2/6.0f);
    float a=(sxx-q)/p,b=(syy-q)/p,c=(szz-q)/p,d=sxy/p,e=sxz/p,f=syz/p;
    float r = (a*(b*c-f*f)-d*(d*c-f*e)+e*(d*f-b*e))/2.0f;
    r = clamp(r,-1.0f,1.0f);
    return q + 2.0f*p*cos(acos(r)/3.0f);
}
kernel void compute_sigmax(device const float* rho[[buffer(0)]], device const float* u[[buffer(1)]],
                           device const int* mat[[buffer(2)]], device const float* Sin[[buffer(3)]],
                           device float* sigmax[[buffer(4)]], constant GMat* mats[[buffer(5)]],
                           constant uint& n[[buffer(6)]], uint i[[thread_position_in_grid]]){
    if(i>=n) return;
    float P = till(rho[i],u[i],mats[mat[i]]);
    sigmax[i] = -P + eig_max(Sin[6*i],Sin[6*i+1],Sin[6*i+2],Sin[6*i+3],Sin[6*i+4],Sin[6*i+5]);
}
kernel void kick(device packed_float3* vel[[buffer(0)]], device const packed_float3* acc[[buffer(1)]],
                 device const uchar* fixed[[buffer(2)]], constant float& cdt[[buffer(3)]],
                 constant uint& n[[buffer(4)]], uint i[[thread_position_in_grid]]){
    if(i>=n || fixed[i]) return;
    vel[i] = packed_float3(float3(vel[i]) + float3(acc[i])*cdt);
}
kernel void drift(device packed_float3* pos[[buffer(0)]], device const packed_float3* vel[[buffer(1)]],
                  device const uchar* fixed[[buffer(2)]], constant float& dt[[buffer(3)]],
                  constant uint& n[[buffer(4)]], uint i[[thread_position_in_grid]]){
    if(i>=n || fixed[i]) return;
    pos[i] = packed_float3(float3(pos[i]) + float3(vel[i])*dt);
}
// adaptive-h density GATHER (rho from own hi, then hh = clamp(eta*(m/rho)^1/3))
kernel void density_adaptive(device const packed_float3* pos[[buffer(0)]], device const float* mass[[buffer(1)]],
                             device const uint* sorted[[buffer(2)]], device const uint* cell_start[[buffer(3)]],
                             device const uint* cell_count[[buffer(4)]], device float* rho[[buffer(5)]],
                             device float* hh[[buffer(6)]],
                             constant float& lox[[buffer(7)]], constant float& loy[[buffer(8)]],
                             constant float& loz[[buffer(9)]], constant float& cw[[buffer(10)]],
                             constant int& ncx[[buffer(11)]], constant int& ncy[[buffer(12)]],
                             constant int& ncz[[buffer(13)]], constant float& eta[[buffer(14)]],
                             constant float& h_init[[buffer(15)]], constant uint& n[[buffer(16)]],
                             uint i[[thread_position_in_grid]]){
    if(i>=n) return;
    float3 ri=float3(pos[i]); float hi=hh[i];
    int bx,by,bz; cell_of(ri,lox,loy,loz,cw,ncx,ncy,ncz,bx,by,bz);
    float s=0.0f;
    for(int dx=-1;dx<=1;dx++){int cx=bx+dx;if(cx<0||cx>=ncx)continue;
      for(int dy=-1;dy<=1;dy++){int cy=by+dy;if(cy<0||cy>=ncy)continue;
        for(int dz=-1;dz<=1;dz++){int cz=bz+dz;if(cz<0||cz>=ncz)continue;
          int c=(cx*ncy+cy)*ncz+cz; uint s0=cell_start[c],cn=cell_count[c];
          for(uint q=0;q<cn;q++){uint j=sorted[s0+q]; float3 dij=ri-float3(pos[j]);
            float r2=dot(dij,dij), th=2.0f*hi; if(r2<th*th) s+=mass[j]*kW(sqrt(r2),hi);}}}}
    float rho_i=max(s,1e-30f); rho[i]=rho_i;
    hh[i]=clamp(eta*pow(mass[i]/rho_i,1.0f/3.0f), 0.5f*h_init, 2.0f*h_init);
}
kernel void grow_damage(device const float* sigmax[[buffer(0)]], device float* D[[buffer(1)]],
                        device const float* eps_act[[buffer(2)]], device const int* mat[[buffer(3)]],
                        device const float* rho[[buffer(4)]], constant GMat* mats[[buffer(5)]],
                        constant float& h_init[[buffer(6)]], constant float& eta[[buffer(7)]],
                        constant float& dt[[buffer(8)]], constant uint& n[[buffer(9)]], uint i[[thread_position_in_grid]]){
    if(i>=n) return; constant GMat& m=mats[mat[i]];
    if(m.Emod<=0.0f || D[i]>=1.0f) return;
    float eps_local=sigmax[i]/m.Emod;
    if(sigmax[i]>0.0f && eps_local>eps_act[i]){
        float cg=0.4f*sqrt(m.Emod/max(rho[i],1e-30f)), Rs=0.5f*h_init/eta;
        float d13=pow(D[i],1.0f/3.0f)+(cg/Rs)*dt;
        D[i]=min(1.0f,d13*d13*d13);
    }
}
// kick2 + trapezoidal u, S + von Mises radial return
kernel void finish(device packed_float3* vel[[buffer(0)]], device const packed_float3* acc[[buffer(1)]],
                   device float* u[[buffer(2)]], device float* Sin[[buffer(3)]],
                   device const float* dudt[[buffer(4)]], device const float* dudt_old[[buffer(5)]],
                   device const float* dSdt[[buffer(6)]], device const float* dSdt_old[[buffer(7)]],
                   device const uchar* fixed[[buffer(8)]], device const int* mat[[buffer(9)]],
                   constant GMat* mats[[buffer(10)]], constant float& dt[[buffer(11)]],
                   constant uint& n[[buffer(12)]], uint i[[thread_position_in_grid]]){
    if(i>=n || fixed[i]) return;
    vel[i]=packed_float3(float3(vel[i])+float3(acc[i])*(0.5f*dt));
    u[i]+=0.5f*dt*(dudt_old[i]+dudt[i]); if(u[i]<1e-12f) u[i]=1e-12f;
    for(int k=0;k<6;k++) Sin[6*i+k]+=0.5f*dt*(dSdt_old[6*i+k]+dSdt[6*i+k]);
    float Y=mats[mat[i]].Y;
    if(Y<=0.0f){ for(int k=0;k<6;k++) Sin[6*i+k]=0.0f; return; }
    float sxx=Sin[6*i],sxy=Sin[6*i+1],sxz=Sin[6*i+2],syy=Sin[6*i+3],syz=Sin[6*i+4],szz=Sin[6*i+5];
    float J2=0.5f*(sxx*sxx+syy*syy+szz*szz)+(sxy*sxy+sxz*sxz+syz*syz);
    float vm=sqrt(3.0f*J2);
    if(vm>Y){ float f=Y/vm; for(int k=0;k<6;k++) Sin[6*i+k]*=f; }
}

// FUSED strain/stress + forces: ONE neighbour traversal computes the velocity
// gradient AND the momentum/energy sums, then post-loop derives dS/dt, eps_work
// (folded into dudt), acc. Replaces the separate strain_stress + forces passes.
kernel void strain_forces(device const packed_float3* pos[[buffer(0)]], device const packed_float3* vel[[buffer(1)]],
        device const float* rho[[buffer(2)]], device const float* hh[[buffer(3)]], device const float* mass[[buffer(4)]],
        device const float* Sin[[buffer(5)]], device const int* mat[[buffer(6)]], device const float* Pf[[buffer(7)]],
        device const float* csig[[buffer(8)]], device const float* dscale[[buffer(9)]], device const float* Rart[[buffer(10)]],
        device const uint* sorted[[buffer(11)]], device const uint* cell_start[[buffer(12)]], device const uint* cell_count[[buffer(13)]],
        constant GMat* mats[[buffer(14)]], device float* dSdt[[buffer(15)]], device packed_float3* acc[[buffer(16)]],
        device float* dudt[[buffer(17)]], constant float& lox[[buffer(18)]], constant float& loy[[buffer(19)]],
        constant float& loz[[buffer(20)]], constant float& cw[[buffer(21)]], constant int& ncx[[buffer(22)]],
        constant int& ncy[[buffer(23)]], constant int& ncz[[buffer(24)]], constant float& eta[[buffer(25)]],
        constant float& alpha[[buffer(26)]], constant float& beta[[buffer(27)]], constant float& eps[[buffer(28)]],
        constant float& eps_as[[buffer(29)]], constant uint& n[[buffer(30)]], uint i[[thread_position_in_grid]]){
    if(i>=n) return;
    float G=mats[mat[i]].G;
    float3 ri=float3(pos[i]), vi=float3(vel[i]);
    float iri2=1.0f/(rho[i]*rho[i]), Pi_over=Pf[i]*iri2, di=dscale[i];
    float sxx=Sin[6*i],sxy=Sin[6*i+1],sxz=Sin[6*i+2],syy=Sin[6*i+3],syz=Sin[6*i+4],szz=Sin[6*i+5];
    int bx,by,bz; cell_of(ri,lox,loy,loz,cw,ncx,ncy,ncz,bx,by,bz);
    float L[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    float3 a=float3(0.0f); float du=0.0f;
    for(int dx=-1;dx<=1;dx++){int cx=bx+dx;if(cx<0||cx>=ncx)continue;
      for(int dy=-1;dy<=1;dy++){int cy=by+dy;if(cy<0||cy>=ncy)continue;
        for(int dz=-1;dz<=1;dz++){int cz=bz+dz;if(cz<0||cz>=ncz)continue;
          int c=(cx*ncy+cy)*ncz+cz; uint s0=cell_start[c],cn=cell_count[c];
          for(uint q=0;q<cn;q++){ uint j=sorted[s0+q]; if(j==i) continue;
            float3 rij=ri-float3(pos[j]); float r2=dot(rij,rij); float hbar=0.5f*(hh[i]+hh[j]);
            float th=2.0f*hbar; if(r2>=th*th||r2==0.0f) continue; float r=sqrt(r2);
            float3 gradW=rij*(kdWdr(r,hbar)/r);
            float Vj=mass[j]/rho[j]; float3 dvj=float3(vel[j])-vi;
            float dd[3]={dvj.x,dvj.y,dvj.z}, gg[3]={gradW.x,gradW.y,gradW.z};
            for(int aa=0;aa<3;aa++)for(int bb=0;bb<3;bb++) L[aa][bb]+=Vj*dd[aa]*gg[bb];
            float3 vij=vi-float3(vel[j]); float rbar=0.5f*(rho[i]+rho[j]); float dvr=dot(vij,rij), Pij=0.0f;
            if(dvr<0.0f){ float mu=hbar*dvr/(r*r+eps*hbar*hbar); float cbar=0.5f*(csig[i]+csig[j]); Pij=(-alpha*cbar*mu+beta*mu*mu)/rbar; }
            float jrj2=1.0f/(rho[j]*rho[j]); float term=Pi_over+Pf[j]*jrj2+Pij;
            du += 0.5f*mass[j]*term*dot(vij,gradW);
            float mterm=term;
            if(eps_as>0.0f){ float wdp=kW(hbar/eta,hbar); float fr=(wdp>0.0f?kW(r,hbar)/wdp:0.0f); float fn=fr*fr; fn*=fn; mterm+=(Rart[i]+Rart[j])*fn; }
            a -= gradW*(mass[j]*mterm);
            float3 si=Sdot(sxx,sxy,sxz,syy,syz,szz,gradW);
            float3 sj=Sdot(Sin[6*j],Sin[6*j+1],Sin[6*j+2],Sin[6*j+3],Sin[6*j+4],Sin[6*j+5],gradW);
            a += (si*(iri2*di)+sj*(jrj2*dscale[j]))*mass[j];
          }}}}
    float epsw=0.0f;
    if(G>0.0f){
        float exx=L[0][0],eyy=L[1][1],ezz=L[2][2];
        float exy=0.5f*(L[0][1]+L[1][0]),exz=0.5f*(L[0][2]+L[2][0]),eyz=0.5f*(L[1][2]+L[2][1]);
        float Rxy=0.5f*(L[0][1]-L[1][0]),Rxz=0.5f*(L[0][2]-L[2][0]),Ryz=0.5f*(L[1][2]-L[2][1]);
        float tr=(exx+eyy+ezz)/3.0f;
        float Rm[3][3]={{0,Rxy,Rxz},{-Rxy,0,Ryz},{-Rxz,-Ryz,0}};
        float Sm[3][3]={{sxx,sxy,sxz},{sxy,syy,syz},{sxz,syz,szz}};
        float J[3][3];
        for(int aa=0;aa<3;aa++)for(int bb=0;bb<3;bb++){ float rs=0,sr=0;
            for(int c=0;c<3;c++){ rs+=Rm[aa][c]*Sm[c][bb]; sr+=Sm[aa][c]*Rm[c][bb]; } J[aa][bb]=rs-sr; }
        dSdt[6*i+0]=2.0f*G*(exx-tr)+J[0][0]; dSdt[6*i+3]=2.0f*G*(eyy-tr)+J[1][1]; dSdt[6*i+5]=2.0f*G*(ezz-tr)+J[2][2];
        dSdt[6*i+1]=2.0f*G*exy+J[0][1]; dSdt[6*i+2]=2.0f*G*exz+J[0][2]; dSdt[6*i+4]=2.0f*G*eyz+J[1][2];
        epsw=(sxx*exx+syy*eyy+szz*ezz+2.0f*(sxy*exy+sxz*exz+syz*eyz))/rho[i];
    } else { for(int k=0;k<6;k++) dSdt[6*i+k]=0.0f; }
    acc[i]=packed_float3(a); dudt[i]=du+epsw*di;
}

// ---- particle reordering (coalesced neighbour reads) ----
// gather 4-byte-word arrays into cell-sorted order: out[k]=in[sorted[k]] (W words/particle)
kernel void gather_w(device const uint* in[[buffer(0)]], device uint* out[[buffer(1)]],
                     device const uint* sorted[[buffer(2)]], constant uint& W[[buffer(3)]],
                     constant uint& n[[buffer(4)]], uint k[[thread_position_in_grid]]){
    if(k>=n) return; uint src=sorted[k]; for(uint w=0;w<W;w++) out[k*W+w]=in[src*W+w];
}
kernel void gather_u8(device const uchar* in[[buffer(0)]], device uchar* out[[buffer(1)]],
                      device const uint* sorted[[buffer(2)]], constant uint& n[[buffer(3)]],
                      uint k[[thread_position_in_grid]]){ if(k<n) out[k]=in[sorted[k]]; }
kernel void set_iota(device uint* sorted[[buffer(0)]], constant uint& n[[buffer(1)]],
                     uint k[[thread_position_in_grid]]){ if(k<n) sorted[k]=k; }
