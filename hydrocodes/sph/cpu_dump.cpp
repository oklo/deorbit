// Dump a realistic SPH state + the CPU's force-evaluation outputs, as the oracle
// for validating the GPU strain/stress and force kernels. Builds a sheared,
// compressed, multi-material cube (frac physics on), runs steps so rho/h/S/D
// develop, then does one consistent density+forces eval at that state and dumps
// inputs (pos,vel,mass,u,rho,hh,mat,S[6],D,eps_act) + outputs (acc,dudt,dSdt[6],
// eps_work). All float64. Header: n, eta,alpha,beta,eps,alpha_u,eps_as.
#include "sph.hpp"
#include <cstdio>
#include <cmath>

int main() {
    System S;
    S.materials.push_back(Material::basalt());   // 0
    S.materials.push_back(Material::iron());      // 1
    S.eta = 1.3; S.eps_as = 0.3; S.damage_on = true;
    const double dx = 30.0; const int L = 22;
    S.h_init = S.eta * dx;                         // <-- the real drivers set this; default 0.01 is wrong here     // 22^3 = 10,648 particles (quick CPU O(N) gates)
    double lo = -L * dx / 2.0, hi = L * dx / 2.0;
    S.set_domain(lo - dx, hi + dx, false, lo - dx, hi + dx, false, lo - dx, hi + dx, false);
    double R = 5.0 * dx;                           // central iron sphere
    for (int i = 0; i < L; i++) for (int j = 0; j < L; j++) for (int k = 0; k < L; k++) {
        double x = lo + (i + 0.5) * dx, y = lo + (j + 0.5) * dx, z = lo + (k + 0.5) * dx;
        bool iron = (x*x + y*y + z*z) < R*R;
        S.cur_mat = iron ? 1 : 0;
        double rho0 = iron ? 7800.0 : 2700.0;
        // shear (develops S_xy) + gentle convergence (develops compression/P) + varied u
        double vx = 80.0 * (y / hi), vz = -40.0 * (z / hi);
        double u0 = 2.0e5 * (1.0 + x / hi);        // 0.1..0.4 MJ/kg gradient
        S.add({x, y, z}, {vx, 0.0, vz}, rho0 * dx*dx*dx, u0, false);
    }
    S.seed_damage(dx);
    S.init();
    auto maxhh = [&]{ double m=0; for(int i=0;i<S.n;i++) m=std::max(m,S.hh[i]); return m; };
    auto maxv  = [&]{ double m=0; for(int i=0;i<S.n;i++) m=std::max(m,norm(S.vel[i])); return m; };
    printf("init: n=%d maxhh=%.1f maxv=%.1f\n", S.n, maxhh(), maxv());
    for (int s = 0; s < 12; s++) {
        double dt = S.timestep(); S.step(dt);
        printf("step %d dt=%.3e maxhh=%.1f maxv=%.1f\n", s, dt, maxhh(), maxv()); fflush(stdout);
    }
    S.compute_density_h(); S.compute_forces();            // consistent eval at the final state

    FILE* f = std::fopen("force_ref.bin", "wb");
    int32_t n = S.n; std::fwrite(&n, 4, 1, f);
    double hdr[6] = { S.eta, S.alpha, S.beta, S.eps, S.alpha_u, S.eps_as };
    std::fwrite(hdr, 8, 6, f);
    for (int i = 0; i < S.n; i++) {
        double rec[29] = {
            S.pos[i].x, S.pos[i].y, S.pos[i].z, S.vel[i].x, S.vel[i].y, S.vel[i].z,
            S.mass[i], S.u[i], S.rho[i], S.hh[i], (double)S.mat[i],
            S.S[i].xx, S.S[i].xy, S.S[i].xz, S.S[i].yy, S.S[i].yz, S.S[i].zz,
            S.D[i], S.eps_act[i],
            S.acc[i].x, S.acc[i].y, S.acc[i].z, S.dudt[i],
            S.dSdt[i].xx, S.dSdt[i].xy, S.dSdt[i].xz, S.dSdt[i].yy, S.dSdt[i].yz, S.dSdt[i].zz };
        std::fwrite(rec, 8, 29, f);
    }
    // one more step at a FIXED dt -> dump the post-step state as the step oracle
    double dtv = S.timestep();
    S.step(dtv);
    std::fwrite(&dtv, 8, 1, f);
    for (int i = 0; i < S.n; i++) {
        double post[13] = { S.pos[i].x, S.pos[i].y, S.pos[i].z, S.vel[i].x, S.vel[i].y, S.vel[i].z,
            S.u[i], S.S[i].xx, S.S[i].xy, S.S[i].xz, S.S[i].yy, S.S[i].yz, S.S[i].zz };
        std::fwrite(post, 8, 13, f);
    }
    std::fclose(f);
    double smax = 0; for (int i = 0; i < S.n; i++) smax = std::max(smax, std::abs(S.S[i].xy));
    printf("dumped force_ref.bin: n=%d  max|S_xy|=%.3e  max hh=%.1f  (frac on)\n",
           S.n, smax, [&]{double m=0;for(int i=0;i<S.n;i++)m=std::max(m,S.hh[i]);return m;}());
    return 0;
}
