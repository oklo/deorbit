// Validation of the strength model (uses the real System methods):
//  A. elastic shear: pure-shear velocity field -> dS_xy/dt = 2 G eps_xy
//  B. von Mises yield: yield_limit caps the von Mises stress at Y (direction kept)
//  C. elastic longitudinal wave speed = sqrt((K + 4G/3)/rho0)  [K = A here]
#include <cstdio>
#include <cmath>
#include "sph.hpp"

static void basalt_block(System& S, double x0, double x1, double dx, double Lyz) {
    S.materials.push_back(Material::basalt());
    S.eta = 1.3; S.h_init = S.eta * dx;
    S.set_domain(x0, x1, false, 0, Lyz, false, 0, Lyz, false);   // non-periodic by default
    double rho0 = Material::basalt().rho0, m = rho0 * dx * dx * dx;
    int nx = (int)std::round((x1 - x0) / dx), ny = (int)std::round(Lyz / dx);
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < ny; k++)
                S.add({x0 + (i + 0.5) * dx, (j + 0.5) * dx, (k + 0.5) * dx},
                      {0, 0, 0}, m, 0.0);
}

static int test_shear() {
    System S; basalt_block(S, -0.2, 0.2, 0.02, 0.16);   // non-periodic (a linear
                                                  // shear field is incompatible with y,z wrap)
    double G = Material::basalt().G, k = 1.0;     // pure shear v=(k y, k x, 0)
    for (int i = 0; i < S.n; i++) S.vel[i] = {k * S.pos[i].y, k * S.pos[i].x, 0};
    S.init();                                     // density + h
    S.compute_strain_stress();                    // fills dSdt from velocity grad
    // pick a fully-interior particle (nearest the block centre 0,0.08,0.08)
    Vec3 c{0, 0.08, 0.08};
    int best = 0; double bd = 1e30;
    for (int i = 0; i < S.n; i++) { double d = norm(S.pos[i] - c); if (d < bd) { bd = d; best = i; } }
    double dSxy = S.dSdt[best].xy, expect = 2 * G * k;   // eps_xy = k
    double ratio = dSxy / expect;
    bool ok = ratio > 0.8 && ratio < 1.2;
    printf("[A shear] dS_xy/dt = %.3e  expect 2G*eps_xy = %.3e  ratio=%.3f -> %s\n",
           dSxy, expect, ratio, ok ? "OK" : "FAIL");
    return ok;
}

static int test_yield() {
    System S; S.materials.push_back(Material::basalt());
    S.add({0, 0, 0}, {0, 0, 0}, 1.0, 0.0);
    double Y = Material::basalt().Y;
    S.S[0] = Sym3{1.0e9, -0.4e9, -0.6e9, 3.0e8, 1.0e8, -2.0e8};  // exceeds Y
    Sym3 before = S.S[0];
    S.yield_limit(0);
    Sym3 a = S.S[0];
    double J2 = 0.5 * (a.xx * a.xx + a.yy * a.yy + a.zz * a.zz)
                + (a.xy * a.xy + a.xz * a.xz + a.yz * a.yz);
    double vm = std::sqrt(3 * J2);
    double dir = a.xx / before.xx;                // scale factor; check uniform
    bool ok = std::abs(vm - Y) / Y < 1e-6
              && std::abs(a.yy / before.yy - dir) < 1e-9
              && std::abs(a.xy / before.xy - dir) < 1e-9;
    printf("[B yield] von Mises after = %.4e  Y = %.4e  uniform-scale=%s -> %s\n",
           vm, Y, (std::abs(a.yy / before.yy - dir) < 1e-9) ? "yes" : "no", ok ? "OK" : "FAIL");
    return ok;
}

static int test_wave() {
    // Well-resolved (w >> dx) small-amplitude Gaussian velocity pulse splits
    // into two elastic waves; track the right-going PEAK (noise-robust) -> c_L.
    System S; basalt_block(S, -0.6, 0.6, 0.01, 0.08);
    S.per[1] = S.per[2] = true;                   // y,z periodic (pulse is y,z-uniform)
    Material b = Material::basalt();
    double cL = std::sqrt((b.A + 4.0 / 3.0 * b.G) / b.rho0);
    double v0 = 1.0, w = 0.06;
    for (int i = 0; i < S.n; i++) {
        double x = S.pos[i].x;
        S.vel[i] = {v0 * std::exp(-(x * x) / (w * w)), 0, 0};
        if (std::abs(x) > 0.55) S.fixed[i] = 1;    // frozen ends (no reflection)
    }
    S.init();
    double t1 = 6.0e-5, t = 0;
    while (t < t1) { double dt = S.timestep(); if (t + dt > t1) dt = t1 - t; S.step(dt); t += dt; }
    // location of the right-going pulse peak (x>0.05)
    double xpk = 0, vpk = 0;
    for (int i = 0; i < S.n; i++)
        if (S.pos[i].x > 0.05 && std::abs(S.vel[i].x) > vpk) { vpk = std::abs(S.vel[i].x); xpk = S.pos[i].x; }
    double speed = xpk / t;
    double ratio = speed / cL;
    bool ok = ratio > 0.85 && ratio < 1.15;
    printf("[C wave] right peak at x=%.3f (v=%.2f) after %.2e s -> speed=%.0f m/s  "
           "c_L=%.0f m/s  ratio=%.3f -> %s\n", xpk, vpk, t, speed, cL, ratio, ok ? "OK" : "FAIL");
    return ok;
}

int main() {
    printf("Strength model validation:\n");
    int ok = 1;
    ok &= test_shear();
    ok &= test_yield();
    ok &= test_wave();
    printf("%s\n", ok ? "\nALL PASS" : "\nSOME CHECKS FAILED");
    return ok ? 0 : 1;
}
