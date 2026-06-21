// Validation gate 1: 3D Sod shock tube (ideal gas, gamma=1.4).
// A thin 3D slab, periodic in y,z, discontinuity at x=0. Compared to the exact
// Riemann solution by analyze_sod.py. This exercises the 3D kernel, neighbour
// search, artificial viscosity and time integration before any impact physics.
#include <cstdio>
#include "sph.hpp"

int main() {
    System S;
    S.eos.gamma = 1.4;
    S.alpha = 1.0; S.beta = 2.0;      // artificial viscosity (standard)
    S.alpha_u = 1.0;                  // artificial conductivity (Price 2008)
    const double aL = 0.005;          // left lattice spacing
    const double aR = 2.0 * aL;       // right spacing (rho ratio 8 -> 2x in 3D)
    const int nxL = 100, nxR = 50;
    const int nyL = 10, nyR = 5;      // y,z layers (right uses 2x spacing)
    S.Ly = S.Lz = nyL * aL;           // = 0.05, periodic
    S.xmin = -0.5; S.xmax = 0.5;
    S.h = 1.6 * aL;

    const double rhoL = 1.0, rhoR = 0.125, PL = 1.0, PR = 0.1;
    const double m = rhoL * aL * aL * aL;
    const double uL = PL / ((S.eos.gamma - 1.0) * rhoL);  // 2.5
    const double uR = PR / ((S.eos.gamma - 1.0) * rhoR);  // 2.0
    const double xfix = 2.0 * S.h;    // freeze particles within 2h of the ends

    // left block
    for (int ix = 0; ix < nxL; ix++)
        for (int iy = 0; iy < nyL; iy++)
            for (int iz = 0; iz < nyL; iz++) {
                double x = -0.5 + (ix + 0.5) * aL;
                double y = (iy + 0.5) * aL, z = (iz + 0.5) * aL;
                S.add({x, y, z}, {0, 0, 0}, m, uL, x < -0.5 + xfix);
            }
    // right block
    for (int ix = 0; ix < nxR; ix++)
        for (int iy = 0; iy < nyR; iy++)
            for (int iz = 0; iz < nyR; iz++) {
                double x = 0.0 + (ix + 0.5) * aR;
                double y = (iy + 0.5) * aR, z = (iz + 0.5) * aR;
                S.add({x, y, z}, {0, 0, 0}, m, uR, x > 0.5 - xfix);
            }

    printf("Sod: %d particles, h=%.4f, m=%.3e\n", S.n, S.h, m);
    S.init();

    const double t_end = 0.2;
    double t = 0; int nstep = 0;
    while (t < t_end) {
        double dt = S.timestep();
        if (t + dt > t_end) dt = t_end - t;
        S.step(dt);
        t += dt; nstep++;
        if (nstep % 50 == 0) printf("  step %d  t=%.4f  dt=%.2e\n", nstep, t, dt);
    }
    printf("done: %d steps to t=%.3f\n", nstep, t);

    FILE* f = std::fopen("sod_out.txt", "w");
    std::fprintf(f, "# x rho vx P\n");
    for (int i = 0; i < S.n; i++)
        if (!S.fixed[i])
            std::fprintf(f, "%.6f %.6f %.6f %.6f\n", S.pos[i].x, S.rho[i], S.vel[i].x, S.P[i]);
    std::fclose(f);
    printf("wrote sod_out.txt\n");
    return 0;
}
