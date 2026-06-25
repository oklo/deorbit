// Gate 2: crater pi-scaling (strength regime, no body-force gravity).
// Vertical iron sphere into a basalt half-space at speed U. Measure the
// transient crater radius via a radial surface profile. Run at several U; the
// crater radius should follow R ~ U^mu (mu ~ 0.4-0.6, the coupling exponent).
//
//   ./pi_scaling <U_km_s> [dx]   ->  prints: U R_crater depth R/a  (+ progress)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "sph.hpp"

int main(int argc, char** argv) {
    double U = (argc > 1 ? atof(argv[1]) : 2.0) * 1000.0;   // m/s
    double dx = (argc > 2 ? atof(argv[2]) : 0.0125);        // spacing
    double a = 0.05;                                        // impactor radius (CPPR=a/dx)
    double Lx = 0.45, Lz = 0.35;                            // target full-width, depth

    System S;
    S.materials.push_back(Material::basalt());   // 0 = target
    S.materials.push_back(Material::iron());      // 1 = impactor
    S.eta = 1.3; S.h_init = S.eta * dx;
    // non-periodic walled container; bounds cover target + impactor + ejecta headroom
    S.set_domain(-Lx / 2 - 0.05, Lx / 2 + 0.05, false, -Lx / 2 - 0.05, Lx / 2 + 0.05, false,
                 -Lz - 0.02, 0.3, false);
    double hb = 2.0 * S.h_init;

    // target half-space: x,y in [-Lx/2,Lx/2], z in [-Lz,0]; walls frozen
    double rho_b = Material::basalt().rho0, m_b = rho_b * dx * dx * dx;
    int nxy = (int)std::round(Lx / dx), nz = (int)std::round(Lz / dx);
    for (int i = 0; i < nxy; i++)
        for (int j = 0; j < nxy; j++)
            for (int k = 0; k < nz; k++) {
                double x = -Lx / 2 + (i + 0.5) * dx, y = -Lx / 2 + (j + 0.5) * dx;
                double z = -Lz + (k + 0.5) * dx;
                bool wall = (z < -Lz + hb) || std::abs(x) > Lx / 2 - hb || std::abs(y) > Lx / 2 - hb;
                S.cur_mat = 0; S.add({x, y, z}, {0, 0, 0}, m_b, 0.0, wall);
            }
    // iron impactor sphere, just above the surface, moving down
    double rho_i = Material::iron().rho0, m_i = rho_i * dx * dx * dx;
    double z0 = a + 0.5 * dx;
    int na = (int)std::ceil(a / dx) + 1;
    S.cur_mat = 1;
    for (int i = -na; i <= na; i++)
        for (int j = -na; j <= na; j++)
            for (int k = -na; k <= na; k++) {
                double x = i * dx, y = j * dx, z = z0 + k * dx;
                if (x * x + y * y + (z - z0) * (z - z0) <= a * a)
                    S.add({x, y, z}, {0, 0, -U}, m_i, 0.0);
            }
    int n_imp = 0; for (int i = 0; i < S.n; i++) if (S.mat[i] == 1) n_imp++;
    printf("pi_scaling: U=%.0f m/s  dx=%.4f (CPPR=%.1f)  N=%d (target+%d impactor)\n",
           U, dx, a / dx, S.n, n_imp);
    S.init();

    // run; track deepest crater floor (min z of target near axis). The crater
    // forms within ~2-3e-5 s and (no gravity, strength regime) then arrests.
    double t_end = 8.0e-5, t = 0;
    int nstep = 0;
    while (t < t_end) {
        double dt = S.timestep();
        if (t + dt > t_end) dt = t_end - t;
        S.step(dt); t += dt; nstep++;
        if (nstep % 100 == 0) {
            double zf = 0;
            for (int i = 0; i < S.n; i++)
                if (S.mat[i] == 0 && std::hypot(S.pos[i].x, S.pos[i].y) < 2 * a)
                    zf = std::min(zf, S.pos[i].z);
            printf("  step %d  t=%.2e  floor z=%.3f\n", nstep, t, zf);
        }
    }

    // radial surface profile of the TARGET: crater radius = outermost radius
    // whose surface is depressed below -0.25 a; depth = -min(surface z).
    int nb = (int)(Lx / 2 / dx) + 1;
    std::vector<double> surf(nb, -1e9);
    for (int i = 0; i < S.n; i++) {
        if (S.mat[i] != 0) continue;
        double r = std::hypot(S.pos[i].x, S.pos[i].y);
        int b = (int)(r / dx);
        if (b < nb) surf[b] = std::max(surf[b], S.pos[i].z);
    }
    double Rc = 0, depth = 0;
    for (int b = 0; b < nb; b++) {
        if (surf[b] < -1e8) continue;
        depth = std::max(depth, -surf[b]);
        if (surf[b] < -0.25 * a) Rc = (b + 0.5) * dx;
    }
    printf("RESULT U=%.0f Rc=%.4f depth=%.4f R/a=%.2f\n", U, Rc, depth, Rc / a);
    return 0;
}
