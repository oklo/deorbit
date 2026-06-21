// Oblique Cayambe impact: 1 km iron monolith grazing a basalt half-space at the
// Phase-1 terminal conditions (v ~ 7.4 km/s, gamma ~ 2 deg). Tests the Phase-2a
// prediction: a grazing plough/ricochet + elongated gouge + leading-contact melt
// (NOT a vertical crater). Full 3D SPH with Tillotson EOS + strength.
//
//   ./cayambe <dx_m> <t_end_s>
// Prints a periodic summary (gouge depth, impactor centroid/speed/ricochet, melt)
// and dumps a downsampled snapshot at the end.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <chrono>
#include "sph.hpp"

static void snapshot(System& S, const char* fn) {
    FILE* f = std::fopen(fn, "w");
    std::fprintf(f, "x,y,z,vx,vy,vz,mat,u\n");
    for (int i = 0; i < S.n; i += 8)
        std::fprintf(f, "%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d,%.3e\n",
                     S.pos[i].x, S.pos[i].y, S.pos[i].z, S.vel[i].x, S.vel[i].y, S.vel[i].z,
                     S.mat[i], S.u[i]);
    std::fclose(f);
}

int main(int argc, char** argv) {
    double dx = (argc > 1 ? atof(argv[1]) : 50.0);
    double t_end = (argc > 2 ? atof(argv[2]) : 0.3);
    double walltime = (argc > 3 ? atof(argv[3]) : 1e30);   // seconds, hard cap
    double U = 7400.0, gam = 2.0 * M_PI / 180.0;     // terminal speed, grazing angle
    double vh = U * std::cos(gam), vn = U * std::sin(gam);
    double a = 500.0;                                // impactor radius (CPPR = a/dx)

    // domain (m): 10 km downrange x 2 km wide x 1.5 km deep rock (+headroom up top)
    double xlo = -1500, xhi = 8500, ylo = -1000, yhi = 1000, zlo = -1500, ztop = 0;
    System S;
    S.materials.push_back(Material::basalt());   // 0 = target rock
    S.materials.push_back(Material::iron());      // 1 = impactor
    S.eta = 1.3; S.h_init = S.eta * dx;
    S.set_domain(xlo, xhi, false, ylo, yhi, false, zlo, 600.0, false);
    double hb = 2.0 * S.h_init;

    double m_b = Material::basalt().rho0 * dx * dx * dx;
    int nx = (int)std::round((xhi - xlo) / dx), ny = (int)std::round((yhi - ylo) / dx),
        nz = (int)std::round((ztop - zlo) / dx);
    S.cur_mat = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                double x = xlo + (i + 0.5) * dx, y = ylo + (j + 0.5) * dx, z = zlo + (k + 0.5) * dx;
                bool wall = z < zlo + hb || std::abs(y) > yhi - hb || x < xlo + hb;
                S.add({x, y, z}, {0, 0, 0}, m_b, 0.0, wall);
            }
    int n_t = S.n;
    // iron impactor: sphere bottom touching z=0 at x=0, moving downrange + down
    double m_i = Material::iron().rho0 * dx * dx * dx, z0 = a + 0.5 * dx;
    int na = (int)std::ceil(a / dx) + 1;
    S.cur_mat = 1;
    for (int i = -na; i <= na; i++)
        for (int j = -na; j <= na; j++)
            for (int k = -na; k <= na; k++) {
                double x = i * dx, y = j * dx, z = z0 + k * dx;
                if (x * x + y * y + (z - z0) * (z - z0) <= a * a)
                    S.add({x, y, z}, {vh, 0, -vn}, m_i, 0.0);
            }
    int n_imp = S.n - n_t;
    printf("cayambe: U=%.0f m/s gamma=2deg  dx=%.1f (CPPR=%.1f)  N=%d (%d target + %d iron)\n",
           U, dx, a / dx, S.n, n_t, n_imp);
    S.init();

    double t = 0; int nstep = 0, M = 25, isnap = 0;
    double u_melt = 1.0e6;                            // ~ basalt/iron melt energy
    double next_snap = 0.1;                           // snapshot every 0.1 s sim time
    auto wall0 = std::chrono::steady_clock::now();
    const char* stop = "t_end";
    while (t < t_end) {
        double dt = S.timestep();
        if (t + dt > t_end) dt = t_end - t;
        S.step(dt); t += dt; nstep++;
        if (nstep % M == 0 || t >= t_end) {
            double cx = 0, cz = 0, vz = 0, sp = 0, pen = 0, umax = 0; int ni = 0, nmelt = 0;
            for (int i = 0; i < S.n; i++) {
                if (S.mat[i] != 1) continue;
                cx += S.pos[i].x; cz += S.pos[i].z; vz += S.vel[i].z; sp += norm(S.vel[i]); ni++;
                pen = std::min(pen, S.pos[i].z);          // impactor penetration (min z)
                if (S.u[i] > u_melt) nmelt++;
                umax = std::max(umax, S.u[i]);
            }
            double wsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
            printf("  t=%.4f step=%d  iron: x=%.0f z=%.0f |v|=%.0f vz=%.0f  pen=%.0f  "
                   "melt=%.0f%% umax=%.1e  [%.0fs wall]\n",
                   t, nstep, cx / ni, cz / ni, sp / ni, vz / ni, pen,
                   100.0 * nmelt / ni, umax, wsec);
            fflush(stdout);
        }
        if (t >= next_snap) {
            char fn[64]; std::snprintf(fn, 64, "cayambe_snap_%03d.csv", isnap++);
            snapshot(S, fn); next_snap += 0.1;
        }
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count() > walltime) {
            stop = "walltime"; break;
        }
    }
    snapshot(S, "cayambe_snap.csv");
    printf("done (%s): %d steps to t=%.4f; wrote cayambe_snap.csv (+%d checkpoints)\n",
           stop, nstep, t, isnap);
    return 0;
}
