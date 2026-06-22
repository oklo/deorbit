// Generic SPH runner: loads a binary particle IC (e.g. from make_cayambe_ic.py)
// and integrates it with Tillotson EOS + strength, walltime cap + checkpoints.
//
//   ./run_ic <ic.bin> <dx_m> <t_end_s> [walltime_s]
// Binary IC (little-endian f8): xlo xhi ylo yhi zlo zhi n, then n*[x y z vx vy
// vz mass u mat fixed].  mat: 0=basalt, 1=iron.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include "sph.hpp"

static void snapshot(System& S, const char* fn) {
    // FULL-resolution binary dump (every particle) so the kernel field can be
    // reconstructed at the simulation's true resolution h.
    // Format: int64 n, then n*[x y z vx vy vz u rho h mat] as float64.
    FILE* f = std::fopen(fn, "wb");
    long n = S.n; std::fwrite(&n, sizeof(long), 1, f);
    std::vector<double> buf((size_t)n * 10);
    for (int i = 0; i < S.n; i++) {
        double* r = &buf[(size_t)i * 10];
        r[0] = S.pos[i].x; r[1] = S.pos[i].y; r[2] = S.pos[i].z;
        r[3] = S.vel[i].x; r[4] = S.vel[i].y; r[5] = S.vel[i].z;
        r[6] = S.u[i]; r[7] = S.rho[i]; r[8] = S.hh[i]; r[9] = (double)S.mat[i];
    }
    std::fwrite(buf.data(), sizeof(double), buf.size(), f);
    std::fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: run_ic ic.bin dx t_end [walltime]\n"); return 1; }
    const char* icf = argv[1];
    double dx = atof(argv[2]), t_end = atof(argv[3]);
    double walltime = (argc > 4 ? atof(argv[4]) : 1e30);

    System S;
    S.materials.push_back(Material::basalt());   // 0
    S.materials.push_back(Material::iron());      // 1
    S.eta = 1.3; S.h_init = S.eta * dx;

    FILE* f = std::fopen(icf, "rb");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", icf); return 1; }
    double hdr[7];
    if (std::fread(hdr, sizeof(double), 7, f) != 7) { std::fprintf(stderr, "bad header\n"); return 1; }
    long N = (long)std::llround(hdr[6]);
    S.set_domain(hdr[0], hdr[1], false, hdr[2], hdr[3], false, hdr[4], hdr[5], false);
    int n_iron = 0;
    for (long i = 0; i < N; i++) {
        double r[10];
        if (std::fread(r, sizeof(double), 10, f) != 10) { std::fprintf(stderr, "short read @%ld\n", i); return 1; }
        S.cur_mat = (int)r[8];
        S.add({r[0], r[1], r[2]}, {r[3], r[4], r[5]}, r[6], r[7], r[9] != 0.0);
        if ((int)r[8] == 1) n_iron++;
    }
    std::fclose(f);
    printf("loaded %s: N=%d (%d iron)  domain x[%.0f,%.0f] y[%.0f,%.0f] z[%.0f,%.0f]  dx=%.1f\n",
           icf, S.n, n_iron, hdr[0], hdr[1], hdr[2], hdr[3], hdr[4], hdr[5], dx);
    auto wall0 = std::chrono::steady_clock::now();
    auto secs = [&]{ return std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count(); };
    fprintf(stderr, "[diag] loaded, starting init...\n"); fflush(stderr);
    S.init();
    { double hlo = 1e30, hhi = 0; for (int i = 0; i < S.n; i++) { hlo = std::min(hlo, S.hh[i]); hhi = std::max(hhi, S.hh[i]); }
      fprintf(stderr, "[diag] init done in %.1fs; h range %.1f..%.1f m; cells nc=%dx%dx%d\n",
              secs(), hlo, hhi, S.nc[0], S.nc[1], S.nc[2]); fflush(stderr); }

    double t = 0, next_snap = 0.05, u_melt = 1.0e6;
    int nstep = 0, isnap = 0;
    const char* stop = "t_end";
    while (t < t_end) {
        double dt = S.timestep();
        if (t + dt > t_end) dt = t_end - t;
        S.step(dt); t += dt; nstep++;
        if (nstep <= 3) { fprintf(stderr, "[diag] step %d dt=%.2e t=%.4f h_max=%.0f wall=%.1fs\n",
              nstep, dt, t, [&]{double m=0;for(int i=0;i<S.n;i++)m=std::max(m,S.hh[i]);return m;}(), secs()); fflush(stderr); }
        if (nstep % 25 == 0 || t >= t_end) {
            double cx = 0, cz = 0, vz = 0, sp = 0, pen = 0, umax = 0; int ni = 0, nmelt = 0;
            for (int i = 0; i < S.n; i++) {
                if (S.mat[i] != 1) continue;
                cx += S.pos[i].x; cz += S.pos[i].z; vz += S.vel[i].z; sp += norm(S.vel[i]); ni++;
                pen = std::min(pen, S.pos[i].z);
                if (S.u[i] > u_melt) nmelt++;
                umax = std::max(umax, S.u[i]);
            }
            double wsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
            printf("  t=%.4f step=%d  iron: x=%.0f z=%.0f |v|=%.0f vz=%.0f  pen_z=%.0f  "
                   "melt=%.0f%% umax=%.1e  [%.0fs]\n",
                   t, nstep, cx / ni, cz / ni, sp / ni, vz / ni, pen,
                   100.0 * nmelt / ni, umax, wsec);
            fflush(stdout);
        }
        if (t >= next_snap) {
            char fn[64]; std::snprintf(fn, 64, "ic_snap_%03d.bin", isnap++);
            snapshot(S, fn); next_snap += 0.05;
        }
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count() > walltime) {
            stop = "walltime"; break;
        }
    }
    snapshot(S, "ic_snap_final.bin");
    printf("done (%s): %d steps to t=%.4f (+%d checkpoints)\n", stop, nstep, t, isnap);
    return 0;
}