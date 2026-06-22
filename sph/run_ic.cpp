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
#include <string>
#include <vector>
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
    // optional --frac flag (anywhere): enable Monaghan artificial stress +
    // Benz-Asphaug brittle damage. Filter it out, then parse positionally.
    bool frac = false;
    std::vector<char*> av;
    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "--frac") frac = true; else av.push_back(argv[i]);
    }
    argc = (int)av.size(); argv = av.data();
    bool resume = (argc > 1 && std::string(argv[1]) == "--resume");
    System S;
    S.materials.push_back(Material::basalt());   // 0
    S.materials.push_back(Material::iron());      // 1
    double dx, t_end, walltime, t = 0;
    int isnap = 0;
    auto wall0 = std::chrono::steady_clock::now();
    auto secs = [&]{ return std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count(); };

    if (resume) {
        // run_ic --resume ic.bin snap.bin dx t0 t_end [walltime]
        // continue a run from a snapshot: domain + frozen flags from the IC (same
        // particle order), current pos/vel/u from the snapshot. NB deviatoric
        // stress S is NOT stored in snapshots -> reset to 0 (fine once the body
        // is largely fragmented/airborne).
        if (argc < 7) { std::fprintf(stderr, "usage: run_ic --resume ic.bin snap.bin dx t0 t_end [walltime]\n"); return 1; }
        const char* icf = argv[2]; const char* snapf = argv[3];
        dx = atof(argv[4]); double t0 = atof(argv[5]); t_end = atof(argv[6]);
        walltime = (argc > 7 ? atof(argv[7]) : 1e30);
        S.eta = 1.3; S.h_init = S.eta * dx;
        FILE* fi = std::fopen(icf, "rb");
        if (!fi) { std::fprintf(stderr, "cannot open %s\n", icf); return 1; }
        double hdr[7]; std::fread(hdr, sizeof(double), 7, fi);
        long Nic = (long)std::llround(hdr[6]);
        S.set_domain(hdr[0], hdr[1], false, hdr[2], hdr[3], false, hdr[4], hdr[5], false);
        std::vector<char> fixedIC(Nic);
        for (long i = 0; i < Nic; i++) { double r[10]; std::fread(r, sizeof(double), 10, fi); fixedIC[i] = (r[9] != 0.0); }
        std::fclose(fi);
        FILE* fs = std::fopen(snapf, "rb");
        if (!fs) { std::fprintf(stderr, "cannot open %s\n", snapf); return 1; }
        long Ns; std::fread(&Ns, sizeof(long), 1, fs);
        int n_iron = 0;
        for (long i = 0; i < Ns; i++) {
            double r[10]; std::fread(r, sizeof(double), 10, fs);
            int mat = (int)r[9]; double mass = (mat == 1 ? 7800.0 : 2700.0) * dx * dx * dx;
            S.cur_mat = mat;
            S.add({r[0], r[1], r[2]}, {r[3], r[4], r[5]}, mass, r[6], (i < Nic && fixedIC[i]));
            if (mat == 1) n_iron++;
        }
        std::fclose(fs);
        t = t0; isnap = (int)std::llround(t0 / 0.05);
        printf("RESUMED from %s at t=%.3f: N=%d (%d iron)  [deviatoric S reset to 0]\n", snapf, t0, S.n, n_iron);
    } else {
        if (argc < 4) { std::fprintf(stderr, "usage: run_ic ic.bin dx t_end [walltime]\n"); return 1; }
        const char* icf = argv[1];
        dx = atof(argv[2]); t_end = atof(argv[3]);
        walltime = (argc > 4 ? atof(argv[4]) : 1e30);
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
    }
    if (frac) {
        S.eps_as = 0.3; S.damage_on = true; S.seed_damage(dx);
        printf("FRAC ON: Monaghan artificial stress (eps=0.3) + Benz-Asphaug brittle damage\n");
    }
    fprintf(stderr, "[diag] loaded, starting init...\n"); fflush(stderr);
    S.init();
    { double hlo = 1e30, hhi = 0; for (int i = 0; i < S.n; i++) { hlo = std::min(hlo, S.hh[i]); hhi = std::max(hhi, S.hh[i]); }
      fprintf(stderr, "[diag] init done in %.1fs; h range %.1f..%.1f m; cells nc=%dx%dx%d\n",
              secs(), hlo, hhi, S.nc[0], S.nc[1], S.nc[2]); fflush(stderr); }

    double next_snap = (std::floor(t / 0.05) + 1) * 0.05, u_melt = 1.0e6;
    int nstep = 0;
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
            char fn[64]; std::snprintf(fn, 64, "%s%03d.bin", frac ? "frac_snap_" : "ic_snap_", isnap++);
            snapshot(S, fn); next_snap += 0.05;
        }
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count() > walltime) {
            stop = "walltime"; break;
        }
    }
    snapshot(S, frac ? "frac_snap_final.bin" : "ic_snap_final.bin");
    printf("done (%s): %d steps to t=%.4f (+%d checkpoints)\n", stop, nstep, t, isnap);
    return 0;
}