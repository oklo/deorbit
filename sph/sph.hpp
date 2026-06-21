// 3D SPH core (header-only): adaptive-h, multi-material EOS (eos.hpp),
// Monaghan artificial viscosity + Price conductivity, KDK leapfrog.
// Parallelised over particles with std::thread (density + forces).
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>
#include "vec3.hpp"
#include "eos.hpp"

inline double kW(double r, double h) {
    double q = r / h, s = 1.0 / (M_PI * h * h * h);
    if (q < 1.0) return s * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    if (q < 2.0) { double t = 2.0 - q; return s * 0.25 * t * t * t; }
    return 0.0;
}
inline double kdWdr(double r, double h) {
    double q = r / h, s = 1.0 / (M_PI * h * h * h);
    if (q < 1.0) return s * (-3.0 * q + 2.25 * q * q) / h;
    if (q < 2.0) { double t = 2.0 - q; return s * (-0.75 * t * t) / h; }
    return 0.0;
}

template <class F>
static void parallel_for(int n, F f) {
    unsigned nt = std::thread::hardware_concurrency();
    if (nt == 0) nt = 1;
    if (nt > 16) nt = 16;
    if (n < 2000 || nt == 1) { for (int i = 0; i < n; i++) f(i); return; }
    std::vector<std::thread> ts;
    int chunk = (n + (int)nt - 1) / (int)nt;
    for (unsigned t = 0; t < nt; t++) {
        int lo = (int)t * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        ts.emplace_back([&f, lo, hi] { for (int i = lo; i < hi; i++) f(i); });
    }
    for (auto& t : ts) t.join();
}

struct System {
    std::vector<Vec3> pos, vel, acc;
    std::vector<double> mass, rho, u, dudt, P, cs, hh;
    std::vector<int> mat;
    std::vector<char> fixed;
    std::vector<Material> materials;
    int cur_mat = 0;
    int n = 0;

    double eta = 1.3;
    double h_init = 0.01, h_max = 0.01;
    double alpha = 1.0, beta = 2.0, eps = 0.01, alpha_u = 0.0;
    double cfl = 0.2;
    double xmin = 0, xmax = 0, Ly = 0, Lz = 0;

    void add(Vec3 p, Vec3 v, double m, double uu, bool fix = false) {
        pos.push_back(p); vel.push_back(v); acc.push_back({});
        mass.push_back(m); rho.push_back(0); u.push_back(uu); dudt.push_back(0);
        P.push_back(0); cs.push_back(0); hh.push_back(h_init);
        mat.push_back(cur_mat); fixed.push_back(fix ? 1 : 0); n++;
    }

    Vec3 sep(const Vec3& a, const Vec3& b) const {
        Vec3 d = a - b;
        if (Ly > 0) d.y -= Ly * std::round(d.y / Ly);
        if (Lz > 0) d.z -= Lz * std::round(d.z / Lz);
        return d;
    }

    int ncx = 1, ncy = 1, ncz = 1;
    double cx = 1, cy = 1, cz = 1;
    std::vector<std::vector<int>> cells;

    int cidx(int ix, int iy, int iz) const {
        ix = std::min(std::max(ix, 0), ncx - 1);
        iy = ((iy % ncy) + ncy) % ncy;
        iz = ((iz % ncz) + ncz) % ncz;
        return (ix * ncy + iy) * ncz + iz;
    }
    void build_cells() {
        h_max = 1e-30;
        for (int i = 0; i < n; i++) h_max = std::max(h_max, hh[i]);
        double rc = 2.0 * h_max;
        ncx = std::max(1, (int)((xmax - xmin) / rc));
        ncy = (int)(Ly / rc); ncz = (int)(Lz / rc);
        if (ncy < 3) ncy = 1;
        if (ncz < 3) ncz = 1;
        cx = (xmax - xmin) / ncx; cy = Ly / std::max(1, ncy); cz = Lz / std::max(1, ncz);
        cells.assign((size_t)ncx * ncy * ncz, {});
        for (int i = 0; i < n; i++) {
            int ix = (int)((pos[i].x - xmin) / cx);
            int iy = ncy > 1 ? (int)(pos[i].y / cy) : 0;
            int iz = ncz > 1 ? (int)(pos[i].z / cz) : 0;
            cells[cidx(ix, iy, iz)].push_back(i);
        }
    }
    template <class F>
    void for_neighbors(int i, F fn) const {
        int ix = (int)((pos[i].x - xmin) / cx);
        int iy = ncy > 1 ? (int)(pos[i].y / cy) : 0;
        int iz = ncz > 1 ? (int)(pos[i].z / cz) : 0;
        int dyl = ncy > 1 ? -1 : 0, dyh = ncy > 1 ? 1 : 0;
        int dzl = ncz > 1 ? -1 : 0, dzh = ncz > 1 ? 1 : 0;
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = dyl; dy <= dyh; dy++)
                for (int dz = dzl; dz <= dzh; dz++)
                    for (int j : cells[cidx(ix + dx, iy + dy, iz + dz)]) fn(j);
    }

    // adaptive density: rebuild cells, iterate rho<->h; if h grew past the cell
    // size, rebuild bigger and repeat. Density gather (uses h_i only) is parallel.
    void compute_density_h(int iters = 2) {
        for (int round = 0; round < 4; round++) {
            build_cells();
            double hmax0 = h_max;
            for (int it = 0; it < iters; it++) {
                parallel_for(n, [&](int i) {
                    double s = 0, hi = hh[i];
                    for_neighbors(i, [&](int j) {
                        double r = norm(sep(pos[i], pos[j]));
                        if (r < 2.0 * hi) s += mass[j] * kW(r, hi);
                    });
                    rho[i] = std::max(s, 1e-30);
                    hh[i] = eta * std::cbrt(mass[i] / rho[i]);
                });
            }
            double hm = 1e-30;
            for (int i = 0; i < n; i++) hm = std::max(hm, hh[i]);
            if (hm <= hmax0 * 1.02) break;       // cell size covered the support
        }
        build_cells();
    }

    void compute_forces() {
        parallel_for(n, [&](int i) {
            P[i] = materials[mat[i]].pressure(rho[i], u[i]);
            cs[i] = materials[mat[i]].sound_speed(rho[i], u[i]);
        });
        parallel_for(n, [&](int i) {
            Vec3 a{}; double du = 0;
            double Pi_over = P[i] / (rho[i] * rho[i]);
            for_neighbors(i, [&](int j) {
                if (j == i) return;
                Vec3 rij = sep(pos[i], pos[j]);
                double r = norm(rij);
                double hbar = 0.5 * (hh[i] + hh[j]);
                if (r >= 2.0 * hbar || r == 0) return;
                Vec3 gradW = rij * (kdWdr(r, hbar) / r);
                Vec3 vij = vel[i] - vel[j];
                double rbar = 0.5 * (rho[i] + rho[j]);
                double Pij = 0, dvr = dot(vij, rij);
                if (dvr < 0) {
                    double mu = hbar * dvr / (r * r + eps * hbar * hbar);
                    double cbar = 0.5 * (cs[i] + cs[j]);
                    Pij = (-alpha * cbar * mu + beta * mu * mu) / rbar;
                }
                double term = Pi_over + P[j] / (rho[j] * rho[j]) + Pij;
                a -= gradW * (mass[j] * term);
                du += 0.5 * mass[j] * term * dot(vij, gradW);
                if (alpha_u > 0) {
                    double vsig = std::sqrt(std::abs(P[i] - P[j]) / rbar);
                    du += alpha_u * mass[j] * vsig * (u[i] - u[j]) * kdWdr(r, hbar) / rbar;
                }
            });
            acc[i] = a; dudt[i] = du;
        });
    }

    double timestep() {
        double dt = 1e30;
        for (int i = 0; i < n; i++)
            dt = std::min(dt, cfl * hh[i] / (cs[i] * (1.0 + alpha) + 1e-30));
        return dt;
    }

    void step(double dt) {
        std::vector<double> dudt_old = dudt;
        for (int i = 0; i < n; i++) if (!fixed[i]) vel[i] += acc[i] * (0.5 * dt);
        for (int i = 0; i < n; i++) if (!fixed[i]) pos[i] += vel[i] * dt;
        compute_density_h();
        compute_forces();
        for (int i = 0; i < n; i++) if (!fixed[i]) {
            vel[i] += acc[i] * (0.5 * dt);
            u[i] += 0.5 * dt * (dudt_old[i] + dudt[i]);
            if (u[i] < 1e-12) u[i] = 1e-12;
        }
    }

    void init() { compute_density_h(); compute_forces(); }
};
