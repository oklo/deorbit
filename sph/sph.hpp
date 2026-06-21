// 3D SPH core (header-only) with ADAPTIVE smoothing length.
//
//   rho_i = sum_j m_j W(r_ij, h_i),   h_i = eta (m_i/rho_i)^(1/3)   [iterated]
//   dv_i/dt = -sum_j m_j (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) grad W_ij(hbar)
//   du_i/dt =  0.5 sum_j m_j (...) (v_i-v_j).grad W_ij(hbar)  + artificial conduction
//
// Density uses gather with the particle's own h_i; forces use a symmetric kernel
// at hbar_ij = 0.5(h_i+h_j) so momentum is conserved exactly. (grad-h correction
// terms are omitted -- a possible later refinement.) Cell list sized to 2*h_max,
// per-dimension, periodic in y,z.
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "vec3.hpp"

// cubic-spline kernel (3D) as free functions of (r,h)
inline double kW(double r, double h) {
    double q = r / h, s = 1.0 / (M_PI * h * h * h);
    if (q < 1.0) return s * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    if (q < 2.0) { double t = 2.0 - q; return s * 0.25 * t * t * t; }
    return 0.0;
}
inline double kdWdr(double r, double h) {   // dW/dr; grad_i W = kdWdr * (r_ij/r)
    double q = r / h, s = 1.0 / (M_PI * h * h * h);
    if (q < 1.0) return s * (-3.0 * q + 2.25 * q * q) / h;
    if (q < 2.0) { double t = 2.0 - q; return s * (-0.75 * t * t) / h; }
    return 0.0;
}

struct IdealGas {
    double gamma = 1.4;
    double pressure(double rho, double u) const { return (gamma - 1.0) * rho * u; }
    double sound_speed(double rho, double u) const {
        return std::sqrt(gamma * pressure(rho, u) / std::max(rho, 1e-30));
    }
};

struct System {
    std::vector<Vec3> pos, vel, acc;
    std::vector<double> mass, rho, u, dudt, P, cs, hh;   // hh = per-particle h
    std::vector<char> fixed;
    int n = 0;

    double eta = 1.3;                 // h = eta (m/rho)^(1/3); ~ neighbour count
    double h_init = 0.01;             // initial guess, seeded into hh on add()
    double h_max = 0.01;
    double alpha = 1.0, beta = 2.0, eps = 0.01;   // artificial viscosity
    double alpha_u = 0.0;             // artificial thermal conductivity
    double cfl = 0.2;
    IdealGas eos;
    double xmin = 0, xmax = 0, Ly = 0, Lz = 0;    // x non-periodic; y,z periodic

    void add(Vec3 p, Vec3 v, double m, double uu, bool fix = false) {
        pos.push_back(p); vel.push_back(v); acc.push_back({});
        mass.push_back(m); rho.push_back(0); u.push_back(uu); dudt.push_back(0);
        P.push_back(0); cs.push_back(0); hh.push_back(h_init); fixed.push_back(fix ? 1 : 0);
        n++;
    }

    Vec3 sep(const Vec3& a, const Vec3& b) const {
        Vec3 d = a - b;
        if (Ly > 0) d.y -= Ly * std::round(d.y / Ly);
        if (Lz > 0) d.z -= Lz * std::round(d.z / Lz);
        return d;
    }

    // ---- per-dimension cell list (size >= 2*h_max) ----
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
        // periodic dims need >=3 cells for the +/-1 stencil; else collapse to 1
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

    // iterate rho <-> h to consistency
    void compute_density_h(int iters = 3) {
        for (int it = 0; it < iters; it++) {
            build_cells();
            for (int i = 0; i < n; i++) {
                double s = 0, hi = hh[i];
                for_neighbors(i, [&](int j) {
                    double r = norm(sep(pos[i], pos[j]));
                    if (r < 2.0 * hi) s += mass[j] * kW(r, hi);
                });
                rho[i] = std::max(s, 1e-30);
                hh[i] = eta * std::cbrt(mass[i] / rho[i]);
            }
        }
        build_cells();   // final cell list consistent with converged h
    }

    void compute_forces() {
        for (int i = 0; i < n; i++) {
            P[i] = eos.pressure(rho[i], u[i]);
            cs[i] = eos.sound_speed(rho[i], u[i]);
        }
        for (int i = 0; i < n; i++) {
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
        }
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
