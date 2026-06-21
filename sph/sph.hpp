// 3D SPH core (header-only): cubic-spline kernel, EOS, cell-list neighbours
// (periodic in y,z), Monaghan artificial viscosity, KDK leapfrog.
//
// Standard weakly/locally-conservative SPH (Monaghan 1992):
//   rho_i = sum_j m_j W_ij
//   dv_i/dt = -sum_j m_j (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) grad_i W_ij
//   du_i/dt = 0.5 sum_j m_j (P_i/rho_i^2 + P_j/rho_j^2 + Pi_ij) (v_i-v_j).grad_i W_ij
//
// This version uses a single global smoothing length h (adaptive-h is the next
// step, needed for the large density contrasts in impacts). EOS is pluggable;
// IdealGas is used for the Sod validation, Tillotson for impacts (later).
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "vec3.hpp"

// ---------------- kernel (cubic spline, 3D) ----------------
struct Kernel {
    double h, sigma;
    explicit Kernel(double h_) : h(h_), sigma(1.0 / (M_PI * h_ * h_ * h_)) {}
    double W(double r) const {
        double q = r / h;
        if (q < 1.0) return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
        if (q < 2.0) { double t = 2.0 - q; return sigma * 0.25 * t * t * t; }
        return 0.0;
    }
    // dW/dr (scalar); grad_i W_ij = dWdr(r) * (r_ij / r)
    double dWdr(double r) const {
        double q = r / h;
        if (q < 1.0) return sigma * (-3.0 * q + 2.25 * q * q) / h;
        if (q < 2.0) { double t = 2.0 - q; return sigma * (-0.75 * t * t) / h; }
        return 0.0;
    }
};

// ---------------- EOS ----------------
struct IdealGas {
    double gamma = 1.4;
    double pressure(double rho, double u) const { return (gamma - 1.0) * rho * u; }
    double sound_speed(double rho, double u) const {
        double P = pressure(rho, u);
        return std::sqrt(gamma * P / std::max(rho, 1e-30));
    }
};

// ---------------- particle system ----------------
struct System {
    std::vector<Vec3> pos, vel, acc;
    std::vector<double> mass, rho, u, dudt, P, cs;
    std::vector<char> fixed;          // frozen boundary particles
    int n = 0;

    double h = 0.01;                  // global smoothing length
    double alpha = 1.0, beta = 2.0, eps = 0.01;   // artificial viscosity
    double alpha_u = 0.0;             // artificial thermal conductivity (Price 2008)
    double cfl = 0.2;
    IdealGas eos;

    // domain: x non-periodic; y,z periodic over [0,Ly],[0,Lz]
    double xmin = 0, xmax = 0, Ly = 0, Lz = 0;

    void add(Vec3 p, Vec3 v, double m, double uu, bool fix = false) {
        pos.push_back(p); vel.push_back(v); acc.push_back({});
        mass.push_back(m); rho.push_back(0); u.push_back(uu); dudt.push_back(0);
        P.push_back(0); cs.push_back(0); fixed.push_back(fix ? 1 : 0); n++;
    }

    // minimum-image separation r_i - r_j (periodic in y,z only)
    Vec3 sep(const Vec3& a, const Vec3& b) const {
        Vec3 d = a - b;
        if (Ly > 0) { d.y -= Ly * std::round(d.y / Ly); }
        if (Lz > 0) { d.z -= Lz * std::round(d.z / Lz); }
        return d;
    }

    // ---- cell list (size = 2h) ----
    int ncx = 0, ncy = 0, ncz = 0;
    double cell = 0;
    std::vector<std::vector<int>> cells;

    int cidx(int ix, int iy, int iz) const {
        // wrap y,z (periodic); clamp x
        ix = std::min(std::max(ix, 0), ncx - 1);
        iy = ((iy % ncy) + ncy) % ncy;
        iz = ((iz % ncz) + ncz) % ncz;
        return (ix * ncy + iy) * ncz + iz;
    }
    void build_cells() {
        cell = 2.0 * h;
        ncx = std::max(1, (int)std::ceil((xmax - xmin) / cell));
        ncy = std::max(1, (int)std::floor(Ly / cell));
        ncz = std::max(1, (int)std::floor(Lz / cell));
        cells.assign(ncx * ncy * ncz, {});
        for (int i = 0; i < n; i++) {
            int ix = (int)((pos[i].x - xmin) / cell);
            int iy = (int)(pos[i].y / cell);
            int iz = (int)(pos[i].z / cell);
            cells[cidx(ix, iy, iz)].push_back(i);
        }
    }
    template <class F>
    void for_neighbors(int i, F fn) const {
        int ix = (int)((pos[i].x - xmin) / cell);
        int iy = (int)(pos[i].y / cell);
        int iz = (int)(pos[i].z / cell);
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                    for (int j : cells[cidx(ix + dx, iy + dy, iz + dz)]) fn(j);
    }

    void compute_density() {
        Kernel k(h);
        for (int i = 0; i < n; i++) {
            double s = 0;
            for_neighbors(i, [&](int j) {
                double r = norm(sep(pos[i], pos[j]));
                if (r < 2.0 * h) s += mass[j] * k.W(r);
            });
            rho[i] = s;
        }
    }

    void compute_forces() {
        Kernel k(h);
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
                if (r >= 2.0 * h || r == 0) return;
                Vec3 gradW = rij * (k.dWdr(r) / r);
                Vec3 vij = vel[i] - vel[j];
                double rbar = 0.5 * (rho[i] + rho[j]);
                double Pij = 0;
                double dvr = dot(vij, rij);
                if (dvr < 0) {
                    double mu = h * dvr / (r * r + eps * h * h);
                    double cbar = 0.5 * (cs[i] + cs[j]);
                    Pij = (-alpha * cbar * mu + beta * mu * mu) / rbar;
                }
                double term = Pi_over + P[j] / (rho[j] * rho[j]) + Pij;
                a -= gradW * (mass[j] * term);
                du += 0.5 * mass[j] * term * dot(vij, gradW);
                // artificial thermal conductivity (Price 2008): diffuses u to
                // remove the contact-discontinuity "pressure blip".
                if (alpha_u > 0) {
                    double vsig_u = std::sqrt(std::abs(P[i] - P[j]) / rbar);
                    du += alpha_u * mass[j] * vsig_u * (u[i] - u[j]) * k.dWdr(r) / rbar;
                }
            });
            acc[i] = a; dudt[i] = du;
        }
    }

    double timestep() {
        double cmax = 1e-30;
        for (int i = 0; i < n; i++) cmax = std::max(cmax, cs[i] + alpha * cs[i]);
        return cfl * h / cmax;
    }

    // KDK leapfrog with trapezoidal energy update; fixed particles are static
    // walls (contribute P, don't move).
    void step(double dt) {
        std::vector<double> dudt_old = dudt;      // from the previous force eval
        for (int i = 0; i < n; i++) if (!fixed[i]) vel[i] += acc[i] * (0.5 * dt);
        for (int i = 0; i < n; i++) if (!fixed[i]) pos[i] += vel[i] * dt;
        build_cells();
        compute_density();
        compute_forces();                         // new acc, dudt
        for (int i = 0; i < n; i++) if (!fixed[i]) {
            vel[i] += acc[i] * (0.5 * dt);
            u[i] += 0.5 * dt * (dudt_old[i] + dudt[i]);
            if (u[i] < 1e-12) u[i] = 1e-12;
        }
    }

    void init() { build_cells(); compute_density(); compute_forces(); }
};
