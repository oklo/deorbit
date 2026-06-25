// 3D SPH core (header-only): adaptive-h, multi-material EOS (eos.hpp), material
// STRENGTH (elastic-perfectly-plastic deviatoric stress, Jaumann rate, von Mises
// yield), Monaghan artificial viscosity + Price conductivity, KDK leapfrog.
// Parallelised over particles with std::thread.
//
// Total stress sigma_ab = -P delta_ab + S_ab (S = deviatoric). Strength refs:
// Libersky & Petschek 1991; Benz & Asphaug 1995; Gray, Monaghan & Swift 2001.
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>
#include <random>
#include "vec3.hpp"
#include "../common/eos.hpp"

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

// symmetric 3x3 (deviatoric stress / strain rate)
struct Sym3 {
    double xx = 0, yy = 0, zz = 0, xy = 0, xz = 0, yz = 0;
};
inline Vec3 Sdot(const Sym3& s, const Vec3& g) {   // S . g
    return {s.xx * g.x + s.xy * g.y + s.xz * g.z,
            s.xy * g.x + s.yy * g.y + s.yz * g.z,
            s.xz * g.x + s.yz * g.y + s.zz * g.z};
}
// largest eigenvalue of a symmetric 3x3 (Smith 1961 trig method) -- used for the
// max principal stress (damage) without needing eigenvectors.
inline double eig_max(const Sym3& s) {
    double p1 = s.xy * s.xy + s.xz * s.xz + s.yz * s.yz;
    double q = (s.xx + s.yy + s.zz) / 3.0;
    if (p1 < 1e-30) return std::max({s.xx, s.yy, s.zz});
    double p2 = (s.xx - q) * (s.xx - q) + (s.yy - q) * (s.yy - q) + (s.zz - q) * (s.zz - q) + 2 * p1;
    double p = std::sqrt(p2 / 6.0);
    // r = det((A-qI)/p)/2
    double a = (s.xx - q) / p, b = (s.yy - q) / p, c = (s.zz - q) / p;
    double d = s.xy / p, e = s.xz / p, f = s.yz / p;
    double r = (a * (b * c - f * f) - d * (d * c - f * e) + e * (d * f - b * e)) / 2.0;
    r = std::max(-1.0, std::min(1.0, r));
    double phi = std::acos(r) / 3.0;
    return q + 2.0 * p * std::cos(phi);            // largest eigenvalue
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
    std::vector<double> mass, rho, u, dudt, P, cs, csig, hh, eps_work;
    // fracture/tensile-stability: D=damage[0,1], eps_act=Weibull flaw strain,
    // sigmax=max principal (tensile) stress, Pf=damaged pressure, Rart=Monaghan
    // artificial-stress coeff, dscale=(1-D) deviatoric scale.
    std::vector<double> D, eps_act, sigmax, Pf, Rart, dscale;
    std::vector<Sym3> S, dSdt;
    std::vector<int> mat;
    std::vector<char> fixed;
    std::vector<Material> materials;
    int cur_mat = 0;
    int n = 0;

    double eta = 1.3;
    double h_init = 0.01, h_max = 0.01;
    double alpha = 1.0, beta = 2.0, eps = 0.01, alpha_u = 0.0;
    double cfl = 0.2;
    double eps_as = 0.0;        // Monaghan artificial-stress coeff (0 = off)
    int    as_n = 4;            // artificial-stress kernel exponent
    bool   damage_on = false;   // Benz-Asphaug brittle damage (off by default)
    // per-axis domain bounds [lo,hi] and periodicity (x,y,z)
    double dlo[3] = {0, 0, 0}, dhi[3] = {0, 0, 0};
    bool per[3] = {false, false, false};
    void set_domain(double xl, double xh, bool px, double yl, double yh, bool py,
                    double zl, double zh, bool pz) {
        dlo[0] = xl; dhi[0] = xh; per[0] = px;
        dlo[1] = yl; dhi[1] = yh; per[1] = py;
        dlo[2] = zl; dhi[2] = zh; per[2] = pz;
    }

    void add(Vec3 p, Vec3 v, double m, double uu, bool fix = false) {
        pos.push_back(p); vel.push_back(v); acc.push_back({});
        mass.push_back(m); rho.push_back(0); u.push_back(uu); dudt.push_back(0);
        P.push_back(0); cs.push_back(0); csig.push_back(0); hh.push_back(h_init);
        eps_work.push_back(0); S.push_back({}); dSdt.push_back({});
        D.push_back(0); eps_act.push_back(1e30); sigmax.push_back(0);
        Pf.push_back(0); Rart.push_back(0); dscale.push_back(1.0);
        mat.push_back(cur_mat); fixed.push_back(fix ? 1 : 0); n++;
    }

    // Benz-Asphaug (1995) flaw seeding: each brittle particle (material wk>0) gets
    // a first-flaw activation strain drawn from the Weibull distribution
    // n(<eps)=k*V*eps^m, i.e. eps_act = (rank/(k*V))^(1/m) with random rank. For
    // one flaw/particle and V_i=dx^3 this is (u/(k*dx^3))^(1/m), u~Uniform(0,1].
    void seed_damage(double dx, unsigned seed = 12345) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> U(1e-9, 1.0);
        double V = dx * dx * dx;
        for (int i = 0; i < n; i++) {
            const Material& m = materials[mat[i]];
            if (m.wk > 0 && m.wm > 0)
                eps_act[i] = std::pow(U(rng) / (m.wk * V), 1.0 / m.wm);
            else
                eps_act[i] = 1e30;     // ductile / non-brittle: never flaw-fractures
        }
    }

    Vec3 sep(const Vec3& a, const Vec3& b) const {
        Vec3 d = a - b;
        if (per[0]) { double L = dhi[0] - dlo[0]; d.x -= L * std::round(d.x / L); }
        if (per[1]) { double L = dhi[1] - dlo[1]; d.y -= L * std::round(d.y / L); }
        if (per[2]) { double L = dhi[2] - dlo[2]; d.z -= L * std::round(d.z / L); }
        return d;
    }

    // proper 3D spatial cell list (cell size 2*h_max), per-axis periodicity
    int nc[3] = {1, 1, 1};
    double cw[3] = {1, 1, 1};
    std::vector<std::vector<int>> cells;
    int binof(int ax, double c) const {
        int i = (int)((c - dlo[ax]) / cw[ax]);
        if (i < 0) i = 0; if (i >= nc[ax]) i = nc[ax] - 1;
        return i;
    }
    void build_cells() {
        h_max = 1e-30;
        for (int i = 0; i < n; i++) h_max = std::max(h_max, hh[i]);
        double rc = 2.0 * h_max;
        for (int a = 0; a < 3; a++) {
            double ext = dhi[a] - dlo[a];
            if (ext <= 0) ext = rc;
            int m = std::max(1, (int)(ext / rc));
            if (per[a] && m < 3) m = 1;        // periodic needs >=3 cells, else collapse
            nc[a] = m; cw[a] = ext / m;
        }
        cells.assign((size_t)nc[0] * nc[1] * nc[2], {});
        for (int i = 0; i < n; i++)
            cells[(binof(0, pos[i].x) * nc[1] + binof(1, pos[i].y)) * nc[2] + binof(2, pos[i].z)]
                .push_back(i);
    }
    template <class F>
    void for_neighbors(int i, F fn) const {
        int b0 = binof(0, pos[i].x), b1 = binof(1, pos[i].y), b2 = binof(2, pos[i].z);
        int lo0 = nc[0] > 1 ? -1 : 0, hi0 = nc[0] > 1 ? 1 : 0;
        int lo1 = nc[1] > 1 ? -1 : 0, hi1 = nc[1] > 1 ? 1 : 0;
        int lo2 = nc[2] > 1 ? -1 : 0, hi2 = nc[2] > 1 ? 1 : 0;
        for (int d0 = lo0; d0 <= hi0; d0++) {
            int j0 = b0 + d0;
            if (per[0]) j0 = ((j0 % nc[0]) + nc[0]) % nc[0]; else if (j0 < 0 || j0 >= nc[0]) continue;
            for (int d1 = lo1; d1 <= hi1; d1++) {
                int j1 = b1 + d1;
                if (per[1]) j1 = ((j1 % nc[1]) + nc[1]) % nc[1]; else if (j1 < 0 || j1 >= nc[1]) continue;
                for (int d2 = lo2; d2 <= hi2; d2++) {
                    int j2 = b2 + d2;
                    if (per[2]) j2 = ((j2 % nc[2]) + nc[2]) % nc[2]; else if (j2 < 0 || j2 >= nc[2]) continue;
                    for (int jj : cells[(j0 * nc[1] + j1) * nc[2] + j2]) fn(jj);
                }
            }
        }
    }

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
                    // clamp h to [0.5, 2]*h_init: keeps the cell size from being
                    // bloated by a few large-h free-surface particles (which made
                    // every particle scan ~7500 candidates -> ~74 s/step at 5.5M).
                    hh[i] = std::min(2.0 * h_init,
                                     std::max(0.5 * h_init, eta * std::cbrt(mass[i] / rho[i])));
                });
            }
            double hm = 1e-30;
            for (int i = 0; i < n; i++) hm = std::max(hm, hh[i]);
            if (hm <= hmax0 * 1.02) break;
        }
        build_cells();
    }

    // velocity gradient -> strain rate + spin -> hypoelastic dS/dt (Jaumann),
    // and the deviatoric power S:eps_dot/rho. Per particle, parallel.
    void compute_strain_stress() {
        parallel_for(n, [&](int i) {
            double G = materials[mat[i]].G;
            if (G <= 0) { dSdt[i] = Sym3{}; eps_work[i] = 0; return; }
            double L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            for_neighbors(i, [&](int j) {
                if (j == i) return;
                Vec3 rij = sep(pos[i], pos[j]);
                double r = norm(rij);
                double hbar = 0.5 * (hh[i] + hh[j]);
                if (r >= 2.0 * hbar || r == 0) return;
                Vec3 gW = rij * (kdWdr(r, hbar) / r);   // grad_i W
                double Vj = mass[j] / rho[j];
                Vec3 dv = vel[j] - vel[i];
                double d[3] = {dv.x, dv.y, dv.z}, g[3] = {gW.x, gW.y, gW.z};
                for (int a = 0; a < 3; a++)
                    for (int bb = 0; bb < 3; bb++) L[a][bb] += Vj * d[a] * g[bb];
            });
            // strain rate (sym) and spin (antisym)
            double exx = L[0][0], eyy = L[1][1], ezz = L[2][2];
            double exy = 0.5 * (L[0][1] + L[1][0]), exz = 0.5 * (L[0][2] + L[2][0]),
                   eyz = 0.5 * (L[1][2] + L[2][1]);
            double Rxy = 0.5 * (L[0][1] - L[1][0]), Rxz = 0.5 * (L[0][2] - L[2][0]),
                   Ryz = 0.5 * (L[1][2] - L[2][1]);
            double tr = (exx + eyy + ezz) / 3.0;
            const Sym3& s = S[i];
            // Jaumann J = R*S - S*R  (R antisym with R[0][1]=Rxy, etc.)
            double R[3][3] = {{0, Rxy, Rxz}, {-Rxy, 0, Ryz}, {-Rxz, -Ryz, 0}};
            double Sm[3][3] = {{s.xx, s.xy, s.xz}, {s.xy, s.yy, s.yz}, {s.xz, s.yz, s.zz}};
            double J[3][3];
            for (int a = 0; a < 3; a++)
                for (int bb = 0; bb < 3; bb++) {
                    double rs = 0, sr = 0;
                    for (int c = 0; c < 3; c++) { rs += R[a][c] * Sm[c][bb]; sr += Sm[a][c] * R[c][bb]; }
                    J[a][bb] = rs - sr;
                }
            dSdt[i].xx = 2 * G * (exx - tr) + J[0][0];
            dSdt[i].yy = 2 * G * (eyy - tr) + J[1][1];
            dSdt[i].zz = 2 * G * (ezz - tr) + J[2][2];
            dSdt[i].xy = 2 * G * exy + J[0][1];
            dSdt[i].xz = 2 * G * exz + J[0][2];
            dSdt[i].yz = 2 * G * eyz + J[1][2];
            // deviatoric power S:eps_dot / rho
            eps_work[i] = (s.xx * exx + s.yy * eyy + s.zz * ezz
                           + 2 * (s.xy * exy + s.xz * exz + s.yz * eyz)) / rho[i];
        });
    }

    void compute_forces() {
        parallel_for(n, [&](int i) {
            const Material& m = materials[mat[i]];
            P[i] = m.pressure(rho[i], u[i]);
            cs[i] = m.sound_speed(rho[i], u[i]);
            csig[i] = std::sqrt(cs[i] * cs[i] + (4.0 / 3.0) * m.G / std::max(rho[i], 1e-30));
            double pf = P[i], dsc = 1.0;
            if (damage_on) {
                sigmax[i] = -P[i] + eig_max(S[i]);      // max principal stress (+ = tension)
                dsc = 1.0 - D[i];
                if (P[i] < 0) pf = dsc * P[i];          // damaged material sheds tension
            }
            Pf[i] = pf; dscale[i] = dsc;
            // Monaghan artificial-stress coeff: repulsive only under tension (pf<0)
            Rart[i] = (eps_as > 0 && pf < 0) ? eps_as * (-pf) / (rho[i] * rho[i]) : 0.0;
        });
        compute_strain_stress();
        parallel_for(n, [&](int i) {
            Vec3 a{}; double du = 0;
            double iri2 = 1.0 / (rho[i] * rho[i]);
            double Pi_over = Pf[i] * iri2;
            double di = dscale[i];
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
                    double cbar = 0.5 * (csig[i] + csig[j]);
                    Pij = (-alpha * cbar * mu + beta * mu * mu) / rbar;
                }
                double jrj2 = 1.0 / (rho[j] * rho[j]);
                double term = Pi_over + Pf[j] * jrj2 + Pij;     // -P (damaged) + AV
                du += 0.5 * mass[j] * term * dot(vij, gradW);    // energy: physical + AV only
                double mterm = term;
                if (eps_as > 0) {       // + Monaghan artificial stress (momentum only, no work)
                    double wdp = kW(hbar / eta, hbar);           // kernel at the particle spacing
                    double fr = (wdp > 0 ? kW(r, hbar) / wdp : 0.0);
                    double fn = fr * fr; fn *= fn;               // (W/W_dp)^4
                    mterm += (Rart[i] + Rart[j]) * fn;
                }
                a -= gradW * (mass[j] * mterm);
                // deviatoric stress term, scaled by (1-D) damage on each side
                Vec3 sterm = Sdot(S[i], gradW) * (iri2 * di) + Sdot(S[j], gradW) * (jrj2 * dscale[j]);
                a += sterm * mass[j];
                if (alpha_u > 0) {
                    double vsig = std::sqrt(std::abs(P[i] - P[j]) / rbar);
                    du += alpha_u * mass[j] * vsig * (u[i] - u[j]) * kdWdr(r, hbar) / rbar;
                }
            });
            acc[i] = a; dudt[i] = du + eps_work[i] * di;
        });
    }

    void yield_limit(int i) {   // von Mises radial return
        double Y = materials[mat[i]].Y;
        Sym3& s = S[i];
        if (Y <= 0) { s = Sym3{}; return; }
        double J2 = 0.5 * (s.xx * s.xx + s.yy * s.yy + s.zz * s.zz)
                    + (s.xy * s.xy + s.xz * s.xz + s.yz * s.yz);
        double vm = std::sqrt(3.0 * J2);
        if (vm > Y) {
            double f = Y / vm;
            s.xx *= f; s.yy *= f; s.zz *= f; s.xy *= f; s.xz *= f; s.yz *= f;
        }
    }

    // Benz-Asphaug crack growth: once the local tensile strain exceeds a flaw's
    // activation strain, damage grows as d(D^1/3)/dt = c_g/R_s with crack speed
    // c_g ~ 0.4 c_long and R_s ~ particle radius (dx/2). D saturates at 1.
    void grow_damage(double dt) {
        parallel_for(n, [&](int i) {
            const Material& m = materials[mat[i]];
            if (m.Emod <= 0 || D[i] >= 1.0) return;
            double eps_local = sigmax[i] / m.Emod;
            if (sigmax[i] > 0 && eps_local > eps_act[i]) {
                double cg = 0.4 * std::sqrt(m.Emod / std::max(rho[i], 1e-30));
                double Rs = 0.5 * h_init / eta;          // ~ dx/2
                double d13 = std::cbrt(D[i]) + (cg / Rs) * dt;
                D[i] = std::min(1.0, d13 * d13 * d13);
            }
        });
    }

    double timestep() {
        double dt = 1e30;
        for (int i = 0; i < n; i++)
            dt = std::min(dt, cfl * hh[i] / (csig[i] * (1.0 + alpha) + 1e-30));
        return dt;
    }

    void step(double dt) {
        std::vector<double> dudt_old = dudt;
        std::vector<Sym3> dSdt_old = dSdt;
        for (int i = 0; i < n; i++) if (!fixed[i]) vel[i] += acc[i] * (0.5 * dt);
        for (int i = 0; i < n; i++) if (!fixed[i]) pos[i] += vel[i] * dt;
        compute_density_h();
        compute_forces();
        if (damage_on) grow_damage(dt);
        for (int i = 0; i < n; i++) if (!fixed[i]) {
            vel[i] += acc[i] * (0.5 * dt);
            u[i] += 0.5 * dt * (dudt_old[i] + dudt[i]);
            if (u[i] < 1e-12) u[i] = 1e-12;
            S[i].xx += 0.5 * dt * (dSdt_old[i].xx + dSdt[i].xx);
            S[i].yy += 0.5 * dt * (dSdt_old[i].yy + dSdt[i].yy);
            S[i].zz += 0.5 * dt * (dSdt_old[i].zz + dSdt[i].zz);
            S[i].xy += 0.5 * dt * (dSdt_old[i].xy + dSdt[i].xy);
            S[i].xz += 0.5 * dt * (dSdt_old[i].xz + dSdt[i].xz);
            S[i].yz += 0.5 * dt * (dSdt_old[i].yz + dSdt[i].yz);
            yield_limit(i);
        }
    }

    void init() { compute_density_h(); compute_forces(); }
};
