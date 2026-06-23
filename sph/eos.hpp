// Equations of state: ideal gas (for the Sod gate) + Tillotson (iron, basalt,
// water) for impacts. Header-only.
//
// Tillotson (1962) parameters from the standard compilations (Melosh 1989
// "Impact Cratering" App. II; Benz & Asphaug 1999). Validated in eos_test.cpp:
// P(rho0, u=0) = 0 and the small-strain sound speed c0 = sqrt(A/rho0) matches
// known bulk sound speeds (iron ~4.0, basalt ~3.1, water ~1.48 km/s). The
// vaporization-branch energies (uiv,ucv) should be re-checked against the
// reference in the stress-test pass; the condensed branch (A,B,a,b) dominates
// the impact contact and is what the sound-speed check pins down.
#pragma once
#include <cmath>
#include <algorithm>

struct Material {
    enum Type { IDEAL, TILLOTSON };
    Type type = IDEAL;
    double gamma = 1.4;
    // Tillotson
    double rho0 = 0, a = 0, b = 0, A = 0, B = 0, alpha = 0, beta = 0;
    double u0 = 0, uiv = 0, ucv = 0;
    // strength
    double G = 0, Y = 0;
    // Young's modulus E = 9KG/(3K+G) with K≈A, and Weibull flaw params for
    // brittle damage (Benz-Asphaug): wm = Weibull modulus, wk = flaw density
    // (m^-3). wk=0 -> no brittle fracture (ductile metals).
    double Emod = 0, wm = 0, wk = 0;

    static Material ideal(double g) { Material m; m.type = IDEAL; m.gamma = g; return m; }
    static Material tillotson(double rho0_, double A_, double B_, double a_, double b_,
                              double al, double be, double u0_, double uiv_, double ucv_,
                              double G_ = 0, double Y_ = 0, double wm_ = 0, double wk_ = 0) {
        Material m; m.type = TILLOTSON;
        m.rho0 = rho0_; m.A = A_; m.B = B_; m.a = a_; m.b = b_;
        m.alpha = al; m.beta = be; m.u0 = u0_; m.uiv = uiv_; m.ucv = ucv_;
        m.G = G_; m.Y = Y_; m.wm = wm_; m.wk = wk_;
        m.Emod = (G_ > 0 ? 9.0 * A_ * G_ / (3.0 * A_ + G_) : 0.0);
        return m;
    }
    // standard parameter sets (SI).  Weibull: basalt m=16, k=1e61 m^-3 (Benz &
    // Asphaug 1995 basalt) -> activation strains ~3-8e-5; iron ductile (wk=0).
    static Material iron()   { return tillotson(7800, 1.28e11, 1.05e11, 0.5, 1.5, 5, 5,
                                                9.5e6, 2.4e6, 8.67e6, 7.75e10, 1.0e9); }
    static Material basalt() { return tillotson(2700, 2.67e10, 2.67e10, 0.5, 1.5, 5, 5,
                                                4.87e8, 4.72e6, 1.82e7, 2.27e10, 3.5e8, 16.0, 1.0e61); }
    static Material water()  { return tillotson(998, 2.2e9, 9.94e9, 0.7, 0.15, 10, 5,
                                                7.0e6, 4.19e5, 2.69e6, 0.0, 0.0); }

    double till(double rho, double u) const {
        if (rho <= 0) return 0.0;
        double eta = rho / rho0, mu = eta - 1.0;
        double w0 = u / (u0 * eta * eta) + 1.0;
        double Pc = (a + b / w0) * rho * u + A * mu + B * mu * mu;   // condensed
        if (rho >= rho0 || u <= uiv) return Pc;
        double z = rho0 / rho - 1.0;
        double Pe = a * rho * u + (b * rho * u / w0 + A * mu * std::exp(-beta * z))
                                  * std::exp(-alpha * z * z);        // expanded hot
        if (u >= ucv) return Pe;
        return ((u - uiv) * Pe + (ucv - u) * Pc) / (ucv - uiv);      // intermediate
    }

    double pressure(double rho, double u) const {
        if (type == IDEAL) return (gamma - 1.0) * rho * u;
        return till(rho, u);
    }
    // ANALYTIC Tillotson derivatives (replaces finite differences, which suffer
    // catastrophic cancellation in dP/du when the density term dwarfs the energy
    // term -- and which break entirely in FP32 on the GPU). Condensed branch:
    void cond_d(double rho, double u, double& P, double& dPdr, double& dPdu) const {
        double eta = rho / rho0, mu = eta - 1.0;
        double w = 1.0 + u / (u0 * eta * eta), aob = a + b / w;
        P = aob * rho * u + A * mu + B * mu * mu;
        dPdu = aob * rho - b * rho * u / (w * w * u0 * eta * eta);
        double dwdr = -2.0 * u / (u0 * eta * eta * eta * rho0);
        dPdr = aob * u + rho * u * (-b / (w * w)) * dwdr + A / rho0 + 2.0 * B * mu / rho0;
    }
    void expd_d(double rho, double u, double& P, double& dPdr, double& dPdu) const {
        double eta = rho / rho0, mu = eta - 1.0;
        double w = 1.0 + u / (u0 * eta * eta);
        double z = rho0 / rho - 1.0, dzdr = -rho0 / (rho * rho);
        double E1 = std::exp(-alpha * z * z), E2 = std::exp(-beta * z);
        double G = b * rho * u / w + A * mu * E2;
        P = a * rho * u + G * E1;
        double dwdr = -2.0 * u / (u0 * eta * eta * eta * rho0);
        double dGu = b * rho / w - b * rho * u / (w * w * u0 * eta * eta);
        dPdu = a * rho + dGu * E1;
        double dbruw_dr = b * u / w + b * rho * u * (-1.0 / (w * w)) * dwdr;
        double dGr = dbruw_dr + A * E2 / rho0 + A * mu * (E2 * (-beta) * dzdr);
        double dE1dr = E1 * (-alpha * 2.0 * z) * dzdr;
        dPdr = a * u + dGr * E1 + G * dE1dr;
    }
    double sound_speed(double rho, double u) const {
        if (type == IDEAL) {
            double P = (gamma - 1.0) * rho * u;
            return std::sqrt(std::max(gamma * P / std::max(rho, 1e-30), 1e-12));
        }
        if (rho <= 0) return 1e-3;
        double P, dPdr, dPdu;
        if (rho >= rho0 || u <= uiv) cond_d(rho, u, P, dPdr, dPdu);
        else if (u >= ucv) expd_d(rho, u, P, dPdr, dPdu);
        else {                                          // intermediate: blend
            double Pc, dPcr, dPcu, Pe, dPer, dPeu;
            cond_d(rho, u, Pc, dPcr, dPcu); expd_d(rho, u, Pe, dPer, dPeu);
            double inv = 1.0 / (ucv - uiv);
            P    = ((u - uiv) * Pe + (ucv - u) * Pc) * inv;
            dPdu = ((Pe + (u - uiv) * dPeu) + (-Pc + (ucv - u) * dPcu)) * inv;
            dPdr = ((u - uiv) * dPer + (ucv - u) * dPcr) * inv;
        }
        return std::sqrt(std::max(dPdr + P / (rho * rho) * dPdu, 1e-6));
    }
};
