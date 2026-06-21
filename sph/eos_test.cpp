// Validation: Tillotson EOS for iron, basalt, water.
//   (1) P(rho0, u=0) = 0    (reference state)
//   (2) small-strain sound speed c0 = sqrt(A/rho0) matches the known bulk sound
//       speed of each material (the strongest independent check: water -> 1481 m/s)
//   (3) compression (eta=1.1) gives a sensible positive pressure
#include <cstdio>
#include <cmath>
#include "eos.hpp"

static int check(const char* name, Material m, double cexp_lo, double cexp_hi) {
    double P0 = m.pressure(m.rho0, 0.0);
    double c0 = m.sound_speed(m.rho0, 0.0);
    double Pc = m.pressure(1.1 * m.rho0, 0.0);     // 10% compression, cold
    bool ok_P0 = std::abs(P0) < 1e3;               // ~0
    bool ok_c = c0 > cexp_lo && c0 < cexp_hi;
    bool ok_comp = Pc > 0;
    printf("  %-7s rho0=%6.0f  P(rho0,0)=%10.3e Pa  c0=%7.1f m/s  "
           "P(1.1 rho0,0)=%8.2f GPa  -> %s\n",
           name, m.rho0, P0, c0, Pc / 1e9,
           (ok_P0 && ok_c && ok_comp) ? "OK" : "FAIL");
    printf("           sqrt(A/rho0)=%7.1f m/s (expect bulk sound speed in [%.0f,%.0f])\n",
           std::sqrt(m.A / m.rho0), cexp_lo, cexp_hi);
    return (ok_P0 && ok_c && ok_comp) ? 1 : 0;
}

int main() {
    printf("Tillotson EOS validation:\n");
    int ok = 1;
    ok &= check("iron",   Material::iron(),   3900, 4200);   // bulk ~4.0 km/s
    ok &= check("basalt", Material::basalt(), 3000, 3300);   // ~3.1 km/s
    ok &= check("water",  Material::water(),  1450, 1520);   // 1481 m/s
    printf("%s\n", ok ? "\nALL PASS" : "\nSOME CHECKS FAILED");
    return ok ? 0 : 1;
}
