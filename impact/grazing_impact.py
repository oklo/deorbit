#!/usr/bin/env python3
"""Phase 2a: reduced-order model of the ULTRA-OBLIQUE terminal impact (stdlib).

The de-orbiting km iron monolith arrives at a flight-path angle gamma ~ 1-4 deg
(from the Phase-1 aerobraking runs), v ~ 7.4-11 km/s. This is NOT vertical
cratering: the normal velocity v*sin(gamma) is small (subsonic in rock) while
the along-track component is hypersonic, so the body ploughs/ricochets and
shock-melts its leading contact rather than excavating a deep crater
(Schultz & Gault 1990: ricochet/decapitation below ~10-15 deg, 60-80% velocity
retained; elongated craters below ~15 deg, extreme below 5 deg).

These are ORDER-OF-MAGNITUDE scaling estimates to bracket the outcome and define
inputs for a 3D hydrocode (Phase 2b) -- not a substitute for it. Each quantity
is tied to a named scaling. Distinguish "scaling estimate" from "verified".

Targets: Cayambe (high equatorial volcanic rock), deep ocean (water/basalt),
Amazon basin (water-saturated sediment).
"""
import math

TNT_J = 4.184e9                         # J per ton TNT


class Body:
    def __init__(self, diameter_m=1000.0, rho=7870.0):
        self.d = diameter_m
        self.rho = rho
        r = diameter_m / 2.0
        self.mass = rho * (4.0 / 3.0) * math.pi * r ** 3
        self.area = math.pi * r ** 2     # frontal area
        self.r = r
        # melt / vaporization specific energies for IRON (J/kg)
        self.E_melt = 0.92e6
        self.E_vap = 6.5e6


# target: (name, density, sound_speed, strength_Pa, E_melt/vap J/kg, layer note)
TARGETS = {
    "cayambe": dict(name="Cayambe (volcanic rock)", rho=2900.0, c=5500.0,
                    Y=2.0e8, E_melt=1.2e6, E_vap=8.0e6, water_depth=0.0,
                    note="hard andesite/basalt edifice, ~5.8 km elevation"),
    "ocean":   dict(name="Deep ocean (water over basalt)", rho=1000.0, c=1500.0,
                    Y=0.0, E_melt=1.0e6, E_vap=2.26e6, water_depth=4000.0,
                    note="~4 km water over basaltic crust; tsunami source"),
    "amazon":  dict(name="Amazon basin (wet sediment)", rho=2000.0, c=2500.0,
                    Y=2.0e7, E_melt=1.0e6, E_vap=2.5e6, water_depth=0.0,
                    note="water-saturated sediment over continental crust"),
}

RICOCHET_ANGLE_DEG = 15.0               # below this: ricochet/decapitation regime


def grazing_impact(v_ms, gamma_deg, body, tgt, cd=1.0):
    g = math.radians(gamma_deg)
    v_n = v_ms * math.sin(g)
    v_h = v_ms * math.cos(g)
    KE = 0.5 * body.mass * v_ms ** 2
    rho_t = tgt["rho"]

    # leading-contact dynamic ("stagnation") pressure from along-track motion
    P_contact = rho_t * v_h ** 2

    # ricochet (Schultz & Gault): below ~15 deg the body skips, retaining a
    # large velocity fraction; crude fit ~ (1 - gamma/RICOCHET_ANGLE)*0.8.
    ricochet = gamma_deg < RICOCHET_ANGLE_DEG
    retained = max(0.0, (1.0 - gamma_deg / RICOCHET_ANGLE_DEG)) * 0.8 if ricochet else 0.0

    # bulldozer deceleration length: a = rho_t*Cd*A*v^2/(2m) -> e-fold length
    # L_e = 2*beta/rho_t, beta = m/(Cd*A)  (Newtonian penetration / impact drag)
    beta = body.mass / (cd * body.area)
    L_efold = 2.0 * beta / rho_t

    # energy deposited (not carried off by the ricochet)
    f_coupled = 1.0 - retained ** 2
    E_dep = f_coupled * KE
    # Melt is NOT energy-limited here (KE >> energy to melt the whole body); it is
    # SHOCK-COUPLING-limited. Report (a) the energy-budget ratio and (b) whether
    # the leading-contact shock pressure exceeds iron's shock-melt threshold.
    E_melt_body = body.mass * body.E_melt
    E_vap_body = body.mass * body.E_vap
    P_IRON_MELT = 220e9                  # ~incipient shock-melt of iron on release
    P_IRON_VAP = 500e9                   # ~incipient vaporization

    out = dict(target=tgt["name"], v_kms=v_ms / 1e3, gamma_deg=gamma_deg,
               KE_J=KE, KE_Gt=KE / TNT_J / 1e9,
               v_normal_ms=v_n, v_horizontal_kms=v_h / 1e3,
               P_contact_GPa=P_contact / 1e9,
               normal_mach=v_n / tgt["c"],
               ricochet=ricochet, retained_v_frac=retained,
               L_efold_km=L_efold / 1e3,
               energy_deposited_frac=f_coupled,
               KE_over_bodymelt=E_dep / E_melt_body,
               KE_over_bodyvap=E_dep / E_vap_body,
               contact_vs_melt=P_contact / P_IRON_MELT,
               contact_melts_iron=P_contact >= P_IRON_MELT,
               contact_vaporizes_iron=P_contact >= P_IRON_VAP)

    # target-specific
    if tgt["water_depth"] > 0:
        # velocity at seafloor after ploughing the water column (e-fold L)
        v_seafloor = v_ms * math.exp(-tgt["water_depth"] / L_efold)
        out["reaches_seafloor"] = True
        out["v_at_seafloor_kms"] = v_seafloor / 1e3
        out["water_efold_km"] = L_efold / 1e3
    return out


def fmt(o):
    L = []
    L.append(f"  {o['target']}")
    L.append(f"    v={o['v_kms']:.2f} km/s  gamma={o['gamma_deg']:.2f} deg  "
             f"KE={o['KE_Gt']:.1f} Gt TNT")
    L.append(f"    normal v={o['v_normal_ms']:.0f} m/s (Mach {o['normal_mach']:.2f})   "
             f"along-track v={o['v_horizontal_kms']:.2f} km/s")
    L.append(f"    leading-contact pressure ~ {o['P_contact_GPa']:.0f} GPa "
             f"(>> rock & iron strength -> fluidized/shocked contact, ploughing)")
    L.append(f"    ricochet: {o['ricochet']}  (retains ~{o['retained_v_frac']*100:.0f}% of v)")
    L.append(f"    bulldozer e-fold length ~ {o['L_efold_km']:.1f} km "
             f"-> elongated gouge of this scale")
    L.append(f"    energy deposited (not ricocheted) ~ {o['energy_deposited_frac']*100:.0f}% of KE "
             f"= {o['KE_over_bodymelt']:.0f}x the energy to melt the whole body "
             f"({o['KE_over_bodyvap']:.0f}x to vaporize it)")
    melt = ("vaporizes" if o['contact_vaporizes_iron'] else
            "melts" if o['contact_melts_iron'] else "does NOT melt")
    L.append(f"    -> melt is shock-COUPLING-limited; leading contact ({o['P_contact_GPa']:.0f} GPa) "
             f"{melt} iron at the contact (threshold ~220 GPa) -> melt-lined gouge, surviving fragment ricochets")
    if o.get("reaches_seafloor"):
        L.append(f"    OCEAN: punches through {4.0:.0f} km water (e-fold {o['water_efold_km']:.1f} km) "
                 f"-> ~{o['v_at_seafloor_kms']:.2f} km/s at seafloor; large water cavity + tsunami")
    return "\n".join(L)


if __name__ == "__main__":
    body = Body()
    print(f"Body: {body.d:.0f} m iron sphere, mass {body.mass:.2e} kg "
          f"(beta {body.mass/body.area:.2e} kg/m^2)\n")
    # nominal terminal conditions from the Phase-1 aerobraked grazing impact:
    print("=== Nominal aerobraked grazing impact: v=7.4 km/s, gamma=1.2 deg ===")
    for key in ("cayambe", "ocean", "amazon"):
        print(fmt(grazing_impact(7390.0, 1.16, body, TARGETS[key])))
        print()
    print("=== Higher-energy / steeper variant: v=11 km/s, gamma=4 deg ===")
    for key in ("cayambe", "ocean", "amazon"):
        print(fmt(grazing_impact(11000.0, 4.0, body, TARGETS[key])))
        print()
