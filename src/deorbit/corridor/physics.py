"""Core flight physics for the deorbit investigation (pure stdlib).

A km-scale monolithic IRON body has an enormous ballistic coefficient
beta = m/(Cd*A) ~ 5e6 kg/m^2, so atmospheric drag barely touches it on any
single pass. The central question of the project is whether such a body,
arriving on a low-energy (low v_inf) Earth encounter, can shed enough energy
in a grazing periapsis passage to be CAPTURED (bound, eps<0) rather than
SKIP back out (eps>=0) or directly IMPACT.

This module integrates a single planar atmospheric passage in an
Earth-centered inertial frame (point-mass gravity + drag) and classifies the
outcome. The capture-corridor sweep in deorbit.py drives it over (v_inf,
periapsis-altitude) grids.

Ground-truthable limits (see validate.py):
  * drag off  -> energy/angular-momentum conserved, Kepler orbit closes.
  * per-pass dv at a given periapsis altitude matches the analytic grazing
    column estimate  dv ~ v * Sigma/beta, Sigma ~ rho(h_p)*sqrt(2*pi*Re*H).
"""
import math

# --- constants (SI) ---
MU = 3.986004418e14      # Earth GM, m^3/s^2
RE = 6.371e6             # Earth mean radius, m
RHO_IRON = 7800.0        # kg/m^3
ENTRY_ALT = 200e3        # atmospheric-interface altitude for passage start/end, m

# Vallado piecewise-exponential atmosphere: (base_alt_km, rho0 kg/m^3, H_km)
_ATM = [
    (0, 1.225, 7.249), (25, 3.899e-2, 6.349), (30, 1.774e-2, 6.682),
    (40, 3.972e-3, 7.554), (50, 1.057e-3, 8.382), (60, 3.206e-4, 7.714),
    (70, 8.770e-5, 6.549), (80, 1.905e-5, 5.799), (90, 3.396e-6, 5.382),
    (100, 5.297e-7, 5.877), (110, 9.661e-8, 7.263), (120, 2.438e-8, 9.473),
    (130, 8.484e-9, 12.636), (140, 3.845e-9, 16.149), (150, 2.070e-9, 22.523),
    (180, 5.464e-10, 29.740), (200, 2.789e-10, 37.105), (250, 7.248e-11, 45.546),
    (300, 2.418e-11, 53.628), (350, 9.518e-12, 53.298), (400, 3.725e-12, 58.515),
    (450, 1.585e-12, 60.828), (500, 6.967e-13, 63.822), (600, 1.454e-13, 71.835),
    (700, 3.614e-14, 88.667), (800, 1.170e-14, 124.64), (900, 5.245e-15, 181.05),
    (1000, 3.019e-15, 268.00),
]


def density(alt_m):
    """Mass density (kg/m^3) at geometric altitude alt_m (m)."""
    if alt_m <= 0.0:
        return _ATM[0][1]
    h = alt_m / 1000.0
    if h >= 1000.0:
        return 0.0
    base = _ATM[0]
    for row in _ATM:
        if row[0] <= h:
            base = row
        else:
            break
    h0, rho0, H = base
    return rho0 * math.exp(-(h - h0) / H)


class Body:
    """Monolithic sphere. beta = m/(Cd*A)."""

    def __init__(self, diameter_m=1000.0, rho=RHO_IRON, cd=1.0):
        self.diameter = diameter_m
        self.rho = rho
        self.cd = cd
        r = diameter_m / 2.0
        self.mass = rho * (4.0 / 3.0) * math.pi * r ** 3
        self.area = math.pi * r ** 2
        self.beta = self.mass / (cd * self.area)
        # drag accel = -k * |v| * v_vec, with k = 0.5*rho*Cd*A/m
        self.k = 0.5 * cd * self.area / self.mass


def _deriv(s, k):
    x, y, vx, vy = s
    r = math.hypot(x, y)
    ag = -MU / (r * r * r)
    ax, ay = ag * x, ag * y
    rho = density(r - RE)
    if rho > 0.0:
        v = math.hypot(vx, vy)
        if v > 0.0:
            fd = k * rho * v
            ax -= fd * vx
            ay -= fd * vy
    return (vx, vy, ax, ay)


def _rk4(s, dt, k):
    a = _deriv(s, k)
    s2 = tuple(s[i] + 0.5 * dt * a[i] for i in range(4))
    b = _deriv(s2, k)
    s3 = tuple(s[i] + 0.5 * dt * b[i] for i in range(4))
    c = _deriv(s3, k)
    s4 = tuple(s[i] + dt * c[i] for i in range(4))
    d = _deriv(s4, k)
    return tuple(s[i] + dt / 6.0 * (a[i] + 2 * b[i] + 2 * c[i] + d[i]) for i in range(4))


def initial_state(v_inf, peri_alt_m, r0=RE + ENTRY_ALT):
    """Inbound state at the interface r0 for an orbit with hyperbolic-excess
    speed v_inf (m/s) and vacuum periapsis radius RE+peri_alt_m."""
    rp = RE + peri_alt_m
    vp = math.sqrt(v_inf * v_inf + 2.0 * MU / rp)
    h_ang = rp * vp                      # specific angular momentum (fpa=0 at peri)
    v0 = math.sqrt(v_inf * v_inf + 2.0 * MU / r0)
    c = max(-1.0, min(1.0, h_ang / (r0 * v0)))
    fpa = math.acos(c)
    vr = -v0 * math.sin(fpa)             # inbound (radial inward)
    vt = v0 * math.cos(fpa)              # prograde transverse
    # at position (r0, 0): radial unit (1,0), transverse unit (0,1)
    return [r0, 0.0, vr, vt], v0


def simulate_passage(v_inf, peri_alt_m, body, max_steps=4_000_000):
    """Integrate one atmospheric passage. Returns a dict with the outcome:
    'impact' | 'capture' | 'skip', plus diagnostics."""
    r0 = RE + ENTRY_ALT
    s, v0 = initial_state(v_inf, peri_alt_m, r0)
    k = body.k
    min_alt = ENTRY_ALT
    peak_q = 0.0
    passed_peri = False
    steps = 0
    outcome = "timeout"
    while steps < max_steps:
        x, y, vx, vy = s
        r = math.hypot(x, y)
        alt = r - RE
        if alt < min_alt:
            min_alt = alt
        v = math.hypot(vx, vy)
        q = 0.5 * density(alt) * v * v
        if q > peak_q:
            peak_q = q
        if alt <= 0.0:
            outcome = "impact"
            break
        vr = (x * vx + y * vy) / r        # radial speed
        if vr > 0.0:
            passed_peri = True
        if passed_peri and r >= r0 and vr > 0.0:
            # exited the interface, climbing
            eps = 0.5 * v * v - MU / r
            outcome = "capture" if eps < 0.0 else "skip"
            break
        # adaptive step: fine in atmosphere, coarse on the way out
        if alt < ENTRY_ALT:
            dt = max(0.02, min(0.5, 200.0 / (abs(vr) + 10.0)))
        else:
            dt = 20.0
        s = _rk4(s, dt, k)
        steps += 1

    x, y, vx, vy = s
    r = math.hypot(x, y)
    v = math.hypot(vx, vy)
    eps = 0.5 * v * v - MU / r
    h_exit = abs(x * vy - y * vx)
    res = {
        "outcome": outcome,
        "v_inf": v_inf,
        "peri_alt_km": peri_alt_m / 1000.0,
        "min_alt_km": min_alt / 1000.0,
        "peak_q_MPa": peak_q / 1e6,
        "dv_ms": v0 - v,                  # speed lost across the passage
        "eps_out": eps,
        "steps": steps,
    }
    if outcome == "capture":
        a = -MU / (2.0 * eps)
        e = math.sqrt(max(0.0, 1.0 + 2.0 * eps * h_exit * h_exit / (MU * MU)))
        res["apoapsis_km"] = (a * (1.0 + e) - RE) / 1000.0
        res["new_peri_km"] = (a * (1.0 - e) - RE) / 1000.0
    return res
