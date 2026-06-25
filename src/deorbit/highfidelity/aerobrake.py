#!/usr/bin/env python3
"""Multi-pass aerobraking sequence at full fidelity (REBOUND + ASSIST + drag).

One continuous geocentric ASSIST integration (Sun/planets/Moon/asteroids +
Earth J2/J3/J4 + Sun J2 + GR) so that lunar + J2 perturbations accumulate
correctly across the long coasts. Drag is applied by operator-splitting only
when the body is below the atmospheric interface. Adaptive stepping: large
steps on the drag-free coast (out to apoapsis and back), fine steps through
each periapsis passage.

Per pass we log osculating elements (apo/peri altitude, e, i, RAAN), the
passage min altitude / peak dynamic pressure / dv, the geocentric Moon
distance, and elapsed time. Checkpointed to a state dir so a long run survives
interruption.

Usage:
  ../.venv/bin/python aerobrake.py --vinf 0.3 --peri 8 --inc 5 --passes 400
"""
import os, sys, json, time, math, argparse
import numpy as np
from . import hifi
from .hifi import AU, DAY, RE, MU_E, ENTRY_ALT


def osculating(r, v):
    rmag = math.sqrt(r @ r)
    vv = v @ v
    energy = 0.5 * vv - MU_E / rmag
    a = -MU_E / (2.0 * energy)
    h = np.cross(r, v)
    hmag = np.linalg.norm(h)
    e = math.sqrt(max(0.0, 1.0 + 2.0 * energy * hmag * hmag / (MU_E * MU_E)))
    inc = math.degrees(math.acos(max(-1, min(1, h[2] / hmag))))
    n = np.array([-h[1], h[0], 0.0])
    raan = math.degrees(math.atan2(n[1], n[0])) if np.linalg.norm(n) > 0 else 0.0
    ra = a * (1 + e); rp = a * (1 - e)
    return dict(a_km=a / 1e3, e=e, inc=inc, raan=raan,
                apo_alt_km=(ra - RE) / 1e3, peri_alt_km=(rp - RE) / 1e3,
                energy=energy)


def moon_geo_dist_m(t_days):
    e = hifi.ephem()
    ea = e.get_particle("earth", t_days)
    mo = e.get_particle("moon", t_days)
    return math.sqrt((mo.x - ea.x) ** 2 + (mo.y - ea.y) ** 2 + (mo.z - ea.z) ** 2) * AU


HILL_M = 1.496e9            # Earth Hill radius ~1.5e6 km
MOON_HILL_M = 6.6e7         # ~66,000 km; flag close lunar approaches
GRAZING_FPA_DEG = 5.0       # |flight-path angle| below this = grazing impact


def classify(outcome, rows):
    """Map a run to a physical regime + the diagnostic row."""
    imp = [r for r in rows if r.get("event") == "impact"]
    if outcome == "impact" and imp:
        m = imp[0]
        regime = "grazing_impact" if abs(m["fpa_deg"]) < GRAZING_FPA_DEG else "steep_impact"
        return regime, m
    if outcome == "escaped":
        return "ejected", None
    if outcome == "lunar_encounter":
        return "lunar_encounter", None
    return "aerobraking", (rows[-1] if rows else None)   # budget/walltime: unresolved


def run(vinf_kms, peri_km, inc_deg, max_passes, statedir,
        diameter=1000.0, cd=1.0, epoch_days=0.0, max_walltime=None):
    from deorbit.corridor import physics
    body = physics.Body(diameter_m=diameter, cd=cd)
    os.makedirs(statedir, exist_ok=True)
    rowsf = os.path.join(statedir, "passes.json")
    statusf = os.path.join(statedir, "status.txt")

    sim, ext = hifi.make_sim(epoch_days=epoch_days)
    r_m, v_m = hifi.initial_state_si(vinf_kms * 1e3, peri_km * 1e3, inc_deg)
    sim.add(x=float(r_m[0] / AU), y=float(r_m[1] / AU), z=float(r_m[2] / AU),
            vx=float(v_m[0] / AU * DAY), vy=float(v_m[1] / AU * DAY),
            vz=float(v_m[2] / AU * DAY))

    rows = []
    meta = dict(vinf_kms=vinf_kms, peri_km=peri_km, inc_deg=inc_deg,
                diameter_m=diameter, cd=cd, beta=body.beta, epoch_days=epoch_days)
    pass_count = 0
    in_atm = False
    pmin = ENTRY_ALT; pq = 0.0; v_entry = 0.0
    min_moon = 1e30
    t0 = epoch_days
    wall0 = time.time()
    outcome = "budget"

    def checkpoint():
        json.dump({"meta": meta, "passes": rows}, open(rowsf + ".tmp", "w"))
        os.replace(rowsf + ".tmp", rowsf)
        last = rows[-1] if rows else {}
        with open(statusf, "w") as fh:
            fh.write(f"aerobrake vinf={vinf_kms} km/s peri={peri_km} km inc={inc_deg} deg\n"
                     f"beta={body.beta:.3e} kg/m^2\n"
                     f"passes={pass_count}/{max_passes}  outcome={outcome}\n"
                     f"last apo_alt={last.get('apo_alt_km',0):,.0f} km  "
                     f"peri_alt={last.get('peri_alt_km',0):.1f} km\n"
                     f"min Moon approach={min_moon/1e3:,.0f} km\n"
                     f"t_elapsed={ (last.get('t_days',t0)-t0):.1f} days\n"
                     f"updated {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    steps = 0
    while True:
        r, v = hifi.geo_state_si(sim)
        rmag = math.sqrt(r @ r)
        alt = rmag - RE
        vr = float(r @ v) / rmag

        if alt <= 0.0:
            outcome = "impact"
            r_unit = r / rmag
            lat = math.degrees(math.asin(r_unit[2]))
            rows.append(dict(pass_n=pass_count + 1, event="impact",
                             speed_kms=math.sqrt(v @ v) / 1e3,
                             fpa_deg=math.degrees(math.asin(max(-1, min(1, vr / math.sqrt(v @ v))))),
                             lat_deg=lat, t_days=sim.t))
            checkpoint()
            break
        if rmag > HILL_M:
            outcome = "escaped"; checkpoint(); break

        md = moon_geo_dist_m(sim.t)
        # body-Moon distance (approx: body geocentric vs Moon geocentric)
        e = hifi.ephem()
        ea = e.get_particle("earth", sim.t); mo = e.get_particle("moon", sim.t)
        moon_geo = np.array([mo.x - ea.x, mo.y - ea.y, mo.z - ea.z]) * AU
        bm = np.linalg.norm(r - moon_geo)
        if bm < min_moon:
            min_moon = bm

        if alt < ENTRY_ALT:
            if not in_atm:
                in_atm = True
                pass_count += 1
                pmin = alt; pq = 0.0; v_entry = math.sqrt(v @ v)
            dt_s = max(0.02, min(0.5, 200.0 / (abs(vr) + 10.0)))
            sim.integrate(sim.t + dt_s / DAY)
            r2, v2 = hifi.geo_state_si(sim)
            a_drag, _, q = hifi._drag_kick_si(r2, v2, body)
            if q > pq:
                pq = q
            a2 = math.sqrt(r2 @ r2) - RE
            if a2 < pmin:
                pmin = a2
            hifi.set_geo_state_si(sim, r2, v2 + a_drag * dt_s)
        else:
            if in_atm:
                in_atm = False
                osc = osculating(r, v)
                dvp = v_entry - math.sqrt(v @ v)
                row = dict(pass_n=pass_count, t_days=sim.t,
                           min_alt_km=pmin / 1e3, peak_q_MPa=pq / 1e6, dv_ms=dvp,
                           moon_dist_km=bm / 1e3, **osc)
                rows.append(row)
                if pass_count % 5 == 0 or pass_count <= 5:
                    checkpoint()
                if bm < MOON_HILL_M:
                    outcome = "lunar_encounter"; checkpoint(); break
                if osc["energy"] > 0:
                    outcome = "escaped"; checkpoint(); break
                if pass_count >= max_passes:
                    outcome = "budget"; checkpoint(); break
            # drag-free coast: large adaptive step, shrink near interface so we
            # never overshoot the atmosphere; big steps near slow apoapsis.
            if vr < 0:   # descending toward interface
                dt_s = min(1200.0, max(5.0, 0.25 * (alt - ENTRY_ALT) / (abs(vr) + 1.0)))
            else:        # ascending / near apoapsis
                dt_s = min(3600.0, max(20.0, 0.05 * alt / (abs(vr) + 1.0)))
            sim.integrate(sim.t + dt_s / DAY)
        steps += 1
        if max_walltime and (time.time() - wall0) > max_walltime:
            outcome = "walltime"; checkpoint(); break

    checkpoint()
    return outcome, rows, meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vinf", type=float, default=0.3, help="km/s")
    ap.add_argument("--peri", type=float, default=8.0, help="periapsis altitude, km")
    ap.add_argument("--inc", type=float, default=5.0, help="inclination, deg")
    ap.add_argument("--passes", type=int, default=400)
    ap.add_argument("--diameter", type=float, default=1000.0)
    ap.add_argument("--walltime", type=float, default=None, help="max wall seconds")
    ap.add_argument("--statedir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "state_aerobrake"))
    a = ap.parse_args()
    t0 = time.time()
    outcome, rows, meta = run(a.vinf, a.peri, a.inc, a.passes, a.statedir,
                              diameter=a.diameter, max_walltime=a.walltime)
    passes = [r for r in rows if "apo_alt_km" in r]
    print(f"\n=== aerobrake outcome: {outcome} ===")
    print(f"IC: v_inf={a.vinf} km/s, peri={a.peri} km, inc={a.inc} deg, "
          f"beta={meta['beta']:.2e} kg/m^2  ({time.time()-t0:.0f}s wall)")
    if passes:
        first, last = passes[0], passes[-1]
        print(f"first apogee alt = {first['apo_alt_km']:,.0f} km "
              f"(Moon orbit = 384,400 km)")
        print(f"completed passes = {len(passes)}  elapsed = {last['t_days']-meta['epoch_days']:.1f} days")
        print(f"min Moon approach = {min(r['moon_dist_km'] for r in passes):,.0f} km")
        print(f"\n{'pass':>5} {'apo_alt_km':>12} {'peri_alt_km':>11} {'dv_ms':>8} "
              f"{'q_MPa':>7} {'inc':>6} {'raan':>7} {'t_days':>8}")
        idx = sorted(set([0, 1, 2] + list(range(0, len(passes), max(1, len(passes)//15))) + [len(passes)-1]))
        idx = [i for i in idx if 0 <= i < len(passes)]
        for i in idx:
            r = passes[i]
            print(f"{r['pass_n']:>5} {r['apo_alt_km']:>12,.0f} {r['peri_alt_km']:>11.1f} "
                  f"{r['dv_ms']:>8.1f} {r['peak_q_MPa']:>7.1f} {r['inc']:>6.2f} "
                  f"{r['raan']:>7.2f} {r['t_days']-meta['epoch_days']:>8.2f}")
    imp = [r for r in rows if r.get("event") == "impact"]
    if imp:
        m = imp[0]
        print(f"\nIMPACT: speed={m['speed_kms']:.2f} km/s  fpa={m['fpa_deg']:.2f} deg  "
              f"lat={m['lat_deg']:.2f} deg")
