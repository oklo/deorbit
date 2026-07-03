"""Arrival prior from heliocentric NEO dynamics -> per-seed catalog weights.

Physics (docs/site_mc_plan.md Phase E): the arrival flux into the capture
corridor is  d(rate) ~ f(U) * U * d(b^2) * dU^3,  where f(U) is the local
velocity-space density of the NEO population at Earth's orbit and b the
impact parameter (gravitational focusing: b^2 = rp^2 + 2 mu rp / U^2).

f(U) is built EMPIRICALLY, not from Opik closed forms: sample orbits
(a,e,i) from a population density, uniform angles/phase, keep samples
inside a thin torus |s-1|<ds, |z|<dz around Earth's (circular, 1 AU) ring,
and record the velocity relative to the local circular velocity. Torus
hits sample the steady-state local density with the correct kinematic
weighting automatically (the Opik 1/|U_x|, 1/sin i pile-ups emerge by
themselves). Gate: |U| of every hit matches the Tisserand form
U^2 = 3 - 1/a - 2 sqrt(a(1-e^2)) cos i.

Population densities (pluggable; the debiased models disagree most in the
Earth-like corner, so prior SENSITIVITY is part of the result):
  flat: constant over the Earth-like box.
  tilt: Granvik-2018-shaped analytic tilt (dN/di ~ i exp(-i/8deg),
        dN/de rising, mild a-slope) — swap for a NEOMOD grid when ingested.

Catalog weighting: each catalog row's arrival direction u_hat is
regenerated from its seed (sample.draw is deterministic); Earth's radial/
transverse/normal frame at the row's epoch comes from the analytic Sun
ephemeris; then
  w  ~  f_KDE(U_local) * U^3 * db^2/drp
(U^3 = flux U times the U^2 Jacobian of the catalog's isotropic-direction
x uniform-speed sampling; db^2/drp = 2 rp + 2 mu/U^2).

Analysis-side tool: numpy allowed (the daemon path stays stdlib).
"""
import csv
import math
import os
import sys

V_EARTH = 29.7846          # km/s, circular at 1 AU
MU_E_KM = 398600.4418      # km^3/s^2
BOX = dict(a=(0.80, 1.30), e=(0.0, 0.35), i=(0.0, 15.0))   # deg
DS, DZ = 0.010, 0.005      # torus half-widths (AU)


def model_density(name):
    if name == "flat":
        return lambda a, e, i: 1.0 + 0.0 * a
    if name == "tilt":
        import numpy as np
        return lambda a, e, i: ((1.0 + 0.5 * (a - 1.0))
                                * (0.05 + e) ** 0.8
                                * (i + 0.5) * np.exp(-i / 8.0))
    raise SystemExit(f"unknown model {name}")


def torus_sample(n, model, seed=1, box=BOX, ds=DS, dz=DZ):
    """Return (U[3,k] km/s, w[k], aei[3,k]) for torus hits."""
    import numpy as np
    rng = np.random.default_rng(seed)
    a = rng.uniform(*box["a"], n)
    e = rng.uniform(*box["e"], n)
    inc = np.radians(rng.uniform(*box["i"], n))
    om = rng.uniform(0.0, 2 * np.pi, n)
    M = rng.uniform(0.0, 2 * np.pi, n)
    E = M.copy()
    for _ in range(6):                       # Kepler (e<=0.35 converges fast)
        E -= (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
    cnu = (np.cos(E) - e) / (1.0 - e * np.cos(E))
    snu = (np.sqrt(1 - e * e) * np.sin(E)) / (1.0 - e * np.cos(E))
    r = a * (1.0 - e * np.cos(E))
    p = a * (1.0 - e * e)
    # perifocal position/velocity (mu=1)
    px, py = r * cnu, r * snu
    vf = 1.0 / np.sqrt(p)
    vx_p, vy_p = -vf * snu, vf * (e + cnu)
    # rotate Rz(om) then Rx(inc)  (node axisymmetry -> Omega omitted)
    co, so = np.cos(om), np.sin(om)
    x1, y1 = co * px - so * py, so * px + co * py
    vx1, vy1 = co * vx_p - so * vy_p, so * vx_p + co * vy_p
    ci, si = np.cos(inc), np.sin(inc)
    x, y, z = x1, ci * y1, si * y1
    vx, vy, vz = vx1, ci * vy1, si * vy1
    s = np.hypot(x, y)
    hit = (np.abs(s - 1.0) < ds) & (np.abs(z) < dz)
    x, y, z, vx, vy, vz = (q[hit] for q in (x, y, z, vx, vy, vz))
    ah, eh, ih = a[hit], e[hit], np.degrees(inc[hit])
    s = np.hypot(x, y)
    cphi, sphi = x / s, y / s
    ur = vx * cphi + vy * sphi           # radial
    ut = -vx * sphi + vy * cphi - 1.0    # transverse minus Earth circular
    U = np.vstack([ur, ut, vz]) * V_EARTH
    w = np.asarray(model(ah, eh, ih), dtype=float)
    return U, w, np.vstack([ah, eh, ih])


def tisserand_U(a, e, i_deg):
    import numpy as np
    p = a * (1 - e * e)
    u2 = 3.0 - 1.0 / a - 2.0 * np.sqrt(p) * np.cos(np.radians(i_deg))
    return np.sqrt(np.maximum(u2, 0.0)) * V_EARTH


class KDE:
    """Equal-weight Gaussian KDE on <=nmax weighted-resampled support."""

    def __init__(self, U, w, h=0.4, nmax=8000, seed=2):
        import numpy as np
        self.np = np
        k = U.shape[1]
        idx = np.random.default_rng(seed).choice(
            k, size=min(nmax, k), p=w / w.sum(), replace=True)
        self.sup = U[:, idx].T.copy()        # (m,3)
        self.h = h

    def __call__(self, u):
        np = self.np
        d2 = ((self.sup - np.asarray(u)) ** 2).sum(axis=1)
        return float(np.exp(-0.5 * d2 / (self.h * self.h)).sum())


def earth_frame(jd):
    """Radial/transverse/normal unit vectors of Earth's heliocentric motion,
    in the geocentric equatorial frame of the catalog."""
    from deorbit.sitemc import ephem
    s0 = ephem.sun_pos(jd - 0.5)
    s1 = ephem.sun_pos(jd + 0.5)
    rE = tuple(-c for c in ephem.sun_pos(jd))          # helio Earth position
    vE = tuple(-(b - a) / 86400.0 for a, b in zip(s0, s1))
    xn = math.sqrt(sum(c * c for c in rE))
    xh = tuple(c / xn for c in rE)
    dot = sum(v * x for v, x in zip(vE, xh))
    ty = tuple(v - dot * x for v, x in zip(vE, xh))
    tn = math.sqrt(sum(c * c for c in ty))
    yh = tuple(c / tn for c in ty)
    zh = (xh[1] * yh[2] - xh[2] * yh[1], xh[2] * yh[0] - xh[0] * yh[2],
          xh[0] * yh[1] - xh[1] * yh[0])
    return xh, yh, zh, math.sqrt(sum(c * c for c in vE))


def weight_catalog(catalog, kde, out_path):
    from deorbit.sitemc import sample as smp
    rows = [r for r in csv.DictReader(open(catalog))
            if r["outcome"] == "impact" and int(r["n_pass"]) >= 2]
    wrows, wsum, w2sum = [], 0.0, 0.0
    for r in rows:
        seed = int(r["seed"])
        d = smp.draw(seed)
        xh, yh, zh, _ = earth_frame(d["epoch_jd"])
        u = d["u_hat"]
        v_kms = d["v_inf"] / 1e3
        Ul = (v_kms * sum(a * b for a, b in zip(u, xh)),
              v_kms * sum(a * b for a, b in zip(u, yh)),
              v_kms * sum(a * b for a, b in zip(u, zh)))
        f = kde(Ul)
        rp_km = d["rp_alt_m"] / 1e3 + 6378.137
        db2 = 2.0 * rp_km + 2.0 * MU_E_KM / (v_kms * v_kms)
        w = f * v_kms ** 3 * db2
        wrows.append((seed, w))
        wsum += w
        w2sum += w * w
    ess = wsum * wsum / w2sum if w2sum else 0.0
    with open(out_path, "w", newline="") as fo:
        wr = csv.writer(fo)
        wr.writerow(["seed", "weight"])
        for s, w in wrows:
            wr.writerow([s, f"{w:.6g}"])
    print(f"wrote {out_path}: {len(wrows)} rows, ESS={ess:.0f} "
          f"({ess / max(1, len(wrows)):.1%} of {len(wrows)})")
    return dict(zip((s for s, _ in wrows), (w for _, w in wrows))), ess


def main():
    import numpy as np
    model = sys.argv[sys.argv.index("--model") + 1] if "--model" in sys.argv else "tilt"
    n = int(float(sys.argv[sys.argv.index("--n") + 1])) if "--n" in sys.argv else 30_000_000
    here = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    catalog = os.path.join(here, "state", "sitemc", "catalog.csv")
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else \
        os.path.join(here, "state", "sitemc", f"weights_{model}.csv")

    U, w, aei = torus_sample(n, model_density(model))
    Umag = np.linalg.norm(U, axis=0)
    print(f"torus hits: {U.shape[1]} of {n:.0f} sampled; "
          f"U<=5 km/s hits: {(Umag <= 5.0).sum()}")
    # Tisserand self-check
    Ut = tisserand_U(*aei)
    err = np.abs(Umag - Ut) / np.maximum(Ut, 0.3)
    print(f"Tisserand check: median rel err {np.median(err):.3f} "
          f"(torus half-widths add O(ds) scatter)")
    keep = Umag <= 6.0
    kde = KDE(U[:, keep], w[keep])
    wmap, ess = weight_catalog(catalog, kde, out)

    # weighted headline numbers
    rows = [r for r in csv.DictReader(open(catalog))
            if r["outcome"] == "impact" and int(r["n_pass"]) >= 2]
    lats = np.array([abs(float(r["lat_deg"])) for r in rows])
    ws = np.array([wmap.get(int(r["seed"]), 0.0) for r in rows])
    order = np.argsort(lats)
    cw = np.cumsum(ws[order]) / ws.sum()
    q = [float(lats[order][np.searchsorted(cw, p)]) for p in (0.25, 0.5, 0.75, 0.9)]
    print(f"weighted |lat| q25/50/75/90 = {q[0]:.1f}/{q[1]:.1f}/{q[2]:.1f}/{q[3]:.1f}")
    hi = [(r, wmap.get(int(r["seed"]), 0.0)) for r in rows
          if r["elev_m"] and float(r["elev_m"]) > 2000]
    def frac(sel):
        tot = sum(w for _, w in hi)
        return sum(w for r, w in hi if sel(r)) / tot if tot else 0.0
    andes = frac(lambda r: -82 <= float(r["lon_e_deg"]) <= -60
                 and -30 <= float(r["lat_deg"]) <= 12)
    tibet = frac(lambda r: 65 <= float(r["lon_e_deg"]) <= 105
                 and 25 <= float(r["lat_deg"]) <= 40)
    print(f"elev>2km weighted share: Andes {andes:.1%}  Tibet/Himalaya {tibet:.1%}")


if __name__ == "__main__":
    main()
