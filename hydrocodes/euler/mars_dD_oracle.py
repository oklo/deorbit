#!/usr/bin/env python3
"""Mars fresh-crater depth-diameter ORACLE for the Phase-3 AF/strength calibration.

The euler `crater` driver is calibrated (acoustic-fluidization + rock-strength params)
to reproduce the observed Mars depth-diameter (d-D) dependence. The METHOD we mirror is
Wuennemann & Ivanov 2003 (Planet. Space Sci. 51:831) -- the AF block model; the DATA we
calibrate against (since our sims use Mars gravity g=3.71) is the MOLA fresh-crater d-D
relation. The canonical published fits:

GARVIN, SAKIMOTO & FRAWLEY (2003), "Craters on Mars: Global Geometric Properties from
Gridded MOLA Topography", Sixth Int'l Conf. on Mars, abstract #3277, Table 1 (d, D in km):
    simple  craters:  d = 0.21 * D**0.81     (fit to MOLA profile data)
    complex craters:  d = 0.36 * D**0.49     (DEM data; complex range ~7-100 km)
    rim height: simple h=0.04 D**0.31, complex h=0.02 D**0.84
  Their worked fresh example: D=6.2 km, d=1385 m -> d/D=0.2234, noted ~400 m DEEPER than
  the GLOBAL fit (which folds in degraded/infilled craters). So a single PRISTINE impact
  (what we simulate) should sit ABOVE the global simple fit -- nearer d/D ~ 0.20-0.23 for
  small simple craters -- not on it.

SIMPLE->COMPLEX TRANSITION on Mars: ~7 km (Garvin et al. 2003); 6.9 km global mean
(Robbins & Hynek 2012, JGR 117, E06001; range ~5-11 km regionally).

Notes / caveats for the comparison:
  - Our impactors a=300..3000 m at 12 km/s give D_app ~ a few -> few tens of km, spanning
    the transition -- exactly the regime W&I/Garvin constrain.
  - The GLOBAL fit is a lower envelope for fresh craters; we also expose a "fresh" target
    band (d/D ~ 0.20-0.23 simple) for the pristine-impact comparison.
  - These are Mars EMPIRICAL fits; W&I 2003's own modelled curve was tuned to the Moon
    (transition ~15 km) and gravity-scaled -- Mars (Garvin) is the right oracle here.

stdlib only.
"""
import csv, os, sys

D_TRANS = 7.0     # km, simple->complex transition (Garvin et al. 2003)
D_TRANS_RH = 6.9  # km, Robbins & Hynek 2012 global mean


def d_simple(D):   # km -> km, Garvin Table 1 simple-crater fit
    return 0.21 * D ** 0.81


def d_complex(D):  # km -> km, Garvin Table 1 complex-crater fit (~7-100 km)
    return 0.36 * D ** 0.49


def d_garvin(D):   # piecewise global fit at the 7 km transition
    return d_simple(D) if D < D_TRANS else d_complex(D)


def dD_garvin(D):  # global depth/diameter ratio
    return d_garvin(D) / D


def dD_fresh_band(D):
    """Fresh (pristine) target band as (low, high) d/D. Fresh simple craters run ~0.20-0.23
    (Garvin's example 0.2234 at D=6.2); complex follow the global complex fit but a bit deeper.
    Low edge = the global Garvin fit; high edge = ~1.2x it (the fresh offset the paper reports)."""
    lo = dD_garvin(D)
    hi = 1.2 * lo if D >= D_TRANS else max(0.21, 1.2 * lo)   # small simple craters cap near ~0.21
    return (lo, hi)


def rim_simple(D):  return 0.04 * D ** 0.31   # km, Garvin rim-height fits (for the rim-uplift cross-check)
def rim_complex(D): return 0.02 * D ** 0.84


def compare(rows):
    """rows: iterable of (label, D_km, d_km). Print sim d/D vs the Garvin oracle + residual."""
    print(f"{'label':<22} {'D_km':>6} {'d_km':>6} {'d/D_sim':>8} {'d/D_oracle':>10} {'fresh_band':>14} {'morph':>8} {'ratio':>6}")
    for label, D, d in rows:
        if D <= 0:
            continue
        dDsim = d / D
        dDor = dD_garvin(D)
        lo, hi = dD_fresh_band(D)
        morph = "simple" if D < D_TRANS else "complex"
        ratio = dDsim / dDor if dDor > 0 else 0.0
        print(f"{label:<22} {D:>6.2f} {d:>6.2f} {dDsim:>8.3f} {dDor:>10.3f} "
              f"{lo:>6.3f}-{hi:<6.3f} {morph:>8} {ratio:>6.2f}")


def _load_analyzer_csv(path):
    """Read analyze_crater.py --csv output (name, depth, D_app in metres) -> (label, D_km, d_km)."""
    out = []
    with open(path) as f:
        for r in csv.DictReader(f):
            label = r["name"].replace(".txt", "")
            D = float(r["D_app"]) / 1e3
            d = float(r["depth"]) / 1e3
            out.append((label, D, d))
    return sorted(out, key=lambda t: t[1])


def main():
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        compare(_load_analyzer_csv(sys.argv[1]))
    else:
        # self-demo: print the oracle curve across the simulated size range
        print("Garvin et al. 2003 Mars fresh-crater d-D oracle (d,D in km):")
        for D in (1, 2, 3, 5, 7, 10, 15, 20, 30, 50):
            lo, hi = dD_fresh_band(D)
            print(f"  D={D:>3} km  d={d_garvin(D):.3f} km  d/D={dD_garvin(D):.3f}  "
                  f"fresh~{lo:.3f}-{hi:.3f}  ({'simple' if D < D_TRANS else 'complex'})")
        print("\nUsage: mars_dD_oracle.py state/dD_settled.csv   # compare sim profiles to the oracle")


if __name__ == "__main__":
    main()
