#!/usr/bin/env python3
"""Robust OFFLINE depth-diameter measurement from dumped crater surface profiles.

The `crater` driver dumps z_s(r)-z0 (surface height relative to the far-field
datum) per radial cell. This reads those profiles and extracts the crater
geometry that the in-loop diagnostic cannot do reliably:
  - floor depth d      = -min(smoothed dz)            (deepest point below datum)
  - apparent diameter  D_app = 2 * r_cross            (where the wall returns to datum)
  - rim crest          = max dz beyond the floor      (ejecta/structural uplift)
  - d/D                = d / D_app                     (the calibration observable)

Robustness (the unsettled profiles showed isolated +350 m spikes on the axis and
at scattered cells, plus concentric ripples): we median-filter (kills isolated
single-cell spikes) then lightly box-smooth before measuring, and locate the
datum crossing by linear interpolation scanning OUTWARD from the floor.

stdlib only. Usage:
    analyze_crater.py [profile.txt ...]      # specific files
    analyze_crater.py                        # all of state/profiles/*.txt
    analyze_crater.py --csv out.csv          # also write a CSV table
"""
import glob, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def read_profile(path):
    r, dz = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            a, b = line.split()[:2]
            r.append(float(a)); dz.append(float(b))
    return r, dz


def median3(y):
    out = list(y)
    for i in range(1, len(y) - 1):
        out[i] = sorted(y[i - 1:i + 2])[1]
    return out


def box3(y):
    out = list(y)
    for i in range(1, len(y) - 1):
        out[i] = (y[i - 1] + y[i] + y[i + 1]) / 3.0
    return out


def analyze(path, datum_tol=None):
    r, dz = read_profile(path)
    if len(r) < 5:
        return None
    dx = r[1] - r[0]
    tol = datum_tol if datum_tol is not None else 0.25 * dx   # "back to datum" band, ~quarter cell
    s = box3(median3(dz))                                     # despike then lightly smooth

    # floor: deepest (most negative) smoothed cell
    ifloor = min(range(len(s)), key=lambda i: s[i])
    depth = -s[ifloor]
    rfloor = r[ifloor]

    # apparent radius: scan OUTWARD from the floor to where the wall rises back through the datum (z=0)
    r_cross = r[-1]
    for i in range(ifloor, len(s) - 1):
        if s[i] < 0.0 <= s[i + 1]:                       # sign change negative -> non-negative
            frac = (0.0 - s[i]) / (s[i + 1] - s[i]) if s[i + 1] != s[i] else 0.0
            r_cross = r[i] + frac * (r[i + 1] - r[i])    # linear-interpolated datum crossing
            break
    else:
        # never returns to datum within the box: fall back to the outermost cell still below -tol
        for i in range(len(s) - 1, ifloor, -1):
            if s[i] < -tol:
                r_cross = r[i]; break

    D_app = 2.0 * r_cross
    dD = depth / D_app if D_app > 0 else 0.0

    # rim crest: highest smoothed cell at r > rfloor (uplift above datum)
    rim_h, rim_r = 0.0, 0.0
    for i in range(ifloor, len(s)):
        if s[i] > rim_h:
            rim_h = s[i]; rim_r = r[i]

    return dict(name=os.path.basename(path), dx=dx, depth=depth, rfloor=rfloor,
                D_app=D_app, dD=dD, rim_h=rim_h, rim_r=rim_r, npts=len(r))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    csv_out = None
    if "--csv" in sys.argv:
        csv_out = sys.argv[sys.argv.index("--csv") + 1]
    files = args if args else sorted(glob.glob(os.path.join(HERE, "state", "profiles", "*.txt")))
    if not files:
        print("no profiles found"); return

    rows = []
    print(f"{'profile':<28} {'dx':>5} {'depth_km':>9} {'D_app_km':>9} {'d/D':>7} {'rim_km':>7} {'rim_r_km':>9}")
    for p in files:
        a = analyze(p)
        if not a:
            print(f"{os.path.basename(p):<28} (too short)"); continue
        rows.append(a)
        print(f"{a['name']:<28} {a['dx']:>5.0f} {a['depth']/1e3:>9.3f} {a['D_app']/1e3:>9.3f} "
              f"{a['dD']:>7.3f} {a['rim_h']/1e3:>7.3f} {a['rim_r']/1e3:>9.3f}")

    if csv_out and rows:
        import csv as _csv
        with open(csv_out, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=["name", "dx", "depth", "rfloor", "D_app", "dD", "rim_h", "rim_r", "npts"])
            w.writeheader()
            for a in rows:
                w.writerow(a)
        print(f"\nwrote {csv_out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
