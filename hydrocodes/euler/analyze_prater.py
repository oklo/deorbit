#!/usr/bin/env python3
"""Compare a prater growth file (t_us R_cm D_cm) against the Prater-1970 Table 4 data.

Usage: analyze_prater.py [growth_file ...]   (default: prater_growth.txt)
For each experimental point (t_exp, X_exp) it linearly interpolates the sim curve at t_exp
and reports mean/RMS relative deviation, split early (<15 us, growth) / late (>=25 us, plateau).
Stdlib only (no numpy on this box). Oracle: prater1970_oracle.md (6061-T6 block).
"""
import sys, re, os

HERE = os.path.dirname(os.path.abspath(__file__))

def read_oracle():
    """Parse the 6061-T6 radius/depth tables out of prater1970_oracle.md."""
    txt = open(os.path.join(HERE, "prater1970_oracle.md")).read()
    m = re.search(r"## Table 4 — Al 6061-T6.*?Radius:\s*t_us R_cm\n(.*?)\n\nDepth:\s*t_us D_cm\n(.*?)\n\n## Table 4 — Al 1100-O",
                  txt, re.S)
    if not m:
        sys.exit("oracle parse failed (6061-T6 block)")
    def parse(block):
        out = []
        for line in block.strip().splitlines():
            p = line.split()
            if len(p) == 2:
                out.append((float(p[0]), float(p[1])))
        return sorted(out)
    return parse(m.group(1)), parse(m.group(2))

def read_growth(path):
    ts, Rs, Ds = [], [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) == 3:
            ts.append(float(p[0])); Rs.append(float(p[1])); Ds.append(float(p[2]))
    return ts, Rs, Ds

def interp(ts, ys, t):
    if t <= ts[0]: return ys[0]
    if t >= ts[-1]: return ys[-1]
    for i in range(1, len(ts)):
        if ts[i] >= t:
            f = (t - ts[i-1]) / (ts[i] - ts[i-1])
            return ys[i-1] + f * (ys[i] - ys[i-1])
    return ys[-1]

def stats(exp, ts, ys, label, tmax):
    rows = [(t, x, interp(ts, ys, t)) for t, x in exp if 3.0 <= t <= tmax]
    if not rows:
        print(f"  {label}: no experimental points in range"); return
    for name, sel in (("early(3-15us)", lambda t: t < 15), ("late(>=25us)", lambda t: t >= 25),
                      ("all(3us-)",     lambda t: True)):
        r = [(s - x) / x for t, x, s in rows if sel(t)]
        if not r: continue
        mean = sum(r) / len(r)
        rms = (sum(v * v for v in r) / len(r)) ** 0.5
        print(f"  {label} {name:14s} n={len(r):2d}  mean {100*mean:+6.1f}%  rms {100*rms:5.1f}%")

def main():
    files = sys.argv[1:] or ["prater_growth.txt"]
    expR, expD = read_oracle()
    for f in files:
        ts, Rs, Ds = read_growth(f)
        if not ts:
            print(f"{f}: empty"); continue
        tmax = ts[-1]
        print(f"{f} (to {tmax:.0f} us):")
        stats(expR, ts, Rs, "R", tmax)
        stats(expD, ts, Ds, "D", tmax)

if __name__ == "__main__":
    main()
