"""Boundary-locating daemon for the §2 bounce-off phase diagram (campaign step 2).

Refines the reduced-order scout boundaries (deorbit.impact.scout) by BISECTION: each
increment picks the widest unresolved bracket, evaluates the regime at its midpoint, and
narrows the interval. Anytime/checkpointed like the corridor daemon (run_corridor.py):
  python3 run_bounce.py --daemon            # 24/7 local, refines all brackets
  python3 run_bounce.py --increment 1800    # bounded step (seconds)
State -> state/bounce/ (boundaries.json + status.txt/json). Resumable.

The regime score is PLUGGABLE:
  --score reduced : the analytic scout classifier (instant; validates the daemon machinery
                    and provides a warm start -- but only reproduces the scout).
  --score sph     : a real GPU-SPH run (generic grazing IC -> gpu_ic -> regime classifier).
                    This is the scientific refinement and the daemon's 24/7 work. The SPH
                    classifier thresholds (see classify_regime) are INITIAL and must be
                    calibrated against the first runs before its boundaries are trusted.

This daemon is LOCAL-only (needs the Mac GPU + Metal), unlike the pure-stdlib corridor daemon.
NumPy is imported lazily (only by the SPH score), so the driver itself stays stdlib.
"""
import os
import sys
import json
import time
import math
import argparse
import subprocess

from . import scout

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
_SPH = os.path.join(_REPO, "hydrocodes", "sph")

# --- SPH run config (per evaluation) ---
SPH = dict(D=1000.0,          # impactor diameter (m): 1 km iron
           cppr=8,            # cells-per-projectile-radius (resolution; >=8 for the boundaries)
           rho_imp=7800.0,    # iron
           imp_mat=1,         # material tag for the impactor (iron) in the solver
           tgt_mat=0,         # basalt target half-space
           v_esc=scout.V_ESC["earth"])
U_MELT_IRON = 1.0e6           # incipient shock-melt specific energy for iron (J/kg; ~grazing.py E_melt)
BISECT_TOL = 0.04             # converged when (hi-lo)/mid < this


# ============================ bisection driver (stdlib) ============================
def _load_brackets(statef):
    """Resume from the checkpoint, else seed from the scout."""
    if os.path.exists(statef):
        return json.load(open(statef))
    res = scout.run()
    bks = []
    for i, b in enumerate(res["brackets"]):
        lo, hi = b["sph_bracket"]
        bks.append({"id": i, "vary": b["vary"], "fixed": b["fixed"], "from": b["from"],
                    "to": b["to"], "lo": lo, "hi": hi, "evals": [], "converged": False})
    return {"brackets": bks, "n_eval": 0, "started": time.time()}


def _pick(st):
    """The non-converged bracket with the widest relative interval."""
    best, bw = None, -1.0
    for b in st["brackets"]:
        if b["converged"]:
            continue
        w = (b["hi"] - b["lo"]) / max(1e-9, 0.5 * (b["hi"] + b["lo"]))
        if w > bw:
            best, bw = b, w
    return best


def _theta_vr(b, x):
    """Map a bracket + sample x -> (theta_deg, v/v_esc)."""
    return (x, b["fixed"]) if b["vary"] == "theta" else (b["fixed"], x)


def _update(b, x, reg):
    """Bisection: keep [lo,hi] straddling the from->(not from) boundary."""
    b["evals"].append([x, reg])
    if reg == b["from"]:
        b["lo"] = x
    else:
        b["to"] = reg            # record the actual adjacent regime (may differ from the scout guess)
        b["hi"] = x
    if (b["hi"] - b["lo"]) / max(1e-9, 0.5 * (b["hi"] + b["lo"])) < BISECT_TOL:
        b["converged"] = True


def _save(st, sd):
    os.makedirs(sd, exist_ok=True)
    json.dump(st, open(os.path.join(sd, "boundaries.json.tmp"), "w"), indent=1)
    os.replace(os.path.join(sd, "boundaries.json.tmp"), os.path.join(sd, "boundaries.json"))
    nconv = sum(1 for b in st["brackets"] if b["converged"])
    with open(os.path.join(sd, "status.txt"), "w") as fh:
        fh.write(f"bounce-off boundary daemon\nbrackets {len(st['brackets'])}  converged {nconv}  "
                 f"evals {st['n_eval']}\nupdated {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        for b in st["brackets"]:
            mid = 0.5 * (b["lo"] + b["hi"]); u = "deg" if b["vary"] == "theta" else "v/v_esc"
            fh.write(f"  [{'OK' if b['converged'] else '..'}] {b['from']}->{b['to']} "
                     f"{b['vary']}={mid:.3f} {u} (+-{0.5*(b['hi']-b['lo']):.3f}) @{('vr' if b['vary']=='theta' else 'th')}={b['fixed']}\n")
    json.dump({"name": "bounce-off boundaries", "metric": "converged_brackets",
               "value": nconv, "iterations": st["n_eval"],
               "detail": f"{nconv}/{len(st['brackets'])} boundaries located", "solved": nconv == len(st['brackets']),
               "updated": time.time()}, open(os.path.join(sd, "status.json"), "w"), indent=1)


# ============================ regime scores ============================
def reduced_evaluate(theta, vr):
    return scout.regime(theta, vr, SPH["v_esc"])


def _arange(start, stop, step):
    out = []; x = start
    while (x < stop) if step > 0 else (x > stop):
        out.append(x); x += step
    return out


def build_grazing_ic(theta_deg, v, path, cppr=8, D=1000.0):
    """Write a generic grazing-impact IC (iron sphere on a basalt half-space) in the gpu_ic
    format: header[xlo xhi ylo yhi zlo zhi N] + N*[x y z vx vy vz mass u mat fixed], f8.
    Pure stdlib (array module). Gravity is omitted in the run (the ~1 s traverse is shock/
    strength-dominated), so no lithostatic pre-stress is needed. Returns (dx, N)."""
    import array
    a = D / 2.0; dx = a / cppr; m_imp = SPH["rho_imp"] * dx ** 3
    th = math.radians(theta_deg); vx = v * math.cos(th); vz = -v * math.sin(th)
    Lx, Ly, Lz = max(20e3, 30 * D), max(4e3, 6 * D), max(6e3, 8 * D); wall = 6.0 * dx
    rho_t = 2700.0; mt = rho_t * dx ** 3; tm = float(SPH["tgt_mat"])
    xs = _arange(-0.4 * Lx + dx / 2, 0.6 * Lx, dx); ys = _arange(-Ly / 2 + dx / 2, Ly / 2, dx)
    zs = _arange(-dx / 2, -Lz, -dx)
    rows = array.array("d")
    for X in xs:
        fx = abs(X - xs[0]) < wall or abs(X - xs[-1]) < wall
        for Y in ys:
            fxy = fx or abs(Y - ys[0]) < wall or abs(Y - ys[-1]) < wall
            for Z in zs:
                fixed = 1.0 if (fxy or Z < zs[-1] + wall) else 0.0
                rows.extend((X, Y, Z, 0.0, 0.0, 0.0, mt, 0.0, tm, fixed))
    nslab = len(xs) * len(ys) * len(zs)
    # iron impactor sphere, entering upper-left, heading +x and down at theta
    x0 = -0.35 * Lx; z0 = a + 3 * dx; na = int(math.ceil(a / dx)) + 1; im = float(SPH["imp_mat"]); ni = 0
    for i in range(-na, na + 1):
        for j in range(-na, na + 1):
            for k in range(-na, na + 1):
                if (i * dx) ** 2 + (j * dx) ** 2 + (k * dx) ** 2 <= a * a:
                    rows.extend((i * dx + x0, j * dx, k * dx + z0, vx, 0.0, vz, m_imp, 0.0, im, 0.0)); ni += 1
    N = nslab + ni
    hdr = array.array("d", [xs[0] - dx, xs[-1] + dx, ys[0] - dx, ys[-1] + dx, zs[-1] - dx, z0 + a + 0.5 * Lz, float(N)])
    with open(path, "wb") as f:
        hdr.tofile(f); rows.tofile(f)
    return dx, N


def regime_from_diag(diag, v_esc, v_imp):
    """Pure regime decision from the projectile diagnostics (so cached runs re-label without the
    snapshot). diag carries v_bulk_kms, vbx_kms, disp_kms, z_mean_m, peak_melt (or final_melt)."""
    speed = diag["v_bulk_kms"] * 1e3
    z_mean = diag["z_mean_m"]
    disp = diag["disp_kms"] * 1e3
    melt = diag.get("peak_melt"); melt = diag.get("final_melt", 0.0) if melt is None else melt
    if speed > v_esc:
        return "E"                                          # bulk exceeds escape speed
    # crater/retention: the projectile ends up BURIED below the surface (z_mean below -0.2*D, i.e.
    # the clump centre sank past ~a fifth of a body diameter) AND decelerated to a small fraction
    # of the impact speed -> it will not re-emerge. (The old test compared a velocity FRACTION to an
    # absolute m/s, so a buried-but-slowly-ploughing body was never caught -> crater never fired.)
    if z_mean < -0.2 * SPH["D"] and speed < 0.3 * v_imp:
        return "C"
    # retained ricochet: INTACT only if both spatially coherent (low velocity dispersion) AND not
    # predominantly shock-melted; a mostly-molten clump is DISPERSED even if moving coherently.
    return "I" if (disp <= 0.3 * max(speed, 1.0) and melt <= 0.5) else "D"


def classify_regime(snap_path, v_esc, imp_mat, peak_melt=None, v_imp=None):
    """Regime from the final snapshot's projectile (mat==imp_mat) particles. Pure stdlib.
    gpu_ic snapshot = [int64 N] then N*[x y z vx vy vz u rho h mat] f8 (mat col 9, equal-mass
    projectile particles -> unweighted means). For the intact/dispersed split, melt is judged on
    PEAK shock-melt (peak_melt, the max over the encounter from the run log) when available --
    NOT the final-snapshot energy, which underreports because a ricocheted clump expands and
    cools below the melt threshold even though it was molten at contact."""
    import array
    with open(snap_path, "rb") as f:
        nb = array.array("q"); nb.fromfile(f, 1); N = nb[0]
        d = array.array("d"); d.fromfile(f, N * 10)
    vbx = vby = vbz = z_acc = 0.0; n_proj = 0
    for p in range(N):
        o = p * 10
        if abs(d[o + 9] - imp_mat) < 0.5:
            vbx += d[o + 3]; vby += d[o + 4]; vbz += d[o + 5]; z_acc += d[o + 2]; n_proj += 1
    if n_proj == 0:
        return "C", {"note": "no projectile particles", "v_bulk_kms": 0.0, "disp_kms": 0.0,
                     "vbx_kms": 0.0, "z_mean_m": 0.0, "n_proj": 0}
    vbx /= n_proj; vby /= n_proj; vbz /= n_proj; z_mean = z_acc / n_proj
    speed = math.sqrt(vbx * vbx + vby * vby + vbz * vbz)
    var = 0.0; n_melt = 0
    for p in range(N):
        o = p * 10
        if abs(d[o + 9] - imp_mat) < 0.5:
            var += (d[o + 3] - vbx) ** 2 + (d[o + 4] - vby) ** 2 + (d[o + 5] - vbz) ** 2
            if d[o + 6] > U_MELT_IRON:                      # specific internal energy past incipient iron melt
                n_melt += 1
    disp = math.sqrt(var / n_proj); final_melt = n_melt / n_proj
    diag = {"v_bulk_kms": speed / 1e3, "disp_kms": disp / 1e3, "vbx_kms": vbx / 1e3,
            "z_mean_m": z_mean, "peak_melt": peak_melt, "final_melt": final_melt, "n_proj": n_proj}
    return regime_from_diag(diag, v_esc, v_imp if v_imp is not None else speed), diag


def sph_evaluate(theta, vr, statedir):
    """One real GPU-SPH grazing run -> regime. Builds the IC, runs gpu_ic, classifies."""
    v = vr * SPH["v_esc"]
    rundir = os.path.join(statedir, "runs", f"th{theta:.2f}_vr{vr:.3f}")
    os.makedirs(rundir, exist_ok=True)
    resf = os.path.join(rundir, "result.json")
    if os.path.exists(resf):                                # cached run: re-derive the regime with the CURRENT
        cached = json.load(open(resf))                      # classifier from stored diagnostics (no GPU re-run)
        reg = regime_from_diag(cached, SPH["v_esc"], v)
        if reg != cached.get("regime"):
            cached["regime"] = reg; json.dump(cached, open(resf, "w"), indent=1)
        return reg
    ic = os.path.join(rundir, "ic.bin")
    dx, N = build_grazing_ic(theta, v, ic, SPH["cppr"], SPH["D"])
    t_end = 25e3 / max(v * math.cos(math.radians(theta)), 1.0)   # ~ downrange traverse time
    gpu = os.path.join(_SPH, "gpu_ic")
    if not os.path.exists(gpu):
        raise RuntimeError(f"gpu_ic not built ({gpu}); build it in hydrocodes/sph (./build.sh)")
    # snapshot sparsely (the classifier only needs the final state) -- a 5.9M-particle run at
    # 0.05 s cadence would write ~30 GB; --snap ~0.6*t_end gives ~3 snapshots, deleted below.
    subprocess.run([gpu, ic, f"{dx}", f"{t_end}", "--snap", f"{0.6 * t_end}"], cwd=rundir, check=True,
                   stdout=open(os.path.join(rundir, "run.log"), "w"), stderr=subprocess.STDOUT)
    snaps = sorted([f for f in os.listdir(rundir) if f.startswith("gpu_snap_") and f.endswith(".bin")])
    if not snaps:
        raise RuntimeError(f"no snapshot from gpu_ic in {rundir}")
    peak_melt = 0.0                                         # max projectile melt fraction over the encounter (gpu_ic "melt=X%")
    with open(os.path.join(rundir, "run.log")) as fh:
        for line in fh:
            i = line.find("melt=")
            if i >= 0:
                try:
                    peak_melt = max(peak_melt, float(line[i + 5:line.index("%", i)]) / 100.0)
                except ValueError:
                    pass
    reg, diag = classify_regime(os.path.join(rundir, snaps[-1]), SPH["v_esc"], SPH["imp_mat"], peak_melt, v)
    json.dump({"theta": theta, "vr": vr, "v": v, "dx": dx, "N": N, "t_end": t_end,
               "regime": reg, **diag}, open(os.path.join(rundir, "result.json"), "w"), indent=1)
    for f in os.listdir(rundir):                            # keep result.json + run.log; drop the heavy binaries
        if f.endswith(".bin"):
            os.remove(os.path.join(rundir, f))
    return reg


# ============================ daemon loop ============================
def run():
    ap = argparse.ArgumentParser(description="bounce-off boundary-locating daemon (SPH bisection)")
    ap.add_argument("--daemon", action="store_true")
    ap.add_argument("--increment", type=float, default=None, help="seconds")
    ap.add_argument("--score", choices=["reduced", "sph"], default="sph")
    ap.add_argument("--cppr", type=int, default=None, help="SPH resolution (cells/projectile-radius); default 8")
    ap.add_argument("--statedir", default=os.path.join(
        os.path.dirname(os.path.abspath(sys.argv[0])), "state", "bounce"))
    a = ap.parse_args()
    if a.cppr:
        SPH["cppr"] = a.cppr
    sd = a.statedir
    statef = os.path.join(sd, "boundaries.json")
    st = _load_brackets(statef)
    _save(st, sd)
    deadline = None if a.daemon else time.time() + (a.increment or 1800)
    while True:
        b = _pick(st)
        if b is None:                                       # all converged
            _save(st, sd)
            if not a.daemon:
                return
            time.sleep(60); continue
        x = 0.5 * (b["lo"] + b["hi"]); theta, vr = _theta_vr(b, x)
        try:
            reg = reduced_evaluate(theta, vr) if a.score == "reduced" else sph_evaluate(theta, vr, sd)
            _update(b, x, reg)
        except Exception as e:
            b["converged"] = True                           # park a failing bracket; don't spin
            b.setdefault("error", str(e))
        st["n_eval"] += 1
        _save(st, sd)
        if deadline is not None and time.time() >= deadline:
            return


if __name__ == "__main__":
    run()
