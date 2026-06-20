# deorbit

**Ultra-oblique de-orbit and grazing impact of a km-scale monolithic iron asteroid.**

A long-horizon background research investigation (see `CLAUDE.md` for the
operating model). The novel problem: the flight regime of a kilometre-scale,
dense, *strong* iron monolith brought down by atmospheric drag at grazing
incidence — a gap between meteor-entry physics (small ablating bodies) and
impact hydrocodes (terminal cratering at steep angles).

## The scientific chain

1. **Capture** — pure gravity can't capture (energy conserved); tidal
   disruption→ring needs a rubble pile (a monolith has strength); binary
   exchange (Triton mechanism) fails by orders of magnitude at km scale
   (v_bin ≪ v∞). The **only viable route for an intact iron monolith is
   aerocapture**: a grazing periapsis pass sheds the small hyperbolic excess.
   But beta = m/(Cd·A) ≈ 5×10⁶ kg/m², so capture is possible only in a narrow
   low-v∞ corridor — *if at all*. **This is the crux open question.**
2. **Circularization** — aerobraking lowers apoapsis over many passes;
   inclination is barely damped; there is no stable ring (periapsis stays in the
   atmosphere), so the endpoint is spiral-in to a grazing impact.
3. **Terminal grazing entry** — ultra-shallow (~1–3°), ~8 km/s, the under-studied
   regime: an elongated gouge / ricochet, not a round crater.
4. **Site selection** — the trajectory is a conic about Earth's *center*, so the
   selecting metric is **geocentric radius**, where the ~21 km equatorial bulge
   dwarfs ~9 km of topography. The body is steered to the point farthest from
   Earth's center: **Chimborazo / Cayambe** (equatorial Andes). Apparently novel.

Three terminal-target cases (Phase 2): **Cayambe** (high equatorial rock),
**deep ocean**, **Amazon basin** (wet sediment). Phobos/Deimos = a later
planet swap.

## Layers (per the budget principle)

- **Free local daemon** — `deorbit.py`, the capture-corridor sweep. Pure stdlib.
- **Sparse cloud routine** — daily bounded `--increment`, commits/pushes state.

## Run

```bash
python3 validate.py                 # ground-truth physics checks
python3 deorbit.py --increment 60   # bounded sweep step (seconds)
python3 deorbit.py --daemon         # 24/7 local, refines the corridor grid
```

Outputs land in `state/`: `corridor.json` (full grid + summary),
`status.txt` / `status.json` (board/dashboard).

## Files

- `physics.py` — atmosphere, single-passage integrator, outcome classifier.
- `deorbit.py` — (v_inf, periapsis) corridor sweep daemon.
- `validate.py` — Kepler closure, per-pass dv vs analytic column, regime sanity.
- `PICKUP.md` — restart prompt for an interactive session after a reboot.

`iSALE2D/` (gitignored) is the gated shock-physics code, kept for Phase-2
terminal-impact runs; the public repo is a stub — the real code is the
iSALE-Dellen release tarball.

## Status

Phase 1a: capture-corridor map (does km iron get captured, and where is the
corridor?). Phase 1b: aerobraking-to-impact + 3D site selection. Phase 2:
terminal oblique-impact hydrocode for the three targets.
