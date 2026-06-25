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

- **Free local daemon** — `run_corridor.py`, the capture-corridor sweep. Pure stdlib.
- **Sparse cloud routine** — daily bounded `--increment`, commits/pushes state.

## Run

```bash
pytest                                  # ground-truth physics checks (stdlib corridor + hi-fi if installed)
python3 run_corridor.py --increment 60  # bounded sweep step (seconds)
python3 run_corridor.py --daemon        # 24/7 local, refines the corridor grid
```

The Python code is the `deorbit` package under `src/`; the core (corridor, site
selection, impact) is **pure stdlib** so the daemon and cloud routine run with no
install. The REBOUND/ASSIST high-fidelity layer is an optional extra:
`pip install -e ".[highfidelity]"`. Outputs land in `state/`: `corridor.json`
(full grid + summary), `status.txt` / `status.json` (board/dashboard).

## Files

- `src/deorbit/corridor/physics.py` — atmosphere, single-passage integrator, outcome classifier.
- `src/deorbit/corridor/sweep.py` — (v_inf, periapsis) corridor sweep; `run_corridor.py` is the root entry point.
- `tests/` — `pytest` ground-truth checks (Kepler closure, per-pass dv vs analytic, regime sanity, hi-fi).
- `PICKUP.md` — restart prompt for an interactive session after a reboot.

`iSALE2D/` (gitignored) is the gated shock-physics code, kept for Phase-2
terminal-impact runs; the public repo is a stub — the real code is the
iSALE-Dellen release tarball.

## High-fidelity dynamics (`highfidelity/`)

The rigorous counterpart to the stdlib sweep: REBOUND + **ASSIST** (JPL DE440
ephemeris, Earth J2/J3/J4, Sun/Moon/planets, GR) in geocentric mode, with drag
by operator-splitting. Validated (Kepler closure, J2 nodal regression vs.
analytic, drag Δv vs. the stdlib model). Local-only (needs a venv + ephemeris
download); the cloud routine stays pure stdlib. See `highfidelity/README.md`.

## Status

Phase 1a: capture-corridor map (does km iron get captured, and where is the
corridor?) — stdlib daemon, running. Phase 1b: high-fidelity engine **built +
validated**; next is the multi-pass aerobraking sequence to test which captures
survive the Moon and the J2-steered terminal great circle (Chimborazo/Cayambe).
Phase 1c: site selection (`src/deorbit/site_selection.py`) — the equatorial bulge dominates
geocentric radius, so a slowly-decaying grazing body is steered to the global
geocentric-radius maximum in its latitude band: **Chimborazo/Cayambe** (validated
against the known "farthest point from Earth's center" ranking). Phase 2:
terminal oblique-impact hydrocode for the three targets.
