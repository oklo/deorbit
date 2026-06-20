# High-fidelity capture dynamics (REBOUND + ASSIST)

Correct capture/aerobraking dynamics for the km iron monolith: JPL DE440
ephemeris (Sun, planets, Moon, 16 asteroids) + Earth **J2/J3/J4** + Sun J2 + GR,
via REBOUND's **ASSIST** extension in **geocentric** mode, with atmospheric
**drag** added by operator-splitting (ASSIST integrates gravity over a small
sub-step; the drag Δv kick is applied between sub-steps near periapsis).

This is the rigorous counterpart to the fast pure-stdlib corridor sweep
(`../physics.py`, `../deorbit.py`). It exists to answer what the 2-body model
can't: **does a "capture" survive the Moon?** (first apogees are typically
~2.4×10⁵ km, deep in the lunar zone), and how J2 steers the long aerobraking.

## Setup (local only — NOT in the cloud routine, which stays pure stdlib)

```bash
# from the repo root
python3 -m venv .venv
.venv/bin/pip install -r highfidelity/requirements.txt

# ephemeris binaries (large; gitignored). Download once into highfidelity/data/:
mkdir -p highfidelity/data && cd highfidelity/data
curl -LO https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440
curl -LO https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp
```

## Validate

```bash
cd highfidelity && ../.venv/bin/python validate_hifi.py
```

Checks: Kepler energy closure (full self-consistent forces), Moon ephemeris
distance, **J2 nodal regression isolated by differencing** vs. the analytic
rate, and per-pass drag Δv vs. the validated stdlib integrator.

## Files

- `hifi.py` — ephemeris/sim setup, geocentric-frame helpers, split-drag, and
  `simulate_passage` (full-forces + drag, comparable to `physics.simulate_passage`).
- `validate_hifi.py` — the ground-truth checks above.
- `data/` — ephemeris binaries (gitignored).

## Notes / gotchas

- ASSIST's geocentric mode is **non-inertial**: Earth's central monopole rides
  on the `PLANETS` force term, so you cannot isolate a "pure two-body" subset by
  dropping forces. Use the full set; isolate an effect (e.g. J2) by *differencing*.
- ASSIST owns rebound's force hook, hence drag-by-splitting rather than a
  rebound `additional_forces` callback.
- Next: promote `simulate_passage` to a multi-pass aerobraking sequence
  (coast→passage→coast) to test lunar survival of captures and the J2-steered
  terminal great circle (the Chimborazo/Cayambe site-selection result).
