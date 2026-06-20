# Pickup prompt — deorbit

Paste the block below into a fresh Claude Code session in this directory to resume.

---

Continue the **deorbit** investigation in this directory (read `CLAUDE.md` and
`README.md` first). Goal: simulate the de-orbit and ultra-oblique grazing impact
of a **km-scale monolithic iron asteroid**, brought down by atmospheric drag.

The physics chain and the reasoning behind it (capture mechanisms ruled in/out,
the geocentric-radius / equatorial-bulge site-selection result, the literature
gap) are summarized in `README.md` and in this session's memory files.

**Where things stand / next steps:**

- **Phase 1a (current):** `deorbit.py` is a pure-stdlib capture-corridor sweep
  over (v_inf, periapsis altitude) for a 1 km iron monolith. Each cell runs
  `physics.simulate_passage` and is classified impact / capture / skip. The crux
  question it answers: **does a km iron body get captured (eps<0) in a low-v∞
  corridor, or does it just skip out / impact directly?** Check `state/corridor.json`
  and `state/status.txt` for the current corridor map.
  - Validate first: `python3 validate.py` (Kepler closure, per-pass dv vs analytic
    column, regime sanity). Then `python3 deorbit.py --increment 60`.
  - The daemon should be running 24/7 via the launchd supervisor (see
    `~/investigations/daemons.json`, entry name `deorbit`). Confirm with
    `~/investigations/bin/board`.

- **Phase 1b (next):** promote the single-passage map into an iterated
  periapsis-passage map → aerobraking circularization → terminal descent, in 3D
  (add inclination + J2), and compute the impact site as the max-geocentric-radius
  intercept along the final great circle (validate the Chimborazo/Cayambe result).

- **Phase 2:** terminal oblique-impact hydrocode (iSALE3D, gated, or SPH) for the
  three targets — Cayambe rock / deep ocean / Amazon sediment — fed by Phase-1b
  terminal conditions. Needs the iSALE-Dellen tarball (not the stub repo).

Conventions: daemon stays pure stdlib (no numpy) so a fresh cloud clone runs with
no setup. Checkpoint frequently. Ground-truth every "solved/captured" claim
against the integrator/oracle, not a diagnostic.
