# deorbit — research project operating guide

This is one of Greg's long-horizon **background research investigations**. Run it
with the same autonomy + "plug away over time" pattern already used by the
`malbolge` (map6 hero run) and `modes` (disk instabilities) projects. This file
is auto-loaded, so you (the Claude Code session running here) already have the
context below — use it.

## What you already have (inherited automatically)

- **Broad tool permissions** + `acceptEdits` default mode (from `~/.claude/settings.json`).
- **`gh` CLI authenticated** as `oklo` (`repo`, `workflow` scopes) → you can
  `gh repo create`, push, open PRs without setup.
- **GitHub App connected** for cloud → scheduled cloud routines can clone/push.
- **Skills**: `schedule` (cloud cron routines), `update-config`, `run`, etc.
- **Shared framework** at `~/investigations/` — reuse it; do NOT reinvent:
  - `template/engine.py` — reusable anytime-search engine (checkpointing, status,
    `--daemon`/`--increment N`, depth growth, solved-detection). You write only
    `propose(rng, depth)` and `score(candidate) -> (value, solved, artifact)`.
  - `bin/new-investigation <slug> "<goal>"` — scaffold + register a study.
  - `bin/board` — terminal dashboard across all studies.
  - `dashboard/server.py` — read-only web dashboard at http://localhost:8787
    (auto-discovers figures/recent files for "static" projects).
  - `daemons.json` — daemons the **launchd supervisor** keeps alive 24/7.
  - `RECOVERY.md` — reboot/OS-upgrade runbook.

## Operating model (the budget principle)

**Local compute is free and unlimited; tokens are the scarce, paced resource.**
So structure the work in two layers:

1. **Free local daemon** — grinds 24/7 on this laptop, self-verifying, checkpointing.
   Zero tokens. This does the heavy lifting.
2. **Sparse cloud routine** — a daily scheduled cloud agent runs a bounded
   increment, commits/pushes progress, and steers when stuck. This is the only
   budgeted cost; keep it short (~25 min) and **stagger its time** away from the
   malbolge routine (which runs 06:00 ET / 10:00 UTC) — e.g. use **07:00 ET**.

Budget governance is heuristic (token meter isn't visible here): keep concurrent
daily routines few and short, treat ~half the monthly allocation as the steady
ceiling, and leave headroom for urgent interactive work. There are already 1–2
active investigations — add this one, but don't pile on more routines without
reason.

## Plug-in checklist (do this once the research goal is defined)

1. **Define the goal** with Greg: what is the win condition, and what's the
   computational core (search? simulation sweep? optimization? symbolic
   derivation)? Pick the layer that fits: an `engine.py`-style `propose/score`
   search, or a custom daemon that periodically advances + checkpoints.
2. **Git**: `git init`; `gh repo create deorbit --private --source=. --push`.
3. **Local daemon**: either
   `~/investigations/bin/new-investigation deorbit "<goal>"` (then edit
   `~/investigations/deorbit/study.py`), or build a project-local daemon here
   that supports `--daemon` and `--increment N` and checkpoints to `state/`.
   Verify a short `--increment` run checkpoints and resumes.
4. **Register for 24/7 supervision**: add an entry to
   `~/investigations/daemons.json` with a unique `match` signature, `cwd`, and
   `cmd` (the supervisor re-reads it every cycle — no restart needed). Confirm
   it appears `live` in `~/investigations/bin/board`.
5. **Dashboard**: add an entry to `~/investigations/dashboard/projects.json`
   (`"investigation"` for a live metric, `"static"` for a passive project) so it
   shows at http://localhost:8787. Config reloads per page view — no restart.
6. **Cloud routine**: use the `/schedule` skill to create a daily routine
   (staggered time, e.g. `0 11 * * *` = 07:00 ET) that pulls this repo, runs one
   bounded `--increment`, commits/pushes state, and opens a PR if it hits the win
   condition. Trigger one test run to confirm the cloud path end-to-end.
7. **Resumability**: write a `PICKUP.md` here so an interactive session can be
   restarted after a reboot (see `modes/PICKUP.md` for the pattern). Daemons +
   cloud routine survive reboot automatically (launchd + server-side); only
   interactive sessions are restarted by hand.

## Conventions

- Keep daemon dependencies to the Python stdlib (+ a local repo on `sys.path`)
  so a fresh cloud clone can run with no setup.
- Checkpoint frequently; never assume the process won't be killed mid-run.
- Ground-truth every claimed result with a real evaluator/oracle before reporting
  it solved — distinguish "diagnostic says X" from "verified X".
- Read-only when observing other projects; never touch another investigation's
  working tree or git index.
