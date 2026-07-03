# Manuscript: The Capture and Deorbit of Kilometer-Scale Iron Asteroids

Target journal: Planetary Science Journal (AAS). First draft 2026-07-03.

- `ms.tex` — AASTeX v6.3.1 manuscript. Scaffolded (unwritten) material is
  marked in-text with *[To be completed: ...]* and `%%TODO` comments:
  the corridor-map figure, the event-rate integral (Sec 6.1), the
  Ordovician application (Sec 6.2), and the second-impact subsection of
  Sec 5 (awaiting the SPH convergence study).
- `references.bib` — entries marked `[CHECK]` were cited from memory and
  must be verified against ADS before submission.
- Figures resolve from `../figures/` (committed): `sitemc_traj`,
  `sitemc_traj_plane`, `sitemc_flashlight_rw`, `bounce_off_phase_diagram`.
- Compile: `tectonic ms.tex` (class + bst included here from CTAN; tectonic
  fetches everything else on demand). Output ms.pdf is untracked (24 MB).

All quantitative statements trace to committed artifacts:
`docs/capture_feasibility.md` (Sec 2), `docs/highfidelity_findings.md`
(Sec 3), `results/sitemc_summary.json` + `docs/site_mc_findings.md`
(Sec 4), `results/bounce_off_phase_diagram.json` (Sec 5).
