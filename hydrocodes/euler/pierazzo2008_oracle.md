# Pierazzo et al. 2008 (MAPS 43:1917) — extracted oracle numbers for the Al benchmark

Source: PierazzoEtAl2008.pdf (Greg's Downloads, 2026-07-02). This file is the digitized
pass/fail oracle for the Al-on-Al peak-shock-pressure head-to-head. Numbers below are
transcribed from Table 1 + the Benchmark Results text (pp. 1922-1925) and the Summary.

## Benchmark setup (what the paper actually ran)
- 1 km **DIAMETER** Al sphere (radius a = 500 m) into Al half-space; **strength neglected** (pure hydro).
- Vertical impacts at **U = 5 km/s and U = 20 km/s** (plus 45-deg oblique in 3D for 5 codes).
- Measured: peak shock pressure at Lagrangian tracers along a line **downward from the impact
  point**; distance d in km (= d in projectile **diameters** D; our axis units r/a = 2*d/D).
- NOTE: **there is NO 10 km/s case in the paper.** Our 2026-07-01/02 ladder (U=10 km/s) is fine
  for internal convergence/boundary studies but is NOT directly comparable to the paper.

## Table 1 — vertical impacts (the primary oracle)
Contact/compression (CC) region = mean peak shock P from tracers over 0-300 m depth (5 km/s)
or 0-600 m (20 km/s), i.e. d/a = 0-0.6 and 0-1.2.
Decay slope n: log-log LSQ of P ~ d^-n over **~2-5 projectile diameters** (2-5 km, d/a = 4-10).

| Code (EoS) | P_CC 5km/s (GPa) | n 5km/s | P_CC 20km/s (GPa) | n 20km/s |
|---|---|---|---|---|
| ALE3D (LEOS) | 39.2 | 1.35±0.02 | 381.1 | 2.38±0.02 |
| AUTODYN-SPH (Tillotson) | 41.3 | 1.307±0.001 | 396.4 | 2.269±0.003 |
| CTH (ANEOS) | 44.5 | 1.31±0.04 | 380.3 | 2.382±0.009 |
| **iSALE (Tillotson)** | **42.7** | **1.13±0.01** | **371.2** | **2.30±0.02** |
| RAGE (SESAME) | 35.5 | 1.13±0.01 | — | — |
| SOVA (ANEOS tables) | 48.0 | 1.207±0.006 | 411.1 | 2.28±0.01 |
| SPH (Tillotson) | 43.5 | 1.41±0.01 | — | — |
| ZEUS-MP (Tillotson) | 28.4 | 1.19±0.01 | 334.7 | 2.53±0.01 |
| **Mean** | **40.4 ± 6.2** | **1.2 ± 0.1** | **379 ± 26** | **2.3 ± 0.1** |

Spot anchor (text, 5 km/s): at 4 projectile diameters (d/a=8), P = **3.2 ± 0.5 GPa**
(min ~2.5 SPH, max ~4 RAGE).

Inter-code shock-pressure variability (Summary): within **15% (5 km/s)** and **8% (20 km/s)**.
Planar impedance-match limits (text, side note): ~60 GPa (5 km/s), 522 GPa (20 km/s) — the CC
tracer averages sit well BELOW these (spherical geometry; no tracer exactly at impact point).

## Table 1 — 45-deg oblique (3D codes; future 3D-oblique validation anchor)
Slope n measured downward (0 deg) and at 45 deg downrange (d = sqrt((x-R_pr)^2+z^2)):
- 5 km/s: mean n = 1.1±0.1 (both directions). iSALE3D: n/a at 5 km/s in 0-deg col; 1.53±0.02 (0deg listed under 20?) — see Table 1 lower block: iSALE3D 0deg 1.53±0.02, 45deg 2.304±0.009 are the 20 km/s values.
- 20 km/s: mean n(0deg) = 1.5±0.3, n(45deg) = 2.1±0.6. iSALE3D: 1.53±0.02 / 2.304±0.009.
- Max peak-shock P occurs ~1 projectile radius DOWNRANGE of the impact point (Fig. 2).

## Resolution guidance (Table 2 + text)
- **>= 20 cppr**: peak shock P within ~10% of converged; 10 cppr "reasonable".
- iSALE 20->10 cppr difference in the decay region: 10.4±3.9%; 40->20: 4.7±5.1%; 80->40: 1.2±3.0%.
- Under-resolution UNDERESTIMATES peak P (systematic sign).

## Pass/fail criterion for our code (what "within the envelope" means)
At each velocity: (1) P_CC within the code min-max range (28.4-48.0 GPa at 5 km/s;
334.7-411.1 at 20 km/s) and ideally within mean±1sigma; (2) n over d/a=4-10 within the code
range (1.13-1.41 at 5; 2.27-2.53 at 20); (3) 5 km/s spot P(d/a=8) in 2.5-4 GPa.
Closest apples-to-apples single code = **iSALE (Tillotson)** — same EOS family as ours.

## Validation tests in the paper (future methods-section anchors, not the current task)
- Validation #1: 2 mm glass sphere -> water at 4.64 km/s (Boeing quarter-space; Table 3 R(t),d(t)).
  Strengthless + gravity; crater growth vs experiment (codes within ~10%, quarter-space ~-10% bias).
- Validation #2: 6.35 mm Al 2017-T4 sphere -> Al 1100-O and 6061-T6 at ~7 km/s (Prater 1970;
  **Table 4 = full R(t), d(t) experimental data, transcribable**). Tests STRENGTH (von Mises /
  JC / SG). iSALE used von Mises Y=414 MPa for 6061-T6. Codes: radius -6% to -13%, depth +4-12%.
  This is the natural EXPERIMENTAL strength validation for our code (we have von Mises + Y knob).
