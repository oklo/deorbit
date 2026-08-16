# PICKUP — shear-strain damage (CMI-2004): breaking the intact-rock arch

Self-contained handoff for a fresh session. **This is the entry point for the
next major effort on the euler code**, following the completed Option C arc
(PICKUP-optionC.md = its full evidence trail + verdict; skim the banner only,
do not re-derive). Repo: github.com/oklo/deorbit, code in `hydrocodes/euler/`
(CPU FP64 oracle `hydro_cpu.cpp` = the engine; GPU port deferred until the
physics closes). Conventions + build lines: bottom of PICKUP-phase3.md.
Gates: `./gates.sh quick` = 44 PASS at the time of writing (July 2026); CI
(`master gates` workflow) green.

## 1. The problem (one paragraph)

Option C (2026-07-07..11) built the full literature AF stack — implicit
Bingham viscosity (eta unbounded, energy-exact), acoustic dP/(rho*c_s) seeding,
W&I overburden relief (P_eff = P - pvib), block-model vibration-energy
REGENERATION (fluidization sustained while flow lives, self-arresting), pvib
saturation — and systematically eliminated the whole family at the a=1000/3000
anchors: no corner of {strength form x seeding x one-shot/regen x amplitude x
TDEC regime x resolution} produces complex-crater FLOOR REBOUND. The craters
are now STABLE at transient depth (d/D ratio 3.5-4x Garvin; wall slumping and
widening reproduced; drain/stew/blowup all fixed), held open by a 2D stress
arch through the off-axis ring. THE DIAGNOSED ROOT CAUSE: our Grady-Kipp
damage grows ONLY UNDER TENSION, so the compressed sub-crater volume never
damages (D ~ 0.04-0.10 in every run) and the arch is INTACT Lundborg rock
(~2e8 Pa at ring pressures). Community codes (iSALE) enter collapse with the
target pervasively D=1 — cohesionless mu_d rubble — via SHEAR/plastic-strain
damage (Collins, Melosh & Ivanov 2004, MAPS 39:217). Their AF only has to beat
damaged friction; ours must beat intact rock. Fix: implement CMI-2004
plastic-strain damage accumulation.

## 2. Evidence trail (all committed/ledgered; do not re-derive)

- Option C verdict + failure-mode matrix (pipe / under-fluidized lock /
  arch-lock / regional meltdown): PICKUP-optionC.md VERDICT banner (2026-07-11).
- Arch evidence: radial af profile at a=3000 — saturated pvib=P lens under the
  floor (af=0.50 cap, SUSTAINED by regen for hours), ring at r~20-30 km /
  z~10-20 km at full lithostatic P with pvib<<P -> friction ~70+ MPa intact.
- D(mean, frac>0.95) = 0.03-0.10, 0.01-0.07 in EVERY complex run (the logs) —
  the smoking gun for tension-only damage.
- Resolution exonerated: cppr8 == cppr5 at a=1000, floors dead flat from
  t=100 s (dynamic phase included). Amplitude exonerated: C_ACT=1.0 melts the
  regional crust (surface sinks 8-10 km domain-wide) while the arch holds.
- Local (gitignored) run evidence: state/optc/ (~25 runs, tags encode params);
  figures/checks/optc_*.png = tracer-frame morphology.

## 3. The fix — design

3.1 PLASTIC-STRAIN ACCUMULATION (the core).
  - The von Mises radial return already computes the excess (vm - Y); the
    plastic strain increment of a return is d_eps_p = (vm - Y)/(3G). Accumulate
    per cell into a new advected field eps_p (advect like D/vib: v.grad upwind
    in Lop, no cylindrical source). Corrector-final only (reg-flag pattern from
    vonmises(g,bool reg) — predictor must not double-count; that machinery
    exists).
  - Damage: D = max(D_tensile_gradykipp, min(1, eps_p/eps_f(P))) with the
    CMI-2004 pressure-dependent failure strain. PULL THE EXACT eps_f(P) FORM
    AND CONSTANTS FROM THE PAPER (Collins, Melosh & Ivanov 2004; basalt-like
    defaults; eps_f is ~0.01 at low P growing with P — do not trust memory,
    read the paper; literature_pdfs/ may have it, else fetch).
  - Keep the existing D semantics downstream (Yb blend intact->damaged in
    vonmises ROCK path) — that is the whole point: D->1 under the crater flips
    the arch from intact Lundborg to cohesionless mu_d friction, WHICH THE
    EXISTING AF RELIEF + REGEN CAN THEN BEAT.

3.2 GATES (before physics runs; extend gates.sh; feature-off = 44 unchanged).
  a. Bit-identity: new path behind a flag (env or arg; crater default ON only
     after the gates pass); all 44 existing PASS with it off.
  b. Unit gate: uniform-shear cell driven past yield accumulates eps_p at the
     analytic radial-return rate; D crosses 1 at eps_f(P) for two P values on
     the CMI curve; strength flips to the damaged branch (extend the friction
     gate pattern).
  c. Crater smoke gate (qualitative): a=1000 short run -> D(mean) under the
     crater PERVASIVE (>0.7 within ~1.5 crater radii) vs the old ~0.1.

3.3 ANCHOR TEST, then calibration.
  - Rerun the two stable anchors at the adopted AF defaults (crater args:
    ROCK=1, Y_D0=5e6, VAMB=0.27, Y_B=3e6, VISC_IMP=1, VOID_CONS=0,
    C_REG=0.1-0.3, AF_SAT=2; TDEC = 2x crater timescale: 900/1400 s;
    eta=1e9/3e9): `run_optc.sh` does all of this (see its header; CREG/SAT/CACT
    env knobs). THE QUESTION: does the floor finally rebound?
  - RISKS to watch: (i) pervasive damage may re-open the deep drain (the
    column that friction held is now rubble) — the guardrails are the pvib
    saturation cap, Y_D0, Y_B, and implicit viscosity; if it drains, that is
    now a CALIBRATABLE regime (raise Y_D0/eta), not a model wall. (ii) the
    a=300 simple-crater guard (target ratio ~1.5) may over-soften — Y_D0 is
    the lever (the "creep knob").
  - If rebound appears: FULL SWEEP vs Garvin (mars_dD_oracle.py; harvest
    analyze_crater.py -> ratio table), TDEC-doubling insensitivity check at one
    anchor (the Option C acceptance criterion, finally testable), anchors
    a=300/1000/2000/3000. Deliver the (eta, C_REG, Y_D0) table to Greg.

## 4. Loose ends / gotchas (inherited)

- VOID_CONS (crater arg17) exists but stays OFF: mass-only transfer is a
  downward pump (a=1000 blowout, 2026-07-09). Momentum-consistent transfer is
  a parked follow-up.
- S-channel energy pump: long full-strength pure-shear-wave setups pump energy
  (pre-existing, found 2026-07-10); keep gate setups yielding-soft. Craters
  unaffected.
- Column-mass depth metric: leak-sensitive when the void churn runs (legacy
  VOID); cross-check depth claims with the geometric axis floor from snapshots
  (contiguous-dense-from-bottom; scratch scripts existed in the scratchpad —
  trivial to rewrite) and tracer frames (CRATER_TRACE=1 CRATER_SNAP_DT=...,
  plot_crater_frames.py).
- Detach long runs with nohup + caffeinate FROM hydrocodes/euler/ (a wrong-cwd
  launch silently no-ops — burned once, 2026-07-11). Verify with ps.
- state/ is gitignored; cite key numbers in ledgers, not paths.
- GPU port of all Option C + damage machinery: AFTER the physics closes.

## 5. Definition of done

1. New gates + all 44 existing PASS; feature-off bit-identical.
2. a=3000: sub-crater D pervasive; floor REBOUNDS from the ~16 km transient;
   endpoint TDEC-insensitive; no drain/stew (stability that Option C bought
   must survive).
3. A parameter point (or narrow region) puts the complex branch in the Garvin
   fresh band (ratio ~1) at a >= 1000 with a=300 staying ~1.5x.
4. Ledgers updated (this file + PICKUP-optionC.md pointer + PICKUP-phase3.md
   banner + memory), then the GPU port as a separate follow-on.
