# PICKUP — Option C: implicit Bingham viscosity for acoustic-fluidization collapse

> ## STATUS (2026-07-09 am; commits e54f12f..HEAD)
> Sections 3.1-3.4 DONE (implicit solver + 4 gates green, 43/43; column-mass metric).
> Section 3.5 calibration BLOCKED by an upstream AF-model finding; diagnostic arc
> (R1-R4, state/optc/) run in its place -- full story + decision options in the
> session consult (2026-07-09) and the banner in PICKUP-phase3.md. Short form:
> 1. STAGE-1 (eta 2.4e8/1e9/3e9, Y_B 3e6): a=3000 sub-crater column DRAINS (axis
>    floor 13->47 km), eta-insensitive; NOT the void-reset mass sink (closing it
>    conservatively -- VOID_CONS, arg17 -- left the drain; naive mass-only transfer
>    also adds a downward-pump blowout at a=1000: keep VOID_CONS OFF for now).
> 2. ROOT CAUSE: vib seeded from 0.5*|flow speed| near the transient (~2 km/s)
>    -> pvib ~ 8 GPa >> overburden -> af 0.96-0.98 down a 13-20 km column for
>    several TDECs (snapshot evidence). R1 (c_vib=0.1) does NOT fix it.
> 3. R2 (env AF_SEEDDP: vib = C_ACT*dP/(rho*cs), acoustic): drain GONE (floor
>    pinned 12.3 km for 4000 s ~ 3 TDECs) but ZERO slump -- (1-af) cut at af~0.5
>    leaves Y ~ 0.5*friction ~ driving stress. Over- vs under-fluidized seesaw.
> 4. R3/R4 (env AF_PREL: W&I/Melosh overburden relief, friction at P_eff=P-pvib,
>    replaces (1-af)): BEST-BEHAVED variant -- a=1000 floor clean+flat 4.50 km all
>    run; a=3000 pinned 12.1 km to t~4800 (3.4 TDECs; slow-stew jump delayed, not
>    gone). But STILL ZERO FLOOR UPLIFT at both anchors.
> 5. WHY (R4 snapshot, t=400): relief active only in a 2-3 km skin under the floor
>    (pvib>P); at 8-14 km depth (the region that must flow for rebound) acoustic
>    vib is 2-9 m/s = pvib 8-42 MPa vs P 50-160 MPa. One-shot seeding is either
>    3 orders too strong (flow-speed) or 1 order too weak at depth (acoustic).
>    The literature block model is a DYNAMIC vibration-energy balance: decay +
>    REGENERATION from the deformation (self-sustained while collapse flows,
>    self-arresting when it stops). Our P_ACT gate deliberately blocks
>    re-excitation -- that amputated the essential feedback.
> DECIDED 2026-07-10 (Greg): (a)+(b) adopted. REGENERATION IMPLEMENTED (AF_REG;
> crater arg18=C_REG, arg19=AF_SAT; commit cc683bb): C_REG fraction of viscous +
> plastic dissipation -> E_vib; decayed E_vib -> heat (fixes the pre-existing
> silent deletion); pvib saturates at AF_SAT*pov. vib_regen gate: sustained 22x
> control while shearing, re-locks after arrest, ledger exact; 44/44 PASS.
> Crater defaults now AF_SEEDDP+AF_PREL (env AF_LEGACY_* for A/B).
>
> ## SWEEP-REDO LADDER (2026-07-10/11, state/optc/ *_cr*/_ca*)
> - C_REG=0.1, C_ACT=0.5: FIRST-EVER structurally clean endpoints at BOTH anchors
>   (a=3000 axis floor 14.90 km flat t=400..5600, af band pvib=P SUSTAINED by
>   regen; a=1000 5.77 km flat). But ZERO floor uplift: ratio ~3.5-4x Garvin =
>   friction-only depth, stable. Cause (radial af profile): the fluidized lens is
>   ARCH-protected by an off-axis ring at full lithostatic P where pvib<<P.
> - C_ACT ladder 1.0/1.5 (SAT=2, C_REG=0.3): arch does NOT yield -- floor pinned
>   16.3 km whole run (deeper transient, more widening). AND the frames show the
>   third failure mode: REGIONAL MELTDOWN -- at C_ACT=1.0 the whole near-field
>   crust (r<70 km, z<10 km) fluidizes, spreads, drains into the sponge; domain
>   surface sinks 8-10 km. Amplitude can NEVER work: relieving the deep ring
>   relieves the shallow crust everywhere more (P grows with z; pvib decays with
>   spherical r). The (seed-amplitude x relief) plane is EXHAUSTED.
> - FINAL DISCRIMINATORS (2026-07-11, both NEGATIVE):
>   (i) RESOLUTION: cppr8 vs cppr5 at a=1000, 100-s snapshots through the dynamic
>       phase: axis floor 6.16 vs 6.10 km, both DEAD FLAT from t=100. Converged;
>       no hidden/starved rebound. Grid exonerated.
>   (ii) SHORT-TDEC (TDEC=200 as regional filter, regen carries the active zone,
>       C_ACT=1.0 SAT=2 C_REG=0.3, a=3000): filter WORKS (no regional meltdown,
>       crater af regen-sustained at 7+ decay times) -- and the floor sits at
>       16.10 km from t=100 to t=2000. No rebound in the dynamic window either.
>
> ## VERDICT (2026-07-11): NEGATIVE RESULT, matrix complete
> The AF model family as implementable in this code -- {(1-af) cut | P_eff
> relief} x {flow-speed | acoustic seeding} x {one-shot | regenerated} x
> amplitude ladder x TDEC regimes x resolution -- does NOT produce complex-crater
> FLOOR REBOUND (stratigraphic uplift). Failure modes span: strengthless pipe /
> under-fluidized lock / arch-lock (stable at transient depth, ratio 3.5-4x
> Garvin) / regional meltdown. Wall slumping, transient craters, simple craters,
> shock physics: all reproduced; stable complex arrest: NEW (regen).
> ROOT CAUSES (vs community codes, iSALE/CMI-2004):
>   1. TENSION-ONLY DAMAGE: grow_damage never damages the compressed sub-crater
>      volume (D~0.04-0.10 in every run) -> the cavity is held open by INTACT
>      Lundborg rock (arch ~2e8 Pa). iSALE targets are pervasively D=1
>      (shear/plastic-strain damage, CMI 2004): cohesionless mu_d rubble BEFORE
>      collapse. Their AF beats damaged friction; ours must beat intact rock.
>   2. STRENGTHLESS VACUUM FACE: the free-surface flux carries no deviatoric
>      stress; incipient uplift columns shed their crown into the void churn
>      (voidrst ~1e9/run; same limitation as the Prater plateau-low).
> RECOMMENDED NEXT (fresh scope, Greg to approve): implement CMI-2004
> shear-strain damage growth (failure strain eps_f(P), compressive branch) --
> well-scoped, directly attacks the arch with standard literature physics; all
> Option C machinery (implicit viscosity, relief, regen, gates) sits ready on
> top. Free-surface strength = the deeper follow-on surgery.

Self-contained handoff for a fresh session. **This is the entry point for the
next major effort on the euler code.** Repo: github.com/oklo/deorbit, code in
`hydrocodes/euler/` (CPU oracle `hydro_cpu.cpp` = FP64 calibration engine;
GPU port follows later per the established CPU-first rhythm). Conventions,
build lines, and gate usage: bottom of PICKUP-phase3.md. All gates green
(39 PASS in `./gates.sh quick`) at the time of writing.

## 1. The problem (one paragraph)

The code reproduces shock physics at benchmark level (Pierazzo-2008 envelope,
vertical + 45-deg oblique, cppr12/16) and simple craters within ~1.5x of the
Garvin (2003) Mars d/D level. It does NOT reproduce gravity-driven
complex-crater collapse: complex craters finish 2.5-5x too deep. A week of
laddering and instrumented diagnosis (2026-07-03..07) established that this
is a RHEOLOGY gap, not a parameter-tuning problem: our fluidized material has
only two states — essentially strengthless (af-scaled strength cut, weak
Newtonian viscosity) or instantly locked (pressure-dependent friction once
vib decays). The literature implementations (Wuennemann & Ivanov 2003;
iSALE block model) live in the missing middle regime: fluidized material
flows as a BINGHAM fluid whose effective viscosity (order 1e9-1e10 Pa s at
these crater scales) regulates the collapse rate, and the collapse
self-terminates when driving stresses fall below the effective yield — NOT
when the vibrational energy happens to decay. Our explicit solver caps the
usable viscosity near 2.4e8 Pa s (viscous dt limit, found empirically
2026-06-28), a factor 10-30 below what the physics needs. Hence Option C:
make the viscous update implicit so eta_eff is unbounded, and give the
fluidized state a proper Bingham constitutive law.

## 2. Evidence trail (all committed; do not re-derive)

- Clean vamb=0.27 baseline: friction-only complex branch 1.9->4.0x too deep,
  growing with size (PICKUP-phase3.md banner 2026-07-03; state/dD_settled.csv).
- Flat AF ladder (Tfrac 60/100/200 x Efrac .003/.01/.03): depth==transient
  everywhere; Efrac irrelevant; all rungs WORSE than friction-only. TDEC
  (impactor-acoustic scaling, 6-200 s) decays before the ~700 s slump time.
- Long-TDEC test (1000/2000 s): no slump either — while fluidized ~strengthless
  the cavity deepens without limit; on decay, friction (Y_d ~ 0.6 rho g z)
  locks any shape. The two-state pathology in its pure form.
- Bingham floor (Y_BINGHAM, crater arg15, committed): FIRST partial slump —
  wall overturn + ~70% widening at a=1000/2000, d/D down 15-25% toward Garvin
  — but Y_B in 1-10 MPa is indistinguishable: arrest is still by vib decay,
  not stress balance. state/bingham/dD.csv.
- Deep diagnostic (counters + field snapshots, state/diag/): reseed hypothesis
  refuted (seed counter frozen at 37,118 for 4,400 s); af>0.9 skin decays to
  ~6 cells; the a=3000 "40 km breakthrough" is the surf() rho>1350 threshold
  crossing a SMOOTHLY dilating sub-floor column (2662->1021 kg/m3 at constant
  positive ~1.4e8 Pa). The floor material itself never collapses (tracer
  frames figures/frames_a3000.png show the sandy/melt floor intact at -15 km
  to the end while the column below drains). ALSO: ~2.6% of domain mass is
  deleted over such a run by RHO_CFL<100 -> VOID_AMB resets (the churn);
  a=1000 control is clean. All of this is DOWNSTREAM of collapse not
  finishing (hours-long stewing of an unloaded ultra-deep transient).

## 3. The fix — design (high-quality, not the cheap variant)

Goal: fluidized cells obey a Bingham law tau = Y_eff + eta_eff * gamma_dot
with eta_eff reachable up to >=1e10 Pa s, unconditionally stable; collapse
endpoint set by stress balance (TDEC-insensitive once decay outlasts collapse).

3.1 IMPLICIT VISCOUS UPDATE (the core).
  - Operator-split, once per step, replacing the explicit af-viscosity term
    currently inside Lop (hydro_cpu.cpp ~line 195: `if(g.af[c]>0&&ETA_AF>0)`
    Cartesian vector Laplacian). REMOVE that explicit term when the implicit
    path is active — no double counting.
  - Solve per velocity component: (I - dt * div(nu grad)) v^{n+1} = v^n with
    nu(c) = af(c)*ETA_AF/rho(c), variable coefficient, FACE-AVERAGED
    (harmonic or arithmetic mean of cell nu) for a symmetric operator.
    Momentum form recommended: solve on velocity, then update momenta
    mu = rho*v (the grid stores momenta; convert in/out).
  - AXISYMMETRY: the r-component Laplacian carries the cylindrical terms
    (d/dr(1/r d(r v_r)/dr) includes -v_r/r^2); z-component is standard.
    Enforce regularity at r=0 (v_r=0 on axis; Neumann for v_z).
  - Solver: the grid is small (~155x115). Red-black SOR or conjugate gradient
    both fine; CG on the SPD operator with Jacobi preconditioner is the
    robust choice. Tolerance ~1e-8 relative; typical iterations O(50-200) at
    the target nu*dt/dx^2 ~ 10-100. Cost is negligible vs the hydro step.
  - Boundary cells / void: nu=0 outside fluidized dense material (rho>RHO_CFL
    and af>1e-3) — operator reduces to identity there; keep it in the system
    (uniform treatment, no masking bugs).
  - ENERGY: viscous dissipation must heat: dE = eta_eff * gamma_dot^2 * dt
    (or compute as the KE lost by the implicit smoothing, added to E in the
    same cells — exact bookkeeping, conserves total energy). The current
    explicit term does NOT heat; fix this in the new path and note it.

3.2 BINGHAM CONSTITUTIVE (keep + document).
  - Retain the Y_BINGHAM floor in vonmises (committed, arg15): while
    fluidized, Y_eff = max((1-af)*Y_base, min(Y_B, Y_base)). Together with
    the implicit eta_eff this IS a (regularized) Bingham material: yield from
    the radial return, plastic viscosity from the implicit solve. This is
    the same decomposition iSALE uses. No Papanastasiou regularization
    needed at these resolutions.

3.3 PARAMETERIZATION + TIMESCALES.
  - eta_eff: sweep DIMENSIONAL eta at the anchors first (2.4e8 [old ceiling,
    sanity overlap], 1e9, 3e9, 1e10 Pa s); once the working range is known,
    re-express as Efrac with a CRATER-scale length (eta = Efrac*rho*c_s*L,
    L ~ transient radius) so one setting spans sizes. Do not re-adopt the
    impactor-scale L=a — that scaling is part of why the ladder was flat.
  - TDEC: keyed to the crater gravitational timescale (2-3x 2pi*sqrt(D/g),
    500-1400 s at the anchors — already the established practice). The
    ACCEPTANCE CRITERION is that the final d/D becomes INSENSITIVE to TDEC
    (say TDEC and 2x TDEC agree within run-to-run scatter): that is the
    signature that stress balance, not decay, terminates collapse.
  - Y_B stays a calibration lever (1-10 MPa) — with real viscosity it should
    finally MATTER; if it still does not, that is diagnostic information.

3.4 GATES (write them BEFORE the calibration sweep; extend gates.sh).
  a. Feature-off bit-identity: eta path disabled -> all 39 existing gates
     PASS unchanged.
  b. Diffusion analytic: a sinusoidal velocity mode in a quiescent fluidized
     slab must decay as exp(-nu k^2 t) with dt >> the explicit stability
     limit (this is THE test that the implicit solve is right; check both
     components, axisym r-term included via a mode at large r).
  c. Explicit-implicit agreement at small eta (both stable): same run, same
     answer to FP tolerance over a short window.
  d. Bingham arrest: a fluidized slope with driving stress < Y_B must be
     static; > Y_B must flow and arrest as it shallows. (Analytic: infinite
     slope stability at tau = rho g h sin(theta).)
  e. Energy: total (KE+IE+PE) drift over a viscous-dominated window < 1e-4.
- Also worth doing in this effort (small, improves calibration reads):
  make the depth/surface diagnostic threshold-robust (column-mass surface
  instead of first rho>1350 cell) so metric jumps cannot recur.

3.5 CALIBRATION SWEEP (after gates).
  - Anchors a = 300/1000/2000/3000, Y_d0 = 5 MPa, vamb = 0.27,
    eta ladder x Y_B {1,3,10 MPa}, TDEC = 2x crater timescale + the TDEC
    doubling check at one anchor. Harvest: analyze_crater.py ->
    mars_dD_oracle.py (the ratio-vs-D table; target ratio -> 1 across the
    complex branch WITHOUT wrecking the a=300 simple guard, currently 1.5).
  - Use the tracer frames for physical inspection of the collapse mode:
    CRATER_TRACE=1 CRATER_SNAP_DT=100 CRATER_SNAP_DIR=... then
    `uv run python hydrocodes/euler/plot_crater_frames.py <dir> --times ...`
    (Collins et al. 2020 Fig. 2 style; see figures/frames_a3000.png for the
    current no-collapse pathology rendered this way).

## 4. Loose ends from the current session (check before starting)

- The a=1000 traced Bingham run may still be finishing (state/trace/a1000/,
  log state/trace/a1000.log). When done: render
  `plot_crater_frames.py ... state/trace/a1000 --times 100 800 1600 3400
   --out figures/frames_a1000` and commit both frames figures.
- state/ is gitignored; the evidence CSVs/logs cited above are local.
- figures/checks/ is gitignored scratch space for quick-look renders.
- Verify processes with `ps aux | grep hydro` before assuming anything runs.

## 5. Definition of done

1. All new gates + the 39 existing PASS; feature off = bit-identical.
2. At eta >= 1e9 Pa s the a=3000 crater collapses from its ~12-16 km
   transient to a settled floor, TDEC-insensitively, with no deep-column
   stewing pathology (it should vanish along with the long-lived cavity).
3. A (eta_eff, Y_B) pair (or narrow region) puts the complex branch within
   the Garvin fresh band (ratio ~1) at a >= 1000 while a=300 stays ~1.5x.
4. Ledger updated (this file + PICKUP-phase3.md banner + memory), then the
   GPU port of the implicit path as a separate follow-on.
