<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 52 — energy-balance transport (proposed 2026-10-29; on `physics/energy-balance`)

**Owner request:** "Start B" (2026-10-29): the non-local transport model that ARCHITECTURE.md 11.3 lists as Planned
(the nanoscale gap found by the FinFET deck; excluded from Unit 23 by T14). The branch starts from `main` (`cb3f9de`).
**Owner rulings (2026-10-29):** (1) energy balance: carrier temperatures with an energy-flux equation, no convective
(momentum) term; not the full hydrodynamic model; (2) both carriers, T_n and T_p, always together when the model is on;
(3) a carrier-temperature mobility that reproduces Canali's v(E) in a uniform field, and impact ionization driven by
the carrier temperature; (4) supported with Boltzmann and Fermi–Dirac statistics, SRH, Auger and radiative
recombination, Masetti and Lombardi mobility, MLDA, transient and small-signal analyses; refused with electrothermal,
thermionic interfaces, gate tunnelling, band-to-band tunnelling, the weighted terminal current and Schottky contacts.
The standing conditions hold: no change to any device without the model switched on, criteria fixed before code,
Release and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without authorization.

**Sources read (before coding).**
- T. Grasser, dissertation, TU Wien (open access, iue.tuwien.ac.at/phd/grasser), §3.2.1–3.2.2: the carrier-temperature
  mobility after Hänsch, μ = μ_LIS / (1 + α ΔT), α = 3 k μ_LIS / (2 q τ_ε v_sat²), and the energy relaxation time
  as a constant or a Monte Carlo fit (the fit's constants are not on the page).
- A. Triebl, dissertation, TU Wien (iue.tuwien.ac.at/phd/triebl, node17): constant energy relaxation times, typically
  0.35 ps for electrons and 0.4 ps for holes in silicon, sometimes used as a fitting parameter.
- A. Aste and R. Vahldieck, "Time-domain simulation of the full hydrodynamic model", arXiv physics/0312021 (open
  access): the energy-balance model as most commercial simulators solve it (current
  j = −(e/m) τ_p (n E + ∇(n k T / e)), energy ∇·((5/2) n k T v) − ∇·(κ ∇T) = −e n v E − (3/2) n k (T − T_L)/τ_w,
  κ = (5/2 + r) n k² μ T / e), and the 0.3 µm silicon ballistic diode benchmark (n⁺ 1e18 / n 2e15 / n⁺ 1e18, 0.1 µm
  each, 1 V, m = 0.26 m₀, v_s = 1.03e7 cm/s, Caughey–Thomas μ₀(N) with μ_min 80, μ_max 1430, N_ref 1.12e17, 0.72,
  β = 2): peak velocity 1.87 v_s with energy balance, 2.61 v_s with the full model. Its τ_w is
  temperature-dependent (Baccarani–Wordeman) and r = −2.1, unlike this unit, so its number is reported, not gated.
- Legacy: `pytcad/hydrodynamic.py` (M29): a local closure after the solve (no energy flux, no Newton unknown, no
  feedback into the current). OLD / NEW / REASON:
  - OLD: T_n from q μ E² = (3/2) k (T_n − T_L)/τ_w evaluated at the solved field, post-processing only.
  - NEW: T_n and T_p as Newton unknowns with the energy flux, self-consistent, exact Jacobian.
  - REASON: ruling 1; without the flux term the closure cannot produce overshoot, the effect the unit exists for.
  Its round trip (temperature to effective field to the field-driven ionization) is kept as ED6's construction.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| ED1 | Unknowns: the scaled temperature rises τ_n = (T_n − T_L)/T_L and τ_p, interleaved after ψ, n, p (stride 5, `DriftDiffusion::dofs()` names them); insulator nodes hold them at 0 like n and p. Ohmic contacts: τ = 0 (carriers at the lattice temperature). Every other boundary, semiconductor–insulator interfaces included: no energy flux. Heterointerfaces without thermionic emission: T continuous (one unknown per node). | τ rather than θ for resolution, as Unit 23's U1. | Stride 5 doubles the edge blocks' size (25 vs 9 entries). |
| ED2 | Particle current (Stratton, energy-independent momentum relaxation): J_n = μ_n(T_n) (n ∇E_c + k ∇(n T_n)), holes mirrored. Discretised as Unit 23's flux with r = 0: Scharfetter–Gummel in w = θ_n n, the edge's mean θ_n,e in δ and in the coefficient. At θ_n = 1 it is the isothermal flux. | The form Aste–Vahldieck call the commercial one; reuses Unit 23's tested scheme. | The edge-mean temperature is second order, not Tang's exponential fit. |
| ED3 | Energy: ∂((3/2) n k T_n)/∂t + ∇·S_n = J_n·E − (3/2) n k (T_n − T_L)/τ_w,n − (3/2) k T_n U_net, S_n = −κ_n ∇T_n − c_n (k T_n / q) J_n, κ_n = c (k²/q) n μ_n T_n, c = 5/2 (r = 0, EQ3; the same energy-independent momentum relaxation as ED2; not the material's thermopower exponent of Unit 23). S_n is one Scharfetter–Gummel flux in n θ_n² with the same δ and the coefficient c V_T. The Joule term on each edge is the edge current times its E_c difference, half to each end's box, so the summed sources equal Σ V_c I_c exactly. | κ by Wiedemann–Franz as all three sources; the edge split makes discrete energy conservation exact (EC7). | r is one number per material, not fitted to Monte Carlo. |
| ED4 | Mobility (ruling 3): an effective field E_eff(T) solves the uniform-field balance q τ_w μ_C(E) E² = (3/2) k (T − T_L), μ_C Canali's law of the edge's low-field mobility (Caughey–Thomas or Masetti, Lombardi by Matthiessen's rule); μ(T) = μ_C(E_eff(T)) at the edge's mean T. β = 1: closed form; β > 1: Newton to round-off; the partial through the implicit-function rule. The Canali parameters therefore enter through T, and `field_mobility = on` with energy balance is refused (counted twice). For β = 1 this is Hänsch's form up to the factor in its denominator. | Exact Canali v(E) in a uniform field for every β (EC3); silicon electrons have β = 2, which Hänsch's form does not reproduce. | One scalar Newton per edge and carrier per evaluation. |
| ED5 | Energy relaxation times: constant per material and carrier, `energy_relaxation_time_n/p` [s]; silicon 0.35 ps (electrons) and 0.4 ps (holes), Triebl. Another semiconductor without them is refused, none invented. Deck overridable per region. | Cited typical values; ruling 3 needs τ_w in ED4 and ED6 too. | No T-dependent τ_w (Aste–Vahldieck's): a follow-up. |
| ED6 | Impact ionization (ruling 3): the existing coefficients α(E) evaluated at E_eff(T_n) (and T_p) of the node, the generation α_n |J_n| + α_p |J_p| as Unit 19, the lattice temperature factor unchanged. In a uniform field E_eff = E, so it equals the field-driven model (EC3). No extra E_g energy sink per ionization event. | The legacy's round trip, now self-consistent; non-local by construction (a carrier ionizes where it is hot, not where the field is). | Without the E_g sink T is overestimated near breakdown; recorded as a limitation. |
| ED7 | Recombination energy: −(3/2) k T_n U_net in the electron energy row and −(3/2) k T_p U_net in the hole's (SRH, Auger, radiative, II as negative U). | The carriers removed take their mean energy. | — |
| ED8 | Fermi–Dirac (ruling 4): pressure n k T F_{3/2}(η)/F_{1/2}(η) in the current (Unit 14's ν-factor generalised: δ gains Δ ln of the ratio), stored energy (3/2) times the pressure, the energy-flux factor (5/2) F_{3/2}/F_{1/2} (r = 0, EQ3: the same ratio as the pressure; normalised integrals). New: F_{3/2} with double-double references (C++). | Every Boltzmann formula is the η → −∞ limit (EC2). | — |
| ED9 | Transient and small-signal: the storage terms of ED3 in the mass matrix; the AC system takes the term's partials. | Ruling 4. | — |
| ED10 | Switch: `[models] energy_balance = off | on` (both carriers, ruling 2). Off: no term is created, the stride stays 3; residuals, Jacobian patterns and values, run identities, device files and the Unit 40 traces are those of `main`. | GC1's pattern. | — |
| ED11 | Refused (invalid_input) with: electrothermal, thermionic interfaces, gate tunnelling, local or nonlocal BTBT, the weighted terminal current, Schottky contacts (ruling 4); `field_mobility = on` (ED4, EQ2); a semiconductor without τ_w (ED5); Lombardi with impact ionization stays refused as on `main`. | Ruling 4. | Follow-ups in ARCHITECTURE 11.3. |
| ED12 | Start and continuation: a sweep's first point starts from the drift-diffusion (field_mobility off) solution with τ = 0; later points as every sweep; bias step cutting (Unit 28) and the line search (Unit 30) apply unchanged. | Equilibrium is exactly τ = 0 (EC5). | — |
| ED13 | Results: T_n and T_p as nodal fields in `BiasPoint`, transient and small-signal results, the result files, `nitcad plot` and the application's field list; the run record names the model and each region's τ_w. | Needed to see overshoot (EC10). | Result files gain two fields only when on. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| EC1 | Model off (and every existing deck): every existing test and pin passes unchanged; residuals, Jacobian patterns and values, run identities, device files and the Unit 40 traces byte-identical to `main`. |
| EC2 | Physics functions: E_eff(T) inverts the uniform-field balance to 1e-13 relative for β = 1 and 2 over 1e2–1e6 V/cm and τ_w 0.1–1 ps; β = 1 closed form against the Newton solve within 1e-14; F_{3/2} against double-double references within 1e-14 relative over η ∈ [−40, 60]; every partial against a central difference within 1e-7 relative. |
| EC3 | Uniform-field limit: 1D silicon resistors (1e16, 1e17, 1e18 cm⁻³, both carriers) at 1e3–1e5 V/cm, at least 50 energy relaxation lengths v_sat τ_w long: in the middle third the velocity J/(q n) equals Canali's v(E) within 1e-3, T equals T_L + (2/3) q τ_w v E / k within 1e-3, and with impact ionization on the generation equals the field-driven model's within 1e-3. |
| EC4 | Drift-diffusion limit: τ_w = 1e-18 s: the diode's and the 2D MOSFET's terminal currents equal drift-diffusion's (field mobility off) within 1e-6 relative over their sweeps. |
| EC5 | Equilibrium: at zero bias |τ_n|, |τ_p| ≤ 1e-12 at every node, and the solution's ψ, n, p equal drift-diffusion's within 1e-12 relative. |
| EC6 | Jacobian: assembled against central differences on 1D, 2D and 3D devices; Boltzmann and Fermi–Dirac; with impact ionization, Lombardi, MLDA, Masetti; the transient mass: every entry within 1e-6 relative (or 1e-12 of the row scale). |
| EC7 | Energy conservation, steady state: Σ_c V_c I_c equals the summed relaxation power plus the recombination energy terms within 1e-9 relative (1D, 2D, 3D; Boltzmann and Fermi–Dirac). |
| EC8 | Current conservation: terminal currents sum to zero within 1e-10 relative of the largest. |
| EC9 | Refinement: the ballistic diode of EC10 on h, h/2, h/4: the peak velocity and peak T_n converge with an observed order ≥ 1.5. |
| EC10 | Overshoot (Aste–Vahldieck's diode, their mobility and v_s, this unit's constant τ_w): at 1 V the peak electron velocity exceeds 1.3 v_s and lies in the n region within 30 nm of the drain junction; drift-diffusion with Canali on the same mesh never exceeds 1.0001 v_s. Their 1.87 v_s is reported beside ours, not gated (different τ_w and r). |
| EC11 | Channel length: a 2D MOSFET at L = 2 µm: energy balance and drift-diffusion with Canali agree within 5% in saturation; at L = 0.1 µm the energy-balance current is larger. Both reported. |
| EC12 | FinFET deck: with energy balance the Id–Vg sweep converges at all 17 points; the on-current change and the wall time against drift-diffusion are reported (cost, no gate on the effect's size). |
| EC13 | Transient: a bias step on the ballistic diode ends within 1e-8 relative of the steady state; over the step, input energy equals relaxation plus recombination energy plus the change in stored energy within 1e-6. |
| EC14 | Small-signal: at 0 Hz the admittance equals the DC sweep's dI/dV (central difference) within 1e-5 relative. |
| EC15 | Every ED11 refusal tested with its message. |
| EC16 | Deck, device file, results and run record: T_n, T_p and τ_w round-trip byte-identical; `nitcad plot` and the application show T_n and T_p. |
| EC17 | Hand mutants: at least one per design row (ED1–ED13), every survivor closed by a new test. |
| EC18 | Release and Debug suites pass at every commit, `/W4 /WX`. |
| EC19 | Cost reported: assembly and solve time per Newton iteration for 1D, 2D, 3D against drift-diffusion. |

**Owner answers (2026-10-29), design and criteria approved:** EQ1 0.35 ps and 0.4 ps (Triebl); EQ2 refused; EQ3
c = 5/2 (r = 0), not the recommended material r (ED3 and ED8 changed to it before code); EQ4 the five slices.

**Found before code, ruled by the owner (2026-10-29): β = 1.** Exact Canali reproduction makes x = μ₀E_eff/v_sat
grow like √a near T_L (a ∝ T − T_L), so for β < 2 dμ/dT is infinite at T = T_L and μ is undefined below it (β = 2 is
analytic through T_L: y = x² = a (a + √(a² + 4))/2). Ruling: regularize. For β ≠ 2, x is exact for a ≥ a₀ = 1e-8
(about 5 µK of hole heating in silicon), the cubic c₂a² + c₃a³ matching the exact value and slope at a₀ below it, and
0 for a ≤ 0, so μ is C¹ and equals μ₀ for cooled carriers; β = 2's mobility takes the analytic form everywhere (cooled
carriers faster, as Hänsch). The effective field takes the regularized x for every β (it feeds only impact
ionization, zero to all orders there). Measured: the regularization moves silicon holes' mobility by at most
3.7e-5 μ₀, only below a₀.

**ED5 as built:** `SemiconductorParameters::energy_relaxation` defaults to 0 (not given) and silicon's set does not
carry it, so device files and run identities without the model are unchanged (EC1); silicon's cited values
(`physics::silicon_energy_relaxation`) apply to regions with silicon's mobility that give none, as Lombardi's and
MLDA's silicon checks.

**Criteria defects found building slice 2, ruled by the owner (2026-10-29).** Measured first; each criterion as written
could not be met for a reason in the criterion, not the model. The rulings replace the rows above where they differ.
- **EC3 at 1e5 V/cm.** The middle third of a 3 µm resistor is not uniform at 1e5 V/cm: the carriers enter cold at the
  contact, and in velocity saturation (dv/dE ≈ 0.6 cm²/Vs) the space charge that adjusts the field is not screened
  (1e17 electrons: E from 2.3e4 to 1.6e5 V/cm across 3 µm; a 30 µm resistor approaches 1.1e5 V/cm but still varies
  6% across its middle third). Misses at 1e5 V/cm: up to 1.4% velocity, 4% temperature (1e16 electrons); 1e3 and 1e4
  V/cm pass (worst 2.4e-4). Drift-diffusion with Canali is exactly uniform because its contact matches the bulk.
  Ruling: at 1e5 V/cm a discrete bulk check: the exact uniform state (neutral densities, E uniform, each carrier's
  rise at the local balance) leaves every interior row zero within 1e-12 of its terms (Σ_c |J_rc x_c|); measured
  2.8e-16. The solved resistors stay at 1e3 and 1e4 V/cm, at 1e-3.
- **EC4.** τ_w → 0 alone tends to drift-diffusion with Canali's μ(E), not constant mobility; against drift-diffusion
  with Canali it still differs by 1–2% in a forward diode (energy balance evaluates Canali at the power-balance field,
  drift-diffusion at its driving field). Ruling: EC4 is the double limit τ_w = 1e-18 s with v_sat = 1e30 cm/s against
  constant-mobility drift-diffusion within 1e-6 (measured 1.3e-8 on the diode); the τ_w-only comparison is printed.
  The MOSFET half (L = 0.6 µm) at 1e-18 misses: 6.8e-7 at Vd = 0.1 V rising to 1.23e-6 at 1.5 V. That is the
  first-order τ_w heating, since with v_sat = 1e30 the channel's drift velocity is unclamped (T_n − T_L about 4e-4 T_L).
  At τ_w = 1e-20 it is exactly 100× smaller, 1.23e-8 (49 Newton iterations, no substeps). Ruled 2026-10-29: the
  MOSFET half uses τ_w = 1e-20 at 1e-6; the diode stays at 1e-18.
  **Unresolved currents (ruled 2026-10-29, "Skip below resolution"):** on the gate ramp the source and drain
  currents lie below drift-diffusion's own `terminal_current_resolution`, so a relative comparison there compares
  rounding. EC4 compares only currents above the drift-diffusion reference's resolution (EC8's resolution) and
  prints how many it skipped. Measured: 12 of 42 currents skipped, worst 1.23e-8 over the other 30.
- **EC8.** The terminal currents' sum is the Newton residual's: on EC4's diode at −0.4 V, 2.1e-10 A/cm² of 4.4e-4
  with energy balance and 1.2e-10 with plain drift-diffusion, against a resolution of 2.3e-7. Ruling:
  |Σ I| ≤ max(1e-10 · largest, terminal_current_resolution).
- **EC10.** Drift-diffusion's J/(q n) reaches 1.34 v_s at the drain junction (diffusion adds to drift), so "never
  exceeds 1.0001 v_s" is wrong. Ruling: drift-diffusion's drift velocity μ_C(E)·E ≤ v_s on every edge, and energy
  balance's peak exceeds drift-diffusion's J/(q n) peak (1.93 against 1.34 v_s).
- **EC7.** Σ V_c I_c is not relaxation plus recombination: the recombined pairs carry their gap (E_c − E_v), and the
  ohmic contacts carry the Peltier term (E_c + V) D_n + (E_v + V) D_p, the energy flux into the contact nodes (τ = 0
  there) and the contact edges' half Joule heat. Ruling: `DriftDiffusion::energy_powers` returns each term and the
  Newton residual's power; the identity with the residual closes within 1e-9 of the largest term (measured ≤ 8e-16,
  1D/2D/3D at −0.3, −0.6, −0.8 V), and |residual| ≤ max(1e-9 · largest term, the current resolution times the largest
  |E_c|, |E_v|) (measured at least 500× inside). Without the residual the miss is ≤ 1.9e-10 at −0.6 and −0.8 V and
  3.4e-6 at −0.3 V. Measured too: the relaxation power is negative in the forward diode (carriers climbing the
  junction cool below T_L and draw heat from the lattice; the Peltier term carries it out).
- **Newton on the EC11 MOSFET (found 2026-10-29).** The energy-balance sweep of the 2 µm MOSFET did not finish
  (one point 210 iterations, the first drain step past 5,000 at residual 1e-11). Traced from the equilibrium start to
  Vg = 0.25 V: (1) the first iterate after a bias step carries spurious currents whose Joule heat throws T_n to about
  1600 K in the drain, and the coupled iteration then diverges (corrections to 1e19); (2) a carrier far below its
  node's majority (holes at 1.6e-5 cm⁻³ against 3e17 electrons) has an energy row that scales with its density,
  so its temperature is numerically undetermined (τ of 30–70 seen) and the update test never passes. At the
  converged first point the Newton step is 5e-15: the system is regular at the solution. Rulings:
  - **ED12 as ruled:** every bias point starts in two stages: ψ, n, p converged with τ held at the point's start
    (`DriftDiffusion::hold_temperatures`: every energy row τ − τ_held), then the coupled system. At the first point
    τ_held = 0, so its first stage is drift-diffusion, ED12's original rule. The first stage's iterations are not
    recorded or reported.
  - **ED14 (new): the relaxation floor.** Each carrier relaxes with density max(c, 1e-10 · max(n, p)) in the
    relaxation term, so a carrier below 1e-10 of its node's majority is held at T_L; elsewhere the physics is
    unchanged. `energy_powers` counts it in the relaxation power.
  - **ED14 amended (ruled 2026-10-29, "Floor on |doping| too"):** the floor is 1e-10 · max(n, p, |N_D − N_A|).
    Found on the 2 µm MOSFET at Vd = 0.6 V: a depleted node under the gate near the drain (y 22.5 nm,
    n 3e12, p 2 cm⁻³) had a floor of 300 cm⁻³, below the hole energy fluxes from its neighbours (whose holes rise
    steeply into the substrate). Newton asked for θ_p < 0 there; the update limit halved θ_p every iteration and
    the update doubled, with the residual at its 1e-10 rounding floor. Referenced to the doping, the floor in
    depletion is about 1e7 cm⁻³ and the node relaxes to T_L; neutral regions, where the majority is about the
    doping, are unchanged. Measured: Vd = 0.6 V in 7 coupled iterations (it diverged), and the whole 21-point
    sweep in 80 s with one substep. The floor's partial in the majority is 0 where the doping sets it.
- **ED11 amended (found building slice 2, not in the row before code; ruled 2026-10-29, "Refuse, add to ED11"):**
  interface charge, traps and surface recombination are refused with energy balance (invalid_input): they act on
  their edges' carrier fluxes, which take no carrier temperature. EC15 tests the refusal with its message; support
  is a follow-up in ARCHITECTURE 11.3.

**Found before slice 3's code, ruled by the owner (2026-10-29).**
- **ED8's flux as ruled ("Derived form").** Read literally (Scharfetter–Gummel in the pressure θ n g, g = F_{3/2}/F_{1/2},
  δ gaining Δ ln g) the edge is not exact at equilibrium: zero flux would need δ/g_e + Δ ln g = δ + Δ ln γ, which an
  edge mean g_e does not give, so EC5 (1e-12) would fail by the discretization error. From the same continuum
  J = μ (n ∇E_c + k ∇(n T g)): with w = θ n, J/(μ k T_L) = ∇w + w (∇e_c/θ − ∇ ln γ + (5/2)(g − 1) ∇ ln θ), so the
  edge is ED2's flux in w with
      d = ±δ/θ_e + Δ ln γ − (5/2)(g_e − 1) Δ ln θ   (+ electrons, − holes in the electron form; Δ = b − a),
  γ Unit 14's degeneracy factor (F_{1/2}(η) e^{−η}), g_e and θ_e the edge's means. The energy flux, with the
  Fermi–Dirac enthalpy (5/2) k T g per particle (ED8's factor (5/2) g in both κ and the convective term), is
  c g_e times the same flux in v = θ² n with the same d (the continuum reduces to it exactly). At θ = 1 it is Unit 14's
  flux exactly; with Boltzmann statistics (γ = g = 1) slice 2's.
- **η of a hot carrier ("Nc at T_n").** n = N_c(T_L) θ^{3/2} F_{1/2}(η): γ and g are evaluated at the density
  n θ^{−3/2} (Unit 14's degeneracy of it), so heating lowers the degeneracy and both carry τ partials; the (5/2) in
  d comes from this choice (N_c at T_L would give 1). Below n θ^{−3/2}/N = 1e-6 g takes its series
  1 + 2^{−5/2} v + (1/8 − 2·3^{−5/2}) v², as Unit 14's γ.
- ED3's relaxation and ED7's recombination energy stay as written ((3/2) k (T − T_L) per carrier and (3/2) k T U), with
  either statistics; EC7's identity is algebraic in them.
- **ED6 as built:** E_eff at a node from the node's own low-field mobility, relaxation time and Canali parameters
  (Lombardi with impact ionization stays refused, so the node's bulk mobility); the generation is Unit 19's
  α(E) |j|²/m per carrier with E_eff(T) in place of the field along the current, and +(3/2) k T G in each energy row
  (ED7's II as negative U).
- **Lombardi** (ruling 4) comes with slice 3, as slice 2's code noted: the edge's low-field mobility takes the
  surface terms in E_perp before μ(T), and the E_perp stencil reaches the energy rows too.

**Criteria defects found building slice 3, ruled by the owner (2026-10-29).**
- **EC3's generation half ("Bulk state, 1e-3").** In the solved resistors α = A e^{−B/E} amplifies E_eff's error by
  B/E (about 120 at 1e4 V/cm): 3.9% at 1e16 electrons, where T is right to 2.4e-4 (inside EC3's 1e-3), 2.6e-4 at
  1e16 holes, the other four ≤ 4.7e-6; at 1e3 V/cm α is 0 in both. α there is about 1e-48 per cm. Ruling: the generation half is gated
  at EC3's 1e5 V/cm discrete bulk state: each interior node's assembled generation equals the field-driven model's
  (Unit 19 with Canali's mobility) at the same ψ, n, p within 1e-3; measured ≤ 2.0e-9 over the six resistors. The
  solved resistors' ratios α(E_eff(T))/α(E) are printed.
- **EC5 at a degenerate n side ("Gate T, ψ, n; print p").** EC5 stays gated in full on its p 1e17 / n 2e17 diode, in
  both statistics (worst 7.5e-15). With the n side at 2e19, T, ψ and n agree within 2e-15, but the minority holes
  (about 1e-21 of the electrons) miss 1e-12 in both statistics (Boltzmann 1.18e-12, Fermi–Dirac 2.64e-12), the
  solve's resolution of them: T, ψ and n are gated there and p printed.

**Slice 3 as built (measured 2026-10-29).** EC6 (FD, impact ionization with either statistics, Lombardi with either,
the held first stage): every column within 1e-6 in 1D, 2D and 3D. At θ = 1 the continuity rows and currents equal
drift-diffusion's within 1e-13 with Fermi–Dirac and with Lombardi. A hot degenerate carrier at constant pressure
across an edge (θ 1.5 to 1.6, η_a = −3, 2, 6) carries flux at most 3.7e-5 of its opposing terms (exactly 0 under
Boltzmann); with N_c at T_L (coefficient 1 for 5/2) it misses by 7.6e-4. EC7 with Fermi–Dirac (n side 2e19, 1D, 2D,
3D) closes within 3e-16 of the largest term with the residual. EC4's diode double limit with Fermi–Dirac (n side 2e19)
against Unit 14's drift-diffusion within 1e-6.

**Slice 4 as built (2026-10-29; owner: "start slice 4: transient and small-signal analyses").**
- **ED9.** Each energy row of a semiconductor node off the contacts stores −V rate (W − c_W), W = (3/2) c θ g the
  carriers' thermal energy (g = F_{3/2}/F_{1/2} at c θ^{−3/2}, ED8 as ruled; 1 with Boltzmann); a step's history holds
  S_n, S_p, W_n, W_p per node (storage width 4). The small-signal matrix takes −V s dW; the small-signal current rows
  take the rise columns and Lombardi's E_perp stencil. `energy_powers` takes an optional time step and returns the
  stored power in two more terms (thermal_storage, carrier_storage; below).
- **Transient** (`solve_transient`): the state at t = 0 in ED12's two stages; the error estimate measures each
  carrier's rise off the contacts (absolute, in units of T_L, as the electrothermal model's); snapshots carry T_n, T_p.
  **Small signal:** the fields carry the carrier temperatures' amplitudes (K/V).
- **ED13 / EC16.** Result files: T_n, T_p under new keys (fields 5, 6; small-signal fields 7, 8), written only when
  present. Device files: the relaxation times under keys 68, 69, written only when given (nonzero), so a file without
  them is the file of before (a test checks the 36-byte difference). Run record: `models.energy_balance` named and each
  region's relaxation times (as given; 0 = silicon's) and the initial carrier temperatures in the identity, only when
  on. Deck: `[models] energy_balance = on | off`; `energy_relaxation.tau_n_s`, `.tau_p_s` per region. `nitcad plot`
  and the application: fields `Tn`, `Tp` (an error naming `energy_balance` when a result has none; the application
  lists them only then).
- **Measured.** EC6 mass: every column within 1e-6 (1D, 1D across a material step, 2D with the oxide and Lombardi, 2D
  with a gate, 3D; Boltzmann and Fermi–Dirac; rates t0/τ_w and 100 t0/τ_w). The small-signal matrix equals the time
  step's Jacobian at s = rate, and at s = iω the steady Jacobian plus iω times the storage's, within 1e-12 of the
  largest entry. The current rows against central differences of the terminal currents within 1e-6 of the row's
  largest entry (beyond the rounding of the current's terms, as EC6's relative steps). EC13 on the ballistic diode
  (0.5 → 1 V at 0.1 ps, backward Euler to 1 ns, 183 steps, 0.4 s): the end state within 9.4e-15 (current), 1.5e-16
  (ψ), 1.0e-15 (n), 2.3e-9 (p), 8.5e-16 (T_n), 1.7e-15 (T_p) of the 1 V steady state; the energy budget below misses
  by 5.5e-14. EC14 at 0 Hz against the DC central difference (±0.1 mV): the ballistic diode at 0.3 V and 1 V within
  6.7e-9 and 7.9e-10; a forward pn diode at 0.6 V within 2.5e-6 with both statistics, which is the central
  difference's own truncation, (ΔV/V_T)²/6 = 2.5e-6. EC15: every ED11 refusal with its message (gate tunnelling,
  nonlocal BTBT, the weighted current and Lombardi with impact ionization added).
- **EC13's energy half: ruled 2026-10-29, accepted as built** (owner, "Accept as built"). Read literally ("input energy equals relaxation plus
  recombination energy plus the change in stored energy") it misses, for EC7's reasons and one more. Built and gated
  at 1e-6: the integral of the conduction power Σ_c V_c I_c equals the integrals of EC7's terms as ruled (relaxation,
  recombination, the energy flux into the contacts, the contact edges' half Joule heat, minus the Peltier term) plus
  the change in stored energy, taken as two terms: the carriers' thermal energy (ΔW, ED9) and the band-edge energy of
  the carriers stored, V rate (E_c ΔS_n − E_v ΔS_p) at each step's band edges (the continuity rows' storage; a step
  that stores carriers moves their potential energy as well as their thermal energy). Measured over the run, in
  units of the input: thermal energy 2.2e-4, band-edge energy 1.5e-4, contact terms 8.9e-2; without the band-edge term
  the miss is 1.5e-4, without the contact terms 8.9e-2. The per-step identity is algebraic (an assembly test closes
  it within 1e-12 at arbitrary states).

**Slice 5 scope, owner ruling (2026-10-29: "Dont take the FinFET case, it will take long").** EC12 (the FinFET deck's
17-point sweep with energy balance, its on-current change and wall time) is not run in Unit 52: not measured, not
claimed. Slice 5 is EC9 (refinement), EC17 (hand mutants), EC19 (cost on 1D, 2D, 3D) and ARCHITECTURE.

**EC9's diode, a criteria defect ruled by the owner (2026-10-29, "Graded 5 nm, gate 1.5").** Measured first on
EC10's diode at 1 V on h = 1 nm, h/2, h/4 and h/8 (nested meshes, the junctions on nodes). The abrupt junctions,
sampled at nodes, make every model first order: energy balance's current converges at order 1.01 (differences 2723,
1349, 671 A/cm²), its peak T_n at 0.79 and then 0.98, its peak velocity at 1.51 and then 1.08 (the 1.51 is
pre-asymptotic); drift-diffusion with Canali's current at 1.0 too (1578, 785, 391). The junction node's doping
averaged cuts the error tenfold but stays first order (the edge mobilities beside it are still wrong). With the
junctions erfc-graded over 5 nm the scheme's own order shows: current 2.00, peak velocity 2.45, peak T_n 2.00
(10 nm: 2.00, 2.83, 1.74). Ruling: EC9 is gated on the diode with its junctions erfc-graded over 5 nm (h, h/2,
h/4, observed order ≥ 1.5 for the current, the peak velocity J/(q n) and the peak T_n over the nodes); the abrupt
diode's orders and drift-diffusion's are printed as the evidence that the first order is the junction's.

**EC19, measured (2026-10-29; Release, P-core pinned 0x3C3; `[.assembly_benchmark]` in
`tests/solve/energy_balance_test.cpp`).** EC7's forward pn diode at −0.6 V; per Newton iteration the assembly
(evaluate, median of 31) and the default backend's (SparseLU) numeric refactorization and solve (best of 3; the
pattern's analysis is once per sweep). Iterations are the sweep's (−0.2, −0.4, −0.6 V) recorded ones: energy
balance's leave out ED12's first stage, as ruled.

| Mesh | Model | Unknowns | Assembly | Refactorization | Solve | Sweep |
|---|---|---|---|---|---|---|
| 1D, 2001 | energy balance | 10,005 | 0.375 ms | 4.17 ms | 0.30 ms | 10 it., 0.16 s |
| | drift-diffusion | 6,003 | 0.094 ms | 1.80 ms | 0.15 ms | 18 it., 0.05 s |
| | drift-diffusion + Canali | 6,003 | 0.219 ms | 1.80 ms | 0.15 ms | 30 it., 0.34 s |
| 2D, 201 × 51 | energy balance | 51,255 | 5.88 ms | 946 ms | 19.8 ms | 10 it., 25.0 s |
| | drift-diffusion | 30,753 | 1.03 ms | 164 ms | 5.2 ms | 18 it., 3.3 s |
| | drift-diffusion + Canali | 30,753 | 1.82 ms | 177 ms | 5.8 ms | 19 it., 273 s |
| 3D, 41 × 11 × 11 | energy balance | 24,805 | 2.83 ms | 2,130 ms | 17.8 ms | 9 it., 56.9 s |
| | drift-diffusion | 14,883 | 0.58 ms | 567 ms | 7.3 ms | 18 it., 10.6 s |
| | drift-diffusion + Canali | 14,883 | 1.11 ms | 564 ms | 6.5 ms | 24 it., 15.8 s |

The factorization dominates in 2D and 3D: energy balance's 5 × 5 blocks make it 5.8× (2D) and 3.8× (3D) drift-
diffusion's; assembly is 4–6× and stays below 1% of an iteration there. The 2D drift-diffusion + Canali sweep's
273 s is its step cutting (Unit 28; the erratic Canali convergence recorded with EC11), not this unit. A first 3D
mesh of 41 × 21 × 21 (92,005 energy-balance unknowns) did not finish its sweep in 35 min on SparseLU and was stopped;
PARDISO or GMRES-ILUT (6.10) is the route for 3D energy balance at that size.

**EC17's hand mutants, as run (2026-10-29).** 19 mutants, at least one per design row (ED1–ED14), each built in
Release and run against the energy-balance tests (physics, F_3/2, assemble, the visible solve cases, deck, io; the
run-record mutants with the whole deck suite), then reverted. 15 were caught at once; 4 survived, and each gap was
closed with a test before the mutant was run again:

| Row | Mutant | Caught by |
|---|---|---|
| ED1 | ohmic contact nodes' temperatures free (relaxation rows) | assemble |
| ED1 | τ_n and τ_p swapped in the unknowns' map | assemble, solve |
| ED2 | d = δ (no 1/θ_e) | assemble (EC6), solve |
| ED3 | c = 2 for 5/2 | survived (EC7's identity is algebraic in c; EC10 gates only > 1.3 v_s); new: a uniform θ carries c θ per particle, both carriers (3.75 within 2e-10) |
| ED3 | Joule heat 0.45 to each end | assemble, solve |
| ED4 | the power balance's 3/2 as 1 | physics, assemble, solve |
| ED4 | the mobility's rise partial doubled | assemble, solve |
| ED5 | electrons given the holes' default relaxation time | assemble, solve |
| ED6 | the ionization energy term without θ_n | assemble |
| ED6 | holes ionizing at the electron temperature | assemble |
| ED7 | no recombination energy in the energy rows | assemble, solve |
| ED8 | (5/2)(g_e − 1) Δ ln θ with 1 for 5/2 | assemble |
| ED9 | W = n θ g for (3/2) n θ g, its partials alike | survived (EC13's budget uses the same storage; the end state does not depend on it); new: the step history's W_n, W_p against (3/2) n θ F_{3/2}/F_{1/2} at hot nodes, both statistics |
| ED10 | `models.energy_balance` named with the model off | deck |
| ED11 | the thermionic-interface refusal removed | assemble |
| ED12 | the two-stage start removed | survived (only the hidden MOSFET sweeps needed it); new: without step cutting, the pn diode to −0.8 V and +5 V in one step and the ballistic diode 0.5 then 1 V (4, 4, 17 iterations; each fails at a singular pivot without the first stage) |
| ED13 | T_n not written to the result file | deck, io |
| ED13 | the identity digesting τ_p for τ_n | survived (the deck test changed both or τ_p); new: changing either time alone changes the identity |
| ED14 | relaxation floor 0 | assemble |

Run again with the new tests: all four caught (ED3a by assemble and solve, ED9 by assemble, ED12 by solve, ED13b by
deck), so 19 of 19; none is recorded as equivalent.

**Slice 5 as built (2026-10-29).** EC9 as ruled (graded 5 nm: orders 2.00, 2.45, 2.00), EC17 (19 of 19), EC19 reported
(above), ARCHITECTURE 6.12.2 and 11.3 as built. EC18: Release (121 tests) and Debug (110) suites pass, zero
warnings under `/W4 /WX`; no existing pin changed. EC12 not run (owner ruling). Unit 52's criteria are then all met
or ruled, for the owner's acceptance.

**Questions for the owner before code (proposal's recommended answer first):**
- EQ1 τ_w for silicon: 0.35 ps / 0.4 ps (Triebl), or another cited pair.
- EQ2 `field_mobility = on` with energy balance: refused (ED4), or silently replaced by the T-driven law.
- EQ3 Energy-flux factor c = 5/2 + r: the material's thermopower r (silicon c = 2, consistent with Unit 23), or c = 5/2
  (r = 0). Aste–Vahldieck note that smaller c suppresses spurious overshoot; c changes EC10's peak.
- EQ4 Build order: slice 1 physics (ED4, ED5, ED8 functions; EC2), slice 2 Boltzmann assembly and steady sweeps
  (EC1, EC3–EC8, EC10, EC11), slice 3 impact ionization and Fermi–Dirac, slice 4 transient, AC, deck, files, app
  (EC13–EC16), slice 5 FinFET, refinement, mutants, ARCHITECTURE (EC9, EC12, EC17, EC19).
